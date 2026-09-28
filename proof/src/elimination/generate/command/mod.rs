//! `rid generate <config.json>`: generates the configured covers, verifies
//! every generated file with the cover loader (pins, tiling, the lemma on
//! every cell) and writes the verified files.
use super::search::Settings;
use super::{generate, GenerateError, Generated, Request};
use crate::elimination::proof::catalogue::{self, LOCAL_COUNT};
use crate::elimination::proof::{Pin, ProvedCover, ProofError, Report};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fmt;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::time::Instant;

/// Exit status: every configured cover was generated, verified and written.
pub const SUCCESS: u8 = 0;
/// Exit status: a refusal or a failure (the reason is on standard error).
pub const FAILURE: u8 = 1;
/// Exit status: a malformed command line.
pub const USAGE_ERROR: u8 = 2;

const USAGE: &str = "usage: rid generate <config.json>";

/// The configuration of `rid generate`: every field required, unknown and
/// repeated fields refused, no defaults.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    /// The data directory: files go to `exotic/<name>.json` and
    /// `local/<number>.json` below it, relative to the working directory.
    pub output: PathBuf,
    /// Worker threads for generation and verification.
    pub threads: NonZeroUsize,
    pub search: Settings,
    /// Exotic covers to generate, by name.
    pub exotic: Vec<Exotic>,
    pub local: Local,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Exotic {
    pub name: String,
    /// Each variable of a zoom's base points takes `samples + 1` values.
    pub samples: NonZeroUsize,
    /// Domain inequalities offered as witnesses (tried first).
    pub inequalities: Vec<usize>,
}

/// Local covers to generate (scope `all`: no domain inequalities).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Local {
    pub numbers: Vec<usize>,
    pub samples: NonZeroUsize,
}

#[derive(Debug)]
pub enum Error {
    Usage(String),
    Read { file: PathBuf, source: io::Error },
    Json(serde_json::Error),
    /// An Exotic cover name that is not in the catalogue.
    UnknownCover(String),
    /// A Local cover number not below `LOCAL_COUNT`.
    UnknownLocal(usize),
    /// A cover listed twice.
    Repeated(String),
    Generate(GenerateError),
    /// A generated file that the cover loader refuses.
    Verify { cover: String, error: ProofError },
    /// The loader's counts differ from the generator's.
    Counts { cover: String },
    Write { file: PathBuf, source: io::Error },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Usage(message) => write!(f, "{message}\n{USAGE}"),
            Error::Read { file, source } => write!(f, "cannot read configuration {}: {source}", file.display()),
            Error::Json(source) => write!(f, "invalid generate configuration: {source}"),
            Error::UnknownCover(name) => write!(f, "no Exotic cover is named {name:?}"),
            Error::UnknownLocal(number) => write!(f, "Local cover {number} is not below {LOCAL_COUNT}"),
            Error::Repeated(name) => write!(f, "cover {name} is listed twice"),
            Error::Generate(error) => write!(f, "{error}"),
            Error::Verify { cover, error } => write!(f, "generated cover {cover} is refused: {error}"),
            Error::Counts { cover } => write!(f, "generated cover {cover}: the verified counts differ"),
            Error::Write { file, source } => write!(f, "cannot write {}: {source}", file.display()),
        }
    }
}

impl std::error::Error for Error {}

impl Configuration {
    pub fn parse(bytes: &[u8]) -> Result<Self, Error> {
        serde_json::from_slice(bytes).map_err(Error::Json)
    }

    /// The requests with their pins, and the path of each cover's file
    /// relative to `output`. Refuses, before any work, unknown and repeated
    /// covers, a sample grid above [`MAX_SAMPLES`](super::MAX_SAMPLES),
    /// unknown domain inequalities or inequalities offered to a cover of
    /// scope `all`, and a `max_depth` above the search's limit.
    pub fn requests(&self) -> Result<Vec<(Request, PathBuf)>, Error> {
        super::check_settings(&self.search).map_err(Error::Generate)?;
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        let mut add = |pin: Pin, samples: NonZeroUsize, inequalities: Vec<usize>, file: PathBuf| {
            if !seen.insert(file.clone()) {
                return Err(Error::Repeated(pin.name.clone()));
            }
            let request = Request { pin, samples, inequalities };
            super::check_request(&request).map_err(Error::Generate)?;
            out.push((request, file));
            Ok(())
        };
        for cover in &self.exotic {
            let pin = catalogue::exotic_pin(&cover.name).ok_or_else(|| Error::UnknownCover(cover.name.clone()))?;
            let file = Path::new("exotic").join(format!("{}.json", cover.name));
            add(pin, cover.samples, cover.inequalities.clone(), file)?;
        }
        for &number in &self.local.numbers {
            let pin = catalogue::local_pin(number).ok_or(Error::UnknownLocal(number))?;
            let file = Path::new("local").join(format!("{number}.json"));
            add(pin, self.local.samples, Vec::new(), file)?;
        }
        Ok(out)
    }
}

/// Writes `bytes` to `file` through a temporary file in the same directory,
/// so that the file is either the old one or the new one; on a failure the
/// temporary file is removed.
fn write_atomically(file: &Path, bytes: &[u8]) -> Result<(), Error> {
    let fail = |source| Error::Write { file: file.to_path_buf(), source };
    if let Some(directory) = file.parent() {
        std::fs::create_dir_all(directory).map_err(fail)?;
    }
    let temporary = file.with_extension("json.partial");
    let written = std::fs::write(&temporary, bytes).and_then(|()| std::fs::rename(&temporary, file));
    if written.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    written.map_err(fail)
}

/// The safety step of the command: verifies the generated cover with the
/// cover loader (pin, tiling, the lemma on every cell), requires the
/// loader's counts to equal the generator's, and only then writes `file`.
/// Returns the verification time and the file's size. A refused file is
/// never written; an interrupted verification is `Generate(Interrupted)`.
pub fn verify_and_write(
    request: &Request,
    generated: &Generated,
    file: &Path,
    threads: NonZeroUsize,
    stop: &AtomicBool,
) -> Result<(f64, usize), Error> {
    let name = || request.pin.name.clone();
    let bytes = generated.file.to_bytes();
    let start = Instant::now();
    let cover = ProvedCover::load(&request.pin, &bytes, threads, stop).map_err(|error| match error {
        ProofError::Interrupted => Error::Generate(GenerateError::Interrupted),
        error => Error::Verify { cover: name(), error },
    })?;
    if !agrees(cover.report(), &generated.statistics) {
        return Err(Error::Counts { cover: name() });
    }
    let seconds = start.elapsed().as_secs_f64();
    write_atomically(file, &bytes)?;
    Ok((seconds, bytes.len()))
}

/// Whether the loader's report agrees with the generator's counts.
fn agrees(report: &Report, statistics: &super::Statistics) -> bool {
    report.zooms == statistics.zooms
        && report.witnesses == statistics.witnesses
        && report.witness_leaves == statistics.witness_leaves
        && report.domain_leaves == statistics.domain_leaves
        && report.delegated_leaves == statistics.delegated_leaves
}

/// Runs the command on the configuration file; the summary goes to `out`.
/// Returns `Ok(true)` when every cover was written.
pub fn run(config: &Path, stop: &AtomicBool, out: &mut dyn Write) -> Result<bool, Error> {
    let bytes = std::fs::read(config).map_err(|source| Error::Read { file: config.to_path_buf(), source })?;
    let configuration = Configuration::parse(&bytes)?;
    let requests = configuration.requests()?;
    let start = Instant::now();
    let plain: Vec<Request> = requests.iter().map(|(r, _)| r.clone()).collect();
    let mut report = |p: &super::Progress| {
        eprintln!("rid generate: depth {}, {} open cells, {} decided, {:.0} s", p.depth, p.cells, p.decided, start.elapsed().as_secs_f64());
    };
    let results = generate(&plain, &configuration.search, configuration.threads, stop, &mut report).map_err(Error::Generate)?;
    let generated = start.elapsed().as_secs_f64();
    let mut summaries = Vec::new();
    let mut complete = true;
    for ((request, file), result) in requests.iter().zip(results) {
        let name = request.pin.name.clone();
        let outcome = result.map_err(Error::Generate).and_then(|g| {
            let (seconds, size) = verify_and_write(request, &g, &configuration.output.join(file), configuration.threads, stop)?;
            Ok((g, seconds, size))
        });
        if matches!(outcome, Err(Error::Generate(GenerateError::Interrupted))) {
            return Err(Error::Generate(GenerateError::Interrupted));
        }
        summaries.push(match outcome {
            Ok((g, verification, size)) => json!({
                "cover": name, "written": true, "statistics": g.statistics,
                "cpu_seconds": g.seconds, "verification_seconds": verification, "bytes": size,
            }),
            Err(error) => {
                complete = false;
                eprintln!("rid generate: {error}");
                json!({"cover": name, "written": false, "error": error.to_string()})
            }
        });
    }
    let summary: Value = json!({
        "command": "generate",
        "covers": summaries,
        "generation_seconds": generated,
        "seconds": start.elapsed().as_secs_f64(),
    });
    writeln!(out, "{summary}").map_err(|source| Error::Write { file: PathBuf::from("<stdout>"), source })?;
    Ok(complete)
}

/// `rid generate` with its arguments after the command name; returns the
/// exit status ([`SUCCESS`], [`FAILURE`] or [`USAGE_ERROR`]).
pub fn entry(args: &[OsString], stop: &AtomicBool, out: &mut dyn Write) -> u8 {
    let [config] = args else {
        eprintln!("rid: {}", Error::Usage(format!("expected one configuration file, got {}", args.len())));
        return USAGE_ERROR;
    };
    match run(Path::new(config), stop, out) {
        Ok(true) => SUCCESS,
        Ok(false) => FAILURE,
        Err(error) => {
            eprintln!("rid generate: {error}");
            FAILURE
        }
    }
}

#[cfg(test)]
mod tests;
