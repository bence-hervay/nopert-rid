//! Pilot searches for estimating a full run: bounded searches below a list
//! of roots, one root per worker, each with the crate's own queue and the
//! collection's verified decisions, reporting records, depths and times.
//!
use crate::config::{self, Roots, Sample};
use crate::experiments::probe::{self, Label, Probes};
use crate::experiments::run::{self, TranscriptError, Writer};
use rid::problem::configuration::{is_path, ConfigurationBox};
use rid::search::certificate::MAX_DEPTH;
use rid::search::queue::{self, Halt, Output};
use rid::search::BoxError;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::AtomicBool;
use std::sync::Mutex;
use std::time::Instant;

pub const FORMAT: &str = "rid-pilot/1";

/// The largest sampling depth: the population `2^depth` fits a `u64`.
pub const MAX_SAMPLE_DEPTH: u32 = 62;

/// The largest sample: drawing `count` distinct paths by rejection then
/// costs at most about `count · ln count` draws.
pub const MAX_SAMPLE_COUNT: u64 = 1 << 16;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    pub format: String,
    pub policy: String,
    pub executable_sha256: String,
    pub components: Vec<String>,
    /// The roots searched, in search order.
    pub roots: Vec<String>,
    /// How the roots were drawn, or `null` for a given list.
    pub sample: Option<Sample>,
    pub depth_limit: usize,
    pub max_decisions: Option<u64>,
}

/// Evaluations with one outcome and the time they took. Times are whole
/// microseconds, so transcript lines have one exact spelling.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Timing {
    pub boxes: u64,
    pub microseconds: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RootResult {
    pub root: String,
    /// `complete`, `decision-limit` or `depth-limit`.
    pub halt: String,
    pub evaluated: u64,
    pub decisions: u64,
    pub records: u64,
    pub splits: u64,
    pub unresolved: u64,
    /// The depth of the deepest record.
    pub deepest: Option<usize>,
    /// Records per label.
    pub labels: BTreeMap<String, u64>,
    /// Records per depth.
    pub depths: BTreeMap<usize, u64>,
    /// Evaluations per outcome (a label, or `refused`) with their time.
    pub outcomes: BTreeMap<String, Timing>,
    /// Wall-clock time of the root's search on its worker.
    pub microseconds: u64,
}

/// A total over the population estimated from a uniform sample without
/// replacement: `population · mean`, with its standard error
/// `population · sd / √n · √(1 - n / population)` (none for one root).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Extrapolation {
    pub mean: f64,
    pub total: f64,
    pub standard_error: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Estimate {
    /// The number of boxes at the sampling depth, `2^depth`.
    pub population: f64,
    pub records: Extrapolation,
    pub evaluated: Extrapolation,
    /// Search time on one thread, in seconds.
    pub seconds: Extrapolation,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub roots: usize,
    /// Roots whose search stopped at a limit: their counts are lower bounds.
    pub incomplete: usize,
    pub records: u64,
    pub evaluated: u64,
    /// The roots' search times added up.
    pub microseconds: u64,
    pub labels: BTreeMap<String, u64>,
    /// For sampled roots only.
    pub estimate: Option<Estimate>,
}

#[derive(Debug)]
pub enum Error {
    Invalid(String),
    Transcript(TranscriptError),
    Search { root: String, source: queue::Error },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Invalid(message) => write!(f, "invalid pilot: {message}"),
            Error::Transcript(e) => e.fmt(f),
            Error::Search { root, source } => write!(f, "root {root:?}: {source}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<TranscriptError> for Error {
    fn from(e: TranscriptError) -> Self {
        Error::Transcript(e)
    }
}

/// The SplitMix64 generator: a fixed, documented sequence for each seed.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// `count` distinct paths of length `depth`, drawn uniformly without
/// replacement from the seeded generator, in search order.
pub fn sample(s: &Sample) -> Result<Vec<String>, Error> {
    if s.depth > MAX_SAMPLE_DEPTH || s.count == 0 || s.count > 1u64 << s.depth || s.count > MAX_SAMPLE_COUNT {
        return Err(Error::Invalid(format!(
            "a sample needs depth at most {MAX_SAMPLE_DEPTH} and 1 ≤ count ≤ min(2^depth, {MAX_SAMPLE_COUNT})"
        )));
    }
    let mut state = s.seed;
    let mut drawn = BTreeSet::new();
    while (drawn.len() as u64) < s.count {
        let bits = if s.depth == 0 { 0 } else { splitmix64(&mut state) >> (64 - s.depth) };
        drawn.insert(bits);
    }
    let width = s.depth as usize;
    Ok(drawn
        .into_iter()
        .map(|bits| if width == 0 { String::new() } else { format!("{bits:0width$b}") })
        .collect())
}

/// The roots of a configuration: distinct binary paths no deeper than the
/// depth limit, none below another (whose subtree would be searched and
/// counted twice), in search order.
pub fn roots(config: &config::Pilot) -> Result<Vec<String>, Error> {
    let mut roots = match &config.roots {
        Roots::Paths(paths) => paths.clone(),
        Roots::Sample(s) => sample(s)?,
    };
    let distinct: BTreeSet<&String> = roots.iter().collect();
    if roots.is_empty() || distinct.len() != roots.len() {
        return Err(Error::Invalid("roots must be nonempty and distinct".into()));
    }
    // In lexicographic order a path is followed directly by its extensions.
    let ordered: Vec<&&String> = distinct.iter().collect();
    if ordered.windows(2).any(|pair| pair[1].starts_with(pair[0].as_str())) {
        return Err(Error::Invalid("no root may lie below another root".into()));
    }
    if config.depth_limit > MAX_DEPTH {
        return Err(Error::Invalid(format!("depth_limit exceeds {MAX_DEPTH}")));
    }
    if roots.iter().any(|r| !is_path(r) || r.len() > config.depth_limit) {
        return Err(Error::Invalid("roots must be binary paths within the depth limit".into()));
    }
    roots.sort_by(|a, b| (a.len(), a).cmp(&(b.len(), b)));
    Ok(roots)
}

/// Counts the records the queue commits.
#[derive(Default)]
struct Tally {
    labels: BTreeMap<String, u64>,
    depths: BTreeMap<usize, u64>,
}

impl Output<Label> for Tally {
    fn append(&mut self, path: &str, label: &Label) -> Result<(), BoxError> {
        *self.labels.entry(label.to_string()).or_default() += 1;
        *self.depths.entry(path.len()).or_default() += 1;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), BoxError> {
        Ok(())
    }
}

/// A bounded search below one root on the calling thread: the queue with one
/// worker, the collection's decision per box verified by its component.
pub fn search<P: Probes>(
    probes: &P,
    names: &[String],
    root: &str,
    depth_limit: usize,
    max_decisions: Option<u64>,
) -> Result<RootResult, Error> {
    let one = NonZeroUsize::new(1).expect("one");
    let options = queue::Options {
        threads: one,
        window: one,
        depth_limit,
        max_decisions,
        after: None,
    };
    let timings: Mutex<BTreeMap<String, Timing>> = Mutex::new(BTreeMap::new());
    let evaluate = |path: &str| -> Result<Option<Label>, BoxError> {
        let b = ConfigurationBox::from_path(path)?;
        // Timed: what the search itself does per box, the collection's check.
        let start = Instant::now();
        let saved = probes.decide(&b)?;
        let elapsed = start.elapsed();
        // Not timed: the experiment's own check of the saved record.
        let label = match saved {
            Some(data) => Some(probe::confirm(probes, names, &b, &data)?.label),
            None => None,
        };
        let outcome = label.as_ref().map_or("refused".into(), Label::to_string);
        let mut timings = timings.lock().unwrap_or_else(|e| e.into_inner());
        let timing = timings.entry(outcome).or_insert(Timing {
            boxes: 0,
            microseconds: 0,
        });
        timing.boxes += 1;
        timing.microseconds += elapsed.as_micros() as u64;
        Ok(label)
    };
    let mut tally = Tally::default();
    let never = AtomicBool::new(false);
    let start = Instant::now();
    let report = queue::run(vec![root.to_owned()], &options, &never, evaluate, &mut tally)
        .map_err(|source| Error::Search {
            root: root.to_owned(),
            source,
        })?;
    let halt = match report.halt {
        Halt::Complete => "complete",
        Halt::Stopped => "stopped",
        Halt::DecisionLimit => "decision-limit",
        Halt::DepthLimit => "depth-limit",
    };
    Ok(RootResult {
        root: root.to_owned(),
        halt: halt.into(),
        evaluated: report.evaluated,
        decisions: report.decisions,
        records: report.records,
        splits: report.splits,
        unresolved: report.unresolved,
        deepest: tally.depths.keys().next_back().copied(),
        labels: tally.labels,
        depths: tally.depths,
        outcomes: timings.into_inner().unwrap_or_else(|e| e.into_inner()),
        microseconds: start.elapsed().as_micros() as u64,
    })
}

/// See [`Extrapolation`].
pub fn extrapolate(values: &[f64], population: f64) -> Extrapolation {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let standard_error = (values.len() > 1).then(|| {
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        population * (variance / n * (1.0 - n / population)).max(0.0).sqrt()
    });
    Extrapolation {
        mean,
        total: population * mean,
        standard_error,
    }
}

fn summarise(results: &[RootResult], sample: Option<&Sample>) -> Summary {
    let mut labels = BTreeMap::new();
    for r in results {
        for (label, count) in &r.labels {
            *labels.entry(label.clone()).or_default() += count;
        }
    }
    let estimate = sample.map(|s| {
        let population = (1u64 << s.depth) as f64;
        let column = |f: fn(&RootResult) -> f64| -> Vec<f64> { results.iter().map(f).collect() };
        Estimate {
            population,
            records: extrapolate(&column(|r| r.records as f64), population),
            evaluated: extrapolate(&column(|r| r.evaluated as f64), population),
            seconds: extrapolate(&column(|r| r.microseconds as f64 / 1e6), population),
        }
    });
    Summary {
        roots: results.len(),
        incomplete: results.iter().filter(|r| r.halt != "complete").count(),
        records: results.iter().map(|r| r.records).sum(),
        evaluated: results.iter().map(|r| r.evaluated).sum(),
        microseconds: results.iter().map(|r| r.microseconds).sum(),
        labels,
        estimate,
    }
}

/// Runs every root's search, one root per worker, writes the transcript and
/// returns the summary.
pub fn run<P: Probes>(probes: &P, config: &config::Pilot) -> Result<Summary, Error> {
    let names = probe::names(probes).map_err(|e| Error::Invalid(e.to_string()))?;
    let roots = roots(config)?;
    let sample = match &config.roots {
        Roots::Sample(s) => Some(s.clone()),
        Roots::Paths(_) => None,
    };
    let header = Header {
        format: FORMAT.into(),
        policy: rid::POLICY.into(),
        executable_sha256: run::executable_sha256()?,
        components: names.clone(),
        roots: roots.clone(),
        sample: sample.clone(),
        depth_limit: config.depth_limit,
        max_decisions: config.max_decisions,
    };
    let mut writer = Writer::create(&config.transcript)?;
    writer.write(&header)?;
    let mut results = Vec::with_capacity(roots.len());
    let total = roots.len();
    run::ordered(
        &roots,
        config.threads,
        |_, root| search(probes, &names, root, config.depth_limit, config.max_decisions),
        |index, result| {
            writer.write(&result)?;
            eprintln!(
                "pilot: {}/{total} root {:?}: {} records, {} boxes, deepest {:?}, {} ({:.1} s)",
                index + 1,
                result.root,
                result.records,
                result.evaluated,
                result.deepest,
                result.halt,
                result.microseconds as f64 / 1e6
            );
            results.push(result);
            Ok(())
        },
    )?;
    let summary = summarise(&results, sample.as_ref());
    writer.write(&summary)?;
    Ok(summary)
}

#[cfg(test)]
mod tests;
