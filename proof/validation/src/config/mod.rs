//! Typed JSON configurations of the commands: one struct per command, every
//! field required, unknown and repeated fields refused, no defaults; `null`
//! only where a selection or limit is explicitly absent.
use serde::{Deserialize, Deserializer, Serialize};
use std::fmt;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// `completeness`: radii around every selected catalogue point.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Completeness {
    /// The point catalogue.
    pub catalogue: PathBuf,
    /// The transcript to create; an existing file is refused.
    pub transcript: PathBuf,
    /// Points evaluated in parallel.
    pub threads: NonZeroUsize,
    /// The finest radius exponent, at most 2048.
    pub max_k: u32,
    /// The identifiers to run, or `null` for every point.
    #[serde(deserialize_with = "required")]
    pub selection: Option<Vec<String>>,
}

/// `irredundancy`: every probe around every selected target at fixed scales.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Irredundancy {
    /// The target catalogue.
    pub catalogue: PathBuf,
    pub transcript: PathBuf,
    /// Targets evaluated in parallel.
    pub threads: NonZeroUsize,
    /// The radius exponents, increasing.
    pub scales: Vec<u32>,
    #[serde(deserialize_with = "required")]
    pub selection: Option<Vec<String>>,
}

/// Roots drawn uniformly without replacement at one depth.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub depth: u32,
    pub count: u64,
    pub seed: u64,
}

/// The roots of a pilot: `{"paths": [...]}` or `{"sample": {...}}`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase", deny_unknown_fields)]
pub enum Roots {
    Paths(Vec<String>),
    Sample(Sample),
}

/// `pilot`: bounded searches below each root.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pilot {
    pub roots: Roots,
    /// Roots searched in parallel, each on one thread.
    pub threads: NonZeroUsize,
    /// Boxes at this depth are evaluated but not split.
    pub depth_limit: usize,
    /// Committed decisions per root after which its search stops, or `null`.
    #[serde(deserialize_with = "required")]
    pub max_decisions: Option<u64>,
    pub transcript: PathBuf,
}

/// `targets`: regenerates, in place, the targets of the Local covers in a
/// target catalogue from the crate's covers.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Targets {
    /// The target catalogue, rewritten in place.
    pub catalogue: PathBuf,
}

/// `check`: the structure and enclosures of one transcript.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Check {
    /// The catalogue the transcript was made from.
    pub catalogue: PathBuf,
    pub transcript: PathBuf,
    /// Entries checked in parallel.
    pub threads: NonZeroUsize,
}

/// `compare`: two completeness transcripts of one catalogue.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Compare {
    pub catalogue: PathBuf,
    pub baseline: PathBuf,
    pub candidate: PathBuf,
    /// Points checked in parallel.
    pub threads: NonZeroUsize,
    /// The least increase of a best exponent that counts as a regression.
    pub worse_by: u32,
    /// Whether a regression makes the command fail with status 3.
    pub fail_on_regression: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Command {
    Completeness(Completeness),
    Irredundancy(Irredundancy),
    Pilot(Pilot),
    Targets(Targets),
    Check(Check),
    Compare(Compare),
}

pub const COMMANDS: [&str; 6] = ["completeness", "irredundancy", "pilot", "targets", "check", "compare"];

#[derive(Debug)]
pub enum Error {
    UnknownCommand(String),
    Read { file: PathBuf, source: std::io::Error },
    Json { command: &'static str, source: serde_json::Error },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::UnknownCommand(name) => {
                write!(f, "unknown command {name:?}; expected one of {}", COMMANDS.join(", "))
            }
            Error::Read { file, source } => {
                write!(f, "cannot read configuration {}: {source}", file.display())
            }
            Error::Json { command, source } => {
                write!(f, "invalid {command} configuration: {source}")
            }
        }
    }
}

impl std::error::Error for Error {}

/// Makes an `Option` field required: only an explicit `null` gives `None`.
pub(crate) fn required<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(deserializer)
}

fn json<'de, T: Deserialize<'de>>(command: &'static str, text: &'de [u8]) -> Result<T, Error> {
    serde_json::from_slice(text).map_err(|source| Error::Json { command, source })
}

impl Command {
    pub fn parse(name: &str, text: &[u8]) -> Result<Self, Error> {
        match name {
            "completeness" => json("completeness", text).map(Command::Completeness),
            "irredundancy" => json("irredundancy", text).map(Command::Irredundancy),
            "pilot" => json("pilot", text).map(Command::Pilot),
            "targets" => json("targets", text).map(Command::Targets),
            "check" => json("check", text).map(Command::Check),
            "compare" => json("compare", text).map(Command::Compare),
            _ => Err(Error::UnknownCommand(name.to_owned())),
        }
    }

    /// Refuses an unknown command before reading the file.
    pub fn load(name: &str, file: &Path) -> Result<Self, Error> {
        if !COMMANDS.contains(&name) {
            return Err(Error::UnknownCommand(name.to_owned()));
        }
        let text = std::fs::read(file).map_err(|source| Error::Read {
            file: file.to_owned(),
            source,
        })?;
        Self::parse(name, &text)
    }
}

#[cfg(test)]
mod tests;
