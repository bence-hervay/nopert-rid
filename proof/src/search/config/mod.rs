//! Typed JSON run configurations: one struct per command, every field
//! required, unknown fields refused, no defaults. `null` is accepted only
//! where a limit may be absent.

use serde::{Deserialize, Deserializer, Serialize};
use std::fmt;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// `rid prepare`: prepare the mathematical components and report on them.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Prepare {
    /// Worker threads for the preparation.
    pub threads: NonZeroUsize,
}

/// `rid search`: create or resume a certificate and extend it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Search {
    /// The certificate file, relative to the working directory.
    pub certificate: PathBuf,
    /// Binary path of the box the certificate covers; `""` is the root box.
    pub root: String,
    /// Worker threads for verification and search.
    pub threads: NonZeroUsize,
    /// Boxes assigned to workers but not yet committed.
    pub window: NonZeroUsize,
    /// The certificate's depth limit, written to its header: no record is
    /// ever deeper. Part of the certificate's identity.
    pub max_depth: usize,
    /// Boxes at this depth are evaluated but not split in this run, or
    /// `null` for `max_depth`. At most `max_depth` and at least the depth of
    /// `root` (checked by the search command). A later run with a deeper
    /// limit continues the same certificate in place.
    #[serde(deserialize_with = "required")]
    pub depth_limit: Option<usize>,
    /// Committed decisions after which this run stops, or `null` for none.
    #[serde(deserialize_with = "required")]
    pub max_decisions: Option<u64>,
}

/// `rid check`: verify a certificate without changing it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Check {
    /// The certificate file, relative to the working directory.
    pub certificate: PathBuf,
    /// The root the certificate must cover.
    pub root: String,
    /// The depth limit the certificate's header must state.
    pub max_depth: usize,
    /// Worker threads for verification.
    pub threads: NonZeroUsize,
}

/// A command with its configuration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Command {
    Prepare(Prepare),
    Search(Search),
    Check(Check),
}

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
                write!(f, "unknown command {name:?}; expected prepare, search or check")
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

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::UnknownCommand(_) => None,
            Error::Read { source, .. } => Some(source),
            Error::Json { source, .. } => Some(source),
        }
    }
}

/// Makes an `Option` field required: serde would otherwise read a missing
/// field as `None`. Only an explicit `null` gives `None`.
fn required<'de, D, T>(deserializer: D) -> Result<T, D::Error>
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
    /// Parses the configuration of the command called `name`.
    pub fn parse(name: &str, text: &[u8]) -> Result<Self, Error> {
        match name {
            "prepare" => json("prepare", text).map(Command::Prepare),
            "search" => json("search", text).map(Command::Search),
            "check" => json("check", text).map(Command::Check),
            _ => Err(Error::UnknownCommand(name.to_owned())),
        }
    }

    /// Reads and parses the configuration file of the command called `name`.
    pub fn load(name: &str, file: &Path) -> Result<Self, Error> {
        if !matches!(name, "prepare" | "search" | "check") {
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
