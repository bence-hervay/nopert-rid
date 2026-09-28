//! The experiments and the checks of their transcripts.
pub mod probe;
pub mod run;
pub mod completeness;
pub mod irredundancy;
pub mod pilot;

use crate::config;
use crate::points::catalogue::{Catalogue, CatalogueError};
use irredundancy::targets::Targets;
use run::TranscriptError;
use serde::Deserialize;
use serde_json::{json, Value};
use std::fmt;

#[derive(Debug)]
pub enum CheckError {
    Transcript(TranscriptError),
    Catalogue(CatalogueError),
    /// The first line names no known transcript format.
    UnknownFormat,
    Completeness(completeness::transcript::CheckError),
    Irredundancy(irredundancy::transcript::CheckError),
    /// Two completeness transcripts that cannot be compared.
    Incompatible(completeness::compare::Incompatible),
}

impl fmt::Display for CheckError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CheckError::Transcript(e) => e.fmt(f),
            CheckError::Catalogue(e) => e.fmt(f),
            CheckError::UnknownFormat => write!(f, "the transcript has no known format"),
            CheckError::Completeness(e) => e.fmt(f),
            CheckError::Irredundancy(e) => e.fmt(f),
            CheckError::Incompatible(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for CheckError {}

/// The `format` field of a transcript's first line.
#[derive(Deserialize)]
struct Format {
    format: String,
}

/// Checks a completeness or irredundancy transcript, whichever its header
/// names, against the catalogue it was made from, with every record of the
/// form `form`; returns a short report. The eliminations themselves are not
/// replayed (`"replayed": false`): every success was verified by its own
/// component when the experiment ran.
pub fn check(config: &config::Check, form: probe::RecordForm) -> Result<Value, CheckError> {
    let bytes = run::read(&config.transcript).map_err(CheckError::Transcript)?;
    let first = bytes.split(|&b| b == b'\n').next().unwrap_or_default();
    let format: Format = serde_json::from_slice(first).map_err(|_| CheckError::UnknownFormat)?;
    match format.format.as_str() {
        completeness::transcript::FORMAT => {
            let catalogue = Catalogue::load(&config.catalogue).map_err(CheckError::Catalogue)?;
            let checked = completeness::transcript::check(&catalogue, &bytes, config.threads, form)
                .map_err(CheckError::Completeness)?;
            let trials: usize = checked
                .points
                .iter()
                .flat_map(|p| std::iter::once(&p.decision).chain(&p.components))
                .map(|schedule| schedule.trials.len())
                .sum();
            Ok(json!({
                "status": "pass",
                "format": format.format,
                "policy": checked.header.policy,
                "replayed": false,
                "trials": trials,
                "summary": checked.summary,
            }))
        }
        irredundancy::transcript::FORMAT => {
            let catalogue = Targets::load(&config.catalogue).map_err(CheckError::Catalogue)?;
            let checked = irredundancy::transcript::check(&catalogue, &bytes, config.threads, form)
                .map_err(CheckError::Irredundancy)?;
            Ok(json!({
                "status": "pass",
                "format": format.format,
                "policy": checked.header.policy,
                "replayed": false,
                "summary": checked.summary,
            }))
        }
        _ => Err(CheckError::UnknownFormat),
    }
}

/// Checks both transcripts against the catalogue (records of the form
/// `form`) and compares them.
pub fn compare(config: &config::Compare, form: probe::RecordForm) -> Result<completeness::compare::Comparison, CheckError> {
    let catalogue = Catalogue::load(&config.catalogue).map_err(CheckError::Catalogue)?;
    let load = |path| -> Result<completeness::transcript::Checked, CheckError> {
        let bytes = run::read(path).map_err(CheckError::Transcript)?;
        completeness::transcript::check(&catalogue, &bytes, config.threads, form)
            .map_err(CheckError::Completeness)
    };
    let (baseline, candidate) = (load(&config.baseline)?, load(&config.candidate)?);
    completeness::compare::compare(&baseline, &candidate, config.worse_by).map_err(CheckError::Incompatible)
}
