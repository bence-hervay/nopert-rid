//! The irredundancy experiment: at the configured scales around each
//! target, every cover and every per-box component is asked on its own
//! whether it eliminates the neighbourhood, and the verdict says whether
//! only the target's own component or cover does.
pub mod targets;
pub mod transcript;

use crate::config;
use crate::experiments::completeness::transcript::spell;
use crate::experiments::probe::{self, consistent, Contradiction, Probe, ProbeError, Probes};
use crate::experiments::run::{self, TranscriptError, Writer};
use crate::points::catalogue::CatalogueError;
use crate::points::neighbourhood::neighbourhood;
use crate::points::relation::{self, Relation};
use std::fmt;
use std::time::Instant;
use targets::{Target, Targets};
use transcript::{verdict, Header, Outcome, ScaleTrial, Summary, TargetResult};

#[derive(Debug)]
pub enum Error {
    Catalogue(CatalogueError),
    Transcript(TranscriptError),
    Invalid(String),
    /// The covers of a component are not records of it.
    Covers(ProbeError),
    Probe { id: String, k: u32, source: ProbeError },
    Contradiction { id: String, k: u32, record: String, source: Contradiction },
    Fit { id: String },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Catalogue(e) => e.fmt(f),
            Error::Transcript(e) => e.fmt(f),
            Error::Invalid(message) => write!(f, "invalid experiment: {message}"),
            Error::Covers(source) => write!(f, "covers: {source}"),
            Error::Probe { id, k, source } => write!(f, "{id} k={k}: {source}"),
            Error::Contradiction {
                id,
                k,
                record,
                source,
            } => write!(f, "{id} k={k}: {source} (record {record})"),
            Error::Fit { id } => write!(f, "{id} is a fit: the catalogue contradicts the theorem"),
        }
    }
}

impl std::error::Error for Error {}

impl From<TranscriptError> for Error {
    fn from(e: TranscriptError) -> Self {
        Error::Transcript(e)
    }
}

/// Every scale of one target.
pub fn evaluate<P: Probes>(
    probes: &P,
    names: &[String],
    list: &[Probe<P::Data>],
    scales: &[u32],
    target: &Target,
) -> Result<TargetResult, Error> {
    let start = Instant::now();
    let in_domain = relation::in_domain(&target.centre);
    let relation = relation::relation(&target.centre);
    if relation == Relation::Fit {
        return Err(Error::Fit { id: target.id.clone() });
    }
    let mut trials = Vec::with_capacity(scales.len());
    for &k in scales {
        let b = neighbourhood(&target.centre, k).expect("scales are valid");
        let mut outcomes = Vec::with_capacity(list.len());
        let mut successes = Vec::new();
        for p in list {
            let asked = Instant::now();
            let found = probe::ask(probes, names, p, &b).map_err(|source| Error::Probe {
                id: target.id.clone(),
                k,
                source,
            })?;
            if let Some(found) = &found {
                consistent(&found.label, in_domain, relation).map_err(|source| {
                    Error::Contradiction {
                        id: target.id.clone(),
                        k,
                        record: found.record.to_string(),
                        source,
                    }
                })?;
                successes.push(found.label.clone());
            }
            outcomes.push(Outcome {
                record: found.map(|f| f.record),
                milliseconds: asked.elapsed().as_millis() as u64,
            });
        }
        trials.push(ScaleTrial {
            k,
            axes: spell(&b),
            outcomes,
            verdict: verdict(&target.target, &successes),
        });
    }
    Ok(TargetResult {
        id: target.id.clone(),
        category: target.category,
        target: target.target.clone(),
        in_domain,
        relation,
        trials,
        milliseconds: start.elapsed().as_millis() as u64,
    })
}

/// Runs the experiment, one target per worker, and writes its transcript;
/// returns the summary.
pub fn run<P: Probes>(probes: &P, config: &config::Irredundancy) -> Result<Summary, Error> {
    if !transcript::valid_scales(&config.scales) {
        return Err(Error::Invalid(
            "scales must be nonempty, increasing and at most the largest exponent".into(),
        ));
    }
    let names = probe::names(probes).map_err(|e| Error::Invalid(e.to_string()))?;
    let list = probe::probes(probes).map_err(Error::Covers)?;
    let catalogue = Targets::load(&config.catalogue).map_err(Error::Catalogue)?;
    let selected = catalogue
        .select(config.selection.as_deref())
        .map_err(Error::Catalogue)?;
    let header = Header {
        format: transcript::FORMAT.into(),
        policy: rid::POLICY.into(),
        executable_sha256: run::executable_sha256()?,
        catalogue_sha256: catalogue.sha256.clone(),
        components: names.clone(),
        probes: list.iter().map(|p| p.label.clone()).collect(),
        scales: config.scales.clone(),
        selection: config.selection.clone(),
    };
    transcript::check_header(&header, &catalogue.sha256)
        .map_err(|e| Error::Invalid(e.to_string()))?;
    let mut writer = Writer::create(&config.transcript)?;
    writer.write(&header)?;
    let mut summary = Summary::new();
    let total = selected.len();
    run::ordered(
        &selected,
        config.threads,
        |_, target| evaluate(probes, &names, &list, &config.scales, target),
        |index, result| {
            summary.add(&result);
            writer.write(&result)?;
            let verdicts: Vec<String> = result
                .trials
                .iter()
                .map(|t| format!("{}:{:?}", t.k, t.verdict))
                .collect();
            eprintln!(
                "irredundancy: {}/{total} {} ({:.1} s) {}",
                index + 1,
                result.id,
                result.milliseconds as f64 / 1000.0,
                verdicts.join(" ")
            );
            Ok(())
        },
    )?;
    writer.write(&summary)?;
    Ok(summary)
}

#[cfg(test)]
mod tests;
