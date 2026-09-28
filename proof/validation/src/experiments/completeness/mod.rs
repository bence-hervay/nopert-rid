//! The completeness experiment: for every catalogue point, the largest
//! tested radius `2^-k` at which the collection eliminates the neighbourhood
//! of the point, the record it uses, the refusals at finer tested radii, and
//! the same for each component on its own.
pub mod compare;
pub mod schedule;
pub mod transcript;

use crate::config;
use crate::experiments::probe::{self, consistent, Contradiction, Found, ProbeError, Probes};
use crate::experiments::run::{self, TranscriptError, Writer};
use crate::points::catalogue::{Catalogue, CatalogueError, Point};
use crate::points::neighbourhood::neighbourhood;
use crate::points::relation::{self, Relation};
use rid::problem::configuration::ConfigurationBox;
use schedule::MAX_K;
use std::collections::BTreeMap;
use std::fmt;
use std::time::Instant;
use transcript::{
    disagreement, spell, Header, Neighbourhood, PointResult, Schedule, Summary, Trial,
};

#[derive(Debug)]
pub enum Error {
    Catalogue(CatalogueError),
    Transcript(TranscriptError),
    /// A configuration or component list the experiment cannot run with.
    Invalid(String),
    Probe { id: String, k: u32, source: ProbeError },
    Contradiction { id: String, k: u32, record: String, source: Contradiction },
    /// A catalogue point that is a fit.
    Fit { id: String },
    /// The decision at radius `2^-k` is not the first component success.
    Disagreement { id: String, k: u32 },
    /// The best record of a schedule is refused by its component on the
    /// nested finer neighbourhood `2^-k`: verification is not inherited by
    /// sub-boxes, which the components' proofs say it must be.
    NotMonotone { id: String, k: u32, record: String, reason: String },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Catalogue(e) => e.fmt(f),
            Error::Transcript(e) => e.fmt(f),
            Error::Invalid(message) => write!(f, "invalid experiment: {message}"),
            Error::Probe { id, k, source } => write!(f, "{id} k={k}: {source}"),
            Error::Contradiction {
                id,
                k,
                record,
                source,
            } => write!(f, "{id} k={k}: {source} (record {record})"),
            Error::Fit { id } => write!(f, "{id} is a fit: the catalogue contradicts the theorem"),
            Error::Disagreement { id, k } => write!(
                f,
                "{id} k={k}: the collection's decision is not the first component success"
            ),
            Error::NotMonotone { id, k, record, reason } => write!(
                f,
                "{id} k={k}: the best record {record} is refused on the finer neighbourhood: {reason}"
            ),
        }
    }
}

impl std::error::Error for Error {}

impl From<TranscriptError> for Error {
    fn from(e: TranscriptError) -> Self {
        Error::Transcript(e)
    }
}

/// Runs one schedule, asking `ask` at every radius it tests. Neighbourhoods
/// are built once per exponent and shared between the schedules of a point.
fn run_schedule(
    point: &Point,
    max_k: u32,
    in_domain: bool,
    relation: Relation,
    boxes: &mut BTreeMap<u32, ConfigurationBox>,
    mut ask: impl FnMut(&ConfigurationBox) -> Result<Option<Found>, ProbeError>,
) -> Result<Schedule, Error> {
    let mut trials = Vec::new();
    let mut outcomes = Vec::new();
    while let Some(k) = schedule::next(max_k, &outcomes).expect("the runner follows the schedule") {
        let b = boxes
            .entry(k)
            .or_insert_with(|| neighbourhood(&point.centre, k).expect("k is at most MAX_K"));
        let start = Instant::now();
        let found = ask(b).map_err(|source| Error::Probe {
            id: point.id.clone(),
            k,
            source,
        })?;
        if let Some(found) = &found {
            consistent(&found.label, in_domain, relation).map_err(|source| {
                Error::Contradiction {
                    id: point.id.clone(),
                    k,
                    record: found.record.to_string(),
                    source,
                }
            })?;
        }
        outcomes.push((k, found.is_some()));
        trials.push(Trial {
            k,
            record: found.map(|f| f.record),
            milliseconds: start.elapsed().as_millis() as u64,
        });
    }
    Ok(Schedule::of(trials))
}

/// Where a schedule was refused at finer radii than its best one, its best
/// record must still verify there: the neighbourhoods are nested and every
/// component's verification is inherited by sub-boxes. A refusal there is
/// then a missed witness of the attempt; a failing verification is a defect.
fn check_inherited<P: Probes>(
    probes: &P,
    names: &[String],
    point: &Point,
    schedule: &Schedule,
    boxes: &BTreeMap<u32, ConfigurationBox>,
) -> Result<(), Error> {
    let (Some(_), Some(record)) = (schedule.best_k, &schedule.record) else {
        return Ok(());
    };
    let data: P::Data = serde_json::from_value(record.clone()).expect("a record the probes wrote");
    let label = probe::Label::of(record).expect("a record the probes wrote");
    let index = names.iter().position(|n| *n == label.component).expect("a checked component");
    for &k in &schedule.finer_refusals {
        probes.verify(index, &boxes[&k], &data).map_err(|reason| Error::NotMonotone {
            id: point.id.clone(),
            k,
            record: record.to_string(),
            reason: reason.to_string(),
        })?;
    }
    Ok(())
}

/// All schedules of one point.
pub fn evaluate<P: Probes>(
    probes: &P,
    names: &[String],
    max_k: u32,
    point: &Point,
) -> Result<PointResult, Error> {
    let start = Instant::now();
    let in_domain = relation::in_domain(&point.centre);
    let relation = relation::relation(&point.centre);
    if relation == Relation::Fit {
        return Err(Error::Fit { id: point.id.clone() });
    }
    let mut boxes = BTreeMap::new();
    let decision = run_schedule(point, max_k, in_domain, relation, &mut boxes, |b| {
        probe::decide(probes, names, b)
    })?;
    let mut components = Vec::with_capacity(names.len());
    for index in 0..names.len() {
        components.push(run_schedule(point, max_k, in_domain, relation, &mut boxes, |b| {
            probe::attempt(probes, names, index, b)
        })?);
    }
    if let Some(k) = disagreement(names, &decision, &components) {
        return Err(Error::Disagreement { id: point.id.clone(), k });
    }
    for schedule in std::iter::once(&decision).chain(&components) {
        check_inherited(probes, names, point, schedule, &boxes)?;
    }
    Ok(PointResult {
        id: point.id.clone(),
        group: point.group.clone(),
        in_domain,
        relation,
        expected: transcript::is_expected(&decision.record, &point.expected),
        neighbourhoods: boxes
            .iter()
            .map(|(&k, b)| Neighbourhood { k, axes: spell(b) })
            .collect(),
        decision,
        components,
        milliseconds: start.elapsed().as_millis() as u64,
    })
}

/// Runs the experiment and writes its transcript; returns the summary.
pub fn run<P: Probes>(probes: &P, config: &config::Completeness) -> Result<Summary, Error> {
    if !(1..=MAX_K).contains(&config.max_k) {
        return Err(Error::Invalid(format!("max_k must lie in 1..={MAX_K}")));
    }
    let names = probe::names(probes).map_err(|e| Error::Invalid(e.to_string()))?;
    let catalogue = Catalogue::load(&config.catalogue).map_err(Error::Catalogue)?;
    let points = catalogue
        .select(config.selection.as_deref())
        .map_err(Error::Catalogue)?;
    let header = Header {
        format: transcript::FORMAT.into(),
        policy: rid::POLICY.into(),
        executable_sha256: run::executable_sha256()?,
        catalogue_sha256: catalogue.sha256.clone(),
        components: names.clone(),
        max_k: config.max_k,
        selection: config.selection.clone(),
    };
    let mut writer = Writer::create(&config.transcript)?;
    writer.write(&header)?;
    let mut summary = Summary::new(&names);
    let total = points.len();
    run::ordered(
        &points,
        config.threads,
        |_, point| evaluate(probes, &names, config.max_k, point),
        |index, result| {
            summary.add(&result, &names);
            writer.write(&result)?;
            let reason = result.decision.record.as_ref().map_or("unresolved".into(), |r| {
                probe::Label::of(r).map_or_else(|_| "?".into(), |l| l.to_string())
            });
            eprintln!(
                "completeness: {}/{total} {}: k={} {reason} ({:.1} s)",
                index + 1,
                result.id,
                result.decision.best_k.map_or("-".into(), |k| k.to_string()),
                result.milliseconds as f64 / 1000.0
            );
            Ok(())
        },
    )?;
    writer.write(&summary)?;
    Ok(summary)
}

#[cfg(test)]
mod tests;
