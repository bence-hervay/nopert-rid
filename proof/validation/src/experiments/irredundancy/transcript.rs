//! The irredundancy transcript: its typed lines, the verdict rule and the
//! check of a whole transcript against its target catalogue.
use super::targets::{Category, Target, Targets};
use crate::experiments::completeness::transcript::{unspell, Axes};
use crate::experiments::probe::{checked_label, consistent, Contradiction, Label, ProbeError, RecordForm};
use crate::experiments::run::{self, TranscriptError};
use crate::points::catalogue::CatalogueError;
use crate::points::neighbourhood::{self, NeighbourhoodError, MAX_EXPONENT};
use crate::points::relation::{self, Relation};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::num::NonZeroUsize;

pub const FORMAT: &str = "rid-irredundancy/1";

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    pub format: String,
    pub policy: String,
    pub executable_sha256: String,
    pub catalogue_sha256: String,
    /// The component names in the collection's order.
    pub components: Vec<String>,
    /// Every probe in the order asked: a cover (component and cover) or
    /// a component's attempt (component alone).
    pub probes: Vec<Label>,
    /// The radius exponents, increasing.
    pub scales: Vec<u32>,
    #[serde(deserialize_with = "crate::config::required")]
    pub selection: Option<Vec<String>>,
}

/// What the probes found at one radius, relative to the target.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Verdict {
    /// Only the target succeeded.
    Exclusive,
    /// The target and something else succeeded.
    Shared,
    /// Something else succeeded, the target did not.
    Other,
    /// Nothing succeeded.
    Unresolved,
}

/// The verdict for the labels of the successful probes.
pub fn verdict(target: &Label, successes: &[Label]) -> Verdict {
    let hit = successes.iter().any(|l| l.meets(target));
    let other = successes.iter().any(|l| !l.meets(target));
    match (hit, other) {
        (true, false) => Verdict::Exclusive,
        (true, true) => Verdict::Shared,
        (false, true) => Verdict::Other,
        (false, false) => Verdict::Unresolved,
    }
}

/// One probe at one radius: the verified record, or `null`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Outcome {
    #[serde(deserialize_with = "crate::config::required")]
    pub record: Option<Value>,
    pub milliseconds: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScaleTrial {
    pub k: u32,
    pub axes: Axes,
    /// One outcome per probe of the header, in its order.
    pub outcomes: Vec<Outcome>,
    pub verdict: Verdict,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetResult {
    pub id: String,
    pub category: Category,
    pub target: Label,
    pub in_domain: bool,
    pub relation: Relation,
    pub trials: Vec<ScaleTrial>,
    pub milliseconds: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub targets: usize,
    pub trials: usize,
    /// Verdict counts per category (`primary`, `control`).
    pub verdicts: BTreeMap<String, BTreeMap<String, usize>>,
    /// Primary targets exclusive at one scale at least, with those scales
    /// (`k` of the radius `2^-k`): their cover or component is needed.
    pub needed: BTreeMap<String, Vec<u32>>,
    /// Primary targets that succeed at some scale but never alone: wherever
    /// they succeed, another probe does too.
    pub covered: Vec<String>,
    /// Primary targets that never succeed.
    pub never: Vec<String>,
}

fn name<T: Serialize>(value: &T) -> String {
    serde_json::to_value(value)
        .ok()
        .and_then(|v| v.as_str().map(String::from))
        .expect("unit variants serialise to strings")
}

impl Summary {
    pub fn new() -> Self {
        Self {
            targets: 0,
            trials: 0,
            verdicts: BTreeMap::new(),
            needed: BTreeMap::new(),
            covered: Vec::new(),
            never: Vec::new(),
        }
    }

    pub fn add(&mut self, result: &TargetResult) {
        self.targets += 1;
        self.trials += result.trials.len();
        let counts = self.verdicts.entry(name(&result.category)).or_default();
        for trial in &result.trials {
            *counts.entry(name(&trial.verdict)).or_default() += 1;
        }
        if result.category == Category::Primary {
            let exclusive: Vec<u32> =
                result.trials.iter().filter(|t| t.verdict == Verdict::Exclusive).map(|t| t.k).collect();
            if !exclusive.is_empty() {
                self.needed.insert(result.id.clone(), exclusive);
            } else if result.trials.iter().any(|t| t.verdict == Verdict::Shared) {
                self.covered.push(result.id.clone());
            } else {
                self.never.push(result.id.clone());
            }
        }
    }
}

impl Default for Summary {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug)]
pub enum Problem {
    Mismatch(&'static str),
    Fit,
    /// The trials are not the header's scales, in order.
    Scales,
    Spelling { k: u32 },
    Neighbourhood { k: u32, source: NeighbourhoodError },
    /// The number of outcomes differs from the number of probes.
    OutcomeCount { k: u32 },
    Record { k: u32, probe: usize, source: ProbeError },
    /// A record that its probe cannot have produced.
    WrongRecord { k: u32, probe: usize },
    Contradiction { k: u32, probe: usize, source: Contradiction },
    Verdict { k: u32 },
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Problem::Mismatch(field) => write!(f, "{field} differs from its recomputation"),
            Problem::Fit => write!(f, "the centre is a fit"),
            Problem::Scales => write!(f, "the trials are not the header's scales"),
            Problem::Spelling { k } => write!(f, "k={k}: malformed box"),
            Problem::Neighbourhood { k, source } => write!(f, "k={k}: {source}"),
            Problem::OutcomeCount { k } => write!(f, "k={k}: one outcome per probe expected"),
            Problem::Record { k, probe, source } => write!(f, "k={k} probe {probe}: {source}"),
            Problem::WrongRecord { k, probe } => {
                write!(f, "k={k} probe {probe}: the record does not belong to the probe")
            }
            Problem::Contradiction { k, probe, source } => {
                write!(f, "k={k} probe {probe}: {source}")
            }
            Problem::Verdict { k } => write!(f, "k={k}: wrong verdict"),
        }
    }
}

#[derive(Debug)]
pub enum CheckError {
    Transcript(TranscriptError),
    Catalogue(CatalogueError),
    Header(String),
    Target { id: String, problem: Problem },
    Summary,
}

impl fmt::Display for CheckError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CheckError::Transcript(e) => e.fmt(f),
            CheckError::Catalogue(e) => e.fmt(f),
            CheckError::Header(message) => write!(f, "header: {message}"),
            CheckError::Target { id, problem } => write!(f, "target {id}: {problem}"),
            CheckError::Summary => write!(f, "the summary differs from its recomputation"),
        }
    }
}

impl std::error::Error for CheckError {}

impl From<TranscriptError> for CheckError {
    fn from(e: TranscriptError) -> Self {
        CheckError::Transcript(e)
    }
}

/// Scales must be nonempty, strictly increasing and at most `MAX_EXPONENT`.
pub fn valid_scales(scales: &[u32]) -> bool {
    !scales.is_empty()
        && scales.windows(2).all(|w| w[0] < w[1])
        && scales.iter().all(|&k| k <= MAX_EXPONENT)
}

pub fn check_header(header: &Header, catalogue_sha256: &str) -> Result<(), CheckError> {
    let fail = |m: &str| Err(CheckError::Header(m.to_owned()));
    let binding = run::Binding {
        format: &header.format,
        policy: &header.policy,
        executable_sha256: &header.executable_sha256,
        catalogue_sha256: &header.catalogue_sha256,
        components: &header.components,
    };
    binding.check(FORMAT, catalogue_sha256).map_err(CheckError::Header)?;
    if !valid_scales(&header.scales) {
        return fail("invalid scales");
    }
    // Probes follow the components' order; a component has either one
    // attempt probe or distinct cover probes.
    let mut position = 0;
    let probes: BTreeSet<&Label> = header.probes.iter().collect();
    if probes.len() != header.probes.len() {
        return fail("a probe is listed twice");
    }
    for probe in &header.probes {
        let Some(index) = header.components.iter().position(|c| *c == probe.component) else {
            return fail("a probe names an unknown component");
        };
        if index < position {
            return fail("probes are not in the components' order");
        }
        position = index;
        let attempt = Label {
            component: probe.component.clone(),
            cover: None,
        };
        if probe.cover.is_some() && probes.contains(&attempt) {
            return fail("a component has both an attempt and covers");
        }
    }
    let covered: BTreeSet<&String> = header.probes.iter().map(|p| &p.component).collect();
    if covered.len() != header.components.len() {
        return fail("a component has no probe");
    }
    Ok(())
}

/// Checks one target's line against its catalogue entry.
pub fn check_target(
    header: &Header,
    target: &Target,
    result: &TargetResult,
    form: RecordForm,
) -> Result<(), Problem> {
    if result.id != target.id {
        return Err(Problem::Mismatch("id"));
    }
    if result.category != target.category {
        return Err(Problem::Mismatch("category"));
    }
    if result.target != target.target {
        return Err(Problem::Mismatch("target"));
    }
    if result.in_domain != relation::in_domain(&target.centre) {
        return Err(Problem::Mismatch("in_domain"));
    }
    if result.relation != relation::relation(&target.centre) {
        return Err(Problem::Mismatch("relation"));
    }
    if result.relation == Relation::Fit {
        return Err(Problem::Fit);
    }
    let ks: Vec<u32> = result.trials.iter().map(|t| t.k).collect();
    if ks != header.scales {
        return Err(Problem::Scales);
    }
    for trial in &result.trials {
        let k = trial.k;
        let b = unspell(&trial.axes).ok_or(Problem::Spelling { k })?;
        neighbourhood::check(&target.centre, k, &b)
            .map_err(|source| Problem::Neighbourhood { k, source })?;
        if trial.outcomes.len() != header.probes.len() {
            return Err(Problem::OutcomeCount { k });
        }
        let mut successes = Vec::new();
        for (index, (probe, outcome)) in header.probes.iter().zip(&trial.outcomes).enumerate() {
            let Some(record) = &outcome.record else {
                continue;
            };
            let label = checked_label(record, form).map_err(|source| Problem::Record {
                k,
                probe: index,
                source,
            })?;
            let belongs = match &probe.cover {
                Some(_) => label == *probe,
                None => label.component == probe.component,
            };
            if !belongs {
                return Err(Problem::WrongRecord { k, probe: index });
            }
            consistent(&label, result.in_domain, result.relation).map_err(|source| {
                Problem::Contradiction {
                    k,
                    probe: index,
                    source,
                }
            })?;
            successes.push(label);
        }
        if verdict(&result.target, &successes) != trial.verdict {
            return Err(Problem::Verdict { k });
        }
    }
    Ok(())
}

/// A transcript that passed [`check`].
#[derive(Clone, Debug)]
pub struct Checked {
    pub header: Header,
    pub targets: Vec<TargetResult>,
    pub summary: Summary,
}

/// Checks a whole transcript against its target catalogue, the targets on
/// `threads` workers.
pub fn check(
    catalogue: &Targets,
    bytes: &[u8],
    threads: NonZeroUsize,
    form: RecordForm,
) -> Result<Checked, CheckError> {
    let (header, targets, summary) = run::check_frame(
        bytes,
        threads,
        "targets",
        |header: &Header| {
            check_header(header, &catalogue.sha256)?;
            let selected = catalogue.select(header.selection.as_deref()).map_err(CheckError::Catalogue)?;
            Ok((selected, Summary::new()))
        },
        |header, target: &Target, result: &TargetResult| {
            check_target(header, target, result, form).map_err(|problem| CheckError::Target { id: target.id.clone(), problem })
        },
        |summary: &mut Summary, _, result| summary.add(result),
        CheckError::Summary,
    )?;
    Ok(Checked { header, targets, summary })
}
