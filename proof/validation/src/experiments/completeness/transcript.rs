//! The completeness transcript: its typed lines and the check of a whole
//! transcript against its catalogue.
use super::schedule::{self, ScheduleError, MAX_K};
use crate::experiments::probe::{checked_label, consistent, Contradiction, Label, ProbeError, RecordForm};
use crate::experiments::run::{self, TranscriptError};
use crate::points::catalogue::{Catalogue, CatalogueError, Point};
use crate::points::neighbourhood::{self, NeighbourhoodError};
use crate::points::relation::{self, Relation};
use rid::arithmetic::exact::{parse_rational, Interval};
use rid::problem::configuration::{ConfigurationBox, AXES};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::num::NonZeroUsize;

pub const FORMAT: &str = "rid-completeness/1";

/// A box as its five intervals `["lo", "hi"]` in canonical rational spelling.
pub type Axes = [[String; 2]; AXES];

pub fn spell(b: &ConfigurationBox) -> Axes {
    b.axes()
        .clone()
        .map(|a| [a.lo().to_string(), a.hi().to_string()])
}

/// The box of a spelling; `None` for a malformed number or reversed interval.
pub fn unspell(axes: &Axes) -> Option<ConfigurationBox> {
    let mut intervals = Vec::with_capacity(AXES);
    for [lo, hi] in axes {
        let (lo, hi) = (parse_rational(lo).ok()?, parse_rational(hi).ok()?);
        intervals.push(Interval::new(lo, hi).ok()?);
    }
    Some(ConfigurationBox::new(intervals.try_into().ok()?))
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    pub format: String,
    /// The `rid` crate's policy fingerprint.
    pub policy: String,
    pub executable_sha256: String,
    pub catalogue_sha256: String,
    /// The component names in the collection's order.
    pub components: Vec<String>,
    /// The finest radius exponent.
    pub max_k: u32,
    /// The selected identifiers, or `null` for the whole catalogue.
    #[serde(deserialize_with = "crate::config::required")]
    pub selection: Option<Vec<String>>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Neighbourhood {
    pub k: u32,
    pub axes: Axes,
}

/// One tested radius `2^-k`: the verified record, or `null` for a refusal.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trial {
    pub k: u32,
    #[serde(deserialize_with = "crate::config::required")]
    pub record: Option<Value>,
    pub milliseconds: u64,
}

/// The trials of one schedule in the order they ran, and what follows from
/// them: the smallest accepted exponent (largest radius), its record, and
/// the refused exponents finer than it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Schedule {
    pub trials: Vec<Trial>,
    #[serde(deserialize_with = "crate::config::required")]
    pub best_k: Option<u32>,
    #[serde(deserialize_with = "crate::config::required")]
    pub record: Option<Value>,
    pub finer_refusals: Vec<u32>,
}

impl Schedule {
    pub fn of(trials: Vec<Trial>) -> Self {
        let outcomes: Vec<(u32, bool)> = trials.iter().map(|t| (t.k, t.record.is_some())).collect();
        let (best_k, finer_refusals) = schedule::summary(&outcomes);
        let record = best_k.and_then(|k| {
            trials
                .iter()
                .find(|t| t.k == k && t.record.is_some())
                .and_then(|t| t.record.clone())
        });
        Self {
            trials,
            best_k,
            record,
            finer_refusals,
        }
    }
}

/// The results for one catalogue point: the collection's schedule and one
/// schedule per component, in the header's order.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointResult {
    pub id: String,
    pub group: String,
    pub in_domain: bool,
    pub relation: Relation,
    /// Every tested neighbourhood, by increasing `k`.
    pub neighbourhoods: Vec<Neighbourhood>,
    pub decision: Schedule,
    pub components: Vec<Schedule>,
    /// Whether the decision's component at the best radius is among the
    /// catalogue's expected components; `null` when unresolved.
    #[serde(deserialize_with = "crate::config::required")]
    pub expected: Option<bool>,
    pub milliseconds: u64,
}

/// Whether the component of `record` is among `expected`.
pub fn is_expected(record: &Option<Value>, expected: &[String]) -> Option<bool> {
    record.as_ref().map(|r| {
        let label = Label::of(r).expect("a checked record has a label");
        expected.contains(&label.component)
    })
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub points: usize,
    /// Points whose decision schedule never succeeded.
    pub unresolved: usize,
    /// Points whose decision schedule has a refusal finer than its best radius.
    pub finer_refusals: usize,
    /// Points decided by a component outside their expected list.
    pub unexpected: usize,
    /// The labels of the decisions at the best radius, with their counts.
    pub reasons: BTreeMap<String, usize>,
    /// For each component, the points at which its own schedule succeeded.
    pub components: BTreeMap<String, usize>,
}

impl Summary {
    pub fn new(components: &[String]) -> Self {
        Self {
            points: 0,
            unresolved: 0,
            finer_refusals: 0,
            unexpected: 0,
            reasons: BTreeMap::new(),
            components: components.iter().map(|c| (c.clone(), 0)).collect(),
        }
    }

    /// Counts one point whose labels were already checked.
    pub fn add(&mut self, point: &PointResult, names: &[String]) {
        self.points += 1;
        match &point.decision.record {
            None => self.unresolved += 1,
            Some(record) => {
                let label = Label::of(record).expect("a checked record has a label");
                *self.reasons.entry(label.to_string()).or_default() += 1;
            }
        }
        if !point.decision.finer_refusals.is_empty() {
            self.finer_refusals += 1;
        }
        if point.expected == Some(false) {
            self.unexpected += 1;
        }
        for (name, schedule) in names.iter().zip(&point.components) {
            if schedule.best_k.is_some() {
                *self.components.entry(name.clone()).or_default() += 1;
            }
        }
    }
}

/// Whether the collection's decision at one radius agrees with the
/// components' own attempts at that radius. The collection returns the
/// first success in its order, so wherever a component was tested at the
/// same radius: each component before the decided one refused, the decided
/// one returned the very same record, and without a decision every tested
/// component refused. `decision` is the decided component's index and
/// record; `tested[i]` is component `i`'s record where it was tested
/// (`Some(None)` for a refusal).
pub fn agrees(decision: Option<(usize, &Value)>, tested: &[Option<Option<&Value>>]) -> bool {
    tested.iter().enumerate().all(|(i, outcome)| match (decision, outcome) {
        (_, None) => true,
        (None, Some(record)) => record.is_none(),
        (Some((d, _)), Some(record)) if i < d => record.is_none(),
        (Some((d, decided)), Some(record)) if i == d => *record == Some(decided),
        (Some(_), Some(_)) => true,
    })
}

/// The first radius exponent of the decision schedule at which it does not
/// agree with the component schedules (see [`agrees`]); a decision naming
/// an unknown component never agrees.
pub fn disagreement(names: &[String], decision: &Schedule, components: &[Schedule]) -> Option<u32> {
    /// The record of `schedule` at `k`, where it was tested.
    fn at(schedule: &Schedule, k: u32) -> Option<Option<&Value>> {
        schedule.trials.iter().find(|t| t.k == k).map(|t| t.record.as_ref())
    }
    decision.trials.iter().map(|t| t.k).find(|&k| {
        let decided = match at(decision, k).flatten() {
            None => None,
            Some(record) => {
                let label = Label::of(record).ok();
                let index = label.and_then(|l| names.iter().position(|n| *n == l.component));
                match index {
                    None => return true,
                    Some(index) => Some((index, record)),
                }
            }
        };
        let tested: Vec<Option<Option<&Value>>> = components.iter().map(|c| at(c, k)).collect();
        !agrees(decided, &tested)
    })
}

/// What is wrong with one point's line.
#[derive(Debug)]
pub enum Problem {
    /// A field differs from the catalogue or from its recomputation.
    Mismatch(&'static str),
    /// A centre that is a fit.
    Fit,
    /// Neighbourhood exponents are not strictly increasing.
    Order,
    Neighbourhood { k: u32, source: NeighbourhoodError },
    /// A malformed interval spelling.
    Spelling { k: u32 },
    /// A trial's neighbourhood is not listed, or a listed one is not used.
    Unlisted { k: u32 },
    Unused { k: u32 },
    /// The number of component schedules differs from the header's list.
    ComponentCount,
    /// The trials of a schedule do not follow the schedule.
    Schedule { schedule: String, source: ScheduleError },
    /// A schedule stops before its end.
    Unfinished { schedule: String },
    /// The best exponent, its record or the finer refusals are wrong.
    Derived { schedule: String },
    Record { schedule: String, k: u32, source: ProbeError },
    /// A record names a component other than its schedule's, or an unknown one.
    WrongComponent { schedule: String, k: u32 },
    Contradiction { schedule: String, k: u32, source: Contradiction },
    /// The decision and the component schedules disagree at radius `2^-k`.
    Disagreement { k: u32 },
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Problem::Mismatch(field) => write!(f, "{field} differs from its recomputation"),
            Problem::Fit => write!(f, "the centre is a fit"),
            Problem::Order => write!(f, "neighbourhoods are not listed by increasing k"),
            Problem::Neighbourhood { k, source } => write!(f, "k={k}: {source}"),
            Problem::Spelling { k } => write!(f, "k={k}: malformed box"),
            Problem::Unlisted { k } => write!(f, "k={k} is tested but has no neighbourhood"),
            Problem::Unused { k } => write!(f, "the neighbourhood k={k} is never tested"),
            Problem::ComponentCount => write!(f, "wrong number of component schedules"),
            Problem::Schedule { schedule, source } => write!(f, "{schedule}: {source}"),
            Problem::Unfinished { schedule } => write!(f, "{schedule}: the schedule is unfinished"),
            Problem::Derived { schedule } => {
                write!(f, "{schedule}: best radius, record or finer refusals are wrong")
            }
            Problem::Record { schedule, k, source } => write!(f, "{schedule} k={k}: {source}"),
            Problem::WrongComponent { schedule, k } => {
                write!(f, "{schedule} k={k}: the record names the wrong component")
            }
            Problem::Contradiction { schedule, k, source } => {
                write!(f, "{schedule} k={k}: {source}")
            }
            Problem::Disagreement { k } => {
                write!(f, "k={k}: the decision is not the first component success")
            }
        }
    }
}

#[derive(Debug)]
pub enum CheckError {
    Transcript(TranscriptError),
    Catalogue(CatalogueError),
    Header(String),
    Point { id: String, problem: Problem },
    Summary,
}

impl fmt::Display for CheckError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CheckError::Transcript(e) => e.fmt(f),
            CheckError::Catalogue(e) => e.fmt(f),
            CheckError::Header(message) => write!(f, "header: {message}"),
            CheckError::Point { id, problem } => write!(f, "point {id}: {problem}"),
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

/// A transcript that passed [`check`].
#[derive(Clone, Debug)]
pub struct Checked {
    pub header: Header,
    pub points: Vec<PointResult>,
    pub summary: Summary,
}

/// Checks the header's own fields.
pub fn check_header(header: &Header, catalogue_sha256: &str) -> Result<(), CheckError> {
    let binding = run::Binding {
        format: &header.format,
        policy: &header.policy,
        executable_sha256: &header.executable_sha256,
        catalogue_sha256: &header.catalogue_sha256,
        components: &header.components,
    };
    binding.check(FORMAT, catalogue_sha256).map_err(CheckError::Header)?;
    if !(1..=MAX_K).contains(&header.max_k) {
        return Err(CheckError::Header("max_k out of range".into()));
    }
    Ok(())
}

/// Checks one schedule's trials against the schedule rule, the listed
/// neighbourhoods, the records' components and the exact properties of the
/// centre.
fn check_schedule(
    name: &str,
    schedule: &Schedule,
    max_k: u32,
    component: Option<usize>,
    names: &[String],
    listed: &BTreeMap<u32, ConfigurationBox>,
    in_domain: bool,
    relation: Relation,
    form: RecordForm,
) -> Result<(), Problem> {
    let named = || name.to_owned();
    let outcomes: Vec<(u32, bool)> = schedule
        .trials
        .iter()
        .map(|t| (t.k, t.record.is_some()))
        .collect();
    match schedule::next(max_k, &outcomes) {
        Ok(None) => {}
        Ok(Some(_)) => return Err(Problem::Unfinished { schedule: named() }),
        Err(source) => return Err(Problem::Schedule { schedule: named(), source }),
    }
    for trial in &schedule.trials {
        let k = trial.k;
        if !listed.contains_key(&k) {
            return Err(Problem::Unlisted { k });
        }
        if let Some(record) = &trial.record {
            let label = checked_label(record, form).map_err(|source| Problem::Record {
                schedule: named(),
                k,
                source,
            })?;
            let index = names.iter().position(|n| *n == label.component);
            if index.is_none() || (component.is_some() && component != index) {
                return Err(Problem::WrongComponent { schedule: named(), k });
            }
            consistent(&label, in_domain, relation).map_err(|source| {
                Problem::Contradiction {
                    schedule: named(),
                    k,
                    source,
                }
            })?;
        }
    }
    if Schedule::of(schedule.trials.clone()) != *schedule {
        return Err(Problem::Derived { schedule: named() });
    }
    Ok(())
}

/// Checks one point's line against its catalogue entry; records must have
/// the form `form`.
pub fn check_point(header: &Header, point: &Point, result: &PointResult, form: RecordForm) -> Result<(), Problem> {
    if result.id != point.id {
        return Err(Problem::Mismatch("id"));
    }
    if result.group != point.group {
        return Err(Problem::Mismatch("group"));
    }
    if result.in_domain != relation::in_domain(&point.centre) {
        return Err(Problem::Mismatch("in_domain"));
    }
    if result.relation != relation::relation(&point.centre) {
        return Err(Problem::Mismatch("relation"));
    }
    if result.relation == Relation::Fit {
        return Err(Problem::Fit);
    }
    let mut listed = BTreeMap::new();
    for n in &result.neighbourhoods {
        if listed.keys().next_back().is_some_and(|&last| last >= n.k) {
            return Err(Problem::Order);
        }
        let b = unspell(&n.axes).ok_or(Problem::Spelling { k: n.k })?;
        neighbourhood::check(&point.centre, n.k, &b)
            .map_err(|source| Problem::Neighbourhood { k: n.k, source })?;
        listed.insert(n.k, b);
    }
    let names = &header.components;
    if result.components.len() != names.len() {
        return Err(Problem::ComponentCount);
    }
    let check = |name: &str, schedule: &Schedule, component: Option<usize>| {
        check_schedule(
            name,
            schedule,
            header.max_k,
            component,
            names,
            &listed,
            result.in_domain,
            result.relation,
            form,
        )
    };
    check("decision", &result.decision, None)?;
    for (i, (name, schedule)) in names.iter().zip(&result.components).enumerate() {
        check(name, schedule, Some(i))?;
    }
    let used: BTreeSet<u32> = std::iter::once(&result.decision)
        .chain(&result.components)
        .flat_map(|s| s.trials.iter().map(|t| t.k))
        .collect();
    if let Some(&k) = listed.keys().find(|k| !used.contains(k)) {
        return Err(Problem::Unused { k });
    }
    if let Some(k) = disagreement(names, &result.decision, &result.components) {
        return Err(Problem::Disagreement { k });
    }
    if result.expected != is_expected(&result.decision.record, &point.expected) {
        return Err(Problem::Mismatch("expected"));
    }
    Ok(())
}

/// Checks a whole transcript against its catalogue: the header binding, one
/// line per selected point in catalogue order (see [`check_point`], run on
/// `threads` workers, records of the form `form`), and the summary.
pub fn check(
    catalogue: &Catalogue,
    bytes: &[u8],
    threads: NonZeroUsize,
    form: RecordForm,
) -> Result<Checked, CheckError> {
    let (header, points, summary) = run::check_frame(
        bytes,
        threads,
        "points",
        |header: &Header| {
            check_header(header, &catalogue.sha256)?;
            let selected = catalogue.select(header.selection.as_deref()).map_err(CheckError::Catalogue)?;
            Ok((selected, Summary::new(&header.components)))
        },
        |header, point: &Point, result: &PointResult| {
            check_point(header, point, result, form).map_err(|problem| CheckError::Point { id: point.id.clone(), problem })
        },
        |summary: &mut Summary, header, result| summary.add(result, &header.components),
        CheckError::Summary,
    )?;
    Ok(Checked { header, points, summary })
}
