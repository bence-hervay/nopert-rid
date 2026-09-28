//! What the experiments ask of the components (the [`Probes`] trait), the
//! labels of records, and the verified probes built on them.
use crate::points::relation::Relation;
use rid::problem::configuration::ConfigurationBox;
use rid::search::certificate::RecordData;
use rid::search::BoxError;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::fmt;

#[cfg(test)]
pub(crate) mod standin;

/// The components as the experiments see them: the collection's decision and
/// each component's own attempt, verification and covers. Components are
/// numbered by their position in [`Probes::components`], the collection's
/// order.
pub trait Probes: Sync {
    /// The record data; its JSON form is an object with a `component` name
    /// and, for a covered record, a `cover` name or number.
    type Data: RecordData + Clone;
    /// The component names, in the order in which the collection tries them.
    fn components(&self) -> Vec<String>;
    /// The collection's decision: the first success in that order. An error
    /// reports a defect of a component (a record its own verification refuses).
    fn decide(&self, b: &ConfigurationBox) -> Result<Option<Self::Data>, BoxError>;
    /// Component `index`'s own attempt on the whole box.
    fn attempt(&self, index: usize, b: &ConfigurationBox) -> Option<Self::Data>;
    /// The records naming component `index`'s covers, in its order; empty
    /// for a component that decides per box.
    fn covers(&self, index: usize) -> Vec<Self::Data>;
    /// Component `index`'s verification of `data` for the box.
    fn verify(&self, index: usize, b: &ConfigurationBox, data: &Self::Data) -> Result<(), BoxError>;
}

/// The component and cover a record names; `cover` is absent for a
/// component that decides per box. Numbers are written in decimal.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Label {
    pub component: String,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub cover: Option<String>,
}

impl Label {
    /// The label of a record's JSON form: its `component` string and its
    /// `cover`, a string or a nonnegative integer, if present.
    pub fn of(record: &Value) -> Result<Self, ProbeError> {
        let invalid = || ProbeError::Label(record.to_string());
        let object = record.as_object().ok_or_else(invalid)?;
        let component = object
            .get("component")
            .and_then(Value::as_str)
            .filter(|c| !c.is_empty())
            .ok_or_else(invalid)?
            .to_owned();
        let cover = match object.get("cover") {
            None => None,
            Some(Value::String(name)) if !name.is_empty() => Some(name.clone()),
            Some(Value::Number(n)) if n.is_u64() => Some(n.to_string()),
            Some(_) => return Err(invalid()),
        };
        Ok(Self { component, cover })
    }

    /// Whether this label meets `target`: the same component and, when the
    /// target names a cover, the same cover.
    pub fn meets(&self, target: &Label) -> bool {
        self.component == target.component
            && (target.cover.is_none() || self.cover == target.cover)
    }
}

impl fmt::Display for Label {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.cover {
            None => f.write_str(&self.component),
            Some(cover) => write!(f, "{} {cover}", self.component),
        }
    }
}

#[derive(Debug)]
pub enum ProbeError {
    /// A record whose JSON form has no valid component name or cover.
    Label(String),
    /// A record the component named by it refuses for the same box: the
    /// component (or the collection) contradicts its own verification.
    Refused { record: String, reason: String },
    /// The collection reported a defect of a component.
    Decision(String),
    /// Component names other than [`ORDER`].
    Names,
    /// A record that is not a canonical record of the components' record type.
    NotCanonical(String),
    /// A decision naming a component that is not in the list.
    UnknownComponent(String),
    /// A component's attempt returned a record of another component.
    Foreign { component: String, record: String },
}

impl fmt::Display for ProbeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProbeError::Label(record) => write!(f, "record {record} has no valid label"),
            ProbeError::Refused { record, reason } => {
                write!(f, "record {record} is refused by its own component: {reason}")
            }
            ProbeError::Decision(reason) => write!(f, "the collection failed: {reason}"),
            ProbeError::Names => write!(f, "the components must be Domain, Exotic, Local, Global in this order"),
            ProbeError::UnknownComponent(name) => write!(f, "unknown component {name:?}"),
            ProbeError::NotCanonical(record) => write!(f, "record {record} is not a canonical record"),
            ProbeError::Foreign { component, record } => {
                write!(f, "{component}'s attempt returned the foreign record {record}")
            }
        }
    }
}

impl std::error::Error for ProbeError {}

/// A record contradicting an exact property of a configuration inside its
/// box: evidence that the component that made it is unsound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Contradiction {
    /// Domain claims that the box misses `D`, but the configuration lies in `D`.
    DomainInD,
    /// Global's gap is negative on the box, so no configuration of it is
    /// weakly contained, but this one is.
    GlobalAtContained,
    /// A configuration in `D` that is a fit was eliminated.
    FitEliminated,
}

impl fmt::Display for Contradiction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Contradiction::DomainInD => "Domain eliminated a box containing a point of D",
            Contradiction::GlobalAtContained => {
                "Global eliminated a box containing a weakly contained configuration"
            }
            Contradiction::FitEliminated => "a box containing a fit in D was eliminated",
        })
    }
}

impl std::error::Error for Contradiction {}

/// Whether a record with this label may hold for a box containing a
/// configuration with these exact properties. The rules follow from what
/// the components claim, whatever their implementation: Domain that the box
/// misses `D`; Global that a support gap is negative, which no weakly
/// contained configuration allows; every component that the box holds no
/// fit in `D`.
pub fn consistent(label: &Label, in_domain: bool, relation: Relation) -> Result<(), Contradiction> {
    if relation == Relation::Fit && in_domain {
        return Err(Contradiction::FitEliminated);
    }
    match label.component.as_str() {
        DOMAIN if in_domain => Err(Contradiction::DomainInD),
        GLOBAL if relation != Relation::Poke => Err(Contradiction::GlobalAtContained),
        _ => Ok(()),
    }
}

/// The names of the two components whose claims [`consistent`] checks.
pub const DOMAIN: &str = "Domain";
pub const GLOBAL: &str = "Global";

/// The components' names in their decided order. The contradiction rules
/// are keyed on `DOMAIN` and `GLOBAL`, so exactly these names are accepted:
/// a renamed component would silently switch a rule off.
pub const ORDER: [&str; 4] = [DOMAIN, "Exotic", "Local", GLOBAL];

/// Whether `names` are the decided component names in order.
pub fn decided(names: &[String]) -> bool {
    names.iter().map(String::as_str).eq(ORDER)
}

/// The component names, which must be [`ORDER`].
pub fn names<P: Probes>(probes: &P) -> Result<Vec<String>, ProbeError> {
    let names = probes.components();
    if !decided(&names) {
        return Err(ProbeError::Names);
    }
    Ok(names)
}

/// How the checkers read a record: whether its JSON form is a canonical
/// record of the components' record type.
pub type RecordForm = fn(&Value) -> bool;

/// The record form of `R`: the value parses as an `R` and serialises back to
/// exactly itself (so cover 9 and cover "9" are not the same record).
pub fn canonical<R: Serialize + DeserializeOwned>(value: &Value) -> bool {
    serde_json::from_value::<R>(value.clone()).ok().and_then(|r| serde_json::to_value(r).ok()).as_ref() == Some(value)
}

/// The label of a record read by a checker: refused unless the record has
/// the form `form`.
pub fn checked_label(record: &Value, form: RecordForm) -> Result<Label, ProbeError> {
    if !form(record) {
        return Err(ProbeError::NotCanonical(record.to_string()));
    }
    Label::of(record)
}

/// A record in JSON form with its label.
#[derive(Clone, Debug, PartialEq)]
pub struct Found {
    pub label: Label,
    pub record: Value,
}

fn found<D: Serialize>(data: &D) -> Result<Found, ProbeError> {
    let record = serde_json::to_value(data).map_err(|e| ProbeError::Label(e.to_string()))?;
    Ok(Found {
        label: Label::of(&record)?,
        record,
    })
}

/// Verifies `data` with component `index`.
fn verified<P: Probes>(
    probes: &P,
    index: usize,
    b: &ConfigurationBox,
    data: &P::Data,
    found: &Found,
) -> Result<(), ProbeError> {
    probes
        .verify(index, b, data)
        .map_err(|reason| ProbeError::Refused {
            record: found.record.to_string(),
            reason: reason.to_string(),
        })
}

/// Component `index`'s attempt, which must name that component and pass
/// its verification.
pub fn attempt<P: Probes>(
    probes: &P,
    names: &[String],
    index: usize,
    b: &ConfigurationBox,
) -> Result<Option<Found>, ProbeError> {
    let Some(data) = probes.attempt(index, b) else {
        return Ok(None);
    };
    let result = found(&data)?;
    if result.label.component != names[index] {
        return Err(ProbeError::Foreign {
            component: names[index].clone(),
            record: result.record.to_string(),
        });
    }
    verified(probes, index, b, &data, &result)?;
    Ok(Some(result))
}

/// The collection's decision, which must pass the verification of the
/// component it names.
pub fn decide<P: Probes>(
    probes: &P,
    names: &[String],
    b: &ConfigurationBox,
) -> Result<Option<Found>, ProbeError> {
    let decided = probes
        .decide(b)
        .map_err(|reason| ProbeError::Decision(reason.to_string()))?;
    let Some(data) = decided else {
        return Ok(None);
    };
    confirm(probes, names, b, &data).map(Some)
}

/// A decided record, which must pass the verification of the component it
/// names.
pub fn confirm<P: Probes>(
    probes: &P,
    names: &[String],
    b: &ConfigurationBox,
    data: &P::Data,
) -> Result<Found, ProbeError> {
    let data = data.clone();
    let result = found(&data)?;
    let index = names
        .iter()
        .position(|n| *n == result.label.component)
        .ok_or_else(|| ProbeError::UnknownComponent(result.label.component.clone()))?;
    verified(probes, index, b, &data, &result)?;
    Ok(result)
}

/// One elementary question an experiment can ask about a box: a component's
/// attempt, or whether one cover of a covered component contains it.
#[derive(Clone, Debug)]
pub struct Probe<D> {
    pub component: usize,
    /// The label of the probe: the component alone for an attempt, the
    /// component and cover for a cover.
    pub label: Label,
    /// The cover's record, or `None` for an attempt.
    pub cover: Option<D>,
}

/// Every probe, in the collection's order: each cover of a covered
/// component, and the attempt of each other component.
pub fn probes<P: Probes>(probes: &P) -> Result<Vec<Probe<P::Data>>, ProbeError> {
    let names = probes.components();
    let mut out = Vec::new();
    for (component, name) in names.iter().enumerate() {
        let covers = probes.covers(component);
        if covers.is_empty() {
            out.push(Probe {
                component,
                label: Label {
                    component: name.clone(),
                    cover: None,
                },
                cover: None,
            });
        }
        for data in covers {
            let result = found(&data)?;
            if result.label.component != *name || result.label.cover.is_none() {
                return Err(ProbeError::Foreign {
                    component: name.clone(),
                    record: result.record.to_string(),
                });
            }
            out.push(Probe {
                component,
                label: result.label,
                cover: Some(data),
            });
        }
    }
    Ok(out)
}

/// Asks one probe about the box: the verified record, or `None`.
pub fn ask<P: Probes>(
    probes: &P,
    names: &[String],
    probe: &Probe<P::Data>,
    b: &ConfigurationBox,
) -> Result<Option<Found>, ProbeError> {
    match &probe.cover {
        None => attempt(probes, names, probe.component, b),
        Some(data) => Ok(probes.verify(probe.component, b, data).is_ok().then(|| Found {
            label: probe.label.clone(),
            record: serde_json::to_value(data).expect("the cover record serialised before"),
        })),
    }
}

#[cfg(test)]
mod tests;
