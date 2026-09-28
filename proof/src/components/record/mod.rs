//! The typed data of a certificate record: which component eliminated a box
//! and what it needs to check that again. One canonical JSON form per value.
//!
use crate::elimination::witness::{DomainInequality, MaximumGap};
use serde::{Deserialize, Serialize};
use std::fmt;

/// The four components; `RecordData` spells them the same way.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ComponentName {
    Domain,
    Exotic,
    Local,
    Global,
}

impl ComponentName {
    /// The order in which a box is checked: Domain, Exotic, Local, Global.
    pub const ORDER: [ComponentName; 4] = [
        ComponentName::Domain,
        ComponentName::Exotic,
        ComponentName::Local,
        ComponentName::Global,
    ];
}

impl fmt::Display for ComponentName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

/// The named covers of the Exotic component. A cover along the arcs is
/// written in the arc-plane coordinates of one sign `σ = ±1`, and its name
/// ends in that sign.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub enum ExoticCover {
    Square,
    Pentagon,
    ArcPlus,
    ArcMinus,
    EndpointPlus,
    EndpointMinus,
    CrossingPlus,
    CrossingMinus,
}

impl ExoticCover {
    /// All covers, in the component's candidate order.
    pub const ALL: [ExoticCover; 8] = [
        ExoticCover::Square,
        ExoticCover::Pentagon,
        ExoticCover::ArcPlus,
        ExoticCover::ArcMinus,
        ExoticCover::EndpointPlus,
        ExoticCover::EndpointMinus,
        ExoticCover::CrossingPlus,
        ExoticCover::CrossingMinus,
    ];

    /// The cover's name, as in records and cover files.
    pub fn name(self) -> &'static str {
        match self {
            ExoticCover::Square => "square",
            ExoticCover::Pentagon => "pentagon",
            ExoticCover::ArcPlus => "arc+",
            ExoticCover::ArcMinus => "arc-",
            ExoticCover::EndpointPlus => "endpoint+",
            ExoticCover::EndpointMinus => "endpoint-",
            ExoticCover::CrossingPlus => "crossing+",
            ExoticCover::CrossingMinus => "crossing-",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|cover| cover.name() == name)
    }
}

impl TryFrom<String> for ExoticCover {
    type Error = String;
    fn try_from(name: String) -> Result<Self, String> {
        Self::from_name(&name).ok_or_else(|| format!("unknown Exotic cover {name:?}"))
    }
}

impl From<ExoticCover> for String {
    fn from(cover: ExoticCover) -> Self {
        cover.name().to_string()
    }
}

impl fmt::Display for ExoticCover {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// The number of a Local cover. Its name is the number in decimal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LocalCover(pub u32);

impl LocalCover {
    /// The cover whose name is `name`: a decimal number without sign or
    /// leading zeros (so every number has one name).
    pub fn from_name(name: &str) -> Option<Self> {
        let number: u32 = name.parse().ok()?;
        (number.to_string() == name).then_some(LocalCover(number))
    }
}

impl fmt::Display for LocalCover {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// The data of one record, internally tagged by the component:
///
/// ```text
/// {"component":"Domain","inequality":13}
/// {"component":"Exotic","cover":"crossing+"}
/// {"component":"Local","cover":9}
/// {"component":"Global","edge":[16,17],"vertex":3}
/// ```
///
/// Parsing checks ranges (inequality below 73, edges of the RID, vertices
/// below 60, known cover names); whether the record is true of its box is
/// decided only by the `holds` of the component it names.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "component", deny_unknown_fields)]
pub enum RecordData {
    /// The box misses `D`: this inequality of `D` is violated throughout it.
    Domain { inequality: DomainInequality },
    /// The box lies in this proved Exotic cover.
    Exotic { cover: ExoticCover },
    /// The box lies in this proved Local cover.
    Local { cover: LocalCover },
    /// The plug vertex lies beyond the hole's full extent in the edge's
    /// direction throughout the box (this maximum form's members certified
    /// negative by the lemma on the identity zoom).
    Global(MaximumGap),
}

impl RecordData {
    pub fn component(&self) -> ComponentName {
        match self {
            RecordData::Domain { .. } => ComponentName::Domain,
            RecordData::Exotic { .. } => ComponentName::Exotic,
            RecordData::Local { .. } => ComponentName::Local,
            RecordData::Global(_) => ComponentName::Global,
        }
    }
}

#[cfg(test)]
mod tests;
