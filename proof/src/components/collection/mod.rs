//! The four components in the order Domain, Exotic, Local, Global: a
//! box's check is the first component's witness, and a record holds when
//! the component it names says so.
use super::covered::{CoversError, Exotic, Local};
use super::domain::Domain;
use super::global::Global;
use super::record::{ComponentName, ExoticCover, RecordData};
use super::{Component, Refusal};
use crate::elimination::proof::ProvedSet;
use crate::problem::configuration::ConfigurationBox;
use serde_json::{json, Value};

/// The four components.
pub struct Collection {
    pub domain: Domain,
    pub exotic: Exotic,
    pub local: Local,
    pub global: Global,
}

impl Collection {
    pub fn new(exotic: Exotic, local: Local) -> Self {
        Self {
            domain: Domain,
            exotic,
            local,
            global: Global::new(),
        }
    }

    /// The collection over proved covers identified by their names: the eight
    /// Exotic names and the decimal numbers of Local's covers. A name of
    /// neither kind, or a name given twice, is refused.
    pub fn from_covers(covers: Vec<Box<dyn ProvedSet>>) -> Result<Self, CoversError> {
        let (exotic, local): (Vec<_>, Vec<_>) = covers
            .into_iter()
            .partition(|cover| ExoticCover::from_name(cover.name()).is_some());
        Ok(Self::new(Exotic::new(exotic)?, Local::new(local)?))
    }

    /// The first component's witness for `b`, in the order Domain,
    /// Exotic, Local, Global, as a record.
    pub fn check(&self, b: &ConfigurationBox) -> Option<RecordData> {
        self.domain
            .check(b)
            .map(|inequality| RecordData::Domain { inequality })
            .or_else(|| self.exotic.check(b).map(|cover| RecordData::Exotic { cover }))
            .or_else(|| self.local.check(b).map(|cover| RecordData::Local { cover }))
            .or_else(|| self.global.check(b).map(RecordData::Global))
    }

    /// Whether the record's witness holds on `b` for the component it names:
    /// what every reading of a certificate runs.
    pub fn holds(&self, b: &ConfigurationBox, data: &RecordData) -> Result<(), Refusal> {
        match data {
            RecordData::Domain { inequality } => self.domain.holds(b, inequality),
            RecordData::Exotic { cover } => self.exotic.holds(b, cover),
            RecordData::Local { cover } => self.local.holds(b, cover),
            RecordData::Global(maximum) => self.global.holds(b, maximum),
        }
    }

    /// A summary for `rid prepare`: the order and the loaded covers.
    pub fn report(&self) -> Value {
        let exotic: Vec<String> = self.exotic.keys().iter().map(|k| k.to_string()).collect();
        let local: Vec<u32> = self.local.keys().iter().map(|k| k.0).collect();
        json!({
            "order": ComponentName::ORDER.map(|name| name.to_string()),
            "exotic_covers": exotic,
            "local_covers": local,
        })
    }
}

#[cfg(test)]
mod tests;
