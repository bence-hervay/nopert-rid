//! The Exotic and Local components: the candidates are the loaded
//! proved covers in key order, and one holds when it contains the box. Both
//! are the same thin component, keyed by cover names or by cover numbers.
//!
use super::record::{ComponentName, ExoticCover, LocalCover};
use super::{Component, Refusal};
use crate::elimination::proof::ProvedSet;
use crate::problem::configuration::ConfigurationBox;
use std::collections::BTreeMap;
use std::fmt;

/// How the covers of one component are identified in records.
pub trait CoverKey: Copy + Ord + fmt::Display {
    const COMPONENT: ComponentName;
    /// The key whose cover is named `name`, if `name` is such a name.
    fn from_name(name: &str) -> Option<Self>;
}

impl CoverKey for ExoticCover {
    const COMPONENT: ComponentName = ComponentName::Exotic;
    fn from_name(name: &str) -> Option<Self> {
        ExoticCover::from_name(name)
    }
}

impl CoverKey for LocalCover {
    const COMPONENT: ComponentName = ComponentName::Local;
    fn from_name(name: &str) -> Option<Self> {
        LocalCover::from_name(name)
    }
}

/// Why a set of covers cannot form a component.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoversError {
    /// A cover's name is not a name of this component's covers.
    Foreign { component: ComponentName, name: String },
    Duplicate { component: ComponentName, name: String },
}

impl fmt::Display for CoversError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CoversError::Foreign { component, name } => {
                write!(f, "{name:?} is not the name of a {component} cover")
            }
            CoversError::Duplicate { component, name } => {
                write!(f, "{component} cover {name} is given twice")
            }
        }
    }
}

impl std::error::Error for CoversError {}

/// A component eliminating the boxes that lie in one of its proved covers.
pub struct Covered<K: CoverKey> {
    covers: BTreeMap<K, Box<dyn ProvedSet>>,
}

/// Exotic: the covers around the square view, the pentagon event, the
/// arcs, their endpoint and crossing, in the order of [`ExoticCover::ALL`].
pub type Exotic = Covered<ExoticCover>;
/// Local: the numbered covers around aligned views, by increasing number.
pub type Local = Covered<LocalCover>;

impl<K: CoverKey> Covered<K> {
    /// The component over `covers`, each identified by its name. Refuses a
    /// name that is not one of this component's and a name given twice.
    pub fn new(covers: Vec<Box<dyn ProvedSet>>) -> Result<Self, CoversError> {
        let mut map = BTreeMap::new();
        for cover in covers {
            let name = cover.name().to_owned();
            let Some(key) = K::from_name(&name) else {
                return Err(CoversError::Foreign { component: K::COMPONENT, name });
            };
            if map.insert(key, cover).is_some() {
                return Err(CoversError::Duplicate { component: K::COMPONENT, name });
            }
        }
        Ok(Self { covers: map })
    }

    /// The keys of the loaded covers, in increasing order.
    pub fn keys(&self) -> Vec<K> {
        self.covers.keys().copied().collect()
    }
}

impl<K: CoverKey> Component for Covered<K> {
    type Witness = K;

    /// The loaded covers, in increasing key order.
    fn candidates(&self, _: &ConfigurationBox) -> Vec<K> {
        self.keys()
    }

    /// The cover of that key is loaded and contains the box.
    fn holds(&self, b: &ConfigurationBox, key: &K) -> Result<(), Refusal> {
        let cover = self.covers.get(key).ok_or_else(|| Refusal::UnknownCover(key.to_string()))?;
        if cover.contains(b) {
            Ok(())
        } else {
            Err(Refusal::NotContained(key.to_string()))
        }
    }
}

#[cfg(test)]
mod tests;
