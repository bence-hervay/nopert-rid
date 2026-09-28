//! The production components of the commands: the collection of the
//! components module over the proved covers of the cover catalogue. Its
//! record data is the certificate's record data, and the collection's
//! `holds` is what every reading of a certificate runs.
use super::Components;
use crate::components::collection::Collection;
use crate::components::record::{ExoticCover, RecordData};
use crate::components::covered::{Exotic, Local, CoversError};
use crate::elimination::proof::catalogue::{self, CatalogueError};
use crate::elimination::proof::ProvedSet;
use crate::problem::configuration::ConfigurationBox;
use crate::search::{certificate, BoxError};
use serde_json::Value;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::AtomicBool;

impl certificate::RecordData for RecordData {}

impl Components for Collection {
    type Data = RecordData;

    fn report(&self) -> Value {
        Collection::report(self)
    }

    /// The collection's first witness for the box of `path`.
    fn check(&self, path: &str) -> Result<Option<RecordData>, BoxError> {
        Ok(Collection::check(self, &ConfigurationBox::from_path(path)?))
    }

    fn holds(&self, path: &str, data: &RecordData) -> Result<(), BoxError> {
        Ok(Collection::holds(self, &ConfigurationBox::from_path(path)?, data)?)
    }
}

/// Why the production components could not be prepared.
#[derive(Debug)]
pub enum PrepareError {
    /// A cover of the catalogue failed to load.
    Catalogue(CatalogueError),
    /// A loaded cover's name is not one of its component's names, or repeats.
    Covers(CoversError),
    /// These Exotic covers are missing from the catalogue.
    MissingExotic(Vec<ExoticCover>),
}

impl fmt::Display for PrepareError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PrepareError::Catalogue(error) => write!(f, "{error}"),
            PrepareError::Covers(error) => write!(f, "{error}"),
            PrepareError::MissingExotic(missing) => {
                let names: Vec<&str> = missing.iter().map(|cover| cover.name()).collect();
                write!(f, "the Exotic covers {} are missing from the catalogue", names.join(", "))
            }
        }
    }
}

impl std::error::Error for PrepareError {}

fn boxed<R: ProvedSet + 'static>(covers: Vec<R>) -> Vec<Box<dyn ProvedSet>> {
    covers.into_iter().map(|cover| Box::new(cover) as Box<dyn ProvedSet>).collect()
}

/// Loads and verifies every cover of the catalogue on `threads` threads and
/// builds the collection, each catalogue's covers for its own component.
/// Refuses unless all eight Exotic covers are present: without them no
/// search of the whole root box can finish.
pub fn prepare(threads: NonZeroUsize, stop: &AtomicBool) -> Result<Collection, BoxError> {
    Ok(prepare_collection(threads, stop)?)
}

/// [`prepare`] with its typed error.
pub fn prepare_collection(threads: NonZeroUsize, stop: &AtomicBool) -> Result<Collection, PrepareError> {
    let exotic = catalogue::exotic(threads, stop).map_err(PrepareError::Catalogue)?;
    let local = catalogue::local(threads, stop).map_err(PrepareError::Catalogue)?;
    let exotic = Exotic::new(boxed(exotic)).map_err(PrepareError::Covers)?;
    let local = Local::new(boxed(local)).map_err(PrepareError::Covers)?;
    let loaded = exotic.keys();
    let missing: Vec<ExoticCover> =
        ExoticCover::ALL.into_iter().filter(|cover| !loaded.contains(cover)).collect();
    if !missing.is_empty() {
        return Err(PrepareError::MissingExotic(missing));
    }
    Ok(Collection::new(exotic, local))
}
