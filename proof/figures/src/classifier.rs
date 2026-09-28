//! The crate's components as a slice classifier: the cover catalogue is
//! loaded and verified once, and each box is checked by the parts of the
//! collection a slice asks for, in the collection's order.
use crate::content::slice::{Category, Classifier, Parts};
use rid::components::collection::Collection;
use rid::components::Component;
use rid::problem::configuration::ConfigurationBox;
use rid::search::command::wiring;
use rid::search::BoxError;
use std::num::NonZeroUsize;
use std::sync::atomic::AtomicBool;

/// The first of `parts` that eliminates `b`: Domain, the Exotic covers,
/// Local, Global. With every part this is the collection's own check.
pub fn classify(collection: &Collection, b: &ConfigurationBox, parts: &Parts) -> Option<Category> {
    if parts.domain && collection.domain.check(b).is_some() {
        return Some(Category::Domain);
    }
    if let Some(&cover) = parts.exotic.iter().find(|cover| collection.exotic.holds(b, cover).is_ok()) {
        return Some(Category::Exotic { cover });
    }
    if parts.local && collection.local.check(b).is_some() {
        return Some(Category::Local);
    }
    if parts.global && collection.global.check(b).is_some() {
        return Some(Category::Global);
    }
    None
}

/// The production classifier, its covers verified on `threads` threads.
pub fn production(threads: NonZeroUsize) -> Result<Box<dyn Classifier>, BoxError> {
    let collection = wiring::prepare(threads, &AtomicBool::new(false))?;
    Ok(Box::new(move |b: &ConfigurationBox, parts: &Parts| -> Result<Option<Category>, BoxError> {
        Ok(classify(&collection, b, parts))
    }))
}
