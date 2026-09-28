//! The components: Domain, Exotic, Local and Global. Each is a fixed,
//! ordered, finite list of candidate witnesses and one exact condition; a
//! box's check is the first candidate that holds.
pub mod record;
pub mod domain;
pub mod covered;
pub mod global;
pub mod proposal;
pub mod collection;

use crate::elimination::lemma::CellRefusal;
use crate::elimination::zoom::{Zoom, ZoomCell};
use crate::problem::configuration::ConfigurationBox;
use std::fmt;
use std::sync::OnceLock;

/// One way of eliminating boxes: a fixed, ordered, finite list of candidate
/// witnesses and one exact condition, `holds`.
///
/// Contract: `holds(b, w)` returns `Ok`
/// only if no configuration of `b` that lies in `D` is a fit.
///
/// `candidates(b)` only proposes what to try, in order; Domain and Global
/// propose in floating point ([`proposal`]). A candidate proposed wrongly or
/// not at all costs a split of the box, never soundness: `check(b)` returns
/// a witness only after `holds` accepted it, and a record is verified by
/// `holds` alone.
pub trait Component {
    type Witness;
    fn candidates(&self, b: &ConfigurationBox) -> Vec<Self::Witness>;
    fn holds(&self, b: &ConfigurationBox, witness: &Self::Witness) -> Result<(), Refusal>;
    /// The first candidate that holds on `b`, if any.
    fn check(&self, b: &ConfigurationBox) -> Option<Self::Witness> {
        self.candidates(b).into_iter().find(|w| self.holds(b, w).is_ok())
    }
}

/// Why a witness does not hold on a box.
#[derive(Debug)]
pub enum Refusal {
    /// Domain or Global: the lemma refuses the box as a cell of the identity zoom.
    Lemma(CellRefusal),
    /// Exotic or Local: no cover of this name is loaded.
    UnknownCover(String),
    /// Exotic or Local: the cover of this name does not contain the box.
    NotContained(String),
}

impl fmt::Display for Refusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Refusal::Lemma(refusal) => write!(f, "the lemma refuses the box: {refusal}"),
            Refusal::UnknownCover(cover) => write!(f, "no cover {cover} is loaded"),
            Refusal::NotContained(cover) => write!(f, "cover {cover} does not contain the box"),
        }
    }
}

impl std::error::Error for Refusal {}

/// The identity zoom on the root box `B₀`: every search box is a cell of
/// it ([`cell`]), with configuration coordinates, view `(s, t, 1)` and no
/// scale variable. Domain and Global apply the lemma on it.
pub(crate) fn identity() -> &'static Zoom {
    static IDENTITY: OnceLock<Zoom> = OnceLock::new();
    IDENTITY.get_or_init(|| Zoom::identity(cell(&ConfigurationBox::root())).expect("B₀ has positive widths"))
}

/// The search box `b` as a cell of the identity zoom.
pub(crate) fn cell(b: &ConfigurationBox) -> ZoomCell {
    ZoomCell::new(b.axes().clone())
}

#[cfg(test)]
mod tests;
