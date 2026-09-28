//! The Domain component: the candidates are the inequalities of `D` that
//! are violated at the box's centre in floating point, and one holds when
//! the lemma certifies it violated throughout the box.
use super::{cell, identity, proposal, Component, Refusal};
use crate::elimination::lemma;
use crate::elimination::witness::{DomainInequality, Witness};
use crate::elimination::zoom::ZoomFactor;
use crate::problem::configuration::ConfigurationBox;

/// Eliminates boxes on which one inequality of `D` is violated throughout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Domain;

impl Component for Domain {
    type Witness = DomainInequality;

    /// The inequalities proposed in floating point ([`proposal::domain`]).
    fn candidates(&self, b: &ConfigurationBox) -> Vec<DomainInequality> {
        proposal::domain(b)
    }

    /// The lemma on the box, as a cell of the identity zoom with the factor
    /// 1, with minus the inequality as the witness.
    fn holds(&self, b: &ConfigurationBox, inequality: &DomainInequality) -> Result<(), Refusal> {
        let witness = Witness::Domain(inequality.clone());
        lemma::check_cell(identity(), &cell(b), &witness, &ZoomFactor::ONE).map_err(Refusal::Lemma)
    }
}

#[cfg(test)]
mod tests;
