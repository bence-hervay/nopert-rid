//! The Global component: the candidate, proposed in floating point, is an
//! edge direction of the RID with a sign and a plug vertex, and it holds
//! when the plug vertex lies strictly beyond the hole's full extent in the direction `u × d`
//! throughout the box (every member of its maximum form certified negative
//! by the lemma).
use super::{cell, identity, proposal, Component, Refusal};
use crate::arithmetic::exact::QSqrt5;
use crate::arithmetic::polynomial::Polynomial;
use crate::elimination::lemma;
use crate::elimination::witness::{Edge, MaximumGap};
use crate::elimination::zoom::ZoomFactor;
use crate::problem::configuration::ConfigurationBox;
use crate::problem::geometry::{self, VERTEX_COUNT};
use std::sync::OnceLock;

/// Eliminates boxes by a maximum form: the fixed candidates are the fifteen
/// edge directions, each with one sign, times the 60 plug vertices.
#[derive(Debug)]
pub struct Global {
    /// One ordered edge per edge direction up to sign, the first of its
    /// class in the order of `geometry::edges()`, oriented as listed there.
    directions: Vec<Edge>,
    /// The 60 members of each candidate of [`Global::all`], in its order, on
    /// the identity zoom: they do not depend on the box, so each is built
    /// the first time any thread needs it and then reused.
    members: Vec<OnceLock<Vec<Polynomial>>>,
}

impl Global {
    pub fn new() -> Self {
        let mut directions: Vec<Edge> = Vec::new();
        for &[from, to] in geometry::edges() {
            let edge = Edge::new(from, to).expect("an edge of the RID");
            let parallel = |other: &Edge| geometry::cross(&edge.direction(), &other.direction()).iter().all(QSqrt5::is_zero);
            if !directions.iter().any(parallel) {
                directions.push(edge);
            }
        }
        let members = (0..directions.len() * VERTEX_COUNT).map(|_| OnceLock::new()).collect();
        Self { directions, members }
    }

    /// The fixed list: each direction, in order, with the plug vertices in
    /// vertex order.
    pub fn all(&self) -> Vec<MaximumGap> {
        self.directions
            .iter()
            .flat_map(|edge| (0..VERTEX_COUNT).map(|p| MaximumGap::new(edge.clone(), p).expect("below 60")))
            .collect()
    }
}

impl Component for Global {
    type Witness = MaximumGap;

    /// The maximum form proposed in floating point ([`proposal::global`]),
    /// if any.
    fn candidates(&self, b: &ConfigurationBox) -> Vec<MaximumGap> {
        proposal::global(&self.directions, b).into_iter().collect()
    }

    /// The lemma for the maximum form on the box, as a cell of the identity
    /// zoom with the factor 1: every one of its 60 members certified negative.
    /// The members come from the table for a candidate of the fixed list,
    /// and are built for any other maximum form (a parallel or reversed edge).
    fn holds(&self, b: &ConfigurationBox, maximum: &MaximumGap) -> Result<(), Refusal> {
        let build = || maximum.members(identity().view(), identity().rotation()).expect("members of degree at most 2");
        let built;
        let members = match self.directions.iter().position(|edge| edge == maximum.edge()) {
            Some(k) => self.members[k * VERTEX_COUNT + maximum.plug()].get_or_init(build),
            None => {
                built = build();
                &built
            }
        };
        lemma::check_members(identity(), &cell(b), members, &ZoomFactor::ONE).map_err(Refusal::Lemma)
    }
}

#[cfg(test)]
mod tests;
