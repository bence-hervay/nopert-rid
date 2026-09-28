//! Which witnesses Domain and Global try on a box, chosen in floating point.
//! A proposal never decides anything: a record is accepted only by the exact
//! `holds`, so a rounding error can cost a split, never soundness.
use crate::arithmetic::polynomial::VARIABLES;
use crate::elimination::float::{self, Dense, Vector};
use crate::elimination::witness::{DomainInequality, Edge, MaximumGap};
use crate::problem::configuration::{ConfigurationBox, AXES, R, S, T};
use crate::problem::domain::{self, CONSTRAINT_COUNT};
use crate::problem::geometry::VERTEX_COUNT;
use std::sync::OnceLock;

/// Domain: the inequalities `c_i` that are positive at the box's centre, in
/// index order. Domain needs `c_i > 0` throughout the box, so an inequality
/// that is not positive at the centre cannot hold.
pub fn domain(b: &ConfigurationBox) -> Vec<DomainInequality> {
    static INEQUALITIES: OnceLock<Vec<Dense>> = OnceLock::new();
    let inequalities = INEQUALITIES.get_or_init(|| domain::affine_polynomials().iter().map(Dense::new).collect());
    let centre: [f64; VARIABLES] = b.midpoint().map(|x| float::rational(&x));
    (0..CONSTRAINT_COUNT)
        .filter(|&i| inequalities[i].evaluate(&centre) > 0.0)
        .map(|i| DomainInequality::new(i).expect("below 73"))
        .collect()
}

/// Global: the maximum form whose plug vertex lies farthest beyond the
/// hole's extent at the box's worst corner, if it lies beyond there at all.
/// The margin of `(d, p)` is the least, over the 32 corners of the box, of
/// the screen distance `(n·R(r)p − max_w n·w) / |n|` with `n = u × d`;
/// ties go to the first in the order of `directions` and plug vertices.
pub fn global(directions: &[Edge], b: &ConfigurationBox) -> Option<MaximumGap> {
    let vertices = float::vertices();
    let axes = b.axes();
    let ends: [[f64; 2]; AXES] = std::array::from_fn(|j| [float::rational(axes[j].lo()), float::rational(axes[j].hi())]);
    let corner = |mask: usize| -> [f64; AXES] { std::array::from_fn(|j| ends[j][mask >> j & 1]) };
    // Corner `mask` has the view of corner `mask & 3`: only S and T move it.
    debug_assert_eq!((S, T), (0, 1));
    let views: Vec<Vector> = (0..4).map(|mask| corner(mask)).map(|c| [c[S], c[T], 1.0]).collect();
    let images: Vec<Vec<Vector>> = (0..1 << AXES)
        .map(|mask| {
            let c = corner(mask);
            let rotation = float::rotation(&R.map(|j| c[j]));
            vertices.iter().map(|p| float::apply(&rotation, p)).collect()
        })
        .collect();
    let mut best: Option<(f64, MaximumGap)> = None;
    for edge in directions {
        let d = float::point(&edge.direction());
        let normals: Vec<Vector> = views.iter().map(|u| float::cross(u, &d)).collect();
        let extents: Vec<f64> = normals
            .iter()
            .map(|n| vertices.iter().map(|w| float::dot(n, w)).fold(f64::NEG_INFINITY, f64::max))
            .collect();
        for plug in 0..VERTEX_COUNT {
            let margin = images
                .iter()
                .enumerate()
                .map(|(mask, image)| {
                    let n = &normals[mask & 3];
                    (float::dot(n, &image[plug]) - extents[mask & 3]) / float::dot(n, n).sqrt()
                })
                .fold(f64::INFINITY, f64::min);
            if best.as_ref().map_or(true, |(most, _)| margin > *most) {
                best = Some((margin, MaximumGap::new(edge.clone(), plug).expect("below 60")));
            }
        }
    }
    best.filter(|(margin, _)| *margin > 0.0).map(|(_, maximum)| maximum)
}

#[cfg(test)]
mod tests;
