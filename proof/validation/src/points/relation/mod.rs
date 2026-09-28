//! Exact properties of one configuration, computed from the solid and the
//! domain alone: whether it lies in `D`, and whether its plug shadow lies in
//! the interior of the hole shadow (fit), in the closed hole shadow only
//! (touch) or not in it (poke). The experiments use them to catch unsound
//! eliminations.
use crate::points::catalogue::Centre;
use rid::arithmetic::exact::QSqrt5;
use rid::problem::domain;
use rid::problem::geometry::{dot, rotate, screen, vertices};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

/// How the plug shadow of a configuration lies relative to the hole shadow.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Relation {
    /// In the interior: the configuration is a fit.
    Fit,
    /// In the closed hole shadow, touching its boundary.
    Touch,
    /// Not in the closed hole shadow.
    Poke,
}

/// Whether the configuration lies in the closed domain `D`: every one of its
/// inequalities `c ≤ 0` holds (the centre is in `B₀` by construction).
pub fn in_domain(centre: &Centre) -> bool {
    domain::affine_polynomials()
        .iter()
        .all(|c| c.evaluate(centre.coordinates()).sign() != Ordering::Greater)
}

/// A point of the screen plane.
pub type Planar = [QSqrt5; 2];

/// `(b - a) × (p - a)`: positive when `p` lies to the left of the line from
/// `a` to `b`.
fn orientation(a: &Planar, b: &Planar, p: &Planar) -> QSqrt5 {
    let (bx, by) = (&b[0] - &a[0], &b[1] - &a[1]);
    let (px, py) = (&p[0] - &a[0], &p[1] - &a[1]);
    &(&bx * &py) - &(&by * &px)
}

/// The screen images of the hole's and the plug's vertices, both multiplied
/// by `1 + |r|² > 0`: the hole's are `(1 + |r|²) T_u(v)`, the plug's are
/// `T_u(R̂(r) v)` with the Cayley numerator `R̂(r) = (1 + |r|²) R(r)`.
fn shadows(centre: &Centre) -> (Vec<Planar>, Vec<Planar>) {
    let [s, t, r1, r2, r3] = centre.coordinates().clone();
    shadows_at([s, t], [r1, r2, r3])
}

/// [`shadows`] at the view `(s, t, 1)` and any rotation vector `r`.
fn shadows_at([s, t]: [QSqrt5; 2], r: [QSqrt5; 3]) -> (Vec<Planar>, Vec<Planar>) {
    let u = [s, t, QSqrt5::one()];
    let scale = &QSqrt5::one() + &dot(&r, &r);
    // rotate((1, r), v) = (1, r)(0, v)(1, r)‾ = (1 + |r|²) R(r) v = R̂(r) v.
    let [x, y, z] = r;
    let cayley = [QSqrt5::one(), x, y, z];
    let hole = vertices()
        .iter()
        .map(|v| screen(&u, v).map(|c| &c * &scale))
        .collect();
    let plug = vertices()
        .iter()
        .map(|v| screen(&u, &rotate(&cayley, v)))
        .collect();
    (hole, plug)
}

/// `Some(σ)` when every hole point `h` has `σ · orientation(a, b, h) ≤ 0`,
/// so that the line through the distinct points `a` and `b` supports the
/// hole shadow with the shadow on its nonpositive side; `None` otherwise.
fn supporting_side(hole: &[Planar], a: usize, b: usize) -> Option<i8> {
    if hole[a] == hole[b] {
        return None;
    }
    let (mut positive, mut negative) = (false, false);
    for h in hole {
        match orientation(&hole[a], &hole[b], h).sign() {
            Ordering::Greater => positive = true,
            Ordering::Less => negative = true,
            Ordering::Equal => {}
        }
        if positive && negative {
            return None;
        }
    }
    match (positive, negative) {
        (false, true) => Some(1),
        (true, false) => Some(-1),
        _ => panic!("the hole points are collinear: the hole has empty interior"),
    }
}

/// The relation of the plug shadow to the hole shadow of the configuration,
/// decided exactly.
pub fn relation(centre: &Centre) -> Relation {
    let (hole, plug) = shadows(centre);
    classify(&hole, &plug)
}

/// The relation of the convex hull of `plug` to the convex hull of `hole`,
/// which must have nonempty interior.
///
/// The edge lines of the hole polygon are exactly the lines through two
/// distinct hole points with all hole points on one closed side. The plug
/// hull lies in the closed (open) polygon iff every plug point lies on the
/// closed (open) inner side of every such line.
pub fn classify(hole: &[Planar], plug: &[Planar]) -> Relation {
    let mut edges = 0;
    let mut strict = true;
    for a in 0..hole.len() {
        for b in (a + 1)..hole.len() {
            let Some(side) = supporting_side(hole, a, b) else {
                continue;
            };
            edges += 1;
            for p in plug {
                let o = orientation(&hole[a], &hole[b], p);
                let o = if side > 0 { o } else { -&o };
                match o.sign() {
                    Ordering::Greater => return Relation::Poke,
                    Ordering::Equal => strict = false,
                    Ordering::Less => {}
                }
            }
        }
    }
    assert!(edges >= 3, "a polygon with nonempty interior has at least three edges");
    if strict {
        Relation::Fit
    } else {
        Relation::Touch
    }
}

#[cfg(test)]
mod tests;
