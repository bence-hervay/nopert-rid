//! The exact decision whether one inequality of `D` is violated throughout
//! a box: a reference for the tests
//! of this module and of the Domain component, which decides through the
//! lemma instead.
use super::super::{constraint, rotation_inequality, triangle, Constraint, DomainError, CONSTRAINT_COUNT};
use crate::arithmetic::exact::{q, Interval, QSqrt5, Q};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::geometry::{self, Point, Quaternion};
use std::cmp::Ordering;

/// `L_g(u, r) = scalar((0, u)(1, r) g) = -g⃗·u - r·(g₀u + g⃗ × u)` at the
/// affine view `u = (s, t, 1)`: its constant `-g⃗·u` and the coefficient
/// vector `-(g₀u + g⃗ × u)` of `r`, with
/// `g₀u + g⃗ × u = (g₀s + g₂ - g₃t, g₀t + g₃s - g₁, g₀ + g₁t - g₂s)`.
fn fold_scalar(g: &Quaternion, s: &Q, t: &Q) -> (QSqrt5, Point) {
    let [g0, g1, g2, g3] = g;
    let constant = -&(&(&g1.scale(s) + &g2.scale(t)) + g3);
    let coefficients = [
        &(&g0.scale(s) + g2) - &g3.scale(t),
        &(&g0.scale(t) + &g3.scale(s)) - g1,
        &(g0 + &g1.scale(t)) - &g2.scale(s),
    ];
    (constant, coefficients.map(|c| -&c))
}

/// Exactly whether inequality `index` is violated (`c > 0`) at every point of
/// the closed box `b`. `false` also covers the case where `b` meets `c ≤ 0`.
pub fn positive_on(b: &ConfigurationBox, index: usize) -> Result<bool, DomainError> {
    let form = match constraint(index)? {
        Constraint::Fold { symmetry } => {
            return Ok(fold_positive(b, &geometry::symmetries()[symmetry]))
        }
        Constraint::Triangle => triangle(),
        Constraint::Rotation {
            axis,
            first,
            second,
        } => rotation_inequality(axis, first, second),
    };
    let axes = b.axes();
    let terms = form
        .view
        .iter()
        .zip(&axes[..2])
        .chain(form.rotation.iter().zip(&axes[2..]));
    let minimum = terms.fold(form.constant.clone(), |sum, (c, axis)| {
        &sum + &affine_minimum(c, axis)
    });
    Ok(minimum.sign() == Ordering::Greater)
}

/// The first inequality violated throughout `b`, if any: then `b` misses `D`.
/// `None` does not assert that `b` meets `D`.
pub fn misses(b: &ConfigurationBox) -> Option<usize> {
    (0..CONSTRAINT_COUNT).find(|&index| positive_on(b, index).expect("index in range"))
}

/// The exact minimum of `c x` over `x ∈ axis`.
fn affine_minimum(c: &QSqrt5, axis: &Interval) -> QSqrt5 {
    c.scale(if c.sign() == Ordering::Less {
        axis.hi()
    } else {
        axis.lo()
    })
}

/// The exact range of `L_g(u, r)` over the rotation box at the view `(s, t, 1)`.
fn fold_range(b: &ConfigurationBox, g: &Quaternion, s: &Q, t: &Q) -> (QSqrt5, QSqrt5) {
    let (constant, coefficients) = fold_scalar(g, s, t);
    let mut lower = constant.clone();
    let mut upper = constant;
    for (c, axis) in coefficients.iter().zip(&b.axes()[2..]) {
        let (low, high) = if c.sign() == Ordering::Less {
            (axis.hi(), axis.lo())
        } else {
            (axis.lo(), axis.hi())
        };
        lower = &lower + &c.scale(low);
        upper = &upper + &c.scale(high);
    }
    (lower, upper)
}

/// `L² - (1 + s² + t²) > 0` on the whole box, by the four-view-corner rule
/// proved in the article: `L` keeps one nonzero sign over the box's corners and
/// the smallest `|L|` at each view corner has square above `1 + s² + t²`.
fn fold_positive(b: &ConfigurationBox, g: &Quaternion) -> bool {
    let axes = b.axes();
    let mut orientation = None;
    for s in [axes[0].lo(), axes[0].hi()] {
        for t in [axes[1].lo(), axes[1].hi()] {
            let (lower, upper) = fold_range(b, g, s, t);
            let positive = if lower.sign() == Ordering::Greater {
                true
            } else if upper.sign() == Ordering::Less {
                false
            } else {
                return false;
            };
            if *orientation.get_or_insert(positive) != positive {
                return false;
            }
            let nearest = if positive { lower } else { upper };
            let norm = QSqrt5::from_rational(q(1) + s * s + t * t);
            if (&(&nearest * &nearest) - &norm).sign() != Ordering::Greater {
                return false;
            }
        }
    }
    true
}

