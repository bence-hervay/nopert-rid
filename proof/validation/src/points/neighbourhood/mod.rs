//! The closed box of radius `2^-k` around an exact centre, intersected with
//! the root box, with rational ends that enclose the exact ends outwards by
//! at most `2^-k / 1024`; and the exact check of that property.
use crate::points::catalogue::Centre;
use num_bigint::BigInt;
use rid::arithmetic::exact::{Interval, QSqrt5, Q};
use rid::problem::configuration::{ConfigurationBox, AXES};
use std::fmt;

/// The largest radius exponent: radius `2^-4096`.
pub const MAX_EXPONENT: u32 = 4096;

/// Irrational ends are rounded outwards to multiples of `2^-(k + 10)`, so
/// the padding beyond the exact neighbourhood is at most `2^-k / 1024`.
const PADDING_BITS: u32 = 10;

/// The radius `2^-k`.
pub fn radius(k: u32) -> Q {
    Q::new(BigInt::from(1), BigInt::from(1) << k)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Lower,
    Upper,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NeighbourhoodError {
    /// The exponent exceeds [`MAX_EXPONENT`].
    Exponent(u32),
    /// The box leaves the root box or is flat on `axis`.
    Shape { axis: usize },
    /// The box omits part of the exact neighbourhood on this side of `axis`.
    Misses { axis: usize, side: Side },
    /// The box reaches more than `2^-k / 1024` beyond the exact neighbourhood.
    Padding { axis: usize, side: Side },
}

impl fmt::Display for NeighbourhoodError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NeighbourhoodError::Exponent(k) => {
                write!(f, "radius exponent {k} exceeds {MAX_EXPONENT}")
            }
            NeighbourhoodError::Shape { axis } => {
                write!(f, "axis {axis} is flat or leaves the root box")
            }
            NeighbourhoodError::Misses { axis, side } => {
                write!(f, "the {side:?} end of axis {axis} omits part of the neighbourhood")
            }
            NeighbourhoodError::Padding { axis, side } => {
                write!(f, "the {side:?} end of axis {axis} is padded by more than radius/1024")
            }
        }
    }
}

impl std::error::Error for NeighbourhoodError {}

/// `⌊x⌋` for `x` in Q(√5). The guess from `⌊a⌋` and `⌊|b|√5⌋ = isqrt(⌊5b²⌋)`
/// is within two of the answer; exact comparisons alone decide it.
fn floor(x: &QSqrt5) -> BigInt {
    let rational_floor = |q: &Q| {
        let (quotient, remainder) = (q.numer() / q.denom(), q.numer() % q.denom());
        if remainder < BigInt::from(0) {
            quotient - 1
        } else {
            quotient
        }
    };
    let (a, b) = (x.rational_part(), x.sqrt5_part());
    let root = rational_floor(&(Q::from_integer(BigInt::from(5)) * b * b)).sqrt();
    let mut n = if *b >= Q::zero() {
        rational_floor(a) + root
    } else {
        rational_floor(a) - root - 1
    };
    let at = |n: &BigInt| QSqrt5::from_rational(Q::from_integer(n.clone()));
    while at(&n) > *x {
        n -= 1;
    }
    while at(&(&n + 1)) <= *x {
        n += 1;
    }
    n
}

/// The largest multiple of `unit` at most `x`, or the smallest at least `x`.
fn round(x: &QSqrt5, unit: &Q, side: Side) -> Q {
    let scaled = x.scale(&(Q::one() / unit));
    let n = match side {
        Side::Lower => floor(&scaled),
        Side::Upper => -floor(&-&scaled),
    };
    Q::from_integer(n) * unit
}

/// The closed box `[c - 2^-k, c + 2^-k] ∩ B₀` around the centre, axis by
/// axis. Rational ends are exact; an irrational end is rounded outwards to a
/// multiple of `2^-(k+10)`.
pub fn neighbourhood(centre: &Centre, k: u32) -> Result<ConfigurationBox, NeighbourhoodError> {
    if k > MAX_EXPONENT {
        return Err(NeighbourhoodError::Exponent(k));
    }
    let r = radius(k);
    let unit = radius(k + PADDING_BITS);
    let root = ConfigurationBox::root();
    let mut axes = Vec::with_capacity(AXES);
    for (c, bounds) in centre.coordinates().iter().zip(root.axes()) {
        let (lo, hi) = if c.sqrt5_part().is_zero() {
            (c.rational_part() - &r, c.rational_part() + &r)
        } else {
            let shift = QSqrt5::from_rational(r.clone());
            (
                round(&(c - &shift), &unit, Side::Lower),
                round(&(c + &shift), &unit, Side::Upper),
            )
        };
        let lo = lo.max(bounds.lo().clone());
        let hi = hi.min(bounds.hi().clone());
        axes.push(Interval::new(lo, hi).expect("a centre in the root box has ordered ends"));
    }
    let b = ConfigurationBox::new(axes.try_into().expect("AXES intervals"));
    // The independent exact check, always (a few exact comparisons per box,
    // negligible next to the components' work on it).
    check(centre, k, &b)?;
    Ok(b)
}

/// Exactly whether `b` is a neighbourhood of radius `2^-k` of the centre as
/// [`neighbourhood`] promises, whatever its construction: on every axis,
/// with root bounds `[L, H]`, centre `c` and radius `ρ`,
/// `L ≤ lo < hi ≤ H`, `c - ρ(1 + 1/1024) ≤ lo ≤ max(L, c - ρ)` and
/// `min(H, c + ρ) ≤ hi ≤ c + ρ(1 + 1/1024)`.
pub fn check(centre: &Centre, k: u32, b: &ConfigurationBox) -> Result<(), NeighbourhoodError> {
    if k > MAX_EXPONENT {
        return Err(NeighbourhoodError::Exponent(k));
    }
    let r = QSqrt5::from_rational(radius(k));
    let reach = QSqrt5::from_rational(radius(k) + radius(k + PADDING_BITS));
    let root = ConfigurationBox::root();
    let exact = |q: &Q| QSqrt5::from_rational(q.clone());
    for (axis, ((c, bounds), interval)) in centre
        .coordinates()
        .iter()
        .zip(root.axes())
        .zip(b.axes())
        .enumerate()
    {
        let (lo, hi) = (exact(interval.lo()), exact(interval.hi()));
        let (bottom, top) = (exact(bounds.lo()), exact(bounds.hi()));
        if lo < bottom || hi > top || lo >= hi {
            return Err(NeighbourhoodError::Shape { axis });
        }
        if lo > (c - &r).max(bottom) {
            return Err(NeighbourhoodError::Misses { axis, side: Side::Lower });
        }
        if hi < (c + &r).min(top) {
            return Err(NeighbourhoodError::Misses { axis, side: Side::Upper });
        }
        if lo < c - &reach {
            return Err(NeighbourhoodError::Padding { axis, side: Side::Lower });
        }
        if hi > c + &reach {
            return Err(NeighbourhoodError::Padding { axis, side: Side::Upper });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
