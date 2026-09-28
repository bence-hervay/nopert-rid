//! The closed representative domain `D ⊂ B₀`: its 73 inequalities as exact
//! polynomials (affine and homogeneous in the view). Every configuration has
//! a representative in `D` (the article, Section 3).
use crate::arithmetic::exact::{frac, q, QSqrt5};
use crate::arithmetic::polynomial::{Polynomial, PolynomialError};
use crate::problem::configuration;
use crate::problem::geometry::{self, Point, View, SYMMETRY_COUNT};
use std::fmt;
use std::sync::OnceLock;

/// Number of inequalities `c ≤ 0` defining `D` inside `B₀`.
pub const CONSTRAINT_COUNT: usize = 73;
/// Index of the triangle inequality `φs + φ²t - 1 ≤ 0`.
pub const TRIANGLE: usize = 0;
/// Index of the first fold; fold `k` belongs to `geometry::symmetries()[k]`.
pub const FOLD_OFFSET: usize = 13;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DomainError {
    IndexOutOfRange { index: usize },
    SymmetryOutOfRange { symmetry: usize },
    /// Building the polynomial overflowed an exponent.
    Polynomial(PolynomialError),
}

impl fmt::Display for DomainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DomainError::IndexOutOfRange { index } => {
                write!(f, "domain inequality {index} is not below {CONSTRAINT_COUNT}")
            }
            DomainError::SymmetryOutOfRange { symmetry } => {
                write!(f, "symmetry {symmetry} is not below {SYMMETRY_COUNT}")
            }
            DomainError::Polynomial(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for DomainError {}

impl From<PolynomialError> for DomainError {
    fn from(error: PolynomialError) -> Self {
        DomainError::Polynomial(error)
    }
}

/// A sign `±1` in a rotation inequality.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Sign {
    Minus,
    Plus,
}

impl Sign {
    /// `-1` or `+1`.
    pub fn value(self) -> i64 {
        match self {
            Sign::Minus => -1,
            Sign::Plus => 1,
        }
    }
}

/// The meaning of one inequality index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Constraint {
    /// Index 0: `φu₁ + φ²u₂ - u₃ ≤ 0`, affinely `φs + φ²t - 1 ≤ 0`.
    Triangle,
    /// Indices 1 to 12: `first·rᵢ + second·(φ - 1)·rᵢ₊₁ + φ - 2 ≤ 0` with
    /// `i = axis` (0, 1 or 2) and indices of `r` mod 3, ordered by `axis`,
    /// then `first`, then `second`, with `Minus` before `Plus`.
    Rotation { axis: usize, first: Sign, second: Sign },
    /// Indices 13 to 72: `L_g(u, r)² - |u|² ≤ 0` for the symmetry
    /// `g = geometry::symmetries()[symmetry]`.
    Fold { symmetry: usize },
}

/// The meaning of inequality `index`.
pub fn constraint(index: usize) -> Result<Constraint, DomainError> {
    match index {
        TRIANGLE => Ok(Constraint::Triangle),
        1..=12 => {
            let k = index - 1;
            let sign = |negative: bool| if negative { Sign::Minus } else { Sign::Plus };
            Ok(Constraint::Rotation {
                axis: k / 4,
                first: sign(k % 4 < 2),
                second: sign(k % 2 == 0),
            })
        }
        FOLD_OFFSET..=72 => Ok(Constraint::Fold {
            symmetry: index - FOLD_OFFSET,
        }),
        _ => Err(DomainError::IndexOutOfRange { index }),
    }
}

/// The inequality index of the fold of symmetry `symmetry`.
pub fn fold_index(symmetry: usize) -> Result<usize, DomainError> {
    if symmetry >= SYMMETRY_COUNT {
        return Err(DomainError::SymmetryOutOfRange { symmetry });
    }
    Ok(FOLD_OFFSET + symmetry)
}

/// The degree in the view vector of the homogeneous form of inequality
/// `index`: at `u = λ(s, t, 1)` it equals `λ^degree` times the affine form.
pub fn view_degree(index: usize) -> Result<u32, DomainError> {
    Ok(match constraint(index)? {
        Constraint::Triangle => 1,
        Constraint::Rotation { .. } => 0,
        Constraint::Fold { .. } => 2,
    })
}

fn number(a: i64, b: i64, denominator: i64) -> QSqrt5 {
    QSqrt5::new(frac(a, denominator), frac(b, denominator))
}

/// An inequality affine in `(s, t, r)`: `view·(s, t) + rotation·r + constant`.
struct Affine {
    view: [QSqrt5; 2],
    rotation: [QSqrt5; 3],
    constant: QSqrt5,
}

fn triangle() -> Affine {
    Affine {
        view: [number(1, 1, 2), number(3, 1, 2)],
        rotation: [QSqrt5::zero(), QSqrt5::zero(), QSqrt5::zero()],
        constant: QSqrt5::integer(-1),
    }
}

fn rotation_inequality(axis: usize, first: Sign, second: Sign) -> Affine {
    let mut rotation = [QSqrt5::zero(), QSqrt5::zero(), QSqrt5::zero()];
    rotation[axis] = QSqrt5::integer(first.value());
    rotation[(axis + 1) % 3] = number(-1, 1, 2).scale(&q(second.value()));
    Affine {
        view: [QSqrt5::zero(), QSqrt5::zero()],
        rotation,
        constant: number(-3, 1, 2),
    }
}

/// The polynomial `c` of inequality `index` (`D` requires `c ≤ 0`) at the
/// view `view` and rotation vector `rotation`. For a homogeneous view
/// `u = λ(s, t, 1)` it is `λ^view_degree(index)` times the affine value.
pub fn polynomial(
    index: usize,
    view: &View,
    rotation: &[Polynomial; 3],
) -> Result<Polynomial, DomainError> {
    let u = view.vector();
    let affine = |form: Affine, constant: Polynomial| {
        let view_part = (0..2).fold(constant, |sum, j| &sum + &u[j].scale(&form.view[j]));
        (0..3).fold(view_part, |sum, j| &sum + &rotation[j].scale(&form.rotation[j]))
    };
    Ok(match constraint(index)? {
        Constraint::Triangle => {
            let form = triangle();
            let constant = u[2].scale(&form.constant);
            affine(form, constant)
        }
        Constraint::Rotation {
            axis,
            first,
            second,
        } => {
            let form = rotation_inequality(axis, first, second);
            let constant = Polynomial::constant(form.constant.clone());
            affine(form, constant)
        }
        Constraint::Fold { symmetry } => {
            let g = &geometry::symmetries()[symmetry];
            let vector: Point = [g[1].clone(), g[2].clone(), g[3].clone()];
            let turn = geometry::polynomial_cross(&u, &vector);
            // -L = g⃗·u + r·(g₀u + g⃗ × u), where g⃗ × u = -(u × g⃗); only L² is used.
            let mut minus_scalar = geometry::polynomial_dot(&u, &vector);
            for j in 0..3 {
                let coefficient = &u[j].scale(&g[0]) - &turn[j];
                minus_scalar = &minus_scalar + &rotation[j].mul(&coefficient)?;
            }
            &minus_scalar.mul(&minus_scalar)? - &geometry::squared_norm(&u)?
        }
    })
}

/// The 73 polynomials in the configuration variables `(s, t, r₁, r₂, r₃)`.
pub fn affine_polynomials() -> &'static [Polynomial] {
    static POLYNOMIALS: OnceLock<Vec<Polynomial>> = OnceLock::new();
    POLYNOMIALS.get_or_init(|| {
        let (view, rotation) = configuration::coordinates();
        (0..CONSTRAINT_COUNT)
            .map(|index| polynomial(index, &view, &rotation).expect("degree at most 4"))
            .collect()
    })
}

#[cfg(test)]
pub(crate) mod tests;
