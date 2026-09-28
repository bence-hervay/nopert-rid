//! Zooms: polynomial maps from closed boxes of zoom variables into
//! configurations, with designated scale variables whose zero sets lie in a
//! no-fit set. Zoom families and their coverage statements are in `cover`,
//! cell trees and their tiling in `tree`.
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::{Exponents, Polynomial, PolynomialError, VARIABLES};
use crate::problem::configuration::{ConfigurationBox, R, S, T};
use crate::problem::geometry::{Point, View};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt;

pub mod cover;
pub mod tree;

/// A point of five zoom variables, or of five coordinates.
pub type Five<T> = [T; VARIABLES];

/// A zoom cell: a closed box of the five zoom variables, such as a zoom's
/// root box or a cell of its tree (each [`Interval`] has `lo ≤ hi`). It
/// reads and writes as its five intervals.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ZoomCell(pub Five<Interval>);

impl ZoomCell {
    pub fn new(intervals: Five<Interval>) -> Self {
        Self(intervals)
    }
    pub fn intervals(&self) -> &Five<Interval> {
        &self.0
    }
    /// The two closed halves at the exact midpoint of `variable`, lower first.
    pub fn bisect(&self, variable: usize) -> [ZoomCell; 2] {
        let [lower, upper] = self.0[variable].bisect();
        let (mut low, mut high) = (self.clone(), self.clone());
        low.0[variable] = lower;
        high.0[variable] = upper;
        [low, high]
    }
    /// Whether `other` lies inside this closed cell.
    pub fn contains_cell(&self, other: &ZoomCell) -> bool {
        (0..VARIABLES).all(|j| self.0[j].lo() <= other.0[j].lo() && other.0[j].hi() <= self.0[j].hi())
    }
}

impl std::ops::Deref for ZoomCell {
    type Target = Five<Interval>;
    fn deref(&self) -> &Five<Interval> {
        &self.0
    }
}

impl std::ops::DerefMut for ZoomCell {
    fn deref_mut(&mut self) -> &mut Five<Interval> {
        &mut self.0
    }
}

impl From<Five<Interval>> for ZoomCell {
    fn from(intervals: Five<Interval>) -> Self {
        Self(intervals)
    }
}

/// The zoom factor `y^a` of a cell: a monomial in the zoom variables, by
/// which the lemma divides a pulled-back witness.
///
/// **Invariant of use:** a factor goes with a zoom only if its positive
/// exponents are on scale variables of that zoom. [`Zoom::factor`] makes
/// factors with this check, and the lemma checks it again for the zoom it
/// is given (`CellRefusal::FactorOutsideScales`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ZoomFactor(Exponents);

impl ZoomFactor {
    /// The factor 1 (every exponent 0), which goes with every zoom.
    pub const ONE: ZoomFactor = ZoomFactor([0; VARIABLES]);
    /// The monomial `y^a` for the exponent vector `a`, unchecked.
    pub fn new(exponents: Exponents) -> Self {
        Self(exponents)
    }
    /// The exponent vector `a`.
    pub fn exponents(&self) -> &Exponents {
        &self.0
    }
}

/// The sign `σ` that labels one of the two arc planes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Sign {
    Plus,
    Minus,
}

impl Sign {
    pub fn value(self) -> i64 {
        match self {
            Sign::Plus => 1,
            Sign::Minus => -1,
        }
    }
}

/// A coordinate system `x = (x₀, …, x₄)` on configurations, with the
/// three coordinates that vanish exactly on its no-fit set.
///
/// - `Configuration`: `x = (s, t, r₁, r₂, r₃)`, the view `(s, t, 1)`; the
///   no-fit set is the aligned set `r = 0`.
/// - `ArcPlane(σ)`: `x = (e, v, θ, ζ₁, ζ₃)` with the homogeneous view
///   `u = (1 − 2e, v, 2 + e)` and `r = σ r_A + θ w_σ + (ζ₁, 0, ζ₃)`; the
///   no-fit set is the arc plane `Π_σ = {v = ζ₁ = ζ₃ = 0}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Coordinates {
    Configuration,
    ArcPlane(Sign),
}

/// `r_A = (a₁, 0, a₃)` with `a₁ = (−5 + 3√5)/10`, `a₃ = (−5 + √5)/10`.
pub fn arc_rotation() -> Point {
    [
        QSqrt5::new(frac(-1, 2), frac(3, 10)),
        QSqrt5::zero(),
        QSqrt5::new(frac(-1, 2), frac(1, 10)),
    ]
}

/// `w_σ = (σa₃, 1, −σa₁)`, so that `(1, θe₂)(1, σr_A) = (1, σr_A + θw_σ)`.
pub fn arc_direction(sign: Sign) -> Point {
    let [a1, _, a3] = arc_rotation();
    let sigma = q(sign.value());
    [a3.scale(&sigma), QSqrt5::one(), a1.scale(&-sigma)]
}

impl Coordinates {
    /// The coordinates that vanish on the no-fit set (and nowhere else).
    pub fn no_fit_coordinates(self) -> [usize; 3] {
        match self {
            Coordinates::Configuration => R,
            Coordinates::ArcPlane(_) => [1, 3, 4],
        }
    }

    /// The view map and rotation vector of the configuration with
    /// coordinates `x` (polynomials in any variables).
    pub fn configuration(self, x: &Five<Polynomial>) -> (View, [Polynomial; 3]) {
        match self {
            Coordinates::Configuration => (
                View::Affine {
                    s: x[S].clone(),
                    t: x[T].clone(),
                },
                R.map(|j| x[j].clone()),
            ),
            Coordinates::ArcPlane(sign) => {
                let constant = |c: QSqrt5| Polynomial::constant(c);
                let two = QSqrt5::integer(2);
                let view = [
                    &constant(QSqrt5::one()) - &x[0].scale(&two),
                    x[1].clone(),
                    &constant(two) + &x[0],
                ];
                let (a, w) = (arc_rotation(), arc_direction(sign));
                let sigma = q(sign.value());
                let rotation = [
                    &(&constant(a[0].scale(&sigma)) + &x[2].scale(&w[0])) + &x[3],
                    x[2].scale(&w[1]),
                    &(&constant(a[2].scale(&sigma)) + &x[2].scale(&w[2])) + &x[4],
                ];
                (View::Homogeneous(view), rotation)
            }
        }
    }

    /// The coordinates of the configuration `c = (s, t, r)`, or `None` for an
    /// arc plane when `s ≤ −2` (there the view has no arc-plane coordinates).
    pub fn of_configuration(self, c: &Five<QSqrt5>) -> Option<Five<QSqrt5>> {
        match self {
            Coordinates::Configuration => Some(c.clone()),
            Coordinates::ArcPlane(sign) => {
                let scale = &QSqrt5::integer(2) + &c[S];
                if scale.sign() != Ordering::Greater {
                    return None;
                }
                let inverse = reciprocal(&scale);
                let e = &(&QSqrt5::one() - &c[S].scale(&q(2))) * &inverse;
                let v = &c[T].scale(&q(5)) * &inverse;
                let (a, sigma, theta) = (arc_rotation(), q(sign.value()), c[R[1]].clone());
                let zeta1 = &c[R[0]] - &(&a[0] + &(&a[2] * &theta)).scale(&sigma);
                let zeta3 = &c[R[2]] - &(&a[2] - &(&a[0] * &theta)).scale(&sigma);
                Some([e, v, theta, zeta1, zeta3])
            }
        }
    }

    /// The exact range of `Σⱼ a[j]·xⱼ` over the configurations of the box `b`,
    /// or `None` for an arc plane when the box reaches `s ≤ −2`; see
    /// [`LinearForm::range`].
    pub fn range(self, a: &Five<QSqrt5>, b: &ConfigurationBox) -> Option<[QSqrt5; 2]> {
        self.linear_form(a).range(b)
    }

    /// The coordinates `(x₀, x₁)` at the four view corners `(s, t)` of the
    /// box, or `None` for an arc plane when the box reaches `s ≤ −2`.
    pub fn view_corners(self, b: &ConfigurationBox) -> Option<[[QSqrt5; 2]; 4]> {
        let axes = b.axes();
        let mut corners: [[QSqrt5; 2]; 4] = std::array::from_fn(|_| [QSqrt5::zero(), QSqrt5::zero()]);
        for (k, corner) in corners.iter_mut().enumerate() {
            let mut c: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
            c[S] = QSqrt5::from_rational(if k & 1 == 0 { axes[S].lo() } else { axes[S].hi() }.clone());
            c[T] = QSqrt5::from_rational(if k & 2 == 0 { axes[T].lo() } else { axes[T].hi() }.clone());
            let x = self.of_configuration(&c)?;
            *corner = [x[0].clone(), x[1].clone()];
        }
        Some(corners)
    }

    /// `Σⱼ a[j]·xⱼ` prepared for ranges over many boxes: the rotation part
    /// folded into a constant and one coefficient per rotation variable.
    pub fn linear_form(self, a: &Five<QSqrt5>) -> LinearForm {
        // x₂..x₄ = constant + Σ_k m[j][k] r_k; collect the coefficient of each r_k.
        let (constant, matrix) = self.rotation_coordinates();
        let mut folded = QSqrt5::zero();
        let mut rotation: [QSqrt5; 3] = std::array::from_fn(|_| QSqrt5::zero());
        for j in 0..3 {
            folded = &folded + &(&a[2 + j] * &constant[j]);
            for k in 0..3 {
                rotation[k] = &rotation[k] + &(&a[2 + j] * &matrix[j][k]);
            }
        }
        LinearForm { coordinates: self, view: [a[0].clone(), a[1].clone()], constant: folded, rotation }
    }

    /// `(x₂, x₃, x₄) = constant + matrix · r` as exact data.
    fn rotation_coordinates(self) -> ([QSqrt5; 3], [[QSqrt5; 3]; 3]) {
        let unit = |i: usize, j: usize| {
            if i == j {
                QSqrt5::one()
            } else {
                QSqrt5::zero()
            }
        };
        match self {
            Coordinates::Configuration => (
                std::array::from_fn(|_| QSqrt5::zero()),
                std::array::from_fn(|i| std::array::from_fn(|j| unit(i, j))),
            ),
            Coordinates::ArcPlane(sign) => {
                // θ = r₂, ζ₁ = r₁ − σa₁ − σa₃ r₂, ζ₃ = r₃ − σa₃ + σa₁ r₂.
                let (a, sigma) = (arc_rotation(), q(sign.value()));
                let zero = QSqrt5::zero;
                (
                    [zero(), a[0].scale(&-sigma.clone()), a[2].scale(&-sigma.clone())],
                    [
                        [zero(), QSqrt5::one(), zero()],
                        [QSqrt5::one(), a[2].scale(&-sigma.clone()), zero()],
                        [zero(), a[0].scale(&sigma), QSqrt5::one()],
                    ],
                )
            }
        }
    }
}

/// A linear form `Σⱼ aⱼxⱼ` in the coordinates of a coordinate system, with
/// its rotation part folded: `a₀x₀ + a₁x₁ + constant + Σ_k rotation[k]·r_k`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearForm {
    coordinates: Coordinates,
    view: [QSqrt5; 2],
    constant: QSqrt5,
    rotation: [QSqrt5; 3],
}

impl LinearForm {
    /// The exact range over the configurations of the box `b`, or `None` for
    /// an arc plane when the box reaches `s ≤ −2`.
    ///
    /// `x₀, x₁` depend on `(s, t)` only and `x₂, x₃, x₄` are affine in `r`
    /// only, so the range is the sum of two ranges. The view part is monotone
    /// in `s` for fixed `t` and affine in `t` for fixed `s`, so its
    /// extremes are among the four `(s, t)` corners; the rotation part is
    /// affine in the independent `r₁, r₂, r₃`, a sum of exact term ranges.
    pub fn range(&self, b: &ConfigurationBox) -> Option<[QSqrt5; 2]> {
        let axes = b.axes();
        match self.coordinates {
            // Affine in (s, t): the corners' values need no field elements.
            Coordinates::Configuration => {
                let term = |a: &QSqrt5, x: &Q| if a.is_zero() { QSqrt5::zero() } else { a.scale(x) };
                let mut values = Vec::with_capacity(4);
                for s in [axes[S].lo(), axes[S].hi()] {
                    for t in [axes[T].lo(), axes[T].hi()] {
                        values.push(&term(&self.view[0], s) + &term(&self.view[1], t));
                    }
                }
                Some(self.finish(&values, b))
            }
            Coordinates::ArcPlane(_) => Some(self.range_at(&self.coordinates.view_corners(b)?, b)),
        }
    }

    /// [`LinearForm::range`] with the box's view corners in these
    /// coordinates (`Coordinates::view_corners`) computed once for several
    /// forms.
    pub fn range_at(&self, corners: &[[QSqrt5; 2]; 4], b: &ConfigurationBox) -> [QSqrt5; 2] {
        let part = |a: &QSqrt5, x: &QSqrt5| if a.is_zero() { QSqrt5::zero() } else { a * x };
        let values = corners.clone().map(|[x0, x1]| &part(&self.view[0], &x0) + &part(&self.view[1], &x1));
        self.finish(&values, b)
    }

    /// The range from the view part's values at the four corners.
    fn finish(&self, values: &[QSqrt5], b: &ConfigurationBox) -> [QSqrt5; 2] {
        let mut lo = &values.iter().min().expect("four corners").clone() + &self.constant;
        let mut hi = &values.iter().max().expect("four corners").clone() + &self.constant;
        for (k, coefficient) in self.rotation.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
            let axis = &b.axes()[R[k]];
            let ends = [coefficient.scale(axis.lo()), coefficient.scale(axis.hi())];
            let (min, max) = if ends[0] <= ends[1] { (0, 1) } else { (1, 0) };
            lo = &lo + &ends[min];
            hi = &hi + &ends[max];
        }
        [lo, hi]
    }
}

/// `1/x` for a nonzero `x` of Q(√5): `x̄ / N(x)` with the rational norm
/// `N(x) = x x̄ ≠ 0`.
pub fn reciprocal(x: &QSqrt5) -> QSqrt5 {
    assert!(!x.is_zero(), "reciprocal of zero");
    x.conjugate().scale(&(Q::one() / x.norm()))
}

/// Why a zoom was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ZoomError {
    /// The root box has zero width in this variable.
    FlatRoot { variable: usize },
    /// A scale variable takes negative values on the root box.
    NegativeScale { variable: usize },
    /// Where the scale variable vanishes, this no-fit coordinate does not
    /// vanish identically.
    ZeroSet { variable: usize, coordinate: usize },
    /// A view component has a term of degree above 1 in some variable.
    ViewNotMultilinear { component: usize },
    /// A zoom factor with a positive exponent on a variable that is not a
    /// scale variable.
    NotScale { variable: usize },
    Polynomial(PolynomialError),
}

impl fmt::Display for ZoomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ZoomError::FlatRoot { variable } => write!(f, "root box is flat in variable {variable}"),
            ZoomError::NegativeScale { variable } => {
                write!(f, "scale variable {variable} is negative on the root box")
            }
            ZoomError::ZeroSet { variable, coordinate } => write!(
                f,
                "where scale variable {variable} vanishes, coordinate {coordinate} does not"
            ),
            ZoomError::ViewNotMultilinear { component } => {
                write!(f, "view component {component} is not multilinear")
            }
            ZoomError::NotScale { variable } => write!(f, "variable {variable} is not a scale variable"),
            ZoomError::Polynomial(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ZoomError {}

/// A zoom: the polynomial map `y ↦ x(y)` from zoom variables to
/// coordinates, the configuration `(view, rotation)` it defines, a closed
/// root box and the designated scale variables.
///
/// **Invariants** (checked by [`Zoom::new`]): the root box has positive
/// widths; every scale variable is nonnegative on it; where a factor
/// variable vanishes, the three no-fit coordinates vanish identically (the
/// zero-set rule), so the configuration lies in the no-fit set; every view
/// component is multilinear in the zoom variables.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Zoom {
    name: String,
    coordinates: Coordinates,
    map: Five<Polynomial>,
    view: View,
    rotation: [Polynomial; 3],
    root: ZoomCell,
    scales: [bool; VARIABLES],
}

impl Zoom {
    pub fn new(
        name: String,
        coordinates: Coordinates,
        map: Five<Polynomial>,
        root: ZoomCell,
        scales: [bool; VARIABLES],
    ) -> Result<Self, ZoomError> {
        for (variable, interval) in root.iter().enumerate() {
            if interval.lo() >= interval.hi() {
                return Err(ZoomError::FlatRoot { variable });
            }
            if scales[variable] && interval.lo() < &Q::zero() {
                return Err(ZoomError::NegativeScale { variable });
            }
        }
        for variable in (0..VARIABLES).filter(|&j| scales[j]) {
            for coordinate in coordinates.no_fit_coordinates() {
                if map[coordinate].terms().any(|(e, _)| e[variable] == 0) {
                    return Err(ZoomError::ZeroSet {
                        variable,
                        coordinate,
                    });
                }
            }
        }
        let (view, rotation) = coordinates.configuration(&map);
        for (component, u) in view.vector().iter().enumerate() {
            if u.terms().any(|(e, _)| e.iter().any(|&k| k > 1)) {
                return Err(ZoomError::ViewNotMultilinear { component });
            }
        }
        Ok(Self {
            name,
            coordinates,
            map,
            view,
            rotation,
            root,
            scales,
        })
    }

    /// The identity zoom `x = y` in configuration coordinates on the root
    /// box `root`, without scale variables: a search box is one cell of it.
    pub fn identity(root: ZoomCell) -> Result<Self, ZoomError> {
        Self::new(
            "identity".into(),
            Coordinates::Configuration,
            Polynomial::variables(),
            root,
            [false; VARIABLES],
        )
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn coordinates(&self) -> Coordinates {
        self.coordinates
    }
    /// The coordinates `x(y)` as polynomials in the zoom variables.
    pub fn map(&self) -> &Five<Polynomial> {
        &self.map
    }
    pub fn view(&self) -> &View {
        &self.view
    }
    pub fn rotation(&self) -> &[Polynomial; 3] {
        &self.rotation
    }
    pub fn root(&self) -> &ZoomCell {
        &self.root
    }
    pub fn is_scale(&self, variable: usize) -> bool {
        self.scales[variable]
    }

    /// The zoom factor `y^a` with these exponents, refused unless every
    /// positive exponent is on a scale variable of this zoom.
    pub fn factor(&self, exponents: Exponents) -> Result<ZoomFactor, ZoomError> {
        match (0..VARIABLES).find(|&j| exponents[j] > 0 && !self.scales[j]) {
            Some(variable) => Err(ZoomError::NotScale { variable }),
            None => Ok(ZoomFactor(exponents)),
        }
    }

    /// The exact views `u(y)` at the corners of `cell` over the variables the
    /// view depends on (at most 32; one view when it is constant).
    pub fn view_corners(&self, cell: &ZoomCell) -> Vec<Point> {
        let u = self.view.vector();
        let used: Vec<usize> = (0..VARIABLES)
            .filter(|&j| u.iter().any(|p| p.degrees()[j] > 0))
            .collect();
        (0..1usize << used.len())
            .map(|mask| {
                let y: Five<QSqrt5> = std::array::from_fn(|j| {
                    let high = used.iter().position(|&v| v == j).is_some_and(|bit| mask >> bit & 1 == 1);
                    QSqrt5::from_rational(if high { cell[j].hi() } else { cell[j].lo() }.clone())
                });
                std::array::from_fn(|i| u[i].evaluate(&y))
            })
            .collect()
    }

    /// The coordinates `x(y)` at an exact point of zoom variables.
    pub fn evaluate(&self, y: &Five<QSqrt5>) -> Five<QSqrt5> {
        std::array::from_fn(|i| self.map[i].evaluate(y))
    }
}

#[cfg(test)]
pub(crate) mod tests;
