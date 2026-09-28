//! Tests of coordinates, no-fit sets and zooms. `Random` is shared with the
//! tests of the other elimination modules.
use super::*;
use crate::arithmetic::exact::frac;
use crate::problem::geometry::{self, rotate, vertices};

/// splitmix64, a small deterministic generator for tests.
pub(crate) struct Random(u64);

impl Random {
    pub(crate) fn new(seed: u64) -> Self {
        Self(seed)
    }
    pub(crate) fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub(crate) fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    /// A random rational of `[lo, hi]` on the grid of step `(hi − lo)/2^bits`
    /// (ends included).
    pub(crate) fn between(&mut self, lo: &Q, hi: &Q, bits: u32) -> Q {
        let steps = 1u64 << bits;
        let k = self.below(steps + 1);
        lo + &((hi - lo) * frac(k as i64, steps as i64))
    }
    /// A random rational `n/d` with `|n| ≤ size`, `1 ≤ d ≤ size`.
    pub(crate) fn small(&mut self, size: u64) -> Q {
        let n = self.below(2 * size + 1) as i64 - size as i64;
        let d = self.below(size) as i64 + 1;
        frac(n, d)
    }
    pub(crate) fn number(&mut self, size: u64) -> QSqrt5 {
        QSqrt5::new(self.small(size), self.small(size))
    }
    /// A random point of a box, on the grid of step `width/2^bits`.
    pub(crate) fn point(&mut self, cell: &Five<Interval>, bits: u32) -> Five<Q> {
        std::array::from_fn(|j| self.between(cell[j].lo(), cell[j].hi(), bits))
    }
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn exact(x: &Five<Q>) -> Five<QSqrt5> {
    std::array::from_fn(|j| QSqrt5::from_rational(x[j].clone()))
}

/// The configuration `(u, r)` of exact coordinates `x`.
fn configuration_at(coordinates: Coordinates, x: &Five<QSqrt5>) -> (Point, Point) {
    let constants: Five<Polynomial> = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (view, rotation) = coordinates.configuration(&constants);
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    let u = view.vector().map(|p| p.evaluate(&origin));
    (u, rotation.map(|p| p.evaluate(&origin)))
}

fn signs() -> [Sign; 2] {
    [Sign::Plus, Sign::Minus]
}

#[test]
fn arc_plane_coordinates_invert_the_configuration() {
    let mut random = Random::new(1);
    for _ in 0..300 {
        let c: Five<Q> = [
            random.between(&frac(-19, 10), &q(3), 12),
            random.between(&q(-1), &q(1), 12),
            random.small(20),
            random.small(20),
            random.small(20),
        ];
        let c = exact(&c);
        for sign in signs() {
            let x = Coordinates::ArcPlane(sign).of_configuration(&c).unwrap();
            let (u, r) = configuration_at(Coordinates::ArcPlane(sign), &x);
            // u = 5/(2 + s) (s, t, 1) and r is the configuration's rotation.
            let lambda = &QSqrt5::integer(5) * &reciprocal(&(&QSqrt5::integer(2) + &c[0]));
            assert_eq!(u, [&lambda * &c[0], &lambda * &c[1], lambda.clone()]);
            assert_eq!(r, [c[2].clone(), c[3].clone(), c[4].clone()]);
            assert_eq!(x[2], c[3], "theta = r2");
        }
        let (u, r) = configuration_at(Coordinates::Configuration, &c);
        assert_eq!(u, [c[0].clone(), c[1].clone(), QSqrt5::one()]);
        assert_eq!(r, [c[2].clone(), c[3].clone(), c[4].clone()]);
    }
    let beyond: Five<QSqrt5> = std::array::from_fn(|j| QSqrt5::integer(if j == 0 { -2 } else { 0 }));
    assert!(Coordinates::ArcPlane(Sign::Plus).of_configuration(&beyond).is_none());
}

fn quaternion(scalar: QSqrt5, v: &Point) -> geometry::Quaternion {
    [scalar, v[0].clone(), v[1].clone(), v[2].clone()]
}

fn highest_second_coordinate(points: impl Iterator<Item = Point>) -> QSqrt5 {
    points.map(|p| p[1].clone()).max().unwrap()
}

#[test]
fn arc_plane_is_a_no_fit_set() {
    let a = arc_rotation();
    let norm = geometry::dot(&a, &a);
    // |r_A|² = 1 − 2/√5 = 1 − (2/5)√5, and r_A is perpendicular to e₂.
    assert_eq!(norm, QSqrt5::new(q(1), frac(-2, 5)));
    let top = highest_second_coordinate(vertices().iter().cloned());
    assert_eq!(top, QSqrt5::new(q(2), q(1)), "the hole's support value in e₂ is 2 + √5");
    let mut random = Random::new(2);
    for sign in signs() {
        let sigma = QSqrt5::integer(sign.value());
        let base: Point = a.clone().map(|c| &c * &sigma);
        let w = arc_direction(sign);
        // R(σ r_A) keeps the support value in e₂ (up to the factor 1 + |r_A|²).
        let turned = vertices().iter().map(|p| rotate(&quaternion(QSqrt5::one(), &base), p));
        let scale = &QSqrt5::one() + &norm;
        assert_eq!(highest_second_coordinate(turned), &scale * &top);
        for _ in 0..40 {
            let theta = QSqrt5::from_rational(random.small(30));
            // (1, θe₂)(1, σr_A) = (1, σr_A + θw_σ).
            let e2 = [QSqrt5::zero(), theta.clone(), QSqrt5::zero()];
            let r: Point = std::array::from_fn(|j| &base[j] + &(&theta * &w[j]));
            let product = geometry::quaternion_product(&quaternion(QSqrt5::one(), &e2), &quaternion(QSqrt5::one(), &base));
            assert_eq!(product, quaternion(QSqrt5::one(), &r));
            // Hence the plug R(r)K has the hole's support value 2 + √5 in the
            // screen direction e₂ (a view with t = 0 is perpendicular to e₂).
            let d = &QSqrt5::one() + &geometry::dot(&r, &r);
            let plug = vertices().iter().map(|p| rotate(&quaternion(QSqrt5::one(), &r), p));
            assert_eq!(highest_second_coordinate(plug), &d * &top);
        }
    }
}

#[test]
fn aligned_rotation_is_the_identity() {
    let zero: [Polynomial; 3] = std::array::from_fn(|_| Polynomial::zero());
    let m = geometry::cayley_matrix(&zero).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { QSqrt5::one() } else { QSqrt5::zero() };
            assert_eq!(m[i][j], Polynomial::constant(expected));
        }
    }
}

fn random_box(random: &mut Random, lo: &Q, hi: &Q) -> ConfigurationBox {
    ConfigurationBox::new(std::array::from_fn(|_| {
        let a = random.between(lo, hi, 6);
        let b = random.between(lo, hi, 6);
        if a <= b {
            interval(a, b)
        } else {
            interval(b, a)
        }
    }))
}

fn corners(b: &ConfigurationBox) -> Vec<Five<QSqrt5>> {
    (0..32u8).map(|mask| exact(&b.corner(mask))).collect()
}

fn combination(a: &Five<QSqrt5>, x: &Five<QSqrt5>) -> QSqrt5 {
    (0..5).fold(QSqrt5::zero(), |sum, j| &sum + &(&a[j] * &x[j]))
}

#[test]
fn range_is_the_exact_range() {
    let mut random = Random::new(3);
    for coordinates in [Coordinates::Configuration, Coordinates::ArcPlane(Sign::Plus), Coordinates::ArcPlane(Sign::Minus)] {
        for _ in 0..150 {
            let b = random_box(&mut random, &frac(-3, 2), &q(2));
            let a: Five<QSqrt5> = std::array::from_fn(|_| random.number(9));
            let [lo, hi] = coordinates.range(&a, &b).unwrap();
            let values: Vec<QSqrt5> = corners(&b)
                .iter()
                .map(|c| combination(&a, &coordinates.of_configuration(c).unwrap()))
                .collect();
            // The ends are attained at corners, and every point lies between them.
            assert_eq!(&lo, values.iter().min().unwrap());
            assert_eq!(&hi, values.iter().max().unwrap());
            let cell = b.axes().clone();
            for _ in 0..10 {
                let c = exact(&random.point(&cell, 5));
                let value = combination(&a, &coordinates.of_configuration(&c).unwrap());
                assert!(lo <= value && value <= hi);
            }
        }
    }
}

#[test]
fn range_refuses_views_at_or_beyond_minus_two() {
    let a: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::one());
    let make = |s_lo: Q| {
        ConfigurationBox::new(std::array::from_fn(|j| if j == 0 { interval(s_lo.clone(), q(0)) } else { interval(q(0), q(1)) }))
    };
    for sign in signs() {
        assert!(Coordinates::ArcPlane(sign).range(&a, &make(q(-2))).is_none());
        assert!(Coordinates::ArcPlane(sign).range(&a, &make(q(-3))).is_none());
        assert!(Coordinates::ArcPlane(sign).range(&a, &make(frac(-1999, 1000))).is_some());
    }
    assert!(Coordinates::Configuration.range(&a, &make(q(-3))).is_some());
}

fn unit_root() -> ZoomCell {
    ZoomCell(std::array::from_fn(|_| interval(q(0), q(1))))
}

#[test]
fn zoom_refusals() {
    let y = Polynomial::variables();
    let zero = Polynomial::zero;
    let configuration = |map: Five<Polynomial>, root: ZoomCell, factors: [bool; 5]| {
        Zoom::new("test".into(), Coordinates::Configuration, map, root, factors)
    };
    // r = y₀ (y₂, y₃, y₄): the factor y₀ vanishes only on the aligned set.
    let good = || -> Five<Polynomial> {
        [y[1].clone(), y[2].clone(), y[0].mul(&y[2]).unwrap(), y[0].mul(&y[3]).unwrap(), y[0].mul(&y[4]).unwrap()]
    };
    let factor0 = [true, false, false, false, false];
    assert!(configuration(good(), unit_root(), factor0).is_ok());
    let mut flat = unit_root();
    flat[3] = Interval::point(q(0));
    assert_eq!(configuration(good(), flat, factor0), Err(ZoomError::FlatRoot { variable: 3 }));
    let mut negative = unit_root();
    negative[0] = interval(q(-1), q(1));
    assert_eq!(configuration(good(), negative, factor0), Err(ZoomError::NegativeScale { variable: 0 }));
    // y₁ is not a factor: at y₁ = 0 the rotation need not vanish.
    assert_eq!(
        configuration(good(), unit_root(), [true, true, false, false, false]),
        Err(ZoomError::ZeroSet { variable: 1, coordinate: 2 })
    );
    // One stray term that survives y₀ = 0.
    let mut stray = good();
    stray[4] = &stray[4] + &y[1].mul(&y[3]).unwrap();
    assert_eq!(configuration(stray, unit_root(), factor0), Err(ZoomError::ZeroSet { variable: 0, coordinate: 4 }));
    let mut square = good();
    square[1] = y[2].mul(&y[2]).unwrap();
    assert_eq!(configuration(square, unit_root(), factor0), Err(ZoomError::ViewNotMultilinear { component: 1 }));
    // Arc-plane coordinates: the no-fit coordinates are v, ζ₁, ζ₃.
    let arc = |map: Five<Polynomial>| Zoom::new("arc".into(), Coordinates::ArcPlane(Sign::Plus), map, unit_root(), factor0);
    assert!(arc([y[1].clone(), y[0].mul(&y[2]).unwrap(), y[3].clone(), y[0].clone(), zero()]).is_ok());
    assert_eq!(
        arc([y[1].clone(), y[0].mul(&y[2]).unwrap(), y[3].clone(), y[4].clone(), zero()]),
        Err(ZoomError::ZeroSet { variable: 0, coordinate: 3 })
    );
    // θ may stay free at the zero set; e in the view must stay multilinear.
    assert_eq!(
        arc([y[1].mul(&y[1]).unwrap(), zero(), zero(), zero(), zero()]),
        Err(ZoomError::ViewNotMultilinear { component: 0 })
    );
}

#[test]
fn identity_zoom_and_view_corners() {
    let root = crate::problem::configuration::ConfigurationBox::root();
    let zoom = Zoom::identity(ZoomCell(root.axes().clone())).unwrap();
    assert_eq!(zoom.map(), &Polynomial::variables());
    assert!((0..5).all(|j| !zoom.is_scale(j)));
    let corners = zoom.view_corners(&ZoomCell(root.axes().clone()));
    assert_eq!(corners.len(), 4);
    let expected: Vec<Point> = root.view_corners().into_iter().collect();
    for u in &corners {
        assert!(expected.contains(u));
    }
    // A homogeneous zoom: e = y₀ + y₁ y₂, v = y₃ y₁; the view uses y₀..y₃.
    let y = Polynomial::variables();
    let map = [&y[0] + &y[1].mul(&y[2]).unwrap(), y[3].mul(&y[1]).unwrap(), Polynomial::zero(), Polynomial::zero(), y[0].clone()];
    let zoom = Zoom::new("arc".into(), Coordinates::ArcPlane(Sign::Minus), map, unit_root(), [false; 5]).unwrap();
    let mut random = Random::new(4);
    let cell = ZoomCell(std::array::from_fn(|_| {
        let a = random.between(&q(0), &q(1), 4);
        interval(a.clone(), &a + &frac(1, 16))
    }));
    let corners = zoom.view_corners(&cell);
    assert_eq!(corners.len(), 16);
    for mask in 0..16usize {
        let y: Five<QSqrt5> = std::array::from_fn(|j| {
            let high = j < 4 && mask >> j & 1 == 1;
            QSqrt5::from_rational(if high { cell[j].hi() } else { cell[j].lo() }.clone())
        });
        let u = zoom.view().vector().map(|p| p.evaluate(&y));
        assert!(corners.contains(&u));
    }
}

#[test]
fn reciprocal_inverts() {
    let mut random = Random::new(5);
    for _ in 0..200 {
        let x = random.number(50);
        if !x.is_zero() {
            assert_eq!(&x * &reciprocal(&x), QSqrt5::one());
        }
    }
}

/// An independent model, written from the definitions.
pub(crate) mod independent;

/// Reparametrisation invariants: both arc-plane labels, homogeneous views
/// and zoom images followed back through the coverage map.
mod reparametrisation;
