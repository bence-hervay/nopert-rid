use super::*;
use crate::arithmetic::exact::tests::{rational_power, sqrt5_enclosure};
use crate::arithmetic::exact::{Interval, Q};
use crate::arithmetic::polynomial::VARIABLES;
use crate::problem::configuration::ConfigurationBox;
use crate::problem::geometry::Quaternion;
use std::cmp::Ordering;
use crate::problem::geometry::tests::symmetry_index;
use crate::problem::geometry::{conjugate, quaternion_product, rotate, symmetries};
use num_bigint::BigInt;
use std::collections::BTreeMap;

mod adversarial;
pub(crate) mod decision;

use decision::{misses, positive_on};

type C = QSqrt5;

fn int(n: i64) -> C {
    C::integer(n)
}

fn rat(x: Q) -> C {
    C::from_rational(x)
}

fn phi() -> C {
    number(1, 1, 2)
}

/// `1/c` in Q(√5): the conjugate divided by the field norm.
fn inverse(c: &C) -> C {
    let norm = c.norm();
    assert_ne!(norm, q(0));
    c.conjugate().scale(&(q(1) / norm))
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    /// A rational in [lo, hi] with a random denominator up to 2^16.
    fn between(&mut self, lo: &Q, hi: &Q) -> Q {
        self.fraction(lo, hi, 1 << 16)
    }
    /// A rational in [lo, hi] with a random denominator up to `bound`.
    fn fraction(&mut self, lo: &Q, hi: &Q, bound: u64) -> Q {
        let denominator = 1 + self.below(bound) as i64;
        let k = self.below(denominator as u64 + 1) as i64;
        lo + (hi - lo) * frac(k, denominator)
    }
    fn rational(&mut self, bound: i64) -> Q {
        self.between(&q(-bound), &q(bound))
    }
    fn number(&mut self, bound: i64) -> C {
        C::new(self.rational(bound), self.rational(bound))
    }
}

// ---- Independent evaluation of the inequalities ------------------------

/// `scalar((0, u)(1, r) g)` by Hamilton multiplication.
fn fold_scalar_by_product(u: &Point, r: &Point, g: &Quaternion) -> C {
    let view = [int(0), u[0].clone(), u[1].clone(), u[2].clone()];
    let rotation = [int(1), r[0].clone(), r[1].clone(), r[2].clone()];
    quaternion_product(&quaternion_product(&view, &rotation), g)[0].clone()
}

/// The value of inequality `index` at the homogeneous view `u` and rotation `r`.
fn independent(index: usize, u: &Point, r: &Point) -> C {
    let phi = phi();
    if index == 0 {
        return &(&(&phi * &u[0]) + &(&(&phi + &int(1)) * &u[1])) - &u[2];
    }
    if index <= 12 {
        let k = index - 1;
        let (i, j) = (k / 4, (k / 4 + 1) % 3);
        let first = if k % 4 < 2 { -1 } else { 1 };
        let second = if k % 2 == 0 { -1 } else { 1 };
        let beta = &phi - &int(1);
        let linear = &r[i].scale(&q(first)) + &(&beta * &r[j]).scale(&q(second));
        return &(&linear + &phi) - &int(2);
    }
    let l = fold_scalar_by_product(u, r, &symmetries()[index - 13]);
    &(&l * &l) - &geometry::dot(u, u)
}

fn split(x: &[C; 5]) -> (Point, Point) {
    (
        [x[0].clone(), x[1].clone(), int(1)],
        [x[2].clone(), x[3].clone(), x[4].clone()],
    )
}

fn independent_affine(index: usize, x: &[C; 5]) -> C {
    let (u, r) = split(x);
    independent(index, &u, &r)
}

fn exact(x: &[Q; 5]) -> [C; 5] {
    x.clone().map(rat)
}

fn at(p: &Polynomial, x: &[Q; VARIABLES]) -> C {
    p.evaluate(&exact(x))
}

fn constants(p: &Point) -> [Polynomial; 3] {
    p.clone().map(Polynomial::constant)
}

fn value_of_constant(p: &Polynomial) -> C {
    p.evaluate(&std::array::from_fn(|_| int(0)))
}

// ---- Layout and polynomial forms ---------------------------------------

#[test]
fn inequality_indices_follow_the_documented_layout() {
    assert_eq!(constraint(0), Ok(Constraint::Triangle));
    let mut index = 1;
    for axis in 0..3 {
        for first in [Sign::Minus, Sign::Plus] {
            for second in [Sign::Minus, Sign::Plus] {
                assert_eq!(
                    constraint(index),
                    Ok(Constraint::Rotation {
                        axis,
                        first,
                        second
                    })
                );
                assert_eq!(view_degree(index), Ok(0));
                index += 1;
            }
        }
    }
    assert_eq!(index, FOLD_OFFSET);
    for symmetry in 0..SYMMETRY_COUNT {
        assert_eq!(fold_index(symmetry), Ok(FOLD_OFFSET + symmetry));
        assert_eq!(
            constraint(FOLD_OFFSET + symmetry),
            Ok(Constraint::Fold { symmetry })
        );
        assert_eq!(view_degree(FOLD_OFFSET + symmetry), Ok(2));
    }
    assert_eq!(view_degree(TRIANGLE), Ok(1));
    assert_eq!(FOLD_OFFSET + SYMMETRY_COUNT, CONSTRAINT_COUNT);
    let (view, rotation) = configuration::coordinates();
    for bad in [CONSTRAINT_COUNT, CONSTRAINT_COUNT + 1, usize::MAX] {
        let error = DomainError::IndexOutOfRange { index: bad };
        assert_eq!(constraint(bad), Err(error.clone()));
        assert_eq!(view_degree(bad), Err(error.clone()));
        assert_eq!(polynomial(bad, &view, &rotation), Err(error.clone()));
        assert_eq!(positive_on(&ConfigurationBox::root(), bad), Err(error));
    }
    assert_eq!(
        fold_index(SYMMETRY_COUNT),
        Err(DomainError::SymmetryOutOfRange {
            symmetry: SYMMETRY_COUNT
        })
    );
}

#[test]
fn affine_polynomials_match_independent_formulas_at_random_points() {
    let mut random = Random(3);
    let polynomials = affine_polynomials();
    assert_eq!(polynomials.len(), CONSTRAINT_COUNT);
    for _ in 0..25 {
        let x: [Q; 5] = std::array::from_fn(|_| random.rational(1));
        for (index, p) in polynomials.iter().enumerate() {
            assert_eq!(at(p, &x), independent_affine(index, &exact(&x)), "index {index}");
        }
    }
}

#[test]
fn homogeneous_forms_scale_with_the_view() {
    let mut random = Random(4);
    let x = Polynomial::variables();
    // Symbolic homogeneous view (x₀, x₁, x₂) with rotation (x₃, x₄, 1/3).
    let view = View::Homogeneous([x[0].clone(), x[1].clone(), x[2].clone()]);
    let rotation = [x[3].clone(), x[4].clone(), Polynomial::constant(rat(frac(1, 3)))];
    let symbolic: Vec<Polynomial> = (0..CONSTRAINT_COUNT)
        .map(|index| polynomial(index, &view, &rotation).unwrap())
        .collect();
    for _ in 0..12 {
        let point: [Q; 5] = std::array::from_fn(|_| random.rational(2));
        if point[2] == q(0) {
            continue;
        }
        let affine = [
            &point[0] / &point[2],
            &point[1] / &point[2],
            point[3].clone(),
            point[4].clone(),
            frac(1, 3),
        ];
        for (index, p) in symbolic.iter().enumerate() {
            let degree = view_degree(index).unwrap();
            let expected = independent_affine(index, &exact(&affine))
                .scale(&rational_power(&point[2], degree));
            assert_eq!(at(p, &point), expected, "index {index}");
        }
    }
    // Irrational configurations and scalings, through constant polynomials.
    for _ in 0..6 {
        let s = random.number(1);
        let t = random.number(1);
        let r: Point = std::array::from_fn(|_| random.number(1));
        let lambda = random.number(3);
        let u: Point = [s.clone(), t.clone(), int(1)];
        let scaled = u.clone().map(|c| &c * &lambda);
        let affine_view = View::Affine {
            s: Polynomial::constant(s.clone()),
            t: Polynomial::constant(t.clone()),
        };
        let homogeneous = View::Homogeneous(constants(&scaled));
        for index in 0..CONSTRAINT_COUNT {
            let expected = independent(index, &u, &r);
            let at_affine = polynomial(index, &affine_view, &constants(&r)).unwrap();
            assert_eq!(value_of_constant(&at_affine), expected);
            let degree = view_degree(index).unwrap();
            let factor = (0..degree).fold(int(1), |f, _| &f * &lambda);
            let at_scaled = polynomial(index, &homogeneous, &constants(&r)).unwrap();
            assert_eq!(value_of_constant(&at_scaled), &expected * &factor);
        }
    }
}

#[test]
fn polynomial_construction_refuses_exponent_overflow() {
    let x = Polynomial::variables();
    let huge = x[0].pow(200).unwrap();
    let view = View::Homogeneous([huge.clone(), x[1].clone(), x[2].clone()]);
    let rotation = [x[2].clone(), x[3].clone(), x[4].clone()];
    let index = fold_index(0).unwrap();
    assert!(matches!(
        polynomial(index, &view, &rotation),
        Err(DomainError::Polynomial(_))
    ));
    // The affine inequalities are linear and never overflow.
    assert!(polynomial(TRIANGLE, &view, &rotation).is_ok());
}

// ---- The exact decision against an exact refutation search --------------

/// A point of the closed box `b` where inequality `index` is not violated,
/// found without the decision procedure: a corner, or for a fold a zero of
/// `L` on a box edge whose end corners give `L` opposite signs.
fn refutation(b: &ConfigurationBox, index: usize) -> Option<[C; 5]> {
    let corners: Vec<[C; 5]> = (0..32).map(|mask| exact(&b.corner(mask))).collect();
    if let Some(corner) = corners
        .iter()
        .find(|x| independent_affine(index, x).sign() != Ordering::Greater)
    {
        return Some(corner.clone());
    }
    if index < FOLD_OFFSET {
        return None;
    }
    let g = &symmetries()[index - FOLD_OFFSET];
    let l = |x: &[C; 5]| {
        let (u, r) = split(x);
        fold_scalar_by_product(&u, &r, g)
    };
    let values: Vec<C> = corners.iter().map(l).collect();
    for mask in 0..32usize {
        for axis in 0..5 {
            let other = mask ^ (1 << axis);
            let (la, lb) = (&values[mask], &values[other]);
            if la.sign() == lb.sign() {
                continue;
            }
            // L is affine along one axis: its zero on this edge.
            let fraction = la * &inverse(&(la - lb));
            let mut point = corners[mask].clone();
            let step = &corners[other][axis] - &corners[mask][axis];
            point[axis] = &corners[mask][axis] + &(&step * &fraction);
            assert!(l(&point).is_zero());
            return Some(point);
        }
    }
    None
}

fn inside(b: &ConfigurationBox, x: &[C; 5]) -> bool {
    b.axes().iter().zip(x).all(|(axis, value)| {
        *value >= rat(axis.lo().clone()) && *value <= rat(axis.hi().clone())
    })
}

fn random_box(random: &mut Random, scale: &Q) -> ConfigurationBox {
    let root = ConfigurationBox::root();
    let axes = std::array::from_fn(|j| {
        let axis = &root.axes()[j];
        let slack = axis.width() / q(8);
        let lo = random.between(&(axis.lo() - &slack), &(axis.hi() + &slack));
        let width = if random.below(6) == 0 {
            q(0)
        } else {
            random.between(&q(0), &(axis.width() * scale))
        };
        interval(lo.clone(), lo + width)
    });
    ConfigurationBox::new(axes)
}

fn check_decisions(b: &ConfigurationBox, random: &mut Random, tally: &mut [[usize; 2]; 2]) {
    for index in 0..CONSTRAINT_COUNT {
        let decided = positive_on(b, index).unwrap();
        let refuted = refutation(b, index);
        if let Some(point) = &refuted {
            assert!(inside(b, point));
            assert_ne!(independent_affine(index, point).sign(), Ordering::Greater);
        }
        assert_eq!(decided, refuted.is_none(), "index {index} on {b:?}");
        tally[usize::from(index >= FOLD_OFFSET)][usize::from(decided)] += 1;
        if decided {
            for _ in 0..3 {
                let x: [Q; 5] =
                    std::array::from_fn(|j| random.between(b.axes()[j].lo(), b.axes()[j].hi()));
                assert_eq!(at(&affine_polynomials()[index], &x).sign(), Ordering::Greater);
            }
        }
    }
    let first = (0..CONSTRAINT_COUNT).find(|&i| positive_on(b, i).unwrap());
    assert_eq!(misses(b), first);
}

#[test]
fn box_decisions_agree_with_exact_refutations_on_random_boxes() {
    let mut random = Random(8);
    let mut tally = [[0usize; 2]; 2];
    for scale in [frac(1, 1), frac(1, 4), frac(1, 32), frac(1, 1024), frac(1, 1 << 30)] {
        for _ in 0..6 {
            check_decisions(&random_box(&mut random, &scale), &mut random, &mut tally);
        }
    }
    for length in [3, 9, 14, 20, 27, 33] {
        let path: String = (0..length)
            .map(|_| if random.below(2) == 0 { '0' } else { '1' })
            .collect();
        check_decisions(&ConfigurationBox::from_path(&path).unwrap(), &mut random, &mut tally);
    }
    // Both outcomes occur for affine inequalities and for folds.
    for row in tally {
        assert!(row[0] > 0 && row[1] > 0, "{tally:?}");
    }
}

// ---- Adversarial cases ------------------------------------------------

fn unit(j: usize) -> Quaternion {
    std::array::from_fn(|i| int(i64::from(i == j)))
}

fn square_fold() -> usize {
    fold_index(symmetry_index(&unit(3)).unwrap()).unwrap()
}

#[test]
fn opposite_signs_of_the_fold_scalar_never_certify() {
    // Half-turn about z at s = 1, t = 0, r = (0, y, 0): L = -1 - y. Both
    // endpoint squared gaps are 7, but L = 0 at y = -1 inside the box.
    let point = |x: i64| Interval::point(q(x));
    let cover =
        ConfigurationBox::new([point(1), point(0), point(0), interval(q(-4), q(2)), point(0)]);
    let index = square_fold();
    for mask in [0u8, 8] {
        let corner = exact(&cover.corner(mask));
        assert_eq!(independent_affine(index, &corner), int(7));
    }
    assert_eq!(positive_on(&cover, index), Ok(false));
    assert!(refutation(&cover, index).is_some());

    // The same trap across view corners: L = -2s and F = 3s² - 1.
    let g = [number(1, 0, 2), number(1, 0, 2), number(-1, 0, 2), number(1, 0, 2)];
    let index = fold_index(symmetry_index(&g).unwrap()).unwrap();
    let trap =
        ConfigurationBox::new([interval(q(-1), q(1)), point(0), point(1), point(1), point(1)]);
    for s in [-1, 1] {
        let mut axes = trap.axes().clone();
        axes[0] = point(s);
        let corner = ConfigurationBox::new(axes);
        assert_eq!(independent_affine(index, &exact(&corner.corner(0))), int(2));
        assert_eq!(positive_on(&corner, index), Ok(true));
    }
    assert_eq!(positive_on(&trap, index), Ok(false));
}

#[test]
fn touching_equalities_and_extremely_small_violations_are_distinguished() {
    let index = square_fold();
    for power in [4usize, 60, 200, 400] {
        let delta = q(1) / Q::from_integer(BigInt::from(1) << power);
        let mirror =
            ConfigurationBox::point(&[delta.clone(), delta.clone(), -&delta, delta.clone(), q(0)]);
        assert_eq!(positive_on(&mirror, index), Ok(true));
        let identity = ConfigurationBox::point(&[delta.clone(), delta, q(0), q(0), q(0)]);
        assert_eq!(misses(&identity), None);
    }
    // The square fold is (1 + s r₂ - t r₁)² ≤ 1 + s² + t²; it is tight at
    // s = 1/3, t = 1/2, r₁ = -1/3, r₂ = 0, where both sides are 49/36.
    let equality =
        ConfigurationBox::point(&[frac(1, 3), frac(1, 2), frac(-1, 3), q(0), frac(1, 5)]);
    assert_eq!(at(&affine_polynomials()[index], &equality.midpoint()), int(0));
    assert_eq!(positive_on(&equality, index), Ok(false));

    // Rational points on both sides of the irrational triangle wall, closer
    // than 2^-200 to it.
    let root = sqrt5_enclosure(220);
    let above = ConfigurationBox::point(&[(root.hi() - q(1)) / q(2), q(0), q(0), q(0), q(0)]);
    let below = ConfigurationBox::point(&[(root.lo() - q(1)) / q(2), q(0), q(0), q(0), q(0)]);
    assert_eq!(positive_on(&above, TRIANGLE), Ok(true));
    assert_eq!(positive_on(&below, TRIANGLE), Ok(false));
    assert_eq!(misses(&above), Some(TRIANGLE));
    // The triangle's corners (1/φ, 0) and (0, 1/φ²) are exactly on the wall.
    for (s, t) in [(number(-1, 1, 2), int(0)), (int(0), number(3, -1, 2))] {
        let x = [s, t, int(0), int(0), int(0)];
        assert!(independent_affine(TRIANGLE, &x).is_zero());
    }
}

#[test]
fn tiny_five_dimensional_fold_boxes_stay_strict_and_seam_crossings_refuse() {
    let index = square_fold();
    for power in [4usize, 7, 17, 31, 63, 100, 180, 250] {
        let delta = q(1) / Q::from_integer(BigInt::from(1) << power);
        let positive = interval(&delta * frac(7, 8), &delta * frac(9, 8));
        let negative = interval(-&delta * frac(9, 8), -&delta * frac(7, 8));
        let thin = interval(-&delta / q(8), &delta / q(8));
        let cover = ConfigurationBox::new([
            positive.clone(),
            positive.clone(),
            negative,
            positive,
            thin.clone(),
        ]);
        assert_eq!(positive_on(&cover, index), Ok(true), "power {power}");
        assert!(refutation(&cover, index).is_none());

        let crossing = ConfigurationBox::new([
            Interval::point(delta.clone()),
            Interval::point(delta.clone()),
            interval(-&delta, q(0)),
            interval(q(0), delta.clone()),
            thin,
        ]);
        assert_eq!(positive_on(&crossing, index), Ok(false));
        assert!(refutation(&crossing, index).is_some());
    }
}

#[test]
fn misses_reports_the_first_violated_inequality_only() {
    assert_eq!(misses(&ConfigurationBox::root()), None);
    let corner = ConfigurationBox::point(&[frac(2, 3), frac(2, 5), q(0), q(0), q(0)]);
    assert_eq!(misses(&corner), Some(TRIANGLE));
    // r₁ = 2/5 violates both rotation inequalities with first = +1 on axis 0
    // (indices 3 and 4); the first is reported.
    let rotated = ConfigurationBox::point(&[q(0), q(0), frac(2, 5), q(0), q(0)]);
    assert_eq!(misses(&rotated), Some(3));
    assert_eq!(positive_on(&rotated, 4), Ok(true));
    // Aligned configurations lie in D.
    for s in [q(0), frac(1, 5), frac(3, 5)] {
        let aligned = ConfigurationBox::point(&[s, q(0), q(0), q(0), q(0)]);
        assert_eq!(misses(&aligned), None);
    }
}

// ---- The finite facts used by the coverage proof ----------------------

#[test]
fn the_fifteen_half_turn_normals_cut_out_the_view_cone() {
    let phi = phi();
    let phi_squared = &phi + &int(1);
    let c: Point = [int(1), int(2), int(9)];
    let walls: [Point; 3] = [
        [int(1), int(0), int(0)],
        [int(0), int(1), int(0)],
        [-&phi, -&phi_squared, int(1)],
    ];
    let rays: [Point; 3] = [
        [int(0), int(0), int(1)],
        [int(1), int(0), phi.clone()],
        [int(0), int(1), phi_squared.clone()],
    ];
    let mut found = [false; 3];
    let mut normals = 0;
    for g in symmetries().iter().filter(|g| g[0].is_zero()) {
        normals += 1;
        let mut n: Point = [g[1].clone(), g[2].clone(), g[3].clone()];
        let orientation = geometry::dot(&n, &c);
        assert_ne!(orientation.sign(), Ordering::Equal);
        if orientation.sign() == Ordering::Less {
            n = n.map(|x| -&x);
        }
        for (wall, seen) in walls.iter().zip(found.iter_mut()) {
            let parallel = geometry::cross(&n, wall).iter().all(C::is_zero);
            if parallel && geometry::dot(&n, wall).sign() == Ordering::Greater {
                *seen = true;
            }
        }
        // Remark: the cone is exactly the intersection of the 15 halfspaces.
        for ray in &rays {
            assert_ne!(geometry::dot(&n, ray).sign(), Ordering::Less);
        }
    }
    assert_eq!(normals, 15);
    assert_eq!(found, [true; 3]);
    // 1/φ = φ - 1 < 2/3 and 1/φ² = 2 - φ < 2/5 put the triangle in B₀.
    assert!(&phi - &int(1) < rat(frac(2, 3)));
    assert!(&int(2) - &phi < rat(frac(2, 5)));
}

#[test]
fn the_rotation_inequalities_are_comparisons_with_group_elements() {
    // Coordinate units give a scalar of absolute value at least 1/2.
    for j in 0..4 {
        assert!(symmetry_index(&unit(j)).is_some());
    }
    let (view, rotation) = configuration::coordinates();
    for index in 1..FOLD_OFFSET {
        let Ok(Constraint::Rotation {
            axis,
            first,
            second,
        }) = constraint(index)
        else {
            panic!("index {index} is a rotation inequality");
        };
        let mut b = [phi().scale(&frac(1, 2)), int(0), int(0), int(0)];
        b[1 + axis] = rat(frac(first.value(), 2));
        b[1 + (axis + 1) % 3] = (&phi() - &int(1)).scale(&frac(second.value(), 2));
        assert!(symmetry_index(&b).is_some(), "index {index}");
        // 2 (b₀ + r·b⃗ - 1) is the inequality's polynomial.
        let comparison = rotation
            .iter()
            .zip(&b[1..])
            .fold(Polynomial::constant(&b[0] - &int(1)), |sum, (r, c)| {
                &sum + &r.scale(c)
            })
            .scale(&int(2));
        assert_eq!(comparison, polynomial(index, &view, &rotation).unwrap());
    }
    // Adding the two choices of the second sign gives |rᵢ| ≤ 2 - φ < 2/5.
    assert!(&int(2) - &phi() < rat(frac(2, 5)));
}

// ---- Executable coverage: every configuration has a representative in D ----

/// A random rational vector of rational length `scale` (stereographic).
fn random_view(random: &mut Random) -> (Point, Q) {
    let a = random.fraction(&q(-3), &q(3), 256);
    let b = random.fraction(&q(-3), &q(3), 256);
    let scale = random.fraction(&frac(1, 4), &q(4), 256);
    let denominator = q(1) + &a * &a + &b * &b;
    let u = [
        q(2) * &a / &denominator,
        q(2) * &b / &denominator,
        (q(1) - &a * &a - &b * &b) / &denominator,
    ];
    (u.map(|x| rat(x * &scale)), scale)
}

/// Invariants of the pair of shadows under isometries of the screen fixing
/// the origin, for the orthogonal projection along `u` (`|u|² = norm`) and
/// plug = rotation of `quaternion`: the multisets of projected squared norms
/// of hole and plug vertices and of hole–plug inner products.
fn shadow_invariants(u: &Point, norm: &Q, quaternion: &Quaternion) -> [BTreeMap<C, usize>; 3] {
    let vertices = geometry::vertices();
    let q_norm = quaternion.iter().fold(int(0), |s, c| &s + &(c * c));
    assert_eq!(*q_norm.sqrt5_part(), q(0));
    let divide = q(1) / q_norm.rational_part();
    let plug: Vec<Point> = vertices
        .iter()
        .map(|v| rotate(quaternion, v).map(|c| c.scale(&divide)))
        .collect();
    let along = |x: &Point| geometry::dot(x, u);
    let hole_along: Vec<C> = vertices.iter().map(along).collect();
    let plug_along: Vec<C> = plug.iter().map(along).collect();
    let inverse_norm = q(1) / norm;
    let inner = |x: &Point, ax: &C, y: &Point, ay: &C| {
        &geometry::dot(x, y) - &(ax * ay).scale(&inverse_norm)
    };
    let mut out: [BTreeMap<C, usize>; 3] = Default::default();
    let mut count = |k: usize, value: C| *out[k].entry(value).or_insert(0) += 1;
    for (h, ah) in vertices.iter().zip(&hole_along) {
        count(0, inner(h, ah, h, ah));
        for (p, ap) in plug.iter().zip(&plug_along) {
            count(2, inner(h, ah, p, ap));
        }
    }
    for (p, ap) in plug.iter().zip(&plug_along) {
        count(1, inner(p, ap, p, ap));
    }
    out
}

fn in_domain(view: &Point, r: &Point) -> bool {
    let homogeneous = View::Homogeneous(constants(view));
    (0..CONSTRAINT_COUNT).all(|index| {
        let value = value_of_constant(&polynomial(index, &homogeneous, &constants(r)).unwrap());
        value.sign() != Ordering::Greater
    })
}

/// Reduces `(u, q)` exactly as in the coverage proof and checks that the
/// representative lies in D ∩ B₀ and has congruent shadows. Returns whether
/// the rotation step was needed: the view-reduced rotation itself, where
/// Cayley-representable, violates an inequality.
fn check_representative(random: &mut Random, compare_shadows: bool) -> bool {
    let (u, scale) = random_view(random);
    let norm = &scale * &scale;
    let q0: Quaternion = std::array::from_fn(|_| rat(random.fraction(&q(-2), &q(2), 256)));
    if q0.iter().all(C::is_zero) {
        return false;
    }
    // View: maximise (1, 2, 9)·(±g u); the rotation becomes g q ḡ.
    let c: Point = [int(1), int(2), int(9)];
    let mut best: Option<(C, Point, Quaternion)> = None;
    for g in symmetries() {
        for sign in [1, -1] {
            let image = rotate(g, &u).map(|x| x.scale(&q(sign)));
            let key = geometry::dot(&c, &image);
            if best.as_ref().map_or(true, |(k, _, _)| key > *k) {
                let turned = quaternion_product(&quaternion_product(g, &q0), &conjugate(g));
                best = Some((key, image, turned));
            }
        }
    }
    let (_, view, q1) = best.unwrap();
    let phi = phi();
    assert_ne!(view[0].sign(), Ordering::Less);
    assert_ne!(view[1].sign(), Ordering::Less);
    let wall = &(&view[2] - &(&phi * &view[0])) - &(&(&phi + &int(1)) * &view[1]);
    assert_ne!(wall.sign(), Ordering::Less);
    assert_eq!(view[2].sign(), Ordering::Greater);
    // Rotation: maximise the scalar over ±q₁b and ±j q₁b, j = (0, view/|u|).
    let j = [int(0), view[0].clone(), view[1].clone(), view[2].clone()]
        .map(|x| x.scale(&(q(1) / &scale)));
    let candidates = symmetries().iter().flat_map(|b| {
        let plain = quaternion_product(&q1, b);
        let turned = quaternion_product(&j, &plain);
        [
            plain.clone(),
            plain.map(|x| -&x),
            turned.clone(),
            turned.map(|x| -&x),
        ]
    });
    let chosen = candidates.max_by(|a, b| a[0].cmp(&b[0])).unwrap();
    let needed = !q1[0].is_zero() && {
        let w = inverse(&q1[0]);
        !in_domain(&view, &[&q1[1] * &w, &q1[2] * &w, &q1[3] * &w])
    };
    assert_eq!(chosen[0].sign(), Ordering::Greater);
    let w = inverse(&chosen[0]);
    let r: Point = [&chosen[1] * &w, &chosen[2] * &w, &chosen[3] * &w];
    // All 73 inequalities hold at the representative.
    assert!(in_domain(&view, &r));
    // The representative lies in B₀.
    let lambda = inverse(&view[2]);
    let (s, t) = (&view[0] * &lambda, &view[1] * &lambda);
    for (value, axis) in [s, t, r[0].clone(), r[1].clone(), r[2].clone()]
        .iter()
        .zip(ConfigurationBox::root().axes())
    {
        assert!(*value >= rat(axis.lo().clone()) && *value <= rat(axis.hi().clone()));
    }
    // The representative's shadows are congruent to the original ones.
    if compare_shadows {
        let invariants = shadow_invariants(&u, &norm, &q0);
        assert_eq!(invariants, shadow_invariants(&view, &norm, &chosen));
        // The invariants distinguish a slightly different rotation.
        let mut other = q0.clone();
        other[1] = &other[1] + &rat(frac(1, 1000));
        assert_ne!(invariants, shadow_invariants(&u, &norm, &other));
    }
    needed
}

#[test]
fn random_configurations_reduce_to_representatives_in_the_domain() {
    let mut random = Random(12);
    let mut needed = 0;
    for sample in 0..16 {
        needed += usize::from(check_representative(&mut random, sample == 0));
    }
    assert!(needed > 0);
}

/// Heavy campaign: `cargo test --release --lib -- --ignored
/// problem::domain::tests::coverage_campaign --test-threads 2`.
#[test]
#[ignore]
fn coverage_campaign() {
    let mut random = Random(1_000_003);
    let mut needed = 0;
    for sample in 0..600 {
        needed += usize::from(check_representative(&mut random, sample % 5 == 0));
    }
    println!("rotation step needed in {needed} of 600 samples");
}

/// Heavy campaign: `cargo test --release --lib -- --ignored
/// problem::domain::tests::decision_campaign --test-threads 2`.
#[test]
#[ignore]
fn decision_campaign() {
    let mut random = Random(2_000_003);
    let mut tally = [[0usize; 2]; 2];
    for round in 0..400 {
        let scale = frac(1, 1 << (round % 24));
        check_decisions(&random_box(&mut random, &scale), &mut random, &mut tally);
    }
    println!("decisions (affine/fold × refused/positive): {tally:?}");
}

/// Costs of the exact box decision, printed:
/// `cargo test --release --lib -- --ignored problem::domain::tests::timing
/// --nocapture --test-threads 1`.
#[test]
#[ignore]
fn timing() {
    use std::time::Instant;
    let mut random = Random(99);
    let boxes: Vec<ConfigurationBox> = (0..200)
        .map(|k| {
            let path: String = (0..10 + k % 40)
                .map(|_| if random.below(2) == 0 { '0' } else { '1' })
                .collect();
            ConfigurationBox::from_path(&path).unwrap()
        })
        .collect();
    let start = Instant::now();
    let missed = boxes.iter().filter(|b| misses(b).is_some()).count();
    println!(
        "misses: {:?} per box ({missed} of {} boxes miss D)",
        start.elapsed() / boxes.len() as u32,
        boxes.len()
    );
    let start = Instant::now();
    for b in &boxes {
        for index in FOLD_OFFSET..CONSTRAINT_COUNT {
            positive_on(b, index).unwrap();
        }
    }
    println!(
        "one fold decision: {:?}",
        start.elapsed() / (boxes.len() * SYMMETRY_COUNT) as u32
    );
    let start = Instant::now();
    let (view, rotation) = configuration::coordinates();
    for index in 0..CONSTRAINT_COUNT {
        polynomial(index, &view, &rotation).unwrap();
    }
    println!("all 73 affine polynomials: {:?}", start.elapsed());
}
