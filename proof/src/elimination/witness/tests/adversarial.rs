//! Independent adversarial tests of support gaps and their validity. Every
//! reference computation here is reimplemented from the definitions (a
//! Rodrigues rotation, its own vector arithmetic and a brute-force validity
//! check), apart from comparing the tests' own `value_at` with them.
use super::value_at;
use crate::arithmetic::exact::{frac, q, QSqrt5, Q};
use crate::arithmetic::polynomial::{Polynomial, PolynomialError};
use crate::elimination::witness::{
    Direction, DomainInequality, Edge, Gap, Support, Witness, WitnessError,
};
use crate::problem::configuration::{self, ConfigurationBox};
use crate::problem::domain::{self, CONSTRAINT_COUNT, FOLD_OFFSET};
use crate::problem::geometry::{self, Point, Quaternion, View, EDGE_COUNT, VERTEX_COUNT};
use std::cmp::Ordering;

type C = QSqrt5;

fn int(n: i64) -> C {
    C::integer(n)
}

fn rat(x: Q) -> C {
    C::from_rational(x)
}

fn inverse(c: &C) -> C {
    let norm = c.norm();
    assert_ne!(norm, q(0));
    c.conjugate().scale(&(q(1) / norm))
}

// ---- Own vector arithmetic (not the geometry module's) -----------------

fn my_dot(a: &Point, b: &Point) -> C {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn my_cross(a: &Point, b: &Point) -> Point {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn my_sub(a: &Point, b: &Point) -> Point {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

fn my_scale(a: &Point, k: &C) -> Point {
    [&a[0] * k, &a[1] * k, &a[2] * k]
}

/// `R(r)x` by the Rodrigues form of the Cayley rotation,
/// `x + 2/(1+|r|²) (r × x + r × (r × x))`.
fn cayley_rotate(r: &Point, x: &Point) -> Point {
    let k = inverse(&(&int(1) + &my_dot(r, r))).scale(&q(2));
    let once = my_cross(r, x);
    let twice = my_cross(r, &once);
    let turn = [&once[0] + &twice[0], &once[1] + &twice[1], &once[2] + &twice[2]];
    [
        &x[0] + &(&k * &turn[0]),
        &x[1] + &(&k * &turn[1]),
        &x[2] + &(&k * &turn[2]),
    ]
}

/// The support normal from the definitions, independently of `witness`.
fn my_normal(support: &Support, u: &Point) -> Point {
    let v = geometry::vertices();
    match support {
        Support::Edge(edge) => my_cross(u, &my_sub(&v[edge.to()], &v[edge.from()])),
        Support::Direction(direction) => {
            let [n1, n2] = direction.screen_normal();
            [&u[2] * n1, &u[2] * n2, -&(&(&u[0] * n1) + &(&u[1] * n2))]
        }
    }
}

/// `(1 + |r|²)(n·h - n·R(r)p)` from the definitions.
fn my_gap(gap: &Gap, u: &Point, r: &Point) -> C {
    let v = geometry::vertices();
    let n = my_normal(gap.support(), u);
    let hole = my_dot(&n, &v[gap.support().contact()]);
    let plug = my_dot(&n, &cayley_rotate(r, &v[gap.plug()]));
    &(&int(1) + &my_dot(r, r)) * &(&hole - &plug)
}

/// Brute-force validity: the contact maximises `n·w` over the vertices, with
/// `n` from the definitions (not the witness module), and `n ≠ 0`.
fn my_valid(support: &Support, u: &Point) -> bool {
    let n = my_normal(support, u);
    if n.iter().all(C::is_zero) {
        return false;
    }
    let v = geometry::vertices();
    let h = my_dot(&n, &v[support.contact()]);
    v.iter().all(|w| my_dot(&n, w) <= h)
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        // xorshift64*, a different generator from the module tests.
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D) >> 8
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn between(&mut self, lo: &Q, hi: &Q) -> Q {
        let d = 1 + self.below(1 << 24) as i64;
        let k = self.below(d as u64 + 1) as i64;
        lo + (hi - lo) * frac(k, d)
    }
    fn rational(&mut self, bound: i64) -> Q {
        self.between(&q(-bound), &q(bound))
    }
    fn number(&mut self, bound: i64) -> C {
        let a = self.rational(bound);
        match self.below(3) {
            0 => rat(a),
            _ => C::new(a / q(2), self.rational(bound) / q(2)),
        }
    }
    fn unit_fraction(&mut self) -> C {
        // A number in [0, 1], irrational two thirds of the time.
        match self.below(3) {
            0 => rat(self.between(&q(0), &q(1))),
            1 => {
                // (√5 - 1)/2 · k/d lies in [0, 0.62].
                let k = self.between(&q(0), &q(1));
                C::new(-&k / q(2), k / q(2))
            }
            _ => {
                // 1 - (√5 - 2)·k/d lies in [0.76, 1].
                let k = self.between(&q(0), &q(1));
                C::new(q(1) + &k * q(2), -k)
            }
        }
    }
    fn path(&mut self, length: usize) -> String {
        (0..length).map(|_| if self.below(2) == 0 { '0' } else { '1' }).collect()
    }
    fn support(&mut self) -> Support {
        if self.below(3) == 0 {
            loop {
                let n = [self.number(3), self.number(3)];
                if let Ok(d) = Direction::new(n, self.below(VERTEX_COUNT as u64) as usize) {
                    return Support::Direction(d);
                }
            }
        }
        let [a, b] = geometry::edges()[self.below(EDGE_COUNT as u64) as usize];
        let (from, to) = if self.below(2) == 0 { (a, b) } else { (b, a) };
        Support::Edge(Edge::new(from, to).unwrap())
    }
}

fn oriented_edges() -> Vec<Support> {
    geometry::edges()
        .iter()
        .flat_map(|[a, b]| [(*a, *b), (*b, *a)])
        .map(|(from, to)| Support::Edge(Edge::new(from, to).unwrap()))
        .collect()
}

fn affine(s: &C, t: &C) -> Point {
    [s.clone(), t.clone(), int(1)]
}

fn value_of_constant(p: &Polynomial) -> C {
    p.evaluate(&std::array::from_fn(|_| int(0)))
}

fn constant_vector(p: &Point) -> [Polynomial; 3] {
    p.clone().map(Polynomial::constant)
}

// ---- Gaps against a third, independent evaluation --------------------

#[test]
fn gap_polynomial_and_value_at_match_the_rodrigues_definition() {
    let mut random = Random(0x9E3779B97F4A7C15);
    let (view, rotation) = configuration::coordinates();
    for _ in 0..40 {
        let gap = Gap::new(random.support(), random.below(60) as usize).unwrap();
        let polynomial = gap.polynomial(&view, &rotation).unwrap();
        for _ in 0..4 {
            let (s, t) = (random.number(1), random.number(1));
            let r: Point = std::array::from_fn(|_| random.number(1));
            let expected = my_gap(&gap, &affine(&s, &t), &r);
            let point = [s.clone(), t.clone(), r[0].clone(), r[1].clone(), r[2].clone()];
            assert_eq!(polynomial.evaluate(&point), expected, "{gap:?}");
            assert_eq!(value_at(&gap, &[s, t], &r), expected, "{gap:?}");
        }
    }
}

#[test]
fn affine_and_homogeneous_views_with_unit_third_coordinate_give_identical_polynomials() {
    let x = Polynomial::variables();
    let one = Polynomial::constant(int(1));
    let affine_view = View::Affine {
        s: x[0].clone(),
        t: x[1].clone(),
    };
    let homogeneous = View::Homogeneous([x[0].clone(), x[1].clone(), one]);
    let rotation = [x[2].clone(), x[3].clone(), x[4].clone()];
    let mut random = Random(77);
    for _ in 0..20 {
        let gap = Gap::new(random.support(), random.below(60) as usize).unwrap();
        assert_eq!(
            gap.polynomial(&affine_view, &rotation),
            gap.polynomial(&homogeneous, &rotation)
        );
    }
    for index in 0..CONSTRAINT_COUNT {
        assert_eq!(
            domain::polynomial(index, &affine_view, &rotation),
            domain::polynomial(index, &homogeneous, &rotation)
        );
    }
}

#[test]
fn a_negative_homogeneous_scale_flips_the_gap_sign() {
    // Documents why a homogeneous zoom must prove u₃ > 0: nothing in the gap
    // polynomial itself guards against a negative multiple of the view.
    let mut random = Random(5);
    for _ in 0..10 {
        let gap = Gap::new(random.support(), random.below(60) as usize).unwrap();
        let (s, t) = (random.number(1), random.number(1));
        let r: Point = std::array::from_fn(|_| random.number(1));
        let negated = affine(&s, &t).map(|c| -&c);
        let view = View::Homogeneous(constant_vector(&negated));
        let value = value_of_constant(&gap.polynomial(&view, &constant_vector(&r)).unwrap());
        assert_eq!(value, -&value_at(&gap, &[s, t], &r));
    }
}

// ---- Tight necessity at configurations with identical shadows -----------

/// A rational unit vector with positive third coordinate (inverse
/// stereographic projection), returned as the affine view `(s, t, 1)` and
/// `|(s, t, 1)|`, which is rational.
fn rational_norm_view(random: &mut Random) -> (Point, Q) {
    loop {
        let a = frac(random.below(17) as i64 - 8, 16);
        let b = frac(random.below(17) as i64 - 8, 16);
        let d = q(1) + &a * &a + &b * &b;
        let z = (q(1) - &a * &a - &b * &b) / &d;
        if z <= q(0) {
            continue;
        }
        let x = q(2) * &a / &d;
        let y = q(2) * &b / &d;
        let norm = q(1) / &z;
        return (affine(&rat(x / &z), &rat(y / &z)), norm);
    }
}

/// Cayley vectors of rotations `R` with `R K = K` or `R K = J_u K`, whose
/// shadows along `u` equal the hole's: `g` and `j g` for group elements `g`.
fn identical_shadow_rotations(u: &Point, norm: &Q) -> Vec<Point> {
    let j: Quaternion =
        [int(0), u[0].clone(), u[1].clone(), u[2].clone()].map(|c| c.scale(&(q(1) / norm)));
    let mut out = vec![[int(0), int(0), int(0)]];
    for g in geometry::symmetries() {
        for candidate in [g.clone(), geometry::quaternion_product(&j, g)] {
            if candidate[0].is_zero() {
                continue;
            }
            let w = inverse(&candidate[0]);
            out.push([&candidate[1] * &w, &candidate[2] * &w, &candidate[3] * &w]);
        }
    }
    out
}

/// Stronger than the module test of necessity: the minimum over the plug
/// vertices is exactly 0 (not only nonnegative), for edges and directions,
/// also for the half-turn class `J_u g`, whose Cayley vectors are irrational.
#[test]
fn valid_gaps_vanish_exactly_at_their_minimum_over_plug_vertices_for_identical_shadows() {
    let mut random = Random(31337);
    let (u, norm) = rational_norm_view(&mut random);
    let (s, t) = (u[0].clone(), u[1].clone());
    let rotations = identical_shadow_rotations(&u, &norm);
    let mut supports: Vec<Support> = oriented_edges()
        .into_iter()
        .filter(|support| support.check_valid(std::slice::from_ref(&u)).is_ok())
        .collect();
    for normal in [[int(1), int(0)], [int(-1), C::new(q(0), frac(1, 3))]] {
        for contact in 0..VERTEX_COUNT {
            let support = Support::Direction(Direction::new(normal.clone(), contact).unwrap());
            if support.check_valid(std::slice::from_ref(&u)).is_ok() {
                supports.push(support);
                break;
            }
        }
    }
    assert!(supports.len() > 10);
    let mut checked = 0;
    for r in rotations.iter().step_by(8) {
        for support in &supports {
            let values: Vec<C> = (0..VERTEX_COUNT)
                .map(|p| {
                    let gap = Gap::new(support.clone(), p).unwrap();
                    let value = value_at(&gap, &[s.clone(), t.clone()], r);
                    assert_eq!(value, my_gap(&gap, &u, r));
                    value
                })
                .collect();
            assert_eq!(values.iter().min().unwrap(), &int(0), "{support:?} at {r:?}");
            checked += 1;
        }
    }
    assert!(checked > 100, "{checked}");
}

// ---- Validity: soundness on view rectangles ---------------------------

#[test]
fn accepted_supports_are_valid_at_irrational_interior_views_of_path_boxes() {
    let mut random = Random(4242);
    let (mut accepted, mut refused) = (0, 0);
    for round in 0..12 {
        let b = ConfigurationBox::from_path(&random.path(8 + 3 * round)).unwrap();
        let corners = b.view_corners();
        let [s, t] = [&b.axes()[0], &b.axes()[1]];
        let mut supports = oriented_edges();
        for _ in 0..40 {
            if let Support::Direction(d) = random.support() {
                supports.push(Support::Direction(d));
            }
        }
        for support in supports {
            match support.check_valid(&corners) {
                Ok(()) => {
                    accepted += 1;
                    for _ in 0..3 {
                        let (a, c) = (random.unit_fraction(), random.unit_fraction());
                        let view = affine(
                            &(&rat(s.lo().clone()) + &(&a * &rat(s.width()))),
                            &(&rat(t.lo().clone()) + &(&c * &rat(t.width()))),
                        );
                        assert!(my_valid(&support, &view), "{support:?} at {view:?}");
                    }
                }
                Err(WitnessError::Invalid { view, vertex }) => {
                    refused += 1;
                    let n = my_normal(&support, &corners[view]);
                    let v = geometry::vertices();
                    assert!(my_dot(&n, &v[support.contact()]) < my_dot(&n, &v[vertex]));
                    for earlier in &corners[..view] {
                        assert!(my_valid(&support, earlier));
                    }
                }
                Err(error) => panic!("unexpected {error:?}"),
            }
        }
    }
    assert!(accepted > 0 && refused > 0);
}

#[test]
fn validity_agrees_with_brute_force_at_views_along_vertex_differences() {
    // Views parallel to a difference of two vertices make projected vertices
    // coincide and outline edges tie: the degenerate cases of validity.
    let v = geometry::vertices();
    let supports = oriented_edges();
    let mut views = 0;
    let mut zero_normals = 0;
    for a in (0..VERTEX_COUNT).step_by(13) {
        for b in 0..VERTEX_COUNT {
            let d = my_sub(&v[b], &v[a]);
            if d[2].sign() != Ordering::Greater {
                continue;
            }
            let u = my_scale(&d, &inverse(&d[2]));
            views += 1;
            for support in &supports {
                match support.check_valid(std::slice::from_ref(&u)) {
                    Ok(()) => assert!(my_valid(support, &u)),
                    Err(WitnessError::ZeroNormal { view: 0 }) => {
                        zero_normals += 1;
                        assert!(my_normal(support, &u).iter().all(C::is_zero));
                    }
                    Err(WitnessError::Invalid { view: 0, .. }) => {
                        assert!(!my_valid(support, &u));
                    }
                    Err(error) => panic!("unexpected {error:?}"),
                }
            }
        }
    }
    assert!(views > 50 && zero_normals > 0, "{views} {zero_normals}");
}

#[test]
fn validity_agrees_with_brute_force_at_face_plane_views() {
    // Views in the plane of a face (the face seen edge-on): ties among outline
    // vertices, where both neighbouring edges can be valid.
    let v = geometry::vertices();
    let supports = oriented_edges();
    let mut both = 0;
    for [a, b] in geometry::edges().iter().step_by(6) {
        for c in 0..VERTEX_COUNT {
            if !geometry::is_edge(*b, c) || c == *a {
                continue;
            }
            let normal = my_cross(&my_sub(&v[*b], &v[*a]), &my_sub(&v[c], &v[*b]));
            // Two views in that plane: along the edge direction rotated a bit.
            for helper in [[int(1), int(0), int(0)], [int(0), int(1), int(3)]] {
                let u = my_cross(&normal, &helper);
                if u[2].sign() != Ordering::Greater {
                    continue;
                }
                let u = my_scale(&u, &inverse(&u[2]));
                let valid: Vec<bool> = supports
                    .iter()
                    .map(|support| {
                        let accepted = support.check_valid(std::slice::from_ref(&u)).is_ok();
                        assert_eq!(accepted, my_valid(support, &u), "{support:?} at {u:?}");
                        accepted
                    })
                    .collect();
                if valid.iter().filter(|x| **x).count() > 0 {
                    both += 1;
                }
            }
        }
    }
    assert!(both > 0);
}

// ---- Robustness ----------------------------------------------------------

#[test]
fn gap_and_domain_polynomials_refuse_exponent_overflow_with_typed_errors() {
    let x = Polynomial::variables();
    let high = x[0].pow(200).unwrap();
    let view = View::Homogeneous([high.clone(), x[1].clone(), Polynomial::constant(int(1))]);
    let rotation = [x[0].pow(100).unwrap(), x[3].clone(), x[4].clone()];
    let mut random = Random(3);
    for _ in 0..6 {
        let gap = Gap::new(random.support(), 0).unwrap();
        assert_eq!(
            gap.polynomial(&view, &rotation),
            Err(WitnessError::Polynomial(PolynomialError::ExponentOverflow { variable: 0 }))
        );
    }
    let fold = Witness::Domain(DomainInequality::new(FOLD_OFFSET).unwrap());
    assert!(matches!(
        fold.polynomial(&view, &rotation),
        Err(WitnessError::Polynomial(PolynomialError::ExponentOverflow { .. }))
    ));
}
