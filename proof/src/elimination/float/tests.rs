//! Tests of the float layer against the exact layers it rounds: evaluation,
//! Bernstein coefficients (against the polynomial module's exact ones), the
//! rotation (against the exact Cayley matrix) and the outline (against a
//! brute-force extremality test and the RID's exactly valid edges).
use super::*;
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5};
use crate::arithmetic::polynomial::{Exponents, Polynomial};
use crate::elimination::zoom::tests::Random;
use crate::elimination::witness::{Edge, Support};
use num_bigint::BigInt;

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-9 * scale.max(1.0)
}

fn random_polynomial(random: &mut Random, degrees: [u64; VARIABLES], terms: usize) -> Polynomial {
    let mut p = Polynomial::zero();
    for _ in 0..terms {
        let e: Exponents = std::array::from_fn(|j| random.below(degrees[j] + 1) as u8);
        p = &p + &Polynomial::monomial(e, random.number(9));
    }
    p
}

fn random_box(random: &mut Random) -> Five<Interval> {
    std::array::from_fn(|_| {
        let lo = random.small(4);
        let width = frac(random.below(64) as i64 + 1, 1 << random.below(12));
        Interval::new(lo.clone(), lo + width).unwrap()
    })
}

fn exact_point(y: &Five<Q>) -> Five<QSqrt5> {
    std::array::from_fn(|j| QSqrt5::from_rational(y[j].clone()))
}

#[test]
fn conversions_round_huge_tiny_and_irrational_numbers() {
    let big = BigInt::from(10).pow(400);
    assert_eq!(rational(&Q::new(big.clone() * 3, big.clone())), 3.0);
    assert_eq!(rational(&Q::new(BigInt::from(1), big.clone())), 0.0);
    assert_eq!(rational(&Q::new(-big.clone() - 1, big * 2)), -0.5);
    assert_eq!(rational(&frac(-7, 8)), -0.875);
    assert!(close(number(&QSqrt5::new(frac(1, 2), frac(1, 2))), (1.0 + 5f64.sqrt()) / 2.0, 1.0));
    assert!(close(number(&QSqrt5::new(q(-2), q(1))), 5f64.sqrt() - 2.0, 1.0));
    let c = cell(&std::array::from_fn(|j| Interval::new(frac(j as i64, 7), q(1)).unwrap()));
    assert_eq!(c[0], [0.0, 1.0]);
    assert!(close(c[1][0], 1.0 / 7.0, 1.0));
}

#[test]
fn dense_evaluation_matches_exact_evaluation() {
    let mut random = Random::new(1);
    for _ in 0..200 {
        let terms = 1 + random.below(20) as usize;
        let p = random_polynomial(&mut random, [3, 2, 4, 1, 2], terms);
        let d = Dense::new(&p);
        assert_eq!(d.shape(), p.degrees().map(|k| usize::from(k) + 1));
        let y: Five<Q> = std::array::from_fn(|_| random.small(5));
        let exact = number(&p.evaluate(&exact_point(&y)));
        let scale: f64 = d.values().iter().map(|c| c.abs()).sum::<f64>() * 100.0;
        assert!(close(d.evaluate(&y.clone().map(|x| rational(&x))), exact, scale), "{p:?}");
    }
    assert_eq!(Dense::new(&Polynomial::zero()).values(), &[0.0]);
}

#[test]
fn float_bernstein_coefficients_match_the_exact_ones() {
    let mut random = Random::new(2);
    for round in 0..150 {
        let terms = 1 + random.below(25) as usize;
        let p = random_polynomial(&mut random, [2, 3, 1, 4, 2], terms);
        let b = random_box(&mut random);
        let exact = p.bernstein_coefficients(&b).unwrap();
        let float = Dense::new(&p).bernstein(&cell(&b));
        assert_eq!(float.values().len(), exact.len(), "round {round}");
        let scale = exact.iter().map(|(_, c)| number(c).abs()).fold(1.0, f64::max);
        let shape = float.shape();
        for (index, c) in exact.iter() {
            let at = (0..VARIABLES).fold(0, |a, j| a * shape[j] + usize::from(index[j]));
            assert!((float.values()[at] - number(c)).abs() <= 1e-9 * scale, "round {round} index {index:?}");
        }
    }
}

#[test]
fn bernstein_coefficients_bound_values_and_equal_corner_values() {
    let mut random = Random::new(3);
    for _ in 0..100 {
        let p = random_polynomial(&mut random, [3, 1, 2, 2, 1], 12);
        let b = random_box(&mut random);
        let d = Dense::new(&p);
        let coefficients = d.bernstein(&cell(&b));
        let (lo, hi) = (
            coefficients.values().iter().copied().fold(f64::INFINITY, f64::min),
            coefficients.max(),
        );
        let slack = 1e-9 * coefficients.magnitude().max(1.0);
        for _ in 0..10 {
            let y = random.point(&b, 8).map(|x| rational(&x));
            let v = d.evaluate(&y);
            assert!(lo - slack <= v && v <= hi + slack);
        }
        // The coefficient at an extreme index is the value at that corner.
        let shape = coefficients.shape();
        for mask in 0..32usize {
            let y: [f64; VARIABLES] = std::array::from_fn(|j| cell(&b)[j][mask >> j & 1]);
            let at = (0..VARIABLES).fold(0, |a, j| a * shape[j] + (mask >> j & 1) * (shape[j] - 1));
            assert!(close(coefficients.values()[at], d.evaluate(&y), coefficients.magnitude() * 1e3));
        }
    }
}

#[test]
fn variation_is_the_largest_neighbouring_difference() {
    let mut random = Random::new(4);
    for _ in 0..50 {
        let d = Dense::new(&random_polynomial(&mut random, [2, 3, 0, 1, 2], 10));
        let shape = d.shape();
        for axis in 0..VARIABLES {
            let mut expected: f64 = 0.0;
            let total = d.values().len();
            for at in 0..total {
                let index: Vec<usize> = (0..VARIABLES)
                    .map(|j| at / shape[j + 1..].iter().product::<usize>() % shape[j])
                    .collect();
                if index[axis] + 1 < shape[axis] {
                    let step: usize = shape[axis + 1..].iter().product();
                    expected = expected.max((d.values()[at + step] - d.values()[at]).abs());
                }
            }
            assert_eq!(d.variation(axis), expected);
        }
        assert_eq!(d.variation(2), 0.0);
    }
}

#[test]
fn float_rotation_is_the_exact_cayley_rotation() {
    let mut random = Random::new(5);
    for _ in 0..100 {
        let r: [Q; 3] = std::array::from_fn(|_| random.small(6));
        let rational_r = r.clone().map(QSqrt5::from_rational);
        let constant = rational_r.clone().map(Polynomial::constant);
        let exact = geometry::cayley_matrix(&constant).unwrap();
        let scale = 1.0 + r.iter().map(|x| rational(x).powi(2)).sum::<f64>();
        let m = rotation(&r.clone().map(|x| rational(&x)));
        for i in 0..3 {
            for j in 0..3 {
                let e = number(&exact[i][j].evaluate(&std::array::from_fn(|_| QSqrt5::zero())));
                assert!(close(m[i][j], e / scale, 1.0));
            }
            for j in 0..3 {
                let product = dot(&m[i], &m[j]);
                assert!(close(product, f64::from(u8::from(i == j)), 1.0));
            }
        }
        let v = vertices()[random.below(60) as usize];
        let turned = apply(&m, &v);
        assert!(close(dot(&turned, &turned), dot(&v, &v), 10.0));
        assert_eq!(cross(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]), [0.0, 0.0, 1.0]);
    }
}

/// Whether the integer point `p` is a strict vertex of the hull of
/// `points`, exactly: it is not in the convex hull of the other (distinct)
/// points, that is (Carathéodory) in no closed triangle or segment of them.
fn extreme(points: &[[f64; 2]], p: [f64; 2]) -> bool {
    let others: Vec<[f64; 2]> = points.iter().copied().filter(|q| *q != p).collect();
    let turn = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let on_segment = |a: [f64; 2], b: [f64; 2]| {
        turn(a, b, p) == 0.0 && (p[0] - a[0]) * (p[0] - b[0]) <= 0.0 && (p[1] - a[1]) * (p[1] - b[1]) <= 0.0
    };
    let n = others.len();
    for i in 0..n {
        for j in i + 1..n {
            if on_segment(others[i], others[j]) {
                return false;
            }
            for k in j + 1..n {
                let (a, b, c) = (others[i], others[j], others[k]);
                let signs = [turn(a, b, p), turn(b, c, p), turn(c, a, p)];
                if turn(a, b, c) != 0.0 && (signs.iter().all(|s| *s >= 0.0) || signs.iter().all(|s| *s <= 0.0)) {
                    return false;
                }
            }
        }
    }
    true
}

#[test]
fn hull_vertices_are_exactly_the_extreme_points_counter_clockwise() {
    let mut random = Random::new(6);
    for round in 0..200 {
        let count = 3 + random.below(20) as usize;
        let mut points: Vec<[f64; 2]> =
            (0..count).map(|_| [random.below(41) as f64 - 20.0, random.below(41) as f64 - 20.0]).collect();
        if round % 3 == 0 {
            // A duplicate and a midpoint of two points.
            points.push(points[0]);
            points.push([(points[1][0] + points[2][0]) / 2.0, (points[1][1] + points[2][1]) / 2.0]);
        }
        let h = hull(&points);
        let mut expected: Vec<[f64; 2]> = points.iter().copied().filter(|&p| extreme(&points, p)).collect();
        expected.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        expected.dedup();
        let mut found: Vec<[f64; 2]> = h.iter().map(|&i| points[i]).collect();
        found.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        if expected.len() >= 3 {
            assert_eq!(found, expected, "round {round}: {points:?}");
            for i in 0..h.len() {
                let (a, b, c) = (points[h[i]], points[h[(i + 1) % h.len()]], points[h[(i + 2) % h.len()]]);
                assert!((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0.0);
            }
        }
    }
    assert_eq!(hull(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.0, 1.0], [1.0, 0.5]]), vec![0, 2, 3]);
}

#[test]
fn hull_of_the_shadow_is_the_outline_of_the_exactly_valid_edges() {
    let mut random = Random::new(7);
    for _ in 0..40 {
        let u: geometry::Point = [
            QSqrt5::from_rational(random.between(&frac(1, 97), &frac(2, 3), 10)),
            QSqrt5::from_rational(random.between(&frac(1, 89), &frac(2, 5), 10)),
            QSqrt5::one(),
        ];
        let shadow: Vec<[f64; 2]> = vertices().iter().map(|v| screen(&point(&u), v)).collect();
        let h = hull(&shadow);
        // Every counter-clockwise outline edge a → b is an edge of the RID whose
        // support (outward normal) is exactly valid at u.
        for i in 0..h.len() {
            let (a, b) = (h[i], h[(i + 1) % h.len()]);
            let edge = Edge::new(a, b).expect("outline edges are RID edges at a generic view");
            let support = Support::Edge(edge);
            let (sa, sb) = (shadow[a], shadow[b]);
            // n(u)·X = d × T_u(X) up to sign; choose the orientation pointing outward.
            let outward = [sb[1] - sa[1], sa[0] - sb[0]];
            let centre = [0.0, 0.0];
            let toward = outward[0] * (centre[0] - sa[0]) + outward[1] * (centre[1] - sa[1]);
            assert!(toward < 0.0);
            let valid = support.check_valid(std::slice::from_ref(&u)).is_ok();
            let reversed = Support::Edge(Edge::new(b, a).unwrap()).check_valid(std::slice::from_ref(&u)).is_ok();
            assert!(valid != reversed, "exactly one orientation of an outline edge is valid");
        }
    }
}
