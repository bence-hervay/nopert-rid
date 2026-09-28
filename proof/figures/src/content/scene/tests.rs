use super::*;
use std::cmp::Ordering;
use crate::random::Random;
use num_traits::ToPrimitive;
use rid::arithmetic::exact::frac;
use rid::problem::geometry::{rotate, Point as Space, VERTEX_COUNT};

fn rational(n: i64, d: i64) -> QSqrt5 {
    QSqrt5::from_rational(frac(n, d))
}

fn float(x: &Q) -> f64 {
    // Exact enough for the comparisons below: both parts are far below 2^1000.
    let scale = 1u64 << 53;
    let n = x.numer() * BigInt::from(scale) / x.denom();
    n.to_f64().expect("finite") / scale as f64
}

fn float_field(x: &QSqrt5) -> f64 {
    float(x.rational_part()) + float(x.sqrt5_part()) * 5f64.sqrt()
}

/// `c·d·√m ≤ p`, decided exactly for rational `c` and `d, m > 0` in Q(√5).
fn scaled_root_at_most(c: &Q, d: &QSqrt5, m: &QSqrt5, p: &QSqrt5) -> bool {
    let g = d.scale(c);
    let zero = QSqrt5::zero();
    match (g >= zero, *p >= zero) {
        (false, true) => true,
        (true, false) => false,
        (true, true) => &(&g * &g) * m <= p * p,
        (false, false) => &(&g * &g) * m >= p * p,
    }
}

/// Whether the rational interval contains `p / (d·√m)` (with `d, m > 0`).
fn contains_exact(interval: &Interval, p: &QSqrt5, d: &QSqrt5, m: &QSqrt5) -> bool {
    let negative = |x: &QSqrt5| -x;
    scaled_root_at_most(interval.lo(), d, m, p)
        && scaled_root_at_most(&-interval.hi(), d, m, &negative(p))
}

/// The exact screen numerators and scales of a configuration, computed with
/// the geometry module's exact rotation instead of this module's polynomials.
struct Independent {
    /// `[body][vertex] = [P₁, P₂]`.
    numerators: [Vec<[QSqrt5; 2]>; 2],
    /// `[body] = D`.
    denominators: [QSqrt5; 2],
    /// `[L, LN]`.
    normalisation: [QSqrt5; 2],
}

fn independent(c: &[QSqrt5; 5]) -> Independent {
    let [s, t, r1, r2, r3] = c;
    let one = QSqrt5::one();
    let l = &one + &(s * s);
    let n = &l + &(t * t);
    // rotate((1, r), p) = (1 + |r|²) R(r) p, the Cayley numerator.
    let rotation = [one.clone(), r1.clone(), r2.clone(), r3.clone()];
    let d = &one + &(&(&(r1 * r1) + &(r2 * r2)) + &(r3 * r3));
    let numerators = |w: &Space| {
        [
            &w[0] - &(s * &w[2]),
            &(&(&l * &w[1]) - &(&(s * t) * &w[0])) - &(t * &w[2]),
        ]
    };
    Independent {
        numerators: [
            vertices().iter().map(numerators).collect(),
            vertices().iter().map(|p| numerators(&rotate(&rotation, p))).collect(),
        ],
        denominators: [one, d],
        normalisation: [l.clone(), &l * &n],
    }
}

fn random_configuration(random: &mut Random, algebraic: bool) -> [QSqrt5; 5] {
    std::array::from_fn(|_| {
        let a = random.rational(&frac(-1, 1), &frac(1, 1), 20);
        let b = if algebraic {
            random.rational(&frac(-1, 3), &frac(1, 3), 20)
        } else {
            Q::zero()
        };
        QSqrt5::new(a, b)
    })
}

#[test]
fn enclosures_contain_the_exact_positions_and_are_narrow() {
    let mut random = Random::new(1);
    let limit = Q::new(BigInt::from(1), BigInt::from(1) << POINT_BITS);
    for case in 0..24 {
        let c = random_configuration(&mut random, case % 2 == 1);
        let scene = ConfigurationScene::new(&c).unwrap();
        let exact = independent(&c);
        for (b, body) in [Body::Hole, Body::Plug].into_iter().enumerate() {
            for v in 0..VERTEX_COUNT {
                let r = scene.enclosure(body, v);
                for j in 0..2 {
                    assert!(r[j].width() <= limit, "case {case} {body:?} {v}");
                    let (p, d, m) = (&exact.numerators[b][v][j], &exact.denominators[b], &exact.normalisation[j]);
                    assert!(contains_exact(&r[j], p, d, m), "case {case} {body:?} {v} {j}");
                }
                // The drawn position is the enclosure's midpoint.
                assert_eq!(scene.position(body, v), &[r[0].midpoint(), r[1].midpoint()]);
            }
        }
    }
}

#[test]
fn positions_agree_with_an_independent_floating_point_projection() {
    let mut random = Random::new(2);
    for case in 0..24 {
        let c = random_configuration(&mut random, case % 2 == 0);
        let scene = ConfigurationScene::new(&c).unwrap();
        let [s, t] = [float_field(&c[0]), float_field(&c[1])];
        let norm = (1.0 + s * s + t * t).sqrt();
        let e1 = [1.0, 0.0, -s].map(|x| x / (1.0 + s * s).sqrt());
        let e2 = [-s * t, 1.0 + s * s, -t].map(|x| x / ((1.0 + s * s).sqrt() * norm));
        // The basis is orthonormal and perpendicular to the view.
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        assert!((dot(e1, e2)).abs() < 1e-12 && (dot(e1, e1) - 1.0).abs() < 1e-12);
        assert!(dot(e1, [s, t, 1.0]).abs() < 1e-12 && dot(e2, [s, t, 1.0]).abs() < 1e-12);
        // The plug by the quaternion (1, r), divided by its squared norm.
        let quaternion = [QSqrt5::one(), c[2].clone(), c[3].clone(), c[4].clone()];
        let norm2: f64 = quaternion.iter().map(|x| float_field(x).powi(2)).sum();
        for (body, turned) in [(Body::Hole, false), (Body::Plug, true)] {
            for (v, p) in vertices().iter().enumerate() {
                let w = if turned { rotate(&quaternion, p) } else { p.clone() };
                let w = w.map(|x| float_field(&x) / if turned { norm2 } else { 1.0 });
                let expected = [dot(e1, w), dot(e2, w)];
                let drawn = scene.position(body, v);
                for j in 0..2 {
                    assert!((float(&drawn[j]) - expected[j]).abs() < 1e-9, "case {case} {body:?} {v}");
                }
            }
        }
    }
}

#[test]
fn aligned_and_symmetric_rotations_give_the_hole_shadow() {
    // r = 0: the plug is the hole.
    let zero = QSqrt5::zero();
    let aligned = [rational(1, 3), rational(1, 7), zero.clone(), zero.clone(), zero];
    let scene = ConfigurationScene::new(&aligned).unwrap();
    for v in 0..VERTEX_COUNT {
        assert_eq!(scene.position(Body::Hole, v), scene.position(Body::Plug, v));
    }
    assert_eq!(scene.silhouette(Body::Hole), scene.silhouette(Body::Plug));
    // r = (1, 1, 1) is the rotation (1, 1, 1, 1)/2 of the solid's group:
    // the plug's vertices are the hole's, permuted, with the same shadow.
    let one = QSqrt5::one();
    let c = [rational(2, 5), rational(-1, 3), one.clone(), one.clone(), one];
    let exact = independent(&c);
    let four = QSqrt5::integer(4);
    let hole: Vec<[QSqrt5; 2]> = exact.numerators[0].clone();
    let mut matched = vec![false; VERTEX_COUNT];
    for plug in &exact.numerators[1] {
        let k = hole
            .iter()
            .position(|h| &(&h[0] * &four) == &plug[0] && &(&h[1] * &four) == &plug[1])
            .expect("every plug vertex is a hole vertex");
        assert!(!matched[k]);
        matched[k] = true;
    }
    let scene = ConfigurationScene::new(&c).unwrap();
    let outline = |body| {
        let mut points: Vec<Point> = scene.silhouette(body).iter().map(|&j| scene.position(body, j).clone()).collect();
        points.sort();
        points
    };
    let (hole, plug) = (outline(Body::Hole), outline(Body::Plug));
    assert_eq!(hole.len(), plug.len());
    let near = Q::new(BigInt::from(1), BigInt::from(1) << 90);
    for (a, b) in hole.iter().zip(&plug) {
        assert!((&a[0] - &b[0]) <= near && (&b[0] - &a[0]) <= near);
        assert!((&a[1] - &b[1]) <= near && (&b[1] - &a[1]) <= near);
    }
}

/// Checks that `hull` of `points` is a strictly convex counterclockwise
/// polygon of input points containing every point, using the lowest index
/// of equal points: the only such polygon is the convex hull.
fn check_hull<T>(points: &[[T; 2]], hull: &[usize])
where
    T: Ord + Clone + std::fmt::Debug,
    for<'a> &'a T: Sub<Output = T> + Mul<Output = T>,
{
    let cross = |a: &[T; 2], b: &[T; 2], c: &[T; 2]| {
        (&(&b[0] - &a[0]) * &(&c[1] - &a[1])).cmp(&(&(&b[1] - &a[1]) * &(&c[0] - &a[0])))
    };
    for &i in hull {
        assert!(points[..i].iter().all(|p| p != &points[i]), "not the lowest index of {i}");
    }
    let mut distinct: Vec<&[T; 2]> = points.iter().collect();
    distinct.sort();
    distinct.dedup();
    match distinct.len() {
        0 => assert!(hull.is_empty()),
        1 => assert_eq!(hull.len(), 1),
        _ => {
            assert!(hull.len() >= 2);
            let n = hull.len();
            if n == 2 {
                // All points on the segment between the two ends.
                let (a, b) = (&points[hull[0]], &points[hull[1]]);
                for p in points {
                    assert_eq!(cross(a, b, p), Ordering::Equal);
                    assert!(p >= a.min(b) && p <= a.max(b));
                }
                return;
            }
            for k in 0..n {
                let (a, b, c) = (&points[hull[k]], &points[hull[(k + 1) % n]], &points[hull[(k + 2) % n]]);
                assert_eq!(cross(a, b, c), Ordering::Greater, "not strictly convex");
                for p in points {
                    assert_ne!(cross(a, b, p), Ordering::Less, "a point outside the hull");
                }
            }
        }
    }
}

#[test]
fn hull_is_the_convex_hull_of_random_small_sets() {
    let mut random = Random::new(3);
    for case in 0..400 {
        let count = random.integer(0, 14) as usize;
        let spread = if case % 4 == 0 { 2 } else { 6 };
        let points: Vec<[Q; 2]> = (0..count)
            .map(|_| {
                if case % 7 == 0 {
                    // Collinear sets.
                    let x = random.integer(-spread, spread);
                    [frac(x, 1), frac(2 * x + 1, 1)]
                } else {
                    [frac(random.integer(-spread, spread), 1), frac(random.integer(-spread, spread), 1)]
                }
            })
            .collect();
        check_hull(&points, &hull(&points));
    }
}

#[test]
fn hull_decides_near_collinear_field_points_exactly() {
    // (0, 0), (1, φ), (2, 2φ + ε) with ε = (φ - 1)^60 ≈ 3·10⁻¹³ > 0: a strict
    // left turn that floating point cannot see.
    let phi = QSqrt5::new(frac(1, 2), frac(1, 2));
    let mut epsilon = QSqrt5::one();
    for _ in 0..60 {
        epsilon = &epsilon * &(&phi - &QSqrt5::one());
    }
    let two = QSqrt5::integer(2);
    let points = vec![
        [QSqrt5::zero(), QSqrt5::zero()],
        [QSqrt5::one(), phi.clone()],
        [two.clone(), &(&two * &phi) + &epsilon],
    ];
    assert_eq!(hull(&points), vec![0, 1, 2]);
    check_hull(&points, &hull(&points));
    let flat = vec![points[0].clone(), points[1].clone(), [two.clone(), &two * &phi]];
    assert_eq!(hull(&flat), vec![0, 2]);
    let below = vec![points[0].clone(), points[1].clone(), [two.clone(), &(&two * &phi) - &epsilon]];
    assert_eq!(hull(&below), vec![0, 2, 1]);
}

#[test]
fn silhouettes_are_the_hulls_of_the_exact_positions() {
    let mut random = Random::new(4);
    for case in 0..12 {
        let c = random_configuration(&mut random, case % 2 == 0);
        let scene = ConfigurationScene::new(&c).unwrap();
        let exact = independent(&c);
        for (b, body) in [Body::Hole, Body::Plug].into_iter().enumerate() {
            check_hull(&exact.numerators[b], scene.silhouette(body));
            match scene.outline(body) {
                Shape::Polygon(points) => assert_eq!(points.len(), scene.silhouette(body).len()),
                _ => panic!("an outline is a polygon"),
            }
        }
        assert_eq!(scene.wireframe(Body::Plug).len(), edges().len());
    }
}

fn random_box(random: &mut Random, width: &Q) -> ConfigurationBox {
    ConfigurationBox::new(std::array::from_fn(|axis| {
        let lo = random.rational(&frac(-1, 2), &frac(1, 2), 16);
        let hi = if axis == 3 { lo.clone() } else { &lo + &random.rational(&Q::zero(), width, 8) };
        Interval::new(lo, hi).unwrap()
    }))
}

#[test]
fn box_envelopes_contain_every_sampled_configuration() {
    let mut random = Random::new(5);
    for case in 0..10 {
        let width = if case % 2 == 0 { frac(1, 8) } else { frac(1, 1000) };
        let b = random_box(&mut random, &width);
        let scene = BoxScene::new(&b).unwrap();
        for sample in 0..6 {
            let x: [Q; 5] = if sample < 2 {
                b.corner(if sample == 0 { 0 } else { 31 })
            } else {
                std::array::from_fn(|j| random.rational(b.axes()[j].lo(), b.axes()[j].hi(), 12))
            };
            let c = x.map(QSqrt5::from_rational);
            let exact = independent(&c);
            for (k, body) in [Body::Hole, Body::Plug].into_iter().enumerate() {
                for v in 0..VERTEX_COUNT {
                    let r = scene.rectangle(body, v);
                    for j in 0..2 {
                        let (p, d, m) = (&exact.numerators[k][v][j], &exact.denominators[k], &exact.normalisation[j]);
                        assert!(contains_exact(&r[j], p, d, m), "case {case} sample {sample} {body:?} {v}");
                    }
                }
            }
        }
    }
}

#[test]
fn edge_envelopes_are_convex_hulls_of_both_rectangles() {
    let mut random = Random::new(6);
    let b = random_box(&mut random, &frac(1, 16));
    let scene = BoxScene::new(&b).unwrap();
    for body in [Body::Hole, Body::Plug] {
        for &edge in edges() {
            let polygon = scene.edge_envelope(body, edge);
            let corners: Vec<Point> = edge.iter().flat_map(|&v| corners(scene.rectangle(body, v))).collect();
            let indices: Vec<usize> = polygon
                .iter()
                .map(|p| corners.iter().position(|c| c == p).expect("a corner"))
                .collect();
            check_hull(&corners, &indices);
        }
        assert_eq!(scene.vertex_envelopes(body).len(), VERTEX_COUNT);
        assert_eq!(scene.edge_envelopes(body).len(), edges().len());
    }
}

#[test]
fn a_point_box_has_tight_envelopes() {
    let x = [frac(1, 5), frac(1, 10), frac(1, 8), frac(1, 16), frac(-1, 32)];
    let scene = BoxScene::new(&ConfigurationBox::point(&x)).unwrap();
    let exact = independent(&x.clone().map(QSqrt5::from_rational));
    let narrow = Q::new(BigInt::from(1), BigInt::from(1) << 120);
    for (k, body) in [Body::Hole, Body::Plug].into_iter().enumerate() {
        for v in 0..VERTEX_COUNT {
            let r = scene.rectangle(body, v);
            for j in 0..2 {
                assert!(r[j].width() < narrow);
                let (p, d, m) = (&exact.numerators[k][v][j], &exact.denominators[k], &exact.normalisation[j]);
                assert!(contains_exact(&r[j], p, d, m));
            }
        }
    }
}

#[test]
fn boxes_with_negative_views_and_large_rotations_are_enclosed() {
    // s ranges over [-3, 2]: the Bernstein lower bound of 1 + s² is below 1
    // and is raised to 1; r is far outside the root box.
    let axes = [(-3, 2), (-1, 1), (1, 2), (-2, -1), (0, 0)]
        .map(|(lo, hi)| Interval::new(frac(lo, 1), frac(hi, 1)).unwrap());
    let b = ConfigurationBox::new(axes);
    let scene = BoxScene::new(&b).unwrap();
    for mask in [0u8, 5, 10, 31] {
        let c = b.corner(mask).map(QSqrt5::from_rational);
        let exact = independent(&c);
        for (k, body) in [Body::Hole, Body::Plug].into_iter().enumerate() {
            for v in 0..VERTEX_COUNT {
                for j in 0..2 {
                    let (p, d, m) = (&exact.numerators[k][v][j], &exact.denominators[k], &exact.normalisation[j]);
                    assert!(contains_exact(&scene.rectangle(body, v)[j], p, d, m));
                }
            }
        }
    }
}

#[test]
fn imprecise_positions_are_refused_rather_than_drawn() {
    // (√5 - 2)^60 ≈ 10⁻³⁸ is written a + b√5 with |a|, |b| ≈ 2¹²⁵: its
    // rational enclosure is far wider than 2^-100.
    let mut s = QSqrt5::one();
    for _ in 0..60 {
        s = &s * &QSqrt5::new(frac(-2, 1), Q::one());
    }
    let zero = QSqrt5::zero();
    let c = [s, zero.clone(), zero.clone(), zero.clone(), zero];
    assert!(matches!(ConfigurationScene::new(&c), Err(SceneError::PrecisionLost { .. })));
    // Large but plainly written coordinates keep their precision.
    let large = QSqrt5::new(Q::new(BigInt::from(1) << 90, BigInt::from(1)), frac(1, 3));
    let c = [large.clone(), large.clone(), large.clone(), large.clone(), large];
    assert!(ConfigurationScene::new(&c).is_ok());
}

#[test]
fn radical_brackets_are_correct() {
    let root = sqrt5();
    assert!(root.lo() * root.lo() < frac(5, 1) && frac(5, 1) < root.hi() * root.hi());
    assert_eq!(root.width(), Q::new(BigInt::from(1), scale()));
    let mut random = Random::new(7);
    for case in 0..300 {
        let lo = random.rational(&Q::zero(), &frac(50, 1), 30);
        let hi = &lo + &random.rational(&Q::zero(), &frac(3, 1), 20);
        let x = Interval::new(lo.clone(), hi.clone()).unwrap();
        let r = square_root(&x);
        assert!(r.lo() * r.lo() <= lo && hi <= r.hi() * r.hi(), "case {case}");
        let step = Q::new(BigInt::from(1), scale());
        let below = r.lo() + &step;
        let above = r.hi() - &step;
        assert!(&below * &below > lo, "the lower bound is the floor");
        assert!(above < Q::zero() || &above * &above < hi, "the upper bound is the ceiling");
        if lo > Q::zero() {
            let inverse = reciprocal(&x);
            assert!(inverse.contains(&(Q::one() / &lo)) && inverse.contains(&(Q::one() / &hi)));
        }
    }
    // Exact squares give exact roots.
    let square = Interval::new(frac(9, 4), frac(49, 16)).unwrap();
    assert_eq!(square_root(&square), Interval::new(frac(3, 2), frac(7, 4)).unwrap());
    // Enclosures of field elements of either sign of the √5 part.
    for (a, b) in [(1, 1), (1, -1), (-3, 2), (0, 0)] {
        let x = QSqrt5::new(frac(a, 1), frac(b, 1));
        let e = enclose(&x);
        assert!(QSqrt5::from_rational(e.lo().clone()) <= x && x <= QSqrt5::from_rational(e.hi().clone()));
    }
}

mod review;

#[test]
fn front_and_back_edges_split_every_body_at_a_generic_view() {
    // A generic configuration: no face is seen edge-on, so the silhouette's
    // edges are exactly the edges with one face towards the viewer. A
    // centrally symmetric solid then has as many back edges as front edges
    // off the silhouette.
    let c = [rational(1, 7), rational(1, 11), rational(1, 13), rational(-1, 17), rational(1, 19)];
    let scene = ConfigurationScene::new(&c).unwrap();
    for body in [Body::Hole, Body::Plug] {
        let (front, back) = (scene.front_wireframe(body).len(), scene.back_wireframe(body).len());
        assert_eq!(front + back, edges().len());
        assert_eq!(front - back, scene.silhouette(body).len(), "{body:?}");
    }
}
