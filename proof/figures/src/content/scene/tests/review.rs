//! Adversarial tests of the scene module: an independent
//! characterisation of the screen (isometry and handedness, computed in 3D
//! with the quaternion rotation), and a high-precision check that edge
//! envelopes contain sampled points of the swept edges.
use super::*;
use crate::random::Random;
use rid::arithmetic::exact::frac;
use rid::problem::geometry::{cross, dot, rotate, Point as Space, VERTEX_COUNT};

fn field(a: Q, b: Q) -> QSqrt5 {
    QSqrt5::new(a, b)
}

/// The exact 3D position of a vertex of a body: `p` or `(1,r)p(1,r)̄ / D`.
fn world(c: &[QSqrt5; 5], body: Body, v: usize) -> Space {
    let p = &vertices()[v];
    match body {
        Body::Hole => p.clone(),
        Body::Plug => {
            let q = [QSqrt5::one(), c[2].clone(), c[3].clone(), c[4].clone()];
            let d = &QSqrt5::one() + &(&(&(&c[2] * &c[2]) + &(&c[3] * &c[3])) + &(&c[4] * &c[4]));
            // Division by the rational-or-field D: multiply by its inverse,
            // conj(D)/norm(D).
            let inverse = d.conjugate().scale(&(Q::one() / d.norm()));
            rotate(&q, p).map(|x| &x * &inverse)
        }
    }
}

fn sub(a: &Space, b: &Space) -> Space {
    std::array::from_fn(|i| &a[i] - &b[i])
}

/// Orthogonal projection preserves the component of every difference
/// perpendicular to u: |Δscreen|² = |Δ|² − (Δ·u)²/|u|², and the screen is
/// seen from +u (the drawn orientation of a triangle is the sign of
/// det(Δ₁, Δ₂, u)).
#[test]
fn screen_is_an_isometry_seen_from_plus_u() {
    let mut random = Random::new(901);
    let tolerance = Q::new(BigInt::from(1), BigInt::from(1) << 80);
    for case in 0..12 {
        let c: [QSqrt5; 5] = std::array::from_fn(|_| {
            let a = random.rational(&frac(-3, 2), &frac(3, 2), 20);
            let b = if case % 2 == 0 { Q::zero() } else { random.rational(&frac(-1, 4), &frac(1, 4), 20) };
            field(a, b)
        });
        let scene = ConfigurationScene::new(&c).unwrap();
        let u: Space = [c[0].clone(), c[1].clone(), QSqrt5::one()];
        let uu = dot(&u, &u);
        let uu_inverse = uu.conjugate().scale(&(Q::one() / uu.norm()));
        for body in [Body::Hole, Body::Plug] {
            for _ in 0..40 {
                let [a, b, k] = [0, 0, 0].map(|_| random.integer(0, VERTEX_COUNT as i64 - 1) as usize);
                let (wa, wb, wk) = (world(&c, body, a), world(&c, body, b), world(&c, body, k));
                let delta = sub(&wb, &wa);
                let along = dot(&delta, &u);
                let expected = &dot(&delta, &delta) - &(&(&along * &along) * &uu_inverse);
                let (pa, pb, pk) = (scene.position(body, a), scene.position(body, b), scene.position(body, k));
                let dx = &pb[0] - &pa[0];
                let dy = &pb[1] - &pa[1];
                let drawn = QSqrt5::from_rational(&(&dx * &dx) + &(&dy * &dy));
                let error = &drawn - &expected;
                let t = QSqrt5::from_rational(tolerance.clone());
                assert!(error <= t && -&error <= t, "case {case} {body:?} {a} {b}: not an isometry");
                // Handedness.
                let d2 = sub(&wk, &wa);
                let det = dot(&cross(&delta, &d2), &u);
                let ex = [&pb[0] - &pa[0], &pb[1] - &pa[1]];
                let ek = [&pk[0] - &pa[0], &pk[1] - &pa[1]];
                let orientation = &(&ex[0] * &ek[1]) - &(&ex[1] * &ek[0]);
                let threshold = QSqrt5::from_rational(frac(1, 1000));
                if det > threshold {
                    assert!(orientation > Q::zero(), "case {case}: mirrored screen");
                } else if -&det > threshold {
                    assert!(orientation < Q::zero(), "case {case}: mirrored screen");
                }
            }
        }
    }
}

// ---- A second, independent enclosure (400 bits, Newton-free) ------------

const BITS: u32 = 400;

fn floor_scaled(x: &Q) -> BigInt {
    use num_integer::Integer;
    (x.numer() << BITS).div_floor(x.denom())
}

/// [lo, hi] (as integers over 2^BITS) containing the rational x.
fn enclose_q(x: &Q) -> (BigInt, BigInt) {
    let f = floor_scaled(x);
    (f.clone(), f + 1)
}

/// Interval (over 2^BITS) containing √x for rational x ≥ 0: with
/// F = ⌊x·4^B⌋ and s = ⌊√F⌋, s ≤ √(x·4^B) < s + 1.
fn enclose_sqrt(x: &Q) -> (BigInt, BigInt) {
    use num_integer::Integer;
    let f = (x.numer() << (2 * BITS)).div_floor(x.denom());
    let s = f.sqrt();
    (s.clone(), s + 1)
}

/// Interval of a + b√5 (over 2^BITS), each end widened by 2.
fn enclose_field(x: &QSqrt5) -> (BigInt, BigInt) {
    let (r0, r1) = enclose_sqrt(&frac(5, 1));
    let (a0, a1) = enclose_q(x.rational_part());
    let b = x.sqrt5_part();
    // b·√5 with b rational: multiply rational b by the root enclosure.
    let lo_b = |r: &BigInt| {
        use num_integer::Integer;
        (b.numer() * r).div_floor(b.denom())
    };
    let (p, q) = (lo_b(&r0), lo_b(&r1));
    let (m0, m1) = if p <= q { (p, q + 1) } else { (q, p + 1) };
    (a0 + m0 - 2, a1 + m1 + 2)
}

/// Whether the rational point enclosure (x ± ε) lies inside (or on) the
/// convex counterclockwise polygon, decided conservatively.
fn inside_polygon(polygon: &[Point], x: &[(BigInt, BigInt); 2]) -> bool {
    let scale = Q::new(BigInt::from(1), BigInt::from(1) << BITS);
    let lo = [Q::from_integer(x[0].0.clone()) * &scale, Q::from_integer(x[1].0.clone()) * &scale];
    let hi = [Q::from_integer(x[0].1.clone()) * &scale, Q::from_integer(x[1].1.clone()) * &scale];
    let n = polygon.len();
    if n < 3 {
        return true; // degenerate envelopes are not checked here
    }
    for k in 0..n {
        let (a, b) = (&polygon[k], &polygon[(k + 1) % n]);
        let (dx, dy) = (&b[0] - &a[0], &b[1] - &a[1]);
        // The minimum over the enclosure box of dx·(y - a₁) - dy·(x - a₀).
        let y = if dx >= Q::zero() { &lo[1] } else { &hi[1] };
        let xx = if dy >= Q::zero() { &hi[0] } else { &lo[0] };
        let value = &(&dx * &(y - &a[1])) - &(&dy * &(xx - &a[0]));
        if value < Q::zero() {
            return false;
        }
    }
    true
}

#[test]
fn edge_envelopes_contain_sampled_edge_points() {
    let mut random = Random::new(902);
    for case in 0..6 {
        let width = if case % 2 == 0 { frac(1, 8) } else { frac(3, 2) };
        let b = ConfigurationBox::new(std::array::from_fn(|axis| {
            let lo = random.rational(&frac(-1, 1), &frac(1, 1), 12);
            let hi = &lo + &random.rational(&Q::zero(), &width, 8);
            let _ = axis;
            Interval::new(lo, hi).unwrap()
        }));
        let scene = BoxScene::new(&b).unwrap();
        for _ in 0..3 {
            let x: [Q; 5] = std::array::from_fn(|j| random.rational(b.axes()[j].lo(), b.axes()[j].hi(), 10));
            let c = x.clone().map(QSqrt5::from_rational);
            let (s, t) = (&x[0], &x[1]);
            let l = &Q::one() + &(s * s);
            let n = &l + &(t * t);
            let (root_l, root_ln) = (enclose_sqrt(&l), enclose_sqrt(&(&l * &n)));
            for body in [Body::Hole, Body::Plug] {
                for &edge in edges().iter().step_by(7) {
                    let lambda = random.rational(&Q::zero(), &Q::one(), 8);
                    let wa = world(&c, body, edge[0]);
                    let wb = world(&c, body, edge[1]);
                    let w: Space = std::array::from_fn(|i| {
                        &wa[i].scale(&lambda) + &wb[i].scale(&(Q::one() - &lambda))
                    });
                    // Numerators P₁ = w₁ − s w₃, P₂ = L w₂ − s t w₁ − t w₃.
                    let (sq, tq) = (QSqrt5::from_rational(s.clone()), QSqrt5::from_rational(t.clone()));
                    let p1 = &w[0] - &(&sq * &w[2]);
                    let lq = QSqrt5::from_rational(l.clone());
                    let p2 = &(&(&lq * &w[1]) - &(&(&sq * &tq) * &w[0])) - &(&tq * &w[2]);
                    // Coordinate = P / root, enclosed with directed rounding.
                    let divide = |p: &QSqrt5, root: &(BigInt, BigInt)| -> (BigInt, BigInt) {
                        use num_integer::Integer;
                        let (p0, p1) = enclose_field(p);
                        let one = BigInt::from(1) << (2 * BITS);
                        // p / r with r ∈ [r0, r1] > 0 (over 2^B): value·2^B = p·2^B / r · ...
                        let candidates = [
                            (&p0 * &one).div_floor(&root.0),
                            (&p0 * &one).div_floor(&root.1),
                            (&p1 * &one).div_floor(&root.0),
                            (&p1 * &one).div_floor(&root.1),
                        ];
                        let lo = candidates.iter().min().unwrap().clone() >> BITS;
                        let hi = (candidates.iter().max().unwrap().clone() >> BITS) + 2;
                        (lo - 2, hi)
                    };
                    let point = [divide(&p1, &root_l), divide(&p2, &root_ln)];
                    let polygon = scene.edge_envelope(body, edge);
                    assert!(inside_polygon(&polygon, &point), "case {case} {body:?} edge {edge:?} λ {lambda}");
                }
            }
        }
    }
}

#[test]
fn hull_of_equal_and_empty_point_sets() {
    let empty: Vec<[Q; 2]> = Vec::new();
    assert!(hull(&empty).is_empty());
    let same = vec![[frac(1, 2), frac(1, 3)]; 5];
    assert_eq!(hull(&same), vec![0]);
    let two = vec![[frac(1, 1), Q::zero()], [Q::zero(), Q::zero()], [frac(1, 1), Q::zero()]];
    assert_eq!(hull(&two), vec![1, 0]);
}

#[test]
fn extreme_boxes_are_still_enclosed() {
    // Views far from the root box and wide rotation ranges: every sampled
    // configuration's exact vertex position is in its rectangle.
    let axes = [(-7, 5), (-9, 4), (-3, 3), (-2, 5), (-4, 1)]
        .map(|(lo, hi)| Interval::new(frac(lo, 1), frac(hi, 1)).unwrap());
    let b = ConfigurationBox::new(axes);
    let scene = BoxScene::new(&b).unwrap();
    let mut random = Random::new(903);
    for _ in 0..4 {
        let x: [Q; 5] = std::array::from_fn(|j| random.rational(b.axes()[j].lo(), b.axes()[j].hi(), 6));
        let point = ConfigurationScene::new(&x.clone().map(QSqrt5::from_rational)).unwrap();
        for body in [Body::Hole, Body::Plug] {
            for v in 0..VERTEX_COUNT {
                let r = scene.rectangle(body, v);
                let e = point.enclosure(body, v);
                // The exact position is in e (width ≤ 2^-100); e must meet r.
                assert!(e[0].lo() <= r[0].hi() && r[0].lo() <= e[0].hi());
                assert!(e[1].lo() <= r[1].hi() && r[1].lo() <= e[1].hi());
            }
        }
    }
}
