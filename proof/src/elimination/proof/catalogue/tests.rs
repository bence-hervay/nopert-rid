//! Tests of the catalogue: pins against data, the Local table, the signs
//! of the arc covers, and the exact statement that the view sets of the
//! Local, square and pentagon covers, shrunk by a positive margin, still
//! cover D's view triangle (every aligned configuration of D has a
//! neighbourhood inside one cover). Loading every cover is tested by the
//! proof module's tests.
use super::*;
use crate::elimination::proof::format::CoverFile;
use crate::elimination::zoom::cover::ZoomCover;
use std::cmp::Ordering;

#[test]
fn pins_and_data_agree() {
    for name in EXOTIC {
        let pin = exotic_pin(name).unwrap();
        let file = CoverFile::parse(exotic_data(name).unwrap()).unwrap();
        assert_eq!((pin.name.as_str(), pin.scope), (file.name.as_str(), file.scope));
        assert_eq!(pin.parameters, file.parameters, "{name}");
    }
    for number in 0..LOCAL_COUNT {
        let pin = local_pin(number).unwrap();
        let file = CoverFile::parse(local_data(number).unwrap()).unwrap();
        assert_eq!((pin.name.as_str(), pin.scope), (file.name.as_str(), file.scope));
        assert_eq!(pin.name, number.to_string());
        assert_eq!(pin.parameters, file.parameters, "{number}");
    }
    assert!(exotic_pin("wall").is_none() && exotic_data("wall").is_none());
    assert!(exotic_pin("arc").is_none() && exotic_pin("second-arc").is_none());
    assert!(local_pin(LOCAL_COUNT).is_none() && local_data(LOCAL_COUNT).is_none());
}

fn rows() -> Vec<[Q; 5]> {
    LOCAL.iter().map(|row| row.map(|x| parse_rational(x).unwrap())).collect()
}

#[test]
fn local_rows_tile_without_overlap_and_the_margin_is_applied() {
    let rows = rows();
    for [s0, s1, t0, t1, radius] in &rows {
        assert!(&q(0) <= s0 && s0 < s1 && s1 <= &frac(2, 3));
        assert!(&q(0) <= t0 && t0 < t1 && t1 <= &frac(2, 5));
        assert!(radius > &q(0) && radius <= &frac(1, 20));
    }
    for (i, a) in rows.iter().enumerate() {
        for b in &rows[i + 1..] {
            let apart = a[1] <= b[0] || b[1] <= a[0] || a[3] <= b[2] || b[3] <= a[2];
            assert!(apart, "{a:?} {b:?}");
        }
    }
    let m = frac(MARGIN.0, MARGIN.1);
    assert!(m > q(0));
    for (n, [s0, s1, t0, t1, _]) in rows.iter().enumerate() {
        let [s, t] = local_rectangle(n).unwrap();
        let expect = |x: Q, lo: bool, end: Q| if lo { (x - &m).max(end) } else { (x + &m).min(end) };
        assert_eq!(s[0], expect(s0.clone(), true, q(0)));
        assert_eq!(s[1], expect(s1.clone(), false, frac(2, 3)));
        assert_eq!(t[0], expect(t0.clone(), true, q(0)));
        assert_eq!(t[1], expect(t1.clone(), false, frac(2, 5)));
    }
}

#[test]
fn the_arc_cover_labels_are_the_signs_of_their_configurations() {
    // r₁ is affine in the cover coordinates, so its extremes over a cover's
    // box are at the box's corners: positive for the `+` covers (near the
    // plane r ≈ +r_A), negative for the `-` covers (r ≈ −r_A).
    for name in EXOTIC.iter().filter(|n| n.ends_with(['+', '-'])) {
        let parameters = exotic_pin(name).unwrap().parameters;
        let sign = if name.ends_with('+') { Ordering::Greater } else { Ordering::Less };
        let expected = if sign == Ordering::Greater { Sign::Plus } else { Sign::Minus };
        assert_eq!(parameters.coordinates, Coordinates::ArcPlane(expected), "{name}");
        let cover = ZoomCover::new(parameters.clone()).unwrap();
        for mask in 0..32u32 {
            let z: [QSqrt5; 5] = std::array::from_fn(|i| cover.bounds()[i][(mask >> i & 1) as usize].clone());
            let x: [QSqrt5; 5] = std::array::from_fn(|i| {
                (0..5).fold(parameters.centre[i].0.clone(), |s, j| &s + &(&parameters.map[i][j].0 * &z[j]))
            });
            let constants = x.clone().map(crate::arithmetic::polynomial::Polynomial::constant);
            let (_, rotation) = parameters.coordinates.configuration(&constants);
            let r1 = rotation[0].evaluate(&x.clone().map(|_| QSqrt5::zero()));
            assert_eq!(r1.sign(), sign, "{name}");
        }
    }
}

// ---- The covers around the aligned set form an open cover ------------------

/// A half-plane `α s + β t ≤ γ` of views.
struct HalfPlane {
    alpha: QSqrt5,
    beta: QSqrt5,
    gamma: QSqrt5,
}

/// The view set of a cover in configuration coordinates whose first two
/// cover coordinates depend on the view `(s, t)` only and the last three on
/// the rotation only: its half-planes, each shrunk by `delta` in the
/// maximum norm unless it is a side of the viewing rectangle
/// `V = [0, 2/3] × [0, 2/5]` (which a neighbourhood relative to V never
/// crosses), and the rotation radius it contains around `r = 0`.
fn shrunk_view_set(cover: &ZoomCover, delta: &Q) -> (Vec<HalfPlane>, Q) {
    let p = cover.parameters();
    assert_eq!(p.coordinates, Coordinates::Configuration);
    let zero: [QSqrt5; 5] = std::array::from_fn(|_| QSqrt5::zero());
    let unit = |j: usize| -> [QSqrt5; 5] { std::array::from_fn(|i| if i == j { QSqrt5::one() } else { QSqrt5::zero() }) };
    let z0 = cover.cover_coordinates(&zero);
    let dz: Vec<[QSqrt5; 5]> = (0..5).map(|j| {
        let z = cover.cover_coordinates(&unit(j));
        std::array::from_fn(|i| &z[i] - &z0[i])
    }).collect();
    // Block structure: z₀, z₁ from (s, t); z₂..z₄ from r (here the identity).
    for i in 0..5 {
        for j in 0..5 {
            if (i < 2) != (j < 2) {
                assert!(dz[j][i].is_zero(), "coordinate {i} depends on variable {j}");
            }
        }
    }
    let mut radius: Option<Q> = None;
    for i in 2..5 {
        assert_eq!(dz[i][i], QSqrt5::one());
        assert!((2..5).all(|j| j == i || dz[j][i].is_zero()));
        let [lo, hi] = &cover.bounds()[i];
        assert!(lo.sqrt5_part().is_zero() && hi.sqrt5_part().is_zero());
        let r = std::cmp::min(-lo.rational_part().clone(), hi.rational_part().clone());
        radius = Some(radius.map_or(r.clone(), |x| std::cmp::min(x, r)));
    }
    let radius = radius.unwrap();
    assert!(radius > q(0));
    let beyond = p.beyond.map(|b| b.axis);
    let side_of_v = |alpha: &QSqrt5, beta: &QSqrt5, gamma: &QSqrt5| {
        // −s ≤ 0, −t ≤ 0, s ≤ 2/3, t ≤ 2/5 up to a positive factor.
        let rational = |x: &QSqrt5| x.sqrt5_part().is_zero();
        if !(rational(alpha) && rational(beta) && rational(gamma)) {
            return false;
        }
        let (a, b, g) = (alpha.rational_part(), beta.rational_part(), gamma.rational_part());
        (b.is_zero() && ((a < &q(0) && g.is_zero()) || (a > &q(0) && g == &(a * &frac(2, 3)))))
            || (a.is_zero() && ((b < &q(0) && g.is_zero()) || (b > &q(0) && g == &(b * &frac(2, 5)))))
    };
    let mut planes = Vec::new();
    for i in 0..2 {
        let [lo, hi] = &cover.bounds()[i];
        // z_i = z0_i + dz[0][i] s + dz[1][i] t.
        let (a, b) = (dz[0][i].clone(), dz[1][i].clone());
        let mut sides = vec![(-&a, -&b, &z0[i] - lo)];
        if beyond != Some(i) {
            sides.push((a, b, hi - &z0[i]));
        }
        for (alpha, beta, gamma) in sides {
            let shrink = if side_of_v(&alpha, &beta, &gamma) {
                QSqrt5::zero()
            } else {
                let norm = |x: &QSqrt5| if x.sign() == Ordering::Less { -x } else { x.clone() };
                (&norm(&alpha) + &norm(&beta)).scale(delta)
            };
            planes.push(HalfPlane { gamma: &gamma - &shrink, alpha, beta });
        }
    }
    (planes, radius)
}

fn inside(planes: &[HalfPlane], s: &Q, t: &Q) -> bool {
    planes.iter().all(|h| &h.alpha.scale(s) + &h.beta.scale(t) <= h.gamma)
}

/// Covers the part of D's view triangle `{s, t ≥ 0, φs + φ²t ≤ 1}` in the
/// box `[s₀, s₁] × [t₀, t₁]` by the shrunk sets, bisecting at most `depth`
/// more times; returns the number of boxes used, or `None` at a point
/// covered by none.
fn cover_triangle(sets: &[Vec<HalfPlane>], s: [Q; 2], t: [Q; 2], depth: u32) -> Option<usize> {
    let phi = QSqrt5::new(frac(1, 2), frac(1, 2));
    let phi2 = &phi * &phi;
    // The box misses the triangle when its lowest corner lies beyond the wall.
    if &phi.scale(&s[0]) + &phi2.scale(&t[0]) > QSqrt5::one() {
        return Some(0);
    }
    let corners_in = |planes: &Vec<HalfPlane>| {
        [(&s[0], &t[0]), (&s[0], &t[1]), (&s[1], &t[0]), (&s[1], &t[1])].iter().all(|(a, b)| inside(planes, a, b))
    };
    if sets.iter().any(corners_in) {
        return Some(1);
    }
    if depth == 0 {
        return None;
    }
    let (sm, tm) = ((&s[0] + &s[1]) * frac(1, 2), (&t[0] + &t[1]) * frac(1, 2));
    let mut used = 0;
    for (s, t) in [
        ([s[0].clone(), sm.clone()], [t[0].clone(), tm.clone()]),
        ([s[0].clone(), sm.clone()], [tm.clone(), t[1].clone()]),
        ([sm.clone(), s[1].clone()], [t[0].clone(), tm.clone()]),
        ([sm, s[1].clone()], [tm, t[1].clone()]),
    ] {
        used += cover_triangle(sets, s, t, depth - 1)?;
    }
    Some(used)
}

#[test]
fn shrunk_view_sets_of_the_aligned_covers_cover_the_view_triangle() {
    // Every aligned configuration (p, 0) of D, p in the triangle, lies with
    // (B∞(p, δ) ∩ V) × [−ρ, ρ]³ in one cover: the covers are an open cover
    // of the aligned part of D relative to B₀, with Lebesgue number at
    // least min(δ, ρ). A quarter of the margin is kept as overlap.
    let delta = frac(MARGIN.0, 4 * MARGIN.1);
    let mut covers: Vec<ZoomCover> =
        ["square", "pentagon"].iter().map(|n| ZoomCover::new(exotic_pin(n).unwrap().parameters).unwrap()).collect();
    covers.extend((0..LOCAL_COUNT).map(|n| ZoomCover::new(local_pin(n).unwrap().parameters).unwrap()));
    let mut sets = Vec::new();
    let mut radius = q(1);
    for cover in &covers {
        let (planes, r) = shrunk_view_set(cover, &delta);
        radius = std::cmp::min(radius, r);
        sets.push(planes);
    }
    assert_eq!(radius, frac(1, 200));
    let used = cover_triangle(&sets, [q(0), frac(2, 3)], [q(0), frac(2, 5)], 14).expect("the shrunk sets cover the triangle");
    assert!(used > 0);
    // The control: with the rows themselves (no margin) the Local covers
    // only abut, and their shrunk sets leave the seams uncovered.
    let mut abutting = Vec::new();
    for cover in &covers[..2] {
        abutting.push(shrunk_view_set(cover, &delta).0);
    }
    for (n, [s0, s1, t0, t1, _]) in rows().into_iter().enumerate() {
        let mut p = local_pin(n).unwrap().parameters;
        let crate::elimination::zoom::cover::Shape::Tube { base, .. } = &mut p.shape else { panic!() };
        *base = vec![[Rational(s0), Rational(s1)], [Rational(t0), Rational(t1)]];
        abutting.push(shrunk_view_set(&ZoomCover::new(p).unwrap(), &delta).0);
    }
    assert!(cover_triangle(&abutting, [q(0), frac(2, 3)], [q(0), frac(2, 5)], 14).is_none());
}

#[test]
fn loading_names_the_failing_cover() {
    let stop = AtomicBool::new(true);
    let threads = NonZeroUsize::new(2).unwrap();
    let error = exotic(threads, &stop).unwrap_err();
    assert_eq!(error.cover, "square");
    assert!(matches!(error.error, ProofError::Interrupted));
    assert_eq!(local(threads, &stop).unwrap_err().cover, "0");
}

/// Loading times on 1 and 4 threads:
/// `cargo test --release --lib catalogue::tests::measure_loading -- --ignored --nocapture --test-threads=1`
#[test]
#[ignore]
fn measure_loading() {
    let stop = AtomicBool::new(false);
    for n in [1, 4] {
        let threads = NonZeroUsize::new(n).unwrap();
        for name in EXOTIC {
            let start = std::time::Instant::now();
            let cover = ProvedCover::load(&exotic_pin(name).unwrap(), exotic_data(name).unwrap(), threads, &stop).unwrap();
            eprintln!("{n} threads: {name}: {:?} in {:.2} s", cover.report(), start.elapsed().as_secs_f64());
        }
        let start = std::time::Instant::now();
        local(threads, &stop).unwrap();
        eprintln!("{n} threads: {LOCAL_COUNT} Local covers in {:.2} s", start.elapsed().as_secs_f64());
    }
}
