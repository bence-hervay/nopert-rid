//! Adversarial tests of Global: an exact
//! reconstruction of the lemma's verdict on the identity zoom by a second
//! implementation (the RID from its defining formula, the plug turned by a
//! hand-written quaternion product, tensor Bernstein coefficients of degree
//! (1, 1, 2, 2, 2) obtained from exact values on a grid instead of from the
//! polynomial) for the maximum form (all 60 members), and an exact check of
//! every candidate the screen at the box's centre drops or keeps. The two
//! exact configurations of D that no valid-edge gap can exclude (the reason
//! Global uses the maximum form, which needs no valid edge) are excluded by
//! the maximum form on small enough boxes.
use super::*;
use crate::components::tests::fixtures::{bare, from_path, search, Random, SMALL_ROOTS};
use crate::problem::configuration::{self, AXES};
use crate::problem::geometry::{self, Point, VERTEX_COUNT};
use std::cmp::Ordering;

// ---- A second exact model ------------------------------------------------

fn n(a: i64, b: i64, d: i64) -> QSqrt5 {
    QSqrt5::new(frac(a, d), frac(b, d))
}

/// The RID from its defining formula (sign changes and cyclic permutations of
/// (1, 1, φ³), (φ², φ, 2φ), (2 + φ, 0, φ²)), in the geometry module's index
/// order, which is found by exact equality.
fn model_vertices() -> Vec<Point> {
    let bases = [
        [n(1, 0, 1), n(1, 0, 1), n(2, 1, 1)],
        [n(3, 1, 2), n(1, 1, 2), n(1, 1, 1)],
        [n(5, 1, 2), n(0, 0, 1), n(3, 1, 2)],
    ];
    let mut own: Vec<Point> = Vec::new();
    for base in &bases {
        for signs in 0..8 {
            let signed: Point = std::array::from_fn(|j| if signs & (1 << j) == 0 { base[j].clone() } else { -&base[j] });
            for shift in 0..3 {
                let v: Point = std::array::from_fn(|j| signed[(j + shift) % 3].clone());
                if !own.contains(&v) {
                    own.push(v);
                }
            }
        }
    }
    assert_eq!(own.len(), VERTEX_COUNT);
    let order: Vec<Point> = geometry::vertices()
        .iter()
        .map(|v| own.iter().find(|w| *w == v).expect("the geometry's vertex is a vertex of the formula").clone())
        .collect();
    order
}

fn vertices() -> &'static [Point] {
    static V: std::sync::OnceLock<Vec<Point>> = std::sync::OnceLock::new();
    V.get_or_init(model_vertices)
}

fn mul(a: &QSqrt5, b: &QSqrt5) -> QSqrt5 {
    a * b
}

fn add(a: &QSqrt5, b: &QSqrt5) -> QSqrt5 {
    a + b
}

fn sub(a: &QSqrt5, b: &QSqrt5) -> QSqrt5 {
    a - b
}

fn inner(a: &Point, b: &Point) -> QSqrt5 {
    add(&add(&mul(&a[0], &b[0]), &mul(&a[1], &b[1])), &mul(&a[2], &b[2]))
}

fn outer(a: &Point, b: &Point) -> Point {
    [
        sub(&mul(&a[1], &b[2]), &mul(&a[2], &b[1])),
        sub(&mul(&a[2], &b[0]), &mul(&a[0], &b[2])),
        sub(&mul(&a[0], &b[1]), &mul(&a[1], &b[0])),
    ]
}

/// The Hamilton product, written out coordinate by coordinate.
fn hamilton(a: &[QSqrt5; 4], b: &[QSqrt5; 4]) -> [QSqrt5; 4] {
    let t = |i: usize, j: usize| mul(&a[i], &b[j]);
    [
        sub(&sub(&sub(&t(0, 0), &t(1, 1)), &t(2, 2)), &t(3, 3)),
        sub(&add(&add(&t(0, 1), &t(1, 0)), &t(2, 3)), &t(3, 2)),
        add(&add(&sub(&t(0, 2), &t(1, 3)), &t(2, 0)), &t(3, 1)),
        add(&sub(&add(&t(0, 3), &t(1, 2)), &t(2, 1)), &t(3, 0)),
    ]
}

/// `(1, r)(0, p)(1, -r)`, which is `(1 + |r|²) R(r) p`.
fn turned(r: &Point, p: &Point) -> Point {
    let one = QSqrt5::one();
    let q = [one.clone(), r[0].clone(), r[1].clone(), r[2].clone()];
    let conj = [one, -&r[0], -&r[1], -&r[2]];
    let pure = [QSqrt5::zero(), p[0].clone(), p[1].clone(), p[2].clone()];
    let x = hamilton(&hamilton(&q, &pure), &conj);
    [x[1].clone(), x[2].clone(), x[3].clone()]
}

fn view(s: &Q, t: &Q) -> Point {
    [QSqrt5::from_rational(s.clone()), QSqrt5::from_rational(t.clone()), QSqrt5::one()]
}

/// The support normal `u × (v_to − v_from)`.
fn normal(u: &Point, from: usize, to: usize) -> Point {
    let v = vertices();
    let d: Point = std::array::from_fn(|j| sub(&v[to][j], &v[from][j]));
    outer(u, &d)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Support {
    Valid,
    /// The normal vanishes (the edge is seen end-on): every condition holds
    /// with equality, and every gap of the edge vanishes at this view.
    EndOn,
    Invalid,
}

fn support(u: &Point, from: usize, to: usize) -> Support {
    let nu = normal(u, from, to);
    if nu.iter().all(QSqrt5::is_zero) {
        return Support::EndOn;
    }
    let h = inner(&nu, &vertices()[from]);
    if vertices().iter().all(|w| inner(&nu, w) <= h) {
        Support::Valid
    } else {
        Support::Invalid
    }
}

/// The gap `(1 + |r|²) n·v_from − n·(1, r)(0, p)(1, −r)` at an exact point.
fn gap(x: &[Q; AXES], from: usize, to: usize, plug: usize) -> QSqrt5 {
    let u = view(&x[0], &x[1]);
    let r: Point = std::array::from_fn(|j| QSqrt5::from_rational(x[2 + j].clone()));
    let nu = normal(&u, from, to);
    let scale = add(&QSqrt5::one(), &inner(&r, &r));
    sub(&mul(&scale, &inner(&nu, &vertices()[from])), &inner(&nu, &turned(&r, &vertices()[plug])))
}

/// Degrees of the gap in each variable are at most (1, 1, 2, 2, 2): `n` is
/// linear in `(s, t)` without an `st` term, and `1 + |r|²` and the quaternion
/// conjugation are quadratic in each `rᵢ`.
const DEGREES: [usize; AXES] = [1, 1, 2, 2, 2];

/// The exact values of `f` on the grid of the box: {lo, hi} along a
/// degree-1 axis, {lo, mid, hi} along a degree-2 axis, row-major with
/// variable 0 slowest; each value is a vector (one entry per polynomial).
fn grid_values(b: &ConfigurationBox, f: impl Fn(&[Q; AXES]) -> Vec<QSqrt5>) -> (Vec<usize>, Vec<Vec<QSqrt5>>) {
    let axes = b.axes();
    let points = |j: usize| -> Vec<Q> {
        let (lo, hi) = (axes[j].lo().clone(), axes[j].hi().clone());
        if DEGREES[j] == 1 {
            vec![lo, hi]
        } else {
            let mid = (&lo + &hi) / frac(2, 1);
            vec![lo, mid, hi]
        }
    };
    let grid: Vec<Vec<Q>> = (0..AXES).map(points).collect();
    let sizes: Vec<usize> = grid.iter().map(Vec::len).collect();
    let total: usize = sizes.iter().product();
    let values = (0..total)
        .map(|mut k| {
            let mut index = [0usize; AXES];
            for j in (0..AXES).rev() {
                index[j] = k % sizes[j];
                k /= sizes[j];
            }
            let x: [Q; AXES] = std::array::from_fn(|j| grid[j][index[j]].clone());
            f(&x)
        })
        .collect();
    (sizes, values)
}

/// The tensor Bernstein coefficients at degrees `DEGREES` from the grid
/// values of one polynomial: along a degree-2 axis
/// `b₁ = 2f(½) − (f(0) + f(1))/2`.
fn bernstein_of(sizes: &[usize], mut values: Vec<QSqrt5>) -> Vec<QSqrt5> {
    let total = values.len();
    let mut stride = 1;
    for j in (0..AXES).rev() {
        if sizes[j] == 3 {
            for k in 0..total {
                if (k / stride) % 3 == 1 {
                    let (f0, f1, f2) = (&values[k - stride], &values[k], &values[k + stride]);
                    let half = QSqrt5::from_rational(frac(1, 2));
                    values[k] = sub(&add(f1, f1), &mul(&half, &add(f0, f2)));
                }
            }
        }
        stride *= sizes[j];
    }
    values
}

/// The tensor Bernstein coefficients of the gap on the box at degrees
/// `DEGREES`, from its exact values on the grid. If the certifier accepts
/// (at the polynomial's own degrees, which are at most these), these are all
/// negative, because raising the degree takes convex combinations of
/// coefficients.
fn bernstein(b: &ConfigurationBox, from: usize, to: usize, plug: usize) -> Vec<QSqrt5> {
    let (sizes, values) = grid_values(b, |x| vec![gap(x, from, to, plug)]);
    bernstein_of(&sizes, values.into_iter().map(|mut v| v.remove(0)).collect())
}

/// The second implementation's verdict on the record `(from → to, plug)` for
/// the box: valid (or end-on) at the four corner views and every Bernstein
/// coefficient of degree `DEGREES` negative.
fn independent_verdict(b: &ConfigurationBox, from: usize, to: usize, plug: usize) -> bool {
    let corners = b.view_corners();
    let views_ok = corners.iter().all(|u| support(u, from, to) != Support::Invalid);
    views_ok && bernstein(b, from, to, plug).iter().all(|c| c.sign() == Ordering::Less)
}

/// The members `(1 + |r|²) n·w − n·(1, r)(0, p)(1, −r)` of the maximum form,
/// one per hole vertex `w`, at an exact point with coordinates in Q(√5).
fn members_exact(x: &[QSqrt5; AXES], from: usize, to: usize, plug: usize) -> Vec<QSqrt5> {
    let u = [x[0].clone(), x[1].clone(), QSqrt5::one()];
    let r = [x[2].clone(), x[3].clone(), x[4].clone()];
    let nu = normal(&u, from, to);
    let scale = add(&QSqrt5::one(), &inner(&r, &r));
    let plug_term = inner(&nu, &turned(&r, &vertices()[plug]));
    vertices().iter().map(|w| sub(&mul(&scale, &inner(&nu, w)), &plug_term)).collect()
}

fn rational_point(x: &[Q; AXES]) -> [QSqrt5; AXES] {
    std::array::from_fn(|j| QSqrt5::from_rational(x[j].clone()))
}

/// The second implementation's verdict on the maximum form: every tensor
/// Bernstein coefficient of degree `DEGREES` of every member negative. (A
/// member nonnegative at a corner of the box refutes it at once: that value
/// is its corner coefficient.)
fn independent_maximum(b: &ConfigurationBox, from: usize, to: usize, plug: usize) -> bool {
    let corner_refutes = (0..1u8 << AXES).any(|m| {
        members_exact(&rational_point(&b.corner(m)), from, to, plug).iter().any(|v| v.sign() != Ordering::Less)
    });
    if corner_refutes {
        return false;
    }
    let (sizes, values) = grid_values(b, |x| members_exact(&rational_point(x), from, to, plug));
    (0..VERTEX_COUNT).all(|w| {
        let member: Vec<QSqrt5> = values.iter().map(|v| v[w].clone()).collect();
        bernstein_of(&sizes, member).iter().all(|c| c.sign() == Ordering::Less)
    })
}

/// Whether every member of the maximum form has the full degrees `DEGREES`
/// (then the certifier's coefficients are exactly those of the grid).
fn members_full_degree(from: usize, to: usize, plug: usize) -> bool {
    let (view, rotation) = configuration::coordinates();
    maximum(from, to, plug)
        .members(&view, &rotation)
        .unwrap()
        .iter()
        .all(|p| p.degrees().iter().zip(DEGREES).all(|(&d, e)| usize::from(d) == e))
}

/// Where the lemma refuses a candidate that the second model accepts, some
/// member has a lower degree (so the certifier works with other
/// coefficients).
fn explained_by_degree(b: &ConfigurationBox, from: usize, to: usize, plug: usize) -> bool {
    independent_maximum(b, from, to, plug) && !members_full_degree(from, to, plug)
}

// ---- Boxes ----------------------------------------------------------------

/// A path whose box has its lower `s` and `t` ends at 0 (the square view is a
/// corner view) and random rotation bits: the views of the root's faces are
/// planes of mirror symmetry, where supports tie.
fn corner_path(random: &mut Random, depth: usize, zero_s: bool, zero_t: bool) -> String {
    (0..depth)
        .map(|d| match d % AXES {
            0 if zero_s => '0',
            1 if zero_t => '0',
            _ => if random.below(2) == 0 { '0' } else { '1' },
        })
        .collect()
}

fn adversarial_boxes(seed: u64, random_count: usize, corner_count: usize) -> Vec<ConfigurationBox> {
    let mut random = Random(seed);
    let mut out = Vec::new();
    for _ in 0..random_count {
        let depth = 10 + random.below(36);
        out.push(from_path(&random.path(depth)));
    }
    for k in 0..corner_count {
        let depth = 10 + random.below(30);
        out.push(from_path(&corner_path(&mut random, depth, k % 3 != 1, k % 3 != 2)));
    }
    // Thin boxes: a random box made degenerate along random axes at a random
    // rational point of that axis.
    for _ in 0..corner_count {
        let depth = 12 + random.below(20);
        let b = from_path(&random.path(depth));
        let axes: [Interval; AXES] = std::array::from_fn(|j| {
            let a = &b.axes()[j];
            if random.below(3) == 0 {
                let k = random.below(9) as i64;
                let x = a.lo() + &((a.hi() - a.lo()) * frac(k, 8));
                Interval::point(x)
            } else {
                a.clone()
            }
        });
        out.push(ConfigurationBox::new(axes));
    }
    out
}

// ---- Tests ----------------------------------------------------------------

#[test]
fn the_second_model_is_the_rid() {
    let v = vertices();
    let four = QSqrt5::integer(4);
    let mut edges = Vec::new();
    for a in 0..VERTEX_COUNT {
        for b in a + 1..VERTEX_COUNT {
            let d: Point = std::array::from_fn(|j| sub(&v[b][j], &v[a][j]));
            if inner(&d, &d) == four {
                edges.push([a, b]);
            }
        }
    }
    assert_eq!(edges, geometry::edges());
    // The quaternion conjugation is (1 + |r|²) times a rotation: it keeps
    // lengths up to that factor.
    let r = [n(1, 0, 3), n(-2, 0, 7), n(0, 1, 5)];
    let scale = add(&QSqrt5::one(), &inner(&r, &r));
    for p in v {
        let x = turned(&r, p);
        assert_eq!(inner(&x, &x), mul(&mul(&scale, &scale), &inner(p, p)));
    }
}

/// Every candidate Global accepts is accepted by the second implementation,
/// and on members of full degree the two verdicts agree exactly: the first
/// kept candidates and random ones on random, search, corner and thin boxes.
#[test]
fn global_verdicts_agree_with_an_independent_exact_reconstruction() {
    let global = Global::new();
    let all = global.all();
    let mut boxes = adversarial_boxes(71, 24, 12);
    let (records, refused) = search(&bare(), SMALL_ROOTS[1], 200).unwrap();
    boxes.extend(records.iter().map(|(p, _)| from_path(p)).take(12));
    boxes.extend(refused.iter().map(|p| from_path(p)).take(6));
    let mut random = Random(72);
    let (mut accepted, mut refused_both, mut lower_degree) = (0, 0, 0);
    for b in &boxes {
        let mut claims: Vec<MaximumGap> = global.candidates(b).into_iter().take(3).collect();
        claims.extend((0..2).map(|_| all[random.below(all.len())].clone()));
        for m in claims {
            let (from, to, plug) = key(&m);
            let lemma = global.holds(b, &m).is_ok();
            let second = independent_maximum(b, from, to, plug);
            assert!(!lemma || second, "Global accepts ({from} → {to}, {plug}) on {b:?}, the second model refuses");
            if lemma {
                accepted += 1;
            } else if second {
                // Only possible when the certifier works at lower degrees.
                assert!(explained_by_degree(b, from, to, plug), "({from} → {to}, {plug}) on {b:?}");
                lower_degree += 1;
            } else {
                refused_both += 1;
            }
        }
    }
    assert!(accepted >= 20 && refused_both >= 20, "{accepted} {refused_both} {lower_degree}");
}

/// The check's witness is the proposed candidate, it survives
/// the record's canonical round trip unchanged, and the second model accepts
/// it.
#[test]
fn the_stored_witness_is_the_checked_witness() {
    let global = Global::new();
    let mut accepted = 0;
    for b in &adversarial_boxes(73, 30, 6) {
        let Some(m) = global.check(b) else { continue };
        accepted += 1;
        let first = global.candidates(b).into_iter().find(|c| global.holds(b, c).is_ok());
        assert_eq!(first.as_ref(), Some(&m));
        let data = RecordData::Global(m.clone());
        let text = serde_json::to_string(&data).unwrap();
        let read: RecordData = serde_json::from_str(&text).unwrap();
        assert_eq!((serde_json::to_string(&read).unwrap(), &read), (text, &data));
        global.holds(b, &m).unwrap();
        let (from, to, plug) = key(&m);
        assert!(independent_maximum(b, from, to, plug));
    }
    assert!(accepted >= 8, "{accepted}");
}

/// Wherever some candidate of the fixed list holds, the proposed one holds
/// too, on real boxes of a small search and on adversarial boxes (the
/// proposal is a heuristic; this measures that it loses nothing here).
fn check_the_proposal(b: &ConfigurationBox) -> (usize, usize) {
    let global = Global::new();
    let some = global.all().into_iter().find(|m| global.holds(b, m).is_ok());
    let proposed = global.check(b);
    if some.is_some() {
        assert!(proposed.is_some(), "a candidate holds on {b:?}, the proposal does not");
    }
    (usize::from(some.is_some()), usize::from(some.is_none()))
}

#[test]
fn the_proposal_holds_wherever_some_candidate_holds() {
    let (records, _) = search(&bare(), SMALL_ROOTS[1], 200).unwrap();
    let mut totals = (0, 0);
    for b in records.iter().map(|(p, _)| from_path(p)).take(2).chain(adversarial_boxes(74, 2, 2)) {
        let (h, n) = check_the_proposal(&b);
        totals = (totals.0 + h, totals.1 + n);
    }
    assert!(totals.0 > 0 && totals.1 > 0, "{totals:?}");
}

/// The same on many more boxes (a few minutes on four threads).
///
/// cargo test --release --lib components::global::tests::adversarial::heavy -- --ignored --nocapture --test-threads=4
#[test]
#[ignore]
fn heavy_the_proposal_holds_wherever_some_candidate_holds() {
    let (records, refused) = search(&bare(), SMALL_ROOTS[1], 2000).unwrap();
    let mut boxes = adversarial_boxes(75, 24, 24);
    boxes.extend(records.iter().map(|(p, _)| from_path(p)).step_by(7).take(40));
    boxes.extend(refused.iter().map(|p| from_path(p)).step_by(7).take(20));
    std::thread::scope(|scope| {
        for chunk in boxes.chunks(boxes.len().div_ceil(4)) {
            scope.spawn(move || {
                for b in chunk {
                    check_the_proposal(b);
                }
            });
        }
    });
}

/// The same agreement as above on many more boxes and the first ten kept
/// candidates (several minutes on four threads).
///
/// cargo test --release --lib components::global::tests::adversarial::heavy -- --ignored --nocapture --test-threads=4
#[test]
#[ignore]
fn heavy_global_verdicts_agree_with_an_independent_exact_reconstruction() {
    let boxes = adversarial_boxes(76, 120, 60);
    let counts = std::sync::Mutex::new([0usize; 3]);
    std::thread::scope(|scope| {
        for chunk in boxes.chunks(boxes.len().div_ceil(4)) {
            let counts = &counts;
            scope.spawn(move || {
                let global = Global::new();
                let mut local = [0usize; 3];
                for b in chunk {
                    for m in global.candidates(b).iter().take(10) {
                        let (from, to, plug) = key(m);
                        let lemma = global.holds(b, m).is_ok();
                        let second = independent_maximum(b, from, to, plug);
                        assert!(!lemma || second, "{m:?} on {b:?}");
                        if !lemma && second {
                            assert!(explained_by_degree(b, from, to, plug), "{m:?} on {b:?}");
                        }
                        local[usize::from(lemma) + usize::from(second)] += 1;
                    }
                }
                let mut counts = counts.lock().unwrap();
                for j in 0..3 {
                    counts[j] += local[j];
                }
            });
        }
    });
    println!("refused by both, by the lemma only, accepted by both: {:?}", counts.lock().unwrap());
}

/// Search boxes whose sibling holds a record: the record moved to the sibling,
/// to the parent and to each child is decided by the second model exactly as
/// by the lemma (on full-degree polynomials), and never accepted alone.
#[test]
fn moved_records_are_decided_like_the_second_model() {
    let global = Global::new();
    for root in SMALL_ROOTS {
        let (records, _) = search(&bare(), root, 200).unwrap();
        for (path, data) in records.iter().filter(|(_, d)| d.component() == ComponentName::Global).take(3) {
            let (from, to, plug) = crate::components::tests::parts(data);
            let parent = &path[..path.len() - 1];
            let sibling = format!("{parent}{}", if path.ends_with('0') { '1' } else { '0' });
            for place in [parent.to_owned(), sibling, format!("{path}0"), format!("{path}1")] {
                let b = from_path(&place);
                let lemma = global.holds(&b, &maximum(from, to, plug)).is_ok();
                let second = independent_maximum(&b, from, to, plug);
                assert!(!lemma || second, "{data:?} at {place}");
                if !lemma && second {
                    assert!(explained_by_degree(&b, from, to, plug), "{data:?} at {place}");
                }
            }
        }
    }
}

/// Diagnosis (ignored): for the box of `RID_PATH`, the proposal checked
/// against the fixed list, and every one of the 900 fixed candidates
/// decided by the lemma and by the second model.
///
/// RID_PATH=<path> cargo test --release --lib components::global::tests::adversarial::every_candidate_on_one_box -- --ignored --nocapture
#[test]
#[ignore]
fn every_candidate_on_one_box() {
    let b = from_path(&std::env::var("RID_PATH").unwrap());
    println!("some candidate holds, none holds: {:?}", check_the_proposal(&b));
    let global = Global::new();
    let mut eliminating = Vec::new();
    for m in global.all() {
        let (from, to, plug) = key(&m);
        let (lemma, second) = (global.holds(&b, &m).is_ok(), independent_maximum(&b, from, to, plug));
        if lemma || second {
            eliminating.push((from, to, plug, lemma, second));
        }
    }
    println!("candidates accepted by the lemma or the second model: {eliminating:?}");
    println!("check: {:?}", global.check(&b));
}

/// Boxes whose view rectangle contains the corner view `(s, t) = (1/φ, 0)` of
/// D (a five-fold axis, where every side of the hole's decagonal shadow is a
/// face seen edge-on, with four tied vertices) keep only the edges of the two
/// squares with normals ±e₂ valid at their four corner views, since every
/// other side's support changes across a line through that view. Two such
/// boxes, found unresolved (before the maximum form) by searches with the
/// cover catalogue at depths 46 and 62 (rotations near (0.28, 0.08, −0.05)
/// and near (0, −0.025, 0)): the second model finds no shortcut record among
/// all 14,400, and the shortcut's ranking proposes nothing the lemma
/// accepts. The maximum form needs no valid edge (see the last test).
#[test]
fn boxes_containing_the_five_fold_corner_view_admit_no_shortcut_record() {
    for path in [
        "1011010101100010011110110100000000010000000001",
        "10100100111001100011100111000100001100010000110001000011000100",
    ] {
        let b = from_path(path);
        let s = b.axes()[0].clone();
        // 1/φ = (√5 − 1)/2 lies strictly inside the s-interval, and t starts at 0.
        let inverse_phi = n(-1, 1, 2);
        assert!(QSqrt5::from_rational(s.lo().clone()) < inverse_phi && inverse_phi < QSqrt5::from_rational(s.hi().clone()));
        assert!(b.axes()[1].lo().is_zero());
        let views = b.view_corners();
        let mut valid = Vec::new();
        for &[x, y] in geometry::edges() {
            for (from, to) in [(x, y), (y, x)] {
                if views.iter().all(|u| support(u, from, to) != Support::Invalid) {
                    valid.push((from, to));
                    for plug in 0..VERTEX_COUNT {
                        // A nonnegative corner value is a nonnegative corner
                        // coefficient; otherwise all coefficients are computed.
                        let corner_refutes = (0..1u8 << AXES).any(|m| gap(&b.corner(m), from, to, plug).sign() != Ordering::Less);
                        assert!(corner_refutes || !independent_verdict(&b, from, to, plug), "({from} → {to}, {plug}) on {path}");
                    }
                }
            }
        }
        // Only the edges of the squares with normals ±e₂ = (0, ±1, 0).
        assert_eq!(valid.len(), 4, "{valid:?}");
        let top = vertices().iter().map(|v| v[1].clone()).max().unwrap();
        for (from, to) in valid {
            let (a, c) = (&vertices()[from], &vertices()[to]);
            assert!(a[1] == c[1] && (a[1] == top || a[1] == -&top), "{from} → {to}");
        }
    }
}

/// The same along the triangle wall `φs + φ²t = 1` (a mirror, irrational, so
/// every box meeting it straddles it): the two squares whose normal is
/// perpendicular to the mirror are seen edge-on on the whole wall, so their
/// four sides never have an edge valid at all four corner views. The box
/// below contains the configuration `s = 0.3044` on the wall,
/// `r = (−0.1087, 0.1244, 0.0028)`, which lies in D and whose plug sticks out
/// by 0.224 through those sides and nowhere else (every other side keeps it
/// inside by at least 0.0046, in floating point); no shortcut record exists
/// for any box containing it (before the maximum form, its subtree had 6,656
/// unresolved boxes at depth 65).
#[test]
fn boxes_containing_a_wall_configuration_hidden_behind_the_mirror_squares_admit_no_shortcut_record() {
    let path = "00011111001101011100011001111000010001101001110011";
    let b = from_path(path);
    // The wall crosses the view rectangle: φs + φ²t − 1 changes sign at its corners.
    let phi = n(1, 1, 2);
    let wall: Vec<Ordering> = b
        .view_corners()
        .iter()
        .map(|u| sub(&add(&mul(&phi, &u[0]), &mul(&mul(&phi, &phi), &u[1])), &QSqrt5::one()).sign())
        .collect();
    assert!(wall.contains(&Ordering::Less) && wall.contains(&Ordering::Greater), "{wall:?}");
    let views = b.view_corners();
    let mut valid = 0;
    for &[x, y] in geometry::edges() {
        for (from, to) in [(x, y), (y, x)] {
            if views.iter().all(|u| support(u, from, to) != Support::Invalid) {
                valid += 1;
                for plug in 0..VERTEX_COUNT {
                    let corner_refutes = (0..1u8 << AXES).any(|m| gap(&b.corner(m), from, to, plug).sign() != Ordering::Less);
                    assert!(corner_refutes || !independent_verdict(&b, from, to, plug), "({from} → {to}, {plug})");
                }
            }
        }
    }
    assert_eq!(valid, 12);
    assert!(crate::components::domain::Domain.check(&b).is_none());
}

/// The gap at an exact configuration with coordinates in Q(√5).
fn gap_exact(x: &[QSqrt5; AXES], from: usize, to: usize, plug: usize) -> QSqrt5 {
    let u = [x[0].clone(), x[1].clone(), QSqrt5::one()];
    let r = [x[2].clone(), x[3].clone(), x[4].clone()];
    let nu = normal(&u, from, to);
    let scale = add(&QSqrt5::one(), &inner(&r, &r));
    sub(&mul(&scale, &inner(&nu, &vertices()[from])), &inner(&nu, &turned(&r, &vertices()[plug])))
}

/// Exactly: `x` lies in the box and in D (all 73 inequalities `c ≤ 0`), and
/// every gap of every ordered edge valid at the box's four corner views is
/// nonnegative at `x`. Then no shortcut record is negative on any box
/// containing `x` (a certified record is negative at every point), however
/// small.
fn hides_from_the_shortcut(b: &ConfigurationBox, x: &[QSqrt5; AXES]) -> bool {
    let views = b.view_corners();
    let hidden = geometry::edges().iter().flat_map(|&[p, q]| [(p, q), (q, p)]).all(|(from, to)| {
        views.iter().any(|u| support(u, from, to) == Support::Invalid)
            || (0..VERTEX_COUNT).all(|plug| gap_exact(x, from, to, plug).sign() != Ordering::Less)
    });
    contains(b, x) && in_d(x) && hidden
}

fn contains(b: &ConfigurationBox, x: &[QSqrt5; AXES]) -> bool {
    (0..AXES).all(|j| {
        let a = &b.axes()[j];
        QSqrt5::from_rational(a.lo().clone()) <= x[j] && x[j] <= QSqrt5::from_rational(a.hi().clone())
    })
}

fn in_d(x: &[QSqrt5; AXES]) -> bool {
    crate::problem::domain::affine_polynomials().iter().all(|p| p.evaluate(x).sign() != Ordering::Greater)
}

/// The two exact configurations of D at the five-fold
/// corner view and on the triangle wall, with the paths of boxes containing
/// them.
fn stuck_configurations() -> [([QSqrt5; AXES], &'static str); 2] {
    let q = |a: i64, d: i64| QSqrt5::from_rational(frac(a, d));
    // The five-fold corner view (1/φ, 0), tilted about e₂: the plug's extent
    // along ±e₂ equals the hole's, and every other side switches at this view.
    let inverse_phi = n(-1, 1, 2);
    let corner = [inverse_phi, QSqrt5::zero(), QSqrt5::zero(), q(-249, 10000), QSqrt5::zero()];
    // On the wall φs + φ²t = 1 at s = 761/2500: t = (1 − φs)/φ² = (1 − φs)(3 − √5)/2.
    let s = q(761, 2500);
    let phi = n(1, 1, 2);
    let t = mul(&sub(&QSqrt5::one(), &mul(&phi, &s)), &n(3, -1, 2));
    let wall = [s, t, q(-1087, 10000), q(1244, 10000), q(28, 10000)];
    [
        (corner, "10100100111001100011100111000100001100010000110001000011000100"),
        (wall, "00011111001101011100011001111000010001101001110011"),
    ]
}

/// The child of `path` whose closed box contains `x` (the lower one on a tie).
fn child_containing(path: &str, x: &[QSqrt5; AXES]) -> String {
    let lower = format!("{path}0");
    if contains(&from_path(&lower), x) {
        lower
    } else {
        format!("{path}1")
    }
}

/// Both configurations hide from every shortcut record, on their boxes and
/// five levels deeper: no gap of an edge valid at the corner views is
/// negative there. Before the maximum form, no search could terminate
/// around them.
#[test]
fn exact_configurations_of_d_that_no_shortcut_record_can_exclude() {
    for (x, path) in stuck_configurations() {
        assert!(hides_from_the_shortcut(&from_path(path), &x), "{path}");
        let mut path = path.to_owned();
        for _ in 0..5 {
            path = child_containing(&path, &x);
            let deeper = from_path(&path);
            assert!(hides_from_the_shortcut(&deeper, &x), "{path}");
        }
    }
}

/// The maximum form excludes both configurations: on the boxes containing
/// them (within a few levels; in fact at once) Global's check finds a
/// witness, which the second model confirms on the box, and whose maximum
/// `G = max_w (1 + |r|²) n·w − n·R̂p` is exactly negative at the
/// configuration (the plug sticks out through a side whose edge is not valid
/// at every corner view). The witness stays valid on deeper boxes.
#[test]
fn the_maximum_form_excludes_the_configurations_on_small_enough_boxes() {
    let global = Global::new();
    for (x, start) in stuck_configurations() {
        let mut path = start.to_owned();
        let found = loop {
            assert!(path.len() < start.len() + 40, "no Global witness within 40 levels below {start}");
            let b = from_path(&path);
            if let Some(m) = global.check(&b) {
                break (path, m);
            }
            path = child_containing(&path, &x);
        };
        let (path, m) = found;
        let (from, to, plug) = key(&m);
        let b = from_path(&path);
        global.holds(&b, &m).unwrap();
        assert!(independent_maximum(&b, from, to, plug), "{path}: {m:?}");
        assert_eq!(members_exact(&x, from, to, plug).into_iter().max().unwrap().sign(), Ordering::Less);
        println!("{start}: the witness {m:?} at depth {} ({} levels deeper)", path.len(), path.len() - start.len());
        let deeper = child_containing(&child_containing(&path, &x), &x);
        global.holds(&from_path(&deeper), &m).unwrap();
    }
}
