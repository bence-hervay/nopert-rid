use super::*;
use crate::arithmetic::exact::{frac, q, Interval, Q};
use crate::components::record::{ComponentName, RecordData};
use crate::components::tests::fixtures::{bare, from_path, search, two_to_minus, Random, SMALL_ROOTS};
use crate::components::tests::oracle::{self, Configuration};
use crate::components::tests::{confirm, oracle_verdict, parts};
use crate::problem::geometry::Point;
use crate::problem::configuration::{self, AXES};
use std::cmp::Ordering;
use std::collections::BTreeSet;

fn maximum(from: usize, to: usize, plug: usize) -> MaximumGap {
    MaximumGap::new(Edge::new(from, to).unwrap(), plug).unwrap()
}

fn key(m: &MaximumGap) -> (usize, usize, usize) {
    (m.edge().from(), m.edge().to(), m.plug())
}

/// The index of `-v_i`, found exactly.
fn antipode(i: usize) -> usize {
    let v = &geometry::vertices()[i];
    let minus = [-&v[0], -&v[1], -&v[2]];
    geometry::vertices().iter().position(|w| *w == minus).unwrap()
}

/// Random boxes and the boxes met by the search of a small subtree.
fn boxes(count: usize, seed: u64) -> Vec<ConfigurationBox> {
    let mut random = Random(seed);
    let mut out: Vec<ConfigurationBox> = (0..count)
        .map(|_| {
            let depth = 12 + random.below(34);
            from_path(&random.path(depth))
        })
        .collect();
    let (records, refused) = search(&bare(), SMALL_ROOTS[0], 200).unwrap();
    out.extend(records.iter().map(|(path, _)| from_path(path)));
    out.extend(refused.iter().map(|path| from_path(path)));
    out
}

/// The members of a candidate at an exact configuration, by the polynomial.
fn member_values(m: &MaximumGap, x: &[QSqrt5; AXES]) -> Vec<QSqrt5> {
    let (view, rotation) = configuration::coordinates();
    m.members(&view, &rotation).unwrap().iter().map(|p| p.evaluate(x)).collect()
}

fn centre(b: &ConfigurationBox) -> [QSqrt5; AXES] {
    b.midpoint().map(QSqrt5::from_rational)
}

#[test]
fn the_fixed_list_is_fifteen_edge_directions_times_sixty_plug_vertices() {
    let global = Global::new();
    let all = global.all();
    assert_eq!(all.len(), 15 * VERTEX_COUNT);
    let directions: Vec<Edge> = all.iter().step_by(VERTEX_COUNT).map(|m| m.edge().clone()).collect();
    assert_eq!(directions, global.directions);
    for (k, m) in all.iter().enumerate() {
        assert_eq!((m.edge(), m.plug()), (&directions[k / VERTEX_COUNT], k % VERTEX_COUNT));
    }
    // Every edge of the RID is parallel to exactly one direction, eight edges
    // per direction, and each direction is the first edge of its class.
    let parallel = |a: &Point, b: &Point| geometry::cross(a, b).iter().all(QSqrt5::is_zero);
    let mut counts = [0; 15];
    for (position, &[a, b]) in geometry::edges().iter().enumerate() {
        let d = Edge::new(a, b).unwrap().direction();
        let classes: Vec<usize> = (0..15).filter(|&k| parallel(&d, &directions[k].direction())).collect();
        assert_eq!(classes.len(), 1, "edge {a} → {b}");
        counts[classes[0]] += 1;
        let first = geometry::edges().iter().position(|&[x, y]| parallel(&Edge::new(x, y).unwrap().direction(), &d)).unwrap();
        assert!(first <= position);
        if first == position {
            assert_eq!(directions[classes[0]], Edge::new(a, b).unwrap());
        }
    }
    assert_eq!(counts, [8; 15]);
    // The directions are in the order of their first edges in the geometry's list.
    let positions: Vec<usize> = directions
        .iter()
        .map(|d| geometry::edges().iter().position(|&[x, y]| [x, y] == [d.from(), d.to()]).unwrap())
        .collect();
    assert!(positions.windows(2).all(|w| w[0] < w[1]), "{positions:?}");
}

/// The fixed list needs one edge per direction and one sign: the members of
/// the maximum form depend only on the edge vector `d`, and the candidates
/// `(d, p)` and `(-d, -p)` have the same members (in another order), while
/// `(-d, p)` differs.
#[test]
fn parallel_edges_and_antipodal_candidates_state_the_same_maximum_form() {
    let (view, rotation) = configuration::coordinates();
    // The members as a multiset: members of hole vertices that differ by a
    // multiple of `d` are equal.
    let members = |from: usize, to: usize, plug: usize| -> Vec<String> {
        let mut all: Vec<String> = maximum(from, to, plug).members(&view, &rotation).unwrap().iter().map(|p| format!("{p:?}")).collect();
        all.sort();
        all
    };
    let mut random = Random(40);
    let global = Global::new();
    for direction in &global.directions {
        let d = direction.direction();
        for _ in 0..3 {
            let plug = random.below(VERTEX_COUNT);
            let own = members(direction.from(), direction.to(), plug);
            assert_eq!(own.len(), VERTEX_COUNT);
            for &[a, b] in geometry::edges() {
                let e = Edge::new(a, b).unwrap().direction();
                if !geometry::cross(&d, &e).iter().all(QSqrt5::is_zero) {
                    continue;
                }
                let (same, opposite) = if e == d { ((a, b), (b, a)) } else { ((b, a), (a, b)) };
                assert_eq!(members(same.0, same.1, plug), own);
                assert_eq!(members(opposite.0, opposite.1, antipode(plug)), own);
                assert_ne!(members(opposite.0, opposite.1, plug), own);
            }
        }
    }
}

#[test]
fn check_and_holds_agree_and_the_oracle_confirms_every_elimination() {
    let global = Global::new();
    let mut eliminated = 0;
    for (n, b) in boxes(100, 41).iter().enumerate() {
        let Some(m) = global.check(b) else { continue };
        eliminated += 1;
        global.holds(b, &m).unwrap();
        confirm(b, &RecordData::Global(m), n as u64);
    }
    assert!(eliminated >= 30, "{eliminated}");
}

/// The proposal is at most one candidate of the fixed list, the same on
/// every call, and never refuted at the box's exact centre when it holds.
#[test]
fn the_proposal_is_one_candidate_of_the_fixed_list() {
    let global = Global::new();
    let all: BTreeSet<(usize, usize, usize)> = global.all().iter().map(key).collect();
    let (mut proposed, mut held) = (0, 0);
    for b in boxes(30, 43).iter() {
        let candidates = global.candidates(b);
        assert!(candidates.len() <= 1, "{candidates:?}");
        assert_eq!(candidates, global.candidates(b));
        for m in &candidates {
            assert!(all.contains(&key(m)));
            proposed += 1;
            if global.holds(b, m).is_ok() {
                held += 1;
                assert!(member_values(m, &centre(b)).iter().all(|v| v.sign() == Ordering::Less));
            }
        }
    }
    assert!(proposed >= 10 && held >= 5, "{proposed} {held}");
}

/// `check(b)` is the first candidate of the whole fixed list that holds: no
/// earlier candidate of the fixed list holds (so the screen dropped none
/// that would), and when it finds none, no candidate at all holds.
#[test]
fn check_is_the_first_candidate_of_the_fixed_list_that_holds() {
    let global = Global::new();
    let all = global.all();
    let (records, refused) = search(&bare(), SMALL_ROOTS[2], 200).unwrap();
    let mut found = 0;
    for (path, _) in records.iter().filter(|(_, d)| d.component() == ComponentName::Global).take(4) {
        let b = from_path(path);
        let m = global.check(&b).unwrap();
        let position = all.iter().position(|x| *x == m).unwrap();
        assert!(all[..position].iter().all(|x| global.holds(&b, x).is_err()), "{path}");
        found += 1;
    }
    for path in refused.iter().take(2) {
        let b = from_path(path);
        assert_eq!(global.check(&b), None);
        assert!(all.iter().all(|x| global.holds(&b, x).is_err()), "{path}");
    }
    assert!(found >= 2, "{found}");
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

/// At `r = 0` the plug is the hole, so some member of every maximum form is
/// at least zero there: a box containing an aligned configuration, at a
/// corner, on a face or inside, is refused, however thin.
#[test]
fn boxes_containing_aligned_configurations_are_refused() {
    let global = Global::new();
    let all = global.all();
    let mut random = Random(45);
    for n in 0..24 {
        let depth = random.below(14);
        let view = from_path(&random.path(depth));
        let width = two_to_minus(random.below(60) as u32 + 1);
        let rotation = |kind: usize| match kind {
            0 => interval(q(0), width.clone()),
            1 => interval(-&width, q(0)),
            2 => interval(-&width, width.clone()),
            _ => interval(q(0), q(0)),
        };
        let mut axes: [Interval; AXES] = view.axes().clone();
        for j in 2..AXES {
            axes[j] = rotation((n + j * random.below(4)) % 4);
        }
        let b = ConfigurationBox::new(axes);
        assert_eq!(global.check(&b), None);
        let mut claims = global.candidates(&b);
        claims.extend((0..12).map(|_| all[random.below(all.len())].clone()));
        for m in claims {
            assert!(matches!(global.holds(&b, &m), Err(Refusal::Lemma(_))));
        }
        let aligned = Configuration { s: oracle::to_f64(b.axes()[0].lo()), t: oracle::to_f64(b.axes()[1].lo()), r: [0.0; 3] };
        assert!(oracle::excursion(&aligned).abs() < 1e-9);
    }
}

/// Changing any part of a found witness, or moving it to its parent or a
/// sibling, is refused unless the oracle confirms the changed claim; a claim
/// the oracle refutes is always refused.
#[test]
fn corrupted_and_moved_witnesses_are_refused_unless_true() {
    let global = Global::new();
    let mut tried = 0;
    // Refuted by the oracle, accepted and confirmed, refused though plausible.
    let mut outcomes = [0; 3];
    for root in SMALL_ROOTS {
        let (records, _) = search(&bare(), root, 200).unwrap();
        let Some((path, data)) = records.iter().find(|(_, d)| d.component() == ComponentName::Global) else {
            panic!("no Global record below {root}")
        };
        let b = from_path(path);
        let (from, to, plug) = parts(data);
        let mut changed: Vec<((usize, usize, usize), ConfigurationBox)> = Vec::new();
        for other in 0..VERTEX_COUNT {
            changed.push(((from, to, other), b.clone()));
        }
        for &[x, y] in geometry::edges() {
            changed.push(((x, y, plug), b.clone()));
            changed.push(((y, x, plug), b.clone()));
        }
        let parent = &path[..path.len() - 1];
        let sibling = format!("{parent}{}", if path.ends_with('0') { '1' } else { '0' });
        changed.push(((from, to, plug), from_path(parent)));
        changed.push(((from, to, plug), from_path(&sibling)));
        for (n, (claim, place)) in changed.into_iter().enumerate() {
            let m = maximum(claim.0, claim.1, claim.2);
            let accepted = global.holds(&place, &m).is_ok();
            match oracle_verdict(&place, claim, n as u64) {
                Some(false) => {
                    assert!(!accepted, "{claim:?} is false but accepted");
                    outcomes[0] += 1;
                }
                _ if accepted => {
                    confirm(&place, &RecordData::Global(m), n as u64);
                    outcomes[1] += 1;
                }
                _ => outcomes[2] += 1,
            }
            tried += 1;
        }
        // The reversed edge points to the other side of the hole.
        assert!(global.holds(&b, &maximum(to, from, plug)).is_err());
        // The antipodal candidate states the same maximum form.
        global.holds(&b, &maximum(antipode(from), antipode(to), antipode(plug))).unwrap();
    }
    assert!(tried > 900);
    assert!(outcomes[0] > 0 && outcomes[1] > 0, "{outcomes:?}");
}

/// The maximum form at a configuration, exactly and without the polynomials:
/// `max_w (1 + |r|²) n·w - n·q(0, p)q̄` with `q = (1, r)`, the plug vertex
/// turned by quaternion conjugation, and `n = u × (v_to - v_from)`.
fn maximum_at(m: &MaximumGap, x: &[Q; AXES]) -> QSqrt5 {
    let v = geometry::vertices();
    let c = |k: usize| QSqrt5::from_rational(x[k].clone());
    let u = [c(0), c(1), QSqrt5::one()];
    let r = [c(2), c(3), c(4)];
    let n = geometry::cross(&u, &m.edge().direction());
    let turned = geometry::rotate(&[QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()], &v[m.plug()]);
    let scale = &QSqrt5::one() + &geometry::dot(&r, &r);
    let plug = geometry::dot(&n, &turned);
    v.iter().map(|w| &(&scale * &geometry::dot(&n, w)) - &plug).max().unwrap()
}

#[test]
fn degenerate_boxes_are_decided_like_any_other() {
    // A configuration as a degenerate box: its Bernstein coefficients are its
    // value, so the lemma accepts exactly when every member is negative there.
    let global = Global::new();
    let mut random = Random(47);
    let mut accepted = 0;
    for n in 0..40 {
        let x: [Q; AXES] = std::array::from_fn(|j| {
            let k = random.below(1001) as i64;
            match j {
                0 => frac(2 * k, 3000),
                1 => frac(2 * k, 5000),
                _ => frac(4 * k - 2000, 5000),
            }
        });
        let b = ConfigurationBox::point(&x);
        let candidates = global.candidates(&b);
        // At a point, the lemma accepts exactly the candidates whose maximum
        // form is negative there, and the proposal finds one when one exists.
        let exists = global.all().iter().any(|m| maximum_at(m, &x).sign() == Ordering::Less);
        assert_eq!(candidates.is_empty(), !exists, "{x:?}");
        for m in &candidates {
            global.holds(&b, m).unwrap();
            assert_eq!(maximum_at(m, &x).sign(), Ordering::Less);
        }
        match global.check(&b) {
            Some(m) => {
                accepted += 1;
                assert_eq!(Some(&m), candidates.first());
                confirm(&b, &RecordData::Global(m), n);
            }
            None => assert!(candidates.is_empty()),
        }
    }
    assert!(accepted > 20, "{accepted}");
}

/// The members table gives the verdicts of members built afresh, for the
/// proposed candidate and random ones of the fixed list (asked twice, the
/// second time from the table), and for reversed and parallel edges, which
/// are not in the table. A filled slot holds exactly the fresh members.
#[test]
fn the_members_table_agrees_with_fresh_members() {
    let global = Global::new();
    let all = global.all();
    let mut random = Random(49);
    let fresh = |b: &ConfigurationBox, m: &MaximumGap| lemma::check_maximum(identity(), &cell(b), m, &ZoomFactor::ONE).is_ok();
    let (mut held, mut filled) = (0, 0);
    for b in boxes(10, 50).iter().step_by(3) {
        let mut claims: Vec<MaximumGap> = global.candidates(b);
        claims.extend((0..2).map(|_| all[random.below(all.len())].clone()));
        let edges = geometry::edges();
        let [from, to] = edges[random.below(edges.len())];
        claims.push(MaximumGap::new(Edge::new(to, from).unwrap(), random.below(VERTEX_COUNT)).unwrap());
        for m in &claims {
            let expected = fresh(b, m);
            assert_eq!(global.holds(b, m).is_ok(), expected, "{m:?} on {b:?}");
            assert_eq!(global.holds(b, m).is_ok(), expected, "{m:?} on {b:?}, from the table");
            held += usize::from(expected);
        }
    }
    for (k, slot) in global.members.iter().enumerate() {
        if let Some(members) = slot.get() {
            filled += 1;
            assert_eq!(members, &all[k].members(identity().view(), identity().rotation()).unwrap());
        }
    }
    assert!(held >= 5 && filled >= 10, "{held} {filled}");
}

mod adversarial;
