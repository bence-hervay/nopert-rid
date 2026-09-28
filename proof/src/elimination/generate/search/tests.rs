//! Tests of the cell decision: the split rule in isolation, accepted leaves
//! re-checked by the lemma, delegation exactly where the window holds,
//! refusal at the depth limit, and decisions that depend only on the cell.
use super::*;
use std::num::{NonZeroU32, NonZeroUsize};
use crate::arithmetic::exact::{frac, q, Interval, Q};
use crate::arithmetic::polynomial::Polynomial;
use crate::elimination::zoom::cover::ZoomCover;
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::tree;
use crate::elimination::proof::catalogue::{exotic_data, exotic_pin};
use crate::elimination::proof::format::{Leaf, CoverFile};

fn settings() -> Settings {
    Settings { max_depth: 60, max_nodes: 100_000, aspect: sixteen(), tried: 6, proposal_depth: 4, proposals: 8, inherited: 16 }
}

fn sixteen() -> NonZeroU32 {
    NonZeroU32::new(16).unwrap()
}

/// Bernstein coefficients on the unit box that vary only along `axis`, by `step`.
fn varying(axis: usize, step: i64) -> Dense {
    let mut unit = [0u8; VARIABLES];
    unit[axis] = 1;
    let p = Polynomial::monomial(unit, crate::arithmetic::exact::QSqrt5::integer(step));
    Dense::new(&p).bernstein(&[[0.0, 1.0]; VARIABLES])
}

#[test]
fn the_split_follows_the_largest_variation_under_the_aspect_cap() {
    let full = [0; VARIABLES];
    // No candidate tried: the relatively widest variable, the first on ties.
    assert_eq!(split_axis(&full, &[], sixteen()), 0);
    let mut narrow = full;
    narrow[0] = 1;
    assert_eq!(split_axis(&narrow, &[], sixteen()), 1);
    // The variation decides among allowed variables.
    let tried = vec![(-1.0, varying(3, 2))];
    assert_eq!(split_axis(&full, &tried, sixteen()), 3);
    // A variable 32 times narrower than the widest one is not split, even
    // where the best candidate varies most.
    let mut thin = full;
    thin[3] = 5;
    let two = vec![(-1.0, varying(3, 2)), (-0.5, varying(1, 1))];
    assert_eq!(split_axis(&thin, &two, sixteen()), 1);
    // ... and with only one candidate varying along it, the widest variable is
    // split (no allowed variable varies).
    assert_eq!(split_axis(&thin, &tried, sixteen()), 0);
    // Exactly 16 times narrower is allowed, at every depth.
    for extra in [0, 7, 40, 500] {
        let mut edge = [extra; VARIABLES];
        edge[3] = extra + 4;
        assert_eq!(split_axis(&edge, &tried, sixteen()), 3);
    }
    // Every variable 32 times narrower than one wide variable: only the wide
    // one is allowed.
    let mut all_thin = [5; VARIABLES];
    all_thin[2] = 0;
    assert_eq!(split_axis(&all_thin, &tried, sixteen()), 2);
    // Only the two best tried candidates count.
    let three = vec![
        (-0.1, varying(4, 100)),
        (-1.0, varying(2, 1)),
        (-0.9, varying(2, 1)),
    ];
    assert_eq!(split_axis(&full, &three, sixteen()), 2);
    // An aspect that is not a power of two: 2^(k − min k) ≤ 20 allows 16, not 32.
    let twenty = NonZeroU32::new(20).unwrap();
    let mut cap = full;
    cap[3] = 4;
    assert_eq!(split_axis(&cap, &tried, twenty), 3);
    cap[3] = 5;
    assert_eq!(split_axis(&cap, &tried, twenty), 0);
    // An aspect of 1 splits only the relatively widest variables; split
    // counts far apart (beyond 64) are handled.
    let mut far = [70; VARIABLES];
    far[1] = 0;
    assert_eq!(split_axis(&far, &tried, NonZeroU32::new(1).unwrap()), 1);
    assert_eq!(split_axis(&far, &tried, NonZeroU32::new(u32::MAX).unwrap()), 1);
}

/// Split counts are read exactly from the cell against the root, also for
/// roots of non-dyadic widths and deep cells.
#[test]
fn split_counts_are_exact() {
    let root: Five<Interval> = std::array::from_fn(|j| Interval::new(q(0), frac(5, 3 + j as i64)).unwrap());
    let mut cell = root.clone();
    let mut expected = [0u32; VARIABLES];
    let mut random = Random::new(5);
    for _ in 0..300 {
        let j = random.below(VARIABLES as u64) as usize;
        let [lower, upper] = cell[j].bisect();
        cell[j] = if random.below(2) == 0 { lower } else { upper };
        expected[j] += 1;
        assert_eq!(splits(&root, &cell), expected);
    }
}

/// The leaves of zoom `zoom` of an Exotic cover of the catalogue's
/// data, with their cells.
fn leaves(name: &str, zoom: usize) -> (ZoomCover, Vec<(ZoomCell, Leaf)>) {
    let file = CoverFile::parse(exotic_data(name).unwrap()).unwrap();
    let a = ZoomCover::new(file.parameters.clone()).unwrap();
    let cells = tree::leaves(a.zooms()[zoom].root(), &file.zooms[zoom].tree).unwrap();
    let list = cells.into_iter().zip(file.zooms[zoom].leaves.clone()).collect();
    (a, list)
}

#[test]
fn accepted_leaves_pass_the_lemma_and_decisions_depend_only_on_the_cell() {
    let (cover, cells) = leaves("arc-", 0);
    let (u, r) = candidates::samples(&cover, 0, NonZeroUsize::new(4).unwrap()).into_iter().next().unwrap();
    let gaps = candidates::active(&u, &r);
    let search = ZoomSearch::new(&cover, 0, &gaps, &[]).unwrap();
    let mut leaves = 0;
    for (cell, _) in cells.iter().take(12) {
        let first = search.decide(cell, 8, &[], &settings()).unwrap();
        match &first {
            Decision::Leaf(c) => {
                leaves += 1;
                lemma::check_cell(search.zoom(), cell, c.witness(), c.factor()).expect("a leaf passes the lemma");
            }
            Decision::Split { axis, .. } => assert!(*axis < VARIABLES),
            Decision::Delegated | Decision::Unresolved => panic!("no window and not at the depth limit"),
        }
        let again = search.decide(cell, 8, &[], &settings()).unwrap();
        assert_eq!(format!("{first:?}"), format!("{again:?}"));
    }
    assert!(leaves > 0);
    // The root is too large for any single witness.
    let root = search.zoom().root().clone();
    assert!(matches!(search.decide(&root, 0, &[], &settings()).unwrap(), Decision::Split { .. }));
    // At the depth limit an undecided cell is unresolved.
    let limited = Settings { max_depth: 0, ..settings() };
    assert!(matches!(search.decide(&root, 0, &[], &limited).unwrap(), Decision::Unresolved));
}

/// A cell of an arc-plane point-cover zoom from its intervals.
fn cell(ends: [(Q, Q); VARIABLES]) -> ZoomCell {
    ZoomCell(ends.map(|(lo, hi)| Interval::new(lo, hi).unwrap()))
}

#[test]
fn cells_are_delegated_exactly_where_the_window_holds() {
    let pin = exotic_pin("crossing+").unwrap();
    let cover = ZoomCover::new(pin.parameters.clone()).unwrap();
    let no_proposals = Settings { tried: 0, proposal_depth: 100, ..settings() };
    let position = |name: &str| cover.zooms().iter().position(|c| c.name() == name).unwrap();
    let (window, other, sheared, b) = (position("A b0+ o1-"), position("A b0- o1-"), position("A' b0+ o1-"), position("B o1-"));
    // Variables of these zooms: μ, ρ, then the offset coordinates v, ζ₁, ζ₃
    // (the θ face is −1). The window vector is (ρv, 1 − ρ, ρζ₁, ρζ₃).
    let near = |rho: (Q, Q), free: Q| {
        cell([(q(0), frac(1, 64)), rho, (q(0), free.clone()), (-free.clone(), free.clone()), (-free.clone(), free)])
    };
    let inside = near((frac(7, 8), frac(9, 8)), frac(1, 8));
    let boundary = near((frac(3, 4), frac(5, 4)), frac(1, 5));
    let outside = near((frac(3, 4), frac(5, 4)), frac(1, 4));
    let far = near((q(0), frac(1, 4)), frac(1, 8));
    let mut random = Random::new(21);
    let mut cases: Vec<(usize, ZoomCell)> = Vec::new();
    for c in [window, other, sheared] {
        for x in [&inside, &boundary, &outside, &far] {
            cases.push((c, x.clone()));
        }
    }
    for c in [window, sheared, b] {
        let root = cover.zooms()[c].root().clone();
        for _ in 0..20 {
            let mut x = root.clone();
            for _ in 0..random.below(14) {
                let axis = random.below(VARIABLES as u64) as usize;
                let [lo, hi] = x[axis].bisect();
                x[axis] = if random.below(2) == 0 { lo } else { hi };
            }
            cases.push((c, x));
        }
    }
    let mut seen = [0; 2];
    for (c, x) in cases {
        let search = ZoomSearch::new(&cover, c, &[], &[]).unwrap();
        let delegated = matches!(search.decide(&x, 60, &[], &no_proposals).unwrap(), Decision::Delegated);
        let holds = c == window && cover.window_holds(c, &x);
        assert_eq!(delegated, holds, "zoom {} cell {x:?}", cover.zooms()[c].name());
        seen[usize::from(delegated)] += 1;
    }
    assert!(seen[1] >= 2 && seen[0] > 0, "{seen:?}");
    let search = ZoomSearch::new(&cover, window, &[], &[]).unwrap();
    for (x, expected) in [(&inside, true), (&boundary, true), (&outside, false), (&far, false)] {
        assert_eq!(matches!(search.decide(x, 60, &[], &no_proposals).unwrap(), Decision::Delegated), expected);
    }
}

#[test]
fn proposals_are_handed_down_and_capped() {
    let pin = exotic_pin("square").unwrap();
    let cover = ZoomCover::new(pin.parameters.clone()).unwrap();
    let search = ZoomSearch::new(&cover, 0, &[], &[]).unwrap();
    let root = search.zoom().root().clone();
    let s = Settings { proposal_depth: 0, proposals: 3, inherited: 4, ..settings() };
    let Decision::Split { proposed, .. } = search.decide(&root, 0, &[], &s).unwrap() else {
        panic!("the root is split");
    };
    assert!(!proposed.is_empty() && proposed.len() <= 4);
    let Decision::Split { proposed: again, .. } = search.decide(&root, 0, &proposed, &s).unwrap() else {
        panic!("the root is split");
    };
    assert!(again.len() <= 4);
    for p in &again {
        assert_eq!(again.iter().filter(|x| x.witness() == p.witness()).count(), 1, "no duplicates");
    }
    // Below the proposal depth nothing is proposed, inherited ones pass on.
    let late = Settings { proposal_depth: 5, ..s };
    let Decision::Split { proposed: passed, .. } = search.decide(&root, 0, &proposed, &late).unwrap() else {
        panic!("the root is split");
    };
    assert_eq!(passed.len(), proposed.len());
    let _ = q(0);
}
