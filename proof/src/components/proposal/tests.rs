use super::*;
use crate::arithmetic::exact::{QSqrt5, Q};
use crate::components::global::Global;
use crate::components::tests::fixtures::{from_path, search, bare, Random, SMALL_ROOTS};
use crate::components::Component;
use crate::problem::geometry;
use std::cmp::Ordering;

fn random_boxes(count: usize, seed: u64) -> Vec<ConfigurationBox> {
    let mut random = Random(seed);
    (0..count)
        .map(|_| {
            let depth = 8 + random.below(40);
            from_path(&random.path(depth))
        })
        .collect()
}

/// The maximum form of `m` at a configuration, exactly: `max_w (1 + |r|²) n·w
/// − n·q(0, p)q̄` with `q = (1, r)` and `n = u × d`, written independently
/// of the polynomials (quaternion conjugation instead of the Cayley matrix).
fn maximum_at(m: &MaximumGap, x: &[Q; AXES]) -> QSqrt5 {
    let v = geometry::vertices();
    let c = |k: usize| QSqrt5::from_rational(x[k].clone());
    let u = [c(S), c(T), QSqrt5::one()];
    let r = R.map(c);
    let n = geometry::cross(&u, &m.edge().direction());
    let turned = geometry::rotate(&[QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()], &v[m.plug()]);
    let scale = &QSqrt5::one() + &geometry::dot(&r, &r);
    let plug = geometry::dot(&n, &turned);
    v.iter().map(|w| &(&scale * &geometry::dot(&n, w)) - &plug).max().unwrap()
}

/// Domain proposes exactly the inequalities positive at the box's exact
/// centre, wherever the float value is clear of rounding.
#[test]
fn domain_proposes_the_inequalities_positive_at_the_exact_centre() {
    let (mut positive, mut other) = (0, 0);
    for b in random_boxes(60, 81) {
        let proposed: Vec<usize> = domain(&b).iter().map(DomainInequality::index).collect();
        let centre = b.midpoint().map(QSqrt5::from_rational);
        for (i, c) in domain::affine_polynomials().iter().enumerate() {
            let exact = c.evaluate(&centre);
            if float::number(&exact).abs() < 1e-9 {
                continue;
            }
            let is_positive = exact.sign() == Ordering::Greater;
            assert_eq!(proposed.contains(&i), is_positive, "inequality {i} on {b:?}");
            if is_positive {
                positive += 1;
            } else {
                other += 1;
            }
        }
    }
    assert!(positive >= 20 && other >= 1000, "{positive} {other}");
}

/// Global's proposal lies beyond the hole's extent at every corner of the
/// box, checked exactly; and at a corner where some other candidate lies
/// beyond by a clear margin while nothing is proposed, the float ranking
/// would have been wrong, which never happens on these boxes.
#[test]
fn global_proposes_a_candidate_beyond_at_every_corner() {
    let global = Global::new();
    let (records, _) = search(&bare(), SMALL_ROOTS[0], 200).unwrap();
    let mut boxes: Vec<ConfigurationBox> = records.iter().map(|(p, _)| from_path(p)).take(20).collect();
    boxes.extend(random_boxes(20, 82));
    let mut proposed = 0;
    for b in &boxes {
        let Some(m) = global.candidates(b).into_iter().next() else { continue };
        proposed += 1;
        for mask in 0..1usize << AXES {
            let x: [Q; AXES] =
                std::array::from_fn(|j| if mask >> j & 1 == 1 { b.axes()[j].hi().clone() } else { b.axes()[j].lo().clone() });
            assert_eq!(maximum_at(&m, &x).sign(), Ordering::Less, "{m:?} at corner {mask} of {b:?}");
        }
    }
    assert!(proposed >= 15, "{proposed}");
}

/// The proposals are pure functions of the box: the same on every call and
/// on a box rebuilt from the same path.
#[test]
fn proposals_are_deterministic() {
    let global = Global::new();
    let mut random = Random(83);
    for _ in 0..20 {
        let depth = 5 + random.below(40);
        let path = random.path(depth);
        let (a, b) = (from_path(&path), from_path(&path));
        assert_eq!(domain(&a), domain(&b));
        assert_eq!(global.candidates(&a), global.candidates(&b));
    }
}
