use super::*;
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::components::tests::fixtures::{from_path, two_to_minus, Random};
use crate::components::tests::oracle::{self, MARGIN};
use crate::problem::configuration::AXES;
use crate::problem::domain::CONSTRAINT_COUNT;
use crate::problem::geometry;
use std::cmp::Ordering;

fn record(index: usize) -> DomainInequality {
    DomainInequality::new(index).unwrap()
}

fn symmetries() -> Vec<[f64; 4]> {
    geometry::symmetries()
        .iter()
        .map(|g| std::array::from_fn(|j| oracle::number(&g[j])))
        .collect()
}

/// Checks every inequality on `b` against the component's check and against
/// the independent oracle: the check names the lowest inequality that holds;
/// an inequality that the oracle shows satisfied somewhere in the box does
/// not hold; one that holds is violated at every sample, and every sample
/// lies outside `D` by the oracle's own group.
fn check_box(b: &ConfigurationBox, seed: u64) -> Option<usize> {
    let accepted: Vec<usize> = (0..CONSTRAINT_COUNT)
        .filter(|&i| match Domain.holds(b, &record(i)) {
            Ok(()) => true,
            Err(Refusal::Lemma(_)) => false,
            Err(other) => panic!("unexpected refusal {other}"),
        })
        .collect();
    assert_eq!(Domain.check(b), accepted.first().map(|&i| record(i)));
    let symmetries = symmetries();
    let samples = oracle::samples(b, seed, 8);
    for sample in &samples {
        let values = oracle::inequalities(sample, &symmetries);
        for &i in &accepted {
            assert!(values[i] > -MARGIN, "inequality {i} accepted but satisfied at {sample:?}");
        }
        if !accepted.is_empty() {
            assert!(oracle::outside_domain(sample) > -MARGIN);
        }
    }
    for i in 0..CONSTRAINT_COUNT {
        let refuted = samples
            .iter()
            .any(|x| oracle::inequalities(x, &symmetries)[i] < -MARGIN);
        if refuted {
            assert!(!accepted.contains(&i), "inequality {i} is refuted by the oracle");
        }
    }
    accepted.first().copied()
}

#[test]
fn check_and_holds_agree_with_each_other_and_the_oracle_on_random_boxes() {
    let mut random = Random(21);
    let mut found = [0; 3];
    for n in 0..240 {
        let depth = random.below(46);
        let b = from_path(&random.path(depth));
        match check_box(&b, n) {
            None => found[0] += 1,
            Some(i) if i < 13 => found[1] += 1,
            Some(_) => found[2] += 1,
        }
    }
    // Refusals, affine and fold inequalities all occur.
    assert!(found.iter().all(|&n| n >= 10), "{found:?}");
}

/// The proposal is in index order and keeps every inequality that holds:
/// holding needs `c_i > 0` throughout the box, so at its centre.
#[test]
fn the_proposal_keeps_every_inequality_that_holds() {
    let mut random = Random(24);
    let mut held = 0;
    for depth in [0, 1, 7, 12, 18, 30, 61].into_iter().cycle().take(70) {
        let b = from_path(&random.path(depth));
        let indices: Vec<usize> = Domain.candidates(&b).iter().map(DomainInequality::index).collect();
        assert!(indices.windows(2).all(|w| w[0] < w[1]));
        for i in 0..CONSTRAINT_COUNT {
            if Domain.holds(&b, &DomainInequality::new(i).unwrap()).is_ok() {
                held += 1;
                assert!(indices.contains(&i), "inequality {i} holds on {b:?} but is not proposed");
            }
        }
    }
    assert!(held >= 5, "{held}");
}

fn phi() -> QSqrt5 {
    QSqrt5::new(frac(1, 2), frac(1, 2))
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

/// Boxes whose view rectangle ends just below, at or across the irrational
/// triangle wall `φs + φ²t = 1`, at widths down to `2⁻²⁰⁰`: inequality 0
/// holds exactly when `φs + φ²t > 1` at the lowest corner, and the check
/// finds it wherever the box is wider than floating-point resolution.
#[test]
fn boxes_at_the_triangle_wall_are_decided_exactly() {
    let root = ConfigurationBox::root();
    let phi = phi();
    let phi2 = &phi * &phi;
    let mut random = Random(22);
    for _ in 0..60 {
        let s = frac(random.below(600) as i64, 1000);
        // t* = (1 - φs)/φ², approximated by rationals on both sides.
        let exact_t = |t: &Q| &(&phi.scale(&s) + &phi2.scale(t)) - &QSqrt5::one();
        let unit = two_to_minus(1 + random.below(200) as u32);
        let mut t = frac(0, 1);
        let mut step = frac(1, 4);
        // Bisection for the largest dyadic t with φs + φ²t < 1.
        for _ in 0..220 {
            let trial = &t + &step;
            if exact_t(&trial).sign() == Ordering::Less {
                t = trial;
            }
            step = &step / &q(2);
        }
        for (lo, hi) in [
            (&t - &unit, t.clone()),
            (t.clone(), &t + &unit),
            (&t + &step * &q(4), &t + &unit),
            (&t + &unit, &t + &unit + &unit),
        ] {
            let mut axes = root.axes().clone();
            axes[0] = interval(s.clone(), s.clone());
            axes[1] = interval(lo.clone(), hi);
            let b = ConfigurationBox::new(axes);
            let violated = exact_t(&lo).sign() == Ordering::Greater;
            let verdict = Domain.holds(&b, &record(0)).is_ok();
            assert_eq!(verdict, violated);
            // The floating-point proposal sees the violation at the centre
            // unless the box is thinner than its resolution; then the box
            // is split instead, never wrongly eliminated.
            let check = Domain.check(&b);
            assert!(check.is_none() || check == Some(record(0)));
            if violated && unit >= two_to_minus(40) {
                assert_eq!(check, Some(record(0)));
            }
        }
    }
}

/// Boxes at the irrational walls of the rotation inequalities
/// `ε r_i + δ(φ - 1) r_{i+1} + φ - 2 ≤ 0`, with the other coordinates at the
/// origin: `r_i` ends just before, at or beyond `2 - φ`.
#[test]
fn boxes_at_the_rotation_walls_are_decided_exactly() {
    let root = ConfigurationBox::root();
    let wall = &QSqrt5::integer(2) - &phi();
    for bits in [8u32, 40, 120, 250] {
        let unit = two_to_minus(bits);
        // The dyadic rational just below 2 - φ at this resolution.
        let mut below = q(0);
        let mut step = frac(1, 2);
        for _ in 0..bits {
            let trial = &below + &step;
            if (&QSqrt5::from_rational(trial.clone()) - &wall).sign() == Ordering::Less {
                below = trial;
            }
            step = &step / &q(2);
        }
        for axis in 0..3 {
            for (first, lo, hi) in [
                (1, below.clone(), &below + &unit),
                (1, &below + &unit, &below + &unit + &unit),
                (-1, -(&below + &unit), -&below),
                (-1, -(&below + &unit + &unit), -(&below + &unit)),
            ] {
                let mut axes: [Interval; AXES] = root.axes().clone();
                for j in 0..3 {
                    axes[2 + j] = interval(q(0), q(0));
                }
                axes[2 + axis] = interval(lo.clone(), hi.clone());
                let b = ConfigurationBox::new(axes);
                // Inequalities 1 + 4·axis + (first + 1) + (second + 1)/2 with
                // r_{i+1} = 0 do not depend on `second`.
                let beyond = if first == 1 { lo.clone() } else { -&hi };
                let violated = (&QSqrt5::from_rational(beyond) - &wall).sign() == Ordering::Greater;
                for second in [0, 1] {
                    let index = 1 + 4 * axis + 2 * usize::from(first == 1) + second;
                    assert_eq!(Domain.holds(&b, &record(index)).is_ok(), violated, "{index}");
                }
                assert_eq!(Domain.check(&b).is_some(), violated);
            }
        }
    }
}

#[test]
fn degenerate_and_deep_boxes_agree_with_the_oracle() {
    let mut random = Random(23);
    for n in 0..40 {
        let depth = 40 + random.below(40);
        let b = from_path(&random.path(depth));
        // The box's corners as degenerate boxes.
        let corner = ConfigurationBox::point(&b.corner(random.below(32) as u8));
        check_box(&corner, n);
        check_box(&b, n);
    }
}

/// Against the exact box decision of the domain module's tests: an
/// inequality that holds is violated throughout the box; for the affine
/// inequalities 0 to 12 (whose Bernstein coefficients are the corner values)
/// the two decisions agree exactly.
#[test]
fn lemma_decisions_agree_with_the_exact_box_decision() {
    use crate::problem::domain::tests::decision::positive_on;
    let mut random = Random(23);
    let (mut certified, mut folds) = (0, 0);
    for _ in 0..160 {
        let depth = random.below(40);
        let b = from_path(&random.path(depth));
        for i in 0..CONSTRAINT_COUNT {
            let lemma = Domain.holds(&b, &record(i)).is_ok();
            let exact = positive_on(&b, i).unwrap();
            assert!(!lemma || exact, "inequality {i} certified on {b:?} but not violated throughout");
            if i < 13 {
                assert_eq!(lemma, exact, "affine inequality {i} on {b:?}");
            }
            if lemma {
                certified += 1;
                folds += usize::from(i >= 13);
            }
        }
    }
    assert!(certified > 50 && folds > 10, "{certified} {folds}");
}
