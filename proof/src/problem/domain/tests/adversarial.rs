//! Independent adversarial tests of the domain decision: exactness under
//! bisection, point boxes, and exact points on the closed boundary of `D`.
use crate::arithmetic::exact::tests::sqrt5_enclosure;
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::problem::configuration::{ConfigurationBox, AXES};
use crate::problem::domain::{self, Constraint, Sign, CONSTRAINT_COUNT, FOLD_OFFSET};
use crate::problem::geometry;
use std::cmp::Ordering;

type C = QSqrt5;

fn int(n: i64) -> C {
    C::integer(n)
}

fn rat(x: Q) -> C {
    C::from_rational(x)
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
}

// ---- Exactness through bisection and at point boxes ----------------------

fn random_box(random: &mut Random) -> ConfigurationBox {
    let root = ConfigurationBox::root();
    ConfigurationBox::new(std::array::from_fn(|j| {
        let axis = &root.axes()[j];
        let slack = axis.width() / q(4);
        let lo = random.between(&(axis.lo() - &slack), &(axis.hi() + &slack));
        let width = match random.below(5) {
            0 => q(0),
            1 => axis.width() / q(1 << random.below(40)),
            _ => {
                let shrink = q(1 << random.below(8));
                random.between(&q(0), &(axis.width() / shrink))
            }
        };
        Interval::new(lo.clone(), lo + width).unwrap()
    }))
}

/// F > 0 on a box iff F > 0 on both closed halves: any inexact decision
/// (in either direction) breaks this for some split. Returns the positive
/// decisions (affine, fold) and the number of splits with different halves.
fn check_bisection(random: &mut Random, boxes: usize) -> ([usize; 2], usize) {
    let mut positives = [0usize; 2];
    let mut mixed = 0;
    for _ in 0..boxes {
        let b = random_box(random);
        for axis in 0..AXES {
            let [lower, upper] = b.split(axis);
            for index in 0..CONSTRAINT_COUNT {
                let parent = super::decision::positive_on(&b, index).unwrap();
                let low = super::decision::positive_on(&lower, index).unwrap();
                let high = super::decision::positive_on(&upper, index).unwrap();
                assert_eq!(parent, low && high, "index {index} axis {axis} on {b:?}");
                positives[usize::from(index >= FOLD_OFFSET)] += usize::from(parent);
                mixed += usize::from(low != high);
            }
        }
    }
    (positives, mixed)
}

#[test]
fn box_decisions_are_exact_under_bisection_of_every_axis() {
    let (positives, mixed) = check_bisection(&mut Random(271828), 40);
    assert!(positives[0] > 0 && positives[1] > 0 && mixed > 0, "{positives:?} {mixed}");
}

/// Heavy campaign: `cargo test --release --lib -- --ignored
/// problem::domain::tests::adversarial::bisection_campaign --nocapture`.
#[test]
#[ignore]
fn bisection_campaign() {
    let (positives, mixed) = check_bisection(&mut Random(314159), 1000);
    println!("positive decisions (affine, fold) {positives:?}; splits with mixed halves {mixed}");
}

#[test]
fn point_box_decisions_equal_the_sign_of_the_polynomial() {
    let mut random = Random(161803);
    let polynomials = domain::affine_polynomials();
    for _ in 0..40 {
        let b = random_box(&mut random);
        let x = b.corner(random.below(32) as u8);
        let point = ConfigurationBox::point(&x);
        let exact = x.clone().map(rat);
        for index in 0..CONSTRAINT_COUNT {
            let positive = polynomials[index].evaluate(&exact).sign() == Ordering::Greater;
            assert_eq!(super::decision::positive_on(&point, index), Ok(positive), "index {index}");
        }
        assert_eq!(
            super::decision::misses(&point),
            (0..CONSTRAINT_COUNT)
                .find(|&i| polynomials[i].evaluate(&exact).sign() == Ordering::Greater)
        );
    }
}

// ---- Exact points on the closed boundary of D ----------------------------

/// Views `(s, t)` with rational `n = |(s, t, 1)|`, from integer solutions of
/// `a² + b² + c² = d²` (then `(s, t, n) = (a/c, b/c, d/c)`).
fn pythagorean_views() -> Vec<(Q, Q, Q)> {
    let views = vec![
        (q(0), q(0), q(1)),
        (q(0), frac(3, 4), frac(5, 4)),
        (frac(3, 4), q(0), frac(5, 4)),
        (frac(1, 3), frac(1, 2), frac(7, 6)),
        (frac(1, 2), frac(1, 3), frac(7, 6)),
        (frac(4, 7), frac(4, 7), frac(9, 7)),
        (frac(1, 8), frac(1, 2), frac(9, 8)),
    ];
    for (s, t, n) in &views {
        assert_eq!(q(1) + s * s + t * t, n * n);
    }
    views
}

#[test]
fn exact_fold_boundary_points_are_not_positive_and_just_outside_is() {
    // For group elements with rational coordinates L is rational; solve
    // L = ±|u| for r₃ and probe the closed boundary of D exactly.
    let mut random = Random(99991);
    let mut probes = 0;
    for (symmetry, g) in geometry::symmetries().iter().enumerate() {
        if g.iter().any(|c| *c.sqrt5_part() != q(0)) {
            continue;
        }
        let index = domain::fold_index(symmetry).unwrap();
        let g: [Q; 4] = std::array::from_fn(|j| g[j].rational_part().clone());
        for (s, t, norm) in pythagorean_views() {
            for _ in 0..3 {
                let (r1, r2) = (random.rational(1) / q(3), random.rational(1) / q(3));
                // L = a + b r₃ with u = (s, t, 1).
                let cross = [
                    &g[2] - &g[3] * &t,
                    &g[3] * &s - &g[1],
                    &g[1] * &t - &g[2] * &s,
                ];
                let a = -(&g[1] * &s + &g[2] * &t + &g[3])
                    - &r1 * (&g[0] * &s + &cross[0])
                    - &r2 * (&g[0] * &t + &cross[1]);
                let b = -(&g[0] + &cross[2]);
                if b == q(0) {
                    continue;
                }
                for sign in [1, -1] {
                    let target = &norm * q(sign);
                    let r3 = (&target - &a) / &b;
                    // Moving r₃ by δ changes L by bδ; outward means |L| grows.
                    let outward = if (&b * q(sign)) > q(0) { q(1) } else { q(-1) };
                    let x = [s.clone(), t.clone(), r1.clone(), r2.clone(), r3.clone()];
                    let value = domain::affine_polynomials()[index].evaluate(&x.clone().map(rat));
                    assert_eq!(value, int(0));
                    assert_eq!(super::decision::positive_on(&ConfigurationBox::point(&x), index), Ok(false));
                    let epsilon = frac(1, 1 << 40);
                    let mut axes: [Interval; AXES] =
                        std::array::from_fn(|j| Interval::point(x[j].clone()));
                    // Closed box touching the boundary from outside: refused.
                    let far = &r3 + &outward * &epsilon;
                    axes[4] =
                        Interval::new(r3.clone().min(far.clone()), r3.clone().max(far.clone()))
                            .unwrap();
                    let touching = ConfigurationBox::new(axes.clone());
                    assert_eq!(super::decision::positive_on(&touching, index), Ok(false));
                    // Strictly outside: positive.
                    let near = &r3 + &outward * &epsilon / q(2);
                    axes[4] = Interval::new(near.clone().min(far.clone()), near.max(far)).unwrap();
                    let outside = ConfigurationBox::new(axes.clone());
                    assert_eq!(super::decision::positive_on(&outside, index), Ok(true));
                    // Straddling the boundary: refused.
                    axes[4] = Interval::new(&r3 - &epsilon, &r3 + &epsilon).unwrap();
                    let straddling = ConfigurationBox::new(axes);
                    assert_eq!(super::decision::positive_on(&straddling, index), Ok(false));
                    probes += 1;
                }
            }
        }
    }
    assert!(probes > 50, "{probes}");
}

#[test]
fn affine_inequalities_are_exact_at_their_irrational_walls() {
    // Inequality 1 is -r₁ - (φ-1) r₂ + φ - 2 ≤ 0; with r₂ = 0 it is violated
    // iff r₁ < φ - 2 = (√5 - 3)/2. Rational r₁ within 2^-300 on either side.
    assert_eq!(
        domain::constraint(1),
        Ok(Constraint::Rotation { axis: 0, first: Sign::Minus, second: Sign::Minus })
    );
    let enclosure = sqrt5_enclosure(320);
    for (root, violated) in [(enclosure.lo(), true), (enclosure.hi(), false)] {
        let r1 = (root - q(3)) / q(2);
        let x = [q(0), q(0), r1.clone(), q(0), q(0)];
        assert_eq!(super::decision::positive_on(&ConfigurationBox::point(&x), 1), Ok(violated));
        // A box reaching across the wall is never positive.
        let across = ConfigurationBox::new([
            Interval::point(q(0)),
            Interval::point(q(0)),
            Interval::new(enclosure.lo() / q(2) - frac(3, 2), enclosure.hi() / q(2) - frac(3, 2))
                .unwrap(),
            Interval::point(q(0)),
            Interval::point(q(0)),
        ]);
        assert_eq!(super::decision::positive_on(&across, 1), Ok(false));
    }
}

/// The builder's suite survived accepting a zero minimum in the affine
/// decision (`!= Less` for `== Greater`); this test refuses that mutant.
#[test]
fn affine_decisions_are_strict_at_the_rational_points_of_the_irrational_walls() {
    // An affine wall with an irrational slope contains exactly one rational
    // point of its plane: (s, t) = (-1, 1) for the triangle (φ² - φ - 1 = 0),
    // and (rᵢ, rᵢ₊₁) = (first, -second) for a rotation inequality. There the
    // exact minimum over a point box is 0, so the box meets the closed wall and
    // must not be reported as missing D.
    let triangle = ConfigurationBox::point(&[q(-1), q(1), q(0), q(0), q(0)]);
    assert_eq!(super::decision::positive_on(&triangle, 0), Ok(false));
    let polynomials = domain::affine_polynomials();
    assert_eq!(polynomials[0].evaluate(&[int(-1), int(1), int(0), int(0), int(0)]), int(0));
    for index in 1..FOLD_OFFSET {
        let Ok(Constraint::Rotation { axis, first, second }) = domain::constraint(index) else {
            panic!("rotation inequality expected at {index}");
        };
        let mut x = [q(0), q(0), q(0), q(0), q(0)];
        x[2 + axis] = q(first.value());
        x[2 + (axis + 1) % 3] = q(-second.value());
        x[2 + (axis + 2) % 3] = frac(1, 3);
        assert_eq!(polynomials[index].evaluate(&x.clone().map(rat)), int(0), "index {index}");
        let point = ConfigurationBox::point(&x);
        assert_eq!(super::decision::positive_on(&point, index), Ok(false), "index {index}");
        // A box with that point as its minimising corner is not positive either.
        let mut axes: [Interval; AXES] = std::array::from_fn(|j| Interval::point(x[j].clone()));
        let lo = x[2 + axis].clone();
        axes[2 + axis] = if first == Sign::Plus {
            Interval::new(lo.clone(), lo + q(1)).unwrap()
        } else {
            Interval::new(lo.clone() - q(1), lo).unwrap()
        };
        let cornered = ConfigurationBox::new(axes);
        assert_eq!(super::decision::positive_on(&cornered, index), Ok(false), "index {index}");
    }
}
