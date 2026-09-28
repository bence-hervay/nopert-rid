use super::*;
use crate::testing::Rng;
use rid::arithmetic::exact::{frac, q};

fn exact(x: &Q) -> QSqrt5 {
    QSqrt5::from_rational(x.clone())
}

fn centre(coordinates: [QSqrt5; AXES]) -> Centre {
    Centre::new(coordinates).unwrap()
}

fn rational(values: [Q; AXES]) -> Centre {
    centre(values.map(QSqrt5::from_rational))
}

/// `x` scaled by `2^-bits`.
fn shrink(x: &QSqrt5, bits: u32) -> QSqrt5 {
    x.scale(&radius(bits))
}

/// A random centre in the root box whose coordinates `a + b√5` have `|b|`
/// up to about `2^bits`: `a = target - b·u` for a rational `u` within
/// `2^-(bits+40)` below √5, so the coordinate lies within `2^-40` of the
/// target.
fn irrational_centre(rng: &mut Rng, bits: u32) -> Centre {
    let root = ConfigurationBox::root();
    let m = bits + 40;
    let u = Q::new((BigInt::from(5) << (2 * m)).sqrt(), BigInt::from(1) << m);
    centre(std::array::from_fn(|j| {
        let axis = &root.axes()[j];
        let margin = radius(30);
        let target = rng.rational(&(axis.lo() + &margin), &(axis.hi() - &margin), 30);
        let size = Q::from_integer(BigInt::from(rng.next() >> (64 - bits)));
        let b = if rng.below(2) == 0 { size + frac(1, 3) } else { -size - frac(1, 3) };
        QSqrt5::new(&target - &(&b * &u), b)
    }))
}

#[test]
fn floor_is_exact_near_integers_and_for_large_coefficients() {
    let mut rng = Rng::new(3);
    for _ in 0..2000 {
        let x = rng.number(24);
        let n = floor(&x);
        let below = QSqrt5::from_rational(Q::from_integer(n.clone()));
        let above = QSqrt5::from_rational(Q::from_integer(n + 1));
        assert!(below <= x && x < above, "{x:?}");
    }
    // Pell units a - b√5 = 1/(a + b√5) approach zero from above; their
    // negatives from below; n plus them has floor n or n - 1.
    let (mut a, mut b) = (BigInt::from(9), BigInt::from(4));
    for _ in 0..40 {
        let unit = QSqrt5::new(Q::from_integer(a.clone()), Q::from_integer(-b.clone()));
        for n in [-3i64, 0, 5] {
            let shifted = &unit + &QSqrt5::integer(n);
            assert_eq!(floor(&shifted), BigInt::from(n));
            assert_eq!(floor(&(&QSqrt5::integer(n) - &unit)), BigInt::from(n - 1));
        }
        let (a2, b2) = (&a * 9 + &b * 20, &a * 4 + &b * 9);
        a = a2;
        b = b2;
    }
    for n in -5..5 {
        assert_eq!(floor(&QSqrt5::integer(n)), BigInt::from(n));
    }
}

#[test]
fn rational_centres_get_the_exact_clipped_cube() {
    let root = ConfigurationBox::root();
    let c = [frac(1, 3), frac(1, 5), frac(1, 7), frac(-2, 9), frac(0, 1)];
    for k in [0, 1, 4, 32, 256, 1024, 2048, MAX_EXPONENT] {
        let b = neighbourhood(&rational(c.clone()), k).unwrap();
        for j in 0..AXES {
            let lo = (&c[j] - &radius(k)).max(root.axes()[j].lo().clone());
            let hi = (&c[j] + &radius(k)).min(root.axes()[j].hi().clone());
            assert_eq!(b.axes()[j], Interval::new(lo, hi).unwrap(), "k={k} axis {j}");
        }
    }
}

#[test]
fn every_root_corner_gets_a_full_dimensional_one_sided_box() {
    let root = ConfigurationBox::root();
    for mask in 0..32u8 {
        let corner = rational(root.corner(mask));
        for k in [0, 8, 256, MAX_EXPONENT] {
            let b = neighbourhood(&corner, k).unwrap();
            assert_eq!(check(&corner, k, &b), Ok(()));
            assert!(root.contains(&b));
            for j in 0..AXES {
                assert!(b.axes()[j].lo() < b.axes()[j].hi());
            }
        }
    }
}

#[test]
fn irrational_centres_are_enclosed_within_the_padding_and_nest() {
    let mut rng = Rng::new(4);
    for round in 0..40 {
        let c = irrational_centre(&mut rng, if round % 2 == 0 { 8 } else { 60 });
        let mut previous: Option<ConfigurationBox> = None;
        for k in [0, 3, 12, 64, 512, 2048] {
            let b = neighbourhood(&c, k).unwrap();
            assert_eq!(check(&c, k, &b), Ok(()), "round {round} k={k}");
            let limit = &(q(2) * radius(k)) + &(q(2) * radius(k + 10));
            for a in b.axes() {
                assert!(a.width() <= limit);
            }
            // Contains the exact centre.
            for (x, a) in c.coordinates().iter().zip(b.axes()) {
                assert!(exact(a.lo()) <= *x && *x <= exact(a.hi()));
            }
            if let Some(p) = &previous {
                // Radii shrink by at least half; the padding is 2^-10 of the
                // radius, so a finer box lies in the coarser one.
                assert!(p.contains(&b), "round {round} k={k}");
            }
            previous = Some(b);
        }
    }
}

#[test]
fn exponents_beyond_the_limit_are_refused() {
    let c = rational([q(0), q(0), q(0), q(0), q(0)]);
    assert_eq!(
        neighbourhood(&c, MAX_EXPONENT + 1),
        Err(NeighbourhoodError::Exponent(MAX_EXPONENT + 1))
    );
    let b = neighbourhood(&c, 4).unwrap();
    assert_eq!(
        check(&c, MAX_EXPONENT + 1, &b),
        Err(NeighbourhoodError::Exponent(MAX_EXPONENT + 1))
    );
}

/// Replaces one end of one axis.
fn with_end(b: &ConfigurationBox, axis: usize, side: Side, value: Q) -> Option<ConfigurationBox> {
    let mut axes = b.axes().clone();
    let (lo, hi) = (axes[axis].lo().clone(), axes[axis].hi().clone());
    axes[axis] = match side {
        Side::Lower => Interval::new(value, hi).ok()?,
        Side::Upper => Interval::new(lo, value).ok()?,
    };
    Some(ConfigurationBox::new(axes))
}

#[test]
fn check_refuses_every_kind_of_wrong_box() {
    let mut rng = Rng::new(5);
    let root = ConfigurationBox::root();
    for round in 0..20 {
        let c = irrational_centre(&mut rng, 20);
        let k = [6, 30, 200][round % 3];
        let b = neighbourhood(&c, k).unwrap();
        let r = QSqrt5::from_rational(radius(k));
        let fine = k + 12;
        for axis in 0..AXES {
            let x = &c.coordinates()[axis];
            let bottom = exact(root.axes()[axis].lo());
            let top = exact(root.axes()[axis].hi());
            let round_up = |y: &QSqrt5| -> Q {
                let unit = radius(fine);
                Q::from_integer(-floor(&(-&y.scale(&(Q::one() / &unit))))) * unit
            };
            let round_down = |y: &QSqrt5| -> Q {
                let unit = radius(fine);
                Q::from_integer(floor(&y.scale(&(Q::one() / &unit)))) * unit
            };
            // An end strictly inside the exact neighbourhood omits part of it.
            let ideal_lo = (x - &r).max(bottom.clone());
            if ideal_lo > bottom {
                let inside = round_up(&(&ideal_lo + &shrink(&r, fine - k + 2)));
                if let Some(bad) = with_end(&b, axis, Side::Lower, inside) {
                    assert_eq!(check(&c, k, &bad), Err(NeighbourhoodError::Misses { axis, side: Side::Lower }));
                }
            }
            let ideal_hi = (x + &r).min(top.clone());
            if ideal_hi < top {
                let inside = round_down(&(&ideal_hi - &shrink(&r, fine - k + 2)));
                if let Some(bad) = with_end(&b, axis, Side::Upper, inside) {
                    assert_eq!(check(&c, k, &bad), Err(NeighbourhoodError::Misses { axis, side: Side::Upper }));
                }
            }
            // An end more than radius/1024 beyond, but inside the root box.
            let far_lo = round_down(&(&(x - &r) - &shrink(&r, 9)));
            if exact(&far_lo) >= bottom {
                let bad = with_end(&b, axis, Side::Lower, far_lo).unwrap();
                assert_eq!(check(&c, k, &bad), Err(NeighbourhoodError::Padding { axis, side: Side::Lower }));
            }
            let far_hi = round_up(&(&(x + &r) + &shrink(&r, 9)));
            if exact(&far_hi) <= top {
                let bad = with_end(&b, axis, Side::Upper, far_hi).unwrap();
                assert_eq!(check(&c, k, &bad), Err(NeighbourhoodError::Padding { axis, side: Side::Upper }));
            }
            // Leaving the root box, or a flat axis.
            let outside = with_end(&b, axis, Side::Lower, root.axes()[axis].lo() - &radius(300)).unwrap();
            assert_eq!(check(&c, k, &outside), Err(NeighbourhoodError::Shape { axis }));
            let lo = b.axes()[axis].lo().clone();
            let flat = with_end(&b, axis, Side::Upper, lo).unwrap();
            assert_eq!(check(&c, k, &flat), Err(NeighbourhoodError::Shape { axis }));
        }
        // The box of another exponent is refused.
        let coarser = neighbourhood(&c, k - 1).unwrap();
        assert!(check(&c, k, &coarser).is_err());
        let finer = neighbourhood(&c, k + 1).unwrap();
        assert!(check(&c, k, &finer).is_err());
    }
}
