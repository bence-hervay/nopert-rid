use super::*;
use num_traits::One;

pub(crate) mod independent;

// ------------------------------------------------------- test-only enclosures

/// Largest precision accepted by [`sqrt5_enclosure`].
pub(crate) const MAX_SQRT5_BITS: u32 = 4096;

/// A rational interval of width `2^-bits` containing √5, for independent
/// cross-checks in tests (no library decision uses an approximation of √5).
/// With `d = 2^bits` and `k = isqrt(5d²)`, the integer check
/// `k² <= 5d² < (k+1)²` proves `k/d <= √5 < (k+1)/d`.
pub(crate) fn sqrt5_enclosure(bits: u32) -> Interval {
    assert!(bits <= MAX_SQRT5_BITS, "sqrt5 precision {bits} exceeds {MAX_SQRT5_BITS} bits");
    let denominator = BigInt::one() << bits;
    let square = BigInt::from(5) * &denominator * &denominator;
    let root = square.sqrt();
    let next = &root + BigInt::one();
    assert!(&root * &root <= square && &next * &next > square, "integer square root");
    Interval::new(Q::new(root, denominator.clone()), Q::new(next, denominator)).unwrap()
}

/// The natural power `x^exponent` (`x^0 = 1`), by repeated multiplication.
pub(crate) fn rational_power(x: &Q, exponent: u32) -> Q {
    (0..exponent).fold(q(1), |product, _| product * x)
}

/// The exact range of `x^exponent` for `x` in the interval (`x^0 = 1`): odd
/// powers are monotone; even powers have their maximum at an end and their
/// minimum at 0 when the interval contains 0, otherwise at an end.
pub(crate) fn interval_power(range: &Interval, exponent: u32) -> Interval {
    if exponent == 0 {
        return Interval::point(q(1));
    }
    let (lo, hi) = (rational_power(range.lo(), exponent), rational_power(range.hi(), exponent));
    if exponent % 2 == 1 {
        return Interval::new(lo, hi).unwrap();
    }
    let low = if range.contains(&q(0)) { q(0) } else { lo.clone().min(hi.clone()) };
    Interval::new(low, lo.max(hi)).unwrap()
}

/// splitmix64: deterministic, well mixed.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
    fn rational(&mut self) -> Q {
        frac(self.int(-1000, 1000), self.int(1, 97))
    }
    fn number(&mut self) -> QSqrt5 {
        QSqrt5::new(self.rational(), self.rational())
    }
    fn interval(&mut self) -> Interval {
        let (a, b) = (self.rational(), self.rational());
        Interval::new(a.clone().min(b.clone()), a.max(b)).unwrap()
    }
}

fn big(n: &BigInt) -> Q {
    Q::from_integer(n.clone())
}

/// Sign decided by a rational enclosure of sqrt5 of the given precision;
/// `None` when the enclosure straddles zero.
fn enclosed_sign(x: &QSqrt5, bits: u32) -> Option<Ordering> {
    let root = sqrt5_enclosure(bits);
    let value = &Interval::point(x.rational_part().clone())
        + &(&Interval::point(x.sqrt5_part().clone()) * &root);
    if value.lo() > &q(0) {
        Some(Ordering::Greater)
    } else if value.hi() < &q(0) {
        Some(Ordering::Less)
    } else if x.is_zero() {
        Some(Ordering::Equal)
    } else {
        None
    }
}

// ------------------------------------------------------------------ rationals

#[test]
fn rational_constructors_reduce_and_normalise() {
    assert_eq!(frac(-3, -6), frac(1, 2));
    assert_eq!(frac(4, -6), frac(-2, 3));
    assert_eq!(frac(0, -5), q(0));
    assert_eq!(frac(10, 5), q(2));
    assert!(frac(7, -3).denom() > &BigInt::from(0));
    assert_eq!(q(i64::MIN) - q(1), big(&(BigInt::from(i64::MIN) - 1)));
}

#[test]
#[should_panic(expected = "zero denominator")]
fn literal_fraction_with_zero_denominator_panics() {
    frac(1, 0);
}

#[test]
fn canonical_rational_spelling_round_trips() {
    let mut rng = Rng(1);
    for _ in 0..2000 {
        let x = rng.rational();
        assert_eq!(parse_rational(&x.to_string()), Ok(x));
    }
    let n = (BigInt::one() << 4096usize) + BigInt::from(7);
    let d = (BigInt::one() << 4000usize) + BigInt::from(3);
    let ratio = Q::new(n.clone(), d.clone());
    assert_eq!(parse_rational(&format!("-{n}/{d}")), Ok(-ratio.clone()));
    assert_eq!(parse_rational(&ratio.to_string()), Ok(ratio));
    for (text, value) in [("0", q(0)), ("-7", q(-7)), ("3/4", frac(3, 4)), ("-3/4", frac(-3, 4))] {
        assert_eq!(parse_rational(text), Ok(value));
    }
}

#[test]
fn rational_parsing_refuses_every_other_spelling() {
    let malformed = [
        "", "/", "1/", "/2", "1.5", "1e3", "a", "1/2/3", " 1", "1 ", "--1", "+1", "1/+2", "0x1",
        "１",
    ];
    for text in malformed {
        assert_eq!(
            parse_rational(text),
            Err(ExactError::MalformedRational(text.to_string())),
            "{text:?}"
        );
    }
    for text in ["1/0", "0/0", "-5/0", "3/00"] {
        assert_eq!(parse_rational(text), Err(ExactError::ZeroDenominator(text.to_string())));
    }
    for text in ["2/4", "-3/-6", "1/1", "4/2", "01", "-0", "0/5", "3/-4", "00", "-0/1", "6/1"] {
        assert_eq!(
            parse_rational(text),
            Err(ExactError::NonCanonicalRational(text.to_string())),
            "{text:?}"
        );
    }
}

// ------------------------------------------------------------------ Q(sqrt5)

#[test]
fn field_identities_hold_exactly() {
    let mut rng = Rng(2);
    let (zero, one, root) = (QSqrt5::zero(), QSqrt5::one(), QSqrt5::sqrt5());
    assert_eq!(&root * &root, QSqrt5::integer(5));
    for _ in 0..500 {
        let (x, y, z) = (rng.number(), rng.number(), rng.number());
        assert_eq!(&x + &y, &y + &x);
        assert_eq!(&x * &y, &y * &x);
        assert_eq!(&(&x + &y) + &z, &x + &(&y + &z));
        assert_eq!(&(&x * &y) * &z, &x * &(&y * &z));
        assert_eq!(&x * &(&y + &z), &(&x * &y) + &(&x * &z));
        assert_eq!(&x - &y, &x + &(-&y));
        assert!((&x - &x).is_zero());
        assert_eq!(&x + &zero, x);
        assert_eq!(&x * &one, x);
        assert_eq!(&x * &x.conjugate(), QSqrt5::from_rational(x.norm()));
        assert_eq!((&x * &y).norm(), x.norm() * y.norm());
        assert_eq!((&x * &y).conjugate(), &x.conjugate() * &y.conjugate());
        let r = rng.rational();
        assert_eq!(x.scale(&r), &x * &QSqrt5::from(r.clone()));
        assert_eq!(QSqrt5::new(r.clone(), q(0)), QSqrt5::from_rational(r));
        // A nonzero element has nonzero norm (sqrt5 is irrational).
        assert_eq!(x.is_zero(), x.norm() == q(0));
    }
}

#[test]
fn integer_sign_rule_matches_enclosures_on_a_full_grid() {
    for a in -60i64..=60 {
        for b in -60i64..=60 {
            let x = QSqrt5::new(q(a), q(b));
            let expected = enclosed_sign(&x, 64).unwrap();
            assert_eq!(integer_sqrt5_sign(&BigInt::from(a), &BigInt::from(b)), expected, "{a} {b}");
            assert_eq!(x.sign(), expected);
        }
    }
}

#[test]
fn sign_survives_arbitrarily_close_cancellation() {
    // 9 + 4·sqrt5 has norm 1, so each power (a + b·sqrt5) has a² - 5b² = 1
    // and a - b·sqrt5 = 1/(a + b·sqrt5) is positive and extremely small.
    let (mut a, mut b) = (BigInt::from(9), BigInt::from(4));
    for _ in 0..200 {
        assert_eq!(&a * &a - BigInt::from(5) * &b * &b, BigInt::one());
        let tiny = QSqrt5::new(big(&a), -big(&b));
        assert_eq!(tiny.sign(), Ordering::Greater);
        assert_eq!((-&tiny).sign(), Ordering::Less);
        assert_eq!(integer_sqrt5_sign(&a, &-&b), Ordering::Greater);
        assert_eq!(integer_sqrt5_sign(&-&a, &b), Ordering::Less);
        // Scaled far down, and shifted by less than its size, it keeps its sign.
        let scale = Q::new(BigInt::one(), BigInt::one() << 5000usize);
        assert_eq!(tiny.scale(&scale).sign(), Ordering::Greater);
        let shifted = &tiny - &QSqrt5::from_rational(Q::new(BigInt::one(), &a * BigInt::from(3)));
        assert_eq!(shifted.sign(), Ordering::Greater);
        let overshot = &tiny - &QSqrt5::from_rational(Q::new(BigInt::one(), a.clone()));
        assert_eq!(overshot.sign(), Ordering::Less);
        (a, b) = (
            BigInt::from(9) * &a + BigInt::from(20) * &b,
            BigInt::from(4) * &a + BigInt::from(9) * &b,
        );
    }
    // Norm -1 units: (2 + sqrt5)^k.
    let unit = QSqrt5::new(q(2), q(1));
    let mut power = QSqrt5::one();
    for k in 1..=150 {
        power = &power * &unit;
        let small = power.conjugate();
        let expected = if k % 2 == 0 { Ordering::Greater } else { Ordering::Less };
        assert_eq!(small.sign(), expected, "k = {k}");
        assert_eq!(power.sign(), Ordering::Greater);
    }
}

#[test]
fn order_is_total_and_consistent_with_enclosures() {
    let mut rng = Rng(3);
    let mut values: Vec<QSqrt5> = (0..300).map(|_| rng.number()).collect();
    values.push(QSqrt5::zero());
    values.push(values[0].clone());
    values.sort();
    for pair in values.windows(2) {
        assert!(pair[0] <= pair[1]);
        let difference = &pair[1] - &pair[0];
        assert_eq!(enclosed_sign(&difference, 128).unwrap() != Ordering::Less, true);
        assert_eq!(pair[0].cmp(&pair[1]) == Ordering::Equal, pair[0] == pair[1]);
    }
    for _ in 0..500 {
        let (x, y) = (rng.number(), rng.number());
        assert_eq!(x.cmp(&y), (&x - &y).sign());
        assert_eq!(x.cmp(&y), y.cmp(&x).reverse());
        assert_eq!(Some(x.cmp(&y)), enclosed_sign(&(&x - &y), 128));
    }
}

#[test]
fn rational_sign_uses_cross_multiplied_denominators() {
    // a = p/r and b = u/v with large coprime denominators on both sides.
    let r = (BigInt::one() << 300usize) + BigInt::from(1);
    let v = (BigInt::one() << 299usize) + BigInt::from(3);
    let (mut p, mut u) = (BigInt::from(9), BigInt::from(4));
    for _ in 0..60 {
        // p/r - (u/v)·sqrt5 with p/u close to sqrt5 and r/v close to 2.
        let x = QSqrt5::new(Q::new(p.clone(), r.clone()), -Q::new(u.clone(), v.clone()));
        let expected = integer_sqrt5_sign(&(&p * &v), &-(&u * &r));
        assert_eq!(x.sign(), expected);
        assert_eq!(enclosed_sign(&x, 4096), Some(expected));
        (p, u) = (
            BigInt::from(9) * &p + BigInt::from(20) * &u,
            BigInt::from(4) * &p + BigInt::from(9) * &u,
        );
    }
}

// ------------------------------------------------------------------ intervals

#[test]
fn intervals_refuse_reversal_and_accept_points() {
    assert_eq!(
        Interval::new(q(1), q(0)),
        Err(ExactError::ReversedInterval { lo: q(1), hi: q(0) })
    );
    let tiny = Q::new(BigInt::one(), BigInt::one() << 3000usize);
    assert!(Interval::new(tiny.clone(), q(0)).is_err());
    assert_eq!(Interval::new(q(3), q(3)), Ok(Interval::point(q(3))));
    let a = Interval::new(q(1), q(1) + &tiny).unwrap();
    assert_eq!(a.width(), tiny);
    assert!(a.contains(&q(1)) && !a.contains(&(q(1) + &tiny + &tiny)));
}

#[test]
fn interval_operations_contain_all_sampled_values() {
    let mut rng = Rng(4);
    for _ in 0..600 {
        let (a, b) = (rng.interval(), rng.interval());
        for i in 0..=4 {
            let x = a.lo() + a.width() * frac(i, 4);
            assert!((-&a).contains(&-&x));
            for exponent in 0..=5 {
                assert!(interval_power(&a, exponent).contains(&rational_power(&x, exponent)));
            }
            for j in 0..=4 {
                let y = b.lo() + b.width() * frac(j, 4);
                assert!((&a + &b).contains(&(&x + &y)));
                assert!((&a - &b).contains(&(&x - &y)));
                assert!((&a * &b).contains(&(&x * &y)));
            }
        }
        // Tightness: every end is attained at an end of the inputs.
        let ends = |i: &Interval| [i.lo().clone(), i.hi().clone()];
        let products: Vec<Q> = ends(&a)
            .iter()
            .flat_map(|x| ends(&b).map(|y| x * &y))
            .collect();
        let product = &a * &b;
        assert!(products.contains(product.lo()) && products.contains(product.hi()));
        assert_eq!((&a + &b).lo(), &(a.lo() + b.lo()));
        assert_eq!((&a - &b).hi(), &(a.hi() - b.lo()));
        let [left, right] = a.bisect();
        assert_eq!((left.lo(), left.hi(), right.hi()), (a.lo(), right.lo(), a.hi()));
        assert_eq!(left.width(), right.width());
        assert_eq!(left.hi(), &a.midpoint());
    }
}

#[test]
fn interval_powers_are_exact_ranges() {
    for (lo, hi) in [(-3, -1), (-2, 3), (-3, 2), (0, 2), (1, 4), (0, 0), (-1, 0), (2, 2)] {
        let axis = Interval::new(q(lo), q(hi)).unwrap();
        for exponent in 0..=9u32 {
            let range = interval_power(&axis, exponent);
            let samples: Vec<Q> = (0..=60)
                .map(|k| rational_power(&(q(lo) + axis.width() * frac(k, 60)), exponent))
                .collect();
            // Both ends are attained by samples (0, the ends and the grid
            // include every extremum of a power).
            assert!(samples.contains(range.lo()), "{lo} {hi} {exponent}");
            assert!(samples.contains(range.hi()), "{lo} {hi} {exponent}");
            assert!(samples.iter().all(|v| range.contains(v)));
        }
    }
}

#[test]
fn sqrt5_enclosures_are_proved_by_rational_squares() {
    assert_eq!(sqrt5_enclosure(0), Interval::new(q(2), q(3)).unwrap());
    for bits in [0, 1, 2, 8, 31, 64, 127, 256, 1024, MAX_SQRT5_BITS] {
        let r = sqrt5_enclosure(bits);
        assert!(r.lo() * r.lo() < q(5), "{bits}");
        assert!(r.hi() * r.hi() > q(5), "{bits}");
        assert_eq!(r.width(), Q::new(BigInt::one(), BigInt::one() << bits));
    }
    // Nested: a finer enclosure lies inside a coarser one.
    let (coarse, fine) = (sqrt5_enclosure(40), sqrt5_enclosure(41));
    assert!(coarse.contains(fine.lo()) && coarse.contains(fine.hi()));
}

#[test]
fn errors_have_readable_messages() {
    for error in [
        ExactError::MalformedRational("x".into()),
        ExactError::ZeroDenominator("1/0".into()),
        ExactError::NonCanonicalRational("2/4".into()),
        ExactError::ReversedInterval { lo: q(1), hi: q(0) },
    ] {
        assert!(!error.to_string().is_empty());
    }
}
