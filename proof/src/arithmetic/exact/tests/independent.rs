//! Checks of `exact` against independent reimplementations
//! of the arithmetic layer: signs in Q(√5) decided by escalating rational
//! enclosures of √5 (never by the crate's integer sign rule), the canonical
//! rational spelling by a separately written grammar, interval powers by
//! brute force, and the normal form of every rational.
use super::{interval_power, rational_power};
use crate::arithmetic::exact::*;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;

/// xorshift64*, a generator unrelated to the other tests' splitmix64.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo + 1) as u64) as i64
    }
    fn chance(&mut self, k: u64, n: u64) -> bool {
        self.below(n) < k
    }
    fn big(&mut self, bits: u64) -> BigInt {
        let mut x = BigInt::zero();
        for _ in 0..bits.div_ceil(32) {
            x = (x << 32usize) + BigInt::from(self.next() >> 32);
        }
        x >> ((bits.div_ceil(32) * 32 - bits) as usize)
    }
}

fn iv(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

/// Independent sign of `a + b sqrt5`: escalating rational enclosures of
/// sqrt5 from an integer square root checked by squares.
pub(crate) fn independent_sign(x: &QSqrt5) -> Ordering {
    let (a, b) = (x.rational_part(), x.sqrt5_part());
    if b.is_zero() {
        return a.cmp(&q(0));
    }
    let mut bits = 64usize;
    loop {
        let d = BigInt::one() << bits;
        let square = BigInt::from(5) * &d * &d;
        let k = square.sqrt();
        assert!(&k * &k <= square && (&k + 1u32) * (&k + 1u32) > square);
        let lo_root = Q::new(k.clone(), d.clone());
        let hi_root = Q::new(k + 1u32, d);
        let (u, v) = (a + b * &lo_root, a + b * &hi_root);
        let (lo, hi) = if u <= v { (u, v) } else { (v, u) };
        if lo > q(0) {
            return Ordering::Greater;
        }
        if hi < q(0) {
            return Ordering::Less;
        }
        bits *= 2;
        assert!(bits <= 1 << 18, "sign undecided at {bits} bits for {x:?}");
    }
}

/// Euclid's algorithm, written here rather than taken from a library.
fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    let (mut a, mut b) = (a.abs(), b.abs());
    while !b.is_zero() {
        let r = &a % &b;
        (a, b) = (b, r);
    }
    a
}

/// Lowest terms with a positive denominator, and the canonical spelling
/// parses back to the same value.
fn assert_normal(x: &Q) {
    assert!(x.denom().is_positive(), "{x:?}");
    assert!(gcd(x.numer(), x.denom()).is_one(), "{x:?}");
    assert_eq!(parse_rational(&x.to_string()).as_ref(), Ok(x));
}

/// The invariant of `Q` cannot be broken from outside: every constructor
/// reduces any integer pair, including negative and unreduced denominators
/// (the spellings that once made the certifier and the sign rule wrong), and
/// every operation keeps the normal form.
#[test]
fn every_rational_is_in_lowest_terms_with_a_positive_denominator() {
    assert_eq!(Q::new(BigInt::from(-1), BigInt::from(-1)), q(1));
    assert_eq!(Q::new(BigInt::from(-1), BigInt::from(-1)).denom(), &BigInt::one());
    assert_eq!(Q::new(BigInt::from(1), BigInt::from(-2)).denom(), &BigInt::from(2));
    assert_eq!(Q::new(BigInt::from(2), BigInt::from(-4)).to_string(), "-1/2");
    let half = QSqrt5::from_rational(Q::new(BigInt::from(1), BigInt::from(-2)));
    assert_eq!(half.sign(), Ordering::Less);
    // -1/2 + √5/5 = -0.052...
    assert_eq!(QSqrt5::new(frac(3, -6), frac(-1, -5)).sign(), Ordering::Less);
    assert_eq!(QSqrt5::new(frac(3, -6), frac(-1, -4)).sign(), Ordering::Greater);
    let mut rng = Rng(0x5eed);
    let mut value = || {
        let n = BigInt::from(rng.int(-10_000, 10_000)) * BigInt::from(rng.int(1, 30));
        let d = BigInt::from(rng.int(1, 5_000)) * BigInt::from(rng.int(1, 30));
        let d = if rng.chance(1, 2) { -d } else { d };
        Q::new(n, d)
    };
    for _ in 0..3000 {
        let (a, b) = (value(), value());
        for x in [&a, &b] {
            assert_normal(x);
        }
        assert_normal(&(&a + &b));
        assert_normal(&(&a - &b));
        assert_normal(&(&a * &b));
        assert_normal(&-&a);
        if !b.is_zero() {
            assert_normal(&(&a / &b));
        }
        let mut c = a.clone();
        c -= &b;
        c *= &b;
        assert_normal(&c);
        let x = QSqrt5::new(a.clone(), b.clone());
        assert_eq!(x.sign(), independent_sign(&x));
        let range = iv(a.clone().min(b.clone()), a.max(b));
        assert_normal(&range.midpoint());
        assert_normal(&range.width());
    }
}

/// Canonical-spelling grammar, written independently: `0`, or an optional
/// `-` then a nonzero digit and digits, optionally `/d` with d > 1 without a
/// leading zero and coprime to the numerator.
fn canonical_by_grammar(text: &str) -> Option<Q> {
    let integer = |s: &str, allow_minus: bool| -> Option<BigInt> {
        let (neg, digits) = match s.strip_prefix('-') {
            Some(rest) if allow_minus => (true, rest),
            Some(_) => return None,
            None => (false, s),
        };
        if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
            return None;
        }
        if digits.len() > 1 && digits.starts_with('0') {
            return None;
        }
        if neg && digits == "0" {
            return None;
        }
        let v: BigInt = digits.parse().ok()?;
        Some(if neg { -v } else { v })
    };
    match text.split_once('/') {
        None => integer(text, true).map(Q::from_integer),
        Some((n, d)) => {
            let (n, d) = (integer(n, true)?, integer(d, false)?);
            if d <= BigInt::one() || n.is_zero() {
                return None;
            }
            let value = Q::new(n, d.clone());
            if value.denom() != &d {
                return None; // not in lowest terms
            }
            Some(value)
        }
    }
}

#[test]
fn parse_rational_agrees_with_an_independent_grammar_on_fuzzed_text() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let alphabet: Vec<char> = "0123456789-/+ 0011".chars().collect();
    let mut accepted = 0;
    for _ in 0..200_000 {
        let length = rng.below(7) as usize;
        let text: String =
            (0..length).map(|_| alphabet[rng.below(alphabet.len() as u64) as usize]).collect();
        let parsed = parse_rational(&text).ok();
        assert_eq!(parsed, canonical_by_grammar(&text), "{text:?}");
        if let Some(v) = parsed {
            assert_eq!(v.to_string(), text);
            accepted += 1;
        }
    }
    assert!(accepted > 1000, "{accepted}");
}

#[test]
fn integer_sign_rule_on_huge_random_and_near_tie_inputs() {
    let mut rng = Rng(77);
    for _ in 0..3000 {
        let bits = 1 + rng.below(700);
        let mut a = rng.big(bits);
        let mut b = rng.big(bits);
        if rng.chance(1, 2) {
            a = -a;
        }
        if rng.chance(1, 2) {
            b = -b;
        }
        let x = QSqrt5::new(Q::from_integer(a.clone()), Q::from_integer(b.clone()));
        assert_eq!(integer_sqrt5_sign(&a, &b), independent_sign(&x));
    }
    // Near ties: a = isqrt(5 b^2) + delta, with b random huge and small delta,
    // rescaled by random rationals (denominators on both parts).
    for _ in 0..2000 {
        let bits = 1 + rng.below(900);
        let b = rng.big(bits) + 1u32;
        let root = (BigInt::from(5) * &b * &b).sqrt();
        let a = &root + BigInt::from(rng.int(-2, 2));
        let (sa, sb) = if rng.chance(1, 2) {
            (-a.clone(), b.clone())
        } else {
            (a.clone(), -b.clone())
        };
        let x = QSqrt5::new(Q::from_integer(sa.clone()), Q::from_integer(sb.clone()));
        assert_eq!(integer_sqrt5_sign(&sa, &sb), independent_sign(&x));
        let r = Q::new(BigInt::from(rng.int(1, 1 << 40)), rng.big(200) + 1u32);
        let scaled = QSqrt5::new(x.rational_part() * &r, x.sqrt5_part() * &r);
        assert_eq!(scaled.sign(), independent_sign(&x));
        // Different denominators on the two parts.
        let odd = QSqrt5::new(
            Q::new(sa.clone() * BigInt::from(7), BigInt::from(21)),
            Q::new(sb.clone() * BigInt::from(11), BigInt::from(33)),
        );
        assert_eq!(odd.sign(), independent_sign(&x));
    }
}

#[test]
fn order_is_transitive_among_near_equal_values() {
    let mut rng = Rng(5);
    let conjugate_phi = QSqrt5::new(frac(1, 2), frac(-1, 2));
    let mut tiny = vec![QSqrt5::one()];
    for _ in 0..80 {
        let next = tiny.last().unwrap() * &conjugate_phi;
        tiny.push(next);
    }
    let values: Vec<QSqrt5> = (0..400)
        .map(|_| {
            let base = QSqrt5::from_rational(frac(rng.int(-3, 3), 7));
            let t = &tiny[rng.below(tiny.len() as u64) as usize];
            let s = frac(rng.int(-5, 5), rng.int(1, 5));
            &base + &t.scale(&s)
        })
        .collect();
    for _ in 0..20_000 {
        let [a, b, c] = [0, 1, 2].map(|_| &values[rng.below(values.len() as u64) as usize]);
        if a <= b && b <= c {
            assert!(a <= c);
        }
        assert_eq!(a.cmp(b), independent_sign(&(a - b)));
    }
}

#[test]
fn interval_powers_by_brute_force_on_fractional_ends() {
    let mut rng = Rng(11);
    for _ in 0..400 {
        let a = frac(rng.int(-40, 40), rng.int(1, 9));
        let b = frac(rng.int(-40, 40), rng.int(1, 9));
        let range = iv(a.clone().min(b.clone()), a.max(b));
        for n in 0..=11u32 {
            let p = interval_power(&range, n);
            // Candidates for the extrema of x^n: the ends and 0 if inside.
            let mut candidates = vec![range.lo().clone(), range.hi().clone()];
            if range.contains(&q(0)) {
                candidates.push(q(0));
            }
            let values: Vec<Q> = candidates.iter().map(|x| rational_power(x, n)).collect();
            assert_eq!(p.lo(), values.iter().min().unwrap(), "{range:?}^{n}");
            assert_eq!(p.hi(), values.iter().max().unwrap(), "{range:?}^{n}");
        }
    }
}
