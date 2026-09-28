//! Exact numbers: validated rationals, the real field Q(√5) with exact sign
//! and order, and closed rational intervals. No floating-point arithmetic.
//!
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, Zero};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// An arbitrary-size rational number.
///
/// **Invariant:** the value is stored in lowest terms with a positive
/// denominator. The field is private, and every value is either a result of
/// num-rational's `Ratio::new`, which reduces any integer pair (the
/// constructors and `+`, `-`, `*`, `/`), or the negation of a value with the
/// invariant (which keeps it). Every value also passes through a private
/// wrapper that asserts the positive denominator. Sign and order decisions in
/// this crate read numerators and denominators directly and rely on this.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Q(BigRational);

impl Q {
    /// Wraps a result of num-rational. The positive denominator, the part of
    /// the invariant that signs rely on, is checked in every build.
    fn checked(value: BigRational) -> Self {
        assert!(value.denom().is_positive(), "rational with a non-positive denominator");
        Self(value)
    }
    /// `numerator / denominator` in lowest terms. Panics when the denominator
    /// is zero, like integer division by zero.
    pub fn new(numerator: BigInt, denominator: BigInt) -> Self {
        assert!(!denominator.is_zero(), "rational with a zero denominator");
        Self::checked(BigRational::new(numerator, denominator))
    }
    pub fn from_integer(n: BigInt) -> Self {
        Self::checked(BigRational::from_integer(n))
    }
    pub fn zero() -> Self {
        q(0)
    }
    pub fn one() -> Self {
        q(1)
    }
    /// The numerator of the lowest-terms representation.
    pub fn numer(&self) -> &BigInt {
        self.0.numer()
    }
    /// The denominator of the lowest-terms representation; always positive.
    pub fn denom(&self) -> &BigInt {
        self.0.denom()
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

/// The integer `n` as a rational.
pub fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

/// The rational `n / d`, reduced. Intended for literal constants: it panics
/// when `d` is zero, like integer division by zero.
pub fn frac(n: i64, d: i64) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

/// The canonical spelling: `n` for an integer, otherwise `n/d` with `d > 1`.
impl fmt::Display for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl fmt::Debug for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// `+`, `-`, `*` and `/` for every combination of owned and borrowed
/// operands. num-rational ends each of them with `Ratio::new`; division by
/// zero panics like integer division by zero.
macro_rules! binary_operation {
    ($trait:ident, $method:ident) => {
        impl $trait<&Q> for &Q {
            type Output = Q;
            fn $method(self, rhs: &Q) -> Q {
                Q::checked($trait::$method(&self.0, &rhs.0))
            }
        }
        impl $trait<Q> for Q {
            type Output = Q;
            fn $method(self, rhs: Q) -> Q {
                Q::checked($trait::$method(self.0, rhs.0))
            }
        }
        impl $trait<&Q> for Q {
            type Output = Q;
            fn $method(self, rhs: &Q) -> Q {
                Q::checked($trait::$method(self.0, &rhs.0))
            }
        }
        impl $trait<Q> for &Q {
            type Output = Q;
            fn $method(self, rhs: Q) -> Q {
                Q::checked($trait::$method(&self.0, rhs.0))
            }
        }
    };
}

binary_operation!(Add, add);
binary_operation!(Sub, sub);
binary_operation!(Mul, mul);
binary_operation!(Div, div);

/// `+=`, `-=`, `*=` and `/=` through the binary operations above.
macro_rules! assign_operation {
    ($trait:ident, $method:ident, $operation:ident) => {
        impl $trait<&Q> for Q {
            fn $method(&mut self, rhs: &Q) {
                *self = std::mem::replace(self, Q::zero()).$operation(rhs);
            }
        }
        impl $trait<Q> for Q {
            fn $method(&mut self, rhs: Q) {
                *self = std::mem::replace(self, Q::zero()).$operation(rhs);
            }
        }
    };
}

assign_operation!(AddAssign, add_assign, add);
assign_operation!(SubAssign, sub_assign, sub);
assign_operation!(MulAssign, mul_assign, mul);
assign_operation!(DivAssign, div_assign, div);

impl Neg for Q {
    type Output = Q;
    fn neg(self) -> Q {
        Q::checked(-self.0)
    }
}

impl Neg for &Q {
    type Output = Q;
    fn neg(self) -> Q {
        Q::checked(-&self.0)
    }
}

/// Why an exact value was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExactError {
    /// The text is not an optionally signed integer, or two such separated
    /// by one `/`.
    MalformedRational(String),
    /// The text names a rational with a zero denominator.
    ZeroDenominator(String),
    /// The text names a rational but is not its canonical spelling
    /// (`n` for integers, otherwise `n/d` in lowest terms with `d > 1`).
    NonCanonicalRational(String),
    /// A closed interval was requested with its lower end above its upper end.
    ReversedInterval { lo: Q, hi: Q },
}

impl fmt::Display for ExactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MalformedRational(text) => write!(f, "malformed rational {text:?}"),
            Self::ZeroDenominator(text) => write!(f, "zero denominator in {text:?}"),
            Self::NonCanonicalRational(text) => write!(f, "non-canonical rational {text:?}"),
            Self::ReversedInterval { lo, hi } => write!(f, "reversed interval [{lo}, {hi}]"),
        }
    }
}

impl std::error::Error for ExactError {}

/// Parse the canonical decimal spelling of a rational: `n` for an integer,
/// otherwise `n/d` in lowest terms with `d > 1`; the only sign allowed is a
/// leading `-` on `n`. Every other spelling of a valid rational is refused,
/// so parsing and printing are inverse to each other.
pub fn parse_rational(text: &str) -> Result<Q, ExactError> {
    let (numerator, denominator) = text.split_once('/').unwrap_or((text, "1"));
    let integer = |part: &str| -> Result<BigInt, ExactError> {
        let digits = part.strip_prefix('-').unwrap_or(part);
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(ExactError::MalformedRational(text.to_string()));
        }
        part.parse()
            .map_err(|_| ExactError::MalformedRational(text.to_string()))
    };
    let (numerator, denominator) = (integer(numerator)?, integer(denominator)?);
    if denominator.is_zero() {
        return Err(ExactError::ZeroDenominator(text.to_string()));
    }
    let value = Q::new(numerator, denominator);
    if value.to_string() != text {
        return Err(ExactError::NonCanonicalRational(text.to_string()));
    }
    Ok(value)
}

/// The real number `a + b·√5` with rational `a` and `b`. Since `√5` is
/// irrational the pair is unique, so structural equality is numeric equality,
/// and the order is exact.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct QSqrt5 {
    rational: Q,
    sqrt5: Q,
}

impl QSqrt5 {
    /// The number `rational + sqrt5·√5`.
    pub fn new(rational: Q, sqrt5: Q) -> Self {
        Self { rational, sqrt5 }
    }
    pub fn from_rational(rational: Q) -> Self {
        Self::new(rational, Q::zero())
    }
    pub fn integer(n: i64) -> Self {
        Self::from_rational(q(n))
    }
    pub fn zero() -> Self {
        Self::integer(0)
    }
    pub fn one() -> Self {
        Self::integer(1)
    }
    /// The number √5 itself.
    pub fn sqrt5() -> Self {
        Self::new(Q::zero(), Q::one())
    }
    /// The rational part `a` of `a + b·√5`.
    pub fn rational_part(&self) -> &Q {
        &self.rational
    }
    /// The coefficient `b` of √5 in `a + b·√5`.
    pub fn sqrt5_part(&self) -> &Q {
        &self.sqrt5
    }
    pub fn is_zero(&self) -> bool {
        self.rational.is_zero() && self.sqrt5.is_zero()
    }
    /// The Galois conjugate `a - b·√5`.
    pub fn conjugate(&self) -> Self {
        Self::new(self.rational.clone(), -&self.sqrt5)
    }
    /// The field norm `a² - 5b²`, the product with the conjugate.
    pub fn norm(&self) -> Q {
        &self.rational * &self.rational - q(5) * &self.sqrt5 * &self.sqrt5
    }
    /// Multiplication by a rational.
    pub fn scale(&self, factor: &Q) -> Self {
        Self::new(&self.rational * factor, &self.sqrt5 * factor)
    }
    /// The exact sign, as the comparison with zero.
    pub fn sign(&self) -> Ordering {
        sign_of_sum(&self.rational, &self.sqrt5)
    }
}

/// Sign of `a + b·√5` for rationals `a = p/r`, `b = u/v`: the denominators
/// `r` and `v` are positive (the invariant of [`Q`]), so multiplying by
/// `rv > 0` gives the integers `pv + ur·√5` with the same sign.
fn sign_of_sum(a: &Q, b: &Q) -> Ordering {
    integer_sqrt5_sign(&(a.numer() * b.denom()), &(b.numer() * a.denom()))
}

/// The exact sign of `a + b·√5` for integers `a` and `b`: the only sign rule
/// of the crate. With equal or zero signs it is immediate; with opposite
/// signs `a` decides exactly when `a² > 5b²` (equality would make √5
/// rational unless both vanish).
pub(crate) fn integer_sqrt5_sign(a: &BigInt, b: &BigInt) -> Ordering {
    let (sa, sb) = (a.cmp(&BigInt::zero()), b.cmp(&BigInt::zero()));
    if sb == Ordering::Equal || sa == sb {
        return sa;
    }
    if sa == Ordering::Equal {
        return sb;
    }
    match (a * a).cmp(&(BigInt::from(5) * b * b)) {
        Ordering::Greater => sa,
        Ordering::Less => sb,
        Ordering::Equal => unreachable!("a² = 5b² with b ≠ 0 would make √5 rational"),
    }
}

impl Ord for QSqrt5 {
    fn cmp(&self, other: &Self) -> Ordering {
        sign_of_sum(
            &(&self.rational - &other.rational),
            &(&self.sqrt5 - &other.sqrt5),
        )
    }
}

impl PartialOrd for QSqrt5 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl From<Q> for QSqrt5 {
    fn from(rational: Q) -> Self {
        Self::from_rational(rational)
    }
}

impl Add for &QSqrt5 {
    type Output = QSqrt5;
    fn add(self, rhs: &QSqrt5) -> QSqrt5 {
        QSqrt5::new(&self.rational + &rhs.rational, &self.sqrt5 + &rhs.sqrt5)
    }
}

impl Sub for &QSqrt5 {
    type Output = QSqrt5;
    fn sub(self, rhs: &QSqrt5) -> QSqrt5 {
        QSqrt5::new(&self.rational - &rhs.rational, &self.sqrt5 - &rhs.sqrt5)
    }
}

impl Mul for &QSqrt5 {
    type Output = QSqrt5;
    fn mul(self, rhs: &QSqrt5) -> QSqrt5 {
        QSqrt5::new(
            &self.rational * &rhs.rational + q(5) * &self.sqrt5 * &rhs.sqrt5,
            &self.rational * &rhs.sqrt5 + &self.sqrt5 * &rhs.rational,
        )
    }
}

impl Neg for &QSqrt5 {
    type Output = QSqrt5;
    fn neg(self) -> QSqrt5 {
        QSqrt5::new(-&self.rational, -&self.sqrt5)
    }
}

/// A closed interval `[lo, hi]` with rational ends, `lo <= hi`. A point is
/// the interval with equal ends.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Interval {
    lo: Q,
    hi: Q,
}

impl Interval {
    pub fn new(lo: Q, hi: Q) -> Result<Self, ExactError> {
        if lo > hi {
            return Err(ExactError::ReversedInterval { lo, hi });
        }
        Ok(Self { lo, hi })
    }
    pub fn point(x: Q) -> Self {
        Self {
            lo: x.clone(),
            hi: x,
        }
    }
    /// Construction from ends already known to be ordered.
    fn ordered(lo: Q, hi: Q) -> Self {
        debug_assert!(lo <= hi);
        Self { lo, hi }
    }
    pub fn lo(&self) -> &Q {
        &self.lo
    }
    pub fn hi(&self) -> &Q {
        &self.hi
    }
    pub fn width(&self) -> Q {
        &self.hi - &self.lo
    }
    pub fn midpoint(&self) -> Q {
        (&self.lo + &self.hi) / q(2)
    }
    pub fn contains(&self, x: &Q) -> bool {
        &self.lo <= x && x <= &self.hi
    }
    /// The two closed halves at the exact midpoint; they share the midpoint
    /// and their union is `self`.
    pub fn bisect(&self) -> [Self; 2] {
        let middle = self.midpoint();
        [
            Self::ordered(self.lo.clone(), middle.clone()),
            Self::ordered(middle, self.hi.clone()),
        ]
    }
}

impl Add for &Interval {
    type Output = Interval;
    fn add(self, rhs: &Interval) -> Interval {
        Interval::ordered(&self.lo + &rhs.lo, &self.hi + &rhs.hi)
    }
}

impl Sub for &Interval {
    type Output = Interval;
    fn sub(self, rhs: &Interval) -> Interval {
        Interval::ordered(&self.lo - &rhs.hi, &self.hi - &rhs.lo)
    }
}

impl Neg for &Interval {
    type Output = Interval;
    fn neg(self) -> Interval {
        Interval::ordered(-&self.hi, -&self.lo)
    }
}

impl Mul for &Interval {
    type Output = Interval;
    fn mul(self, rhs: &Interval) -> Interval {
        let products = [
            &self.lo * &rhs.lo,
            &self.lo * &rhs.hi,
            &self.hi * &rhs.lo,
            &self.hi * &rhs.hi,
        ];
        let lo = products.iter().min().expect("four products").clone();
        let hi = products.iter().max().expect("four products").clone();
        Interval::ordered(lo, hi)
    }
}

#[cfg(test)]
pub(crate) mod tests;
