//! Polynomials over Q(sqrt5) in a fixed number of variables: exact
//! arithmetic, evaluation, substitution, exact monomial division, tensor
//! Bernstein coefficients on closed rational boxes and the certifier of
//! strict negativity on such a box.
use crate::arithmetic::exact::QSqrt5;
use std::collections::{btree_map::Entry, BTreeMap};
use std::fmt;
use std::ops::{Add, Neg, Sub};

mod bernstein;
pub use bernstein::{BernsteinCoefficients, Refusal, MAX_BERNSTEIN_COEFFICIENTS};

/// The number of variables of every polynomial.
pub const VARIABLES: usize = 5;

/// The exponent of each variable in a monomial. Also used for multi-indices
/// of Bernstein coefficients, which are bounded by the degrees.
pub type Exponents = [u8; VARIABLES];

/// Why a polynomial operation was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PolynomialError {
    /// A product would raise the exponent of `variable` above `u8::MAX`.
    ExponentOverflow { variable: usize },
    /// The monomial `term` of the dividend is not a multiple of `factor`.
    NotDivisible {
        term: Exponents,
        factor: Exponents,
    },
    /// The tensor of Bernstein coefficients for these per-variable degrees
    /// would exceed [`MAX_BERNSTEIN_COEFFICIENTS`] entries.
    TooManyBernsteinCoefficients { degrees: Exponents },
}

impl fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ExponentOverflow { variable } => {
                write!(f, "exponent of variable {variable} exceeds {}", u8::MAX)
            }
            Self::NotDivisible { term, factor } => {
                write!(f, "monomial {term:?} is not divisible by {factor:?}")
            }
            Self::TooManyBernsteinCoefficients { degrees } => write!(
                f,
                "degrees {degrees:?} need more than {MAX_BERNSTEIN_COEFFICIENTS} \
                 Bernstein coefficients"
            ),
        }
    }
}

impl std::error::Error for PolynomialError {}

/// A polynomial with coefficients in Q(sqrt5), stored as its nonzero terms.
/// Every term is kept: there is no truncation of any degree.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Polynomial {
    /// Invariant: no stored coefficient is zero.
    terms: BTreeMap<Exponents, QSqrt5>,
}

impl Polynomial {
    pub fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }
    pub fn constant(value: QSqrt5) -> Self {
        Self::monomial([0; VARIABLES], value)
    }
    /// The single term `coefficient · x^exponents` (zero if the coefficient is).
    pub fn monomial(exponents: Exponents, coefficient: QSqrt5) -> Self {
        let mut p = Self::zero();
        p.accumulate(exponents, coefficient);
        p
    }
    /// The coordinate polynomials `x0, …, x4`.
    pub fn variables() -> [Self; VARIABLES] {
        std::array::from_fn(|j| {
            let mut exponents = [0; VARIABLES];
            exponents[j] = 1;
            Self::monomial(exponents, QSqrt5::one())
        })
    }
    /// The nonzero terms in increasing lexicographic order of exponents.
    pub fn terms(&self) -> impl Iterator<Item = (&Exponents, &QSqrt5)> {
        self.terms.iter()
    }
    /// The coefficient of one monomial (zero when absent).
    pub fn coefficient(&self, exponents: &Exponents) -> QSqrt5 {
        self.terms.get(exponents).cloned().unwrap_or_else(QSqrt5::zero)
    }
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    /// The degree in each variable separately (all zero for a constant or the
    /// zero polynomial).
    pub fn degrees(&self) -> Exponents {
        let mut degrees = [0; VARIABLES];
        for exponents in self.terms.keys() {
            for j in 0..VARIABLES {
                degrees[j] = degrees[j].max(exponents[j]);
            }
        }
        degrees
    }
    /// Add one term, collecting like terms and dropping exact cancellations.
    fn accumulate(&mut self, exponents: Exponents, coefficient: QSqrt5) {
        if coefficient.is_zero() {
            return;
        }
        match self.terms.entry(exponents) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                let sum = entry.get() + &coefficient;
                if sum.is_zero() {
                    entry.remove();
                } else {
                    *entry.get_mut() = sum;
                }
            }
        }
    }
    /// Multiplication by a constant.
    pub fn scale(&self, factor: &QSqrt5) -> Self {
        let mut p = Self::zero();
        for (exponents, coefficient) in &self.terms {
            p.accumulate(*exponents, coefficient * factor);
        }
        p
    }
    pub fn mul(&self, rhs: &Self) -> Result<Self, PolynomialError> {
        let mut p = Self::zero();
        for (m, c) in &self.terms {
            for (n, d) in &rhs.terms {
                p.accumulate(multiply_monomials(m, n)?, c * d);
            }
        }
        Ok(p)
    }
    /// The `exponent`-th power; the zeroth power of every polynomial is 1.
    /// A `u8` exponent suffices: any higher power of a nonconstant polynomial
    /// overflows the exponents, and it bounds the work for constants.
    pub fn pow(&self, exponent: u8) -> Result<Self, PolynomialError> {
        let mut result = Self::constant(QSqrt5::one());
        for _ in 0..exponent {
            result = result.mul(self)?;
        }
        Ok(result)
    }
    /// The exact value at a point of Q(sqrt5)^5.
    pub fn evaluate(&self, point: &[QSqrt5; VARIABLES]) -> QSqrt5 {
        let degrees = self.degrees();
        let powers: [Vec<QSqrt5>; VARIABLES] = std::array::from_fn(|j| {
            let mut powers = vec![QSqrt5::one()];
            for k in 1..=usize::from(degrees[j]) {
                let next = &powers[k - 1] * &point[j];
                powers.push(next);
            }
            powers
        });
        let mut sum = QSqrt5::zero();
        for (exponents, coefficient) in &self.terms {
            let mut term = coefficient.clone();
            for j in 0..VARIABLES {
                if exponents[j] > 0 {
                    term = &term * &powers[j][usize::from(exponents[j])];
                }
            }
            sum = &sum + &term;
        }
        sum
    }
    /// The composition `self(images[0], …, images[4])`: every variable is
    /// replaced by a polynomial (the pull-back along a polynomial map).
    pub fn substitute(&self, images: &[Self; VARIABLES]) -> Result<Self, PolynomialError> {
        let degrees = self.degrees();
        let mut powers: [Vec<Self>; VARIABLES] = std::array::from_fn(|_| Vec::new());
        for j in 0..VARIABLES {
            powers[j].push(Self::constant(QSqrt5::one()));
            for k in 1..=usize::from(degrees[j]) {
                let next = powers[j][k - 1].mul(&images[j])?;
                powers[j].push(next);
            }
        }
        let mut result = Self::zero();
        for (exponents, coefficient) in &self.terms {
            let mut term = Self::constant(coefficient.clone());
            for j in 0..VARIABLES {
                if exponents[j] > 0 {
                    term = term.mul(&powers[j][usize::from(exponents[j])])?;
                }
            }
            for (exponents, coefficient) in term.terms {
                result.accumulate(exponents, coefficient);
            }
        }
        Ok(result)
    }
    /// The exact quotient by the monomial `x^factor`. Refused, naming the
    /// first offending term, unless every term is a multiple of it.
    pub fn divide_by_monomial(&self, factor: &Exponents) -> Result<Self, PolynomialError> {
        let mut terms = BTreeMap::new();
        for (exponents, coefficient) in &self.terms {
            let mut quotient = [0; VARIABLES];
            for j in 0..VARIABLES {
                quotient[j] = exponents[j].checked_sub(factor[j]).ok_or(
                    PolynomialError::NotDivisible {
                        term: *exponents,
                        factor: *factor,
                    },
                )?;
            }
            terms.insert(quotient, coefficient.clone());
        }
        Ok(Self { terms })
    }
}

fn multiply_monomials(m: &Exponents, n: &Exponents) -> Result<Exponents, PolynomialError> {
    let mut product = [0; VARIABLES];
    for j in 0..VARIABLES {
        product[j] = m[j]
            .checked_add(n[j])
            .ok_or(PolynomialError::ExponentOverflow { variable: j })?;
    }
    Ok(product)
}

impl Add for &Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: &Polynomial) -> Polynomial {
        let mut p = self.clone();
        for (exponents, coefficient) in &rhs.terms {
            p.accumulate(*exponents, coefficient.clone());
        }
        p
    }
}

impl Sub for &Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: &Polynomial) -> Polynomial {
        let mut p = self.clone();
        for (exponents, coefficient) in &rhs.terms {
            p.accumulate(*exponents, -coefficient);
        }
        p
    }
}

impl Neg for &Polynomial {
    type Output = Polynomial;
    fn neg(self) -> Polynomial {
        Polynomial {
            terms: self
                .terms
                .iter()
                .map(|(exponents, coefficient)| (*exponents, -coefficient))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests;
