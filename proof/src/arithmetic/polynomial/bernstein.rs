//! Tensor Bernstein coefficients on a closed rational box, computed in
//! integers with an exactly known positive scale, and the certifier of strict
//! negativity.
use super::{Exponents, Polynomial, PolynomialError, VARIABLES};
use crate::arithmetic::exact::{integer_sqrt5_sign, Interval, QSqrt5, Q};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::fmt;

/// Largest tensor of Bernstein coefficients computed for one polynomial,
/// `2^16`. It bounds the worst case of one call.
pub const MAX_BERNSTEIN_COEFFICIENTS: usize = 1 << 16;

/// Why [`Polynomial::certify_negative`] did not certify.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Refusal {
    /// The Bernstein coefficient at this multi-index is zero or positive
    /// (the first such in the order of [`BernsteinCoefficients::iter`]).
    NonNegativeCoefficient { index: Exponents },
    /// The tensor for these per-variable degrees would exceed
    /// [`MAX_BERNSTEIN_COEFFICIENTS`] entries.
    TooManyCoefficients { degrees: Exponents },
}

impl fmt::Display for Refusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonNegativeCoefficient { index } => {
                write!(f, "the Bernstein coefficient {index:?} is not negative")
            }
            // The same condition as PolynomialError::TooManyBernsteinCoefficients.
            Self::TooManyCoefficients { degrees } => write!(
                f,
                "degrees {degrees:?} need more than {MAX_BERNSTEIN_COEFFICIENTS} \
                 Bernstein coefficients"
            ),
        }
    }
}

impl std::error::Error for Refusal {}

/// The tensor Bernstein coefficients of a polynomial on a closed box, of
/// degree `degrees[j]` (the polynomial's own degree) in variable `j`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BernsteinCoefficients {
    layout: Layout,
    values: Vec<QSqrt5>,
}

impl BernsteinCoefficients {
    pub fn degrees(&self) -> Exponents {
        self.layout.degrees
    }
    /// The number of coefficients, the product of `degrees[j] + 1`.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// The coefficient with multi-index `index` (`index[j] <= degrees[j]`).
    pub fn get(&self, index: &Exponents) -> Option<&QSqrt5> {
        self.layout.offset(index).map(|offset| &self.values[offset])
    }
    /// All coefficients, in lexicographic order of their multi-indices.
    pub fn iter(&self) -> impl Iterator<Item = (Exponents, &QSqrt5)> {
        self.values
            .iter()
            .enumerate()
            .map(|(offset, value)| (self.layout.index(offset), value))
    }
}

impl Polynomial {
    /// The exact tensor Bernstein coefficients on the closed box `cell`
    /// (intervals may have zero width).
    pub fn bernstein_coefficients(
        &self,
        cell: &[Interval; VARIABLES],
    ) -> Result<BernsteinCoefficients, PolynomialError> {
        let scaled = Scaled::new(self, cell)
            .ok_or(PolynomialError::TooManyBernsteinCoefficients { degrees: self.degrees() })?;
        let binomials = pascal(scaled.layout.max_degree());
        let values = scaled
            .values
            .iter()
            .enumerate()
            .map(|(offset, value)| {
                let index = scaled.layout.index(offset);
                let mut scale = scaled.common.clone();
                for j in 0..VARIABLES {
                    let n = usize::from(scaled.layout.degrees[j]);
                    scale *= &binomials[n][usize::from(index[j])];
                    scale *= &scaled.axis_scales[j];
                }
                QSqrt5::new(Q::new(value.a.clone(), scale.clone()), Q::new(value.b.clone(), scale))
            })
            .collect();
        Ok(BernsteinCoefficients {
            layout: scaled.layout,
            values,
        })
    }

    /// Certify that the polynomial is strictly negative on the closed box
    /// `cell`: accepted exactly when every tensor Bernstein coefficient is
    /// strictly negative. A zero or positive coefficient is a refusal.
    pub fn certify_negative(&self, cell: &[Interval; VARIABLES]) -> Result<(), Refusal> {
        let scaled = Scaled::new(self, cell).ok_or(Refusal::TooManyCoefficients {
            degrees: self.degrees(),
        })?;
        match scaled.values.iter().position(|v| v.sign() != Ordering::Less) {
            None => Ok(()),
            Some(offset) => Err(Refusal::NonNegativeCoefficient {
                index: scaled.layout.index(offset),
            }),
        }
    }
}

/// Dense row-major multi-index layout: variable 0 varies slowest.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Layout {
    degrees: Exponents,
    strides: [usize; VARIABLES],
    size: usize,
}

impl Layout {
    /// `None` when the tensor would exceed the size limit.
    fn new(degrees: Exponents) -> Option<Self> {
        let mut strides = [0; VARIABLES];
        let mut size = 1usize;
        for j in (0..VARIABLES).rev() {
            strides[j] = size;
            size = size.checked_mul(usize::from(degrees[j]) + 1)?;
        }
        (size <= MAX_BERNSTEIN_COEFFICIENTS).then_some(Self {
            degrees,
            strides,
            size,
        })
    }
    fn max_degree(&self) -> usize {
        usize::from(self.degrees.into_iter().max().unwrap_or(0))
    }
    fn offset(&self, index: &Exponents) -> Option<usize> {
        let mut offset = 0;
        for j in 0..VARIABLES {
            if index[j] > self.degrees[j] {
                return None;
            }
            offset += usize::from(index[j]) * self.strides[j];
        }
        Some(offset)
    }
    fn index(&self, mut offset: usize) -> Exponents {
        let mut index = [0; VARIABLES];
        for j in 0..VARIABLES {
            let digit = offset / self.strides[j];
            offset %= self.strides[j];
            index[j] = u8::try_from(digit).expect("layout digit bounded by a u8 degree");
        }
        index
    }
}

/// The integer element `a + b·sqrt5`.
#[derive(Clone, Debug, Default)]
struct Pair {
    a: BigInt,
    b: BigInt,
}

impl Pair {
    fn sign(&self) -> Ordering {
        integer_sqrt5_sign(&self.a, &self.b)
    }
}

/// Bernstein coefficients up to positive integer factors: the coefficient
/// with multi-index `I` is `values[I] / (common · Π_j C(n_j, I_j) · axis_scales[j])`.
struct Scaled {
    layout: Layout,
    values: Vec<Pair>,
    common: BigInt,
    axis_scales: [BigInt; VARIABLES],
}

impl Scaled {
    /// `None` when the tensor would exceed the size limit.
    fn new(p: &Polynomial, cell: &[Interval; VARIABLES]) -> Option<Self> {
        let layout = Layout::new(p.degrees())?;
        // A common denominator of all coefficient parts turns every
        // coefficient into an element of Z[sqrt5].
        let mut common = BigInt::one();
        for coefficient in p.terms.values() {
            common = lcm_with_denominator(&common, coefficient.rational_part());
            common = lcm_with_denominator(&common, coefficient.sqrt5_part());
        }
        let mut values = vec![Pair::default(); layout.size];
        for (exponents, coefficient) in &p.terms {
            let offset = layout.offset(exponents).expect("term within its own degrees");
            values[offset] = Pair {
                a: scaled_integer(coefficient.rational_part(), &common),
                b: scaled_integer(coefficient.sqrt5_part(), &common),
            };
        }
        let binomials = pascal(layout.max_degree());
        let mut axis_scales: [BigInt; VARIABLES] = std::array::from_fn(|_| BigInt::one());
        for j in 0..VARIABLES {
            let n = usize::from(layout.degrees[j]);
            if n == 0 {
                continue;
            }
            let axis = AxisMap::new(&cell[j], n, &binomials);
            axis.apply(&mut values, &layout, j);
            axis_scales[j] = axis.scale;
        }
        Some(Self {
            layout,
            values,
            common,
            axis_scales,
        })
    }
}

/// `lcm(accumulated, denominator of x)` for a positive `accumulated`. With
/// `g = gcd(accumulated, den)`, the reduced fraction `accumulated/den` has
/// the positive denominator `den/g` (both invariants of `Q`), and
/// `accumulated · den / g` is the least common multiple, again positive.
fn lcm_with_denominator(accumulated: &BigInt, x: &Q) -> BigInt {
    let reduced = Q::new(accumulated.clone(), x.denom().clone());
    accumulated * reduced.denom()
}

/// `x · common` for a multiple `common` of the denominator of `x`.
fn scaled_integer(x: &Q, common: &BigInt) -> BigInt {
    x.numer() * (common / x.denom())
}

/// Binomial coefficients `C(n, k)` for `n <= max`.
fn pascal(max: usize) -> Vec<Vec<BigInt>> {
    let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
    for n in 1..=max {
        let previous = &rows[n - 1];
        let row = (0..=n)
            .map(|k| {
                let left = if k > 0 { previous[k - 1].clone() } else { BigInt::zero() };
                let right = previous.get(k).cloned().unwrap_or_else(BigInt::zero);
                left + right
            })
            .collect();
        rows.push(row);
    }
    rows
}

/// The integer matrix taking the power coefficients `c_i` of a degree-`n`
/// offset in one variable to its scaled Bernstein coefficients on `[lo, hi]`:
/// row `r` gives `C(n, r) · d^n · b_r`, where `d` is the common denominator of
/// `lo` and `hi` and `b_r` the exact Bernstein coefficient.
struct AxisMap {
    degree: usize,
    /// Row-major `(n+1) × (n+1)`.
    matrix: Vec<BigInt>,
    /// `d^n`.
    scale: BigInt,
}

impl AxisMap {
    /// With `lo = l/d`, `hi = h/d`, `m = h - l`, substituting
    /// `x = (l + m·y)/d` and multiplying by `d^n` gives the power
    /// coefficients `c'_k = Σ_i C(i,k) l^(i-k) m^k d^(n-i) c_i` in `y`. The
    /// Bernstein coefficients on `[0,1]` then satisfy
    /// `C(n,r) b_r = Σ_{k<=r} C(n-k, r-k) c'_k`, all in integers.
    fn new(interval: &Interval, n: usize, binomials: &[Vec<BigInt>]) -> Self {
        // Positive by construction, independently of how the ends are
        // spelled: an lcm accumulated from 1 over positive denominators.
        let d = lcm_with_denominator(&BigInt::one(), interval.lo());
        let d = lcm_with_denominator(&d, interval.hi());
        let l = scaled_integer(interval.lo(), &d);
        let h = scaled_integer(interval.hi(), &d);
        let m = &h - &l;
        let powers = |x: &BigInt| -> Vec<BigInt> {
            let mut out = vec![BigInt::one()];
            for k in 1..=n {
                let next = &out[k - 1] * x;
                out.push(next);
            }
            out
        };
        let (lp, mp, dp) = (powers(&l), powers(&m), powers(&d));
        // substitution[k][i] = C(i,k) l^(i-k) m^k d^(n-i) for k <= i.
        let substitution = |k: usize, i: usize| -> BigInt {
            &binomials[i][k] * &lp[i - k] * &mp[k] * &dp[n - i]
        };
        let mut matrix = vec![BigInt::zero(); (n + 1) * (n + 1)];
        for r in 0..=n {
            for i in 0..=n {
                let mut entry = BigInt::zero();
                for k in 0..=r.min(i) {
                    entry += &binomials[n - k][r - k] * substitution(k, i);
                }
                matrix[r * (n + 1) + i] = entry;
            }
        }
        Self {
            degree: n,
            matrix,
            scale: dp[n].clone(),
        }
    }

    /// Replace every offset along `axis` by its image under the matrix.
    fn apply(&self, values: &mut [Pair], layout: &Layout, axis: usize) {
        let n = self.degree;
        let stride = layout.strides[axis];
        let block = stride * (n + 1);
        let mut offset: Vec<Pair> = vec![Pair::default(); n + 1];
        for start in (0..layout.size).step_by(block) {
            for base in start..start + stride {
                for (i, slot) in offset.iter_mut().enumerate() {
                    *slot = std::mem::take(&mut values[base + i * stride]);
                }
                for r in 0..=n {
                    let row = &self.matrix[r * (n + 1)..(r + 1) * (n + 1)];
                    let out = &mut values[base + r * stride];
                    for (entry, input) in row.iter().zip(&offset) {
                        if entry.is_zero() {
                            continue;
                        }
                        if !input.a.is_zero() {
                            out.a += entry * &input.a;
                        }
                        if !input.b.is_zero() {
                            out.b += entry * &input.b;
                        }
                    }
                }
            }
        }
    }
}
