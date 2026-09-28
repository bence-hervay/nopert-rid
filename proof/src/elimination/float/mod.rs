//! Floating-point proposals: exact numbers rounded to `f64`, polynomials as
//! dense tensors with their Bernstein coefficients on boxes, and the float
//! geometry of views, rotations and shadow outlines. Nothing computed here
//! decides anything.
use crate::arithmetic::exact::{Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::{Polynomial, VARIABLES};
use crate::elimination::zoom::Five;
use crate::problem::geometry;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use std::sync::OnceLock;

/// A float vector of R³.
pub type Vector = [f64; 3];

/// A nearby `f64` of a rational: the quotient of the rounded numerator and
/// denominator when both are finite, otherwise num-rational's scaled
/// conversion (which handles huge numerators and denominators; NaN only if
/// the value itself is out of range).
pub fn rational(x: &Q) -> f64 {
    match (x.numer().to_f64(), x.denom().to_f64()) {
        (Some(n), Some(d)) if n.is_finite() && d.is_finite() => n / d,
        _ => BigRational::new(x.numer().clone(), x.denom().clone()).to_f64().unwrap_or(f64::NAN),
    }
}

/// `a + b√5` rounded to `f64`.
pub fn number(x: &QSqrt5) -> f64 {
    rational(x.rational_part()) + rational(x.sqrt5_part()) * 5f64.sqrt()
}

pub fn point(x: &geometry::Point) -> Vector {
    std::array::from_fn(|j| number(&x[j]))
}

/// The ends of every interval of a box.
pub fn cell(c: &Five<Interval>) -> Five<[f64; 2]> {
    std::array::from_fn(|j| [rational(c[j].lo()), rational(c[j].hi())])
}

/// A polynomial in the five variables as a dense tensor of `f64`
/// coefficients, indexed by exponents (variable 0 slowest).
///
/// **Invariant:** `values.len()` is the product of `shape`, and `shape[j]` is
/// one more than the degree in variable `j`.
#[derive(Clone, Debug, PartialEq)]
pub struct Dense {
    shape: [usize; VARIABLES],
    values: Vec<f64>,
}

fn strides(shape: &[usize; VARIABLES]) -> [usize; VARIABLES] {
    let mut out = [1; VARIABLES];
    for j in (0..VARIABLES - 1).rev() {
        out[j] = out[j + 1] * shape[j + 1];
    }
    out
}

/// `C(n, k)` as a float (small arguments only).
fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |c, i| c * (n - i) as f64 / (i + 1) as f64)
}

impl Dense {
    pub fn new(p: &Polynomial) -> Self {
        let shape = p.degrees().map(|d| usize::from(d) + 1);
        let stride = strides(&shape);
        let mut values = vec![0.0; shape.iter().product()];
        for (exponents, coefficient) in p.terms() {
            let at: usize = (0..VARIABLES).map(|j| usize::from(exponents[j]) * stride[j]).sum();
            values[at] += number(coefficient);
        }
        Self { shape, values }
    }

    pub fn shape(&self) -> [usize; VARIABLES] {
        self.shape
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The value at a point, as `Σ c_e Π y_j^{e_j}`.
    pub fn evaluate(&self, y: &[f64; VARIABLES]) -> f64 {
        let powers: Vec<Vec<f64>> = (0..VARIABLES)
            .map(|j| (0..self.shape[j]).scan(1.0, |p, _| { let v = *p; *p *= y[j]; Some(v) }).collect())
            .collect();
        let stride = strides(&self.shape);
        self.values
            .iter()
            .enumerate()
            .map(|(at, c)| (0..VARIABLES).fold(*c, |v, j| v * powers[j][at / stride[j] % self.shape[j]]))
            .sum()
    }

    /// The tensor Bernstein coefficients on the box `cell` (lower end first
    /// in each pair), in the same shape: along each variable `x = l + w y`
    /// gives power coefficients `c'_k = Σ_{i≥k} C(i,k) l^{i−k} w^k c_i`, and
    /// `b_r = Σ_{k≤r} C(r,k)/C(n,k) c'_k`.
    pub fn bernstein(&self, cell: &Five<[f64; 2]>) -> Dense {
        let stride = strides(&self.shape);
        let mut values = self.values.clone();
        for axis in 0..VARIABLES {
            let n = self.shape[axis] - 1;
            if n == 0 {
                continue;
            }
            let [lo, hi] = cell[axis];
            let width = hi - lo;
            let matrix: Vec<Vec<f64>> = (0..=n)
                .map(|r| {
                    (0..=n)
                        .map(|i| {
                            (0..=r.min(i))
                                .map(|k| {
                                    binomial(r, k) / binomial(n, k)
                                        * binomial(i, k)
                                        * lo.powi((i - k) as i32)
                                        * width.powi(k as i32)
                                })
                                .sum()
                        })
                        .collect()
                })
                .collect();
            let block = stride[axis] * self.shape[axis];
            let mut offset = vec![0.0; n + 1];
            for outer in (0..values.len()).step_by(block) {
                for inner in 0..stride[axis] {
                    let base = outer + inner;
                    for (i, f) in offset.iter_mut().enumerate() {
                        *f = values[base + i * stride[axis]];
                    }
                    for (r, row) in matrix.iter().enumerate() {
                        values[base + r * stride[axis]] = row.iter().zip(&offset).map(|(m, f)| m * f).sum();
                    }
                }
            }
        }
        Dense { shape: self.shape, values }
    }

    /// The largest entry.
    pub fn max(&self) -> f64 {
        self.values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// The largest absolute value of an entry.
    pub fn magnitude(&self) -> f64 {
        self.values.iter().fold(0.0, |m, v| m.max(v.abs()))
    }

    /// The largest absolute difference of neighbouring entries along `axis`
    /// (0 when the degree in `axis` is 0).
    pub fn variation(&self, axis: usize) -> f64 {
        let stride = strides(&self.shape)[axis];
        (0..self.values.len())
            .filter(|at| at / stride % self.shape[axis] + 1 < self.shape[axis])
            .fold(0.0, |m, at| m.max((self.values[at + stride] - self.values[at]).abs()))
    }
}

/// The 60 vertices as float points.
pub fn vertices() -> &'static [Vector] {
    static VERTICES: OnceLock<Vec<Vector>> = OnceLock::new();
    VERTICES.get_or_init(|| geometry::vertices().iter().map(point).collect())
}

pub fn dot(a: &Vector, b: &Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: &Vector, b: &Vector) -> Vector {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// The rotation `R(r) = R̂(r)/(1 + |r|²)` of the quaternion `(1, r)`, with
/// `R̂ = (1 − |r|²) I + 2 r rᵀ + 2 [r]×`.
pub fn rotation(r: &Vector) -> [Vector; 3] {
    let n = dot(r, r);
    let skew = [[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]];
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let diagonal = if i == j { 1.0 - n } else { 0.0 };
            (diagonal + 2.0 * r[i] * r[j] + 2.0 * skew[i][j]) / (1.0 + n)
        })
    })
}

pub fn apply(m: &[Vector; 3], x: &Vector) -> Vector {
    std::array::from_fn(|i| dot(&m[i], x))
}

/// The screen map `T_u(X) = (u₃X₁ − u₁X₃, u₃X₂ − u₂X₃)`.
pub fn screen(u: &Vector, x: &Vector) -> [f64; 2] {
    [u[2] * x[0] - u[0] * x[2], u[2] * x[1] - u[1] * x[2]]
}

/// The vertices of the convex hull of `points` in counter-clockwise order,
/// starting from the lowest-leftmost, without collinear points (Andrew's
/// monotone chain). At least three points in general position are assumed.
pub fn hull(points: &[[f64; 2]]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&a, &b| points[a][0].total_cmp(&points[b][0]).then(points[a][1].total_cmp(&points[b][1])));
    let turn = |o: usize, a: usize, b: usize| {
        let (p, q, r) = (points[o], points[a], points[b]);
        (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    };
    let mut chain: Vec<usize> = Vec::new();
    for pass in 0..2 {
        let start = chain.len();
        for &i in order.iter() {
            while chain.len() >= start + 2 && turn(chain[chain.len() - 2], chain[chain.len() - 1], i) <= 0.0 {
                chain.pop();
            }
            chain.push(i);
        }
        chain.pop();
        if pass == 0 {
            order.reverse();
        }
    }
    chain
}

#[cfg(test)]
mod tests;
