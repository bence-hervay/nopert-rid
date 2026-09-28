//! The graph of a polynomial in one variable on a closed interval, and its
//! Bernstein coefficients there, computed by the crate's own exact
//! conversion.
use super::{Point, Shape};
use rid::arithmetic::exact::{q, Interval, QSqrt5, Q};
use rid::arithmetic::polynomial::{Polynomial, VARIABLES};
use std::fmt;

/// A polynomial `Σ cₖ xᵏ` on `[l, h]`, drawn in the plane of `(x, p(x))`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    /// The graph through `samples + 1` equally spaced exact points.
    pub curve: Vec<Shape>,
    /// The control polygon: the Bernstein coefficient `b_k` at
    /// `x = l + k(h − l)/n`, joined in order.
    pub polygon: Vec<Shape>,
    /// The control points of the polygon.
    pub points: Vec<Point>,
    /// The segment of the horizontal axis over `[l, h]`.
    pub axis: Vec<Shape>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphError {
    /// No coefficients, or a degree the crate's polynomials cannot hold.
    Degree,
    /// `h ≤ l`.
    Interval,
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GraphError::Degree => write!(f, "a polynomial needs between 1 and 255 coefficients"),
            GraphError::Interval => write!(f, "the interval must have positive width"),
        }
    }
}

impl std::error::Error for GraphError {}

fn value(coefficients: &[Q], x: &Q) -> Q {
    coefficients.iter().rev().fold(Q::zero(), |sum, c| &(&sum * x) + c)
}

/// The graph of `Σ cₖ xᵏ` (`coefficients[k] = cₖ`) on `[l, h]`.
pub fn graph(coefficients: &[Q], interval: [Q; 2], samples: usize) -> Result<Graph, GraphError> {
    if coefficients.is_empty() || coefficients.len() > 255 {
        return Err(GraphError::Degree);
    }
    let [l, h] = interval;
    if h <= l {
        return Err(GraphError::Interval);
    }
    let width = &h - &l;
    let samples = samples.max(1);
    let at = |i: usize, of: usize| &l + &(&width * &Q::new((i as i64).into(), (of as i64).into()));
    let curve_points: Vec<Point> = (0..=samples)
        .map(|i| {
            let x = at(i, samples);
            let y = value(coefficients, &x);
            [x, y]
        })
        .collect();
    let curve = curve_points
        .windows(2)
        .map(|pair| Shape::Segment([pair[0].clone(), pair[1].clone()]))
        .collect();
    // The polynomial in the crate's first variable, and its tensor Bernstein
    // coefficients on [l, h] × {0}⁴.
    let mut p = Polynomial::zero();
    for (k, c) in coefficients.iter().enumerate() {
        let mut exponents = [0u8; VARIABLES];
        exponents[0] = u8::try_from(k).expect("at most 255 coefficients");
        p = &p + &Polynomial::monomial(exponents, QSqrt5::from_rational(c.clone()));
    }
    let degree = usize::from(p.degrees()[0]);
    let mut cell: [Interval; VARIABLES] = std::array::from_fn(|_| Interval::point(q(0)));
    cell[0] = Interval::new(l.clone(), h.clone()).expect("l < h");
    let bernstein = p.bernstein_coefficients(&cell).map_err(|_| GraphError::Degree)?;
    let points: Vec<Point> = (0..=degree)
        .map(|k| {
            let mut index = [0u8; VARIABLES];
            index[0] = u8::try_from(k).expect("degree below 256");
            let b = bernstein.get(&index).expect("an index within the degree");
            [at(k, degree.max(1)), b.rational_part().clone()]
        })
        .collect();
    let polygon = points
        .windows(2)
        .map(|pair| Shape::Segment([pair[0].clone(), pair[1].clone()]))
        .collect();
    Ok(Graph {
        curve,
        polygon,
        points,
        axis: vec![Shape::Segment([[l.clone(), Q::zero()], [h.clone(), Q::zero()]])],
    })
}

/// Diamonds of half-diagonal `size` around points.
pub fn markers(points: &[Point], size: &Q) -> Vec<Shape> {
    points
        .iter()
        .map(|p| {
            Shape::Polygon(vec![
                [&p[0] + size, p[1].clone()],
                [p[0].clone(), &p[1] + size],
                [&p[0] - size, p[1].clone()],
                [p[0].clone(), &p[1] - size],
            ])
        })
        .collect()
}

#[cfg(test)]
mod tests;
