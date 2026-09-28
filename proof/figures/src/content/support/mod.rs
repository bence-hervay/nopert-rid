//! A support drawn at one configuration: the supporting line of the hole's
//! shadow with the outward normal `u × (v_b − v_a)` of a hole edge, the part
//! of the hole's silhouette on that line, the normal as an arrow from the
//! line towards a plug vertex, and a marker at that vertex.
use super::scene::{Body, ConfigurationScene};
use super::{Point, Shape};
use num_bigint::BigInt;
use rid::arithmetic::exact::{frac, q, Q};
use rid::problem::geometry::VERTEX_COUNT;
use std::fmt;

/// The four parts of a drawn support, in world (screen) coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Support {
    /// The part of the hole's silhouette on the supporting line: the segment
    /// between the extreme projected hole vertices on it (empty when a single
    /// vertex attains the extent).
    pub edge: Vec<Shape>,
    /// The supporting line, long enough to cross any panel of its scale.
    pub line: Vec<Shape>,
    /// The outward normal: a shaft and two barbs.
    pub normal: Vec<Shape>,
    /// A dot at the projected plug vertex.
    pub vertex: Vec<Shape>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SupportError {
    /// A vertex index of 60 or more.
    Vertex(usize),
    /// The edge's two ends are the same vertex.
    SameEnds,
    /// The edge projects to a point: seen end-on, it has no normal.
    EndOn,
    /// A length that is not positive.
    Length,
}

impl fmt::Display for SupportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SupportError::Vertex(v) => write!(f, "vertex {v} is not below {VERTEX_COUNT}"),
            SupportError::SameEnds => write!(f, "the edge's ends are the same vertex"),
            SupportError::EndOn => write!(f, "the edge projects to a point"),
            SupportError::Length => write!(f, "the normal length and the marker must be positive"),
        }
    }
}

impl std::error::Error for SupportError {}

/// Bits of the rational square roots used to normalise directions; far
/// below any visible scale.
const ROOT_BITS: u32 = 64;

/// A rational within `2^-ROOT_BITS` of `√x` for `x ≥ 0`.
fn root(x: &Q) -> Q {
    let scaled: BigInt = (x.numer() << (2 * ROOT_BITS)) / x.denom();
    Q::new(scaled.sqrt(), BigInt::from(1) << ROOT_BITS)
}

fn dot(a: &Point, b: &Point) -> Q {
    &(&a[0] * &b[0]) + &(&a[1] * &b[1])
}

fn along(p: &Point, d: &Point, k: &Q) -> Point {
    [&p[0] + &(k * &d[0]), &p[1] + &(k * &d[1])]
}

/// The support of the ordered hole edge `edge = [a, b]` with the plug vertex
/// `plug`, at the configuration of `scene`. The normal `n` is the edge's
/// screen direction turned by a quarter turn counterclockwise, which is the
/// screen image of `u × (v_b − v_a)` in the orthonormal screen; the line is
/// `{x : n·x = max_w n·w}` over the hole's vertices; the silhouette segment
/// joins the extreme hole vertices whose `n·w` is within `2^-80·|n|` of the
/// extent; the arrow of length `length` starts at the foot of the
/// perpendicular from the plug vertex to the line; the dot has radius
/// `marker`.
pub fn support(
    scene: &ConfigurationScene,
    edge: [usize; 2],
    plug: usize,
    length: &Q,
    marker: &Q,
) -> Result<Support, SupportError> {
    for v in [edge[0], edge[1], plug] {
        if v >= VERTEX_COUNT {
            return Err(SupportError::Vertex(v));
        }
    }
    if edge[0] == edge[1] {
        return Err(SupportError::SameEnds);
    }
    if *length <= q(0) || *marker <= q(0) {
        return Err(SupportError::Length);
    }
    let (a, b) = (scene.position(Body::Hole, edge[0]), scene.position(Body::Hole, edge[1]));
    let d: Point = [&b[0] - &a[0], &b[1] - &a[1]];
    let n: Point = [-&d[1], d[0].clone()];
    let squared = dot(&n, &n);
    if squared.is_zero() {
        return Err(SupportError::EndOn);
    }
    let extent = (0..VERTEX_COUNT)
        .map(|w| dot(&n, scene.position(Body::Hole, w)))
        .max()
        .expect("sixty vertices");
    let p = scene.position(Body::Plug, plug);
    let base = along(p, &n, &(&(&extent - &dot(&n, p)) / &squared));
    let norm = root(&squared);
    // The hole vertices on the line, up to the drawn positions' precision,
    // ordered along it.
    let tolerance = &norm * &Q::new(BigInt::from(1), BigInt::from(1) << 80);
    let mut on_line: Vec<(Q, &Point)> = (0..VERTEX_COUNT)
        .map(|w| scene.position(Body::Hole, w))
        .filter(|w| &extent - &dot(&n, w) <= tolerance)
        .map(|w| (dot(&d, w), w))
        .collect();
    on_line.sort_by(|x, y| x.0.cmp(&y.0));
    let edge_shapes = match (on_line.first(), on_line.last()) {
        (Some(first), Some(last)) if first.0 < last.0 => vec![Shape::Segment([first.1.clone(), last.1.clone()])],
        _ => Vec::new(),
    };
    let unit_n: Point = [&n[0] / &norm, &n[1] / &norm];
    let unit_d: Point = [&d[0] / &norm, &d[1] / &norm];
    let far = length * &q(1000);
    let line = Shape::Segment([along(&base, &unit_d, &-&far), along(&base, &unit_d, &far)]);
    let tip = along(&base, &unit_n, length);
    let back = along(&tip, &unit_n, &-(length * &frac(3, 10)));
    let barb = length * &frac(3, 20);
    let normal = vec![
        Shape::Segment([base.clone(), tip.clone()]),
        Shape::Segment([tip.clone(), along(&back, &unit_d, &barb)]),
        Shape::Segment([tip, along(&back, &unit_d, &-&barb)]),
    ];
    let dot = Shape::Circle { centre: p.clone(), radius: marker.clone() };
    Ok(Support {
        edge: edge_shapes,
        line: vec![line],
        normal,
        vertex: vec![dot],
    })
}

#[cfg(test)]
mod tests;
