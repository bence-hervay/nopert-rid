//! The shadows of the hole and the plug in the orthonormal screen: exact
//! positions and silhouettes at one configuration, and outer envelopes of
//! every vertex and edge over a box of configurations.
use super::{Point, Shape};
use num_bigint::BigInt;
use rid::arithmetic::exact::{q, Interval, QSqrt5, Q};
use rid::arithmetic::polynomial::{Polynomial, PolynomialError};
use rid::problem::configuration::{coordinates, ConfigurationBox, AXES};
use rid::problem::geometry::{cayley_matrix, edges, polynomial_matrix_point, squared_norm, vertices, Point as SpacePoint};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::ops::{Mul, Sub};
use std::sync::OnceLock;

/// Bits of the rational brackets of √5 and of square roots.
pub const RADICAL_BITS: u32 = 160;
/// Every drawn position of a configuration's vertex is within `2^-POINT_BITS`
/// of the exact screen position in each coordinate; this is checked.
pub const POINT_BITS: u32 = 100;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Body {
    /// The RID itself.
    Hole,
    /// The RID rotated by the configuration's rotation.
    Plug,
}

impl Body {
    fn index(self) -> usize {
        match self {
            Body::Hole => 0,
            Body::Plug => 1,
        }
    }
}

/// A closed axis-parallel rectangle `[x] × [y]` of the screen.
pub type Rectangle = [Interval; 2];

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SceneError {
    Polynomial(PolynomialError),
    /// The rational enclosure of a configuration's vertex is wider than
    /// `2^-POINT_BITS`: the configuration's coordinates are too large.
    PrecisionLost { body: Body, vertex: usize },
}

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SceneError::Polynomial(e) => write!(f, "screen polynomial: {e}"),
            SceneError::PrecisionLost { body, vertex } => write!(
                f,
                "the {body:?} vertex {vertex} is not located to 2^-{POINT_BITS}; \
                 the configuration's coordinates are too large"
            ),
        }
    }
}

impl std::error::Error for SceneError {}

/// The screen map as polynomials in the configuration `(s, t, r)`.
///
/// With `L = 1 + s²` and `N = 1 + s² + t²`, the orthonormal screen of the
/// view `u = (s, t, 1)` has the rows `e₁ = (1, 0, -s)/√L` and
/// `e₂ = (-st, L, -t)/√(LN)`. A point `w` of a body has the screen position
/// `(P₁(w)/√L, P₂(w)/√(LN))` with `P₁ = w₁ - s·w₃` and
/// `P₂ = L·w₂ - st·w₁ - t·w₃`. The plug's vertex `p` is `w = R̂(r)p / D` with
/// `D = 1 + |r|²`; the hole's is `w = p` with `D = 1`.
struct ScreenPolynomials {
    /// `[body][vertex] = [P₁(R̂p), P₂(R̂p)]` (with `R̂ = I` for the hole).
    numerators: [Vec<[Polynomial; 2]>; 2],
    /// `[body] = D`.
    denominators: [Polynomial; 2],
    /// `[L, LN]`.
    normalisation: [Polynomial; 2],
}

fn build() -> Result<ScreenPolynomials, PolynomialError> {
    let (view, r) = coordinates();
    let [s, t, _] = view.vector();
    let one = Polynomial::constant(QSqrt5::one());
    let l = &one + &s.mul(&s)?;
    let n = &l + &t.mul(&t)?;
    let st = s.mul(&t)?;
    let numerator = |w: [Polynomial; 3]| -> Result<[Polynomial; 2], PolynomialError> {
        let first = &w[0] - &s.mul(&w[2])?;
        let second = &(&l.mul(&w[1])? - &st.mul(&w[0])?) - &t.mul(&w[2])?;
        Ok([first, second])
    };
    let rotation = cayley_matrix(&r)?;
    let hole = vertices()
        .iter()
        .map(|p| numerator(p.clone().map(Polynomial::constant)))
        .collect::<Result<_, _>>()?;
    let plug = vertices()
        .iter()
        .map(|p| numerator(polynomial_matrix_point(&rotation, p)))
        .collect::<Result<_, _>>()?;
    Ok(ScreenPolynomials {
        numerators: [hole, plug],
        denominators: [one.clone(), &one + &squared_norm(&r)?],
        normalisation: [l.clone(), l.mul(&n)?],
    })
}

fn screen_polynomials() -> Result<&'static ScreenPolynomials, SceneError> {
    static POLYNOMIALS: OnceLock<Result<ScreenPolynomials, PolynomialError>> = OnceLock::new();
    POLYNOMIALS
        .get_or_init(build)
        .as_ref()
        .map_err(|e| SceneError::Polynomial(e.clone()))
}

/// The screen rectangles `P_j · 1/(D·√norm_j)` from enclosures of the
/// numerators, the denominators and the normalisations. `D`, `L` and `LN`
/// are at least 1 at every real configuration, so their lower bounds are
/// raised to 1: this keeps the enclosures sound and the divisions defined.
fn rectangles(
    numerators: [Vec<[Interval; 2]>; 2],
    denominators: [Interval; 2],
    normalisation: [Interval; 2],
) -> [Vec<Rectangle>; 2] {
    let roots = normalisation.map(|n| square_root(&at_least_one(n)));
    let [hole, plug] = numerators;
    let [hole_denominator, plug_denominator] = denominators;
    [(hole, hole_denominator), (plug, plug_denominator)].map(|(numerators, denominator)| {
        let d = at_least_one(denominator);
        let scale = [0, 1].map(|j| reciprocal(&(&d * &roots[j])));
        numerators
            .iter()
            .map(|p| [&p[0] * &scale[0], &p[1] * &scale[1]])
            .collect()
    })
}

/// The shadows at one configuration `(s, t, r₁, r₂, r₃)` of Q(√5)⁵.
#[derive(Clone, Debug)]
pub struct ConfigurationScene {
    /// Rational enclosures of the exact screen positions, `[body][vertex]`.
    rectangles: [Vec<Rectangle>; 2],
    /// The drawn positions: the enclosures' midpoints.
    positions: [Vec<Point>; 2],
    /// Hull vertex indices, decided exactly.
    silhouettes: [Vec<usize>; 2],
    /// `[body][edge]`: whether both faces at the edge face away from the
    /// viewer, so that the body itself hides it; decided exactly.
    back: [Vec<bool>; 2],
}

impl ConfigurationScene {
    pub fn new(configuration: &[QSqrt5; AXES]) -> Result<Self, SceneError> {
        let polynomials = screen_polynomials()?;
        let exact: [Vec<[QSqrt5; 2]>; 2] = std::array::from_fn(|body| {
            polynomials.numerators[body]
                .iter()
                .map(|p| [p[0].evaluate(configuration), p[1].evaluate(configuration)])
                .collect()
        });
        let value = |p: &Polynomial| enclose(&p.evaluate(configuration));
        let rectangles = rectangles(
            std::array::from_fn(|body| {
                exact[body].iter().map(|x| [enclose(&x[0]), enclose(&x[1])]).collect()
            }),
            std::array::from_fn(|body| value(&polynomials.denominators[body])),
            std::array::from_fn(|j| value(&polynomials.normalisation[j])),
        );
        let limit = Q::new(BigInt::from(1), BigInt::from(1) << POINT_BITS);
        for (body, name) in [Body::Hole, Body::Plug].into_iter().enumerate() {
            if let Some(vertex) = rectangles[body]
                .iter()
                .position(|r| r[0].width() > limit || r[1].width() > limit)
            {
                return Err(SceneError::PrecisionLost { body: name, vertex });
            }
        }
        // A positive diagonal scaling of the numerators gives the screen
        // positions, and it preserves orientations: the hull of the exact
        // numerators is the silhouette.
        let silhouettes = std::array::from_fn(|body| hull(&exact[body]));
        let positions = std::array::from_fn(|body| {
            rectangles[body].iter().map(|r| [r[0].midpoint(), r[1].midpoint()]).collect()
        });
        let back = back_edges(configuration)?;
        Ok(Self {
            rectangles,
            positions,
            silhouettes,
            back,
        })
    }

    /// A rational rectangle containing the exact screen position of a vertex.
    pub fn enclosure(&self, body: Body, vertex: usize) -> &Rectangle {
        &self.rectangles[body.index()][vertex]
    }

    /// The drawn screen position of a vertex (within `2^-POINT_BITS`).
    pub fn position(&self, body: Body, vertex: usize) -> &Point {
        &self.positions[body.index()][vertex]
    }

    /// The silhouette's vertices, counterclockwise.
    pub fn silhouette(&self, body: Body) -> &[usize] {
        &self.silhouettes[body.index()]
    }

    /// All 120 edges of the body as segments.
    pub fn wireframe(&self, body: Body) -> Vec<Shape> {
        let positions = &self.positions[body.index()];
        edges()
            .iter()
            .map(|&[a, b]| Shape::Segment([positions[a].clone(), positions[b].clone()]))
            .collect()
    }

    /// The edges on the side of the body facing the viewer (or on its
    /// silhouette), as segments.
    pub fn front_wireframe(&self, body: Body) -> Vec<Shape> {
        self.wireframe_where(body, false)
    }

    /// The edges on the far side of the body, hidden by the body itself.
    pub fn back_wireframe(&self, body: Body) -> Vec<Shape> {
        self.wireframe_where(body, true)
    }

    fn wireframe_where(&self, body: Body, back: bool) -> Vec<Shape> {
        let positions = &self.positions[body.index()];
        edges()
            .iter()
            .zip(&self.back[body.index()])
            .filter(|(_, &b)| b == back)
            .map(|(&[a, b], _)| Shape::Segment([positions[a].clone(), positions[b].clone()]))
            .collect()
    }

    /// The boundary of the body's shadow as one polygon.
    pub fn outline(&self, body: Body) -> Shape {
        let positions = &self.positions[body.index()];
        Shape::Polygon(self.silhouette(body).iter().map(|&j| positions[j].clone()).collect())
    }
}

/// The outward normals (not normalised) of the two faces at each edge of
/// the RID, found exactly: a plane through the edge and a third vertex with
/// every vertex on its inner side.
fn edge_normals() -> &'static [[SpacePoint; 2]] {
    static NORMALS: OnceLock<Vec<[SpacePoint; 2]>> = OnceLock::new();
    NORMALS.get_or_init(|| {
        let v = vertices();
        let minus = |a: &SpacePoint, b: &SpacePoint| -> SpacePoint { std::array::from_fn(|i| &a[i] - &b[i]) };
        let dot = |a: &SpacePoint, b: &SpacePoint| -> QSqrt5 { &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2]) };
        let cross = |a: &SpacePoint, b: &SpacePoint| -> SpacePoint {
            [
                &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
                &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
                &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
            ]
        };
        let parallel = |a: &SpacePoint, b: &SpacePoint| cross(a, b).iter().all(|x| x.sign() == std::cmp::Ordering::Equal);
        edges()
            .iter()
            .map(|&[a, b]| {
                let along = minus(&v[b], &v[a]);
                let mut found: Vec<SpacePoint> = Vec::new();
                for (c, w) in v.iter().enumerate() {
                    if c == a || c == b || found.len() == 2 {
                        continue;
                    }
                    let mut n = cross(&along, &minus(w, &v[a]));
                    if n.iter().all(|x| x.sign() == std::cmp::Ordering::Equal) {
                        continue;
                    }
                    // Outward: the centre, the origin, is on the inner side.
                    if dot(&n, &v[a]).sign() == std::cmp::Ordering::Less {
                        n = n.map(|x| -&x);
                    }
                    let supporting = v.iter().all(|u| dot(&n, &minus(u, &v[a])).sign() != std::cmp::Ordering::Greater);
                    if supporting && !found.iter().any(|m| parallel(m, &n)) {
                        found.push(n);
                    }
                }
                let [first, second]: [SpacePoint; 2] = found.try_into().expect("two faces at every edge");
                [first, second]
            })
            .collect()
    })
}

/// `[body][edge]`: whether both faces at the edge face away from the viewer
/// at `u = (s, t, 1)`, which the screen map orients towards the viewer. The
/// plug's normals are turned by the rotation's numerator `R̂(r)`, a positive
/// multiple of the rotation.
fn back_edges(configuration: &[QSqrt5; AXES]) -> Result<[Vec<bool>; 2], SceneError> {
    let (_, r) = coordinates();
    let rotation = cayley_matrix(&r).map_err(SceneError::Polynomial)?;
    let view: SpacePoint = [configuration[0].clone(), configuration[1].clone(), QSqrt5::one()];
    let facing = |n: &SpacePoint| -> QSqrt5 { &(&(&n[0] * &view[0]) + &(&n[1] * &view[1])) + &(&n[2] * &view[2]) };
    let away = |n: &SpacePoint| facing(n).sign() == std::cmp::Ordering::Less;
    let turned = |n: &SpacePoint| -> SpacePoint {
        polynomial_matrix_point(&rotation, n).map(|p| p.evaluate(configuration))
    };
    let hole = edge_normals().iter().map(|[m, n]| away(m) && away(n)).collect();
    let plug = edge_normals().iter().map(|[m, n]| away(&turned(m)) && away(&turned(n))).collect();
    Ok([hole, plug])
}

/// Outer envelopes of the shadows over every configuration of a closed box.
#[derive(Clone, Debug)]
pub struct BoxScene {
    /// `[body][vertex]`: every screen position of the vertex over the box.
    rectangles: [Vec<Rectangle>; 2],
}

impl BoxScene {
    pub fn new(b: &ConfigurationBox) -> Result<Self, SceneError> {
        let polynomials = screen_polynomials()?;
        let cell = b.axes();
        let range = |p: &Polynomial| range_on(p, cell).map_err(SceneError::Polynomial);
        let mut numerators = [Vec::new(), Vec::new()];
        for (body, target) in numerators.iter_mut().enumerate() {
            for p in &polynomials.numerators[body] {
                target.push([range(&p[0])?, range(&p[1])?]);
            }
        }
        let [d0, d1] = &polynomials.denominators;
        let [l, ln] = &polynomials.normalisation;
        Ok(Self {
            rectangles: rectangles(numerators, [range(d0)?, range(d1)?], [range(l)?, range(ln)?]),
        })
    }

    /// Every screen position of the vertex over the box lies in this rectangle.
    pub fn rectangle(&self, body: Body, vertex: usize) -> &Rectangle {
        &self.rectangles[body.index()][vertex]
    }

    /// Every screen position of the edge `[a, b]` over the box lies in this
    /// convex polygon, the hull of the two endpoints' rectangles: at each
    /// configuration the edge is the segment between its endpoints, and both
    /// lie in the hull.
    pub fn edge_envelope(&self, body: Body, edge: [usize; 2]) -> Vec<Point> {
        let corners: Vec<Point> = edge
            .iter()
            .flat_map(|&v| corners(self.rectangle(body, v)))
            .collect();
        hull(&corners).into_iter().map(|j| corners[j].clone()).collect()
    }

    pub fn vertex_envelopes(&self, body: Body) -> Vec<Shape> {
        self.rectangles[body.index()]
            .iter()
            .map(|r| Shape::Polygon(corners(r).to_vec()))
            .collect()
    }

    pub fn edge_envelopes(&self, body: Body) -> Vec<Shape> {
        edges()
            .iter()
            .map(|&edge| Shape::Polygon(self.edge_envelope(body, edge)))
            .collect()
    }
}

/// The corners of a rectangle, counterclockwise from the lower left.
fn corners(r: &Rectangle) -> [Point; 4] {
    let [x, y] = r;
    [
        [x.lo().clone(), y.lo().clone()],
        [x.hi().clone(), y.lo().clone()],
        [x.hi().clone(), y.hi().clone()],
        [x.lo().clone(), y.hi().clone()],
    ]
}

/// The vertices of the convex hull of `points`, as indices, counterclockwise
/// from the smallest point in lexicographic order. Points inside the hull or
/// on an edge between two hull vertices are left out, and of equal points
/// only the lowest index is kept, so a point or a segment gives one or two
/// indices. Orientation is decided exactly.
pub fn hull<T>(points: &[[T; 2]]) -> Vec<usize>
where
    T: Ord,
    for<'a> &'a T: Sub<Output = T> + Mul<Output = T>,
{
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&i, &j| points[i].cmp(&points[j]).then(i.cmp(&j)));
    order.dedup_by(|later, earlier| points[*later] == points[*earlier]);
    if order.len() <= 2 {
        return order;
    }
    // Whether c lies strictly to the left of the directed line a → b.
    let left = |a: usize, b: usize, c: usize| {
        let [a, b, c] = [&points[a], &points[b], &points[c]];
        let along = &(&b[0] - &a[0]) * &(&c[1] - &a[1]);
        let across = &(&b[1] - &a[1]) * &(&c[0] - &a[0]);
        along > across
    };
    // Andrew's monotone chain: the lower hull left to right, then the upper
    // hull right to left, each keeping only strict left turns.
    let mut chain: Vec<usize> = Vec::with_capacity(order.len() + 1);
    for &i in &order {
        while chain.len() >= 2 && !left(chain[chain.len() - 2], chain[chain.len() - 1], i) {
            chain.pop();
        }
        chain.push(i);
    }
    let lower = chain.len();
    for &i in order.iter().rev().skip(1) {
        while chain.len() > lower && !left(chain[chain.len() - 2], chain[chain.len() - 1], i) {
            chain.pop();
        }
        chain.push(i);
    }
    chain.pop(); // the upper hull ends where the lower hull began
    chain
}

// ---- Rational enclosures ----------------------------------------------------

fn scale() -> BigInt {
    BigInt::from(1) << RADICAL_BITS
}

/// `[k, k+1]/2^B` with `k = ⌊√5·2^B⌋`; `√5` is irrational, so it lies strictly inside.
fn sqrt5() -> &'static Interval {
    static SQRT5: OnceLock<Interval> = OnceLock::new();
    SQRT5.get_or_init(|| {
        let k = (BigInt::from(5) << (2 * RADICAL_BITS)).sqrt();
        Interval::new(Q::new(k.clone(), scale()), Q::new(k + 1, scale())).expect("k < k + 1")
    })
}

/// A rational interval containing `a + b√5`.
fn enclose(x: &QSqrt5) -> Interval {
    let (a, b) = (x.rational_part(), x.sqrt5_part());
    let root = sqrt5();
    let (lo, hi) = if *b >= q(0) {
        (root.lo(), root.hi())
    } else {
        (root.hi(), root.lo())
    };
    Interval::new(a + b * lo, a + b * hi).expect("ordered by the sign of b")
}

/// The range of `p` on the closed box is inside the hull of its Bernstein
/// coefficients (each value is their convex combination).
fn range_on(p: &Polynomial, cell: &[Interval; AXES]) -> Result<Interval, PolynomialError> {
    let coefficients = p.bernstein_coefficients(cell)?;
    let mut values = coefficients.iter().map(|(_, value)| value);
    let first = values.next().expect("at least one Bernstein coefficient");
    let (lo, hi) = values.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v)));
    Ok(Interval::new(enclose(lo).lo().clone(), enclose(hi).hi().clone())
        .expect("the enclosure of the minimum starts below that of the maximum"))
}

/// `[max(lo, 1), hi]` for an enclosure of a quantity that is at least 1.
fn at_least_one(x: Interval) -> Interval {
    let one = q(1);
    if *x.lo() >= one {
        return x;
    }
    Interval::new(one, x.hi().clone()).expect("an enclosure of a value ≥ 1 ends at or above 1")
}

/// A rational interval containing `√x` for every `x` of a nonnegative interval.
fn square_root(x: &Interval) -> Interval {
    assert!(*x.lo() >= q(0), "square root of a negative interval");
    // ⌊√(n/d)·2^B⌋ = ⌊√⌊n·4^B/d⌋⌋ and ⌈√(n/d)·2^B⌉ = ⌈√⌈n·4^B/d⌉⌉.
    let scaled = |v: &Q| -> BigInt { v.numer() << (2 * RADICAL_BITS) };
    let floor: BigInt = (scaled(x.lo()) / x.lo().denom()).sqrt();
    let ceiling_square: BigInt = (scaled(x.hi()) + x.hi().denom() - 1) / x.hi().denom();
    let mut ceiling = ceiling_square.sqrt();
    if &ceiling * &ceiling < ceiling_square {
        ceiling += 1;
    }
    Interval::new(Q::new(floor, scale()), Q::new(ceiling, scale())).expect("floor ≤ ceiling")
}

/// `[1/hi, 1/lo]` for a positive interval.
fn reciprocal(x: &Interval) -> Interval {
    assert!(*x.lo() > q(0), "reciprocal of a non-positive interval");
    Interval::new(q(1) / x.hi(), q(1) / x.lo()).expect("1/hi ≤ 1/lo")
}

#[cfg(test)]
mod tests;

