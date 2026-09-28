//! Candidate witnesses for a zoom: which supports and plug vertices to try,
//! found in floating point at sample configurations and at cells, built
//! exactly, pulled back along the zoom and divided by their factor
//! monomial. A candidate is only a proposal.
use crate::elimination::float::{self, Dense, Vector};
use crate::arithmetic::exact::{frac, QSqrt5};
use crate::arithmetic::polynomial::{Exponents, PolynomialError, VARIABLES};
use crate::elimination::zoom::cover::{ZoomCover, Kind};
use crate::elimination::zoom::{Five, Zoom, ZoomFactor};
use crate::elimination::witness::{Direction, DomainInequality, Edge, Gap, Support, Witness, WitnessError};
use crate::problem::geometry::{self, Point, VERTEX_COUNT};
use std::collections::BTreeSet;
use std::num::NonZeroUsize;

/// Relative tolerance of the float filters: a condition that holds up to
/// this much is kept, so that exact ties (touching vertices, supports at
/// silhouette events) are proposed.
pub const TOLERANCE: f64 = 1e-9;

/// `separating`: a plug vertex closer than this to the hole's outline in the
/// screen (absolute) takes the normal of the outline side it is most beyond
/// instead of the direction to its nearest outline point.
pub const ON_OUTLINE: f64 = 1e-9;

/// `separating`: an edge whose normal is shorter than this (seen nearly
/// end-on) is not proposed.
pub const END_ON: f64 = 1e-12;

/// Denominator of the dyadic grid to which proposed planar directions are
/// rounded.
pub const DIRECTION_GRID: i64 = 4096;

/// A witness pulled back along one zoom: the exact witness and factor
/// that are written into cover files, and float data used only to rank and
/// filter it.
///
/// **Invariant:** `factor` is positive only on scale variables of the
/// zoom, and `quotient` is the float tensor of `g ∘ zoom / y^factor`,
/// which is exact as a polynomial (checked by the division).
#[derive(Clone, Debug)]
pub struct Candidate {
    witness: Witness,
    factor: ZoomFactor,
    quotient: Dense,
    /// `n(u)·(h − w) = Σ_k u_k row_k[w]` for every vertex `w`; empty for a
    /// domain inequality.
    rows: Vec<Vector>,
}

impl Candidate {
    /// The candidate of `witness` on `zoom`: the factor is the largest
    /// monomial in the zoom's scale variables dividing `g ∘ zoom`. `None`
    /// when `g ∘ zoom` vanishes identically.
    pub fn new(zoom: &Zoom, witness: Witness) -> Result<Option<Self>, PolynomialError> {
        let pulled = match witness.polynomial(zoom.view(), zoom.rotation()) {
            Ok(p) => p,
            Err(WitnessError::Polynomial(error)) => return Err(error),
            Err(error) => unreachable!("a constructed witness pulls back: {error}"),
        };
        if pulled.is_zero() {
            return Ok(None);
        }
        let factor: Exponents = std::array::from_fn(|j| {
            if zoom.is_scale(j) {
                pulled.terms().map(|(e, _)| e[j]).min().expect("a nonzero polynomial")
            } else {
                0
            }
        });
        let quotient = Dense::new(&pulled.divide_by_monomial(&factor)?);
        let factor = zoom.factor(factor).expect("exponents only on scale variables");
        let rows = match &witness {
            Witness::Gap(gap) => validity_rows(gap.support()),
            Witness::Domain(_) => Vec::new(),
        };
        Ok(Some(Self { witness, factor, quotient, rows }))
    }

    pub fn witness(&self) -> &Witness {
        &self.witness
    }
    pub fn factor(&self) -> &ZoomFactor {
        &self.factor
    }
    pub fn quotient(&self) -> &Dense {
        &self.quotient
    }
    pub fn is_domain(&self) -> bool {
        matches!(self.witness, Witness::Domain(_))
    }

    /// Whether the support is not clearly invalid at any of the float
    /// `views`: no validity condition is negative beyond the tolerance.
    pub fn plausibly_valid(&self, views: &[Vector]) -> bool {
        views.iter().all(|u| {
            let size = float::dot(u, u).sqrt();
            self.rows
                .iter()
                .all(|row| float::dot(row, u) >= -TOLERANCE * size * float::dot(row, row).sqrt())
        })
    }
}

/// The rows of the validity conditions: `n(u)` is linear in `u`, so
/// `n(u) = Σ_k u_k n(e_k)` and `n(u)·(h − w) = Σ_k u_k n(e_k)·(h − w)`; the
/// three normals are exact, the products are formed in floating point.
fn validity_rows(support: &Support) -> Vec<Vector> {
    let unit = |k: usize| -> Point { std::array::from_fn(|j| QSqrt5::integer(i64::from(j == k))) };
    let basis: [Vector; 3] = std::array::from_fn(|k| float::point(&support.normal_at(&unit(k))));
    let vertices = float::vertices();
    let h = vertices[support.contact()];
    vertices
        .iter()
        .map(|w| {
            let d = [h[0] - w[0], h[1] - w[1], h[2] - w[2]];
            std::array::from_fn(|k| float::dot(&basis[k], &d))
        })
        .collect()
}

/// The domain inequality `index` as a witness.
pub fn inequality(index: usize) -> Result<Witness, WitnessError> {
    Ok(Witness::Domain(DomainInequality::new(index)?))
}

/// The supports proposed at the exact view `u`: every oriented edge, the
/// radial directions `T_u(c)` and the plain directions `(c₁, c₂)` of every
/// vertex `c` with contact `c`, and directions inside the normal cone of
/// every vertex of the hole's shadow outline (the sum of its two outline
/// edge normals, and the 2:1 and 1:2 mixtures). The outline is found in
/// floating point; the directions are exact.
pub fn supports(u: &Point) -> Vec<Support> {
    let vertices = geometry::vertices();
    let mut out: Vec<Support> = Vec::new();
    for &[a, b] in geometry::edges() {
        for (from, to) in [(a, b), (b, a)] {
            out.push(Support::Edge(Edge::new(from, to).expect("an edge of the RID")));
        }
    }
    let projected: Vec<[QSqrt5; 2]> = vertices.iter().map(|v| geometry::screen(u, v)).collect();
    for (c, v) in vertices.iter().enumerate() {
        for normal in [projected[c].clone(), [v[0].clone(), v[1].clone()]] {
            if let Ok(direction) = Direction::new(normal, c) {
                out.push(Support::Direction(direction));
            }
        }
    }
    let outline = float::hull(&projected.iter().map(|p| [float::number(&p[0]), float::number(&p[1])]).collect::<Vec<_>>());
    let m = outline.len();
    for i in 0..m {
        let (a, c, b) = (outline[(i + m - 1) % m], outline[i], outline[(i + 1) % m]);
        // Outward normals (d₂, −d₁) of the counter-clockwise edges a→c and c→b.
        let normal = |p: usize, q: usize| {
            let d: [QSqrt5; 2] = std::array::from_fn(|j| &projected[q][j] - &projected[p][j]);
            [d[1].clone(), -&d[0]]
        };
        let (first, second) = (normal(a, c), normal(c, b));
        for (x, y) in [(1, 1), (2, 1), (1, 2)] {
            let mix: [QSqrt5; 2] =
                std::array::from_fn(|j| &first[j].scale(&frac(x, 1)) + &second[j].scale(&frac(y, 1)));
            if let Ok(direction) = Direction::new(mix, c) {
                out.push(Support::Direction(direction));
            }
        }
    }
    out
}

/// The float gap `((1 + |r|²) n·h − n·R̂(r)p)/(1 + |r|²) = n·h − n·R(r)p`
/// relative to `|n|` and the vertex size, for every plug vertex `p`.
fn relative_gaps(n: &Vector, h: &Vector, turned: &[Vector]) -> Vec<f64> {
    let size = float::dot(n, n).sqrt() * float::dot(h, h).sqrt();
    turned.iter().map(|p| (float::dot(n, h) - float::dot(n, p)) / size).collect()
}

/// The gaps proposed at the exact configuration `(u, r)`: supports of
/// [`supports`] that are valid at `u` up to the tolerance, with every plug
/// vertex whose gap is at most the tolerance (touching or protruding).
pub fn active(u: &Point, r: &Point) -> Vec<Gap> {
    let rf = float::point(r);
    let rotation = float::rotation(&rf);
    let turned: Vec<Vector> = float::vertices().iter().map(|v| float::apply(&rotation, v)).collect();
    let mut out = Vec::new();
    for support in supports(u) {
        let n = float::point(&support.normal_at(u));
        if float::dot(&n, &n) == 0.0 {
            continue;
        }
        let h = float::vertices()[support.contact()];
        if relative_gaps(&n, &h, float::vertices()).iter().any(|g| *g < -TOLERANCE) {
            continue;
        }
        let gaps = relative_gaps(&n, &h, &turned);
        for (plug, g) in gaps.iter().enumerate() {
            if *g <= TOLERANCE {
                out.push(Gap::new(support.clone(), plug).expect("a vertex"));
            }
        }
    }
    out
}

/// A float direction rounded to the dyadic grid, as an exact planar
/// direction (`None` when it rounds to zero).
fn rounded(n: [f64; 2]) -> Option<[QSqrt5; 2]> {
    let size = n[0].abs().max(n[1].abs());
    if !(size > 0.0) {
        return None;
    }
    let grid = |x: f64| (x / size * DIRECTION_GRID as f64).round() as i64;
    let (a, b) = (grid(n[0]), grid(n[1]));
    if a == 0 && b == 0 {
        return None;
    }
    Some([QSqrt5::from_rational(frac(a, DIRECTION_GRID)), QSqrt5::from_rational(frac(b, DIRECTION_GRID))])
}

/// Gaps proposed at a float configuration `(u, r)` of a cell:
///
/// - for the `count` plug vertices whose shadows protrude furthest beyond
///   (or come closest to) the outline of the hole's shadow, the planar
///   direction from the nearest outline point to the plug vertex, rounded
///   to the dyadic grid, with the float-extreme hole vertex as contact;
/// - the `count` oriented edges valid at `u` (in floating point) whose gap
///   is most negative, with the plug vertex attaining it.
///
/// Everything is decided later, exactly, by the lemma.
pub fn separating(u: &Vector, r: &Vector, count: usize) -> Vec<Gap> {
    let vertices = float::vertices();
    let rotation = float::rotation(r);
    let hole: Vec<[f64; 2]> = vertices.iter().map(|v| float::screen(u, v)).collect();
    let plug: Vec<[f64; 2]> = vertices.iter().map(|v| float::screen(u, &float::apply(&rotation, v))).collect();
    let outline = float::hull(&hole);
    let m = outline.len();
    // Outward edge lines n·x ≤ c of the counter-clockwise outline.
    let lines: Vec<([f64; 2], f64)> = (0..m)
        .map(|i| {
            let (p, q) = (hole[outline[i]], hole[outline[(i + 1) % m]]);
            let d = [q[0] - p[0], q[1] - p[1]];
            let length = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let n = [d[1] / length, -d[0] / length];
            (n, n[0] * p[0] + n[1] * p[1])
        })
        .collect();
    let excess = |x: &[f64; 2]| lines.iter().map(|(n, c)| n[0] * x[0] + n[1] * x[1] - c).fold(f64::NEG_INFINITY, f64::max);
    let mut ranked: Vec<(f64, usize)> = plug.iter().enumerate().map(|(p, x)| (excess(x), p)).collect();
    ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut out = Vec::new();
    for &(_, p) in ranked.iter().take(count) {
        let x = plug[p];
        let mut best = (f64::INFINITY, x);
        for i in 0..m {
            let (a, b) = (hole[outline[i]], hole[outline[(i + 1) % m]]);
            let d = [b[0] - a[0], b[1] - a[1]];
            let t = (((x[0] - a[0]) * d[0] + (x[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1])).clamp(0.0, 1.0);
            let c = [a[0] + t * d[0], a[1] + t * d[1]];
            let distance = ((x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2)).sqrt();
            if distance < best.0 {
                best = (distance, c);
            }
        }
        let mut n = [x[0] - best.1[0], x[1] - best.1[1]];
        if best.0 < ON_OUTLINE {
            let (i, _) = lines
                .iter()
                .enumerate()
                .map(|(i, (l, c))| (i, l[0] * x[0] + l[1] * x[1] - c))
                .fold((0, f64::NEG_INFINITY), |a, b| if b.1 > a.1 { b } else { a });
            n = lines[i].0;
        }
        let Some(direction) = rounded(n) else { continue };
        let nf = [float::number(&direction[0]), float::number(&direction[1])];
        let contact = (0..VERTEX_COUNT)
            .fold(0, |best, c| if nf[0] * hole[c][0] + nf[1] * hole[c][1] > nf[0] * hole[best][0] + nf[1] * hole[best][1] { c } else { best });
        let support = Support::Direction(Direction::new(direction, contact).expect("a nonzero direction"));
        out.push(Gap::new(support, p).expect("a vertex"));
    }
    let turned: Vec<Vector> = vertices.iter().map(|v| float::apply(&rotation, v)).collect();
    let mut edges: Vec<(f64, usize, usize, usize)> = Vec::new();
    for &[a, b] in geometry::edges() {
        for (from, to) in [(a, b), (b, a)] {
            let d = [0, 1, 2].map(|j| vertices[to][j] - vertices[from][j]);
            let n = float::cross(u, &d);
            if float::dot(&n, &n).sqrt() < END_ON {
                continue;
            }
            if relative_gaps(&n, &vertices[from], vertices).iter().any(|g| *g < -TOLERANCE) {
                continue;
            }
            let gaps = relative_gaps(&n, &vertices[from], &turned);
            let (plug, g) = gaps.iter().enumerate().fold((0, f64::INFINITY), |a, (p, g)| if *g < a.1 { (p, *g) } else { a });
            if g < 0.0 {
                edges.push((g, from, to, plug));
            }
        }
    }
    edges.sort_by(|a, b| a.0.total_cmp(&b.0).then((a.1, a.2, a.3).cmp(&(b.1, b.2, b.3))));
    for &(_, from, to, plug) in edges.iter().take(count) {
        out.push(Gap::new(Support::Edge(Edge::new(from, to).expect("an edge")), plug).expect("a vertex"));
    }
    out
}

/// The zoom variable that the offset is proportional to, so that setting it
/// to 0 gives the zoom's base points (in the no-fit set): `δ` for tubes and
/// B zooms, `ρ` for A zooms.
fn offset_scale(kind: Kind, base_dimension: usize) -> usize {
    match kind {
        Kind::Tube { .. } => base_dimension,
        Kind::A { .. } | Kind::Sheared { .. } => 1,
        Kind::B { .. } => 0,
    }
}

/// The exact configurations `(u, r)` at the base points of zoom `zoom`
/// of `cover` on a grid: the offset scale is 0 and every other variable the
/// image depends on takes `grid + 1` equally spaced values of its root
/// interval. These configurations lie in the zoom's no-fit set.
pub fn samples(cover: &ZoomCover, zoom: usize, grid: NonZeroUsize) -> Vec<(Point, Point)> {
    let c = &cover.zooms()[zoom];
    let zero = offset_scale(cover.kind(zoom), cover.parameters().shape.base_dimension());
    let used: Vec<usize> = (0..VARIABLES)
        .filter(|&j| j != zero && c.map().iter().any(|x| x.terms().any(|(e, _)| e[j] > 0 && e[zero] == 0)))
        .collect();
    let mut out = BTreeSet::new();
    let steps = grid.get();
    let count = (steps + 1).checked_pow(used.len() as u32).expect("the sample grid is bounded by the configuration");
    for index in 0..count {
        let mut y: Five<QSqrt5> = std::array::from_fn(|j| QSqrt5::from_rational(c.root()[j].lo().clone()));
        y[zero] = QSqrt5::zero();
        for (position, &j) in used.iter().enumerate() {
            let k = index / (steps + 1).pow(position as u32) % (steps + 1);
            let (lo, width) = (c.root()[j].lo(), c.root()[j].width());
            y[j] = QSqrt5::from_rational(lo + &(width * frac(k as i64, steps as i64)));
        }
        let u: Point = std::array::from_fn(|i| c.view().vector()[i].evaluate(&y));
        let r: Point = std::array::from_fn(|i| c.rotation()[i].evaluate(&y));
        out.insert((u, r));
    }
    out.into_iter().collect()
}

#[cfg(test)]
mod tests;
