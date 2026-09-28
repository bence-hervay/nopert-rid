//! An independent floating-point model of the problem, used only to flag
//! claims that are clearly false: the RID from its defining formula (not the
//! geometry module's list), rotations by quaternion products, orthogonal
//! projection onto the screen plane, planar convex hulls, and the
//! inequalities of `D` from their formulas with the rotation group generated
//! by closure. It never accepts anything on the components' behalf.
use crate::arithmetic::exact::{QSqrt5, Q};
use crate::problem::configuration::{ConfigurationBox, AXES, R, S, T};
use num_rational::BigRational;
use num_traits::ToPrimitive;
use std::sync::OnceLock;

pub type V3 = [f64; 3];
type V2 = [f64; 2];
type Quaternion = [f64; 4];

/// A claim is "clearly false" only when the model contradicts it by more
/// than this; the model's own rounding error is below `10⁻¹²`.
pub const MARGIN: f64 = 1e-7;

pub fn phi() -> f64 {
    (1.0 + 5f64.sqrt()) / 2.0
}

pub fn to_f64(x: &Q) -> f64 {
    BigRational::new(x.numer().clone(), x.denom().clone())
        .to_f64()
        .unwrap()
}

pub fn number(x: &QSqrt5) -> f64 {
    to_f64(x.rational_part()) + to_f64(x.sqrt5_part()) * 5f64.sqrt()
}

pub fn dot(a: &V3, b: &V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: &V3, b: &V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn scaled(a: &V3, k: f64) -> V3 {
    a.map(|x| x * k)
}

fn close(a: &V3, b: &V3) -> bool {
    (0..3).all(|j| (a[j] - b[j]).abs() < 1e-9)
}

/// The 60 vertices of the edge-length-2 RID from the defining formula: sign
/// changes and cyclic permutations of three points, duplicates removed.
pub fn vertices() -> &'static [V3] {
    static VERTICES: OnceLock<Vec<V3>> = OnceLock::new();
    VERTICES.get_or_init(|| {
        let f = phi();
        let mut out: Vec<V3> = Vec::new();
        for base in [[1.0, 1.0, f * f * f], [f * f, f, 2.0 * f], [f + 2.0, 0.0, f * f]] {
            for signs in 0..8 {
                let signed: V3 = std::array::from_fn(|j| {
                    if signs & (1 << j) == 0 {
                        base[j]
                    } else {
                        -base[j]
                    }
                });
                for shift in 0..3 {
                    let v: V3 = std::array::from_fn(|j| signed[(j + shift) % 3]);
                    if !out.iter().any(|w| close(w, &v)) {
                        out.push(v);
                    }
                }
            }
        }
        assert_eq!(out.len(), 60);
        out
    })
}

pub fn product(a: &Quaternion, b: &Quaternion) -> Quaternion {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// The 120 unit quaternions of the RID's rotations, by closure under products
/// of three generators.
pub fn group() -> &'static [Quaternion] {
    static GROUP: OnceLock<Vec<Quaternion>> = OnceLock::new();
    GROUP.get_or_init(|| {
        let f = phi();
        let generators = [[0.0, 1.0, 0.0, 0.0], [0.5; 4], [0.0, 0.5, f / 2.0, 0.5 / f]];
        let mut out: Vec<Quaternion> = vec![[1.0, 0.0, 0.0, 0.0]];
        let mut i = 0;
        while i < out.len() {
            for g in &generators {
                let p = product(&out[i], g);
                if !out.iter().any(|q| (0..4).all(|j| (q[j] - p[j]).abs() < 1e-9)) {
                    out.push(p);
                }
            }
            i += 1;
            assert!(out.len() <= 120, "the generators do not generate the RID's group");
        }
        assert_eq!(out.len(), 120);
        out
    })
}

/// `R(r)p` by the quaternion product `(1, r)(0, p)(1, -r) / (1 + |r|²)`.
pub fn rotate(r: &V3, p: &V3) -> V3 {
    let q = [1.0, r[0], r[1], r[2]];
    let conjugate = [1.0, -r[0], -r[1], -r[2]];
    let x = product(&product(&q, &[0.0, p[0], p[1], p[2]]), &conjugate);
    let scale = 1.0 + dot(r, r);
    [x[1] / scale, x[2] / scale, x[3] / scale]
}

/// A configuration `(s, t, r)` in floating point.
#[derive(Clone, Copy, Debug)]
pub struct Configuration {
    pub s: f64,
    pub t: f64,
    pub r: V3,
}

impl Configuration {
    pub fn view(&self) -> V3 {
        [self.s, self.t, 1.0]
    }
}

/// An orthonormal basis of the plane perpendicular to `u`.
fn screen_basis(u: &V3) -> (V3, V3) {
    let axis = (0..3)
        .min_by(|&i, &j| u[i].abs().total_cmp(&u[j].abs()))
        .unwrap();
    let mut a = [0.0; 3];
    a[axis] = 1.0;
    let e1 = cross(u, &a);
    let e1 = scaled(&e1, 1.0 / dot(&e1, &e1).sqrt());
    let unit = scaled(u, 1.0 / dot(u, u).sqrt());
    (e1, cross(&unit, &e1))
}

fn cross2(o: &V2, a: &V2, b: &V2) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

/// The convex hull, counter-clockwise, by Andrew's monotone chain.
fn hull(points: &[V2]) -> Vec<V2> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    let mut out: Vec<V2> = Vec::new();
    for pass in 0..2 {
        let start = out.len();
        for x in &p {
            while out.len() >= start + 2 && cross2(&out[out.len() - 2], &out[out.len() - 1], x) <= 0.0 {
                out.pop();
            }
            out.push(*x);
        }
        out.pop();
        if pass == 0 {
            p.reverse();
        }
    }
    out
}

/// How far `x` lies outside the convex polygon: positive outside, minus the
/// distance to the boundary inside.
fn beyond(polygon: &[V2], x: &V2) -> f64 {
    (0..polygon.len())
        .map(|i| {
            let (a, b) = (&polygon[i], &polygon[(i + 1) % polygon.len()]);
            let length = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
            -cross2(a, b, x) / length
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

/// The largest distance by which a plug vertex's shadow lies outside the
/// hole's shadow (orthogonal projection along `u`). Positive: the
/// configuration is not a fit. Negative: every plug vertex is inside.
pub fn excursion(c: &Configuration) -> f64 {
    let (e1, e2) = screen_basis(&c.view());
    let project = |x: &V3| [dot(x, &e1), dot(x, &e2)];
    let hole: Vec<V2> = vertices().iter().map(project).collect();
    let polygon = hull(&hole);
    vertices()
        .iter()
        .map(|p| beyond(&polygon, &project(&rotate(&c.r, p))))
        .fold(f64::NEG_INFINITY, f64::max)
}

/// `scalar((0, u)(1, r)g)`.
fn fold_scalar(c: &Configuration, g: &Quaternion) -> f64 {
    let u = c.view();
    let x = product(&product(&[0.0, u[0], u[1], u[2]], &[1.0, c.r[0], c.r[1], c.r[2]]), g);
    x[0]
}

/// The values of the 73 inequalities `c_i ≤ 0` of `D`, from their stated
/// formulas, in the documented index order; fold `13 + k` uses the group
/// element `symmetries[k]` (given in floating point by the caller).
pub fn inequalities(c: &Configuration, symmetries: &[Quaternion]) -> Vec<f64> {
    let f = phi();
    let mut out = vec![f * c.s + f * f * c.t - 1.0];
    for axis in 0..3 {
        for first in [-1.0, 1.0] {
            for second in [-1.0, 1.0] {
                out.push(first * c.r[axis] + second * (f - 1.0) * c.r[(axis + 1) % 3] + f - 2.0);
            }
        }
    }
    let norm = 1.0 + c.s * c.s + c.t * c.t;
    out.extend(symmetries.iter().map(|g| fold_scalar(c, g).powi(2) - norm));
    out
}

/// How far the configuration lies outside `D`, from the independent group:
/// the largest value of the triangle, rotation and fold inequalities (and of
/// the root box's bounds). Positive: outside `D`.
pub fn outside_domain(c: &Configuration) -> f64 {
    let f = phi();
    let mut worst = f * c.s + f * f * c.t - 1.0;
    for axis in 0..3 {
        let (a, b) = (c.r[axis], c.r[(axis + 1) % 3]);
        worst = worst.max(a.abs() + (f - 1.0) * b.abs() + f - 2.0);
    }
    let norm = 1.0 + c.s * c.s + c.t * c.t;
    for g in group() {
        worst = worst.max(fold_scalar(c, g).powi(2) - norm);
    }
    let bounds = [-c.s, c.s - 2.0 / 3.0, -c.t, c.t - 0.4];
    bounds.iter().chain(c.r.map(|x| x.abs() - 0.4).iter()).fold(worst, |w, &x| w.max(x))
}

/// The signed distance of the rotated plug vertex beyond the line of the
/// ordered edge `from → to` of the hole, on the side of `n = u × (to - from)`.
pub fn plug_beyond(c: &Configuration, from: &V3, to: &V3, plug: &V3) -> f64 {
    let d: V3 = std::array::from_fn(|j| to[j] - from[j]);
    let n = cross(&c.view(), &d);
    let length = dot(&n, &n).sqrt();
    (dot(&n, &rotate(&c.r, plug)) - dot(&n, from)) / length
}

/// How far the hole extends beyond the line of the ordered edge on the side
/// of `n = u × (to - from)`: at most about zero exactly when the edge supports
/// the hole's shadow with outer normal `n`.
pub fn hole_beyond(c: &Configuration, from: &V3, to: &V3) -> f64 {
    let d: V3 = std::array::from_fn(|j| to[j] - from[j]);
    let n = cross(&c.view(), &d);
    let length = dot(&n, &n).sqrt();
    vertices()
        .iter()
        .map(|w| (dot(&n, w) - dot(&n, from)) / length)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// The 32 corners of the box, its midpoint and `interior` pseudo-random points.
pub fn samples(b: &ConfigurationBox, seed: u64, interior: usize) -> Vec<Configuration> {
    let point = |x: [f64; AXES]| Configuration {
        s: x[S],
        t: x[T],
        r: R.map(|j| x[j]),
    };
    let lo: [f64; AXES] = std::array::from_fn(|j| to_f64(b.axes()[j].lo()));
    let hi: [f64; AXES] = std::array::from_fn(|j| to_f64(b.axes()[j].hi()));
    let mut out: Vec<Configuration> = (0..1u8 << AXES)
        .map(|mask| point(std::array::from_fn(|j| if mask & (1 << j) == 0 { lo[j] } else { hi[j] })))
        .collect();
    out.push(point(std::array::from_fn(|j| (lo[j] + hi[j]) / 2.0)));
    let mut state = seed;
    for _ in 0..interior {
        out.push(point(std::array::from_fn(|j| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
            lo[j] + (hi[j] - lo[j]) * unit
        })));
    }
    out
}
