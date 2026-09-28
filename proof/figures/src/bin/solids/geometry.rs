//! Vectors, rotations, planar hulls and convex polyhedra in floating point.
//! These figures only illustrate; nothing here decides anything.

pub type V2 = [f64; 2];
pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale(k: f64, a: V3) -> V3 {
    [k * a[0], k * a[1], k * a[2]]
}

pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

pub fn unit(a: V3) -> V3 {
    scale(1.0 / dot(a, a).sqrt(), a)
}

pub fn apply(m: &M3, v: V3) -> V3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// The rotation of a nonzero quaternion `(w, x, y, z)`.
pub fn quaternion(q: [f64; 4]) -> M3 {
    let [w, x, y, z] = q;
    let n = w * w + x * x + y * y + z * z;
    [
        [(w * w + x * x - y * y - z * z) / n, 2.0 * (x * y - z * w) / n, 2.0 * (z * x + y * w) / n],
        [2.0 * (x * y + z * w) / n, (w * w - x * x + y * y - z * z) / n, 2.0 * (y * z - x * w) / n],
        [2.0 * (z * x - y * w) / n, 2.0 * (y * z + x * w) / n, (w * w - x * x - y * y + z * z) / n],
    ]
}

/// The rotation with rotation vector `r` in the paper's sense: the
/// quaternion `(1, r)`.
pub fn rotation(r: V3) -> M3 {
    quaternion([1.0, r[0], r[1], r[2]])
}

/// The convex hull of planar points, counter-clockwise, without collinear
/// points.
pub fn hull(points: &[V2]) -> Vec<V2> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    p.dedup();
    if p.len() < 3 {
        return p;
    }
    let turn = |o: V2, a: V2, b: V2| (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    let mut h: Vec<V2> = Vec::new();
    for pass in 0..2 {
        let start = h.len();
        let ordered: Box<dyn Iterator<Item = &V2>> = if pass == 0 { Box::new(p.iter()) } else { Box::new(p.iter().rev()) };
        for &q in ordered {
            while h.len() >= start + 2 && turn(h[h.len() - 2], h[h.len() - 1], q) <= 0.0 {
                h.pop();
            }
            h.push(q);
        }
        h.pop();
    }
    h
}

/// The supporting lines `n·x = h` of a counter-clockwise convex polygon,
/// with unit outward normals `n`.
pub fn sides(polygon: &[V2]) -> Vec<(V2, f64)> {
    (0..polygon.len())
        .map(|i| {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let l = dx.hypot(dy);
            let n = [dy / l, -dx / l];
            (n, n[0] * a[0] + n[1] * a[1])
        })
        .collect()
}

/// The largest ratio `n·p / h` over the sides of `hole` (which must contain
/// the origin in its interior) and the points `p` of `plug`: the least
/// scaling of the hole about the origin that contains the plug.
pub fn containment_scale(hole: &[V2], plug: &[V2]) -> f64 {
    let sides = sides(hole);
    plug.iter()
        .flat_map(|p| sides.iter().map(move |(n, h)| (n[0] * p[0] + n[1] * p[1]) / h))
        .fold(f64::NEG_INFINITY, f64::max)
}

/// The distance from `p` to the boundary of the convex polygon, positive
/// inside.
pub fn depth(polygon: &[V2], p: V2) -> f64 {
    sides(polygon).iter().map(|(n, h)| h - n[0] * p[0] - n[1] * p[1]).fold(f64::INFINITY, f64::min)
}

/// A convex polyhedron: its vertices and its faces, each face counter-
/// clockwise seen from outside.
pub struct Polyhedron {
    pub vertices: Vec<V3>,
    pub faces: Vec<Vec<usize>>,
}

impl Polyhedron {
    /// The convex polyhedron spanned by `vertices`, all of which must be
    /// extreme points; faces are found by brute force, merging vertices
    /// within `tolerance` of a supporting plane.
    pub fn new(vertices: Vec<V3>, tolerance: f64) -> Self {
        let n = vertices.len();
        let mut faces: Vec<Vec<usize>> = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        for i in 0..n {
            for j in i + 1..n {
                for k in j + 1..n {
                    let normal = cross(sub(vertices[j], vertices[i]), sub(vertices[k], vertices[i]));
                    if dot(normal, normal) < 1e-18 {
                        continue;
                    }
                    let normal = unit(normal);
                    let d: Vec<f64> = vertices.iter().map(|&v| dot(normal, sub(v, vertices[i]))).collect();
                    let outward = if d.iter().all(|&x| x <= tolerance) {
                        normal
                    } else if d.iter().all(|&x| x >= -tolerance) {
                        scale(-1.0, normal)
                    } else {
                        continue;
                    };
                    let on: Vec<usize> = (0..n).filter(|&m| d[m].abs() <= tolerance).collect();
                    if seen.insert(on.clone()) {
                        faces.push(order(&vertices, on, outward));
                    }
                }
            }
        }
        Polyhedron { vertices, faces }
    }

    /// The unit outward normal of face `f`.
    pub fn normal(&self, f: usize) -> V3 {
        let face = &self.faces[f];
        let v = |i: usize| self.vertices[face[i]];
        let mut n = [0.0; 3];
        for i in 1..face.len() - 1 {
            n = add(n, cross(sub(v(i), v(0)), sub(v(i + 1), v(0))));
        }
        unit(n)
    }

    /// Every edge `(a, b)` with `a < b`, with the two faces it borders.
    pub fn edges(&self) -> Vec<(usize, usize, [usize; 2])> {
        let mut map: std::collections::BTreeMap<(usize, usize), Vec<usize>> = Default::default();
        for (f, face) in self.faces.iter().enumerate() {
            for i in 0..face.len() {
                let (a, b) = (face[i], face[(i + 1) % face.len()]);
                map.entry((a.min(b), a.max(b))).or_default().push(f);
            }
        }
        map.into_iter().filter(|(_, f)| f.len() == 2).map(|((a, b), f)| (a, b, [f[0], f[1]])).collect()
    }

    /// The same polyhedron with every vertex mapped by `f`.
    pub fn map(&self, f: impl Fn(V3) -> V3) -> Self {
        Polyhedron { vertices: self.vertices.iter().map(|&v| f(v)).collect(), faces: self.faces.clone() }
    }
}

/// The face's vertices in counter-clockwise order about `outward`.
fn order(vertices: &[V3], face: Vec<usize>, outward: V3) -> Vec<usize> {
    let c = scale(1.0 / face.len() as f64, face.iter().fold([0.0; 3], |a, &i| add(a, vertices[i])));
    let x = unit(sub(vertices[face[0]], c));
    let y = cross(outward, x);
    let mut keyed: Vec<(f64, usize)> = face
        .into_iter()
        .map(|i| {
            let d = sub(vertices[i], c);
            (dot(d, y).atan2(dot(d, x)), i)
        })
        .collect();
    keyed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    keyed.into_iter().map(|(_, i)| i).collect()
}

/// All sign changes of `(a, b, c)` and their cyclic shifts, without
/// repetitions.
pub fn cyclic_signed(a: f64, b: f64, c: f64) -> Vec<V3> {
    let mut out: Vec<V3> = Vec::new();
    for signs in 0..8 {
        let s = |k: usize, x: f64| if signs >> k & 1 == 1 { -x } else { x };
        let p = [s(0, a), s(1, b), s(2, c)];
        for shift in 0..3 {
            let q = [p[shift], p[(shift + 1) % 3], p[(shift + 2) % 3]];
            if !out.iter().any(|o| sub(*o, q).iter().all(|x| x.abs() < 1e-12)) {
                out.push(q);
            }
        }
    }
    out
}

pub const PHI: f64 = 1.618_033_988_749_895;

/// The rhombicosidodecahedron with edge length 2.
pub fn rid() -> Polyhedron {
    let p = PHI;
    let mut v = cyclic_signed(1.0, 1.0, p * p * p);
    v.extend(cyclic_signed(p * p, p, 2.0 * p));
    v.extend(cyclic_signed(2.0 + p, 0.0, p * p));
    assert_eq!(v.len(), 60);
    Polyhedron::new(v, 1e-9)
}

/// Orthographic view: `right` and `up` span the drawing plane, `toward`
/// points to the viewer.
#[derive(Clone, Copy)]
pub struct Camera {
    pub right: V3,
    pub up: V3,
    pub toward: V3,
}

impl Camera {
    /// Looking at the origin from azimuth `azimuth` (about `z`, from `x`) and
    /// elevation `elevation`, with `z` drawn upwards.
    pub fn new(azimuth: f64, elevation: f64) -> Self {
        let toward = [elevation.cos() * azimuth.cos(), elevation.cos() * azimuth.sin(), elevation.sin()];
        let right = unit(cross([0.0, 0.0, 1.0], toward));
        let up = cross(toward, right);
        Camera { right, up, toward }
    }

    /// Looking along `-toward`, with `up_hint` drawn as nearly upwards as
    /// possible.
    pub fn facing(toward: V3, up_hint: V3) -> Self {
        let toward = unit(toward);
        let right = unit(cross(up_hint, toward));
        Camera { right, up: cross(toward, right), toward }
    }

    pub fn project(&self, v: V3) -> V2 {
        [dot(v, self.right), dot(v, self.up)]
    }
}
