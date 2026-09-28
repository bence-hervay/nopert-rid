//! The rhombicosidodecahedron over Q(√5): its vertices in a fixed order, its
//! edges and rotation group, quaternion and Cayley rotations, the screen
//! projection and the view maps. Everything is exact.
use crate::arithmetic::exact::{frac, QSqrt5};
use crate::arithmetic::polynomial::{Polynomial, PolynomialError};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::OnceLock;

/// An exact vector of Q(√5)³.
pub type Point = [QSqrt5; 3];
/// A quaternion `(w, x, y, z)`, scalar first, over Q(√5).
pub type Quaternion = [QSqrt5; 4];
/// A 3×3 matrix, indexed `[row][column]`.
pub type Matrix<T> = [[T; 3]; 3];

pub const VERTEX_COUNT: usize = 60;
pub const EDGE_COUNT: usize = 120;
pub const SYMMETRY_COUNT: usize = 60;

/// Literal vertex coordinates: the pair `[a, b]` denotes `(a + b√5)/2`.
type VertexCode = [[i64; 2]; 3];
/// Literal quaternion coordinates: the pair `[a, b]` denotes `(a + b√5)/4`.
type SymmetryCode = [[i64; 2]; 4];

fn number(a: i64, b: i64, denominator: i64) -> QSqrt5 {
    QSqrt5::new(frac(a, denominator), frac(b, denominator))
}

/// The literal coordinates of the 60 vertices; their lexicographic order is
/// the fixed vertex order.
fn vertex_codes() -> &'static [VertexCode] {
    static CODES: OnceLock<Vec<VertexCode>> = OnceLock::new();
    CODES.get_or_init(|| {
        let mut out = BTreeSet::new();
        for base in [
            [[2, 0], [2, 0], [4, 2]],
            [[3, 1], [1, 1], [2, 2]],
            [[5, 1], [0, 0], [3, 1]],
        ] {
            for signs in 0..8 {
                let signed: VertexCode = std::array::from_fn(|j| {
                    let sign = if signs & (1 << j) == 0 { 1 } else { -1 };
                    [sign * base[j][0], sign * base[j][1]]
                });
                for shift in 0..3 {
                    out.insert(std::array::from_fn(|j| signed[(j + shift) % 3]));
                }
            }
        }
        assert_eq!(out.len(), VERTEX_COUNT);
        out.into_iter().collect()
    })
}

/// The 60 vertices of the edge-length-2 RID as exact points, in the fixed order.
pub fn vertices() -> &'static [Point] {
    static VERTICES: OnceLock<Vec<Point>> = OnceLock::new();
    VERTICES.get_or_init(|| {
        vertex_codes()
            .iter()
            .map(|v| std::array::from_fn(|j| number(v[j][0], v[j][1], 2)))
            .collect()
    })
}

/// The 120 edges as index pairs `[a, b]` with `a < b`, sorted: exactly the
/// vertex pairs at distance 2.
pub fn edges() -> &'static [[usize; 2]] {
    static EDGES: OnceLock<Vec<[usize; 2]>> = OnceLock::new();
    EDGES.get_or_init(|| {
        let v = vertices();
        let four = QSqrt5::integer(4);
        let mut out = Vec::new();
        for a in 0..VERTEX_COUNT {
            for b in (a + 1)..VERTEX_COUNT {
                let d = difference(&v[b], &v[a]);
                if dot(&d, &d) == four {
                    out.push([a, b]);
                }
            }
        }
        assert_eq!(out.len(), EDGE_COUNT);
        out
    })
}

/// Whether the unordered pair `{a, b}` is an edge. Indices out of range are not.
pub fn is_edge(a: usize, b: usize) -> bool {
    edges().binary_search(&[a.min(b), a.max(b)]).is_ok()
}

/// The 60 rotations of the RID as unit quaternions, one of each antipodal
/// pair: the one whose first nonzero coordinate is positive. Their order is
/// the lexicographic order of the literal pairs `[a, b]` denoting `(a + b√5)/4`.
pub fn symmetries() -> &'static [Quaternion] {
    static GROUP: OnceLock<Vec<Quaternion>> = OnceLock::new();
    GROUP.get_or_init(|| {
        let decode = |v: &SymmetryCode| -> Quaternion {
            std::array::from_fn(|j| number(v[j][0], v[j][1], 4))
        };
        let mut codes = BTreeSet::new();
        let mut insert = |mut v: SymmetryCode| {
            let first = decode(&v)
                .into_iter()
                .find(|x| !x.is_zero())
                .expect("a group element is nonzero");
            if first.sign() == Ordering::Less {
                v = v.map(|[a, b]| [-a, -b]);
            }
            codes.insert(v);
        };
        for j in 0..4 {
            let mut v = [[0, 0]; 4];
            v[j] = [4, 0];
            insert(v);
        }
        for signs in 0..16 {
            insert(std::array::from_fn(|j| {
                [if signs & (1 << j) == 0 { 2 } else { -2 }, 0]
            }));
        }
        // (0, 1/2, φ/2, (φ-1)/2) under even permutations and all signs.
        let seed = [[0, 0], [2, 0], [1, 1], [-1, 1]];
        for permutation in even_permutations() {
            for signs in 0..16 {
                insert(std::array::from_fn(|j| {
                    let sign = if signs & (1 << j) == 0 { 1 } else { -1 };
                    let [a, b] = seed[permutation[j]];
                    [sign * a, sign * b]
                }));
            }
        }
        assert_eq!(codes.len(), SYMMETRY_COUNT);
        codes.iter().map(decode).collect()
    })
}

fn even_permutations() -> Vec<[usize; 4]> {
    let mut out = Vec::new();
    for code in 0..256usize {
        let p: [usize; 4] = std::array::from_fn(|j| (code >> (2 * j)) & 3);
        let distinct = (0..4).all(|i| ((i + 1)..4).all(|j| p[i] != p[j]));
        let inversions = (0..4)
            .flat_map(|i| ((i + 1)..4).map(move |j| (i, j)))
            .filter(|&(i, j)| p[i] > p[j])
            .count();
        if distinct && inversions % 2 == 0 {
            out.push(p);
        }
    }
    out
}

// ---- Exact vectors and quaternions -------------------------------------

pub fn dot(a: &Point, b: &Point) -> QSqrt5 {
    (0..3).fold(QSqrt5::zero(), |sum, j| &sum + &(&a[j] * &b[j]))
}

pub fn cross(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&a[j] * &b[k]) - &(&a[k] * &b[j])
    })
}

/// `a - b`.
pub fn difference(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|j| &a[j] - &b[j])
}

/// The Hamilton product `a b`.
pub fn quaternion_product(a: &Quaternion, b: &Quaternion) -> Quaternion {
    let (u, v) = (vector_part(a), vector_part(b));
    let scalar = &(&a[0] * &b[0]) - &dot(&u, &v);
    let w = cross(&u, &v);
    let [x, y, z]: Point =
        std::array::from_fn(|j| &(&(&a[0] * &v[j]) + &(&b[0] * &u[j])) + &w[j]);
    [scalar, x, y, z]
}

/// The quaternion conjugate `q̄ = (w, -v)` of `q = (w, v)`.
pub fn conjugate(quaternion: &Quaternion) -> Quaternion {
    let [w, x, y, z] = quaternion;
    [w.clone(), -x, -y, -z]
}

fn vector_part(quaternion: &Quaternion) -> Point {
    let [_, x, y, z] = quaternion;
    [x.clone(), y.clone(), z.clone()]
}

/// `q (0, p) q̄` for the quaternion `q`: this is `|q|² R(q) p`, where `R(q)`
/// is the rotation of `q ≠ 0`.
pub fn rotate(quaternion: &Quaternion, p: &Point) -> Point {
    let pure = [QSqrt5::zero(), p[0].clone(), p[1].clone(), p[2].clone()];
    let turned = quaternion_product(&quaternion_product(quaternion, &pure), &conjugate(quaternion));
    vector_part(&turned)
}

/// The sign `±1` of the cross-product matrix `[v]×` at `(i, j)`, `i ≠ j`:
/// its entry there is `sign · v_k` with `k` the third index.
fn cross_sign(i: usize, j: usize) -> i64 {
    if (j + 3 - i) % 3 == 1 {
        -1
    } else {
        1
    }
}

// ---- Polynomial vectors ------------------------------------------------

/// `a · p` for a polynomial vector `a` and an exact point `p`.
pub fn polynomial_dot(a: &[Polynomial; 3], p: &Point) -> Polynomial {
    (0..3).fold(Polynomial::zero(), |sum, j| &sum + &a[j].scale(&p[j]))
}

/// `a × p` for a polynomial vector `a` and an exact point `p`.
pub fn polynomial_cross(a: &[Polynomial; 3], p: &Point) -> [Polynomial; 3] {
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &a[j].scale(&p[k]) - &a[k].scale(&p[j])
    })
}

/// `|r|²` for a polynomial vector.
pub fn squared_norm(r: &[Polynomial; 3]) -> Result<Polynomial, PolynomialError> {
    let mut sum = Polynomial::zero();
    for x in r {
        sum = &sum + &x.mul(x)?;
    }
    Ok(sum)
}

/// The Cayley rotation numerator `R̂(r) = (1 + |r|²) R(r)` of the quaternion
/// `(1, r)`: `(1 - |r|²) I + 2 r rᵀ + 2 [r]×`, for polynomial `r`.
pub fn cayley_matrix(r: &[Polynomial; 3]) -> Result<Matrix<Polynomial>, PolynomialError> {
    let diagonal = &Polynomial::constant(QSqrt5::one()) - &squared_norm(r)?;
    let entry = |i: usize, j: usize| -> Result<Polynomial, PolynomialError> {
        let outer = r[i].mul(&r[j])?.scale(&QSqrt5::integer(2));
        Ok(if i == j {
            &outer + &diagonal
        } else {
            &outer + &r[3 - i - j].scale(&QSqrt5::integer(2 * cross_sign(i, j)))
        })
    };
    Ok([
        [entry(0, 0)?, entry(0, 1)?, entry(0, 2)?],
        [entry(1, 0)?, entry(1, 1)?, entry(1, 2)?],
        [entry(2, 0)?, entry(2, 1)?, entry(2, 2)?],
    ])
}

/// `m p` for a polynomial matrix and an exact point.
pub fn polynomial_matrix_point(m: &Matrix<Polynomial>, p: &Point) -> [Polynomial; 3] {
    std::array::from_fn(|i| polynomial_dot(&m[i], p))
}

// ---- Screen and views --------------------------------------------------

/// The screen map `T_u(X) = (u₃X₁ - u₁X₃, u₃X₂ - u₂X₃)`. Its kernel contains
/// `u`; for `u₃ ≠ 0` it is `u₃` times the oblique projection
/// `(X₁ - sX₃, X₂ - tX₃)` of the affine view `(s, t, 1) = u/u₃`.
pub fn screen(u: &Point, x: &Point) -> [QSqrt5; 2] {
    [
        &(&u[2] * &x[0]) - &(&u[0] * &x[2]),
        &(&u[2] * &x[1]) - &(&u[1] * &x[2]),
    ]
}

/// How a zoom presents the viewing line as a polynomial vector `u`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum View {
    /// `u = (s, t, 1)`.
    Affine { s: Polynomial, t: Polynomial },
    /// `u` given homogeneously. The zoom must prove `u₃ > 0` on its domain;
    /// the viewing line is then that of the affine view `(u₁/u₃, u₂/u₃, 1)`.
    Homogeneous([Polynomial; 3]),
}

impl View {
    /// The view vector `u`.
    pub fn vector(&self) -> [Polynomial; 3] {
        match self {
            View::Affine { s, t } => [s.clone(), t.clone(), Polynomial::constant(QSqrt5::one())],
            View::Homogeneous(u) => u.clone(),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
