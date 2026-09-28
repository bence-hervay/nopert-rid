//! An independent model of zooms, covers and cover files,
//! written from the definitions without calling the code under test (only the
//! base layers: exact numbers, polynomial arithmetic, the RID's vertices and
//! D's polynomials). Used by the adversarial tests of `cover` and
//! `cover`:
//!
//! - coordinates and their inverse (the arc coordinates);
//! - the zoom families of a cover, their roots, factors and maps
//!   (tube and point zooms and the hand-over);
//! - the covered set and the exact inclusion decision by the 32 corners;
//! - witness polynomials from the definition `g = (1+|r|²) n·h − n·R̂(r)p`
//!   with `R̂(r)p = (1−|r|²)p + 2(r·p)r + 2 r×p`;
//! - exact monomial division, tensor Bernstein coefficients by per-axis
//!   shift and basis change, and a recursive reading of cell trees.
//!
//! **Why 32 corners decide inclusion.** Each cover coordinate is
//! `zᵢ = α e + β v + γ(r)` with `γ` affine (arc planes) or `zᵢ` affine
//! (configuration coordinates). With `e = (1−2s)/(2+s)`, `v = 5t/(2+s)` and
//! `2 + s > 0`, `α e + β v ≥ L` is `α(1−2s) + 5βt − L(2+s) ≥ 0`, affine in
//! `(s, t)`; so the minimum and maximum over a box of each `zᵢ` are attained
//! at corners.
use crate::arithmetic::exact::{frac, q, QSqrt5, Q};
use crate::arithmetic::polynomial::{Exponents, Polynomial};
use crate::elimination::zoom::cover::{Parameters, Shape};
use crate::elimination::zoom::Coordinates;
use crate::elimination::witness::{Support, Witness};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::domain;
use crate::problem::geometry::{self, View};

pub(crate) type P5 = [QSqrt5; 5];
pub(crate) type Poly3 = [Polynomial; 3];

pub(crate) fn rat(x: &Q) -> QSqrt5 {
    QSqrt5::from_rational(x.clone())
}

pub(crate) fn int(n: i64) -> QSqrt5 {
    QSqrt5::integer(n)
}

/// `1/x` as `x̄/(a² − 5b²)`.
pub(crate) fn inverse_of(x: &QSqrt5) -> QSqrt5 {
    let (a, b) = (x.rational_part(), x.sqrt5_part());
    let norm = a * a - &(q(5) * b * b);
    assert!(!norm.is_zero());
    QSqrt5::new(a / &norm, -(b / &norm))
}

/// `a₁ = (−5 + 3√5)/10`, `a₃ = (−5 + √5)/10`.
pub(crate) fn arc_constants() -> (QSqrt5, QSqrt5) {
    (QSqrt5::new(frac(-1, 2), frac(3, 10)), QSqrt5::new(frac(-1, 2), frac(1, 10)))
}

fn sigma(coordinates: Coordinates) -> i64 {
    match coordinates {
        Coordinates::ArcPlane(crate::elimination::zoom::Sign::Plus) => 1,
        Coordinates::ArcPlane(crate::elimination::zoom::Sign::Minus) => -1,
        Coordinates::Configuration => 0,
    }
}

/// The view vector and rotation vector of coordinates `x` (polynomials).
pub(crate) fn view_rotation(coordinates: Coordinates, x: &[Polynomial; 5]) -> (Poly3, Poly3) {
    let k = |c: QSqrt5| Polynomial::constant(c);
    match coordinates {
        Coordinates::Configuration => (
            [x[0].clone(), x[1].clone(), k(int(1))],
            [x[2].clone(), x[3].clone(), x[4].clone()],
        ),
        Coordinates::ArcPlane(_) => {
            let s = int(sigma(coordinates));
            let (a1, a3) = arc_constants();
            let u = [&k(int(1)) - &x[0].scale(&int(2)), x[1].clone(), &k(int(2)) + &x[0]];
            // r = σ r_A + θ (σa₃, 1, −σa₁) + (ζ₁, 0, ζ₃).
            let r = [
                &(&k(&s * &a1) + &x[2].scale(&(&s * &a3))) + &x[3],
                x[2].clone(),
                &(&k(&s * &a3) - &x[2].scale(&(&s * &a1))) + &x[4],
            ];
            (u, r)
        }
    }
}

/// The coordinates of the configuration `c = (s, t, r)`; `None` for an arc
/// plane when `2 + s ≤ 0`.
pub(crate) fn coordinates_of(coordinates: Coordinates, c: &P5) -> Option<P5> {
    match coordinates {
        Coordinates::Configuration => Some(c.clone()),
        Coordinates::ArcPlane(_) => {
            let d = &int(2) + &c[0];
            if d.sign() != std::cmp::Ordering::Greater {
                return None;
            }
            let inv = inverse_of(&d);
            let s = int(sigma(coordinates));
            let (a1, a3) = arc_constants();
            let theta = c[3].clone();
            let e = &(&int(1) - &(&int(2) * &c[0])) * &inv;
            let v = &(&int(5) * &c[1]) * &inv;
            let z1 = &c[2] - &(&s * &(&a1 + &(&a3 * &theta)));
            let z3 = &c[4] - &(&s * &(&a3 - &(&a1 * &theta)));
            Some([e, v, theta, z1, z3])
        }
    }
}

/// The configuration `(s, t, r)` of exact coordinates `x` (`u₃ ≠ 0`).
pub(crate) fn configuration_of(coordinates: Coordinates, x: &P5) -> P5 {
    let (u, r) = exact_view_rotation(coordinates, x);
    let inv = inverse_of(&u[2]);
    [&u[0] * &inv, &u[1] * &inv, r[0].clone(), r[1].clone(), r[2].clone()]
}

pub(crate) fn exact_view_rotation(coordinates: Coordinates, x: &P5) -> ([QSqrt5; 3], [QSqrt5; 3]) {
    let constants: [Polynomial; 5] = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (u, r) = view_rotation(coordinates, &constants);
    let origin: P5 = std::array::from_fn(|_| QSqrt5::zero());
    (u.map(|p| p.evaluate(&origin)), r.map(|p| p.evaluate(&origin)))
}

pub(crate) fn no_fit_coordinates(coordinates: Coordinates) -> [usize; 3] {
    match coordinates {
        Coordinates::Configuration => [2, 3, 4],
        Coordinates::ArcPlane(_) => [1, 3, 4],
    }
}

// ---- Covers ---------------------------------------------------------------

pub(crate) type FaceOf = (usize, i8);

/// The faces of a unit box: per axis the lower face (end −1) then the upper
/// face (end +1).
pub(crate) fn faces_of(unit_box: &[[i8; 2]]) -> Vec<FaceOf> {
    let mut out = Vec::new();
    for (axis, &[lo, hi]) in unit_box.iter().enumerate() {
        if lo == -1 {
            out.push((axis, -1));
        }
        if hi == 1 {
            out.push((axis, 1));
        }
    }
    out
}

fn face_name((axis, side): FaceOf) -> String {
    format!("{axis}{}", if side > 0 { '+' } else { '-' })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Family {
    Tube,
    A,
    B,
    Sheared,
}

#[derive(Clone, Debug)]
pub(crate) struct ModelZoom {
    pub name: String,
    pub family: Family,
    pub base_face: Option<FaceOf>,
    pub offset_face: FaceOf,
    pub root: Vec<(Q, Q)>,
    pub factors: Vec<usize>,
    /// Cover coordinates `z(y)` and coordinates `x(y) = c₀ + M z(y)`.
    pub z: Vec<Polynomial>,
    pub x: [Polynomial; 5],
}

#[derive(Clone, Debug)]
pub(crate) struct ModelZoomCover {
    pub parameters: Parameters,
    pub zooms: Vec<ModelZoom>,
    /// The covered set in cover coordinates.
    pub bounds: Vec<(QSqrt5, QSqrt5)>,
    pub inverse: [[QSqrt5; 5]; 5],
    pub base_dimension: usize,
}

struct Vars(Vec<(Q, Q)>);

impl Vars {
    fn take(&mut self, lo: Q, hi: Q) -> Polynomial {
        let mut e = [0u8; 5];
        e[self.0.len()] = 1;
        self.0.push((lo, hi));
        Polynomial::monomial(e, int(1))
    }
    fn unit(&mut self, unit_box: &[[i8; 2]], face: FaceOf) -> Vec<Polynomial> {
        (0..unit_box.len())
            .map(|j| {
                if j == face.0 {
                    Polynomial::constant(int(i64::from(face.1)))
                } else {
                    self.take(q(i64::from(unit_box[j][0])), q(i64::from(unit_box[j][1])))
                }
            })
            .collect()
    }
}

fn times(f: &Polynomial, v: &[Polynomial]) -> Vec<Polynomial> {
    v.iter().map(|p| f.mul(p).unwrap()).collect()
}

/// Gauss–Jordan inverse over Q(√5), written independently.
pub(crate) fn invert5(m: &[[QSqrt5; 5]; 5]) -> [[QSqrt5; 5]; 5] {
    let mut a: Vec<Vec<QSqrt5>> = (0..5)
        .map(|i| (0..10).map(|j| if j < 5 { m[i][j].clone() } else if j - 5 == i { int(1) } else { QSqrt5::zero() }).collect())
        .collect();
    for col in 0..5 {
        let pivot = (col..5).find(|&r| !a[r][col].is_zero()).expect("invertible map");
        a.swap(col, pivot);
        let inv = inverse_of(&a[col][col]);
        a[col] = a[col].iter().map(|v| v * &inv).collect();
        for r in 0..5 {
            if r != col && !a[r][col].is_zero() {
                let f = a[r][col].clone();
                a[r] = (0..10).map(|j| &a[r][j] - &(&f * &a[col][j])).collect();
            }
        }
    }
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j + 5].clone()))
}

impl ModelZoomCover {
    pub(crate) fn new(p: &Parameters) -> Self {
        let centre: P5 = std::array::from_fn(|i| p.centre[i].0.clone());
        let map: [[QSqrt5; 5]; 5] = std::array::from_fn(|i| std::array::from_fn(|j| p.map[i][j].0.clone()));
        let mut zooms = Vec::new();
        let mut bounds = Vec::new();
        let base_dimension;
        let mut push = |name: String, family, base_face, offset_face, vars: Vars, factors: Vec<usize>, z: Vec<Polynomial>| {
            assert_eq!(z.len(), 5);
            assert_eq!(vars.0.len(), 5, "{name}: five zoom variables");
            let x: [Polynomial; 5] = std::array::from_fn(|i| {
                (0..5).fold(Polynomial::constant(centre[i].clone()), |s, j| &s + &z[j].scale(&map[i][j]))
            });
            zooms.push(ModelZoom { name, family, base_face, offset_face, root: vars.0, factors, z, x });
        };
        match &p.shape {
            Shape::Tube { base, offset, radius } => {
                let k = base.len();
                base_dimension = k;
                for [lo, hi] in base {
                    bounds.push((rat(&lo.0), rat(&hi.0)));
                }
                for &[lo, hi] in offset {
                    bounds.push((rat(&(&radius.0 * q(i64::from(lo)))), rat(&(&radius.0 * q(i64::from(hi))))));
                }
                for f in faces_of(offset) {
                    let mut v = Vars(Vec::new());
                    let b: Vec<Polynomial> = base.iter().map(|[lo, hi]| v.take(lo.0.clone(), hi.0.clone())).collect();
                    let delta = v.take(Q::zero(), radius.0.clone());
                    let unit = v.unit(offset, f);
                    let mut z = b;
                    z.extend(times(&delta, &unit));
                    push(format!("T o{}", face_name(f)), Family::Tube, None, f, v, vec![k], z);
                }
            }
            Shape::Point { base, offset, radii, radius, ratio, window } => {
                let k = base.len();
                base_dimension = k;
                let base_faces = faces_of(base);
                assert_eq!(radii.len(), base_faces.len());
                let delta_of = |g: FaceOf| radii[base_faces.iter().position(|&h| h == g).unwrap()].0.clone();
                for (j, &[lo, hi]) in base.iter().enumerate() {
                    let l = if lo == -1 { -delta_of((j, -1)) } else { Q::zero() };
                    let h = if hi == 1 { delta_of((j, 1)) } else { Q::zero() };
                    bounds.push((rat(&l), rat(&h)));
                }
                for &[lo, hi] in offset {
                    bounds.push((rat(&(&radius.0 * q(i64::from(lo)))), rat(&(&radius.0 * q(i64::from(hi))))));
                }
                let a_zoom = |g: FaceOf, f: FaceOf, limit: &Q, shear: Option<&Vec<Vec<crate::elimination::zoom::cover::Rational>>>| {
                    let mut v = Vars(Vec::new());
                    let mu = v.take(Q::zero(), delta_of(g));
                    let rho = v.take(Q::zero(), limit.clone());
                    let bh = v.unit(base, g);
                    let fh = v.unit(offset, f);
                    let zb = times(&mu, &bh);
                    let mut zf = times(&mu.mul(&rho).unwrap(), &fh);
                    if let Some(s) = shear {
                        // z_offset = z'_offset − shear · z_base.
                        for (l, zl) in zf.iter_mut().enumerate() {
                            for j in 0..k {
                                *zl = &*zl - &zb[j].scale(&rat(&s[l][j].0));
                            }
                        }
                    }
                    let mut z = zb;
                    z.extend(zf);
                    (v, z)
                };
                for &g in &base_faces {
                    for f in faces_of(offset) {
                        let (v, z) = a_zoom(g, f, &ratio.0, None);
                        push(format!("A b{} o{}", face_name(g), face_name(f)), Family::A, Some(g), f, v, vec![0, 1], z);
                    }
                }
                for f in faces_of(offset) {
                    let mut v = Vars(Vec::new());
                    let delta = v.take(Q::zero(), radius.0.clone());
                    let b: Vec<Polynomial> = base
                        .iter()
                        .map(|&[lo, hi]| v.take(q(i64::from(lo)) / &ratio.0, q(i64::from(hi)) / &ratio.0))
                        .collect();
                    let fh = v.unit(offset, f);
                    let mut z = times(&delta, &b);
                    z.extend(times(&delta, &fh));
                    push(format!("B o{}", face_name(f)), Family::B, None, f, v, vec![0], z);
                }
                if let Some(w) = window {
                    for g in w.faces.iter().map(|g| (g.axis, g.side)) {
                        for f in faces_of(offset) {
                            let (v, z) = a_zoom(g, f, &w.radius.0, Some(&w.shear));
                            push(format!("A' b{} o{}", face_name(g), face_name(f)), Family::Sheared, Some(g), f, v, vec![0, 1], z);
                        }
                    }
                }
            }
        }
        Self { parameters: p.clone(), zooms, bounds, inverse: invert5(&map), base_dimension }
    }

    pub(crate) fn centre(&self) -> P5 {
        std::array::from_fn(|i| self.parameters.centre[i].0.clone())
    }

    pub(crate) fn cover_coordinates(&self, x: &P5) -> P5 {
        let c = self.centre();
        std::array::from_fn(|i| (0..5).fold(QSqrt5::zero(), |s, j| &s + &(&self.inverse[i][j] * &(&x[j] - &c[j]))))
    }

    pub(crate) fn image(&self, z: &P5) -> P5 {
        let c = self.centre();
        std::array::from_fn(|i| (0..5).fold(c[i].clone(), |s, j| &s + &(&self.parameters.map[i][j].0 * &z[j])))
    }

    fn beyond_axis(&self) -> Option<usize> {
        self.parameters.beyond.map(|b| b.axis)
    }

    /// Whether cover coordinates `z` lie in the (extended) cover.
    pub(crate) fn z_inside(&self, z: &P5) -> bool {
        (0..5).all(|i| self.bounds[i].0 <= z[i] && (self.beyond_axis() == Some(i) || z[i] <= self.bounds[i].1))
    }

    /// Whether the configuration `c` lies in the (extended) cover.
    pub(crate) fn configuration_inside(&self, c: &P5) -> bool {
        coordinates_of(self.parameters.coordinates, c).is_some_and(|x| self.z_inside(&self.cover_coordinates(&x)))
    }

    /// Exactly whether every configuration of `b` lies in the cover: all 32
    /// corners do (module comment).
    pub(crate) fn box_inside(&self, b: &ConfigurationBox) -> bool {
        (0..32u8).all(|m| self.configuration_inside(&b.corner(m).map(|x| rat(&x))))
    }
}

// ---- Witnesses, division, Bernstein ------------------------------------------

fn dot3(a: &[Polynomial; 3], p: &[QSqrt5; 3]) -> Polynomial {
    (0..3).fold(Polynomial::zero(), |s, j| &s + &a[j].scale(&p[j]))
}

fn pdot(a: &Poly3, b: &Poly3) -> Polynomial {
    (0..3).fold(Polynomial::zero(), |s, j| &s + &a[j].mul(&b[j]).unwrap())
}

/// The witness polynomial from its definition, at a polynomial view `u`
/// (homogeneous) and rotation vector `r`.
pub(crate) fn witness_polynomial(w: &Witness, u: &Poly3, r: &Poly3) -> Polynomial {
    match w {
        Witness::Gap(gap) => {
            let v = geometry::vertices();
            let (n, h): (Poly3, &[QSqrt5; 3]) = match gap.support() {
                Support::Edge(e) => {
                    let d: [QSqrt5; 3] = std::array::from_fn(|j| &v[e.to()][j] - &v[e.from()][j]);
                    // n = u × d.
                    let n = [
                        &u[1].scale(&d[2]) - &u[2].scale(&d[1]),
                        &u[2].scale(&d[0]) - &u[0].scale(&d[2]),
                        &u[0].scale(&d[1]) - &u[1].scale(&d[0]),
                    ];
                    (n, &v[e.from()])
                }
                Support::Direction(d) => {
                    let [n1, n2] = d.screen_normal();
                    let n = [u[2].scale(n1), u[2].scale(n2), -&(&u[0].scale(n1) + &u[1].scale(n2))];
                    (n, &v[d.contact()])
                }
            };
            let p = &v[gap.plug()];
            let rr = pdot(r, r);
            let rp = dot3(r, p);
            let one = Polynomial::constant(int(1));
            // R̂(r)p = (1 − |r|²)p + 2(r·p)r + 2 r × p.
            let cross = [
                &r[1].scale(&p[2]) - &r[2].scale(&p[1]),
                &r[2].scale(&p[0]) - &r[0].scale(&p[2]),
                &r[0].scale(&p[1]) - &r[1].scale(&p[0]),
            ];
            let turned: Poly3 = std::array::from_fn(|j| {
                &(&(&one - &rr).scale(&p[j]) + &rp.mul(&r[j]).unwrap().scale(&int(2))) + &cross[j].scale(&int(2))
            });
            &(&one + &rr).mul(&dot3(&n, h)).unwrap() - &pdot(&n, &turned)
        }
        Witness::Domain(d) => -&domain::polynomial(d.index(), &View::Homogeneous(u.clone()), r).unwrap(),
    }
}

/// The exact value of a witness at the configuration of coordinates `x`
/// (its homogeneous view; `u₃ > 0` gives it the configuration's sign).
pub(crate) fn witness_value(w: &Witness, coordinates: Coordinates, x: &P5) -> QSqrt5 {
    let constants: [Polynomial; 5] = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (u, r) = view_rotation(coordinates, &constants);
    let origin: P5 = std::array::from_fn(|_| QSqrt5::zero());
    witness_polynomial(w, &u, &r).evaluate(&origin)
}

/// Validity of a gap's support at an exact view, by brute force.
pub(crate) fn valid_at(w: &Witness, u: &[QSqrt5; 3]) -> bool {
    let Witness::Gap(gap) = w else { return true };
    let v = geometry::vertices();
    let (n, h): ([QSqrt5; 3], usize) = match gap.support() {
        Support::Edge(e) => {
            let d: [QSqrt5; 3] = std::array::from_fn(|j| &v[e.to()][j] - &v[e.from()][j]);
            ([&(&u[1] * &d[2]) - &(&u[2] * &d[1]), &(&u[2] * &d[0]) - &(&u[0] * &d[2]), &(&u[0] * &d[1]) - &(&u[1] * &d[0])], e.from())
        }
        Support::Direction(d) => {
            let [n1, n2] = d.screen_normal();
            ([&u[2] * n1, &u[2] * n2, -&(&(&u[0] * n1) + &(&u[1] * n2))], d.contact())
        }
    };
    let dot = |a: &[QSqrt5; 3], b: &[QSqrt5; 3]| (0..3).fold(QSqrt5::zero(), |s, j| &s + &(&a[j] * &b[j]));
    let top = dot(&n, &v[h]);
    v.iter().all(|w| dot(&n, w) <= top)
}

/// The exact quotient by `y^a`, or `None` when some term is not a multiple.
pub(crate) fn divide(p: &Polynomial, a: &Exponents) -> Option<Polynomial> {
    let mut out = Polynomial::zero();
    for (e, c) in p.terms() {
        if (0..5).any(|j| e[j] < a[j]) {
            return None;
        }
        let e2: Exponents = std::array::from_fn(|j| e[j] - a[j]);
        out = &out + &Polynomial::monomial(e2, c.clone());
    }
    Some(out)
}

fn binomial(n: usize, k: usize) -> Q {
    (0..k).fold(Q::one(), |c, i| c * frac((n - i) as i64, (i + 1) as i64))
}

/// Whether every tensor Bernstein coefficient of `p` on the closed `cell`
/// (degree = `p`'s degree per variable) is strictly negative. Dense tensor;
/// per axis the shift `y = lo + w x` and then the change to the Bernstein
/// basis `b_i = Σ_{k≤i} C(i,k)/C(n,k) a_k`.
pub(crate) fn bernstein_negative(p: &Polynomial, cell: &[(Q, Q)]) -> bool {
    let mut n = [0usize; 5];
    for (e, _) in p.terms() {
        for j in 0..5 {
            n[j] = n[j].max(usize::from(e[j]));
        }
    }
    let mut stride = [1usize; 5];
    for j in (0..4).rev() {
        stride[j] = stride[j + 1] * (n[j + 1] + 1);
    }
    let size = stride[0] * (n[0] + 1);
    let mut t = vec![QSqrt5::zero(); size];
    for (e, c) in p.terms() {
        let offset: usize = (0..5).map(|j| usize::from(e[j]) * stride[j]).sum();
        t[offset] = c.clone();
    }
    for j in 0..5 {
        let (lo, hi) = &cell[j];
        let w = hi - lo;
        let nj = n[j];
        if nj == 0 {
            continue;
        }
        let lo_pow: Vec<Q> = (0..=nj).scan(Q::one(), |acc, i| {
            let v = acc.clone();
            if i < nj {
                *acc = &*acc * lo;
            }
            Some(v)
        }).collect();
        let w_pow: Vec<Q> = (0..=nj).scan(Q::one(), |acc, i| {
            let v = acc.clone();
            if i < nj {
                *acc = &*acc * &w;
            }
            Some(v)
        }).collect();
        for start in 0..size {
            if (start / stride[j]) % (nj + 1) != 0 {
                continue;
            }
            let a: Vec<QSqrt5> = (0..=nj).map(|m| t[start + m * stride[j]].clone()).collect();
            // Shift and scale: a'_k = Σ_{m≥k} a_m C(m,k) lo^{m−k} w^k.
            let shifted: Vec<QSqrt5> = (0..=nj)
                .map(|k| {
                    (k..=nj).fold(QSqrt5::zero(), |s, m| &s + &a[m].scale(&(binomial(m, k) * &lo_pow[m - k] * &w_pow[k])))
                })
                .collect();
            for i in 0..=nj {
                let b = (0..=i).fold(QSqrt5::zero(), |s, k| &s + &shifted[k].scale(&(binomial(i, k) / binomial(nj, k))));
                t[start + i * stride[j]] = b;
            }
        }
    }
    t.iter().all(|b| b.sign() == std::cmp::Ordering::Less)
}

/// The leaf cells of a tree, read recursively; `None` for a malformed tree.
pub(crate) fn tree_leaves(root: &[(Q, Q)], tree: &str) -> Option<Vec<Vec<(Q, Q)>>> {
    fn walk(cell: Vec<(Q, Q)>, bytes: &[u8], at: &mut usize, out: &mut Vec<Vec<(Q, Q)>>) -> Option<()> {
        let b = *bytes.get(*at)?;
        *at += 1;
        match b {
            b'.' => out.push(cell),
            b'0'..=b'4' => {
                let j = usize::from(b - b'0');
                let mid = (&cell[j].0 + &cell[j].1) / q(2);
                let mut lower = cell.clone();
                lower[j].1 = mid.clone();
                let mut upper = cell;
                upper[j].0 = mid;
                walk(lower, bytes, at, out)?;
                walk(upper, bytes, at, out)?;
            }
            _ => return None,
        }
        Some(())
    }
    let mut out = Vec::new();
    let mut at = 0;
    walk(root.to_vec(), tree.as_bytes(), &mut at, &mut out)?;
    (at == tree.len()).then_some(out)
}

pub(crate) fn corner(cell: &[(Q, Q)], mask: u32) -> P5 {
    std::array::from_fn(|j| rat(if mask >> j & 1 == 1 { &cell[j].1 } else { &cell[j].0 }))
}

pub(crate) fn evaluate_all(x: &[Polynomial], y: &P5) -> Vec<QSqrt5> {
    x.iter().map(|p| p.evaluate(y)).collect()
}

pub(crate) fn is_multilinear(p: &Polynomial) -> bool {
    p.terms().all(|(e, _)| e.iter().all(|&k| k <= 1))
}
