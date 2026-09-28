use super::*;
use crate::arithmetic::exact::{q, Q};
use crate::arithmetic::polynomial::VARIABLES;
use std::collections::{BTreeMap, BTreeSet};

type C = QSqrt5;

fn int(n: i64) -> C {
    C::integer(n)
}

fn rat(x: Q) -> C {
    C::from_rational(x)
}

fn var(j: usize) -> Polynomial {
    Polynomial::variables()[j].clone()
}

fn at(p: &Polynomial, x: &[Q; VARIABLES]) -> C {
    p.evaluate(&x.clone().map(C::from_rational))
}

/// Deterministic generator of test data.
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    /// A rational in [-bound, bound] with a random denominator up to 2^20.
    fn rational(&mut self, bound: i64) -> Q {
        let denominator = 1 + self.below(1 << 20) as i64;
        let numerator =
            self.below((2 * bound * denominator + 1) as u64) as i64 - bound * denominator;
        frac(numerator, denominator)
    }
    /// A number of Q(√5) with rational and radical parts in [-bound, bound].
    fn number(&mut self, bound: i64) -> C {
        C::new(self.rational(bound), self.rational(bound))
    }
    fn point(&mut self, bound: i64) -> Point {
        std::array::from_fn(|_| self.number(bound))
    }
    fn quaternion(&mut self, bound: i64) -> Quaternion {
        std::array::from_fn(|_| self.number(bound))
    }
}

// ---- Test-only helpers, also used by the tests of other modules ----------

/// The index of a rotation of the group given by either of its two unit quaternions.
pub(crate) fn symmetry_index(quaternion: &Quaternion) -> Option<usize> {
    let negated = quaternion.clone().map(|c| -&c);
    symmetries()
        .iter()
        .position(|g| g == quaternion || *g == negated)
}

/// The rotation numerator `|q|² R(q) = (w² - |v|²) I + 2 v vᵀ + 2 w [v]×` of
/// `q = (w, v)`, the matrix form that the tests compare `rotate` and
/// `cayley_matrix` with.
fn rotation_matrix(quaternion: &Quaternion) -> Matrix<QSqrt5> {
    let (w, v) = (&quaternion[0], vector_part(quaternion));
    let diagonal = &(w * w) - &dot(&v, &v);
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let outer = (&v[i] * &v[j]).scale(&q(2));
            if i == j {
                &outer + &diagonal
            } else {
                &outer + &(w * &v[3 - i - j]).scale(&q(2 * cross_sign(i, j)))
            }
        })
    })
}

fn matrix_point(m: &Matrix<QSqrt5>, p: &Point) -> Point {
    std::array::from_fn(|i| dot(&m[i], p))
}

fn vertex_lookup() -> BTreeMap<Point, usize> {
    vertices()
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect()
}

fn transpose<T: Clone>(m: &Matrix<T>) -> Matrix<T> {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i].clone()))
}

fn exact_product(a: &Matrix<C>, b: &Matrix<C>) -> Matrix<C> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).fold(int(0), |s, k| &s + &(&a[i][k] * &b[k][j])))
    })
}

fn exact_determinant(m: &Matrix<C>) -> C {
    dot(&m[0], &cross(&m[1], &m[2]))
}

fn polynomial_product(a: &Matrix<Polynomial>, b: &Matrix<Polynomial>) -> Matrix<Polynomial> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..3).fold(Polynomial::zero(), |s, k| &s + &a[i][k].mul(&b[k][j]).unwrap())
        })
    })
}

fn polynomial_determinant(m: &Matrix<Polynomial>) -> Polynomial {
    let product = |a: &Polynomial, b: &Polynomial| a.mul(b).unwrap();
    let minor = |i: usize, j: usize, k: usize, l: usize| {
        &product(&m[1][i], &m[2][j]) - &product(&m[1][k], &m[2][l])
    };
    &(&product(&m[0][0], &minor(1, 2, 2, 1)) - &product(&m[0][1], &minor(0, 2, 2, 0)))
        + &product(&m[0][2], &minor(0, 1, 1, 0))
}

#[test]
fn vertex_codes_are_sorted_and_decode_to_sixty_centred_vertices_on_one_sphere() {
    let codes = vertex_codes();
    let points = vertices();
    assert_eq!(codes.len(), VERTEX_COUNT);
    assert_eq!(points.len(), VERTEX_COUNT);
    assert!(codes.windows(2).all(|pair| pair[0] < pair[1]));
    let lookup = vertex_lookup();
    assert_eq!(lookup.len(), VERTEX_COUNT);
    for (code, point) in codes.iter().zip(points) {
        for j in 0..3 {
            assert_eq!(point[j], C::new(frac(code[j][0], 2), frac(code[j][1], 2)));
        }
        assert_eq!(dot(point, point), C::new(q(11), q(4)));
        assert!(lookup.contains_key(&point.clone().map(|c| -&c)));
    }
    // The three generating vertices of the vertex list.
    for literal in [
        [int(1), int(1), C::new(q(2), q(1))],
        [
            C::new(frac(3, 2), frac(1, 2)),
            C::new(frac(1, 2), frac(1, 2)),
            C::new(q(1), q(1)),
        ],
        [
            C::new(frac(5, 2), frac(1, 2)),
            int(0),
            C::new(frac(3, 2), frac(1, 2)),
        ],
    ] {
        assert!(lookup.contains_key(&literal));
    }
    // The fixed order starts and ends with these vertices.
    assert_eq!(codes[0], [[-5, -1], [0, 0], [-3, -1]]);
    assert_eq!(codes[59], [[5, 1], [0, 0], [3, 1]]);
}

#[test]
fn edges_are_exactly_the_shortest_distances_and_every_vertex_has_degree_four() {
    let v = vertices();
    let list = edges();
    assert_eq!(list.len(), EDGE_COUNT);
    assert!(list.windows(2).all(|pair| pair[0] < pair[1]));
    let mut degree = [0usize; VERTEX_COUNT];
    for a in 0..VERTEX_COUNT {
        for b in 0..VERTEX_COUNT {
            if a == b {
                assert!(!is_edge(a, b));
                continue;
            }
            let d = difference(&v[a], &v[b]);
            let comparison = dot(&d, &d).cmp(&int(4));
            assert_ne!(comparison, Ordering::Less, "pair {a} {b}");
            assert_eq!(is_edge(a, b), comparison == Ordering::Equal);
            assert_eq!(is_edge(a, b), is_edge(b, a));
            if a < b && comparison == Ordering::Equal {
                assert!(list.contains(&[a, b]));
                degree[a] += 1;
                degree[b] += 1;
            }
        }
    }
    assert!(degree.iter().all(|d| *d == 4));
    assert!(!is_edge(0, VERTEX_COUNT));
    assert!(!is_edge(VERTEX_COUNT, usize::MAX));
}

fn coplanar(points: &[&Point]) -> bool {
    let base = points[0];
    let normal = cross(&difference(points[1], base), &difference(points[2], base));
    points[3..]
        .iter()
        .all(|p| dot(&normal, &difference(p, base)).is_zero())
}

/// Cycles of the edge graph with `length` vertices, as sorted vertex sets.
fn cycles(length: usize) -> BTreeSet<Vec<usize>> {
    fn extend(path: &mut Vec<usize>, length: usize, out: &mut BTreeSet<Vec<usize>>) {
        let last = *path.last().unwrap();
        if path.len() == length {
            if is_edge(last, path[0]) {
                let mut set = path.clone();
                set.sort();
                out.insert(set);
            }
            return;
        }
        for next in 0..VERTEX_COUNT {
            if next > path[0] && is_edge(last, next) && !path.contains(&next) {
                path.push(next);
                extend(path, length, out);
                path.pop();
            }
        }
    }
    let mut out = BTreeSet::new();
    for start in 0..VERTEX_COUNT {
        extend(&mut vec![start], length, &mut out);
    }
    out
}

#[test]
fn faces_are_twenty_triangles_thirty_squares_and_twelve_pentagons() {
    let v = vertices();
    let mut faces = 0;
    for (length, count, per_vertex) in [(3, 20, 1), (4, 30, 2), (5, 12, 1)] {
        let planar: Vec<Vec<usize>> = cycles(length)
            .into_iter()
            .filter(|c| coplanar(&c.iter().map(|i| &v[*i]).collect::<Vec<_>>()))
            .collect();
        assert_eq!(planar.len(), count, "faces with {length} vertices");
        let mut membership = [0usize; VERTEX_COUNT];
        for face in &planar {
            for i in face {
                membership[*i] += 1;
            }
        }
        assert!(membership.iter().all(|m| *m == per_vertex));
        faces += planar.len();
    }
    // Euler's formula for the closed surface.
    assert_eq!(VERTEX_COUNT + faces, EDGE_COUNT + 2);
}

fn code_of(g: &Quaternion) -> SymmetryCode {
    g.clone().map(|c| {
        let a = c.rational_part() * q(4);
        let b = c.sqrt5_part() * q(4);
        let one = num_bigint::BigInt::from(1);
        assert!(a.denom() == &one && b.denom() == &one);
        [
            i64::try_from(a.numer()).unwrap(),
            i64::try_from(b.numer()).unwrap(),
        ]
    })
}

#[test]
fn symmetries_are_sixty_normalized_unit_quaternions_closed_under_products() {
    let group = symmetries();
    assert_eq!(group.len(), SYMMETRY_COUNT);
    let codes: Vec<_> = group.iter().map(code_of).collect();
    assert!(codes.windows(2).all(|pair| pair[0] < pair[1]));
    for g in group {
        assert_eq!(g.iter().fold(int(0), |s, c| &s + &(c * c)), int(1));
        let first = g.iter().find(|c| !c.is_zero()).unwrap();
        assert_eq!(first.sign(), Ordering::Greater);
    }
    assert!(symmetry_index(&[int(1), int(0), int(0), int(0)]).is_some());
    for a in group {
        assert!(symmetry_index(&conjugate(a)).is_some());
        for b in group {
            assert!(symmetry_index(&quaternion_product(a, b)).is_some());
        }
    }
}

#[test]
fn symmetry_index_accepts_either_sign_and_refuses_non_members() {
    for (i, g) in symmetries().iter().enumerate() {
        assert_eq!(symmetry_index(g), Some(i));
        assert_eq!(symmetry_index(&g.clone().map(|c| -&c)), Some(i));
        assert_eq!(symmetry_index(&g.clone().map(|c| c.scale(&q(2)))), None);
    }
    let not_a_symmetry = [rat(frac(3, 5)), rat(frac(4, 5)), int(0), int(0)];
    assert_eq!(symmetry_index(&not_a_symmetry), None);
}

#[test]
fn every_symmetry_permutes_vertices_and_edges_and_is_a_proper_rotation() {
    let lookup = vertex_lookup();
    let identity: Matrix<C> =
        std::array::from_fn(|i| std::array::from_fn(|j| int(i64::from(i == j))));
    let mut permutations = BTreeSet::new();
    for g in symmetries() {
        let m = rotation_matrix(g);
        assert_eq!(exact_product(&m, &transpose(&m)), identity);
        assert_eq!(exact_determinant(&m), int(1));
        let image: Vec<usize> = vertices()
            .iter()
            .map(|v| {
                let turned = rotate(g, v);
                assert_eq!(turned, matrix_point(&m, v));
                *lookup.get(&turned).expect("image is a vertex")
            })
            .collect();
        assert_eq!(image.iter().collect::<BTreeSet<_>>().len(), VERTEX_COUNT);
        for [a, b] in edges() {
            assert!(is_edge(image[*a], image[*b]));
        }
        permutations.insert(image);
    }
    assert_eq!(permutations.len(), SYMMETRY_COUNT);
}

#[test]
fn quaternion_rotation_numerator_matches_conjugation_and_composes() {
    let mut random = Random(17);
    for _ in 0..40 {
        let a = random.quaternion(3);
        let b = random.quaternion(3);
        let p = random.point(5);
        let norm = a.iter().fold(int(0), |s, c| &s + &(c * c));
        let m = rotation_matrix(&a);
        assert_eq!(matrix_point(&m, &p), rotate(&a, &p));
        assert_eq!(rotate(&quaternion_product(&a, &b), &p), rotate(&a, &rotate(&b, &p)));
        assert_eq!(
            rotation_matrix(&quaternion_product(&a, &b)),
            exact_product(&m, &rotation_matrix(&b))
        );
        let square = &norm * &norm;
        let expected: Matrix<C> = std::array::from_fn(|i| {
            std::array::from_fn(|j| if i == j { square.clone() } else { int(0) })
        });
        assert_eq!(exact_product(&m, &transpose(&m)), expected);
        assert_eq!(exact_determinant(&m), &square * &norm);
        // The rotation fixes its axis: |q|² R(q) v = |q|² v for q = (w, v).
        let axis: Point = [a[1].clone(), a[2].clone(), a[3].clone()];
        assert_eq!(rotate(&a, &axis), axis.map(|c| &c * &norm));
    }
    // A half-turn about the z axis negates x and y.
    let half_turn = [int(0), int(0), int(0), int(1)];
    assert_eq!(
        rotate(&half_turn, &[int(1), int(2), int(3)]),
        [int(-1), int(-2), int(3)]
    );
}

fn rotation_variables() -> [Polynomial; 3] {
    [2, 3, 4].map(var)
}

#[test]
fn cayley_matrix_is_symbolically_orthogonal_with_the_cube_of_its_scale_as_determinant() {
    let r = rotation_variables();
    let m = cayley_matrix(&r).unwrap();
    let scale = &Polynomial::constant(int(1)) + &squared_norm(&r).unwrap();
    let square = scale.mul(&scale).unwrap();
    let product = polynomial_product(&m, &transpose(&m));
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { square.clone() } else { Polynomial::zero() };
            assert_eq!(product[i][j], expected);
        }
    }
    assert_eq!(polynomial_determinant(&m), square.mul(&scale).unwrap());
    // R̂(r) r = (1 + |r|²) r: the rotation vector is the axis.
    for i in 0..3 {
        let image = (0..3).fold(Polynomial::zero(), |s, k| &s + &m[i][k].mul(&r[k]).unwrap());
        assert_eq!(image, scale.mul(&r[i]).unwrap());
    }
}

#[test]
fn cayley_matrix_equals_the_quaternion_numerator_at_points_and_for_polynomial_rotations() {
    let mut random = Random(29);
    let m = cayley_matrix(&rotation_variables()).unwrap();
    // A rotation vector that is itself a nonlinear polynomial.
    let composite = [
        var(0).mul(&var(1)).unwrap(),
        &var(2) - &Polynomial::constant(rat(frac(1, 3))),
        &var(4).scale(&int(2)) + &var(3),
    ];
    let n = cayley_matrix(&composite).unwrap();
    for _ in 0..60 {
        let x: [Q; 5] = std::array::from_fn(|_| random.rational(1));
        let point = [rat(x[2].clone()), rat(x[3].clone()), rat(x[4].clone())];
        let quaternion = [int(1), point[0].clone(), point[1].clone(), point[2].clone()];
        let numerator = rotation_matrix(&quaternion);
        let values = [
            rat(&x[0] * &x[1]),
            rat(&x[2] - frac(1, 3)),
            rat(&x[4] * q(2) + &x[3]),
        ];
        let composite_numerator =
            rotation_matrix(&[int(1), values[0].clone(), values[1].clone(), values[2].clone()]);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(at(&m[i][j], &x), numerator[i][j]);
                assert_eq!(at(&n[i][j], &x), composite_numerator[i][j]);
            }
        }
        let p = random.point(3);
        let image = polynomial_matrix_point(&m, &p);
        let turned = rotate(&quaternion, &p);
        for i in 0..3 {
            assert_eq!(at(&image[i], &x), turned[i]);
        }
    }
}

#[test]
fn cayley_matrix_refuses_exponent_overflow() {
    let huge = var(2).pow(200).unwrap();
    let r = [huge, Polynomial::zero(), Polynomial::zero()];
    assert!(cayley_matrix(&r).is_err());
    assert!(squared_norm(&r).is_err());
}

#[test]
fn polynomial_vector_helpers_agree_with_exact_arithmetic() {
    let mut random = Random(31);
    let a: [Polynomial; 3] = [
        &var(0) + &var(1).mul(&var(2)).unwrap(),
        var(3).scale(&int(-3)),
        &var(4).mul(&var(4)).unwrap() + &Polynomial::constant(int(2)),
    ];
    for _ in 0..60 {
        let x: [Q; 5] = std::array::from_fn(|_| random.rational(2));
        let exact: Point = [
            rat(&x[0] + &x[1] * &x[2]),
            rat(&x[3] * q(-3)),
            rat(&x[4] * &x[4] + q(2)),
        ];
        let p = random.point(4);
        assert_eq!(at(&polynomial_dot(&a, &p), &x), dot(&exact, &p));
        let c = polynomial_cross(&a, &p);
        let expected = cross(&exact, &p);
        for j in 0..3 {
            assert_eq!(at(&c[j], &x), expected[j]);
        }
        assert_eq!(at(&squared_norm(&a).unwrap(), &x), dot(&exact, &exact));
    }
}

#[test]
fn screen_contains_the_view_in_its_kernel_and_scales_linearly_with_the_view() {
    let mut random = Random(37);
    for _ in 0..80 {
        let s = random.number(1);
        let t = random.number(1);
        let lambda = random.number(3);
        let u = [s.clone(), t.clone(), int(1)];
        let scaled = u.clone().map(|c| &c * &lambda);
        let x = random.point(4);
        assert_eq!(screen(&u, &u), [int(0), int(0)]);
        assert_eq!(screen(&scaled, &scaled), [int(0), int(0)]);
        let affine = screen(&u, &x);
        assert_eq!(affine, [&x[0] - &(&s * &x[2]), &x[1] - &(&t * &x[2])]);
        assert_eq!(screen(&scaled, &x), affine.clone().map(|c| &c * &lambda));
        // Adding a multiple of the view does not move the image.
        let shifted: Point = std::array::from_fn(|j| &x[j] + &(&u[j] * &lambda));
        assert_eq!(screen(&u, &shifted), affine);
    }
}

#[test]
fn views_give_the_affine_and_homogeneous_view_vectors() {
    let s = var(0);
    let t = var(1).scale(&int(3));
    let affine = View::Affine {
        s: s.clone(),
        t: t.clone(),
    };
    assert_eq!(affine.vector(), [s.clone(), t.clone(), Polynomial::constant(int(1))]);
    let u = [s, t, &var(2) + &Polynomial::constant(int(1))];
    assert_eq!(View::Homogeneous(u.clone()).vector(), u);
}
