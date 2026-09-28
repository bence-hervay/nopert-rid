//! Tests of the polynomial layer. `reference` is an independent, deliberately
//! simple second implementation used to cross-check results that could be
//! subtly wrong (Bernstein coefficients and every certified verdict).
use super::*;
use crate::arithmetic::exact::tests::{rational_power, sqrt5_enclosure};
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use num_bigint::BigInt;
use num_traits::One;
use std::cmp::Ordering;
use std::time::Instant;

mod interpolation;
mod reference;

type Cell = [Interval; VARIABLES];
type Point = [Q; VARIABLES];

/// splitmix64: deterministic, well mixed.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    /// Uniform in `lo..=hi`.
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo + 1) as u64) as i64
    }
    fn chance(&mut self, numerator: u64, denominator: u64) -> bool {
        self.below(denominator) < numerator
    }
    fn rational(&mut self, magnitude: i64) -> Q {
        let denominators = [1, 2, 3, 4, 5, 7, 8, 16, 60];
        let d = denominators[self.below(denominators.len() as u64) as usize];
        frac(self.int(-magnitude, magnitude), d)
    }
    fn coefficient(&mut self) -> QSqrt5 {
        let b = if self.chance(1, 2) { self.rational(6) } else { q(0) };
        QSqrt5::new(self.rational(9), b)
    }
    /// Up to `terms` random terms with exponents bounded by `max`.
    fn polynomial(&mut self, max: [u8; VARIABLES], terms: usize) -> Polynomial {
        let mut p = Polynomial::zero();
        for _ in 0..terms {
            let exponents = std::array::from_fn(|j| self.below(u64::from(max[j]) + 1) as u8);
            p = &p + &Polynomial::monomial(exponents, self.coefficient());
        }
        p
    }
    /// A box with a random scale: most axes of moderate size, some of zero
    /// width, some tiny, some huge and far from the origin.
    fn cell(&mut self) -> Cell {
        std::array::from_fn(|_| {
            let lo = self.rational(4);
            let kind = self.below(10);
            let width = match kind {
                0 => q(0),
                1 => Q::new(BigInt::one(), BigInt::one() << (40 + self.below(80) as usize)),
                2 => q(self.int(1, 1 << 20)),
                _ => frac(self.int(1, 32), 16),
            };
            let offset = if kind == 2 { q(self.int(-(1 << 30), 1 << 30)) } else { q(0) };
            let lo = lo + offset;
            Interval::new(lo.clone(), lo + width).unwrap()
        })
    }
    /// A sub-box of the root search box reached by `depth` random bisections.
    fn search_cell(&mut self, depth: usize) -> Cell {
        let mut cell = root_cell();
        for level in 0..depth {
            let halves = cell[level % VARIABLES].bisect();
            cell[level % VARIABLES] = halves[self.below(2) as usize].clone();
        }
        cell
    }
    /// A sub-box of the root search box at a random depth below `max_depth`.
    fn random_search_cell(&mut self, max_depth: u64) -> Cell {
        let depth = self.below(max_depth) as usize;
        self.search_cell(depth)
    }
    /// Between 1 and `max_terms` random terms.
    fn polynomial_up_to(&mut self, max: [u8; VARIABLES], max_terms: u64) -> Polynomial {
        let terms = 1 + self.below(max_terms) as usize;
        self.polynomial(max, terms)
    }
    /// A point of the closed box; each coordinate is a random end with
    /// probability 1/3, otherwise a random interior rational.
    fn point_in(&mut self, cell: &Cell) -> Point {
        std::array::from_fn(|j| match self.below(6) {
            0 => cell[j].lo().clone(),
            1 => cell[j].hi().clone(),
            _ => {
                let t = frac(self.int(1, 999), 1000);
                cell[j].lo() + cell[j].width() * t
            }
        })
    }
    fn point(&mut self) -> Point {
        std::array::from_fn(|_| self.rational(5))
    }
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn root_cell() -> Cell {
    [
        interval(q(0), frac(2, 3)),
        interval(q(0), frac(2, 5)),
        interval(frac(-2, 5), frac(2, 5)),
        interval(frac(-2, 5), frac(2, 5)),
        interval(frac(-2, 5), frac(2, 5)),
    ]
}

fn unit_cell() -> Cell {
    std::array::from_fn(|_| interval(q(0), q(1)))
}

fn algebraic(point: &Point) -> [QSqrt5; VARIABLES] {
    std::array::from_fn(|j| QSqrt5::from_rational(point[j].clone()))
}

fn value(p: &Polynomial, point: &Point) -> QSqrt5 {
    p.evaluate(&algebraic(point))
}

fn constant(c: QSqrt5) -> Polynomial {
    Polynomial::constant(c)
}

fn integer(n: i64) -> Polynomial {
    constant(QSqrt5::integer(n))
}

fn qs(a: Q, b: Q) -> QSqrt5 {
    QSqrt5::new(a, b)
}

fn certified(p: &Polynomial, cell: &Cell) -> bool {
    p.certify_negative(cell).is_ok()
}

fn exponents_of(index: &reference::Index) -> Exponents {
    std::array::from_fn(|j| u8::try_from(index[j]).unwrap())
}

/// Checks everything that must hold for a certified verdict and the verdict
/// itself against the independent implementation: the verdict equals
/// "every reference coefficient is negative", and a certified polynomial is
/// negative on an exact grid and at random points of the box, and is proved
/// negative independently by subdivision.
fn cross_check(p: &Polynomial, cell: &Cell, rng: &mut Rng, grid: usize, budget: usize) -> bool {
    let verdict = p.certify_negative(cell);
    let coefficients = reference::bernstein(p, cell);
    let all_negative = coefficients.values().all(|b| b.sign() == Ordering::Less);
    assert_eq!(verdict.is_ok(), all_negative, "{p:?} on {cell:?}");
    match &verdict {
        Err(Refusal::NonNegativeCoefficient { index }) => {
            let first = coefficients
                .iter()
                .find(|(_, b)| b.sign() != Ordering::Less)
                .map(|(i, _)| exponents_of(i))
                .unwrap();
            assert_eq!(*index, first);
        }
        Err(Refusal::TooManyCoefficients { .. }) => panic!("unexpected size refusal"),
        Ok(()) => {
            for point in reference::grid(cell, grid) {
                assert_eq!(value(p, &point).sign(), Ordering::Less, "{p:?} at {point:?}");
            }
            for _ in 0..20 {
                let point = rng.point_in(cell);
                assert_eq!(reference::evaluate(p, &point).sign(), Ordering::Less);
            }
            assert!(
                reference::proves_negative(p, cell, budget),
                "no independent proof for {p:?} on {cell:?}"
            );
        }
    }
    verdict.is_ok()
}

// ---------------------------------------------------------------- arithmetic

#[test]
fn ring_operations_agree_with_pointwise_values() {
    let mut rng = Rng(1);
    for _ in 0..150 {
        let p = rng.polynomial([3, 2, 2, 1, 2], 8);
        let r = rng.polynomial([2, 2, 1, 3, 2], 8);
        let c = rng.coefficient();
        let point = rng.point();
        let (pv, rv) = (reference::evaluate(&p, &point), reference::evaluate(&r, &point));
        assert_eq!(value(&p, &point), pv);
        assert_eq!(value(&(&p + &r), &point), &pv + &rv);
        assert_eq!(value(&(&p - &r), &point), &pv - &rv);
        assert_eq!(value(&(-&p), &point), -&pv);
        assert_eq!(value(&p.mul(&r).unwrap(), &point), &pv * &rv);
        assert_eq!(value(&p.scale(&c), &point), &pv * &c);
        assert_eq!(value(&p.pow(3).unwrap(), &point), &(&pv * &pv) * &pv);
    }
}

#[test]
fn evaluation_at_algebraic_points_matches_binomial_expansion() {
    // For a polynomial with rational coefficients, expand each power of
    // a + b·sqrt5 by the binomial theorem independently.
    let mut rng = Rng(2);
    for _ in 0..80 {
        let mut p = Polynomial::zero();
        for _ in 0..6 {
            let exponents = std::array::from_fn(|_| rng.below(4) as u8);
            p = &p + &Polynomial::monomial(exponents, QSqrt5::from_rational(rng.rational(9)));
        }
        let point: [QSqrt5; VARIABLES] =
            std::array::from_fn(|_| qs(rng.rational(3), rng.rational(3)));
        let power = |x: &QSqrt5, k: usize| -> (Q, Q) {
            let (a, b) = (x.rational_part(), x.sqrt5_part());
            let (mut rational, mut radical) = (q(0), q(0));
            for i in 0..=k {
                let term = reference::binomial(k, i)
                    * rational_power(a, (k - i) as u32)
                    * rational_power(b, i as u32);
                let five = rational_power(&q(5), (i / 2) as u32);
                if i % 2 == 0 {
                    rational += term * five;
                } else {
                    radical += term * five;
                }
            }
            (rational, radical)
        };
        let mut expected = QSqrt5::zero();
        for (exponents, coefficient) in p.terms() {
            let mut term = coefficient.clone();
            for j in 0..VARIABLES {
                let (a, b) = power(&point[j], usize::from(exponents[j]));
                term = &term * &qs(a, b);
            }
            expected = &expected + &term;
        }
        assert_eq!(p.evaluate(&point), expected);
        // Galois symmetry for rational coefficients.
        let conjugate: [QSqrt5; VARIABLES] = std::array::from_fn(|j| point[j].conjugate());
        assert_eq!(p.evaluate(&conjugate), expected.conjugate());
    }
}

#[test]
fn ring_identities_hold_exactly() {
    let mut rng = Rng(3);
    let one = integer(1);
    for _ in 0..60 {
        let [a, b, c] = [0, 1, 2].map(|_| rng.polynomial([2, 1, 2, 1, 1], 6));
        assert_eq!(a.mul(&b).unwrap(), b.mul(&a).unwrap());
        assert_eq!(
            a.mul(&b).unwrap().mul(&c).unwrap(),
            a.mul(&b.mul(&c).unwrap()).unwrap()
        );
        assert_eq!(
            a.mul(&(&b + &c)).unwrap(),
            &a.mul(&b).unwrap() + &a.mul(&c).unwrap()
        );
        assert_eq!(&(&a + &b) + &c, &a + &(&b + &c));
        assert_eq!(&a - &b, &a + &(-&b));
        assert!((&a - &a).is_zero());
        assert_eq!(a.mul(&one).unwrap(), a);
        assert!(a.mul(&Polynomial::zero()).unwrap().is_zero());
        assert_eq!(&a + &Polynomial::zero(), a);
        assert_eq!(a.pow(2).unwrap(), a.mul(&a).unwrap());
        assert_eq!(a.pow(0).unwrap(), one);
        assert!(a.terms().all(|(_, c)| !c.is_zero()));
    }
}

#[test]
fn no_term_is_truncated_and_cancellation_is_exact() {
    let [x, y, z, v, w] = Polynomial::variables();
    let cubic = &(&(&integer(1) + &x) + &x.pow(2).unwrap()) + &x.pow(3).unwrap();
    let squared = cubic.mul(&cubic).unwrap();
    for (degree, expected) in [1, 2, 3, 4, 3, 2, 1].into_iter().enumerate() {
        assert_eq!(
            squared.coefficient(&[degree as u8, 0, 0, 0, 0]),
            QSqrt5::integer(expected)
        );
    }
    assert_eq!(squared.terms().count(), 7);
    // (1 + x + y + z + v + w)^4 has every one of the C(9,5) = 126 monomials of
    // total degree at most four, with multinomial coefficients.
    let sum = [&x, &y, &z, &v, &w].into_iter().fold(integer(1), |s, t| &s + t);
    let fourth = sum.pow(4).unwrap();
    assert_eq!(fourth.terms().count(), 126);
    assert_eq!(fourth.coefficient(&[1, 1, 1, 1, 0]), QSqrt5::integer(24));
    assert_eq!(fourth.coefficient(&[0, 0, 0, 0, 0]), QSqrt5::integer(1));
    assert_eq!(fourth.coefficient(&[0, 0, 0, 0, 4]), QSqrt5::integer(1));
    // (x + sqrt5)(x - sqrt5) = x² - 5: the linear sqrt5 terms cancel exactly.
    let root = constant(QSqrt5::sqrt5());
    let product = (&x + &root).mul(&(&x - &root)).unwrap();
    assert_eq!(product, &x.pow(2).unwrap() - &integer(5));
    assert_eq!(product.terms().count(), 2);
    assert!(Polynomial::monomial([1, 2, 3, 4, 5], QSqrt5::zero()).is_zero());
    assert_eq!(Polynomial::zero().degrees(), [0; VARIABLES]);
    assert_eq!(product.coefficient(&[1, 0, 0, 0, 0]), QSqrt5::zero());
    assert_eq!(
        fourth.mul(&v.pow(3).unwrap()).unwrap().degrees(),
        [4, 4, 4, 7, 4]
    );
}

#[test]
fn exponent_overflow_is_a_typed_error() {
    let [x, _, _, v, _] = Polynomial::variables();
    let top = x.pow(255).unwrap();
    assert_eq!(top.degrees(), [255, 0, 0, 0, 0]);
    assert_eq!(top.mul(&x), Err(PolynomialError::ExponentOverflow { variable: 0 }));
    let square = v.mul(&v).unwrap();
    assert_eq!(square.pow(128), Err(PolynomialError::ExponentOverflow { variable: 3 }));
    assert_eq!(square.pow(127).unwrap().degrees(), [0, 0, 0, 254, 0]);
    // Constants never overflow; the u8 exponent bounds their work.
    let three = Polynomial::constant(QSqrt5::integer(3));
    let power = three.pow(u8::MAX).unwrap();
    assert_eq!(power, Polynomial::constant(QSqrt5::from_rational(rational_power(&q(3), 255))));
    let images = [x.pow(2).unwrap(), x.clone(), x.clone(), x.clone(), x.clone()];
    assert_eq!(
        x.pow(128).unwrap().substitute(&images),
        Err(PolynomialError::ExponentOverflow { variable: 0 })
    );
    assert!(x.pow(127).unwrap().substitute(&images).is_ok());
}

// -------------------------------------------------------------- substitution

#[test]
fn substitution_is_pointwise_composition() {
    let mut rng = Rng(4);
    for _ in 0..60 {
        let p = rng.polynomial([2, 2, 2, 1, 1], 6);
        let images: [Polynomial; VARIABLES] =
            std::array::from_fn(|_| rng.polynomial([1, 1, 2, 1, 1], 3));
        let composed = p.substitute(&images).unwrap();
        for _ in 0..4 {
            let point = rng.point();
            let inner: [QSqrt5; VARIABLES] = std::array::from_fn(|j| value(&images[j], &point));
            assert_eq!(value(&composed, &point), p.evaluate(&inner));
        }
    }
}

#[test]
fn substitution_identities() {
    let mut rng = Rng(5);
    let variables = Polynomial::variables();
    for _ in 0..30 {
        let p = rng.polynomial([2, 2, 2, 2, 2], 8);
        assert_eq!(p.substitute(&variables).unwrap(), p);
        let point = rng.point();
        let constants: [Polynomial; VARIABLES] =
            std::array::from_fn(|j| constant(QSqrt5::from_rational(point[j].clone())));
        assert_eq!(p.substitute(&constants).unwrap(), constant(value(&p, &point)));
        assert!(Polynomial::zero().substitute(&constants).unwrap().is_zero());
        // Associativity with zoom-like images: a variable times an affine
        // function, or a quadratic in few variables.
        let small = rng.polynomial([2, 1, 1, 0, 0], 4);
        let zoom = |rng: &mut Rng| -> [Polynomial; VARIABLES] {
            std::array::from_fn(|j| {
                let affine = &rng.polynomial([1, 1, 0, 0, 1], 1) + &integer(1);
                &variables[j].mul(&affine).unwrap() + &rng.polynomial([0, 0, 2, 1, 0], 1)
            })
        };
        let (f, g) = (zoom(&mut rng), zoom(&mut rng));
        let f_then_g: [Polynomial; VARIABLES] =
            std::array::from_fn(|j| f[j].substitute(&g).unwrap());
        assert_eq!(
            small.substitute(&f).unwrap().substitute(&g).unwrap(),
            small.substitute(&f_then_g).unwrap()
        );
    }
}

// ------------------------------------------------------------ monomial division

#[test]
fn monomial_division_inverts_multiplication_exactly() {
    let mut rng = Rng(6);
    for _ in 0..100 {
        let quotient = rng.polynomial([3, 3, 3, 3, 3], 8);
        let factor: Exponents = std::array::from_fn(|_| rng.below(4) as u8);
        let product = quotient
            .mul(&Polynomial::monomial(factor, QSqrt5::one()))
            .unwrap();
        assert_eq!(product.divide_by_monomial(&factor).unwrap(), quotient);
        let point = rng.point();
        let monomial_value =
            (0..VARIABLES).fold(q(1), |a, j| a * rational_power(&point[j], u32::from(factor[j])));
        assert_eq!(
            value(&product, &point),
            value(&quotient, &point).scale(&monomial_value)
        );
        assert_eq!(product.divide_by_monomial(&[0; VARIABLES]).unwrap(), product);
    }
}

#[test]
fn monomial_division_refuses_every_non_multiple() {
    let mut rng = Rng(7);
    for _ in 0..100 {
        let quotient = rng.polynomial([2, 2, 2, 2, 2], 6);
        let mut factor: Exponents = std::array::from_fn(|_| rng.below(3) as u8);
        let variable = rng.below(VARIABLES as u64) as usize;
        factor[variable] = factor[variable].max(1);
        let product = quotient
            .mul(&Polynomial::monomial(factor, QSqrt5::one()))
            .unwrap();
        // One extra term lacking a single power of one variable.
        let mut short = factor;
        short[variable] -= 1;
        let corrupted = &product + &Polynomial::monomial(short, rng.coefficient());
        if corrupted.coefficient(&short).is_zero() {
            continue;
        }
        let first_bad = *corrupted
            .terms()
            .map(|(e, _)| e)
            .find(|e| (0..VARIABLES).any(|j| e[j] < factor[j]))
            .unwrap();
        assert_eq!(
            corrupted.divide_by_monomial(&factor),
            Err(PolynomialError::NotDivisible {
                term: first_bad,
                factor
            })
        );
    }
    let [x, y, ..] = Polynomial::variables();
    assert!(Polynomial::zero().divide_by_monomial(&[9, 9, 9, 9, 9]).unwrap().is_zero());
    assert_eq!(x.divide_by_monomial(&[1, 0, 0, 0, 0]).unwrap(), integer(1));
    assert_eq!(
        x.mul(&y).unwrap().divide_by_monomial(&[1, 0, 1, 0, 0]),
        Err(PolynomialError::NotDivisible {
            term: [1, 1, 0, 0, 0],
            factor: [1, 0, 1, 0, 0]
        })
    );
    assert!(integer(3).divide_by_monomial(&[0, 0, 0, 0, 1]).is_err());
}

// ---------------------------------------------------- Bernstein coefficients

fn as_tensor(b: &BernsteinCoefficients) -> reference::Tensor {
    b.iter()
        .map(|(index, value)| (std::array::from_fn(|j| usize::from(index[j])), value.clone()))
        .collect()
}

#[test]
fn bernstein_coefficients_match_the_reference_implementation() {
    let mut rng = Rng(8);
    for case in 0..250 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial_up_to(max, 8);
        let cell = if case % 3 == 0 {
            rng.random_search_cell(30)
        } else {
            rng.cell()
        };
        let b = p.bernstein_coefficients(&cell).unwrap();
        assert_eq!(b.degrees(), p.degrees());
        assert_eq!(
            b.len(),
            p.degrees().iter().map(|n| usize::from(*n) + 1).product::<usize>()
        );
        assert_eq!(as_tensor(&b), reference::bernstein(&p, &cell), "{p:?} on {cell:?}");
    }
}

#[test]
fn bernstein_form_reproduces_the_polynomial_exactly() {
    let mut rng = Rng(9);
    for _ in 0..120 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial(max, 6);
        let cell = rng.cell();
        let b = as_tensor(&p.bernstein_coefficients(&cell).unwrap());
        let n = reference::degrees(&p);
        for sample in 0..6 {
            // Inside or outside the box: the identity is one of polynomials.
            let x: Point = if sample < 4 { rng.point_in(&cell) } else { rng.point() };
            let y: Point = std::array::from_fn(|j| {
                if cell[j].width() == q(0) {
                    // Every coefficient along a zero-width axis is the same.
                    rng.rational(3)
                } else {
                    (&x[j] - cell[j].lo()) / cell[j].width()
                }
            });
            // Along a zero-width axis the form is the value at its point.
            let projected: Point = std::array::from_fn(|j| {
                if cell[j].width() == q(0) { cell[j].lo().clone() } else { x[j].clone() }
            });
            assert_eq!(reference::bernstein_form(&b, &n, &y), reference::evaluate(&p, &projected));
        }
    }
}

#[test]
fn corner_coefficients_are_corner_values_and_bound_all_values() {
    let mut rng = Rng(10);
    for _ in 0..120 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial(max, 6);
        let cell = rng.cell();
        let b = p.bernstein_coefficients(&cell).unwrap();
        let n = p.degrees();
        for corner in 0..(1 << VARIABLES) {
            let upper = |j: usize| corner & (1 << j) != 0;
            let index: Exponents = std::array::from_fn(|j| if upper(j) { n[j] } else { 0 });
            let point: Point = std::array::from_fn(|j| {
                if upper(j) { cell[j].hi().clone() } else { cell[j].lo().clone() }
            });
            assert_eq!(b.get(&index).unwrap(), &value(&p, &point));
        }
        let lowest = b.iter().map(|(_, v)| v).min().unwrap();
        let highest = b.iter().map(|(_, v)| v).max().unwrap();
        for _ in 0..10 {
            let v = value(&p, &rng.point_in(&cell));
            assert!(lowest <= &v && &v <= highest);
        }
    }
}

#[test]
fn de_casteljau_subdivision_matches_direct_coefficients() {
    let mut rng = Rng(11);
    for case in 0..120 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial(max, 6);
        let cell = if case % 2 == 0 {
            rng.random_search_cell(20)
        } else {
            rng.cell()
        };
        let n = reference::degrees(&p);
        let parent = as_tensor(&p.bernstein_coefficients(&cell).unwrap());
        let axis = rng.below(VARIABLES as u64) as usize;
        // The midpoint, then an arbitrary rational split parameter.
        for tau in [frac(1, 2), frac(rng.int(1, 12), 13)] {
            let split = cell[axis].lo() + cell[axis].width() * &tau;
            let (mut left_cell, mut right_cell) = (cell.clone(), cell.clone());
            left_cell[axis] = interval(cell[axis].lo().clone(), split.clone());
            right_cell[axis] = interval(split, cell[axis].hi().clone());
            if tau == frac(1, 2) {
                let halves = [left_cell[axis].clone(), right_cell[axis].clone()];
                assert_eq!(halves, cell[axis].bisect());
            }
            let (left, right) = reference::de_casteljau(&parent, &n, axis, &tau);
            assert_eq!(as_tensor(&p.bernstein_coefficients(&left_cell).unwrap()), left);
            assert_eq!(as_tensor(&p.bernstein_coefficients(&right_cell).unwrap()), right);
        }
    }
}

#[test]
fn bernstein_coefficients_are_invariant_under_affine_reparametrisation() {
    // p(a ⊙ x + c) on the preimage box has the same coefficients, with the
    // index reversed along every axis where a is negative.
    let mut rng = Rng(12);
    for _ in 0..80 {
        let max = std::array::from_fn(|_| rng.below(3) as u8);
        let p = rng.polynomial(max, 6);
        let cell = rng.cell();
        let a: Point = std::array::from_fn(|_| {
            let magnitude = frac(rng.int(1, 9), rng.int(1, 9));
            if rng.chance(1, 2) { -magnitude } else { magnitude }
        });
        let c: Point = std::array::from_fn(|_| rng.rational(5));
        let variables = Polynomial::variables();
        let images: [Polynomial; VARIABLES] = std::array::from_fn(|j| {
            &variables[j].scale(&QSqrt5::from_rational(a[j].clone()))
                + &constant(QSqrt5::from_rational(c[j].clone()))
        });
        let pulled = p.substitute(&images).unwrap();
        let preimage: Cell = std::array::from_fn(|j| {
            let ends = [(cell[j].lo() - &c[j]) / &a[j], (cell[j].hi() - &c[j]) / &a[j]];
            let [u, v] = ends;
            interval(u.clone().min(v.clone()), u.max(v))
        });
        let original = p.bernstein_coefficients(&cell).unwrap();
        let transformed = pulled.bernstein_coefficients(&preimage).unwrap();
        assert_eq!(original.degrees(), transformed.degrees());
        let n = original.degrees();
        for (index, b) in original.iter() {
            let mirrored: Exponents = std::array::from_fn(|j| {
                if a[j] < q(0) { n[j] - index[j] } else { index[j] }
            });
            assert_eq!(transformed.get(&mirrored).unwrap(), b);
        }
    }
}

#[test]
fn constant_and_multilinear_cases() {
    let mut rng = Rng(13);
    for _ in 0..40 {
        let c = rng.coefficient();
        let cell = rng.cell();
        let b = constant(c.clone()).bernstein_coefficients(&cell).unwrap();
        assert_eq!(b.len(), 1);
        assert_eq!(b.get(&[0; VARIABLES]), Some(&c));
        assert_eq!(b.get(&[0, 0, 1, 0, 0]), None);
        // Multilinear: every coefficient is the value at a corner.
        let p = rng.polynomial([1, 1, 1, 1, 1], 12);
        let b = p.bernstein_coefficients(&cell).unwrap();
        let n = p.degrees();
        for (index, v) in b.iter() {
            let corner: Point = std::array::from_fn(|j| {
                if n[j] == 1 && index[j] == 1 { cell[j].hi().clone() } else { cell[j].lo().clone() }
            });
            assert_eq!(v, &value(&p, &corner));
        }
        let listed: Vec<Exponents> = b.iter().map(|(i, _)| i).collect();
        assert!(listed.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(listed.len(), b.len());
    }
    let zero = Polynomial::zero().bernstein_coefficients(&unit_cell()).unwrap();
    assert!(zero.get(&[0; VARIABLES]).unwrap().is_zero());
}

#[test]
fn tensor_size_limit_is_a_typed_error_before_any_work() {
    let [x, y, z, v, w] = Polynomial::variables();
    let product = [&x, &y, &z, &v, &w]
        .into_iter()
        .fold(integer(1), |s, t| s.mul(&t.pow(16).unwrap()).unwrap());
    let degrees = [16; VARIABLES];
    assert_eq!(
        product.bernstein_coefficients(&unit_cell()),
        Err(PolynomialError::TooManyBernsteinCoefficients { degrees })
    );
    assert_eq!(
        product.certify_negative(&unit_cell()),
        Err(Refusal::TooManyCoefficients { degrees })
    );
    let wide = [&x, &y, &z, &v, &w]
        .into_iter()
        .fold(integer(1), |s, t| s.mul(&t.pow(40).unwrap()).unwrap());
    // 41^5 ≈ 1.2e8 coefficients would need more than 5 GB before the first
    // sign could be read: the refusal comes from the size alone.
    assert_eq!(
        wide.certify_negative(&unit_cell()),
        Err(Refusal::TooManyCoefficients { degrees: [40; VARIABLES] })
    );
    // A single high degree is fine.
    let tall = &(-&x.pow(60).unwrap()) - &integer(1);
    assert!(certified(&tall, &unit_cell()));
    assert_eq!(MAX_BERNSTEIN_COEFFICIENTS, 1 << 16);
    // One factor above the limit is refused before any work. (Exactly at the
    // limit, a 256 x 256 tensor is certified in `interpolation::degree_255_axes`.)
    let square = x.pow(255).unwrap().mul(&y.pow(255).unwrap()).unwrap();
    assert_eq!(
        square.mul(&z).unwrap().certify_negative(&unit_cell()),
        Err(Refusal::TooManyCoefficients { degrees: [255, 255, 1, 0, 0] })
    );
    let quartic = [&x, &y, &z, &v]
        .into_iter()
        .fold(integer(1), |s, t| s.mul(&t.pow(15).unwrap()).unwrap());
    assert_eq!(
        quartic.mul(&w).unwrap().bernstein_coefficients(&unit_cell()),
        Err(PolynomialError::TooManyBernsteinCoefficients { degrees: [15, 15, 15, 15, 1] })
    );
}

// --------------------------------------------------------------- certifier

#[test]
fn random_certified_verdicts_are_independently_confirmed() {
    let mut rng = Rng(14);
    let (mut accepted, mut refused) = (0, 0);
    for case in 0..500 {
        let cell = if case % 2 == 0 {
            rng.random_search_cell(40)
        } else {
            std::array::from_fn(|_| {
                let lo = rng.rational(2);
                interval(lo.clone(), lo + frac(rng.int(0, 8), 8))
            })
        };
        let p = centred_random(&mut rng, &cell);
        if cross_check(&p, &cell, &mut rng, 2, 4000) {
            // Coefficients on a piece are convex combinations of those on
            // the whole box, so both halves of any bisection are certified.
            let axis = rng.below(VARIABLES as u64) as usize;
            for half in cell[axis].bisect() {
                let mut piece = cell.clone();
                piece[axis] = half;
                assert!(certified(&p, &piece));
            }
            accepted += 1;
        } else {
            refused += 1;
        }
    }
    println!("accepted {accepted}, refused {refused}");
    assert!(accepted > 100 && refused > 100, "{accepted} {refused}");
}

/// A random polynomial shifted so that its value at the midpoint of `cell`
/// is a random number in [-2, 1/8]: verdicts are mixed and many are close to
/// the boundary between certifiable and not.
fn centred_random(rng: &mut Rng, cell: &Cell) -> Polynomial {
    let max = std::array::from_fn(|_| rng.below(3) as u8);
    let noise = rng.polynomial_up_to(max, 6);
    let middle: Point = std::array::from_fn(|j| cell[j].midpoint());
    let target = QSqrt5::from_rational(frac(-rng.int(-2, 32), 16));
    &(&noise - &constant(value(&noise, &middle))) + &constant(target)
}

/// `-Σ_j a_j (x_j - c_j)²`: nonpositive, zero only at `c`.
fn bowl(centre: &Point, weights: [i64; VARIABLES]) -> Polynomial {
    let variables = Polynomial::variables();
    let mut p = Polynomial::zero();
    for j in 0..VARIABLES {
        let shifted = &variables[j] - &constant(QSqrt5::from_rational(centre[j].clone()));
        p = &p - &shifted.pow(2).unwrap().scale(&QSqrt5::integer(weights[j]));
    }
    p
}

/// The points of the box with each coordinate at its lower end, upper end
/// or a fixed interior value: all corners, edge, face and interior points.
fn structured_points(cell: &Cell) -> Vec<Point> {
    reference::indices_below(&[2; VARIABLES])
        .into_iter()
        .map(|choice| {
            std::array::from_fn(|j| match choice[j] {
                0 => cell[j].lo().clone(),
                1 => cell[j].hi().clone(),
                _ => cell[j].lo() + cell[j].width() * frac(3, 7),
            })
        })
        .collect()
}

#[test]
fn refuses_zero_or_tiny_positive_values_at_interior_face_edge_and_corner_points() {
    let mut rng = Rng(15);
    // Tiny positive values: a power of two and a Pell-type number
    // conj(phi)^80 = (L_80 - F_80·sqrt5)/2, about 2^-55, with parts near 2^55.
    let conjugate_phi = qs(frac(1, 2), frac(-1, 2));
    let pell = (0..80).fold(QSqrt5::one(), |x, _| &x * &conjugate_phi);
    assert_eq!(pell.sign(), Ordering::Greater);
    let bumps = [
        QSqrt5::from_rational(Q::new(BigInt::one(), BigInt::one() << 300usize)),
        pell,
    ];
    for cell in [unit_cell(), root_cell(), rng.search_cell(33), rng.search_cell(7)] {
        let mut centres = structured_points(&cell);
        // Interior points with large denominators, off every grid.
        for _ in 0..20 {
            centres.push(std::array::from_fn(|j| {
                let t = Q::new(BigInt::from(rng.int(1, (1 << 40) - 1)), BigInt::one() << 40usize);
                cell[j].lo() + cell[j].width() * t
            }));
        }
        for centre in centres {
            let weights = std::array::from_fn(|_| rng.int(1, 5));
            let zero_at_centre = bowl(&centre, weights);
            assert_eq!(value(&zero_at_centre, &centre).sign(), Ordering::Equal);
            assert!(!certified(&zero_at_centre, &cell));
            for bump in &bumps {
                assert!(!certified(&(&zero_at_centre + &constant(bump.clone())), &cell));
            }
            // The same shape strictly below zero is certified where the
            // reference agrees, and then confirmed.
            let lowered = &zero_at_centre - &integer(1);
            cross_check(&lowered, &cell, &mut rng, 1, 2000);
        }
    }
}

#[test]
fn refuses_polynomials_vanishing_on_a_face_or_an_edge() {
    let mut rng = Rng(16);
    let variables = Polynomial::variables();
    for _ in 0..40 {
        let cell = rng.random_search_cell(35);
        let j = rng.below(VARIABLES as u64) as usize;
        let k = (j + 1 + rng.below(4) as usize) % VARIABLES;
        let lower = &variables[j] - &constant(QSqrt5::from_rational(cell[j].lo().clone()));
        let upper = &constant(QSqrt5::from_rational(cell[k].hi().clone())) - &variables[k];
        // Nonpositive on the box; zero on the face x_j = lo_j.
        let positive = &integer(1) - &bowl(&rng.point(), [1; VARIABLES]);
        let face = (-&lower).mul(&positive).unwrap();
        assert!(!certified(&face, &cell));
        // Nonpositive; zero on the edge x_j = lo_j, x_k = hi_k (and faces).
        let edge = (-&lower.pow(2).unwrap()).mul(&upper).unwrap();
        assert!(!certified(&edge, &cell));
        let edge_only = &(-&lower.pow(2).unwrap()) - &upper.pow(2).unwrap();
        assert!(!certified(&edge_only, &cell));
        // Lowered by a positive constant they are strictly negative.
        cross_check(&(&face - &integer(1)), &cell, &mut rng, 1, 1000);
        cross_check(&(&edge_only - &integer(1)), &cell, &mut rng, 1, 1000);
    }
}

#[test]
fn refuses_zero_sets_through_irrational_points_only() {
    let [x, y, ..] = Polynomial::variables();
    let root = constant(QSqrt5::sqrt5());
    // -(x - sqrt5·y)² vanishes on an irrational line through the box but at
    // no rational point of it except where x = y = 0.
    let line = -&(&x - &root.mul(&y).unwrap()).pow(2).unwrap();
    let cell: Cell = [
        interval(q(2), q(3)),
        interval(q(1), q(2)),
        interval(q(0), q(0)),
        interval(q(0), q(0)),
        interval(q(0), q(0)),
    ];
    assert!(!certified(&line, &cell));
    let mut rng = Rng(17);
    cross_check(&(&line - &integer(1)), &cell, &mut rng, 2, 4000);
    // 5 - x² changes sign at sqrt5; boxes one 2^-200 step from it.
    let parabola = &integer(5) - &x.pow(2).unwrap();
    let enclosure = sqrt5_enclosure(200);
    let step = enclosure.width();
    let at = |lo: &Q, hi: &Q| -> Cell {
        let mut cell = unit_cell();
        cell[0] = interval(lo.clone(), hi.clone());
        cell
    };
    // Contains sqrt5: refused.
    assert!(!certified(&parabola, &at(enclosure.lo(), enclosure.hi())));
    // Just above sqrt5: negative, certified and confirmed.
    let above = at(enclosure.hi(), &(enclosure.hi() + &step));
    assert!(certified(&parabola, &above));
    cross_check(&parabola, &above, &mut rng, 3, 1000);
    // Just below: positive, refused; negated, certified.
    let below = at(&(enclosure.lo() - &step), enclosure.lo());
    assert!(!certified(&parabola, &below));
    assert!(certified(&(-&parabola), &below));
    assert!(!certified(&(-&parabola), &above));
}

#[test]
fn pell_type_near_cancellation_decides_exactly() {
    // conj(phi)^k = (L_k - F_k sqrt5)/2 has absolute value phi^-k and sign
    // (-1)^k, with rational and sqrt5 parts nearly cancelling.
    let [x, ..] = Polynomial::variables();
    let conjugate_phi = qs(frac(1, 2), frac(-1, 2));
    let mut power = QSqrt5::one();
    let mut rng = Rng(18);
    for k in 1..=160 {
        power = &power * &conjugate_phi;
        assert_eq!(power.sign(), if k % 2 == 0 { Ordering::Greater } else { Ordering::Less });
        let negative = if k % 2 == 0 { -&power } else { power.clone() };
        // A tiny negative constant, alone or times a positive polynomial.
        let shapes = [
            constant(negative.clone()),
            constant(negative.clone()).mul(&(&integer(1) + &x.pow(2).unwrap())).unwrap(),
        ];
        for p in shapes {
            for cell in [unit_cell(), rng.search_cell(k % 40)] {
                assert!(certified(&p, &cell), "k = {k}");
                assert!(!certified(&(-&p), &cell), "k = {k}");
                // A + B·sqrt5 is tiny with A, B of opposite signs; flipping
                // one part leaves a value of about 2A (flipped sqrt5 part)
                // or -2A (flipped rational part), so the sign of A decides.
                let rational_negative = negative.rational_part() < &q(0);
                let flip = |rational: i64, radical: i64| {
                    p.terms().fold(Polynomial::zero(), |s, (e, c)| {
                        let c = qs(c.rational_part() * q(rational), c.sqrt5_part() * q(radical));
                        &s + &Polynomial::monomial(*e, c)
                    })
                };
                assert_eq!(certified(&flip(1, -1), &cell), rational_negative, "k = {k}");
                assert_eq!(certified(&flip(-1, 1), &cell), !rational_negative, "k = {k}");
            }
        }
        // The same near-cancellation as the largest value of a non-constant
        // polynomial: x + b·sqrt5 on [a - 1, a] and on the point a has
        // maximum a + b·sqrt5 = conj(phi)^k.
        let (a, b) = (power.rational_part().clone(), power.sqrt5_part().clone());
        let linear = &x + &constant(qs(q(0), b));
        for axis in [interval(&a - q(1), a.clone()), Interval::point(a.clone())] {
            let mut cell = unit_cell();
            cell[0] = axis;
            assert_eq!(certified(&linear, &cell), power.sign() == Ordering::Less, "k = {k}");
        }
    }
}

#[test]
fn degree_zero_variables_ignore_their_intervals() {
    let mut rng = Rng(19);
    for _ in 0..60 {
        let p = &rng.polynomial([2, 0, 2, 0, 1], 5) - &integer(rng.int(0, 6));
        let cell = rng.random_search_cell(25);
        let verdict = p.certify_negative(&cell);
        let coefficients = p.bernstein_coefficients(&cell).unwrap();
        for replacement in [
            Interval::point(q(7)),
            interval(q(-(1 << 40)), q(1 << 40)),
            interval(frac(1, 3), frac(1, 3) + Q::new(BigInt::one(), BigInt::one() << 200usize)),
        ] {
            let mut other = cell.clone();
            other[1] = replacement.clone();
            other[3] = replacement;
            assert_eq!(p.certify_negative(&other), verdict);
            assert_eq!(p.bernstein_coefficients(&other).unwrap(), coefficients);
        }
        cross_check(&p, &cell, &mut rng, 2, 2000);
    }
}

#[test]
fn zero_width_boxes_are_decided_on_the_face_they_describe() {
    let [x, y, ..] = Polynomial::variables();
    let third = constant(QSqrt5::from_rational(frac(1, 3)));
    let p = &(-&(&x - &third).pow(2).unwrap()) - &y.pow(2).unwrap();
    let mut cell = unit_cell();
    cell[0] = Interval::point(frac(1, 3));
    cell[1] = Interval::point(q(0));
    assert_eq!(
        p.certify_negative(&cell),
        Err(Refusal::NonNegativeCoefficient { index: [0; VARIABLES] })
    );
    cell[1] = interval(q(0), frac(1, 5));
    assert!(!certified(&p, &cell));
    cell[1] = interval(frac(1, 10), frac(1, 5));
    assert!(certified(&p, &cell));
    cell[0] = Interval::point(frac(1, 2));
    cell[1] = Interval::point(q(0));
    assert!(certified(&p, &cell));
    let all_points: Cell = std::array::from_fn(|_| Interval::point(frac(1, 7)));
    assert!(certified(&p, &all_points));
    let mut rng = Rng(20);
    for _ in 0..60 {
        let p = &rng.polynomial([2, 2, 2, 1, 1], 6) - &integer(rng.int(0, 4));
        let point = rng.point();
        let degenerate: Cell = std::array::from_fn(|j| Interval::point(point[j].clone()));
        // On a single point the certifier decides the sign of the value.
        assert_eq!(
            certified(&p, &degenerate),
            value(&p, &point).sign() == Ordering::Less
        );
        let mut mixed = rng.cell();
        mixed[rng.below(VARIABLES as u64) as usize] = Interval::point(rng.rational(3));
        cross_check(&p, &mixed, &mut rng, 2, 2000);
    }
}

#[test]
fn boxes_of_widely_different_scales() {
    let mut rng = Rng(21);
    let huge = |bits: usize| Q::from_integer(BigInt::one() << bits);
    let small = |bits: usize| Q::new(BigInt::one(), BigInt::one() << bits);
    for _ in 0..40 {
        let p = &rng.polynomial([2, 2, 1, 2, 1], 6) - &integer(rng.int(0, 5));
        let centre = huge(100) * rng.rational(1);
        let cell: Cell = [
            interval(small(120), small(120) * q(2)),
            interval(-huge(90), huge(90)),
            interval(centre.clone(), centre + small(60)),
            interval(frac(1, 3), frac(1, 3) + small(200)),
            Interval::point(small(500)),
        ];
        cross_check(&p, &cell, &mut rng, 1, 4000);
    }
    // x·y - 1 with x tiny and y huge: the product is exactly representable.
    let [x, y, ..] = Polynomial::variables();
    let p = &x.mul(&y).unwrap() - &integer(1);
    let mut cell = unit_cell();
    cell[0] = interval(small(200), small(200) * q(2));
    cell[1] = interval(huge(199), huge(199) * frac(3, 2));
    // x·y ranges over [1/2, 3/2]: not negative.
    assert!(!certified(&p, &cell));
    // x·y at most 2^-199·(2^199 - 1) = 1 - 2^-199: negative.
    cell[1] = interval(huge(198), huge(199) - q(1));
    assert!(certified(&p, &cell));
    // x·y reaches 1 at a corner: refused.
    cell[1] = interval(huge(198), huge(199));
    assert!(!certified(&p, &cell));
}

#[test]
fn mutations_change_the_verdict_exactly_where_they_should() {
    let mut rng = Rng(22);
    let (mut changed, mut unchanged) = (0, 0);
    for _ in 0..60 {
        let max = std::array::from_fn(|_| rng.below(3) as u8);
        let p = &rng.polynomial(max, 5).scale(&QSqrt5::from_rational(frac(1, 16)))
            - &integer(rng.int(1, 4));
        let cell = rng.random_search_cell(20);
        let original = cross_check(&p, &cell, &mut rng, 1, 2000);
        for (exponents, coefficient) in p.terms() {
            for mutated_coefficient in [
                -coefficient,
                qs(coefficient.rational_part().clone(), -coefficient.sqrt5_part()),
                qs(-coefficient.rational_part(), coefficient.sqrt5_part().clone()),
            ] {
                if &mutated_coefficient == coefficient {
                    continue;
                }
                let mutated = &(&p - &Polynomial::monomial(*exponents, coefficient.clone()))
                    + &Polynomial::monomial(*exponents, mutated_coefficient);
                let verdict = cross_check(&mutated, &cell, &mut rng, 1, 2000);
                if verdict != original {
                    changed += 1;
                } else {
                    unchanged += 1;
                }
                // A nonnegative value anywhere forces a refusal.
                if structured_points(&cell)
                    .iter()
                    .any(|x| value(&mutated, x).sign() != Ordering::Less)
                {
                    assert!(!verdict);
                }
            }
        }
    }
    println!("changed {changed}, unchanged {unchanged}");
    assert!(changed > 20 && unchanged > 20, "{changed} {unchanged}");
    // Targeted: -1 - x² is certified on [0,1]; flipping the x² sign gives a
    // zero at x = 1; flipping the constant gives positive values.
    let [x, ..] = Polynomial::variables();
    let square = x.pow(2).unwrap();
    assert!(certified(&(&integer(-1) - &square), &unit_cell()));
    assert!(!certified(&(&integer(-1) + &square), &unit_cell()));
    assert!(!certified(&(&integer(1) - &square), &unit_cell()));
}

#[test]
fn zero_polynomial_and_constants() {
    let mut rng = Rng(23);
    for _ in 0..30 {
        let cell = rng.cell();
        assert_eq!(
            Polynomial::zero().certify_negative(&cell),
            Err(Refusal::NonNegativeCoefficient { index: [0; VARIABLES] })
        );
        let c = rng.coefficient();
        assert_eq!(certified(&constant(c.clone()), &cell), c.sign() == Ordering::Less);
    }
}

/// A support-gap-shaped polynomial `(1+|r|²) n·h - n·N(r)p` with
/// `n = u × e`, `u = (s,t,1)`, `N(r) = (1-|r|²)I + 2rrᵀ + 2[r]×`: degree one in
/// s and in t, two in each rotation coordinate.
fn support_gap(edge: &[QSqrt5; 3], hole: &[QSqrt5; 3], plug: &[QSqrt5; 3]) -> Polynomial {
    let [s, t, x, y, z] = Polynomial::variables();
    let one = integer(1);
    let u = [s, t, one.clone()];
    let e: [Polynomial; 3] = std::array::from_fn(|j| constant(edge[j].clone()));
    let n: [Polynomial; 3] = std::array::from_fn(|i| {
        let (k, l) = ((i + 1) % 3, (i + 2) % 3);
        &u[k].mul(&e[l]).unwrap() - &u[l].mul(&e[k]).unwrap()
    });
    let r = [x, y, z];
    let rr = r.iter().fold(Polynomial::zero(), |sum, c| &sum + &c.pow(2).unwrap());
    let dot = |v: &[Polynomial; 3], w: &[QSqrt5; 3]| {
        (0..3).fold(Polynomial::zero(), |sum, j| &sum + &v[j].scale(&w[j]))
    };
    let rp = dot(&r, plug);
    let two = QSqrt5::integer(2);
    let rotated: [Polynomial; 3] = std::array::from_fn(|i| {
        let (k, l) = ((i + 1) % 3, (i + 2) % 3);
        let cross = &r[k].scale(&plug[l]) - &r[l].scale(&plug[k]);
        &(&(&one - &rr).scale(&plug[i]) + &r[i].mul(&rp).unwrap().scale(&two)) + &cross.scale(&two)
    });
    let mut gap = (&one + &rr).mul(&dot(&n, hole)).unwrap();
    for i in 0..3 {
        gap = &gap - &n[i].mul(&rotated[i]).unwrap();
    }
    gap
}

/// A few vertices of the edge-length-2 solid, coordinates `(a + b·sqrt5)/2`.
fn sample_vertices() -> Vec<[QSqrt5; 3]> {
    let half = |a: i64, b: i64| qs(frac(a, 2), frac(b, 2));
    let bases = [[(2, 0), (2, 0), (4, 2)], [(3, 1), (1, 1), (2, 2)], [(5, 1), (0, 0), (3, 1)]];
    let mut out = Vec::new();
    for base in bases {
        for signs in [[1, 1, 1], [-1, 1, 1], [1, -1, -1], [-1, -1, 1]] {
            out.push(std::array::from_fn(|j| half(signs[j] * base[j].0, signs[j] * base[j].1)));
        }
    }
    out
}

fn support_gap_cases(rng: &mut Rng, count: usize) -> Vec<(Polynomial, Cell)> {
    let vertices = sample_vertices();
    (0..count)
        .map(|_| {
            let pick = |rng: &mut Rng| vertices[rng.below(vertices.len() as u64) as usize].clone();
            let (a, b) = (pick(rng), pick(rng));
            let edge: [QSqrt5; 3] = std::array::from_fn(|j| &a[j] - &b[j]);
            let p = support_gap(&edge, &pick(rng), &pick(rng));
            let depth = 5 + rng.below(40) as usize;
            (p, rng.search_cell(depth))
        })
        .collect()
}

#[test]
fn support_gap_shaped_polynomials_on_search_boxes() {
    let mut rng = Rng(24);
    let (mut accepted, mut refused) = (0, 0);
    for (p, cell) in support_gap_cases(&mut rng, 150) {
        let degrees = p.degrees();
        assert!(degrees[0] <= 1 && degrees[1] <= 1 && degrees[2..].iter().all(|d| *d <= 2));
        if cross_check(&p, &cell, &mut rng, 1, 3000) {
            accepted += 1;
        } else {
            refused += 1;
        }
    }
    println!("accepted {accepted}, refused {refused}");
    assert!(accepted > 10 && refused > 10, "{accepted} {refused}");
}

/// One share of the heavy randomised campaign: the checks of the default
/// campaigns with more cases, finer grids and a larger subdivision budget.
fn heavy_random_campaign(seed: u64) {
    let mut rng = Rng(seed);
    let start = Instant::now();
    let (mut accepted, mut refused) = (0, 0);
    for case in 0..6000 {
        let cell = match case % 3 {
            0 => rng.random_search_cell(50),
            1 => rng.cell(),
            _ => std::array::from_fn(|_| {
                let lo = rng.rational(2);
                interval(lo.clone(), lo + frac(rng.int(0, 8), 8))
            }),
        };
        let p = centred_random(&mut rng, &cell);
        if cross_check(&p, &cell, &mut rng, 2, 20000) {
            accepted += 1;
        } else {
            refused += 1;
        }
    }
    for (p, cell) in support_gap_cases(&mut rng, 1500) {
        if cross_check(&p, &cell, &mut rng, 2, 20000) {
            accepted += 1;
        } else {
            refused += 1;
        }
    }
    println!(
        "seed {seed}: accepted {accepted}, refused {refused}, {:.1} s",
        start.elapsed().as_secs_f64()
    );
}

/// cargo test --release --lib arithmetic::polynomial::tests::heavy_random_campaign --
///     --ignored --nocapture --test-threads=3
#[test]
#[ignore]
fn heavy_random_campaign_1() {
    heavy_random_campaign(101);
}

#[test]
#[ignore]
fn heavy_random_campaign_2() {
    heavy_random_campaign(202);
}

#[test]
#[ignore]
fn heavy_random_campaign_3() {
    heavy_random_campaign(303);
}

/// cargo test --release --lib arithmetic::polynomial::tests::measure_certifier --
///     --ignored --nocapture
#[test]
#[ignore]
fn measure_certifier() {
    let mut rng = Rng(26);
    let vertices = sample_vertices();
    let gap = (0..vertices.len())
        .flat_map(|a| (0..vertices.len()).map(move |b| (a, b)))
        .map(|(a, b)| {
            let edge: [QSqrt5; 3] = std::array::from_fn(|j| &vertices[a][j] - &vertices[b][j]);
            support_gap(&edge, &vertices[(a + 3) % 12], &vertices[(b + 5) % 12])
        })
        .find(|p| p.degrees() == [1, 1, 2, 2, 2])
        .unwrap();
    println!("support gap: {} terms, degrees {:?}", gap.terms().count(), gap.degrees());
    for depth in [10, 25, 45, 60] {
        let cells: Vec<Cell> = (0..4000).map(|_| rng.search_cell(depth)).collect();
        let start = Instant::now();
        let accepted = cells.iter().filter(|cell| certified(&gap, cell)).count();
        let micros = start.elapsed().as_secs_f64() * 1e6 / cells.len() as f64;
        let total = cells.len();
        println!("depth {depth}: {micros:.1} us per certification, {accepted}/{total} certified");
    }
    let cells: Vec<Cell> = (0..2000).map(|_| rng.search_cell(45)).collect();
    let start = Instant::now();
    for cell in &cells {
        gap.bernstein_coefficients(cell).unwrap();
    }
    println!(
        "exact coefficients of the support gap: {:.1} us per box",
        start.elapsed().as_secs_f64() * 1e6 / cells.len() as f64
    );
    // Pulled back along a point blow-up in the first rotation coordinate.
    let [s, t, x, y, z] = Polynomial::variables();
    let images = [
        s.mul(&x).unwrap(),
        t.mul(&x).unwrap(),
        x.pow(2).unwrap(),
        x.mul(&y).unwrap(),
        x.mul(&z).unwrap(),
    ];
    let pulled = gap.substitute(&images).unwrap();
    // Dense polynomials of higher degree.
    let dense = |degree: u8, rng: &mut Rng| {
        let mut p = rng.polynomial([degree; VARIABLES], 60);
        for index in reference::indices_below(&[usize::from(degree); VARIABLES]) {
            p = &p + &Polynomial::monomial(exponents_of(&index), QSqrt5::integer(-1));
        }
        p
    };
    let cases = [
        ("pulled back", pulled),
        ("dense degree 3", dense(3, &mut rng)),
        ("dense degree 4", dense(4, &mut rng)),
    ];
    for (name, p) in cases {
        let cells: Vec<Cell> = (0..100).map(|_| rng.search_cell(40)).collect();
        let start = Instant::now();
        for cell in &cells {
            let _ = p.certify_negative(cell);
        }
        println!(
            "{name}: {} terms, degrees {:?}, {:.1} us per certification",
            p.terms().count(),
            p.degrees(),
            start.elapsed().as_secs_f64() * 1e6 / cells.len() as f64
        );
    }
}
