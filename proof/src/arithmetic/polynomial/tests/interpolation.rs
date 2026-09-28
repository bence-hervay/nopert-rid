//! Checks of the Bernstein coefficients and the certifier against a second
//! derivation: coefficients are
//! recovered by *interpolation* (evaluating on a grid of `n_j + 1` points per
//! axis and inverting the Bernstein collocation matrix by rational
//! Gauss-Jordan elimination), and signs in Q(√5) are decided by escalating
//! rational enclosures of √5, never by the crate's integer sign rule. Also
//! metamorphic tests: permuted variables, positive scaling, sub-boxes.
use crate::arithmetic::exact::tests::independent::independent_sign;
use crate::arithmetic::exact::tests::rational_power;
use crate::arithmetic::exact::*;
use crate::arithmetic::polynomial::*;
use num_bigint::BigInt;
use num_traits::One;
use std::cmp::Ordering;
use std::time::Instant;

type Cell = [Interval; VARIABLES];
type Point = [Q; VARIABLES];

/// xorshift64*, a generator unrelated to the other tests' splitmix64.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo + 1) as u64) as i64
    }
    fn chance(&mut self, k: u64, n: u64) -> bool {
        self.below(n) < k
    }
    /// Rationals with awkward denominators: powers of 3, 7, 11, large primes,
    /// mixed with plain small fractions.
    fn rational(&mut self, magnitude: i64) -> Q {
        let d = match self.below(6) {
            0 => BigInt::from(3).pow(self.below(40) as u32),
            1 => BigInt::from(7).pow(self.below(25) as u32),
            2 => BigInt::from(1_000_000_007i64) * BigInt::from(self.int(1, 97)),
            3 => BigInt::one() << (self.below(90) as usize),
            _ => BigInt::from(self.int(1, 60)),
        };
        let n = BigInt::from(self.int(-magnitude, magnitude)) * &d
            + BigInt::from(self.int(-1000, 1000));
        Q::new(n, d)
    }
    fn coefficient(&mut self) -> QSqrt5 {
        match self.below(8) {
            // A near-cancelling Pell-type value conj(phi)^k scaled by +-2^e.
            0 => {
                let k = self.int(1, 60) as usize;
                let conjugate_phi = QSqrt5::new(frac(1, 2), frac(-1, 2));
                let mut x = QSqrt5::one();
                for _ in 0..k {
                    x = &x * &conjugate_phi;
                }
                let sign = if self.chance(1, 2) { q(1) } else { q(-1) };
                x.scale(&(sign * Q::from_integer(BigInt::one() << (self.below(30) as usize))))
            }
            1 | 2 => QSqrt5::from_rational(self.rational(9)),
            _ => QSqrt5::new(self.rational(9), self.rational(5)),
        }
    }
    fn polynomial(&mut self, max: [u8; VARIABLES], terms: usize) -> Polynomial {
        let mut p = Polynomial::zero();
        for _ in 0..terms {
            let e = std::array::from_fn(|j| self.below(u64::from(max[j]) + 1) as u8);
            p = &p + &Polynomial::monomial(e, self.coefficient());
        }
        p
    }
    fn interval(&mut self) -> Interval {
        match self.below(7) {
            0 => Interval::point(self.rational(3)),
            1 => {
                // Both ends with unrelated large denominators.
                let (a, b) = (self.rational(3), self.rational(3));
                Interval::new(a.clone().min(b.clone()), a.max(b)).unwrap()
            }
            2 => {
                // Straddling zero.
                let magnitude = |x: Q| if x < q(0) { -x } else { x };
                let a = magnitude(self.rational(2)) + frac(1, 9);
                let b = magnitude(self.rational(2)) + frac(1, 11);
                Interval::new(-a, b).unwrap()
            }
            3 => {
                // Tiny width at an awkward place.
                let lo = self.rational(2);
                let w = Q::new(BigInt::one(), BigInt::from(3).pow(self.int(20, 120) as u32));
                Interval::new(lo.clone(), lo + w).unwrap()
            }
            _ => {
                let lo = self.rational(1);
                Interval::new(lo.clone(), lo + frac(self.int(1, 40), 17)).unwrap()
            }
        }
    }
    fn cell(&mut self) -> Cell {
        std::array::from_fn(|_| self.interval())
    }
    fn inside(&mut self, i: &Interval) -> Q {
        match self.below(5) {
            0 => i.lo().clone(),
            1 => i.hi().clone(),
            _ => {
                let t = Q::new(BigInt::from(self.int(0, 1 << 30)), BigInt::from(1 << 30));
                i.lo() + i.width() * t
            }
        }
    }
    fn sub_interval(&mut self, i: &Interval) -> Interval {
        let (a, b) = (self.inside(i), self.inside(i));
        Interval::new(a.clone().min(b.clone()), a.max(b)).unwrap()
    }
}

fn iv(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn unit() -> Cell {
    std::array::from_fn(|_| iv(q(0), q(1)))
}

fn constant(c: QSqrt5) -> Polynomial {
    Polynomial::constant(c)
}

fn integer(n: i64) -> Polynomial {
    constant(QSqrt5::integer(n))
}

fn binomial(n: usize, k: usize) -> Q {
    (0..k).fold(q(1), |v, i| v * q((n - i) as i64) / q((i + 1) as i64))
}

/// Inverse of the collocation matrix `M[s][r] = B_{n,r}(s/n)` by Gauss-Jordan.
fn inverse_collocation(n: usize) -> Vec<Vec<Q>> {
    if n == 0 {
        return vec![vec![q(1)]];
    }
    let size = n + 1;
    let mut m: Vec<Vec<Q>> = (0..size)
        .map(|s| {
            let y = Q::new(BigInt::from(s), BigInt::from(n));
            let one_minus = q(1) - &y;
            let mut row: Vec<Q> = (0..size)
                .map(|r| {
                    binomial(n, r)
                        * rational_power(&y, r as u32)
                        * rational_power(&one_minus, (n - r) as u32)
                })
                .collect();
            row.extend((0..size).map(|c| if c == s { q(1) } else { q(0) }));
            row
        })
        .collect();
    for col in 0..size {
        let pivot = (col..size).find(|&r| !m[r][col].is_zero()).expect("invertible");
        m.swap(col, pivot);
        let inv = q(1) / &m[col][col];
        for c in 0..2 * size {
            m[col][c] = &m[col][c] * &inv;
        }
        for r in 0..size {
            if r != col && !m[r][col].is_zero() {
                let factor = m[r][col].clone();
                for c in 0..2 * size {
                    let delta = &factor * &m[col][c];
                    m[r][c] -= delta;
                }
            }
        }
    }
    m.into_iter().map(|row| row[size..].to_vec()).collect()
}

fn own_evaluate(p: &Polynomial, x: &Point) -> QSqrt5 {
    p.terms().fold(QSqrt5::zero(), |sum, (e, c)| {
        let f = (0..VARIABLES).fold(q(1), |f, j| f * rational_power(&x[j], u32::from(e[j])));
        &sum + &c.scale(&f)
    })
}

fn own_degrees(p: &Polynomial) -> [usize; VARIABLES] {
    let mut n = [0; VARIABLES];
    for (e, _) in p.terms() {
        for j in 0..VARIABLES {
            n[j] = n[j].max(usize::from(e[j]));
        }
    }
    n
}

/// All multi-indices below `n` (inclusive), variable 0 slowest.
fn indices(n: &[usize; VARIABLES]) -> Vec<[usize; VARIABLES]> {
    let mut out = vec![[0; VARIABLES]];
    for j in 0..VARIABLES {
        out = out
            .into_iter()
            .flat_map(|p| {
                (0..=n[j]).map(move |k| {
                    let mut i = p;
                    i[j] = k;
                    i
                })
            })
            .collect();
    }
    out
}

/// Bernstein coefficients by interpolation, in lexicographic order.
fn interpolated(p: &Polynomial, cell: &Cell) -> Vec<([usize; VARIABLES], QSqrt5)> {
    let n = own_degrees(p);
    let all = indices(&n);
    let position = |i: &[usize; VARIABLES]| {
        (0..VARIABLES).fold(0usize, |o, j| o * (n[j] + 1) + i[j])
    };
    let mut values: Vec<QSqrt5> = all
        .iter()
        .map(|s| {
            let x: Point = std::array::from_fn(|j| {
                if n[j] == 0 {
                    cell[j].lo().clone()
                } else {
                    cell[j].lo() + cell[j].width() * Q::new(BigInt::from(s[j]), BigInt::from(n[j]))
                }
            });
            own_evaluate(p, &x)
        })
        .collect();
    for j in 0..VARIABLES {
        if n[j] == 0 {
            continue;
        }
        let inverse = inverse_collocation(n[j]);
        let mut next = vec![QSqrt5::zero(); values.len()];
        for i in &all {
            let mut sum = QSqrt5::zero();
            for s in 0..=n[j] {
                let mut from = *i;
                from[j] = s;
                sum = &sum + &values[position(&from)].scale(&inverse[i[j]][s]);
            }
            next[position(i)] = sum;
        }
        values = next;
    }
    all.into_iter().zip(values).collect()
}

fn as_exponents(i: &[usize; VARIABLES]) -> Exponents {
    std::array::from_fn(|j| u8::try_from(i[j]).unwrap())
}

/// Verdict and first nonnegative index against the interpolation reference
/// with independent signs; also exact equality of `bernstein_coefficients`.
fn differential(p: &Polynomial, cell: &Cell) -> bool {
    let reference = interpolated(p, cell);
    let exact = p.bernstein_coefficients(cell).unwrap();
    assert_eq!(exact.len(), reference.len());
    for (i, b) in &reference {
        assert_eq!(exact.get(&as_exponents(i)), Some(b), "{p:?} on {cell:?} at {i:?}");
    }
    let first_bad = reference
        .iter()
        .find(|(_, b)| independent_sign(b) != Ordering::Less)
        .map(|(i, _)| as_exponents(i));
    match (p.certify_negative(cell), first_bad) {
        (Ok(()), None) => true,
        (Err(Refusal::NonNegativeCoefficient { index }), Some(bad)) => {
            assert_eq!(index, bad);
            false
        }
        (verdict, bad) => panic!("verdict {verdict:?} vs reference {bad:?}: {p:?} on {cell:?}"),
    }
}

/// x0 is positive on [1, 2]; the lower end spelled -1/-1 (once a way to make
/// the certifier read a negative scale) is the normal 1 by construction.
#[test]
fn a_positive_polynomial_is_refused_however_its_box_ends_are_spelled() {
    let [x, ..] = Polynomial::variables();
    let one = Q::new(BigInt::from(-1), BigInt::from(-1));
    let mut cell = unit();
    cell[0] = iv(one.clone(), q(2));
    assert!(x.certify_negative(&cell).is_err());
    let p = &x - &constant(QSqrt5::from_rational(frac(1, 2)));
    cell[0] = iv(one, Q::new(BigInt::from(-3), BigInt::from(-2)));
    let coefficients = p.bernstein_coefficients(&cell).unwrap();
    assert!(coefficients.iter().all(|(_, b)| b.sign() == Ordering::Greater));
    assert!(p.certify_negative(&cell).is_err());
    assert!((-&p).certify_negative(&cell).is_ok());
}

#[test]
fn differential_against_interpolation_on_awkward_boxes() {
    let mut rng = Rng(2024);
    let (mut yes, mut no) = (0, 0);
    for case in 0..500 {
        let bound = if case % 5 == 0 { 5 } else { 3 };
        let max: [u8; VARIABLES] = std::array::from_fn(|_| rng.below(bound) as u8);
        let terms = 1 + rng.below(7) as usize;
        let mut p = rng.polynomial(max, terms);
        let cell = rng.cell();
        // Shift so that the value at a random point of the box is slightly
        // negative: mixed verdicts close to the boundary.
        let point: Point = std::array::from_fn(|j| rng.inside(&cell[j]));
        let shift = &own_evaluate(&p, &point) + &QSqrt5::from_rational(frac(rng.int(0, 40), 32));
        p = &p - &constant(shift);
        if differential(&p, &cell) {
            yes += 1;
        } else {
            no += 1;
        }
    }
    println!("certified {yes}, refused {no}");
    assert!(yes > 40 && no > 40, "{yes} {no}");
}

#[test]
fn verdict_is_invariant_under_variable_permutation() {
    let mut rng = Rng(31);
    let variables = Polynomial::variables();
    let permutations = [[1, 2, 3, 4, 0], [4, 3, 2, 1, 0], [0, 2, 1, 4, 3], [3, 0, 4, 2, 1]];
    for _ in 0..150 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial(max, 6);
        let cell = rng.cell();
        let point: Point = std::array::from_fn(|j| rng.inside(&cell[j]));
        let margin = QSqrt5::from_rational(frac(rng.int(0, 8), 4));
        let p = &p - &constant(&own_evaluate(&p, &point) + &margin);
        let base = p.bernstein_coefficients(&cell).unwrap();
        for sigma in permutations {
            // pq(x) = p(x_sigma(0), ..., x_sigma(4)) on the box with
            // cell'[sigma(j)] = cell[j].
            let images: [Polynomial; VARIABLES] =
                std::array::from_fn(|j| variables[sigma[j]].clone());
            let pq = p.substitute(&images).unwrap();
            let mut moved = cell.clone();
            for j in 0..VARIABLES {
                moved[sigma[j]] = cell[j].clone();
            }
            assert_eq!(pq.certify_negative(&moved).is_ok(), p.certify_negative(&cell).is_ok());
            let permuted = pq.bernstein_coefficients(&moved).unwrap();
            for (i, b) in base.iter() {
                let mut k = [0u8; VARIABLES];
                for j in 0..VARIABLES {
                    k[sigma[j]] = i[j];
                }
                assert_eq!(permuted.get(&k), Some(b));
            }
        }
    }
}

#[test]
fn verdict_is_invariant_under_positive_algebraic_scaling() {
    let mut rng = Rng(32);
    let conjugate_phi = QSqrt5::new(frac(1, 2), frac(-1, 2));
    let phi = QSqrt5::new(frac(1, 2), frac(1, 2));
    for _ in 0..150 {
        let max = std::array::from_fn(|_| rng.below(3) as u8);
        let p = &rng.polynomial(max, 5) - &integer(rng.int(0, 6));
        let cell = rng.cell();
        let verdict = p.certify_negative(&cell);
        let k = rng.int(1, 120);
        let factor = if rng.chance(1, 2) {
            let sign = q(if k % 2 == 0 { 1 } else { -1 });
            (0..k).fold(QSqrt5::one(), |x, _| &x * &conjugate_phi).scale(&sign)
        } else {
            (0..k).fold(QSqrt5::one(), |x, _| &x * &phi)
        };
        assert_eq!(independent_sign(&factor), Ordering::Greater);
        assert_eq!(p.scale(&factor).certify_negative(&cell), verdict);
        let negative = -&factor;
        if verdict.is_ok() {
            assert!(p.scale(&negative).certify_negative(&cell).is_err());
        }
    }
}

/// Coefficients on any sub-box are convex combinations of the parent's, so a
/// certified box certifies every sub-box with arbitrary rational ends.
#[test]
fn certified_boxes_certify_every_random_sub_box() {
    let mut rng = Rng(33);
    let mut checked = 0;
    for _ in 0..400 {
        let max = std::array::from_fn(|_| rng.below(3) as u8);
        let p = rng.polynomial(max, 5);
        let cell = rng.cell();
        let point: Point = std::array::from_fn(|j| rng.inside(&cell[j]));
        let margin = QSqrt5::from_rational(frac(rng.int(1, 16), 4));
        let p = &p - &constant(&own_evaluate(&p, &point) + &margin);
        if p.certify_negative(&cell).is_err() {
            continue;
        }
        for _ in 0..8 {
            let sub: Cell = std::array::from_fn(|j| rng.sub_interval(&cell[j]));
            assert!(p.certify_negative(&sub).is_ok(), "{p:?} on {sub:?} inside {cell:?}");
        }
        checked += 1;
    }
    assert!(checked > 50, "{checked}");
}

/// A nonnegative value anywhere in the closed box forces a refusal.
#[test]
fn any_nonnegative_sample_forces_refusal() {
    let mut rng = Rng(34);
    for _ in 0..300 {
        let max = std::array::from_fn(|_| rng.below(4) as u8);
        let p = rng.polynomial(max, 6);
        let cell = rng.cell();
        // Make the value at one sampled point exactly 0 or a tiny positive.
        let point: Point = std::array::from_fn(|j| rng.inside(&cell[j]));
        let bump = if rng.chance(1, 2) {
            QSqrt5::zero()
        } else {
            QSqrt5::from_rational(Q::new(BigInt::one(), BigInt::one() << 400usize))
        };
        let p = &(&p - &constant(own_evaluate(&p, &point))) + &constant(bump);
        assert!(p.certify_negative(&cell).is_err(), "{p:?} on {cell:?}, value >= 0 at {point:?}");
    }
}

#[test]
fn degree_255_axes() {
    let [x, y, ..] = Polynomial::variables();
    let top = x.pow(255).unwrap();
    let mut cell = unit();
    // -x^255 - 1 on [0,1]: certified.
    assert!((&(-&top) - &integer(1)).certify_negative(&cell).is_ok());
    // x^255 - 1 on [0,1]: zero at x = 1, refused at the last index.
    assert_eq!(
        (&top - &integer(1)).certify_negative(&cell),
        Err(Refusal::NonNegativeCoefficient { index: [255, 0, 0, 0, 0] })
    );
    // On [-1, 1] the coefficients of x^255 are (-1)^(255-r): x^255 - 1 has
    // zero coefficients at odd r; the first is r = 1.
    cell[0] = iv(q(-1), q(1));
    let b = top.bernstein_coefficients(&cell).unwrap();
    for (i, v) in b.iter() {
        let expected = if (255 - i[0]) % 2 == 0 { 1 } else { -1 };
        assert_eq!(v, &QSqrt5::integer(expected));
    }
    assert_eq!(
        (&top - &integer(1)).certify_negative(&cell),
        Err(Refusal::NonNegativeCoefficient { index: [1, 0, 0, 0, 0] })
    );
    assert!((&top - &integer(2)).certify_negative(&cell).is_ok());
    // A 256 x 256 tensor, exactly at the size limit: x^255 y^255 - 1 on
    // [0, 1 - 2^-8]^2.
    assert_eq!(256 * 256, MAX_BERNSTEIN_COEFFICIENTS);
    let almost = q(1) - frac(1, 256);
    cell[0] = iv(q(0), almost.clone());
    cell[1] = iv(q(0), almost);
    let start = Instant::now();
    let p = &top.mul(&y.pow(255).unwrap()).unwrap() - &integer(1);
    assert!(p.certify_negative(&cell).is_ok());
    cell[1] = iv(q(0), q(1));
    cell[0] = iv(q(0), q(1));
    assert_eq!(
        p.certify_negative(&cell),
        Err(Refusal::NonNegativeCoefficient { index: [255, 255, 0, 0, 0] })
    );
    println!("256x256 tensor: {:?}", start.elapsed());
}

/// Zero maximum attained only at an irrational interior point.
#[test]
fn zero_maximum_at_irrational_interior_point_is_refused() {
    let [x, y, z, ..] = Polynomial::variables();
    let five = integer(5);
    let well = -&(&x.pow(2).unwrap() - &five).pow(2).unwrap();
    let positive = &(&integer(1) + &y.pow(2).unwrap()) + &z.pow(4).unwrap();
    let p = well.mul(&positive).unwrap();
    let mut cell = unit();
    cell[0] = iv(q(2), q(3));
    cell[1] = iv(q(-1), q(2));
    assert!(p.certify_negative(&cell).is_err());
    // Plus a positive constant far below any rational sample's resolution.
    let tiny = constant(QSqrt5::from_rational(Q::new(BigInt::one(), BigInt::one() << 1000usize)));
    assert!((&p + &tiny).certify_negative(&cell).is_err());
    // phi-type zero: x^2 - x - 1 = 0 at x = phi; -(x^2-x-1)^2 - 2^-1000 is negative
    // everywhere but its maximum is -2^-1000; certify refuses or accepts, but
    // then every exact sample must be negative (already implied).
    let golden = -&(&(&x.pow(2).unwrap() - &x) - &integer(1)).pow(2).unwrap();
    cell[0] = iv(frac(8, 5), frac(17, 10));
    assert!(golden.certify_negative(&cell).is_err());
}

#[test]
fn substitute_and_divide_round_trip_on_zoom_like_maps() {
    let mut rng = Rng(35);
    let [s, t, x, y, z] = Polynomial::variables();
    for _ in 0..60 {
        let g = rng.polynomial([1, 1, 2, 2, 2], 10);
        // Point blow-up in r1: (s, t, r1, r2, r3) -> (s, t, r1, r1 r2, r1 r3).
        let images = [s.clone(), t.clone(), x.clone(), x.mul(&y).unwrap(), x.mul(&z).unwrap()];
        let pulled = g.substitute(&images).unwrap();
        // Remove the part of g not vanishing at r = 0 so the pull-back is
        // divisible by r1.
        let origin = [s.clone(), t.clone(), integer(0), integer(0), integer(0)];
        let at_zero = g.substitute(&origin).unwrap();
        let vanishing = &g - &at_zero;
        let pulled_vanishing = vanishing.substitute(&images).unwrap();
        let quotient = pulled_vanishing.divide_by_monomial(&[0, 0, 1, 0, 0]).unwrap();
        assert_eq!(quotient.mul(&x).unwrap(), pulled_vanishing);
        // Pointwise.
        for _ in 0..3 {
            let pt: Point = std::array::from_fn(|_| rng.rational(2));
            let mapped: Point = [
                pt[0].clone(),
                pt[1].clone(),
                pt[2].clone(),
                &pt[2] * &pt[3],
                &pt[2] * &pt[4],
            ];
            assert_eq!(own_evaluate(&pulled, &pt), own_evaluate(&g, &mapped));
        }
        if !at_zero.is_zero() {
            assert!(pulled.divide_by_monomial(&[0, 0, 1, 0, 0]).is_err());
        }
    }
}

// ============================================================== resources

/// Cost of the certifier for large tensors that pass the size guard.
/// cargo test --release --lib
///     arithmetic::polynomial::tests::interpolation::measure_large_tensors -- --ignored --nocapture
#[test]
#[ignore]
fn measure_large_tensors() {
    let mut rng = Rng(36);
    let cell: Cell = [
        iv(frac(123, 3 * 512), frac(124, 3 * 512)),
        iv(frac(77, 5 * 512), frac(78, 5 * 512)),
        iv(frac(-201, 5 * 512), frac(-200, 5 * 512)),
        iv(frac(3, 5 * 512), frac(4, 5 * 512)),
        iv(frac(-9, 5 * 512), frac(-8, 5 * 512)),
    ];
    for degrees in [[7u8, 7, 7, 0, 0], [15, 15, 15, 0, 0], [31, 31, 31, 0, 0], [39, 39, 39, 0, 0]] {
        let mut p = rng.polynomial(degrees, 40);
        p = &p + &Polynomial::monomial(degrees, QSqrt5::one());
        let size: usize = degrees.iter().map(|d| usize::from(*d) + 1).product();
        let start = Instant::now();
        let verdict = p.certify_negative(&cell);
        let elapsed = start.elapsed();
        println!("degrees {degrees:?}: {size} coefficients, {elapsed:?}, {:?}", verdict.is_ok());
    }
}

/// Near the worst case the size limit admits: 2^16 coefficients with the
/// largest per-axis work, degrees (255, 255, 0, 0, 0), and for comparison
/// (15, 15, 15, 15, 0); 41 terms with awkward coefficient denominators (up to
/// 3^39, 7^24, 2^89) on a box whose ends have 11- to 13-bit denominators.
/// cargo test --release --lib
///     arithmetic::polynomial::tests::interpolation::measure_at_cap -- --ignored --nocapture
#[test]
#[ignore]
fn measure_at_cap() {
    let mut rng = Rng(37);
    let cell: Cell = [
        iv(frac(123, 3 * 512), frac(124, 3 * 512)),
        iv(frac(77, 5 * 512), frac(78, 5 * 512)),
        iv(frac(-201, 5 * 512), frac(-200, 5 * 512)),
        iv(frac(3, 5 * 512), frac(4, 5 * 512)),
        iv(frac(-9, 5 * 512), frac(-8, 5 * 512)),
    ];
    for degrees in [[255u8, 255, 0, 0, 0], [15, 15, 15, 15, 0]] {
        let p = &rng.polynomial(degrees, 40) + &Polynomial::monomial(degrees, QSqrt5::one());
        let size: usize = degrees.iter().map(|d| usize::from(*d) + 1).product();
        assert_eq!(size, MAX_BERNSTEIN_COEFFICIENTS);
        let start = Instant::now();
        let verdict = p.certify_negative(&cell);
        println!("degrees {degrees:?}: {:?}, {:?}", start.elapsed(), verdict.is_ok());
        let start = Instant::now();
        p.bernstein_coefficients(&cell).unwrap();
        println!("  exact coefficients: {:?}", start.elapsed());
    }
}
