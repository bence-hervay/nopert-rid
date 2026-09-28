//! A second, deliberately simple implementation used only to cross-check the
//! production code. It reads a polynomial only through its list of terms and
//! uses plain rational arithmetic and textbook formulas:
//!
//! - Bernstein coefficients: binomial expansion of `x = lo + w·y` term by term,
//!   then `b_I = Σ_{K<=I} Π_j C(i_j,k_j)/C(n_j,k_j) · a_K`;
//! - the Bernstein form evaluated at a point, and de Casteljau subdivision;
//! - negativity proved independently by bisection and a centred interval
//!   enclosure with a rational enclosure of sqrt5;
//! - exact grids of points in a box.
use crate::arithmetic::exact::tests::sqrt5_enclosure;
use crate::arithmetic::exact::{q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::{Polynomial, VARIABLES};
use std::collections::BTreeMap;

pub type Index = [usize; VARIABLES];
pub type Tensor = BTreeMap<Index, QSqrt5>;

/// Every multi-index `K` with `K_j <= bounds_j`, lexicographically.
pub fn indices_below(bounds: &Index) -> Vec<Index> {
    let mut out = vec![[0; VARIABLES]];
    for j in 0..VARIABLES {
        let mut next = Vec::new();
        for prefix in &out {
            for k in 0..=bounds[j] {
                let mut index = *prefix;
                index[j] = k;
                next.push(index);
            }
        }
        out = next;
    }
    out
}

/// `C(n, k)` by the product formula.
pub fn binomial(n: usize, k: usize) -> Q {
    assert!(k <= n);
    let mut value = q(1);
    for i in 0..k {
        value = value * q((n - i) as i64) / q((i + 1) as i64);
    }
    value
}

fn power(x: &Q, k: usize) -> Q {
    let mut value = q(1);
    for _ in 0..k {
        value = &value * x;
    }
    value
}

/// The degree of `p` in each variable, from its terms.
pub fn degrees(p: &Polynomial) -> Index {
    let mut degrees = [0; VARIABLES];
    for (exponents, _) in p.terms() {
        for j in 0..VARIABLES {
            degrees[j] = degrees[j].max(usize::from(exponents[j]));
        }
    }
    degrees
}

/// Plain evaluation at a rational point, with powers by repeated product.
pub fn evaluate(p: &Polynomial, point: &[Q; VARIABLES]) -> QSqrt5 {
    let mut sum = QSqrt5::zero();
    for (exponents, coefficient) in p.terms() {
        let mut factor = q(1);
        for j in 0..VARIABLES {
            factor *= power(&point[j], usize::from(exponents[j]));
        }
        sum = &sum + &coefficient.scale(&factor);
    }
    sum
}

/// Coefficients of `p(origin + scale ⊙ y)` as a polynomial in `y`, by the
/// binomial theorem applied to every term.
pub fn affine_expansion(p: &Polynomial, origin: &[Q; VARIABLES], scale: &[Q; VARIABLES]) -> Tensor {
    let mut out = Tensor::new();
    for (exponents, coefficient) in p.terms() {
        let e: Index = std::array::from_fn(|j| usize::from(exponents[j]));
        for k in indices_below(&e) {
            let mut factor = q(1);
            for j in 0..VARIABLES {
                factor *= binomial(e[j], k[j])
                    * power(&origin[j], e[j] - k[j])
                    * power(&scale[j], k[j]);
            }
            let entry = out.entry(k).or_insert_with(QSqrt5::zero);
            *entry = &*entry + &coefficient.scale(&factor);
        }
    }
    out
}

/// The tensor Bernstein coefficients of `p` on `cell`, of `p`'s own degrees.
pub fn bernstein(p: &Polynomial, cell: &[Interval; VARIABLES]) -> Tensor {
    let n = degrees(p);
    let lo: [Q; VARIABLES] = std::array::from_fn(|j| cell[j].lo().clone());
    let width: [Q; VARIABLES] = std::array::from_fn(|j| cell[j].width());
    let a = affine_expansion(p, &lo, &width);
    let mut out = Tensor::new();
    for i in indices_below(&n) {
        let mut b = QSqrt5::zero();
        for k in indices_below(&i) {
            if let Some(a_k) = a.get(&k) {
                let mut weight = q(1);
                for j in 0..VARIABLES {
                    weight *= binomial(i[j], k[j]) / binomial(n[j], k[j]);
                }
                b = &b + &a_k.scale(&weight);
            }
        }
        out.insert(i, b);
    }
    out
}

/// `Σ_I b_I Π_j C(n_j,i_j) y_j^i_j (1-y_j)^(n_j-i_j)`.
pub fn bernstein_form(b: &Tensor, n: &Index, y: &[Q; VARIABLES]) -> QSqrt5 {
    let mut sum = QSqrt5::zero();
    for (i, value) in b {
        let mut basis = q(1);
        for j in 0..VARIABLES {
            let complement = q(1) - &y[j];
            basis *= binomial(n[j], i[j]) * power(&y[j], i[j]) * power(&complement, n[j] - i[j]);
        }
        sum = &sum + &value.scale(&basis);
    }
    sum
}

/// de Casteljau subdivision of a Bernstein tensor along `axis` at the
/// parameter `tau`: the coefficients on the two pieces `[0,tau]`, `[tau,1]`.
pub fn de_casteljau(b: &Tensor, n: &Index, axis: usize, tau: &Q) -> (Tensor, Tensor) {
    let (mut left, mut right) = (Tensor::new(), Tensor::new());
    let one_minus = q(1) - tau;
    let mut offset_bounds = *n;
    offset_bounds[axis] = 0;
    for base in indices_below(&offset_bounds) {
        let mut level: Vec<QSqrt5> = (0..=n[axis])
            .map(|i| {
                let mut index = base;
                index[axis] = i;
                b[&index].clone()
            })
            .collect();
        let degree = n[axis];
        let mut firsts = vec![level[0].clone()];
        let mut lasts = vec![level[degree].clone()];
        for _ in 0..degree {
            level = level
                .windows(2)
                .map(|w| &w[0].scale(&one_minus) + &w[1].scale(tau))
                .collect();
            firsts.push(level[0].clone());
            lasts.push(level[level.len() - 1].clone());
        }
        for i in 0..=degree {
            let mut index = base;
            index[axis] = i;
            left.insert(index, firsts[i].clone());
            right.insert(index, lasts[degree - i].clone());
        }
    }
    (left, right)
}

/// Every point `lo + width·(i/k)` of a `(k+1)^5` grid: corners, edges,
/// faces and interior points of the closed box.
pub fn grid(cell: &[Interval; VARIABLES], k: usize) -> Vec<[Q; VARIABLES]> {
    indices_below(&[k; VARIABLES])
        .into_iter()
        .map(|i| {
            std::array::from_fn(|j| {
                cell[j].lo() + cell[j].width() * Q::new((i[j] as i64).into(), (k as i64).into())
            })
        })
        .collect()
}

/// An upper bound of `p` on `cell` by the centred form: expand about the
/// midpoint and bound every term with interval arithmetic, using a rational
/// interval around sqrt5. Also returns, per variable, the total width of the
/// ranges of the terms involving it (a guide for where to subdivide).
pub fn centred_upper_bound(
    p: &Polynomial,
    cell: &[Interval; VARIABLES],
    root: &Interval,
) -> (Q, [Q; VARIABLES]) {
    let middle: [Q; VARIABLES] = std::array::from_fn(|j| cell[j].midpoint());
    let radius: [Q; VARIABLES] = std::array::from_fn(|j| cell[j].width() / q(2));
    let expansion = affine_expansion(p, &middle, &[q(1), q(1), q(1), q(1), q(1)]);
    let mut upper = q(0);
    let mut spread: [Q; VARIABLES] = std::array::from_fn(|_| q(0));
    for (k, coefficient) in &expansion {
        // Range of Π_j δ_j^k_j over |δ_j| <= radius_j.
        let mut magnitude = q(1);
        let mut signed = false;
        for j in 0..VARIABLES {
            magnitude *= power(&radius[j], k[j]);
            if k[j] % 2 == 1 {
                signed = true;
            }
        }
        let range = if signed {
            Interval::new(-magnitude.clone(), magnitude).unwrap()
        } else if k.iter().all(|e| *e == 0) {
            Interval::point(q(1))
        } else {
            Interval::new(q(0), magnitude).unwrap()
        };
        let value = &Interval::point(coefficient.rational_part().clone())
            + &(&Interval::point(coefficient.sqrt5_part().clone()) * root);
        let term = &value * &range;
        upper += term.hi().clone();
        for j in 0..VARIABLES {
            if k[j] > 0 {
                spread[j] += term.width();
            }
        }
    }
    (upper, spread)
}

/// Independent proof that `p < 0` on `cell`: bisect until the centred upper
/// bound is negative on every piece, splitting the axis whose terms spread
/// the bound most. `false` means only that no proof was found within
/// `budget` pieces.
pub fn proves_negative(p: &Polynomial, cell: &[Interval; VARIABLES], budget: usize) -> bool {
    let root = sqrt5_enclosure(200);
    let mut pending = vec![cell.clone()];
    let mut examined = 0;
    while let Some(piece) = pending.pop() {
        examined += 1;
        if examined > budget {
            return false;
        }
        let (upper, spread) = centred_upper_bound(p, &piece, &root);
        if upper < q(0) {
            continue;
        }
        let axis = (0..VARIABLES)
            .filter(|j| spread[*j] > q(0))
            .max_by(|a, b| spread[*a].cmp(&spread[*b]));
        match axis {
            Some(axis) => {
                for half in piece[axis].bisect() {
                    let mut next = piece.clone();
                    next[axis] = half;
                    pending.push(next);
                }
            }
            // Every range is a point: the bound is the enclosed exact value.
            None => return false,
        }
    }
    true
}
