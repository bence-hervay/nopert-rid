//! Tests of the lemma on one cell: every refusal path, and the soundness of
//! accepted cells at exact points against an independent evaluation of the
//! witness (direct projection for gaps, the affine domain polynomials).
use super::*;
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5};
use crate::arithmetic::polynomial::Exponents;
use crate::elimination::zoom::{Five, ZoomCell, ZoomFactor};
use crate::elimination::zoom::cover::ZoomCover;
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::{tree, Coordinates, Sign};
use crate::elimination::proof::catalogue::{exotic_data, local_data};
use crate::elimination::proof::format::{Leaf, CoverFile};
use crate::elimination::witness::{Edge, Gap, MaximumGap, Support, WitnessError};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::{domain, geometry};

/// A witness leaf of a cover file: its zoom, cell, witness and factor.
struct Fixture {
    cover: ZoomCover,
    zoom: usize,
    cell: ZoomCell,
    witness: Witness,
    factor: ZoomFactor,
}

/// Every `step`-th witness leaf of a cover file.
fn fixtures(bytes: &[u8], step: usize) -> Vec<Fixture> {
    let file = CoverFile::parse(bytes).unwrap();
    let mut out = Vec::new();
    let mut count = 0;
    {
        let cover = ZoomCover::new(file.parameters.clone()).unwrap();
        for (c, data) in file.zooms.iter().enumerate() {
            let cells = tree::leaves(cover.zooms()[c].root(), &data.tree).unwrap();
            for (cell, leaf) in cells.into_iter().zip(&data.leaves) {
                if let Leaf::Witness { index, factor } = leaf {
                    count += 1;
                    if count % step == 0 {
                        out.push(Fixture {
                            cover: cover.clone(),
                            zoom: c,
                            cell,
                            witness: file.witnesses[*index].clone(),
                            factor: ZoomFactor::new(*factor),
                        });
                    }
                }
            }
        }
    }
    out
}

fn exotic(name: &str, step: usize) -> Vec<Fixture> {
    fixtures(exotic_data(name).unwrap(), step)
}

impl Fixture {
    fn zoom(&self) -> &Zoom {
        &self.cover.zooms()[self.zoom]
    }
}

fn exact(x: &Five<Q>) -> Five<QSqrt5> {
    std::array::from_fn(|j| QSqrt5::from_rational(x[j].clone()))
}

/// The homogeneous view and rotation of the zoom at `y`, computed from the
/// coordinates `x(y)` by the coordinate formulas.
fn configuration(zoom: &Zoom, y: &Five<QSqrt5>) -> (Point, Point) {
    let x = zoom.evaluate(y);
    match zoom.coordinates() {
        Coordinates::Configuration => (
            [x[0].clone(), x[1].clone(), QSqrt5::one()],
            [x[2].clone(), x[3].clone(), x[4].clone()],
        ),
        Coordinates::ArcPlane(sign) => {
            let two = QSqrt5::integer(2);
            let u = [&QSqrt5::one() - &(&two * &x[0]), x[1].clone(), &two + &x[0]];
            let (a, w) = (crate::elimination::zoom::arc_rotation(), crate::elimination::zoom::arc_direction(sign));
            let sigma = QSqrt5::integer(sign.value());
            let r = [
                &(&(&sigma * &a[0]) + &(&x[2] * &w[0])) + &x[3],
                &x[2] * &w[1],
                &(&(&sigma * &a[2]) + &(&x[2] * &w[2])) + &x[4],
            ];
            (u, r)
        }
    }
}

/// The witness at the homogeneous configuration `(u, r)`, evaluated
/// directly: for a gap, the plug vertex is turned by quaternion conjugation,
/// `(1, r)(0, p)(1, r)‾ = R̂(r)p`, and read against the support normal,
/// `g = (1 + |r|²) n(u)·h − n(u)·R̂(r)p`; for a domain inequality, D's affine
/// polynomial at `(u₁/u₃, u₂/u₃, r)` times `u₃^k`, negated.
pub(crate) fn independent_value(witness: &Witness, u: &Point, r: &Point) -> QSqrt5 {
    match witness {
        Witness::Gap(gap) => {
            let v = geometry::vertices();
            let n = gap.support().normal_at(u);
            let turned = geometry::rotate(&[QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()], &v[gap.plug()]);
            let scale = &QSqrt5::one() + &geometry::dot(r, r);
            &(&scale * &geometry::dot(&n, &v[gap.support().contact()])) - &geometry::dot(&n, &turned)
        }
        Witness::Domain(inequality) => {
            let inverse = crate::elimination::zoom::reciprocal(&u[2]);
            let (s, t) = (&u[0] * &inverse, &u[1] * &inverse);
            let c = &domain::affine_polynomials()[inequality.index()];
            let value = c.evaluate(&[s, t, r[0].clone(), r[1].clone(), r[2].clone()]);
            let k = domain::view_degree(inequality.index()).unwrap();
            let lambda = (0..k).fold(QSqrt5::one(), |p, _| &p * &u[2]);
            -&(&lambda * &value)
        }
    }
}

/// Validity at one view, by brute force over all vertices; a zero normal
/// makes every condition an equality.
fn brute_force_valid(witness: &Witness, u: &Point) -> bool {
    let Witness::Gap(gap) = witness else { return true };
    let n = gap.support().normal_at(u);
    let v = geometry::vertices();
    let h = geometry::dot(&n, &v[gap.support().contact()]);
    v.iter().all(|w| geometry::dot(&n, w) <= h)
}

fn monomial(y: &Five<QSqrt5>, a: &Exponents) -> QSqrt5 {
    (0..5).fold(QSqrt5::one(), |p, j| (0..a[j]).fold(p, |p, _| &p * &y[j]))
}

fn assert_sound(f: &Fixture, random: &mut Random, points: usize) {
    check_cell(f.zoom(), &f.cell, &f.witness, &f.factor).unwrap();
    let quotient = Quotient::new(f.zoom(), &f.witness, &f.factor).unwrap();
    for i in 0..points {
        // Corners first (where the margins are smallest), then random points.
        let y = if i < 4 {
            let mask = random.below(32);
            exact(&std::array::from_fn(|j| if mask >> j & 1 == 1 { f.cell[j].hi() } else { f.cell[j].lo() }.clone()))
        } else {
            exact(&random.point(&f.cell, 20))
        };
        let (u, r) = configuration(f.zoom(), &y);
        let n = quotient.polynomial().evaluate(&y);
        assert_eq!(independent_value(&f.witness, &u, &r), &monomial(&y, f.factor.exponents()) * &n);
        assert_eq!(n.sign(), std::cmp::Ordering::Less);
        assert!(u[2].sign() == std::cmp::Ordering::Greater);
        assert!(brute_force_valid(&f.witness, &u));
    }
}

#[test]
fn accepted_cells_are_sound_at_exact_points() {
    let mut random = Random::new(31);
    let mut sets = vec![
        exotic("square", 9),
        exotic("pentagon", 100),
        exotic("arc+", 8),
        exotic("arc-", 4),
        exotic("endpoint+", 9),
        exotic("endpoint-", 9),
        exotic("crossing+", 40),
        exotic("crossing-", 20),
    ];
    for n in [0, 16, 29] {
        sets.push(fixtures(local_data(n).unwrap(), 3));
    }
    let mut checked = 0;
    for set in sets {
        assert!(!set.is_empty());
        for f in &set {
            assert_sound(f, &mut random, 6);
            checked += 1;
        }
    }
    assert!(checked > 400, "{checked}");
}

fn first_with(name: &str, pick: impl Fn(&Fixture) -> bool) -> Fixture {
    exotic(name, 1).into_iter().find(|f| pick(f)).expect("a fixture")
}

#[test]
fn wrong_factors_are_refused() {
    let f = first_with("square", |f| f.factor.exponents()[0] > 0 && f.factor.exponents()[1] > 0 && f.cell[0].lo() == &Q::zero());
    // A divided variable that is not a scale variable of the zoom.
    for j in 2..5 {
        let mut factor = *f.factor.exponents();
        factor[j] = 1;
        assert_eq!(
            check_cell(f.zoom(), &f.cell, &f.witness, &ZoomFactor::new(factor)),
            Err(CellRefusal::FactorOutsideScales { variable: j })
        );
    }
    // Raising an exponent past the pulled-back witness's valuation.
    for j in 0..2 {
        let mut factor = *f.factor.exponents();
        loop {
            factor[j] += 1;
            match Quotient::new(f.zoom(), &f.witness, &ZoomFactor::new(factor)) {
                Ok(_) => continue,
                Err(refusal) => {
                    assert!(matches!(refusal, CellRefusal::NotDivisible(_)));
                    break;
                }
            }
        }
    }
    // Lowering the exponent of μ leaves a factor μ in the quotient, which
    // vanishes on the cell's face μ = 0: not strictly negative.
    let mut factor = *f.factor.exponents();
    factor[0] -= 1;
    assert!(matches!(
        check_cell(f.zoom(), &f.cell, &f.witness, &ZoomFactor::new(factor)),
        Err(CellRefusal::NotNegative(_))
    ));
}

#[test]
fn negative_divided_variables_are_refused() {
    let f = first_with("arc+", |f| f.factor.exponents()[1] > 0);
    let mut cell = f.cell.clone();
    let width = cell[1].width();
    cell[1] = Interval::new(-&width, cell[1].hi().clone()).unwrap();
    assert_eq!(
        check_cell(f.zoom(), &cell, &f.witness, &f.factor),
        Err(CellRefusal::NegativeScale { variable: 1 })
    );
}

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

#[test]
fn support_invalid_at_one_corner_is_refused() {
    // Search boxes of the identity zoom and edge supports that are valid
    // at exactly three of the four corner views.
    let root = ConfigurationBox::root();
    let identity = Zoom::identity(ZoomCell(root.axes().clone())).unwrap();
    let mut random = Random::new(32);
    let mut found = 0;
    while found < 20 {
        let s = random.between(&q(0), &frac(2, 3), 8);
        let t = random.between(&q(0), &frac(2, 5), 8);
        let w = frac(1, 1 << (3 + random.below(4)));
        let cell = ZoomCell(std::array::from_fn(|j| match j {
            0 => interval(s.clone(), &s + &w),
            1 => interval(t.clone(), &t + &w),
            _ => interval(frac(-1, 100), frac(1, 100)),
        }));
        let [a, b] = geometry::edges()[random.below(120) as usize];
        let (from, to) = if random.below(2) == 0 { (a, b) } else { (b, a) };
        let support = Support::Edge(Edge::new(from, to).unwrap());
        let views = identity.view_corners(&cell);
        let valid = views.iter().filter(|u| support.check_valid(std::slice::from_ref(u)).is_ok()).count();
        if valid != 3 {
            continue;
        }
        let witness = Witness::Gap(Gap::new(support, from).unwrap());
        assert!(matches!(
            check_cell(&identity, &cell, &witness, &ZoomFactor::ONE),
            Err(CellRefusal::InvalidSupport(crate::elimination::witness::WitnessError::Invalid { .. }))
        ));
        found += 1;
    }
    // A real arc leaf with its edge reversed: the support is on the far side.
    let f = first_with("endpoint+", |f| {
        matches!(&f.witness, Witness::Gap(g) if matches!(g.support(), Support::Edge(_)))
    });
    let Witness::Gap(gap) = &f.witness else { unreachable!() };
    let Support::Edge(edge) = gap.support() else { unreachable!() };
    let reversed = Witness::Gap(Gap::new(Support::Edge(Edge::new(edge.to(), edge.from()).unwrap()), gap.plug()).unwrap());
    assert!(matches!(
        check_cell(f.zoom(), &f.cell, &reversed, &f.factor),
        Err(CellRefusal::InvalidSupport(_)) | Err(CellRefusal::NotDivisible(_))
    ));
}

#[test]
fn views_without_positive_third_coordinate_are_refused() {
    // The identity map in arc-plane coordinates: u₃ = 2 + e.
    let root = ZoomCell(std::array::from_fn(|j| if j == 0 { interval(q(-3), q(1)) } else { interval(q(0), q(1)) }));
    let zoom = Zoom::new("arc".into(), Coordinates::ArcPlane(Sign::Plus), Polynomial::variables(), root, [false; 5]).unwrap();
    let witness = Witness::Domain(crate::elimination::witness::DomainInequality::new(10).unwrap());
    let cell = |lo: Q, hi: Q| -> ZoomCell {
        ZoomCell(std::array::from_fn(|j| if j == 0 { interval(lo.clone(), hi.clone()) } else { interval(q(0), frac(1, 8)) }))
    };
    for (lo, hi) in [(q(-3), q(0)), (q(-2), q(-1)), (q(-3), q(-2)), (frac(-5, 2), frac(-9, 4))] {
        assert!(matches!(check_cell(&zoom, &cell(lo, hi), &witness, &ZoomFactor::ONE), Err(CellRefusal::ViewNotPositive(_))));
    }
    let refusal = check_cell(&zoom, &cell(frac(-1, 2), q(0)), &witness, &ZoomFactor::ONE);
    assert!(!matches!(refusal, Err(CellRefusal::ViewNotPositive(_))));
}

#[test]
fn nonnegative_quotients_are_refused() {
    let root = ConfigurationBox::root();
    let identity = Zoom::identity(ZoomCell(root.axes().clone())).unwrap();
    // The triangle inequality of D holds at the root's origin: −c₀ is not
    // negative there.
    let witness = Witness::Domain(crate::elimination::witness::DomainInequality::new(0).unwrap());
    assert!(matches!(check_cell(&identity, &ZoomCell(root.axes().clone()), &witness, &ZoomFactor::ONE), Err(CellRefusal::NotNegative(_))));
    // Aligned configurations: a gap with the plug vertex at its contact
    // vanishes at r = 0, so a cell containing r = 0 is refused.
    let f = first_with("square", |f| f.factor.exponents() == &[1, 1, 0, 0, 0]);
    let mut factor = *f.factor.exponents();
    factor[1] = 0;
    factor[0] = 0;
    let refused = check_cell(f.zoom(), &f.cell, &f.witness, &ZoomFactor::new(factor));
    assert!(matches!(refused, Err(CellRefusal::NotNegative(_))), "{refused:?}");
    // A witness taken from another leaf of a different zoom.
    let others = exotic("square", 1);
    let mut refusals = 0;
    for g in others.iter().step_by(7) {
        if g.witness != f.witness && check_cell(f.zoom(), &f.cell, &g.witness, &f.factor).is_err() {
            refusals += 1;
        }
    }
    assert!(refusals > 50, "{refusals}");
}

#[test]
fn zero_normal_views_are_dropped_only_where_the_normal_vanishes() {
    // An edge parallel to e₃ is seen end-on at the square view u = (0, 0, 1),
    // a corner view of cells of the identity zoom at s = t = 0.
    let v = geometry::vertices();
    let [from, to] = geometry::edges()
        .iter()
        .copied()
        .find(|&[a, b]| v[a][0] == v[b][0] && v[a][1] == v[b][1])
        .expect("an edge along e₃");
    let witness = Witness::Gap(Gap::new(Support::Edge(Edge::new(from, to).unwrap()), 0).unwrap());
    let identity = Zoom::identity(ZoomCell(ConfigurationBox::root().axes().clone())).unwrap();
    let cell =
        ZoomCell(std::array::from_fn(|j| if j < 2 { interval(q(0), frac(1, 8)) } else { interval(frac(-1, 100), frac(1, 100)) }));
    let views = identity.view_corners(&cell);
    let square = [QSqrt5::zero(), QSqrt5::zero(), QSqrt5::one()];
    let kept = valid_views(&witness, views.clone());
    assert_eq!(kept.len(), 3);
    assert!(views.contains(&square) && !kept.contains(&square));
    // Every validity condition is an equality at the dropped view, and the
    // lemma never asks the witness module about it.
    assert!(brute_force_valid(&witness, &square));
    let refusal = check_cell(&identity, &cell, &witness, &ZoomFactor::ONE);
    assert!(!matches!(refusal, Err(CellRefusal::InvalidSupport(WitnessError::ZeroNormal { .. }))), "{refusal:?}");
    // In the cover data, every dropped view is the square view.
    for f in &exotic("square", 1) {
        let Witness::Gap(gap) = &f.witness else { continue };
        for u in f.zoom().view_corners(&f.cell) {
            if gap.support().normal_at(&u).iter().all(|x| x.is_zero()) {
                assert!(brute_force_valid(&f.witness, &u));
                assert_eq!(u, square);
            }
        }
    }
    // Domain inequalities keep all views.
    let views = vec![square];
    let domain = Witness::Domain(crate::elimination::witness::DomainInequality::new(29).unwrap());
    assert_eq!(valid_views(&domain, views.clone()), views);
}

// ---- The maximum form --------------------------------------------------------

/// Member `w` of the maximum form at the homogeneous configuration `(u, r)`,
/// evaluated directly (quaternion conjugation, as [`independent_value`]).
fn member_value(maximum: &MaximumGap, w: usize, u: &Point, r: &Point) -> QSqrt5 {
    let v = geometry::vertices();
    let n = Support::Edge(maximum.edge().clone()).normal_at(u);
    let turned = geometry::rotate(&[QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()], &v[maximum.plug()]);
    let scale = &QSqrt5::one() + &geometry::dot(r, r);
    &(&scale * &geometry::dot(&n, &v[w])) - &geometry::dot(&n, &turned)
}

fn to_f64(x: &Q) -> f64 {
    use num_traits::ToPrimitive;
    num_rational::BigRational::new(x.numer().clone(), x.denom().clone()).to_f64().unwrap()
}

/// The oriented edge and plug vertex whose maximum form is most negative at
/// the configuration `x` in floating point (a proposal only).
fn best_maximum(x: &Five<Q>) -> MaximumGap {
    let v: Vec<[f64; 3]> = geometry::vertices().iter().map(|p| p.clone().map(|c| {
        to_f64(c.rational_part()) + to_f64(c.sqrt5_part()) * 5f64.sqrt()
    })).collect();
    let x: Vec<f64> = x.iter().map(to_f64).collect();
    let (u, r) = ([x[0], x[1], 1.0], [x[2], x[3], x[4]]);
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: &[f64; 3], b: &[f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let turn = |p: &[f64; 3]| {
        let (rr, rp, c) = (dot(&r, &r), dot(&r, p), cross(&r, p));
        [0, 1, 2].map(|j| (1.0 - rr) * p[j] + 2.0 * rp * r[j] + 2.0 * c[j])
    };
    let mut best = (f64::INFINITY, 0, 0, 0);
    for &[a, b] in geometry::edges() {
        for (from, to) in [(a, b), (b, a)] {
            let n = cross(&u, &[0, 1, 2].map(|j| v[to][j] - v[from][j]));
            let support = v.iter().map(|w| dot(&n, w)).fold(f64::NEG_INFINITY, f64::max) * (1.0 + dot(&r, &r));
            for plug in 0..v.len() {
                let value = (support - dot(&n, &turn(&v[plug]))) / dot(&n, &n).sqrt();
                if value < best.0 {
                    best = (value, from, to, plug);
                }
            }
        }
    }
    MaximumGap::new(Edge::new(best.1, best.2).unwrap(), best.3).unwrap()
}

/// Small cells of the identity zoom around random configurations whose
/// plug clearly sticks out: the maximum form's best proposal is accepted,
/// and at exact corners and random points of the cell every member is
/// negative by direct evaluation (so the maximum is).
#[test]
fn accepted_maximum_forms_are_sound_at_exact_points() {
    let identity = Zoom::identity(ZoomCell(ConfigurationBox::root().axes().clone())).unwrap();
    let mut random = Random::new(34);
    let mut accepted = 0;
    for _ in 0..24 {
        let centre: Five<Q> = std::array::from_fn(|j| match j {
            0 => random.between(&q(0), &frac(2, 3), 10),
            1 => random.between(&q(0), &frac(2, 5), 10),
            _ => random.between(&frac(-2, 5), &frac(2, 5), 10),
        });
        let width = frac(1, 1 << (8 + random.below(4)));
        let cell = ZoomCell(std::array::from_fn(|j| interval(&centre[j] - &width, &centre[j] + &width)));
        let maximum = best_maximum(&centre);
        if check_maximum(&identity, &cell, &maximum, &ZoomFactor::ONE).is_err() {
            continue;
        }
        accepted += 1;
        for i in 0..6 {
            let y = if i < 3 {
                let mask = random.below(32);
                exact(&std::array::from_fn(|j| if mask >> j & 1 == 1 { cell[j].hi() } else { cell[j].lo() }.clone()))
            } else {
                exact(&random.point(&cell, 20))
            };
            let (u, r) = configuration(&identity, &y);
            for w in 0..geometry::VERTEX_COUNT {
                assert_eq!(member_value(&maximum, w, &u, &r).sign(), std::cmp::Ordering::Less, "member {w}");
            }
        }
    }
    assert!(accepted >= 12, "{accepted}");
}

/// Every refusal of the maximum form: a factor outside the zoom, a negative
/// divided variable, a view without positive third coordinate, a factor that
/// does not divide a member, and a cell containing an aligned configuration
/// (where the member of the plug vertex vanishes).
#[test]
fn maximum_form_refusals() {
    let identity = Zoom::identity(ZoomCell(ConfigurationBox::root().axes().clone())).unwrap();
    let [a, b] = geometry::edges()[5];
    let maximum = MaximumGap::new(Edge::new(a, b).unwrap(), 7).unwrap();
    let small = ZoomCell(std::array::from_fn(|_| interval(frac(1, 8), frac(9, 64))));
    assert_eq!(
        check_maximum(&identity, &small, &maximum, &ZoomFactor::new([0, 0, 1, 0, 0])),
        Err(CellRefusal::FactorOutsideScales { variable: 2 })
    );
    // A cover zoom with scale variables: a divided variable made negative.
    let f = first_with("arc+", |f| f.factor.exponents()[1] > 0);
    let mut cell = f.cell.clone();
    let width = cell[1].width();
    cell[1] = Interval::new(-&width, cell[1].hi().clone()).unwrap();
    assert_eq!(check_maximum(f.zoom(), &cell, &maximum, &f.factor), Err(CellRefusal::NegativeScale { variable: 1 }));
    // The factor μ of a square zoom divides a member only if it vanishes at
    // μ = 0, the aligned square view u₀ = e₃, where member 0 is n(u₀)·(v₀ − p):
    // with a plug vertex making that nonzero, member 0 is not divisible.
    let f = first_with("square", |f| f.factor.exponents()[0] > 0);
    let mut factor = [0; 5];
    factor[0] = 1;
    let v = geometry::vertices();
    let square = [QSqrt5::zero(), QSqrt5::zero(), QSqrt5::one()];
    let n = Support::Edge(Edge::new(a, b).unwrap()).normal_at(&square);
    let plug = (0..geometry::VERTEX_COUNT)
        .find(|&p| !(&geometry::dot(&n, &v[0]) - &geometry::dot(&n, &v[p])).is_zero())
        .unwrap();
    let off = MaximumGap::new(Edge::new(a, b).unwrap(), plug).unwrap();
    assert!(matches!(
        check_maximum(f.zoom(), &f.cell, &off, &ZoomFactor::new(factor)),
        Err(CellRefusal::MemberNotDivisible { vertex: 0, .. })
    ));
    // u₃ ≤ 0 on an arc-plane zoom.
    let root = ZoomCell(std::array::from_fn(|j| if j == 0 { interval(q(-3), q(1)) } else { interval(q(0), q(1)) }));
    let zoom = Zoom::new("arc".into(), Coordinates::ArcPlane(Sign::Plus), Polynomial::variables(), root, [false; 5]).unwrap();
    let low = ZoomCell(std::array::from_fn(|j| if j == 0 { interval(q(-3), q(-2)) } else { interval(q(0), frac(1, 8)) }));
    assert!(matches!(check_maximum(&zoom, &low, &maximum, &ZoomFactor::ONE), Err(CellRefusal::ViewNotPositive(_))));
    // Cells containing r = 0 (at a corner, on a face, inside): the member of
    // the plug vertex is 0 there, so no member set is certified.
    for rotation in [interval(q(0), frac(1, 64)), interval(frac(-1, 64), q(0)), interval(frac(-1, 64), frac(1, 64))] {
        let cell = ZoomCell(std::array::from_fn(|j| if j < 2 { interval(frac(1, 8), frac(9, 64)) } else { rotation.clone() }));
        for plug in [0, 7, 33] {
            let maximum = MaximumGap::new(Edge::new(a, b).unwrap(), plug).unwrap();
            assert!(matches!(
                check_maximum(&identity, &cell, &maximum, &ZoomFactor::ONE),
                Err(CellRefusal::MemberNotNegative { .. })
            ));
        }
    }
}

#[test]
fn domain_witnesses_are_checked_on_views_only() {
    let f = first_with("crossing+", |f| matches!(f.witness, Witness::Domain(_)));
    check_cell(f.zoom(), &f.cell, &f.witness, &f.factor).unwrap();
    let mut random = Random::new(33);
    assert_sound(&f, &mut random, 10);
}

#[test]
fn check_all_agrees_with_check_and_names_the_first_refused_cell() {
    let all = exotic("arc+", 1);
    let f = &all[0];
    let same: Vec<&Fixture> = all
        .iter()
        .filter(|g| g.zoom == f.zoom && g.witness == f.witness && g.factor == f.factor && g.cover.parameters() == f.cover.parameters())
        .collect();
    assert!(same.len() > 1);
    let quotient = Quotient::new(f.zoom(), &f.witness, &f.factor).unwrap();
    let mut cells: Vec<ZoomCell> = same.iter().map(|g| g.cell.clone()).collect();
    assert_eq!(quotient.check_all(&cells), Ok(()));
    for cell in &cells {
        assert_eq!(quotient.check(cell), Ok(()));
    }
    // A cell reaching negative values of the divided variable, placed second.
    let mut bad = cells[0].clone();
    let j = (0..5).find(|&j| f.factor.exponents()[j] > 0).unwrap();
    bad[j] = Interval::new(-cells[0][j].width(), cells[0][j].hi().clone()).unwrap();
    cells.insert(1, bad.clone());
    assert_eq!(quotient.check_all(&cells), Err((1, CellRefusal::NegativeScale { variable: j })));
    assert_eq!(quotient.check(&bad), Err(CellRefusal::NegativeScale { variable: j }));
    assert_eq!(quotient.check_all(&cells[..1]), Ok(()));
}
