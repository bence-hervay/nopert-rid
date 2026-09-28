//! Tests of candidate construction: exact properties of the proposed
//! supports and gaps, completeness of the proposals against brute force,
//! samples in the no-fit sets, factors that are exactly maximal, the float
//! validity filter never excluding an exactly valid support, and separating
//! proposals that really separate.
use super::*;
use std::num::NonZeroUsize;
use crate::arithmetic::exact::{q, Interval};
use crate::arithmetic::polynomial::Polynomial;
use crate::elimination::zoom::cover::ZoomCover;
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::{Coordinates, ZoomCell};
use crate::elimination::proof::catalogue::{exotic_pin, local_pin};
use std::cmp::Ordering;

fn view(random: &mut Random) -> Point {
    [
        QSqrt5::from_rational(random.between(&frac(1, 97), &frac(2, 3), 10)),
        QSqrt5::from_rational(random.between(&frac(1, 89), &frac(2, 5), 10)),
        QSqrt5::one(),
    ]
}

fn zero() -> Point {
    std::array::from_fn(|_| QSqrt5::zero())
}

/// The exact gap `(1 + |r|²) n·h − n·R̂(r)p` at a constant configuration,
/// with `R̂(r)p` computed by quaternion conjugation (not the Cayley matrix).
fn exact_gap(gap: &Gap, u: &Point, r: &Point) -> QSqrt5 {
    let n = gap.support().normal_at(u);
    let v = geometry::vertices();
    let quaternion = [QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()];
    let turned = geometry::rotate(&quaternion, &v[gap.plug()]);
    let scale = &QSqrt5::one() + &geometry::dot(r, r);
    &(&scale * &geometry::dot(&n, &v[gap.support().contact()])) - &geometry::dot(&n, &turned)
}

#[test]
fn cone_directions_are_exactly_valid_and_supports_are_complete() {
    let mut random = Random::new(11);
    for _ in 0..30 {
        let u = view(&mut random);
        let list = supports(&u);
        let edges = list.iter().filter(|s| matches!(s, Support::Edge(_))).count();
        assert_eq!(edges, 240);
        let outline = float::hull(
            &geometry::vertices().iter().map(|v| float::screen(&float::point(&u), &float::point(v))).collect::<Vec<_>>(),
        );
        let directions = &list[240..];
        // radial and plain directions: two per vertex unless a direction is zero
        assert!(directions.len() >= 3 * outline.len() + 100);
        for cone in &directions[directions.len() - 3 * outline.len()..] {
            cone.check_valid(std::slice::from_ref(&u)).expect("a direction in a vertex's normal cone is valid");
        }
        // The radial and plain directions come vertex by vertex (at most two
        // each, fewer where a direction is zero), with that vertex as contact.
        let contacts: Vec<usize> = directions[..directions.len() - 3 * outline.len()].iter().map(|s| s.contact()).collect();
        assert!(contacts.windows(2).all(|w| w[0] <= w[1]), "{contacts:?}");
        for c in 0..geometry::VERTEX_COUNT {
            assert!(contacts.iter().filter(|&&k| k == c).count() <= 2);
        }
    }
}

#[test]
fn active_gaps_at_the_aligned_set_touch_and_are_complete() {
    let mut random = Random::new(12);
    for _ in 0..12 {
        let u = view(&mut random);
        let gaps = active(&u, &zero());
        assert!(!gaps.is_empty());
        for gap in &gaps {
            gap.support().check_valid(std::slice::from_ref(&u)).expect("valid at the sample");
            assert_eq!(exact_gap(gap, &u, &zero()).sign(), Ordering::Equal, "{gap:?}");
        }
        // Brute force: every exactly valid oriented edge appears with both of
        // its end points as plug vertices.
        for &[a, b] in geometry::edges() {
            for (from, to) in [(a, b), (b, a)] {
                let support = Support::Edge(Edge::new(from, to).unwrap());
                if support.check_valid(std::slice::from_ref(&u)).is_ok() {
                    for plug in [from, to] {
                        assert!(gaps.contains(&Gap::new(support.clone(), plug).unwrap()), "{from} {to} {plug}");
                    }
                }
            }
        }
    }
}

#[test]
fn active_gaps_on_an_arc_plane_are_not_positive() {
    // A point of the arc plane Π₊: v = ζ₁ = ζ₃ = 0.
    let x: Five<QSqrt5> = [QSqrt5::from_rational(frac(1, 10)), QSqrt5::zero(), QSqrt5::zero(), QSqrt5::zero(), QSqrt5::zero()];
    let coordinates = Coordinates::ArcPlane(crate::elimination::zoom::Sign::Plus);
    let constant = x.clone().map(Polynomial::constant);
    let (view, rotation) = coordinates.configuration(&constant);
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    let u: Point = std::array::from_fn(|i| view.vector()[i].evaluate(&origin));
    let r: Point = std::array::from_fn(|i| rotation[i].evaluate(&origin));
    let gaps = active(&u, &r);
    assert!(!gaps.is_empty());
    for gap in &gaps {
        assert_ne!(exact_gap(gap, &u, &r).sign(), Ordering::Greater, "{gap:?}");
    }
}

#[test]
fn samples_lie_in_the_no_fit_sets() {
    let covers: Vec<ZoomCover> = ["square", "pentagon", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]
        .iter()
        .map(|name| exotic_pin(name).unwrap().parameters)
        .chain([local_pin(7).unwrap().parameters])
        .map(|p| ZoomCover::new(p).unwrap())
        .collect();
    for cover in &covers {
        let shape = &cover.parameters().shape;
        let base = shape.base_dimension();
        for zoom in 0..cover.zooms().len() {
            let list = samples(cover, zoom, NonZeroUsize::new(2).unwrap());
            assert!(!list.is_empty());
            for (u, r) in &list {
                assert_eq!(u[2].sign(), Ordering::Greater);
                // Back to coordinates: the offset part of the cover coordinates
                // vanishes, in the zoom's own (possibly sheared) family.
                let c: Five<QSqrt5> = {
                    let inverse = crate::elimination::zoom::reciprocal(&u[2]);
                    [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()]
                };
                let x = cover.parameters().coordinates.of_configuration(&c).unwrap();
                let no_fit = cover.parameters().coordinates.no_fit_coordinates();
                let sheared = matches!(cover.kind(zoom), Kind::Sheared { .. });
                let z = cover.cover_coordinates(&x);
                if sheared {
                    // On the second arc θ = −e: the offset coordinates other than θ vanish.
                    assert!(no_fit.iter().all(|&i| x[i].is_zero()));
                } else {
                    assert!(z[base..].iter().all(QSqrt5::is_zero), "{z:?}");
                }
            }
        }
    }
}

#[test]
fn factors_are_exactly_maximal_and_quotients_exact() {
    let mut random = Random::new(13);
    for name in ["square", "endpoint+"] {
        let cover = ZoomCover::new(exotic_pin(name).unwrap().parameters).unwrap();
        for zoom in [0, cover.zooms().len() - 1] {
            let c = &cover.zooms()[zoom];
            let (u, r) = samples(&cover, zoom, NonZeroUsize::new(1).unwrap()).into_iter().last().unwrap();
            for gap in active(&u, &r).into_iter().step_by(7) {
                let witness = Witness::Gap(gap);
                let Some(candidate) = Candidate::new(c, witness.clone()).unwrap() else { continue };
                let pulled = witness.polynomial(c.view(), c.rotation()).unwrap();
                let quotient = pulled.divide_by_monomial(candidate.factor().exponents()).unwrap();
                for j in 0..VARIABLES {
                    let mut more = *candidate.factor().exponents();
                    more[j] += 1;
                    if c.is_scale(j) {
                        assert!(pulled.divide_by_monomial(&more).is_err(), "factor not maximal");
                    } else {
                        assert_eq!(candidate.factor().exponents()[j], 0);
                    }
                }
                for _ in 0..5 {
                    let y = random.point(c.root(), 6);
                    let exact = float::number(&quotient.evaluate(&y.clone().map(QSqrt5::from_rational)));
                    let value = candidate.quotient().evaluate(&y.map(|x| float::rational(&x)));
                    assert!((exact - value).abs() <= 1e-9 * exact.abs().max(1.0));
                }
            }
        }
    }
}

#[test]
fn a_vanishing_pull_back_is_no_candidate() {
    // A zoom without rotation: a gap whose plug is its contact vanishes.
    let variables = Polynomial::variables();
    let map = [variables[0].clone(), variables[1].clone(), Polynomial::zero(), Polynomial::zero(), Polynomial::zero()];
    let root = ZoomCell(std::array::from_fn(|_| Interval::new(q(0), q(1)).unwrap()));
    let zoom = Zoom::new("flat".into(), Coordinates::Configuration, map, root, [false; VARIABLES]).unwrap();
    let [a, b] = geometry::edges()[0];
    let gap = Gap::new(Support::Edge(Edge::new(a, b).unwrap()), a).unwrap();
    assert!(Candidate::new(&zoom, Witness::Gap(gap)).unwrap().is_none());
    assert!(inequality(73).is_err());
    assert!(inequality(72).is_ok());
}

#[test]
fn the_validity_filter_never_excludes_an_exactly_valid_support() {
    let mut random = Random::new(14);
    let zoom = Zoom::identity(ZoomCell(std::array::from_fn(|_| Interval::new(q(0), q(1)).unwrap()))).unwrap();
    let mut kept_invalid = 0;
    for _ in 0..20 {
        let views: Vec<Point> = (0..1 + random.below(3)).map(|_| view(&mut random)).collect();
        let float_views: Vec<Vector> = views.iter().map(float::point).collect();
        for support in supports(&views[0]).into_iter().step_by(3) {
            let gap = Witness::Gap(Gap::new(support.clone(), 0).unwrap());
            let candidate = Candidate::new(&zoom, gap).unwrap().unwrap();
            let exact = support.check_valid(&views).is_ok();
            let plausible = candidate.plausibly_valid(&float_views);
            assert!(!exact || plausible, "an exactly valid support was filtered out");
            kept_invalid += usize::from(plausible && !exact);
        }
    }
    assert!(kept_invalid < 5, "the filter keeps {kept_invalid} clearly invalid supports");
}

#[test]
fn separating_proposals_separate_protruding_plugs() {
    let mut random = Random::new(15);
    let mut separated = 0;
    for _ in 0..40 {
        let u = view(&mut random);
        let r: Point = std::array::from_fn(|_| QSqrt5::from_rational(random.between(&frac(-1, 20), &frac(1, 20), 8)));
        let gaps = separating(&float::point(&u), &float::point(&r), 4);
        for gap in &gaps {
            if let Support::Direction(d) = gap.support() {
                for x in d.screen_normal() {
                    assert!(x.sqrt5_part().is_zero());
                    let scaled = x.rational_part() * &q(DIRECTION_GRID);
                    assert!(scaled.denom() == &num_bigint::BigInt::from(1), "on the dyadic grid");
                }
            }
        }
        // The first proposal is exactly valid with a negative gap whenever the
        // most protruding plug vertex protrudes clearly (float margin).
        let first = &gaps[0];
        let valid = first.support().check_valid(std::slice::from_ref(&u)).is_ok();
        if valid && exact_gap(first, &u, &r).sign() == Ordering::Less {
            separated += 1;
        }
    }
    assert!(separated >= 30, "only {separated} of 40 proposals separate");
}

