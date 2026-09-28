//! Reparametrisation invariants of zooms and their coordinates:
//!
//! - the arc-plane coordinates of both signs `σ = ±1` are coordinates of the
//!   same configurations: they agree in `e, v, θ`, differ in `ζ` by an
//!   explicit shift, and map back to the same view and rotation;
//! - a witness at the homogeneous view `λ(s, t, 1)` of arc-plane
//!   coordinates is `λ^k` times its value at the configuration, `λ > 0`, so
//!   the two have the same sign (gaps, maximum-form members, domain
//!   inequalities of every view degree);
//! - every point of a zoom's root whose image lies in the cover is found
//!   again by the coverage map, in a zoom whose image at the found
//!   variables is the same configuration.
use super::*;
use crate::elimination::proof::catalogue::{exotic_pin, local_pin, EXOTIC, LOCAL_COUNT};
use crate::elimination::witness::{DomainInequality, Edge, Gap, MaximumGap, Support, Witness};
use crate::elimination::zoom::cover::{Location, ZoomCover};
use crate::problem::domain;

fn constant(x: &Five<QSqrt5>) -> Five<Polynomial> {
    std::array::from_fn(|j| Polynomial::constant(x[j].clone()))
}

fn value(p: &Polynomial) -> QSqrt5 {
    p.evaluate(&std::array::from_fn(|_| QSqrt5::zero()))
}

fn random_configuration(random: &mut Random) -> Five<QSqrt5> {
    exact(&[
        random.between(&frac(-19, 10), &q(3), 12),
        random.between(&q(-1), &q(1), 12),
        random.small(12),
        random.small(12),
        random.small(12),
    ])
}

#[test]
fn both_arc_plane_labels_describe_the_same_configurations() {
    let mut random = Random::new(61);
    let [a1, _, a3] = arc_rotation();
    for _ in 0..300 {
        let c = random_configuration(&mut random);
        let plus = Coordinates::ArcPlane(Sign::Plus).of_configuration(&c).unwrap();
        let minus = Coordinates::ArcPlane(Sign::Minus).of_configuration(&c).unwrap();
        // The same view and rotation from either label.
        assert_eq!(configuration_at(Coordinates::ArcPlane(Sign::Plus), &plus), configuration_at(Coordinates::ArcPlane(Sign::Minus), &minus));
        // e, v, θ agree; ζ₁, ζ₃ differ by 2(a₁ + a₃θ), 2(a₃ − a₁θ).
        assert_eq!(plus[..3], minus[..3]);
        let theta = &plus[2];
        let two = QSqrt5::integer(2);
        assert_eq!(&minus[3] - &plus[3], &two * &(&a1 + &(&a3 * theta)));
        assert_eq!(&minus[4] - &plus[4], &two * &(&a3 - &(&a1 * theta)));
        // A point of one arc plane is never on the other: r₁ = ±(a₁ + a₃θ)
        // with a₁ + a₃θ ≠ 0 for |θ| < 1 (a₁ > 0 > a₃, |a₃| < a₁).
        if plus[1].is_zero() && plus[3].is_zero() && plus[4].is_zero() && theta.sign() != std::cmp::Ordering::Less {
            assert!(!(minus[3].is_zero() && minus[4].is_zero()));
        }
    }
    assert_eq!(a1.sign(), std::cmp::Ordering::Greater);
    assert_eq!(a3.sign(), std::cmp::Ordering::Less);
}

#[test]
fn witnesses_at_homogeneous_zoom_views_have_the_sign_of_the_configuration() {
    let mut random = Random::new(62);
    let edges = geometry::edges();
    for round in 0..120 {
        let c = random_configuration(&mut random);
        let (affine_u, r) = configuration_at(Coordinates::Configuration, &c);
        let rotation: [Polynomial; 3] = std::array::from_fn(|j| Polynomial::constant(r[j].clone()));
        let affine = View::Affine { s: Polynomial::constant(c[0].clone()), t: Polynomial::constant(c[1].clone()) };
        let [a, b] = edges[random.below(edges.len() as u64) as usize];
        let plug = random.below(60) as usize;
        // Gaps with edge and direction supports and domain inequalities of
        // view degree 0, 1 and 2, with their degree in the view.
        let mut cases: Vec<(Witness, u32)> = vec![
            (Witness::Gap(Gap::new(Support::Edge(Edge::new(a, b).unwrap()), plug).unwrap()), 1),
            (Witness::Domain(DomainInequality::new(round % 73).unwrap()), domain::view_degree(round % 73).unwrap()),
        ];
        let normal = [random.number(4), random.number(4)];
        if let Ok(support) = crate::elimination::witness::Direction::new(normal, random.below(60) as usize) {
            cases.push((Witness::Gap(Gap::new(Support::Direction(support), plug).unwrap()), 1));
        }
        for sign in signs() {
            let coordinates = Coordinates::ArcPlane(sign);
            let x = coordinates.of_configuration(&c).unwrap();
            let (view, rotation_x) = coordinates.configuration(&constant(&x));
            let u = view.vector().map(|p| value(&p));
            // λ = u₃ = 5/(2 + s) > 0 and u = λ(s, t, 1).
            let lambda = u[2].clone();
            assert_eq!(lambda.sign(), std::cmp::Ordering::Greater);
            assert_eq!(u, affine_u.clone().map(|k| &k * &lambda));
            assert_eq!(rotation_x.clone().map(|p| value(&p)), r);
            let power = |k: u32| (0..k).fold(QSqrt5::one(), |p, _| &p * &lambda);
            for (witness, degree) in &cases {
                let at_zoom = value(&witness.polynomial(&view, &rotation_x).unwrap());
                let at_configuration = value(&witness.polynomial(&affine, &rotation).unwrap());
                assert_eq!(at_zoom, &power(*degree) * &at_configuration, "{witness:?}");
                assert_eq!(at_zoom.sign(), at_configuration.sign());
            }
            let maximum = MaximumGap::new(Edge::new(a, b).unwrap(), plug).unwrap();
            let at_zoom = maximum.members(&view, &rotation_x).unwrap();
            let at_configuration = maximum.members(&affine, &rotation).unwrap();
            for (m, n) in at_zoom.iter().zip(&at_configuration) {
                assert_eq!(value(m), &lambda * &value(n));
            }
        }
    }
}

#[test]
fn zoom_images_in_the_cover_are_found_again_by_the_coverage_map() {
    let mut random = Random::new(63);
    let mut pins: Vec<_> = EXOTIC.iter().map(|n| (n.to_string(), exotic_pin(n).unwrap().parameters)).collect();
    pins.extend((0..LOCAL_COUNT).step_by(5).map(|n| (n.to_string(), local_pin(n).unwrap().parameters)));
    let mut found = 0;
    for (name, parameters) in pins {
        let cover = ZoomCover::new(parameters).unwrap();
        let coordinates = cover.parameters().coordinates;
        for (k, zoom) in cover.zooms().iter().enumerate() {
            for round in 0..6 {
                // Dyadic points of the root, on faces now and then.
                let y = exact(&random.point(zoom.root(), if round % 3 == 0 { 1 } else { 10 }));
                let x = zoom.evaluate(&y);
                let (u, r) = configuration_at(coordinates, &x);
                assert_eq!(u[2].sign(), std::cmp::Ordering::Greater, "{name} {}", zoom.name());
                let inverse = reciprocal(&u[2]);
                let c: Five<QSqrt5> = [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()];
                // The configuration's coordinates are the zoom's image again.
                assert_eq!(coordinates.of_configuration(&c).unwrap(), x, "{name} {}", zoom.name());
                match cover.locate(&x) {
                    None => {} // The zoom's image reaches beyond the cover's set.
                    Some(Location::NoFit) => {
                        assert!(coordinates.no_fit_coordinates().iter().all(|&i| x[i].is_zero()));
                    }
                    Some(Location::OutsideDomain) => assert!(cover.parameters().beyond.is_some()),
                    Some(Location::Zoom { zoom: other, y: found_y }) => {
                        let back = cover.zooms()[other].evaluate(&found_y);
                        assert_eq!(back, x, "{name}: zoom {k} → {other}");
                        found += 1;
                    }
                }
            }
        }
    }
    assert!(found > 200, "{found}");
}

#[test]
fn zoom_cells_bisect_into_closed_halves_and_factors_stay_on_scale_variables() {
    let mut random = Random::new(64);
    for _ in 0..50 {
        let cell = ZoomCell(std::array::from_fn(|_| {
            let a = random.small(50);
            let b = &a + &frac(1 + random.below(9) as i64, 1 << random.below(40));
            interval(a, b)
        }));
        let j = random.below(5) as usize;
        let [low, high] = cell.bisect(j);
        assert!(cell.contains_cell(&low) && cell.contains_cell(&high) && cell.contains_cell(&cell));
        assert!(!low.contains_cell(&cell) && !high.contains_cell(&cell));
        assert_eq!(low[j].hi(), high[j].lo());
        assert_eq!((low[j].lo(), high[j].hi()), (cell[j].lo(), cell[j].hi()));
        assert_eq!(&low[j].width() + &high[j].width(), cell[j].width());
        assert!((0..5).filter(|&i| i != j).all(|i| low[i] == cell[i] && high[i] == cell[i]));
    }
    // A point cell bisects into two equal point cells.
    let point = ZoomCell(std::array::from_fn(|_| Interval::point(q(1))));
    assert_eq!(point.bisect(2), [point.clone(), point.clone()]);
    // Factors: only on scale variables; 1 goes with every zoom.
    let covers = ["square", "arc+", "crossing-"].map(|n| ZoomCover::new(exotic_pin(n).unwrap().parameters).unwrap());
    for cover in &covers {
        for zoom in cover.zooms() {
            let scales: Vec<usize> = (0..5).filter(|&j| zoom.is_scale(j)).collect();
            assert!(!scales.is_empty());
            assert_eq!(zoom.factor([0; 5]), Ok(ZoomFactor::ONE));
            for j in 0..5 {
                let mut a = [0u8; 5];
                a[j] = 1 + random.below(3) as u8;
                match zoom.factor(a) {
                    Ok(f) => assert!(scales.contains(&j) && f.exponents() == &a),
                    Err(ZoomError::NotScale { variable }) => assert!(!scales.contains(&j) && variable == j),
                    Err(other) => panic!("{other}"),
                }
            }
        }
    }
    let identity = Zoom::identity(unit_root()).unwrap();
    assert_eq!(identity.factor([1, 0, 0, 0, 0]), Err(ZoomError::NotScale { variable: 0 }));
}
