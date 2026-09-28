//! Tests of the whole generator on small pinned covers (catalogue covers
//! and reduced pins of the same kinds): the loader accepts every generated
//! file, the files do not depend on the thread count, every witness leaf is
//! checked again by an independent exact evaluation at random points of its
//! cell, and every failure is typed and confined to its cover.
use super::*;
use crate::arithmetic::exact::{frac, QSqrt5, Q};
use crate::arithmetic::polynomial::{Polynomial, VARIABLES};
use crate::elimination::zoom::cover::{Parameters, Shape};
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::{tree, Five};
use crate::elimination::proof::catalogue::{exotic_pin, local_pin};
use crate::elimination::proof::ProvedCover;
use crate::problem::{domain, geometry};
use std::cmp::Ordering as Order;

fn settings() -> Settings {
    Settings { max_depth: 60, max_nodes: 100_000, aspect: std::num::NonZeroU32::new(16).unwrap(), tried: 6, proposal_depth: 4, proposals: 8, inherited: 16 }
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn request(pin: Pin, samples: NonZeroUsize, inequalities: Vec<usize>) -> Request {
    Request { pin, samples, inequalities }
}

/// The arc cover's parameters with the base reduced to `[3/64, 1/16]`.
fn short_arc() -> Pin {
    let mut pin = exotic_pin("arc+").unwrap();
    pin.name = "short-arc".into();
    if let Shape::Tube { base, .. } = &mut pin.parameters.shape {
        base[0][1].0 = frac(1, 16);
    }
    pin
}

/// The crossing's parameters for `σ = +1` with radii `1/64` (the window is
/// kept).
fn small_crossing() -> Pin {
    let mut pin = exotic_pin("crossing+").unwrap();
    pin.name = "small-crossing".into();
    let parameters: &mut Parameters = &mut pin.parameters;
    if let Shape::Point { radii, radius, .. } = &mut parameters.shape {
        for r in radii.iter_mut() {
            r.0 = frac(1, 64);
        }
        radius.0 = frac(1, 64);
    }
    pin
}

fn run(requests: &[Request], s: &Settings, n: usize) -> Vec<Result<Generated, GenerateError>> {
    generate(requests, s, threads(n), &AtomicBool::new(false), &mut |_| {}).unwrap()
}

/// The configuration `(u, r)` at coordinates `x` (constant polynomials).
fn configuration(coordinates: crate::elimination::zoom::Coordinates, x: &Five<QSqrt5>) -> (Point, Point) {
    let (view, rotation) = coordinates.configuration(&x.clone().map(Polynomial::constant));
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    (
        std::array::from_fn(|i| view.vector()[i].evaluate(&origin)),
        std::array::from_fn(|i| rotation[i].evaluate(&origin)),
    )
}

/// Checks every witness leaf of a generated file at `points` random exact
/// points of its cell, independently of the pull-back and the certifier:
/// where the zoom factor is positive the witness is negative at the
/// configuration (a gap by quaternion conjugation, a domain inequality by
/// D's affine polynomial at the affine view), and a gap's support is valid
/// at the point's view by brute force over all vertices. Returns the number
/// of points checked.
fn check_independently(pin: &Pin, file: &CoverFile, points: usize, random: &mut Random) -> usize {
    let vertices = geometry::vertices();
    let mut checked = 0;
    {
        let (entry, parameters) = (file, &pin.parameters);
        let cover = ZoomCover::new(parameters.clone()).unwrap();
        for (zoom, data) in cover.zooms().iter().zip(&entry.zooms) {
            let cells = tree::leaves(zoom.root(), &data.tree).unwrap();
            for (cell, leaf) in cells.iter().zip(&data.leaves) {
                let Leaf::Witness { index, factor } = leaf else { continue };
                for _ in 0..points {
                    let y: Five<QSqrt5> = random.point(cell, 12).map(QSqrt5::from_rational);
                    let monomial =
                        (0..VARIABLES).fold(QSqrt5::one(), |m, j| (0..factor[j]).fold(m, |m, _| &m * &y[j]));
                    let (u, r) = configuration(parameters.coordinates, &zoom.evaluate(&y));
                    assert_eq!(u[2].sign(), Order::Greater);
                    let value = match &file.witnesses[*index] {
                        Witness::Gap(gap) => {
                            let n = gap.support().normal_at(&u);
                            let support = geometry::dot(&n, &vertices[gap.support().contact()]);
                            assert!(vertices.iter().all(|w| geometry::dot(&n, w) <= support), "invalid support");
                            let quaternion = [QSqrt5::one(), r[0].clone(), r[1].clone(), r[2].clone()];
                            let turned = geometry::rotate(&quaternion, &vertices[gap.plug()]);
                            let scale = &QSqrt5::one() + &geometry::dot(&r, &r);
                            &(&scale * &support) - &geometry::dot(&n, &turned)
                        }
                        Witness::Domain(inequality) => {
                            let inverse = crate::elimination::zoom::reciprocal(&u[2]);
                            let affine: Five<QSqrt5> =
                                [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()];
                            -&domain::affine_polynomials()[inequality.index()].evaluate(&affine)
                        }
                    };
                    if monomial.sign() == Order::Greater {
                        assert_eq!(value.sign(), Order::Less, "{} {}", pin.name, zoom.name());
                    }
                    checked += 1;
                }
            }
        }
    }
    checked
}

/// Generates `requests` with 1 and `threads` threads (or only `threads`),
/// requires identical files, verifies them with the loader and checks every
/// witness leaf independently. Returns the delegated leaves.
fn generate_and_check(requests: &[Request], compare: bool, n: usize) -> usize {
    let many = run(requests, &settings(), n);
    let one = if compare { Some(run(requests, &settings(), 1)) } else { None };
    let mut random = Random::new(31);
    let mut delegated = 0;
    for (k, (request, a)) in requests.iter().zip(many).enumerate() {
        let a = a.unwrap();
        if let Some(one) = &one {
            let b = one[k].as_ref().unwrap();
            assert_eq!(a.file.to_bytes(), b.file.to_bytes(), "{}", request.pin.name);
            assert_eq!(a.statistics.nodes, b.statistics.nodes);
        }
        let bytes = a.file.to_bytes();
        let cover = ProvedCover::load(&request.pin, &bytes, threads(2), &AtomicBool::new(false)).unwrap();
        let report = cover.report();
        assert_eq!(report.witness_leaves, a.statistics.witness_leaves);
        assert_eq!(report.delegated_leaves, a.statistics.delegated_leaves);
        assert_eq!(report.domain_leaves, a.statistics.domain_leaves);
        assert_eq!(report.witnesses, a.file.witnesses.len());
        let leaves = a.statistics.witness_leaves + a.statistics.delegated_leaves;
        assert_eq!(a.statistics.nodes, 2 * leaves - a.statistics.zooms);
        delegated += a.statistics.delegated_leaves;
        // Witnesses are numbered by first use, each used, none repeated.
        let mut first_use = Vec::new();
        {
            for zoom in &a.file.zooms {
                for leaf in &zoom.leaves {
                    if let Leaf::Witness { index, .. } = leaf {
                        if !first_use.contains(index) {
                            first_use.push(*index);
                        }
                    }
                }
            }
        }
        assert_eq!(first_use, (0..a.file.witnesses.len()).collect::<Vec<_>>());
        assert!(check_independently(&request.pin, &a.file, 2, &mut random) > 0);
    }
    delegated
}

#[test]
fn generated_covers_verify_and_do_not_depend_on_the_thread_count() {
    let requests = vec![
        request(local_pin(19).unwrap(), NonZeroUsize::new(2).unwrap(), vec![]),
        request(local_pin(21).unwrap(), NonZeroUsize::new(1).unwrap(), vec![]),
        request(short_arc(), NonZeroUsize::new(2).unwrap(), vec![10, 11]),
    ];
    generate_and_check(&requests, true, 3);
}

#[test]
fn a_generated_window_cover_verifies() {
    let requests = [request(small_crossing(), NonZeroUsize::new(2).unwrap(), vec![10])];
    assert!(generate_and_check(&requests, false, 4) > 0, "the window is exercised");
}

/// Heavy: every cover of the catalogue, generated with 4 and 1 threads,
/// identical and verified (about 17 minutes on one core plus 5 on four):
///
/// ```text
/// cargo test --release --lib generate::tests::catalogue_campaign -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn catalogue_campaign() {
    let two = NonZeroUsize::new(2).unwrap();
    let eight = NonZeroUsize::new(8).unwrap();
    let mut requests = vec![
        request(exotic_pin("square").unwrap(), NonZeroUsize::new(1).unwrap(), vec![29]),
        request(exotic_pin("pentagon").unwrap(), two, vec![]),
    ];
    for name in ["arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"] {
        requests.push(request(exotic_pin(name).unwrap(), eight, vec![10, 11]));
    }
    for number in 0..crate::elimination::proof::catalogue::LOCAL_COUNT {
        requests.push(request(local_pin(number).unwrap(), two, vec![]));
    }
    generate_and_check(&requests, true, 4);
}

#[test]
fn failures_are_typed_and_confined_to_their_cover() {
    let good = request(local_pin(21).unwrap(), NonZeroUsize::new(1).unwrap(), vec![]);
    let results = run(
        &[
            request(local_pin(19).unwrap(), NonZeroUsize::new(1).unwrap(), vec![29]),
            request(exotic_pin("square").unwrap(), NonZeroUsize::new(1).unwrap(), vec![73]),
            good.clone(),
        ],
        &settings(),
        2,
    );
    assert!(matches!(&results[0], Err(GenerateError::Scope { cover }) if cover == "19"));
    assert!(matches!(&results[1], Err(GenerateError::Inequality { cover, .. }) if cover == "square"));
    assert!(results[2].is_ok());
    // The depth limit names the zoom and the cell.
    let shallow = Settings { max_depth: 0, ..settings() };
    let results = run(&[good.clone()], &shallow, 2);
    let Err(GenerateError::Depth { cover, zoom, cell }) = &results[0] else {
        panic!("{:?}", results[0].as_ref().err())
    };
    assert_eq!(cover, "21");
    let cover = ZoomCover::new(good.pin.parameters.clone()).unwrap();
    assert_eq!(cell, cover.zooms().iter().find(|c| c.name() == zoom).unwrap().root());
    // The node limit.
    let few = Settings { max_nodes: 2, ..settings() };
    assert!(matches!(&run(&[good.clone()], &few, 2)[0], Err(GenerateError::Nodes { .. })));
    // A broken pin is refused by the cover's own checks.
    let mut broken = good.clone();
    if let Shape::Tube { radius, .. } = &mut broken.pin.parameters.shape {
        radius.0 = Q::zero();
    }
    assert!(matches!(&run(&[broken], &settings(), 1)[0], Err(GenerateError::Cover { .. })));
    // A stop request ends the generation.
    let result = generate(&[good], &settings(), threads(2), &AtomicBool::new(true), &mut |_| {});
    assert!(matches!(result, Err(GenerateError::Interrupted)));
}

#[test]
fn progress_is_reported_level_by_level() {
    let mut seen: Vec<Progress> = Vec::new();
    let result = generate(
        &[request(local_pin(21).unwrap(), NonZeroUsize::new(1).unwrap(), vec![])],
        &settings(),
        threads(2),
        &AtomicBool::new(false),
        &mut |p| seen.push(*p),
    )
    .unwrap();
    let generated = result.into_iter().next().unwrap().unwrap();
    assert_eq!(seen[0], Progress { depth: 0, cells: generated.statistics.zooms, decided: 0 });
    for pair in seen.windows(2) {
        assert_eq!(pair[1].depth, pair[0].depth + 1);
        assert_eq!(pair[1].decided, pair[0].decided + pair[0].cells);
    }
    let last = seen.last().unwrap();
    assert_eq!(last.decided + last.cells, generated.statistics.nodes);
}

#[test]
fn parallel_keeps_the_order_and_stops() {
    let items: Vec<usize> = (0..1000).collect();
    let squares = parallel(&items, threads(4), &AtomicBool::new(false), |x| x * x).unwrap();
    assert_eq!(squares, items.iter().map(|x| x * x).collect::<Vec<_>>());
    assert!(parallel(&items, threads(4), &AtomicBool::new(true), |x| *x).is_err());
    assert!(parallel(&Vec::<usize>::new(), threads(4), &AtomicBool::new(false), |x| *x).unwrap().is_empty());
}

mod review;
