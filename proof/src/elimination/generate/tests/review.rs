//! Adversarial tests of the generator: settings
//! that the configuration accepts but the search cannot handle, failures and
//! files that must not depend on the thread count (also with proposals and
//! inherited candidates), mutated generated files that the loader must refuse
//! unless an independent exact evaluation confirms them, and degenerate float
//! inputs.
use super::check_independently;
use crate::arithmetic::exact::{frac, Interval, QSqrt5};
use crate::arithmetic::polynomial::VARIABLES;
use crate::elimination::zoom::cover::{ZoomCover, Kind};
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::{Five, Zoom, ZoomCell};
use crate::elimination::generate::candidates::{self, Candidate};
use crate::elimination::generate::command::Configuration;
use crate::elimination::float::{self, Vector};
use crate::elimination::generate::search::{split_axis, ZoomSearch, Settings};
use crate::elimination::generate::{generate, GenerateError, Generated, Request};
use crate::elimination::proof::catalogue::{exotic_pin, local_pin, LOCAL_COUNT};
use crate::elimination::proof::format::{Leaf, CoverFile};
use crate::elimination::proof::ProvedCover;
use crate::elimination::witness::{Gap, Support, Witness};
use crate::problem::geometry::{self, Point};
use std::num::{NonZeroU32, NonZeroUsize};
use std::sync::atomic::AtomicBool;

fn settings() -> Settings {
    Settings { max_depth: 60, max_nodes: 200_000, aspect: NonZeroU32::new(16).unwrap(), tried: 6, proposal_depth: 4, proposals: 8, inherited: 16 }
}

fn one(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn local(number: usize, samples: usize) -> Request {
    Request { pin: local_pin(number).unwrap(), samples: one(samples), inequalities: Vec::new() }
}

fn run(requests: &[Request], s: &Settings, threads: usize) -> Vec<Result<Generated, GenerateError>> {
    generate(requests, s, one(threads), &AtomicBool::new(false), &mut |_| {}).unwrap()
}

/// A configuration text with the given search settings and Local samples.
fn configuration(aspect: u64, samples: u64, max_depth: u64) -> String {
    format!(
        r#"{{"output":"x","threads":1,"search":{{"max_depth":{max_depth},"max_nodes":1000,"aspect":{aspect},"tried":6,"proposal_depth":4,"proposals":8,"inherited":16}},"exotic":[],"local":{{"numbers":[21],"samples":{samples}}}}}"#
    )
}

// ---- Settings the search cannot handle are refused ------------------------
//
// (`aspect: 0`, `samples: 2^32` and a very large `max_depth`, which would
// make float cell widths vanish, are refused with typed errors before
// any work, and the split rule uses exact split counts.)

#[test]
fn an_aspect_of_zero_is_refused_by_the_configuration() {
    let error = Configuration::parse(configuration(0, 1, 60).as_bytes()).unwrap_err();
    assert!(error.to_string().contains("nonzero") || error.to_string().contains("zero"), "{error}");
}

#[test]
fn a_huge_sample_grid_is_refused_before_any_work() {
    let huge = Configuration::parse(configuration(16, 1 << 32, 60).as_bytes()).unwrap();
    assert!(matches!(
        huge.requests(),
        Err(crate::elimination::generate::command::Error::Generate(GenerateError::Samples { samples, .. })) if samples == 1 << 32
    ));
    // The library refuses it too, for its cover only.
    let results = run(&[local(21, 1 << 32)], &settings(), 1);
    assert!(matches!(results[0], Err(GenerateError::Samples { .. })));
    // The largest accepted grid is still enumerable.
    assert!(Configuration::parse(configuration(16, 64, 60).as_bytes()).unwrap().requests().is_ok());
}

#[test]
fn a_huge_max_depth_is_refused_and_deep_cells_split_exactly() {
    let deep_configuration = Configuration::parse(configuration(16, 1, 100_000).as_bytes()).unwrap();
    assert!(matches!(
        deep_configuration.requests(),
        Err(crate::elimination::generate::command::Error::Generate(GenerateError::MaxDepth { max_depth: 100_000 }))
    ));
    let deep = Settings { max_depth: 100_000, ..settings() };
    assert!(matches!(
        generate(&[local(21, 1)], &deep, one(1), &AtomicBool::new(false), &mut |_| {}),
        Err(GenerateError::MaxDepth { .. })
    ));
    // A cell of positive exact width below float resolution: the split rule
    // no longer depends on float widths, so the decision is a split.
    let cover = ZoomCover::new(local_pin(21).unwrap().parameters.clone()).unwrap();
    let search = ZoomSearch::new(&cover, 0, &[], &[]).unwrap();
    let mut tiny = search.zoom().root().clone();
    for j in 0..VARIABLES {
        for _ in 0..80 {
            tiny[j] = tiny[j].bisect()[1].clone();
        }
    }
    assert!(float::cell(&tiny).iter().all(|[lo, hi]| lo == hi), "the float ends coincide");
    let decision = search.decide(&tiny, 400, &[], &Settings { max_depth: 1000, proposal_depth: 1000, ..settings() }).unwrap();
    assert!(matches!(decision, crate::elimination::generate::search::Decision::Split { axis: 0, .. }), "{decision:?}");
}

// ---- Determinism --------------------------------------------------------

#[test]
fn failures_do_not_depend_on_the_thread_count() {
    let requests: Vec<Request> = [19, 22, 14, 12].into_iter().map(|n| local(n, 2)).collect();
    for s in [
        Settings { max_depth: 6, ..settings() },
        Settings { max_nodes: 30, ..settings() },
        Settings { max_depth: 9, max_nodes: 60, ..settings() },
    ] {
        let text = |results: Vec<Result<Generated, GenerateError>>| -> Vec<String> {
            results
                .into_iter()
                .map(|r| match r {
                    Ok(g) => String::from_utf8(g.file.to_bytes()).unwrap(),
                    Err(e) => e.to_string(),
                })
                .collect()
        };
        let a = text(run(&requests, &s, 1));
        let b = text(run(&requests, &s, 4));
        assert_eq!(a, b, "{s:?}");
        assert!(a.iter().any(|t| t.contains("cover ")), "a failure is exercised: {s:?}");
    }
}

#[test]
fn proposals_and_inherited_candidates_do_not_depend_on_the_thread_count() {
    // Local 22 and 14 lie next to the pentagon view and need proposals.
    let requests = vec![local(22, 2), local(14, 2)];
    let a = run(&requests, &settings(), 1);
    let b = run(&requests, &settings(), 3);
    let mut proposed = 0;
    for (k, (x, y)) in a.iter().zip(&b).enumerate() {
        let (x, y) = (x.as_ref().unwrap(), y.as_ref().unwrap());
        assert_eq!(x.file.to_bytes(), y.file.to_bytes());
        assert_eq!(x.statistics, y.statistics);
        proposed += x.statistics.proposed_leaves + x.statistics.inherited_leaves;
        ProvedCover::load(&requests[k].pin, &x.file.to_bytes(), one(2), &AtomicBool::new(false)).unwrap();
    }
    assert!(proposed > 0, "proposals are exercised");
}

// ---- Verification is independent of generation -------------------------

/// Every mutated file that the loader accepts must pass the independent
/// exact evaluation at random points of every witness leaf.
fn accepted_mutants_are_confirmed(pin: &crate::elimination::proof::Pin, mutants: Vec<(String, CoverFile)>) -> [usize; 2] {
    let mut random = Random::new(71);
    let mut counts = [0; 2];
    for (what, file) in mutants {
        match ProvedCover::load(pin, &file.to_bytes(), one(2), &AtomicBool::new(false)) {
            Ok(_) => {
                counts[1] += 1;
                assert!(check_independently(pin, &file, 3, &mut random) > 0, "{what}");
            }
            Err(_) => counts[0] += 1,
        }
    }
    counts
}

#[test]
fn mutated_generated_files_are_refused_or_independently_confirmed() {
    let request = local(19, 2);
    let generated = run(&[request.clone()], &settings(), 2).remove(0).unwrap();
    let file = generated.file;
    let mut mutants: Vec<(String, CoverFile)> = Vec::new();
    for (c, zoom) in file.zooms.iter().enumerate() {
        for (l, leaf) in zoom.leaves.iter().enumerate() {
            let Leaf::Witness { index, factor } = leaf else { continue };
            // Another witness on the same cell.
            for other in [(index + 1) % file.witnesses.len(), (index + 7) % file.witnesses.len()] {
                let mut m = file.clone();
                m.zooms[c].leaves[l] = Leaf::Witness { index: other, factor: *factor };
                mutants.push((format!("zoom {c} leaf {l} witness {other}"), m));
            }
            // A larger and a smaller factor.
            for j in 0..VARIABLES {
                let mut more = *factor;
                more[j] += 1;
                let mut m = file.clone();
                m.zooms[c].leaves[l] = Leaf::Witness { index: *index, factor: more };
                mutants.push((format!("zoom {c} leaf {l} factor +{j}"), m));
                if factor[j] > 0 {
                    let mut less = *factor;
                    less[j] -= 1;
                    let mut m = file.clone();
                    m.zooms[c].leaves[l] = Leaf::Witness { index: *index, factor: less };
                    mutants.push((format!("zoom {c} leaf {l} factor -{j}"), m));
                }
            }
            // A delegation without a window.
            let mut m = file.clone();
            m.zooms[c].leaves[l] = Leaf::Delegated;
            mutants.push((format!("zoom {c} leaf {l} delegated"), m));
        }
        // Two sibling leaves merged into their parent (keeping the first
        // leaf's entry), and a split along another variable.
        if let Some(at) = zoom.tree.find("..") {
            if at > 0 {
                let mut m = file.clone();
                let tree = &mut m.zooms[c].tree;
                let leaf = tree[..at].matches('.').count();
                tree.replace_range(at - 1..at + 2, ".");
                m.zooms[c].leaves.remove(leaf + 1);
                mutants.push((format!("zoom {c} merged at {at}"), m));
                let mut m = file.clone();
                let tree = &mut m.zooms[c].tree;
                let axis = tree.as_bytes()[at - 1] - b'0';
                let other = char::from(b'0' + (axis + 1) % VARIABLES as u8);
                tree.replace_range(at - 1..at, &other.to_string());
                mutants.push((format!("zoom {c} split {other} at {at}"), m));
            }
        }
    }
    // Every witness with another plug vertex or another contact.
    for (w, witness) in file.witnesses.iter().enumerate() {
        let Witness::Gap(gap) = witness else { continue };
        let mut m = file.clone();
        m.witnesses[w] = Witness::Gap(Gap::new(gap.support().clone(), (gap.plug() + 1) % 60).unwrap());
        mutants.push((format!("witness {w} plug"), m));
        if let Support::Edge(edge) = gap.support() {
            let reversed = crate::elimination::witness::Edge::new(edge.to(), edge.from()).unwrap();
            let mut m = file.clone();
            m.witnesses[w] = Witness::Gap(Gap::new(Support::Edge(reversed), gap.plug()).unwrap());
            mutants.push((format!("witness {w} reversed"), m));
        }
    }
    let total = mutants.len();
    let [refused, accepted] = accepted_mutants_are_confirmed(&request.pin, mutants);
    assert_eq!(refused + accepted, total);
    assert!(refused * 2 > total, "most mutants are refused: {refused} of {total}");
}

// ---- Proposals and float filters ------------------------------------------

#[test]
fn samples_of_every_cover_lie_in_their_no_fit_sets() {
    let mut pins: Vec<_> = crate::elimination::proof::catalogue::EXOTIC
        .iter()
        .map(|n| exotic_pin(n).unwrap())
        .collect();
    pins.extend((0..LOCAL_COUNT).map(|n| local_pin(n).unwrap()));
    for pin in pins {
        {
            let cover = ZoomCover::new(pin.parameters.clone()).unwrap();
            let k = cover.parameters().shape.base_dimension();
            for zoom in 0..cover.zooms().len() {
                for (u, r) in candidates::samples(&cover, zoom, one(1)) {
                    assert_eq!(u[2].sign(), std::cmp::Ordering::Greater);
                    let inverse = crate::elimination::zoom::reciprocal(&u[2]);
                    let c: Five<QSqrt5> =
                        [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()];
                    let x = cover.parameters().coordinates.of_configuration(&c).unwrap();
                    if matches!(cover.kind(zoom), Kind::Sheared { .. }) {
                        let no_fit = cover.parameters().coordinates.no_fit_coordinates();
                        assert!(no_fit.iter().all(|&i| x[i].is_zero()), "{} {}", pin.name, zoom);
                    } else {
                        let z = cover.cover_coordinates(&x);
                        assert!(z[k..].iter().all(QSqrt5::is_zero), "{} {}", pin.name, zoom);
                    }
                }
            }
        }
    }
}

#[test]
fn the_validity_filter_keeps_exactly_valid_supports_at_homogeneous_corner_views() {
    // Corner views of arc-plane zooms are homogeneous (u₃ = 2 + e ≠ 1), and
    // cells are tried with all their corner views at once. Whatever the
    // float filter drops must be exactly invalid at some corner view.
    let mut random = Random::new(72);
    let identity = Zoom::identity(ZoomCell(std::array::from_fn(|_| Interval::new(frac(0, 1), frac(1, 1)).unwrap()))).unwrap();
    let exactly_valid = |support: &Support, views: &[Point]| {
        views.iter().all(|u| {
            let n = support.normal_at(u);
            let h = geometry::dot(&n, &geometry::vertices()[support.contact()]);
            geometry::vertices().iter().all(|w| geometry::dot(&n, w) <= h)
        })
    };
    let (mut dropped, mut kept_valid) = (0, 0);
    for name in ["arc+", "crossing-", "square"] {
        let cover = ZoomCover::new(exotic_pin(name).unwrap().parameters.clone()).unwrap();
        for zoom in [0, cover.zooms().len() - 1] {
            let c = &cover.zooms()[zoom];
            for _ in 0..2 {
                let mut cell = c.root().clone();
                for _ in 0..4 + random.below(12) {
                    let axis = random.below(VARIABLES as u64) as usize;
                    let [lo, hi] = cell[axis].bisect();
                    cell[axis] = if random.below(2) == 0 { lo } else { hi };
                }
                let views: Vec<Point> = c.view_corners(&cell);
                let float_views: Vec<Vector> = views.iter().map(float::point).collect();
                for support in candidates::supports(&views[0]).into_iter().step_by(2) {
                    let gap = Witness::Gap(Gap::new(support.clone(), 0).unwrap());
                    let candidate = Candidate::new(&identity, gap).unwrap().unwrap();
                    if candidate.plausibly_valid(&float_views) {
                        if kept_valid < 50 && exactly_valid(&support, &views) {
                            kept_valid += 1;
                        }
                    } else {
                        assert!(!exactly_valid(&support, &views), "{name} zoom {zoom}: a valid support dropped");
                        dropped += 1;
                    }
                }
            }
        }
    }
    assert!(dropped > 100 && kept_valid > 10, "{dropped} {kept_valid}");
}

#[test]
fn degenerate_float_inputs_do_not_panic() {
    for u in [[0.0; 3], [f64::NAN; 3], [f64::INFINITY, 0.0, 1.0], [0.0, 0.0, 1e-300]] {
        for r in [[0.0; 3], [f64::NAN; 3], [1e300, 0.0, 0.0]] {
            let _ = candidates::separating(&u, &r, 8);
        }
    }
    assert!(float::hull(&[]).is_empty());
    // All points equal: documented as outside the assumptions; only no panic.
    let _ = float::hull(&[[1.0, 1.0]; 5]);
    assert_eq!(float::hull(&[[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]).len(), 2);
    let _ = float::hull(&[[f64::NAN, 0.0], [1.0, 0.0], [0.0, 1.0]]);
}

#[test]
fn the_aspect_cap_is_applied_to_exact_split_counts() {
    // The documented rule allows a variable exactly `aspect` times narrower
    // (relative to its root) than the widest one. Applied to rounded float
    // widths it would exclude the variable at 5 of 16 positions on this
    // root; with exact split counts it is allowed at every
    // position.
    use crate::arithmetic::exact::{q, QSqrt5};
    use crate::arithmetic::polynomial::Polynomial;
    use crate::elimination::generate::search::splits;
    let cover = ZoomCover::new(local_pin(21).unwrap().parameters.clone()).unwrap();
    let root = cover.zooms()[0].root().clone();
    let mut unit = [0u8; VARIABLES];
    unit[0] = 1;
    let along_0 = float::Dense::new(&Polynomial::monomial(unit, QSqrt5::integer(2))).bernstein(&[[0.0, 1.0]; VARIABLES]);
    let tried = vec![(-1.0, along_0)];
    for i in 0..16 {
        let mut cell = root.clone();
        let sixteenth = root[0].width() * frac(1, 16);
        let lo = root[0].lo() + &(&sixteenth * &q(i));
        cell[0] = Interval::new(lo.clone(), &lo + &sixteenth).unwrap();
        let axis = split_axis(&splits(&root, &cell), &tried, NonZeroU32::new(16).unwrap());
        assert_eq!(axis, 0, "position {i}");
    }
}

/// Once `data/` holds generated covers: the shipped configuration
/// reproduces them byte for byte (also across build profiles, since the data
/// is generated by a release build and this test runs in the test profile).
/// Ignored while `data/` holds covers from elsewhere:
///
/// ```text
/// cargo test --lib review_the_shipped_configuration_reproduces_the_data -- --ignored
/// ```
#[test]
#[ignore]
fn the_shipped_configuration_reproduces_the_data() {
    use crate::elimination::proof::catalogue::{exotic_data, local_data};
    let configuration = Configuration::parse(include_bytes!("../command/covers.json")).unwrap();
    let requests: Vec<Request> = configuration
        .requests()
        .unwrap()
        .into_iter()
        .map(|(r, _)| r)
        .filter(|r| ["arc-", "endpoint+", "15", "19", "29"].contains(&r.pin.name.as_str()))
        .collect();
    assert_eq!(requests.len(), 5);
    for (request, result) in requests.iter().zip(run(&requests, &configuration.search, 4)) {
        let bytes = result.unwrap().file.to_bytes();
        let data = match request.pin.name.parse::<usize>() {
            Ok(n) => local_data(n).unwrap(),
            Err(_) => exotic_data(&request.pin.name).unwrap(),
        };
        assert!(bytes == data, "cover {} differs from data/", request.pin.name);
    }
}

/// Every witness leaf of every cover in `data/` checked at random exact
/// points of its cell by the independent evaluation (no pull-back, no
/// Bernstein coefficients), on four threads:
///
/// ```text
/// cargo test --lib review_every_leaf_of_the_data_is_confirmed_at_random_points -- --ignored
/// ```
#[test]
#[ignore]
fn every_leaf_of_the_data_is_confirmed_at_random_points() {
    use crate::elimination::proof::catalogue::{exotic_data, local_data, EXOTIC};
    let mut covers: Vec<(crate::elimination::proof::Pin, &'static [u8])> =
        EXOTIC.iter().map(|n| (exotic_pin(n).unwrap(), exotic_data(n).unwrap())).collect();
    covers.extend((0..LOCAL_COUNT).map(|n| (local_pin(n).unwrap(), local_data(n).unwrap())));
    let next = std::sync::atomic::AtomicUsize::new(0);
    let checked = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| loop {
                let k = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                let Some((pin, bytes)) = covers.get(k) else { return };
                let file = CoverFile::parse(bytes).unwrap();
                let mut random = Random::new(1000 + k as u64);
                let n = check_independently(pin, &file, 3, &mut random);
                checked.fetch_add(n, std::sync::atomic::Ordering::SeqCst);
            });
        }
    });
    eprintln!("{} points checked", checked.into_inner());
}
