//! Test support shared by the components' tests, and the tests of the
//! independent oracle itself against the crate's exact model.
pub(crate) mod fixtures;
pub(crate) mod oracle;

use crate::arithmetic::exact::{frac, QSqrt5};
use crate::arithmetic::polynomial::VARIABLES;
use crate::problem::domain;
use crate::problem::geometry;
use crate::components::record::RecordData;
use crate::problem::configuration::ConfigurationBox;
use fixtures::Random;
use oracle::{Configuration, MARGIN};

fn geometry_vertices() -> Vec<oracle::V3> {
    geometry::vertices()
        .iter()
        .map(|v| std::array::from_fn(|j| oracle::number(&v[j])))
        .collect()
}

fn symmetries() -> Vec<[f64; 4]> {
    geometry::symmetries()
        .iter()
        .map(|g| std::array::from_fn(|j| oracle::number(&g[j])))
        .collect()
}

#[test]
fn the_oracle_and_the_geometry_have_the_same_solid_and_group() {
    let own = oracle::vertices();
    for v in geometry_vertices() {
        assert!(own.iter().any(|w| (0..3).all(|j| (v[j] - w[j]).abs() < 1e-12)));
    }
    for g in oracle::group() {
        for v in own {
            let x = oracle::product(&oracle::product(g, &[0.0, v[0], v[1], v[2]]), &[g[0], -g[1], -g[2], -g[3]]);
            assert!(own.iter().any(|w| (0..3).all(|j| (x[j + 1] - w[j]).abs() < 1e-9)));
        }
    }
    for g in symmetries() {
        let found = oracle::group().iter().any(|h| (0..4).all(|j| (g[j] - h[j]).abs() < 1e-12));
        assert!(found, "{g:?}");
    }
}

#[test]
fn the_oracles_inequalities_equal_the_domain_polynomials() {
    let mut random = Random(5);
    let symmetries = symmetries();
    for _ in 0..40 {
        let exact: [QSqrt5; VARIABLES] = std::array::from_fn(|j| {
            let bound = if j < 2 { 2 } else { 1 };
            let k = random.below(2001) as i64 - 1000;
            QSqrt5::from_rational(frac(k * bound, 2500))
        });
        let c = Configuration {
            s: oracle::number(&exact[0]),
            t: oracle::number(&exact[1]),
            r: std::array::from_fn(|j| oracle::number(&exact[j + 2])),
        };
        let values = oracle::inequalities(&c, &symmetries);
        for (i, p) in domain::affine_polynomials().iter().enumerate() {
            let expected = oracle::number(&p.evaluate(&exact));
            assert!((values[i] - expected).abs() < 1e-9, "inequality {i}");
        }
    }
}

#[test]
fn the_oracle_sees_touching_shadows_at_aligned_and_symmetric_rotations() {
    let mut random = Random(6);
    for _ in 0..50 {
        let s = random.below(1000) as f64 / 1500.0;
        let t = random.below(1000) as f64 / 2500.0;
        let aligned = Configuration { s, t, r: [0.0; 3] };
        assert!(oracle::excursion(&aligned).abs() < 1e-9);
        for g in oracle::group().iter().filter(|g| g[0].abs() > 0.5) {
            let r = [g[1] / g[0], g[2] / g[0], g[3] / g[0]];
            let symmetric = Configuration { s, t, r };
            assert!(oracle::excursion(&symmetric).abs() < 1e-9);
        }
        let turned = Configuration { s, t, r: [0.3, -0.2, 0.1] };
        assert!(oracle::excursion(&turned) > 1e-3);
    }
}

/// The ordered edge and plug vertex of a Global record.
pub(crate) fn parts(data: &RecordData) -> (usize, usize, usize) {
    let RecordData::Global(maximum) = data else { panic!("not a Global record") };
    (maximum.edge().from(), maximum.edge().to(), maximum.plug())
}

pub(crate) fn vertex(i: usize) -> oracle::V3 {
    std::array::from_fn(|j| oracle::number(&geometry::vertices()[i][j]))
}

/// What the oracle says about the claim of a Global record, that the plug
/// vertex lies beyond the hole's support in the direction `n(u)` of the
/// ordered edge throughout the box (the maximum form's claim): `Some(true)`
/// if it holds at every sample within the margin, `Some(false)` if some
/// sample clearly contradicts it, `None` in between.
pub(crate) fn oracle_verdict(b: &ConfigurationBox, (from, to, plug): (usize, usize, usize), seed: u64) -> Option<bool> {
    let (f, t, p) = (vertex(from), vertex(to), vertex(plug));
    let samples = oracle::samples(b, seed, 12);
    let beyond = |c: &Configuration| oracle::plug_beyond(c, &f, &t, &p) - oracle::hole_beyond(c, &f, &t);
    let clearly_false = samples.iter().any(|c| beyond(c) < -MARGIN);
    let holds = samples.iter().all(|c| beyond(c) >= -MARGIN);
    match (clearly_false, holds) {
        (true, _) => Some(false),
        (false, true) => Some(true),
        (false, false) => None,
    }
}

/// The oracle's confirmation of an accepted Domain or Global record: a Domain
/// inequality is violated at every sample, which also lies outside `D` by the
/// oracle's own group; a Global record's plug vertex lies beyond the hole's
/// support in the edge's normal direction at every sample, and every sample
/// is a non-fit by the oracle's convex hulls.
pub(crate) fn confirm(b: &ConfigurationBox, data: &RecordData, seed: u64) {
    let samples = oracle::samples(b, seed, 12);
    match data {
        RecordData::Domain { inequality } => {
            let symmetries = symmetries();
            for c in &samples {
                assert!(oracle::inequalities(c, &symmetries)[inequality.index()] > -MARGIN, "{data:?} at {c:?}");
                assert!(oracle::outside_domain(c) > -MARGIN, "{c:?} is in D");
            }
        }
        RecordData::Global(_) => {
            assert_eq!(oracle_verdict(b, parts(data), seed), Some(true), "{data:?}");
            for c in &samples {
                assert!(oracle::excursion(c) > -MARGIN, "a fit at {c:?}");
            }
        }
        _ => panic!("the oracle knows no covers: {data:?}"),
    }
}

/// 12,000 random boxes of depths 10 to 60, most of them meeting `D`, checked
/// by Domain and Global on four threads, every witness confirmed by the
/// oracle.
///
/// cargo test --release --lib components::tests::heavy_random_campaign -- --ignored --nocapture
#[test]
#[ignore]
fn heavy_random_campaign() {
    let collection = fixtures::bare();
    let counts = std::sync::Mutex::new([0usize; 3]);
    std::thread::scope(|scope| {
        for seed in 0..4u64 {
            let (collection, counts) = (&collection, &counts);
            scope.spawn(move || {
                let mut random = Random(1000 + seed);
                let mut local = [0usize; 3];
                for n in 0..3000 {
                    // A random descent that prefers children whose midpoint
                    // lies in D, so that most boxes meet D.
                    let depth = 10 + random.below(51);
                    let mut path = String::new();
                    while path.len() < depth {
                        let inside: Vec<char> = ['0', '1']
                            .into_iter()
                            .filter(|bit| {
                                let child = fixtures::from_path(&format!("{path}{bit}"));
                                oracle::outside_domain(&oracle::samples(&child, 0, 0)[32]) < 0.0
                            })
                            .collect();
                        let choice = if inside.is_empty() { ['0', '1'][random.below(2)] } else { inside[random.below(inside.len())] };
                        path.push(choice);
                    }
                    let b = fixtures::from_path(&path);
                    match collection.check(&b) {
                        None => local[0] += 1,
                        Some(data) => {
                            confirm(&b, &data, n);
                            local[1 + usize::from(matches!(data, RecordData::Global(_)))] += 1;
                        }
                    }
                }
                let mut counts = counts.lock().unwrap();
                for j in 0..3 {
                    counts[j] += local[j];
                }
            });
        }
    });
    let [refused, domain, global] = *counts.lock().unwrap();
    println!("{refused} refused, {domain} Domain and {global} Global records, all confirmed");
    assert!(domain > 0 && global > 0);
}
