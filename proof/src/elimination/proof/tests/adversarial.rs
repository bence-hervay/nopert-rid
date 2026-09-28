//! Adversarial tests of proved covers:
//!
//! - a second verification of cover files written from the definitions only
//!   (zooms, pull-backs, division, Bernstein coefficients, validity, trees,
//!   window condition and `beyond` identity are an independent code; see
//!   `zoom::tests::independent`), sampled by default and complete under
//!   `#[ignore]`:
//!   `cargo test --release --lib cover::tests::adversarial::every_cell -- --ignored --nocapture --test-threads=4`
//! - end to end: every sampled configuration of every box accepted by
//!   `contains` (real search boxes of B₀ and small boxes around cover
//!   points, faces and ties) is proved at the exact point, through every
//!   zoom and every leaf whose closed cell contains its preimage (including
//!   hand-overs through the window);
//! - a cover whose cover has no zoom.
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::{Exponents, Polynomial};
use crate::elimination::zoom::cover::{Number, Parameters, Rational, Shape};
use crate::elimination::zoom::tests::independent::{self as model, rat, Family, ModelZoomCover, ModelZoom, P5};
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::Coordinates;
use crate::elimination::proof::catalogue::{exotic_data, exotic_pin, local_data, local_pin, EXOTIC, LOCAL_COUNT};
use crate::elimination::proof::format::{Leaf, CoverFile, Scope};
use crate::elimination::proof::{Pin, ProvedSet, ProvedCover};
use crate::elimination::witness::{DomainInequality, Witness};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::domain;
use crate::problem::geometry::View;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::num::NonZeroUsize;
use std::sync::atomic::AtomicBool;

fn covers() -> Vec<(Pin, &'static [u8])> {
    let mut out: Vec<(Pin, &'static [u8])> =
        EXOTIC.iter().map(|n| (exotic_pin(n).unwrap(), exotic_data(n).unwrap())).collect();
    out.extend((0..LOCAL_COUNT).map(|n| (local_pin(n).unwrap(), local_data(n).unwrap())));
    out
}

#[derive(Debug, Default, PartialEq, Eq)]
struct Counts {
    zooms: usize,
    witness_leaves: usize,
    domain_leaves: usize,
    delegated: usize,
    checked: usize,
}

/// The `beyond` identity from an independent coordinates: D's inequality
/// at `x(z) = c₀ + Mz` is `κ(z_axis − end)` with a constant `κ > 0`.
fn beyond_holds(m: &ModelZoomCover) -> bool {
    let Some(b) = m.parameters.beyond else { return true };
    let z = Polynomial::variables();
    let c = m.centre();
    let x: [Polynomial; 5] = std::array::from_fn(|i| {
        (0..5).fold(Polynomial::constant(c[i].clone()), |s, j| &s + &z[j].scale(&m.parameters.map[i][j].0))
    });
    let (u, r) = model::view_rotation(m.parameters.coordinates, &x);
    let ck = domain::polynomial(b.inequality, &View::Homogeneous(u), &r).unwrap();
    let mut unit = [0u8; 5];
    unit[b.axis] = 1;
    let kappa = ck.coefficient(&unit);
    let expected = &z[b.axis].scale(&kappa) - &Polynomial::constant(&kappa * &m.bounds[b.axis].1);
    kappa.sign() == Ordering::Greater && ck == expected
}

/// The window vector `w = (z_f + shear·z_b)/μ` of an A zoom, derived from
/// the zoom's own cover coordinates.
fn window_vector(m: &ModelZoomCover, zoom: &ModelZoom) -> Option<(Vec<Polynomial>, Vec<[i8; 2]>, Q)> {
    let Shape::Point { offset, window: Some(w), .. } = &m.parameters.shape else { return None };
    let k = m.base_dimension;
    let v: Vec<Polynomial> = (0..offset.len())
        .map(|l| {
            let sum = (0..k).fold(zoom.z[k + l].clone(), |s, j| &s + &zoom.z[j].scale(&rat(&w.shear[l][j].0)));
            model::divide(&sum, &[1, 0, 0, 0, 0]).expect("z is a multiple of μ")
        })
        .collect();
    Some((v, offset.clone(), w.radius.0.clone()))
}

fn window_faces(m: &ModelZoomCover) -> Vec<(usize, i8)> {
    match &m.parameters.shape {
        Shape::Point { window: Some(w), .. } => w.faces.iter().map(|f| (f.axis, f.side)).collect(),
        _ => Vec::new(),
    }
}

/// The second verification: returns the counts, panicking on the first
/// defect. Only witness leaves with `leaf % step == offset` are checked by
/// the lemma, and only jobs with `job % parts == part`.
fn verify_independently(pin: &Pin, bytes: &[u8], step: usize, offset: usize, part: (usize, usize)) -> Counts {
    let file = CoverFile::parse(bytes).unwrap();
    let name = &pin.name;
    assert_eq!(&file.name, name);
    assert_eq!(file.scope, pin.scope, "{name}");
    let mut counts = Counts::default();
    let mut models = Vec::new();
    // (cover, zoom, witness, factor) → cells.
    let mut jobs: BTreeMap<(usize, usize, usize, Exponents), Vec<Vec<(Q, Q)>>> = BTreeMap::new();
    {
        let (a, entry, pinned) = (0usize, &file, &pin.parameters);
        assert_eq!(&entry.parameters, pinned, "{name}");
        let m = ModelZoomCover::new(&entry.parameters);
        assert!(!m.zooms.is_empty(), "{name}: a cover without zooms covers nothing");
        assert!(pin.scope == Scope::Domain || m.parameters.beyond.is_none(), "{name}");
        assert!(beyond_holds(&m), "{name}: beyond");
        let names: Vec<&str> = m.zooms.iter().map(|c| c.name.as_str()).collect();
        let listed: Vec<&str> = entry.zooms.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, listed, "{name}: zooms");
        for (c, (zoom, data)) in m.zooms.iter().zip(&entry.zooms).enumerate() {
            counts.zooms += 1;
            let cells = model::tree_leaves(&zoom.root, &data.tree).expect("a complete tree");
            assert_eq!(cells.len(), data.leaves.len(), "{name} {}", zoom.name);
            for (leaf, (cell, e)) in cells.into_iter().zip(&data.leaves).enumerate() {
                match e {
                    Leaf::Witness { index, factor } => {
                        let w = &file.witnesses[*index];
                        if matches!(w, Witness::Domain(_)) {
                            assert_eq!(pin.scope, Scope::Domain, "{name}");
                            counts.domain_leaves += 1;
                        }
                        counts.witness_leaves += 1;
                        if leaf % step == offset {
                            jobs.entry((a, c, *index, *factor)).or_default().push(cell);
                        }
                    }
                    Leaf::Delegated => {
                        assert_eq!(zoom.family, Family::A, "{name} {}: only A zooms delegate", zoom.name);
                        assert!(window_faces(&m).contains(&zoom.base_face.unwrap()), "{name} {}", zoom.name);
                        let (w, offset, radius) = window_vector(&m, zoom).unwrap();
                        assert!(w.iter().all(model::is_multilinear));
                        for mask in 0..32 {
                            let y = model::corner(&cell, mask);
                            for (l, p) in w.iter().enumerate() {
                                let v = p.evaluate(&y);
                                let [lo, hi] = offset[l];
                                assert!(rat(&(&radius * q(i64::from(lo)))) <= v && v <= rat(&(&radius * q(i64::from(hi)))), "{name} {} leaf {leaf}", zoom.name);
                            }
                        }
                        counts.delegated += 1;
                    }
                }
            }
        }
        models.push(m);
    }
    for (job, ((a, c, index, factor), cells)) in jobs.into_iter().enumerate() {
        if job % part.1 != part.0 {
            continue;
        }
        let m = &models[a];
        let zoom = &m.zooms[c];
        let w = &file.witnesses[index];
        assert!((0..5).all(|j| factor[j] == 0 || zoom.factors.contains(&j)), "{name} {}", zoom.name);
        let (u, r) = model::view_rotation(m.parameters.coordinates, &zoom.x);
        assert!(u.iter().all(model::is_multilinear));
        let g = model::witness_polynomial(w, &u, &r);
        let n = model::divide(&g, &factor).unwrap_or_else(|| panic!("{name} {}: not divisible", zoom.name));
        let mut valid: HashSet<[QSqrt5; 3]> = HashSet::new();
        for cell in cells {
            assert!((0..5).all(|j| factor[j] == 0 || cell[j].0 >= Q::zero()), "{name} {}", zoom.name);
            for mask in 0..32 {
                let y = model::corner(&cell, mask);
                let view: [QSqrt5; 3] = std::array::from_fn(|i| u[i].evaluate(&y));
                assert_eq!(view[2].sign(), Ordering::Greater, "{name} {}", zoom.name);
                if !valid.contains(&view) {
                    assert!(model::valid_at(w, &view), "{name} {}: invalid at a corner", zoom.name);
                    valid.insert(view);
                }
            }
            assert!(model::bernstein_negative(&n, &cell), "{name} {}: {cell:?}", zoom.name);
            counts.checked += 1;
        }
    }
    counts
}

#[test]
fn a_sample_of_every_cover_verifies_independently() {
    let mut checked = 0;
    for (i, (pin, bytes)) in covers().into_iter().enumerate() {
        // Every 40th leaf (a different residue per cover), every 5th for
        // the small arc covers; the structure is checked completely.
        let step = if pin.name.contains("arc") || pin.name == "endpoint" { 5 } else { 40 };
        let counts = verify_independently(&pin, bytes, step, i % step, (0, 1));
        checked += counts.checked;
    }
    assert!(checked > 400, "{checked}");
}

#[test]
fn independent_counts_equal_the_loaded_reports() {
    let (exotic, local) = super::loaded();
    for (cover, (pin, bytes)) in exotic.iter().chain(local).zip(covers()) {
        assert_eq!(cover.name(), pin.name);
        let counts = verify_independently(&pin, bytes, usize::MAX, 1, (0, 1));
        let r = cover.report();
        assert_eq!(
            (r.zooms, r.witness_leaves, r.domain_leaves, r.delegated_leaves),
            (counts.zooms, counts.witness_leaves, counts.domain_leaves, counts.delegated),
            "{}",
            pin.name
        );
    }
}

fn every_cell(part: usize) {
    let start = std::time::Instant::now();
    let mut total = 0;
    for (pin, bytes) in covers() {
        let counts = verify_independently(&pin, bytes, 1, 0, (part, 4));
        total += counts.checked;
        eprintln!("part {part}: {} {} cells ({:.0?})", pin.name, counts.checked, start.elapsed());
    }
    eprintln!("part {part}: {total} cells verified independently in {:.0?}", start.elapsed());
}

#[test]
#[ignore]
fn every_cell_part_0() {
    every_cell(0);
}
#[test]
#[ignore]
fn every_cell_part_1() {
    every_cell(1);
}
#[test]
#[ignore]
fn every_cell_part_2() {
    every_cell(2);
}
#[test]
#[ignore]
fn every_cell_part_3() {
    every_cell(3);
}

// ---- End to end: points of accepted boxes are proved ----------------------------

/// The zoom variables of cover coordinates `z` in `zoom`, when they lie in
/// its root and map back exactly.
fn preimage(m: &ModelZoomCover, zoom: &ModelZoom, z: &P5) -> Option<P5> {
    let k = m.base_dimension;
    let (zb, zf) = (&z[..k], &z[k..]);
    let side = |s: i8| QSqrt5::integer(i64::from(s));
    fn free(v: &[QSqrt5], axis: usize) -> Vec<QSqrt5> {
        v.iter().enumerate().filter(|(j, _)| *j != axis).map(|(_, x)| x.clone()).collect()
    }
    let over = |v: &[QSqrt5], d: &QSqrt5| v.iter().map(|x| x * &model::inverse_of(d)).collect::<Vec<_>>();
    let (fa, fs) = zoom.offset_face;
    let delta = &side(fs) * &zf[fa];
    if delta.sign() != Ordering::Greater {
        return None;
    }
    let fhat = over(zf, &delta);
    let mut y: Vec<QSqrt5> = Vec::new();
    match zoom.family {
        Family::Tube => {
            y.extend(zb.iter().cloned());
            y.push(delta.clone());
            y.extend(free(&fhat, fa));
        }
        Family::A => {
            let (ga, gs) = zoom.base_face.unwrap();
            let mu = &side(gs) * &zb[ga];
            if mu.sign() != Ordering::Greater {
                return None;
            }
            y.push(mu.clone());
            y.push(&delta * &model::inverse_of(&mu));
            y.extend(free(&over(zb, &mu), ga));
            y.extend(free(&fhat, fa));
        }
        Family::B => {
            y.push(delta.clone());
            y.extend(over(zb, &delta));
            y.extend(free(&fhat, fa));
        }
        Family::Sheared => return None,
    }
    let y: P5 = y.try_into().ok()?;
    let inside = (0..5).all(|j| rat(&zoom.root[j].0) <= y[j] && y[j] <= rat(&zoom.root[j].1));
    (inside && model::evaluate_all(&zoom.z, &y) == z.to_vec()).then_some(y)
}

fn cell_contains(cell: &[(Q, Q)], y: &P5) -> bool {
    (0..5).all(|j| rat(&cell[j].0) <= y[j] && y[j] <= rat(&cell[j].1))
}

fn in_no_fit_set(coordinates: Coordinates, x: &[QSqrt5]) -> bool {
    model::no_fit_coordinates(coordinates).iter().all(|&i| x[i].is_zero())
}

/// Proves the configuration with zoom variables `y` of `zoom` at the
/// exact point, through every leaf containing it. Returns the number of
/// leaves used.
fn prove_in_zoom(file: &CoverFile, a: usize, m: &ModelZoomCover, c: usize, y: &P5, depth: usize) -> usize {
    let zoom = &m.zooms[c];
    let coordinates = m.parameters.coordinates;
    let x = model::evaluate_all(&zoom.x, y);
    let x5: P5 = x.clone().try_into().unwrap();
    debug_assert_eq!(a, 0);
    let data = &file.zooms[c];
    let cells = model::tree_leaves(&zoom.root, &data.tree).unwrap();
    let mut used = 0;
    for (cell, entry) in cells.iter().zip(&data.leaves) {
        if !cell_contains(cell, y) {
            continue;
        }
        used += 1;
        match entry {
            Leaf::Witness { index, factor } => {
                let monomial = (0..5).fold(QSqrt5::one(), |p, j| (0..factor[j]).fold(p, |p, _| &p * &y[j]));
                let (u, _) = model::exact_view_rotation(coordinates, &x5);
                assert_eq!(u[2].sign(), Ordering::Greater);
                if monomial.is_zero() {
                    assert!(in_no_fit_set(coordinates, &x), "{}: a vanishing factor off the no-fit set", zoom.name);
                } else {
                    let w = &file.witnesses[*index];
                    assert_eq!(model::witness_value(w, coordinates, &x5).sign(), Ordering::Less, "{} at {y:?}", zoom.name);
                    assert!(model::valid_at(w, &u), "{}", zoom.name);
                }
            }
            Leaf::Delegated => {
                assert_eq!(depth, 0, "a sheared zoom delegates");
                let (w, offset, radius) = window_vector(m, zoom).unwrap();
                let w: Vec<QSqrt5> = w.iter().map(|p| p.evaluate(y)).collect();
                // The sheared family's gauge and face.
                let (mut size, mut face) = (QSqrt5::zero(), None);
                for (l, v) in w.iter().enumerate() {
                    let [lo, hi] = offset[l];
                    assert!(rat(&(&radius * q(i64::from(lo)))) <= *v && *v <= rat(&(&radius * q(i64::from(hi)))));
                    for s in [-1i8, 1] {
                        let a = v * &QSqrt5::integer(i64::from(s));
                        if a > size {
                            size = a;
                            face = Some((l, s));
                        }
                    }
                }
                let Some((fa, fs)) = face else {
                    assert!(in_no_fit_set(coordinates, &x), "{}: w = 0 off the no-fit set", zoom.name);
                    continue;
                };
                let target = m
                    .zooms
                    .iter()
                    .position(|t| t.family == Family::Sheared && t.base_face == zoom.base_face && t.offset_face == (fa, fs))
                    .expect("a sheared zoom");
                let k = m.base_dimension;
                let mut y2: Vec<QSqrt5> = vec![y[0].clone(), size.clone()];
                y2.extend(y[2..1 + k].iter().cloned());
                y2.extend(w.iter().enumerate().filter(|(l, _)| *l != fa).map(|(_, v)| v * &model::inverse_of(&size)));
                let y2: P5 = y2.try_into().unwrap();
                let t = &m.zooms[target];
                assert!((0..5).all(|j| rat(&t.root[j].0) <= y2[j] && y2[j] <= rat(&t.root[j].1)), "hand-over leaves the sheared root");
                assert_eq!(model::evaluate_all(&t.x, &y2), x, "hand-over moves the configuration");
                assert!(prove_in_zoom(file, a, m, target, &y2, depth + 1) > 0);
            }
        }
    }
    used
}

/// Proves configuration `c` (exact) of a cover; `false` when `c` is
/// not in this covered set.
fn prove_point(file: &CoverFile, a: usize, m: &ModelZoomCover, c: &P5) -> bool {
    let coordinates = m.parameters.coordinates;
    let Some(x) = model::coordinates_of(coordinates, c) else { return false };
    let z = m.cover_coordinates(&x);
    if !m.z_inside(&z) {
        return false;
    }
    if let Some(b) = m.parameters.beyond {
        if z[b.axis] > m.bounds[b.axis].1 {
            // Outside D: its inequality is violated at c.
            let w = Witness::Domain(DomainInequality::new(b.inequality).unwrap());
            assert_eq!(model::witness_value(&w, coordinates, &x).sign(), Ordering::Less);
            return true;
        }
    }
    let k = m.base_dimension;
    if z[k..].iter().all(QSqrt5::is_zero) {
        assert!(in_no_fit_set(coordinates, &x));
        return true;
    }
    let mut zooms = 0;
    for (ci, zoom) in m.zooms.iter().enumerate() {
        if let Some(y) = preimage(m, zoom, &z) {
            zooms += 1;
            assert!(prove_in_zoom(file, a, m, ci, &y, 0) > 0, "no leaf contains the preimage");
        }
    }
    assert!(zooms > 0, "a cover point in no zoom root: z = {z:?}");
    true
}

fn near(x: &QSqrt5) -> Q {
    x.rational_part() + &(x.sqrt5_part() * frac(2_236_067_977, 1_000_000_000))
}

/// The cyclic-bisection search box of depth `depth` containing `c`.
fn search_box(c: &[Q; 5], depth: usize, random: &mut Random) -> ConfigurationBox {
    let mut b = ConfigurationBox::root();
    for d in 0..depth {
        let [lower, upper] = b.split(d % 5);
        let in_lower = lower.axes()[d % 5].contains(&c[d % 5]);
        let in_upper = upper.axes()[d % 5].contains(&c[d % 5]);
        b = match (in_lower, in_upper) {
            (true, true) => if random.below(2) == 0 { lower } else { upper },
            (true, false) => lower,
            (false, true) => upper,
            (false, false) => return b,
        };
    }
    b
}

#[test]
fn every_point_of_accepted_boxes_is_proved() {
    let (exotic, local) = super::loaded();
    let mut random = Random::new(0x5eed);
    let all: Vec<(&ProvedCover, (Pin, &'static [u8]))> = exotic.iter().chain(local).zip(covers()).collect();
    for (cover, (pin, bytes)) in all.into_iter().filter(|(r, _)| r.name().parse::<usize>().map_or(true, |n| n % 8 == 0)) {
        let file = CoverFile::parse(bytes).unwrap();
        let models: Vec<ModelZoomCover> = vec![ModelZoomCover::new(&pin.parameters)];
        let (mut accepted, mut points) = (0, 0);
        for round in 0..60 {
            let a = random.below(models.len() as u64) as usize;
            let m = &models[a];
            let z: P5 = std::array::from_fn(|i| {
                let (lo, hi) = &m.bounds[i];
                let t = match random.below(3) {
                    0 => frac(random.below(5) as i64, 4),
                    _ => random.between(&q(0), &q(1), 8),
                };
                lo + &(&(hi - lo) * &rat(&t))
            });
            let c = model::configuration_of(m.parameters.coordinates, &m.image(&z));
            let c: [Q; 5] = std::array::from_fn(|j| near(&c[j]));
            let b = if round % 2 == 0 {
                search_box(&c, 10 + random.below(30) as usize, &mut random)
            } else {
                let w = frac(1, 1 << (6 + random.below(10)));
                ConfigurationBox::new(std::array::from_fn(|j| {
                    let (w1, w2) = (&w * random.between(&q(0), &q(1), 3), &w * random.between(&q(0), &q(1), 3));
                    Interval::new(&c[j] - &w1, &c[j] + &w2).unwrap()
                }))
            };
            if !cover.contains(&b) {
                continue;
            }
            accepted += 1;
            let mut samples: Vec<P5> = (0..32u8).map(|mask| b.corner(mask).map(|v| rat(&v))).collect();
            let axes: [Interval; 5] = b.axes().clone();
            for _ in 0..6 {
                samples.push(random.point(&axes, 8).map(|v| rat(&v)));
            }
            for s in &samples {
                let proved = models.iter().enumerate().filter(|(a, m)| prove_point(&file, *a, m, s)).count();
                assert!(proved > 0, "{}: a point of an accepted box is in no cover", pin.name);
                points += 1;
            }
        }
        eprintln!("{}: {accepted} accepted boxes, {points} points proved", pin.name);
        assert!(accepted > 3, "{}: {accepted}", pin.name);
    }
}

// ---- Degenerate covers -----------------------------------------------------

/// DEFECT (soundness, latent): a pinned cover with an empty offset has no
/// zoom; its cover file lists none, verification succeeds and `contains`
/// accepts boxes that nothing proved.
#[test]
fn a_cover_whose_cover_has_no_zoom_is_refused() {
    let unit = || [Rational(q(0)), Rational(q(1))];
    let parameters = Parameters {
        coordinates: Coordinates::Configuration,
        centre: std::array::from_fn(|_| Number(QSqrt5::zero())),
        map: std::array::from_fn(|i| std::array::from_fn(|j| Number(QSqrt5::integer(i64::from(i == j))))),
        shape: Shape::Tube { base: vec![unit(), unit(), unit(), unit(), unit()], offset: vec![], radius: Rational(q(1)) },
        beyond: None,
    };
    let pin = Pin { name: "empty".into(), scope: Scope::All, parameters: parameters.clone() };
    let file = CoverFile {
        format: crate::elimination::proof::format::Format::Version1,
        name: "empty".into(),
        scope: Scope::All,
        parameters,
        witnesses: vec![],
        zooms: vec![],
    };
    let stop = AtomicBool::new(false);
    if let Ok(cover) = ProvedCover::load(&pin, &file.to_bytes(), NonZeroUsize::new(1).unwrap(), &stop) {
        let b = ConfigurationBox::new(std::array::from_fn(|j| {
            if j < 2 { Interval::new(q(0), frac(1, 4)).unwrap() } else { Interval::new(frac(1, 4), frac(2, 5)).unwrap() }
        }));
        panic!("a cover without zooms loaded; contains(B) = {}", cover.contains(&b));
    }
}

// ---- Measurements -----------------------------------------------------------

/// Cost of `contains` over every cover for random search boxes of B₀, and
/// of the tiling check:
/// `cargo test --release --lib cover::tests::adversarial::measure -- --ignored --nocapture --test-threads=1`
#[test]
#[ignore]
fn measure_contains_and_tiling() {
    let (exotic, local) = super::loaded();
    let mut random = Random::new(99);
    let root = ConfigurationBox::root();
    let boxes: Vec<ConfigurationBox> = (0..2000)
        .map(|_| {
            let c: [Q; 5] = std::array::from_fn(|j| random.between(root.axes()[j].lo(), root.axes()[j].hi(), 20));
            search_box(&c, 10 + random.below(40) as usize, &mut random)
        })
        .collect();
    for (label, covers) in [("exotic", exotic), ("local", local)] {
        let start = std::time::Instant::now();
        let mut hits = 0;
        for b in &boxes {
            hits += covers.iter().filter(|r| r.contains(b)).count();
        }
        let per_box = start.elapsed().as_secs_f64() / boxes.len() as f64;
        eprintln!("{label}: {} covers, {:.1} µs per box for all of them ({hits} hits)", covers.len(), per_box * 1e6);
    }
    for (pin, bytes) in covers().into_iter().take(EXOTIC.len()) {
        let file = CoverFile::parse(bytes).unwrap();
        let start = std::time::Instant::now();
        let (mut cells, mut tiling) = (0, 0.0);
        {
            let entry = &file;
            let cover = crate::elimination::zoom::cover::ZoomCover::new(entry.parameters.clone()).unwrap();
            for (zoom, data) in cover.zooms().iter().zip(&entry.zooms) {
                let leaves = crate::elimination::zoom::tree::leaves(zoom.root(), &data.tree).unwrap();
                cells += leaves.len();
                let t = std::time::Instant::now();
                crate::elimination::zoom::tree::check_tiling(zoom.root(), &leaves).unwrap();
                tiling += t.elapsed().as_secs_f64();
            }
        }
        eprintln!("{}: covers, trees and tiling of {cells} cells in {:.3} s, of which tiling {tiling:.3} s", pin.name, start.elapsed().as_secs_f64());
    }
}
