//! Tests of proved covers: randomised coverage through verified cells, the
//! inclusion test against exact sampling of search boxes, refusal of
//! corrupted files and pins, determinism across thread counts, and an
//! independent second verification of the cells.
use super::catalogue::{self, exotic_data, exotic_pin, local_data, local_pin, EXOTIC, LOCAL_COUNT};
use super::format::{Leaf, CoverFile, Scope};
use super::*;
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::Polynomial;
use crate::elimination::zoom::cover::{Location, Number, Rational, Shape};
use crate::elimination::zoom::tests::Random;
use crate::arithmetic::polynomial::Exponents;
use crate::elimination::zoom::{reciprocal, tree, Five, Zoom, ZoomCell, ZoomError, ZoomFactor};
use crate::elimination::witness::DomainInequality;
use crate::problem::configuration::{split_axis, ConfigurationBox};
use crate::problem::geometry::{self, Point};
use crate::problem::domain;
use std::collections::HashMap;
use std::sync::OnceLock;

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn running() -> AtomicBool {
    AtomicBool::new(false)
}

/// Every cover of the catalogue, loaded once for all tests.
fn loaded() -> &'static (Vec<ProvedCover>, Vec<ProvedCover>) {
    static LOADED: OnceLock<(Vec<ProvedCover>, Vec<ProvedCover>)> = OnceLock::new();
    LOADED.get_or_init(|| {
        let stop = running();
        (catalogue::exotic(threads(4), &stop).unwrap(), catalogue::local(threads(4), &stop).unwrap())
    })
}

#[test]
fn every_cover_of_the_catalogue_verifies() {
    let (exotic, local) = loaded();
    // The counts of the data that `rid generate` writes from the shipped
    // configuration.
    let expected = [
        ("square", 877, 95, 0),
        ("pentagon", 6173, 0, 0),
        ("arc+", 343, 32, 0),
        ("arc-", 288, 72, 0),
        ("endpoint+", 305, 10, 0),
        ("endpoint-", 246, 50, 0),
        ("crossing+", 1431, 31, 34),
        ("crossing-", 687, 160, 14),
    ];
    for (cover, (name, leaves, domain, delegated)) in exotic.iter().zip(expected) {
        assert_eq!(cover.name(), name);
        assert_eq!(cover.scope(), Scope::Domain);
        let r = cover.report();
        assert_eq!((r.witness_leaves, r.domain_leaves, r.delegated_leaves), (leaves, domain, delegated), "{name}");
    }
    assert_eq!(local.len(), LOCAL_COUNT);
    let mut leaves = 0;
    for (n, cover) in local.iter().enumerate() {
        assert_eq!(cover.name(), n.to_string());
        assert_eq!(cover.scope(), Scope::All);
        assert_eq!(cover.report().domain_leaves, 0);
        assert_eq!(cover.report().zooms, 6);
        leaves += cover.report().witness_leaves;
    }
    assert_eq!(leaves, 833);
}

// ---- Coverage through verified cells ---------------------------------------

/// The leaves of every zoom of a cover, with their entries.
struct Cells {
    file: CoverFile,
    leaves: Vec<Vec<ZoomCell>>,
}

impl Cells {
    fn new(cover: &ProvedCover, bytes: &[u8]) -> Self {
        let file = CoverFile::parse(bytes).unwrap();
        let leaves =
            cover.cover().zooms().iter().zip(&file.zooms).map(|(z, e)| tree::leaves(z.root(), &e.tree).unwrap()).collect();
        Self { file, leaves }
    }
}

fn contains_point(cell: &ZoomCell, y: &Five<QSqrt5>) -> bool {
    (0..5).all(|j| QSqrt5::from_rational(cell[j].lo().clone()) <= y[j] && y[j] <= QSqrt5::from_rational(cell[j].hi().clone()))
}

fn exact(x: &Five<Q>) -> Five<QSqrt5> {
    std::array::from_fn(|j| QSqrt5::from_rational(x[j].clone()))
}

/// The configuration `(s, t, r)` of zoom coordinates `x`.
fn configuration(zoom: &Zoom, x: &Five<QSqrt5>) -> Five<QSqrt5> {
    let constants: Five<Polynomial> = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (view, rotation) = zoom.coordinates().configuration(&constants);
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    let u = view.vector().map(|p| p.evaluate(&origin));
    let r = rotation.map(|p| p.evaluate(&origin));
    let inverse = reciprocal(&u[2]);
    [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()]
}

/// Follows a configuration of the cover into a verified cell and checks
/// the lemma's conclusion there at the exact point: the divided monomial
/// vanishes (and the configuration is in the no-fit set) or the quotient is
/// negative, and the support is valid at the exact view.
fn follow(cover: &ProvedCover, cells: &Cells, c: &Five<QSqrt5>, quotients: &mut HashMap<(usize, usize, Exponents), Polynomial>) -> &'static str {
    let a = cover.cover();
    let x = a.parameters().coordinates.of_configuration(c).unwrap();
    let mut location = a.locate(&x).expect("a cover configuration is located");
    loop {
        match location {
            Location::NoFit => {
                let nofit = a.parameters().coordinates.no_fit_coordinates();
                assert!(nofit.iter().all(|&i| x[i].is_zero()));
                return "no-fit";
            }
            Location::OutsideDomain => {
                let inequality = a.parameters().beyond.unwrap().inequality;
                let value = domain::affine_polynomials()[inequality].evaluate(c);
                assert_eq!(value.sign(), std::cmp::Ordering::Greater);
                return "outside";
            }
            Location::Zoom { zoom, y } => {
                let ch = &a.zooms()[zoom];
                assert_eq!(configuration(ch, &ch.evaluate(&y)), *c);
                let leaf = cells.leaves[zoom].iter().position(|cell| contains_point(cell, &y)).expect("in a leaf");
                match &cells.file.zooms[zoom].leaves[leaf] {
                    Leaf::Delegated => {
                        location = a.hand_over(zoom, &y).expect("the window holds");
                        continue;
                    }
                    Leaf::Witness { index, factor } => {
                        let witness = &cells.file.witnesses[*index];
                        let n = quotients
                            .entry((zoom, *index, *factor))
                            .or_insert_with(|| {
                                let factor = ZoomFactor::new(*factor);
                                crate::elimination::lemma::Quotient::new(ch, witness, &factor).unwrap().polynomial().clone()
                            });
                        let monomial = (0..5).fold(QSqrt5::one(), |p, j| (0..factor[j]).fold(p, |p, _| &p * &y[j]));
                        if monomial.is_zero() {
                            let nofit = ch.coordinates().no_fit_coordinates();
                            let xy = ch.evaluate(&y);
                            assert!(nofit.iter().all(|&i| xy[i].is_zero()));
                            return "factor";
                        }
                        assert_eq!(n.evaluate(&y).sign(), std::cmp::Ordering::Less);
                        if let Witness::Gap(gap) = witness {
                            let u = ch.view().vector().map(|p| p.evaluate(&y));
                            let normal = gap.support().normal_at(&u);
                            let h = geometry::dot(&normal, &geometry::vertices()[gap.support().contact()]);
                            assert!(geometry::vertices().iter().all(|w| geometry::dot(&normal, w) <= h));
                        }
                        return "cell";
                    }
                }
            }
        }
    }
}

fn cover_data(name: &str) -> &'static [u8] {
    match name.parse::<usize>() {
        Ok(n) => local_data(n).unwrap(),
        Err(_) => exotic_data(name).unwrap(),
    }
}

/// A random configuration of a cover (often on faces and ties).
fn random_configuration(random: &mut Random, cover: &ProvedCover) -> Five<QSqrt5> {
    let a = cover.cover();
    let z: Five<QSqrt5> = std::array::from_fn(|i| {
        let [lo, hi] = &a.bounds()[i];
        let t = if random.below(4) == 0 { frac(random.below(5) as i64, 4) } else { random.between(&q(0), &q(1), 16) };
        lo + &(&(hi - lo) * &QSqrt5::from_rational(t))
    });
    let p = a.parameters();
    let x: Five<QSqrt5> = std::array::from_fn(|i| (0..5).fold(p.centre[i].0.clone(), |s, j| &s + &(&p.map[i][j].0 * &z[j])));
    configuration(&a.zooms()[0], &x)
}

#[test]
fn random_configurations_of_a_cover_lie_in_proved_cells() {
    let (exotic, local) = loaded();
    let mut random = Random::new(41);
    let chosen: Vec<&ProvedCover> = exotic.iter().chain([&local[0], &local[14], &local[29]]).collect();
    for cover in chosen {
        let cells = Cells::new(cover, cover_data(cover.name()));
        let mut quotients = HashMap::new();
        let mut outcomes: HashMap<&str, usize> = HashMap::new();
        let rounds = if cover.name().starts_with("crossing") { 300 } else { 75 };
        for _ in 0..rounds {
            let c = random_configuration(&mut random, cover);
            *outcomes.entry(follow(cover, &cells, &c, &mut quotients)).or_default() += 1;
        }
        assert!(outcomes.get("cell").copied().unwrap_or(0) > rounds / 2, "{}: {outcomes:?}", cover.name());
    }
}

#[test]
fn configurations_near_the_second_arc_are_handed_over() {
    let (exotic, _) = loaded();
    let mut random = Random::new(42);
    let mut handed = 0;
    for round in 0..300 {
        let crossing = &exotic[6 + round % 2];
        let cells = Cells::new(crossing, exotic_data(crossing.name()).unwrap());
        let mut quotients = HashMap::new();
        // θ close to −e with e > 0: the second arc, where zoom A delegates.
        let a = crossing.cover();
        let sign = match a.parameters().coordinates {
            crate::elimination::zoom::Coordinates::ArcPlane(s) => s,
            _ => unreachable!(),
        };
        // Offset offsets relative to e keep the zoom point near ρ = 1, f̂ = −e_θ.
        let e = random.between(&frac(1, 256), &frac(1, 32), 12);
        let theta = -&e * (q(1) + random.between(&frac(-1, 32), &frac(1, 32), 8));
        let small = |r: &mut Random| &e * r.between(&frac(-1, 32), &frac(1, 32), 8);
        let (zeta1, zeta3) = (small(&mut random), small(&mut random));
        let v = &e * random.between(&q(0), &frac(1, 32), 8);
        let x: Five<QSqrt5> = exact(&[e.clone(), v, theta, zeta1, zeta3]);
        let c = configuration(&a.zooms()[0], &x);
        assert_eq!(crate::elimination::zoom::Coordinates::ArcPlane(sign).of_configuration(&c).unwrap(), x);
        if let Some(Location::Zoom { zoom, y }) = a.locate(&x) {
            let leaf = cells.leaves[zoom].iter().position(|cell| contains_point(cell, &y)).unwrap();
            if cells.file.zooms[zoom].leaves[leaf] == Leaf::Delegated {
                handed += 1;
            }
        }
        follow(crossing, &cells, &c, &mut quotients);
    }
    assert!(handed > 50, "{handed}");
}

// ---- The inclusion test ------------------------------------------------------

/// The search box of depth `depth` (cyclic midpoint bisection of B₀) that
/// contains the configuration `c`, preferring the lower half on ties.
fn search_box(c: &Five<QSqrt5>, depth: usize) -> ConfigurationBox {
    let mut b = ConfigurationBox::root();
    for d in 0..depth {
        let axis = split_axis(d);
        let [lower, upper] = b.split(axis);
        b = if c[axis] <= QSqrt5::from_rational(lower.axes()[axis].hi().clone()) { lower } else { upper };
    }
    b
}

fn in_root(c: &Five<QSqrt5>) -> bool {
    let root = ConfigurationBox::root();
    (0..5).all(|j| {
        QSqrt5::from_rational(root.axes()[j].lo().clone()) <= c[j] && c[j] <= QSqrt5::from_rational(root.axes()[j].hi().clone())
    })
}

/// Whether a configuration lies in the cover, by its cover coordinates.
fn point_inside(cover: &ProvedCover, c: &Five<QSqrt5>) -> bool {
    let a = cover.cover();
    let Some(x) = a.parameters().coordinates.of_configuration(c) else { return false };
    let z = a.cover_coordinates(&x);
    (0..5).all(|i| {
        let open = a.parameters().beyond.is_some_and(|b| b.axis == i);
        a.bounds()[i][0] <= z[i] && (open || z[i] <= a.bounds()[i][1])
    })
}

#[test]
fn inclusion_test_agrees_with_exact_sampling_of_search_boxes() {
    let (exotic, local) = loaded();
    let mut random = Random::new(43);
    for cover in exotic.iter().chain([&local[5], &local[19]]) {
        let (mut accepted, mut refused) = (0, 0);
        for _ in 0..200 {
            let c = random_configuration(&mut random, cover);
            if !in_root(&c) {
                continue;
            }
            let b = search_box(&c, 10 + random.below(40) as usize);
            let corners: Vec<Five<QSqrt5>> = (0..32u8).map(|m| exact(&b.corner(m))).collect();
            let by_cover = corners.iter().all(|k| point_inside(cover, k));
            let contains = cover.contains(&b);
            assert_eq!(contains, by_cover, "{}", cover.name());
            if contains {
                accepted += 1;
                for _ in 0..6 {
                    assert!(point_inside(cover, &exact(&random.point(b.axes(), 10))));
                }
            } else {
                // Refusal is exact: some corner lies outside.
                refused += 1;
            }
        }
        assert!(accepted > 10 && refused > 10, "{}: {accepted} {refused}", cover.name());
    }
}

// ---- Refusals ---------------------------------------------------------------

fn load(pin: &Pin, bytes: &[u8]) -> Result<ProvedCover, ProofError> {
    ProvedCover::load(pin, bytes, threads(4), &running())
}

fn edited(name: &str, change: impl Fn(&mut CoverFile)) -> Vec<u8> {
    let mut file = CoverFile::parse(cover_data(name)).unwrap();
    change(&mut file);
    file.to_bytes()
}

fn pin(name: &str) -> Pin {
    match name.parse::<usize>() {
        Ok(n) => local_pin(n).unwrap(),
        Err(_) => exotic_pin(name).unwrap(),
    }
}

#[test]
fn malformed_and_noncanonical_files_are_refused() {
    let bytes = cover_data("arc-");
    let p = pin("arc-");
    assert!(load(&p, bytes).is_ok());
    let text = std::str::from_utf8(bytes).unwrap();
    let variants = [
        text.replacen("\"name\"", " \"name\"", 1),
        text.trim_end().to_string(),
        format!("{text}\n"),
        text.replacen("{\"format\":\"rid-cover/1\",\"name\":\"arc-\"", "{\"name\":\"arc-\",\"format\":\"rid-cover/1\"", 1),
    ];
    for v in &variants {
        assert!(matches!(load(&p, v.as_bytes()), Err(ProofError::Format(FormatError::NonCanonical))), "{}", &v[..60]);
    }
    let json = [
        text.replacen("rid-cover/1", "rid-cover/2", 1),
        text.replacen("\"scope\":\"domain\"", "\"scope\":\"domain\",\"extra\":1", 1),
        text.replacen("\"beyond\":null", "\"beyond\":null,\"note\":\"x\"", 1),
        text.replacen(",\"beyond\":null", "", 1),
        text.replacen("\"radius\":\"1/32\"", "\"radius\":\"2/64\"", 1),
        text.replacen("\"radius\":\"1/32\"", "\"radius\":0.03125", 1),
        text.replacen("\"tree\":\"", "\"tree\":1,\"x\":\"", 1),
        text[..text.len() / 2].to_string(),
        String::new(),
    ];
    for v in &json {
        assert!(matches!(load(&p, v.as_bytes()), Err(ProofError::Format(FormatError::Json(_)))), "{}", &v[..v.len().min(60)]);
    }
}

#[test]
fn files_that_differ_from_the_pin_are_refused() {
    let p = pin("crossing+");
    let bytes = cover_data("crossing+");
    let mut other = p.clone();
    other.name = "crossing-".into();
    assert!(matches!(load(&other, bytes), Err(ProofError::Pin("name"))));
    let mut other = p.clone();
    other.scope = Scope::All;
    assert!(matches!(load(&other, bytes), Err(ProofError::Pin("scope"))));
    // The file of the other sign under this pin.
    assert!(matches!(load(&pin("crossing-"), bytes), Err(ProofError::Pin("name"))));
    let mut renamed = CoverFile::parse(cover_data("crossing-")).unwrap();
    renamed.name = "crossing+".into();
    assert!(matches!(load(&p, &renamed.to_bytes()), Err(ProofError::Pin("parameters"))));
    // Any parameter changed in the file alone.
    let changes: Vec<Box<dyn Fn(&mut CoverFile)>> = vec![
        Box::new(|f| f.parameters.centre[0] = Number(QSqrt5::from_rational(frac(1, 1000)))),
        Box::new(|f| f.parameters.map[1][1] = Number(QSqrt5::integer(2))),
        Box::new(|f| if let Shape::Point { radius, .. } = &mut f.parameters.shape { radius.0 = frac(1, 8) }),
        Box::new(|f| if let Shape::Point { ratio, .. } = &mut f.parameters.shape { ratio.0 = frac(3, 2) }),
        Box::new(|f| if let Shape::Point { window: Some(w), .. } = &mut f.parameters.shape { w.radius.0 = frac(1, 2) }),
        Box::new(|f| if let Shape::Point { radii, .. } = &mut f.parameters.shape { radii[1].0 = frac(1, 8) }),
        Box::new(|f| f.parameters.coordinates = crate::elimination::zoom::Coordinates::ArcPlane(crate::elimination::zoom::Sign::Minus)),
    ];
    for change in changes {
        assert!(matches!(load(&p, &edited("crossing+", change)), Err(ProofError::Pin("parameters"))));
    }
    let mut scope_all = pin("pentagon");
    scope_all.scope = Scope::All;
    let bytes = edited("pentagon", |f| f.scope = Scope::All);
    assert!(matches!(load(&scope_all, &bytes), Err(ProofError::Pin("beyond"))));
}

/// Loads a file whose parameters were changed together with the pin: the
/// cells no longer verify, or the zooms are refused.
fn load_moved(name: &str, change: impl Fn(&mut Parameters)) -> Result<ProvedCover, ProofError> {
    let mut p = pin(name);
    change(&mut p.parameters);
    let bytes = edited(name, |f| change(&mut f.parameters));
    load(&p, &bytes)
}

#[test]
fn moved_scaled_or_enlarged_zooms_are_refused() {
    // Square: the centre moved along s, and off the aligned set.
    assert!(matches!(load_moved("square", |p| p.centre[0] = Number(QSqrt5::from_rational(frac(1, 64)))), Err(ProofError::Cell { .. })));
    assert!(matches!(load_moved("square", |p| p.centre[3] = Number(QSqrt5::from_rational(frac(1, 64)))), Err(ProofError::Cover(_))));
    // Square and a Local tube: enlarged radius.
    assert!(matches!(load_moved("square", |p| if let Shape::Point { radii, radius, .. } = &mut p.shape {
        radii[0].0 = frac(1, 4);
        radius.0 = frac(1, 4);
    }), Err(ProofError::Cell { .. })));
    assert!(matches!(load_moved("7", |p| if let Shape::Tube { radius, .. } = &mut p.shape { radius.0 = &radius.0 * q(8) }), Err(ProofError::Cell { .. })));
    // A Local tube moved to another view rectangle.
    assert!(matches!(load_moved("7", |p| if let Shape::Tube { base, .. } = &mut p.shape {
        base[0] = [Rational(&base[0][0].0 + frac(1, 4)), Rational(&base[0][1].0 + frac(1, 4))];
    }), Err(ProofError::Cell { .. })));
    // The first arc's tube: offset scaling doubled, base moved to the crossing.
    assert!(matches!(load_moved("arc+", |p| p.map[1][1] = Number(QSqrt5::integer(4))), Err(ProofError::Cell { .. })));
    assert!(matches!(load_moved("arc+", |p| if let Shape::Tube { base, .. } = &mut p.shape {
        base[0] = [Rational(q(0)), Rational(frac(9, 64))];
    }), Err(ProofError::Cell { .. })));
    // The crossing: an enlarged window, and a changed shear.
    assert!(matches!(load_moved("crossing+", |p| if let Shape::Point { window: Some(w), .. } = &mut p.shape { w.radius.0 = frac(1, 2) }), Err(ProofError::Cell { .. })));
    assert!(matches!(load_moved("crossing+", |p| if let Shape::Point { window: Some(w), .. } = &mut p.shape { w.shear[1][0].0 = frac(1, 2) }), Err(ProofError::Cell { .. } | ProofError::Delegation(_) | ProofError::Cover(_))));
    // The endpoint moved to a rational centre.
    assert!(matches!(load_moved("endpoint-", |p| p.centre[0] = Number(QSqrt5::from_rational(frac(1, 4)))), Err(ProofError::Cell { .. })));
}

#[test]
fn corrupted_trees_and_leaves_are_refused() {
    let p = pin("arc-");
    fn zoom(f: &mut CoverFile) -> &mut crate::elimination::proof::format::ZoomEntry {
        &mut f.zooms[3]
    }
    // Dropped and appended cells.
    let bytes = edited("arc-", |f| { zoom(f).tree.pop(); });
    assert!(matches!(load(&p, &bytes), Err(ProofError::Tree { .. })));
    let bytes = edited("arc-", |f| zoom(f).tree.push('.'));
    assert!(matches!(load(&p, &bytes), Err(ProofError::Tree { .. })));
    // A leaf split in two without a leaf entry for the new cell.
    let bytes = edited("arc-", |f| {
        let t = &mut zoom(f).tree;
        let at = t.find('.').unwrap();
        t.replace_range(at..at + 1, "0..");
    });
    assert!(matches!(load(&p, &bytes), Err(ProofError::LeafCount { .. })));
    let bytes = edited("arc-", |f| { zoom(f).leaves.pop(); });
    assert!(matches!(load(&p, &bytes), Err(ProofError::LeafCount { .. })));
    // Splits along other variables: the trees still tile, but the cells
    // change; a changed tree is only accepted if its cells verify, which
    // fails for some of these.
    let mut refused = 0;
    for axis in ["0", "2", "3", "4"] {
        let bytes = edited("arc-", |f| {
            let t = &mut zoom(f).tree;
            let at = t.find(|c: char| c.is_ascii_digit()).unwrap();
            t.replace_range(at..at + 1, axis);
        });
        match load(&p, &bytes) {
            Ok(_) => {}
            Err(ProofError::Cell { .. }) => refused += 1,
            Err(other) => panic!("{other}"),
        }
    }
    assert!(refused >= 2, "{refused}");
    // Zooms renamed, reordered or dropped.
    let bytes = edited("arc-", |f| f.zooms.swap(0, 1));
    assert!(matches!(load(&p, &bytes), Err(ProofError::Zooms)));
    let bytes = edited("arc-", |f| { f.zooms.pop(); });
    assert!(matches!(load(&p, &bytes), Err(ProofError::Zooms)));
    let bytes = edited("arc-", |f| f.zooms[2].name = "T o1".into());
    assert!(matches!(load(&p, &bytes), Err(ProofError::Zooms)));
    let bytes = edited("arc-", |f| f.zooms.clear());
    assert!(matches!(load(&p, &bytes), Err(ProofError::Zooms)));
    // A witness index out of range.
    let bytes = edited("arc-", |f| {
        let n = f.witnesses.len();
        f.zooms[0].leaves[0] = Leaf::Witness { index: n, factor: [0, 1, 0, 0, 0] };
    });
    assert!(matches!(load(&p, &bytes), Err(ProofError::WitnessIndex(_))));
    // Wrong factors: outside the scale variables, too high, too low.
    for factor in [[1, 1, 0, 0, 0], [0, 3, 0, 0, 0], [0, 0, 0, 0, 0]] {
        let bytes = edited("arc-", |f| {
            for leaf in &mut f.zooms[1].leaves {
                if let Leaf::Witness { factor: d, .. } = leaf {
                    *d = factor;
                }
            }
        });
        let outside = factor[0] > 0;
        match load(&p, &bytes) {
            Err(ProofError::Factor { error: ZoomError::NotScale { variable: 0 }, .. }) => assert!(outside),
            Err(ProofError::Cell { .. }) => assert!(!outside),
            other => panic!("{factor:?}: {:?}", other.map(|_| ())),
        }
    }
    // Every leaf of a zoom reassigned to one witness.
    let bytes = edited("arc-", |f| {
        for leaf in &mut f.zooms[2].leaves {
            if let Leaf::Witness { index, .. } = leaf {
                *index = 0;
            }
        }
    });
    assert!(matches!(load(&p, &bytes), Err(ProofError::Cell { .. })));
}

#[test]
fn domain_inequalities_outside_the_scope_are_refused() {
    // A Local cover's statement holds without D: a domain inequality there
    // is refused even though it would prove the cell on D.
    let p = pin("3");
    let bytes = edited("3", |f| {
        f.witnesses.push(Witness::Domain(DomainInequality::new(13).unwrap()));
        let index = f.witnesses.len() - 1;
        f.zooms[4].leaves[0] = Leaf::Witness { index, factor: [0, 0, 0, 0, 0] };
    });
    assert!(matches!(load(&p, &bytes), Err(ProofError::Scope(_))));
    // The square uses D's fold: claiming the statement without D is refused.
    let mut all = pin("square");
    all.scope = Scope::All;
    let bytes = edited("square", |f| f.scope = Scope::All);
    assert!(matches!(load(&all, &bytes), Err(ProofError::Scope(_))));
}

#[test]
fn delegations_outside_the_window_are_refused() {
    let p = pin("crossing+");
    let file = CoverFile::parse(cover_data("crossing+")).unwrap();
    let names: Vec<&str> = file.zooms.iter().map(|c| c.name.as_str()).collect();
    let find = |name: &str| names.iter().position(|n| *n == name).unwrap();
    // A B zoom, an A zoom on the base face e < 0 and a sheared zoom
    // never delegate.
    for zoom in [find("B o1-"), find("A b0- o1-"), find("A' b0+ o0+")] {
        let bytes = edited("crossing+", |f| f.zooms[zoom].leaves[0] = Leaf::Delegated);
        assert!(matches!(load(&p, &bytes), Err(ProofError::Delegation(_))), "{}", names[zoom]);
    }
    // An A zoom on the window face, at a leaf far from the second arc.
    let zoom = find("A b0+ o1+");
    let bytes = edited("crossing+", |f| {
        let leaves = &mut f.zooms[zoom].leaves;
        let at = leaves.iter().position(|l| *l != Leaf::Delegated).unwrap();
        leaves[at] = Leaf::Delegated;
    });
    assert!(matches!(load(&p, &bytes), Err(ProofError::Delegation(_))));
    // The delegated cell at the second arc's factor point (μ = 0, ρ = 1,
    // f̂ = −e_θ) proved by a witness instead: no witness is negative there.
    let zoom = find("A b0- o1-") + 7;
    assert_eq!(names[zoom], "A b0+ o1-");
    let cover = crate::elimination::zoom::cover::ZoomCover::new(file.parameters.clone()).unwrap();
    let cells = tree::leaves(cover.zooms()[zoom].root(), &file.zooms[zoom].tree).unwrap();
    let point = exact(&[q(0), q(1), q(0), q(0), q(0)]);
    let at = cells
        .iter()
        .zip(&file.zooms[zoom].leaves)
        .position(|(cell, leaf)| *leaf == Leaf::Delegated && contains_point(cell, &point))
        .unwrap();
    for index in 0..file.witnesses.len() {
        for factor in [[1, 0, 0, 0, 0], [1, 1, 0, 0, 0]] {
            let factor = ZoomFactor::new(factor);
            let refusal = crate::elimination::lemma::check_cell(&cover.zooms()[zoom], &cells[at], &file.witnesses[index], &factor);
            assert!(refusal.is_err(), "witness {index}");
        }
    }
    for index in [0, file.witnesses.len() - 1] {
        let bytes = edited("crossing+", |f| f.zooms[zoom].leaves[at] = Leaf::Witness { index, factor: [1, 0, 0, 0, 0] });
        assert!(matches!(load(&p, &bytes), Err(ProofError::Cell { .. })), "witness {index}");
    }
}

#[test]
fn verification_is_deterministic_and_interruptible() {
    let p = pin("endpoint+");
    let bytes = cover_data("endpoint+");
    let one = ProvedCover::load(&p, bytes, threads(1), &running()).unwrap();
    let four = ProvedCover::load(&p, bytes, threads(4), &running()).unwrap();
    assert_eq!(one.report(), four.report());
    // Several corrupted leaves: the reported refusal is the same for every
    // thread count.
    let corrupted = edited("endpoint+", |f| {
        for zoom in f.zooms.iter_mut().skip(3) {
            for leaf in zoom.leaves.iter_mut().step_by(5) {
                if let Leaf::Witness { index, .. } = leaf {
                    *index = (*index + 1) % 7;
                }
            }
        }
    });
    let messages: Vec<String> = [1, 2, 3, 4]
        .iter()
        .map(|&n| ProvedCover::load(&p, &corrupted, threads(n), &running()).unwrap_err().to_string())
        .collect();
    assert!(messages.windows(2).all(|w| w[0] == w[1]), "{messages:?}");
    let stop = AtomicBool::new(true);
    assert!(matches!(ProvedCover::load(&p, bytes, threads(2), &stop), Err(ProofError::Interrupted)));
}

// ---- An independent second verification -------------------------------------

fn binomial(n: u32, k: u32) -> Q {
    (0..k).fold(Q::one(), |c, i| c * frac(i64::from(n - i), i64::from(i + 1)))
}

/// The textbook tensor Bernstein coefficients of `p` on `cell`: with
/// `p(l + w∘x) = Σ_J c_J x^J`, `b_I = Σ_{J ≤ I} Π_j C(I_j, J_j)/C(n_j, J_j) c_J`.
fn textbook_bernstein(p: &Polynomial, cell: &Five<Interval>) -> Vec<QSqrt5> {
    let x = Polynomial::variables();
    let affine: Five<Polynomial> = std::array::from_fn(|j| {
        &Polynomial::constant(QSqrt5::from_rational(cell[j].lo().clone())) + &x[j].scale(&QSqrt5::from_rational(cell[j].width()))
    });
    let shifted = p.substitute(&affine).unwrap();
    let n = p.degrees();
    let mut out = Vec::new();
    let count: usize = n.iter().map(|&d| usize::from(d) + 1).product();
    for flat in 0..count {
        // Mixed radix, variable 4 fastest.
        let mut rest = flat;
        let mut index = [0u8; 5];
        for j in (0..5).rev() {
            let base = usize::from(n[j]) + 1;
            index[j] = (rest % base) as u8;
            rest /= base;
        }
        let mut b = QSqrt5::zero();
        for (e, c) in shifted.terms() {
            if (0..5).all(|j| e[j] <= index[j]) {
                let weight = (0..5).fold(Q::one(), |w, j| {
                    w * binomial(u32::from(index[j]), u32::from(e[j])) / binomial(u32::from(n[j]), u32::from(e[j]))
                });
                b = &b + &c.scale(&weight);
            }
        }
        out.push(b);
    }
    out
}

/// The witness at zoom variables `y`, evaluated directly at the zoom's
/// configuration (see `lemma::tests::independent_value`).
fn witness_value(witness: &Witness, zoom: &Zoom, y: &Five<QSqrt5>) -> QSqrt5 {
    let x = zoom.evaluate(y);
    let constants: Five<Polynomial> = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (view, rotation) = zoom.coordinates().configuration(&constants);
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    let u: Point = view.vector().map(|p| p.evaluate(&origin));
    let r: Point = rotation.map(|p| p.evaluate(&origin));
    crate::elimination::lemma::tests::independent_value(witness, &u, &r)
}

/// Checks `g ∘ zoom = y^a N` exactly on a tensor grid of `degree + 1`
/// points per variable: both sides have degree at most `degree` in each
/// variable (the zoom coordinates are multilinear and `g` has total degree
/// at most 4), so agreement on the grid is the polynomial identity.
fn identity_on_grid(witness: &Witness, zoom: &Zoom, factor: &Exponents, n: &Polynomial) {
    assert!(zoom.map().iter().all(|p| p.terms().all(|(e, _)| e.iter().all(|&k| k <= 1))));
    let degree = 4usize;
    for (j, d) in n.degrees().iter().enumerate() {
        assert!(usize::from(*d) + usize::from(factor[j]) <= degree);
    }
    // Grid points k/8, k = 0..4: small enough that u₃ ≠ 0, which the affine
    // evaluation needs.
    let points = degree + 1;
    for flat in 0..points.pow(5) {
        let y: Five<QSqrt5> = std::array::from_fn(|j| QSqrt5::from_rational(frac((flat / points.pow(j as u32) % points) as i64, 8)));
        let monomial = (0..5).fold(QSqrt5::one(), |p, j| (0..factor[j]).fold(p, |p, _| &p * &y[j]));
        assert!(!zoom.view().vector()[2].evaluate(&y).is_zero());
        assert_eq!(witness_value(witness, zoom, &y), &monomial * &n.evaluate(&y));
    }
}

/// The second verification of every `step`-th witness leaf of a cover.
fn independent_verification(name: &str, step: usize, grid: bool) -> usize {
    let file = CoverFile::parse(cover_data(name)).unwrap();
    let mut checked = 0;
    let mut identities: HashMap<(usize, usize, Exponents), Polynomial> = HashMap::new();
    {
        let cover = crate::elimination::zoom::cover::ZoomCover::new(file.parameters.clone()).unwrap();
        for (c, data) in file.zooms.iter().enumerate() {
            let zoom = &cover.zooms()[c];
            let cells = tree::leaves(zoom.root(), &data.tree).unwrap();
            for (leaf, (cell, e)) in cells.iter().zip(&data.leaves).enumerate() {
                let Leaf::Witness { index, factor } = e else { continue };
                if leaf % step != 0 {
                    continue;
                }
                let witness = &file.witnesses[*index];
                let n = identities.entry((c, *index, *factor)).or_insert_with(|| {
                    let n = crate::elimination::lemma::Quotient::new(zoom, witness, &ZoomFactor::new(*factor)).unwrap().polynomial().clone();
                    if grid {
                        identity_on_grid(witness, zoom, factor, &n);
                    }
                    n
                });
                // Divided variables nonnegative, u₃ > 0 at the corners (u₃ is
                // multilinear), validity by brute force at every corner view.
                assert!((0..5).all(|j| factor[j] == 0 || (zoom.is_scale(j) && cell[j].lo() >= &Q::zero())));
                for mask in 0..32u32 {
                    let y = exact(&std::array::from_fn(|j| if mask >> j & 1 == 1 { cell[j].hi() } else { cell[j].lo() }.clone()));
                    let u = zoom.view().vector().map(|p| p.evaluate(&y));
                    assert_eq!(u[2].sign(), std::cmp::Ordering::Greater);
                    if let Witness::Gap(gap) = witness {
                        let normal = gap.support().normal_at(&u);
                        let h = geometry::dot(&normal, &geometry::vertices()[gap.support().contact()]);
                        assert!(geometry::vertices().iter().all(|w| geometry::dot(&normal, w) <= h), "{name} {} leaf {leaf}", zoom.name());
                    }
                }
                let coefficients = textbook_bernstein(n, cell);
                assert!(coefficients.iter().all(|b| b.sign() == std::cmp::Ordering::Less), "{name} {} leaf {leaf}", zoom.name());
                checked += 1;
            }
        }
    }
    checked
}

#[test]
fn a_sample_of_cells_verifies_independently() {
    let mut checked = 0;
    for (name, step) in [("square", 15), ("pentagon", 200), ("arc+", 12), ("arc-", 6), ("endpoint+", 12), ("endpoint-", 12), ("crossing+", 50), ("crossing-", 25), ("0", 4), ("28", 4)] {
        checked += independent_verification(name, step, false);
    }
    assert!(checked > 200, "{checked}");
    // The exact identity of pull-back and division on full grids, for a
    // few leaves of each kind of zoom.
    for name in ["arc-", "11"] {
        let file = CoverFile::parse(cover_data(name)).unwrap();
        let cover = crate::elimination::zoom::cover::ZoomCover::new(file.parameters.clone()).unwrap();
        for (c, data) in file.zooms.iter().enumerate().step_by(3) {
            if let Some(Leaf::Witness { index, factor }) = data.leaves.first() {
                let zoom = &cover.zooms()[c];
                let witness = &file.witnesses[*index];
                let n = crate::elimination::lemma::Quotient::new(zoom, witness, &ZoomFactor::new(*factor)).unwrap();
                identity_on_grid(witness, zoom, factor, n.polynomial());
            }
        }
    }
}

/// The full campaign: every witness leaf of every cover, with the grid
/// identity for every zoom, witness and factor. About an hour on one
/// core; run on several with
/// `cargo test --release --lib proof::tests::independent_verification_of_every_cell -- --ignored --test-threads=4`
/// (the covers are split over four tests).
fn campaign(part: usize) {
    let names: Vec<String> = EXOTIC.iter().map(|s| s.to_string()).chain((0..LOCAL_COUNT).map(|n| n.to_string())).collect();
    for (i, name) in names.iter().enumerate() {
        // The pentagon, the largest cover, is split across all parts.
        if name == "pentagon" || i % 4 != part {
            continue;
        }
        let n = independent_verification(name, 1, true);
        eprintln!("{name}: {n} cells verified independently");
    }
}

#[test]
#[ignore]
fn independent_verification_of_every_cell_part_0() {
    campaign(0);
}
#[test]
#[ignore]
fn independent_verification_of_every_cell_part_1() {
    campaign(1);
}
#[test]
#[ignore]
fn independent_verification_of_every_cell_part_2() {
    campaign(2);
}
#[test]
#[ignore]
fn independent_verification_of_every_cell_part_3() {
    campaign(3);
    let n = independent_verification("pentagon", 1, true);
    eprintln!("pentagon: {n} cells verified independently");
}

/// Adversarial tests and a full re-verification.
mod adversarial;
