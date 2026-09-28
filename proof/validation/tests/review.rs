//! Independent adversarial review tests of the validation package, written
//! against its public API only, with a toy collection of components.
use rid::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use rid::problem::configuration::{ConfigurationBox, AXES, R};
use rid::components::domain::Domain;
use rid::components::Component;
use rid::problem::geometry::{rotate, vertices};
use rid::search::certificate::RecordData;
use rid::search::BoxError;
use rid_validation::config;
use rid_validation::experiments::completeness::{self, compare::compare, transcript};
use rid_validation::experiments::probe::{canonical, Probes, RecordForm};
use rid_validation::experiments::run;
use rid_validation::points::catalogue::{Catalogue, Centre};
use rid_validation::points::neighbourhood::{check, neighbourhood, radius};
use rid_validation::points::relation::{relation, Relation};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

// ---- A toy collection with configurable names and faults --------------

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
struct Rec(Value);
impl RecordData for Rec {}

/// The checkers read the toy's records as the crate's record data, as the
/// binary does: the toy writes canonical crate records.
const FORM: RecordForm = canonical::<rid::components::record::RecordData>;

/// The toy's Global record for axis `w`: a well-formed Global witness whose
/// plug vertex encodes the axis.
fn global_witness(w: usize) -> Value {
    json!({ "edge": [16, 17], "vertex": w })
}

struct Toy {
    names: [&'static str; 4],
    square: ConfigurationBox,
    /// Global: some rotation coordinate keeps this distance from zero.
    global: Q,
    /// Global also accepts every box whose rotation part contains zero.
    unsound_global: bool,
    /// Global's attempt refuses boxes whose `s` width is `2^(1-k)` for these
    /// `k` (a missed witness: its verification still accepts them).
    refuse: Vec<u32>,
    /// Global's verification refuses those boxes too (not inherited by
    /// sub-boxes: a defect the runner must catch).
    nonmonotone: bool,
}

impl Toy {
    fn new() -> Self {
        let r = || Interval::new(frac(-1, 64), frac(1, 64)).unwrap();
        Self {
            names: ["Domain", "Exotic", "Local", "Global"],
            square: ConfigurationBox::new([
                Interval::new(q(0), frac(1, 16)).unwrap(),
                Interval::new(q(0), frac(1, 16)).unwrap(),
                r(),
                r(),
                r(),
            ]),
            global: frac(1, 8),
            unsound_global: false,
            refuse: Vec::new(),
            nonmonotone: false,
        }
    }

    fn clone_settings(&self) -> Toy {
        Toy {
            names: self.names,
            square: self.square.clone(),
            global: self.global.clone(),
            unsound_global: self.unsound_global,
            refuse: self.refuse.clone(),
            nonmonotone: self.nonmonotone,
        }
    }

    fn record(&self, index: usize, fields: Value) -> Rec {
        let mut object = fields.as_object().unwrap().clone();
        object.insert("component".into(), json!(self.names[index]));
        Rec(Value::Object(object))
    }

    fn global_axis(&self, b: &ConfigurationBox) -> Option<usize> {
        let width = b.axes()[0].width();
        if self.refuse.iter().any(|&k| width == radius(k - 1)) {
            return None;
        }
        let far = |axis: usize| {
            let a = &b.axes()[axis];
            *a.lo() >= self.global || *a.hi() <= -&self.global
        };
        if let Some(axis) = R.into_iter().find(|&a| far(a)) {
            return Some(axis);
        }
        let zero = Q::zero();
        (self.unsound_global && R.into_iter().all(|a| b.axes()[a].contains(&zero))).then_some(R[0])
    }
}

impl Probes for Toy {
    type Data = Rec;
    fn components(&self) -> Vec<String> {
        self.names.iter().map(|s| s.to_string()).collect()
    }
    fn decide(&self, b: &ConfigurationBox) -> Result<Option<Rec>, BoxError> {
        Ok((0..4).find_map(|i| self.attempt(i, b)))
    }
    fn attempt(&self, index: usize, b: &ConfigurationBox) -> Option<Rec> {
        match index {
            0 => Domain.check(b).map(|inequality| self.record(0, json!({ "inequality": inequality.index() }))),
            1 => self.square.contains(b).then(|| self.record(1, json!({"cover": "square"}))),
            2 => None,
            _ => self.global_axis(b).map(|w| self.record(3, global_witness(w))),
        }
    }
    fn covers(&self, index: usize) -> Vec<Rec> {
        if index == 1 {
            vec![self.record(1, json!({"cover": "square"}))]
        } else {
            Vec::new()
        }
    }
    fn verify(&self, index: usize, b: &ConfigurationBox, data: &Rec) -> Result<(), BoxError> {
        let expected = if index == 3 && !self.nonmonotone {
            let honest = Toy { refuse: Vec::new(), ..self.clone_settings() };
            honest.attempt(index, b)
        } else {
            self.attempt(index, b)
        };
        if expected.as_ref() == Some(data) {
            Ok(())
        } else {
            Err("refused".into())
        }
    }
}

// ---- Helpers ------------------------------------------------------------

struct Dir(PathBuf);
impl Dir {
    fn new(name: &str) -> Self {
        let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("review-{name}"));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    fn join(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

fn two() -> NonZeroUsize {
    NonZeroUsize::new(2).unwrap()
}

fn rational_centre(values: [&str; 5]) -> String {
    let pairs: Vec<String> = values.iter().map(|v| format!(r#"["{v}","0"]"#)).collect();
    format!("[{}]", pairs.join(","))
}

fn catalogue(dir: &Dir) -> PathBuf {
    let points = [
        ("square", rational_centre(["0", "0", "0", "0", "0"])),
        ("poke", rational_centre(["1/4", "1/8", "1/4", "0", "0"])),
        ("aligned", rational_centre(["1/5", "1/10", "0", "0", "0"])),
        ("outside", rational_centre(["2/3", "2/5", "0", "0", "0"])),
        ("inner-poke", rational_centre(["1/5", "1/10", "1/32", "0", "0"])),
    ];
    let lines: Vec<String> = points
        .iter()
        .map(|(id, c)| format!(r#"{{"id":"{id}","group":"g","centre":{c},"expected":["Global"]}}"#))
        .collect();
    let path = dir.join("points.json");
    std::fs::write(
        &path,
        format!(r#"{{"format":"rid-validation-points/1","points":[{}]}}"#, lines.join(",")),
    )
    .unwrap();
    path
}

fn completeness_config(dir: &Dir, name: &str, max_k: u32, selection: Option<Vec<&str>>) -> config::Completeness {
    config::Completeness {
        catalogue: catalogue(dir),
        transcript: dir.join(name),
        threads: two(),
        max_k,
        selection: selection.map(|s| s.into_iter().map(String::from).collect()),
    }
}

fn checked(config: &config::Completeness) -> transcript::Checked {
    let catalogue = Catalogue::load(&config.catalogue).unwrap();
    transcript::check(&catalogue, &run::read(&config.transcript).unwrap(), two(), FORM).unwrap()
}

// ---- Compare: regressions -------------

/// A candidate transcript that omits points the baseline has counts each
/// missing point as a regression, so `fail_on_regression` fails such a run.
#[test]
fn compare_counts_points_missing_from_the_candidate_as_regressions() {
    let dir = Dir::new("missing");
    let base = completeness_config(&dir, "base.jsonl", 40, None);
    completeness::run(&Toy::new(), &base).unwrap();
    let part = completeness_config(&dir, "part.jsonl", 40, Some(vec!["aligned"]));
    completeness::run(&Toy::new(), &part).unwrap();
    let c = compare(&checked(&base), &checked(&part), 1).unwrap();
    assert_eq!(c.only_baseline, ["square", "poke", "outside", "inner-poke"]);
    assert_eq!(c.regressions, 4);
}

/// A finer refusal the baseline did not have is flagged and counted, even
/// when the baseline already had another one; both sets are listed.
#[test]
fn compare_counts_additional_finer_refusals() {
    let dir = Dir::new("finer");
    let base = completeness_config(&dir, "base.jsonl", 40, Some(vec!["inner-poke"]));
    completeness::run(&Toy { refuse: vec![14], global: frac(1, 64), ..Toy::new() }, &base).unwrap();
    let worse = completeness_config(&dir, "worse.jsonl", 40, Some(vec!["inner-poke"]));
    completeness::run(&Toy { refuse: vec![8, 14], global: frac(1, 64), ..Toy::new() }, &worse).unwrap();
    let (b, w) = (checked(&base), checked(&worse));
    assert!(b.points[0].in_domain && b.points[0].relation == Relation::Poke);
    assert_eq!(b.points[0].decision.best_k, Some(6));
    assert_eq!(b.points[0].decision.finer_refusals, [14]);
    assert_eq!(w.points[0].decision.finer_refusals, [8, 14]);
    let c = compare(&b, &w, 1).unwrap();
    assert_eq!(c.regressions, 1);
    assert_eq!(c.changes[0].finer_refusals, [vec![14], vec![8, 14]]);
    assert!(c.changes[0].flags.contains(&completeness::compare::Flag::NewFinerRefusals));
    // Fewer finer refusals: flagged as cleared, not a regression.
    let back = compare(&w, &b, 1).unwrap();
    assert_eq!(back.regressions, 0);
    assert!(back.changes[0].flags.contains(&completeness::compare::Flag::FinerRefusalsCleared));
}

/// A verification that a finer, nested neighbourhood refuses although it
/// accepted the coarser one stops the run: verification must be inherited.
#[test]
fn a_verification_not_inherited_by_sub_boxes_stops_the_run() {
    let dir = Dir::new("monotone");
    let config = completeness_config(&dir, "t.jsonl", 40, Some(vec!["inner-poke"]));
    let toy = Toy { refuse: vec![8], global: frac(1, 64), nonmonotone: true, ..Toy::new() };
    assert!(matches!(
        completeness::run(&toy, &config),
        Err(completeness::Error::NotMonotone { k: 8, .. })
    ));
    // The same refusal by the attempt alone (a missed witness) runs.
    let config = completeness_config(&dir, "u.jsonl", 40, Some(vec!["inner-poke"]));
    let toy = Toy { refuse: vec![8], global: frac(1, 64), ..Toy::new() };
    completeness::run(&toy, &config).unwrap();
}

/// Transcripts with different finest exponents are not compared.
#[test]
fn compare_refuses_transcripts_with_different_finest_exponents() {
    let dir = Dir::new("maxk");
    let base = completeness_config(&dir, "base.jsonl", 2, Some(vec!["poke"]));
    completeness::run(&Toy::new(), &base).unwrap();
    let finer = completeness_config(&dir, "finer.jsonl", 40, Some(vec!["poke"]));
    completeness::run(&Toy::new(), &finer).unwrap();
    assert_eq!(
        compare(&checked(&base), &checked(&finer), 1).unwrap_err(),
        completeness::compare::Incompatible::MaxK([2, 40])
    );
}

// ---- The checker reads records as typed records ----------------------------

/// A record replaced by one without its data, or with a cover of the wrong
/// type, is refused by the checker, which parses every record as the
/// components' record type. (A well-formed record that its component would
/// refuse still passes: the checker does not replay eliminations, which it
/// reports as `"replayed": false`.)
#[test]
fn the_checker_accepts_records_their_component_would_refuse() {
    let dir = Dir::new("records");
    let config = completeness_config(&dir, "t.jsonl", 40, Some(vec!["inner-poke"]));
    let toy = Toy { global: frac(1, 64), ..Toy::new() };
    completeness::run(&toy, &config).unwrap();
    let catalogue = Catalogue::load(&config.catalogue).unwrap();
    let bytes = run::read(&config.transcript).unwrap();
    let lines = run::lines(&bytes).unwrap();
    for (replacement, typed) in [
        (json!({"component": "Global", "edge": [16, 17], "vertex": 4}), true),
        (json!({"component": "Global"}), false),
        (json!({"component": "Global", "edge": [16, 58], "vertex": 4}), false),
        (json!({"component": "Global", "witness": global_witness(4)}), false),
        (json!({"component": "Global", "witness": 4}), false),
    ] {
        let mut point: transcript::PointResult = run::parse(lines[1], 2).unwrap();
        let replace = |s: &mut transcript::Schedule| {
            for t in &mut s.trials {
                if t.record.is_some() {
                    t.record = Some(replacement.clone());
                }
            }
            if s.record.is_some() {
                s.record = Some(replacement.clone());
            }
        };
        replace(&mut point.decision);
        replace(&mut point.components[3]);
        let mut forged = lines[0].to_vec();
        forged.push(b'\n');
        forged.extend(serde_json::to_vec(&point).unwrap());
        forged.push(b'\n');
        forged.extend_from_slice(lines[2]);
        forged.push(b'\n');
        // Toy's Global verifies only witness 2 (axis r1) for this point.
        let b = neighbourhood(&catalogue.points[4].centre, 6).unwrap();
        assert!(toy.verify(3, &b, &Rec(replacement.clone())).is_err());
        let outcome = transcript::check(&catalogue, &forged, two(), FORM);
        assert_eq!(outcome.is_ok(), typed, "{replacement}: {:?}", outcome.err());
    }
}

// ---- The contradiction rules need the exact component names ------------

/// A collection whose Global is named differently would escape the Global
/// rule; such names are refused before any work.
#[test]
fn renamed_components_are_refused() {
    let dir = Dir::new("names");
    let unsound = Toy { unsound_global: true, ..Toy::new() };
    let named = completeness_config(&dir, "named.jsonl", 8, Some(vec!["aligned"]));
    assert!(matches!(
        completeness::run(&unsound, &named),
        Err(completeness::Error::Contradiction { .. })
    ));
    let renamed = Toy {
        names: ["Domain", "Exotic", "Local", "global"],
        ..unsound
    };
    let config = completeness_config(&dir, "renamed.jsonl", 8, Some(vec!["aligned"]));
    assert!(completeness::run(&renamed, &config).is_err());
    assert!(!config.transcript.exists());
}

// ---- Neighbourhoods of irrational centres near the root faces ---------

/// `(161 - 72√5)·2^-shift`: a tiny positive irrational number.
fn tiny(shift: u32) -> QSqrt5 {
    QSqrt5::new(q(161), q(-72)).scale(&radius(shift))
}

#[test]
fn neighbourhoods_near_faces_contain_their_centre_and_check() {
    let root = ConfigurationBox::root();
    for shift in [0, 20, 60, 300] {
        let d = tiny(shift);
        for face in 0..2 {
            let coordinates: [QSqrt5; AXES] = std::array::from_fn(|j| {
                let a = &root.axes()[j];
                let end = QSqrt5::from_rational(if face == 0 { a.lo().clone() } else { a.hi().clone() });
                if face == 0 { &end + &d } else { &end - &d }
            });
            let c = Centre::new(coordinates.clone()).unwrap();
            for k in [0, 1, 9, 10, 11, 58, 59, 60, 61, 62, 299, 300, 301, 310, 2048] {
                let b = neighbourhood(&c, k).unwrap();
                assert_eq!(check(&c, k, &b), Ok(()), "shift {shift} face {face} k {k}");
                for (x, a) in coordinates.iter().zip(b.axes()) {
                    let (lo, hi) = (QSqrt5::from_rational(a.lo().clone()), QSqrt5::from_rational(a.hi().clone()));
                    assert!(lo <= *x && *x <= hi && lo < hi);
                }
                assert!(root.contains(&b));
            }
        }
    }
}

// ---- The relation against an independent floating-point hull ----------

fn project(u: [f64; 2], p: [f64; 3]) -> [f64; 2] {
    [p[0] - u[0] * p[2], p[1] - u[1] * p[2]]
}

fn f(x: &QSqrt5) -> f64 {
    let (a, b) = (x.rational_part(), x.sqrt5_part());
    let v = |q: &Q| q.numer().to_string().parse::<f64>().unwrap() / q.denom().to_string().parse::<f64>().unwrap();
    v(a) + v(b) * 5f64.sqrt()
}

/// Andrew's monotone chain, counter-clockwise, collinear points dropped.
fn hull(mut points: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    let mut lower: Vec<[f64; 2]> = Vec::new();
    for &p in &points {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 1e-12 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<[f64; 2]> = Vec::new();
    for &p in points.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 1e-12 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// The largest signed distance of a plug point outside the hole hull.
fn excess(c: &[f64; 5]) -> f64 {
    let u = [c[0], c[1]];
    let r = [c[2], c[3], c[4]];
    let n2 = r.iter().map(|x| x * x).sum::<f64>();
    // Cayley rotation R(r) = ((1-|r|²)I + 2rrᵀ + 2[r]×)/(1+|r|²).
    let m = |i: usize, j: usize| {
        let mut v = 2.0 * r[i] * r[j];
        if i == j {
            v += 1.0 - n2;
        }
        let cross = [[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]];
        (v + 2.0 * cross[i][j]) / (1.0 + n2)
    };
    let verts: Vec<[f64; 3]> = vertices().iter().map(|v| [f(&v[0]), f(&v[1]), f(&v[2])]).collect();
    let hole = hull(verts.iter().map(|&v| project(u, v)).collect());
    let plug: Vec<[f64; 2]> = verts
        .iter()
        .map(|v| {
            let p: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| m(i, j) * v[j]).sum());
            project(u, p)
        })
        .collect();
    let mut worst = f64::NEG_INFINITY;
    for p in &plug {
        let mut inside = f64::NEG_INFINITY;
        for i in 0..hole.len() {
            let (a, b) = (hole[i], hole[(i + 1) % hole.len()]);
            let (ex, ey) = (b[0] - a[0], b[1] - a[1]);
            let len = (ex * ex + ey * ey).sqrt();
            // Outward normal of a counter-clockwise polygon: (ey, -ex).
            let d = ((p[0] - a[0]) * ey - (p[1] - a[1]) * ex) / len;
            inside = inside.max(d);
        }
        worst = worst.max(inside);
    }
    worst
}

#[test]
fn the_relation_agrees_with_a_floating_point_hull_where_the_margin_is_clear() {
    // The crate's rotation matrix is the Cayley numerator of (1, r).
    let r = [frac(1, 7), frac(-1, 5), frac(2, 9)];
    let m = [QSqrt5::one(), QSqrt5::from_rational(r[0].clone()), QSqrt5::from_rational(r[1].clone()), QSqrt5::from_rational(r[2].clone())];
    let v = &vertices()[5];
    let exact = rotate(&m, v);
    let rf: Vec<f64> = r.iter().map(|x| f(&QSqrt5::from_rational(x.clone()))).collect();
    let n2: f64 = rf.iter().map(|x| x * x).sum();
    let vf = [f(&v[0]), f(&v[1]), f(&v[2])];
    let cross = [rf[1] * vf[2] - rf[2] * vf[1], rf[2] * vf[0] - rf[0] * vf[2], rf[0] * vf[1] - rf[1] * vf[0]];
    let dot: f64 = (0..3).map(|i| rf[i] * vf[i]).sum();
    for i in 0..3 {
        let expected = (1.0 - n2) * vf[i] + 2.0 * dot * rf[i] + 2.0 * cross[i];
        assert!((f(&exact[i]) - expected).abs() < 1e-9, "rotation convention");
    }
    let root = ConfigurationBox::root();
    let mut state = 12345u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as i64
    };
    let (mut pokes, mut touches) = (0, 0);
    for round in 0..60 {
        let scale = [1, 8, 64][round % 3];
        let values: [Q; 5] = std::array::from_fn(|j| {
            let a = &root.axes()[j];
            let t = frac(next() % 1001, 1000);
            let full = a.lo() + &((a.hi() - a.lo()) * t);
            if j >= 2 { full / q(scale) } else { full }
        });
        let c = Centre::new(values.clone().map(QSqrt5::from_rational)).unwrap();
        let exact = relation(&c);
        let floats: [f64; 5] = std::array::from_fn(|j| f(&QSqrt5::from_rational(values[j].clone())));
        let e = excess(&floats);
        if e > 1e-9 {
            assert_eq!(exact, Relation::Poke, "{values:?} excess {e}");
            pokes += 1;
        } else if e < -1e-9 {
            assert_eq!(exact, Relation::Fit, "{values:?} excess {e}");
        } else {
            touches += 1;
        }
    }
    // r = 0 gives identical shadows: a touch, with zero excess.
    for view in [[frac(1, 5), frac(1, 10)], [frac(2, 3), frac(0, 1)], [frac(1, 3), frac(1, 7)]] {
        let values = [view[0].clone(), view[1].clone(), q(0), q(0), q(0)];
        let c = Centre::new(values.clone().map(QSqrt5::from_rational)).unwrap();
        assert_eq!(relation(&c), Relation::Touch);
        let floats: [f64; 5] = std::array::from_fn(|j| f(&QSqrt5::from_rational(values[j].clone())));
        assert!(excess(&floats).abs() < 1e-9);
    }
    assert!(pokes > 40, "{pokes} pokes, {touches} near touches");
}

/// Whether the documented `cargo test --release` keeps debug assertions
/// (and so the `debug_assert` in `neighbourhood`).
#[test]
fn report_the_profile_of_the_documented_test_command() {
    eprintln!("debug_assertions = {}", cfg!(debug_assertions));
}

// ---- Pilot roots ---------------------------------------------------------

/// Nested explicit roots are refused (their subtree would be searched and
/// counted twice), and so are samples above the documented size.
#[test]
fn nested_pilot_roots_and_huge_samples_are_refused() {
    use rid_validation::experiments::pilot;
    let config = |roots| config::Pilot { roots, threads: two(), depth_limit: 8, max_decisions: None, transcript: "unused".into() };
    for nested in [vec!["0", "01"], vec!["01", "0"], vec!["", "1"], vec!["10", "1011", "11"]] {
        let paths = nested.into_iter().map(String::from).collect();
        assert!(pilot::roots(&config(config::Roots::Paths(paths))).is_err());
    }
    let side_by_side = vec!["00".to_owned(), "01".to_owned(), "1".to_owned()];
    assert_eq!(pilot::roots(&config(config::Roots::Paths(side_by_side))).unwrap(), ["1", "00", "01"]);
    let huge = config::Sample { depth: 62, count: 1 << 62, seed: 1 };
    assert!(pilot::sample(&huge).is_err());
    let largest = config::Sample { depth: 20, count: pilot::MAX_SAMPLE_COUNT, seed: 1 };
    assert_eq!(pilot::sample(&largest).unwrap().len() as u64, pilot::MAX_SAMPLE_COUNT);
}

#[test]
fn irredundancy_summary_tells_small_covers_from_redundant_ones() {
    use rid_validation::experiments::irredundancy::{self, transcript::Verdict};
    let dir = Dir::new("irredundancy");
    let targets = dir.join("targets.json");
    std::fs::write(
        &targets,
        format!(
            r#"{{"format":"rid-validation-targets/1","targets":[{{"id":"square","category":"primary","target":{{"component":"Exotic","cover":"square"}},"centre":{}}}]}}"#,
            rational_centre(["0", "0", "0", "0", "0"])
        ),
    )
    .unwrap();
    let config = config::Irredundancy {
        catalogue: targets.clone(),
        transcript: dir.join("t.jsonl"),
        threads: two(),
        scales: vec![4, 8, 12],
        selection: None,
    };
    let summary = irredundancy::run(&Toy::new(), &config).unwrap();
    let bytes = run::read(&config.transcript).unwrap();
    let lines = run::lines(&bytes).unwrap();
    let result: irredundancy::transcript::TargetResult = run::parse(lines[1], 2).unwrap();
    let verdicts: Vec<Verdict> = result.trials.iter().map(|t| t.verdict).collect();
    assert_eq!(verdicts, [Verdict::Unresolved, Verdict::Exclusive, Verdict::Exclusive]);
    assert_eq!(summary.needed.get("square"), Some(&vec![8, 12]));
    assert!(summary.covered.is_empty() && summary.never.is_empty());
}

/// Labels write covers in decimal, so `9` and `"9"` have one label; the
/// checkers therefore read every record as a typed record first, which
/// refuses the wrong type.
#[test]
fn records_with_covers_of_the_wrong_type_are_refused_before_labels() {
    use rid_validation::experiments::probe::{checked_label, Label};
    assert_eq!(
        Label::of(&json!({"component": "Local", "cover": 9})).unwrap(),
        Label::of(&json!({"component": "Local", "cover": "9"})).unwrap()
    );
    assert!(checked_label(&json!({"component": "Local", "cover": 9}), FORM).is_ok());
    assert!(checked_label(&json!({"component": "Local", "cover": "9"}), FORM).is_err());
    assert!(checked_label(&json!({"component": "Exotic", "cover": 9}), FORM).is_err());
    assert!(checked_label(&json!({"component": "Exotic", "cover": "square"}), FORM).is_ok());
}

// ---- The command line ---------------------------------------------------

/// `compare` exits with status 3 on a regression when asked to, and 0
/// otherwise (the command tests do not cover status 3).
#[test]
fn compare_exits_with_status_3_on_a_regression() {
    let dir = Dir::new("exit");
    let base = completeness_config(&dir, "base.jsonl", 40, Some(vec!["inner-poke"]));
    completeness::run(&Toy { global: frac(1, 64), ..Toy::new() }, &base).unwrap();
    let worse = completeness_config(&dir, "worse.jsonl", 40, Some(vec!["inner-poke"]));
    completeness::run(&Toy { global: frac(1, 33), ..Toy::new() }, &worse).unwrap();
    let status = |fail: bool, baseline: &Path, candidate: &Path| {
        let file = dir.join("compare.json");
        std::fs::write(
            &file,
            json!({"catalogue": base.catalogue, "baseline": baseline, "candidate": candidate,
                   "threads": 2, "worse_by": 1, "fail_on_regression": fail})
            .to_string(),
        )
        .unwrap();
        std::process::Command::new(env!("CARGO_BIN_EXE_rid-validation"))
            .args(["compare", file.to_str().unwrap()])
            .output()
            .unwrap()
            .status
            .code()
            .unwrap()
    };
    // "worse" needs 2^-11 instead of 2^-6: a regression from base to worse.
    assert_eq!(status(true, &base.transcript, &worse.transcript), 3);
    assert_eq!(status(false, &base.transcript, &worse.transcript), 0);
    assert_eq!(status(true, &worse.transcript, &base.transcript), 0);
}
