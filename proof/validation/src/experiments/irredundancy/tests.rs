use super::targets::Category;
use super::transcript::{check, CheckError, Problem, Verdict};
use super::*;
use crate::experiments::probe::standin::{Fault, StandIn, FORM};
use crate::experiments::probe::Label;
use crate::experiments::run::{lines, parse, read};
use crate::testing::{shipped, Scratch};
use rid::elimination::proof::catalogue::{local_rectangle, EXOTIC, LOCAL_COUNT};
use serde_json::json;
use std::num::NonZeroUsize;
use std::path::Path;

fn label(component: &str, cover: Option<&str>) -> Label {
    Label {
        component: component.into(),
        cover: cover.map(String::from),
    }
}

#[test]
fn the_shipped_targets_name_every_known_cover_once() {
    let targets = Targets::load(&shipped("targets.json")).unwrap();
    assert_eq!(targets.targets.len(), 15 + LOCAL_COUNT);
    let primary: Vec<String> = targets
        .targets
        .iter()
        .filter(|t| t.category == Category::Primary)
        .map(|t| t.target.to_string())
        .collect();
    // Every component, then every cover of the crate's catalogue, once.
    let mut expected: Vec<String> = ["Global", "Domain", "Local"].map(String::from).to_vec();
    expected.extend(EXOTIC.iter().map(|name| format!("Exotic {name}")));
    expected.extend((0..LOCAL_COUNT).map(|n| format!("Local {n}")));
    assert_eq!(primary, expected);
    // The exact properties the targets were chosen for.
    for t in &targets.targets {
        let relation = relation::relation(&t.centre);
        let in_domain = relation::in_domain(&t.centre);
        match t.id.as_str() {
            "global" => assert_eq!((in_domain, relation), (true, Relation::Poke)),
            "domain" | "second-sheet" => assert_eq!((in_domain, relation), (false, Relation::Touch)),
            "outside-poke" => assert_eq!((in_domain, relation), (false, Relation::Poke)),
            "pentagon-rotation" => assert_eq!((in_domain, relation), (true, Relation::Poke)),
            _ => assert_eq!((in_domain, relation), (true, Relation::Touch), "{}", t.id),
        }
    }
}

#[test]
fn the_local_cover_targets_are_generated_from_the_crate_catalogue() {
    // The shipped catalogue is the generator's fixed point.
    let bytes = std::fs::read(shipped("targets.json")).unwrap();
    assert_eq!(targets::regenerate(&bytes).unwrap(), bytes);
    // Each representative lies in D at r = 0, strictly inside its own view
    // rectangle and outside every other Local cover's (closed) rectangle.
    let shipped = Targets::parse(&bytes).unwrap();
    for n in 0..LOCAL_COUNT {
        let centre = targets::local_representative(n).unwrap();
        let t = shipped.targets.iter().find(|t| t.id == targets::local_id(n)).unwrap();
        assert_eq!((&t.centre, t.category), (&centre, Category::Primary));
        assert_eq!(t.target, label("Local", Some(&n.to_string())));
        let [s, tt] = local_rectangle(n).unwrap();
        let x = centre.coordinates();
        assert!(x[2..].iter().all(|r| r.is_zero()));
        assert!(relation::in_domain(&centre));
        let inside = |v: &rid::arithmetic::exact::QSqrt5, [lo, hi]: &[rid::arithmetic::exact::Q; 2]| {
            v.sqrt5_part().is_zero() && v.rational_part() > lo && v.rational_part() < hi
        };
        assert!(inside(&x[0], &s) && inside(&x[1], &tt), "{n}");
        for m in (0..LOCAL_COUNT).filter(|&m| m != n) {
            let [s, t] = local_rectangle(m).unwrap();
            let within = |v: &rid::arithmetic::exact::QSqrt5, [lo, hi]: &[rid::arithmetic::exact::Q; 2]| {
                v.rational_part() >= lo && v.rational_part() <= hi
            };
            assert!(!(within(&x[0], &s) && within(&x[1], &t)), "{n} in {m}");
        }
    }
    // Stale or missing Local targets are replaced, the others kept in order.
    let stale = catalogue(&[
        target("a", "control", r#"{"component":"Local"}"#, ["1/5", "1/10", "0", "0", "0"]),
        target("old", "primary", r#"{"component":"Local","cover":3}"#, ["1/5", "1/10", "0", "0", "0"]),
        target("b", "primary", r#"{"component":"Domain"}"#, ["2/3", "2/5", "0", "0", "0"]),
    ]);
    let fresh = Targets::parse(&targets::regenerate(&stale).unwrap()).unwrap();
    let ids: Vec<&str> = fresh.targets.iter().map(|t| t.id.as_str()).collect();
    assert_eq!(ids[..2], ["a", "b"]);
    assert_eq!(ids.len(), 2 + LOCAL_COUNT);
    assert!(ids[2..].iter().zip(0..).all(|(id, n)| *id == targets::local_id(n)));
    assert!(targets::regenerate(b"{}").is_err());
}

fn catalogue(targets: &[String]) -> Vec<u8> {
    format!(
        r#"{{"format":"rid-validation-targets/1","targets":[{}]}}"#,
        targets.join(",")
    )
    .into_bytes()
}

fn target(id: &str, category: &str, target: &str, centre: [&str; 5]) -> String {
    let pairs: Vec<String> = centre.iter().map(|v| format!(r#"["{v}","0"]"#)).collect();
    format!(
        r#"{{"id":"{id}","category":"{category}","target":{target},"centre":[{}]}}"#,
        pairs.join(",")
    )
}

#[test]
fn malformed_target_catalogues_are_refused() {
    let origin = ["0", "0", "0", "0", "0"];
    let good = target("a", "primary", r#"{"component":"Local","cover":3}"#, origin);
    let parsed = Targets::parse(&catalogue(&[good.clone()])).unwrap();
    assert_eq!(parsed.targets[0].target, label("Local", Some("3")));
    for (bad, what) in [
        (target("a", "primary", r#"{"component":"Local","extra":1}"#, origin), "extra key"),
        (target("a", "primary", r#"{"component":""}"#, origin), "empty component"),
        (target("a", "primary", r#"{"component":"Local","cover":{}}"#, origin), "cover object"),
        (target("a", "primary", r#""Local""#, origin), "not an object"),
        (target("a", "secondary", r#"{"component":"Local"}"#, origin), "category"),
        (target("a", "primary", r#"{"component":"Local"}"#, ["1", "0", "0", "0", "0"]), "outside the root"),
    ] {
        assert!(Targets::parse(&catalogue(&[bad])).is_err(), "{what} accepted");
    }
    assert!(Targets::parse(&catalogue(&[good.clone(), good.clone()])).is_err());
    assert!(Targets::parse(&catalogue(&[])).is_err());
    let wrong = String::from_utf8(catalogue(&[good])).unwrap().replace("targets/1", "targets/2");
    assert!(Targets::parse(wrong.as_bytes()).is_err());
}

#[test]
fn verdicts_compare_the_successes_with_the_target() {
    let local = label("Local", None);
    let local0 = label("Local", Some("0"));
    let local1 = label("Local", Some("1"));
    let global = label("Global", None);
    assert_eq!(verdict(&local, &[local0.clone(), local1.clone()]), Verdict::Exclusive);
    assert_eq!(verdict(&local0, &[local0.clone()]), Verdict::Exclusive);
    assert_eq!(verdict(&local0, &[local0.clone(), local1.clone()]), Verdict::Shared);
    assert_eq!(verdict(&local0, &[local1.clone()]), Verdict::Other);
    assert_eq!(verdict(&local, &[global.clone()]), Verdict::Other);
    assert_eq!(verdict(&local, &[]), Verdict::Unresolved);
}

const SCALES: [u32; 3] = [12, 16, 32];

fn test_targets() -> Vec<u8> {
    catalogue(&[
        target("square", "primary", r#"{"component":"Exotic","cover":"square"}"#, ["0", "0", "0", "0", "0"]),
        target("overlap-any", "primary", r#"{"component":"Local"}"#, ["9/32", "1/8", "0", "0", "0"]),
        target("overlap-0", "control", r#"{"component":"Local","cover":0}"#, ["9/32", "1/8", "0", "0", "0"]),
        target("outside", "primary", r#"{"component":"Domain"}"#, ["2/3", "2/5", "0", "0", "0"]),
        target("poke", "control", r#"{"component":"Global"}"#, ["1/4", "1/8", "1/4", "0", "0"]),
        target("nothing", "primary", r#"{"component":"Global"}"#, ["1/5", "1/5", "1/20", "1/20", "1/20"]),
    ])
}

fn configuration(dir: &Scratch, name: &str, threads: usize) -> config::Irredundancy {
    let catalogue = dir.join("targets.json");
    if !catalogue.exists() {
        std::fs::write(&catalogue, test_targets()).unwrap();
    }
    config::Irredundancy {
        catalogue,
        transcript: dir.join(name),
        threads: NonZeroUsize::new(threads).unwrap(),
        scales: SCALES.to_vec(),
        selection: None,
    }
}

fn checked(config: &config::Irredundancy) -> Result<transcript::Checked, CheckError> {
    let targets = Targets::load(&config.catalogue).unwrap();
    check(&targets, &read(&config.transcript).unwrap(), NonZeroUsize::new(2).unwrap(), FORM)
}

#[test]
fn the_experiment_asks_every_probe_and_finds_the_predicted_verdicts() {
    let dir = Scratch::new("irredundancy");
    let config = configuration(&dir, "t.jsonl", 3);
    let summary = run(&StandIn::typical(), &config).unwrap();
    let result = checked(&config).unwrap();
    assert_eq!(result.summary, summary);
    let probes: Vec<String> = result.header.probes.iter().map(Label::to_string).collect();
    assert_eq!(probes, ["Domain", "Exotic square", "Exotic pentagon", "Local 0", "Local 1", "Global"]);
    let verdicts = |id: &str| -> Vec<Verdict> {
        result.targets.iter().find(|t| t.id == id).unwrap().trials.iter().map(|t| t.verdict).collect()
    };
    assert_eq!(verdicts("square"), [Verdict::Exclusive; 3]);
    assert_eq!(verdicts("overlap-any"), [Verdict::Exclusive; 3]);
    assert_eq!(verdicts("overlap-0"), [Verdict::Shared; 3]);
    assert_eq!(verdicts("outside"), [Verdict::Shared; 3]);
    assert_eq!(verdicts("poke"), [Verdict::Shared; 3]);
    assert_eq!(verdicts("nothing"), [Verdict::Unresolved; 3]);
    let square = &result.targets[0].trials[0];
    let found: Vec<bool> = square.outcomes.iter().map(|o| o.record.is_some()).collect();
    assert_eq!(found, [false, true, false, false, false, false]);
    assert_eq!(square.outcomes[1].record, Some(json!({"component": "Exotic", "cover": "square"})));
    let needed: Vec<(&str, Vec<u32>)> = summary.needed.iter().map(|(id, k)| (id.as_str(), k.clone())).collect();
    let scales = config.scales.clone();
    assert_eq!(needed, [("overlap-any", scales.clone()), ("square", scales)]);
    assert_eq!(summary.covered, ["outside"]);
    assert_eq!(summary.never, ["nothing"]);
    assert_eq!(summary.trials, 18);
    assert_eq!(summary.verdicts["control"]["shared"], 6);
    // Every thread count writes the same transcript up to timings.
    let other = configuration(&dir, "one.jsonl", 1);
    run(&StandIn::typical(), &other).unwrap();
    let strip = |path: &Path| -> Vec<transcript::TargetResult> {
        checked_at(&config.catalogue, path)
            .targets
            .into_iter()
            .map(|mut t| {
                t.milliseconds = 0;
                for trial in &mut t.trials {
                    for o in &mut trial.outcomes {
                        o.milliseconds = 0;
                    }
                }
                t
            })
            .collect()
    };
    assert_eq!(strip(&config.transcript), strip(&other.transcript));
}

fn checked_at(catalogue: &Path, transcript: &Path) -> transcript::Checked {
    let targets = Targets::load(catalogue).unwrap();
    check(&targets, &read(transcript).unwrap(), NonZeroUsize::new(2).unwrap(), FORM).unwrap()
}

#[test]
fn unsound_components_stop_the_experiment() {
    let dir = Scratch::new("irredundancy-faults");
    let unsound = StandIn {
        fault: Fault::UnsoundGlobal,
        ..StandIn::typical()
    };
    let config = configuration(&dir, "t.jsonl", 2);
    match run(&unsound, &config) {
        Err(Error::Contradiction { id, source, .. }) => {
            assert_eq!(id, "square");
            assert_eq!(source, Contradiction::GlobalAtContained);
        }
        other => panic!("{other:?}"),
    }
    let liar = StandIn {
        fault: Fault::LyingAttempt(3),
        ..StandIn::typical()
    };
    let config = configuration(&dir, "u.jsonl", 2);
    assert!(matches!(run(&liar, &config), Err(Error::Probe { .. })));
    let mut bad = configuration(&dir, "v.jsonl", 2);
    bad.scales = vec![16, 12];
    assert!(matches!(run(&StandIn::typical(), &bad), Err(Error::Invalid(_))));
    assert!(!bad.transcript.exists());
}

type Lines = (transcript::Header, Vec<transcript::TargetResult>, transcript::Summary);

fn parts(path: &Path) -> Lines {
    let bytes = read(path).unwrap();
    let lines = lines(&bytes).unwrap();
    (
        parse(lines[0], 1).unwrap(),
        lines[1..lines.len() - 1].iter().map(|l| parse(l, 2).unwrap()).collect(),
        parse(lines[lines.len() - 1], lines.len()).unwrap(),
    )
}

fn write_parts(path: &Path, (header, targets, summary): &Lines) {
    let mut out = Vec::new();
    out.extend(serde_json::to_vec(header).unwrap());
    out.push(b'\n');
    for t in targets {
        out.extend(serde_json::to_vec(t).unwrap());
        out.push(b'\n');
    }
    out.extend(serde_json::to_vec(summary).unwrap());
    out.push(b'\n');
    std::fs::write(path, out).unwrap();
}

#[test]
fn every_corruption_of_a_transcript_is_refused() {
    let dir = Scratch::new("irredundancy-corrupt");
    let config = configuration(&dir, "t.jsonl", 2);
    run(&StandIn::typical(), &config).unwrap();
    let target = dir.join("mutated.jsonl");
    write_parts(&target, &parts(&config.transcript));
    assert!(checked_at(&config.catalogue, &target).summary.targets == 6);
    type Mutation = Box<dyn Fn(&mut Lines)>;
    let mutations: Vec<(&str, Mutation)> = vec![
        ("a verdict", Box::new(|l| l.1[0].trials[0].verdict = Verdict::Shared)),
        ("a record moved to another cover", Box::new(|l| {
            let outcomes = &mut l.1[0].trials[0].outcomes;
            outcomes[2].record = outcomes[1].record.take();
        })),
        ("an added success", Box::new(|l| {
            l.1[0].trials[0].outcomes[3].record = Some(json!({"component": "Local", "cover": 0}));
        })),
        ("a Domain success in D", Box::new(|l| {
            l.1[0].trials[0].outcomes[0].record = Some(json!({"component": "Domain", "inequality": 0}));
            l.1[0].trials[0].verdict = Verdict::Shared;
        })),
        ("a Global success at a touch", Box::new(|l| {
            l.1[0].trials[0].outcomes[5].record = Some(json!({"component": "Global", "witness": 2}));
            l.1[0].trials[0].verdict = Verdict::Shared;
        })),
        ("a missing outcome", Box::new(|l| { l.1[0].trials[0].outcomes.pop(); })),
        ("a missing scale", Box::new(|l| { l.1[0].trials.pop(); })),
        ("a box", Box::new(|l| l.1[0].trials[1].axes = l.1[0].trials[0].axes.clone())),
        ("the target", Box::new(|l| l.1[1].target = label("Local", Some("1")))),
        ("the category", Box::new(|l| l.1[1].category = Category::Control)),
        ("the relation", Box::new(|l| l.1[0].relation = Relation::Poke)),
        ("domain membership", Box::new(|l| l.1[0].in_domain = false)),
        ("target order", Box::new(|l| l.1.swap(0, 1))),
        ("reordered probes", Box::new(|l| l.0.probes.swap(0, 5))),
        ("a repeated probe", Box::new(|l| l.0.probes[2] = l.0.probes[1].clone())),
        ("an attempt beside covers", Box::new(|l| l.0.probes[1] = label("Exotic", None))),
        ("an unknown probe", Box::new(|l| l.0.probes[5] = label("Other", None))),
        ("the scales", Box::new(|l| l.0.scales = vec![12, 16, 33])),
        ("decreasing scales", Box::new(|l| l.0.scales = vec![16, 12, 32])),
        ("the summary", Box::new(|l| { l.2.covered.pop(); })),
        ("the catalogue binding", Box::new(|l| l.0.catalogue_sha256 = "1".repeat(64))),
        ("the policy", Box::new(|l| l.0.policy = "A".repeat(64))),
    ];
    for (what, mutate) in &mutations {
        let mut l = parts(&config.transcript);
        mutate(&mut l);
        write_parts(&target, &l);
        let targets = Targets::load(&config.catalogue).unwrap();
        let outcome = check(&targets, &read(&target).unwrap(), NonZeroUsize::new(2).unwrap(), FORM);
        assert!(outcome.is_err(), "{what} accepted");
    }
    let mut l = parts(&config.transcript);
    l.1[0].trials[0].outcomes[5].record = Some(json!({"component": "Global", "witness": 2}));
    l.1[0].trials[0].verdict = Verdict::Shared;
    write_parts(&target, &l);
    match checked(&config::Irredundancy { transcript: target, ..configuration(&dir, "unused", 1) }) {
        Err(CheckError::Target { problem: Problem::Contradiction { source, .. }, .. }) => {
            assert_eq!(source, Contradiction::GlobalAtContained)
        }
        other => panic!("{other:?}"),
    }
}

/// The experiment's own cost on the shipped targets at the nine scales up to
/// 2^-2048, with the stand-in (its Global disabled, as it would eliminate
/// the arc touches).
/// `cargo test --release measure_the_harness_on_the_shipped_targets -- --ignored --nocapture`
#[test]
#[ignore]
fn measure_the_harness_on_the_shipped_targets() {
    let dir = Scratch::new("measure-targets");
    let threads = std::thread::available_parallelism().unwrap();
    let config = config::Irredundancy {
        catalogue: shipped("targets.json"),
        transcript: dir.join("t.jsonl"),
        threads,
        scales: vec![12, 16, 32, 64, 128, 256, 512, 1024, 2048],
        selection: None,
    };
    let start = std::time::Instant::now();
    let standin = StandIn {
        global: rid::arithmetic::exact::q(1),
        ..StandIn::typical()
    };
    let summary = run(&standin, &config).unwrap();
    let ran = start.elapsed();
    let start = std::time::Instant::now();
    checked(&config).unwrap();
    eprintln!(
        "{} targets on {threads} threads: run {:.1} s, check {:.1} s; {summary:?}",
        summary.targets,
        ran.as_secs_f64(),
        start.elapsed().as_secs_f64()
    );
}
