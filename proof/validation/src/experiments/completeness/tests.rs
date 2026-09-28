use super::compare::{compare, Flag};
use super::schedule::{next, summary, ScheduleError};
use super::transcript::{check, CheckError, Checked, Problem};
use super::*;
use crate::experiments::probe::standin::{boxed, Fault, StandIn, FORM};
use crate::experiments::probe::Label;
use crate::points::catalogue::Centre;
use crate::experiments::run::{lines, parse, read};
use crate::testing::{Rng, Scratch};
use rid::arithmetic::exact::{frac, q};
use serde_json::{json, Value};
use std::num::NonZeroUsize;
use std::path::Path;

// ---- The schedule ------------------------------------------------------

/// The trials the schedule makes for an outcome rule.
fn run_rule(max_k: u32, rule: impl Fn(u32) -> bool) -> Vec<(u32, bool)> {
    let mut trials = Vec::new();
    while let Some(k) = next(max_k, &trials).unwrap() {
        trials.push((k, rule(k)));
    }
    trials
}

#[test]
fn the_schedule_finds_the_threshold_of_monotone_outcomes() {
    for max_k in 1..=40 {
        for threshold in 0..=max_k + 1 {
            let trials = run_rule(max_k, |k| k >= threshold);
            let ks: Vec<u32> = trials.iter().map(|t| t.0).collect();
            let distinct: std::collections::BTreeSet<u32> = ks.iter().copied().collect();
            assert_eq!(distinct.len(), ks.len(), "repeated exponent");
            assert_eq!(ks[0], max_k);
            let (best, finer) = summary(&trials);
            if threshold > max_k {
                assert_eq!(trials, vec![(max_k, false)]);
                assert_eq!(best, None);
            } else {
                assert_eq!(best, Some(threshold), "max_k {max_k}");
                assert!(finer.is_empty());
                for offset in [1, 2, 8] {
                    if threshold + offset <= max_k {
                        assert!(ks.contains(&(threshold + offset)));
                    }
                }
            }
        }
    }
}

#[test]
fn replays_accept_exactly_the_schedules_own_trials() {
    let mut rng = Rng::new(21);
    for _ in 0..500 {
        let max_k = 1 + rng.below(300) as u32;
        let seed = rng.next();
        // An arbitrary, not necessarily monotone, outcome rule.
        let rule = |k: u32| (seed >> (k % 64)) & 1 == 1 || k > max_k / 2;
        let trials = run_rule(max_k, rule);
        assert_eq!(next(max_k, &trials), Ok(None));
        let (best, finer) = summary(&trials);
        if let Some(best) = best {
            assert!(finer.iter().all(|&k| k > best && trials.contains(&(k, false))));
        }
        // Every prefix asks for exactly the next recorded exponent.
        for i in 0..trials.len() {
            assert_eq!(next(max_k, &trials[..i]), Ok(Some(trials[i].0)));
        }
        // Changing an exponent, adding a trial or repeating one is refused.
        let mut changed = trials.clone();
        let i = rng.below(changed.len() as u64) as usize;
        changed[i].0 = if changed[i].0 == 0 { 1 } else { changed[i].0 - 1 };
        assert!(matches!(next(max_k, &changed), Err(ScheduleError::Mismatch { .. }) | Ok(Some(_)) | Err(ScheduleError::Trailing { .. })));
        assert_ne!(next(max_k, &changed), Ok(None));
        let mut longer = trials.clone();
        longer.push(*trials.last().unwrap());
        assert_eq!(next(max_k, &longer), Err(ScheduleError::Trailing { position: trials.len() }));
    }
}

#[test]
fn a_refused_finest_radius_ends_the_schedule() {
    assert_eq!(run_rule(256, |_| false), vec![(256, false)]);
    assert_eq!(summary(&[(256, false)]), (None, vec![]));
    // Refused finer controls are reported, not hidden.
    let trials = run_rule(64, |k| k >= 10 && k != 11);
    assert_eq!(summary(&trials), (Some(10), vec![11]));
}

// ---- Agreement between the decision and the components ---------------

#[test]
fn a_decision_agrees_only_with_the_first_success_and_its_record() {
    let a = json!({"component": "Local", "cover": 1});
    let b = json!({"component": "Local", "cover": 2});
    use transcript::agrees;
    assert!(agrees(None, &[None, Some(None), None]));
    assert!(!agrees(None, &[Some(Some(&a))]));
    assert!(agrees(Some((1, &a)), &[Some(None), Some(Some(&a)), Some(Some(&b))]));
    assert!(agrees(Some((1, &a)), &[None, None, Some(None)]));
    assert!(!agrees(Some((1, &a)), &[Some(Some(&b)), Some(Some(&a))]));
    assert!(!agrees(Some((1, &a)), &[Some(None), Some(None)]));
    assert!(!agrees(Some((1, &a)), &[Some(None), Some(Some(&b))]));
}

// ---- The experiment with the stand-in ----------------------------------

const MAX_K: u32 = 40;

fn centre(values: [&str; 5]) -> String {
    let pairs: Vec<String> = values.iter().map(|v| format!(r#"["{v}","0"]"#)).collect();
    format!("[{}]", pairs.join(","))
}

/// Seven points whose outcomes the stand-in's rules predict.
fn test_catalogue() -> String {
    let points = [
        ("square", "g1", centre(["0", "0", "0", "0", "0"])),
        ("aligned", "g1", centre(["1/5", "1/10", "0", "0", "0"])),
        ("overlap", "g2", centre(["9/32", "1/8", "0", "0", "0"])),
        ("outside", "g2", centre(["2/3", "2/5", "0", "0", "0"])),
        ("poke", "g3", centre(["1/4", "1/8", "1/4", "0", "0"])),
        ("unresolved", "g3", centre(["1/5", "1/5", "1/20", "1/20", "1/20"])),
        (
            "pentagon",
            "g1",
            r#"[["-1/2","3/10"],["1/2","-1/10"],["0","0"],["0","0"],["0","0"]]"#.to_string(),
        ),
    ];
    let lines: Vec<String> = points
        .iter()
        .map(|(id, group, c)| format!(r#"{{"id":"{id}","group":"{group}","centre":{c},"expected":["Local"]}}"#))
        .collect();
    format!(
        r#"{{"format":"rid-validation-points/1","points":[{}]}}"#,
        lines.join(",")
    )
}

fn configuration(dir: &Scratch, name: &str, threads: usize, selection: Option<Vec<&str>>) -> config::Completeness {
    let catalogue = dir.join("points.json");
    if !catalogue.exists() {
        std::fs::write(&catalogue, test_catalogue()).unwrap();
    }
    config::Completeness {
        catalogue,
        transcript: dir.join(name),
        threads: NonZeroUsize::new(threads).unwrap(),
        max_k: MAX_K,
        selection: selection.map(|s| s.into_iter().map(String::from).collect()),
    }
}

fn checked(config: &config::Completeness) -> Result<Checked, CheckError> {
    let catalogue = Catalogue::load(&config.catalogue).unwrap();
    check(&catalogue, &read(&config.transcript).unwrap(), NonZeroUsize::new(2).unwrap(), FORM)
}

fn point<'a>(checked: &'a Checked, id: &str) -> &'a PointResult {
    checked.points.iter().find(|p| p.id == id).unwrap()
}

#[test]
fn the_experiment_finds_the_predicted_radii_and_its_transcript_checks() {
    let dir = Scratch::new("completeness");
    let config = configuration(&dir, "t.jsonl", 3, None);
    let summary = run(&StandIn::typical(), &config).unwrap();
    let result = checked(&config).unwrap();
    assert_eq!(result.summary, summary);
    assert_eq!(summary.points, 7);
    assert_eq!(summary.unresolved, 1);
    assert_eq!(summary.finer_refusals, 0);
    // The square cover holds the clipped cube exactly when 2^-k ≤ 1/64.
    let square = point(&result, "square");
    assert_eq!(square.expected, Some(false), "the catalogue expects Local");
    assert_eq!(point(&result, "aligned").expected, Some(true));
    assert_eq!(point(&result, "unresolved").expected, None);
    assert_eq!(summary.unexpected, 4);
    assert_eq!(square.decision.best_k, Some(6));
    assert_eq!(square.decision.record, Some(json!({"component": "Exotic", "cover": "square"})));
    assert_eq!(square.relation, Relation::Touch);
    assert!(square.in_domain);
    let ks: Vec<u32> = square.decision.trials.iter().map(|t| t.k).collect();
    assert_eq!(ks, [40, 0, 20, 10, 5, 7, 6, 8, 14]);
    let own: Vec<Option<u32>> = square.components.iter().map(|s| s.best_k).collect();
    assert_eq!(own, [None, Some(6), None, None]);
    assert_eq!(square.components[0].trials.len(), 1, "Domain refuses the finest radius");
    // Local's first cover, at 2^-k ≤ 1/128.
    let aligned = point(&result, "aligned");
    assert_eq!(aligned.decision.best_k, Some(7));
    assert_eq!(aligned.decision.record, Some(json!({"component": "Local", "cover": 0})));
    assert_eq!(point(&result, "overlap").decision.record, Some(json!({"component": "Local", "cover": 0})));
    // The pentagon event is irrational; its cover needs 2^-k ≤ 1/64.
    let pentagon = point(&result, "pentagon");
    assert_eq!(pentagon.decision.best_k, Some(6));
    // Outside D: Domain first; Local 1 also holds finer boxes on its own.
    let outside = point(&result, "outside");
    assert!(!outside.in_domain);
    assert_eq!(Label::of(outside.decision.record.as_ref().unwrap()).unwrap().component, "Domain");
    assert_eq!(outside.components[2].best_k, Some(7));
    // A poke with a large rotation: Global on its own from 2^-k ≤ 1/8.
    let poke = point(&result, "poke");
    assert_eq!(poke.relation, Relation::Poke);
    assert_eq!(poke.components[3].best_k, Some(3));
    let unresolved = point(&result, "unresolved");
    assert_eq!(unresolved.decision.best_k, None);
    assert_eq!(unresolved.decision.trials.len(), 1);
    assert_eq!(summary.components["Exotic"], 2);
}

/// The transcript without timings.
fn timeless(path: &Path) -> Vec<Value> {
    fn strip(value: &mut Value) {
        match value {
            Value::Object(map) => {
                map.remove("milliseconds");
                map.values_mut().for_each(strip);
            }
            Value::Array(items) => items.iter_mut().for_each(strip),
            _ => {}
        }
    }
    let bytes = read(path).unwrap();
    lines(&bytes)
        .unwrap()
        .into_iter()
        .map(|line| {
            let mut value: Value = serde_json::from_slice(line).unwrap();
            strip(&mut value);
            value
        })
        .collect()
}

#[test]
fn transcripts_are_the_same_for_every_thread_count() {
    let dir = Scratch::new("threads");
    let reference = configuration(&dir, "1.jsonl", 1, None);
    run(&StandIn::typical(), &reference).unwrap();
    for threads in [2, 5] {
        let config = configuration(&dir, &format!("{threads}.jsonl"), threads, None);
        run(&StandIn::typical(), &config).unwrap();
        assert_eq!(timeless(&config.transcript), timeless(&reference.transcript));
    }
}

#[test]
fn selections_run_in_catalogue_order_and_existing_transcripts_are_kept() {
    let dir = Scratch::new("selection");
    let config = configuration(&dir, "t.jsonl", 2, Some(vec!["pentagon", "square"]));
    run(&StandIn::typical(), &config).unwrap();
    let result = checked(&config).unwrap();
    let ids: Vec<&str> = result.points.iter().map(|p| p.id.as_str()).collect();
    assert_eq!(ids, ["square", "pentagon"]);
    let before = read(&config.transcript).unwrap();
    assert!(matches!(run(&StandIn::typical(), &config), Err(Error::Transcript(_))));
    assert_eq!(read(&config.transcript).unwrap(), before);
    let unknown = configuration(&dir, "u.jsonl", 2, Some(vec!["nowhere"]));
    assert!(matches!(run(&StandIn::typical(), &unknown), Err(Error::Catalogue(_))));
    assert!(!unknown.transcript.exists());
    let mut invalid = configuration(&dir, "k.jsonl", 2, None);
    invalid.max_k = 0;
    assert!(matches!(run(&StandIn::typical(), &invalid), Err(Error::Invalid(_))));
}

#[test]
fn unsound_or_inconsistent_components_stop_the_experiment() {
    let dir = Scratch::new("faults");
    // Global accepting boxes around r = 0 contradicts the touch at the square view.
    let unsound = StandIn {
        fault: Fault::UnsoundGlobal,
        ..StandIn::typical()
    };
    let config = configuration(&dir, "a.jsonl", 2, Some(vec!["square"]));
    match run(&unsound, &config) {
        Err(Error::Contradiction { id, source, .. }) => {
            assert_eq!(id, "square");
            assert_eq!(source, Contradiction::GlobalAtContained);
        }
        other => panic!("{other:?}"),
    }
    // The partial transcript has no summary and does not check.
    assert!(checked(&config).is_err());
    // An attempt its own verification refuses.
    let liar = StandIn {
        fault: Fault::LyingAttempt(2),
        ..StandIn::typical()
    };
    let config = configuration(&dir, "b.jsonl", 2, Some(vec!["square"]));
    assert!(matches!(run(&liar, &config), Err(Error::Probe { source: ProbeError::Refused { .. }, .. })));
    // A collection that does not take the first success: outside D both
    // Domain and Local succeed at fine radii, and the reversed order picks Local.
    let reversed = StandIn {
        fault: Fault::ReversedDecision,
        ..StandIn::typical()
    };
    let config = configuration(&dir, "c.jsonl", 2, Some(vec!["outside"]));
    assert!(matches!(run(&reversed, &config), Err(Error::Disagreement { .. })));
}

// ---- The checker against corrupted transcripts ------------------------

/// The lines of a transcript as typed values.
struct Parts {
    header: Header,
    points: Vec<PointResult>,
    summary: Summary,
}

fn parts(path: &Path) -> Parts {
    let bytes = read(path).unwrap();
    let lines = lines(&bytes).unwrap();
    Parts {
        header: parse(lines[0], 1).unwrap(),
        points: lines[1..lines.len() - 1]
            .iter()
            .enumerate()
            .map(|(i, l)| parse(l, i + 2).unwrap())
            .collect(),
        summary: parse(lines[lines.len() - 1], lines.len()).unwrap(),
    }
}

fn square(p: &mut Parts) -> &mut PointResult {
    p.points.iter_mut().find(|p| p.id == "square").unwrap()
}

fn write_parts(path: &Path, parts: &Parts) {
    let mut bytes = serde_json::to_vec(&parts.header).unwrap();
    bytes.push(b'\n');
    for p in &parts.points {
        bytes.extend(serde_json::to_vec(p).unwrap());
        bytes.push(b'\n');
    }
    bytes.extend(serde_json::to_vec(&parts.summary).unwrap());
    bytes.push(b'\n');
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn every_corruption_of_a_transcript_is_refused() {
    let dir = Scratch::new("corrupt");
    let config = configuration(&dir, "t.jsonl", 2, None);
    run(&StandIn::typical(), &config).unwrap();
    let catalogue = Catalogue::load(&config.catalogue).unwrap();
    let two = NonZeroUsize::new(2).unwrap();
    type Mutation = Box<dyn Fn(&mut Parts)>;
    let mutations: Vec<(&str, Mutation)> = vec![
        ("an end of a box", Box::new(|p| square(p).neighbourhoods[2].axes[0][1] = "1/3".into())),
        ("a malformed end", Box::new(|p| square(p).neighbourhoods[2].axes[0][1] = "2/4".into())),
        ("neighbourhood order", Box::new(|p| square(p).neighbourhoods.reverse())),
        ("an unused neighbourhood", Box::new(|p| {
            let c = Centre::parse(&[["0", "0"], ["0", "0"], ["0", "0"], ["0", "0"], ["0", "0"]].map(|a| a.map(String::from))).unwrap();
            let b = neighbourhood(&c, 39).unwrap();
            let s = square(p);
            s.neighbourhoods.push(Neighbourhood { k: 39, axes: spell(&b) });
            s.neighbourhoods.sort_by_key(|n| n.k);
        })),
        ("a trial's exponent", Box::new(|p| square(p).decision.trials[2].k = 21)),
        ("a dropped trial", Box::new(|p| { square(p).decision.trials.pop(); })),
        ("an added trial", Box::new(|p| {
            let s = square(p);
            let extra = s.decision.trials[1].clone();
            s.decision.trials.push(extra);
        })),
        ("the best exponent", Box::new(|p| square(p).decision.best_k = Some(7))),
        ("the finer refusals", Box::new(|p| square(p).decision.finer_refusals = vec![40])),
        ("the best record", Box::new(|p| square(p).decision.record = None)),
        ("a record of another component", Box::new(|p| {
            square(p).components[1].trials[0].record = Some(json!({"component": "Local", "cover": 0}));
        })),
        ("another cover at the same radius", Box::new(|p| {
            let s = square(p);
            for t in s.decision.trials.iter_mut().filter(|t| t.record.is_some()) {
                t.record = Some(json!({"component": "Exotic", "cover": "pentagon"}));
            }
            s.decision.record = Some(json!({"component": "Exotic", "cover": "pentagon"}));
        })),
        ("a Domain record in D", Box::new(|p| {
            square(p).decision.trials[0].record = Some(json!({"component": "Domain", "inequality": 0}));
        })),
        ("an unknown component", Box::new(|p| {
            square(p).decision.trials[0].record = Some(json!({"component": "Other"}));
        })),
        ("an invalid label", Box::new(|p| square(p).decision.trials[0].record = Some(json!(3)))),
        ("swapped component schedules", Box::new(|p| square(p).components.swap(0, 1))),
        ("a missing component schedule", Box::new(|p| { square(p).components.pop(); })),
        ("the relation", Box::new(|p| square(p).relation = Relation::Poke)),
        ("the expected flag", Box::new(|p| square(p).expected = Some(true))),
        ("a fit", Box::new(|p| square(p).relation = Relation::Fit)),
        ("domain membership", Box::new(|p| square(p).in_domain = false)),
        ("the group", Box::new(|p| square(p).group = "g2".into())),
        ("the identifier", Box::new(|p| square(p).id = "Square".into())),
        ("point order", Box::new(|p| p.points.swap(0, 1))),
        ("a dropped point", Box::new(|p| { p.points.pop(); })),
        ("the summary", Box::new(|p| p.summary.unresolved += 1)),
        ("the summary's reasons", Box::new(|p| { p.summary.reasons.insert("Nothing".into(), 1); })),
        ("the catalogue binding", Box::new(|p| p.header.catalogue_sha256 = "0".repeat(64))),
        ("the finest exponent", Box::new(|p| p.header.max_k = 41)),
        ("an impossible finest exponent", Box::new(|p| p.header.max_k = 0)),
        ("the component names", Box::new(|p| p.header.components[0] = "Exotic".into())),
        ("the format", Box::new(|p| p.header.format = "rid-completeness/2".into())),
        ("the policy", Box::new(|p| p.header.policy = "policy".into())),
        ("the executable digest", Box::new(|p| p.header.executable_sha256.push('0'))),
        ("the selection", Box::new(|p| p.header.selection = Some(vec!["square".into()]))),
    ];
    let original = parts(&config.transcript);
    let target = dir.join("mutated.jsonl");
    write_parts(&target, &original);
    assert!(check(&catalogue, &read(&target).unwrap(), two, FORM).is_ok(), "the rewrite is exact");
    for (what, mutate) in &mutations {
        let mut p = parts(&config.transcript);
        mutate(&mut p);
        write_parts(&target, &p);
        let outcome = check(&catalogue, &read(&target).unwrap(), two, FORM);
        assert!(outcome.is_err(), "{what} accepted");
    }
    // Byte-level damage: a truncated line, extra spacing, a trailing line.
    let bytes = read(&config.transcript).unwrap();
    for damaged in [
        bytes[..bytes.len() - 1].to_vec(),
        String::from_utf8(bytes.clone()).unwrap().replacen(",\"", ", \"", 1).into_bytes(),
        [bytes.clone(), b"{}\n".to_vec()].concat(),
        Vec::new(),
    ] {
        assert!(check(&catalogue, &damaged, two, FORM).is_err());
    }
    // The problems are reported precisely.
    let mut p = parts(&config.transcript);
    square(&mut p).decision.trials[0].record = Some(json!({"component": "Domain", "inequality": 0}));
    write_parts(&target, &p);
    match check(&catalogue, &read(&target).unwrap(), two, FORM) {
        Err(CheckError::Point { id, problem: Problem::Contradiction { source, .. } }) => {
            assert_eq!(id, "square");
            assert_eq!(source, Contradiction::DomainInD);
        }
        other => panic!("{other:?}"),
    }
}

// ---- Comparisons -------------------------------------------------------

fn standin_with_square(side: rid::arithmetic::exact::Q, rotation: rid::arithmetic::exact::Q) -> StandIn {
    let mut standin = StandIn::typical();
    let r = || (-rotation.clone(), rotation.clone());
    standin.exotic[0].1 = boxed([(q(0), side.clone()), (q(0), side.clone()), r(), r(), r()]);
    standin
}

#[test]
fn comparisons_report_changes_per_group_and_regressions() {
    let dir = Scratch::new("compare");
    let base = configuration(&dir, "base.jsonl", 2, None);
    run(&StandIn::typical(), &base).unwrap();
    let base = checked(&base).unwrap();
    // A larger square cover: a larger radius at the square view only.
    let larger = configuration(&dir, "larger.jsonl", 2, None);
    run(&standin_with_square(frac(1, 8), frac(1, 32)), &larger).unwrap();
    let larger = checked(&larger).unwrap();
    let c = compare(&base, &larger, 1).unwrap();
    assert_eq!(c.points, 7);
    assert_eq!(c.regressions, 0);
    assert_eq!(c.changes.len(), 1);
    assert_eq!(c.changes[0].id, "square");
    assert_eq!(c.changes[0].best_k, [Some(6), Some(5)]);
    assert_eq!(c.changes[0].flags, [Flag::LargerRadius]);
    assert_eq!(c.changes[0].components[0].component, "Exotic");
    let g1 = c.groups.iter().find(|g| g.group == "g1").unwrap();
    assert_eq!((g1.points, g1.larger_radius, g1.smaller_radius), (3, 1, 0));
    assert_eq!(c.groups.iter().map(|g| g.group.as_str()).collect::<Vec<_>>(), ["g1", "g2", "g3"]);
    // The reverse direction is a regression by one exponent.
    let back = compare(&larger, &base, 1).unwrap();
    assert_eq!(back.regressions, 1);
    assert_eq!(back.changes[0].flags, [Flag::SmallerRadius]);
    assert_eq!(compare(&larger, &base, 2).unwrap().regressions, 0);
    // Without Local's first cover the aligned point becomes unresolved.
    let fewer = configuration(&dir, "fewer.jsonl", 2, None);
    let mut standin = StandIn::typical();
    standin.local.remove(0);
    run(&standin, &fewer).unwrap();
    let fewer = checked(&fewer).unwrap();
    let c = compare(&base, &fewer, 1).unwrap();
    let aligned = c.changes.iter().find(|p| p.id == "aligned").unwrap();
    assert!(aligned.flags.contains(&Flag::NewlyUnresolved) && aligned.regression);
    let overlap = c.changes.iter().find(|p| p.id == "overlap").unwrap();
    assert_eq!(overlap.reason, [Some("Local 0".into()), Some("Local 1".into())]);
    assert!(overlap.flags.contains(&Flag::ReasonChanged));
    // Partial selections compare their common points; the points missing
    // from the candidate are regressions.
    let part = configuration(&dir, "part.jsonl", 2, Some(vec!["square", "poke"]));
    run(&StandIn::typical(), &part).unwrap();
    let part = checked(&part).unwrap();
    let c = compare(&base, &part, 1).unwrap();
    assert_eq!(c.points, 2);
    assert_eq!(c.only_baseline.len(), 5);
    assert_eq!(c.regressions, 5);
    assert!(c.only_candidate.is_empty() && c.changes.is_empty());
    // The other way round, the extra points are listed and nothing regresses.
    let c = compare(&part, &base, 1).unwrap();
    assert_eq!((c.only_candidate.len(), c.regressions), (5, 0));
}

#[test]
fn the_compare_command_checks_both_transcripts_first() {
    let dir = Scratch::new("compare-command");
    let base = configuration(&dir, "base.jsonl", 2, None);
    run(&StandIn::typical(), &base).unwrap();
    let mut command = config::Compare {
        catalogue: base.catalogue.clone(),
        baseline: base.transcript.clone(),
        candidate: base.transcript.clone(),
        threads: NonZeroUsize::new(2).unwrap(),
        worse_by: 1,
        fail_on_regression: true,
    };
    let same = crate::experiments::compare(&command, FORM).unwrap();
    assert!(same.changes.is_empty() && same.regressions == 0);
    // A damaged candidate is refused before any comparison.
    let damaged = dir.join("damaged.jsonl");
    let mut p = parts(&base.transcript);
    p.summary.points += 1;
    write_parts(&damaged, &p);
    command.candidate = damaged;
    assert!(crate::experiments::compare(&command, FORM).is_err());
    // The check command reports a pass.
    let report = crate::experiments::check(&config::Check {
        catalogue: base.catalogue.clone(),
        transcript: base.transcript.clone(),
        threads: NonZeroUsize::new(2).unwrap(),
    }, FORM)
    .unwrap();
    assert_eq!(report["status"], "pass");
    assert_eq!(report["summary"]["points"], 7);
}

/// The experiment's own cost on the shipped catalogue with the stand-in
/// (whose probes are nearly free): exact properties, neighbourhoods and
/// schedules for 683 points, then the check. The stand-in's Global is
/// disabled: its crude rule would eliminate arc points, which are touches,
/// and the experiment rightly stops on that contradiction.
/// `cargo test --release measure_the_harness_on_the_shipped_catalogue -- --ignored --nocapture`
#[test]
#[ignore]
fn measure_the_harness_on_the_shipped_catalogue() {
    let dir = Scratch::new("measure");
    let threads = std::thread::available_parallelism().unwrap();
    let config = config::Completeness {
        catalogue: crate::testing::shipped("points.json"),
        transcript: dir.join("t.jsonl"),
        threads,
        max_k: 256,
        selection: None,
    };
    let start = std::time::Instant::now();
    let standin = StandIn {
        global: q(1),
        ..StandIn::typical()
    };
    let summary = run(&standin, &config).unwrap();
    let ran = start.elapsed();
    let bytes = read(&config.transcript).unwrap();
    let catalogue = Catalogue::load(&config.catalogue).unwrap();
    let start = std::time::Instant::now();
    check(&catalogue, &bytes, threads, FORM).unwrap();
    eprintln!(
        "{} points on {threads} threads: run {:.1} s, check {:.1} s, transcript {} bytes; {summary:?}",
        summary.points,
        ran.as_secs_f64(),
        start.elapsed().as_secs_f64(),
        bytes.len()
    );
}
