//! End-to-end tests of the commands, run in this process with the stand-in
//! components (and, for the production binary, as a real process).

mod support;

use rid::search::certificate::RecordData;
use rid::search::command::{entry, Components, FAILURE, SUCCESS, USAGE_ERROR};
use rid::search::BoxError;
use rid::POLICY;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::Digest;
use std::ffi::OsString;
use std::fs;
use std::process::Command as Process;
use std::sync::atomic::AtomicBool;
use std::thread;
use std::time::{Duration, Instant};
use support::*;

#[test]
fn prepare_reports_the_policy_and_the_components() {
    let _files = files();
    let dir = Folder::new();
    let summary = succeed(&dir.args("prepare", json!({ "threads": 2 })), 7);
    assert_eq!(
        summary,
        json!({ "command": "prepare", "policy": POLICY, "components": { "stand_in": 7 } })
    );
}

#[test]
fn a_search_completes_and_checks_complete() {
    let _files = files();
    let dir = Folder::new();
    let summary = succeed(&dir.search(3, 8, 64, None), 1);
    assert_eq!(summary["command"], "search");
    assert_eq!(summary["halt"], "complete");
    assert_eq!(summary["unresolved"], 0);
    assert_eq!((&summary["max_depth"], &summary["depth_limit"]), (&json!(64), &json!(64)));
    assert_eq!(summary["records"], summary["new_records"]);
    assert_eq!(summary["policy"], POLICY);
    let before = snapshot(&dir.certificate());
    assert_eq!(summary["certificate_bytes"], before.0.len());
    let checked = succeed(&dir.check(64), 1);
    assert_eq!(
        (checked["complete"].clone(), checked["records"].clone(), checked["frontier"].clone()),
        (json!(true), summary["records"].clone(), json!(0))
    );
    assert_eq!(snapshot(&dir.certificate()), before);
    // A completed search is a no-op when repeated.
    let again = succeed(&dir.search(2, 3, 64, None), 1);
    assert_eq!((again["halt"].clone(), again["decisions"].clone()), (json!("complete"), json!(0)));
    assert_eq!(fs::read(dir.certificate()).unwrap(), before.0);
}

#[test]
fn every_thread_count_writes_the_same_certificate() {
    let _files = files();
    for seed in 0..4 {
        let expected = reference(seed, 64);
        for (threads, window) in [(1, 4), (2, 2), (3, 1), (3, 9), (4, 64), (8, 3)] {
            let dir = Folder::new();
            succeed(&dir.search(threads, window, 64, None), seed);
            let bytes = fs::read(dir.certificate()).unwrap();
            assert_eq!(bytes, expected, "seed {seed} threads {threads}");
        }
    }
}

#[test]
fn bounded_searches_resume_to_the_same_certificate() {
    let _files = files();
    let expected = reference(3, 64);
    let dir = Folder::new();
    for attempt in 1.. {
        assert!(attempt < 200, "no progress");
        let budget = Some(10 * attempt as u64);
        let summary = succeed(&dir.search(1 + attempt % 3, 1 + attempt % 5, 64, budget), 3);
        let bytes = fs::read(dir.certificate()).unwrap();
        assert!(expected.starts_with(&bytes));
        if summary["halt"] == "complete" {
            break;
        }
        assert_eq!(summary["halt"], "decision-limit");
        let checked = succeed(&dir.check(64), 3);
        assert_eq!(checked["complete"], false);
    }
    assert_eq!(fs::read(dir.certificate()).unwrap(), expected);
}

#[test]
fn the_header_limit_is_part_of_the_certificate() {
    let _files = files();
    let dir = Folder::new();
    let summary = succeed(&dir.search(3, 4, 5, None), 2);
    assert_eq!(summary["halt"], "depth-limit");
    assert!(summary["unresolved"].as_u64().unwrap() > 0);
    assert_eq!((&summary["max_depth"], &summary["depth_limit"]), (&json!(5), &json!(5)));
    let checked = succeed(&dir.check(5), 2);
    assert_eq!(checked["complete"], false);
    // The search counts its unresolved leaves; the check, the minimal boxes
    // they make up (refused siblings merge).
    let frontier = checked["frontier"].as_u64().unwrap();
    assert!(frontier > 0 && frontier <= summary["unresolved"].as_u64().unwrap());
    // The header fixes the limit: another one refuses without a change.
    let before = snapshot(&dir.certificate());
    let deeper = dir.search(3, 4, 64, None);
    let (status, summary) = run(&deeper, &AtomicBool::new(false), StandIn::new(2));
    assert_eq!((status, summary), (FAILURE, None));
    assert_eq!(run(&dir.check(64), &AtomicBool::new(false), StandIn::new(2)).0, FAILURE);
    assert_eq!(snapshot(&dir.certificate()), before);
}

/// A pilot search with a small `depth_limit` writes exactly the records of
/// the full search up to that depth; later runs deepen the same file in
/// place, with any limits, budgets and thread counts, and end with exactly
/// the certificate of one uninterrupted search.
#[test]
fn a_depth_limited_pilot_is_deepened_in_place() {
    let _files = files();
    for seed in [5, 8] {
        let expected = reference(seed, 64);
        let dir = Folder::new();
        let pilot = succeed(&dir.search_to(3, 4, 64, Some(4), None), seed);
        assert_eq!((&pilot["halt"], &pilot["depth_limit"]), (&json!("depth-limit"), &json!(4)));
        let bytes = fs::read(dir.certificate()).unwrap();
        assert!(expected.starts_with(&bytes) && bytes.len() < expected.len());
        // Every record of the pilot is at most 4 deep.
        let text = String::from_utf8(bytes).unwrap();
        for line in text.lines().skip(1) {
            let payload: Value = serde_json::from_str(line.split_once(' ').unwrap().1).unwrap();
            assert!(payload["path"].as_str().unwrap().len() <= 4, "{line}");
        }
        assert_eq!(succeed(&dir.check(64), seed)["complete"], false);
        // Deeper, shallower again, budgeted, and finally unlimited.
        let runs = [(Some(7), None), (Some(3), None), (Some(9), Some(15)), (None, Some(40))];
        for (attempt, (limit, budget)) in runs.into_iter().enumerate() {
            let summary = succeed(&dir.search_to(1 + attempt % 3, 2, 64, limit, budget), seed);
            if budget.is_none() {
                // The stand-in accepts every box at depth 10, not before.
                assert_eq!(summary["halt"], "depth-limit", "{summary}");
            }
            let bytes = fs::read(dir.certificate()).unwrap();
            assert!(expected.starts_with(&bytes), "seed {seed} run {attempt}");
        }
        loop {
            let summary = succeed(&dir.search_to(2, 3, 64, None, Some(50)), seed);
            if summary["halt"] == "complete" {
                break;
            }
        }
        assert_eq!(fs::read(dir.certificate()).unwrap(), expected);
        assert_eq!(succeed(&dir.check(64), seed)["complete"], true);
    }
}

#[test]
fn false_or_foreign_certificates_are_refused_without_change() {
    let _files = files();
    let dir = Folder::new();
    succeed(&dir.search(3, 4, 64, Some(40)), 5);
    let before = snapshot(&dir.certificate());
    // Another decision rule (seed) does not accept these records.
    for args in [dir.check(64), dir.search(3, 4, 64, None)] {
        assert_eq!(run(&args, &AtomicBool::new(false), StandIn::new(6)).0, FAILURE);
    }
    assert_eq!(snapshot(&dir.certificate()), before);
    // A record whose data was replaced (with a matching checksum) is false.
    let text = String::from_utf8(before.0.clone()).unwrap();
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let (_, payload) = lines[2].split_once(' ').unwrap();
    let (head, tail) = payload.split_once("\"weight\":").unwrap();
    let digits = tail.bytes().take_while(u8::is_ascii_digit).count();
    let weight: u64 = tail[..digits].parse().unwrap();
    let payload = format!("{head}\"weight\":{}{}", weight + 1, &tail[digits..]);
    let digest = sha2::Sha256::digest(payload.as_bytes());
    let checksum: String = digest[..4].iter().map(|b| format!("{b:02x}")).collect();
    lines[2] = format!("{checksum} {payload}");
    fs::write(dir.certificate(), lines.join("\n") + "\n").unwrap();
    let forged = snapshot(&dir.certificate());
    for args in [dir.check(64), dir.search(3, 4, 64, None)] {
        assert_eq!(run(&args, &AtomicBool::new(false), StandIn::new(5)).0, FAILURE);
    }
    assert_eq!(snapshot(&dir.certificate()), forged);
}

/// A component whose `check` disagrees with its `holds` cannot make a false
/// certificate pass: the search writes what `check` returns (the store does
/// not verify appends), but every reading runs `holds` on every record, so
/// `check` and a resumed search refuse the file without changing it, and a
/// corrected component refuses it as well.
#[test]
fn a_record_its_own_holds_refuses_is_never_accepted() {
    let _files = files();
    let expected = reference(4, 64);
    let text = String::from_utf8(expected.clone()).unwrap();
    // A path that the stand-in accepts, somewhere in the middle.
    let lines: Vec<&str> = text.lines().collect();
    let line = lines[lines.len() / 2];
    let payload: Value = serde_json::from_str(line.split_once(' ').unwrap().1).unwrap();
    let path = payload["path"].as_str().unwrap().to_owned();
    for threads in [1, 3] {
        let dir = Folder::new();
        let mut stand_in = StandIn::new(4);
        stand_in.false_record = Some(path.clone());
        let args = dir.search(threads, 4, 64, None);
        let (status, summary) = run(&args, &AtomicBool::new(false), stand_in);
        assert_eq!(status, SUCCESS);
        assert_eq!(summary.unwrap()["halt"], "complete");
        let written = snapshot(&dir.certificate());
        assert!(String::from_utf8(fs::read(dir.certificate()).unwrap()).unwrap().contains(&format!("\"{path}\"")));
        assert_ne!(fs::read(dir.certificate()).unwrap(), expected);
        for args in [dir.check(64), dir.search(threads, 4, 64, None)] {
            assert_eq!(run(&args, &AtomicBool::new(false), StandIn::new(4)).0, FAILURE);
        }
        assert_eq!(snapshot(&dir.certificate()), written);
    }
}

/// A search configuration with the given header fields.
fn search_with(
    dir: &Folder,
    root: &str,
    max_depth: u64,
    depth_limit: Option<u64>,
) -> Vec<OsString> {
    dir.args(
        "search",
        json!({
            "certificate": dir.certificate(),
            "root": root,
            "threads": 1,
            "window": 1,
            "max_depth": max_depth,
            "depth_limit": depth_limit,
            "max_decisions": null,
        }),
    )
}

#[test]
fn invalid_invocations_are_refused_before_any_file_is_created() {
    let _files = files();
    let dir = Folder::new();
    let stop = AtomicBool::new(false);
    let good = dir.search(2, 2, 64, None);
    let usage: [Vec<OsString>; 5] = [
        vec![],
        vec!["search".into()],
        vec![good[0].clone(), good[1].clone(), "extra".into()],
        vec!["run".into(), good[1].clone()],
        vec!["Search".into(), good[1].clone()],
    ];
    for args in usage {
        assert_eq!(run(&args, &stop, StandIn::new(0)), (USAGE_ERROR, None), "{args:?}");
    }
    let mut missing_limit = json!({
        "certificate": dir.certificate(), "root": "", "threads": 1, "window": 1,
        "max_depth": 64, "depth_limit": null, "max_decisions": null,
    });
    missing_limit.as_object_mut().unwrap().remove("depth_limit");
    let failures = [
        vec!["search".into(), dir.0.join("missing.json").into()],
        dir.args("search", json!({ "certificate": dir.certificate(), "root": "" })),
        dir.args("search", missing_limit),
        dir.args(
            "check",
            json!({ "certificate": dir.certificate(), "root": "", "max_depth": 64, "threads": 0 }),
        ),
        search_with(&dir, "2", 64, None),
        search_with(&dir, "0101", 3, None),
        search_with(&dir, "", 4097, None),
        search_with(&dir, "", 64, Some(65)),
        search_with(&dir, "", 64, Some(4097)),
        search_with(&dir, "0101", 64, Some(2)),
        search_with(&dir, "0101", 64, Some(0)),
        dir.check(64),
    ];
    for args in failures {
        assert_eq!(run(&args, &stop, StandIn::new(0)), (FAILURE, None), "{args:?}");
    }
    assert!(!dir.certificate().exists());
    // A stop requested before the command starts changes nothing either.
    let stopped = AtomicBool::new(true);
    assert_eq!(run(&dir.search(1, 1, 64, None), &stopped, StandIn::new(0)), (FAILURE, None));
    assert!(!dir.certificate().exists());
}

/// The production binary refuses malformed invocations before it prepares
/// its components or creates any file. (Its runs with the components are in
/// components.rs.)
#[test]
fn the_production_binary_refuses_malformed_invocations_before_any_file_exists() {
    let _processes = processes();
    let dir = Folder::new();
    let rid = env!("CARGO_BIN_EXE_rid");
    let wrong = Process::new(rid).args(["search"]).output().unwrap();
    assert_eq!(wrong.status.code(), Some(i32::from(USAGE_ERROR)));
    let invalid = dir.args("check", json!({ "certificate": dir.certificate() }));
    let output = Process::new(rid).args(&invalid).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(FAILURE)));
    assert!(!dir.certificate().exists());
}

/// `rid generate` through the binary: a usage error without a
/// configuration, a failure for an unknown cover (before any file exists),
/// and a small run (Local cover 21) that writes a file the cover loader
/// accepts.
#[test]
fn the_production_binary_generates_a_cover() {
    let _processes = processes();
    let dir = Folder::new();
    let rid = env!("CARGO_BIN_EXE_rid");
    let usage = Process::new(rid).args(["generate"]).output().unwrap();
    assert_eq!(usage.status.code(), Some(i32::from(USAGE_ERROR)));
    let output = dir.0.join("generated");
    let config = |numbers: Value| {
        json!({
            "output": output,
            "threads": 2,
            "search": {"max_depth": 60, "max_nodes": 200000, "aspect": 16, "tried": 6,
                       "proposal_depth": 4, "proposals": 8, "inherited": 16},
            "exotic": [],
            "local": {"numbers": numbers, "samples": 1}
        })
    };
    let file = dir.0.join("generate.json");
    std::fs::write(&file, config(json!([49])).to_string()).unwrap();
    let unknown = Process::new(rid).arg("generate").arg(&file).output().unwrap();
    assert_eq!(unknown.status.code(), Some(i32::from(FAILURE)));
    assert!(!output.exists());
    std::fs::write(&file, config(json!([21])).to_string()).unwrap();
    let run = Process::new(rid).arg("generate").arg(&file).output().unwrap();
    assert_eq!(run.status.code(), Some(0), "{}", String::from_utf8_lossy(&run.stderr));
    let summary: Value = serde_json::from_slice(&run.stdout).unwrap();
    assert_eq!(summary["covers"][0]["written"], true);
    let bytes = std::fs::read(output.join("local/21.json")).unwrap();
    let pin = rid::elimination::proof::catalogue::local_pin(21).unwrap();
    let stop = AtomicBool::new(false);
    rid::elimination::proof::ProvedCover::load(&pin, &bytes, std::num::NonZeroUsize::new(2).unwrap(), &stop).unwrap();
}

#[test]
fn a_subtree_root_is_searched_and_checked_on_its_own() {
    let _files = files();
    let dir = Folder::new();
    let config_to = |command: &str, depth_limit: Option<u64>| {
        let mut config = json!({
            "certificate": dir.certificate(),
            "root": "0110",
            "max_depth": 64,
            "threads": 2,
        });
        if command == "search" {
            config["window"] = json!(4);
            config["depth_limit"] = json!(depth_limit);
            config["max_decisions"] = json!(null);
        }
        dir.args(command, config)
    };
    let config = |command: &str| config_to(command, None);
    // A depth limit equal to the root's depth evaluates the root and nothing
    // else. (A shallower one is refused: see the invalid invocations.)
    let root_only = succeed(&config_to("search", Some(4)), 4);
    assert_eq!(root_only["evaluated"], 1);
    let recorded = root_only["records"].as_u64().unwrap();
    assert_eq!(recorded + root_only["unresolved"].as_u64().unwrap(), 1);
    assert_eq!(root_only["halt"], if recorded == 1 { "complete" } else { "depth-limit" });
    let summary = succeed(&config("search"), 4);
    assert_eq!((&summary["halt"], &summary["root"]), (&json!("complete"), &json!("0110")));
    let checked = succeed(&config("check"), 4);
    assert_eq!(checked["complete"], true);
    // The same file does not pass as a certificate for the whole root.
    assert_eq!(run(&dir.check(64), &AtomicBool::new(false), StandIn::new(4)).0, FAILURE);
}

/// A check of a valid but incomplete certificate exits with status 0, like a
/// complete one: only the summary's `complete` field tells them apart.
#[test]
fn an_incomplete_check_succeeds_and_says_so() {
    let _files = files();
    let dir = Folder::new();
    succeed(&dir.search(1, 1, 64, Some(3)), 9);
    let (status, summary) = run(&dir.check(64), &AtomicBool::new(false), StandIn::new(9));
    assert_eq!(status, SUCCESS);
    let summary = summary.unwrap();
    assert_eq!((&summary["complete"], &summary["header_complete"]), (&json!(false), &json!(true)));
    assert!(summary["frontier"].as_u64().unwrap() > 0);
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Depth {
    depth: usize,
}

impl RecordData for Depth {}

/// Accepts the boxes at depth 3 and deeper, slowly.
struct Slow;

impl Components for Slow {
    type Data = Depth;
    fn report(&self) -> Value {
        json!({})
    }
    fn check(&self, path: &str) -> Result<Option<Depth>, BoxError> {
        thread::sleep(Duration::from_millis(20));
        Ok((path.len() >= 3).then(|| Depth { depth: path.len() }))
    }
    fn holds(&self, path: &str, data: &Depth) -> Result<(), BoxError> {
        if path.len() >= 3 && data.depth == path.len() {
            Ok(())
        } else {
            Err("false".into())
        }
    }
}

/// Two invocations on one certificate: the second is refused (locked) and
/// changes nothing; the first then finishes normally.
#[test]
fn a_second_invocation_on_a_running_certificate_is_refused() {
    let _files = files();
    let dir = Folder::new();
    let stop = AtomicBool::new(false);
    let first = search_with(&dir, "", 16, None);
    let file = dir.certificate();
    thread::scope(|scope| {
        let search = scope.spawn(|| {
            let mut out = Vec::new();
            let status = entry(&first, &stop, |_, _| Ok(Slow), &mut out);
            (status, out)
        });
        let started = Instant::now();
        while fs::metadata(&file).map(|m| m.len()).unwrap_or(0) == 0 {
            assert!(started.elapsed() < Duration::from_secs(10));
            thread::sleep(Duration::from_millis(2));
        }
        let check = dir.args(
            "check",
            json!({ "certificate": dir.certificate(), "root": "", "max_depth": 16, "threads": 1 }),
        );
        for args in [search_with(&dir, "", 16, None), check] {
            let (status, summary) = run_with(&args, &stop, Slow);
            assert_eq!((status, summary), (FAILURE, None));
        }
        let (status, out) = search.join().unwrap();
        assert_eq!(status, SUCCESS);
        let summary: Value = serde_json::from_slice(&out).unwrap();
        assert_eq!(summary["halt"], "complete");
    });
    let check = dir.args(
        "check",
        json!({ "certificate": dir.certificate(), "root": "", "max_depth": 16, "threads": 1 }),
    );
    let (status, summary) = run_with(&check, &stop, Slow);
    assert_eq!(status, SUCCESS);
    let summary = summary.unwrap();
    assert_eq!((summary["complete"].clone(), summary["records"].clone()), (json!(true), json!(8)));
}
