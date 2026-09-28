//! End-to-end tests of the commands with the real components: searches and
//! checks of small subtrees by Domain and Global, and by Exotic and
//! Local over small proved stand-in covers.

mod support;

use rid::components::collection::Collection;
use rid::components::record::{ComponentName, ExoticCover, RecordData};
use rid::elimination::proof::{catalogue, ProvedSet};
use rid::problem::configuration::ConfigurationBox;
use rid::search::command::{entry, FAILURE, SUCCESS};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::ffi::OsString;
use std::fs;
use std::process::{Command as Process, Stdio};
use std::sync::atomic::AtomicBool;
use support::{files, processes, snapshot, Folder};

/// Subtree roots whose search by Domain and Global alone finishes after a few
/// dozen boxes, with records of both.
const ROOTS: [&str; 3] = ["01110101000", "0100110011011", "10010001110"];
const MAX_DEPTH: usize = 64;

fn from_path(path: &str) -> ConfigurationBox {
    ConfigurationBox::from_path(path).unwrap()
}

/// A box cover whose proof is the completed search of its box (asserted by
/// [`proved_covers`]); for these tests only.
struct Proved {
    name: &'static str,
    cover: ConfigurationBox,
}

impl ProvedSet for Proved {
    fn name(&self) -> &str {
        self.name
    }
    fn contains(&self, b: &ConfigurationBox) -> bool {
        self.cover.contains(b)
    }
}

fn collection(covers: &[(&'static str, &str)]) -> Collection {
    let covers = covers
        .iter()
        .map(|&(name, path)| Box::new(Proved { name, cover: from_path(path) }) as Box<dyn ProvedSet>)
        .collect();
    Collection::from_covers(covers).unwrap()
}

/// The records of the search of `root` by `collection`, box by box in search
/// order, independently of the queue and the certificate.
fn sequential(collection: &Collection, root: &str) -> Vec<(String, RecordData)> {
    let mut records = Vec::new();
    let mut layer = vec![root.to_owned()];
    while !layer.is_empty() {
        let mut next = Vec::new();
        for path in layer {
            assert!(path.len() <= MAX_DEPTH);
            match collection.check(&from_path(&path)) {
                Some(data) => records.push((path, data)),
                None => next.extend([format!("{path}0"), format!("{path}1")]),
            }
        }
        layer = next;
    }
    records
}

fn search_args(dir: &Folder, root: &str, threads: usize, window: usize, max_decisions: Option<u64>) -> Vec<OsString> {
    dir.args(
        "search",
        json!({
            "certificate": dir.certificate(),
            "root": root,
            "threads": threads,
            "window": window,
            "max_depth": MAX_DEPTH,
            "depth_limit": null,
            "max_decisions": max_decisions,
        }),
    )
}

fn check_args(dir: &Folder, root: &str, threads: usize) -> Vec<OsString> {
    dir.args(
        "check",
        json!({ "certificate": dir.certificate(), "root": root, "max_depth": MAX_DEPTH, "threads": threads }),
    )
}

fn run(args: &[OsString], collection: Collection) -> (u8, Option<Value>) {
    let mut out = Vec::new();
    let status = entry(args, &AtomicBool::new(false), |_, _| Ok(collection), &mut out);
    (status, (!out.is_empty()).then(|| serde_json::from_slice(&out).unwrap()))
}

fn succeed(args: &[OsString], collection: Collection) -> Value {
    let (status, summary) = run(args, collection);
    assert_eq!(status, SUCCESS);
    summary.unwrap()
}

/// The records of a certificate file: `(path, data)` from each line after the header.
fn records(bytes: &[u8]) -> Vec<(String, RecordData)> {
    let text = std::str::from_utf8(bytes).unwrap();
    text.lines()
        .skip(1)
        .map(|line| {
            let mut value: Value = serde_json::from_str(&line[9..]).unwrap();
            let path = value["path"].as_str().unwrap().to_owned();
            value.as_object_mut().unwrap().remove("path");
            (path, serde_json::from_value(value).unwrap())
        })
        .collect()
}

/// A line with a correct checksum for `payload`.
fn line(payload: &str) -> String {
    let hash = Sha256::digest(payload.as_bytes());
    format!("{} {payload}\n", hash[..4].iter().map(|b| format!("{b:02x}")).collect::<String>())
}

#[test]
fn every_thread_count_writes_the_certificate_of_the_sequential_search() {
    let _files = files();
    for root in ROOTS {
        let expected = sequential(&collection(&[]), root);
        let kinds: Vec<ComponentName> = expected.iter().map(|(_, d)| d.component()).collect();
        assert!(kinds.contains(&ComponentName::Domain) && kinds.contains(&ComponentName::Global));
        let mut first: Option<Vec<u8>> = None;
        for (threads, window) in [(1, 1), (2, 3), (3, 16), (4, 64)] {
            let dir = Folder::new();
            let summary = succeed(&search_args(&dir, root, threads, window, None), collection(&[]));
            assert_eq!((&summary["halt"], &summary["unresolved"]), (&json!("complete"), &json!(0)));
            let bytes = fs::read(dir.certificate()).unwrap();
            assert_eq!(records(&bytes), expected, "{root} on {threads} threads");
            match &first {
                None => first = Some(bytes.clone()),
                Some(bytes_before) => assert_eq!(&bytes, bytes_before),
            }
            let checked = succeed(&check_args(&dir, root, threads), collection(&[]));
            assert_eq!((&checked["complete"], &checked["records"]), (&json!(true), &json!(expected.len())));
        }
    }
}

#[test]
fn proved_covers_give_exotic_and_local_records() {
    let _files = files();
    let root = ROOTS[2];
    let bare = sequential(&collection(&[]), root);
    // Two boxes that the bare search splits and whose subtrees it completes:
    // covers proved by that search.
    let split: Vec<String> = {
        let mut parents: Vec<String> = bare.iter().map(|(p, _)| p[..p.len() - 1].to_owned()).collect();
        parents.sort_by_key(|p| (p.len(), p.clone()));
        parents.dedup();
        parents.into_iter().filter(|p| p.len() > root.len()).collect()
    };
    let (square, local) = (split[0].clone(), split[split.len() - 1].clone());
    assert!(!square.starts_with(&local) && !local.starts_with(&square));
    for cover in [&square, &local] {
        let below: Vec<_> = sequential(&collection(&[]), cover);
        assert!(below.len() >= 2 && bare.iter().filter(|(p, _)| p.starts_with(cover.as_str())).count() == below.len());
    }
    let covers = [("square", square.as_str()), ("4", local.as_str())];
    let expected = sequential(&collection(&covers), root);
    let find = |name| expected.iter().find(|(_, d)| d.component() == name).map(|(p, _)| p.clone());
    assert_eq!(find(ComponentName::Exotic), Some(square.clone()));
    assert_eq!(find(ComponentName::Local), Some(local.clone()));
    let dir = Folder::new();
    let mut certificate: Option<Vec<u8>> = None;
    for threads in [1, 3] {
        let _ = fs::remove_file(dir.certificate());
        succeed(&search_args(&dir, root, threads, 5, None), collection(&covers));
        let bytes = fs::read(dir.certificate()).unwrap();
        assert_eq!(records(&bytes), expected);
        assert!(certificate.as_ref().map_or(true, |before| *before == bytes));
        certificate = Some(bytes);
    }
    let checked = succeed(&check_args(&dir, root, 2), collection(&covers));
    assert_eq!(checked["complete"], true);
    // Without the covers, or with them swapped, the same file is refused.
    let before = snapshot(&dir.certificate());
    for other in [vec![], vec![("square", local.as_str()), ("4", square.as_str())], vec![("arc+", square.as_str()), ("4", local.as_str())]] {
        assert_eq!(run(&check_args(&dir, root, 2), collection(&other)), (FAILURE, None));
        assert_eq!(run(&search_args(&dir, root, 2, 5, None), collection(&other)), (FAILURE, None));
    }
    assert_eq!(snapshot(&dir.certificate()), before);
}

#[test]
fn false_records_with_correct_checksums_are_refused() {
    let _files = files();
    let root = ROOTS[0];
    let dir = Folder::new();
    succeed(&search_args(&dir, root, 2, 4, None), collection(&[]));
    let original = fs::read_to_string(dir.certificate()).unwrap();
    let lines: Vec<&str> = original.lines().collect();
    let checker = collection(&[]);
    let mut refused = 0;
    for (k, (path, data)) in records(original.as_bytes()).into_iter().enumerate() {
        let b = from_path(&path);
        // A changed claim for the same box, false by the component's own rule.
        let mut value = serde_json::to_value(&data).unwrap();
        match data.component() {
            ComponentName::Domain => {
                let i = value["inequality"].as_u64().unwrap();
                let Some(other) = (0..73).map(|j| (i + j) % 73).find(|&j| {
                    let changed: RecordData = serde_json::from_value(json!({"component": "Domain", "inequality": j})).unwrap();
                    checker.holds(&b, &changed).is_err()
                }) else { continue };
                value["inequality"] = json!(other);
            }
            ComponentName::Global => {
                // The reversed edge points to the opposite side of the hole.
                let edge = value["edge"].as_array().unwrap().clone();
                value["edge"] = json!([edge[1], edge[0]]);
            }
            other => panic!("unexpected {other}"),
        }
        let changed: RecordData = serde_json::from_value(value).unwrap();
        assert!(checker.holds(&b, &changed).is_err());
        let mut payload = serde_json::to_value(&changed).unwrap();
        let object = payload.as_object_mut().unwrap();
        let mut ordered = serde_json::Map::new();
        ordered.insert("path".into(), json!(path));
        ordered.append(object);
        let payload = serde_json::to_string(&Value::Object(ordered)).unwrap();
        let mut text = String::new();
        for (j, l) in lines.iter().enumerate() {
            if j == k + 1 {
                text.push_str(&line(&payload));
            } else {
                text.push_str(l);
                text.push('\n');
            }
        }
        fs::write(dir.certificate(), &text).unwrap();
        assert_eq!(run(&check_args(&dir, root, 2), collection(&[])), (FAILURE, None), "{payload}");
        assert_eq!(fs::read_to_string(dir.certificate()).unwrap(), text);
        refused += 1;
    }
    assert!(refused >= 10, "{refused}");
    // The original, rewritten line by line with fresh checksums, still checks.
    let rebuilt: String = lines.iter().map(|l| line(&l[9..])).collect();
    assert_eq!(rebuilt, original);
}

#[test]
fn bounded_searches_resume_to_the_same_certificate() {
    let _files = files();
    let root = ROOTS[1];
    let reference = Folder::new();
    succeed(&search_args(&reference, root, 1, 1, None), collection(&[]));
    let expected = fs::read(reference.certificate()).unwrap();
    let dir = Folder::new();
    let mut runs = 0;
    for attempt in 1.. {
        assert!(attempt < 100, "no progress");
        runs += 1;
        // Refusals after the last saved record are evaluated again by the
        // next run, so the budget grows until a record follows them.
        let budget = Some(5 * attempt as u64);
        let summary = succeed(&search_args(&dir, root, 1 + attempt % 4, 1 + attempt % 7, budget), collection(&[]));
        let bytes = fs::read(dir.certificate()).unwrap();
        assert!(expected.starts_with(&bytes));
        if summary["halt"] == "complete" {
            break;
        }
    }
    assert_eq!(fs::read(dir.certificate()).unwrap(), expected);
    assert!(runs > 2);
}

/// The production binary runs with the cover catalogue's components: when
/// the catalogue holds the eight Exotic covers, it prepares them and
/// searches and checks a small subtree to completion; otherwise every
/// command refuses, naming the missing covers, before any file exists.
#[test]
fn the_production_binary_uses_the_catalogues_covers() {
    let _processes = processes();
    let dir = Folder::new();
    let rid = env!("CARGO_BIN_EXE_rid");
    let run = |args: &[OsString]| Process::new(rid).args(args).stdin(Stdio::null()).output().unwrap();
    let complete = ExoticCover::ALL
        .iter()
        .all(|cover| catalogue::exotic_data(cover.name()).is_some());
    let prepare = dir.args("prepare", json!({ "threads": 4 }));
    let (search, check) = (search_args(&dir, ROOTS[0], 4, 16, None), check_args(&dir, ROOTS[0], 4));
    if !complete {
        for args in [&prepare, &search, &check] {
            let output = run(args);
            assert_eq!(output.status.code(), Some(i32::from(FAILURE)), "{args:?}");
            assert!(output.stdout.is_empty());
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains("Exotic covers") && stderr.contains("missing"), "{stderr}");
        }
        assert!(!dir.certificate().exists());
        return;
    }
    let output = run(&prepare);
    assert_eq!(output.status.code(), Some(i32::from(SUCCESS)));
    let summary: Value = serde_json::from_slice(&output.stdout).unwrap();
    let names: Vec<&str> = ExoticCover::ALL.iter().map(|r| r.name()).collect();
    assert_eq!(summary["components"]["exotic_covers"], json!(names));
    assert!(!summary["components"]["local_covers"].as_array().unwrap().is_empty());
    for (args, field) in [(&search, "halt"), (&check, "complete")] {
        let output = run(args);
        assert_eq!(output.status.code(), Some(i32::from(SUCCESS)), "{}", String::from_utf8_lossy(&output.stderr));
        let summary: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert!(summary[field] == json!("complete") || summary[field] == json!(true), "{summary}");
    }
}

/// Boxes per second of the real search command with Domain and Global: eight
/// subtrees at depths 9 to 12, each searched 20 levels deep on 4 threads
/// (evaluation, the reader's verification of every record and the
/// certificate's appends).
///
/// cargo test --release --test components measure_search -- --ignored --nocapture
#[test]
#[ignore]
fn measure_search() {
    let _files = files();
    let mut state = 12345u64;
    let mut bit = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        if state >> 63 == 0 { '0' } else { '1' }
    };
    let (mut evaluated, mut records, mut milliseconds) = (0u64, 0u64, 0u64);
    for k in 0..8 {
        let root: String = (0..9 + k % 4).map(|_| bit()).collect();
        let dir = Folder::new();
        let args = dir.args(
            "search",
            json!({
                "certificate": dir.certificate(),
                "root": root,
                "threads": 4,
                "window": 64,
                "max_depth": MAX_DEPTH,
                "depth_limit": root.len() + 20,
                "max_decisions": null,
            }),
        );
        let summary = succeed(&args, collection(&[]));
        evaluated += summary["evaluated"].as_u64().unwrap();
        records += summary["new_records"].as_u64().unwrap();
        milliseconds += summary["milliseconds"].as_u64().unwrap();
        println!("{root}: {}", summary);
    }
    println!(
        "{evaluated} boxes, {records} records in {milliseconds} ms: {} boxes per second on 4 threads",
        1000 * evaluated / milliseconds.max(1)
    );
}

/// The time of each component's check with the catalogue's covers, on real
/// boxes: evenly spaced records of each component in the committed full
/// certificate (`results/full/search.cert`) and their parents (boxes every
/// component refused). For Global, also the time of the same check without
/// the floating-point proposal (whether any of all 900 fixed candidates
/// holds), on fewer boxes, since that is what the proposal saves. Prints
/// a note when the catalogue or the certificate is missing.
///
/// MEASURE_BOXES=48 cargo test --release --test components measure_catalogue_components -- --ignored --nocapture
#[test]
#[ignore]
fn measure_catalogue_components() {
    use rid::components::Component;
    use rid::search::command::wiring;
    use std::num::NonZeroUsize;
    use std::time::Instant;
    let certificate = concat!(env!("CARGO_MANIFEST_DIR"), "/results/full/search.cert");
    let Ok(text) = fs::read_to_string(certificate) else {
        return println!("no certificate at {certificate}");
    };
    let stop = AtomicBool::new(false);
    let started = Instant::now();
    let collection = match wiring::prepare(NonZeroUsize::new(4).unwrap(), &stop) {
        Ok(collection) => collection,
        Err(error) => return println!("no catalogue: {error}"),
    };
    println!("covers loaded and verified in {} ms on 4 threads", started.elapsed().as_millis());
    let count: usize = std::env::var("MEASURE_BOXES").map_or(24, |n| n.parse().unwrap());
    let recorded: Vec<(String, String)> = text
        .lines()
        .skip(1)
        .map(|line| {
            let value: Value = serde_json::from_str(&line[9..]).unwrap();
            (value["path"].as_str().unwrap().to_owned(), value["component"].as_str().unwrap().to_owned())
        })
        .collect();
    let spread = |paths: Vec<String>| -> Vec<String> {
        let step = (paths.len() / count).max(1);
        paths.into_iter().step_by(step).take(count).collect()
    };
    let mut sets: Vec<(String, Vec<String>)> = ["Domain", "Exotic", "Local", "Global"]
        .iter()
        .map(|name| {
            let paths = recorded.iter().filter(|(_, c)| c == name).map(|(p, _)| p.clone()).collect();
            (format!("{name} records"), spread(paths))
        })
        .collect();
    let mut parents: Vec<String> = recorded.iter().map(|(p, _)| p[..p.len() - 1].to_owned()).collect();
    parents.sort_by(|a, b| (a.len(), a).cmp(&(b.len(), b)));
    parents.dedup();
    sets.push(("parents".into(), spread(parents)));
    let all = collection.global.all();
    for (label, paths) in sets {
        let mut micros = [0u128; 5];
        let mut outcome = [0usize; 5];
        let mut unproposed = (0u128, 0usize, 0usize);
        for (n, path) in paths.iter().enumerate() {
            let b = from_path(path);
            let timed = |f: &dyn Fn() -> bool| {
                let started = Instant::now();
                let held = f();
                (started.elapsed().as_micros(), held)
            };
            let results = [
                timed(&|| collection.domain.check(&b).is_some()),
                timed(&|| collection.exotic.check(&b).is_some()),
                timed(&|| collection.local.check(&b).is_some()),
                timed(&|| collection.global.check(&b).is_some()),
                timed(&|| collection.check(&b).is_some()),
            ];
            for (k, (time, held)) in results.into_iter().enumerate() {
                micros[k] += time;
                outcome[k] += usize::from(held);
            }
            if n < count / 4 {
                let (time, held) = timed(&|| all.iter().any(|m| collection.global.holds(&b, m).is_ok()));
                unproposed = (unproposed.0 + time, unproposed.1 + 1, unproposed.2 + usize::from(held));
            }
        }
        let boxes = paths.len().max(1) as u128;
        println!(
            "{label} ({} boxes): Domain {} µs ({} held), Exotic {} µs ({}), Local {} µs ({}), Global {} µs ({}); \
             the collection {} µs per box; Global without the proposal {} µs ({} of {} held)",
            paths.len(),
            micros[0] / boxes,
            outcome[0],
            micros[1] / boxes,
            outcome[1],
            micros[2] / boxes,
            outcome[2],
            micros[3] / boxes,
            outcome[3],
            micros[4] / boxes,
            unproposed.0 / (unproposed.1.max(1) as u128),
            unproposed.2,
            unproposed.1,
        );
    }
}
