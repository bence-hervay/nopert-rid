//! Tests of the generate command: the typed configuration (every field
//! required, unknown, repeated and mistyped fields refused), the expansion
//! into pinned requests, usage and read errors, and an end-to-end run that
//! writes verified files, reports failures and leaves no partial file.
use super::*;
use crate::elimination::proof::format::CoverFile;
use serde_json::Value;
use std::sync::atomic::AtomicBool;

fn example() -> Value {
    json!({
        "output": "data",
        "threads": 4,
        "search": {"max_depth": 60, "max_nodes": 200000, "aspect": 16, "tried": 6,
                   "proposal_depth": 4, "proposals": 8, "inherited": 16},
        "exotic": [{"name": "square", "samples": 1, "inequalities": [29]},
                        {"name": "crossing+", "samples": 8, "inequalities": [10, 11]}],
        "local": {"numbers": [0, 29], "samples": 2}
    })
}

fn parse(value: &Value) -> Result<Configuration, Error> {
    Configuration::parse(value.to_string().as_bytes())
}

/// A scratch directory of this test process.
fn scratch(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("rid-generate-test-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn the_documented_configuration_parses_and_expands() {
    let configuration = parse(&example()).unwrap();
    assert_eq!(configuration.threads.get(), 4);
    assert_eq!(configuration.search.aspect.get(), 16);
    let requests = configuration.requests().unwrap();
    let files: Vec<PathBuf> = requests.iter().map(|(_, f)| f.clone()).collect();
    assert_eq!(
        files,
        ["exotic/square.json", "exotic/crossing+.json", "local/0.json", "local/29.json"]
            .map(PathBuf::from)
            .to_vec()
    );
    assert_eq!(requests[0].0.pin, catalogue::exotic_pin("square").unwrap());
    assert_eq!(requests[1].0.inequalities, vec![10, 11]);
    assert_eq!(requests[3].0.pin, catalogue::local_pin(29).unwrap());
    assert!(requests[2].0.inequalities.is_empty() && requests[2].0.samples.get() == 2);
    // Round trip.
    let text = serde_json::to_string(&configuration).unwrap();
    assert_eq!(Configuration::parse(text.as_bytes()).unwrap(), configuration);
}

/// Every path to a field of an object in `value`.
fn fields(value: &Value, prefix: Vec<String>, out: &mut Vec<Vec<String>>) {
    if let Value::Object(map) = value {
        for (key, child) in map {
            let mut path = prefix.clone();
            path.push(key.clone());
            out.push(path.clone());
            fields(child, path, out);
        }
    }
    if let Value::Array(items) = value {
        if let Some(first) = items.first() {
            let mut path = prefix.clone();
            path.push("0".into());
            fields(first, path, out);
        }
    }
}

fn at<'v>(value: &'v mut Value, path: &[String]) -> &'v mut Value {
    path.iter().fold(value, |v, key| match v {
        Value::Array(items) => &mut items[key.parse::<usize>().unwrap()],
        _ => &mut v[key.as_str()],
    })
}

#[test]
fn every_missing_unknown_or_mistyped_field_is_refused() {
    let mut paths = Vec::new();
    fields(&example(), Vec::new(), &mut paths);
    assert_eq!(paths.len(), 17);
    for path in &paths {
        let (last, parent) = path.split_last().unwrap();
        let mut missing = example();
        at(&mut missing, parent).as_object_mut().unwrap().remove(last);
        assert!(matches!(parse(&missing), Err(Error::Json(_))), "missing {path:?} accepted");
        let mut extra = example();
        at(&mut extra, parent).as_object_mut().unwrap().insert("extra".into(), json!(1));
        assert!(matches!(parse(&extra), Err(Error::Json(_))), "unknown field beside {path:?} accepted");
        let mut null = example();
        *at(&mut null, path) = Value::Null;
        assert!(matches!(parse(&null), Err(Error::Json(_))), "null {path:?} accepted");
        let mut text = example();
        *at(&mut text, path) = json!("text");
        if path != &["output"] && path.last().map(String::as_str) != Some("name") {
            assert!(matches!(parse(&text), Err(Error::Json(_))), "string {path:?} accepted");
        }
    }
    for (path, value) in [
        (vec!["threads"], json!(0)),
        (vec!["threads"], json!(-1)),
        (vec!["threads"], json!(2.0)),
        (vec!["search", "aspect"], json!(1.5)),
        (vec!["search", "max_depth"], json!(-3)),
        (vec!["local", "samples"], json!(-1)),
        (vec!["local", "samples"], json!(0)),
        (vec!["exotic", "0", "samples"], json!(0)),
        (vec!["exotic", "0", "inequalities"], json!([-1])),
    ] {
        let mut bad = example();
        *at(&mut bad, &path.iter().map(|s| s.to_string()).collect::<Vec<_>>()) = value.clone();
        assert!(matches!(parse(&bad), Err(Error::Json(_))), "{path:?} = {value} accepted");
    }
    let repeated = r#"{"output":"a","output":"b","threads":1,"search":{"max_depth":1,"max_nodes":1,"aspect":1,"tried":1,"proposal_depth":1,"proposals":1,"inherited":1},"exotic":[],"local":{"numbers":[],"samples":1}}"#;
    assert!(matches!(Configuration::parse(repeated.as_bytes()), Err(Error::Json(_))));
    assert!(matches!(Configuration::parse(b"{"), Err(Error::Json(_))));
    assert!(matches!(Configuration::parse(b"[]"), Err(Error::Json(_))));
}

#[test]
fn unknown_and_repeated_covers_are_refused() {
    let mut unknown = example();
    unknown["exotic"][0]["name"] = json!("wall");
    assert!(matches!(parse(&unknown).unwrap().requests(), Err(Error::UnknownCover(n)) if n == "wall"));
    let mut far = example();
    far["local"]["numbers"] = json!([32]);
    assert!(matches!(parse(&far).unwrap().requests(), Err(Error::UnknownLocal(32))));
    let mut twice = example();
    twice["exotic"][1]["name"] = json!("square");
    assert!(matches!(parse(&twice).unwrap().requests(), Err(Error::Repeated(n)) if n == "square"));
    let mut twice = example();
    twice["local"]["numbers"] = json!([3, 3]);
    assert!(matches!(parse(&twice).unwrap().requests(), Err(Error::Repeated(n)) if n == "3"));
}

#[test]
fn usage_and_read_errors_have_their_exit_status() {
    let stop = AtomicBool::new(false);
    let mut out = Vec::new();
    assert_eq!(entry(&[], &stop, &mut out), USAGE_ERROR);
    assert_eq!(entry(&["a".into(), "b".into()], &stop, &mut out), USAGE_ERROR);
    let missing = scratch("missing").join("absent.json");
    assert_eq!(entry(&[missing.clone().into()], &stop, &mut out), FAILURE);
    assert!(matches!(run(&missing, &stop, &mut out), Err(Error::Read { .. })));
    assert!(out.is_empty());
}

#[test]
fn a_run_writes_verified_files_and_reports_failures() {
    let directory = scratch("run");
    let output = directory.join("data");
    let mut value = example();
    value["output"] = json!(output);
    value["threads"] = json!(2);
    value["exotic"] = json!([]);
    value["local"] = json!({"numbers": [13, 19], "samples": 1});
    let config = directory.join("generate.json");
    std::fs::write(&config, value.to_string()).unwrap();
    let stop = AtomicBool::new(false);
    let mut out = Vec::new();
    assert_eq!(entry(&[config.clone().into()], &stop, &mut out), SUCCESS);
    let summary: Value = serde_json::from_slice(&out).unwrap();
    assert_eq!(summary["command"], "generate");
    assert_eq!(summary["covers"][0]["cover"], "13");
    assert_eq!(summary["covers"][1]["written"], true);
    let first = std::fs::read(output.join("local/13.json")).unwrap();
    let file = CoverFile::parse(&first).unwrap();
    let pin = catalogue::local_pin(13).unwrap();
    let cover = ProvedCover::load(&pin, &first, NonZeroUsize::new(2).unwrap(), &stop).unwrap();
    assert_eq!(serde_json::to_value(cover.report()).unwrap()["witness_leaves"], summary["covers"][0]["statistics"]["witness_leaves"]);
    assert_eq!(file.name, "13");
    // A second run writes the same bytes.
    let mut again = Vec::new();
    assert_eq!(entry(&[config.clone().into()], &stop, &mut again), SUCCESS);
    assert_eq!(std::fs::read(output.join("local/13.json")).unwrap(), first);
    // A cover that cannot be generated is reported, not written, and fails
    // the command; the other cover is still written.
    std::fs::remove_dir_all(&output).unwrap();
    value["search"]["max_depth"] = json!(3);
    value["local"]["numbers"] = json!([14, 13]);
    std::fs::write(&config, value.to_string()).unwrap();
    let mut failed = Vec::new();
    assert_eq!(entry(&[config.into()], &stop, &mut failed), FAILURE);
    let summary: Value = serde_json::from_slice(&failed).unwrap();
    assert_eq!(summary["covers"][0]["written"], false);
    assert!(summary["covers"][0]["error"].as_str().unwrap().contains("depth limit"));
    assert!(!output.join("local/14.json").exists());
    assert_eq!(summary["covers"][1]["written"], true);
    assert_eq!(std::fs::read(output.join("local/13.json")).unwrap(), first);
    let leftovers: Vec<_> = walk(&directory).into_iter().filter(|p| p.to_string_lossy().ends_with(".partial")).collect();
    assert!(leftovers.is_empty());
    std::fs::remove_dir_all(&directory).unwrap();
}

fn walk(directory: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for entry in std::fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            out.extend(walk(&path));
        } else {
            out.push(path);
        }
    }
    out
}

#[test]
fn counts_must_agree() {
    let mut statistics = super::super::Statistics { zooms: 1, witnesses: 2, witness_leaves: 3, ..Default::default() };
    let report = Report { zooms: 1, witnesses: 2, witness_leaves: 3, domain_leaves: 0, delegated_leaves: 0 };
    assert!(agrees(&report, &statistics));
    statistics.delegated_leaves = 1;
    assert!(!agrees(&report, &statistics));
}

#[test]
fn the_shipped_configuration_generates_every_cover_once_into_data() {
    let configuration = Configuration::parse(include_bytes!("covers.json")).unwrap();
    assert_eq!(configuration.output, PathBuf::from("data"));
    let requests = configuration.requests().unwrap();
    let names: Vec<String> = requests.iter().map(|(r, _)| r.pin.name.clone()).collect();
    let mut expected: Vec<String> = catalogue::EXOTIC.iter().map(|n| n.to_string()).collect();
    expected.extend((0..LOCAL_COUNT).map(|n| n.to_string()));
    assert_eq!(names, expected);
    for (request, file) in &requests {
        let data = match file.parent().unwrap().to_str().unwrap() {
            "exotic" => catalogue::exotic_data(&request.pin.name),
            _ => catalogue::local_data(request.pin.name.parse().unwrap()),
        };
        // The catalogue's data file of the same name has the same pin.
        let parsed = CoverFile::parse(data.unwrap()).unwrap();
        assert_eq!(parsed.name, request.pin.name);
        assert!(request.inequalities.iter().all(|i| *i < crate::problem::domain::CONSTRAINT_COUNT));
    }
}

/// The safety step in isolation: a generated cover whose file the loader
/// refuses, or whose counts disagree, is never written (no file, no partial
/// file), and an interrupted verification is reported as an interruption.
#[test]
fn verify_and_write_writes_only_what_the_loader_accepts() {
    use crate::elimination::proof::format::Leaf;
    let directory = scratch("verify");
    let threads = NonZeroUsize::new(2).unwrap();
    let stop = AtomicBool::new(false);
    let request = Request { pin: catalogue::local_pin(19).unwrap(), samples: NonZeroUsize::new(1).unwrap(), inequalities: Vec::new() };
    let settings = Configuration::parse(include_bytes!("covers.json")).unwrap().search;
    let generated = super::super::generate(std::slice::from_ref(&request), &settings, threads, &stop, &mut |_| {})
        .unwrap()
        .remove(0)
        .unwrap();
    let file = directory.join("local/19.json");
    let nothing_written = |file: &Path| {
        assert!(!file.exists());
        assert!(walk(&directory).iter().all(|p| !p.to_string_lossy().ends_with(".partial")));
    };
    // A leaf's witness moved to another witness, a leaf dropped, a count off.
    let mut moved = generated.clone();
    let leaves = &mut moved.file.zooms[0].leaves;
    let at = leaves.iter().position(|l| matches!(l, Leaf::Witness { .. })).unwrap();
    if let Leaf::Witness { index, .. } = &mut leaves[at] {
        *index = (*index + 1) % moved.file.witnesses.len();
    }
    let mut dropped = generated.clone();
    dropped.file.zooms[0].leaves.pop();
    let mut miscounted = generated.clone();
    miscounted.statistics.witness_leaves += 1;
    for (corrupted, counts) in [(moved, false), (dropped, false), (miscounted, true)] {
        let refusal = verify_and_write(&request, &corrupted, &file, threads, &stop).unwrap_err();
        if counts {
            assert!(matches!(refusal, Error::Counts { .. }), "{refusal}");
        } else {
            assert!(matches!(refusal, Error::Verify { .. }), "{refusal}");
        }
        nothing_written(&file);
    }
    // Interrupted: a typed interruption, nothing written.
    let stopped = AtomicBool::new(true);
    assert!(matches!(
        verify_and_write(&request, &generated, &file, threads, &stopped),
        Err(Error::Generate(GenerateError::Interrupted))
    ));
    nothing_written(&file);
    // The genuine file is written.
    verify_and_write(&request, &generated, &file, threads, &stop).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), generated.file.to_bytes());
    std::fs::remove_dir_all(&directory).unwrap();
}

/// Refused before any work: domain inequalities that do not exist or are
/// offered to a cover of scope `all`, a sample grid above the limit and a
/// `max_depth` above the limit.
#[test]
fn requests_are_checked_before_any_work() {
    let mut unknown = example();
    unknown["exotic"][0]["inequalities"] = json!([73]);
    assert!(matches!(parse(&unknown).unwrap().requests(), Err(Error::Generate(GenerateError::Inequality { .. }))));
    let mut samples = example();
    samples["local"]["samples"] = json!(65);
    assert!(matches!(parse(&samples).unwrap().requests(), Err(Error::Generate(GenerateError::Samples { samples: 65, .. }))));
    let mut depth = example();
    depth["search"]["max_depth"] = json!(201);
    assert!(matches!(parse(&depth).unwrap().requests(), Err(Error::Generate(GenerateError::MaxDepth { max_depth: 201 }))));
    depth["search"]["max_depth"] = json!(200);
    assert!(parse(&depth).unwrap().requests().is_ok());
    let mut aspect = example();
    aspect["search"]["aspect"] = json!(0);
    assert!(matches!(parse(&aspect), Err(Error::Json(_))));
}
