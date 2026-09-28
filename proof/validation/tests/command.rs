//! The command line: exit statuses, and the experiments run end to end with
//! the crate's components (its cover catalogue loaded and verified).
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_rid-validation")
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("command-{name}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn status(args: &[&str]) -> (i32, String) {
    let output = Command::new(binary()).args(args).output().unwrap();
    (
        output.status.code().unwrap(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

#[test]
fn malformed_command_lines_exit_with_status_2() {
    assert_eq!(status(&[]).0, 2);
    assert_eq!(status(&["check"]).0, 2);
    assert_eq!(status(&["check", "a.json", "b.json"]).0, 2);
    let (code, message) = status(&["search", "a.json"]);
    assert_eq!(code, 2);
    assert!(message.contains("unknown command"));
}

#[test]
fn unreadable_and_invalid_configurations_exit_with_status_1() {
    let dir = scratch("invalid");
    assert_eq!(status(&["check", dir.join("missing.json").to_str().unwrap()]).0, 1);
    let config = dir.join("check.json");
    std::fs::write(&config, r#"{"catalogue": "c", "transcript": "t"}"#).unwrap();
    let (code, message) = status(&["check", config.to_str().unwrap()]);
    assert_eq!(code, 1);
    assert!(message.contains("threads"), "{message}");
    // A valid configuration naming a missing transcript.
    std::fs::write(&config, r#"{"catalogue": "c", "transcript": "t", "threads": 1}"#).unwrap();
    assert_eq!(status(&["check", config.to_str().unwrap()]).0, 1);
}

/// Completeness on three points (one each for Local, Domain and Global)
/// with the crate's components, its check and a comparison with itself, and
/// a small pilot search.
#[test]
fn experiments_run_with_the_crate_components() {
    let dir = scratch("components");
    let catalogue = Path::new(env!("CARGO_MANIFEST_DIR")).join("catalogue").join("points.json");
    let transcript = dir.join("completeness.jsonl");
    let write = |name: &str, text: String| {
        let file = dir.join(name);
        std::fs::write(&file, text).unwrap();
        file
    };
    let completeness = write(
        "completeness.json",
        format!(
            r#"{{"catalogue": {catalogue:?}, "transcript": {transcript:?}, "threads": 4, "max_k": 12,
                "selection": ["aligned-generic-rational", "outside-root-corner", "global-generic-mixed"]}}"#
        ),
    );
    let output = Command::new(binary()).args(["completeness", completeness.to_str().unwrap()]).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", String::from_utf8_lossy(&output.stderr));
    let summary: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(summary["points"], 3);
    let check = write("check.json", format!(r#"{{"catalogue": {catalogue:?}, "transcript": {transcript:?}, "threads": 4}}"#));
    let output = Command::new(binary()).args(["check", check.to_str().unwrap()]).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", String::from_utf8_lossy(&output.stderr));
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!((report["status"].as_str(), report["replayed"].as_bool()), (Some("pass"), Some(false)));
    let compare = write(
        "compare.json",
        format!(
            r#"{{"catalogue": {catalogue:?}, "baseline": {transcript:?}, "candidate": {transcript:?}, "threads": 4,
                "worse_by": 1, "fail_on_regression": true}}"#
        ),
    );
    assert_eq!(status(&["compare", compare.to_str().unwrap()]).0, 0);
    let pilot_transcript = dir.join("pilot.jsonl");
    let pilot = write(
        "pilot.json",
        format!(
            r#"{{"roots": {{"paths": ["0110", "10010"]}}, "threads": 2, "depth_limit": 14, "max_decisions": 400,
                "transcript": {pilot_transcript:?}}}"#
        ),
    );
    let output = Command::new(binary()).args(["pilot", pilot.to_str().unwrap()]).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", String::from_utf8_lossy(&output.stderr));
    assert!(pilot_transcript.exists());
    // `targets` regenerates the Local covers' targets in place: removed, they
    // come back as the shipped file; a second run changes nothing.
    let shipped = std::fs::read(Path::new(env!("CARGO_MANIFEST_DIR")).join("catalogue").join("targets.json")).unwrap();
    let mut spelling: serde_json::Value = serde_json::from_slice(&shipped).unwrap();
    spelling["targets"].as_array_mut().unwrap().retain(|t| !t["id"].as_str().unwrap().starts_with("local-") || t["id"] == "local-overlap");
    let targets = write("targets-catalogue.json", spelling.to_string());
    let config = write("targets.json", format!(r#"{{"catalogue": {targets:?}}}"#));
    for changed in [true, false] {
        let output = Command::new(binary()).args(["targets", config.to_str().unwrap()]).output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", String::from_utf8_lossy(&output.stderr));
        let summary: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!((summary["changed"].as_bool(), summary["targets"].as_u64()), (Some(changed), Some(45)));
        assert_eq!(std::fs::read(&targets).unwrap(), shipped);
    }
}
