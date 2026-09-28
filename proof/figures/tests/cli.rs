//! End-to-end tests of `rid-figures`: refusals leave no output behind, and
//! (with a real converter, see `renders_the_examples`) complete renders.
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BINARY: &str = env!("CARGO_BIN_EXE_rid-figures");

fn scratch(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join("cli").join(name);
    if directory.exists() {
        fs::remove_dir_all(&directory).unwrap();
    }
    fs::create_dir_all(&directory).unwrap();
    directory
}

fn run(arguments: &[&Path]) -> Output {
    Command::new(BINARY).args(arguments).output().unwrap()
}

fn example(name: &str) -> serde_json::Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples").join(format!("{name}.json"));
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

/// The configuration example with another converter, written to `directory`.
fn description(directory: &Path, converter: &Path) -> PathBuf {
    let mut value = example("configuration");
    value["png"]["converter"] = converter.to_str().unwrap().into();
    let path = directory.join("figure.json");
    fs::write(&path, value.to_string()).unwrap();
    path
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn wrong_arguments_are_a_usage_error() {
    let directory = scratch("usage");
    for arguments in [vec![], vec![directory.as_path()], vec![directory.as_path(); 3]] {
        let output = run(&arguments);
        assert_eq!(output.status.code(), Some(2));
        assert!(stderr(&output).starts_with("usage: rid-figures"));
    }
}

#[test]
fn existing_outputs_are_never_touched() {
    let directory = scratch("existing");
    let figure = description(&directory, Path::new("/nonexistent/converter"));
    let output_directory = directory.join("out");
    fs::create_dir(&output_directory).unwrap();
    fs::write(output_directory.join("keep"), "kept").unwrap();
    let output = run(&[&figure, &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).contains("already exists"));
    assert_eq!(fs::read_to_string(output_directory.join("keep")).unwrap(), "kept");
    // A leftover staging directory is refused as well.
    let other = directory.join("other");
    fs::create_dir(directory.join("other.incomplete")).unwrap();
    let output = run(&[&figure, &other]);
    assert_eq!(output.status.code(), Some(1));
    assert!(!other.exists());
}

#[test]
fn refused_renders_leave_no_directory() {
    let directory = scratch("refused");
    let output_directory = directory.join("out");
    let staging = directory.join("out.incomplete");
    // An invalid description.
    let invalid = directory.join("invalid.json");
    fs::write(&invalid, "{\"png\": {}}").unwrap();
    let output = run(&[&invalid, &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).contains("invalid description"));
    assert!(!output_directory.exists() && !staging.exists());
    // A missing description.
    let output = run(&[&directory.join("missing.json"), &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    // A converter that does not exist.
    let figure = description(&directory, Path::new("/nonexistent/converter"));
    let output = run(&[&figure, &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).contains("converter"), "{}", stderr(&output));
    assert!(!output_directory.exists() && !staging.exists());
    // A converter whose `--version` fails (this program itself, which
    // refuses the argument) is refused before rendering.
    let figure = description(&directory, Path::new(BINARY));
    let output = run(&[&figure, &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).contains("--version failed"), "{}", stderr(&output));
    assert!(!output_directory.exists() && !staging.exists());
    // A converter that answers `--version` but writes no PNG (`true`): the
    // staging directory is removed. With a trailing separator on the output
    // path the staging directory is still a sibling.
    let figure = description(&directory, Path::new("true"));
    let slashed = PathBuf::from(format!("{}/", output_directory.display()));
    for target in [&output_directory, &slashed] {
        let output = run(&[&figure, target]);
        assert_eq!(output.status.code(), Some(1));
        assert!(stderr(&output).contains("no PNG was written"), "{}", stderr(&output));
        assert!(!output_directory.exists() && !staging.exists());
    }
    // A slice example with a converter that is not there is refused before
    // any subdivision (quickly, and without loading the components).
    let mut slice = example("slice");
    slice["png"]["converter"] = "/nonexistent/converter".into();
    let file = directory.join("slice.json");
    fs::write(&file, slice.to_string()).unwrap();
    let started = std::time::Instant::now();
    let output = run(&[&file, &output_directory]);
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).contains("converter"), "{}", stderr(&output));
    assert!(started.elapsed().as_secs() < 5);
    assert!(!output_directory.exists() && !staging.exists());
}

/// Renders the configuration and box examples with the converter named by
/// `RID_FIGURES_CONVERTER` and compares their SVGs with the committed ones:
///
/// ```sh
/// RID_FIGURES_CONVERTER=/path/to/rsvg-convert cargo test -- --ignored
/// ```
#[test]
#[ignore = "needs RID_FIGURES_CONVERTER, the path of rsvg-convert"]
fn renders_the_examples() {
    let converter = PathBuf::from(std::env::var_os("RID_FIGURES_CONVERTER").expect("RID_FIGURES_CONVERTER"));
    let directory = scratch("examples");
    for name in ["configuration", "box"] {
        let mut value = example(name);
        value["png"]["converter"] = converter.to_str().unwrap().into();
        let figure = directory.join(format!("{name}.json"));
        fs::write(&figure, value.to_string()).unwrap();
        let output_directory = directory.join(name);
        let output = run(&[&figure, &output_directory]);
        assert!(output.status.success(), "{}", stderr(&output));
        let mut files: Vec<String> = fs::read_dir(&output_directory)
            .unwrap()
            .map(|e| e.unwrap().file_name().into_string().unwrap())
            .collect();
        files.sort();
        assert_eq!(files, ["figure.png", "figure.svg"]);
        let png = fs::read(output_directory.join("figure.png")).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        let committed = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples").join(name).join("figure.svg");
        assert_eq!(
            fs::read_to_string(output_directory.join("figure.svg")).unwrap(),
            fs::read_to_string(committed).unwrap(),
            "{name}: the rendering differs from the committed example"
        );
    }
}
