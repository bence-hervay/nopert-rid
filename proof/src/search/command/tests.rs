//! Unit tests of the command wiring that need no certificate file; the
//! end-to-end and interruption tests are in the crate's `tests/` directory.
use super::*;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Tagged {
    tag: u32,
}

impl RecordData for Tagged {}

/// Accepts the paths of length at least 2 with `tag` = path length.
struct Rule;

impl Components for Rule {
    type Data = Tagged;
    fn report(&self) -> Value {
        json!({})
    }
    fn check(&self, path: &str) -> Result<Option<Tagged>, BoxError> {
        if path.len() < 2 {
            return Ok(None);
        }
        Ok(Some(Tagged { tag: path.len() as u32 }))
    }
    fn holds(&self, path: &str, data: &Tagged) -> Result<(), BoxError> {
        if path.len() >= 2 && data.tag == path.len() as u32 {
            Ok(())
        } else {
            Err("false record".into())
        }
    }
}

#[test]
fn command_lines_are_usage_errors_before_any_file_is_read() {
    let usage = |args: &[&str]| {
        let args: Vec<OsString> = args.iter().map(OsString::from).collect();
        matches!(command(&args), Err(Error::Usage(_)))
    };
    assert!(usage(&[]));
    assert!(usage(&["search"]));
    assert!(usage(&["search", "a.json", "b.json"]));
    assert!(usage(&["run", "/nonexistent/config.json"]));
    assert!(usage(&["Check", "/nonexistent/config.json"]));
    let args = [OsString::from("check"), OsString::from("/nonexistent/config.json")];
    assert!(matches!(command(&args), Err(Error::Config(config::Error::Read { .. }))));
}

#[test]
fn a_depth_limit_beyond_the_header_limit_is_refused_before_the_components() {
    let config = config::Search {
        certificate: PathBuf::from("/nonexistent/never-created.cert"),
        root: String::new(),
        threads: NonZeroUsize::MIN,
        window: NonZeroUsize::MIN,
        max_depth: 10,
        depth_limit: Some(11),
        max_decisions: None,
    };
    let make = |_: NonZeroUsize, _: &AtomicBool| -> Result<Rule, BoxError> {
        panic!("the components must not be built")
    };
    let error = search(&config, &AtomicBool::new(false), make).unwrap_err();
    assert!(matches!(error, Error::Invalid(_)), "{error}");
    assert!(error.to_string().contains("depth_limit 11 exceeds"), "{error}");
}

#[test]
fn a_depth_limit_shallower_than_the_root_is_refused_before_the_components() {
    let config = |depth_limit| config::Search {
        certificate: PathBuf::from("/nonexistent/never-created.cert"),
        root: "0101".into(),
        threads: NonZeroUsize::MIN,
        window: NonZeroUsize::MIN,
        max_depth: 10,
        depth_limit,
        max_decisions: None,
    };
    let make = |_: NonZeroUsize, _: &AtomicBool| -> Result<Rule, BoxError> {
        panic!("the components must not be built")
    };
    for shallow in [0, 2, 3] {
        let error = search(&config(Some(shallow)), &AtomicBool::new(false), make).unwrap_err();
        assert!(matches!(error, Error::Invalid(_)), "{error}");
        let expected = format!("depth_limit {shallow} is shallower than the root");
        assert!(error.to_string().contains(&expected), "{error}");
        assert!(error.to_string().contains("\"0101\" (depth 4)"), "{error}");
    }
    // At the root's own depth the root is evaluated but not split: accepted
    // as far as the configuration goes (here the missing folder stops it).
    let built = |_: NonZeroUsize, _: &AtomicBool| Ok(Rule);
    for fine in [Some(4), None] {
        let error = search(&config(fine), &AtomicBool::new(false), built).unwrap_err();
        assert!(matches!(error, Error::Certificate(_)), "{error}");
    }
}

#[test]
fn a_stop_before_or_during_preparation_is_reported_as_such() {
    let prepared = |_: NonZeroUsize, _: &AtomicBool| Ok(Rule);
    let stopped = AtomicBool::new(true);
    let config = config::Prepare {
        threads: NonZeroUsize::MIN,
    };
    let error = prepare(&config, &stopped, prepared).unwrap_err();
    assert!(matches!(error, Error::Stopped("during preparation")), "{error}");
    let failing = |_: NonZeroUsize, _: &AtomicBool| -> Result<Rule, BoxError> { Err("no".into()) };
    let error = prepare(&config, &AtomicBool::new(false), failing).unwrap_err();
    assert!(matches!(error, Error::Components(_)), "{error}");
}
