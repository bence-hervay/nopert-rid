//! Every record shape shown in this module's doc comments is itself a
//! canonical record, byte for byte.
use super::*;

/// The lines of `text` that begin (after `///` and spaces) with a record.
fn documented(text: &str) -> Vec<String> {
    text.lines()
        .map(|l| l.trim_start().trim_start_matches("///").trim())
        .filter(|l| l.starts_with(r#"{"component""#))
        .map(str::to_owned)
        .collect()
}

#[test]
fn every_documented_record_is_a_canonical_record() {
    let shown = documented(include_str!("../mod.rs"));
    assert!(!shown.is_empty(), "{shown:?}");
    for text in &shown {
        assert!(canonical(text).is_some(), "the documented record {text} does not parse canonically");
    }
}
