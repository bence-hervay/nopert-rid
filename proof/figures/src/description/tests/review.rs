//! Adversarial tests of descriptions.
use super::*;
use crate::content::slice::{Category, Parts};
use serde_json::{json, Value};
use std::sync::atomic::{AtomicUsize, Ordering};

static CALLS: AtomicUsize = AtomicUsize::new(0);

fn counting() -> Result<Box<dyn Classifier>, BoxError> {
    Ok(Box::new(|b: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        CALLS.fetch_add(1, Ordering::SeqCst);
        Ok(if *b.axes()[0].lo() >= rid::arithmetic::exact::frac(1, 2) { Some(Category::Global) } else { None })
    }))
}

fn cells_styles() -> Value {
    json!({"Domain": "c", "Global": "c", "Local": "c",
           "Exotic": {"square": "c", "pentagon": "c", "arc+": "c", "arc-": "c", "endpoint+": "c", "endpoint-": "c", "crossing+": "c", "crossing-": "c"},
           "unresolved": "c", "outside": "c"})
}

fn slice_source(name: &str, max_depth: usize) -> Value {
    let mut sources = serde_json::Map::new();
    sources.insert(name.to_owned(), json!({"kind": "slice", "threads": 2, "slice": {
        "origin": ["0", "1/10", "0", "0", "0"],
        "horizontal": ["1", "0", "0", "0", "0"],
        "vertical": ["0", "0", "0", "0", "1"],
        "plot": [["0", "2/3"], ["-2/5", "2/5"]],
        "thickness": ["0", "0", "0", "0", "0"],
        "max_depth": max_depth,
        "max_cells": 65536, "components": {"Domain": true, "Global": true, "Local": true, "Exotic": ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]}}}));
    Value::Object(sources)
}

fn description(name: &str, max_depth: usize) -> Value {
    json!({
        "png": {"converter": "/nonexistent", "dpi": 100},
        "canvas": {"width_mm": "100", "height_mm": "100", "background": null},
        "sources": slice_source(name, max_depth),
        "styles": {"c": {"fill": "#123456", "stroke": null, "opacity": "1"},
                   "line": {"fill": null, "stroke": {"color": "#000000", "width_mm": "1/10"}, "opacity": "1"}},
        "panels": {"plot": {"viewport_mm": [["10", "60"], ["10", "70"]], "world": [["0", "2/3"], ["-2/5", "2/5"]]},
                   "inset": {"viewport_mm": [["70", "90"], ["10", "30"]], "world": [["0", "1/10"], ["0", "1/10"]]}},
        "layers": [{"kind": "cells", "panel": "plot", "source": name, "styles": cells_styles()},
                   {"kind": "inset", "overview": "plot", "detail": "inset", "connectors": [], "mark_mm": "0", "style": "line"}]
    })
}

fn render(v: &Value) -> Result<Rendered, DescriptionError> {
    Description::parse(v.to_string().as_bytes())?.render(Path::new("."), "policy", &counting)
}

/// Source names become file names (`<name>.partition.json` in the output
/// directory), so separators, `..` and absolute paths are refused: names are
/// `[a-z0-9][a-z0-9-]*`, checked while parsing, for
/// sources, styles and panels alike.
#[test]
fn names_that_are_not_file_name_safe_are_refused() {
    for name in ["../escaped", "/tmp/absolute", "a/b", "", ".hidden", "x\ny", "Upper", "-dash", "a.b", "a_b"] {
        assert!(matches!(render(&description(name, 2)), Err(DescriptionError::Json(_))), "{name:?}");
        let mut styles = description("slice", 2);
        styles["styles"][name] = json!({"fill": "#000000", "stroke": null, "opacity": "1"});
        assert!(matches!(render(&styles), Err(DescriptionError::Json(_))), "style {name:?}");
    }
    for name in ["slice", "a", "0x", "slice-2"] {
        assert!(render(&description(name, 2)).unwrap().partitions.contains_key(name), "{name:?}");
    }
}

/// Cheap, purely structural refusals (an inset outside its overview, a filled
/// outline) are found before any source is prepared: the classifier is
/// never called.
#[test]
fn structural_refusals_come_before_any_subdivision() {
    // A counter of this test's own: other tests classify concurrently.
    let calls = std::sync::Arc::new(AtomicUsize::new(0));
    let counter = calls.clone();
    let counted = move || -> Result<Box<dyn Classifier>, BoxError> {
        let counter = counter.clone();
        Ok(Box::new(move |_: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
            counter.fetch_add(1, Ordering::SeqCst);
            Ok(None)
        }))
    };
    let render = |v: &Value| Description::parse(v.to_string().as_bytes())?.render(Path::new("."), "policy", &counted);
    let mut v = description("slice", 8);
    // An inset whose detail world is outside the overview's.
    v["panels"]["inset"]["world"] = json!([["1", "11/10"], ["0", "1/10"]]);
    let error = render(&v).err().unwrap();
    assert!(matches!(error, DescriptionError::Inset { .. }), "{error:?}");
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    // The same for a filled inset style.
    let mut v = description("slice", 8);
    v["layers"][1]["style"] = json!("c");
    let error = render(&v).err().unwrap();
    assert!(matches!(error, DescriptionError::FilledOutline { .. }), "{error:?}");
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    // `check` alone never prepares a source.
    Description::parse(description("slice", 8).to_string().as_bytes()).unwrap().check().unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

/// Every field of slice sources, partition sources and cells layers is
/// required, and unknown fields are refused there too (the builder's test
/// covers only configuration and box sources and scene/inset layers).
#[test]
fn slice_descriptions_need_every_field() {
    let good = description("slice", 2);
    assert!(render(&good).is_ok());
    let mut paths: Vec<Vec<String>> = Vec::new();
    fn walk(v: &Value, at: Vec<String>, out: &mut Vec<Vec<String>>) {
        match v {
            Value::Object(map) => {
                let named = matches!(at.as_slice(), [one] if one == "sources" || one == "styles" || one == "panels");
                for (k, inner) in map {
                    let mut p = at.clone();
                    p.push(k.clone());
                    if !named {
                        out.push(p.clone());
                    }
                    walk(inner, p, out);
                }
            }
            Value::Array(items) => {
                for (k, inner) in items.iter().enumerate() {
                    let mut p = at.clone();
                    p.push(k.to_string());
                    walk(inner, p, out);
                }
            }
            _ => {}
        }
    }
    walk(&good, Vec::new(), &mut paths);
    let locate = |v: &mut Value, path: &[String]| -> *mut Value {
        let mut cur = v;
        for key in path {
            cur = match cur {
                Value::Array(items) => &mut items[key.parse::<usize>().unwrap()],
                other => &mut other[key.as_str()],
            };
        }
        cur as *mut Value
    };
    for path in &paths {
        let mut missing = good.clone();
        let (last, parent) = path.split_last().unwrap();
        unsafe { (*locate(&mut missing, parent)).as_object_mut().unwrap().remove(last) };
        assert!(matches!(render(&missing), Err(DescriptionError::Json(_))), "without {path:?}");
        let mut extra = good.clone();
        let target = unsafe { &mut *locate(&mut extra, path) };
        if target.is_object() {
            target["unexpected"] = json!(1);
            assert!(matches!(render(&extra), Err(DescriptionError::Json(_))), "extra in {path:?}");
        }
    }
    // The partition source kind.
    let mut v = description("slice", 2);
    v["sources"]["slice"] = json!({"kind": "partition", "file": "x.json", "extra": 1});
    assert!(matches!(render(&v), Err(DescriptionError::Json(_))));
    v["sources"]["slice"] = json!({"kind": "partition"});
    assert!(matches!(render(&v), Err(DescriptionError::Json(_))));
}

/// An inset between panels that show different worlds (plot parameters of a
/// slice and the screen of a configuration), or a panel showing both, is
/// refused.
#[test]
fn inset_between_unrelated_worlds_is_refused() {
    let mut v = description("slice", 2);
    v["sources"]["point"] = json!({"kind": "configuration", "coordinates": [["1/5", "0"], ["1/10", "0"], ["1/8", "0"], ["1/16", "0"], ["-1/32", "0"]]});
    v["layers"].as_array_mut().unwrap().push(json!({"kind": "scene", "panel": "inset", "source": "point", "body": "hole", "geometry": "wireframe", "style": "line"}));
    assert!(matches!(render(&v), Err(DescriptionError::WorldMismatch { panel, .. }) if panel == "inset"));
    let mut v = description("slice", 2);
    v["sources"]["point"] = json!({"kind": "configuration", "coordinates": [["1/5", "0"], ["1/10", "0"], ["1/8", "0"], ["1/16", "0"], ["-1/32", "0"]]});
    v["layers"].as_array_mut().unwrap().push(json!({"kind": "scene", "panel": "plot", "source": "point", "body": "hole", "geometry": "wireframe", "style": "line"}));
    assert!(matches!(render(&v), Err(DescriptionError::WorldMismatch { panel, .. }) if panel == "plot"));
}

/// Unused panels and styles are refused, like unused sources.
#[test]
fn unused_panels_and_styles_are_refused() {
    let mut v = description("slice", 2);
    v["styles"]["never"] = json!({"fill": "#000000", "stroke": null, "opacity": "1"});
    assert!(matches!(render(&v), Err(DescriptionError::UnusedStyle(name)) if name == "never"));
    let mut v = description("slice", 2);
    v["panels"]["spare"] = json!({"viewport_mm": [["70", "90"], ["40", "60"]], "world": [["0", "1"], ["0", "1"]]});
    assert!(matches!(render(&v), Err(DescriptionError::UnusedPanel(name)) if name == "spare"));
}

/// Huge PNG resolutions and thread counts are accepted without bound.
#[test]
fn resource_fields_are_unbounded() {
    let mut v = description("slice", 2);
    v["png"]["dpi"] = json!(4294967295u64);
    v["sources"]["slice"]["threads"] = json!(1_000_000);
    let parsed = Description::parse(v.to_string().as_bytes()).unwrap();
    assert_eq!(parsed.png.dpi.get(), u32::MAX);
}
