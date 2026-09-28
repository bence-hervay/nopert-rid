use super::*;
use crate::content::slice::{Category, Outcome, Parts};
use rid::components::record::ExoticCover;
use serde_json::{json, Value};
use std::cell::Cell;

fn base() -> Value {
    json!({
        "png": {"converter": "/nonexistent/rsvg-convert", "dpi": 150},
        "canvas": {"width_mm": "180", "height_mm": "100", "background": "#ffffff"},
        "sources": {
            "point": {"kind": "configuration", "coordinates": [["1/5", "0"], ["1/10", "0"], ["1/8", "0"], ["1/16", "0"], ["-1/32", "0"]]},
            "box": {"kind": "box", "axes": [["3/16", "7/32"], ["3/32", "7/64"], ["1/8", "9/64"], ["1/16", "5/64"], ["-1/32", "-1/64"]]}
        },
        "styles": {
            "hole": {"fill": null, "stroke": {"color": "#ff8000", "width_mm": "7/50"}, "opacity": "1"},
            "plug": {"fill": null, "stroke": {"color": "#0080ff", "width_mm": "7/100"}, "opacity": "1"},
            "shadow": {"fill": "#0080ff", "stroke": null, "opacity": "1/5"},
            "outline": {"fill": null, "stroke": {"color": "#404040", "width_mm": "1/5"}, "opacity": "1"}
        },
        "panels": {
            "overview": {"viewport_mm": [["5", "95"], ["5", "95"]], "world": [["-5", "5"], ["-5", "5"]]},
            "detail": {"viewport_mm": [["110", "170"], ["20", "80"]], "world": [["37/50", "43/50"], ["-441/100", "-429/100"]]}
        },
        "layers": [
            {"kind": "scene", "panel": "overview", "source": "box", "body": "plug", "geometry": "edge-envelopes", "style": "shadow"},
            {"kind": "scene", "panel": "overview", "source": "point", "body": "hole", "geometry": "wireframe", "style": "hole"},
            {"kind": "scene", "panel": "overview", "source": "point", "body": "plug", "geometry": "silhouette", "style": "plug"},
            {"kind": "scene", "panel": "detail", "source": "point", "body": "plug", "geometry": "wireframe", "style": "plug"},
            {"kind": "inset", "overview": "overview", "detail": "detail",
             "connectors": [{"from": "top-right", "to": "top-left"}, {"from": "bottom-right", "to": "bottom-left"}], "mark_mm": "0", "style": "outline"}
        ]
    })
}

fn no_components() -> Result<Box<dyn Classifier>, BoxError> {
    Err("no components in this test".into())
}

fn render(value: &Value) -> Result<Rendered, DescriptionError> {
    Description::parse(value.to_string().as_bytes())?.render(Path::new("."), "policy", &no_components)
}

/// The groups of an SVG: each group's attributes and its paths' data.
fn groups(svg: &str) -> Vec<(String, Vec<String>)> {
    svg.split("<g ")
        .skip(1)
        .map(|group| {
            let (attributes, rest) = group.split_once('>').unwrap();
            let body = rest.split("</g>").next().unwrap();
            let paths = body
                .split("<path d=\"")
                .skip(1)
                .map(|p| p.split('"').next().unwrap().to_owned())
                .collect();
            (attributes.to_owned(), paths)
        })
        .collect()
}

#[test]
fn a_description_renders_its_layers_in_order() {
    let rendered = render(&base()).unwrap();
    assert!(rendered.partitions.is_empty());
    let groups = groups(&rendered.svg);
    assert_eq!(groups.len(), 5);
    let expected = [
        ("fill=\"#0080ff\"", "clip-1"),
        ("stroke=\"#ff8000\"", "clip-2"),
        ("stroke=\"#0080ff\"", "clip-3"),
        ("stroke=\"#0080ff\"", "clip-4"),
        ("stroke=\"#404040\"", ""),
    ];
    for ((attributes, _), (style, clip)) in groups.iter().zip(expected) {
        assert!(attributes.contains(style), "{attributes}");
        assert_eq!(attributes.contains("clip-path"), !clip.is_empty());
        assert!(clip.is_empty() || attributes.contains(clip));
    }
    // 120 edge envelopes and 120 edges in the overview, one silhouette.
    assert_eq!(groups[0].1.len(), 120);
    assert_eq!(groups[1].1.len(), 120);
    assert_eq!(groups[2].1.len(), 1);
    // The detail shows the few edges near its window.
    assert!(!groups[3].1.is_empty() && groups[3].1.len() < 20);
    // The inset: the window x ∈ [0.74, 0.86], y ∈ [-4.41, -4.29] at 9 mm per
    // unit around the overview's centre (50, 50), the frame, two connectors.
    assert_eq!(
        groups[4].1,
        vec![
            "M56.66 88.61L57.74 88.61L57.74 89.69L56.66 89.69Z",
            "M110 20L170 20L170 80L110 80Z",
            "M57.74 88.61L110 20",
            "M57.74 89.69L110 80",
        ]
    );
}

/// Every object path of a JSON value, excluding the entries of named maps
/// (removing one is a missing name, not a missing field).
fn fields(value: &Value, at: Vec<String>, out: &mut Vec<Vec<String>>) {
    match value {
        Value::Object(map) => {
            let named = matches!(at.as_slice(), [one] if one == "sources" || one == "styles" || one == "panels");
            for (key, inner) in map {
                let mut path = at.clone();
                path.push(key.clone());
                if !named {
                    out.push(path.clone());
                }
                fields(inner, path, out);
            }
        }
        Value::Array(items) => {
            for (k, inner) in items.iter().enumerate() {
                let mut path = at.clone();
                path.push(k.to_string());
                fields(inner, path, out);
            }
        }
        _ => {}
    }
}

fn locate<'a>(value: &'a mut Value, path: &[String]) -> &'a mut Value {
    path.iter().fold(value, |v, key| match v {
        Value::Array(items) => &mut items[key.parse::<usize>().unwrap()],
        other => &mut other[key.as_str()],
    })
}

#[test]
fn every_field_is_required_and_unknown_fields_are_refused() {
    let good = base();
    let mut paths = Vec::new();
    fields(&good, Vec::new(), &mut paths);
    assert!(paths.len() > 60, "{}", paths.len());
    for path in &paths {
        let mut missing = good.clone();
        let (last, parent) = path.split_last().unwrap();
        locate(&mut missing, parent).as_object_mut().unwrap().remove(last);
        assert!(matches!(render(&missing), Err(DescriptionError::Json(_))), "without {path:?}");
        if locate(&mut good.clone(), path).is_object() {
            let mut extra = good.clone();
            locate(&mut extra, path)["unexpected"] = json!(1);
            assert!(matches!(render(&extra), Err(DescriptionError::Json(_))), "extra field in {path:?}");
        }
    }
    let mut extra = good.clone();
    extra["unexpected"] = json!(1);
    assert!(matches!(render(&extra), Err(DescriptionError::Json(_))));
}

#[test]
fn numbers_names_and_kinds_are_strict() {
    let text = base().to_string();
    let duplicated = text.replacen("\"point\":{", "\"box\":{\"kind\":\"configuration\",\"coordinates\":[[\"0\",\"0\"],[\"0\",\"0\"],[\"0\",\"0\"],[\"0\",\"0\"],[\"0\",\"0\"]]},\"point\":{", 1);
    assert_ne!(duplicated, text);
    assert!(matches!(Description::parse(duplicated.as_bytes()), Err(DescriptionError::Json(_))));
    for (path, value) in [
        (vec!["canvas", "width_mm"], json!(180)),
        (vec!["canvas", "width_mm"], json!("180.0")),
        (vec!["canvas", "background"], json!("white")),
        (vec!["png", "dpi"], json!(0)),
        (vec!["png", "dpi"], json!("150")),
        (vec!["sources", "point", "kind"], json!("point")),
        (vec!["sources", "point", "coordinates", "0"], json!(["1/5"])),
        (vec!["layers", "0", "body"], json!("Plug")),
        (vec!["layers", "0", "geometry"], json!("edges")),
        (vec!["layers", "4", "connectors", "0", "from"], json!("top")),
    ] {
        let mut bad = base();
        let path: Vec<String> = path.iter().map(|s| s.to_string()).collect();
        *locate(&mut bad, &path) = value;
        assert!(matches!(render(&bad), Err(DescriptionError::Json(_))), "{path:?}");
    }
}

#[test]
fn names_must_resolve_and_every_source_must_be_drawn() {
    let check = |edit: &dyn Fn(&mut Value), expected: &dyn Fn(&DescriptionError) -> bool| {
        let mut v = base();
        edit(&mut v);
        let error = render(&v).err().expect("refused");
        assert!(expected(&error), "{error:?}");
    };
    check(&|v| v["layers"][1]["panel"] = json!("missing"), &|e| matches!(e, DescriptionError::Unknown { layer: 1, what: "panel", .. }));
    check(&|v| v["layers"][4]["detail"] = json!("missing"), &|e| matches!(e, DescriptionError::Unknown { layer: 4, what: "panel", .. }));
    check(&|v| v["layers"][2]["source"] = json!("missing"), &|e| matches!(e, DescriptionError::Unknown { layer: 2, what: "source", .. }));
    check(&|v| v["layers"][3]["style"] = json!("missing"), &|e| matches!(e, DescriptionError::Unknown { layer: 3, what: "style", .. }));
    check(&|v| {
        v["layers"].as_array_mut().unwrap().remove(0);
    }, &|e| matches!(e, DescriptionError::UnusedSource(name) if name == "box"));
    check(&|v| v["layers"][0]["geometry"] = json!("wireframe"), &|e| matches!(e, DescriptionError::Mismatch { layer: 0, .. }));
    check(&|v| v["layers"][1]["geometry"] = json!("vertex-envelopes"), &|e| matches!(e, DescriptionError::Mismatch { layer: 1, .. }));
    check(&|v| v["layers"][1]["style"] = json!("shadow"), &|e| matches!(e, DescriptionError::FilledOutline { layer: 1, .. }));
    check(&|v| v["layers"][4]["style"] = json!("shadow"), &|e| matches!(e, DescriptionError::FilledOutline { layer: 4, .. }));
    check(&|v| v["panels"]["detail"]["world"][0] = json!(["6", "153/25"]), &|e| matches!(e, DescriptionError::Inset { layer: 4, .. }));
    check(&|v| v["panels"]["detail"]["world"][0] = json!(["37/50", "1"]), &|e| matches!(e, DescriptionError::Panel { name, .. } if name == "detail"));
    check(&|v| v["panels"]["detail"]["world"][0] = json!(["1", "37/50"]), &|e| matches!(e, DescriptionError::Panel { .. }));
    check(&|v| v["canvas"]["height_mm"] = json!("0"), &|e| matches!(e, DescriptionError::Canvas(_)));
    check(&|v| v["sources"]["box"]["axes"][0] = json!(["1", "0"]), &|e| matches!(e, DescriptionError::Source { name, .. } if name == "box"));
    // A filled silhouette is a shadow.
    let mut v = base();
    v["layers"][2]["style"] = json!("shadow");
    assert!(render(&v).is_ok());
}

fn slice_description(sources: Value) -> Value {
    let mut styles = serde_json::Map::new();
    for k in 0..13 {
        styles.insert(format!("c{k}"), json!({"fill": format!("#0000{k:02x}"), "stroke": null, "opacity": "1"}));
    }
    let cells = |source: &str| json!({"kind": "cells", "panel": "plot", "source": source, "styles": {
        "Domain": "c0", "Global": "c1", "Local": "c2",
        "Exotic": {"square": "c3", "pentagon": "c4", "arc+": "c5", "arc-": "c6", "endpoint+": "c7", "endpoint-": "c8", "crossing+": "c9", "crossing-": "c10"},
        "unresolved": "c11", "outside": "c12"}});
    let layers: Vec<Value> = sources.as_object().unwrap().keys().map(|name| cells(name)).collect();
    json!({
        "png": {"converter": "/nonexistent/rsvg-convert", "dpi": 150},
        "canvas": {"width_mm": "100", "height_mm": "100", "background": null},
        "sources": sources,
        "styles": styles,
        "panels": {"plot": {"viewport_mm": [["10", "60"], ["10", "70"]], "world": [["0", "2/3"], ["-2/5", "2/5"]]}},
        "layers": layers
    })
}

fn slice_source(max_depth: usize) -> Value {
    json!({"kind": "slice", "threads": 3, "slice": {
        "origin": ["0", "1/10", "0", "0", "0"],
        "horizontal": ["1", "0", "0", "0", "0"],
        "vertical": ["0", "0", "0", "0", "1"],
        "plot": [["0", "2/3"], ["-2/5", "2/5"]],
        "thickness": ["0", "0", "0", "0", "0"],
        "max_depth": max_depth,
        "max_cells": 65536, "components": {"Domain": true, "Global": true, "Local": true, "Exotic": ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]}}})
}

/// Global for s ≥ 1/2, Exotic crossing for r₃ ≥ 1/5, otherwise none.
fn stand_in() -> Result<Box<dyn Classifier>, BoxError> {
    Ok(Box::new(|b: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        let axes = b.axes();
        Ok(if *axes[0].lo() >= rid::arithmetic::exact::frac(1, 2) {
            Some(Category::Global)
        } else if *axes[4].lo() >= rid::arithmetic::exact::frac(1, 5) {
            Some(Category::Exotic { cover: ExoticCover::CrossingPlus })
        } else {
            None
        })
    }))
}

#[test]
fn slices_are_computed_once_with_the_components_and_drawn_by_outcome() {
    let v = slice_description(json!({"first": slice_source(4), "second": slice_source(3)}));
    let calls = Cell::new(0);
    let factory = || {
        calls.set(calls.get() + 1);
        stand_in()
    };
    let description = Description::parse(v.to_string().as_bytes()).unwrap();
    let rendered = description.render(Path::new("."), "policy", &factory).unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(rendered.partitions.len(), 2);
    let first = &rendered.partitions["first"];
    // Groups follow the fixed outcome order with each outcome's own style.
    let groups = groups(&rendered.svg);
    let present: Vec<Outcome> = outcomes().into_iter().filter(|o| first.cells().iter().any(|c| c.outcome == *o)).collect();
    assert_eq!(present, vec![Outcome::Eliminated(Category::Global), Outcome::Eliminated(Category::Exotic { cover: ExoticCover::CrossingPlus }), Outcome::Unresolved]);
    for ((attributes, paths), (colour, outcome)) in groups.iter().zip([("#000001", present[0]), ("#000009", present[1]), ("#00000b", present[2])]) {
        assert!(attributes.contains(colour), "{attributes}");
        assert_eq!(paths.len(), first.cells().iter().filter(|c| c.outcome == outcome).count());
    }
    // Without components a slice cannot be computed.
    let error = description.render(Path::new("."), "policy", &no_components).err().unwrap();
    assert!(matches!(error, DescriptionError::Source { error: SourceError::Components(_), .. }));
}

#[test]
fn saved_partitions_draw_the_same_cells_without_the_components() {
    let computed = slice_description(json!({"slice": slice_source(4)}));
    let description = Description::parse(computed.to_string().as_bytes()).unwrap();
    let rendered = description.render(Path::new("."), "policy", &stand_in).unwrap();
    let directory = std::env::temp_dir().join(format!("rid-figures-description-test-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("cells.json"), rendered.partitions["slice"].to_json("policy")).unwrap();
    let saved = slice_description(json!({"slice": {"kind": "partition", "file": "cells.json"}}));
    let description = Description::parse(saved.to_string().as_bytes()).unwrap();
    let again = description.render(&directory, "policy", &no_components).unwrap();
    assert_eq!(again.svg, rendered.svg);
    assert!(again.partitions.is_empty());
    let error = description.render(&directory, "another policy", &no_components).err().unwrap();
    assert!(matches!(error, DescriptionError::Source { error: SourceError::Slice(SliceError::Policy { .. }), .. }));
    let error = description.render(Path::new("/nonexistent"), "policy", &no_components).err().unwrap();
    assert!(matches!(error, DescriptionError::Source { error: SourceError::Read { .. }, .. }));
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn cells_need_a_slice_and_scenes_refuse_slices() {
    let mut v = slice_description(json!({"slice": slice_source(2)}));
    v["sources"]["point"] = base()["sources"]["point"].clone();
    v["layers"][0]["source"] = json!("point");
    v["layers"].as_array_mut().unwrap().push(json!({"kind": "scene", "panel": "plot", "source": "slice", "body": "hole", "geometry": "wireframe", "style": "c0"}));
    let description = Description::parse(v.to_string().as_bytes()).unwrap();
    let error = description.render(Path::new("."), "policy", &stand_in).err().unwrap();
    assert!(matches!(error, DescriptionError::Mismatch { layer: 0, .. }), "{error:?}");
}

mod review;
