use super::*;
use crate::elimination::witness::{Edge, Gap, WitnessError};
use crate::problem::geometry::{self, VERTEX_COUNT};

fn global(from: usize, to: usize, plug: usize) -> RecordData {
    RecordData::Global(MaximumGap::new(Edge::new(from, to).unwrap(), plug).unwrap())
}

fn json(data: &RecordData) -> String {
    serde_json::to_string(data).unwrap()
}

fn parse(text: &str) -> Result<RecordData, serde_json::Error> {
    serde_json::from_str(text)
}

/// A payload is accepted as a record only if it parses and serialises back to
/// the same bytes (the certificate's canonical-form rule).
fn canonical(text: &str) -> Option<RecordData> {
    let data = parse(text).ok()?;
    (json(&data) == text).then_some(data)
}

fn examples() -> Vec<(RecordData, &'static str)> {
    let [a, b] = geometry::edges()[0];
    let edge_json = Box::leak(format!(
        r#"{{"component":"Global","edge":[{b},{a}],"vertex":3}}"#
    ).into_boxed_str());
    vec![
        (
            RecordData::Domain { inequality: DomainInequality::new(13).unwrap() },
            r#"{"component":"Domain","inequality":13}"#,
        ),
        (
            RecordData::Exotic { cover: ExoticCover::CrossingPlus },
            r#"{"component":"Exotic","cover":"crossing+"}"#,
        ),
        (
            RecordData::Exotic { cover: ExoticCover::ArcMinus },
            r#"{"component":"Exotic","cover":"arc-"}"#,
        ),
        (RecordData::Local { cover: LocalCover(9) }, r#"{"component":"Local","cover":9}"#),
        (global(b, a, 3), edge_json),
    ]
}

#[test]
fn every_record_shape_has_its_documented_canonical_form() {
    for (data, text) in examples() {
        assert_eq!(json(&data), text);
        assert!(text.starts_with(&format!(r#"{{"component":"{}","#, data.component())));
        assert_eq!(canonical(text), Some(data.clone()));
    }
}

#[test]
fn the_design_examples_are_records() {
    let [a, b] = geometry::edges()[17];
    for text in [
        r#"{"component":"Domain","inequality":0}"#,
        r#"{"component":"Domain","inequality":72}"#,
        r#"{"component":"Local","cover":0}"#,
        r#"{"component":"Local","cover":4294967295}"#,
        &format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":59}}"#),
        &format!(r#"{{"component":"Global","edge":[{b},{a}],"vertex":0}}"#),
    ] {
        assert!(canonical(text).is_some(), "{text}");
    }
    for cover in ExoticCover::ALL {
        let text = format!(r#"{{"component":"Exotic","cover":"{cover}"}}"#);
        assert_eq!(canonical(&text), Some(RecordData::Exotic { cover }));
    }
}

#[test]
fn malformed_and_non_canonical_records_are_refused() {
    let [a, b] = geometry::edges()[0];
    let refused = [
        String::new(),
        "{}".into(),
        "null".into(),
        r#"{"component":"Domain"}"#.into(),
        r#"{"component":"Domain","inequality":73}"#.into(),
        r#"{"component":"Domain","inequality":-1}"#.into(),
        r#"{"component":"Domain","inequality":13.0}"#.into(),
        r#"{"component":"Domain","inequality":"13"}"#.into(),
        r#"{"component":"Domain","inequality":13,"inequality":13}"#.into(),
        r#"{"component":"Domain","inequality":13,"cover":9}"#.into(),
        r#"{"inequality":13,"component":"Domain"}"#.into(),
        r#"{"component":"Domain","inequality":013}"#.into(),
        r#"{"component":"Domain", "inequality":13}"#.into(),
        r#"{"component":"domain","inequality":13}"#.into(),
        r#"{"component":"LocalCover","cover":9}"#.into(),
        r#"{"component":"Local","cover":-9}"#.into(),
        r#"{"component":"Local","cover":9.0}"#.into(),
        r#"{"component":"Local","cover":"9"}"#.into(),
        r#"{"component":"Local","cover":4294967296}"#.into(),
        r#"{"component":"Local","cover":1e1}"#.into(),
        r#"{"component":"Exotic","cover":"Square"}"#.into(),
        r#"{"component":"Exotic","cover":"second_arc"}"#.into(),
        r#"{"component":"Exotic","cover":"second-arc"}"#.into(),
        r#"{"component":"Exotic","cover":"arc"}"#.into(),
        r#"{"component":"Exotic","cover":"arc±"}"#.into(),
        r#"{"component":"Exotic","cover":"crossing +"}"#.into(),
        r#"{"component":"Exotic","cover":"ArcPlus"}"#.into(),
        r#"{"component":"Exotic","cover":"secondArc"}"#.into(),
        r#"{"component":"Exotic","cover":0}"#.into(),
        r#"{"component":"Exotic","cover":"squ\u0061re"}"#.into(),
        r#"{"component":"Exotic","cover":"square","path":""}"#.into(),
        r#"{"component":"Global","edge":[0,1],"vertex":3}"#.into(),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":60}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{a}],"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b},{a}],"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a}],"vertex":3}}"#),
        format!(r#"{{"component":"Global","vertex":3,"edge":[{a},{b}]}}"#),
        format!(r#"{{"vertex":3,"edge":[{a},{b}],"component":"Global"}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}]}}"#),
        format!(r#"{{"component":"Global","vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":3,"sign":1}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":3,"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":3,"contact":0}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":null}}"#),
        format!(r#"{{"component":"Global","edge":null,"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":3.0}}"#),
        format!(r#"{{"component":"Global","edge":["{a}",{b}],"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":{{"from":{a},"to":{b}}},"vertex":3}}"#),
        format!(r#"{{"component":"Global","witness":{{"edge":[{a},{b}],"vertex":3}}}}"#),
        format!(r#"{{"component":"Global","inequality":13}}"#),
        format!(r#"{{"component":"Global","cover":9}}"#),
        r#"{"component":"Global","direction":[["1","0"],["0","0"]],"contact":0,"vertex":3}"#.into(),
        format!(r#"{{"component":"Global","inequality":13,"edge":[{a},{b}],"vertex":3}}"#),
        format!(r#"{{"component":"Domain","edge":[{a},{b}],"vertex":3}}"#),
        format!(r#"{{"component":"Global","edge":[{a},{b}],"vertex":3,"path":"0"}}"#),
    ];
    for text in &refused {
        assert_eq!(canonical(text), None, "{text}");
    }
    // A planar direction is a valid gap in a cover file but never Global's.
    let direction = r#"{"direction":[["1","0"],["0","0"]],"contact":0,"vertex":3}"#;
    assert!(serde_json::from_str::<Gap>(direction).is_ok());
    assert!(serde_json::from_str::<MaximumGap>(direction).is_err());
}

#[test]
fn global_witnesses_are_exactly_the_ordered_edges_with_a_vertex() {
    let mut edges = 0;
    for from in 0..VERTEX_COUNT {
        for to in 0..VERTEX_COUNT {
            let plug = (from + to) % 60;
            let text = format!(r#"{{"component":"Global","edge":[{from},{to}],"vertex":{plug}}}"#);
            let parsed = canonical(&text);
            assert_eq!(parsed.is_some(), geometry::is_edge(from, to), "{text}");
            if let Some(RecordData::Global(maximum)) = parsed {
                edges += 1;
                assert_eq!((maximum.edge().from(), maximum.edge().to(), maximum.plug()), (from, to, plug));
                assert_eq!(serde_json::to_string(&maximum).unwrap(), format!(r#"{{"edge":[{from},{to}],"vertex":{plug}}}"#));
            }
        }
    }
    assert_eq!(edges, 240);
    let error = serde_json::from_str::<MaximumGap>(r#"{"edge":[0,1],"vertex":3}"#).unwrap_err().to_string();
    assert!(error.contains(&WitnessError::NotAnEdge { from: 0, to: 1 }.to_string()), "{error}");
    let [a, b] = geometry::edges()[0];
    let error = serde_json::from_str::<MaximumGap>(&format!(r#"{{"edge":[{a},{b}],"vertex":60}}"#)).unwrap_err().to_string();
    assert!(error.contains(&WitnessError::VertexOutOfRange { index: 60 }.to_string()), "{error}");
}

#[test]
fn every_single_byte_change_of_a_record_is_refused_or_a_different_valid_record() {
    // The certificate refuses a changed line by its checksum; here the
    // payload alone: a changed payload never parses back to the same record,
    // and it is accepted only in canonical form.
    let bytes: Vec<u8> = (0x20..0x7f).collect();
    for (data, text) in examples() {
        for position in 0..text.len() {
            for &byte in &bytes {
                let mut changed = text.as_bytes().to_vec();
                if changed[position] == byte {
                    continue;
                }
                changed[position] = byte;
                let Ok(changed) = String::from_utf8(changed) else { continue };
                if let Some(other) = canonical(&changed) {
                    assert_ne!(other, data, "{changed}");
                    assert_eq!(json(&other), changed);
                }
            }
        }
    }
}

#[test]
fn names_and_orders_are_fixed() {
    assert_eq!(
        ExoticCover::ALL.map(ExoticCover::name),
        ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]
    );
    let mut sorted = ExoticCover::ALL;
    sorted.sort();
    assert_eq!(sorted, ExoticCover::ALL, "key order is the documented order");
    for cover in ExoticCover::ALL {
        assert_eq!(ExoticCover::from_name(cover.name()), Some(cover));
        assert_eq!(cover.to_string(), cover.name());
    }
    for name in ["Square", "second_arc", "", "square ", "local", "0"] {
        assert_eq!(ExoticCover::from_name(name), None, "{name}");
    }
    for (name, number) in [("0", 0), ("9", 9), ("10", 10), ("4294967295", u32::MAX)] {
        assert_eq!(LocalCover::from_name(name), Some(LocalCover(number)));
        assert_eq!(LocalCover(number).to_string(), name);
    }
    for name in ["", "09", "00", "-1", "+1", " 1", "1 ", "1.0", "4294967296", "square", "٣"] {
        assert_eq!(LocalCover::from_name(name), None, "{name}");
    }
    assert_eq!(
        ComponentName::ORDER,
        [ComponentName::Domain, ComponentName::Exotic, ComponentName::Local, ComponentName::Global]
    );
    let names = examples().into_iter().map(|(data, _)| data.component().to_string());
    assert_eq!(names.collect::<Vec<_>>(), ["Domain", "Exotic", "Exotic", "Local", "Global"]);
}

mod adversarial;
