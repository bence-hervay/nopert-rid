use super::*;
use crate::testing::Scratch;

const COMPLETENESS: &str = r#"{
  "catalogue": "validation/catalogue/points.json",
  "transcript": "runs/completeness.jsonl",
  "threads": 28,
  "max_k": 256,
  "selection": null
}"#;

const IRREDUNDANCY: &str = r#"{
  "catalogue": "validation/catalogue/targets.json",
  "transcript": "runs/irredundancy.jsonl",
  "threads": 12,
  "scales": [12, 16, 32, 64, 128, 256, 512, 1024, 2048],
  "selection": null
}"#;

const PILOT: &str = r#"{
  "roots": {"sample": {"depth": 20, "count": 256, "seed": 1}},
  "threads": 28,
  "depth_limit": 200,
  "max_decisions": 100000,
  "transcript": "runs/pilot.jsonl"
}"#;

const TARGETS: &str = r#"{
  "catalogue": "validation/catalogue/targets.json"
}"#;

const CHECK: &str = r#"{
  "catalogue": "validation/catalogue/points.json",
  "transcript": "runs/completeness.jsonl",
  "threads": 28
}"#;

const COMPARE: &str = r#"{
  "catalogue": "validation/catalogue/points.json",
  "baseline": "runs/before.jsonl",
  "candidate": "runs/after.jsonl",
  "threads": 28,
  "worse_by": 1,
  "fail_on_regression": true
}"#;

const DOCUMENTED: [(&str, &str); 6] = [
    ("completeness", COMPLETENESS),
    ("irredundancy", IRREDUNDANCY),
    ("pilot", PILOT),
    ("targets", TARGETS),
    ("check", CHECK),
    ("compare", COMPARE),
];

#[test]
fn the_documented_configurations_parse_and_round_trip() {
    for (name, text) in DOCUMENTED {
        let command = Command::parse(name, text.as_bytes()).unwrap();
        let again = match &command {
            Command::Completeness(c) => serde_json::to_vec(c),
            Command::Irredundancy(c) => serde_json::to_vec(c),
            Command::Pilot(c) => serde_json::to_vec(c),
            Command::Targets(c) => serde_json::to_vec(c),
            Command::Check(c) => serde_json::to_vec(c),
            Command::Compare(c) => serde_json::to_vec(c),
        }
        .unwrap();
        assert_eq!(Command::parse(name, &again).unwrap(), command);
    }
    let paths = r#"{"roots": {"paths": ["0", "10"]}, "threads": 1, "depth_limit": 5, "max_decisions": null, "transcript": "t"}"#;
    match Command::parse("pilot", paths.as_bytes()).unwrap() {
        Command::Pilot(p) => {
            assert_eq!(p.roots, Roots::Paths(vec!["0".into(), "10".into()]));
            assert_eq!(p.max_decisions, None);
        }
        other => panic!("{other:?}"),
    }
}

/// The documented configuration as an object, to remove or change fields.
fn object(text: &str) -> serde_json::Map<String, serde_json::Value> {
    serde_json::from_str(text).unwrap()
}

#[test]
fn every_field_is_required_and_no_other_is_accepted() {
    for (name, text) in DOCUMENTED {
        let fields = object(text);
        for field in fields.keys() {
            let mut missing = fields.clone();
            missing.remove(field);
            let bytes = serde_json::to_vec(&missing).unwrap();
            assert!(Command::parse(name, &bytes).is_err(), "{name} without {field}");
        }
        let mut extra = fields.clone();
        extra.insert("default".into(), serde_json::Value::Null);
        assert!(Command::parse(name, &serde_json::to_vec(&extra).unwrap()).is_err());
        // Repeated fields.
        let repeated = text.replacen('{', r#"{"threads": 1, "#, 1);
        assert!(Command::parse(name, repeated.as_bytes()).is_err(), "{name} repeated");
        // A configuration for another command.
        for (other, _) in DOCUMENTED {
            if other != name && !(name == "check" && other == "check") {
                let foreign = DOCUMENTED.iter().find(|(n, _)| *n == other).unwrap().1;
                assert!(Command::parse(name, foreign.as_bytes()).is_err(), "{other} as {name}");
            }
        }
    }
}

#[test]
fn values_of_the_wrong_type_are_refused() {
    let with = |text: &str, field: &str, value: serde_json::Value| {
        let mut fields = object(text);
        fields.insert(field.into(), value);
        serde_json::to_vec(&fields).unwrap()
    };
    use serde_json::json;
    let cases = [
        ("completeness", with(COMPLETENESS, "threads", json!(0))),
        ("completeness", with(COMPLETENESS, "threads", json!(-1))),
        ("completeness", with(COMPLETENESS, "threads", json!(2.0))),
        ("completeness", with(COMPLETENESS, "threads", json!("2"))),
        ("completeness", with(COMPLETENESS, "max_k", json!(null))),
        ("completeness", with(COMPLETENESS, "max_k", json!(4294967296u64))),
        ("completeness", with(COMPLETENESS, "selection", json!("a"))),
        ("completeness", with(COMPLETENESS, "catalogue", json!(null))),
        ("irredundancy", with(IRREDUNDANCY, "scales", json!(null))),
        ("irredundancy", with(IRREDUNDANCY, "scales", json!([-1]))),
        ("pilot", with(PILOT, "roots", json!({"paths": ["0"], "sample": {"depth": 1, "count": 1, "seed": 0}}))),
        ("pilot", with(PILOT, "roots", json!({}))),
        ("pilot", with(PILOT, "roots", json!({"list": ["0"]}))),
        ("pilot", with(PILOT, "roots", json!({"sample": {"depth": 1, "count": 1}}))),
        ("pilot", with(PILOT, "roots", json!({"sample": {"depth": 1, "count": 1, "seed": 0, "x": 1}}))),
        ("pilot", with(PILOT, "depth_limit", json!(null))),
        ("compare", with(COMPARE, "fail_on_regression", json!(1))),
        ("compare", with(COMPARE, "worse_by", json!(null))),
    ];
    for (name, bytes) in cases {
        assert!(Command::parse(name, &bytes).is_err(), "{name}: {}", String::from_utf8_lossy(&bytes));
    }
    assert!(Command::parse("completeness", b"").is_err());
    assert!(Command::parse("completeness", b"[]").is_err());
}

#[test]
fn unknown_commands_are_refused_before_reading_and_unreadable_files_after() {
    let dir = Scratch::new("config");
    let missing = dir.join("missing.json");
    assert!(matches!(Command::load("search", &missing), Err(Error::UnknownCommand(_))));
    assert!(matches!(Command::load("check", &missing), Err(Error::Read { .. })));
    let file = dir.join("check.json");
    std::fs::write(&file, CHECK).unwrap();
    assert!(matches!(Command::load("check", &file), Ok(Command::Check(_))));
    assert!(matches!(Command::parse("Check", CHECK.as_bytes()), Err(Error::UnknownCommand(_))));
}
