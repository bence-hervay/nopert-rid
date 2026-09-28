use super::*;
use serde_json::{json, Map, Value};

fn search() -> Value {
    json!({
        "certificate": "runs/full.cert",
        "root": "",
        "threads": 28,
        "window": 4096,
        "max_depth": 4096,
        "depth_limit": null,
        "max_decisions": null
    })
}

fn check() -> Value {
    json!({
        "certificate": "runs/full.cert",
        "root": "0101",
        "max_depth": 4096,
        "threads": 3
    })
}

fn prepare() -> Value {
    json!({ "threads": 2 })
}

fn all() -> [(&'static str, Value); 3] {
    [("prepare", prepare()), ("search", search()), ("check", check())]
}

fn parse(name: &str, value: &Value) -> Result<Command, Error> {
    Command::parse(name, value.to_string().as_bytes())
}

fn fields(value: &Value) -> &Map<String, Value> {
    value.as_object().unwrap()
}

#[test]
fn the_documented_configurations_parse() {
    let nz = |n| NonZeroUsize::new(n).unwrap();
    assert_eq!(
        parse("search", &search()).unwrap(),
        Command::Search(Search {
            certificate: PathBuf::from("runs/full.cert"),
            root: String::new(),
            threads: nz(28),
            window: nz(4096),
            max_depth: 4096,
            depth_limit: None,
            max_decisions: None,
        })
    );
    let mut limited = search();
    limited["max_decisions"] = json!(15);
    limited["depth_limit"] = json!(40);
    let Command::Search(config) = parse("search", &limited).unwrap() else {
        panic!("not a search")
    };
    assert_eq!((config.max_decisions, config.depth_limit), (Some(15), Some(40)));
    assert_eq!(
        parse("check", &check()).unwrap(),
        Command::Check(Check {
            certificate: PathBuf::from("runs/full.cert"),
            root: "0101".into(),
            max_depth: 4096,
            threads: nz(3),
        })
    );
    assert_eq!(
        parse("prepare", &prepare()).unwrap(),
        Command::Prepare(Prepare { threads: nz(2) })
    );
}

#[test]
fn configurations_round_trip_through_serialisation() {
    for (name, value) in all() {
        let command = parse(name, &value).unwrap();
        let text = match &command {
            Command::Prepare(c) => serde_json::to_vec(c),
            Command::Search(c) => serde_json::to_vec(c),
            Command::Check(c) => serde_json::to_vec(c),
        }
        .unwrap();
        assert_eq!(Command::parse(name, &text).unwrap(), command);
        assert_eq!(serde_json::from_slice::<Value>(&text).unwrap(), value);
    }
}

#[test]
fn every_missing_field_is_refused() {
    for (name, value) in all() {
        for field in fields(&value).keys() {
            let mut missing = value.clone();
            missing.as_object_mut().unwrap().remove(field);
            let error = parse(name, &missing).unwrap_err();
            assert!(
                error.to_string().contains(&format!("missing field `{field}`")),
                "{name}.{field}: {error}"
            );
        }
    }
}

#[test]
fn unknown_and_duplicate_fields_are_refused() {
    for (name, value) in all() {
        for extra in ["extra", "Threads", "max-depth", ""] {
            let mut more = value.clone();
            more[extra] = json!(1);
            assert!(parse(name, &more).is_err(), "{name} with {extra:?}");
        }
        let first = fields(&value).keys().next().unwrap();
        let text = value.to_string();
        let duplicated = format!(
            "{{\"{first}\":{},{}",
            fields(&value)[first],
            &text[1..]
        );
        assert!(Command::parse(name, duplicated.as_bytes()).is_err(), "{duplicated}");
    }
}

#[test]
fn mistyped_values_are_refused() {
    let wrong: [(&str, Vec<Value>); 7] = [
        ("certificate", vec![json!(null), json!(1), json!(["a"]), json!({})]),
        ("root", vec![json!(null), json!(0), json!(false)]),
        (
            "threads",
            vec![
                json!(0),
                json!(-1),
                json!(2.0),
                json!(2.5),
                json!("2"),
                json!(null),
                json!(1e3),
                json!(1e40),
            ],
        ),
        ("window", vec![json!(0), json!(-4), json!("4096"), json!(null), json!(true)]),
        ("max_depth", vec![json!(-1), json!(4096.0), json!("4096"), json!(null)]),
        ("depth_limit", vec![json!(-1), json!(40.0), json!("40"), json!(false), json!([])]),
        ("max_decisions", vec![json!(-1), json!(1.5), json!("10"), json!(false), json!([])]),
    ];
    for (name, value) in all() {
        for (field, values) in &wrong {
            if !fields(&value).contains_key(*field) {
                continue;
            }
            for bad in values {
                let mut changed = value.clone();
                changed[*field] = bad.clone();
                assert!(parse(name, &changed).is_err(), "{name}.{field} = {bad}");
            }
        }
    }
    // Explicit null is accepted only for the absent limits.
    let mut value = search();
    value["max_decisions"] = json!(null);
    value["depth_limit"] = json!(null);
    assert!(parse("search", &value).is_ok());
    // Integers beyond the field's range are refused rather than wrapped.
    value["max_decisions"] = json!(u64::MAX);
    assert!(parse("search", &value).is_ok());
    let text = value.to_string().replace(&u64::MAX.to_string(), "18446744073709551616");
    assert!(Command::parse("search", text.as_bytes()).is_err());
}

#[test]
fn malformed_documents_are_refused() {
    for text in [
        "",
        "null",
        "[]",
        "{}",
        "\"search\"",
        r#"{"threads": 2} {"threads": 2}"#,
        r#"{"threads": 2},"#,
        r#"{"threads": 2"#,
        r#"{threads: 2}"#,
        r#"{"threads": 2,}"#,
        r#"{'threads': 2}"#,
    ] {
        let parsed = Command::parse("prepare", text.as_bytes());
        assert!(matches!(parsed, Err(Error::Json { .. })), "{text}");
    }
    // Surrounding whitespace, including a final newline, is ordinary JSON.
    assert!(Command::parse("prepare", b" {\"threads\": 2}\n").is_ok());
}

#[test]
fn configurations_of_one_command_do_not_parse_as_another() {
    for (name, value) in all() {
        for (other, _) in all() {
            if other != name {
                assert!(parse(other, &value).is_err(), "{name} as {other}");
            }
        }
    }
}

#[test]
fn unknown_commands_and_unreadable_files_are_refused() {
    for name in ["", "Search", "run", "search ", "help"] {
        assert!(matches!(
            Command::parse(name, search().to_string().as_bytes()),
            Err(Error::UnknownCommand(_))
        ));
        // An unknown command is reported before the file is read.
        assert!(matches!(
            Command::load(name, Path::new("/nonexistent/config.json")),
            Err(Error::UnknownCommand(_))
        ));
    }
    assert!(matches!(
        Command::load("search", Path::new("/nonexistent/config.json")),
        Err(Error::Read { .. })
    ));
    let file = std::env::temp_dir().join(format!("rid-config-test-{}.json", std::process::id()));
    std::fs::write(&file, search().to_string()).unwrap();
    let loaded = Command::load("search", &file);
    std::fs::remove_file(&file).unwrap();
    assert_eq!(loaded.unwrap(), parse("search", &search()).unwrap());
}

#[test]
fn the_example_configurations_parse_as_their_commands() {
    let examples: [(&str, &[u8]); 4] = [
        ("prepare", include_bytes!("../../../examples/prepare.json")),
        ("search", include_bytes!("../../../examples/pilot.json")),
        ("search", include_bytes!("../../../examples/search.json")),
        ("check", include_bytes!("../../../examples/check.json")),
    ];
    for (name, bytes) in examples {
        let command = Command::parse(name, bytes).unwrap();
        match (name, command) {
            ("prepare", Command::Prepare(_)) | ("search", Command::Search(_)) | ("check", Command::Check(_)) => {}
            (name, other) => panic!("{name}: {other:?}"),
        }
    }
}
