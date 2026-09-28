use super::*;
use rid::arithmetic::exact::frac;

#[test]
fn colours_are_six_hex_digits() {
    assert_eq!(Color::parse("#FF8000").unwrap().hex(), "#ff8000");
    assert_eq!(Color::parse("#0a0B0c").unwrap(), Color([10, 11, 12]));
    for bad in ["ff8000", "#ff800", "#ff80000", "#gg8000", "#ff 800", "", "#", "#+f8000", "red"] {
        assert_eq!(Color::parse(bad), Err(StyleError::Color(bad.to_owned())), "{bad}");
    }
}

#[test]
fn styles_are_validated() {
    let c = Color::parse("#123456").unwrap();
    let stroke = || Stroke::new(c, frac(1, 10)).unwrap();
    assert!(Style::new(Some(c), Some(stroke()), Q::one()).is_ok());
    assert!(Style::new(None, Some(stroke()), frac(1, 5)).is_ok());
    assert!(Style::new(Some(c), None, frac(1, 5)).is_ok());
    assert_eq!(Style::new(None, None, Q::one()), Err(StyleError::Invisible));
    assert_eq!(Style::new(Some(c), None, Q::zero()), Err(StyleError::Opacity(Q::zero())));
    assert_eq!(Style::new(Some(c), None, frac(3, 2)), Err(StyleError::Opacity(frac(3, 2))));
    assert_eq!(Stroke::new(c, Q::zero()), Err(StyleError::Width(Q::zero())));
    assert_eq!(Stroke::new(c, frac(-1, 10)), Err(StyleError::Width(frac(-1, 10))));
}

#[test]
fn styles_read_from_json_need_every_field() {
    let parse = |text: &str| serde_json::from_str::<Style>(text);
    let style = parse(r##"{"fill": null, "stroke": {"color": "#ff8000", "width_mm": "7/50"}, "opacity": "1"}"##).unwrap();
    assert_eq!(style.fill(), None);
    assert_eq!(style.stroke().unwrap().width_mm(), &frac(7, 50));
    assert_eq!(style.opacity(), &Q::one());
    for bad in [
        r##"{"stroke": null, "opacity": "1", "fill": "#000000", "extra": 1}"##,
        r##"{"stroke": null, "opacity": "1"}"##,
        r##"{"fill": "#000000", "opacity": "1"}"##,
        r##"{"fill": "#000000", "stroke": null}"##,
        r##"{"fill": "#000000", "stroke": null, "opacity": 1}"##,
        r##"{"fill": "#000000", "stroke": null, "opacity": "0.5"}"##,
        r##"{"fill": null, "stroke": null, "opacity": "1"}"##,
        r##"{"fill": null, "stroke": {"color": "#000000"}, "opacity": "1"}"##,
        r##"{"fill": null, "stroke": {"color": "#000000", "width_mm": "0"}, "opacity": "1"}"##,
        r##"{"fill": null, "stroke": {"color": "#000000", "width_mm": "1", "dash": []}, "opacity": "1"}"##,
    ] {
        assert!(parse(bad).is_err(), "{bad}");
    }
}

fn numbered() -> PerOutcome<usize> {
    let covers = r#"{"square": 3, "pentagon": 4, "arc+": 5, "arc-": 6, "endpoint+": 7, "endpoint-": 8, "crossing+": 9, "crossing-": 10}"#;
    PerOutcome {
        domain: 0,
        global: 1,
        local: 2,
        exotic: serde_json::from_str::<PerCover<usize>>(covers).unwrap(),
        unresolved: 11,
        outside: 12,
    }
}

#[test]
fn every_outcome_has_its_own_entry() {
    let styles = numbered();
    for (k, outcome) in outcomes().iter().enumerate() {
        assert_eq!(*styles.get(outcome), k);
    }
    let mut distinct = outcomes().to_vec();
    distinct.dedup();
    assert_eq!(distinct.len(), 13);
}

#[test]
fn outcome_styles_read_from_json_need_every_entry() {
    let full = r#"{"Domain": 0, "Global": 1, "Local": 2,
        "Exotic": {"square": 3, "pentagon": 4, "arc+": 5, "arc-": 6, "endpoint+": 7, "endpoint-": 8, "crossing+": 9, "crossing-": 10},
        "unresolved": 11, "outside": 12}"#;
    assert_eq!(serde_json::from_str::<PerOutcome<usize>>(full).unwrap(), numbered());
    let missing = full.replace(r#""outside": 12"#, r#""unused": 12"#);
    assert!(serde_json::from_str::<PerOutcome<usize>>(&missing).is_err());
    let no_cover = full.replace(r#""crossing-": 10"#, r#""cross": 10"#);
    assert!(serde_json::from_str::<PerOutcome<usize>>(&no_cover).is_err());
    let dropped = full.replace(r#", "crossing-": 10"#, "");
    assert!(serde_json::from_str::<PerOutcome<usize>>(&dropped).is_err());
    let repeated = full.replace(r#""crossing-": 10"#, r#""crossing+": 10"#);
    assert!(serde_json::from_str::<PerOutcome<usize>>(&repeated).is_err());
    let old = full.replace(r#""arc-": 6"#, r#""second-arc": 6"#);
    assert!(serde_json::from_str::<PerOutcome<usize>>(&old).is_err());
    let lower = full.replace(r#""Domain""#, r#""domain""#);
    assert!(serde_json::from_str::<PerOutcome<usize>>(&lower).is_err());
}
