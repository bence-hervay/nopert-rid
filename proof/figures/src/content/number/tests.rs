//! Tests of the exact numbers of figure files: canonical rational strings
//! only, pairs as elements of Q(√5) and closed intervals, and required
//! optional fields.
use super::*;
use rid::arithmetic::exact::frac;
use serde_json::json;

fn rational(value: serde_json::Value) -> Result<Rational, serde_json::Error> {
    serde_json::from_value(value)
}

#[test]
fn rationals_are_canonical_strings() {
    for (text, expected) in [("3", frac(3, 1)), ("-7/12", frac(-7, 12)), ("0", frac(0, 1))] {
        assert_eq!(rational(json!(text)).unwrap(), Rational(expected.clone()));
        assert_eq!(serde_json::to_value(Rational(expected)).unwrap(), json!(text));
    }
    for bad in [json!("3/1"), json!("6/4"), json!("-0"), json!("+3"), json!("1/-2"), json!("0.5"), json!(" 3"), json!(""), json!(3), json!(0.5), json!(null)] {
        assert!(rational(bad.clone()).is_err(), "{bad}");
    }
}

#[test]
fn pairs_are_field_elements_and_closed_intervals() {
    let pair = |a: i64, b: i64| [Rational(frac(a, 1)), Rational(frac(b, 1))];
    assert_eq!(field(&pair(1, 2)), QSqrt5::new(frac(1, 1), frac(2, 1)));
    let point = interval(&pair(2, 2)).unwrap();
    assert_eq!((point.lo(), point.hi()), (&frac(2, 1), &frac(2, 1)));
    assert!(interval(&pair(1, 3)).is_ok());
    assert!(interval(&pair(3, 1)).is_err(), "a reversed interval is refused");
}

#[derive(Debug, serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct WithOption {
    #[serde(deserialize_with = "required")]
    value: Option<Rational>,
}

#[test]
fn optional_fields_are_required_and_only_null_is_none() {
    let parse = |v: serde_json::Value| serde_json::from_value::<WithOption>(v);
    assert!(parse(json!({"value": null})).unwrap().value.is_none());
    assert_eq!(parse(json!({"value": "1/2"})).unwrap().value, Some(Rational(frac(1, 2))));
    assert!(parse(json!({})).is_err(), "a missing field is not None");
}
