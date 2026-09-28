use super::standin::{boxed, Data, Fault, StandIn};
use super::*;
use rid::arithmetic::exact::{frac, q};
use serde_json::json;

fn label(component: &str, cover: Option<&str>) -> Label {
    Label {
        component: component.into(),
        cover: cover.map(String::from),
    }
}

/// A box of side `2^-k` at the given corner.
fn cube(corner: [rid::arithmetic::exact::Q; 5], k: u32) -> ConfigurationBox {
    let side = crate::points::neighbourhood::radius(k);
    boxed(corner.map(|c| (c.clone(), &c + &side)))
}

#[test]
fn labels_are_read_from_record_objects() {
    assert_eq!(Label::of(&json!({"component": "Global", "witness": {"edge": [0, 2], "vertex": 3}})).unwrap(), label("Global", None));
    assert_eq!(Label::of(&json!({"component": "Local", "cover": 9})).unwrap(), label("Local", Some("9")));
    assert_eq!(Label::of(&json!({"component": "Exotic", "cover": "arc-"})).unwrap(), label("Exotic", Some("arc-")));
    for bad in [
        json!(null),
        json!("Local"),
        json!([{"component": "Local"}]),
        json!({}),
        json!({"component": ""}),
        json!({"component": 3}),
        json!({"component": "Local", "cover": null}),
        json!({"component": "Local", "cover": -1}),
        json!({"component": "Local", "cover": 1.5}),
        json!({"component": "Local", "cover": ""}),
        json!({"component": "Local", "cover": true}),
        json!({"component": "Local", "cover": {"n": 1}}),
    ] {
        assert!(matches!(Label::of(&bad), Err(ProbeError::Label(_))), "{bad}");
    }
}

#[test]
fn labels_meet_targets_by_component_and_named_cover() {
    let local9 = label("Local", Some("9"));
    assert!(local9.meets(&label("Local", None)));
    assert!(local9.meets(&label("Local", Some("9"))));
    assert!(!local9.meets(&label("Local", Some("10"))));
    assert!(!local9.meets(&label("Exotic", None)));
    assert!(!label("Local", None).meets(&local9));
    assert_eq!(local9.to_string(), "Local 9");
    assert_eq!(label("Global", None).to_string(), "Global");
    // The transcript spelling omits an absent cover and refuses unknown fields.
    assert_eq!(serde_json::to_string(&label("Global", None)).unwrap(), r#"{"component":"Global"}"#);
    assert!(serde_json::from_str::<Label>(r#"{"component":"Global","x":1}"#).is_err());
}

#[test]
fn contradictions_follow_the_components_claims() {
    use Relation::*;
    let (domain, global, local) = (label("Domain", None), label("Global", None), label("Local", Some("1")));
    for (l, in_d, relation, expected) in [
        (&domain, true, Poke, Err(Contradiction::DomainInD)),
        (&domain, false, Touch, Ok(())),
        (&domain, false, Fit, Ok(())),
        (&global, false, Touch, Err(Contradiction::GlobalAtContained)),
        (&global, true, Touch, Err(Contradiction::GlobalAtContained)),
        (&global, false, Fit, Err(Contradiction::GlobalAtContained)),
        (&global, true, Poke, Ok(())),
        (&local, true, Touch, Ok(())),
        (&local, true, Fit, Err(Contradiction::FitEliminated)),
        (&local, false, Fit, Ok(())),
        (&domain, true, Fit, Err(Contradiction::FitEliminated)),
    ] {
        assert_eq!(consistent(l, in_d, relation), expected, "{l} {in_d} {relation:?}");
    }
}

#[test]
fn attempts_and_decisions_are_verified_by_their_components() {
    let standin = StandIn::typical();
    let names = standin.components();
    // Inside the square cover only: Exotic, first in order after Domain.
    let b = cube([q(0), q(0), q(0), q(0), q(0)], 8);
    let found = decide(&standin, &names, &b).unwrap().unwrap();
    assert_eq!(found.label, label("Exotic", Some("square")));
    assert_eq!(found.record, json!({"component": "Exotic", "cover": "square"}));
    assert!(attempt(&standin, &names, 0, &b).unwrap().is_none());
    assert_eq!(attempt(&standin, &names, 1, &b).unwrap().unwrap(), found);
    assert!(attempt(&standin, &names, 3, &b).unwrap().is_none());
    // A lying attempt is refused by its own verification.
    for index in 0..4 {
        let liar = StandIn {
            fault: Fault::LyingAttempt(index),
            ..StandIn::typical()
        };
        let far = cube([frac(1, 5), frac(1, 10), frac(1, 20), frac(1, 20), frac(1, 20)], 10);
        assert!(matches!(attempt(&liar, &names, index, &far), Err(ProbeError::Refused { .. })), "{index}");
    }
}

/// Probes whose decision and attempts return records of other components.
struct Confused;

impl Probes for Confused {
    type Data = Data;
    fn components(&self) -> Vec<String> {
        vec!["Domain".into(), "Local".into()]
    }
    fn decide(&self, _: &ConfigurationBox) -> Result<Option<Data>, BoxError> {
        Ok(Some(Data::Global { witness: 2 }))
    }
    fn attempt(&self, _: usize, _: &ConfigurationBox) -> Option<Data> {
        Some(Data::Local { cover: 0 })
    }
    fn covers(&self, index: usize) -> Vec<Data> {
        if index == 1 {
            vec![Data::Exotic { cover: "square".into() }]
        } else {
            Vec::new()
        }
    }
    fn verify(&self, _: usize, _: &ConfigurationBox, _: &Data) -> Result<(), BoxError> {
        Ok(())
    }
}

/// Probes whose collection reports a defect.
struct Broken;

impl Probes for Broken {
    type Data = Data;
    fn components(&self) -> Vec<String> {
        vec!["Domain".into()]
    }
    fn decide(&self, _: &ConfigurationBox) -> Result<Option<Data>, BoxError> {
        Err("Domain refuses its own record".into())
    }
    fn attempt(&self, _: usize, _: &ConfigurationBox) -> Option<Data> {
        None
    }
    fn covers(&self, _: usize) -> Vec<Data> {
        Vec::new()
    }
    fn verify(&self, _: usize, _: &ConfigurationBox, _: &Data) -> Result<(), BoxError> {
        Err("never".into())
    }
}

/// Probes with repeated component names.
struct Twice;

impl Probes for Twice {
    type Data = Data;
    fn components(&self) -> Vec<String> {
        vec!["Local".into(), "Local".into()]
    }
    fn decide(&self, _: &ConfigurationBox) -> Result<Option<Data>, BoxError> {
        Ok(None)
    }
    fn attempt(&self, _: usize, _: &ConfigurationBox) -> Option<Data> {
        None
    }
    fn covers(&self, _: usize) -> Vec<Data> {
        Vec::new()
    }
    fn verify(&self, _: usize, _: &ConfigurationBox, _: &Data) -> Result<(), BoxError> {
        Ok(())
    }
}

#[test]
fn component_names_must_be_distinct() {
    assert!(matches!(names(&Twice), Err(ProbeError::Names)));
}

#[test]
fn records_of_the_wrong_component_are_refused() {
    let b = ConfigurationBox::root();
    assert!(matches!(decide(&Broken, &Broken.components(), &b), Err(ProbeError::Decision(_))));
    let names = Confused.components();
    assert!(matches!(decide(&Confused, &names, &b), Err(ProbeError::UnknownComponent(c)) if c == "Global"));
    assert!(matches!(attempt(&Confused, &names, 0, &b), Err(ProbeError::Foreign { .. })));
    assert!(attempt(&Confused, &names, 1, &b).unwrap().is_some());
    assert!(matches!(probes(&Confused), Err(ProbeError::Foreign { .. })));
    // Only the decided names in their order are accepted.
    assert!(matches!(super::names(&Confused), Err(ProbeError::Names)));
    assert_eq!(super::names(&StandIn::typical()).unwrap(), ORDER);
}

#[test]
fn probes_list_covers_and_attempts_in_the_collections_order() {
    let standin = StandIn::typical();
    let names = standin.components();
    let list = probes(&standin).unwrap();
    let labels: Vec<String> = list.iter().map(|p| p.label.to_string()).collect();
    assert_eq!(
        labels,
        ["Domain", "Exotic square", "Exotic pentagon", "Local 0", "Local 1", "Global"]
    );
    // A box in the overlap of both Local covers: both cover probes succeed,
    // and the attempt of Local would name only the first.
    let b = cube([frac(1, 4) + frac(1, 64), frac(1, 8), frac(-1, 256), frac(-1, 256), frac(-1, 256)], 9);
    let asked: Vec<bool> = list
        .iter()
        .map(|p| ask(&standin, &names, p, &b).unwrap().is_some())
        .collect();
    assert_eq!(asked, [false, false, false, true, true, false]);
    let found = ask(&standin, &names, &list[4], &b).unwrap().unwrap();
    assert_eq!(found.record, json!({"component": "Local", "cover": 1}));
    assert_eq!(found.label, label("Local", Some("1")));
}
