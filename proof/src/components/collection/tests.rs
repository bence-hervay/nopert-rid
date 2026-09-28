use super::*;
use crate::components::record::{ExoticCover, LocalCover};
use crate::components::tests::fixtures::{bare, from_path, search, BoxCover, Random, SMALL_ROOTS};
use crate::problem::configuration::ConfigurationBox;

fn collection(covers: &[(&str, &str)]) -> Result<Collection, CoversError> {
    let covers = covers.iter().map(|(name, path)| BoxCover::boxed(name, path)).collect();
    Collection::from_covers(covers)
}

/// The collection's check, component by component in the documented order.
fn first_witness(collection: &Collection, b: &ConfigurationBox) -> Option<RecordData> {
    if let Some(inequality) = collection.domain.check(b) {
        return Some(RecordData::Domain { inequality });
    }
    if let Some(cover) = collection.exotic.check(b) {
        return Some(RecordData::Exotic { cover });
    }
    if let Some(cover) = collection.local.check(b) {
        return Some(RecordData::Local { cover });
    }
    collection.global.check(b).map(RecordData::Global)
}

#[test]
fn the_first_components_witness_is_the_record() {
    // Below SMALL_ROOTS[0], find a box outside D and a box only Global eliminates.
    let (records, _) = search(&bare(), SMALL_ROOTS[0], 200).unwrap();
    let find = |name| records.iter().find(|(_, d)| d.component() == name).unwrap().0.clone();
    let outside = find(ComponentName::Domain);
    let global = find(ComponentName::Global);
    // Covers containing both boxes, of both kinds.
    let collection = collection(&[("7", &outside), ("3", &global), ("square", &global), ("arc+", &outside)]).unwrap();
    let check = |path: &str| collection.check(&from_path(path)).unwrap();
    assert_eq!(check(&outside).component(), ComponentName::Domain);
    assert_eq!(check(&global), RecordData::Exotic { cover: ExoticCover::Square });
    let collection = self::collection(&[("7", &outside), ("3", &global)]).unwrap();
    assert_eq!(collection.check(&from_path(&global)), Some(RecordData::Local { cover: LocalCover(3) }));
    assert_eq!(bare().check(&from_path(&global)).unwrap().component(), ComponentName::Global);
    // Every record of the bare search also holds in the collection with covers.
    for (path, data) in &records {
        let b = from_path(path);
        collection.holds(&b, data).unwrap();
        assert!(collection.check(&b).is_some());
    }
    assert_eq!(bare().check(&ConfigurationBox::root()), None);
}

#[test]
fn a_record_holds_by_the_component_it_names() {
    let collection = collection(&[("square", "0110"), ("4", "10")]).unwrap();
    let b = from_path("01101");
    collection.holds(&b, &RecordData::Exotic { cover: ExoticCover::Square }).unwrap();
    let refusal = collection.holds(&b, &RecordData::Local { cover: LocalCover(4) }).unwrap_err();
    assert!(matches!(refusal, Refusal::NotContained(ref cover) if cover == "4"));
    let refusal = collection.holds(&b, &RecordData::Exotic { cover: ExoticCover::ArcPlus }).unwrap_err();
    assert!(matches!(refusal, Refusal::UnknownCover(ref cover) if cover == "arc+"));
    let refusal = collection.holds(&b, &RecordData::Local { cover: LocalCover(5) }).unwrap_err();
    assert!(matches!(refusal, Refusal::UnknownCover(ref cover) if cover == "5"));
    // The box of an Exotic record, with the same cover named as Local
    // by number, is refused: the cover's key is part of the record.
    let refusal = collection.holds(&b, &RecordData::Local { cover: LocalCover(0) }).unwrap_err();
    assert!(matches!(refusal, Refusal::UnknownCover(_)));
    let domain = RecordData::Domain { inequality: crate::elimination::witness::DomainInequality::new(0).unwrap() };
    assert!(matches!(collection.holds(&b, &domain), Err(Refusal::Lemma(_))));
}

#[test]
fn covers_are_assigned_by_their_names() {
    let collection = collection(&[("12", ""), ("crossing-", ""), ("0", ""), ("square", "")]).unwrap();
    assert_eq!(collection.exotic.keys(), [ExoticCover::Square, ExoticCover::CrossingMinus]);
    assert_eq!(collection.local.keys(), [LocalCover(0), LocalCover(12)]);
    let refused = |covers: &[(&str, &str)]| self::collection(covers).err().unwrap();
    for name in ["Square", "01", "local-3", "", "second_arc"] {
        assert_eq!(refused(&[(name, "")]), CoversError::Foreign { component: ComponentName::Local, name: name.into() });
    }
    assert!(matches!(refused(&[("arc+", ""), ("arc+", "0")]), CoversError::Duplicate { component: ComponentName::Exotic, .. }));
    assert!(matches!(refused(&[("3", ""), ("3", "0")]), CoversError::Duplicate { component: ComponentName::Local, .. }));
}

/// On random boxes, the collection's check is the first component's own
/// check, and it holds; when it is `None`, no component's witness holds.
#[test]
fn check_is_the_first_components_witness_and_holds_on_random_boxes() {
    let collection = collection(&[("3", "0111"), ("endpoint+", "1010"), ("11", "00")]).unwrap();
    let mut random = Random(51);
    let mut kinds = std::collections::BTreeSet::new();
    for _ in 0..60 {
        let depth = random.below(30);
        let b = from_path(&random.path(depth));
        let checked = collection.check(&b);
        assert_eq!(checked, first_witness(&collection, &b));
        match &checked {
            Some(data) => {
                collection.holds(&b, data).unwrap();
                kinds.insert(data.component().to_string());
            }
            None => {
                assert!(collection.exotic.candidates(&b).iter().all(|k| collection.exotic.holds(&b, k).is_err()));
                assert!(collection.local.candidates(&b).iter().all(|k| collection.local.holds(&b, k).is_err()));
            }
        }
    }
    assert!(kinds.len() >= 3, "{kinds:?}");
}

#[test]
fn the_report_lists_the_order_and_the_covers() {
    let collection = collection(&[("12", ""), ("arc+", ""), ("0", "")]).unwrap();
    assert_eq!(
        collection.report(),
        serde_json::json!({
            "order": ["Domain", "Exotic", "Local", "Global"],
            "exotic_covers": ["arc+"],
            "local_covers": [0, 12],
        })
    );
}
