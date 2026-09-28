use super::*;
use crate::components::tests::fixtures::{from_path, two_to_minus, BoxCover, Random};
use crate::problem::configuration::AXES;

fn exotic(covers: &[(&str, &str)]) -> Result<Exotic, CoversError> {
    Exotic::new(covers.iter().map(|(name, path)| BoxCover::boxed(name, path)).collect())
}

fn local(covers: &[(&str, &str)]) -> Result<Local, CoversError> {
    Local::new(covers.iter().map(|(name, path)| BoxCover::boxed(name, path)).collect())
}

#[test]
fn only_the_components_own_names_are_accepted_once() {
    let foreign = |component, name: &str| CoversError::Foreign { component, name: name.into() };
    let duplicate = |component, name: &str| CoversError::Duplicate { component, name: name.into() };
    let e = ComponentName::Exotic;
    let l = ComponentName::Local;
    for name in ["Square", "second_arc", "0", "", "local"] {
        assert_eq!(exotic(&[("arc+", ""), (name, "")]).err(), Some(foreign(e, name)));
    }
    for name in ["square", "01", "-1", "+1", "1.0", "", "4294967296"] {
        assert_eq!(local(&[("1", ""), (name, "")]).err(), Some(foreign(l, name)));
    }
    assert_eq!(exotic(&[("arc+", "0"), ("arc+", "1")]).err(), Some(duplicate(e, "arc+")));
    assert_eq!(local(&[("7", "0"), ("3", ""), ("7", "1")]).err(), Some(duplicate(l, "7")));
    let all: Vec<_> = ExoticCover::ALL.iter().rev().map(|r| (r.name(), "")).collect();
    let component = exotic(&all).unwrap();
    assert_eq!(component.keys(), ExoticCover::ALL);
    assert_eq!(component.candidates(&from_path("0110")), ExoticCover::ALL);
    let component = local(&[("10", ""), ("2", ""), ("0", "")]).unwrap();
    assert_eq!(component.keys(), [LocalCover(0), LocalCover(2), LocalCover(10)]);
    assert_eq!(component.candidates(&from_path("1")), [LocalCover(0), LocalCover(2), LocalCover(10)]);
    assert!(exotic(&[]).unwrap().keys().is_empty());
}

#[test]
fn the_first_containing_cover_in_key_order_is_the_witness() {
    // Nested covers: every one of them contains the deepest box.
    let component = local(&[("10", "0110"), ("2", "01"), ("33", "011")]).unwrap();
    let deep = from_path("01101");
    assert_eq!(component.check(&deep), Some(LocalCover(2)));
    assert_eq!(component.check(&from_path("0111")), Some(LocalCover(2)));
    assert_eq!(component.check(&from_path("0")), None);
    for number in [2, 10, 33] {
        assert!(component.holds(&deep, &LocalCover(number)).is_ok());
    }
    let component = exotic(&[("crossing-", "1"), ("pentagon", "10"), ("square", "0")]).unwrap();
    assert_eq!(component.check(&from_path("101")), Some(ExoticCover::Pentagon));
    assert_eq!(component.check(&from_path("11")), Some(ExoticCover::CrossingMinus));
    assert_eq!(component.check(&from_path("00")), Some(ExoticCover::Square));
    assert_eq!(component.check(&from_path("")), None);
}

/// Boxes on a cover's boundary: the cover itself, boxes sharing a face from
/// inside and outside, straddling boxes and points on and next to the faces.
#[test]
fn boxes_at_cover_boundaries_are_decided_by_closed_containment() {
    let path = "0110100";
    let component = local(&[("5", path)]).unwrap();
    let cover = from_path(path);
    let accepted = |b: &ConfigurationBox| {
        let verdict = component.holds(b, &LocalCover(5));
        assert_eq!(component.check(b).is_some(), verdict.is_ok());
        match verdict {
            Ok(()) => true,
            Err(Refusal::NotContained(cover)) => {
                assert_eq!(cover, "5");
                false
            }
            Err(other) => panic!("{other}"),
        }
    };
    assert!(accepted(&cover));
    for child in ["0", "1", "00", "11111", "0101010101"] {
        assert!(accepted(&from_path(&format!("{path}{child}"))));
    }
    for outside in ["011010", "0110101", "0110110", "", "1"] {
        assert!(!accepted(&from_path(outside)), "{outside}");
    }
    // Every corner of the cover, as a point, is inside; moved outward by a
    // tiny amount along any axis where it is on the upper face, it is not.
    let tiny = two_to_minus(300);
    for mask in 0..32u8 {
        let corner = cover.corner(mask);
        assert!(accepted(&ConfigurationBox::point(&corner)));
        for axis in 0..AXES {
            let mut moved = corner.clone();
            moved[axis] = if mask & (1 << axis) == 0 { &moved[axis] - &tiny } else { &moved[axis] + &tiny };
            assert!(!accepted(&ConfigurationBox::point(&moved)));
            // A box from the corner to the moved point straddles the face.
            let axes = std::array::from_fn(|j| {
                let (a, b) = (corner[j].clone(), moved[j].clone());
                crate::arithmetic::exact::Interval::new(a.clone().min(b.clone()), a.max(b)).unwrap()
            });
            assert!(!accepted(&ConfigurationBox::new(axes)));
        }
    }
}

#[test]
fn random_covers_agree_with_the_tree_structure() {
    // In the bisection tree, a box lies in another exactly when the other's
    // path is a prefix of its path.
    let mut random = Random(31);
    for _ in 0..30 {
        let covers: Vec<(String, String)> = (0..4)
            .map(|k| {
                let (number, depth) = (random.below(40) * 4 + k, random.below(4) + 1);
                (number.to_string(), random.path(depth))
            })
            .collect();
        let named: Vec<(&str, &str)> = covers.iter().map(|(n, p)| (n.as_str(), p.as_str())).collect();
        let component = local(&named).unwrap();
        for _ in 0..20 {
            let depth = random.below(8);
            let path = random.path(depth);
            let b = from_path(&path);
            let mut containing: Vec<u32> = covers
                .iter()
                .filter(|(_, p)| path.starts_with(p.as_str()))
                .map(|(n, _)| n.parse().unwrap())
                .collect();
            containing.sort();
            let expected = containing.first().map(|&n| LocalCover(n));
            assert_eq!(component.check(&b), expected, "{path} {covers:?}");
            for (name, _) in &covers {
                let number = name.parse().unwrap();
                let verdict = component.holds(&b, &LocalCover(number));
                assert_eq!(verdict.is_ok(), containing.contains(&number));
            }
        }
    }
}

#[test]
fn unknown_covers_are_refused() {
    let component = exotic(&[("square", "")]).unwrap();
    let b = from_path("0101");
    assert!(component.holds(&b, &ExoticCover::Square).is_ok());
    let unknown = component.holds(&b, &ExoticCover::ArcPlus);
    assert!(matches!(unknown, Err(Refusal::UnknownCover(ref cover)) if cover == "arc+"));
    let local_component = local(&[("0", "")]).unwrap();
    let unknown = local_component.holds(&b, &LocalCover(1));
    assert!(matches!(unknown, Err(Refusal::UnknownCover(ref cover)) if cover == "1"));
    let unknown = local_component.holds(&b, &LocalCover(u32::MAX));
    assert!(matches!(unknown, Err(Refusal::UnknownCover(ref cover)) if cover == "4294967295"));
}

/// A witness names one cover, and only that cover's containment counts: a
/// box lying only in another cover of the same component does not hold for
/// the named one, although the component's check finds the other cover.
#[test]
fn a_witness_naming_a_cover_does_not_hold_for_a_box_only_in_another_cover() {
    let component = local(&[("3", "0110"), ("7", "1001")]).unwrap();
    for (inside, other) in [("0110", 7), ("01101", 7), ("1001", 3), ("10011010", 3)] {
        let b = from_path(inside);
        assert!(component.check(&b).is_some());
        let refused = component.holds(&b, &LocalCover(other));
        assert!(matches!(&refused, Err(Refusal::NotContained(cover)) if *cover == other.to_string()), "{inside}: {refused:?}");
    }
    let component = exotic(&[("square", "00"), ("crossing-", "11")]).unwrap();
    let refused = component.holds(&from_path("001"), &ExoticCover::CrossingMinus);
    assert!(matches!(refused, Err(Refusal::NotContained(_))));
    component.holds(&from_path("001"), &ExoticCover::Square).unwrap();
}
