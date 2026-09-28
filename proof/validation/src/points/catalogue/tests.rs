use super::*;
use crate::testing::{shipped, Rng};
use num_bigint::BigInt;
use rid::arithmetic::exact::{q, Q};

fn point(id: &str, centre: &str) -> String {
    format!(r#"{{"id":"{id}","group":"g","centre":{centre},"expected":["Local"]}}"#)
}

fn catalogue(points: &[String]) -> Vec<u8> {
    format!(r#"{{"format":"{FORMAT}","points":[{}]}}"#, points.join(",")).into_bytes()
}

const ORIGIN: &str = r#"[["0","0"],["0","0"],["0","0"],["0","0"],["0","0"]]"#;

fn with_axis(axis: usize, pair: [&str; 2]) -> String {
    let mut spelling: Vec<String> = (0..AXES).map(|_| r#"["0","0"]"#.to_string()).collect();
    spelling[axis] = format!(r#"["{}","{}"]"#, pair[0], pair[1]);
    format!("[{}]", spelling.join(","))
}

#[test]
fn the_shipped_catalogue_has_683_points_with_the_85_original_points_first() {
    let catalogue = Catalogue::load(&shipped("points.json")).unwrap();
    assert_eq!(catalogue.points.len(), 683);
    assert_eq!(catalogue.points[0].id, "aligned-generic-rational");
    assert!(catalogue.points[..85].iter().all(|p| !p.group.starts_with("dense-")));
    assert!(catalogue.points[85..].iter().all(|p| p.group.starts_with("dense-")));
    let names = ["Domain", "Exotic", "Local", "Global"];
    for p in &catalogue.points {
        assert!(!p.expected.is_empty());
        assert!(p.expected.iter().all(|e| names.contains(&e.as_str())), "{}", p.id);
        assert_eq!(Centre::parse(&p.centre.spelling()).unwrap(), p.centre);
    }
    assert_eq!(catalogue.sha256, sha256_hex(&read(&shipped("points.json")).unwrap()));
}

#[test]
fn sha256_matches_the_published_vectors() {
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

#[test]
fn centres_round_trip_through_their_spelling() {
    let mut rng = Rng::new(1);
    let root = ConfigurationBox::root();
    for _ in 0..200 {
        let coordinates: [QSqrt5; AXES] = std::array::from_fn(|j| {
            let a = &root.axes()[j];
            QSqrt5::from_rational(rng.rational(a.lo(), a.hi(), 20))
        });
        let centre = Centre::new(coordinates).unwrap();
        assert_eq!(Centre::parse(&centre.spelling()).unwrap(), centre);
    }
}

#[test]
fn every_root_face_is_accepted_and_the_outside_refused_exactly() {
    let root = ConfigurationBox::root();
    // 161 - 72√5 ≈ 0.0031 > 0 and its negative: a Pell unit, decided exactly.
    let tiny = QSqrt5::new(q(161), q(-72));
    assert_eq!(tiny.sign(), std::cmp::Ordering::Greater);
    for axis in 0..AXES {
        for (end, outward) in [(root.axes()[axis].lo(), -1), (root.axes()[axis].hi(), 1)] {
            let at = |x: QSqrt5| -> [QSqrt5; AXES] {
                std::array::from_fn(|j| if j == axis { x.clone() } else { QSqrt5::zero() })
            };
            let face = QSqrt5::from_rational(end.clone());
            // The origin lies on the lower faces of s and t and inside the rest.
            assert!(Centre::new(at(face.clone())).is_ok());
            let scale = Q::new(BigInt::from(1), BigInt::from(1) << 200);
            let step = tiny.scale(&scale).scale(&q(outward));
            assert_eq!(
                Centre::new(at(&face + &step)),
                Err(CentreError::OutsideRoot { axis })
            );
            assert!(Centre::new(at(&face - &step)).is_ok());
        }
    }
}

#[test]
fn well_formed_catalogues_parse() {
    let bytes = catalogue(&[point("a", ORIGIN), point("b", &with_axis(0, ["-2", "1"]))]);
    let parsed = Catalogue::parse(&bytes).unwrap();
    assert_eq!(parsed.points.len(), 2);
    assert_eq!(
        parsed.points[1].centre.coordinates()[0],
        QSqrt5::new(q(-2), q(1))
    );
    assert_eq!(parsed.sha256, sha256_hex(&bytes));
}

#[test]
fn malformed_catalogues_are_refused() {
    let good = point("a", ORIGIN);
    let cases: Vec<(Vec<u8>, &str)> = vec![
        (b"{}".to_vec(), "missing fields"),
        (catalogue(&[]), "empty"),
        (
            format!(r#"{{"format":"other/1","points":[{good}]}}"#).into_bytes(),
            "format",
        ),
        (
            format!(r#"{{"format":"{FORMAT}","points":[{good}],"extra":1}}"#).into_bytes(),
            "unknown top-level field",
        ),
        (
            catalogue(&[good.replace(r#""expected""#, r#""note":1,"expected""#)]),
            "unknown point field",
        ),
        (catalogue(&[good.replace(r#","expected":["Local"]"#, "")]), "missing expected"),
        (catalogue(&[good.clone(), good.clone()]), "duplicate id"),
        (catalogue(&[point("", ORIGIN)]), "empty id"),
        (catalogue(&[good.replace(r#""group":"g""#, r#""group":"""#)]), "empty group"),
        (catalogue(&[good.replace(r#"["Local"]"#, r#"["Local","Local"]"#)]), "repeated expected"),
        (catalogue(&[good.replace(r#"["Local"]"#, r#"[""]"#)]), "empty expected"),
        (
            catalogue(&[point("a", r#"[["0","0"],["0","0"],["0","0"],["0","0"]]"#)]),
            "four coordinates",
        ),
        (
            catalogue(&[point("a", r#"[["0","0"],["0","0"],["0","0"],["0","0"],["0","0"],["0","0"]]"#)]),
            "six coordinates",
        ),
        (catalogue(&[point("a", &with_axis(0, ["1", "0"]))]), "outside the root"),
        (catalogue(&[point("a", &with_axis(2, ["0", "-1/5"]))]), "irrational outside"),
        (catalogue(&[point("a", &with_axis(0, ["0", "1/0"]))]), "zero denominator"),
        (catalogue(&[point("a", &with_axis(0, ["0.5", "0"]))]), "decimal"),
    ];
    for spelling in ["2/4", "+1", "01", "-0", "1/1", "3/-4", " 1", "1 "] {
        let bytes = catalogue(&[point("a", &with_axis(0, [spelling, "0"]))]);
        assert!(Catalogue::parse(&bytes).is_err(), "{spelling:?} accepted");
    }
    for (bytes, what) in cases {
        assert!(Catalogue::parse(&bytes).is_err(), "{what} accepted");
    }
}

#[test]
fn selections_keep_catalogue_order_and_refuse_bad_names() {
    let ids = ["a", "b", "c", "d"];
    assert_eq!(select(&ids, None).unwrap(), vec![0, 1, 2, 3]);
    let pick = |names: &[&str]| {
        let owned: Vec<String> = names.iter().map(|s| s.to_string()).collect();
        select(&ids, Some(&owned))
    };
    assert_eq!(pick(&["d", "b"]).unwrap(), vec![1, 3]);
    assert!(pick(&[]).is_err());
    assert!(pick(&["e"]).is_err());
    assert!(pick(&["a", "a"]).is_err());
}

#[test]
fn names_must_be_nonempty_and_distinct() {
    assert!(check_names([("a", "g"), ("b", "g")].into_iter()).is_ok());
    assert!(matches!(
        check_names(std::iter::empty()),
        Err(CatalogueError::Empty)
    ));
    assert!(matches!(
        check_names([("a", "g"), ("a", "h")].into_iter()),
        Err(CatalogueError::DuplicateId(_))
    ));
    assert!(matches!(
        check_names([("a", "")].into_iter()),
        Err(CatalogueError::EmptyName { index: 0 })
    ));
}
