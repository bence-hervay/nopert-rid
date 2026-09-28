use super::*;
use crate::random::Random;
use num_bigint::BigInt;
use rid::arithmetic::exact::frac;
use rid::arithmetic::polynomial::Polynomial;
use rid::problem::configuration::coordinates;
use rid::components::domain::Domain;
use rid::components::record::RecordData;
use rid::components::Component;
use rid::problem::domain;
use rid::arithmetic::exact::QSqrt5;

fn r(n: i64, d: i64) -> Rational {
    Rational(frac(n, d))
}

fn vector(values: [(i64, i64); 5]) -> [Rational; 5] {
    values.map(|(n, d)| r(n, d))
}

/// `s ∈ [0, 2/3]` across, `r₃ ∈ [-2/5, 2/5]` up, at `t = 1/10`, `r₁ = r₂ = 0`.
fn definition(max_depth: usize) -> Definition {
    Definition {
        origin: vector([(0, 1), (1, 10), (0, 1), (0, 1), (0, 1)]),
        horizontal: vector([(1, 1), (0, 1), (0, 1), (0, 1), (0, 1)]),
        vertical: vector([(0, 1), (0, 1), (0, 1), (0, 1), (1, 1)]),
        plot: [[r(0, 1), r(2, 3)], [r(-2, 5), r(2, 5)]],
        thickness: vector([(0, 1), (0, 1), (1, 100), (1, 100), (0, 1)]),
        max_depth,
        max_cells: NonZeroUsize::new(1 << 16).unwrap(),
        components: Parts::all(),
    }
}

/// A stand-in for the components: a fixed rule on the box, covering every
/// category, that declines boxes crossing its boundaries.
fn stand_in(b: &ConfigurationBox, _: &Parts) -> Result<Option<Category>, BoxError> {
    let [s, _, _, _, r3] = b.axes();
    let exotic = |cover| Some(Category::Exotic { cover });
    Ok(if *s.lo() >= frac(1, 2) {
        Some(Category::Domain)
    } else if *s.lo() >= frac(1, 3) && *r3.lo() >= frac(1, 5) {
        Some(Category::Global)
    } else if *r3.hi() <= frac(-1, 5) {
        Some(Category::Local)
    } else if *s.hi() <= frac(1, 10) {
        exotic(if *r3.lo() >= Q::zero() { ExoticCover::Square } else if *r3.hi() <= Q::zero() { ExoticCover::CrossingPlus } else { return Ok(None) })
    } else if *r3.lo() >= frac(3, 10) {
        exotic(ExoticCover::Pentagon)
    } else if *s.lo() >= frac(1, 5) && *s.hi() <= frac(1, 4) {
        exotic(ExoticCover::ArcPlus)
    } else if *s.lo() >= frac(1, 4) && *s.hi() <= frac(3, 10) && *r3.hi() <= Q::zero() {
        exotic(ExoticCover::ArcMinus)
    } else if *s.lo() >= frac(5, 12) && *s.hi() <= frac(1, 2) && *r3.hi() <= Q::zero() {
        exotic(ExoticCover::EndpointPlus)
    } else {
        None
    })
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

#[test]
fn enclosures_are_the_exact_bounding_boxes_of_the_thickened_cells() {
    let mut random = Random::new(11);
    for case in 0..40 {
        let mut d = definition(6);
        for axis in 0..AXES {
            d.origin[axis] = Rational(random.rational(&frac(-1, 2), &frac(1, 2), 10));
            d.horizontal[axis] = Rational(random.rational(&frac(-1, 1), &frac(1, 1), 6));
            d.vertical[axis] = Rational(random.rational(&frac(-1, 1), &frac(1, 1), 6));
            d.thickness[axis] = Rational(if case % 3 == 0 { Q::zero() } else { random.rational(&Q::zero(), &frac(1, 50), 6) });
        }
        let Ok(slice) = Slice::new(d.clone()) else { continue };
        let path: String = (0..random.integer(0, 6)).map(|_| if random.next() % 2 == 0 { '0' } else { '1' }).collect();
        let rectangle = slice.rectangle(&path).unwrap();
        let enclosure = slice.enclosure(&rectangle);
        // Every thickened image of a corner lies inside, and each bound is
        // attained at some corner with an extreme thickness.
        for i in 0..AXES {
            let mut values = Vec::new();
            for a in [rectangle[0].lo(), rectangle[0].hi()] {
                for b in [rectangle[1].lo(), rectangle[1].hi()] {
                    let centre = &(&d.origin[i].0 + &(&d.horizontal[i].0 * a)) + &(&d.vertical[i].0 * b);
                    values.push(&centre - &d.thickness[i].0);
                    values.push(&centre + &d.thickness[i].0);
                }
            }
            assert_eq!(enclosure[i].lo(), values.iter().min().unwrap(), "case {case} axis {i}");
            assert_eq!(enclosure[i].hi(), values.iter().max().unwrap(), "case {case} axis {i}");
        }
        // Random interior points map inside.
        for _ in 0..8 {
            let a = random.rational(rectangle[0].lo(), rectangle[0].hi(), 12);
            let b = random.rational(rectangle[1].lo(), rectangle[1].hi(), 12);
            for i in 0..AXES {
                let x = &(&d.origin[i].0 + &(&d.horizontal[i].0 * &a)) + &(&d.vertical[i].0 * &b);
                assert!(enclosure[i].contains(&x));
            }
        }
    }
}

#[test]
fn clipping_keeps_boundary_contact_and_drops_disjoint_cells() {
    let mut d = definition(4);
    // A slice along s from -1 to 0: only its right edge s = 0 touches B₀.
    d.plot = [[r(-1, 1), r(0, 1)], [r(-1, 10), r(1, 10)]];
    let slice = Slice::new(d).unwrap();
    let touching = slice.clipped(&slice.rectangle("1").unwrap()).expect("s = 0 is in B₀");
    assert_eq!(touching.axes()[0], Interval::point(Q::zero()));
    assert!(slice.clipped(&slice.rectangle("0").unwrap()).is_none());
    // Clipping intersects every axis with the root box's.
    let full = definition(4);
    let slice = Slice::new(full).unwrap();
    let b = slice.clipped(slice.plot()).unwrap();
    assert_eq!(b.axes()[2], Interval::new(frac(-1, 100), frac(1, 100)).unwrap());
    assert!(ConfigurationBox::root().contains(&b));
}

#[test]
fn rectangles_halve_alternate_axes_and_refuse_bad_paths() {
    let slice = Slice::new(definition(3)).unwrap();
    let [a, b] = slice.rectangle("10").unwrap();
    assert_eq!(a, Interval::new(frac(1, 3), frac(2, 3)).unwrap());
    assert_eq!(b, Interval::new(frac(-2, 5), Q::zero()).unwrap());
    assert!(matches!(slice.rectangle("0101"), Err(SliceError::InvalidPath(_))));
    assert!(matches!(slice.rectangle("0a"), Err(SliceError::InvalidPath(_))));
}

#[test]
fn definitions_are_validated() {
    let mut d = definition(3);
    d.plot[1] = [r(1, 5), r(1, 5)];
    assert!(matches!(Slice::new(d), Err(SliceError::EmptyPlot { axis: 1 })));
    let mut d = definition(3);
    d.plot[0] = [r(1, 5), r(1, 7)];
    assert!(matches!(Slice::new(d), Err(SliceError::Number(_))));
    let mut d = definition(3);
    d.vertical = vector([(2, 1), (0, 1), (0, 1), (0, 1), (0, 1)]);
    assert!(matches!(Slice::new(d), Err(SliceError::DependentDirections)));
    let mut d = definition(3);
    d.thickness[4] = r(-1, 100);
    assert!(matches!(Slice::new(d), Err(SliceError::NegativeThickness { axis: 4 })));
    assert!(matches!(Slice::new(definition(MAX_DEPTH + 1)), Err(SliceError::DepthLimit { .. })));
    assert!(Slice::new(definition(MAX_DEPTH)).is_ok());
}

/// The exact area of a cell as a fraction of the plot: `2^-depth`.
fn area(path: &str) -> Q {
    Q::new(BigInt::from(1), BigInt::from(1) << path.len())
}

#[test]
fn computed_partitions_agree_with_an_independent_reconstruction() {
    let slice = Slice::new(definition(7)).unwrap();
    let partition = Partition::compute(&slice, &stand_in, threads(4)).unwrap();
    let cells = partition.cells();
    // The cells tile the plot: prefix-free paths whose areas sum to one.
    let total = cells.iter().fold(Q::zero(), |sum, c| sum + area(&c.path));
    assert_eq!(total, Q::one());
    for (k, a) in cells.iter().enumerate() {
        for b in &cells[k + 1..] {
            assert!(!b.path.starts_with(&a.path) && !a.path.starts_with(&b.path));
        }
    }
    let mut seen = std::collections::HashSet::new();
    for cell in cells {
        seen.insert(cell.outcome);
        let b = slice.clipped(&slice.rectangle(&cell.path).unwrap());
        // Each outcome is what the rule says about the cell's own box.
        match (&b, cell.outcome) {
            (None, Outcome::Outside) => {}
            (Some(b), Outcome::Eliminated(category)) => assert_eq!(stand_in(b, &Parts::all()).unwrap(), Some(category)),
            (Some(b), Outcome::Unresolved) => {
                assert_eq!(cell.path.len(), 7);
                assert_eq!(stand_in(b, &Parts::all()).unwrap(), None);
            }
            _ => panic!("cell {} has outcome {:?}", cell.path, cell.outcome),
        }
        // Every ancestor was split because the rule declined it.
        for k in 0..cell.path.len() {
            let parent = slice.clipped(&slice.rectangle(&cell.path[..k]).unwrap()).expect("split cells meet B₀");
            assert_eq!(stand_in(&parent, &Parts::all()).unwrap(), None, "ancestor {}", &cell.path[..k]);
        }
    }
    // The stand-in reaches every eliminated category and unresolved cells.
    assert!(seen.len() >= 10, "{seen:?}");
}

#[test]
fn partitions_do_not_depend_on_the_thread_count() {
    let slice = Slice::new(definition(6)).unwrap();
    let one = Partition::compute(&slice, &stand_in, threads(1)).unwrap();
    for n in [2, 3, 4, 16] {
        assert_eq!(Partition::compute(&slice, &stand_in, threads(n)).unwrap(), one);
    }
}

#[test]
fn classifier_errors_stop_the_subdivision() {
    let slice = Slice::new(definition(4)).unwrap();
    let failing = |b: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        if *b.axes()[0].hi() <= frac(1, 6) {
            Err("broken".into())
        } else {
            Ok(None)
        }
    };
    assert!(matches!(Partition::compute(&slice, &failing, threads(3)), Err(SliceError::Classifier(_))));
}

#[test]
fn domain_cells_of_the_crate_s_domain_test_hold_at_sampled_points() {
    // A real classifier from the crate: its Domain component. At random
    // exact points of each Domain cell some inequality of D fails.
    let domain_only = |b: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        Ok(Domain.check(b).map(|inequality| Category::of(&RecordData::Domain { inequality })))
    };
    let slice = Slice::new(definition(6)).unwrap();
    let partition = Partition::compute(&slice, &domain_only, threads(4)).unwrap();
    let (view, rotation) = coordinates();
    let inequalities: Vec<Polynomial> = (0..domain::CONSTRAINT_COUNT)
        .map(|i| domain::polynomial(i, &view, &rotation).unwrap())
        .collect();
    let mut random = Random::new(12);
    let mut domain_cells = 0;
    for cell in partition.cells().iter().filter(|c| c.outcome == Outcome::Eliminated(Category::Domain)) {
        domain_cells += 1;
        let b = slice.clipped(&slice.rectangle(&cell.path).unwrap()).unwrap();
        for _ in 0..3 {
            let x: [QSqrt5; 5] = std::array::from_fn(|j| {
                QSqrt5::from_rational(random.rational(b.axes()[j].lo(), b.axes()[j].hi(), 10))
            });
            assert!(
                inequalities.iter().any(|p| p.evaluate(&x) > QSqrt5::zero()),
                "a point of Domain cell {} lies in D",
                cell.path
            );
        }
    }
    assert!(domain_cells > 0);
}

#[test]
fn records_map_to_categories_by_type() {
    let record = |text: &str| Category::of(&serde_json::from_str::<RecordData>(text).unwrap());
    assert_eq!(record(r#"{"component":"Domain","inequality":13}"#), Category::Domain);
    assert_eq!(record(r#"{"component":"Global","edge":[16,17],"vertex":3}"#), Category::Global);
    assert_eq!(record(r#"{"component":"Local","cover":9}"#), Category::Local);
    for (name, cover) in [
        ("square", ExoticCover::Square),
        ("pentagon", ExoticCover::Pentagon),
        ("arc+", ExoticCover::ArcPlus),
        ("arc-", ExoticCover::ArcMinus),
        ("endpoint+", ExoticCover::EndpointPlus),
        ("endpoint-", ExoticCover::EndpointMinus),
        ("crossing+", ExoticCover::CrossingPlus),
        ("crossing-", ExoticCover::CrossingMinus),
    ] {
        let text = format!(r#"{{"component":"Exotic","cover":"{name}"}}"#);
        assert_eq!(record(&text), Category::Exotic { cover });
        // The category's own spelling of the cover is the record's.
        assert_eq!(serde_json::to_value(Category::Exotic { cover }).unwrap()["cover"], name);
    }
}

const POLICY: &str = "0123abcd";

fn file_of(partition: &Partition) -> serde_json::Value {
    serde_json::from_str(&partition.to_json(POLICY)).unwrap()
}

fn read(value: &serde_json::Value) -> Result<Partition, SliceError> {
    Partition::from_json(value.to_string().as_bytes(), POLICY)
}

#[test]
fn partition_files_round_trip_one_cell_per_line() {
    let slice = Slice::new(definition(5)).unwrap();
    let partition = Partition::compute(&slice, &stand_in, threads(4)).unwrap();
    let text = partition.to_json(POLICY);
    assert_eq!(text.lines().count(), partition.cells().len() + 2);
    assert_eq!(Partition::from_json(text.as_bytes(), POLICY).unwrap(), partition);
}

#[test]
fn corrupted_partition_files_are_refused() {
    let slice = Slice::new(definition(5)).unwrap();
    let partition = Partition::compute(&slice, &stand_in, threads(4)).unwrap();
    let good = file_of(&partition);
    let cells = good["cells"].as_array().unwrap().clone();
    let with_cells = |cells: Vec<serde_json::Value>| {
        let mut v = good.clone();
        v["cells"] = serde_json::Value::Array(cells);
        v
    };
    assert!(read(&good).is_ok());
    assert!(matches!(Partition::from_json(partition.to_json(POLICY).as_bytes(), "other"), Err(SliceError::Policy { .. })));
    let mut v = good.clone();
    v["format"] = "rid-figures-partition/0".into();
    assert!(matches!(read(&v), Err(SliceError::Format(_))));
    let mut v = good.clone();
    v["extra"] = 1.into();
    assert!(matches!(read(&v), Err(SliceError::Json(_))));
    let mut v = good.clone();
    v["slice"]["max_depth"] = 6.into();
    assert!(matches!(read(&v), Err(SliceError::Gap { .. }) | Err(SliceError::EarlyUnresolved { .. })));
    let mut v = good.clone();
    v["slice"]["origin"][0] = serde_json::json!(0);
    assert!(matches!(read(&v), Err(SliceError::Json(_))), "JSON numbers are not rationals");
    let mut v = good.clone();
    v["slice"]["origin"][1] = "2/20".into();
    assert!(matches!(read(&v), Err(SliceError::Json(_))), "non-canonical rational");
    let mut v = good.clone();
    v["slice"]["thickness"][0] = "-1".into();
    assert!(matches!(read(&v), Err(SliceError::NegativeThickness { .. })));
    // A missing cell leaves a gap.
    for k in [0, cells.len() / 2, cells.len() - 1] {
        let mut fewer = cells.clone();
        fewer.remove(k);
        assert!(matches!(read(&with_cells(fewer)), Err(SliceError::Gap { .. })), "cell {k}");
    }
    // A repeated or reordered cell.
    let mut repeated = cells.clone();
    repeated.insert(3, cells[3].clone());
    assert!(matches!(read(&with_cells(repeated)), Err(SliceError::Order { .. })));
    let mut swapped = cells.clone();
    swapped.swap(1, 2);
    assert!(matches!(read(&with_cells(swapped)), Err(SliceError::Order { .. })));
    // A cell inside another cell, valid in every other respect: at the
    // depth limit, unresolved, in order.
    let leaf = cells.iter().position(|c| c["path"].as_str().unwrap().len() < 5).unwrap();
    let mut inner = cells.clone();
    let path = format!("{:0<5}", cells[leaf]["path"].as_str().unwrap());
    inner.push(serde_json::json!({"path": path, "outcome": "unresolved"}));
    inner.sort_by_key(|c| {
        let p = c["path"].as_str().unwrap().to_owned();
        (p.len(), p)
    });
    assert!(matches!(read(&with_cells(inner)), Err(SliceError::Overlap { .. })));
    // Outcomes that contradict the boxes or the depth.
    let unresolved = cells.iter().position(|c| c["outcome"] == "unresolved").unwrap();
    let mut outside = cells.clone();
    outside[unresolved]["outcome"] = "outside".into();
    assert!(matches!(read(&with_cells(outside)), Err(SliceError::WrongOutside { .. })));
    let mut early = cells.clone();
    early[leaf]["outcome"] = "unresolved".into();
    assert!(matches!(read(&with_cells(early)), Err(SliceError::EarlyUnresolved { .. })));
    let mut unknown = cells.clone();
    unknown[0]["outcome"] = serde_json::json!({"eliminated": {"component": "Domain", "inequality": 3}});
    assert!(matches!(read(&with_cells(unknown)), Err(SliceError::Json(_))));
    let mut bad_path = cells.clone();
    bad_path[0]["path"] = "2".into();
    assert!(matches!(read(&with_cells(bad_path)), Err(SliceError::InvalidPath(_))));
}

#[test]
fn outside_cells_are_checked_against_the_root_box() {
    // A slice mostly outside B₀: s ∈ [-2, 2/3].
    let mut d = definition(3);
    d.plot[0] = [r(-2, 1), r(2, 3)];
    let slice = Slice::new(d).unwrap();
    let partition = Partition::compute(&slice, &stand_in, threads(2)).unwrap();
    assert!(partition.cells().iter().any(|c| c.outcome == Outcome::Outside));
    let good = file_of(&partition);
    assert!(read(&good).is_ok());
    // Declaring an outside cell eliminated contradicts its box.
    let mut v = good.clone();
    let k = partition.cells().iter().position(|c| c.outcome == Outcome::Outside).unwrap();
    v["cells"][k]["outcome"] = serde_json::json!({"eliminated": {"component": "Global"}});
    assert!(matches!(read(&v), Err(SliceError::WrongOutside { .. })));
    // Splitting an outside cell is refused too.
    let mut v = good.clone();
    let path = partition.cells()[k].path.clone();
    let mut cells = v["cells"].as_array().unwrap().clone();
    cells.remove(k);
    for bit in ["0", "1"] {
        cells.push(serde_json::json!({"path": format!("{path}{bit}"), "outcome": "outside"}));
    }
    cells.sort_by_key(|c| {
        let p = c["path"].as_str().unwrap().to_owned();
        (p.len(), p)
    });
    v["cells"] = serde_json::Value::Array(cells);
    assert!(matches!(read(&v), Err(SliceError::WrongOutside { .. })));
}

#[test]
fn categories_are_written_strictly() {
    let text = |c: Category| serde_json::to_string(&c).unwrap();
    assert_eq!(text(Category::Local), r#"{"component":"Local"}"#);
    assert_eq!(
        text(Category::Exotic { cover: ExoticCover::ArcMinus }),
        r#"{"component":"Exotic","cover":"arc-"}"#
    );
    let parse = |t: &str| serde_json::from_str::<Category>(t);
    assert_eq!(parse(r#"{"component":"Global"}"#).unwrap(), Category::Global);
    for bad in [
        r#"{"component":"Global","cover":"arc+"}"#,
        r#"{"component":"Exotic"}"#,
        r#"{"component":"Domain","inequality":3}"#,
        r#"{"component":"domain"}"#,
        r#"{}"#,
    ] {
        assert!(parse(bad).is_err(), "{bad}");
    }
}

mod review;
