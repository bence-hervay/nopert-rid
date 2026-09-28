//! Adversarial tests of the slice module.
use super::*;
use rid::arithmetic::exact::frac;
use std::sync::atomic::AtomicUsize;
use std::time::Instant;

fn r(n: i64, d: i64) -> Rational {
    Rational(frac(n, d))
}

fn definition(max_depth: usize) -> Definition {
    let v = |x: [(i64, i64); 5]| x.map(|(n, d)| r(n, d));
    Definition {
        origin: v([(0, 1), (1, 10), (0, 1), (0, 1), (0, 1)]),
        horizontal: v([(1, 1), (0, 1), (0, 1), (0, 1), (0, 1)]),
        vertical: v([(0, 1), (0, 1), (0, 1), (0, 1), (1, 1)]),
        plot: [[r(0, 1), r(2, 3)], [r(-2, 5), r(2, 5)]],
        thickness: v([(0, 1); 5]),
        max_depth,
        max_cells: NonZeroUsize::new(1 << 16).unwrap(),
        components: Parts::all(),
    }
}

#[test]
fn categories_refuse_an_explicit_null_cover() {
    // One spelling per category: `cover: null` is refused; absence still
    // means no cover.
    for text in [r#"{"component":"Domain","cover":null}"#, r#"{"component":"Exotic","cover":null}"#] {
        assert!(serde_json::from_str::<Category>(text).is_err(), "{text}");
    }
    assert_eq!(serde_json::from_str::<Category>(r#"{"component":"Domain"}"#).unwrap(), Category::Domain);
}

#[test]
fn inconsistent_records_have_no_category() {
    // Records are typed: JSON with inconsistent fields is no record at all
    // and has no category.
    use rid::components::record::RecordData;
    for v in [
        serde_json::json!({"component": "Domain", "cover": "crossing", "witness": 7}),
        serde_json::json!({"component": "Local", "cover": "crossing"}),
    ] {
        assert!(serde_json::from_value::<RecordData>(v.clone()).is_err(), "{v}");
    }
}

#[test]
fn duplicate_keys_in_cells_are_refused() {
    let slice = Slice::new(definition(0)).unwrap();
    let text = format!(
        "{{\"format\":\"{FORMAT}\",\"policy\":\"p\",\"slice\":{},\"cells\":[{{\"path\":\"\",\"path\":\"\",\"outcome\":\"unresolved\"}}]}}",
        serde_json::to_string(slice.definition()).unwrap()
    );
    assert!(matches!(Partition::from_json(text.as_bytes(), "p"), Err(SliceError::Json(_))));
}

#[test]
fn an_empty_cell_list_is_a_gap() {
    let slice = Slice::new(definition(3)).unwrap();
    let text = format!(
        "{{\"format\":\"{FORMAT}\",\"policy\":\"p\",\"slice\":{},\"cells\":[]}}",
        serde_json::to_string(slice.definition()).unwrap()
    );
    assert!(matches!(Partition::from_json(text.as_bytes(), "p"), Err(SliceError::Gap { .. })));
}

/// The subdivision is bounded by `max_cells`, not only by `max_depth ≤ 40`
/// (up to 2^40 cells): a classifier that declines
/// everywhere is stopped before the layer that would exceed the limit is
/// classified.
#[test]
fn subdivision_size_is_bounded_by_max_cells() {
    let calls = AtomicUsize::new(0);
    let declining = |_: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        Ok(None)
    };
    let mut d = definition(MAX_DEPTH);
    d.max_cells = NonZeroUsize::new(1000).unwrap();
    let slice = Slice::new(d).unwrap();
    let start = Instant::now();
    let refused = Partition::compute(&slice, &declining, NonZeroUsize::new(4).unwrap());
    assert!(matches!(refused, Err(SliceError::TooManyCells { max_cells: 1000 })));
    // Layers of 1, 2, …, 512 nodes were classified; the next (1024) was not.
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 1023);
    assert!(start.elapsed().as_secs() < 10);
    // Exactly enough cells is accepted.
    let mut d = definition(10);
    d.max_cells = NonZeroUsize::new(1 << 10).unwrap();
    let full = Partition::compute(&Slice::new(d).unwrap(), &declining, NonZeroUsize::new(4).unwrap()).unwrap();
    assert_eq!(full.cells().len(), 1 << 10);
}

/// A slab touching B₀ only at a corner point is classified on a degenerate
/// (point) box; compute and from_json agree on it.
#[test]
fn corner_contact_is_consistent() {
    let mut d = definition(2);
    // s ∈ [-1, 0], r₃ ∈ [2/5, 1]: only (0, 2/5) is in B₀'s projection.
    d.plot = [[r(-1, 1), r(0, 1)], [r(2, 5), r(1, 1)]];
    let slice = Slice::new(d).unwrap();
    let boxes = std::sync::Mutex::new(Vec::new());
    let recorder = |b: &ConfigurationBox, _: &Parts| -> Result<Option<Category>, BoxError> {
        boxes.lock().unwrap().push(b.clone());
        Ok(Some(Category::Global))
    };
    let partition = Partition::compute(&slice, &recorder, NonZeroUsize::new(2).unwrap()).unwrap();
    let text = partition.to_json("p");
    assert_eq!(Partition::from_json(text.as_bytes(), "p").unwrap(), partition);
    let boxes = boxes.into_inner().unwrap();
    assert!(boxes.iter().all(|b| ConfigurationBox::root().contains(b)));
    println!("corner contact: {:?}", partition.cells());
}
