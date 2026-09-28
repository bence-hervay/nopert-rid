use super::*;
use crate::experiments::probe::standin::StandIn;
use crate::experiments::run::{lines, parse, read};
use crate::testing::Scratch;

#[test]
fn the_generator_is_splitmix64() {
    let mut state = 0;
    assert_eq!(splitmix64(&mut state), 0xE220_A839_7B1D_CDAF);
    assert_eq!(splitmix64(&mut state), 0x6E78_9E6A_A1B9_65F4);
    assert_eq!(splitmix64(&mut state), 0x06C4_5D18_8009_454F);
}

fn sampled(depth: u32, count: u64, seed: u64) -> Result<Vec<String>, Error> {
    sample(&Sample { depth, count, seed })
}

#[test]
fn samples_are_distinct_paths_of_one_depth_in_search_order() {
    let a = sampled(12, 300, 5).unwrap();
    assert_eq!(a, sampled(12, 300, 5).unwrap());
    assert_ne!(a, sampled(12, 300, 6).unwrap());
    assert_eq!(a.len(), 300);
    assert!(a.iter().all(|p| p.len() == 12 && is_path(p)));
    assert!(a.windows(2).all(|w| w[0] < w[1]));
    // A whole population, and the root alone.
    assert_eq!(sampled(3, 8, 1).unwrap(), ["000", "001", "010", "011", "100", "101", "110", "111"]);
    assert_eq!(sampled(0, 1, 9).unwrap(), [""]);
    for (depth, count) in [(3, 9), (3, 0), (0, 2), (63, 1)] {
        assert!(matches!(sampled(depth, count, 1), Err(Error::Invalid(_))), "{depth} {count}");
    }
    let deep = sampled(62, 5, 1).unwrap();
    assert!(deep.iter().all(|p| p.len() == 62));
}

fn pilot(roots: Roots, depth_limit: usize, max_decisions: Option<u64>, transcript: std::path::PathBuf) -> config::Pilot {
    config::Pilot {
        roots,
        threads: NonZeroUsize::new(2).unwrap(),
        depth_limit,
        max_decisions,
        transcript,
    }
}

#[test]
fn roots_are_checked_and_put_in_search_order() {
    let paths = |p: &[&str]| Roots::Paths(p.iter().map(|s| s.to_string()).collect());
    let make = |roots, limit| pilot(roots, limit, None, "unused".into());
    assert_eq!(roots(&make(paths(&["11", "0", "10"]), 4)).unwrap(), ["0", "10", "11"]);
    for (bad, limit) in [
        (paths(&[]), 4),
        (paths(&["0", "0"]), 4),
        (paths(&["0a"]), 4),
        (paths(&["00000"]), 4),
        (paths(&["0"]), MAX_DEPTH + 1),
    ] {
        assert!(matches!(roots(&make(bad, limit)), Err(Error::Invalid(_))));
    }
}

/// An independent recursive search: records per label and depth, refused
/// boxes at the depth limit, and evaluations.
#[derive(Default, Debug, PartialEq)]
struct Reference {
    labels: BTreeMap<String, u64>,
    depths: BTreeMap<usize, u64>,
    unresolved: u64,
    evaluated: u64,
    splits: u64,
}

fn reference(probes: &StandIn, path: &str, limit: usize, out: &mut Reference) {
    let b = ConfigurationBox::from_path(path).unwrap();
    out.evaluated += 1;
    match probes.decide(&b).unwrap() {
        Some(data) => {
            let label = Label::of(&serde_json::to_value(data).unwrap()).unwrap();
            *out.labels.entry(label.to_string()).or_default() += 1;
            *out.depths.entry(path.len()).or_default() += 1;
        }
        None if path.len() == limit => out.unresolved += 1,
        None => {
            out.splits += 1;
            reference(probes, &format!("{path}0"), limit, out);
            reference(probes, &format!("{path}1"), limit, out);
        }
    }
}

#[test]
fn bounded_searches_match_an_independent_recursive_search() {
    let standin = StandIn::typical();
    let names = standin.components();
    for (root, limit) in [("", 9), ("0", 11), ("00000", 16), ("1101", 12)] {
        let result = search(&standin, &names, root, limit, None).unwrap();
        let mut expected = Reference::default();
        reference(&standin, root, limit, &mut expected);
        let found = Reference {
            labels: result.labels.clone(),
            depths: result.depths.clone(),
            unresolved: result.unresolved,
            evaluated: result.evaluated,
            splits: result.splits,
        };
        assert_eq!(found, expected, "root {root:?}");
        assert_eq!(result.records, expected.labels.values().sum::<u64>());
        assert_eq!(result.deepest, expected.depths.keys().next_back().copied());
        let halt = if expected.unresolved == 0 { "complete" } else { "depth-limit" };
        assert_eq!(result.halt, halt);
        let timed: u64 = result.outcomes.values().map(|t| t.boxes).sum();
        assert_eq!(timed, result.evaluated);
    }
    let limited = search(&standin, &names, "", 30, Some(10)).unwrap();
    assert_eq!((limited.halt.as_str(), limited.decisions), ("decision-limit", 10));
}

#[test]
fn extrapolations_follow_the_sampling_formula() {
    let e = extrapolate(&[1.0, 2.0, 3.0], 10.0);
    assert_eq!((e.mean, e.total), (2.0, 20.0));
    let expected = 10.0 * (1.0f64 / 3.0 * (1.0 - 0.3)).sqrt();
    assert!((e.standard_error.unwrap() - expected).abs() < 1e-12);
    assert_eq!(extrapolate(&[4.0], 8.0).standard_error, None);
    assert_eq!(extrapolate(&[1.0, 5.0], 2.0).standard_error, Some(0.0));
}

#[test]
fn a_sampled_pilot_writes_its_roots_and_an_estimate() {
    let dir = Scratch::new("pilot");
    let config = pilot(Roots::Sample(Sample { depth: 5, count: 6, seed: 3 }), 14, None, dir.join("p.jsonl"));
    let summary = run(&StandIn::typical(), &config).unwrap();
    let bytes = read(&config.transcript).unwrap();
    let lines = lines(&bytes).unwrap();
    assert_eq!(lines.len(), 8);
    let header: Header = parse(lines[0], 1).unwrap();
    assert_eq!(header.roots, sampled(5, 6, 3).unwrap());
    let results: Vec<RootResult> = lines[1..7].iter().map(|l| parse(l, 2).unwrap()).collect();
    assert_eq!(results.iter().map(|r| r.root.clone()).collect::<Vec<_>>(), header.roots);
    // The summary's estimates are floating point: read without the exact
    // spelling check.
    let written: Summary = serde_json::from_slice(lines[7]).unwrap();
    assert_eq!((written.records, written.evaluated, written.microseconds), (summary.records, summary.evaluated, summary.microseconds));
    assert_eq!(summary.records, results.iter().map(|r| r.records).sum::<u64>());
    let estimate = summary.estimate.unwrap();
    assert_eq!(estimate.population, 32.0);
    assert_eq!(estimate.records.total, 32.0 * summary.records as f64 / 6.0);
    // Each root's result is its own search, whatever the thread count.
    let names = StandIn::typical().components();
    for r in &results {
        let alone = search(&StandIn::typical(), &names, &r.root, 14, None).unwrap();
        assert_eq!((alone.labels, alone.depths, alone.evaluated), (r.labels.clone(), r.depths.clone(), r.evaluated));
    }
    assert!(matches!(run(&StandIn::typical(), &config), Err(Error::Transcript(_))));
}
