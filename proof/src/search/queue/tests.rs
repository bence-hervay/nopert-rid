use super::*;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::SeqCst};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Records appends and the number of appends at each flush.
#[derive(Default)]
struct Memory {
    records: Vec<(String, String)>,
    flushes: Vec<usize>,
}

impl Output<String> for Memory {
    fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
        assert_eq!(path, data);
        if let Some((last, _)) = self.records.last() {
            assert!(search_order(last) < search_order(path), "{last} then {path}");
        }
        self.records.push((path.into(), data.clone()));
        Ok(())
    }
    fn flush(&mut self) -> Result<(), BoxError> {
        self.flushes.push(self.records.len());
        Ok(())
    }
}

impl Memory {
    fn paths(&self) -> Vec<&str> {
        self.records.iter().map(|(p, _)| p.as_str()).collect()
    }
}

fn nonzero(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn options(threads: usize, window: usize) -> Options {
    Options {
        threads: nonzero(threads),
        window: nonzero(window),
        depth_limit: 64,
        max_decisions: None,
        after: None,
    }
}

fn paths(depth: usize) -> Vec<String> {
    (0..(1usize << depth)).map(|n| format!("{n:0depth$b}")).collect()
}

fn go() -> AtomicBool {
    AtomicBool::new(false)
}

fn wait_until(mut condition: impl FnMut() -> bool) {
    let deadline = Instant::now() + Duration::from_secs(20);
    while !condition() {
        assert!(Instant::now() < deadline, "concurrency condition timed out");
        thread::sleep(Duration::from_millis(1));
    }
}

fn hash(path: &str, seed: u64) -> u64 {
    path.bytes().fold(seed.wrapping_add(0x9e3779b97f4a7c15), |n, b| {
        n.wrapping_mul(6364136223846793005).wrapping_add(u64::from(b) + 1)
    })
}

/// A deterministic decision rule with trees of a few hundred boxes.
fn accept(path: &str, seed: u64) -> bool {
    path.len() >= 8 || (path.len() >= 2 && hash(path, seed) % 4 == 0)
}

/// Sequential breadth-first reference: the records in order and the count
/// of decisions.
fn reference(seed: u64) -> (Vec<(String, String)>, u64) {
    let mut queue = VecDeque::from([String::new()]);
    let mut records = Vec::new();
    let mut decisions = 0;
    while let Some(path) = queue.pop_front() {
        decisions += 1;
        if accept(&path, seed) {
            records.push((path.clone(), path));
        } else {
            queue.push_back(format!("{path}0"));
            queue.push_back(format!("{path}1"));
        }
    }
    (records, decisions)
}

/// Independent (recursive) minimal complement of a set of records.
fn complement(records: &[(String, String)]) -> Vec<String> {
    fn visit(path: String, records: &[(String, String)], out: &mut Vec<String>) {
        if records.iter().any(|(p, _)| *p == path) {
            return;
        }
        if records.iter().any(|(p, _)| p.starts_with(&path)) {
            visit(format!("{path}0"), records, out);
            visit(format!("{path}1"), records, out);
        } else {
            out.push(path);
        }
    }
    let mut out = Vec::new();
    visit(String::new(), records, &mut out);
    out.sort_by(|a, b| search_order(a).cmp(&search_order(b)));
    out
}

/// Random delays that differ between calls, not only between paths.
fn jitter(counter: &AtomicU64, seed: u64) {
    let n = hash(&counter.fetch_add(1, SeqCst).to_string(), seed);
    match n % 8 {
        0 => thread::sleep(Duration::from_micros(50 + n % 400)),
        1 | 2 => thread::sleep(Duration::from_micros(n % 40)),
        _ => thread::yield_now(),
    }
}

// ---------------------------------------------------------------------------
// Order, determinism and results
// ---------------------------------------------------------------------------

#[test]
fn empty_frontier_and_a_single_record() {
    let mut out = Memory::default();
    let r = run(Vec::new(), &options(3, 4), &go(), |_| panic!("nothing to do"), &mut out).unwrap();
    assert_eq!((r.halt, r.evaluated, r.unresolved), (Halt::Complete, 0, 0));
    assert_eq!(out.flushes, [0]);
    let accept_all = |p: &str| Ok(Some(p.to_string()));
    let r = run(vec![String::new()], &options(3, 4), &go(), accept_all, &mut out).unwrap();
    assert_eq!(
        (r.halt, r.evaluated, r.decisions, r.records, r.unresolved),
        (Halt::Complete, 1, 1, 1, 0)
    );
    assert_eq!(out.paths(), [""]);
}

#[test]
fn invalid_inputs_start_no_work_and_no_output() {
    for (frontier, overlap) in [
        (vec!["2"], false),
        (vec!["0", "a"], false),
        (vec!["0", "0"], true),
        (vec!["", "01"], true),
        (vec!["10", "101"], true),
        (vec!["11", "0", "1"], true),
    ] {
        let mut out = Memory::default();
        let error = run(
            frontier.iter().map(|p| p.to_string()).collect(),
            &options(2, 2),
            &go(),
            |_| panic!("invalid input"),
            &mut out,
        )
        .unwrap_err();
        assert_eq!(matches!(error, Error::Overlap(..)), overlap, "{frontier:?}");
        assert_eq!(matches!(error, Error::InvalidPath(_)), !overlap, "{frontier:?}");
        assert!(out.records.is_empty() && out.flushes.is_empty());
    }
    let error = run::<String, _, _>(
        vec![String::new()],
        &Options {
            after: Some("01x".into()),
            ..options(1, 1)
        },
        &go(),
        |_| panic!("invalid after"),
        &mut Memory::default(),
    )
    .unwrap_err();
    assert!(matches!(error, Error::InvalidPath(_)));
}

#[test]
fn every_thread_count_and_window_commits_the_sequential_reference() {
    determinism_campaign(0..48);
}

/// The heavier campaign: 1000 further trees. Run with `cargo test --release
/// --lib search::queue::tests::determinism_campaign_heavy -- --ignored`.
#[test]
#[ignore]
fn determinism_campaign_heavy() {
    determinism_campaign(48..1048);
}

fn determinism_campaign(seeds: std::ops::Range<u64>) {
    for seed in seeds {
        let (expected, decisions) = reference(seed);
        for threads in [1, 2, 3, 5, 8] {
            for window in [1, threads, 3 * threads + 1] {
                let counter = AtomicU64::new(0);
                let mut out = Memory::default();
                let r = run(
                    vec![String::new()],
                    &options(threads, window),
                    &go(),
                    |p| {
                        jitter(&counter, seed);
                        Ok(accept(p, seed).then(|| p.to_string()))
                    },
                    &mut out,
                )
                .unwrap();
                assert_eq!(out.records, expected, "seed {seed} threads {threads} window {window}");
                assert_eq!(
                    (r.halt, r.decisions, r.evaluated, r.records, r.unresolved),
                    (Halt::Complete, decisions, decisions, expected.len() as u64, 0)
                );
                assert_eq!(r.splits, decisions - r.records);
                assert_eq!(out.flushes.last(), Some(&expected.len()));
            }
        }
    }
}

#[test]
fn forced_reverse_completion_still_commits_in_order() {
    let next = AtomicUsize::new(7);
    let mut out = Memory::default();
    run(
        paths(3),
        &options(8, 8),
        &go(),
        |p| {
            let index = usize::from_str_radix(p, 2).unwrap();
            wait_until(|| next.load(SeqCst) == index);
            if index > 0 {
                next.fetch_sub(1, SeqCst);
            }
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!(out.paths(), paths(3));
}

#[test]
fn mixed_depth_frontiers_merge_with_new_children() {
    let mut out = Memory::default();
    let frontier = ["11", "00", "10", "010", "011"].map(String::from).to_vec();
    run(
        frontier,
        &options(4, 4),
        &go(),
        |p| Ok((p.len() == 3).then(|| p.to_string())),
        &mut out,
    )
    .unwrap();
    assert_eq!(out.paths(), paths(3));
}

#[test]
fn a_single_worker_takes_tickets_in_order() {
    let seen = Mutex::new(Vec::new());
    run(
        vec![String::new()],
        &options(1, 5),
        &go(),
        |p| {
            seen.lock().unwrap().push(p.to_string());
            Ok((p.len() == 3).then(|| p.to_string()))
        },
        &mut Memory::default(),
    )
    .unwrap();
    let seen = seen.into_inner().unwrap();
    let mut sorted = seen.clone();
    sorted.sort_by(|a, b| search_order(a).cmp(&search_order(b)));
    assert_eq!(seen, sorted);
    assert_eq!(seen.len(), 15);
}

// ---------------------------------------------------------------------------
// Scheduling: window, refill and barrier
// ---------------------------------------------------------------------------

/// Counts appends so that workers can see how far commits have progressed.
struct Counting<'a> {
    memory: Memory,
    appended: &'a AtomicUsize,
}

impl Output<String> for Counting<'_> {
    fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
        self.memory.append(path, data)?;
        self.appended.fetch_add(1, SeqCst);
        Ok(())
    }
    fn flush(&mut self) -> Result<(), BoxError> {
        self.memory.flush()
    }
}

#[test]
fn the_window_bounds_assigned_but_uncommitted_tickets() {
    for threads in [1, 2, 4, 8] {
        for window in [1, 2, 3, 5, 8, 13] {
            let appended = AtomicUsize::new(0);
            let running = AtomicUsize::new(0);
            let peak = AtomicUsize::new(0);
            let counter = AtomicU64::new(0);
            let mut out = Counting {
                memory: Memory::default(),
                appended: &appended,
            };
            // One layer of accepted boxes: every commit is an append, and a
            // ticket is its path's index.
            run(
                paths(6),
                &options(threads, window),
                &go(),
                |p| {
                    let ticket = usize::from_str_radix(p, 2).unwrap();
                    assert!(ticket < appended.load(SeqCst) + window, "ticket {ticket}");
                    let now = running.fetch_add(1, SeqCst) + 1;
                    peak.fetch_max(now, SeqCst);
                    jitter(&counter, (threads * 100 + window) as u64);
                    running.fetch_sub(1, SeqCst);
                    Ok(Some(p.into()))
                },
                &mut out,
            )
            .unwrap();
            assert_eq!(out.memory.records.len(), 64);
            assert!(peak.load(SeqCst) <= threads.min(window));
        }
    }
}

#[test]
fn a_stalled_head_fills_exactly_the_window() {
    let started = AtomicUsize::new(0);
    let later = AtomicUsize::new(0);
    let mut out = Memory::default();
    let r = run(
        paths(5),
        &options(4, 4),
        &go(),
        |p| {
            started.fetch_add(1, SeqCst);
            if p == "00000" {
                wait_until(|| later.load(SeqCst) == 3);
                thread::sleep(Duration::from_millis(30));
                assert_eq!(started.load(SeqCst), 4);
            } else {
                later.fetch_add(1, SeqCst);
            }
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!((out.records.len(), r.evaluated), (32, 32));
}

#[test]
fn workers_refill_without_waiting_for_a_slow_head() {
    let later = AtomicUsize::new(0);
    let mut out = Memory::default();
    run(
        paths(3),
        &options(2, 8),
        &go(),
        |p| {
            if p == "000" {
                // Only the other worker can make this progress.
                wait_until(|| later.load(SeqCst) >= 5);
            } else {
                later.fetch_add(1, SeqCst);
            }
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!(out.records.len(), 8);
}

#[test]
fn workers_refill_while_the_coordinator_is_stuck_in_output() {
    struct Slow<'a> {
        memory: Memory,
        done: &'a AtomicUsize,
    }
    impl Output<String> for Slow<'_> {
        fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
            if self.memory.records.is_empty() {
                // The first commit waits until workers filled the window.
                wait_until(|| self.done.load(SeqCst) >= 6);
            }
            self.memory.append(path, data)
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            self.memory.flush()
        }
    }
    let done = AtomicUsize::new(0);
    let mut out = Slow {
        memory: Memory::default(),
        done: &done,
    };
    run(
        paths(4),
        &options(2, 6),
        &go(),
        |p| {
            done.fetch_add(1, SeqCst);
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!(out.memory.records.len(), 16);
}

#[test]
fn no_deeper_box_starts_before_the_layer_is_durable() {
    struct Durable<'a> {
        flushed: &'a AtomicUsize,
        memory: Memory,
    }
    impl Output<String> for Durable<'_> {
        fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
            self.memory.append(path, data)
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            thread::sleep(Duration::from_millis(3));
            self.flushed.fetch_add(1, SeqCst);
            self.memory.flush()
        }
    }
    let flushed = AtomicUsize::new(0);
    let mut out = Durable {
        flushed: &flushed,
        memory: Memory::default(),
    };
    run(
        vec![String::new()],
        &options(4, 16),
        &go(),
        |p| {
            // Layers 0..depth have each been flushed exactly once.
            assert_eq!(flushed.load(SeqCst), p.len());
            Ok((p.len() == 4 || p.ends_with("01")).then(|| p.to_string()))
        },
        &mut out,
    )
    .unwrap();
    // One flush per layer, the last also closing the run.
    assert_eq!(out.memory.flushes.len(), 6);
    assert_eq!(out.memory.flushes.last(), Some(&out.memory.records.len()));
}

#[test]
fn a_failed_barrier_flush_starts_nothing_deeper() {
    struct Fail;
    impl Output<String> for Fail {
        fn append(&mut self, _: &str, _: &String) -> Result<(), BoxError> {
            panic!("no record is produced")
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            Err("disk full".into())
        }
    }
    let started = AtomicUsize::new(0);
    let r = run(
        vec![String::new()],
        &options(8, 8),
        &go(),
        |p| {
            started.fetch_add(1, SeqCst);
            assert!(p.is_empty());
            Ok(None)
        },
        &mut Fail,
    );
    assert!(matches!(r, Err(Error::Output(_))));
    assert_eq!(started.load(SeqCst), 1);
}

// ---------------------------------------------------------------------------
// Failures and cancellation
// ---------------------------------------------------------------------------

#[test]
fn evaluator_errors_and_panics_stop_commits_at_their_ticket() {
    for panics in [false, true] {
        for threads in [1, 4] {
            let mut out = Memory::default();
            let error = run(
                paths(3),
                &options(threads, 8),
                &go(),
                |p| {
                    if p == "010" {
                        thread::sleep(Duration::from_millis(10));
                        if panics {
                            panic!("injected evaluator panic");
                        }
                        return Err("injected evaluator error".into());
                    }
                    Ok(Some(p.into()))
                },
                &mut out,
            )
            .unwrap_err();
            match error {
                Error::Evaluator { path, .. } if !panics => assert_eq!(path, "010"),
                Error::EvaluatorPanic { path, message } if panics => {
                    assert_eq!(path, "010");
                    assert!(message.contains("injected evaluator panic"));
                }
                other => panic!("{other}"),
            }
            assert_eq!(out.paths(), ["000", "001"]);
        }
    }
}

#[test]
fn an_early_failure_behind_a_full_window_commits_nothing() {
    for panics in [false, true] {
        let later = AtomicUsize::new(0);
        let mut out = Memory::default();
        let result = run(
            paths(5),
            &options(4, 4),
            &go(),
            |p| {
                if p == "00000" {
                    wait_until(|| later.load(SeqCst) == 3);
                    if panics {
                        panic!("injected head panic");
                    }
                    return Err("injected head error".into());
                }
                later.fetch_add(1, SeqCst);
                Ok(Some(p.into()))
            },
            &mut out,
        );
        assert!(result.is_err());
        assert_eq!(later.load(SeqCst), 3);
        assert!(out.records.is_empty());
    }
}

#[test]
fn output_errors_and_panics_retire_every_worker() {
    struct Broken {
        calls: usize,
        panics: bool,
    }
    impl Output<String> for Broken {
        fn append(&mut self, _: &str, _: &String) -> Result<(), BoxError> {
            self.calls += 1;
            if self.calls == 2 {
                if self.panics {
                    panic!("injected output panic");
                }
                return Err("injected output error".into());
            }
            Ok(())
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            Ok(())
        }
    }
    for panics in [false, true] {
        let live = AtomicUsize::new(0);
        let mut out = Broken { calls: 0, panics };
        let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
            run(
                paths(4),
                &options(4, 8),
                &go(),
                |p| {
                    live.fetch_add(1, SeqCst);
                    thread::sleep(Duration::from_millis(2));
                    live.fetch_sub(1, SeqCst);
                    Ok(Some(p.to_string()))
                },
                &mut out,
            )
        }));
        match result {
            Ok(Err(Error::Output(_))) => assert!(!panics),
            Err(_) => assert!(panics, "an output panic reaches the caller"),
            Ok(other) => panic!("{:?}", other.map(|r| r.halt)),
        }
        assert_eq!(live.load(SeqCst), 0);
        assert_eq!(out.calls, 2);
    }
}

#[test]
fn a_stop_before_the_start_or_a_zero_budget_evaluates_nothing() {
    let stopped = AtomicBool::new(true);
    let r = run::<String, _, _>(
        vec![String::new()],
        &options(4, 4),
        &stopped,
        |_| panic!("already stopped"),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!((r.halt, r.evaluated, r.unresolved), (Halt::Stopped, 0, 1));
    let r = run::<String, _, _>(
        vec!["0".into()],
        &Options {
            after: Some("111".into()),
            ..options(4, 4)
        },
        &stopped,
        |_| panic!("stopped before reconstruction"),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!(
        (r.halt, r.evaluated, r.reconstructed, r.unresolved),
        (Halt::Stopped, 0, 0, 1)
    );
    let saved = "1".repeat(40);
    let frontier = complement(&[(saved.clone(), saved.clone())]);
    let r = run::<String, _, _>(
        frontier,
        &Options {
            after: Some(saved),
            max_decisions: Some(0),
            ..options(4, 4)
        },
        &go(),
        |_| panic!("zero budget"),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!(
        (r.halt, r.evaluated, r.reconstructed, r.unresolved),
        (Halt::DecisionLimit, 0, 0, 40)
    );
}

#[test]
fn a_stop_during_a_huge_reconstruction_is_observed_promptly() {
    let stop = go();
    thread::scope(|scope| {
        scope.spawn(|| {
            thread::sleep(Duration::from_millis(20));
            stop.store(true, SeqCst);
        });
        let started = Instant::now();
        let r = run::<String, _, _>(
            vec!["0".into()],
            &Options {
                after: Some("1".repeat(40)),
                ..options(2, 2)
            },
            &stop,
            |_| panic!("reconstruction has not reached new work"),
            &mut Memory::default(),
        )
        .unwrap();
        assert!(started.elapsed() < Duration::from_secs(10));
        assert_eq!((r.halt, r.evaluated), (Halt::Stopped, 0));
        assert!(r.reconstructed > 0);
        assert_eq!(r.unresolved, r.reconstructed + 1);
    });
}

#[test]
fn a_stop_commits_nothing_more_and_waits_for_running_evaluations() {
    struct Interrupt<'a> {
        stop: &'a AtomicBool,
        memory: Memory,
    }
    impl Output<String> for Interrupt<'_> {
        fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
            self.memory.append(path, data)?;
            if self.memory.records.len() == 3 {
                self.stop.store(true, SeqCst);
            }
            Ok(())
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            self.memory.flush()
        }
    }
    let stop = go();
    let live = AtomicUsize::new(0);
    let mut out = Interrupt {
        stop: &stop,
        memory: Memory::default(),
    };
    let r = run(
        paths(5),
        &options(4, 8),
        &stop,
        |p| {
            live.fetch_add(1, SeqCst);
            thread::sleep(Duration::from_millis(1));
            live.fetch_sub(1, SeqCst);
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!(live.load(SeqCst), 0);
    assert_eq!((r.halt, r.decisions, r.unresolved), (Halt::Stopped, 3, 29));
    assert_eq!(out.memory.records.len(), 3);
    assert_eq!(out.memory.flushes.last(), Some(&3));
    assert!(r.evaluated >= 3 && r.evaluated <= 3 + 8);
}

#[test]
fn a_stop_from_another_thread_ends_a_slow_search() {
    let stop = go();
    let evaluated = AtomicUsize::new(0);
    thread::scope(|scope| {
        scope.spawn(|| {
            wait_until(|| evaluated.load(SeqCst) >= 20);
            stop.store(true, SeqCst);
        });
        let mut out = Memory::default();
        let r = run(
            vec![String::new()],
            &options(3, 6),
            &stop,
            |p| {
                thread::sleep(Duration::from_millis(2));
                evaluated.fetch_add(1, SeqCst);
                Ok(accept(p, 7).then(|| p.to_string()))
            },
            &mut out,
        )
        .unwrap();
        assert_eq!(r.halt, Halt::Stopped);
        assert_eq!(r.evaluated as usize, evaluated.load(SeqCst));
        // The committed records are a prefix of the full search.
        let (expected, _) = reference(7);
        assert_eq!(out.records[..], expected[..out.records.len()]);
    });
}

// ---------------------------------------------------------------------------
// Limits and recovery
// ---------------------------------------------------------------------------

#[test]
fn decision_and_depth_limits_leave_refused_boxes_unresolved() {
    for budget in 1..20 {
        let r = run::<String, _, _>(
            vec![String::new()],
            &Options {
                max_decisions: Some(budget),
                ..options(4, 8)
            },
            &go(),
            |_| Ok(None),
            &mut Memory::default(),
        )
        .unwrap();
        assert_eq!(
            (r.halt, r.decisions, r.unresolved),
            (Halt::DecisionLimit, budget, budget + 1)
        );
        assert!(r.evaluated <= budget + 8);
    }
    for depth in 0..6 {
        let r = run::<String, _, _>(
            vec![String::new()],
            &Options {
                depth_limit: depth,
                ..options(4, 8)
            },
            &go(),
            |_| Ok(None),
            &mut Memory::default(),
        )
        .unwrap();
        assert_eq!(
            (r.halt, r.unresolved, r.decisions),
            (Halt::DepthLimit, 1 << depth, (1 << (depth + 1)) - 1)
        );
    }
    // Frontier paths beyond the limit are never evaluated.
    let r = run::<String, _, _>(
        vec!["0".into(), "10".into(), "11".into()],
        &Options {
            depth_limit: 1,
            ..options(2, 2)
        },
        &go(),
        |p| Ok(Some(p.to_string())),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!((r.halt, r.records, r.unresolved), (Halt::DepthLimit, 1, 2));
}

#[test]
fn a_budget_that_ends_a_layer_starts_no_further_box() {
    let r = run::<String, _, _>(
        vec![String::new()],
        &Options {
            max_decisions: Some(3),
            ..options(4, 8)
        },
        &go(),
        |_| Ok(None),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!((r.halt, r.evaluated, r.decisions, r.unresolved), (Halt::DecisionLimit, 3, 3, 4));
}

#[test]
fn resuming_from_saved_records_alone_reproduces_the_uninterrupted_search() {
    for seed in 0..24 {
        let mut out = Memory::default();
        // Refusals are not saved, so a budget must grow to reach new records.
        for attempt in 1..200 {
            let frontier = complement(&out.records);
            let config = Options {
                max_decisions: Some(2 * attempt + attempt % 3),
                after: out.records.last().map(|(p, _)| p.clone()),
                ..options(1 + attempt as usize % 4, 1 + attempt as usize % 5)
            };
            let r = run(
                frontier,
                &config,
                &go(),
                |p| Ok(accept(p, seed).then(|| p.to_string())),
                &mut out,
            )
            .unwrap();
            if r.halt == Halt::Complete {
                break;
            }
        }
        assert_eq!(out.records, reference(seed).0, "seed {seed}");
        assert!(complement(&out.records).is_empty());
    }
}

#[test]
fn recovery_splits_earlier_boxes_without_evaluating_them() {
    let mut out = Memory::default();
    let r = run(
        vec!["0".into()],
        &Options {
            after: Some("111".into()),
            ..options(2, 4)
        },
        &go(),
        |p| {
            assert!(search_order(p) > search_order("111"));
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!((r.reconstructed, r.evaluated, r.records), (7, 8, 8));
    assert!(out.records.iter().all(|(p, _)| p.len() == 4));
}

#[test]
fn recovery_examples_from_minimal_complements() {
    for (saved, expected, count) in [
        ("00", vec!["01", "10", "11"], 1),
        ("10", vec!["11", "000", "001", "010", "011"], 3),
    ] {
        let old = vec![(saved.to_string(), saved.to_string())];
        let mut out = Memory::default();
        let r = run(
            complement(&old),
            &Options {
                after: Some(saved.into()),
                ..options(3, 3)
            },
            &go(),
            |p| Ok(Some(p.into())),
            &mut out,
        )
        .unwrap();
        assert_eq!(out.paths(), expected);
        assert_eq!(r.reconstructed, count);
    }
}

#[test]
fn recovery_leaves_earlier_boxes_at_the_depth_limit_unresolved() {
    let mut out = Memory::default();
    let r = run(
        vec!["00".into(), "1".into()],
        &Options {
            after: Some("01".into()),
            depth_limit: 2,
            ..options(3, 3)
        },
        &go(),
        |p| {
            assert!(["10", "11"].contains(&p));
            Ok(Some(p.into()))
        },
        &mut out,
    )
    .unwrap();
    assert_eq!(out.paths(), ["10", "11"]);
    assert_eq!(
        (r.halt, r.reconstructed, r.evaluated, r.unresolved),
        (Halt::DepthLimit, 1, 2, 1)
    );
    let r = run::<String, _, _>(
        vec!["00".into(), "1".into()],
        &Options {
            after: Some("01".into()),
            depth_limit: 1,
            ..options(3, 3)
        },
        &go(),
        |_| panic!("nothing may be evaluated at this limit"),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!(
        (r.halt, r.reconstructed, r.evaluated, r.unresolved),
        (Halt::DepthLimit, 0, 0, 2)
    );
    let r = run::<String, _, _>(
        vec!["0".repeat(4096)],
        &Options {
            after: Some("1".repeat(4096)),
            depth_limit: 4096,
            ..options(2, 2)
        },
        &go(),
        |_| panic!("an earlier box at the limit stays unresolved"),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!(
        (r.halt, r.reconstructed, r.evaluated, r.unresolved),
        (Halt::DepthLimit, 0, 0, 1)
    );
}

/// Measures the per-box overhead of the queue with a trivial evaluator. Run
/// with `cargo test --release --lib search::queue::tests::measure -- --ignored
/// --nocapture`.
#[test]
#[ignore]
fn measure() {
    struct Sink(u64);
    impl Output<()> for Sink {
        fn append(&mut self, _: &str, _: &()) -> Result<(), BoxError> {
            self.0 += 1;
            Ok(())
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            Ok(())
        }
    }
    for threads in [1, 3] {
        let started = Instant::now();
        let r = run(
            vec![String::new()],
            &options(threads, 64 * threads),
            &go(),
            |p| Ok((p.len() >= 20).then_some(())),
            &mut Sink(0),
        )
        .unwrap();
        let seconds = started.elapsed().as_secs_f64();
        println!(
            "{threads} threads: {} boxes in {seconds:.2} s, {:.2} us per box",
            r.decisions,
            seconds * 1e6 / r.decisions as f64
        );
    }
}

mod adversarial;

/// Utilisation and ordering with manually timed evaluators.
mod utilisation;
