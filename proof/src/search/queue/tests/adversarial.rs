//! Adversarial tests of the queue: randomised
//! failures, stops and depth-limited resumptions.

use super::*;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicUsize, Ordering::SeqCst};
use std::sync::mpsc;
use std::time::Instant;

fn nonzero(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn mix(path: &str, seed: u64) -> u64 {
    path.bytes().fold(seed ^ 0x2545_f491_4f6c_dd1d, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
    })
}

/// Accepts a box with probability about 1/3 from depth 2 on.
fn accepts(path: &str, seed: u64) -> bool {
    path.len() >= 2 && mix(path, seed) % 3 == 0
}

/// Sequential breadth-first reference with a depth limit: every decision in
/// order (path, accepted), and the unresolved leaves at the limit.
fn reference(seed: u64, depth_limit: usize) -> (Vec<(String, bool)>, u64) {
    let mut queue = VecDeque::from([String::new()]);
    let mut decisions = Vec::new();
    let mut unresolved = 0;
    while let Some(path) = queue.pop_front() {
        let accepted = accepts(&path, seed);
        decisions.push((path.clone(), accepted));
        if !accepted {
            if path.len() < depth_limit {
                queue.push_back(format!("{path}0"));
                queue.push_back(format!("{path}1"));
            } else {
                unresolved += 1;
            }
        }
    }
    (decisions, unresolved)
}

fn records_of(decisions: &[(String, bool)]) -> Vec<String> {
    decisions.iter().filter(|(_, a)| *a).map(|(p, _)| p.clone()).collect()
}

/// Independent minimal complement of a prefix-free record set (recursive).
fn complement(records: &[String]) -> Vec<String> {
    fn visit(path: String, records: &[String], out: &mut Vec<String>) {
        if records.iter().any(|p| *p == path) {
            return;
        }
        if records.iter().any(|p| p.starts_with(&path)) {
            visit(format!("{path}0"), records, out);
            visit(format!("{path}1"), records, out);
        } else {
            out.push(path);
        }
    }
    let mut out = Vec::new();
    visit(String::new(), records, &mut out);
    out
}

#[derive(Default)]
struct Memory {
    records: Vec<String>,
    stop_after: Option<(usize, &'static AtomicBool)>,
}

impl Output<String> for Memory {
    fn append(&mut self, path: &str, data: &String) -> Result<(), BoxError> {
        assert_eq!(path, data);
        if let Some(last) = self.records.last() {
            assert!(search_order(last) < search_order(path));
        }
        self.records.push(path.to_owned());
        if let Some((n, stop)) = self.stop_after {
            if self.records.len() == n {
                stop.store(true, SeqCst);
            }
        }
        Ok(())
    }
    fn flush(&mut self) -> Result<(), BoxError> {
        Ok(())
    }
}

fn jitter(state: u64) {
    match state % 6 {
        0 => thread::sleep(Duration::from_micros(state % 300)),
        1 => thread::yield_now(),
        _ => {}
    }
}

#[test]
fn an_evaluator_error_anywhere_commits_exactly_the_decisions_before_it() {
    let mut state = 0xfeed_beef_0123_4567;
    for round in 0..150 {
        let seed = random(&mut state) % 1000;
        let depth_limit = 3 + (random(&mut state) as usize % 6);
        let (decisions, _) = reference(seed, depth_limit);
        let at = random(&mut state) as usize % decisions.len();
        let bad = decisions[at].0.clone();
        let panics = round % 2 == 1;
        let threads = 1 + random(&mut state) as usize % 6;
        let window = 1 + random(&mut state) as usize % 12;
        let mut out = Memory::default();
        let calls = AtomicUsize::new(0);
        let error = run(
            vec![String::new()],
            &Options {
                threads: nonzero(threads),
                window: nonzero(window),
                depth_limit,
                max_decisions: None,
                after: None,
            },
            &AtomicBool::new(false),
            |p: &str| {
                jitter(mix(p, calls.fetch_add(1, SeqCst) as u64));
                if p == bad {
                    if panics {
                        panic!("review panic");
                    }
                    return Err("review error".into());
                }
                Ok(accepts(p, seed).then(|| p.to_owned()))
            },
            &mut out,
        )
        .unwrap_err();
        match error {
            Error::Evaluator { path, .. } | Error::EvaluatorPanic { path, .. } => {
                assert_eq!(path, bad)
            }
            other => panic!("{other}"),
        }
        assert_eq!(out.records, records_of(&decisions[..at]), "round {round}");
    }
}

#[test]
fn depth_limited_resumptions_reproduce_the_uninterrupted_search() {
    let mut state = 0x0dd0_c0de_7777_1111;
    for round in 0..60 {
        let seed = random(&mut state) % 10_000;
        let depth_limit = 2 + (random(&mut state) as usize % 7);
        let (decisions, unresolved) = reference(seed, depth_limit);
        let expected = records_of(&decisions);
        let mut out = Memory::default();
        let mut last = None;
        for attempt in 1..500u64 {
            let frontier = complement(&out.records);
            let budget = 1 + random(&mut state) % (3 * attempt);
            let r = run(
                frontier,
                &Options {
                    threads: nonzero(1 + random(&mut state) as usize % 4),
                    window: nonzero(1 + random(&mut state) as usize % 6),
                    depth_limit,
                    max_decisions: Some(budget),
                    after: out.records.last().cloned(),
                },
                &AtomicBool::new(false),
                |p: &str| Ok(accepts(p, seed).then(|| p.to_owned())),
                &mut out,
            )
            .unwrap();
            if matches!(r.halt, Halt::Complete | Halt::DepthLimit) {
                last = Some(r);
                break;
            }
            assert_eq!(r.halt, Halt::DecisionLimit);
        }
        let last = last.expect("no attempt finished");
        assert_eq!(out.records, expected, "round {round} seed {seed} depth {depth_limit}");
        assert_eq!(last.unresolved, unresolved, "round {round}");
        assert_eq!(
            last.halt,
            if unresolved == 0 { Halt::Complete } else { Halt::DepthLimit }
        );
        // A further run finds nothing new and evaluates only the boxes after
        // the last record.
        let again = run(
            complement(&out.records),
            &Options {
                threads: nonzero(2),
                window: nonzero(3),
                depth_limit,
                max_decisions: None,
                after: out.records.last().cloned(),
            },
            &AtomicBool::new(false),
            |p: &str| {
                if let Some(last) = expected.last() {
                    assert!(search_order(p) > search_order(last), "{p} precedes {last}");
                }
                Ok(accepts(p, seed).then(|| p.to_owned()))
            },
            &mut Memory::default(),
        )
        .unwrap();
        assert_eq!((again.records, again.unresolved), (0, unresolved));
    }
}

#[test]
fn random_stops_commit_a_prefix_and_resume_to_the_reference() {
    let mut state = 0x5151_7272_9393_abab;
    for round in 0..40 {
        let seed = random(&mut state) % 10_000;
        let (decisions, _) = reference(seed, 10);
        let expected = records_of(&decisions);
        let mut out = Memory::default();
        let mut runs = 0;
        loop {
            runs += 1;
            assert!(runs < 100);
            let stop: &'static AtomicBool = Box::leak(Box::new(AtomicBool::new(false)));
            out.stop_after = Some((out.records.len() + 1 + random(&mut state) as usize % 20, stop));
            let calls = AtomicUsize::new(0);
            let r = run(
                complement(&out.records),
                &Options {
                    threads: nonzero(1 + random(&mut state) as usize % 5),
                    window: nonzero(1 + random(&mut state) as usize % 9),
                    depth_limit: 10,
                    max_decisions: None,
                    after: out.records.last().cloned(),
                },
                stop,
                |p: &str| {
                    jitter(mix(p, calls.fetch_add(1, SeqCst) as u64));
                    Ok(accepts(p, seed).then(|| p.to_owned()))
                },
                &mut out,
            )
            .unwrap();
            assert!(expected.starts_with(&out.records), "round {round}");
            if r.halt != Halt::Stopped {
                break;
            }
        }
        assert_eq!(out.records, expected, "round {round}");
    }
}

/// Every evaluation fails: the run must end promptly with the first box's
/// error, for every thread count and window.
#[test]
fn failures_on_every_box_end_the_run_promptly() {
    for threads in [1, 2, 8] {
        for window in [1, 3, 64] {
            for panics in [false, true] {
                let (send, receive) = mpsc::channel();
                thread::spawn(move || {
                    let frontier: Vec<String> =
                        (0..64).map(|n| format!("{n:06b}")).collect();
                    let result = run::<String, _, _>(
                        frontier,
                        &Options {
                            threads: nonzero(threads),
                            window: nonzero(window),
                            depth_limit: 64,
                            max_decisions: None,
                            after: None,
                        },
                        &AtomicBool::new(false),
                        |p: &str| {
                            if panics {
                                panic!("every box panics");
                            }
                            Err(format!("fail {p}").into())
                        },
                        &mut Memory::default(),
                    );
                    let _ = send.send(match result {
                        Err(Error::Evaluator { path, .. } | Error::EvaluatorPanic { path, .. }) => {
                            path
                        }
                        _ => "unexpected".into(),
                    });
                });
                let path = receive
                    .recv_timeout(Duration::from_secs(20))
                    .expect("the run did not end");
                assert_eq!(path, "000000");
            }
        }
    }
}

#[test]
fn a_failing_final_flush_is_an_error_even_without_work() {
    struct Refuse;
    impl Output<String> for Refuse {
        fn append(&mut self, _: &str, _: &String) -> Result<(), BoxError> {
            Ok(())
        }
        fn flush(&mut self) -> Result<(), BoxError> {
            Err("no flush".into())
        }
    }
    let options = Options {
        threads: nonzero(2),
        window: nonzero(2),
        depth_limit: 4,
        max_decisions: None,
        after: None,
    };
    for (frontier, stopped) in [(vec![], false), (vec![String::new()], true)] {
        let r = run::<String, _, _>(
            frontier,
            &options,
            &AtomicBool::new(stopped),
            |_| panic!("no work"),
            &mut Refuse,
        );
        assert!(matches!(r, Err(Error::Output(_))));
    }
}

#[test]
fn a_budget_spent_on_the_last_box_reports_complete() {
    // Root refused, both children accepted: exactly three decisions.
    let r = run(
        vec![String::new()],
        &Options {
            threads: nonzero(3),
            window: nonzero(3),
            depth_limit: 4,
            max_decisions: Some(3),
            after: None,
        },
        &AtomicBool::new(false),
        |p: &str| Ok((!p.is_empty()).then(|| p.to_owned())),
        &mut Memory::default(),
    )
    .unwrap();
    assert_eq!((r.halt, r.decisions, r.records, r.unresolved), (Halt::Complete, 3, 2, 0));
}

/// A documented resource limit: a certificate that an ordered search could
/// not have written (one deep record) makes reconstruction exponential in the
/// record's depth; only a stop ends it early. Counted here at a modest depth.
#[test]
fn reconstruction_after_a_lone_deep_record_splits_every_earlier_box() {
    let saved = "1".repeat(18);
    let started = Instant::now();
    let r = run(
        complement(&[saved.clone()]),
        &Options {
            threads: nonzero(1),
            window: nonzero(1),
            depth_limit: 64,
            max_decisions: Some(1),
            after: Some(saved),
        },
        &AtomicBool::new(false),
        |p: &str| Ok(Some(p.to_owned())),
        &mut Memory::default(),
    )
    .unwrap();
    eprintln!(
        "one record at depth 18: {} reconstructed splits in {:.2} s",
        r.reconstructed,
        started.elapsed().as_secs_f64()
    );
    // Every box of depths 1 to 18 except the record's ancestors and itself.
    assert_eq!(r.reconstructed, (1 << 19) - 20);
}
