//! Utilisation and ordering with manually timed evaluators.
//!
//! Every box sleeps for a duration fixed by its path (slow heads of layers,
//! bursts of slow boxes, many short boxes), and the decisions come from a
//! fixed rule. The queue must commit the records of the sequential
//! reference for every thread count and window, and its wall-clock time
//! must match the *greedy schedule*: workers take the tickets of a layer in
//! order, each as soon as some worker is free and the window allows (ticket
//! `i` may start once ticket `i − window` is committed), and a layer ends
//! when its last ticket is committed. The greedy schedule is the earliest
//! any queue with in-order tickets and this window can finish, so the
//! measured time must lie between it and it plus a small overhead. A queue
//! that let a worker idle while work and window allowed (batches, waiting
//! for commits, dispatch by the coordinator) is slower by far more than the
//! overhead: each profile is checked to separate the greedy schedule from
//! the schedule in batches of `threads` tickets by a wide margin.
use super::*;
use std::time::Duration;

/// Durations in milliseconds of the boxes of a layer, by ticket.
type Profile = fn(ticket: usize, depth: usize) -> u64;

fn slow_heads(ticket: usize, _depth: usize) -> u64 {
    if ticket == 0 {
        60
    } else {
        8
    }
}

fn bursts(ticket: usize, depth: usize) -> u64 {
    if (ticket + depth) % 8 == 0 {
        24
    } else {
        2
    }
}

fn many_short(ticket: usize, depth: usize) -> u64 {
    1 + ((ticket * 7 + depth) % 3) as u64
}

/// The decision rule: a third of the boxes eliminated, all at depth 8.
fn eliminated(path: &str, seed: u64) -> bool {
    path.len() >= 8 || hash(path, seed) % 3 == 0
}

/// The layers of the sequential reference from the 16 boxes of depth 4:
/// the paths of every depth in search order.
fn layers(seed: u64) -> Vec<Vec<String>> {
    let mut layers = Vec::new();
    let mut layer = paths(4);
    while !layer.is_empty() {
        let next = layer.iter().filter(|p| !eliminated(p, seed)).flat_map(|p| [format!("{p}0"), format!("{p}1")]).collect();
        layers.push(layer);
        layer = next;
    }
    layers
}

/// The greedy schedule's time for one layer of these durations.
fn greedy(durations: &[u64], threads: usize, window: usize) -> u64 {
    let mut free = vec![0u64; threads];
    let mut commit: Vec<u64> = Vec::with_capacity(durations.len());
    for (i, d) in durations.iter().enumerate() {
        let worker = (0..threads).min_by_key(|&w| free[w]).unwrap();
        let allowed = if i >= window { commit[i - window] } else { 0 };
        let finish = free[worker].max(allowed) + d;
        free[worker] = finish;
        commit.push(finish.max(commit.last().copied().unwrap_or(0)));
    }
    commit.last().copied().unwrap_or(0)
}

/// The time in batches of `threads` tickets, each batch waiting for the
/// previous one to finish: what an idling queue would need.
fn batched(durations: &[u64], threads: usize) -> u64 {
    durations.chunks(threads).map(|batch| batch.iter().max().copied().unwrap_or(0)).sum()
}

fn run_timed(seed: u64, profile: Profile, threads: usize, window: usize) -> (Vec<(String, String)>, Duration) {
    let depth_of_first: Vec<String> = layers(seed).iter().map(|l| l[0].clone()).collect();
    let tickets: std::collections::HashMap<String, usize> = layers(seed)
        .iter()
        .flat_map(|layer| layer.iter().enumerate().map(|(i, p)| (p.clone(), i)))
        .collect();
    assert!(!depth_of_first.is_empty());
    let mut out = Memory::default();
    let start = Instant::now();
    let report = run(
        paths(4),
        &options(threads, window),
        &go(),
        |p| {
            thread::sleep(Duration::from_millis(profile(tickets[p], p.len())));
            Ok(eliminated(p, seed).then(|| p.to_string()))
        },
        &mut out,
    )
    .unwrap();
    let elapsed = start.elapsed();
    assert_eq!(report.halt, Halt::Complete);
    (out.records, elapsed)
}

#[test]
fn timed_boxes_keep_every_worker_busy_and_commit_the_reference() {
    let seed = 5;
    let layers = layers(seed);
    let expected: Vec<(String, String)> =
        layers.iter().flatten().filter(|p| eliminated(p, seed)).map(|p| (p.clone(), p.clone())).collect();
    assert!(layers.len() >= 5 && layers.iter().map(Vec::len).sum::<usize>() > 100, "{:?}", layers.iter().map(Vec::len).collect::<Vec<_>>());
    let profiles: [(&str, Profile); 3] = [("slow heads", slow_heads), ("bursts", bursts), ("many short", many_short)];
    for (name, profile) in profiles {
        for (threads, window) in [(1, 4), (4, 32), (8, 64), (16, 64), (32, 128)] {
            let per_layer: Vec<Vec<u64>> =
                layers.iter().map(|l| (0..l.len()).map(|i| profile(i, l[0].len())).collect()).collect();
            let ideal: u64 = per_layer.iter().map(|d| greedy(d, threads, window)).sum();
            // Two idling queues: tickets in batches of `threads`, and workers
            // that wait for commits (a window of only `threads` tickets).
            let batches: u64 = per_layer.iter().map(|d| batched(d, threads)).sum();
            let waiting: u64 = per_layer.iter().map(|d| greedy(d, threads, threads)).sum();
            let (records, elapsed) = run_timed(seed, profile, threads, window);
            assert_eq!(records, expected, "{name}, {threads} threads, window {window}");
            let ms = elapsed.as_secs_f64() * 1000.0;
            // Sleeping never ends early; allow 2 ms of overhead per layer and
            // 10% on top of the greedy schedule.
            let bound = ideal as f64 * 1.1 + 2.0 * layers.len() as f64 + 20.0;
            eprintln!(
                "{name}, {threads} threads, window {window}: {ms:.0} ms; greedy {ideal}, bound {bound:.0}, batches {batches}, waiting {waiting} ms"
            );
            assert!(ms >= ideal as f64 * 0.99, "{name}: {ms:.1} ms is faster than the greedy schedule {ideal} ms");
            assert!(ms <= bound, "{name}, {threads} threads, window {window}: {ms:.1} ms, greedy {ideal} ms");
            // The profiles tell a busy queue from the idling ones.
            match (name, threads) {
                ("bursts", 4 | 8 | 16) => assert!(bound < batches as f64, "{name}: batches {batches} ms within {bound:.0}"),
                ("slow heads", 4 | 8) => assert!(bound < waiting as f64, "{name}: waiting {waiting} ms within {bound:.0}"),
                _ => {}
            }
        }
    }
}

#[test]
fn the_greedy_schedule_respects_the_window_and_the_order() {
    // A stalled head fills exactly the window; later tickets wait for its
    // commit even with free workers.
    assert_eq!(greedy(&[10, 1, 1, 1, 1], 4, 3), 11);
    assert_eq!(greedy(&[10, 1, 1, 1], 4, 4), 10);
    // One worker: the sum.
    assert_eq!(greedy(&[3, 4, 5], 1, 8), 12);
    // Enough workers and window: the longest box, committed in order.
    assert_eq!(greedy(&[1, 9, 1, 1], 8, 8), 9);
    assert_eq!(batched(&[1, 9, 1, 1], 2), 10);
}
