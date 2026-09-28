use super::*;
use crate::testing::{Rng, Scratch};
use serde::Deserialize;
use std::sync::atomic::AtomicU64;
use std::time::Duration;

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

#[test]
fn results_are_committed_in_item_order_whatever_the_delays() {
    let mut rng = Rng::new(11);
    for round in 0..24 {
        let n = rng.below(40) as usize;
        let delays: Vec<u64> = (0..n).map(|_| rng.below(4)).collect();
        let evaluated: Vec<AtomicU64> = (0..n).map(|_| AtomicU64::new(0)).collect();
        let mut committed = Vec::new();
        let result: Result<(), ()> = ordered(
            &delays,
            threads(1 + round % 5),
            |index, delay| {
                evaluated[index].fetch_add(1, Ordering::Relaxed);
                thread::sleep(Duration::from_millis(*delay));
                Ok(index * 10)
            },
            |index, value| {
                assert_eq!(value, index * 10);
                committed.push(index);
                Ok(())
            },
        );
        assert_eq!(result, Ok(()));
        assert_eq!(committed, (0..n).collect::<Vec<_>>());
        assert!(evaluated.iter().all(|c| c.load(Ordering::Relaxed) == 1));
    }
}

#[test]
fn the_earliest_failing_item_stops_the_run_after_its_predecessors() {
    for workers in 1..=4 {
        for failing in [0usize, 3, 17] {
            let items: Vec<usize> = (0..30).collect();
            let mut committed = Vec::new();
            let result = ordered(
                &items,
                threads(workers),
                |index, _| {
                    // A later item fails first in time; the earlier one wins.
                    if index == failing + 2 {
                        return Err(index);
                    }
                    thread::sleep(Duration::from_millis(if index == failing { 20 } else { 0 }));
                    if index == failing {
                        Err(index)
                    } else {
                        Ok(index)
                    }
                },
                |index, _| {
                    committed.push(index);
                    Ok(())
                },
            );
            assert_eq!(result, Err(failing));
            assert_eq!(committed, (0..failing).collect::<Vec<_>>());
        }
    }
}

#[test]
fn a_failing_commit_stops_the_run() {
    let items: Vec<u32> = (0..50).collect();
    let evaluated = AtomicU64::new(0);
    let mut committed = 0;
    let result = ordered(
        &items,
        threads(3),
        |_, _| {
            evaluated.fetch_add(1, Ordering::Relaxed);
            thread::sleep(Duration::from_millis(2));
            Ok(())
        },
        |index, ()| {
            if index == 5 {
                return Err("full disk");
            }
            committed += 1;
            Ok(())
        },
    );
    assert_eq!(result, Err("full disk"));
    assert_eq!(committed, 5);
    assert!(evaluated.load(Ordering::Relaxed) < 50, "the workers stopped early");
}

#[test]
fn a_panic_propagates_and_stops_the_other_workers() {
    let items: Vec<u32> = (0..100).collect();
    let evaluated = AtomicU64::new(0);
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        ordered(
            &items,
            threads(2),
            |index, _| {
                evaluated.fetch_add(1, Ordering::Relaxed);
                thread::sleep(Duration::from_millis(1));
                if index == 3 {
                    panic!("evaluator bug");
                }
                Ok::<_, ()>(())
            },
            |_, ()| Ok(()),
        )
    }));
    assert!(outcome.is_err());
    assert!(evaluated.load(Ordering::Relaxed) < 100);
}

#[test]
fn nothing_to_do_is_done_at_once() {
    let items: Vec<u8> = Vec::new();
    let result: Result<(), ()> = ordered(&items, threads(4), |_, _| Ok(()), |_, ()| panic!("no commit"));
    assert_eq!(result, Ok(()));
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Line {
    name: String,
    count: u32,
}

#[test]
fn transcripts_are_created_once_and_read_back_exactly() {
    let dir = Scratch::new("writer");
    let path = dir.join("t.jsonl");
    let mut writer = Writer::create(&path).unwrap();
    writer.write(&Line { name: "a".into(), count: 1 }).unwrap();
    writer.write(&Line { name: "b".into(), count: 2 }).unwrap();
    assert!(matches!(Writer::create(&path), Err(TranscriptError::Io { .. })));
    let bytes = read(&path).unwrap();
    assert_eq!(bytes, b"{\"name\":\"a\",\"count\":1}\n{\"name\":\"b\",\"count\":2}\n");
    let lines = lines(&bytes).unwrap();
    assert_eq!(parse::<Line>(lines[1], 2).unwrap(), Line { name: "b".into(), count: 2 });
}

#[test]
fn incomplete_and_non_canonical_lines_are_refused() {
    for bad in [&b""[..], b"{}", b"{}\n\n{}\n", b"\n"] {
        assert!(matches!(lines(bad), Err(TranscriptError::Incomplete)), "{bad:?}");
    }
    for (bad, canonical) in [
        (&br#"{"name":"a","count":1}"#[..], true),
        (br#"{"count":1,"name":"a"}"#, false),
        (br#"{"name": "a","count":1}"#, false),
        (br#"{"name":"a","count":1.0}"#, false),
        (b"{\"name\":\"\\u0061\",\"count\":1}", false),
        (br#"{"name":"a","count":1,"extra":0}"#, false),
        (br#"{"name":"a"}"#, false),
    ] {
        assert_eq!(parse::<Line>(bad, 7).is_ok(), canonical, "{}", String::from_utf8_lossy(bad));
    }
}
