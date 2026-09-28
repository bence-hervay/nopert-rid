//! The radii a completeness schedule tests. One rule serves the runner,
//! which asks for the next exponent, and the checker, which replays a
//! recorded schedule.
use std::fmt;

/// The largest finest exponent a schedule may use.
pub const MAX_K: u32 = 2048;

/// Recorded trials that do not follow the schedule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ScheduleError {
    /// Trial `position` tests `found` where the schedule tests `expected`.
    Mismatch { position: usize, expected: u32, found: u32 },
    /// The schedule ends before trial `position`.
    Trailing { position: usize },
}

impl fmt::Display for ScheduleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScheduleError::Mismatch {
                position,
                expected,
                found,
            } => write!(f, "trial {position} tests k={found}, the schedule tests k={expected}"),
            ScheduleError::Trailing { position } => {
                write!(f, "the schedule ends before trial {position}")
            }
        }
    }
}

/// Why a replay stopped early.
enum Stop {
    /// The recorded trials end; this exponent is tested next.
    Next(u32),
    Mismatch(ScheduleError),
}

/// Replays recorded `(k, success)` trials against the schedule.
struct Replay<'a> {
    trials: &'a [(u32, bool)],
    position: usize,
}

impl Replay<'_> {
    /// The outcome of the schedule's next trial, which must test `k`.
    fn test(&mut self, k: u32) -> Result<bool, Stop> {
        let Some(&(found, success)) = self.trials.get(self.position) else {
            return Err(Stop::Next(k));
        };
        if found != k {
            return Err(Stop::Mismatch(ScheduleError::Mismatch {
                position: self.position,
                expected: k,
                found,
            }));
        }
        self.position += 1;
        Ok(success)
    }

    fn tested(&self, k: u32) -> bool {
        self.trials[..self.position].iter().any(|&(found, _)| found == k)
    }
}

/// The schedule, as a program over the outcomes: test the finest exponent
/// `max_k` and stop if it fails; test `0`; bisect between the largest
/// refused and the smallest accepted exponent; then test the untested
/// exponents 1, 2 and 8 finer than the accepted one. Nothing is inferred
/// about untested exponents.
fn program(max_k: u32, replay: &mut Replay) -> Result<(), Stop> {
    if !replay.test(max_k)? {
        return Ok(());
    }
    let mut accepted = if replay.test(0)? { 0 } else { max_k };
    let mut refused = 0;
    while accepted - refused > 1 {
        let middle = refused + (accepted - refused) / 2;
        if replay.test(middle)? {
            accepted = middle;
        } else {
            refused = middle;
        }
    }
    for offset in [1, 2, 8] {
        let k = accepted + offset;
        if k <= max_k && !replay.tested(k) {
            replay.test(k)?;
        }
    }
    Ok(())
}

/// The next exponent to test after `trials` (`None` when the schedule is
/// complete), or the first place where `trials` leave the schedule.
pub fn next(max_k: u32, trials: &[(u32, bool)]) -> Result<Option<u32>, ScheduleError> {
    assert!((1..=MAX_K).contains(&max_k), "max_k {max_k} out of range");
    let mut replay = Replay { trials, position: 0 };
    match program(max_k, &mut replay) {
        Err(Stop::Next(k)) => Ok(Some(k)),
        Err(Stop::Mismatch(error)) => Err(error),
        Ok(()) if replay.position == trials.len() => Ok(None),
        Ok(()) => Err(ScheduleError::Trailing {
            position: replay.position,
        }),
    }
}

/// The smallest accepted exponent (the largest radius) and the refused
/// exponents finer than it, in increasing order.
pub fn summary(trials: &[(u32, bool)]) -> (Option<u32>, Vec<u32>) {
    let best = trials.iter().filter(|t| t.1).map(|t| t.0).min();
    let mut finer: Vec<u32> = match best {
        None => Vec::new(),
        Some(best) => trials
            .iter()
            .filter(|&&(k, success)| !success && k > best)
            .map(|t| t.0)
            .collect(),
    };
    finer.sort_unstable();
    (best, finer)
}
