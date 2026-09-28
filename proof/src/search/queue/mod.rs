//! Parallel evaluation of search boxes with commits in a fixed order:
//! increasing depth, then path. Knows binary paths, an evaluator and an
//! output; no geometry and no file format.

use super::panic_message;
pub use super::BoxError;
use crate::problem::configuration::{is_path, search_order};
use std::collections::BTreeMap;
use std::fmt;
use std::num::NonZeroUsize;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread;
use std::time::Duration;

/// How often the coordinator looks at the stop flag while it waits.
const POLL: Duration = Duration::from_millis(10);

/// Receives committed records, only on the calling thread and in strictly
/// increasing (depth, path) order. `flush` makes every earlier append durable.
pub trait Output<D> {
    fn append(&mut self, path: &str, data: &D) -> Result<(), BoxError>;
    fn flush(&mut self) -> Result<(), BoxError>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    /// Worker threads evaluating boxes.
    pub threads: NonZeroUsize,
    /// Tickets assigned but not yet committed, running or finished.
    pub window: NonZeroUsize,
    /// Boxes at this depth are evaluated but not split; deeper ones are not
    /// evaluated. Both stay unresolved.
    pub depth_limit: usize,
    /// Commit at most this many decisions in this run.
    pub max_decisions: Option<u64>,
    /// The last record already saved: unresolved boxes at or before it in
    /// the search order are split without evaluation.
    pub after: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Halt {
    /// No unresolved box remains.
    Complete,
    /// The stop flag was observed.
    Stopped,
    /// `max_decisions` decisions were committed.
    DecisionLimit,
    /// Only boxes at or beyond `depth_limit` remain unresolved.
    DepthLimit,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Report {
    pub halt: Halt,
    /// Tickets handed to workers, including results never committed.
    pub evaluated: u64,
    /// Committed decisions: records plus refusals.
    pub decisions: u64,
    pub records: u64,
    /// Refusals that replaced a box by its two children.
    pub splits: u64,
    /// Recovery splits at or before `after`, made without evaluation.
    pub reconstructed: u64,
    /// Unresolved boxes after the committed decisions.
    pub unresolved: u64,
}

#[derive(Debug)]
pub enum Error {
    /// A frontier path or `after` is not a binary string.
    InvalidPath(String),
    /// Two frontier paths are equal or one contains the other.
    Overlap(String, String),
    /// The evaluator returned an error.
    Evaluator { path: String, source: BoxError },
    /// The evaluator panicked.
    EvaluatorPanic { path: String, message: String },
    /// The output refused an append or a flush.
    Output(BoxError),
    /// A worker thread could not be started.
    Spawn(std::io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidPath(path) => write!(f, "path {path:?} is not binary"),
            Error::Overlap(a, b) => write!(f, "frontier paths {a:?} and {b:?} overlap"),
            Error::Evaluator { path, source } => write!(f, "evaluating {path:?} failed: {source}"),
            Error::EvaluatorPanic { path, message } => {
                write!(f, "evaluating {path:?} panicked: {message}")
            }
            Error::Output(source) => write!(f, "output failed: {source}"),
            Error::Spawn(source) => write!(f, "cannot start a worker: {source}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Evaluator { source, .. } | Error::Output(source) => Some(source.as_ref()),
            Error::Spawn(source) => Some(source),
            _ => None,
        }
    }
}

/// Groups a prefix-free frontier by depth.
fn layers(mut frontier: Vec<String>) -> Result<BTreeMap<usize, Vec<String>>, Error> {
    if let Some(path) = frontier.iter().find(|p| !is_path(p)) {
        return Err(Error::InvalidPath(path.clone()));
    }
    // In lexical order a path is immediately followed by its extensions.
    frontier.sort_unstable();
    if let Some(pair) = frontier.windows(2).find(|p| p[1].starts_with(p[0].as_str())) {
        return Err(Error::Overlap(pair[0].clone(), pair[1].clone()));
    }
    let mut layers = BTreeMap::<usize, Vec<String>>::new();
    for path in frontier {
        layers.entry(path.len()).or_default().push(path);
    }
    Ok(layers)
}

/// Splits unresolved boxes at or before `after` without evaluating them,
/// except at `depth_limit`, where they stay unresolved. Returns false if the
/// stop flag was observed.
fn reconstruct(
    pending: &mut BTreeMap<usize, Vec<String>>,
    after: &str,
    depth_limit: usize,
    stop: &AtomicBool,
    report: &mut Report,
) -> bool {
    while let Some((depth, mut layer)) = pending.pop_first() {
        if stop.load(Ordering::Acquire) {
            pending.insert(depth, layer);
            return false;
        }
        if depth > depth_limit {
            pending.insert(depth, layer);
            return true;
        }
        layer.sort_unstable();
        let cut = layer.partition_point(|p| search_order(p) <= search_order(after));
        let later = layer.split_off(cut);
        if !later.is_empty() {
            pending.insert(depth, later);
        }
        if layer.is_empty() {
            // Deeper layers hold only longer paths, which follow `after` too.
            return true;
        }
        if depth == depth_limit {
            // These boxes were refused at the limit; they stay unresolved.
            continue;
        }
        let children = pending.entry(depth + 1).or_default();
        for path in layer {
            // A lone deep record can make this expansion very large.
            if stop.load(Ordering::Acquire) {
                return false;
            }
            children.push(format!("{path}0"));
            children.push(format!("{path}1"));
            report.reconstructed += 1;
            report.unresolved += 1;
        }
    }
    true
}

type Outcome<D> = Result<Option<D>, Error>;

struct State<D> {
    /// The layer being evaluated, sorted; ticket `i` is `layer[i]`.
    layer: Arc<[String]>,
    /// Next ticket to hand out.
    next: usize,
    /// Tickets committed by the coordinator.
    committed: usize,
    /// Finished, not yet committed outcomes (fewer than `window`).
    done: BTreeMap<usize, Outcome<D>>,
    shutdown: bool,
    evaluated: u64,
}

struct Shared<D> {
    state: Mutex<State<D>>,
    /// Workers wait here for a ticket.
    work: Condvar,
    /// The coordinator waits here for the outcome it must commit next.
    results: Condvar,
    window: usize,
}

impl<D> Shared<D> {
    fn lock(&self) -> MutexGuard<'_, State<D>> {
        // Critical sections never panic, but a poisoned lock must not hide
        // the original panic behind a second one.
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn start(&self, layer: Arc<[String]>) {
        let mut state = self.lock();
        debug_assert!(state.done.is_empty() && state.committed == state.layer.len());
        state.layer = layer;
        state.next = 0;
        state.committed = 0;
        drop(state);
        self.work.notify_all();
    }

    /// The outcome of `ticket`, or `None` once the stop flag is observed.
    fn wait_for(&self, ticket: usize, stop: &AtomicBool) -> Option<Outcome<D>> {
        let mut state = self.lock();
        loop {
            if stop.load(Ordering::Acquire) {
                return None;
            }
            if let Some(outcome) = state.done.remove(&ticket) {
                return Some(outcome);
            }
            state = self
                .results
                .wait_timeout(state, POLL)
                .unwrap_or_else(|e| e.into_inner())
                .0;
        }
    }

    fn commit(&self, committed: usize) {
        self.lock().committed = committed;
        self.work.notify_all();
    }

    fn work<F>(&self, evaluate: &F)
    where
        F: Fn(&str) -> Result<Option<D>, BoxError>,
    {
        loop {
            let (ticket, path) = {
                let mut state = self.lock();
                while !state.shutdown
                    && !(state.next < state.layer.len()
                        && state.next - state.committed < self.window)
                {
                    state = self.work.wait(state).unwrap_or_else(|e| e.into_inner());
                }
                if state.shutdown {
                    return;
                }
                let ticket = state.next;
                state.next += 1;
                state.evaluated += 1;
                (ticket, state.layer[ticket].clone())
            };
            let outcome = match catch_unwind(AssertUnwindSafe(|| evaluate(&path))) {
                Ok(Ok(decision)) => Ok(decision),
                Ok(Err(source)) => Err(Error::Evaluator { path, source }),
                Err(panic) => Err(Error::EvaluatorPanic {
                    path,
                    message: panic_message(panic.as_ref()),
                }),
            };
            let failed = outcome.is_err();
            self.lock().done.insert(ticket, outcome);
            self.results.notify_one();
            if failed {
                // The coordinator stops at this ticket; later work is wasted.
                return;
            }
        }
    }
}

/// Releases the workers however the coordinator leaves, including a panic.
struct Shutdown<'a, D>(&'a Shared<D>);

impl<D> Drop for Shutdown<'_, D> {
    fn drop(&mut self) {
        self.0.lock().shutdown = true;
        self.0.work.notify_all();
    }
}

/// Evaluates the unresolved boxes of `frontier`, layer by layer.
///
/// `evaluate` returns `Some(data)` for a box it eliminates (the record) and
/// `None` for a box to split; it must be deterministic and must return. Only
/// this thread calls `output`. Every `Ok` return follows a successful flush.
/// The caller guarantees that `frontier` and the saved records cover the root.
pub fn run<D, F, O>(
    frontier: Vec<String>,
    options: &Options,
    stop: &AtomicBool,
    evaluate: F,
    output: &mut O,
) -> Result<Report, Error>
where
    D: Send,
    F: Fn(&str) -> Result<Option<D>, BoxError> + Sync,
    O: Output<D> + ?Sized,
{
    if let Some(after) = options.after.as_deref().filter(|p| !is_path(p)) {
        return Err(Error::InvalidPath(after.to_owned()));
    }
    let mut report = Report {
        halt: Halt::Complete,
        evaluated: 0,
        decisions: 0,
        records: 0,
        splits: 0,
        reconstructed: 0,
        unresolved: frontier.len() as u64,
    };
    let mut pending = layers(frontier)?;
    let early = if report.unresolved == 0 {
        Some(Halt::Complete)
    } else if stop.load(Ordering::Acquire) {
        Some(Halt::Stopped)
    } else if options.max_decisions == Some(0) {
        Some(Halt::DecisionLimit)
    } else if let Some(after) = options.after.as_deref() {
        (!reconstruct(&mut pending, after, options.depth_limit, stop, &mut report))
            .then_some(Halt::Stopped)
    } else {
        None
    };
    report.halt = match early {
        Some(halt) => halt,
        None => {
            let shared = Shared {
                state: Mutex::new(State {
                    layer: Arc::from(Vec::new()),
                    next: 0,
                    committed: 0,
                    done: BTreeMap::new(),
                    shutdown: false,
                    evaluated: 0,
                }),
                work: Condvar::new(),
                results: Condvar::new(),
                window: options.window.get(),
            };
            let halt = thread::scope(|scope| {
                let _shutdown = Shutdown(&shared);
                for n in 0..options.threads.get() {
                    let (shared, evaluate) = (&shared, &evaluate);
                    thread::Builder::new()
                        .name(format!("box-{n}"))
                        .spawn_scoped(scope, move || shared.work(evaluate))
                        .map_err(Error::Spawn)?;
                }
                coordinate(&shared, pending, options, stop, output, &mut report)
            });
            report.evaluated = shared.lock().evaluated;
            halt?
        }
    };
    output.flush().map_err(Error::Output)?;
    Ok(report)
}

/// The coordinator: publishes each layer, commits outcomes strictly in ticket
/// order and flushes the output before the next layer.
fn coordinate<D, O>(
    shared: &Shared<D>,
    mut pending: BTreeMap<usize, Vec<String>>,
    options: &Options,
    stop: &AtomicBool,
    output: &mut O,
    report: &mut Report,
) -> Result<Halt, Error>
where
    O: Output<D> + ?Sized,
{
    let exhausted = |report: &Report| {
        options
            .max_decisions
            .is_some_and(|cap| report.decisions >= cap)
    };
    while let Some((depth, layer)) = pending.pop_first() {
        if depth > options.depth_limit {
            break;
        }
        if stop.load(Ordering::Acquire) {
            return Ok(Halt::Stopped);
        }
        if exhausted(report) {
            return Ok(Halt::DecisionLimit);
        }
        let mut layer = layer;
        // Recovered paths and new children of one depth arrive unsorted.
        layer.sort_unstable();
        let layer: Arc<[String]> = layer.into();
        shared.start(Arc::clone(&layer));
        for (ticket, path) in layer.iter().enumerate() {
            if exhausted(report) {
                return Ok(Halt::DecisionLimit);
            }
            let Some(outcome) = shared.wait_for(ticket, stop) else {
                return Ok(Halt::Stopped);
            };
            match outcome? {
                Some(data) => {
                    output.append(path, &data).map_err(Error::Output)?;
                    report.records += 1;
                    report.unresolved -= 1;
                }
                None if depth >= options.depth_limit => {}
                None => {
                    let children = pending.entry(depth + 1).or_default();
                    children.push(format!("{path}0"));
                    children.push(format!("{path}1"));
                    report.splits += 1;
                    report.unresolved += 1;
                }
            }
            report.decisions += 1;
            shared.commit(ticket + 1);
        }
        // The depth barrier: nothing deeper starts before this succeeds.
        output.flush().map_err(Error::Output)?;
    }
    Ok(if report.unresolved == 0 {
        Halt::Complete
    } else {
        Halt::DepthLimit
    })
}

#[cfg(test)]
mod tests;
