//! The commands `prepare`, `search` and `check` of `rid <command>
//! <config.json>`: they connect the certificate and the queue to the
//! mathematical components through the [`Components`] trait. `main.rs`
//! installs the signal handlers and calls [`entry`] with the components.
//!

use super::certificate::{self, Header, RecordData, Store};
use super::config::{self, Command};
use super::{queue, BoxError};
use crate::POLICY;
use serde_json::{json, Value};
use std::ffi::OsString;
use std::fmt;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::{Duration, Instant};

/// During a search, an append at least this long after the previous flush
/// flushes the certificate first (besides the flush at every depth barrier).
const FLUSH_INTERVAL: Duration = Duration::from_secs(1);
/// Interval between progress lines on standard error during a search.
const PROGRESS_INTERVAL: Duration = Duration::from_secs(10);

/// Exit status: the summary was printed.
pub const SUCCESS: u8 = 0;
/// Exit status: a refusal or a failure (the reason is on standard error).
pub const FAILURE: u8 = 1;
/// Exit status: a malformed command line or an unknown command.
pub const USAGE_ERROR: u8 = 2;

const USAGE: &str = "usage: rid <prepare|search|check> <config.json>";

/// What the commands need from the mathematical components.
pub trait Components: Sync {
    type Data: RecordData;
    /// A summary of the prepared components, printed by `prepare`.
    fn report(&self) -> Value;
    /// `Some(data)` eliminates the box at `path` and is appended as its
    /// record; `None` refuses it, so that it is split.
    fn check(&self, path: &str) -> Result<Option<Self::Data>, BoxError>;
    /// Whether a record holds: run on every saved record whenever a
    /// certificate is opened or checked.
    fn holds(&self, path: &str, data: &Self::Data) -> Result<(), BoxError>;
}

/// The verifier observed the stop flag.
#[derive(Debug)]
struct Stopped;

impl fmt::Display for Stopped {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "stopped by a signal")
    }
}

impl std::error::Error for Stopped {}

#[derive(Debug)]
enum Error {
    Usage(String),
    Config(config::Error),
    /// A configuration whose fields contradict each other.
    Invalid(String),
    Certificate(certificate::Error),
    Queue(queue::Error),
    Components(BoxError),
    /// A signal arrived before the command could finish; nothing was changed.
    Stopped(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Usage(message) => write!(f, "{message}\n{USAGE}"),
            Error::Config(e) => e.fmt(f),
            Error::Invalid(message) => write!(f, "invalid configuration: {message}"),
            Error::Certificate(e) => e.fmt(f),
            Error::Queue(e) => e.fmt(f),
            Error::Components(e) => write!(f, "cannot prepare the components: {e}"),
            Error::Stopped(when) => {
                write!(f, "stopped by a signal {when}; the certificate is unchanged")
            }
        }
    }
}

impl From<certificate::Error> for Error {
    fn from(error: certificate::Error) -> Self {
        match &error {
            certificate::Error::Refuted { source, .. } if source.is::<Stopped>() => {
                Error::Stopped("during certificate verification")
            }
            _ => Error::Certificate(error),
        }
    }
}

impl From<queue::Error> for Error {
    fn from(error: queue::Error) -> Self {
        Error::Queue(error)
    }
}

/// Verification of saved records, which a stop request interrupts.
fn verify_saved<C: Components>(
    components: &C,
    stop: &AtomicBool,
    path: &str,
    data: &C::Data,
) -> Result<(), BoxError> {
    if stop.load(Ordering::Acquire) {
        return Err(Box::new(Stopped));
    }
    components.holds(path, data)
}

fn build<C>(
    threads: NonZeroUsize,
    stop: &AtomicBool,
    build: impl FnOnce(NonZeroUsize, &AtomicBool) -> Result<C, BoxError>,
) -> Result<C, Error> {
    let components = build(threads, stop).map_err(Error::Components)?;
    if stop.load(Ordering::Acquire) {
        return Err(Error::Stopped("during preparation"));
    }
    Ok(components)
}

/// Connects the queue to the store. Besides the queue's flush at every depth
/// barrier, an append at least `FLUSH_INTERVAL` after the previous flush
/// flushes, so that a crash inside a long layer loses little work. (A record
/// followed only by refusals waits for the next append or the barrier.)
struct Durable<'a, D> {
    store: &'a mut Store<D>,
    flushed: Instant,
}

impl<D: RecordData> queue::Output<D> for Durable<'_, D> {
    fn append(&mut self, path: &str, data: &D) -> Result<(), BoxError> {
        self.store.append(path, data)?;
        if self.flushed.elapsed() >= FLUSH_INTERVAL {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), BoxError> {
        self.store.flush()?;
        self.flushed = Instant::now();
        Ok(())
    }
}

/// Runs `work` while printing the evaluation count every `PROGRESS_INTERVAL`.
fn with_progress<T>(evaluated: &AtomicU64, work: impl FnOnce() -> T) -> T {
    let started = Instant::now();
    let (finished, waiting) = mpsc::channel::<()>();
    thread::scope(|scope| {
        scope.spawn(move || {
            while let Err(mpsc::RecvTimeoutError::Timeout) =
                waiting.recv_timeout(PROGRESS_INTERVAL)
            {
                let count = evaluated.load(Ordering::Relaxed);
                let seconds = started.elapsed().as_secs();
                eprintln!("rid: evaluated {count} boxes in {seconds} s");
            }
        });
        let result = work();
        drop(finished);
        result
    })
}

fn prepare<C: Components>(
    config: &config::Prepare,
    stop: &AtomicBool,
    make: impl FnOnce(NonZeroUsize, &AtomicBool) -> Result<C, BoxError>,
) -> Result<Value, Error> {
    let components = build(config.threads, stop, make)?;
    Ok(json!({
        "command": "prepare",
        "policy": POLICY,
        "components": components.report(),
    }))
}

fn search<C: Components>(
    config: &config::Search,
    stop: &AtomicBool,
    make: impl FnOnce(NonZeroUsize, &AtomicBool) -> Result<C, BoxError>,
) -> Result<Value, Error> {
    let header = Header::new(POLICY, &config.root, config.max_depth)?;
    let depth_limit = config.depth_limit.unwrap_or(config.max_depth);
    if depth_limit > config.max_depth {
        return Err(Error::Invalid(format!(
            "depth_limit {depth_limit} exceeds the certificate's max_depth {}",
            config.max_depth
        )));
    }
    // Such a search would evaluate nothing and report a depth-limit halt.
    if depth_limit < config.root.len() {
        return Err(Error::Invalid(format!(
            "depth_limit {depth_limit} is shallower than the root {:?} (depth {})",
            config.root,
            config.root.len()
        )));
    }
    let components = build(config.threads, stop, make)?;
    let started = Instant::now();
    let verifier = |path: &str, data: &C::Data| verify_saved(&components, stop, path, data);
    let (mut store, found) = Store::open(&config.certificate, &header, config.threads, verifier)?;
    eprintln!(
        "rid: {} verified records, a frontier of {} boxes, {} unterminated bytes removed",
        found.records, found.frontier, found.unterminated_bytes
    );
    let options = queue::Options {
        threads: config.threads,
        window: config.window,
        depth_limit,
        max_decisions: config.max_decisions,
        after: found.last.clone(),
    };
    let frontier = store.frontier();
    let evaluated = AtomicU64::new(0);
    let result = with_progress(&evaluated, || {
        let mut output = Durable {
            store: &mut store,
            flushed: Instant::now(),
        };
        let evaluate = |path: &str| {
            let decision = components.check(path);
            evaluated.fetch_add(1, Ordering::Relaxed);
            decision
        };
        queue::run(frontier, &options, stop, evaluate, &mut output)
    });
    // Keep what was committed even when the run failed.
    let flushed = store.flush();
    let report = result?;
    flushed?;
    let saved = store.report();
    Ok(json!({
        "command": "search",
        "policy": POLICY,
        "root": header.root(),
        "max_depth": header.max_depth(),
        "depth_limit": depth_limit,
        "halt": match report.halt {
            queue::Halt::Complete => "complete",
            queue::Halt::Stopped => "stopped",
            queue::Halt::DecisionLimit => "decision-limit",
            queue::Halt::DepthLimit => "depth-limit",
        },
        "records": saved.records,
        "new_records": report.records,
        "decisions": report.decisions,
        "evaluated": report.evaluated,
        "splits": report.splits,
        "reconstructed": report.reconstructed,
        "unresolved": report.unresolved,
        "removed_unterminated_bytes": found.unterminated_bytes,
        "certificate_bytes": saved.bytes,
        "milliseconds": u64::try_from(started.elapsed().as_millis()).unwrap_or(u64::MAX),
    }))
}

fn check<C: Components>(
    config: &config::Check,
    stop: &AtomicBool,
    make: impl FnOnce(NonZeroUsize, &AtomicBool) -> Result<C, BoxError>,
) -> Result<Value, Error> {
    let header = Header::new(POLICY, &config.root, config.max_depth)?;
    let components = build(config.threads, stop, make)?;
    let found = certificate::check(&config.certificate, &header, config.threads, |path, data| {
        verify_saved(&components, stop, path, data)
    })?;
    Ok(json!({
        "command": "check",
        "policy": POLICY,
        "root": header.root(),
        "max_depth": header.max_depth(),
        "complete": found.complete(),
        "header_complete": found.header_complete,
        "records": found.records,
        "frontier": found.frontier,
        "unterminated_bytes": found.unterminated_bytes,
    }))
}

fn command(args: &[OsString]) -> Result<Command, Error> {
    let [name, file] = args else {
        return Err(Error::Usage("expected a command and a configuration file".into()));
    };
    let name = name
        .to_str()
        .ok_or_else(|| Error::Usage("the command is not valid text".into()))?;
    Command::load(name, Path::new(file)).map_err(|error| match error {
        config::Error::UnknownCommand(_) => Error::Usage(error.to_string()),
        error => Error::Config(error),
    })
}

/// Runs one invocation `<command> <config.json>` (`args` without the program
/// name): prints the summary as one JSON line on `out` and returns the exit
/// status ([`SUCCESS`], [`FAILURE`] or [`USAGE_ERROR`]). `make` builds the
/// components from the configured thread count; `stop` requests a
/// cooperative stop.
pub fn entry<C: Components>(
    args: &[OsString],
    stop: &AtomicBool,
    make: impl FnOnce(NonZeroUsize, &AtomicBool) -> Result<C, BoxError>,
    out: &mut dyn Write,
) -> u8 {
    let result = command(args).and_then(|command| match &command {
        Command::Prepare(config) => prepare(config, stop, make),
        Command::Search(config) => search(config, stop, make),
        Command::Check(config) => check(config, stop, make),
    });
    match result {
        Ok(summary) => match writeln!(out, "{summary}").and_then(|()| out.flush()) {
            Ok(()) => SUCCESS,
            Err(e) => {
                eprintln!("rid: cannot write the summary: {e}");
                FAILURE
            }
        },
        Err(error @ Error::Usage(_)) => {
            eprintln!("rid: {error}");
            USAGE_ERROR
        }
        Err(error) => {
            eprintln!("rid: {error}");
            FAILURE
        }
    }
}

/// The first SIGINT or SIGTERM sets `stop`, requesting a cooperative stop; a
/// second one ends the process at once with status 1 (recovery from an
/// abrupt end is the same as after SIGKILL).
pub fn install_signal_handlers(stop: &Arc<AtomicBool>) -> io::Result<()> {
    for signal in [signal_hook::consts::SIGINT, signal_hook::consts::SIGTERM] {
        // Registered first, so it sees the flag as the previous signal left it.
        let status = i32::from(FAILURE);
        signal_hook::flag::register_conditional_shutdown(signal, status, Arc::clone(stop))?;
        signal_hook::flag::register(signal, Arc::clone(stop))?;
    }
    Ok(())
}

pub mod wiring;

#[cfg(test)]
mod tests;
