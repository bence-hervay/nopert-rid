//! Fixtures shared by the end-to-end tests: a stand-in for the mathematical
//! components, temporary folders with configuration files, and in-process
//! runs of the commands.
#![allow(dead_code)]

use rid::search::certificate::RecordData;
use rid::search::command::{entry, Components, SUCCESS};
use rid::search::BoxError;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::thread;
use std::time::{Duration, SystemTime};

// ---------------------------------------------------------------------------
// A stand-in for the mathematical components
// ---------------------------------------------------------------------------

/// Record data shaped like the components' (an internally tagged enum).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "component", deny_unknown_fields)]
pub enum Mark {
    Even { weight: u32 },
    Odd { weight: u32, pair: [u8; 2] },
}

impl RecordData for Mark {}

pub fn hash(path: &str, seed: u64) -> u64 {
    path.bytes().fold(seed ^ 0xcbf29ce484222325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
    })
}

/// Trees of a few hundred boxes, about half of them records.
fn accepts(path: &str, seed: u64) -> bool {
    path.len() >= 10 || (path.len() >= 2 && hash(path, seed) % 4 == 0)
}

pub fn mark(path: &str, seed: u64) -> Mark {
    let h = hash(path, seed.wrapping_add(1));
    if h % 2 == 0 {
        Mark::Even {
            weight: (h >> 8) as u32 % 1000,
        }
    } else {
        Mark::Odd {
            weight: (h >> 8) as u32 % 1000,
            pair: [(h >> 40) as u8, (h >> 48) as u8],
        }
    }
}

/// Deterministic decisions; optional delays, self-raised signals and a
/// deliberate disagreement between `check` and `holds`.
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StandIn {
    pub seed: u64,
    /// Each decision sleeps 1 to 3 times this long.
    pub delay_us: u64,
    /// Each `holds` sleeps this long.
    pub holds_delay_us: u64,
    /// Raise this signal at the n-th decision.
    pub raise_in_check: Option<(u64, i32)>,
    /// Raise this signal at the n-th `holds`.
    pub raise_in_holds: Option<(u64, i32)>,
    /// Instead of the random rule, accept exactly the boxes this deep or deeper.
    pub accept_from_depth: Option<usize>,
    /// `check` returns a record for this path that `holds` refuses.
    pub false_record: Option<String>,
    /// Each decision first writes its path to this file, so that another
    /// process can see that an evaluation is running.
    pub check_marker: Option<PathBuf>,
    #[serde(skip)]
    pub checked: AtomicU64,
    #[serde(skip)]
    pub held: AtomicU64,
}

impl StandIn {
    pub fn new(seed: u64) -> Self {
        Self {
            seed,
            delay_us: 0,
            holds_delay_us: 0,
            raise_in_check: None,
            raise_in_holds: None,
            accept_from_depth: None,
            false_record: None,
            check_marker: None,
            checked: AtomicU64::new(0),
            held: AtomicU64::new(0),
        }
    }

    fn accepts(&self, path: &str) -> bool {
        match self.accept_from_depth {
            Some(depth) => path.len() >= depth,
            None => accepts(path, self.seed),
        }
    }

    /// The same configuration with fresh counters.
    pub fn clone_config(&self) -> Self {
        Self {
            seed: self.seed,
            delay_us: self.delay_us,
            holds_delay_us: self.holds_delay_us,
            raise_in_check: self.raise_in_check,
            raise_in_holds: self.raise_in_holds,
            accept_from_depth: self.accept_from_depth,
            false_record: self.false_record.clone(),
            check_marker: self.check_marker.clone(),
            checked: AtomicU64::new(0),
            held: AtomicU64::new(0),
        }
    }
}

fn raise_at(count: &AtomicU64, raise: Option<(u64, i32)>) -> Result<(), BoxError> {
    let n = count.fetch_add(1, Ordering::SeqCst) + 1;
    if let Some((at, signal)) = raise {
        if n == at {
            signal_hook::low_level::raise(signal)?;
        }
    }
    Ok(())
}

impl Components for StandIn {
    type Data = Mark;
    fn report(&self) -> Value {
        json!({ "stand_in": self.seed })
    }
    fn check(&self, path: &str) -> Result<Option<Mark>, BoxError> {
        if let Some(marker) = &self.check_marker {
            fs::write(marker, path)?;
        }
        if self.delay_us > 0 {
            thread::sleep(Duration::from_micros(self.delay_us * (1 + hash(path, 99) % 3)));
        }
        raise_at(&self.checked, self.raise_in_check)?;
        if self.false_record.as_deref() == Some(path) {
            return Ok(Some(mark(path, self.seed + 1)));
        }
        Ok(self.accepts(path).then(|| mark(path, self.seed)))
    }
    fn holds(&self, path: &str, data: &Mark) -> Result<(), BoxError> {
        if self.holds_delay_us > 0 {
            thread::sleep(Duration::from_micros(self.holds_delay_us));
        }
        raise_at(&self.held, self.raise_in_holds)?;
        if self.accepts(path) && *data == mark(path, self.seed) {
            Ok(())
        } else {
            Err("not a stand-in record".into())
        }
    }
}

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Spawning a process (fork, then exec) lends every open file of this process
/// to the child until its exec, including a certificate's lock: a lock can
/// seem held just after its handle was dropped. Tests that spawn processes
/// hold this guard exclusively; tests that use certificate files, shared.
static PROCESSES: RwLock<()> = RwLock::new(());

pub fn files() -> RwLockReadGuard<'static, ()> {
    PROCESSES.read().unwrap_or_else(|e| e.into_inner())
}

pub fn processes() -> RwLockWriteGuard<'static, ()> {
    PROCESSES.write().unwrap_or_else(|e| e.into_inner())
}

static SERIAL: AtomicU64 = AtomicU64::new(0);

/// A fresh directory under the system temporary directory, removed on drop.
pub struct Folder(pub PathBuf);

impl Folder {
    pub fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "rid-end-to-end-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    pub fn certificate(&self) -> PathBuf {
        self.0.join("run.cert")
    }
    /// Writes a configuration file and returns the command line for it.
    pub fn args(&self, command: &str, config: Value) -> Vec<OsString> {
        let serial = SERIAL.fetch_add(1, Ordering::Relaxed);
        let file = self.0.join(format!("{command}-{serial}.json"));
        fs::write(&file, config.to_string()).unwrap();
        vec![command.into(), file.into()]
    }
    /// A search of the whole root with a certificate limit of `max_depth`
    /// and no depth limit of its own.
    pub fn search(
        &self,
        threads: usize,
        window: usize,
        max_depth: usize,
        max_decisions: Option<u64>,
    ) -> Vec<OsString> {
        self.search_to(threads, window, max_depth, None, max_decisions)
    }
    pub fn search_to(
        &self,
        threads: usize,
        window: usize,
        max_depth: usize,
        depth_limit: Option<usize>,
        max_decisions: Option<u64>,
    ) -> Vec<OsString> {
        self.args(
            "search",
            json!({
                "certificate": self.certificate(),
                "root": "",
                "threads": threads,
                "window": window,
                "max_depth": max_depth,
                "depth_limit": depth_limit,
                "max_decisions": max_decisions,
            }),
        )
    }
    pub fn check(&self, max_depth: usize) -> Vec<OsString> {
        self.args(
            "check",
            json!({
                "certificate": self.certificate(),
                "root": "",
                "max_depth": max_depth,
                "threads": 3,
            }),
        )
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Runs one invocation in this process with the given components and
/// returns its status and summary.
pub fn run_with<C: Components>(
    args: &[OsString],
    stop: &AtomicBool,
    components: C,
) -> (u8, Option<Value>) {
    let mut out = Vec::new();
    let status = entry(args, stop, |_, _| Ok(components), &mut out);
    let summary = (!out.is_empty()).then(|| serde_json::from_slice(&out).unwrap());
    (status, summary)
}

pub fn run(args: &[OsString], stop: &AtomicBool, stand_in: StandIn) -> (u8, Option<Value>) {
    run_with(args, stop, stand_in)
}

pub fn succeed(args: &[OsString], seed: u64) -> Value {
    let (status, summary) = run(args, &AtomicBool::new(false), StandIn::new(seed));
    assert_eq!(status, SUCCESS);
    summary.unwrap()
}

pub fn snapshot(file: &Path) -> (Vec<u8>, SystemTime) {
    (
        fs::read(file).unwrap(),
        fs::metadata(file).unwrap().modified().unwrap(),
    )
}

/// The certificate of an uninterrupted single-threaded search.
pub fn reference(seed: u64, max_depth: usize) -> Vec<u8> {
    reference_with_summary(seed, max_depth).0
}

pub fn reference_with_summary(seed: u64, max_depth: usize) -> (Vec<u8>, Value) {
    let dir = Folder::new();
    let summary = succeed(&dir.search(1, 1, max_depth, None), seed);
    assert_eq!(summary["unresolved"], 0);
    (fs::read(dir.certificate()).unwrap(), summary)
}
