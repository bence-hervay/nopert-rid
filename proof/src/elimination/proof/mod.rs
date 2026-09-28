//! Proved zoom covers: a cover file is parsed, pinned against the
//! parameters in code, every cell is verified by the lemma (in parallel),
//! every tree is checked to tile its zoom's root, and the result answers the
//! exact per-box inclusion question.
use crate::elimination::zoom::cover::{CoverError, Kind, Parameters, ZoomCover};
use crate::elimination::zoom::tree::{self, TilingError, TreeError};
use crate::elimination::zoom::{ZoomCell, ZoomError, ZoomFactor};
use crate::elimination::lemma::{CellRefusal, Quotient};
use crate::elimination::witness::Witness;
use crate::problem::configuration::ConfigurationBox;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;

pub mod catalogue;
pub mod format;

use format::{CoverFile, FormatError, Leaf, Scope};

/// A set of configurations proved once to contain no fit (in D, for scope
/// `Domain`), answering whether a box lies inside it.
pub trait ProvedSet: Send + Sync {
    fn name(&self) -> &str;
    /// Exactly whether every configuration of `b` lies in the cover.
    fn contains(&self, b: &ConfigurationBox) -> bool;
}

/// What the code expects of a cover file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pin {
    pub name: String,
    pub scope: Scope,
    pub parameters: Parameters,
}

/// Counts of a verified cover, for reports.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Report {
    pub zooms: usize,
    pub witnesses: usize,
    pub witness_leaves: usize,
    pub domain_leaves: usize,
    pub delegated_leaves: usize,
}

/// Where in a cover file something failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Place {
    pub zoom: String,
    pub leaf: usize,
}

#[derive(Debug)]
pub enum ProofError {
    Format(FormatError),
    /// The file differs from the pin in this field.
    Pin(&'static str),
    Cover(CoverError),
    /// The file's zooms are not the cover's zooms in order.
    Zooms,
    Tree { zoom: String, error: TreeError },
    Tiling { zoom: String, error: TilingError },
    /// The number of leaf entries differs from the tree's leaves.
    LeafCount { zoom: String },
    WitnessIndex(Place),
    /// A domain inequality in a cover whose statement holds everywhere.
    Scope(Place),
    /// A delegated leaf outside an A zoom of a window face, or whose window
    /// condition fails.
    Delegation(Place),
    /// A zoom factor with a positive exponent on a variable that is not a
    /// scale variable of the zoom.
    Factor { place: Place, error: ZoomError },
    Cell { place: Place, refusal: CellRefusal },
    Interrupted,
}

impl fmt::Display for ProofError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let at = |p: &Place| format!("zoom {}, leaf {}", p.zoom, p.leaf);
        match self {
            ProofError::Format(error) => write!(f, "{error}"),
            ProofError::Pin(field) => write!(f, "the cover file differs from the code in {field}"),
            ProofError::Cover(error) => write!(f, "{error}"),
            ProofError::Zooms => write!(f, "the file's zooms differ from the cover's"),
            ProofError::Tree { zoom, error } => write!(f, "zoom {zoom}: {error}"),
            ProofError::Tiling { zoom, error } => write!(f, "zoom {zoom}: {error}"),
            ProofError::LeafCount { zoom } => write!(f, "zoom {zoom}: leaf entries differ from the tree"),
            ProofError::WitnessIndex(p) => write!(f, "{}: no such witness", at(p)),
            ProofError::Scope(p) => write!(f, "{}: domain inequality outside the cover's scope", at(p)),
            ProofError::Delegation(p) => write!(f, "{}: invalid delegation", at(p)),
            ProofError::Factor { place, error } => write!(f, "{}: {error}", at(place)),
            ProofError::Cell { place, refusal } => write!(f, "{}: {refusal}", at(place)),
            ProofError::Interrupted => write!(f, "verification interrupted"),
        }
    }
}

impl std::error::Error for ProofError {}

/// A zoom cover whose file (cell trees and witnesses) has been verified
/// completely.
#[derive(Debug)]
pub struct ProvedCover {
    name: String,
    scope: Scope,
    cover: ZoomCover,
    report: Report,
}

impl ProvedSet for ProvedCover {
    fn name(&self) -> &str {
        &self.name
    }
    fn contains(&self, b: &ConfigurationBox) -> bool {
        self.cover.contains(b)
    }
}

/// The witness leaves of one zoom that share a witness and a zoom factor.
struct Job {
    zoom: usize,
    witness: usize,
    factor: ZoomFactor,
    cells: Vec<(usize, ZoomCell)>,
}

/// The work of verification, done in parallel: the tiling check of one
/// zoom's leaves, or the lemma on one job.
enum Task {
    Tiling { zoom: usize, cells: Vec<ZoomCell> },
    Lemma(Job),
}

impl ProvedCover {
    pub fn scope(&self) -> Scope {
        self.scope
    }
    pub fn cover(&self) -> &ZoomCover {
        &self.cover
    }
    pub fn report(&self) -> &Report {
        &self.report
    }

    /// Parses, pins and verifies a cover file on `threads` threads. A set
    /// `stop` flag ends the verification with `Interrupted`.
    pub fn load(pin: &Pin, bytes: &[u8], threads: NonZeroUsize, stop: &AtomicBool) -> Result<Self, ProofError> {
        let file = CoverFile::parse(bytes).map_err(ProofError::Format)?;
        if file.name != pin.name {
            return Err(ProofError::Pin("name"));
        }
        if file.scope != pin.scope {
            return Err(ProofError::Pin("scope"));
        }
        if file.parameters != pin.parameters {
            return Err(ProofError::Pin("parameters"));
        }
        let cover = ZoomCover::new(file.parameters.clone()).map_err(ProofError::Cover)?;
        if pin.scope == Scope::All && cover.parameters().beyond.is_some() {
            return Err(ProofError::Pin("beyond"));
        }
        let names: Vec<&str> = cover.zooms().iter().map(|z| z.name()).collect();
        let listed: Vec<&str> = file.zooms.iter().map(|z| z.name.as_str()).collect();
        // A cover without zooms would prove nothing (`ZoomCover::new` already
        // refuses one; this keeps the invariant at the place that relies on it).
        if names.is_empty() || names != listed {
            return Err(ProofError::Zooms);
        }
        let mut report = Report {
            witnesses: file.witnesses.len(),
            ..Report::default()
        };
        let mut jobs: BTreeMap<(usize, usize, ZoomFactor), Vec<(usize, ZoomCell)>> = BTreeMap::new();
        let mut tasks = Vec::new();
        for (z, (zoom, data)) in cover.zooms().iter().zip(&file.zooms).enumerate() {
            let name = || data.name.clone();
            let cells = tree::leaves(zoom.root(), &data.tree).map_err(|error| ProofError::Tree { zoom: name(), error })?;
            if cells.len() != data.leaves.len() {
                return Err(ProofError::LeafCount { zoom: name() });
            }
            tasks.push(Task::Tiling { zoom: z, cells: cells.clone() });
            report.zooms += 1;
            for (leaf, (cell, entry)) in cells.into_iter().zip(&data.leaves).enumerate() {
                let place = || Place { zoom: name(), leaf };
                match entry {
                    Leaf::Witness { index, factor } => {
                        let witness = file.witnesses.get(*index).ok_or_else(|| ProofError::WitnessIndex(place()))?;
                        if matches!(witness, Witness::Domain(_)) {
                            if pin.scope == Scope::All {
                                return Err(ProofError::Scope(place()));
                            }
                            report.domain_leaves += 1;
                        }
                        // The lemma's first hypothesis: only scale variables are divided.
                        let factor = zoom.factor(*factor).map_err(|error| ProofError::Factor { place: place(), error })?;
                        report.witness_leaves += 1;
                        jobs.entry((z, *index, factor)).or_default().push((leaf, cell));
                    }
                    Leaf::Delegated => {
                        // Only unsheared A zooms delegate (so the sheared family,
                        // whose trees are complete, never does), and only when the
                        // window holds on the whole cell.
                        if !matches!(cover.kind(z), Kind::A { .. }) || !cover.window_holds(z, &cell) {
                            return Err(ProofError::Delegation(place()));
                        }
                        report.delegated_leaves += 1;
                    }
                }
            }
        }
        tasks.extend(
            jobs.into_iter()
                .map(|((zoom, witness, factor), cells)| Task::Lemma(Job { zoom, witness, factor, cells })),
        );
        verify(&tasks, &cover, &file, threads, stop)?;
        Ok(Self {
            name: file.name,
            scope: file.scope,
            cover,
            report,
        })
    }
}

/// Runs every task (the tiling checks, then the lemma on every job) on
/// `threads` threads. Tasks are taken in order; after a refusal no new task
/// starts, and the refusal of the first failing task is reported, so the
/// error does not depend on the timing. A worker that sees `stop` takes no
/// new task and ends its job early; that it did is recorded here, so the
/// result never depends on reading the caller's flag again: a verification
/// that skipped any cell is `Interrupted` (unless a refusal was found).
fn verify(
    tasks: &[Task],
    cover: &ZoomCover,
    file: &CoverFile,
    threads: NonZeroUsize,
    stop: &AtomicBool,
) -> Result<(), ProofError> {
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let interrupted = AtomicBool::new(false);
    let first: Mutex<Option<(usize, ProofError)>> = Mutex::new(None);
    // `Ok(false)`: the job was ended early by `stop`.
    let run = |index: usize| -> Result<bool, ProofError> {
        let job = match &tasks[index] {
            Task::Tiling { zoom, cells } => {
                let root = cover.zooms()[*zoom].root();
                return tree::check_tiling(root, cells).map(|()| true).map_err(|error| ProofError::Tiling {
                    zoom: cover.zooms()[*zoom].name().to_string(),
                    error,
                });
            }
            Task::Lemma(job) => job,
        };
        let zoom = &cover.zooms()[job.zoom];
        let place = |leaf: usize| Place { zoom: zoom.name().to_string(), leaf };
        let witness = &file.witnesses[job.witness];
        let quotient = Quotient::new(zoom, witness, &job.factor)
            .map_err(|refusal| ProofError::Cell { place: place(job.cells[0].0), refusal })?;
        let mut complete = true;
        let cells = job.cells.iter().map(|(_, cell)| cell).take_while(|_| {
            complete = !stop.load(Ordering::Relaxed);
            complete
        });
        quotient
            .check_all(cells)
            .map_err(|(position, refusal)| ProofError::Cell { place: place(job.cells[position].0), refusal })?;
        Ok(complete)
    };
    std::thread::scope(|scope| {
        for _ in 0..threads.get() {
            scope.spawn(|| loop {
                if failed.load(Ordering::SeqCst) {
                    return;
                }
                let index = next.fetch_add(1, Ordering::SeqCst);
                if index >= tasks.len() {
                    return;
                }
                // A job left out, or ended early, because of `stop`.
                if stop.load(Ordering::Relaxed) {
                    interrupted.store(true, Ordering::SeqCst);
                    return;
                }
                match run(index) {
                    Ok(true) => {}
                    Ok(false) => interrupted.store(true, Ordering::SeqCst),
                    Err(error) => {
                        failed.store(true, Ordering::SeqCst);
                        let mut first = first.lock().expect("no panics while holding the lock");
                        if first.as_ref().map_or(true, |(i, _)| index < *i) {
                            *first = Some((index, error));
                        }
                    }
                }
            });
        }
    });
    match first.into_inner().expect("no panics while holding the lock") {
        Some((_, error)) => Err(error),
        None if interrupted.into_inner() => Err(ProofError::Interrupted),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests;
