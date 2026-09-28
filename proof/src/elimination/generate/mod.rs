//! The cover generator: for pinned zoom covers, finds cell trees and
//! witnesses for every zoom of every cover, in parallel over all open cells
//! of all zooms, with a result independent of the thread count. Every leaf is
//! accepted by `lemma::check_cell`; the generated file is then verified as
//! a whole by the cover loader before it is written.
use crate::arithmetic::polynomial::PolynomialError;
use crate::elimination::proof::format::{CoverFile, Format, Leaf, Scope, ZoomEntry};
use crate::elimination::zoom::cover::{CoverError, ZoomCover};
use crate::elimination::zoom::ZoomCell;
use crate::elimination::proof::Pin;
use crate::elimination::witness::{Gap, Witness, WitnessError};
use crate::problem::geometry::Point;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

pub mod candidates;
pub mod command;
pub mod search;

use candidates::Candidate;
use search::{ZoomSearch, Decision, Settings};

/// One cover to generate: its pin (name, scope, parameters), the
/// sample grid and the domain inequalities offered as witnesses.
#[derive(Clone, Debug)]
pub struct Request {
    pub pin: Pin,
    /// Each variable of a zoom's base points takes `samples + 1` values; at
    /// most [`MAX_SAMPLES`].
    pub samples: NonZeroUsize,
    pub inequalities: Vec<usize>,
}

/// The largest sample grid accepted: with at most four sampled variables a
/// zoom then has at most 65⁴ ≈ 1.8·10⁷ base points.
pub const MAX_SAMPLES: usize = 64;

/// Counts of a generated cover.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Statistics {
    pub zooms: usize,
    /// Witnesses offered to the zooms from the samples (summed over zooms).
    pub candidates: usize,
    pub nodes: usize,
    pub witnesses: usize,
    pub witness_leaves: usize,
    pub domain_leaves: usize,
    pub delegated_leaves: usize,
    pub depth: usize,
    /// Float-negative candidates refused by the lemma.
    pub refused: usize,
    /// Witness leaves by the origin of their witness: sampled, inherited
    /// from an ancestor's proposals, proposed by the cell itself.
    pub sampled_leaves: usize,
    pub inherited_leaves: usize,
    pub proposed_leaves: usize,
}

/// Reported before each level of the search: the depth of its cells, how
/// many there are, and how many cells were decided before.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Progress {
    pub depth: usize,
    pub cells: usize,
    pub decided: usize,
}

/// A generated cover file and its counts; `seconds` is the processor time
/// spent on its cells (summed over threads).
#[derive(Clone, Debug)]
pub struct Generated {
    pub file: CoverFile,
    pub statistics: Statistics,
    pub seconds: f64,
}

/// Why a cover was not generated.
#[derive(Debug)]
pub enum GenerateError {
    Cover { cover: String, error: CoverError },
    /// A domain inequality index that does not exist.
    Inequality { cover: String, error: WitnessError },
    /// Domain inequalities offered to a cover whose statement holds without D.
    Scope { cover: String },
    Polynomial { cover: String, error: PolynomialError },
    /// A zoom needed a deeper tree than `max_depth` at this cell.
    Depth { cover: String, zoom: String, cell: ZoomCell },
    /// A zoom needed more than `max_nodes` nodes.
    Nodes { cover: String, zoom: String },
    /// A sample grid above [`MAX_SAMPLES`].
    Samples { cover: String, samples: usize },
    /// A `max_depth` above [`search::MAX_DEPTH`] (for all covers).
    MaxDepth { max_depth: usize },
    Interrupted,
}

impl fmt::Display for GenerateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GenerateError::Cover { cover, error } => write!(f, "cover {cover}: {error}"),
            GenerateError::Inequality { cover, error } => write!(f, "cover {cover}: {error}"),
            GenerateError::Scope { cover } => {
                write!(f, "cover {cover}: domain inequalities offered to a cover of scope all")
            }
            GenerateError::Polynomial { cover, error } => write!(f, "cover {cover}: {error}"),
            GenerateError::Depth { cover, zoom, cell } => {
                write!(f, "cover {cover}, zoom {zoom}: no witness found at the depth limit for the cell {cell:?}")
            }
            GenerateError::Nodes { cover, zoom } => write!(f, "cover {cover}, zoom {zoom}: node limit reached"),
            GenerateError::Samples { cover, samples } => {
                write!(f, "cover {cover}: {samples} samples, above the limit {MAX_SAMPLES}")
            }
            GenerateError::MaxDepth { max_depth } => {
                write!(f, "max_depth {max_depth} is above the limit {}", search::MAX_DEPTH)
            }
            GenerateError::Interrupted => write!(f, "generation interrupted"),
        }
    }
}

impl std::error::Error for GenerateError {}

/// `f` on every item on `threads` threads, results in item order. A set
/// `stop` ends the work with `Interrupted`.
fn parallel<T: Sync, R: Send>(
    items: &[T],
    threads: NonZeroUsize,
    stop: &AtomicBool,
    f: impl Fn(&T) -> R + Sync,
) -> Result<Vec<R>, GenerateError> {
    let next = AtomicUsize::new(0);
    let done: Mutex<Vec<(usize, R)>> = Mutex::new(Vec::with_capacity(items.len()));
    std::thread::scope(|scope| {
        for _ in 0..threads.get().min(items.len()) {
            scope.spawn(|| {
                let mut local = Vec::new();
                loop {
                    let index = next.fetch_add(1, Ordering::SeqCst);
                    if index >= items.len() || stop.load(Ordering::Relaxed) {
                        break;
                    }
                    local.push((index, f(&items[index])));
                }
                done.lock().expect("no panics while holding the lock").extend(local);
            });
        }
    });
    if stop.load(Ordering::Relaxed) {
        return Err(GenerateError::Interrupted);
    }
    let mut done = done.into_inner().expect("no panics while holding the lock");
    done.sort_by_key(|(index, _)| *index);
    Ok(done.into_iter().map(|(_, r)| r).collect())
}

/// A node of a zoom's cell tree.
struct Node {
    cell: ZoomCell,
    depth: usize,
    /// Candidates proposed by ancestors (dropped once decided).
    inherited: Vec<Arc<Candidate>>,
    outcome: Option<Outcome>,
}

enum Outcome {
    Leaf(Arc<Candidate>),
    Delegated,
    Split { axis: usize, lower: usize, upper: usize },
}

/// One zoom being searched.
struct Search<'a> {
    cover: usize,
    prepared: ZoomSearch<'a>,
    nodes: Vec<Node>,
    failure: Option<GenerateError>,
    seconds: f64,
}

/// The checks of a request that need no computation: the sample grid is
/// bounded, and domain inequalities exist and are offered only to a cover
/// whose statement holds on D.
pub fn check_request(request: &Request) -> Result<(), GenerateError> {
    let cover = || request.pin.name.clone();
    if request.samples.get() > MAX_SAMPLES {
        return Err(GenerateError::Samples { cover: cover(), samples: request.samples.get() });
    }
    if request.pin.scope == Scope::All && !request.inequalities.is_empty() {
        return Err(GenerateError::Scope { cover: cover() });
    }
    for &i in &request.inequalities {
        candidates::inequality(i).map_err(|error| GenerateError::Inequality { cover: cover(), error })?;
    }
    Ok(())
}

/// The checks of the settings: `max_depth` at most [`search::MAX_DEPTH`].
pub fn check_settings(settings: &Settings) -> Result<(), GenerateError> {
    if settings.max_depth > search::MAX_DEPTH {
        return Err(GenerateError::MaxDepth { max_depth: settings.max_depth });
    }
    Ok(())
}

/// The cover of a request and the domain inequalities it offers.
fn prepare(request: &Request) -> Result<(ZoomCover, Vec<Witness>), GenerateError> {
    let cover = || request.pin.name.clone();
    check_request(request)?;
    let inequalities = request
        .inequalities
        .iter()
        .map(|&i| candidates::inequality(i))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| GenerateError::Inequality { cover: cover(), error })?;
    let zoom_cover =
        ZoomCover::new(request.pin.parameters.clone()).map_err(|error| GenerateError::Cover { cover: cover(), error })?;
    Ok((zoom_cover, inequalities))
}

/// Generates the requested covers on `threads` threads. The result for
/// each request is its cover file or the reason it failed; the files do not
/// depend on the thread count. Every leaf was accepted by the lemma; the
/// file as a whole is verified by the caller before it is written (see
/// [`command`]). The outer error is an interruption or a `max_depth` above
/// [`search::MAX_DEPTH`]; once a zoom of a cover fails, no zoom of that
/// cover is expanded further.
pub fn generate(
    requests: &[Request],
    settings: &Settings,
    threads: NonZeroUsize,
    stop: &AtomicBool,
    progress: &mut dyn FnMut(&Progress),
) -> Result<Vec<Result<Generated, GenerateError>>, GenerateError> {
    check_settings(settings)?;
    let mut errors: Vec<Option<GenerateError>> = Vec::new();
    let mut covers: Vec<Option<ZoomCover>> = Vec::new();
    let mut offered: Vec<Vec<Witness>> = Vec::new();
    for request in requests {
        let (c, w, e) = match prepare(request) {
            Ok((c, w)) => (Some(c), w, None),
            Err(e) => (None, Vec::new(), Some(e)),
        };
        covers.push(c);
        offered.push(w);
        errors.push(e);
    }
    // Every zoom, in cover and zoom order.
    let zooms: Vec<(usize, usize)> = covers
        .iter()
        .enumerate()
        .flat_map(|(r, cover)| (0..cover.as_ref().map_or(0, |c| c.zooms().len())).map(move |z| (r, z)))
        .collect();
    let cover = |r: usize| covers[r].as_ref().expect("zooms only of built covers");
    // Samples at the zooms' base points, and the gaps active there, per cover.
    let sampled = parallel(&zooms, threads, stop, |&(r, z)| candidates::samples(cover(r), z, requests[r].samples))?;
    let mut samples: BTreeMap<usize, BTreeSet<(Point, Point)>> = BTreeMap::new();
    for (&(r, _), points) in zooms.iter().zip(sampled) {
        samples.entry(r).or_default().extend(points);
    }
    let points: Vec<(usize, (Point, Point))> =
        samples.into_iter().flat_map(|(key, set)| set.into_iter().map(move |p| (key, p))).collect();
    let active = parallel(&points, threads, stop, |(_, (u, r))| candidates::active(u, r))?;
    let mut gaps: BTreeMap<usize, BTreeMap<String, Gap>> = BTreeMap::new();
    for ((key, _), list) in points.iter().zip(active) {
        let entry = gaps.entry(*key).or_default();
        for gap in list {
            entry.insert(serde_json::to_string(&gap).expect("gaps serialise"), gap);
        }
    }
    let gaps: BTreeMap<usize, Vec<Gap>> = gaps.into_iter().map(|(k, m)| (k, m.into_values().collect())).collect();
    // Every zoom prepared with its candidates.
    let prepared = parallel(&zooms, threads, stop, |&(r, z)| {
        let list = gaps.get(&r).map_or(&[][..], |v| &v[..]);
        ZoomSearch::new(cover(r), z, list, &offered[r])
    })?;
    let mut searches: Vec<Search> = Vec::new();
    for (&(r, _), p) in zooms.iter().zip(prepared) {
        let cover = || requests[r].pin.name.clone();
        match p {
            Ok(prepared) => {
                let root = prepared.zoom().root().clone();
                let node = Node { cell: root, depth: 0, inherited: Vec::new(), outcome: None };
                searches.push(Search { cover: r, prepared, nodes: vec![node], failure: None, seconds: 0.0 });
            }
            Err(error) => {
                errors[r].get_or_insert(GenerateError::Polynomial { cover: cover(), error });
            }
        }
    }
    // Level by level: every open node of every zoom is decided in parallel.
    let mut failed = vec![false; requests.len()];
    let mut frontier: Vec<(usize, usize)> =
        (0..searches.len()).filter(|&s| errors[searches[s].cover].is_none()).map(|s| (s, 0)).collect();
    let mut level = Progress { depth: 0, cells: 0, decided: 0 };
    while !frontier.is_empty() {
        level.cells = frontier.len();
        progress(&level);
        let decided = parallel(&frontier, threads, stop, |&(s, n)| {
            let (search, start) = (&searches[s], Instant::now());
            let node = &search.nodes[n];
            let decision = search.prepared.decide(&node.cell, node.depth, &node.inherited, settings);
            (decision, start.elapsed().as_secs_f64())
        })?;
        let mut next = Vec::new();
        for (&(s, n), (decision, seconds)) in frontier.iter().zip(decided) {
            let search = &mut searches[s];
            search.seconds += seconds;
            search.nodes[n].inherited = Vec::new();
            if failed[search.cover] {
                continue;
            }
            let cover = requests[search.cover].pin.name.clone();
            let zoom = search.prepared.zoom().name().to_string();
            let outcome = match decision {
                Err(error) => Err(GenerateError::Polynomial { cover, error }),
                Ok(Decision::Unresolved) => Err(GenerateError::Depth { cover, zoom, cell: search.nodes[n].cell.clone() }),
                Ok(Decision::Leaf(candidate)) => Ok(Outcome::Leaf(candidate)),
                Ok(Decision::Delegated) => Ok(Outcome::Delegated),
                Ok(Decision::Split { .. }) if search.nodes.len() + 2 > settings.max_nodes => {
                    Err(GenerateError::Nodes { cover, zoom })
                }
                Ok(Decision::Split { axis, proposed }) => {
                    let depth = search.nodes[n].depth + 1;
                    let halves = search.nodes[n].cell.bisect(axis);
                    let mut ids = [0; 2];
                    for (k, cell) in halves.into_iter().enumerate() {
                        ids[k] = search.nodes.len();
                        search.nodes.push(Node { cell, depth, inherited: proposed.clone(), outcome: None });
                        next.push((s, ids[k]));
                    }
                    Ok(Outcome::Split { axis, lower: ids[0], upper: ids[1] })
                }
            };
            match outcome {
                Ok(outcome) => search.nodes[n].outcome = Some(outcome),
                Err(error) => {
                    // The first failure of the cover in search order.
                    failed[search.cover] = true;
                    search.failure = Some(error);
                }
            }
        }
        // Nodes of covers that failed in this level are not expanded.
        next.retain(|&(s, _)| !failed[searches[s].cover]);
        level.decided += frontier.len();
        level.depth += 1;
        frontier = next;
    }
    // Assemble every cover whose zooms all succeeded.
    let mut results: Vec<Result<Generated, GenerateError>> = Vec::new();
    for (r, request) in requests.iter().enumerate() {
        let own: Vec<&mut Search> = searches.iter_mut().filter(|s| s.cover == r).collect();
        let failure = own.into_iter().find_map(|s| s.failure.take());
        results.push(match errors[r].take().or(failure) {
            Some(error) => Err(error),
            None => Ok(assemble(request, searches.iter().filter(|s| s.cover == r))),
        });
    }
    Ok(results)
}

/// The cover file of one request from its searched zooms (in zoom order,
/// all complete). Witnesses are numbered in the order of their
/// first leaf.
fn assemble<'s, 'a: 's>(request: &Request, searches: impl Iterator<Item = &'s Search<'a>>) -> Generated {
    let pin = &request.pin;
    let mut zooms: Vec<ZoomEntry> = Vec::new();
    let mut index: BTreeMap<String, usize> = BTreeMap::new();
    let mut witnesses: Vec<Witness> = Vec::new();
    let mut statistics = Statistics::default();
    let mut seconds = 0.0;
    for search in searches {
        statistics.zooms += 1;
        statistics.candidates += search.prepared.candidate_count();
        statistics.nodes += search.nodes.len();
        let counts = search.prepared.counts();
        let read = |c: &std::sync::atomic::AtomicUsize| c.load(Ordering::Relaxed);
        statistics.refused += read(&counts.refused);
        statistics.sampled_leaves += read(&counts.sampled);
        statistics.inherited_leaves += read(&counts.inherited);
        statistics.proposed_leaves += read(&counts.proposed);
        seconds += search.seconds;
        let (mut tree, mut leaves) = (String::new(), Vec::new());
        let mut pending = vec![0];
        while let Some(n) = pending.pop() {
            let node = &search.nodes[n];
            statistics.depth = statistics.depth.max(node.depth);
            match node.outcome.as_ref().expect("a complete tree") {
                Outcome::Split { axis, lower, upper } => {
                    tree.push(char::from(b'0' + *axis as u8));
                    pending.push(*upper);
                    pending.push(*lower);
                }
                Outcome::Delegated => {
                    tree.push('.');
                    leaves.push(Leaf::Delegated);
                    statistics.delegated_leaves += 1;
                }
                Outcome::Leaf(candidate) => {
                    tree.push('.');
                    let witness = candidate.witness();
                    let key = serde_json::to_string(witness).expect("witnesses serialise");
                    let next = witnesses.len();
                    let i = *index.entry(key).or_insert(next);
                    if i == next {
                        witnesses.push(witness.clone());
                    }
                    leaves.push(Leaf::Witness { index: i, factor: *candidate.factor().exponents() });
                    statistics.witness_leaves += 1;
                    statistics.domain_leaves += usize::from(candidate.is_domain());
                }
            }
        }
        let name = search.prepared.zoom().name().to_string();
        zooms.push(ZoomEntry { name, tree, leaves });
    }
    statistics.witnesses = witnesses.len();
    let file = CoverFile {
        format: Format::Version1,
        name: pin.name.clone(),
        scope: pin.scope,
        parameters: pin.parameters.clone(),
        witnesses,
        zooms,
    };
    Generated { file, statistics, seconds }
}

#[cfg(test)]
mod tests;
