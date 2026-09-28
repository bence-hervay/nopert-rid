//! The decision for one cell of one zoom: a witness leaf accepted by the
//! lemma, a delegated leaf, or a split axis. Floating point ranks and
//! filters the candidates and chooses the axis; only `lemma::check_cell`
//! accepts a leaf and only the cover's exact window test delegates.
use super::candidates::{self, Candidate};
use crate::elimination::float::{self, Dense, Vector};
use crate::arithmetic::exact::{Interval, Q};
use crate::arithmetic::polynomial::{PolynomialError, VARIABLES};
use crate::elimination::zoom::cover::{ZoomCover, Kind, Shape};
use crate::elimination::zoom::{Five, Zoom, ZoomCell};
use crate::elimination::lemma;
use crate::elimination::witness::{Gap, Witness};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::num::NonZeroU32;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

/// How cells are searched. Every field is required in the configuration.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Settings {
    /// A cell at this depth that is not a leaf makes its zoom fail; at
    /// most [`MAX_DEPTH`].
    pub max_depth: usize,
    /// A zoom whose tree would need more nodes fails.
    pub max_nodes: usize,
    /// The aspect cap (positive): only variables whose width relative to the
    /// root is at least `1/aspect` of the relatively widest variable's are
    /// split, decided exactly from the split counts.
    pub aspect: NonZeroU32,
    /// At most this many sampled or inherited candidates get their Bernstein
    /// coefficients computed per cell (domain inequalities, tried first, and
    /// the cell's own proposals, up to twice `proposals`, not counted).
    pub tried: usize,
    /// Cells at this depth or deeper propose candidates of their own.
    pub proposal_depth: usize,
    /// Separating directions, and edges, proposed per cell (each).
    pub proposals: usize,
    /// At most this many proposed candidates are handed down to a cell's
    /// descendants.
    pub inherited: usize,
}

/// A float Bernstein maximum below `−MARGIN · max(1, max |b|)` is worth an
/// exact check.
const MARGIN: f64 = 1e-10;

/// The deepest `max_depth` accepted: below it every variable of a cell is
/// split at most 200 times, far from where the cell's floating-point ends
/// would coincide.
pub const MAX_DEPTH: usize = 200;

/// The decision for one cell.
#[derive(Clone, Debug)]
pub enum Decision {
    /// The lemma accepted the cell with this candidate.
    Leaf(Arc<Candidate>),
    /// The cover's window condition holds on the whole cell.
    Delegated,
    /// Split along `axis`; `proposed` are the candidates the cell proposed
    /// (with those it inherited), for its descendants.
    Split { axis: usize, proposed: Vec<Arc<Candidate>> },
    /// The cell is at the depth limit and no candidate was accepted.
    Unresolved,
}

/// One zoom prepared for the search: its candidates and float maps.
pub struct ZoomSearch<'a> {
    cover: &'a ZoomCover,
    zoom: usize,
    /// Domain inequalities, tried first where negative at the midpoint.
    inequalities: Vec<Arc<Candidate>>,
    /// Gap candidates from the samples, in a fixed order.
    candidates: Vec<Arc<Candidate>>,
    view: [Dense; 3],
    rotation: [Dense; 3],
    /// Whether cells may be delegated through the cover's window.
    delegates: bool,
    /// Proposed candidates already built, by their witness's JSON form: a
    /// cache of a pure function, so it does not affect any decision.
    built: Mutex<HashMap<String, Option<Arc<Candidate>>>>,
    counts: Counts,
}

/// Counters of what the lemma was asked and what it accepted, by the origin
/// of the candidate (diagnostics only: they do not affect any decision).
#[derive(Debug, Default)]
pub struct Counts {
    /// Float-negative candidates that the lemma refused.
    pub refused: AtomicUsize,
    /// Leaves accepted with a domain inequality, a sampled candidate, an
    /// inherited proposal and a proposal of the cell itself.
    pub domain: AtomicUsize,
    pub sampled: AtomicUsize,
    pub inherited: AtomicUsize,
    pub proposed: AtomicUsize,
}

impl<'a> ZoomSearch<'a> {
    /// Prepares zoom `zoom` of `cover` with the given witnesses (gaps
    /// from samples, domain inequalities). Witnesses whose pull-back
    /// vanishes are dropped.
    pub fn new(cover: &'a ZoomCover, zoom: usize, gaps: &[Gap], inequalities: &[Witness]) -> Result<Self, PolynomialError> {
        let c = &cover.zooms()[zoom];
        let build = |w: Witness| Candidate::new(c, w).map(|o| o.map(Arc::new));
        let mut candidates = Vec::new();
        for gap in gaps {
            candidates.extend(build(Witness::Gap(gap.clone()))?);
        }
        let mut domain = Vec::new();
        for w in inequalities {
            domain.extend(build(w.clone())?);
        }
        let delegates = matches!(
            (&cover.parameters().shape, cover.kind(zoom)),
            (Shape::Point { window: Some(window), .. }, Kind::A { base, .. }) if window.faces.contains(&base)
        );
        Ok(Self {
            cover,
            zoom,
            inequalities: domain,
            candidates,
            view: c.view().vector().map(|p| Dense::new(&p)),
            rotation: c.rotation().clone().map(|p| Dense::new(&p)),
            delegates,
            built: Mutex::new(HashMap::new()),
            counts: Counts::default(),
        })
    }

    pub fn zoom(&self) -> &Zoom {
        &self.cover.zooms()[self.zoom]
    }
    pub fn counts(&self) -> &Counts {
        &self.counts
    }
    pub fn candidate_count(&self) -> usize {
        self.candidates.len() + self.inequalities.len()
    }

    fn view_at(&self, y: &[f64; VARIABLES]) -> Vector {
        std::array::from_fn(|i| self.view[i].evaluate(y))
    }

    /// The candidate of a proposed gap on this zoom, built once.
    fn proposal(&self, gap: Gap) -> Result<Option<Arc<Candidate>>, PolynomialError> {
        let witness = Witness::Gap(gap);
        let key = serde_json::to_string(&witness).expect("witnesses serialise");
        if let Some(found) = self.built.lock().expect("no panics while holding the lock").get(&key) {
            return Ok(found.clone());
        }
        let made = Candidate::new(self.zoom(), witness)?.map(Arc::new);
        self.built.lock().expect("no panics while holding the lock").insert(key, made.clone());
        Ok(made)
    }

    /// The decision for the closed `cell` at `depth`, given the candidates
    /// its ancestors proposed. It depends only on these arguments and the
    /// prepared zoom, never on other cells or on timing.
    pub fn decide(
        &self,
        cell: &ZoomCell,
        depth: usize,
        inherited: &[Arc<Candidate>],
        settings: &Settings,
    ) -> Result<Decision, PolynomialError> {
        if self.delegates && self.cover.window_holds(self.zoom, cell) {
            return Ok(Decision::Delegated);
        }
        let bounds = float::cell(cell);
        let middle: [f64; VARIABLES] = bounds.map(|[lo, hi]| (lo + hi) / 2.0);
        let corners: Vec<Vector> = (0..1usize << VARIABLES)
            .map(|mask| self.view_at(&std::array::from_fn(|j| bounds[j][mask >> j & 1])))
            .collect();
        let mut attempt = Attempt { search: self, cell, bounds: &bounds, tried: Vec::new() };
        let leaf = |c: &Arc<Candidate>, counter: &AtomicUsize| {
            counter.fetch_add(1, Ordering::Relaxed);
            Ok(Decision::Leaf(c.clone()))
        };
        for c in &self.inequalities {
            if c.quotient().evaluate(&middle) < 0.0 && attempt.accepts(c) {
                return leaf(c, &self.counts.domain);
            }
        }
        let mut ranked: Vec<(f64, usize, &Arc<Candidate>)> = self
            .candidates
            .iter()
            .chain(inherited)
            .enumerate()
            .map(|(i, c)| (c.quotient().evaluate(&middle), i, c))
            .collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let before = attempt.tried.len();
        for (value, i, c) in ranked {
            let tried = attempt.tried.len() - before;
            if tried >= settings.tried || (value >= 0.0 && tried > 0) {
                break;
            }
            if c.plausibly_valid(&corners) && attempt.accepts(c) {
                let counter = if i < self.candidates.len() { &self.counts.sampled } else { &self.counts.inherited };
                return leaf(c, counter);
            }
        }
        let mut proposed: Vec<Arc<Candidate>> = Vec::new();
        if depth >= settings.proposal_depth {
            let (u, r) = (self.view_at(&middle), std::array::from_fn(|i| self.rotation[i].evaluate(&middle)));
            for gap in candidates::separating(&u, &r, settings.proposals) {
                let Some(c) = self.proposal(gap)? else { continue };
                if proposed.iter().chain(inherited).any(|p| p.witness() == c.witness()) {
                    continue;
                }
                if c.plausibly_valid(&corners) && attempt.accepts(&c) {
                    return leaf(&c, &self.counts.proposed);
                }
                proposed.push(c);
            }
        }
        if depth >= settings.max_depth {
            return Ok(Decision::Unresolved);
        }
        let axis = split_axis(&splits(self.zoom().root(), cell), &attempt.tried, settings.aspect);
        proposed.extend(inherited.iter().cloned());
        proposed.truncate(settings.inherited);
        Ok(Decision::Split { axis, proposed })
    }
}

/// The candidates tried on one cell, with their float Bernstein
/// coefficients and relative margins.
struct Attempt<'s, 'a> {
    search: &'s ZoomSearch<'a>,
    cell: &'s ZoomCell,
    bounds: &'s Five<[f64; 2]>,
    tried: Vec<(f64, Dense)>,
}

impl Attempt<'_, '_> {
    /// Computes the candidate's float Bernstein coefficients on the cell and,
    /// when they are clearly negative, asks the lemma. (Callers skip gaps
    /// whose support is clearly invalid at the corner views.)
    fn accepts(&mut self, c: &Candidate) -> bool {
        let b = c.quotient().bernstein(self.bounds);
        let (max, scale) = (b.max(), b.magnitude().max(1.0));
        self.tried.push((max / scale, b));
        if max >= -MARGIN * scale {
            return false;
        }
        let accepted = lemma::check_cell(self.search.zoom(), self.cell, c.witness(), c.factor()).is_ok();
        if !accepted {
            self.search.counts.refused.fetch_add(1, Ordering::Relaxed);
        }
        accepted
    }
}

/// How often each variable of `cell` was halved from the zoom's `root`: the
/// cell's width is exactly `2^-k` times the root's.
pub fn splits(root: &Five<Interval>, cell: &Five<Interval>) -> [u32; VARIABLES] {
    std::array::from_fn(|j| {
        let ratio: Q = root[j].width() / cell[j].width();
        debug_assert!(ratio.denom() == &1.into(), "a cell of the tree");
        u32::try_from(ratio.numer().bits() - 1).expect("fewer than 2³² splits")
    })
}

/// The split axis: where the Bernstein coefficients of the two best tried
/// candidates vary most (domain inequalities included), among the variables
/// allowed by the aspect cap: those whose width relative to the root is at
/// least `1/aspect` of the relatively widest variable's. Relative widths are
/// `2^-k` for the split counts `k`, so the cap is the exact integer test
/// `2^(k_j − min k) ≤ aspect`. Without variation, the relatively widest
/// variable (fewest splits) is split; ties go to the lower variable. Never
/// driven by a parent's witness.
pub fn split_axis(splits: &[u32; VARIABLES], tried: &[(f64, Dense)], aspect: NonZeroU32) -> usize {
    let fewest = *splits.iter().min().expect("five variables");
    let allowed: Vec<usize> = (0..VARIABLES)
        .filter(|&j| 1u64.checked_shl(splits[j] - fewest).is_some_and(|ratio| ratio <= u64::from(aspect.get())))
        .collect();
    let mut best: Vec<&(f64, Dense)> = tried.iter().collect();
    best.sort_by(|a, b| a.0.total_cmp(&b.0));
    let variation: [f64; VARIABLES] =
        std::array::from_fn(|j| best.iter().take(2).fold(0.0_f64, |m, (_, b)| m.max(b.variation(j))));
    let by_variation = allowed.iter().copied().fold(allowed[0], |b, j| if variation[j] > variation[b] { j } else { b });
    if variation[by_variation] > 0.0 {
        by_variation
    } else {
        allowed.iter().copied().fold(allowed[0], |b, j| if splits[j] < splits[b] { j } else { b })
    }
}

#[cfg(test)]
mod tests;
