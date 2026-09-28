//! Comparison of two checked completeness transcripts of one catalogue:
//! per-group changes and the points that changed, with regressions marked.
//! Transcripts with different finest radii or component lists are not
//! compared; a point of the baseline missing from the candidate is a
//! regression.
use super::transcript::{Checked, PointResult};
use crate::experiments::probe::Label;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Why two transcripts are not compared.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Incompatible {
    /// Different finest exponents `max_k`: the same point would be judged on
    /// different schedules.
    MaxK(Pair<u32>),
    /// Different component lists.
    Components(Pair<Vec<String>>),
}

impl fmt::Display for Incompatible {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Incompatible::MaxK([b, c]) => write!(f, "the transcripts have different max_k ({b} and {c})"),
            Incompatible::Components([b, c]) => {
                write!(f, "the transcripts have different components ({b:?} and {c:?})")
            }
        }
    }
}

impl std::error::Error for Incompatible {}

/// How a point changed from the baseline to the candidate. The radius is
/// `2^-k`, so a smaller best `k` is a larger radius.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Flag {
    NewlyUnresolved,
    NewlyResolved,
    LargerRadius,
    SmallerRadius,
    ReasonChanged,
    NewFinerRefusals,
    FinerRefusalsCleared,
}

/// A pair of values: baseline first, candidate second.
pub type Pair<T> = [T; 2];

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ComponentChange {
    pub component: String,
    pub best_k: Pair<Option<u32>>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PointChange {
    pub id: String,
    pub group: String,
    pub best_k: Pair<Option<u32>>,
    pub reason: Pair<Option<String>>,
    /// The finer refusals of both transcripts (listed when they differ).
    pub finer_refusals: Pair<Vec<u32>>,
    pub flags: Vec<Flag>,
    /// The components, present in both transcripts, whose own best radius changed.
    pub components: Vec<ComponentChange>,
    /// Newly unresolved, a best `k` larger by at least the threshold, or a
    /// finer refusal the baseline did not have.
    pub regression: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct GroupChange {
    pub group: String,
    pub points: usize,
    pub unresolved: Pair<usize>,
    pub finer_refusals: Pair<usize>,
    pub median_k: Pair<Option<u32>>,
    pub max_k: Pair<Option<u32>>,
    pub larger_radius: usize,
    pub smaller_radius: usize,
    pub reason_changed: usize,
    pub regressions: usize,
    pub reasons: Pair<BTreeMap<String, usize>>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Comparison {
    pub policy: Pair<String>,
    pub max_k: Pair<u32>,
    pub components: Pair<Vec<String>>,
    /// Points of both transcripts, compared.
    pub points: usize,
    pub only_baseline: Vec<String>,
    pub only_candidate: Vec<String>,
    /// Groups in catalogue order.
    pub groups: Vec<GroupChange>,
    /// Changed points in catalogue order.
    pub changes: Vec<PointChange>,
    /// Regressed points plus the points missing from the candidate.
    pub regressions: usize,
}

fn reason(point: &PointResult) -> Option<String> {
    point
        .decision
        .record
        .as_ref()
        .map(|r| Label::of(r).expect("checked").to_string())
}

/// The lower median.
fn median(mut values: Vec<u32>) -> Option<u32> {
    values.sort_unstable();
    values.get(values.len().saturating_sub(1) / 2).copied()
}

fn change(
    baseline: &PointResult,
    candidate: &PointResult,
    worse_by: u32,
    names: &Pair<Vec<String>>,
) -> PointChange {
    let best_k = [baseline.decision.best_k, candidate.decision.best_k];
    let reasons = [reason(baseline), reason(candidate)];
    let mut flags = Vec::new();
    let mut regression = false;
    match best_k {
        [Some(_), None] => {
            flags.push(Flag::NewlyUnresolved);
            regression = true;
        }
        [None, Some(_)] => flags.push(Flag::NewlyResolved),
        [Some(b), Some(c)] if c < b => flags.push(Flag::LargerRadius),
        [Some(b), Some(c)] if c > b => {
            flags.push(Flag::SmallerRadius);
            regression |= c - b >= worse_by;
        }
        _ => {}
    }
    if reasons[0] != reasons[1] {
        flags.push(Flag::ReasonChanged);
    }
    let finer = [&baseline.decision.finer_refusals, &candidate.decision.finer_refusals];
    if finer[1].iter().any(|k| !finer[0].contains(k)) {
        flags.push(Flag::NewFinerRefusals);
        regression = true;
    }
    if finer[0].iter().any(|k| !finer[1].contains(k)) {
        flags.push(Flag::FinerRefusalsCleared);
    }
    let mut components = Vec::new();
    for (i, name) in names[0].iter().enumerate() {
        if let Some(j) = names[1].iter().position(|n| n == name) {
            let pair = [baseline.components[i].best_k, candidate.components[j].best_k];
            if pair[0] != pair[1] {
                components.push(ComponentChange {
                    component: name.clone(),
                    best_k: pair,
                });
            }
        }
    }
    PointChange {
        id: baseline.id.clone(),
        group: baseline.group.clone(),
        best_k,
        reason: reasons,
        finer_refusals: [finer[0].clone(), finer[1].clone()],
        flags,
        components,
        regression,
    }
}

/// Compares two transcripts that were checked against the same catalogue,
/// with the same `max_k` and components. Points are matched by identifier;
/// `worse_by` is the least increase of the best `k` that counts as a
/// regression; a baseline point missing from the candidate counts as one.
pub fn compare(baseline: &Checked, candidate: &Checked, worse_by: u32) -> Result<Comparison, Incompatible> {
    let names = [
        baseline.header.components.clone(),
        candidate.header.components.clone(),
    ];
    if baseline.header.max_k != candidate.header.max_k {
        return Err(Incompatible::MaxK([baseline.header.max_k, candidate.header.max_k]));
    }
    if names[0] != names[1] {
        return Err(Incompatible::Components(names));
    }
    let by_id: BTreeMap<&str, &PointResult> =
        candidate.points.iter().map(|p| (p.id.as_str(), p)).collect();
    let in_baseline: BTreeSet<&str> = baseline.points.iter().map(|p| p.id.as_str()).collect();
    let mut groups: Vec<GroupChange> = Vec::new();
    let mut changes = Vec::new();
    let mut ks: BTreeMap<String, Pair<Vec<u32>>> = BTreeMap::new();
    let mut compared = 0;
    for b in &baseline.points {
        let Some(c) = by_id.get(b.id.as_str()) else {
            continue;
        };
        compared += 1;
        let point = change(b, c, worse_by, &names);
        let index = match groups.iter().position(|g| g.group == b.group) {
            Some(index) => index,
            None => {
                groups.push(GroupChange {
                    group: b.group.clone(),
                    points: 0,
                    unresolved: [0, 0],
                    finer_refusals: [0, 0],
                    median_k: [None, None],
                    max_k: [None, None],
                    larger_radius: 0,
                    smaller_radius: 0,
                    reason_changed: 0,
                    regressions: 0,
                    reasons: [BTreeMap::new(), BTreeMap::new()],
                });
                groups.len() - 1
            }
        };
        let group = &mut groups[index];
        group.points += 1;
        let k = ks.entry(b.group.clone()).or_default();
        for (side, p) in [(0, b), (1, *c)] {
            match p.decision.best_k {
                None => group.unresolved[side] += 1,
                Some(best) => k[side].push(best),
            }
            if !p.decision.finer_refusals.is_empty() {
                group.finer_refusals[side] += 1;
            }
            let r = reason(p).unwrap_or_else(|| "unresolved".into());
            *group.reasons[side].entry(r).or_default() += 1;
        }
        let has = |f: Flag| point.flags.contains(&f);
        group.larger_radius += usize::from(has(Flag::LargerRadius) || has(Flag::NewlyResolved));
        group.smaller_radius += usize::from(has(Flag::SmallerRadius) || has(Flag::NewlyUnresolved));
        group.reason_changed += usize::from(has(Flag::ReasonChanged));
        group.regressions += usize::from(point.regression);
        if !point.flags.is_empty() || !point.components.is_empty() {
            changes.push(point);
        }
    }
    for group in &mut groups {
        let k = &ks[&group.group];
        group.median_k = [median(k[0].clone()), median(k[1].clone())];
        group.max_k = [k[0].iter().max().copied(), k[1].iter().max().copied()];
    }
    let only_baseline: Vec<String> = baseline
        .points
        .iter()
        .filter(|p| !by_id.contains_key(p.id.as_str()))
        .map(|p| p.id.clone())
        .collect();
    Ok(Comparison {
        policy: [baseline.header.policy.clone(), candidate.header.policy.clone()],
        max_k: [baseline.header.max_k, candidate.header.max_k],
        components: names,
        points: compared,
        only_candidate: candidate
            .points
            .iter()
            .filter(|p| !in_baseline.contains(p.id.as_str()))
            .map(|p| p.id.clone())
            .collect(),
        regressions: changes.iter().filter(|c| c.regression).count() + only_baseline.len(),
        only_baseline,
        groups,
        changes,
    })
}
