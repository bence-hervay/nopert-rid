//! The target catalogue of the irredundancy experiment: exact centres, each
//! with the component (and cover) it is meant to need; the targets of the
//! Local covers are generated from the crate's catalogue of covers.
use crate::config;
use crate::experiments::probe::Label;
use crate::points::catalogue::{self, CatalogueError, Centre, Spelling};
use crate::points::relation;
use rid::arithmetic::exact::{frac, QSqrt5, Q};
use rid::elimination::proof::catalogue as covers;
use rid::elimination::zoom::cover::ZoomCover;
use rid::problem::configuration::ConfigurationBox;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::path::Path;

pub const FORMAT: &str = "rid-validation-targets/1";

/// A primary target should be eliminated only by its own component or
/// cover; a control demonstrates an expected overlap.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Category {
    Primary,
    Control,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Target {
    pub id: String,
    pub category: Category,
    /// The component, and optionally its cover, the target is meant to need.
    pub target: Label,
    pub centre: Centre,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetSpelling {
    id: String,
    category: Category,
    /// Spelt like a record: `{"component": …, "cover": …}`.
    target: Value,
    centre: Spelling,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetsSpelling {
    format: String,
    targets: Vec<TargetSpelling>,
}

#[derive(Clone, Debug)]
pub struct Targets {
    pub targets: Vec<Target>,
    pub sha256: String,
}

impl Targets {
    pub fn load(file: &Path) -> Result<Self, CatalogueError> {
        Self::parse(&catalogue::read(file)?)
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, CatalogueError> {
        let spelling: TargetsSpelling =
            serde_json::from_slice(bytes).map_err(CatalogueError::Json)?;
        if spelling.format != FORMAT {
            return Err(CatalogueError::Format {
                found: spelling.format,
                expected: FORMAT,
            });
        }
        catalogue::check_names(spelling.targets.iter().map(|t| (t.id.as_str(), "target")))?;
        let mut targets = Vec::with_capacity(spelling.targets.len());
        for t in spelling.targets {
            let target = Label::of(&t.target).map_err(|_| CatalogueError::Target {
                id: t.id.clone(),
            })?;
            let object = t.target.as_object().expect("a label is read from an object");
            if object.keys().any(|k| k != "component" && k != "cover") {
                return Err(CatalogueError::Target { id: t.id });
            }
            let centre = Centre::parse(&t.centre).map_err(|source| CatalogueError::Centre {
                id: t.id.clone(),
                source,
            })?;
            targets.push(Target {
                id: t.id,
                category: t.category,
                target,
                centre,
            });
        }
        Ok(Self {
            targets,
            sha256: catalogue::sha256_hex(bytes),
        })
    }

    /// The targets named by `selection`, in catalogue order; all for `None`.
    pub fn select(&self, selection: Option<&[String]>) -> Result<Vec<&Target>, CatalogueError> {
        let ids: Vec<&str> = self.targets.iter().map(|t| t.id.as_str()).collect();
        Ok(catalogue::select(&ids, selection)?
            .into_iter()
            .map(|i| &self.targets[i])
            .collect())
    }
}

// ---- Targets generated from the crate's covers ------------------------------

/// How often the view rectangle is halved, at most, while looking for a
/// representative configuration.
const LEVELS: u32 = 6;

/// The identifier of the generated target of Local cover `number`.
pub fn local_id(number: usize) -> String {
    format!("local-{number}")
}

/// The covers of the crate's catalogue, Exotic in their order, then
/// Local cover `i` at index `8 + i`.
fn all_covers() -> Vec<ZoomCover> {
    let pins = covers::EXOTIC
        .iter()
        .map(|name| covers::exotic_pin(name).expect("pinned"))
        .chain((0..covers::LOCAL_COUNT).map(|n| covers::local_pin(n).expect("pinned")));
    pins.map(|pin| ZoomCover::new(pin.parameters).expect("a pinned cover builds"))
        .collect()
}

/// The representative configuration of Local cover `number`: the first view
/// `(s, t)`, with the aligned rotation `r = 0`, that lies in D, in the cover
/// and in no other cover of the catalogue, among the centre of the cover's
/// view rectangle, then the centres of its 2 × 2, 4 × 4, … equal parts (by
/// increasing `s`, then `t`), down to `2^LEVELS` parts per side. No other
/// cover can eliminate a neighbourhood of such a configuration (they are
/// closed and miss it), Domain cannot (it lies in D) and Global cannot (it
/// is a touch), while its own cover eliminates every small one.
pub fn local_representative(number: usize) -> Option<Centre> {
    representative(number, &all_covers())
}

fn representative(number: usize, all: &[ZoomCover]) -> Option<Centre> {
    let own = covers::EXOTIC.len() + number;
    let [s, t] = covers::local_rectangle(number)?;
    for level in 0..=LEVELS {
        let parts = 1i64 << level;
        for i in 0..parts {
            for j in 0..parts {
                let at = |[lo, hi]: &[Q; 2], k: i64| lo + &((hi - lo) * frac(2 * k + 1, 2 * parts));
                let x: [Q; 5] = [at(&s, i), at(&t, j), Q::zero(), Q::zero(), Q::zero()];
                let point = ConfigurationBox::point(&x);
                let Ok(centre) = Centre::new(x.clone().map(QSqrt5::from_rational)) else {
                    continue;
                };
                let alone = all
                    .iter()
                    .enumerate()
                    .all(|(k, cover)| cover.contains(&point) == (k == own));
                if alone && relation::in_domain(&centre) {
                    return Some(centre);
                }
            }
        }
    }
    None
}

/// The generated targets, one primary target per Local cover, in cover
/// order: `local-<n>`, the cover `{"component": "Local", "cover": n}` and its
/// representative configuration.
fn generated() -> Result<Vec<TargetSpelling>, CatalogueError> {
    let all = all_covers();
    (0..covers::LOCAL_COUNT)
        .map(|n| {
            let centre = representative(n, &all).ok_or(CatalogueError::Target { id: local_id(n) })?;
            Ok(TargetSpelling {
                id: local_id(n),
                category: Category::Primary,
                target: json!({"component": "Local", "cover": n}),
                centre: centre.spelling(),
            })
        })
        .collect()
}

/// Whether a spelt target names one Local cover (and is therefore generated).
fn names_a_local_cover(target: &Value) -> bool {
    Label::of(target).is_ok_and(|l| l.component == "Local" && l.cover.is_some())
}

/// The target catalogue `bytes` with the targets that name a Local cover
/// replaced by the generated ones, appended in cover order, and spelt
/// canonically: one target per line.
pub fn regenerate(bytes: &[u8]) -> Result<Vec<u8>, CatalogueError> {
    Targets::parse(bytes)?;
    let spelling: TargetsSpelling = serde_json::from_slice(bytes).map_err(CatalogueError::Json)?;
    let mut targets: Vec<TargetSpelling> =
        spelling.targets.into_iter().filter(|t| !names_a_local_cover(&t.target)).collect();
    targets.extend(generated()?);
    let lines: Vec<String> = targets
        .iter()
        .map(|t| serde_json::to_string(t).map(|line| format!("  {line}")))
        .collect::<Result<_, _>>()
        .map_err(CatalogueError::Json)?;
    let body = lines.join(",\n");
    let text = format!("{{\n \"format\": \"{FORMAT}\",\n \"targets\": [\n{body}\n ]\n}}\n");
    Targets::parse(text.as_bytes())?;
    Ok(text.into_bytes())
}

/// The summary of the `targets` command.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Written {
    pub targets: usize,
    pub generated: usize,
    /// Whether the file changed.
    pub changed: bool,
    pub sha256: String,
}

/// `targets`: regenerates the Local covers' targets of the catalogue file in
/// place.
pub fn write(config: &config::Targets) -> Result<Written, CatalogueError> {
    let bytes = catalogue::read(&config.catalogue)?;
    let new = regenerate(&bytes)?;
    let changed = new != bytes;
    if changed {
        std::fs::write(&config.catalogue, &new)
            .map_err(|source| CatalogueError::Write { file: config.catalogue.clone(), source })?;
    }
    let targets = Targets::parse(&new)?;
    let generated = targets.targets.iter().filter(|t| t.target.component == "Local" && t.target.cover.is_some());
    Ok(Written {
        generated: generated.count(),
        targets: targets.targets.len(),
        changed,
        sha256: targets.sha256,
    })
}
