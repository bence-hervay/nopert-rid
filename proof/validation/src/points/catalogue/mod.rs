//! Exact centres of configurations and the point catalogue of the
//! completeness experiment: its format, loading and refusals.
use rid::arithmetic::exact::{parse_rational, ExactError, QSqrt5};
use rid::problem::configuration::{ConfigurationBox, AXES};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fmt;
use std::path::{Path, PathBuf};

/// The format name of a point catalogue.
pub const FORMAT: &str = "rid-validation-points/1";

/// The spelling of an exact coordinate `a + b√5`: `["a", "b"]`, both rationals
/// in canonical form.
pub type Spelling = [[String; 2]; AXES];

/// An exact configuration `(s, t, r₁, r₂, r₃)` in the closed root box `B₀`,
/// every coordinate in Q(√5).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Centre {
    coordinates: [QSqrt5; AXES],
}

impl Centre {
    /// Refuses a configuration outside the closed root box.
    pub fn new(coordinates: [QSqrt5; AXES]) -> Result<Self, CentreError> {
        let root = ConfigurationBox::root();
        for (axis, (c, interval)) in coordinates.iter().zip(root.axes()).enumerate() {
            let lo = QSqrt5::from_rational(interval.lo().clone());
            let hi = QSqrt5::from_rational(interval.hi().clone());
            if *c < lo || *c > hi {
                return Err(CentreError::OutsideRoot { axis });
            }
        }
        Ok(Self { coordinates })
    }

    /// Parses `[["a","b"], …]`, meaning the coordinates `a + b√5`.
    pub fn parse(spelling: &Spelling) -> Result<Self, CentreError> {
        let mut coordinates = Vec::with_capacity(AXES);
        for [a, b] in spelling {
            coordinates.push(QSqrt5::new(
                parse_rational(a).map_err(CentreError::Number)?,
                parse_rational(b).map_err(CentreError::Number)?,
            ));
        }
        let coordinates: [QSqrt5; AXES] = coordinates.try_into().expect("AXES coordinates");
        Self::new(coordinates)
    }

    pub fn coordinates(&self) -> &[QSqrt5; AXES] {
        &self.coordinates
    }

    /// The canonical spelling, inverse to [`Centre::parse`].
    pub fn spelling(&self) -> Spelling {
        self.coordinates.clone().map(|c| {
            [c.rational_part().to_string(), c.sqrt5_part().to_string()]
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CentreError {
    /// A coordinate part is not a rational in canonical spelling.
    Number(ExactError),
    /// The coordinate on `axis` lies outside the closed root box.
    OutsideRoot { axis: usize },
}

impl fmt::Display for CentreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CentreError::Number(e) => write!(f, "invalid coordinate: {e}"),
            CentreError::OutsideRoot { axis } => {
                write!(f, "coordinate {axis} lies outside the root box")
            }
        }
    }
}

impl std::error::Error for CentreError {}

/// One point of the catalogue. `expected` names the components that are
/// expected to eliminate small neighbourhoods of it; it is reported, never
/// enforced.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point {
    pub id: String,
    pub group: String,
    pub centre: Centre,
    pub expected: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PointSpelling {
    id: String,
    group: String,
    centre: Spelling,
    expected: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CatalogueSpelling {
    format: String,
    points: Vec<PointSpelling>,
}

/// A point catalogue: nonempty, with unique nonempty identifiers and groups,
/// every centre exact and in the closed root box.
#[derive(Clone, Debug)]
pub struct Catalogue {
    pub points: Vec<Point>,
    /// SHA-256 of the file's bytes, in lowercase hexadecimal: transcripts
    /// are bound to it.
    pub sha256: String,
}

#[derive(Debug)]
pub enum CatalogueError {
    Read { file: PathBuf, source: std::io::Error },
    Write { file: PathBuf, source: std::io::Error },
    Json(serde_json::Error),
    /// The format field is not the expected one.
    Format { found: String, expected: &'static str },
    Empty,
    EmptyName { index: usize },
    DuplicateId(String),
    Centre { id: String, source: CentreError },
    /// An empty or repeated name among a point's expected components.
    Expected { id: String },
    /// A target that is not a component name with an optional cover.
    Target { id: String },
    /// A selection that is empty, names an unknown point or names one twice.
    Selection(String),
}

impl fmt::Display for CatalogueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CatalogueError::Read { file, source } => {
                write!(f, "cannot read catalogue {}: {source}", file.display())
            }
            CatalogueError::Write { file, source } => {
                write!(f, "cannot write catalogue {}: {source}", file.display())
            }
            CatalogueError::Json(e) => write!(f, "invalid catalogue: {e}"),
            CatalogueError::Format { found, expected } => {
                write!(f, "catalogue format {found:?}, expected {expected:?}")
            }
            CatalogueError::Empty => write!(f, "the catalogue has no entries"),
            CatalogueError::EmptyName { index } => {
                write!(f, "entry {index} has an empty identifier or group")
            }
            CatalogueError::DuplicateId(id) => write!(f, "identifier {id:?} appears twice"),
            CatalogueError::Centre { id, source } => write!(f, "{id}: {source}"),
            CatalogueError::Expected { id } => {
                write!(f, "{id}: empty or repeated expected component")
            }
            CatalogueError::Target { id } => write!(f, "{id}: invalid target"),
            CatalogueError::Selection(message) => write!(f, "invalid selection: {message}"),
        }
    }
}

impl std::error::Error for CatalogueError {}

/// SHA-256 of `bytes` in lowercase hexadecimal.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Reads a whole catalogue file.
pub fn read(file: &Path) -> Result<Vec<u8>, CatalogueError> {
    std::fs::read(file).map_err(|source| CatalogueError::Read {
        file: file.to_owned(),
        source,
    })
}

/// Refuses an empty identifier or group and a repeated identifier.
pub fn check_names<'a>(
    entries: impl Iterator<Item = (&'a str, &'a str)>,
) -> Result<(), CatalogueError> {
    let mut seen = BTreeSet::new();
    let mut count = 0;
    for (index, (id, group)) in entries.enumerate() {
        if id.is_empty() || group.is_empty() {
            return Err(CatalogueError::EmptyName { index });
        }
        if !seen.insert(id) {
            return Err(CatalogueError::DuplicateId(id.to_owned()));
        }
        count += 1;
    }
    if count == 0 {
        return Err(CatalogueError::Empty);
    }
    Ok(())
}

/// The indices of the entries named by `selection`, in catalogue order; all
/// entries for `None`. A selection must be nonempty and name each known
/// identifier at most once.
pub fn select(ids: &[&str], selection: Option<&[String]>) -> Result<Vec<usize>, CatalogueError> {
    let Some(selection) = selection else {
        return Ok((0..ids.len()).collect());
    };
    if selection.is_empty() {
        return Err(CatalogueError::Selection("no identifier selected".into()));
    }
    let mut wanted = BTreeSet::new();
    for id in selection {
        if !ids.contains(&id.as_str()) {
            return Err(CatalogueError::Selection(format!("unknown identifier {id:?}")));
        }
        if !wanted.insert(id.as_str()) {
            return Err(CatalogueError::Selection(format!("{id:?} selected twice")));
        }
    }
    Ok((0..ids.len()).filter(|&i| wanted.contains(ids[i])).collect())
}

impl Catalogue {
    pub fn load(file: &Path) -> Result<Self, CatalogueError> {
        Self::parse(&read(file)?)
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, CatalogueError> {
        let spelling: CatalogueSpelling =
            serde_json::from_slice(bytes).map_err(CatalogueError::Json)?;
        if spelling.format != FORMAT {
            return Err(CatalogueError::Format {
                found: spelling.format,
                expected: FORMAT,
            });
        }
        check_names(spelling.points.iter().map(|p| (p.id.as_str(), p.group.as_str())))?;
        let mut points = Vec::with_capacity(spelling.points.len());
        for p in spelling.points {
            let centre = Centre::parse(&p.centre).map_err(|source| CatalogueError::Centre {
                id: p.id.clone(),
                source,
            })?;
            let distinct: BTreeSet<&String> = p.expected.iter().collect();
            if distinct.len() != p.expected.len() || p.expected.iter().any(String::is_empty) {
                return Err(CatalogueError::Expected { id: p.id });
            }
            points.push(Point {
                id: p.id,
                group: p.group,
                centre,
                expected: p.expected,
            });
        }
        Ok(Self {
            points,
            sha256: sha256_hex(bytes),
        })
    }

    /// The points named by `selection`, in catalogue order; all for `None`.
    pub fn select(&self, selection: Option<&[String]>) -> Result<Vec<&Point>, CatalogueError> {
        let ids: Vec<&str> = self.points.iter().map(|p| p.id.as_str()).collect();
        Ok(select(&ids, selection)?
            .into_iter()
            .map(|i| &self.points[i])
            .collect())
    }
}

#[cfg(test)]
mod tests;
