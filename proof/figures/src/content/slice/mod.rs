//! Bounded subdivision of a two-dimensional slice of configuration space,
//! each cell classified by the component that eliminates it; the complete
//! partition as data.
use super::number::{interval, Rational};
use super::Point;
use rid::arithmetic::exact::{ExactError, Interval, Q};
use rid::components::record::{ExoticCover, RecordData};
use rid::problem::configuration::{ConfigurationBox, AXES};
use rid::search::BoxError;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// The deepest subdivision a slice may ask for.
pub const MAX_DEPTH: usize = 40;
/// The first field of every partition file.
pub const FORMAT: &str = "rid-figures-partition/1";

/// A slice as written in figure descriptions and partition files. The plot
/// parameters `(a, b)` range over `plot = [[a₀, a₁], [b₀, b₁]]` and stand
/// for the configurations `origin + a·horizontal + b·vertical + δ` with
/// `|δ_i| ≤ thickness_i` in each configuration coordinate.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub origin: [Rational; AXES],
    pub horizontal: [Rational; AXES],
    pub vertical: [Rational; AXES],
    pub plot: [[Rational; 2]; 2],
    pub thickness: [Rational; AXES],
    /// Cells at this depth are not split further.
    pub max_depth: usize,
    /// The largest number of cells: a subdivision that would need more is
    /// refused before the layer that would exceed it is classified, so at
    /// most about `2 · max_cells` boxes are ever classified.
    pub max_cells: NonZeroUsize,
    /// The parts of the collection the slice is subdivided with.
    pub components: Parts,
}

/// Parts of the collection: some of the components and some of the
/// Exotic covers, tried in the collection's order. With every part,
/// a cell's outcome is the collection's own decision. Written
/// `{"Domain": true, "Global": true, "Local": false, "Exotic": ["square"]}`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parts {
    #[serde(rename = "Domain")]
    pub domain: bool,
    #[serde(rename = "Global")]
    pub global: bool,
    #[serde(rename = "Local")]
    pub local: bool,
    /// In the collection's order; repeats are refused.
    #[serde(rename = "Exotic", deserialize_with = "covers")]
    pub exotic: Vec<ExoticCover>,
}

impl Parts {
    /// The whole collection.
    pub fn all() -> Self {
        Self { domain: true, global: true, local: true, exotic: ExoticCover::ALL.to_vec() }
    }
}

/// Exotic covers in the collection's order, without repeats.
fn covers<'de, D: serde::Deserializer<'de>>(deserializer: D) -> Result<Vec<ExoticCover>, D::Error> {
    let covers = Vec::<ExoticCover>::deserialize(deserializer)?;
    if covers.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(serde::de::Error::custom("Exotic covers must be in the collection's order, without repeats"));
    }
    Ok(covers)
}

/// What eliminated a cell: the component, and the cover for Exotic.
/// Written `{"component": "Global"}` or
/// `{"component": "Exotic", "cover": "crossing+"}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "CategoryFields", into = "CategoryFields")]
pub enum Category {
    Domain,
    Global,
    Local,
    Exotic { cover: ExoticCover },
}

/// The fields of a category: `cover` exactly for Exotic. (serde's
/// tagged unit variants would accept unknown fields.)
#[derive(Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CategoryFields {
    component: Component,
    /// Absent, or a cover: an explicit `null` is refused, so a category has
    /// one spelling.
    #[serde(default, skip_serializing_if = "Option::is_none", deserialize_with = "present")]
    cover: Option<ExoticCover>,
}

/// A present field: its value, never `null`.
fn present<'de, D: serde::Deserializer<'de>>(deserializer: D) -> Result<Option<ExoticCover>, D::Error> {
    ExoticCover::deserialize(deserializer).map(Some)
}

#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
enum Component {
    Domain,
    Global,
    Local,
    Exotic,
}

impl TryFrom<CategoryFields> for Category {
    type Error = String;
    fn try_from(fields: CategoryFields) -> Result<Self, String> {
        match (fields.component, fields.cover) {
            (Component::Domain, None) => Ok(Category::Domain),
            (Component::Global, None) => Ok(Category::Global),
            (Component::Local, None) => Ok(Category::Local),
            (Component::Exotic, Some(cover)) => Ok(Category::Exotic { cover }),
            (Component::Exotic, None) => {
                Err("an Exotic category needs its cover".into())
            }
            (_, Some(_)) => Err("only Exotic categories have a cover".into()),
        }
    }
}

impl From<Category> for CategoryFields {
    fn from(category: Category) -> Self {
        let (component, cover) = match category {
            Category::Domain => (Component::Domain, None),
            Category::Global => (Component::Global, None),
            Category::Local => (Component::Local, None),
            Category::Exotic { cover } => (Component::Exotic, Some(cover)),
        };
        Self { component, cover }
    }
}

/// The final state of a cell of the partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Outcome {
    /// The component eliminated the cell's box.
    Eliminated(Category),
    /// At the depth limit and not eliminated.
    Unresolved,
    /// The cell's box misses the root box `B₀`, so no configuration of the
    /// cell has to be searched.
    Outside,
}

/// A leaf of the subdivision. Bit `d` of the path keeps the lower (`0`) or
/// upper (`1`) closed half of plot axis `d mod 2` (`a` first).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Cell {
    pub path: String,
    pub outcome: Outcome,
}

/// The components as the slice sees them.
pub trait Classifier: Sync {
    /// The category of the first of the `parts` that eliminates every
    /// configuration of the closed box `b ⊂ B₀`, or `None`.
    fn classify(&self, b: &ConfigurationBox, parts: &Parts) -> Result<Option<Category>, BoxError>;
}

impl<F> Classifier for F
where
    F: Fn(&ConfigurationBox, &Parts) -> Result<Option<Category>, BoxError> + Sync,
{
    fn classify(&self, b: &ConfigurationBox, parts: &Parts) -> Result<Option<Category>, BoxError> {
        self(b, parts)
    }
}

#[derive(Debug)]
pub enum SliceError {
    Number(ExactError),
    /// A plot range of zero width (`axis` 0 is `a`).
    EmptyPlot { axis: usize },
    /// The two directions are linearly dependent.
    DependentDirections,
    NegativeThickness { axis: usize },
    DepthLimit { max_depth: usize },
    /// The subdivision needs more than `max_cells` cells.
    TooManyCells { max_cells: usize },
    /// A path with a byte other than `0`/`1`, or deeper than the slice's limit.
    InvalidPath(String),
    Classifier(BoxError),
    Json(serde_json::Error),
    Format(String),
    Policy { expected: String, found: String },
    /// Cells not in strictly increasing order of depth, then path.
    Order { path: String },
    /// No cell covers this node at the depth limit.
    Gap { path: String },
    /// A cell lies inside another cell.
    Overlap { path: String },
    /// A split node whose box misses `B₀`, or an `outside` cell whose box
    /// meets it, or an eliminated cell whose box misses it.
    WrongOutside { path: String },
    /// An unresolved cell above the depth limit.
    EarlyUnresolved { path: String },
}

impl fmt::Display for SliceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SliceError::Number(e) => write!(f, "{e}"),
            SliceError::EmptyPlot { axis } => write!(f, "plot range {axis} has zero width"),
            SliceError::DependentDirections => {
                write!(f, "the horizontal and vertical directions are linearly dependent")
            }
            SliceError::NegativeThickness { axis } => write!(f, "thickness {axis} is negative"),
            SliceError::DepthLimit { max_depth } => {
                write!(f, "max_depth {max_depth} exceeds {MAX_DEPTH}")
            }
            SliceError::TooManyCells { max_cells } => {
                write!(f, "the subdivision needs more than max_cells = {max_cells} cells")
            }
            SliceError::InvalidPath(path) => write!(f, "invalid slice path {path:?}"),
            SliceError::Classifier(e) => write!(f, "classification failed: {e}"),
            SliceError::Json(e) => write!(f, "partition file: {e}"),
            SliceError::Format(found) => {
                write!(f, "partition format {found:?}, expected {FORMAT:?}")
            }
            SliceError::Policy { expected, found } => write!(
                f,
                "partition computed under policy {found}, the crate's policy is {expected}"
            ),
            SliceError::Order { path } => write!(f, "cell {path:?} is out of order"),
            SliceError::Gap { path } => write!(f, "no cell covers {path:?}"),
            SliceError::Overlap { path } => write!(f, "cell {path:?} lies inside another cell"),
            SliceError::WrongOutside { path } => {
                write!(f, "cell {path:?} contradicts its box's intersection with the root box")
            }
            SliceError::EarlyUnresolved { path } => {
                write!(f, "cell {path:?} is unresolved above the depth limit")
            }
        }
    }
}

impl std::error::Error for SliceError {}

impl Category {
    /// The category of a record of the crate's components, by its type (an
    /// exhaustive match: a new component is a compile error here; the
    /// Exotic covers are the crate's own type).
    pub fn of(record: &RecordData) -> Self {
        match record {
            RecordData::Domain { .. } => Category::Domain,
            RecordData::Global(_) => Category::Global,
            RecordData::Local { .. } => Category::Local,
            RecordData::Exotic { cover } => Category::Exotic { cover: *cover },
        }
    }
}

/// A validated slice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Slice {
    definition: Definition,
    plot: [Interval; 2],
}

impl Slice {
    pub fn new(definition: Definition) -> Result<Self, SliceError> {
        let plot = [
            interval(&definition.plot[0]).map_err(SliceError::Number)?,
            interval(&definition.plot[1]).map_err(SliceError::Number)?,
        ];
        if let Some(axis) = plot.iter().position(|p| p.lo() == p.hi()) {
            return Err(SliceError::EmptyPlot { axis });
        }
        let (h, v) = (&definition.horizontal, &definition.vertical);
        let independent = (0..AXES)
            .any(|i| (0..i).any(|j| &h[i].0 * &v[j].0 != &h[j].0 * &v[i].0));
        if !independent {
            return Err(SliceError::DependentDirections);
        }
        if let Some(axis) = definition.thickness.iter().position(|w| w.0 < Q::zero()) {
            return Err(SliceError::NegativeThickness { axis });
        }
        if definition.max_depth > MAX_DEPTH {
            return Err(SliceError::DepthLimit {
                max_depth: definition.max_depth,
            });
        }
        Ok(Self { definition, plot })
    }

    pub fn definition(&self) -> &Definition {
        &self.definition
    }

    /// The plot rectangle `[a] × [b]`: the world rectangle of the slice.
    pub fn plot(&self) -> &[Interval; 2] {
        &self.plot
    }

    /// The closed plot rectangle of a path.
    pub fn rectangle(&self, path: &str) -> Result<[Interval; 2], SliceError> {
        if path.len() > self.definition.max_depth {
            return Err(SliceError::InvalidPath(path.to_owned()));
        }
        let mut rectangle = self.plot.clone();
        for (depth, byte) in path.bytes().enumerate() {
            let half = match byte {
                b'0' => 0,
                b'1' => 1,
                _ => return Err(SliceError::InvalidPath(path.to_owned())),
            };
            let axis = depth % 2;
            rectangle[axis] = rectangle[axis].bisect()[half].clone();
        }
        Ok(rectangle)
    }

    /// The smallest box containing every configuration of the rectangle's
    /// thickened image: coordinate `i` is affine in the independent `a` and
    /// `b`, so its range is `oᵢ + hᵢ·[a] + vᵢ·[b] + [-wᵢ, wᵢ]` exactly.
    pub fn enclosure(&self, rectangle: &[Interval; 2]) -> [Interval; AXES] {
        let d = &self.definition;
        std::array::from_fn(|i| {
            let w = &d.thickness[i].0;
            let thickness = Interval::new(-w, w.clone()).expect("nonnegative thickness");
            let h = &Interval::point(d.horizontal[i].0.clone()) * &rectangle[0];
            let v = &Interval::point(d.vertical[i].0.clone()) * &rectangle[1];
            &(&(&Interval::point(d.origin[i].0.clone()) + &h) + &v) + &thickness
        })
    }

    /// The enclosure intersected with the root box `B₀`, or `None` when they
    /// are disjoint. Every configuration of the cell that lies in the domain
    /// `D ⊂ B₀` lies in this box.
    pub fn clipped(&self, rectangle: &[Interval; 2]) -> Option<ConfigurationBox> {
        let enclosure = self.enclosure(rectangle);
        let root = ConfigurationBox::root();
        let mut axes = Vec::with_capacity(AXES);
        for (e, r) in enclosure.iter().zip(root.axes()) {
            let lo = e.lo().max(r.lo()).clone();
            let hi = e.hi().min(r.hi()).clone();
            axes.push(Interval::new(lo, hi).ok()?);
        }
        Some(ConfigurationBox::new(axes.try_into().expect("five axes")))
    }

    /// The outcome of the cell at `path`, or `None` when it is to be split.
    fn decide(
        &self,
        path: &str,
        classifier: &dyn Classifier,
    ) -> Result<Option<Outcome>, SliceError> {
        let Some(b) = self.clipped(&self.rectangle(path)?) else {
            return Ok(Some(Outcome::Outside));
        };
        match classifier.classify(&b, &self.definition.components).map_err(SliceError::Classifier)? {
            Some(category) => Ok(Some(Outcome::Eliminated(category))),
            None if path.len() == self.definition.max_depth => Ok(Some(Outcome::Unresolved)),
            None => Ok(None),
        }
    }
}

/// A complete subdivision of a slice's plot rectangle into cells.
///
/// **Invariant:** the cells' paths are in strictly increasing order of depth
/// then path, and they are the leaves of a binary tree whose root is the
/// empty path: their closed rectangles cover the plot rectangle and have
/// disjoint interiors. Each outcome is the collection's decision for its box.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Partition {
    slice: Slice,
    cells: Vec<Cell>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct File {
    format: String,
    policy: String,
    slice: Definition,
    cells: Vec<Cell>,
}

impl Partition {
    /// Subdivides the slice depth by depth: a cell whose box misses `B₀` is
    /// `outside`, one the classifier eliminates is `eliminated`, one at the
    /// depth limit is `unresolved`, and every other cell is split at the
    /// midpoint of plot axis `depth mod 2`. Each depth is classified by
    /// `threads` workers; the result does not depend on their number.
    pub fn compute(
        slice: &Slice,
        classifier: &dyn Classifier,
        threads: NonZeroUsize,
    ) -> Result<Self, SliceError> {
        let mut cells = Vec::new();
        let mut layer = vec![String::new()];
        let max_cells = slice.definition.max_cells.get();
        while !layer.is_empty() {
            // Every node of the layer ends in at least one cell.
            if cells.len() + layer.len() > max_cells {
                return Err(SliceError::TooManyCells { max_cells });
            }
            let decisions = parallel(&layer, threads, |path| slice.decide(path, classifier))?;
            let mut next = Vec::new();
            for (path, decision) in layer.into_iter().zip(decisions) {
                match decision {
                    Some(outcome) => cells.push(Cell { path, outcome }),
                    None => {
                        next.push(format!("{path}0"));
                        next.push(format!("{path}1"));
                    }
                }
            }
            layer = next;
        }
        Ok(Self {
            slice: slice.clone(),
            cells,
        })
    }

    pub fn slice(&self) -> &Slice {
        &self.slice
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// The partition file: one line per cell.
    pub fn to_json(&self, policy: &str) -> String {
        let cells: Vec<String> = self.cells.iter().map(|c| format!("  {}", json(c))).collect();
        format!(
            "{{\"format\":{},\"policy\":{},\"slice\":{},\"cells\":[\n{}\n]}}\n",
            json(FORMAT),
            json(policy),
            json(self.slice.definition()),
            cells.join(",\n")
        )
    }

    /// Reads a partition file written under `policy` and checks its
    /// structure: the format and policy, the slice, the order and tiling of
    /// the cells, and each outcome against its recomputed box. The
    /// classifications themselves are not repeated.
    pub fn from_json(bytes: &[u8], policy: &str) -> Result<Self, SliceError> {
        let file: File = serde_json::from_slice(bytes).map_err(SliceError::Json)?;
        if file.format != FORMAT {
            return Err(SliceError::Format(file.format));
        }
        if file.policy != policy {
            return Err(SliceError::Policy {
                expected: policy.to_owned(),
                found: file.policy,
            });
        }
        let slice = Slice::new(file.slice)?;
        let max_depth = slice.definition.max_depth;
        let max_cells = slice.definition.max_cells.get();
        if file.cells.len() > max_cells {
            return Err(SliceError::TooManyCells { max_cells });
        }
        for pair in file.cells.windows(2) {
            let key = |c: &Cell| (c.path.len(), c.path.clone());
            if key(&pair[0]) >= key(&pair[1]) {
                return Err(SliceError::Order {
                    path: pair[1].path.clone(),
                });
            }
        }
        for cell in &file.cells {
            let inside = slice.clipped(&slice.rectangle(&cell.path)?).is_some();
            let wrong = match cell.outcome {
                Outcome::Outside => inside,
                Outcome::Eliminated(_) => !inside,
                Outcome::Unresolved if cell.path.len() < max_depth => {
                    return Err(SliceError::EarlyUnresolved {
                        path: cell.path.clone(),
                    })
                }
                Outcome::Unresolved => !inside,
            };
            if wrong {
                return Err(SliceError::WrongOutside {
                    path: cell.path.clone(),
                });
            }
        }
        let leaves: BTreeSet<&str> = file.cells.iter().map(|c| c.path.as_str()).collect();
        let mut reached = 0;
        tile(&slice, &leaves, &mut String::new(), &mut reached)?;
        if reached != leaves.len() {
            let inner = file
                .cells
                .iter()
                .find(|c| (0..c.path.len()).any(|k| leaves.contains(&c.path[..k])))
                .expect("an unreached leaf lies below another leaf");
            return Err(SliceError::Overlap {
                path: inner.path.clone(),
            });
        }
        Ok(Self {
            slice,
            cells: file.cells,
        })
    }
}

/// Walks the tree from `node`: a leaf ends the walk; any other node must be
/// above the depth limit with a box meeting `B₀`, and both children must be
/// covered. Counts the leaves reached.
fn tile(
    slice: &Slice,
    leaves: &BTreeSet<&str>,
    node: &mut String,
    reached: &mut usize,
) -> Result<(), SliceError> {
    if leaves.contains(node.as_str()) {
        *reached += 1;
        return Ok(());
    }
    if node.len() == slice.definition.max_depth {
        return Err(SliceError::Gap { path: node.clone() });
    }
    if slice.clipped(&slice.rectangle(node)?).is_none() {
        return Err(SliceError::WrongOutside { path: node.clone() });
    }
    for bit in ['0', '1'] {
        node.push(bit);
        tile(slice, leaves, node, reached)?;
        node.pop();
    }
    Ok(())
}

/// `f` of every item on `threads` workers, in item order; the first error
/// in item order is returned.
fn parallel<T: Send>(
    items: &[String],
    threads: NonZeroUsize,
    f: impl Fn(&str) -> Result<T, SliceError> + Sync,
) -> Result<Vec<T>, SliceError> {
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<T, SliceError>>>> =
        items.iter().map(|_| Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads.get().min(items.len()) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                let Some(item) = items.get(i) else { break };
                let result = f(item);
                *results[i].lock().expect("no worker panicked") = Some(result);
            });
        }
    });
    results
        .into_iter()
        .map(|slot| slot.into_inner().expect("no worker panicked").expect("every item evaluated"))
        .collect()
}

/// The corners of a plot rectangle as a polygon.
pub fn polygon(rectangle: &[Interval; 2]) -> Vec<Point> {
    let [a, b] = rectangle;
    vec![
        [a.lo().clone(), b.lo().clone()],
        [a.hi().clone(), b.lo().clone()],
        [a.hi().clone(), b.hi().clone()],
        [a.lo().clone(), b.hi().clone()],
    ]
}

/// Compact JSON of a part of a partition file.
fn json<T: Serialize + ?Sized>(value: &T) -> String {
    serde_json::to_string(value).expect("plain data serialises")
}

#[cfg(test)]
mod tests;

