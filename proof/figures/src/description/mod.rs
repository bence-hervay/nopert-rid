//! Typed JSON figure descriptions: parsing (every field required, unknown
//! fields and duplicate names refused) and rendering into one SVG.
use crate::content::graph::{graph, markers, Graph, GraphError};
use crate::content::number::{field, interval, required, Rational};
use crate::content::scene::{Body, BoxScene, ConfigurationScene, SceneError};
use crate::content::slice::{polygon, Classifier, Definition, Partition, Slice, SliceError};
use crate::content::support::{support, SupportError};
use crate::content::Shape;
use crate::drawing::layout::{self, Canvas, Corner, LayoutError, Panel};
use crate::drawing::style::{outcomes, Color, PerOutcome, Style};
use crate::drawing::svg::Document;
use rid::arithmetic::exact::{frac, ExactError, Q};
use rid::problem::configuration::{ConfigurationBox, AXES};
use rid::search::BoxError;
use serde::de::{Error as _, MapAccess, Visitor};
use serde::{Deserialize, Deserializer};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::marker::PhantomData;
use std::num::{NonZeroU32, NonZeroUsize};
use std::path::{Path, PathBuf};

/// A figure: one canvas with named sources, styles and panels, drawn by an
/// ordered list of layers.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Description {
    pub png: Png,
    pub canvas: CanvasSpec,
    pub sources: Named<Source>,
    pub styles: Named<Style>,
    pub panels: Named<PanelSpec>,
    pub layers: Vec<Layer>,
}

/// The PNG preview: the converter's path and the resolution.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Png {
    pub converter: PathBuf,
    pub dpi: NonZeroU32,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CanvasSpec {
    pub width_mm: Rational,
    pub height_mm: Rational,
    /// `null` for a transparent page.
    #[serde(deserialize_with = "required")]
    pub background: Option<Color>,
}

/// What a figure shows.
#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Source {
    /// One configuration `(s, t, r₁, r₂, r₃)` of Q(√5)⁵.
    Configuration { coordinates: [[Rational; 2]; AXES] },
    /// Every configuration of a closed box.
    Box { axes: [[Rational; 2]; AXES] },
    /// A slice subdivided now, with the crate's components, on `threads`
    /// workers; its partition is written next to the figure.
    Slice { slice: Definition, threads: NonZeroUsize },
    /// A partition written by an earlier rendering, relative to the
    /// description's directory.
    Partition { file: PathBuf },
    /// A polynomial `Σ cₖ xᵏ` in one variable on `[l, h]`, drawn in the plane
    /// of `(x, p(x))` through `samples + 1` exact points.
    Polynomial {
        coefficients: Vec<Rational>,
        interval: [Rational; 2],
        samples: NonZeroUsize,
    },
    /// Explicit segments and polygons in the plane of `world`: guide lines
    /// over a scene or a slice, or a diagram of its own.
    Shapes { world: World, shapes: Vec<ShapeSpec> },
}

/// An explicit shape: `{"segment": [p, q]}`, `{"polygon": [p, q, …]}` or
/// `{"circle": {"centre": p, "radius": r}}`, each point a pair of rationals.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ShapeSpec {
    Segment([[Rational; 2]; 2]),
    Polygon(Vec<[Rational; 2]>),
    Circle { centre: [Rational; 2], radius: Rational },
}

impl ShapeSpec {
    fn shape(&self) -> Shape {
        let point = |p: &[Rational; 2]| [p[0].0.clone(), p[1].0.clone()];
        match self {
            ShapeSpec::Segment([p, q]) => Shape::Segment([point(p), point(q)]),
            ShapeSpec::Polygon(points) => Shape::Polygon(points.iter().map(point).collect()),
            ShapeSpec::Circle { centre, radius } => Shape::Circle { centre: point(centre), radius: radius.0.clone() },
        }
    }
}

/// The styles of the four parts of a support layer.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SupportStyles {
    /// The projected hole edge (an outline).
    pub edge: String,
    /// The supporting line (an outline).
    pub line: String,
    /// The outward normal (an outline).
    pub normal: String,
    /// The marker at the plug vertex.
    pub vertex: String,
}

/// A part of a polynomial's graph drawn as lines.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum GraphPart {
    /// The graph of the polynomial.
    Curve,
    /// The polygon through the Bernstein coefficients.
    ControlPolygon,
    /// The horizontal axis over the interval.
    Axis,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PanelSpec {
    /// `[[x₀, x₁], [y₀, y₁]]` in page millimetres, `y` pointing down.
    pub viewport_mm: [[Rational; 2]; 2],
    /// `[[x₀, x₁], [y₀, y₁]]` in world units, `y` pointing up.
    pub world: [[Rational; 2]; 2],
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Geometry {
    /// All edges at a configuration.
    Wireframe,
    /// The edges on the side of the body facing the viewer.
    FrontWireframe,
    /// The edges on the far side of the body, hidden by the body itself.
    BackWireframe,
    /// The shadow's boundary at a configuration.
    Silhouette,
    /// Over a box: a rectangle containing every position of each vertex.
    VertexEnvelopes,
    /// Over a box: a polygon containing every position of each edge.
    EdgeEnvelopes,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Connector {
    /// A corner of the window marker in the overview.
    pub from: Corner,
    /// A corner of the detail panel.
    pub to: Corner,
}

/// One drawing step, in order.
#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Layer {
    Scene {
        panel: String,
        source: String,
        body: Body,
        geometry: Geometry,
        style: String,
    },
    Cells {
        panel: String,
        source: String,
        styles: PerOutcome<String>,
    },
    Inset {
        overview: String,
        detail: String,
        connectors: Vec<Connector>,
        /// The smallest side of the window marker in millimetres: a
        /// smaller window is marked by a square of this side.
        mark_mm: Rational,
        style: String,
    },
    /// A support at a configuration: the ordered hole edge `edge`, the
    /// supporting line of the hole's shadow at the hole's extent in the
    /// direction `u × (v_b − v_a)`, the outward normal of length `length`
    /// and a dot of radius `marker` at the plug vertex `vertex`.
    Support {
        panel: String,
        source: String,
        edge: [usize; 2],
        vertex: usize,
        length: Rational,
        marker: Rational,
        styles: SupportStyles,
    },
    /// Lines of a polynomial's graph.
    Graph {
        panel: String,
        source: String,
        part: GraphPart,
        style: String,
    },
    /// Diamonds of half-diagonal `marker` at a polynomial's Bernstein
    /// coefficients.
    ControlPoints {
        panel: String,
        source: String,
        marker: Rational,
        style: String,
    },
    /// The shapes of a `shapes` source.
    Shapes {
        panel: String,
        source: String,
        style: String,
    },
    /// A small text in the top left corner of a panel, such as an inset's
    /// magnification.
    Label { panel: String, text: String },
}

/// Whether `name` is a name of a source, style or panel: `[a-z0-9][a-z0-9-]*`.
/// Source names become file names (`<name>.partition.json`), so a name is
/// never a path: no separator, no dot, nothing empty or hidden.
pub fn is_name(name: &str) -> bool {
    let mut bytes = name.bytes();
    bytes.next().is_some_and(|b| b.is_ascii_lowercase() || b.is_ascii_digit())
        && bytes.all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
}

/// Named entries of a JSON object; a repeated name, or one that is not
/// [`is_name`], is refused.
#[derive(Debug)]
pub struct Named<T>(pub BTreeMap<String, T>);

impl<'de, T: Deserialize<'de>> Deserialize<'de> for Named<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Entries<T>(PhantomData<T>);
        impl<'de, T: Deserialize<'de>> Visitor<'de> for Entries<T> {
            type Value = Named<T>;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("an object of named entries")
            }
            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Named<T>, A::Error> {
                let mut entries = BTreeMap::new();
                while let Some((name, value)) = map.next_entry::<String, T>()? {
                    if !is_name(&name) {
                        return Err(A::Error::custom(format!(
                            "name {name:?} is not a lowercase letter or digit followed by lowercase letters, digits and hyphens"
                        )));
                    }
                    if entries.contains_key(&name) {
                        return Err(A::Error::custom(format!("duplicate name {name:?}")));
                    }
                    entries.insert(name, value);
                }
                Ok(Named(entries))
            }
        }
        deserializer.deserialize_map(Entries(PhantomData))
    }
}

#[derive(Debug)]
pub enum DescriptionError {
    Json(serde_json::Error),
    Canvas(LayoutError),
    Panel { name: String, error: PanelError },
    /// A layer names a panel, source or style that is not defined.
    Unknown { layer: usize, what: &'static str, name: String },
    /// A source no layer draws: it would be computed for nothing.
    UnusedSource(String),
    /// A style no layer uses.
    UnusedStyle(String),
    /// A panel no layer uses.
    UnusedPanel(String),
    /// A panel that would show two different worlds (the screen of
    /// configurations and boxes, and the plot parameters of slices), or an
    /// inset between panels showing different worlds.
    WorldMismatch { layer: usize, panel: String },
    Source { name: String, error: SourceError },
    /// A layer draws a geometry its source does not have.
    Mismatch { layer: usize, source: String },
    /// A filled style on outlines (a wireframe or an inset).
    FilledOutline { layer: usize, style: String },
    Inset { layer: usize, error: LayoutError },
    /// A support layer that cannot be drawn at its configuration.
    Support { layer: usize, error: SupportError },
}

#[derive(Debug)]
pub enum PanelError {
    Number(ExactError),
    Layout(LayoutError),
}

#[derive(Debug)]
pub enum SourceError {
    Number(ExactError),
    Scene(SceneError),
    Slice(SliceError),
    Read { file: PathBuf, error: std::io::Error },
    Components(BoxError),
    Graph(GraphError),
}

impl fmt::Display for DescriptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DescriptionError::Json(e) => write!(f, "invalid description: {e}"),
            DescriptionError::Canvas(e) => write!(f, "canvas: {e}"),
            DescriptionError::Panel { name, error } => write!(f, "panel {name:?}: {error}"),
            DescriptionError::Unknown { layer, what, name } => {
                write!(f, "layer {layer}: no {what} named {name:?}")
            }
            DescriptionError::UnusedSource(name) => write!(f, "source {name:?} is not drawn"),
            DescriptionError::UnusedStyle(name) => write!(f, "style {name:?} is not used"),
            DescriptionError::UnusedPanel(name) => write!(f, "panel {name:?} is not used"),
            DescriptionError::WorldMismatch { layer, panel } => {
                write!(f, "layer {layer}: panel {panel:?} would show two different worlds")
            }
            DescriptionError::Source { name, error } => write!(f, "source {name:?}: {error}"),
            DescriptionError::Mismatch { layer, source } => {
                write!(f, "layer {layer}: source {source:?} has no such geometry")
            }
            DescriptionError::FilledOutline { layer, style } => {
                write!(f, "layer {layer}: style {style:?} has a fill, but the layer draws outlines")
            }
            DescriptionError::Inset { layer, error } => write!(f, "layer {layer}: {error}"),
            DescriptionError::Support { layer, error } => write!(f, "layer {layer}: {error}"),
        }
    }
}

impl fmt::Display for PanelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PanelError::Number(e) => write!(f, "{e}"),
            PanelError::Layout(e) => write!(f, "{e}"),
        }
    }
}

impl fmt::Display for SourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SourceError::Number(e) => write!(f, "{e}"),
            SourceError::Scene(e) => write!(f, "{e}"),
            SourceError::Slice(e) => write!(f, "{e}"),
            SourceError::Read { file, error } => {
                write!(f, "cannot read {}: {error}", file.display())
            }
            SourceError::Components(e) => write!(f, "components unavailable: {e}"),
            SourceError::Graph(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DescriptionError {}

/// A rendered figure.
pub struct Rendered {
    pub svg: String,
    /// The partitions of the `slice` sources, by source name.
    pub partitions: BTreeMap<String, Partition>,
}

/// The components' classifier, built only when a slice is subdivided.
pub type Components<'a> = &'a dyn Fn() -> Result<Box<dyn Classifier>, BoxError>;

impl Description {
    pub fn parse(text: &[u8]) -> Result<Self, DescriptionError> {
        serde_json::from_slice(text).map_err(DescriptionError::Json)
    }

    /// Every check that needs only the description, before any source is
    /// prepared: the canvas and panels, the names (every name defined and
    /// used), each layer's geometry against its source's kind, filled styles
    /// on outlines, inset windows inside their overviews, and one world per
    /// panel.
    pub fn check(&self) -> Result<(), DescriptionError> {
        self.canvas_and_panels().map(|_| ())
    }

    fn canvas_and_panels(&self) -> Result<(Canvas, BTreeMap<&str, Panel>), DescriptionError> {
        let (width, height) = (&self.canvas.width_mm.0, &self.canvas.height_mm.0);
        let canvas = Canvas::new(width.clone(), height.clone()).map_err(DescriptionError::Canvas)?;
        let mut panels = BTreeMap::new();
        for (name, spec) in &self.panels.0 {
            let panel = panel(&canvas, spec).map_err(|error| DescriptionError::Panel {
                name: name.clone(),
                error,
            })?;
            panels.insert(name.as_str(), panel);
        }
        self.check_names(&panels)?;
        self.check_layers(&panels)?;
        self.check_used()?;
        Ok((canvas, panels))
    }

    /// Checks the description ([`Description::check`]), prepares every
    /// source once, then draws the layers in order. Partition files are read
    /// relative to `base`; they must have been written under `policy`.
    pub fn render(
        &self,
        base: &Path,
        policy: &str,
        components: Components,
    ) -> Result<Rendered, DescriptionError> {
        let (canvas, panels) = self.canvas_and_panels()?;
        let mut classifier = None;
        let mut sources = BTreeMap::new();
        let mut partitions = BTreeMap::new();
        for (name, source) in &self.sources.0 {
            let prepared = prepare(source, base, policy, components, &mut classifier)
                .map_err(|error| DescriptionError::Source {
                    name: name.clone(),
                    error,
                })?;
            if let (Source::Slice { .. }, Prepared::Partition(p)) = (source, &prepared) {
                partitions.insert(name.clone(), p.clone());
            }
            sources.insert(name.as_str(), prepared);
        }
        let mut document = Document::new(&canvas, self.canvas.background);
        for (index, layer) in self.layers.iter().enumerate() {
            self.draw(index, layer, &panels, &sources, &mut document)?;
        }
        Ok(Rendered {
            svg: document.finish(),
            partitions,
        })
    }

    /// Every name a layer uses is defined.
    fn check_names(&self, panels: &BTreeMap<&str, Panel>) -> Result<(), DescriptionError> {
        for (layer, entry) in self.layers.iter().enumerate() {
            let unknown = |what, name: &String| DescriptionError::Unknown {
                layer,
                what,
                name: name.clone(),
            };
            let (panel_names, source, styles) = entry.names();
            if let Some(name) = panel_names.iter().find(|n| !panels.contains_key(n.as_str())) {
                return Err(unknown("panel", name));
            }
            if let Some(name) = styles.iter().find(|n| !self.styles.0.contains_key(**n)) {
                return Err(unknown("style", name));
            }
            if let Some(name) = source.filter(|name| !self.sources.0.contains_key(*name)) {
                return Err(unknown("source", name));
            }
        }
        Ok(())
    }

    /// Every source, style and panel is used by some layer.
    fn check_used(&self) -> Result<(), DescriptionError> {
        let (mut drawn, mut used_panels, mut used_styles) = (BTreeSet::new(), BTreeSet::new(), BTreeSet::new());
        for layer in &self.layers {
            let (panel_names, source, styles) = layer.names();
            used_panels.extend(panel_names);
            used_styles.extend(styles);
            drawn.extend(source);
        }
        if let Some(name) = self.sources.0.keys().find(|name| !drawn.contains(name)) {
            return Err(DescriptionError::UnusedSource(name.clone()));
        }
        if let Some(name) = self.styles.0.keys().find(|name| !used_styles.contains(name)) {
            return Err(DescriptionError::UnusedStyle(name.clone()));
        }
        match self.panels.0.keys().find(|name| !used_panels.contains(name)) {
            Some(name) => Err(DescriptionError::UnusedPanel(name.clone())),
            None => Ok(()),
        }
    }

    /// Each layer against its source's kind and its styles, inset windows,
    /// and the world of each panel: the screen of configurations and boxes,
    /// or the plot parameters of slices.
    fn check_layers(&self, panels: &BTreeMap<&str, Panel>) -> Result<(), DescriptionError> {
        let mut worlds: BTreeMap<&str, World> = BTreeMap::new();
        for (index, layer) in self.layers.iter().enumerate() {
            let mismatch = |source: &String| DescriptionError::Mismatch { layer: index, source: source.clone() };
            let filled = |name: &String| DescriptionError::FilledOutline { layer: index, style: name.clone() };
            match layer {
                Layer::Scene { panel, source, geometry, style, .. } => {
                    let fits = matches!(
                        (&self.sources.0[source], geometry),
                        (
                            Source::Configuration { .. },
                            Geometry::Wireframe | Geometry::FrontWireframe | Geometry::BackWireframe | Geometry::Silhouette
                        )
                            | (Source::Box { .. }, Geometry::VertexEnvelopes | Geometry::EdgeEnvelopes)
                    );
                    if !fits {
                        return Err(mismatch(source));
                    }
                    let edges = matches!(geometry, Geometry::Wireframe | Geometry::FrontWireframe | Geometry::BackWireframe);
                    if edges && self.styles.0[style].fill().is_some() {
                        return Err(filled(style));
                    }
                    shows(&mut worlds, index, panel, World::Screen)?;
                }
                Layer::Cells { panel, source, .. } => {
                    if !matches!(self.sources.0[source], Source::Slice { .. } | Source::Partition { .. }) {
                        return Err(mismatch(source));
                    }
                    shows(&mut worlds, index, panel, World::Plot)?;
                }
                Layer::Support { panel, source, styles, .. } => {
                    if !matches!(self.sources.0[source], Source::Configuration { .. }) {
                        return Err(mismatch(source));
                    }
                    for name in [&styles.edge, &styles.line, &styles.normal] {
                        if self.styles.0[name].fill().is_some() {
                            return Err(filled(name));
                        }
                    }
                    shows(&mut worlds, index, panel, World::Screen)?;
                }
                Layer::Graph { panel, source, style, .. } => {
                    if !matches!(self.sources.0[source], Source::Polynomial { .. }) {
                        return Err(mismatch(source));
                    }
                    if self.styles.0[style].fill().is_some() {
                        return Err(filled(style));
                    }
                    shows(&mut worlds, index, panel, World::Graph)?;
                }
                Layer::ControlPoints { panel, source, .. } => {
                    if !matches!(self.sources.0[source], Source::Polynomial { .. }) {
                        return Err(mismatch(source));
                    }
                    shows(&mut worlds, index, panel, World::Graph)?;
                }
                Layer::Shapes { panel, source, .. } => {
                    let Source::Shapes { world, .. } = &self.sources.0[source] else {
                        return Err(mismatch(source));
                    };
                    shows(&mut worlds, index, panel, *world)?;
                }
                Layer::Label { .. } => {}
                Layer::Inset { overview, detail, connectors, mark_mm, style } => {
                    if self.styles.0[style].fill().is_some() {
                        return Err(filled(style));
                    }
                    let pairs: Vec<(Corner, Corner)> = connectors.iter().map(|c| (c.from, c.to)).collect();
                    layout::inset(&panels[overview.as_str()], &panels[detail.as_str()], &pairs, &mark_mm.0)
                        .map_err(|error| DescriptionError::Inset { layer: index, error })?;
                }
            }
        }
        for (index, layer) in self.layers.iter().enumerate() {
            if let Layer::Inset { overview, detail, .. } = layer {
                let (a, b) = (worlds.get(overview.as_str()), worlds.get(detail.as_str()));
                if a.is_some() && b.is_some() && a != b {
                    return Err(DescriptionError::WorldMismatch { layer: index, panel: detail.clone() });
                }
            }
        }
        Ok(())
    }

    fn draw(
        &self,
        index: usize,
        layer: &Layer,
        panels: &BTreeMap<&str, Panel>,
        sources: &BTreeMap<&str, Prepared>,
        document: &mut Document,
    ) -> Result<(), DescriptionError> {
        // The kinds, styles and insets were checked by `check_layers`.
        let style = |name: &String| &self.styles.0[name];
        match layer {
            Layer::Scene {
                panel,
                source,
                body,
                geometry,
                style: name,
            } => {
                let shapes = match (&sources[source.as_str()], geometry) {
                    (Prepared::Configuration(c), Geometry::Wireframe) => c.wireframe(*body),
                    (Prepared::Configuration(c), Geometry::FrontWireframe) => c.front_wireframe(*body),
                    (Prepared::Configuration(c), Geometry::BackWireframe) => c.back_wireframe(*body),
                    (Prepared::Configuration(c), Geometry::Silhouette) => vec![c.outline(*body)],
                    (Prepared::Box(b), Geometry::VertexEnvelopes) => b.vertex_envelopes(*body),
                    (Prepared::Box(b), Geometry::EdgeEnvelopes) => b.edge_envelopes(*body),
                    _ => unreachable!("the geometry was checked against the source's kind"),
                };
                draw_in(document, &panels[panel.as_str()], &shapes, style(name));
            }
            Layer::Cells {
                panel,
                source,
                styles,
            } => {
                let Prepared::Partition(partition) = &sources[source.as_str()] else {
                    unreachable!("the source's kind was checked")
                };
                for outcome in outcomes() {
                    let shapes: Vec<Shape> = partition
                        .cells()
                        .iter()
                        .filter(|cell| cell.outcome == outcome)
                        .map(|cell| {
                            let rectangle = partition.slice().rectangle(&cell.path);
                            Shape::Polygon(polygon(&rectangle.expect("a valid path")))
                        })
                        .collect();
                    if !shapes.is_empty() {
                        let panel = &panels[panel.as_str()];
                        draw_in(document, panel, &shapes, style(styles.get(&outcome)));
                    }
                }
            }
            Layer::Support {
                panel,
                source,
                edge,
                vertex,
                length,
                marker,
                styles,
            } => {
                let Prepared::Configuration(scene) = &sources[source.as_str()] else {
                    unreachable!("the source's kind was checked")
                };
                let drawn = support(scene, *edge, *vertex, &length.0, &marker.0)
                    .map_err(|error| DescriptionError::Support { layer: index, error })?;
                let panel = &panels[panel.as_str()];
                draw_in(document, panel, &drawn.line, style(&styles.line));
                draw_in(document, panel, &drawn.edge, style(&styles.edge));
                draw_in(document, panel, &drawn.normal, style(&styles.normal));
                draw_in(document, panel, &drawn.vertex, style(&styles.vertex));
            }
            Layer::Graph {
                panel,
                source,
                part,
                style: name,
            } => {
                let Prepared::Graph(g) = &sources[source.as_str()] else {
                    unreachable!("the source's kind was checked")
                };
                let shapes = match part {
                    GraphPart::Curve => &g.curve,
                    GraphPart::ControlPolygon => &g.polygon,
                    GraphPart::Axis => &g.axis,
                };
                draw_in(document, &panels[panel.as_str()], shapes, style(name));
            }
            Layer::ControlPoints {
                panel,
                source,
                marker,
                style: name,
            } => {
                let Prepared::Graph(g) = &sources[source.as_str()] else {
                    unreachable!("the source's kind was checked")
                };
                draw_in(document, &panels[panel.as_str()], &markers(&g.points, &marker.0), style(name));
            }
            Layer::Shapes {
                panel,
                source,
                style: name,
            } => {
                let Prepared::Shapes(shapes) = &sources[source.as_str()] else {
                    unreachable!("the source's kind was checked")
                };
                draw_in(document, &panels[panel.as_str()], shapes, style(name));
            }
            Layer::Label { panel, text } => {
                let [x, y] = panels[panel.as_str()].viewport();
                document.label([x.lo() + &frac(9, 10), y.lo() + &frac(12, 5)], &frac(19, 10), text);
            }
            Layer::Inset {
                overview,
                detail,
                connectors,
                mark_mm,
                style: name,
            } => {
                let pairs: Vec<(Corner, Corner)> =
                    connectors.iter().map(|c| (c.from, c.to)).collect();
                let (overview, detail) = (&panels[overview.as_str()], &panels[detail.as_str()]);
                let inset = layout::inset(overview, detail, &pairs, &mark_mm.0)
                    .map_err(|error| DescriptionError::Inset { layer: index, error })?;
                let mut shapes = vec![inset.window, inset.frame];
                shapes.extend(inset.connectors);
                document.layer(&shapes, style(name), None);
            }
        }
        Ok(())
    }
}

impl Layer {
    /// The panels, the source and the styles a layer names.
    fn names(&self) -> (Vec<&String>, Option<&String>, Vec<&String>) {
        match self {
            Layer::Scene {
                panel,
                source,
                style,
                ..
            } => (vec![panel], Some(source), vec![style]),
            Layer::Cells {
                panel,
                source,
                styles,
            } => {
                let names = outcomes().iter().map(|o| styles.get(o)).collect();
                (vec![panel], Some(source), names)
            }
            Layer::Inset {
                overview,
                detail,
                style,
                ..
            } => (vec![overview, detail], None, vec![style]),
            Layer::Support {
                panel,
                source,
                styles,
                ..
            } => (vec![panel], Some(source), vec![&styles.edge, &styles.line, &styles.normal, &styles.vertex]),
            Layer::Graph {
                panel,
                source,
                style,
                ..
            }
            | Layer::ControlPoints {
                panel,
                source,
                style,
                ..
            }
            | Layer::Shapes {
                panel,
                source,
                style,
            } => (vec![panel], Some(source), vec![style]),
            Layer::Label { panel, .. } => (vec![panel], None, vec![]),
        }
    }
}

/// Records that layer `layer` draws `world` in `panel`; a panel keeps one world.
fn shows<'a>(
    worlds: &mut BTreeMap<&'a str, World>,
    layer: usize,
    panel: &'a str,
    world: World,
) -> Result<(), DescriptionError> {
    if *worlds.entry(panel).or_insert(world) != world {
        return Err(DescriptionError::WorldMismatch { layer, panel: panel.to_owned() });
    }
    Ok(())
}

/// What a panel shows: the screen of configurations and boxes, or the plot
/// parameters `(a, b)` of slices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum World {
    Screen,
    Plot,
    /// The plane `(x, p(x))` of a polynomial's graph.
    Graph,
    /// A plane of its own, drawn only by `shapes` layers.
    Diagram,
}

/// A source ready to be drawn.
enum Prepared {
    Configuration(ConfigurationScene),
    Box(BoxScene),
    Partition(Partition),
    Graph(Graph),
    Shapes(Vec<Shape>),
}

fn prepare(
    source: &Source,
    base: &Path,
    policy: &str,
    components: Components,
    classifier: &mut Option<Box<dyn Classifier>>,
) -> Result<Prepared, SourceError> {
    Ok(match source {
        Source::Configuration { coordinates } => {
            let point = std::array::from_fn(|j| field(&coordinates[j]));
            Prepared::Configuration(ConfigurationScene::new(&point).map_err(SourceError::Scene)?)
        }
        Source::Box { axes } => {
            let axes: Vec<_> = axes
                .iter()
                .map(interval)
                .collect::<Result<_, _>>()
                .map_err(SourceError::Number)?;
            let b = ConfigurationBox::new(axes.try_into().expect("five axes"));
            Prepared::Box(BoxScene::new(&b).map_err(SourceError::Scene)?)
        }
        Source::Slice { slice, threads } => {
            let slice = Slice::new(slice.clone()).map_err(SourceError::Slice)?;
            if classifier.is_none() {
                *classifier = Some(components().map_err(SourceError::Components)?);
            }
            let classifier = classifier.as_deref().expect("built above");
            let partition = Partition::compute(&slice, classifier, *threads);
            Prepared::Partition(partition.map_err(SourceError::Slice)?)
        }
        Source::Partition { file } => {
            let path = base.join(file);
            let bytes = std::fs::read(&path).map_err(|error| SourceError::Read {
                file: path,
                error,
            })?;
            Prepared::Partition(Partition::from_json(&bytes, policy).map_err(SourceError::Slice)?)
        }
        Source::Polynomial {
            coefficients,
            interval: [l, h],
            samples,
        } => {
            let c: Vec<Q> = coefficients.iter().map(|x| x.0.clone()).collect();
            Prepared::Graph(graph(&c, [l.0.clone(), h.0.clone()], samples.get()).map_err(SourceError::Graph)?)
        }
        Source::Shapes { shapes, .. } => Prepared::Shapes(shapes.iter().map(ShapeSpec::shape).collect()),
    })
}

fn panel(canvas: &Canvas, spec: &PanelSpec) -> Result<Panel, PanelError> {
    let rectangle = |pairs: &[[Rational; 2]; 2]| -> Result<layout::Rectangle, PanelError> {
        Ok([
            interval(&pairs[0]).map_err(PanelError::Number)?,
            interval(&pairs[1]).map_err(PanelError::Number)?,
        ])
    };
    let (viewport, world) = (rectangle(&spec.viewport_mm)?, rectangle(&spec.world)?);
    Panel::new(canvas, viewport, world).map_err(PanelError::Layout)
}

/// Draws world shapes in a panel: cut one millimetre beyond the reach of the
/// stroke (see the layout's clipping) and clipped to the viewport.
fn draw_in(document: &mut Document, panel: &Panel, shapes: &[Shape], style: &Style) {
    let reach = style.stroke().map_or(Q::zero(), |s| s.width_mm().clone());
    let placed = panel.place(shapes, &(reach + Q::one()));
    document.layer(&placed, style, Some(panel.viewport()));
}

#[cfg(test)]
mod tests;

