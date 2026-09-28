//! `paper-figures <figures-directory> [<name>…]`: the article's figures.
//! For each figure (all, or those named) it writes the description
//! `source/<name>.json`, subdivides the slices it needs into
//! `source/partitions/` (a partition file whose slice is unchanged is
//! reused), renders `<name>.svg` and converts it to `<name>.pdf` with
//! `rsvg-convert`, found through `PATH`. With `RID_FIGURES_PREVIEW=<dir>`,
//! a PNG preview of each figure is written there too.
//!
//! Floating point only chooses what to show (inset windows, examples); every
//! drawn position is computed exactly by the library.
use num_traits::ToPrimitive;
use rid::arithmetic::exact::{frac, q, QSqrt5, Q};
use rid::problem::configuration::{ConfigurationBox, AXES};
use rid::problem::geometry;
use rid::search::BoxError;
use rid_figures::content::slice::{Classifier, Definition, Parts, Partition, Slice};
use rid_figures::description::Description;
use serde_json::{json, Value};
use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

mod diagrams;
mod shadows;
mod slices;

/// Worker threads for subdivisions and cover verification.
const THREADS: usize = 8;

/// A figure: its description, and the slices its partition sources name.
pub struct Figure {
    pub description: Value,
    pub slices: Vec<(String, Definition)>,
}

impl Figure {
    pub fn plain(description: Value) -> Self {
        Self { description, slices: Vec::new() }
    }
}

// ---- Exact numbers as description values ----------------------------------

pub fn rat(x: &Q) -> Value {
    Value::String(x.to_string())
}

pub fn fr(n: i64, d: i64) -> Q {
    frac(n, d)
}

/// `[a, b]` for `a + b√5`.
pub fn pair(x: &QSqrt5) -> Value {
    json!([rat(x.rational_part()), rat(x.sqrt5_part())])
}

pub fn interval(lo: &Q, hi: &Q) -> Value {
    json!([rat(lo), rat(hi)])
}

pub fn point(p: &[Q; 2]) -> Value {
    json!([rat(&p[0]), rat(&p[1])])
}

pub fn to_f64(x: &Q) -> f64 {
    x.numer().to_f64().expect("finite") / x.denom().to_f64().expect("finite")
}

pub fn field_f64(x: &QSqrt5) -> f64 {
    to_f64(x.rational_part()) + to_f64(x.sqrt5_part()) * 5f64.sqrt()
}

/// The nearest multiple of `1/d` to `x`.
pub fn grid(x: f64, d: i64) -> Q {
    frac((x * d as f64).round() as i64, d)
}

/// A configuration source at an exact point of Q(√5)⁵.
pub fn configuration(c: &[QSqrt5; AXES]) -> Value {
    json!({"kind": "configuration", "coordinates": c.iter().map(pair).collect::<Vec<_>>()})
}

pub fn rational_configuration(c: &[Q; AXES]) -> [QSqrt5; AXES] {
    std::array::from_fn(|j| QSqrt5::from_rational(c[j].clone()))
}

/// The midpoint of the search box of a certificate path.
pub fn midpoint(path: &str) -> [Q; AXES] {
    ConfigurationBox::from_path(path).expect("a valid path").midpoint()
}

pub fn style(fill: Option<&str>, stroke: Option<(&str, Q)>, opacity: Q) -> Value {
    json!({
        "fill": fill,
        "stroke": stroke.map(|(color, width)| json!({"color": color, "width_mm": rat(&width)})),
        "opacity": rat(&opacity),
    })
}

pub fn line(color: &str, width: Q) -> Value {
    style(None, Some((color, width)), q(1))
}

/// A line of dashes and gaps `dash` millimetres long.
pub fn dashed(color: &str, width: Q, dash: Q) -> Value {
    json!({
        "fill": null,
        "stroke": {"color": color, "width_mm": rat(&width), "dash_mm": rat(&dash)},
        "opacity": "1",
    })
}

pub fn fill(color: &str) -> Value {
    style(Some(color), None, q(1))
}

pub fn panel(x: [Q; 2], y: [Q; 2], world: [[Q; 2]; 2]) -> Value {
    json!({
        "viewport_mm": [interval(&x[0], &x[1]), interval(&y[0], &y[1])],
        "world": [interval(&world[0][0], &world[0][1]), interval(&world[1][0], &world[1][1])],
    })
}

pub fn description(width: Q, height: Q, sources: Value, styles: Value, panels: Value, layers: Vec<Value>) -> Value {
    json!({
        "png": {"converter": "rsvg-convert", "dpi": 150},
        "canvas": {"width_mm": rat(&width), "height_mm": rat(&height), "background": "#ffffff"},
        "sources": sources,
        "styles": styles,
        "panels": panels,
        "layers": layers,
    })
}

/// A `shapes` source.
pub fn shapes(world: &str, shapes: Vec<Value>) -> Value {
    json!({"kind": "shapes", "world": world, "shapes": shapes})
}

pub fn segment(p: [Q; 2], q: [Q; 2]) -> Value {
    json!({"segment": [point(&p), point(&q)]})
}

pub fn polygon(points: &[[Q; 2]]) -> Value {
    json!({"polygon": points.iter().map(point).collect::<Vec<_>>()})
}

pub fn shapes_layer(panel: &str, source: &str, style: &str) -> Value {
    json!({"kind": "shapes", "panel": panel, "source": source, "style": style})
}

pub fn all_parts() -> Parts {
    Parts::all()
}

// ---- The figures -----------------------------------------------------------

type Builder = fn() -> Figure;

fn figures() -> Vec<(&'static str, Builder)> {
    vec![
        ("global", shadows::global as Builder),
        ("local", shadows::local),
        ("square-pentagon", shadows::square_pentagon),
        ("arcs", shadows::arcs),
        ("endpoint-crossing", shadows::endpoint_crossing),
        ("bernstein", diagrams::bernstein),
        ("slice-aligned", slices::slice_aligned),
        ("slice-arcs", slices::slice_arcs),
        ("hierarchy", slices::hierarchy),
        ("global-refinement", slices::global_refinement),
        ("local-refinement", slices::local_refinement),
        ("exotic-directions", slices::exotic_directions),
        ("arcs-remaining", slices::arcs_remaining),
        ("arcs-residual", slices::arcs_residual),
        ("toy-zoom", diagrams::toy_zoom),
        ("pentagon-support", shadows::pentagon_support),
        ("crossing-plane", diagrams::crossing_plane),
        ("local-tiling", diagrams::local_tiling),
    ]
}

// ---- Running ---------------------------------------------------------------

fn fail(what: impl std::fmt::Display) -> String {
    what.to_string()
}

/// The production classifier, built at the first subdivision.
struct Components(OnceCell<Box<dyn Classifier>>);

impl Components {
    fn get(&self) -> Result<&dyn Classifier, BoxError> {
        if self.0.get().is_none() {
            eprintln!("verifying the zoom covers");
            let classifier = rid_figures::classifier::production(threads())?;
            let _ = self.0.set(classifier);
        }
        Ok(self.0.get().expect("set above").as_ref())
    }
}

fn threads() -> NonZeroUsize {
    NonZeroUsize::new(THREADS).expect("positive")
}

/// Writes `source/partitions/<name>.partition.json` unless a file with the
/// same slice, written under the crate's policy, is already there.
fn partition(directory: &Path, name: &str, definition: &Definition, components: &Components) -> Result<(), String> {
    let file = directory.join("source/partitions").join(format!("{name}.partition.json"));
    if let Ok(bytes) = fs::read(&file) {
        if let Ok(existing) = Partition::from_json(&bytes, rid::POLICY) {
            if existing.slice().definition() == definition {
                return Ok(());
            }
        }
    }
    eprintln!("subdividing {name}");
    let slice = Slice::new(definition.clone()).map_err(fail)?;
    let classifier = components.get().map_err(fail)?;
    let computed = Partition::compute(&slice, classifier, threads()).map_err(|e| format!("{name}: {e}"))?;
    fs::create_dir_all(file.parent().expect("a parent")).map_err(fail)?;
    fs::write(&file, computed.to_json(rid::POLICY)).map_err(fail)
}

fn convert(format: &str, svg: &Path, output: &Path) -> Result<(), String> {
    let status = Command::new("rsvg-convert")
        .args(["-f", format, "--dpi-x", "150", "--dpi-y", "150", "-o"])
        .arg(output)
        .arg(svg)
        .status()
        .map_err(|e| format!("rsvg-convert: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("rsvg-convert failed on {}", svg.display()))
    }
}

fn render(directory: &Path, name: &str, figure: Figure, components: &Components) -> Result<(), String> {
    let source = directory.join("source");
    for (slice, definition) in &figure.slices {
        partition(directory, slice, definition, components)?;
    }
    let text = serde_json::to_string_pretty(&figure.description).map_err(fail)? + "\n";
    let file = source.join(format!("{name}.json"));
    fs::write(&file, &text).map_err(fail)?;
    let parsed = Description::parse(text.as_bytes()).map_err(|e| format!("{name}: {e}"))?;
    let no_slices = || -> Result<Box<dyn Classifier>, BoxError> {
        unreachable!("the figures read their slices from partition files")
    };
    let rendered = parsed.render(&source, rid::POLICY, &no_slices).map_err(|e| format!("{name}: {e}"))?;
    let svg = directory.join(format!("{name}.svg"));
    let pdf = directory.join(format!("{name}.pdf"));
    let unchanged = fs::read_to_string(&svg).is_ok_and(|old| old == rendered.svg) && pdf.exists();
    if !unchanged {
        fs::write(&svg, &rendered.svg).map_err(fail)?;
        convert("pdf", &svg, &pdf)?;
    }
    if let Some(preview) = std::env::var_os("RID_FIGURES_PREVIEW") {
        let preview = PathBuf::from(preview);
        fs::create_dir_all(&preview).map_err(fail)?;
        convert("png", &svg, &preview.join(format!("{name}.png")))?;
    }
    println!("{name}{}", if unchanged { " (unchanged)" } else { "" });
    Ok(())
}

fn run(arguments: &[String]) -> Result<(), String> {
    let Some((directory, names)) = arguments.split_first() else {
        return Err("usage: paper-figures <figures-directory> [<name>…]".into());
    };
    let directory = Path::new(directory);
    let known: BTreeMap<&str, Builder> = figures().into_iter().collect();
    if let Some(unknown) = names.iter().find(|n| !known.contains_key(n.as_str())) {
        return Err(format!("no figure named {unknown:?}"));
    }
    let components = Components(OnceCell::new());
    for (name, build) in figures() {
        if names.is_empty() || names.iter().any(|n| n == name) {
            render(directory, name, build(), &components)?;
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    match run(&arguments) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("paper-figures: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Whether a configuration lies in the reduced domain `D`: in `B₀`, and every
/// inequality of `D` holds (exactly).
pub fn in_domain(c: &[QSqrt5; AXES]) -> bool {
    let root = ConfigurationBox::root();
    let inside = root.axes().iter().zip(c).all(|(axis, x)| {
        QSqrt5::from_rational(axis.lo().clone()) <= *x && *x <= QSqrt5::from_rational(axis.hi().clone())
    });
    inside && rid::problem::domain::affine_polynomials().iter().all(|p| p.evaluate(c).sign() != std::cmp::Ordering::Greater)
}

/// The RID's vertices in floating point, in the crate's order.
pub fn vertices_f64() -> Vec<[f64; 3]> {
    geometry::vertices().iter().map(|v| std::array::from_fn(|j| field_f64(&v[j]))).collect()
}
