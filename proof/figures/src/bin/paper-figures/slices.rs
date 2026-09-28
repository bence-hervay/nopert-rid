//! Slices of configuration space, each cell coloured by the part of the
//! collection that eliminates it; unresolved cells at the depth limit dark.
use super::*;
use rid::components::record::ExoticCover;

/// The colours of the components and covers, and of unresolved cells.
pub fn slice_styles(edge: Q) -> Value {
    let cell = |color: &str| style(Some(color), Some(("#ffffff", edge.clone())), q(1));
    json!({
        "domain": cell("#bdbdbd"),
        "global": cell("#9ecae1"),
        "local": cell("#74c476"),
        "square": cell("#fdae6b"),
        "pentagon": cell("#e6550d"),
        "arc": cell("#9e9ac8"),
        "endpoint": cell("#fa9fb5"),
        "crossing": cell("#c51b8a"),
        // Grey edges show the grid of unresolved cells without washing out
        // the thinnest ones.
        "unresolved": style(Some("#252525"), Some(("#8c8c8c", edge.clone() / q(2))), q(1)),
        "outside": cell("#f7f7f7"),
    })
}

pub fn cells_layer(panel: &str, source: &str) -> Value {
    json!({
        "kind": "cells", "panel": panel, "source": source,
        "styles": {
            "Domain": "domain", "Global": "global", "Local": "local",
            "Exotic": {"square": "square", "pentagon": "pentagon", "arc+": "arc", "arc-": "arc",
                "endpoint+": "endpoint", "endpoint-": "endpoint", "crossing+": "crossing", "crossing-": "crossing"},
            "unresolved": "unresolved", "outside": "outside",
        },
    })
}

pub fn definition(origin: [Q; AXES], horizontal: [Q; AXES], vertical: [Q; AXES], plot: [[Q; 2]; 2], thickness: [Q; AXES], max_depth: usize, max_cells: usize, parts: Parts) -> Definition {
    let value = json!({
        "origin": origin.iter().map(rat).collect::<Vec<_>>(),
        "horizontal": horizontal.iter().map(rat).collect::<Vec<_>>(),
        "vertical": vertical.iter().map(rat).collect::<Vec<_>>(),
        "plot": [interval(&plot[0][0], &plot[0][1]), interval(&plot[1][0], &plot[1][1])],
        "thickness": thickness.iter().map(rat).collect::<Vec<_>>(),
        "max_depth": max_depth,
        "max_cells": max_cells,
        "components": parts,
    });
    serde_json::from_value(value).expect("a valid slice definition")
}

fn partition_source(name: &str) -> Value {
    json!({"kind": "partition", "file": format!("partitions/{name}.partition.json")})
}

/// The final collection's slice through `t = 1/10`, `r₁ = r₂ = 0`.
fn aligned_definition() -> Definition {
    let zero = || q(0);
    definition(
        [zero(), fr(1, 10), zero(), zero(), zero()],
        [q(1), zero(), zero(), zero(), zero()],
        [zero(), zero(), zero(), zero(), q(1)],
        [[q(0), fr(2, 3)], [fr(-2, 5), fr(2, 5)]],
        std::array::from_fn(|_| zero()),
        16,
        65536,
        Parts::all(),
    )
}

/// `r_A` rounded to rationals within `10⁻¹⁴`: `(a₁, 0, a₃)`.
pub fn arc_rationals() -> (Q, Q) {
    (fr(726103, 4250681), fr(-2178309, 7881196))
}

/// The final collection's slice through the arc plane `Π₊`: `s` and `θ`.
fn arcs_definition() -> Definition {
    let (a1, a3) = arc_rationals();
    let zero = || q(0);
    definition(
        [zero(), zero(), a1.clone(), zero(), a3.clone()],
        [q(1), zero(), zero(), zero(), zero()],
        [zero(), zero(), a3, q(1), -a1],
        [[fr(3, 16), fr(5, 8)], [fr(-5, 16), fr(1, 8)]],
        std::array::from_fn(|_| zero()),
        18,
        65536,
        Parts::all(),
    )
}

fn single(name: &str, definition: Definition, width: Q, height: Q) -> Figure {
    let plot = definition.plot.clone();
    let world = [[plot[0][0].0.clone(), plot[0][1].0.clone()], [plot[1][0].0.clone(), plot[1][1].0.clone()]];
    let panels = json!({"plot": panel([q(5), &width - q(5)], [q(5), &height - q(5)], world)});
    let description = description(width, height, json!({"slice": partition_source(name)}), slice_styles(fr(1, 40)), panels, vec![cells_layer("plot", "slice")]);
    Figure { description, slices: vec![(name.into(), definition)] }
}

pub fn slice_aligned() -> Figure {
    single("slice-aligned", aligned_definition(), q(100), q(118))
}

pub fn slice_arcs() -> Figure {
    single("slice-arcs", arcs_definition(), q(100), q(100))
}

/// The Domain and Global components, then Local, then everything.
pub fn stages() -> [Parts; 3] {
    let dg = Parts { domain: true, global: true, local: false, exotic: vec![] };
    let dgl = Parts { local: true, ..dg.clone() };
    [dg, dgl, Parts::all()]
}

/// Domain, Global, Local and some Exotic covers.
fn with_covers(covers: &[ExoticCover]) -> Parts {
    Parts { exotic: covers.to_vec(), ..stages()[1].clone() }
}

/// The slice with other parts of the collection.
fn with(definition: &Definition, parts: Parts) -> Definition {
    Definition { components: parts, ..definition.clone() }
}

fn world(definition: &Definition) -> [[Q; 2]; 2] {
    let p = &definition.plot;
    [[p[0][0].0.clone(), p[0][1].0.clone()], [p[1][0].0.clone(), p[1][1].0.clone()]]
}

/// The slice restricted to a window of its plot, `depth` halvings deep.
fn window(definition: &Definition, plot: [[Q; 2]; 2], depth: usize) -> Definition {
    let value = json!([interval(&plot[0][0], &plot[0][1]), interval(&plot[1][0], &plot[1][1])]);
    Definition { plot: serde_json::from_value(value).expect("a plot"), max_depth: depth, ..definition.clone() }
}

/// Square windows of half-widths `halves` around `centre`.
fn zooms(definition: &Definition, centre: [Q; 2], halves: &[Q], depth: usize) -> Vec<Definition> {
    halves
        .iter()
        .map(|h| window(definition, [[&centre[0] - h, &centre[0] + h], [&centre[1] - h, &centre[1] + h]], depth))
        .collect()
}

/// A panel of a grid: a named slice and guide shapes drawn over it.
pub struct Cell {
    pub name: String,
    pub definition: Definition,
    pub guides: Vec<(Value, &'static str)>,
}

pub fn cell(name: String, definition: Definition) -> Cell {
    Cell { name, definition, guides: Vec::new() }
}

/// Rows of slices, each row `height` millimetres high with equal scales on
/// both axes of every panel, panels `gap` apart and rows centred; `insets`
/// are `(row, overview, detail, connected)` panel links.
pub fn grid_figure(rows: Vec<Vec<Cell>>, height: Q, gap: Q, insets: &[(usize, usize, usize, bool)], guide_styles: Value) -> Figure {
    let mut sources = serde_json::Map::new();
    let mut panels = serde_json::Map::new();
    let mut layers = Vec::new();
    let mut slices = Vec::new();
    let widths: Vec<Vec<Q>> = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|cell| {
                    let w = world(&cell.definition);
                    &height * &(&w[0][1] - &w[0][0]) / (&w[1][1] - &w[1][0])
                })
                .collect()
        })
        .collect();
    let row_width = |row: &Vec<Q>| row.iter().fold(Q::zero(), |sum, w| &sum + w) + &gap * &q(row.len() as i64 - 1);
    let total = widths.iter().map(row_width).max().expect("a row") + &gap * &q(2);
    let mut guides = 0;
    for (i, row) in rows.iter().enumerate() {
        let y = &gap + &(&(&height + &gap) * &q(i as i64));
        let mut x = (&total - &row_width(&widths[i])) / q(2);
        for (j, cell) in row.iter().enumerate() {
            let panel_name = format!("p{i}-{j}");
            let w = &widths[i][j];
            panels.insert(panel_name.clone(), panel([x.clone(), &x + w], [y.clone(), &y + &height], world(&cell.definition)));
            x = &(&x + w) + &gap;
            sources.insert(cell.name.clone(), partition_source(&cell.name));
            layers.push(cells_layer(&panel_name, &cell.name));
            for (source, style) in &cell.guides {
                let name = format!("guide-{guides}");
                guides += 1;
                sources.insert(name.clone(), source.clone());
                layers.push(shapes_layer(&panel_name, &name, style));
            }
            slices.push((cell.name.clone(), cell.definition.clone()));
        }
    }
    for &(i, overview, detail, connected) in insets {
        let connectors = if connected {
            json!([{"from": "top-right", "to": "top-left"}, {"from": "bottom-right", "to": "bottom-left"}])
        } else {
            json!([])
        };
        layers.push(json!({
            "kind": "inset", "overview": format!("p{i}-{overview}"), "detail": format!("p{i}-{detail}"),
            "connectors": connectors, "mark_mm": "6/5", "style": "inset",
        }));
    }
    let total_height = &gap + &(&(&height + &gap) * &q(rows.len() as i64));
    let mut styles = slice_styles(fr(1, 40));
    let extra = styles.as_object_mut().expect("an object");
    if !insets.is_empty() {
        extra.insert("inset".into(), line("#404040", fr(1, 6)));
    }
    extra.extend(guide_styles.as_object().expect("an object").clone());
    Figure { description: description(total, total_height, Value::Object(sources), styles, Value::Object(panels), layers), slices }
}

/// The plane `t = 0`, `r₁ = r₂ = 0` through the square configuration:
/// `s` and `r₃`.
fn square_plane() -> Definition {
    let zero = || q(0);
    definition(
        std::array::from_fn(|_| zero()),
        [q(1), zero(), zero(), zero(), zero()],
        [zero(), zero(), zero(), zero(), q(1)],
        [[q(0), fr(2, 3)], [fr(-2, 5), fr(2, 5)]],
        std::array::from_fn(|_| zero()),
        16,
        1 << 17,
        Parts::all(),
    )
}

pub fn hierarchy() -> Figure {
    let rows = [("aligned", square_plane()), ("arcs", arcs_definition())]
        .into_iter()
        .map(|(name, definition)| {
            stages().into_iter().enumerate().map(|(k, parts)| cell(format!("hierarchy-{name}-{k}"), with(&definition, parts))).collect()
        })
        .collect();
    grid_figure(rows, q(46), q(4), &[], json!({}))
}

/// Two ordinary aligned touch configurations, each in a plane of one view
/// and one rotation coordinate, in three successive windows.
fn ordinary(name: &str, parts: Parts) -> Figure {
    let zero = || q(0);
    let first = definition(
        [zero(), fr(1, 10), zero(), zero(), zero()],
        [q(1), zero(), zero(), zero(), zero()],
        [zero(), zero(), zero(), zero(), q(1)],
        [[q(0), q(1)], [q(0), q(1)]],
        std::array::from_fn(|_| zero()),
        14,
        1 << 17,
        parts.clone(),
    );
    let second = definition(
        [fr(1, 5), zero(), zero(), zero(), zero()],
        [zero(), q(1), zero(), zero(), zero()],
        [zero(), zero(), q(1), zero(), zero()],
        [[q(0), q(1)], [q(0), q(1)]],
        std::array::from_fn(|_| zero()),
        14,
        1 << 17,
        parts,
    );
    let halves = [fr(1, 8), fr(1, 32), fr(1, 128)];
    let rows = [(first, [fr(1, 4), q(0)]), (second, [fr(1, 8), q(0)])]
        .into_iter()
        .enumerate()
        .map(|(i, (d, centre))| {
            zooms(&d, centre, &halves, 14).into_iter().enumerate().map(|(j, d)| cell(format!("{name}-{i}-{j}"), d)).collect()
        })
        .collect();
    let insets = [(0, 0, 1, true), (0, 1, 2, true), (1, 0, 1, true), (1, 1, 2, true)];
    grid_figure(rows, q(40), q(6), &insets, json!({}))
}

pub fn global_refinement() -> Figure {
    ordinary("global-refinement", stages()[0].clone())
}

pub fn local_refinement() -> Figure {
    ordinary("local-refinement", stages()[1].clone())
}

/// The plane `r = 0` of aligned configurations: `s` and `t`.
fn aligned_plane(parts: Parts) -> Definition {
    let zero = || q(0);
    definition(
        std::array::from_fn(|_| zero()),
        [q(1), zero(), zero(), zero(), zero()],
        [zero(), q(1), zero(), zero(), zero()],
        [[q(0), fr(2, 3)], [q(0), fr(2, 5)]],
        std::array::from_fn(|_| zero()),
        18,
        1 << 17,
        parts,
    )
}

pub fn exotic_directions() -> Figure {
    // In the plane r = 0 only covers eliminate anything, so an unresolved
    // region is the same at every depth: a moderate depth draws it.
    let d = Definition { max_depth: 14, ..aligned_plane(stages()[1].clone()) };
    let (pentagon, h) = ([fr(1708, 10000), fr(2764, 10000)], fr(1, 48));
    let zoom = window(&d, [[&pentagon[0] - &h, &pentagon[0] + &h], [&pentagon[1] - &h, &pentagon[1] + &h]], 14);
    let row = vec![cell("exotic-directions-0".into(), d), cell("exotic-directions-1".into(), zoom)];
    grid_figure(vec![row], q(50), q(5), &[(0, 0, 1, true)], json!({}))
}

pub fn arcs_remaining() -> Figure {
    let parts = with_covers(&[ExoticCover::Square, ExoticCover::Pentagon]);
    let rows = vec![vec![
        cell("arcs-remaining-0".into(), aligned_plane(parts.clone())),
        cell("arcs-remaining-1".into(), with(&arcs_definition(), parts)),
    ]];
    grid_figure(rows, q(50), q(5), &[], json!({}))
}

pub fn arcs_residual() -> Figure {
    use ExoticCover::*;
    let parts = with_covers(&[Square, Pentagon, ArcPlus, ArcMinus]);
    let d = with(&arcs_definition(), parts);
    let endpoint = [fr(236068, 1000000), q(0)];
    let crossing = [fr(1, 2), q(0)];
    let half = fr(1, 32);
    let rows = vec![vec![
        cell("arcs-residual-0".into(), d.clone()),
        cell("arcs-residual-1".into(), window(&d, [[&endpoint[0] - &half, &endpoint[0] + &half], [-half.clone(), half.clone()]], 16)),
        cell("arcs-residual-2".into(), window(&d, [[&crossing[0] - &half, &crossing[0] + &half], [-half.clone(), half.clone()]], 16)),
    ]];
    grid_figure(rows, q(44), q(5), &[(0, 0, 1, false), (0, 0, 2, false)], json!({}))
}
