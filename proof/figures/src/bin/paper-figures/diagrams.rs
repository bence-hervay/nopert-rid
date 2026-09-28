//! Figures drawn from explicit exact geometry: the Bernstein coefficients,
//! the model of a point zoom, and the arc plane near the crossing.
use super::*;
use std::fs;

pub fn bernstein() -> Figure {
    let polynomial = |lo: Q, hi: Q, samples: u32| {
        json!({"kind": "polynomial", "coefficients": ["-1", "3", "-3"], "interval": interval(&lo, &hi), "samples": samples})
    };
    let sources = json!({
        "whole": polynomial(q(0), q(1), 240),
        "left": polynomial(q(0), fr(1, 2), 120),
        "right": polynomial(fr(1, 2), q(1), 120),
    });
    let styles = json!({
        "axis": line("#808080", fr(1, 8)),
        "curve": line("#2166ac", fr(2, 5)),
        "faint": style(None, Some(("#2166ac", fr(1, 5))), fr(1, 4)),
        "polygon": line("#e08214", fr(1, 5)),
        "points": fill("#e08214"),
    });
    let world = [[fr(-1, 10), fr(11, 10)], [fr(-23, 20), fr(13, 20)]];
    let mut panels = serde_json::Map::new();
    let mut layers = Vec::new();
    for (k, (name, source)) in [("p-whole", "whole"), ("p-left", "left"), ("p-right", "right")].into_iter().enumerate() {
        let x0 = q(3 + 46 * k as i64);
        panels.insert(name.into(), panel([x0.clone(), x0 + q(40)], [q(3), q(63)], world.clone()));
        let graph = |source: &str, part: &str, style: &str| {
            json!({"kind": "graph", "panel": name, "source": source, "part": part, "style": style})
        };
        layers.push(graph("whole", "axis", "axis"));
        if source != "whole" {
            layers.push(graph("whole", "curve", "faint"));
        }
        layers.push(graph(source, "control-polygon", "polygon"));
        layers.push(graph(source, "curve", "curve"));
        layers.push(json!({"kind": "control-points", "panel": name, "source": source, "marker": "1/40", "style": "points"}));
    }
    Figure::plain(description(q(138), q(66), sources, styles, Value::Object(panels), layers))
}

/// A sequential blue: white at 0, light blue at 1/2, dark blue at 1
/// (clamped), for values `t = value/2` of the model's quotients.
fn ramp(t: f64) -> String {
    let stops = [(0.0, [255.0, 255.0, 255.0]), (0.5, [158.0, 202.0, 225.0]), (1.0, [8.0, 81.0, 156.0])];
    let t = t.clamp(0.0, 1.0);
    let k = if t <= 0.5 { 0 } else { 1 };
    let (t0, a) = stops[k];
    let (t1, b) = stops[k + 1];
    let u = (t - t0) / (t1 - t0);
    let c: Vec<u8> = (0..3).map(|i| (a[i] + u * (b[i] - a[i])).round() as u8).collect();
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

/// Cells of the model `p = y(x + y)` coloured by the smallest value of
/// the quotient on them: grouped by that value (in eighths).
struct Quotients(BTreeMap<i64, Vec<Value>>);

impl Quotients {
    fn add(&mut self, eighths: i64, corners: &[[Q; 2]]) {
        self.0.entry(eighths).or_default().push(polygon(corners));
    }
}

/// The model of a point zoom: `x` the view distance, `y` the size of the
/// rotation, `p = y(x + y)`. Left, the cells of a face zoom, which divides
/// by `y`: the quotient `x + y` has smallest value 0 on the cell at the
/// special configuration. Right, point zooms with `ρ₀ = 1`: family A
/// (`x = μ`, `y = ρμ`, quotient `1 + ρ`) below the diagonal and family B
/// (`y = δ`, `x = bδ`, quotient `1 + b`) above it.
pub fn toy_zoom() -> Figure {
    let n = 8;
    let at = |i: i64| fr(i, n);
    let mut face = Quotients(BTreeMap::new());
    for i in 0..n {
        for j in 0..n {
            face.add(i + j, &[[at(i), at(j)], [at(i + 1), at(j)], [at(i + 1), at(j + 1)], [at(i), at(j + 1)]]);
        }
    }
    let mut point = Quotients(BTreeMap::new());
    let m = 4;
    let ratio = |j: i64| fr(j, m);
    for i in 0..n {
        for j in 0..m {
            // Family A: μ ∈ [i/n, (i+1)/n], ρ ∈ [j/m, (j+1)/m].
            let (a, b) = (at(i), at(i + 1));
            let (r0, r1) = (ratio(j), ratio(j + 1));
            point.add(8 + 8 * j / m, &[[a.clone(), &r0 * &a], [b.clone(), &r0 * &b], [b.clone(), &r1 * &b], [a.clone(), &r1 * &a]]);
            // Family B: δ ∈ [i/n, (i+1)/n], b ∈ [j/m, (j+1)/m].
            point.add(8 + 8 * j / m, &[[&r0 * &a, a.clone()], [&r1 * &a, a.clone()], [&r1 * &b, b.clone()], [&r0 * &b, b.clone()]]);
        }
    }
    let (size, gap, between) = (q(60), q(4), q(10));
    let world = [[fr(-1, 20), fr(21, 20)], [fr(-1, 20), fr(21, 20)]];
    let mut sources = serde_json::Map::new();
    let mut styles = serde_json::Map::new();
    let mut panels = serde_json::Map::new();
    let mut layers = Vec::new();
    for (k, (name, cells)) in [("face", face), ("point", point)].into_iter().enumerate() {
        let x = &gap + &(&(&size + &between) * &q(k as i64));
        panels.insert(name.into(), panel([x.clone(), &x + &size], [gap.clone(), &gap + &size], world.clone()));
        for (eighths, shapes_of) in cells.0 {
            let source = format!("{name}-{eighths}");
            let style_name = format!("quotient-{eighths}");
            sources.insert(source.clone(), shapes("diagram", shapes_of));
            styles.insert(style_name.clone(), style(Some(&ramp(eighths as f64 / 16.0)), Some(("#ffffff", fr(1, 8))), q(1)));
            layers.push(shapes_layer(name, &source, &style_name));
        }
        // The no-fit set y = 0 and the special configuration.
        let no_fit = format!("{name}-no-fit");
        sources.insert(no_fit.clone(), shapes("diagram", vec![segment([q(0), q(0)], [q(1), q(0)])]));
        layers.push(shapes_layer(name, &no_fit, "no-fit"));
        let special = format!("{name}-special");
        sources.insert(special.clone(), shapes("diagram", vec![json!({"circle": {"centre": ["0", "0"], "radius": "1/60"}})]));
        layers.push(shapes_layer(name, &special, "special"));
    }
    sources.insert("diagonal".into(), shapes("diagram", vec![segment([q(0), q(0)], [q(1), q(1)])]));
    layers.insert(layers.len() - 2, shapes_layer("point", "diagonal", "diagonal"));
    styles.insert("no-fit".into(), line("#1a1a1a", fr(1, 2)));
    styles.insert("special".into(), fill("#1a1a1a"));
    styles.insert("diagonal".into(), line("#e08214", fr(2, 5)));
    let width = &(&size * &q(2)) + &(&gap * &q(2)) + &between;
    Figure::plain(description(width, &size + &(&gap * &q(2)), Value::Object(sources), Value::Object(styles), Value::Object(panels), layers))
}

/// A rational within `2⁻⁶⁰` of `a + b√5`, for drawing only.
fn near(x: &QSqrt5) -> Q {
    use num_bigint::BigInt;
    let scale = BigInt::from(1) << 60u32;
    let root5 = Q::new((BigInt::from(5) * &scale * &scale).sqrt(), scale);
    x.rational_part() + &(x.sqrt5_part() * &root5)
}

/// The views of the covers around aligned configurations: the viewing
/// triangle, the 30 rows of the Local covers (shaded by their radius), the
/// square cover's views `[0, 1/8]²` and the pentagon cover's parallelogram
/// `|η| ≤ 1/32, −1/32 ≤ ξ ≤ 0`.
pub fn local_tiling() -> Figure {
    let rows: [[&str; 5]; 30] = [
        ["1/12", "1/6", "0", "1/10", "1/96"], ["0", "1/6", "1/10", "1/5", "1/80"], ["1/6", "1/3", "0", "1/10", "1/50"],
        ["1/6", "1/4", "1/10", "1/5", "1/50"], ["1/4", "1/3", "1/10", "3/20", "1/50"], ["1/4", "1/3", "3/20", "1/5", "1/50"],
        ["0", "1/12", "1/5", "3/10", "1/50"], ["1/12", "1/8", "1/5", "9/40", "1/50"], ["1/12", "1/8", "9/40", "1/4", "1/50"],
        ["1/8", "1/6", "1/5", "1/4", "1/200"], ["1/12", "1/8", "1/4", "3/10", "1/50"], ["1/8", "7/48", "1/4", "11/40", "1/50"],
        ["7/48", "1/6", "1/4", "21/80", "1/200"], ["7/48", "5/32", "21/80", "11/40", "1/50"], ["5/32", "31/192", "21/80", "43/160", "1/200"],
        ["1/8", "1/6", "11/40", "3/10", "1/50"], ["0", "1/6", "3/10", "2/5", "1/50"], ["1/6", "5/24", "1/5", "9/40", "1/50"],
        ["1/6", "3/16", "9/40", "1/4", "1/200"], ["3/16", "5/24", "9/40", "1/4", "1/50"], ["5/24", "1/4", "1/5", "9/40", "1/50"],
        ["5/24", "1/4", "9/40", "1/4", "1/200"], ["1/6", "3/16", "1/4", "21/80", "1/200"], ["3/16", "5/24", "1/4", "21/80", "1/200"],
        ["5/24", "11/48", "1/4", "11/40", "1/200"], ["1/4", "1/3", "1/5", "1/4", "1/50"], ["1/3", "1/2", "0", "1/10", "1/50"],
        ["1/3", "1/2", "1/10", "1/5", "1/50"], ["1/2", "7/12", "0", "1/10", "1/50"], ["7/12", "5/8", "0", "1/20", "1/20"],
    ];
    let r = |text: &str| rid::arithmetic::exact::parse_rational(text).expect("a rational");
    let mut by_radius: BTreeMap<String, Vec<Value>> = BTreeMap::new();
    for row in rows {
        let [s0, s1, t0, t1] = [r(row[0]), r(row[1]), r(row[2]), r(row[3])];
        by_radius.entry(row[4].to_owned()).or_default().push(polygon(&[[s0.clone(), t0.clone()], [s1.clone(), t0.clone()], [s1, t1.clone()], [s0, t1]]));
    }
    let mut sources = serde_json::Map::new();
    let mut styles = serde_json::Map::new();
    let mut order: Vec<(String, String)> = Vec::new();
    let mut draw = |source: &str, style: &str| order.push((source.to_owned(), style.to_owned()));
    // The viewing triangle s, t ≥ 0, φs + φ²t ≤ 1: corners (0, 0), (1/φ, 0), (0, 1/φ²).
    let inverse_phi = near(&QSqrt5::new(fr(-1, 2), fr(1, 2)));
    let inverse_phi2 = near(&QSqrt5::new(fr(3, 2), fr(-1, 2)));
    sources.insert("triangle".into(), shapes("diagram", vec![polygon(&[[q(0), q(0)], [inverse_phi, q(0)], [q(0), inverse_phi2]])]));
    styles.insert("triangle".into(), fill("#f7f7f7"));
    draw("triangle", "triangle");
    // The square and pentagon covers beneath the Local rows, and outlined
    // on top: the rows overlap them.
    let square = vec![polygon(&[[q(0), q(0)], [fr(1, 8), q(0)], [fr(1, 8), fr(1, 8)], [q(0), fr(1, 8)]])];
    let view = |eta: Q, xi: Q| {
        let c = super::shadows::pentagon_configuration(eta, xi, [q(0), q(0), q(0)]);
        [near(&c[0]), near(&c[1])]
    };
    let c = fr(1, 32);
    let pentagon = vec![polygon(&[view(-c.clone(), -c.clone()), view(c.clone(), -c.clone()), view(c.clone(), q(0)), view(-c.clone(), q(0))])];
    for (name, color, shapes_of) in [("square", "#fdae6b", square), ("pentagon", "#e6550d", pentagon)] {
        sources.insert(name.into(), shapes("diagram", shapes_of));
        // A light tint: each is one region with no grid of its own.
        styles.insert(name.into(), style(Some(color), None, fr(2, 5)));
        styles.insert(format!("{name}-edge"), line(color, fr(1, 4)));
        draw(name, name);
    }
    // Lighter green for larger radii: 1/20 and 1/50 light, 1/200 dark.
    for (radius, color) in [("1/20", "#c7e9c0"), ("1/50", "#a1d99b"), ("1/80", "#74c476"), ("1/96", "#74c476"), ("1/200", "#31a354")] {
        let shapes_of = by_radius.remove(radius).expect("a radius of the table");
        let name = format!("local-{}", radius.replace('/', "-"));
        sources.insert(name.clone(), shapes("diagram", shapes_of));
        styles.insert(name.clone(), style(Some(color), Some(("#ffffff", fr(1, 6))), q(1)));
        draw(&name, &name);
    }
    assert!(by_radius.is_empty(), "every radius has a colour");
    for name in ["square", "pentagon"] {
        draw(name, &format!("{name}-edge"));
    }
    styles.insert("triangle-edge".into(), line("#404040", fr(1, 5)));
    draw("triangle", "triangle-edge");
    // The whole triangle, and the pentagon's neighbourhood enlarged.
    let pentagon_view = [fr(1708, 10000), fr(2764, 10000)];
    let h = fr(1, 32);
    let zoom_world = [[&pentagon_view[0] - &h, &pentagon_view[0] + &h], [&pentagon_view[1] - &h, &pentagon_view[1] + &h]];
    let panels = json!({
        "views": panel([q(4), q(88)], [q(4), fr(272, 5)], [[q(0), fr(2, 3)], [q(0), fr(2, 5)]]),
        "pentagon-views": panel([q(94), q(144)], [fr(11, 5), fr(261, 5)], zoom_world),
    });
    let mut layers = Vec::new();
    for panel_name in ["views", "pentagon-views"] {
        layers.extend(order.iter().map(|(source, style)| shapes_layer(panel_name, source, style)));
    }
    layers.push(json!({
        "kind": "inset", "overview": "views", "detail": "pentagon-views",
        "connectors": [{"from": "top-right", "to": "top-left"}, {"from": "bottom-right", "to": "bottom-left"}],
        "mark_mm": "6/5", "style": "inset",
    }));
    styles.insert("inset".into(), line("#404040", fr(1, 6)));
    Figure::plain(description(q(148), fr(292, 5), Value::Object(sources), Value::Object(styles), panels, layers))
}

/// The leaves of a zoom's tree (preorder: a digit splits that axis at its
/// midpoint, `.` is a leaf) over the root box.
fn leaves(tree: &str, root: Vec<[Q; 2]>) -> Vec<Vec<[Q; 2]>> {
    fn walk(tree: &mut std::str::Chars, cell: Vec<[Q; 2]>, out: &mut Vec<Vec<[Q; 2]>>) {
        match tree.next().expect("a complete tree") {
            '.' => out.push(cell),
            digit => {
                let axis = digit.to_digit(10).expect("an axis") as usize;
                let middle = (&cell[axis][0] + &cell[axis][1]) / q(2);
                let mut lower = cell.clone();
                lower[axis][1] = middle.clone();
                let mut upper = cell;
                upper[axis][0] = middle;
                walk(tree, lower, out);
                walk(tree, upper, out);
            }
        }
    }
    let mut out = Vec::new();
    walk(&mut tree.chars(), root, &mut out);
    out
}

/// The crossing cover of the plane `Π₊` where it meets the plane itself
/// (`v = ζ₁ = ζ₃ = 0`), in the coordinates `(e, θ)`: the cells of family A
/// (`e = ±μ`, `θ = ±ρμ`), of family B (`θ = ±δ`, `e = bδ`) and of the
/// sheared family around the second arc (`e = μ`, `θ + e = ±ρ'μ`), the cells
/// of family A handed over to it outlined, and the two arcs.
pub fn crossing_plane() -> Figure {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../data/exotic/crossing+.json");
    let file: Value = serde_json::from_slice(&fs::read(&path).expect("the crossing cover")).expect("JSON");
    let point = &file["parameters"]["shape"]["point"];
    let r = |v: &Value| rid::arithmetic::exact::parse_rational(v.as_str().expect("a rational")).expect("a rational");
    let (radii, radius, ratio) = ([r(&point["radii"][0]), r(&point["radii"][1])], r(&point["radius"]), r(&point["ratio"]));
    let window_radius = r(&point["window"]["radius"]);
    // Offset axes: 4v/5 in [0, 1], θ, ζ₁, ζ₃ in [-1, 1]; the plane is the
    // zero of every offset axis except θ, which is offset axis 1.
    let unit = |axis: usize| if axis == 0 { [q(0), q(1)] } else { [q(-1), q(1)] };
    let mut groups: BTreeMap<&str, Vec<Value>> = BTreeMap::new();
    // The zoom through which the second arc runs, in its own coordinates
    // (20μ, ρ): its cells in the plane, and those handed over.
    let mut zoom_cells: BTreeMap<&str, Vec<Value>> = BTreeMap::new();
    for zoom in file["zooms"].as_array().expect("zooms") {
        let name = zoom["name"].as_str().expect("a name");
        let parts: Vec<&str> = name.split(' ').collect();
        let offset = parts.last().expect("an offset face");
        let (axis, sign) = (offset[1..2].parse::<usize>().expect("an axis"), if offset.ends_with('+') { 1 } else { -1 });
        if axis != 1 {
            continue; // Meets the plane only on the first arc.
        }
        let free: Vec<[Q; 2]> = (0..4).filter(|&a| a != axis).map(unit).collect();
        let (root, family) = match parts[0] {
            "A" | "A'" => {
                let g = if parts[1].ends_with('+') { 1 } else { 0 };
                let limit = if parts[0] == "A" { ratio.clone() } else { window_radius.clone() };
                let mut root = vec![[q(0), radii[g].clone()], [q(0), limit]];
                root.extend(free);
                (root, parts[0])
            }
            _ => {
                let mut root = vec![[q(0), radius.clone()], [-(&q(1) / &ratio), &q(1) / &ratio]];
                root.extend(free);
                (root, "B")
            }
        };
        let base_sign = if parts[1].ends_with('-') { -1 } else { 1 };
        let cells = leaves(zoom["tree"].as_str().expect("a tree"), root);
        let marks = zoom["leaves"].as_array().expect("leaves");
        assert_eq!(cells.len(), marks.len(), "{name}: one mark per leaf");
        for (cell, mark) in cells.iter().zip(marks) {
            if !cell[2..].iter().all(|[lo, hi]| *lo <= q(0) && q(0) <= *hi) {
                continue; // Misses the plane.
            }
            let corners: Vec<[Q; 2]> = [(0, 0), (1, 0), (1, 1), (0, 1)]
                .iter()
                .map(|&(i, j)| {
                    let (d, x) = (&cell[0][i], &cell[1][j]);
                    let (sd, sx) = (q(sign), q(base_sign));
                    match family {
                        "A" => [&sx * d, &(&sd * x) * d],
                        "A'" => [d.clone(), &(&(&sd * x) * d) - d],
                        _ => [x * d, &sd * d],
                    }
                })
                .collect();
            let group = if mark.as_str() == Some("delegated") { "delegated" } else { family };
            groups.entry(group).or_default().push(polygon(&corners));
            if name == "A b0+ o1-" {
                let [m0, m1] = [&cell[0][0] * &q(20), &cell[0][1] * &q(20)];
                let [r0, r1] = cell[1].clone();
                let rectangle = polygon(&[[m0.clone(), r0.clone()], [m1.clone(), r0], [m1, r1.clone()], [m0, r1]]);
                zoom_cells.entry(if group == "delegated" { "delegated" } else { "A" }).or_default().push(rectangle);
            }
        }
    }
    let colours = [("A", "#9ecae1"), ("B", "#fdd49e"), ("A'", "#a1d99b")];
    let mut sources = serde_json::Map::new();
    let mut styles = serde_json::Map::new();
    let mut order = Vec::new();
    for (family, color) in colours {
        let name = format!("family-{}", if family == "A'" { "sheared" } else { family }).to_lowercase();
        sources.insert(name.clone(), shapes("diagram", groups.remove(family).unwrap_or_default()));
        styles.insert(name.clone(), style(Some(color), Some(("#ffffff", fr(1, 12))), q(1)));
        order.push((name.clone(), name));
    }
    sources.insert("delegated".into(), shapes("diagram", groups.remove("delegated").unwrap_or_default()));
    styles.insert("delegated".into(), line("#c51b8a", fr(1, 4)));
    order.push(("delegated".into(), "delegated".into()));
    let reach = &radii[1] + &fr(1, 64);
    sources.insert("first-arc".into(), shapes("diagram", vec![segment([q(0), q(0)], [reach.clone(), q(0)])]));
    sources.insert("second-arc".into(), shapes("diagram", vec![segment([q(0), q(0)], [reach.clone(), -reach.clone()])]));
    styles.insert("arc".into(), line("#1a1a1a", fr(2, 5)));
    order.push(("first-arc".into(), "arc".into()));
    order.push(("second-arc".into(), "arc".into()));
    let h = &radii[1] + &fr(1, 128);
    sources.insert("zoom-a".into(), shapes("diagram", zoom_cells.remove("A").unwrap_or_default()));
    sources.insert("zoom-delegated".into(), shapes("diagram", zoom_cells.remove("delegated").unwrap_or_default()));
    sources.insert("zoom-arcs".into(), shapes("diagram", vec![
        segment([q(0), q(0)], [&radii[1] * &q(20), q(0)]),
        segment([q(0), q(1)], [&radii[1] * &q(20), q(1)]),
    ]));
    let panels = json!({
        "plane": panel([q(4), q(68)], [q(4), q(68)], [[-h.clone(), h.clone()], [-h.clone(), h.clone()]]),
        "zoom": panel([q(80), q(144)], [q(4), q(68)], [[fr(-1, 16), &h * &q(20)], [fr(-1, 16), &h * &q(20)]]),
    });
    let mut layers: Vec<Value> = order.iter().map(|(source, style)| shapes_layer("plane", source, style)).collect();
    for (source, style) in [("zoom-a", "family-a"), ("zoom-delegated", "family-sheared"), ("zoom-delegated", "delegated"), ("zoom-arcs", "arc")] {
        layers.push(shapes_layer("zoom", source, style));
    }
    Figure::plain(description(q(148), q(72), Value::Object(sources), Value::Object(styles), panels, layers))
}
