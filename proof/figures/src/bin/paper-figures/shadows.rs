//! Shadows of the two copies at single configurations, with enlarged insets
//! (and, for some, a support): the Global, Local, square and pentagon, arc,
//! endpoint and crossing examples.
use super::*;

/// Floating-point shadows, only to place inset windows: the orthonormal
/// screen `e₁ = (1, 0, -s)/√L`, `e₂ = (-st, L, -t)/√(LN)` of the library.
fn shadows(c: &[f64; AXES]) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
    let (s, t, r) = (c[0], c[1], [c[2], c[3], c[4]]);
    let l = 1.0 + s * s;
    let n = l + t * t;
    let e1 = [1.0 / l.sqrt(), 0.0, -s / l.sqrt()];
    let d = (l * n).sqrt();
    let e2 = [-s * t / d, l / d, -t / d];
    let rr: f64 = r.iter().map(|x| x * x).sum();
    let k = [[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]];
    let rotation: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| (((1.0 - rr) * if i == j { 1.0 } else { 0.0 }) + 2.0 * r[i] * r[j] + 2.0 * k[i][j]) / (1.0 + rr))
    });
    let project = |w: [f64; 3]| {
        [(0..3).map(|k| w[k] * e1[k]).sum::<f64>(), (0..3).map(|k| w[k] * e2[k]).sum::<f64>()]
    };
    let vertices = vertices_f64();
    let hole = vertices.iter().map(|v| project(*v)).collect();
    let plug = vertices
        .iter()
        .map(|v| project(std::array::from_fn(|i| (0..3).map(|k| rotation[i][k] * v[k]).sum())))
        .collect();
    (hole, plug)
}

/// Indices of the convex hull, counter-clockwise (monotone chain).
fn hull(points: &[[f64; 2]]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&a, &b| points[a].partial_cmp(&points[b]).expect("finite"));
    let cross = |o: usize, a: usize, b: usize| {
        let (o, a, b) = (points[o], points[a], points[b]);
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    };
    let chain = |indices: &mut dyn Iterator<Item = usize>| {
        let mut out: Vec<usize> = Vec::new();
        for i in indices {
            while out.len() >= 2 && cross(out[out.len() - 2], out[out.len() - 1], i) <= 0.0 {
                out.pop();
            }
            out.push(i);
        }
        out.pop();
        out
    };
    let mut lower = chain(&mut order.iter().copied());
    lower.extend(chain(&mut order.iter().rev().copied()));
    lower
}

/// The plug vertex farthest outside the hole's shadow: its index, position
/// and distance outside (negative inside).
fn protrusion(c: &[f64; AXES]) -> (usize, [f64; 2], f64) {
    let (hole, plug) = shadows(c);
    let h = hull(&hole);
    let outside = |p: [f64; 2]| {
        (0..h.len())
            .map(|k| {
                let (a, b) = (hole[h[k]], hole[h[(k + 1) % h.len()]]);
                let normal = [b[1] - a[1], a[0] - b[0]];
                ((p[0] - a[0]) * normal[0] + (p[1] - a[1]) * normal[1]) / normal[0].hypot(normal[1])
            })
            .fold(f64::NEG_INFINITY, f64::max)
    };
    // The first farthest vertex: its antipode protrudes equally far.
    let best = (1..plug.len()).fold(0, |best, i| if outside(plug[i]) > outside(plug[best]) { i } else { best });
    (best, plug[best], outside(plug[best]))
}

fn floats(c: &[QSqrt5; AXES]) -> [f64; AXES] {
    std::array::from_fn(|j| field_f64(&c[j]))
}

/// The hole's shadow is tinted to show its interior; the plug is not. Each
/// body's edges on its far side, which the body itself hides, are dashed.
pub fn scene_styles() -> Value {
    json!({
        "hole-fill": fill("#fcefe0"),
        "hole-edges": line("#e08214", fr(1, 10)),
        "hole-back-edges": dashed("#e08214", fr(1, 10), fr(2, 5)),
        "plug-edges": line("#2166ac", fr(1, 10)),
        "plug-back-edges": dashed("#2166ac", fr(1, 10), fr(2, 5)),
        "inset": line("#404040", fr(1, 6)),
    })
}

pub fn support_styles() -> Value {
    json!({
        "support-line": line("#1a1a1a", fr(1, 7)),
        "support-edge": line("#b35806", fr(1, 2)),
        "support-normal": line("#1b7837", fr(1, 4)),
        "support-vertex": fill("#08306b"),
    })
}

pub fn scene_layers(panel: &str, source: &str) -> Vec<Value> {
    [
        ("hole", "silhouette", "hole-fill"),
        ("hole", "back-wireframe", "hole-back-edges"),
        ("hole", "front-wireframe", "hole-edges"),
        ("plug", "back-wireframe", "plug-back-edges"),
        ("plug", "front-wireframe", "plug-edges"),
    ]
        .iter()
        .map(|(body, geometry, style)| {
            json!({"kind": "scene", "panel": panel, "source": source, "body": body, "geometry": geometry, "style": style})
        })
        .collect()
}

/// A square window of half-width `half` around `centre`, snapped to 2⁻²⁴.
fn window(centre: [f64; 2], half: &Q) -> [[Q; 2]; 2] {
    let (x, y) = (grid(centre[0], 1 << 24), grid(centre[1], 1 << 24));
    [[&x - half, &x + half], [&y - half, &y + half]]
}

/// The largest of 1, 2, 5, 10, 20, 50, … not above `m` (at least 1).
fn round_down(m: f64) -> i64 {
    let mut best = 1;
    let mut unit = 1i64;
    while unit as f64 <= m {
        for lead in [1, 2, 5] {
            if (lead * unit) as f64 <= m {
                best = lead * unit;
            }
        }
        unit *= 10;
    }
    best
}

/// The magnification at which an inset as tall as its overview (which
/// shows the world `[-5, 5]²`) shows about `across` world units above and
/// below its centre.
fn magnification(across: f64) -> i64 {
    round_down(5.0 / across)
}

/// One example: its name, configuration, inset centre and magnification,
/// the length in the screen of what its inset must show (a protrusion, the
/// distance between two edge ends; 0 where they coincide), and optionally
/// a support `(edge, vertex, normal length, dot radius)`.
pub struct Example {
    pub name: String,
    pub configuration: [QSqrt5; AXES],
    pub centre: [f64; 2],
    pub magnification: i64,
    pub feature: f64,
    pub support: Option<([usize; 2], usize, Q, Q)>,
}

/// The smallest printed size of an inset's feature, in millimetres.
const VISIBLE_MM: f64 = 0.8;

/// One example per row, printed `printed_mm` wide: the overview, a square
/// of side `over` showing the world `[-5, 5]²`, and to its right an inset
/// as tall as the overview filling the rest of the width, labelled with its
/// magnification. Every inset's feature must print at least
/// [`VISIBLE_MM`] long (or be zero).
pub fn stack(examples: &[Example], over: i64, gap: i64, printed_mm: f64) -> Value {
    let rows = examples.len() as i64;
    let width = 150;
    let zoom = width - over - 3 * gap;
    let print = printed_mm / width as f64;
    for e in examples {
        // An inset shows `M · over/10` page millimetres per screen unit.
        let printed = e.feature * (e.magnification * over) as f64 / 10.0 * print;
        eprintln!("  {}: ×{}, feature {:.3e}, printed {printed:.2} mm", e.name, e.magnification, e.feature);
        assert!(e.feature == 0.0 || printed >= VISIBLE_MM, "{}: the feature is printed {printed:.2} mm long", e.name);
    }
    let height = q(rows * over + (rows + 1) * gap);
    let mut sources = serde_json::Map::new();
    let mut panels = serde_json::Map::new();
    let mut styles = scene_styles();
    let mut layers = Vec::new();
    if examples.iter().any(|e| e.support.is_some()) {
        styles.as_object_mut().expect("an object").extend(support_styles().as_object().expect("an object").clone());
    }
    for (k, e) in examples.iter().enumerate() {
        let name = e.name.as_str();
        let zoom_name = format!("{name}-zoom");
        let y0 = q(gap + k as i64 * (over + gap));
        let y1 = &y0 + q(over);
        sources.insert(name.into(), configuration(&e.configuration));
        panels.insert(name.into(), panel([q(gap), q(gap + over)], [y0.clone(), y1.clone()], [[q(-5), q(5)], [q(-5), q(5)]]));
        // Vertically `over / 2h = M · over / 10`; horizontally in proportion.
        let half = frac(5, e.magnification);
        let (x, y) = (grid(e.centre[0], 1 << 24), grid(e.centre[1], 1 << 24));
        let across = &half * &frac(zoom, over);
        let window = [[&x - &across, &x + &across], [&y - &half, &y + &half]];
        let zx = q(2 * gap + over);
        panels.insert(zoom_name.clone(), panel([zx.clone(), &zx + q(zoom)], [y0, y1], window));
        layers.extend(scene_layers(name, name));
        layers.extend(scene_layers(&zoom_name, name));
        if let Some((edge, vertex, length, marker)) = &e.support {
            // Three times larger in the overview: legible, yet inside it.
            for (p, l, m) in [(name, length * q(3), marker * q(3)), (zoom_name.as_str(), length.clone(), marker.clone())] {
                layers.push(json!({
                    "kind": "support", "panel": p, "source": name, "edge": edge, "vertex": vertex,
                    "length": rat(&l), "marker": rat(&m),
                    "styles": {"edge": "support-edge", "line": "support-line", "normal": "support-normal", "vertex": "support-vertex"},
                }));
            }
        }
        layers.push(json!({
            "kind": "inset", "overview": name, "detail": zoom_name,
            "connectors": [{"from": "top-right", "to": "top-left"}, {"from": "bottom-right", "to": "bottom-left"}],
            "mark_mm": "6/5", "style": "inset",
        }));
        layers.push(json!({"kind": "label", "panel": zoom_name, "text": format!("×{}", e.magnification)}));
    }
    description(q(width), height, Value::Object(sources), styles, Value::Object(panels), layers)
}

/// A poke configuration, its inset on the farthest protruding plug vertex
/// at a magnification showing about three times the protrusion.
fn poke(name: &str, configuration: [QSqrt5; AXES]) -> Example {
    let (vertex, centre, distance) = protrusion(&floats(&configuration));
    eprintln!("  {name}: plug vertex {vertex} protrudes {distance:.3e}");
    assert!(distance > 0.0, "{name} is not a poke configuration");
    assert!(in_domain(&configuration), "{name} is not in D");
    let magnification = magnification(3.0 * distance);
    Example { name: name.into(), configuration, centre, magnification, feature: distance, support: None }
}

pub fn global() -> Figure {
    let examples: Vec<Example> = [("global-a", "00101000111101011", [1, 10], 11), ("global-b", "01000101010110001111", [0, 9], 10)]
        .into_iter()
        .map(|(name, path, edge, vertex)| {
            let configuration = rational_configuration(&midpoint(path));
            let (_, plug) = shadows(&floats(&configuration));
            let (_, _, feature) = protrusion(&floats(&configuration));
            let support = Some((edge, vertex, fr(3, 20), fr(3, 200)));
            Example { name: name.into(), configuration, centre: plug[vertex], magnification: 10, feature, support }
        })
        .collect();
    Figure::plain(stack(&examples, 65, 4, 150.0))
}

/// `[s, t, r]` with rational coordinates.
fn at(c: [Q; AXES]) -> [QSqrt5; AXES] {
    rational_configuration(&c)
}

pub fn local() -> Figure {
    let examples = [
        poke("local-a", at([fr(1, 4), fr(1, 20), fr(1, 1000), fr(1, 2000), fr(-1, 1000)])),
        poke("local-b", at([fr(2, 5), fr(1, 10), fr(-1, 1000), fr(1, 1000), fr(1, 2000)])),
    ];
    Figure::plain(stack(&examples, 65, 4, 150.0))
}

/// The pentagon view `(η, ξ)` (exact) with the rotation `r`.
pub fn pentagon_configuration(eta: Q, xi: Q, r: [Q; 3]) -> [QSqrt5; AXES] {
    // s = (−5 + 3√5)/10 + (φη − φ̄ξ)/√5 and t = (5 − √5)/10 + (ξ − η)/√5,
    // with φ/√5 = (5 + √5)/10, φ̄/√5 = (−5 + √5)/10 and 1/√5 = √5/5.
    let phi_over_root5 = QSqrt5::new(fr(1, 2), fr(1, 10));
    let conjugate_over_root5 = QSqrt5::new(fr(-1, 2), fr(1, 10));
    let inverse_root5 = QSqrt5::new(q(0), fr(1, 5));
    let s = &(&QSqrt5::new(fr(-1, 2), fr(3, 10)) + &phi_over_root5.scale(&eta)) - &conjugate_over_root5.scale(&xi);
    let t = &QSqrt5::new(fr(1, 2), fr(-1, 10)) + &inverse_root5.scale(&(&xi - &eta));
    let [r1, r2, r3] = r.map(QSqrt5::from_rational);
    [s, t, r1, r2, r3]
}

/// The two pentagon examples, on the two sides of the edge-on line `η = 0`.
pub fn pentagon_examples() -> [[QSqrt5; AXES]; 2] {
    let r = [fr(1, 1000), q(0), q(0)];
    [pentagon_configuration(fr(1, 200), fr(-1, 400), r.clone()), pentagon_configuration(fr(-1, 200), fr(-1, 400), r)]
}

pub fn square_pentagon() -> Figure {
    let half_turn = rational_configuration(&[fr(1, 16), fr(1, 16), fr(-1, 16), fr(1, 16), q(0)]);
    // The plug is the hole turned half about the viewing line, which maps
    // the projected wireframe to itself: the two coincide edge for edge.
    // Show it at the upper right corner of the top edge of the silhouette.
    let (hole, plug) = shadows(&floats(&half_turn));
    let segments = |points: &[[f64; 2]]| -> Vec<[[f64; 2]; 2]> {
        geometry::edges().iter().map(|&[a, b]| [points[a], points[b]]).collect()
    };
    let same = |x: [[f64; 2]; 2], y: [[f64; 2]; 2]| {
        let close = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]) < 1e-12;
        (close(x[0], y[0]) && close(x[1], y[1])) || (close(x[0], y[1]) && close(x[1], y[0]))
    };
    let (hole_edges, plug_edges) = (segments(&hole), segments(&plug));
    assert!(plug_edges.iter().all(|x| hole_edges.iter().any(|y| same(*x, *y))), "the wireframes coincide");
    let corner = hull(&hole).into_iter().map(|i| hole[i]).fold([f64::NEG_INFINITY; 2], |best, p| {
        if p[1] + p[0] / 10.0 > best[1] + best[0] / 10.0 { p } else { best }
    });
    let [pentagon_a, pentagon_b] = pentagon_examples();
    let examples = [
        Example { name: "square-a".into(), configuration: half_turn, centre: corner, magnification: 100, feature: 0.0, support: None },
        poke("square-b", at([fr(1, 64), fr(1, 128), fr(1, 1000), fr(-1, 1000), fr(-1, 1000)])),
        poke("pentagon-a", pentagon_a),
        poke("pentagon-b", pentagon_b),
    ];
    Figure::plain(stack(&examples, 42, 4, 150.0))
}
/// `a₁ = (-5 + 3√5)/10` and `a₃ = (-5 + √5)/10`: `r_A = (a₁, 0, a₃)`.
pub fn arc_constants() -> (QSqrt5, QSqrt5) {
    (QSqrt5::new(fr(-1, 2), fr(3, 10)), QSqrt5::new(fr(-1, 2), fr(1, 10)))
}

/// The configuration `[s, 0, σr_A + θw_σ]` of the arc plane of sign `σ` at
/// arc coordinate `e = (1 - 2s)/(2 + s)`, with `w_σ = (σa₃, 1, -σa₁)`.
pub fn arc_configuration(e: &QSqrt5, theta: &QSqrt5, sigma: i64) -> [QSqrt5; AXES] {
    let (a1, a3) = arc_constants();
    let sigma = QSqrt5::integer(sigma);
    let numerator = &QSqrt5::one() - &(&QSqrt5::integer(2) * e);
    let denominator = &QSqrt5::integer(2) + e;
    let inverse = denominator.conjugate().scale(&(Q::one() / denominator.norm()));
    let s = &numerator * &inverse;
    let st = &sigma * theta;
    [
        s,
        QSqrt5::zero(),
        &(&sigma * &a1) + &(&st * &a3),
        theta.clone(),
        &(&sigma * &a3) - &(&st * &a1),
    ]
}

/// The endpoint's poke neighbour `e = e_* + 1/256`, with `e_* = √5 − 2`.
fn beyond_endpoint() -> [QSqrt5; AXES] {
    arc_configuration(&QSqrt5::new(fr(-511, 256), q(1)), &QSqrt5::zero(), 1)
}

/// The right ends of the plug's and the hole's top edges: the rightmost
/// vertices at the top of the hole's shadow.
fn top_ends(c: &[QSqrt5; AXES]) -> ([f64; 2], [f64; 2]) {
    let (hole, plug) = shadows(&floats(c));
    let top = hole.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
    let end = |points: &[[f64; 2]]| {
        points.iter().filter(|p| p[1] > top - 1e-9).fold([f64::NEG_INFINITY, top], |best, p| if p[0] > best[0] { *p } else { best })
    };
    (end(&plug), end(&hole))
}

/// A configuration of the plane `Π₊`, its inset at the right ends of the
/// two top edges: midway between them, with the distance between them
/// across the middle half of the inset, or, where they coincide, at them
/// with the given magnification.
fn arc_example(name: &str, e: QSqrt5, theta: QSqrt5, coincident: i64) -> Example {
    let configuration = arc_configuration(&e, &theta, 1);
    let (plug, hole) = top_ends(&configuration);
    let distance = (plug[0] - hole[0]).hypot(plug[1] - hole[1]);
    let (magnification, feature) = if distance < 1e-9 { (coincident, 0.0) } else { (magnification(distance), distance) };
    // The top edges a third of the way down: nothing lies above them.
    let half = 5.0 / magnification as f64;
    let centre = [(plug[0] + hole[0]) / 2.0, (plug[1] + hole[1]) / 2.0 - half / 3.0];
    Example { name: name.into(), configuration, centre, magnification, feature, support: None }
}

pub fn arcs() -> Figure {
    let examples = [
        arc_example("arc-a", QSqrt5::from_rational(fr(1, 10)), QSqrt5::zero(), 0),
        arc_example("arc-b", QSqrt5::from_rational(fr(3, 20)), QSqrt5::zero(), 0),
    ];
    Figure::plain(stack(&examples, 65, 4, 150.0))
}

pub fn endpoint_crossing() -> Figure {
    assert!(in_domain(&beyond_endpoint()), "the endpoint's poke neighbour is not in D");
    let (_, _, protruding) = protrusion(&floats(&beyond_endpoint()));
    let mut beyond = arc_example("endpoint-b", QSqrt5::new(fr(-511, 256), q(1)), QSqrt5::zero(), 0);
    // What must be seen is the protrusion of the plug's corner.
    beyond.feature = protruding;
    let examples = [
        arc_example("endpoint-a", QSqrt5::new(q(-2), q(1)), QSqrt5::zero(), beyond.magnification),
        beyond,
        arc_example("crossing-a", QSqrt5::zero(), QSqrt5::zero(), 0),
        arc_example("crossing-b", QSqrt5::from_rational(fr(1, 200)), QSqrt5::from_rational(fr(-1, 200)), 0),
    ];
    Figure::plain(stack(&examples, 42, 4, 150.0))
}

/// The ordered silhouette edge of the hole whose supporting line the plug
/// vertex `vertex` lies farthest beyond (the library's normal is the edge
/// direction turned counterclockwise).
fn support_edge(c: &[f64; AXES], vertex: usize) -> [usize; 2] {
    let (hole, plug) = shadows(c);
    let p = plug[vertex];
    let mut best = ([0, 0], f64::NEG_INFINITY);
    for &[a, b] in geometry::edges() {
        for [a, b] in [[a, b], [b, a]] {
            let d = [hole[b][0] - hole[a][0], hole[b][1] - hole[a][1]];
            let n = [-d[1], d[0]];
            let length = n[0].hypot(n[1]);
            let height = |w: [f64; 2]| (n[0] * w[0] + n[1] * w[1]) / length;
            let extent = hole.iter().map(|w| height(*w)).fold(f64::NEG_INFINITY, f64::max);
            let on_line = extent - height(hole[a]) < 1e-9 && extent - height(hole[b]) < 1e-9;
            if on_line && height(p) - extent > best.1 {
                best = ([a, b], height(p) - extent);
            }
        }
    }
    best.0
}

pub fn pentagon_support() -> Figure {
    let (size, gap) = (q(44), q(4));
    let mut sources = serde_json::Map::new();
    let mut panels = serde_json::Map::new();
    let mut layers = Vec::new();
    let mut styles = scene_styles();
    let extra = styles.as_object_mut().expect("an object");
    extra.extend(support_styles().as_object().expect("an object").clone());
    extra.remove("inset");
    extra.insert("outside".into(), fill("#eeeeee"));
    extra.insert("cover".into(), style(Some("#e6550d"), None, fr(1, 4)));
    extra.insert("edge-on".into(), line("#1a1a1a", fr(1, 5)));
    extra.insert("example".into(), fill("#08306b"));
    extra.insert("frame".into(), line("#404040", fr(1, 6)));
    let frame = |world: &[[Q; 2]; 2], kind: &str| {
        let [[x0, x1], [y0, y1]] = world;
        shapes(kind, vec![polygon(&[[x0.clone(), y0.clone()], [x1.clone(), y0.clone()], [x1.clone(), y1.clone()], [x0.clone(), y1.clone()]])])
    };
    // The pentagon cover's views in the coordinates (η, ξ): the covered
    // rectangle, the views beyond the viewing triangle (ξ > 0), the edge-on
    // line η = 0 and the two examples.
    let h = fr(3, 64);
    let world = [[-h.clone(), h.clone()], [-&h - &fr(1, 64), &h - &fr(1, 64)]];
    panels.insert("views".into(), panel([gap.clone(), &gap + &size], [gap.clone(), &gap + &size], world.clone()));
    let (lo, hi) = (world[0][0].clone(), world[0][1].clone());
    sources.insert("outside".into(), shapes("diagram", vec![polygon(&[[lo.clone(), q(0)], [hi.clone(), q(0)], [hi.clone(), world[1][1].clone()], [lo.clone(), world[1][1].clone()]])]));
    let c = fr(1, 32);
    sources.insert("cover".into(), shapes("diagram", vec![polygon(&[[-c.clone(), -c.clone()], [c.clone(), -c.clone()], [c.clone(), q(0)], [-c.clone(), q(0)]])]));
    sources.insert("edge-on".into(), shapes("diagram", vec![segment([q(0), world[1][0].clone()], [q(0), world[1][1].clone()])]));
    sources.insert("views-frame".into(), frame(&world, "diagram"));
    let dots: Vec<Value> = [fr(1, 200), fr(-1, 200)]
        .iter()
        .map(|eta| json!({"circle": {"centre": point(&[eta.clone(), fr(-1, 400)]), "radius": rat(&fr(1, 1000))}}))
        .collect();
    sources.insert("examples".into(), shapes("diagram", dots));
    for (source, style) in [("outside", "outside"), ("cover", "cover"), ("edge-on", "edge-on"), ("examples", "example"), ("views-frame", "frame")] {
        layers.push(shapes_layer("views", source, style));
    }
    // The two examples' shadows, enlarged at their protruding vertices, with
    // the support beyond which that vertex lies.
    // Left to right as in the diagram: η < 0 first.
    for (k, example) in pentagon_examples().into_iter().rev().enumerate() {
        let name = format!("side-{k}");
        let (vertex, centre, distance) = protrusion(&floats(&example));
        let edge = support_edge(&floats(&example), vertex);
        eprintln!("  {name}: plug vertex {vertex} beyond hole edge {edge:?} by {distance:.3e}");
        let magnification = magnification(3.0 * distance * 40.0 / 34.0);
        let half = frac(5 * 34, magnification * 40);
        // Printed at the canvas width: 44 mm across `2·half`.
        let printed = distance * 44.0 / (2.0 * to_f64(&half));
        eprintln!("  {name}: ×{magnification}, printed {printed:.2} mm");
        assert!(printed >= VISIBLE_MM, "{name}: the protrusion is printed {printed:.2} mm long");
        let x = &gap + &(&(&size + &gap) * &q(k as i64 + 1)) + &q(2);
        let world = window(centre, &half);
        panels.insert(name.clone(), panel([x.clone(), &x + &size], [gap.clone(), &gap + &size], world.clone()));
        sources.insert(format!("{name}-frame"), frame(&world, "screen"));
        sources.insert(name.clone(), configuration(&example));
        layers.extend(scene_layers(&name, &name));
        layers.push(json!({
            "kind": "support", "panel": name, "source": name, "edge": edge, "vertex": vertex,
            "length": rat(&(&half / &q(2))), "marker": rat(&(&half / &q(30))),
            "styles": {"edge": "support-edge", "line": "support-line", "normal": "support-normal", "vertex": "support-vertex"},
        }));
        layers.push(shapes_layer(&name, &format!("{name}-frame"), "frame"));
        layers.push(json!({"kind": "label", "panel": name, "text": format!("×{magnification}")}));
    }
    let width = &(&(&size + &gap) * &q(3)) + &gap + &q(2);
    Figure::plain(description(width, &size + &(&gap * &q(2)), Value::Object(sources), styles, Value::Object(panels), layers))
}
