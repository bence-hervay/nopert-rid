//! The figures: passages of other solids, the margin, and the
//! rhombicosidodecahedron's parameters.
use super::geometry::*;
use super::optimise::{minimise, Random};
use super::passage::*;
use super::svg::{Panel, Style, Svg};
use std::f64::consts::PI;

/// A panel in `rect` showing all of `points` with a relative margin.
fn fit(svg: &mut Svg, rect: [f64; 4], points: &[V2], margin: f64, clip: bool) -> Panel {
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in points {
        for k in 0..2 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let centre = [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0];
    let scale = (rect[2] / (hi[0] - lo[0])).min(rect[3] / (hi[1] - lo[1])) * (1.0 - 2.0 * margin);
    svg.panel(rect, centre, scale, clip)
}

fn scene_points(passage: &Passage, camera: &Camera) -> Vec<V2> {
    passage.hole.vertices.iter().chain(&passage.plug.vertices).map(|&v| camera.project(v)).collect()
}

/// Where an enlargement at a tight place looks: between the plug's vertex
/// and the nearest point of the hole's shadow.
fn focus(passage: &Passage, p: V2) -> V2 {
    let hole = passage.hole_shadow();
    let (n, h) = sides(&hole)
        .into_iter()
        .min_by(|a, b| (a.1 - a.0[0] * p[0] - a.0[1] * p[1]).total_cmp(&(b.1 - b.0[0] * p[0] - b.0[1] * p[1])))
        .unwrap();
    let gap = h - n[0] * p[0] - n[1] * p[1];
    [p[0] + n[0] * gap / 2.0, p[1] + n[1] * gap / 2.0]
}

/// The magnification that shows a gap `gap` as about a fifth of a detail
/// panel `side` millimetres wide, or `fallback` for a gap of zero.
fn magnification(overview: &Panel, side: f64, gap: f64, fallback: f64) -> f64 {
    if gap <= 0.0 {
        return fallback;
    }
    round_down(side / 5.0 / (gap * overview.scale))
}

/// Extra drawing on the three-dimensional view.
type Decoration<'a> = &'a dyn Fn(&mut Svg, &Panel, &Camera);

/// A three-part figure: the solids, the shadows, and enlargements at the
/// `places` stacked on the right.
fn three_parts(passage: &Passage, look: &Look, places: &[(V2, f64)], fallback: f64) -> String {
    three_parts_with(passage, look, places, fallback, &|_, _, _| {}, &[])
}

fn three_parts_with(
    passage: &Passage,
    look: &Look,
    places: &[(V2, f64)],
    fallback: f64,
    decorate: Decoration,
    reach_also: &[V3],
) -> String {
    let (h, gap, w2) = (48.0, 5.0, 48.0);
    let mut seen = scene_points(passage, &look.camera);
    seen.extend(reach_also.iter().map(|&v| look.camera.project(v)));
    // The solids' panel is as wide as their picture needs, at most 52 mm.
    let span = |k: usize| seen.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max) - seen.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min);
    let w1 = (h * span(0) / span(1) * 1.04).min(52.0);
    let n = places.len() as f64;
    let side = ((h - (n - 1.0) * 3.0) / n).min(34.0);
    let width = w1 + w2 + side + 2.0 * gap + 2.0;
    let mut svg = Svg::new(width, h + 2.0);
    let solids = fit(&mut svg, [1.0, 1.0, w1, h], &seen, 0.02, false);
    draw_solids(&mut svg, &solids, passage, look);
    decorate(&mut svg, &solids, &look.camera);
    let shadows = fit(&mut svg, [1.0 + w1 + gap, 1.0, w2, h], &passage.hole_shadow(), 0.03, false);
    draw_shadows(&mut svg, &shadows, passage);
    for (k, &(p, gap_size)) in places.iter().enumerate() {
        let top = 1.0 + (h - n * side - (n - 1.0) * 3.0) / 2.0;
        let rect = [1.0 + w1 + w2 + 2.0 * gap, top + k as f64 * (side + 3.0), side, side];
        let m = magnification(&shadows, side, gap_size.abs(), fallback);
        let detail = enlarge(&mut svg, &shadows, rect, focus(passage, p), m, 1.2);
        draw_shadows(&mut svg, &detail, passage);
        label_corner(&mut svg, rect, m);
    }
    svg.finish()
}

/// The screen axes and direction of the projection `M_{θ,φ}` of
/// Fredriksson, whose rows are the screen axes.
fn fredriksson_frame(theta: f64, phi: f64) -> [V3; 3] {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    let e1 = [-st, ct, 0.0];
    let e2 = [-ct * cp, -st * cp, sp];
    [e1, e2, cross(e1, e2)]
}

/// The passage of Fredriksson's parameters `[θp, φp, θq, φq, α, u, v]`: the
/// hole seen along `(θq, φq)`, the plug along `(θp, φp)`, turned by `α` in
/// the screen and moved by `(u, v)`.
fn fredriksson(solid: &Polyhedron, x: &[f64]) -> Passage {
    let p = fredriksson_frame(x[0], x[1]);
    let q = fredriksson_frame(x[2], x[3]);
    let (sa, ca) = x[4].sin_cos();
    let (u, v) = (x.get(5).copied().unwrap_or(0.0), x.get(6).copied().unwrap_or(0.0));
    // In p-coordinates, turn the screen by α; then read them as q-coordinates.
    let turn = |c: V3| [ca * c[0] - sa * c[1], sa * c[0] + ca * c[1], c[2]];
    let plug = solid.map(|w| {
        let c = turn([dot(p[0], w), dot(p[1], w), dot(p[2], w)]);
        add(add(scale(c[0] + u, q[0]), scale(c[1] + v, q[1])), scale(c[2], q[2]))
    });
    Passage { hole: solid.map(|w| w), plug, frame: q }
}

fn containment(passage: &Passage) -> f64 {
    containment_scale(&passage.hole_shadow(), &passage.plug_shadow())
}

pub fn cube(half: f64) -> Polyhedron {
    Polyhedron::new(cyclic_signed(half, half, half), 1e-9)
}

/// Nieuwland's passage: the unit cube, and a cube of edge `3√2/4` through
/// it along `(2, 2, -1)`; its square cross-section has its corners on the
/// unit cube's edges at distances `1/4` and `3/4` from the vertices.
fn nieuwland_passage() -> Passage {
    let d = unit([2.0, 2.0, -1.0]);
    let square = [[0.25, -0.5, -0.5], [-0.5, 0.25, -0.5], [-0.25, 0.5, 0.5], [0.5, -0.25, 0.5]];
    let edge = 3.0 * 2f64.sqrt() / 4.0;
    let mut plug = Vec::new();
    for sign in [-0.5, 0.5] {
        plug.extend(square.iter().map(|&q| add(q, scale(sign * edge, d))));
    }
    let e1 = unit(sub(square[1], square[0]));
    Passage { hole: cube(0.5), plug: Polyhedron::new(plug, 1e-9), frame: [e1, cross(d, e1), d] }
}

pub fn nieuwland_cube() -> String {
    let passage = nieuwland_passage();
    let scale_needed = containment(&passage);
    assert!((scale_needed - 1.0).abs() < 1e-12, "the square touches: {scale_needed}");
    let places = passage.tight_places(1, 0.1);
    let look = Look { camera: exit_view(&passage.frame, 0.6, -2.2), tunnel: true, hole_fill: 0.07, opaque: false };
    let r = reach(&passage);
    three_parts(&passage.advance(1.6 * r), &look, &places, 20.0)
}

/// The five Platonic solids with circumradius about one.
fn platonic() -> Vec<(&'static str, Polyhedron)> {
    let p = PHI;
    let tetra = vec![[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];
    let mut dodeca = cyclic_signed(1.0, 1.0, 1.0);
    dodeca.extend(cyclic_signed(0.0, 1.0 / p, p));
    let normalise = |v: Vec<V3>| {
        let r = dot(v[0], v[0]).sqrt();
        Polyhedron::new(v.into_iter().map(|w| scale(1.0 / r, w)).collect(), 1e-9)
    };
    vec![
        ("tetrahedron", normalise(tetra)),
        ("cube", normalise(cyclic_signed(1.0, 1.0, 1.0))),
        ("octahedron", normalise(cyclic_signed(1.0, 0.0, 0.0))),
        ("dodecahedron", normalise(dodeca)),
        ("icosahedron", normalise(cyclic_signed(0.0, 1.0, p))),
    ]
}

/// The best passage found for `solid`: least containment scale over the
/// viewing directions, the screen turn and (if `translate`) a shift.
fn best_passage(solid: &Polyhedron, translate: bool, seed: u64) -> (Vec<f64>, f64) {
    let f = |x: &[f64]| containment(&fredriksson(solid, x));
    minimise(&f, 300, seed, |r: &mut Random| {
        let mut x = vec![2.0 * PI * r.next(), PI * r.next(), 2.0 * PI * r.next(), PI * r.next(), 2.0 * PI * r.next()];
        if translate {
            x.extend([0.0, 0.0]);
        }
        x
    })
}

pub fn platonic_passages() -> Vec<(&'static str, Passage, f64)> {
    platonic()
        .into_iter()
        .enumerate()
        .map(|(k, (name, solid))| {
            let (x, s) = best_passage(&solid, name == "tetrahedron", 11 + k as u64);
            (name, fredriksson(&solid, &x), s)
        })
        .collect()
}

pub fn platonic_figure() -> String {
    let passages = platonic_passages();
    let (cell, gap, row_gap) = (27.0, 3.0, 4.0);
    let width = 5.0 * cell + 4.0 * gap + 2.0;
    let mut svg = Svg::new(width, 3.0 * cell + 2.0 * row_gap + 2.0);
    for (k, (name, passage, s)) in passages.iter().enumerate() {
        eprintln!("{name}: containment scale {s:.6}, Nieuwland number at least {:.6}", 1.0 / s);
        let x = 1.0 + k as f64 * (cell + gap);
        let passage_view = Passage {
            hole: passage.hole.map(|v| v),
            plug: passage.plug.map(|v| v),
            frame: passage.frame,
        }
        .advance(1.6 * reach(passage));
        let look = Look { camera: exit_view(&passage.frame, 0.6, -2.2), tunnel: true, hole_fill: 0.07, opaque: false };
        let solids = fit(&mut svg, [x, 1.0, cell, cell], &scene_points(&passage_view, &look.camera), 0.02, false);
        draw_solids(&mut svg, &solids, &passage_view, &look);
        let shadows = fit(&mut svg, [x, 1.0 + cell + row_gap, cell, cell], &passage.hole_shadow(), 0.03, false);
        draw_shadows(&mut svg, &shadows, passage);
        let (p, g) = passage.tight_places(1, 0.1)[0];
        let rect = [x, 1.0 + 2.0 * (cell + row_gap), cell, cell];
        let m = magnification(&shadows, cell, g, 10.0);
        let detail = enlarge(&mut svg, &shadows, rect, focus(passage, p), m, 1.2);
        draw_shadows(&mut svg, &detail, passage);
        label_corner(&mut svg, rect, m);
    }
    svg.finish()
}

/// A camera facing the exit of the tunnel, tilted `tilt` radians from the
/// direction of motion towards the screen direction at angle `towards`:
/// with no tilt it sees exactly the shadows' picture.
fn exit_view(frame: &[V3; 3], tilt: f64, towards: f64) -> Camera {
    let [e1, e2, d] = *frame;
    let side = add(scale(towards.cos(), e1), scale(towards.sin(), e2));
    Camera::facing(add(scale(tilt.cos(), d), scale(tilt.sin(), side)), e2)
}

/// How far the solids reach along the direction of motion.
fn reach(passage: &Passage) -> f64 {
    passage.hole.vertices.iter().map(|&v| dot(v, passage.frame[2]).abs()).fold(0.0, f64::max)
}

/// The vertices of Tom 7's candidate #214 (Renshaw's Lean file).
const NOPERT214: [V3; 20] = [
    [0.5542570167628148, 0.13498214234883502, 0.5670539264866502],
    [0.839503072954794, 0.4526456900329921, 0.3005768949769993],
    [0.7619849874129984, 0.5429603653859462, -0.04074876955046458],
    [0.591853727475924, 0.13110932665305655, -0.7653829489805983],
    [0.04289919136693016, 0.5688415234175154, 0.5670539264866502],
    [-0.17107091670577093, 0.9382900786342319, 0.3005768949769993],
    [-0.2809196830211049, 0.8924747677745009, -0.04074876955046458],
    [0.05820048051375984, 0.6034013542664042, -0.7653829489805983],
    [-0.5277438584081657, 0.2165812533354586, 0.5670539264866502],
    [-0.9452307139655625, 0.1272494698697746, 0.3005768949769993],
    [-0.9356028996288878, 0.008619375200364626, -0.04074876955046458],
    [-0.5558838523568445, 0.2418132191412975, -0.7653829489805983],
    [-0.36906283321718863, -0.4349869475301504, 0.5670539264866502],
    [-0.41311379173527674, -0.8596455812043055, 0.3005768949769993],
    [-0.2973147089225043, -0.8871477009388876, -0.04074876955046458],
    [-0.40175559506751807, -0.45395256590805566, -0.7653829489805983],
    [0.29965048349560947, -0.4854179715716587, 0.5670539264866502],
    [0.6899123494518161, -0.6585396573326932, 0.3005768949769993],
    [0.7518523041594987, -0.5569068074219243, -0.04074876955046458],
    [0.3075852394346786, -0.5223713341527026, -0.7653829489805983],
];

/// Renshaw's passage for #214: the outer rotation turns the hole, the inner
/// one the plug, which is then shifted; both are seen along `z`.
pub fn nopert214_passage() -> Passage {
    let solid = Polyhedron::new(NOPERT214.to_vec(), 1e-6);
    let outer = quaternion([0.513670889522, -0.357385492934, -0.445478171674, -0.640286674279]);
    let inner = quaternion([0.513409524171, -0.357274843179, -0.445838043814, -0.640307571101]);
    let offset = [-0.000068618674499, 0.00004658859372, 0.0];
    Passage {
        hole: solid.map(|v| apply(&outer, v)),
        plug: solid.map(|v| add(apply(&inner, v), offset)),
        frame: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    }
}

pub fn nopert214() -> String {
    let passage = nopert214_passage();
    let places = passage.tight_places(2, 0.2);
    eprintln!("#214 clearances: {:?}", places.iter().map(|p| p.1).collect::<Vec<_>>());
    let look = Look { camera: exit_view(&passage.frame, 0.6, -2.2), tunnel: true, hole_fill: 0.07, opaque: false };
    let r = reach(&passage);
    three_parts(&passage.advance(1.6 * r), &look, &places, 10.0)
}

/// The triakis tetrahedron (coordinates of McCooey), `parity` choosing
/// which tetrahedron carries the six-valent vertices.
pub fn triakis(parity: bool) -> Polyhedron {
    let (c0, c1) = (3.0 * 2f64.sqrt() / 8.0, 5.0 * 2f64.sqrt() / 8.0);
    let even = [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];
    let (small, large) = if parity { (c0, c1) } else { (c1, c0) };
    let mut v: Vec<V3> = even.iter().map(|&e| scale(small, e)).collect();
    v.extend(even.iter().map(|&e| scale(-large, e)));
    Polyhedron::new(v, 1e-9)
}

/// Fredriksson's parameters for the triakis tetrahedron.
pub const TRIAKIS: [f64; 7] = [6.283130, 0.817234, 1.548107, 2.356150, 6.2671031, 0.0001408, -0.0000022];

/// Fredriksson's passage, polished by a local minimisation: as printed,
/// the parameters overshoot by about `2·10⁻⁵` in this reading, while a
/// passage with clearance `4·10⁻⁶` lies within `2·10⁻⁴` of them.
pub fn triakis_passage() -> Passage {
    fredriksson(&triakis(false), &triakis_parameters())
}

fn triakis_parameters() -> Vec<f64> {
    let solid = triakis(false);
    let g = |x: &[f64]| containment(&fredriksson(&solid, x));
    let (x, _) = super::optimise::nelder_mead(&g, &TRIAKIS, 1e-5, 4000);
    super::optimise::nelder_mead(&g, &x, 1e-6, 4000).0
}

pub fn triakis_figure() -> String {
    let x = triakis_parameters();
    let change = x.iter().zip(TRIAKIS).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    eprintln!("triakis parameters {x:?}, at most {change:.1e} from the printed ones");
    eprintln!("printed parameters: containment scale {:.7}", containment(&fredriksson(&triakis(false), &TRIAKIS)));
    let passage = triakis_passage();
    eprintln!("polished: containment scale {:.9}", containment(&passage));
    let places = passage.tight_places(2, 0.2);
    eprintln!("triakis clearances: {:?}", places.iter().map(|p| p.1).collect::<Vec<_>>());
    let look = Look { camera: exit_view(&passage.frame, 0.6, -2.2), tunnel: true, hole_fill: 0.07, opaque: false };
    let r = reach(&passage);
    three_parts(&passage.advance(1.6 * r), &look, &places, 10.0)
}

const AXIS: &str = "#1b7837";

/// An arrow in space, seen by `camera`.
fn arrow3(svg: &mut Svg, panel: &Panel, camera: &Camera, from: V3, to: V3, color: &'static str, width: f64) {
    svg.arrow(panel, camera.project(from), camera.project(to), color, width, 1.6);
}

/// A configuration of the reduced domain: `u = (s, t, 1)` and the rotation
/// vector `r`.
pub const VIEW: [f64; 2] = [0.2, 0.1];
pub const ROTATION: V3 = [0.25, -1.0 / 12.0, 5.0 / 24.0];

/// The camera of the rhombicosidodecahedron's figures: `z` upwards, seen
/// from the side of negative `y`, a little from above.
fn rid_camera() -> Camera {
    Camera::new(-1.17, 0.35)
}

/// The hole, and the plug turned by `ROTATION` and moved along `u`.
fn rid_passage() -> Passage {
    let d = unit([VIEW[0], VIEW[1], 1.0]);
    let e2 = unit(sub([0.0, 1.0, 0.0], scale(d[1], d)));
    let turn = rotation(ROTATION);
    let solid = rid();
    let plug = solid.map(|v| apply(&turn, v));
    Passage { hole: solid, plug, frame: [cross(e2, d), e2, d] }
}

pub fn rid_parametrisation() -> String {
    let passage = rid_passage();
    let places = passage.tight_places(1, 1.0);
    eprintln!("RID: plug vertex furthest outside by {:.4}", -places[0].1);
    let r = reach(&passage);
    let shift = 2.0 * r;
    let passage = passage.advance(shift);
    let (axis, d) = (unit(ROTATION), passage.frame[2]);
    let camera = rid_camera();
    let look = Look { camera, tunnel: false, hole_fill: 0.0, opaque: true };
    let centre = scale(shift, d);
    let angle = 2.0 * dot(ROTATION, ROTATION).sqrt().atan();
    eprintln!("RID: rotation angle {:.2}°", angle.to_degrees());
    // The circle showing the angle sits on the axis just outside the plug.
    let plug = &passage.plug;
    let leave = (0..plug.faces.len())
        .map(|f| (plug.normal(f), plug.vertices[plug.faces[f][0]]))
        .filter(|(n, _)| dot(*n, axis) > 1e-9)
        .map(|(n, v)| dot(n, sub(v, centre)) / dot(n, axis))
        .fold(f64::INFINITY, f64::min);
    let hub = add(centre, scale(leave + 1.0 * r, axis));
    let tip = add(hub, scale(0.8 * r, axis));
    let radius = 0.55 * r;
    // The circle is seen as an ellipse; the angle is drawn half-way between
    // the top of the ellipse and its end away from the plug, where neither
    // its sides nor its arc are foreshortened much.
    let across = unit(cross(axis, camera.toward));
    let a = if dot(across, camera.up) > 0.0 { across } else { scale(-1.0, across) };
    let b = cross(axis, a);
    let away = unit(sub(scale(dot(camera.toward, axis), axis), camera.toward));
    let middle = add(a, away);
    let centre_angle = dot(middle, b).atan2(dot(middle, a));
    let around = move |phi: f64| add(hub, add(scale(radius * phi.cos(), a), scale(radius * phi.sin(), b)));
    let circle: Vec<V3> = (0..=120).map(|k| around(2.0 * PI * k as f64 / 120.0)).collect();
    let (from, to) = (scale(-1.2 * r, d), scale(shift + 1.45 * r, d));
    let mut ends = circle.clone();
    ends.extend([from, to, tip]);
    let bodies = (passage.hole.map(|v| v), passage.plug.map(|v| v));
    let decorate = move |svg: &mut Svg, panel: &Panel, camera: &Camera| {
        let solids = [&bodies.0, &bodies.1];
        let faint = |color: &'static str| Some(Style::stroke(color, 0.15).dashed(0.5).faded(0.7));
        // The viewing direction, through both centres.
        curve3(svg, panel, camera, &solids, &[from, to], Style::stroke(INK, 0.22), faint(INK));
        arrow3(svg, panel, camera, sub(to, scale(0.1, d)), to, INK, 0.22);
        // The rotation axis from the plug's centre, and the angle about it.
        curve3(svg, panel, camera, &solids, &[centre, tip], Style::stroke(AXIS, 0.25), faint(AXIS));
        arrow3(svg, panel, camera, sub(tip, scale(0.1, axis)), tip, AXIS, 0.25);
        curve3(svg, panel, camera, &solids, &circle, Style::stroke(AXIS, 0.15), None);
        let (start, end) = (centre_angle - angle / 2.0, centre_angle + angle / 2.0);
        let arc: Vec<V3> = (0..=40).map(|k| around(start + (end - start) * k as f64 / 40.0)).collect();
        let mut sector: Vec<V2> = vec![camera.project(hub)];
        sector.extend(arc.iter().map(|&p| camera.project(p)));
        svg.polygon(panel, &sector, Style::fill(AXIS, 0.3));
        for p in [arc[0], arc[40]] {
            svg.line(panel, camera.project(hub), camera.project(p), Style::stroke(AXIS, 0.15));
        }
        svg.polyline(panel, &arc[..37].iter().map(|&p| camera.project(p)).collect::<Vec<_>>(), Style::stroke(AXIS, 0.25));
        svg.arrow(panel, camera.project(arc[36]), camera.project(arc[40]), AXIS, 0.25, 1.3);
    };
    three_parts_with(&passage, &look, &places, 10.0, &decorate, &ends)
}

/// The radial extent of a convex polygon containing the origin, in the
/// direction `w` (a unit vector).
fn extent(polygon: &[V2], w: V2) -> f64 {
    sides(polygon)
        .into_iter()
        .filter(|(n, _)| n[0] * w[0] + n[1] * w[1] > 1e-12)
        .map(|(n, h)| h / (n[0] * w[0] + n[1] * w[1]))
        .fold(f64::INFINITY, f64::min)
}

/// Two congruent cubes in the orientations of Nieuwland's passage: the
/// hole seen along `(2, 2, -1)`, the plug with a face across it.
fn margin_passage() -> Passage {
    let n = nieuwland_passage();
    let [e1, e2, d] = n.frame;
    let mut plug = Vec::new();
    for (a, b, c) in [(1.0, 1.0, 1.0), (1.0, 1.0, -1.0), (1.0, -1.0, 1.0), (1.0, -1.0, -1.0)] {
        for sign in [1.0, -1.0] {
            plug.push(scale(0.5 * sign, add(add(scale(a, e1), scale(b, e2)), scale(c, d))));
        }
    }
    Passage { hole: cube(0.5), plug: Polyhedron::new(plug, 1e-9), frame: n.frame }
}

pub fn margin() -> String {
    let passage = margin_passage();
    let (hole, plug) = (passage.hole_shadow(), passage.plug_shadow());
    let samples = 3600;
    let ratio: Vec<(f64, f64)> = (0..samples)
        .map(|k| {
            let a = 2.0 * PI * k as f64 / samples as f64;
            let w = [a.cos(), a.sin()];
            (a, extent(&plug, w) / extent(&hole, w))
        })
        .collect();
    let (best, top) = ratio.iter().copied().fold((0.0, 0.0), |b, r| if r.1 > b.1 + 1e-12 { r } else { b });
    eprintln!("margin of the congruent cubes: {top:.6} (1/margin {:.6})", 1.0 / top);
    let (w, h, gap) = (52.0, 52.0, 8.0);
    let graph_w = 62.0;
    let mut svg = Svg::new(w + gap + graph_w + 2.0, h + 2.0);
    let panel = fit(&mut svg, [1.0, 1.0, w, h], &hole, 0.03, false);
    // The two shadows, and faint rays every fifteen degrees.
    draw_shadows(&mut svg, &panel, &passage);
    for k in 0..24 {
        let a = PI * k as f64 / 12.0;
        let wv = [a.cos(), a.sin()];
        let e = extent(&hole, wv) * 1.08;
        svg.line(&panel, [0.0, 0.0], [e * wv[0], e * wv[1]], Style::stroke(INK, 0.08).faded(0.45));
    }
    // The ray where the ratio is largest, with both extents along it.
    let wv = [best.cos(), best.sin()];
    let (eh, ep) = (extent(&hole, wv), extent(&plug, wv));
    svg.line(&panel, [0.0, 0.0], [1.12 * eh * wv[0], 1.12 * eh * wv[1]], Style::stroke(INK, 0.15));
    // The two extents side by side, half a millimetre either side of it.
    let off = 0.55 / panel.scale;
    let side = [-wv[1] * off, wv[0] * off];
    for (e, sign, color) in [(eh, 1.0, HOLE), (ep, -1.0, PLUG)] {
        let o = [sign * side[0], sign * side[1]];
        svg.line(&panel, o, [o[0] + e * wv[0], o[1] + e * wv[1]], Style::stroke(color, 0.5));
    }
    svg.dot(&panel, [0.0, 0.0], 0.45, INK);
    // The ratio over a full turn of directions, with the level 1 dashed
    // and the maximum marked.
    let (lo, hi) = (0.6, 1.05);
    let gx = |a: f64| 1.0 + w + gap + graph_w * a / (2.0 * PI);
    let gy = |v: f64| 1.0 + h * (hi - v) / (hi - lo);
    svg.canvas_line([gx(0.0), gy(lo)], [gx(2.0 * PI), gy(lo)], Style::stroke(INK, 0.15));
    svg.canvas_line([gx(0.0), gy(lo)], [gx(0.0), gy(hi)], Style::stroke(INK, 0.15));
    svg.canvas_line([gx(0.0), gy(1.0)], [gx(2.0 * PI), gy(1.0)], Style::stroke(INK, 0.12).dashed(0.6));
    svg.canvas_line([gx(0.0), gy(top)], [gx(2.0 * PI), gy(top)], Style::stroke(PLUG, 0.1).faded(0.6));
    for k in 0..24 {
        let a = PI * k as f64 / 12.0;
        let wv = [a.cos(), a.sin()];
        let v = extent(&plug, wv) / extent(&hole, wv);
        svg.canvas_line([gx(a), gy(lo)], [gx(a), gy(v)], Style::stroke(INK, 0.08).faded(0.35));
    }
    let curve: Vec<V2> = ratio.iter().map(|&(a, v)| [gx(a), gy(v)]).chain([[gx(2.0 * PI), gy(ratio[0].1)]]).collect();
    let d: Vec<String> = curve.iter().enumerate().map(|(i, p)| format!("{}{:.3} {:.3}", if i == 0 { "M" } else { "L" }, p[0], p[1])).collect();
    svg.raw(&format!(r#"<path d="{}" fill="none" stroke="{PLUG}" stroke-width="0.25" stroke-linejoin="round"/>"#, d.concat()));
    svg.canvas_line([gx(best), gy(lo)], [gx(best), gy(top)], Style::stroke(INK, 0.15));
    svg.raw(&format!(r#"<circle cx="{:.3}" cy="{:.3}" r="0.45" fill="{PLUG}"/>"#, gx(best), gy(top)));
    svg.finish()
}

/// The rotation vectors allowed by the twelve comparisons
/// `|rᵢ| + (φ - 1)|rᵢ₊₁| ≤ 2 - φ`.
fn rotation_polytope() -> Polyhedron {
    let mut planes: Vec<(V3, f64)> = Vec::new();
    for i in 0..3 {
        for (a, b) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
            let mut n = [0.0; 3];
            n[i] = a;
            n[(i + 1) % 3] = b * (PHI - 1.0);
            planes.push((n, 2.0 - PHI));
        }
    }
    let mut corners: Vec<V3> = Vec::new();
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let [(a, p), (b, q), (c, r)] = [planes[i], planes[j], planes[k]];
                let det = dot(a, cross(b, c));
                if det.abs() < 1e-9 {
                    continue;
                }
                let x = scale(1.0 / det, add(add(scale(p, cross(b, c)), scale(q, cross(c, a))), scale(r, cross(a, b))));
                let inside = planes.iter().all(|&(n, h)| dot(n, x) <= h + 1e-9);
                if inside && !corners.iter().any(|&y| dot(sub(x, y), sub(x, y)) < 1e-18) {
                    corners.push(x);
                }
            }
        }
    }
    Polyhedron::new(corners, 1e-9)
}

/// Edges of a solid seen by `camera`: solid on the side facing the viewer,
/// dashed on the far side.
fn draw_edges(svg: &mut Svg, panel: &Panel, camera: &Camera, body: &Polyhedron, color: &'static str, width: f64) {
    for (a, b, faces) in body.edges() {
        let style = Style::stroke(color, width);
        let style = if back_edge(body, faces, camera.toward) { style.dashed(DASH) } else { style };
        svg.line(panel, camera.project(body.vertices[a]), camera.project(body.vertices[b]), style);
    }
}

pub fn reduced_ranges() -> String {
    let camera = rid_camera();
    let (w, h, gap) = (60.0, 52.0, 8.0);
    let mut svg = Svg::new(2.0 * w + gap + 2.0, h + 2.0);
    // The viewing cone from the origin through the viewing triangle in the
    // plane z = 1, inside the rectangle of B₀'s first two intervals.
    let corners = [[0.0, 0.0, 1.0], [1.0 / PHI, 0.0, 1.0], [0.0, 1.0 / (PHI * PHI), 1.0]];
    let far = corners.map(|c| scale(1.35, c));
    let rectangle = [[0.0, 0.0, 1.0], [2.0 / 3.0, 0.0, 1.0], [2.0 / 3.0, 0.4, 1.0], [0.0, 0.4, 1.0]];
    let axes = [([-0.25, 0.0, 0.0], [0.95, 0.0, 0.0]), ([0.0, -0.25, 0.0], [0.0, 0.7, 0.0]), ([0.0, 0.0, -0.2], [0.0, 0.0, 1.5])];
    let mut seen: Vec<V2> = far.iter().chain(&rectangle).map(|&v| camera.project(v)).collect();
    seen.extend(axes.iter().flat_map(|&(a, b)| [camera.project(a), camera.project(b)]));
    let left = fit(&mut svg, [1.0, 1.0, w, h], &seen, 0.04, false);
    for (a, b) in axes {
        svg.line(&left, camera.project(a), camera.project(b), Style::stroke(INK, 0.1).faded(0.5));
    }
    // The cone's sides, and its edges: dashed where both sides at an edge
    // face away from the viewer.
    let centroid = scale(1.0 / 3.0, add(add(far[0], far[1]), far[2]));
    let outward = |k: usize| {
        let n = cross(far[k], far[(k + 1) % 3]);
        if dot(n, centroid) > 0.0 { scale(-1.0, n) } else { n }
    };
    for k in 0..3 {
        let side = [[0.0; 3], far[k], far[(k + 1) % 3]].map(|v| camera.project(v));
        svg.polygon(&left, &side, Style::fill(INK, 0.07));
    }
    for k in 0..3 {
        let back = dot(outward((k + 2) % 3), camera.toward) < 0.0 && dot(outward(k), camera.toward) < 0.0;
        let style = if back { Style::stroke(INK, 0.15).dashed(DASH) } else { Style::stroke(INK, 0.15) };
        svg.line(&left, camera.project([0.0; 3]), camera.project(far[k]), style);
    }
    svg.polygon(&left, &rectangle.map(|v| camera.project(v)), Style::stroke(INK, 0.12).faded(0.6));
    svg.polygon(&left, &corners.map(|v| camera.project(v)), Style::fill(INK, 0.3).with_stroke(INK, 0.2));
    svg.dot(&left, camera.project([0.0; 3]), 0.4, INK);
    // The rotation box and, inside it, the rotation comparisons.
    let polytope = rotation_polytope();
    let bound = 0.4;
    let cube_box = cube(bound);
    let seen: Vec<V2> = cube_box.vertices.iter().map(|&v| camera.project(v)).collect();
    let right = fit(&mut svg, [1.0 + w + gap, 1.0, w, h], &seen, 0.12, false);
    for f in 0..polytope.faces.len() {
        let ring: Vec<V2> = polytope.faces[f].iter().map(|&i| camera.project(polytope.vertices[i])).collect();
        svg.polygon(&right, &ring, Style::fill(PLUG, 0.06));
    }
    draw_edges(&mut svg, &right, &camera, &polytope, PLUG, 0.18);
    draw_edges(&mut svg, &right, &camera, &cube_box, INK, 0.15);
    for axis in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
        svg.line(&right, camera.project(scale(-0.55, axis)), camera.project(scale(0.55, axis)), Style::stroke(INK, 0.1).faded(0.5));
    }
    svg.dot(&right, camera.project([0.0; 3]), 0.4, INK);
    svg.finish()
}

/// The rhombicosidodecahedron alone, opaque and shaded, for the title page.
pub fn rid_emblem() -> String {
    let camera = rid_camera();
    let solid = rid();
    let seen: Vec<V2> = solid.vertices.iter().map(|&v| camera.project(v)).collect();
    let size = 30.0;
    let mut svg = Svg::new(size, size);
    let panel = fit(&mut svg, [0.0, 0.0, size, size], &seen, 0.02, false);
    draw_opaque_with(&mut svg, &panel, &camera, &[(&solid, PLUG)], false);
    svg.finish()
}
