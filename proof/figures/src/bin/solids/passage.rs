//! A passage of one solid (the plug) through a copy of it (the hole), and
//! the three drawings of it: the hole with the tunnel cut out and the plug
//! in it, the two shadows, and enlargements where the clearance is least.
use super::geometry::*;
use super::svg::{Panel, Style, Svg};

pub const HOLE: &str = "#e08214";
pub const PLUG: &str = "#2166ac";
pub const INK: &str = "#404040";
/// Dashes and gaps of edges on the far side of their solid, in millimetres.
pub const DASH: f64 = 0.4;

/// Whether the edge between faces `f` and `g` of a convex body lies on its
/// far side, hidden by the body itself from a viewer in the direction
/// `toward`.
pub fn back_edge(body: &Polyhedron, [f, g]: [usize; 2], toward: V3) -> bool {
    dot(body.normal(f), toward) < 0.0 && dot(body.normal(g), toward) < 0.0
}

/// A style's hidden version: the same stroke, dashed.
fn hidden_style(style: Style) -> Style {
    style.dashed(DASH)
}

pub struct Passage {
    pub hole: Polyhedron,
    pub plug: Polyhedron,
    /// Two screen axes and the direction of motion, orthonormal and
    /// right-handed.
    pub frame: [V3; 3],
}

impl Passage {
    pub fn screen(&self, v: V3) -> V2 {
        [dot(v, self.frame[0]), dot(v, self.frame[1])]
    }

    fn shadow(&self, body: &Polyhedron) -> Vec<V2> {
        hull(&body.vertices.iter().map(|&v| self.screen(v)).collect::<Vec<_>>())
    }

    pub fn hole_shadow(&self) -> Vec<V2> {
        self.shadow(&self.hole)
    }

    pub fn plug_shadow(&self) -> Vec<V2> {
        self.shadow(&self.plug)
    }

    /// The vertices of the plug's shadow with their clearance, the distance
    /// to the boundary of the hole's shadow (negative outside), least first.
    pub fn clearances(&self) -> Vec<(V2, f64)> {
        let hole = self.hole_shadow();
        let mut c: Vec<(V2, f64)> = self.plug_shadow().into_iter().map(|p| (p, depth(&hole, p))).collect();
        c.sort_by(|a, b| a.1.total_cmp(&b.1));
        c
    }

    /// The plug moved by `mu` along the direction of motion.
    pub fn advance(mut self, mu: f64) -> Self {
        let d = scale(mu, self.frame[2]);
        self.plug = self.plug.map(|v| add(v, d));
        self
    }

    /// Places for enlargements: the `count` tightest plug vertices, skipping
    /// any within `apart` of one already taken or of its reflection
    /// through the origin (for centrally symmetric passages).
    pub fn tight_places(&self, count: usize, apart: f64) -> Vec<(V2, f64)> {
        let mut out: Vec<(V2, f64)> = Vec::new();
        for (p, c) in self.clearances() {
            let near = |q: V2| (p[0] - q[0]).hypot(p[1] - q[1]) < apart || (p[0] + q[0]).hypot(p[1] + q[1]) < apart;
            if out.len() < count && !out.iter().any(|&(q, _)| near(q)) {
                out.push((p, c));
            }
        }
        out
    }
}

/// The part of a polygon (in space) where `n·x ≤ h`.
fn clip(polygon: &[V3], n: V3, h: f64) -> Vec<V3> {
    let mut out = Vec::new();
    for i in 0..polygon.len() {
        let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        let (da, db) = (dot(n, a) - h, dot(n, b) - h);
        if da <= 0.0 {
            out.push(a);
        }
        if (da <= 0.0) != (db <= 0.0) {
            out.push(add(a, scale(da / (da - db), sub(b, a))));
        }
    }
    out
}

fn area(polygon: &[V3]) -> f64 {
    let mut n = [0.0; 3];
    for i in 1..polygon.len().saturating_sub(1) {
        n = add(n, cross(sub(polygon[i], polygon[0]), sub(polygon[i + 1], polygon[0])));
    }
    dot(n, n).sqrt() / 2.0
}

/// How the three-dimensional view is drawn.
pub struct Look {
    pub camera: Camera,
    /// Cut the tunnel out of the hole.
    pub tunnel: bool,
    /// The opacity of the hole's faces.
    pub hole_fill: f64,
    /// Draw both solids opaque, with hidden edges omitted, instead.
    pub opaque: bool,
}

/// Whether `p` is hidden from the camera looking along `-toward` by one of
/// the convex `bodies`: whether the ray from `p` towards the viewer meets
/// one of them (which it does from inside).
pub fn hidden(bodies: &[&Polyhedron], toward: V3, p: V3) -> bool {
    bodies.iter().any(|body| {
        let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
        for f in 0..body.faces.len() {
            let n = body.normal(f);
            let (a, b) = (dot(n, toward), dot(n, sub(body.vertices[body.faces[f][0]], p)));
            if a > 1e-12 {
                hi = hi.min(b / a);
            } else if a < -1e-12 {
                lo = lo.max(b / a);
            } else if b < 0.0 {
                return false;
            }
        }
        lo < hi
    })
}

/// A curve in space: its visible pieces in `style`, its hidden pieces in
/// `behind` (or not at all).
pub fn curve3(svg: &mut Svg, panel: &Panel, camera: &Camera, bodies: &[&Polyhedron], points: &[V3], style: Style, behind: Option<Style>) {
    let mut fine: Vec<V3> = Vec::new();
    for w in points.windows(2) {
        for k in 0..40 {
            fine.push(add(w[0], scale(k as f64 / 40.0, sub(w[1], w[0]))));
        }
    }
    fine.extend(points.last());
    let mut run: Vec<V2> = Vec::new();
    let mut state = None;
    for &p in &fine {
        let h = hidden(bodies, camera.toward, p);
        if state.is_some_and(|s| s != h) {
            run.push(camera.project(p));
            draw_run(svg, panel, &run, if state == Some(true) { behind } else { Some(style) });
            run.clear();
        }
        state = Some(h);
        run.push(camera.project(p));
    }
    draw_run(svg, panel, &run, if state == Some(true) { behind } else { Some(style) });
}

fn draw_run(svg: &mut Svg, panel: &Panel, run: &[V2], style: Option<Style>) {
    if let Some(style) = style {
        svg.polyline(panel, run, style);
    }
}

/// Convex solids drawn opaque, back to front: their faces seen from the
/// front, shaded, and the edges of those faces.
pub fn draw_opaque(svg: &mut Svg, panel: &Panel, camera: &Camera, bodies: &[(&Polyhedron, &'static str)]) {
    draw_opaque_with(svg, panel, camera, bodies, true)
}

/// As [`draw_opaque`], with the edges on each solid's far side dashed or
/// left out.
pub fn draw_opaque_with(svg: &mut Svg, panel: &Panel, camera: &Camera, bodies: &[(&Polyhedron, &'static str)], back_edges: bool) {
    let mut order: Vec<&(&Polyhedron, &'static str)> = bodies.iter().collect();
    let depth = |b: &Polyhedron| {
        let c = b.vertices.iter().fold([0.0; 3], |a, &v| add(a, v));
        dot(c, camera.toward) / b.vertices.len() as f64
    };
    order.sort_by(|a, b| depth(a.0).total_cmp(&depth(b.0)));
    let light = unit(add(add(camera.toward, scale(0.6, camera.up)), scale(-0.4, camera.right)));
    for &&(body, color) in &order {
        for f in 0..body.faces.len() {
            let n = body.normal(f);
            if dot(n, camera.toward) <= 0.0 {
                continue;
            }
            let ring: Vec<V2> = body.faces[f].iter().map(|&i| camera.project(body.vertices[i])).collect();
            svg.polygon(panel, &ring, Style::fill("#ffffff", 1.0));
            svg.polygon(panel, &ring, Style::fill(color, 0.08 + 0.2 * (1.0 - dot(n, light).max(0.0))));
        }
        for (a, b, [f, g]) in body.edges() {
            let (ff, fg) = (dot(body.normal(f), camera.toward) > 0.0, dot(body.normal(g), camera.toward) > 0.0);
            let (p, q) = (camera.project(body.vertices[a]), camera.project(body.vertices[b]));
            if ff || fg {
                let width = if ff && fg { 0.1 } else { 0.18 };
                svg.line(panel, p, q, Style::stroke(color, width));
            } else if back_edges {
                svg.line(panel, p, q, hidden_style(Style::stroke(color, 0.1)));
            }
        }
    }
}

/// The hole (faint, with the tunnel cut out) and the plug (transparent,
/// every edge drawn) seen by the camera.
pub fn draw_solids(svg: &mut Svg, panel: &Panel, passage: &Passage, look: &Look) {
    if look.opaque {
        return draw_opaque(svg, panel, &look.camera, &[(&passage.hole, HOLE), (&passage.plug, PLUG)]);
    }
    let cam = look.camera;
    let pr = |v: V3| cam.project(v);
    let [e1, e2, d] = passage.frame;
    // The prism swept by the plug's shadow, as half-spaces.
    let prism: Vec<(V3, f64)> = if look.tunnel {
        sides(&passage.plug_shadow()).into_iter().map(|(m, c)| (add(scale(m[0], e1), scale(m[1], e2)), c)).collect()
    } else {
        Vec::new()
    };
    let hole = &passage.hole;
    let speck = hole.vertices.iter().map(|&v| dot(v, v).sqrt()).fold(0.0, f64::max) * 2e-3;
    // The hole's faces outside the prism.
    for face in &hole.faces {
        let polygon: Vec<V3> = face.iter().map(|&i| hole.vertices[i]).collect();
        let mut cut = polygon.clone();
        for &(n, h) in &prism {
            cut = clip(&cut, n, h);
        }
        let whole = area(&polygon);
        let removed = if prism.is_empty() || cut.len() < 3 { 0.0 } else { area(&cut) };
        if removed > whole * (1.0 - 1e-6) {
            continue;
        }
        let mut rings = vec![polygon.iter().map(|&v| pr(v)).collect::<Vec<_>>()];
        if removed > 0.0 {
            rings.push(cut.iter().map(|&v| pr(v)).collect());
        }
        svg.region(panel, &rings, Style::fill(HOLE, look.hole_fill));
    }
    // The hole's edges outside the prism.
    for (a, b, faces) in hole.edges() {
        let (p, q) = (hole.vertices[a], hole.vertices[b]);
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for &(n, h) in &prism {
            let (dp, dq) = (dot(n, p) - h, dot(n, q) - h);
            if dp > 0.0 && dq > 0.0 {
                t0 = 1.0;
                t1 = 0.0;
            } else if dp > 0.0 {
                t0 = t0.max(dp / (dp - dq));
            } else if dq > 0.0 {
                t1 = t1.min(dp / (dp - dq));
            }
        }
        let at = |t: f64| pr(add(p, scale(t, sub(q, p))));
        let style = Style::stroke(HOLE, 0.12).faded(0.55);
        let style = if back_edge(hole, faces, cam.toward) { hidden_style(style) } else { style };
        let pieces = if prism.is_empty() || t0 >= t1 { vec![(0.0, 1.0)] } else { vec![(0.0, t0), (t1, 1.0)] };
        for (a, b) in pieces {
            // Pieces far too short to see would only leave specks.
            if (b - a) * dot(sub(q, p), sub(q, p)).sqrt() > speck {
                svg.line(panel, at(a), at(b), style);
            }
        }
    }
    // The tunnel's walls: each side of the prism, inside the hole.
    if !prism.is_empty() {
        let shadow = passage.plug_shadow();
        let reach = hole.vertices.iter().map(|&v| dot(v, v).sqrt()).fold(0.0, f64::max) * 4.0;
        let lift = |p: V2| add(scale(p[0], e1), scale(p[1], e2));
        for i in 0..shadow.len() {
            let (a, b) = (lift(shadow[i]), lift(shadow[(i + 1) % shadow.len()]));
            let mut wall = vec![sub(a, scale(reach, d)), sub(b, scale(reach, d)), add(b, scale(reach, d)), add(a, scale(reach, d))];
            for f in 0..hole.faces.len() {
                let n = hole.normal(f);
                wall = clip(&wall, n, dot(n, hole.vertices[hole.faces[f][0]]));
            }
            let size = wall.iter().map(|&v| dot(sub(v, wall[0]), sub(v, wall[0])).sqrt()).fold(0.0, f64::max);
            if wall.len() >= 3 && size > speck {
                let ring: Vec<V2> = wall.iter().map(|&v| pr(v)).collect();
                svg.polygon(panel, &ring, Style::fill(HOLE, 0.14));
                // The wall's material lies outside the tunnel, so its outward
                // normal points into the tunnel, towards the prism's axis.
                let (m, _) = prism[i];
                let outline = Style::stroke(HOLE, 0.12);
                let outline = if dot(m, cam.toward) > 0.0 { hidden_style(outline) } else { outline };
                let mut closed = ring.clone();
                closed.push(ring[0]);
                svg.polyline(panel, &closed, outline);
            }
        }
    }
    // The plug, transparent.
    let plug = &passage.plug;
    for face in &plug.faces {
        let ring: Vec<V2> = face.iter().map(|&i| pr(plug.vertices[i])).collect();
        svg.polygon(panel, &ring, Style::fill(PLUG, 0.08));
    }
    for (a, b, faces) in plug.edges() {
        let style = Style::stroke(PLUG, 0.15);
        let style = if back_edge(plug, faces, cam.toward) { hidden_style(style) } else { style };
        svg.line(panel, pr(plug.vertices[a]), pr(plug.vertices[b]), style);
    }
}

/// A light opaque tint of the hole's colour.
pub const HOLE_TINT: &str = "#fbefe1";

/// Both shadows as wireframes of every projected edge, over a light tint of
/// the hole's shadow.
pub fn draw_shadows(svg: &mut Svg, panel: &Panel, passage: &Passage) {
    svg.polygon(panel, &passage.hole_shadow(), Style::fill(HOLE_TINT, 1.0));
    // The screen is right-handed, so the viewer lies along the direction of
    // motion; each body's edges on its far side are dashed.
    for (body, color) in [(&passage.hole, HOLE), (&passage.plug, PLUG)] {
        for (a, b, faces) in body.edges() {
            let style = Style::stroke(color, 0.1);
            let style = if back_edge(body, faces, passage.frame[2]) { hidden_style(style) } else { style };
            svg.line(panel, passage.screen(body.vertices[a]), passage.screen(body.vertices[b]), style);
        }
    }
}

/// The largest of 1, 2, 5 times a power of ten not above `x`.
pub fn round_down(x: f64) -> f64 {
    let p = 10f64.powf(x.log10().floor());
    [5.0, 2.0, 1.0].into_iter().map(|m| m * p).find(|&m| m <= x * (1.0 + 1e-9)).unwrap_or(p)
}

/// `×M` with large powers of ten raised.
pub fn magnification_label(m: f64) -> String {
    if m < 1e5 {
        format!("×{}", m.round())
    } else {
        let k = m.log10().floor();
        let lead = (m / 10f64.powf(k)).round();
        if lead == 1.0 { format!("×10^{}", k) } else { format!("×{}·10^{}", lead, k) }
    }
}

/// An enlargement of `overview` in `rect` around `centre`, magnified `m`
/// times: the window is outlined on the overview (a small circle if it
/// would be tinier than `min_mark` millimetres), joined to the
/// enlargement, whose magnification is written in its corner.
pub fn enlarge(svg: &mut Svg, overview: &Panel, rect: [f64; 4], centre: V2, m: f64, min_mark: f64) -> Panel {
    let detail = svg.panel(rect, centre, overview.scale * m, true);
    let [hx, hy] = detail.half();
    let a = overview.map([centre[0] - hx, centre[1] + hy]);
    let b = overview.map([centre[0] + hx, centre[1] - hy]);
    let ink = Style::stroke(INK, 0.15);
    let thin = Style::stroke(INK, 0.1).faded(0.6);
    let (w, h) = (b[0] - a[0], b[1] - a[1]);
    let mark = if w.max(h) >= min_mark {
        svg.frame([a[0], a[1], w, h], ink);
        [a[0], a[1], b[0], b[1]]
    } else {
        let c = overview.map(centre);
        let r = min_mark / 2.0;
        svg.frame([c[0] - r, c[1] - r, 2.0 * r, 2.0 * r], ink);
        [c[0] - r, c[1] - r, c[0] + r, c[1] + r]
    };
    // Join the two nearest corner pairs on the side facing the enlargement.
    let [x, y, rw, rh] = rect;
    let (ox, oy) = ((mark[0] + mark[2]) / 2.0, (mark[1] + mark[3]) / 2.0);
    let (dx, dy) = (x + rw / 2.0 - ox, y + rh / 2.0 - oy);
    let pairs = if dx.abs() * rh >= dy.abs() * rw {
        let (mx, px) = if dx > 0.0 { (mark[2], x) } else { (mark[0], x + rw) };
        [([mx, mark[1]], [px, y]), ([mx, mark[3]], [px, y + rh])]
    } else {
        let (my, py) = if dy > 0.0 { (mark[3], y) } else { (mark[1], y + rh) };
        [([mark[0], my], [x, py]), ([mark[2], my], [x + rw, py])]
    };
    for (p, q) in pairs {
        svg.canvas_line(p, q, thin);
    }
    svg.frame(rect, ink);
    detail
}

/// Writes the magnification in the top-left corner of `rect`.
pub fn label_corner(svg: &mut Svg, rect: [f64; 4], m: f64) {
    svg.label([rect[0] + 0.9, rect[1] + 2.4], 1.9, &magnification_label(m));
}
