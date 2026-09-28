//! A minimal SVG writer in millimetres: panels map a world window to a
//! rectangle of the canvas with one scale for both axes, `y` pointing up.
use super::geometry::V2;
use std::fmt::Write;

/// How a shape is painted.
#[derive(Clone, Copy)]
pub struct Style {
    pub fill: Option<(&'static str, f64)>,
    pub stroke: Option<(&'static str, f64)>,
    pub dash: Option<f64>,
    pub opacity: f64,
}

impl Style {
    pub fn fill(color: &'static str, opacity: f64) -> Self {
        Style { fill: Some((color, opacity)), stroke: None, dash: None, opacity: 1.0 }
    }

    pub fn stroke(color: &'static str, width: f64) -> Self {
        Style { fill: None, stroke: Some((color, width)), dash: None, opacity: 1.0 }
    }

    pub fn dashed(mut self, length: f64) -> Self {
        self.dash = Some(length);
        self
    }

    pub fn faded(mut self, opacity: f64) -> Self {
        self.opacity = opacity;
        self
    }

    pub fn with_stroke(mut self, color: &'static str, width: f64) -> Self {
        self.stroke = Some((color, width));
        self
    }

    fn attributes(&self) -> String {
        let mut a = String::new();
        match self.fill {
            Some((c, o)) => write!(a, r#" fill="{c}" fill-opacity="{}""#, num(o)).unwrap(),
            None => a.push_str(r#" fill="none""#),
        }
        if let Some((c, w)) = self.stroke {
            write!(a, r#" stroke="{c}" stroke-width="{}" stroke-linecap="round" stroke-linejoin="round""#, num(w)).unwrap();
        }
        if let Some(d) = self.dash {
            write!(a, r#" stroke-dasharray="{} {}""#, num(d), num(d)).unwrap();
        }
        if self.opacity < 1.0 {
            write!(a, r#" opacity="{}""#, num(self.opacity)).unwrap();
        }
        a
    }
}

/// A rectangle of the canvas showing the world window centred at `centre`
/// with `scale` millimetres per world unit.
#[derive(Clone, Copy)]
pub struct Panel {
    pub rect: [f64; 4],
    pub centre: V2,
    pub scale: f64,
    clip: Option<usize>,
}

impl Panel {
    /// The canvas position of a world point.
    pub fn map(&self, p: V2) -> V2 {
        let [x, y, w, h] = self.rect;
        [x + w / 2.0 + (p[0] - self.centre[0]) * self.scale, y + h / 2.0 - (p[1] - self.centre[1]) * self.scale]
    }

    /// Half the width and height of the world window.
    pub fn half(&self) -> V2 {
        [self.rect[2] / 2.0 / self.scale, self.rect[3] / 2.0 / self.scale]
    }
}

pub struct Svg {
    width: f64,
    height: f64,
    body: String,
    clips: usize,
}

fn num(x: f64) -> String {
    let s = format!("{:.4}", x);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.into() }
}

/// The part of a polygon inside the panel's window, enlarged by a tenth
/// (Sutherland–Hodgman), in coordinates relative to the window's centre:
/// far-away points never reach the renderer, whatever the magnification.
fn clip_polygon(panel: &Panel, pts: &[V2]) -> Vec<V2> {
    let [hx, hy] = panel.half();
    let mut ring: Vec<V2> = pts.iter().map(|p| [p[0] - panel.centre[0], p[1] - panel.centre[1]]).collect();
    for (axis, sign, bound) in [(0, 1.0, hx), (0, -1.0, hx), (1, 1.0, hy), (1, -1.0, hy)] {
        let limit = bound * 1.1;
        let inside = |p: &V2| sign * p[axis] <= limit;
        let mut next = Vec::new();
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if inside(&a) {
                next.push(a);
            }
            if inside(&a) != inside(&b) {
                let t = (limit - sign * a[axis]) / (sign * (b[axis] - a[axis]));
                next.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        }
        ring = next;
        if ring.is_empty() {
            break;
        }
    }
    ring
}

/// The part of a segment inside the enlarged window, relative to its centre.
fn clip_segment(panel: &Panel, a: V2, b: V2) -> Option<[V2; 2]> {
    let [hx, hy] = panel.half();
    let a = [a[0] - panel.centre[0], a[1] - panel.centre[1]];
    let b = [b[0] - panel.centre[0], b[1] - panel.centre[1]];
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (axis, bound) in [(0, hx * 1.1), (1, hy * 1.1)] {
        let d = b[axis] - a[axis];
        for (p, q) in [(-d, a[axis] + bound), (d, bound - a[axis])] {
            if p == 0.0 {
                if q < 0.0 {
                    return None;
                }
            } else {
                let r = q / p;
                if p < 0.0 { t0 = t0.max(r) } else { t1 = t1.min(r) }
            }
        }
    }
    (t0 <= t1).then(|| [0, 1].map(|k| {
        let t = if k == 0 { t0 } else { t1 };
        [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
    }))
}

/// Path data for points given relative to the window's centre.
fn relative(panel: &Panel, pts: &[V2]) -> String {
    let [x, y, w, h] = panel.rect;
    let mut d = String::new();
    for (i, p) in pts.iter().enumerate() {
        let q = [x + w / 2.0 + p[0] * panel.scale, y + h / 2.0 - p[1] * panel.scale];
        write!(d, "{}{} {}", if i == 0 { "M" } else { "L" }, num(q[0]), num(q[1])).unwrap();
    }
    d
}

impl Svg {
    pub fn new(width: f64, height: f64) -> Self {
        Svg { width, height, body: String::new(), clips: 0 }
    }

    /// A panel filling `rect = [x, y, width, height]` (mm) around `centre`;
    /// with `clip`, drawing outside the rectangle is cut off.
    pub fn panel(&mut self, rect: [f64; 4], centre: V2, scale: f64, clip: bool) -> Panel {
        let clip = clip.then(|| {
            self.clips += 1;
            writeln!(
                self.body,
                r#"<clipPath id="c{}"><rect x="{}" y="{}" width="{}" height="{}"/></clipPath>"#,
                self.clips,
                num(rect[0]),
                num(rect[1]),
                num(rect[2]),
                num(rect[3])
            )
            .unwrap();
            self.clips
        });
        Panel { rect, centre, scale, clip }
    }

    fn element(&mut self, panel: &Panel, body: String, style: &Style) {
        let clip = panel.clip.map(|c| format!(r#" clip-path="url(#c{c})""#)).unwrap_or_default();
        writeln!(self.body, "<path d=\"{body}\"{}{clip}/>", style.attributes()).unwrap();
    }

    pub fn polygon(&mut self, panel: &Panel, pts: &[V2], style: Style) {
        let ring = clip_polygon(panel, pts);
        if ring.len() >= 2 {
            self.element(panel, relative(panel, &ring) + "Z", &style);
        }
    }

    /// Polygons filled with the even-odd rule, so that an inner polygon
    /// cuts a hole in an outer one.
    pub fn region(&mut self, panel: &Panel, rings: &[Vec<V2>], style: Style) {
        let d: String = rings
            .iter()
            .map(|r| clip_polygon(panel, r))
            .filter(|r| r.len() >= 3)
            .map(|r| relative(panel, &r) + "Z")
            .collect();
        if !d.is_empty() {
            let clip = panel.clip.map(|c| format!(r#" clip-path="url(#c{c})""#)).unwrap_or_default();
            writeln!(self.body, "<path d=\"{d}\" fill-rule=\"evenodd\"{}{clip}/>", style.attributes()).unwrap();
        }
    }

    /// A polyline as one path (so that dashes run on across its corners),
    /// broken only where it leaves the panel's window.
    pub fn polyline(&mut self, panel: &Panel, pts: &[V2], style: Style) {
        let mut pieces: Vec<Vec<V2>> = Vec::new();
        for w in pts.windows(2) {
            if let Some([a, b]) = clip_segment(panel, w[0], w[1]) {
                match pieces.last_mut() {
                    Some(last) if last.last().is_some_and(|p| (p[0] - a[0]).abs() + (p[1] - a[1]).abs() < 1e-12) => last.push(b),
                    _ => pieces.push(vec![a, b]),
                }
            }
        }
        let d: String = pieces.iter().map(|piece| relative(panel, piece)).collect();
        if !d.is_empty() {
            self.element(panel, d, &style);
        }
    }

    pub fn line(&mut self, panel: &Panel, a: V2, b: V2, style: Style) {
        self.polyline(panel, &[a, b], style);
    }

    /// A disc of radius `r` millimetres.
    pub fn dot(&mut self, panel: &Panel, p: V2, r: f64, color: &'static str) {
        let q = panel.map(p);
        writeln!(self.body, r#"<circle cx="{}" cy="{}" r="{}" fill="{color}"/>"#, num(q[0]), num(q[1]), num(r)).unwrap();
    }

    /// An arrow from `a` to `b` with a filled head `head` millimetres long.
    pub fn arrow(&mut self, panel: &Panel, a: V2, b: V2, color: &'static str, width: f64, head: f64) {
        let (p, q) = (panel.map(a), panel.map(b));
        let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
        let l = dx.hypot(dy);
        let (ux, uy) = (dx / l, dy / l);
        let base = [q[0] - ux * head, q[1] - uy * head];
        let side = [-uy * head * 0.35, ux * head * 0.35];
        writeln!(
            self.body,
            r#"<path d="M{} {}L{} {}" stroke="{color}" stroke-width="{}" stroke-linecap="round" fill="none"/>"#,
            num(p[0]),
            num(p[1]),
            num(base[0] + ux * head * 0.3),
            num(base[1] + uy * head * 0.3),
            num(width)
        )
        .unwrap();
        writeln!(
            self.body,
            r#"<path d="M{} {}L{} {}L{} {}Z" fill="{color}"/>"#,
            num(q[0]),
            num(q[1]),
            num(base[0] + side[0]),
            num(base[1] + side[1]),
            num(base[0] - side[0]),
            num(base[1] - side[1])
        )
        .unwrap();
    }

    /// A rectangle given in canvas millimetres.
    pub fn frame(&mut self, rect: [f64; 4], style: Style) {
        writeln!(
            self.body,
            r#"<rect x="{}" y="{}" width="{}" height="{}"{}/>"#,
            num(rect[0]),
            num(rect[1]),
            num(rect[2]),
            num(rect[3]),
            style.attributes()
        )
        .unwrap();
    }

    /// A segment given in canvas millimetres.
    pub fn canvas_line(&mut self, a: V2, b: V2, style: Style) {
        writeln!(
            self.body,
            "<path d=\"M{} {}L{} {}\"{}/>",
            num(a[0]),
            num(a[1]),
            num(b[0]),
            num(b[1]),
            style.attributes()
        )
        .unwrap();
    }

    /// Small text at a canvas position, anchored at its start; `^` starts a
    /// raised exponent up to the next space or the end.
    pub fn label(&mut self, at: V2, size: f64, text: &str) {
        let mut spans = String::new();
        for (i, part) in text.split('^').enumerate() {
            if i == 0 {
                spans.push_str(part);
            } else {
                let (exp, rest) = part.split_once(' ').map(|(e, r)| (e, format!(" {r}"))).unwrap_or((part, String::new()));
                write!(
                    spans,
                    r#"<tspan dy="{}" font-size="{}">{exp}</tspan><tspan dy="{}">{rest}</tspan>"#,
                    num(-size * 0.38),
                    num(size * 0.7),
                    num(size * 0.38)
                )
                .unwrap();
            }
        }
        writeln!(
            self.body,
            r##"<text x="{}" y="{}" font-family="Latin Modern Roman" font-size="{}" fill="#404040">{spans}</text>"##,
            num(at[0]),
            num(at[1]),
            num(size)
        )
        .unwrap();
    }

    /// An element written as given, in canvas millimetres.
    pub fn raw(&mut self, element: &str) {
        self.body.push_str(element);
        self.body.push('\n');
    }

    pub fn finish(self) -> String {
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}mm\" height=\"{h}mm\" viewBox=\"0 0 {w} {h}\">\n<rect x=\"0\" y=\"0\" width=\"{w}\" height=\"{h}\" fill=\"#ffffff\"/>\n{}</svg>\n",
            self.body,
            w = num(self.width),
            h = num(self.height)
        )
    }
}
