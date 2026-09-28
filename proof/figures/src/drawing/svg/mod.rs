//! The SVG writer: shapes in page millimetres, one group per layer with its
//! explicit style and optional clip rectangle, exact decimal output.
use super::layout::{Canvas, Rectangle};
use super::style::{Color, Style};
use crate::content::{Point, Shape};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Signed;
use rid::arithmetic::exact::Q;
use std::fmt::Write as _;

/// Decimal places of every written number: 0.1 µm on the page.
pub const DECIMALS: usize = 4;

/// An SVG document under construction; one SVG unit is one millimetre.
pub struct Document {
    text: String,
    clips: usize,
}

impl Document {
    /// The page of the canvas, painted with the background unless it is `None`.
    pub fn new(canvas: &Canvas, background: Option<Color>) -> Self {
        let (w, h) = (decimal(canvas.width_mm()), decimal(canvas.height_mm()));
        let mut text = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" \
             width=\"{w}mm\" height=\"{h}mm\" viewBox=\"0 0 {w} {h}\">\n"
        );
        if let Some(color) = background {
            let fill = color.hex();
            writeln!(text, "<rect x=\"0\" y=\"0\" width=\"{w}\" height=\"{h}\" fill=\"{fill}\"/>")
                .expect("writing to a string");
        }
        Self { text, clips: 0 }
    }

    /// Adds one layer: the shapes, in page coordinates, drawn with the style
    /// and, when `clip` is given, clipped to that page rectangle.
    pub fn layer(&mut self, shapes: &[Shape], style: &Style, clip: Option<&Rectangle>) {
        let mut group = String::from("<g");
        if let Some([x, y]) = clip {
            self.clips += 1;
            let (left, top) = (decimal(x.lo()), decimal(y.lo()));
            let (width, height) = (decimal(&x.width()), decimal(&y.width()));
            writeln!(
                self.text,
                "<clipPath id=\"clip-{}\"><rect x=\"{left}\" y=\"{top}\" \
                 width=\"{width}\" height=\"{height}\"/></clipPath>",
                self.clips
            )
            .expect("writing to a string");
            write!(group, " clip-path=\"url(#clip-{})\"", self.clips).expect("writing to a string");
        }
        let fill = style.fill().map_or("none".to_owned(), |c| c.hex());
        let opacity = decimal(style.opacity());
        write!(group, " opacity=\"{opacity}\" fill=\"{fill}\" fill-rule=\"nonzero\"")
            .expect("writing to a string");
        match style.stroke() {
            Some(stroke) => write!(
                group,
                " stroke=\"{}\" stroke-width=\"{}\" \
                 stroke-linecap=\"round\" stroke-linejoin=\"round\"{}",
                stroke.color().hex(),
                decimal(stroke.width_mm()),
                stroke.dash_mm().map_or(String::new(), |d| {
                    let d = decimal(d);
                    format!(" stroke-dasharray=\"{d} {d}\"")
                })
            ),
            None => write!(group, " stroke=\"none\""),
        }
        .expect("writing to a string");
        self.text.push_str(&group);
        self.text.push_str(">\n");
        for shape in shapes {
            writeln!(self.text, "<path d=\"{}\"/>", path(shape)).expect("writing to a string");
        }
        self.text.push_str("</g>\n");
    }

    /// Small grey text starting at a page position (its baseline), in the
    /// article's typeface, on a thin white halo so that lines behind it do
    /// not cross it.
    pub fn label(&mut self, at: Point, size_mm: &Q, text: &str) {
        let escaped = text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
        writeln!(
            self.text,
            "<text x=\"{}\" y=\"{}\" font-family=\"Latin Modern Roman, LM Roman 10, serif\" \
             font-size=\"{}\" fill=\"#404040\" stroke=\"#ffffff\" stroke-width=\"0.5\" \
             stroke-linejoin=\"round\" paint-order=\"stroke\">{escaped}</text>",
            decimal(&at[0]),
            decimal(&at[1]),
            decimal(size_mm)
        )
        .expect("writing to a string");
    }

    pub fn finish(mut self) -> String {
        self.text.push_str("</svg>\n");
        self.text
    }
}

/// Path data: `M x y L x y …`, closed with `Z` for polygons; a circle as
/// two half-circle arcs.
fn path(shape: &Shape) -> String {
    let (points, closed): (&[Point], bool) = match shape {
        Shape::Segment(points) => (points, false),
        Shape::Polygon(points) => (points, true),
        Shape::Circle { centre, radius } => {
            let (x, y, r) = (&centre[0], decimal(&centre[1]), decimal(radius));
            let (left, right) = (decimal(&(x - radius)), decimal(&(x + radius)));
            return format!("M{left} {y}A{r} {r} 0 1 0 {right} {y}A{r} {r} 0 1 0 {left} {y}Z");
        }
    };
    let mut d = String::new();
    for (k, p) in points.iter().enumerate() {
        let command = if k == 0 { "M" } else { "L" };
        write!(d, "{command}{} {}", decimal(&p[0]), decimal(&p[1])).expect("writing to a string");
    }
    if closed {
        d.push('Z');
    }
    d
}

/// The rational rounded to `DECIMALS` places, halves away from zero,
/// without trailing zeros.
pub fn decimal(x: &Q) -> String {
    let unit = BigInt::from(10).pow(DECIMALS as u32);
    let (whole, rest) = (x.numer().abs() * &unit).div_rem(x.denom());
    let rounded = if rest * 2 >= *x.denom() { whole + 1 } else { whole };
    let digits = format!("{:0>width$}", rounded.to_string(), width = DECIMALS + 1);
    let (int, frac) = digits.split_at(digits.len() - DECIMALS);
    let frac = frac.trim_end_matches('0');
    let sign = if x.numer().is_negative() && rounded != BigInt::from(0) { "-" } else { "" };
    if frac.is_empty() {
        format!("{sign}{int}")
    } else {
        format!("{sign}{int}.{frac}")
    }
}

#[cfg(test)]
mod tests;
