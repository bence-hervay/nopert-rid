//! Placement of world content on the page: panels with one scale for both
//! axes, exact clipping, and inset windows with their true-size markers,
//! frames and connecting lines.
use crate::content::{Point, Shape};
use num_bigint::BigInt;
use rid::arithmetic::exact::{Interval, Q};
use serde::Deserialize;
use std::fmt;

/// Panels narrower than `2^-MIN_WORLD_BITS` world units are refused, so
/// that a drawn configuration point (within `2^-100` of its exact position)
/// is exact to far below the resolution of any panel.
pub const MIN_WORLD_BITS: u32 = 60;

/// A closed rectangle `[x] × [y]`.
pub type Rectangle = [Interval; 2];

/// The page, in millimetres from its upper left corner, `y` pointing down.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Canvas {
    width_mm: Q,
    height_mm: Q,
}

/// A viewport of the page showing a world rectangle at one scale: world
/// `x` maps to page `x`, world `y` (pointing up) to page `y` (pointing down).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Panel {
    viewport: Rectangle,
    world: Rectangle,
    /// Millimetres per world unit, the same on both axes.
    scale: Q,
}

/// A corner as seen on the page.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Corner {
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

/// An inset in page coordinates: the detail's world rectangle marked at true
/// size in the overview, the detail panel's frame, and connecting lines.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Inset {
    pub window: Shape,
    pub frame: Shape,
    pub connectors: Vec<Shape>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LayoutError {
    EmptyCanvas,
    EmptyViewport,
    EmptyWorld,
    OutsideCanvas,
    /// Viewport and world rectangles of different aspect ratios.
    UnequalScale,
    /// A world rectangle narrower than `2^-MIN_WORLD_BITS`.
    TooNarrow,
    /// An inset whose detail world rectangle is not inside the overview's.
    WindowOutsideOverview,
}

impl fmt::Display for LayoutError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            LayoutError::EmptyCanvas => "the canvas has no area",
            LayoutError::EmptyViewport => "the viewport has no area",
            LayoutError::EmptyWorld => "the world rectangle has no area",
            LayoutError::OutsideCanvas => "the viewport is not inside the canvas",
            LayoutError::UnequalScale => {
                "the viewport and world rectangles have different aspect ratios"
            }
            LayoutError::TooNarrow => "the world rectangle is narrower than 2^-60",
            LayoutError::WindowOutsideOverview => {
                "the detail's world rectangle is not inside the overview's"
            }
        };
        f.write_str(text)
    }
}

impl std::error::Error for LayoutError {}

impl Canvas {
    pub fn new(width_mm: Q, height_mm: Q) -> Result<Self, LayoutError> {
        if width_mm <= Q::zero() || height_mm <= Q::zero() {
            return Err(LayoutError::EmptyCanvas);
        }
        Ok(Self {
            width_mm,
            height_mm,
        })
    }
    pub fn width_mm(&self) -> &Q {
        &self.width_mm
    }
    pub fn height_mm(&self) -> &Q {
        &self.height_mm
    }
}

impl Panel {
    pub fn new(
        canvas: &Canvas,
        viewport: Rectangle,
        world: Rectangle,
    ) -> Result<Self, LayoutError> {
        let zero = Q::zero();
        if viewport.iter().any(|i| i.width() == zero) {
            return Err(LayoutError::EmptyViewport);
        }
        if world.iter().any(|i| i.width() == zero) {
            return Err(LayoutError::EmptyWorld);
        }
        let [x, y] = &viewport;
        let on_canvas = *x.lo() >= zero
            && *y.lo() >= zero
            && x.hi() <= canvas.width_mm()
            && y.hi() <= canvas.height_mm();
        if !on_canvas {
            return Err(LayoutError::OutsideCanvas);
        }
        let scale = x.width() / world[0].width();
        if scale != y.width() / world[1].width() {
            return Err(LayoutError::UnequalScale);
        }
        let narrowest = Q::new(BigInt::from(1), BigInt::from(1) << MIN_WORLD_BITS);
        if world.iter().any(|i| i.width() < narrowest) {
            return Err(LayoutError::TooNarrow);
        }
        Ok(Self {
            viewport,
            world,
            scale,
        })
    }

    pub fn viewport(&self) -> &Rectangle {
        &self.viewport
    }

    pub fn world(&self) -> &Rectangle {
        &self.world
    }

    /// The page position of a world point (exact).
    pub fn to_page(&self, p: &Point) -> Point {
        let [x, y] = &self.viewport;
        let [wx, wy] = &self.world;
        [
            x.lo() + &(&(&p[0] - wx.lo()) * &self.scale),
            y.lo() + &(&(wy.hi() - &p[1]) * &self.scale),
        ]
    }

    /// The world shapes in page coordinates, cut exactly to the viewport
    /// grown by `padding_mm` on every side; shapes outside it are dropped.
    /// With a padding above half the stroke width, every cut lies outside
    /// the viewport, so clipping the drawing to the viewport shows exactly
    /// the uncut shapes while keeping every page coordinate bounded.
    pub fn place(&self, shapes: &[Shape], padding_mm: &Q) -> Vec<Shape> {
        let window = self.viewport.clone().map(|i| {
            Interval::new(i.lo() - padding_mm, i.hi() + padding_mm).expect("nonnegative padding")
        });
        shapes
            .iter()
            .filter_map(|shape| match shape {
                Shape::Segment([p, q]) => {
                    clip_segment([self.to_page(p), self.to_page(q)], &window).map(Shape::Segment)
                }
                Shape::Polygon(points) => {
                    let page: Vec<Point> = points.iter().map(|p| self.to_page(p)).collect();
                    let clipped = clip_polygon(page, &window);
                    (!clipped.is_empty()).then_some(Shape::Polygon(clipped))
                }
                // Kept whole when it meets the window, so its page
                // coordinates stay bounded; the viewport clip cuts it.
                Shape::Circle { centre, radius } => {
                    let (centre, radius) = (self.to_page(centre), radius * &self.scale);
                    let meets = (0..2).all(|axis| {
                        &centre[axis] + &radius >= *window[axis].lo() && &centre[axis] - &radius <= *window[axis].hi()
                    });
                    meets.then_some(Shape::Circle { centre, radius })
                }
            })
            .collect()
    }
}

/// The window marker, frame and connecting lines of an inset from `overview`
/// to `detail`. Each connector joins a corner of the marker to a corner of
/// the detail's viewport. A window whose sides are both shorter than
/// `mark_mm` is marked by a square of that side around its centre, so that
/// it stays visible.
pub fn inset(
    overview: &Panel,
    detail: &Panel,
    connectors: &[(Corner, Corner)],
    mark_mm: &Q,
) -> Result<Inset, LayoutError> {
    let inside = overview
        .world
        .iter()
        .zip(&detail.world)
        .all(|(outer, inner)| outer.lo() <= inner.lo() && inner.hi() <= outer.hi());
    if !inside {
        return Err(LayoutError::WindowOutsideOverview);
    }
    let [wx, wy] = &detail.world;
    let [a, b] = [
        overview.to_page(&[wx.lo().clone(), wy.hi().clone()]),
        overview.to_page(&[wx.hi().clone(), wy.lo().clone()]),
    ];
    let mut window: Rectangle = [
        Interval::new(a[0].clone(), b[0].clone()).expect("x keeps its order"),
        Interval::new(a[1].clone(), b[1].clone()).expect("the top maps above the bottom"),
    ];
    if window.iter().all(|side| side.width() < *mark_mm) {
        let half = mark_mm / &Q::from_integer(2.into());
        window = window.map(|side| {
            let centre = side.midpoint();
            Interval::new(&centre - &half, &centre + &half).expect("a nonnegative mark")
        });
    }
    Ok(Inset {
        window: Shape::Polygon(rectangle_polygon(&window)),
        frame: Shape::Polygon(rectangle_polygon(&detail.viewport)),
        connectors: connectors
            .iter()
            .map(|&(from, to)| {
                Shape::Segment([corner(&window, from), corner(&detail.viewport, to)])
            })
            .collect(),
    })
}

/// A page rectangle's corner.
fn corner(r: &Rectangle, c: Corner) -> Point {
    let [x, y] = r;
    let (left, top) = match c {
        Corner::TopLeft => (true, true),
        Corner::TopRight => (false, true),
        Corner::BottomLeft => (true, false),
        Corner::BottomRight => (false, false),
    };
    [
        (if left { x.lo() } else { x.hi() }).clone(),
        (if top { y.lo() } else { y.hi() }).clone(),
    ]
}

fn rectangle_polygon(r: &Rectangle) -> Vec<Point> {
    [Corner::TopLeft, Corner::TopRight, Corner::BottomRight, Corner::BottomLeft]
        .map(|c| corner(r, c))
        .to_vec()
}

/// The part of the closed segment inside the closed rectangle (Liang–Barsky,
/// exact): `p + t(q - p)` for `t` in the intersection of the parameter
/// ranges of the four half-planes.
fn clip_segment([p, q]: [Point; 2], r: &Rectangle) -> Option<[Point; 2]> {
    let d = [&q[0] - &p[0], &q[1] - &p[1]];
    let (mut start, mut end) = (Q::zero(), Q::one());
    for axis in 0..2 {
        // lo ≤ p + t·d  and  p + t·d ≤ hi
        for (direction, room) in [
            (-&d[axis], &p[axis] - r[axis].lo()),
            (d[axis].clone(), r[axis].hi() - &p[axis]),
        ] {
            if direction == Q::zero() {
                if room < Q::zero() {
                    return None;
                }
                continue;
            }
            let t = room / &direction;
            if direction < Q::zero() {
                start = start.max(t);
            } else {
                end = end.min(t);
            }
        }
    }
    if start > end {
        return None;
    }
    let at = |t: &Q| [&p[0] + &(&d[0] * t), &p[1] + &(&d[1] * t)];
    Some([at(&start), at(&end)])
}

/// The polygon cut to the closed rectangle by one half-plane after another
/// (Sutherland–Hodgman, exact). The result can have edges along the
/// rectangle's boundary; it is empty when nothing is inside.
fn clip_polygon(mut points: Vec<Point>, r: &Rectangle) -> Vec<Point> {
    for axis in 0..2 {
        for (bound, below) in [(r[axis].lo(), false), (r[axis].hi(), true)] {
            let inside = |p: &Point| if below { &p[axis] <= bound } else { &p[axis] >= bound };
            let cut = |a: &Point, b: &Point| -> Point {
                let t = (bound - &a[axis]) / (&b[axis] - &a[axis]);
                let mut p = [&a[0] + &(&(&b[0] - &a[0]) * &t), &a[1] + &(&(&b[1] - &a[1]) * &t)];
                p[axis] = bound.clone();
                p
            };
            let mut kept = Vec::with_capacity(points.len() + 2);
            for (k, current) in points.iter().enumerate() {
                let previous = &points[(k + points.len() - 1) % points.len()];
                match (inside(previous), inside(current)) {
                    (true, true) => kept.push(current.clone()),
                    (true, false) => kept.push(cut(previous, current)),
                    (false, true) => {
                        kept.push(cut(previous, current));
                        kept.push(current.clone());
                    }
                    (false, false) => {}
                }
            }
            points = kept;
        }
    }
    points
}

#[cfg(test)]
mod tests;

