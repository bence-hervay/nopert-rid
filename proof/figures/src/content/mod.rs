//! What is drawn: exact geometric content in world coordinates, with no
//! colours, panels or output format.
use rid::arithmetic::exact::Q;

pub mod graph;
pub mod number;
pub mod scene;
pub mod slice;
pub mod support;

/// An exact point of a world plane: screen coordinates for shadows, the two
/// plot parameters for slices. The second coordinate points up.
pub type Point = [Q; 2];

/// A drawable piece of content in world coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Shape {
    /// A closed segment; it has no interior.
    Segment([Point; 2]),
    /// A closed polygon given by its vertices in order, without repeating
    /// the first. Fewer than three vertices (a degenerate enclosure) are
    /// allowed and draw only a stroke.
    Polygon(Vec<Point>),
    /// A circle, filled as a dot by a style with a fill.
    Circle { centre: Point, radius: Q },
}
