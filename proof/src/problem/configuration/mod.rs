//! Closed boxes of configurations `(s, t, r₁, r₂, r₃)`, binary paths and
//! cyclic midpoint bisection from the root box `B₀`.
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::Polynomial;
use crate::problem::geometry::{Point, View};
use std::fmt;

/// Axis index of each configuration coordinate; it is also the index of the
/// polynomial variable standing for that coordinate.
pub const S: usize = 0;
pub const T: usize = 1;
pub const R: [usize; 3] = [2, 3, 4];
pub const AXES: usize = 5;
const _: () = assert!(AXES == crate::arithmetic::polynomial::VARIABLES);

/// The largest depth (path length) of a search box, and so the largest
/// `max_depth` a certificate header may state. [`ConfigurationBox::from_path`]
/// refuses longer paths before any arithmetic: its exact midpoints make the
/// cost grow faster than the square of the length.
pub const MAX_DEPTH: usize = 4096;

/// A closed box `[lo₀, hi₀] × … × [lo₄, hi₄]` of configurations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConfigurationBox {
    axes: [Interval; AXES],
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PathError {
    /// A byte other than `0` or `1`.
    InvalidByte { position: usize, byte: u8 },
    /// A path longer than [`MAX_DEPTH`].
    TooDeep { length: usize },
}

impl fmt::Display for PathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PathError::InvalidByte { position, byte } => {
                write!(f, "path byte {position} is {byte:#04x}, not '0' or '1'")
            }
            PathError::TooDeep { length } => {
                write!(f, "a path of {length} bytes is deeper than {MAX_DEPTH}")
            }
        }
    }
}

impl std::error::Error for PathError {}

/// Whether `path` is a binary path: only the bytes `0` and `1`. The empty
/// path names the root box. The length is not checked here (see
/// [`ConfigurationBox::from_path`]).
pub fn is_path(path: &str) -> bool {
    path.bytes().all(|b| b == b'0' || b == b'1')
}

/// The search order of paths as a sort key: increasing depth, then lexical
/// order within a depth.
pub fn search_order(path: &str) -> (usize, &str) {
    (path.len(), path)
}

/// The axis bisected at depth `depth`: `depth mod 5`.
pub fn split_axis(depth: usize) -> usize {
    depth % AXES
}

/// The configuration coordinates as polynomial variables: the affine view
/// `(s, t, 1)` and the rotation vector `(r₁, r₂, r₃)`.
pub fn coordinates() -> (View, [Polynomial; 3]) {
    let x = Polynomial::variables();
    (
        View::Affine {
            s: x[S].clone(),
            t: x[T].clone(),
        },
        R.map(|j| x[j].clone()),
    )
}

impl ConfigurationBox {
    pub fn new(axes: [Interval; AXES]) -> Self {
        Self { axes }
    }

    /// The box containing exactly one configuration.
    pub fn point(x: &[Q; AXES]) -> Self {
        Self {
            axes: std::array::from_fn(|j| Interval::point(x[j].clone())),
        }
    }

    /// `B₀ = [0, 2/3] × [0, 2/5] × [-2/5, 2/5]³`.
    pub fn root() -> Self {
        let interval = |lo: Q, hi: Q| Interval::new(lo, hi).expect("ordered root interval");
        let rotation = interval(frac(-2, 5), frac(2, 5));
        Self {
            axes: [
                interval(q(0), frac(2, 3)),
                interval(q(0), frac(2, 5)),
                rotation.clone(),
                rotation.clone(),
                rotation,
            ],
        }
    }

    pub fn axes(&self) -> &[Interval; AXES] {
        &self.axes
    }

    /// The search box of a binary path: starting from `B₀`, the bit at depth
    /// `d` keeps the lower (`0`) or upper (`1`) closed half of axis `d mod 5`.
    /// A path longer than [`MAX_DEPTH`] is refused before any arithmetic, and
    /// otherwise the first byte other than `0` or `1`.
    pub fn from_path(path: &str) -> Result<Self, PathError> {
        if path.len() > MAX_DEPTH {
            return Err(PathError::TooDeep { length: path.len() });
        }
        let mut cover = Self::root();
        for (depth, byte) in path.bytes().enumerate() {
            let child = match byte {
                b'0' => 0,
                b'1' => 1,
                _ => {
                    return Err(PathError::InvalidByte {
                        position: depth,
                        byte,
                    })
                }
            };
            let [lower, upper] = cover.split(split_axis(depth));
            cover = if child == 0 { lower } else { upper };
        }
        Ok(cover)
    }

    /// The two closed halves at the exact midpoint of `axis`; they share that face.
    ///
    /// # Panics
    ///
    /// If `axis` is not below [`AXES`] (a programming error: axes come from
    /// [`split_axis`] or the constants above).
    pub fn split(&self, axis: usize) -> [Self; 2] {
        assert!(axis < AXES, "axis {axis} out of range");
        self.axes[axis].bisect().map(|half| {
            let mut b = self.clone();
            b.axes[axis] = half;
            b
        })
    }

    pub fn contains(&self, other: &ConfigurationBox) -> bool {
        self.axes
            .iter()
            .zip(&other.axes)
            .all(|(a, b)| a.lo() <= b.lo() && b.hi() <= a.hi())
    }

    pub fn contains_point(&self, x: &[Q; AXES]) -> bool {
        self.axes.iter().zip(x).all(|(a, x)| a.contains(x))
    }

    pub fn midpoint(&self) -> [Q; AXES] {
        std::array::from_fn(|j| self.axes[j].midpoint())
    }

    /// The view vectors `(s, t, 1)` at the four corners of the box's view
    /// rectangle, `s` varying fastest: the view set of the box.
    pub fn view_corners(&self) -> [Point; 4] {
        std::array::from_fn(|k| {
            let s = if k & 1 == 0 { self.axes[S].lo() } else { self.axes[S].hi() };
            let t = if k & 2 == 0 { self.axes[T].lo() } else { self.axes[T].hi() };
            [
                QSqrt5::from_rational(s.clone()),
                QSqrt5::from_rational(t.clone()),
                QSqrt5::one(),
            ]
        })
    }

    /// The corner taking the upper endpoint on exactly the axes whose bit is set in `mask`.
    ///
    /// # Panics
    ///
    /// If `mask` is not below `2^AXES = 32` (a programming error).
    pub fn corner(&self, mask: u8) -> [Q; AXES] {
        assert!(mask < 1 << AXES, "corner mask {mask} out of range");
        std::array::from_fn(|j| {
            if mask & (1 << j) == 0 {
                self.axes[j].lo().clone()
            } else {
                self.axes[j].hi().clone()
            }
        })
    }
}

#[cfg(test)]
mod tests;
