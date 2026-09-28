//! The constructive coverage map and the window hand-over: the proofs of the
//! coverage lemmas and of the hand-over lemma as computations.
use super::parameters::{small, Face, UnitBox, Shape};
use super::{ZoomCover, Kind};
use crate::arithmetic::exact::{Interval, QSqrt5};
use crate::arithmetic::polynomial::VARIABLES;
use crate::elimination::zoom::{reciprocal, Five, ZoomCell};
use std::cmp::Ordering;

/// Where a point of a cover lies.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Location {
    /// In the no-fit set.
    NoFit,
    /// Beyond the cover's end on the `beyond` axis: outside D.
    OutsideDomain,
    /// The image of the zoom variables `y` of zoom `zoom`.
    Zoom { zoom: usize, y: Five<QSqrt5> },
}

/// `max_j |z_j|` and the first face attaining it, for `z` in the cone of the
/// unit_box; `None` outside the cone. The face is `None` when `z = 0`.
fn gauge(z: &[QSqrt5], unit_box: &UnitBox) -> Option<(QSqrt5, Option<Face>)> {
    let mut best = (QSqrt5::zero(), None);
    for (axis, (x, &[lo, hi])) in z.iter().zip(unit_box).enumerate() {
        let side = match x.sign() {
            Ordering::Greater if hi == 1 => 1,
            Ordering::Less if lo == -1 => -1,
            Ordering::Equal => continue,
            _ => return None,
        };
        let size = if side > 0 { x.clone() } else { -x };
        if size > best.0 {
            best = (size, Some(Face { axis, side }));
        }
    }
    Some(best)
}

/// The coordinates of `v` other than the face coordinate.
fn free(v: &[QSqrt5], face: Face) -> impl Iterator<Item = QSqrt5> + '_ {
    v.iter().enumerate().filter(move |(j, _)| *j != face.axis).map(|(_, x)| x.clone())
}

/// The unit vector of length `len` on `face`, its free coordinates taken in
/// order from `values`.
fn unit<T>(len: usize, face: Face, values: &mut impl Iterator<Item = T>, side: impl Fn(i8) -> T) -> Vec<T> {
    (0..len)
        .map(|j| if j == face.axis { side(face.side) } else { values.next().expect("a free variable") })
        .collect()
}

impl ZoomCover {
    /// For coordinates `x` in the cover: a zoom and zoom variables in its
    /// root with image `x`, or the no-fit set, or outside D beyond the
    /// `beyond` face. `None` when `x` is not in the cover.
    pub fn locate(&self, x: &Five<QSqrt5>) -> Option<Location> {
        let z = self.cover_coordinates(x);
        let bounds = self.bounds();
        let inside = |i: usize| bounds[i][0] <= z[i] && z[i] <= bounds[i][1];
        if let Some(i) = (0..VARIABLES).find(|&i| !inside(i)) {
            let beyond = self.parameters().beyond.is_some_and(|b| b.axis == i) && z[i] > bounds[i][1];
            let rest = (0..VARIABLES).filter(|&j| j != i).all(inside);
            return (beyond && rest).then_some(Location::OutsideDomain);
        }
        let shape = &self.parameters().shape;
        let k = shape.base_dimension();
        let (delta, offset_face) = gauge(&z[k..], shape.offset())?;
        let Some(offset_face) = offset_face else {
            return Some(Location::NoFit);
        };
        let scaled = |v: &[QSqrt5], by: &QSqrt5| -> Vec<QSqrt5> { v.iter().map(|x| x * &reciprocal(by)).collect() };
        let unit_offset = scaled(&z[k..], &delta);
        let mut y: Vec<QSqrt5> = Vec::new();
        let kind = match shape {
            Shape::Tube { .. } => {
                y.extend(z[..k].iter().cloned());
                y.push(delta);
                y.extend(free(&unit_offset, offset_face));
                Kind::Tube { offset: offset_face }
            }
            Shape::Point { base, ratio, .. } => {
                let (mu, base_face) = gauge(&z[..k], base)?;
                if delta <= &mu * &QSqrt5::from_rational(ratio.0.clone()) {
                    let base_face = base_face.expect("μ > 0 when δ > 0");
                    y.push(mu.clone());
                    y.push(&delta * &reciprocal(&mu));
                    y.extend(free(&scaled(&z[..k], &mu), base_face));
                    y.extend(free(&unit_offset, offset_face));
                    Kind::A { base: base_face, offset: offset_face }
                } else {
                    y.push(delta.clone());
                    y.extend(scaled(&z[..k], &delta));
                    y.extend(free(&unit_offset, offset_face));
                    Kind::B { offset: offset_face }
                }
            }
        };
        let zoom = (0..self.zooms().len()).find(|&c| self.kind(c) == kind)?;
        Some(Location::Zoom { zoom, y: y.try_into().ok()? })
    }

    /// The window vector `w = ρf̂ + shear · b̂` of an A zoom at `y`, with the
    /// zoom's base face.
    pub(super) fn window_vector(&self, zoom: usize, y: &Five<QSqrt5>) -> Option<(Face, Vec<QSqrt5>)> {
        let Kind::A { base: g, offset: f } = self.kind(zoom) else { return None };
        let Shape::Point { base, offset, window: Some(window), .. } = &self.parameters().shape else {
            return None;
        };
        let integer = |x: i8| QSqrt5::integer(i64::from(x));
        let mut values = y[2..].iter().cloned();
        let b = unit(base.len(), g, &mut values, integer);
        let f = unit(offset.len(), f, &mut values, integer);
        let w = (0..offset.len())
            .map(|l| {
                (0..base.len()).fold(&y[1] * &f[l], |sum, j| {
                    &sum + &(&b[j] * &QSqrt5::from_rational(window.shear[l][j].0.clone()))
                })
            })
            .collect();
        Some((g, w))
    }

    /// For zoom variables `y` of an A zoom whose window vector lies in
    /// `radius · offset`: the same configuration in the sheared family, or
    /// its no-fit set (window lemma).
    pub fn hand_over(&self, zoom: usize, y: &Five<QSqrt5>) -> Option<Location> {
        let (g, w) = self.window_vector(zoom, y)?;
        let Shape::Point { offset, window: Some(window), .. } = &self.parameters().shape else {
            return None;
        };
        let (size, face) = gauge(&w, offset)?;
        if size > QSqrt5::from_rational(window.radius.0.clone()) {
            return None;
        }
        let Some(face) = face else {
            return Some(Location::NoFit);
        };
        let target = (0..self.zooms().len()).find(|&c| self.kind(c) == Kind::Sheared { base: g, offset: face })?;
        let k = VARIABLES - offset.len();
        let mut image: Vec<QSqrt5> = vec![y[0].clone(), size.clone()];
        image.extend(y[2..1 + k].iter().cloned());
        image.extend(free(&w, face).map(|x| &x * &reciprocal(&size)));
        Some(Location::Zoom { zoom: target, y: image.try_into().ok()? })
    }

    /// Exactly whether the window vector of the A zoom `zoom` lies in
    /// `radius · offset` throughout `cell` (only for window faces): each
    /// coordinate is a product of two independent zoom variables plus terms
    /// in further independent variables, so the sum of the exact interval
    /// ranges is its exact range.
    pub fn window_holds(&self, zoom: usize, cell: &ZoomCell) -> bool {
        let Kind::A { base: g, offset: f } = self.kind(zoom) else { return false };
        let Shape::Point { base, offset, window: Some(window), .. } = &self.parameters().shape else {
            return false;
        };
        if !window.faces.contains(&g) {
            return false;
        }
        let point = |x: i8| Interval::point(small(x));
        let mut values = cell[2..].iter().cloned();
        let b = unit(base.len(), g, &mut values, point);
        let f = unit(offset.len(), f, &mut values, point);
        (0..offset.len()).all(|l| {
            let range = (0..base.len()).fold(&cell[1] * &f[l], |sum, j| {
                &sum + &(&Interval::point(window.shear[l][j].0.clone()) * &b[j])
            });
            let [lo, hi] = offset[l];
            range.lo() >= &(small(lo) * &window.radius.0) && range.hi() <= &(small(hi) * &window.radius.0)
        })
    }
}
