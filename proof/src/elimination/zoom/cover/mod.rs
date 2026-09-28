//! Zoom covers: families of zooms built from pinned parameters, with the
//! set of configurations they cover (a box in cover coordinates), the exact
//! per-box inclusion test, the window hand-over to a sheared family and a
//! constructive coverage map. The coverage lemmas are proved in the article
//! (point, tube and hand-over coverage, Sections 8 to 10).
use super::{reciprocal, Coordinates, Five, LinearForm, Zoom, ZoomCell, ZoomError};
use crate::arithmetic::exact::{Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::{Polynomial, PolynomialError, VARIABLES};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::domain::{self, CONSTRAINT_COUNT};
use std::cmp::Ordering;
use std::fmt;

mod coverage;
mod parameters;

pub use coverage::Location;
pub use parameters::{faces, Beyond, Face, Number, Parameters, Rational, UnitBox, Shape, Window};
use parameters::{small, valid_unit_box};

/// What a zoom of a cover is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Tube { offset: Face },
    A { base: Face, offset: Face },
    B { offset: Face },
    /// A zoom of the sheared family (receives delegations, never delegates).
    Sheared { base: Face, offset: Face },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoverError {
    /// A unit box with an end outside `{−1, 0, 1}` or not around 0.
    UnitBox,
    /// The base and offset dimensions are not both positive with sum 5, or
    /// the shear or radii have the wrong shape.
    Dimensions,
    /// A radius, ratio or base interval that must be positive is not.
    NotPositive,
    /// The map is singular.
    Singular,
    /// A window face is not a base face, is repeated, or the window has no
    /// face.
    WindowFace,
    /// The beyond axis or inequality is out of range, or D's inequality is
    /// not a positive multiple of the distance beyond the cover's end.
    Beyond,
    Zoom { zoom: String, error: ZoomError },
    Polynomial(PolynomialError),
}

impl fmt::Display for CoverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CoverError::UnitBox => write!(f, "invalid unit_box"),
            CoverError::Dimensions => write!(f, "inconsistent dimensions"),
            CoverError::NotPositive => write!(f, "a radius, ratio or width is not positive"),
            CoverError::Singular => write!(f, "the cover map is singular"),
            CoverError::WindowFace => write!(f, "the window faces are not distinct base faces"),
            CoverError::Beyond => write!(f, "the beyond statement does not hold"),
            CoverError::Zoom { zoom, error } => write!(f, "zoom {zoom}: {error}"),
            CoverError::Polynomial(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for CoverError {}

impl From<PolynomialError> for CoverError {
    fn from(error: PolynomialError) -> Self {
        CoverError::Polynomial(error)
    }
}

/// A zoom cover: its parameters, its zooms (covering zooms first, then the
/// sheared family), and the covered set in cover coordinates.
#[derive(Clone, Debug)]
pub struct ZoomCover {
    parameters: Parameters,
    zooms: Vec<Zoom>,
    kinds: Vec<Kind>,
    inverse: [[QSqrt5; VARIABLES]; VARIABLES],
    bounds: Five<[QSqrt5; 2]>,
    /// Cover coordinate `i` as a linear form in the coordinates, `(M⁻¹x)ᵢ`,
    /// and its offset `(M⁻¹c₀)ᵢ`, prepared for `contains`.
    forms: Five<(LinearForm, QSqrt5)>,
}

fn rational(x: &Q) -> QSqrt5 {
    QSqrt5::from_rational(x.clone())
}

fn positive(x: &Q) -> bool {
    x > &Q::zero()
}

/// The exact inverse of a matrix over Q(√5) by Gauss–Jordan elimination.
fn invert(m: &[[QSqrt5; VARIABLES]; VARIABLES]) -> Option<[[QSqrt5; VARIABLES]; VARIABLES]> {
    let n = VARIABLES;
    let mut a: Vec<Vec<QSqrt5>> = m.iter().map(|row| row.to_vec()).collect();
    let mut b: Vec<Vec<QSqrt5>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { QSqrt5::one() } else { QSqrt5::zero() }).collect())
        .collect();
    for col in 0..n {
        let pivot = (col..n).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, pivot);
        b.swap(col, pivot);
        let inv = reciprocal(&a[col][col]);
        for j in 0..n {
            a[col][j] = &a[col][j] * &inv;
            b[col][j] = &b[col][j] * &inv;
        }
        for r in (0..n).filter(|&r| r != col) {
            let factor = a[r][col].clone();
            if factor.is_zero() {
                continue;
            }
            for j in 0..n {
                a[r][j] = &a[r][j] - &(&factor * &a[col][j]);
                b[r][j] = &b[r][j] - &(&factor * &b[col][j]);
            }
        }
    }
    Some(std::array::from_fn(|i| std::array::from_fn(|j| b[i][j].clone())))
}

impl ZoomCover {
    pub fn new(parameters: Parameters) -> Result<Self, CoverError> {
        let map: [[QSqrt5; VARIABLES]; VARIABLES] =
            std::array::from_fn(|i| std::array::from_fn(|j| parameters.map[i][j].0.clone()));
        let inverse = invert(&map).ok_or(CoverError::Singular)?;
        let shape = &parameters.shape;
        let (k, offset) = (shape.base_dimension(), shape.offset());
        // Every family has one zoom per offset face, which the coverage
        // lemmas use: without an offset there would be no zoom to verify.
        if k == 0 || offset.is_empty() || k + offset.len() != VARIABLES {
            return Err(CoverError::Dimensions);
        }
        if !valid_unit_box(offset) {
            return Err(CoverError::UnitBox);
        }
        if !positive(shape.radius()) {
            return Err(CoverError::NotPositive);
        }
        let mut bounds: Five<[QSqrt5; 2]> = std::array::from_fn(|_| [QSqrt5::zero(), QSqrt5::zero()]);
        for (l, &[lo, hi]) in offset.iter().enumerate() {
            bounds[k + l] = [rational(&(shape.radius() * small(lo))), rational(&(shape.radius() * small(hi)))];
        }
        let offset_faces = faces(offset);
        let mut kinds = Vec::new();
        let mut sheared = None;
        match shape {
            Shape::Tube { base, .. } => {
                if base.iter().any(|[lo, hi]| lo.0 >= hi.0) {
                    return Err(CoverError::NotPositive);
                }
                for (j, [lo, hi]) in base.iter().enumerate() {
                    bounds[j] = [rational(&lo.0), rational(&hi.0)];
                }
                kinds.extend(offset_faces.iter().map(|&f| Kind::Tube { offset: f }));
            }
            Shape::Point { base, radii, ratio, window, .. } => {
                if !valid_unit_box(base) {
                    return Err(CoverError::UnitBox);
                }
                let base_faces = faces(base);
                if radii.len() != base_faces.len() {
                    return Err(CoverError::Dimensions);
                }
                if !positive(&ratio.0) || radii.iter().any(|r| !positive(&r.0)) {
                    return Err(CoverError::NotPositive);
                }
                for (face, r) in base_faces.iter().zip(radii) {
                    bounds[face.axis][usize::from(face.side > 0)] = rational(&(&r.0 * small(face.side)));
                }
                for &b in &base_faces {
                    kinds.extend(offset_faces.iter().map(|&f| Kind::A { base: b, offset: f }));
                }
                kinds.extend(offset_faces.iter().map(|&f| Kind::B { offset: f }));
                if let Some(window) = window {
                    if window.shear.len() != offset.len() || window.shear.iter().any(|row| row.len() != k) {
                        return Err(CoverError::Dimensions);
                    }
                    if !positive(&window.radius.0) {
                        return Err(CoverError::NotPositive);
                    }
                    let repeated = window.faces.iter().enumerate().any(|(i, g)| window.faces[..i].contains(g));
                    if window.faces.is_empty() || repeated || window.faces.iter().any(|g| !base_faces.contains(g)) {
                        return Err(CoverError::WindowFace);
                    }
                    for &b in &window.faces {
                        kinds.extend(offset_faces.iter().map(|&f| Kind::Sheared { base: b, offset: f }));
                    }
                    // M' = M S with S = [[I, 0], [−shear, I]].
                    sheared = Some(std::array::from_fn(|i| {
                        std::array::from_fn(|j| {
                            (0..offset.len()).filter(|_| j < k).fold(map[i][j].clone(), |entry, l| {
                                &entry - &(&map[i][k + l] * &rational(&window.shear[l][j].0))
                            })
                        })
                    }));
                }
            }
        }
        let forms = std::array::from_fn(|i| {
            let offset = (0..VARIABLES)
                .fold(QSqrt5::zero(), |sum, j| &sum + &(&inverse[i][j] * &parameters.centre[j].0));
            (parameters.coordinates.linear_form(&inverse[i]), offset)
        });
        let mut cover = Self {
            parameters,
            zooms: Vec::new(),
            kinds: Vec::new(),
            inverse,
            bounds,
            forms,
        };
        for kind in kinds {
            let m = if matches!(kind, Kind::Sheared { .. }) { sheared.as_ref().expect("window") } else { &map };
            let zoom = cover.zoom(m, kind)?;
            cover.zooms.push(zoom);
            cover.kinds.push(kind);
        }
        if let Some(beyond) = cover.parameters.beyond {
            cover.check_beyond(&beyond, &map)?;
        }
        Ok(cover)
    }

    /// The zoom of `kind` with the map `map` (the cover map, or the sheared
    /// map for the sheared family).
    fn zoom(&self, map: &[[QSqrt5; VARIABLES]; VARIABLES], kind: Kind) -> Result<Zoom, CoverError> {
        let shape = &self.parameters.shape;
        let offset = shape.offset();
        let k = shape.base_dimension();
        let mut v = Variables(Vec::new());
        let (z_base, z_offset, factors, name) = match (kind, shape) {
            (Kind::Tube { offset: f }, Shape::Tube { base, .. }) => {
                let b = base
                    .iter()
                    .map(|[lo, hi]| v.take(lo.0.clone(), hi.0.clone()))
                    .collect::<Result<Vec<_>, _>>()?;
                let delta = v.take(Q::zero(), shape.radius().clone())?;
                let unit = v.unit(offset, f)?;
                (b, times(&delta, &unit)?, vec![k], format!("T o{f}"))
            }
            (
                Kind::A { base: g, offset: f } | Kind::Sheared { base: g, offset: f },
                Shape::Point { base, radii, ratio, window, .. },
            ) => {
                let position = faces(base).iter().position(|b| *b == g).expect("a base face");
                let (limit, prefix) = match (kind, window) {
                    (Kind::Sheared { .. }, Some(window)) => (&window.radius.0, "A'"),
                    _ => (&ratio.0, "A"),
                };
                let mu = v.take(Q::zero(), radii[position].0.clone())?;
                let rho = v.take(Q::zero(), limit.clone())?;
                let b = v.unit(base, g)?;
                let unit = v.unit(offset, f)?;
                (times(&mu, &b)?, times(&rho.mul(&mu)?, &unit)?, vec![0, 1], format!("{prefix} b{g} o{f}"))
            }
            (Kind::B { offset: f }, Shape::Point { base, ratio, .. }) => {
                let delta = v.take(Q::zero(), shape.radius().clone())?;
                let b = base
                    .iter()
                    .map(|&[lo, hi]| v.take(small(lo) / &ratio.0, small(hi) / &ratio.0))
                    .collect::<Result<Vec<_>, _>>()?;
                let unit = v.unit(offset, f)?;
                (times(&delta, &b)?, times(&delta, &unit)?, vec![0], format!("B o{f}"))
            }
            _ => unreachable!("zoom kinds follow the shape"),
        };
        debug_assert_eq!(z_base.len(), k);
        let z: Vec<Polynomial> = z_base.into_iter().chain(z_offset).collect();
        let x: Five<Polynomial> = std::array::from_fn(|i| {
            (0..VARIABLES).fold(Polynomial::constant(self.parameters.centre[i].0.clone()), |sum, j| {
                &sum + &z[j].scale(&map[i][j])
            })
        });
        let root: Five<Interval> = v.0.try_into().map_err(|_| CoverError::Dimensions)?;
        let root = ZoomCell::new(root);
        let mut flags = [false; VARIABLES];
        for d in factors {
            flags[d] = true;
        }
        Zoom::new(name.clone(), self.parameters.coordinates, x, root, flags)
            .map_err(|error| CoverError::Zoom { zoom: name, error })
    }

    /// Checks that D's inequality, written in cover coordinates, is
    /// `κ (z_axis − end)` with a constant `κ > 0`, `end` the cover's upper
    /// end on that axis: beyond the end, the inequality is violated.
    fn check_beyond(&self, beyond: &Beyond, map: &[[QSqrt5; VARIABLES]; VARIABLES]) -> Result<(), CoverError> {
        if beyond.axis >= VARIABLES || beyond.inequality >= CONSTRAINT_COUNT {
            return Err(CoverError::Beyond);
        }
        let z = Polynomial::variables();
        let x: Five<Polynomial> = std::array::from_fn(|i| {
            (0..VARIABLES).fold(
                Polynomial::constant(self.parameters.centre[i].0.clone()),
                |sum, j| &sum + &z[j].scale(&map[i][j]),
            )
        });
        let (view, rotation) = self.parameters.coordinates.configuration(&x);
        let c = domain::polynomial(beyond.inequality, &view, &rotation).map_err(|_| CoverError::Beyond)?;
        let mut unit = [0u8; VARIABLES];
        unit[beyond.axis] = 1;
        let kappa = c.coefficient(&unit);
        let end = &self.bounds[beyond.axis][1];
        let expected = &z[beyond.axis].scale(&kappa) - &Polynomial::constant(&kappa * end);
        if kappa.sign() != Ordering::Greater || c != expected {
            return Err(CoverError::Beyond);
        }
        Ok(())
    }

    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    pub fn zooms(&self) -> &[Zoom] {
        &self.zooms
    }
    pub fn kind(&self, zoom: usize) -> Kind {
        self.kinds[zoom]
    }
    /// The covered set in cover coordinates `z` (closed; on the `beyond`
    /// axis the cover also contains everything above the upper end).
    pub fn bounds(&self) -> &Five<[QSqrt5; 2]> {
        &self.bounds
    }

    /// The cover coordinates `z = map⁻¹ (x − centre)` of coordinates `x`.
    pub fn cover_coordinates(&self, x: &Five<QSqrt5>) -> Five<QSqrt5> {
        std::array::from_fn(|i| {
            (0..VARIABLES).fold(QSqrt5::zero(), |sum, j| {
                &sum + &(&self.inverse[i][j] * &(&x[j] - &self.parameters.centre[j].0))
            })
        })
    }

    /// Exactly whether every configuration of `b` lies in the covered set:
    /// for each cover coordinate, its exact range over `b` lies in the
    /// covered set's interval (above the upper end on the `beyond` axis).
    pub fn contains(&self, b: &ConfigurationBox) -> bool {
        // Arc-plane coordinates of the four view corners, shared by the five
        // coordinates (configuration coordinates need none).
        let corners = match self.parameters.coordinates {
            Coordinates::Configuration => None,
            Coordinates::ArcPlane(_) => match self.parameters.coordinates.view_corners(b) {
                Some(corners) => Some(corners),
                None => return false,
            },
        };
        (0..VARIABLES).all(|i| {
            let (form, offset) = &self.forms[i];
            let [lo, hi] = match &corners {
                Some(corners) => form.range_at(corners, b),
                None => form.range(b).expect("configuration coordinates have a range"),
            };
            let open_above = self.parameters.beyond.is_some_and(|beyond| beyond.axis == i);
            &lo - offset >= self.bounds[i][0] && (open_above || &hi - offset <= self.bounds[i][1])
        })
    }
}

/// Zoom variables handed out in order, with their root intervals.
struct Variables(Vec<Interval>);

impl Variables {
    /// The next variable, ranging over `[lo, hi]` in the root.
    fn take(&mut self, lo: Q, hi: Q) -> Result<Polynomial, CoverError> {
        self.0.push(Interval::new(lo, hi).map_err(|_| CoverError::NotPositive)?);
        Ok(Polynomial::variables()[self.0.len() - 1].clone())
    }
    /// The unit vector on `face` of `unit_box`: `±1` on the face coordinate, a
    /// new variable ranging over the unit box's interval elsewhere.
    fn unit(&mut self, unit_box: &UnitBox, face: Face) -> Result<Vec<Polynomial>, CoverError> {
        (0..unit_box.len())
            .map(|j| {
                if j == face.axis {
                    Ok(Polynomial::constant(QSqrt5::integer(i64::from(face.side))))
                } else {
                    self.take(small(unit_box[j][0]), small(unit_box[j][1]))
                }
            })
            .collect()
    }
}

/// `factor · v` for a polynomial vector.
fn times(factor: &Polynomial, v: &[Polynomial]) -> Result<Vec<Polynomial>, CoverError> {
    v.iter().map(|p| factor.mul(p).map_err(CoverError::from)).collect()
}

#[cfg(test)]
mod tests;
