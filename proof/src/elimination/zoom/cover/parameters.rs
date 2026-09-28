//! The pinned parameters of a zoom cover, in their canonical JSON form.
use crate::arithmetic::exact::{parse_rational, q, ExactError, QSqrt5, Q};
use crate::arithmetic::polynomial::VARIABLES;
use crate::elimination::zoom::Coordinates;
use serde::{Deserialize, Deserializer, Serialize};
use std::fmt;

/// A number of Q(√5), `a + b√5`, written `["a", "b"]` with canonical rationals.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "[String; 2]", into = "[String; 2]")]
pub struct Number(pub QSqrt5);

/// A rational number written as its canonical string.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct Rational(pub Q);

impl TryFrom<[String; 2]> for Number {
    type Error = ExactError;
    fn try_from([a, b]: [String; 2]) -> Result<Self, ExactError> {
        Ok(Number(QSqrt5::new(parse_rational(&a)?, parse_rational(&b)?)))
    }
}
impl From<Number> for [String; 2] {
    fn from(Number(x): Number) -> Self {
        [x.rational_part().to_string(), x.sqrt5_part().to_string()]
    }
}
impl TryFrom<String> for Rational {
    type Error = ExactError;
    fn try_from(text: String) -> Result<Self, ExactError> {
        parse_rational(&text).map(Rational)
    }
}
impl From<Rational> for String {
    fn from(Rational(x): Rational) -> Self {
        x.to_string()
    }
}

/// A unit box: per axis `[lo, hi]` with `lo ∈ {−1, 0}`, `hi ∈ {0, 1}`, `lo < hi`.
pub type UnitBox = Vec<[i8; 2]>;

/// The face `{z_axis = side}` of a unit box (`side` is `−1` or `+1`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Face {
    pub axis: usize,
    pub side: i8,
}

impl fmt::Display for Face {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.axis, if self.side > 0 { '+' } else { '-' })
    }
}

/// The pinned parameters of a cover: coordinates `x = centre + map · z`
/// with cover coordinates `z = (z_base, z_offset)` and the shape of the family.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    pub coordinates: Coordinates,
    pub centre: [Number; VARIABLES],
    pub map: [[Number; VARIABLES]; VARIABLES],
    pub shape: Shape,
    /// Beyond the upper end of this cover coordinate, an inequality of D is
    /// violated, so the cover extends there (on D nothing is added).
    #[serde(deserialize_with = "required")]
    pub beyond: Option<Beyond>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase", deny_unknown_fields)]
pub enum Shape {
    /// A tube over the base box: `z_base ∈ base`, `z_offset ∈ radius · offset`.
    Tube {
        base: Vec<[Rational; 2]>,
        offset: UnitBox,
        radius: Rational,
    },
    /// A point blow-up at `z = 0`: `z_base` in the base unit_box with one radius
    /// per base face (in face order), `z_offset ∈ radius · offset`, zooms A for
    /// offset gauge at most `ratio` times the base gauge and B beyond.
    Point {
        base: UnitBox,
        offset: UnitBox,
        radii: Vec<Rational>,
        radius: Rational,
        ratio: Rational,
        #[serde(deserialize_with = "required")]
        window: Option<Window>,
    },
}

impl Shape {
    pub fn base_dimension(&self) -> usize {
        match self {
            Shape::Tube { base, .. } => base.len(),
            Shape::Point { base, .. } => base.len(),
        }
    }
    pub fn offset(&self) -> &UnitBox {
        match self {
            Shape::Tube { offset, .. } | Shape::Point { offset, .. } => offset,
        }
    }
    /// The offset radius `ε`.
    pub fn radius(&self) -> &Q {
        match self {
            Shape::Tube { radius, .. } | Shape::Point { radius, .. } => &radius.0,
        }
    }
}

/// The sheared family `z_offset = z'_offset − shear · z_base` on the given base
/// faces with offset ratio at most `radius`, which receives delegated cells of
/// the unsheared A zooms.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Window {
    pub faces: Vec<Face>,
    /// One row per offset coordinate, one entry per base coordinate.
    pub shear: Vec<Vec<Rational>>,
    pub radius: Rational,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Beyond {
    pub axis: usize,
    pub inequality: usize,
}

/// Makes an `Option` field required: only an explicit `null` is `None`.
fn required<'de, D: Deserializer<'de>, T: Deserialize<'de>>(d: D) -> Result<T, D::Error> {
    T::deserialize(d)
}

/// The faces of a unit box: per axis, the lower then the upper face if present.
pub fn faces(unit_box: &UnitBox) -> Vec<Face> {
    let mut out = Vec::new();
    for (axis, ends) in unit_box.iter().enumerate() {
        for side in *ends {
            if side != 0 {
                out.push(Face { axis, side });
            }
        }
    }
    out
}

pub(super) fn valid_unit_box(unit_box: &UnitBox) -> bool {
    unit_box.iter().all(|&[lo, hi]| (lo == -1 || lo == 0) && (hi == 0 || hi == 1) && lo < hi)
}

/// A small integer as a rational.
pub(super) fn small(x: i8) -> Q {
    q(i64::from(x))
}
