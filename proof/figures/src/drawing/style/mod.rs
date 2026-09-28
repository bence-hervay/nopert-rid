//! Explicit styles: colours, strokes and opacity, and one style for every
//! outcome of a slice cell. Nothing has a default.
use crate::content::number::{required, Rational};
use crate::content::slice::{Category, Outcome};
use rid::components::record::ExoticCover;
use std::collections::BTreeMap;
use rid::arithmetic::exact::Q;
use serde::de::Error as _;
use serde::{Deserialize, Deserializer};
use std::fmt;

/// An sRGB colour, written `#rrggbb`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color([u8; 3]);

/// How lines are drawn. Caps and joins are always round, so a stroke
/// reaches at most half its width beyond the drawn geometry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stroke {
    color: Color,
    width_mm: Q,
    /// Dashes and gaps of this length, or a solid line.
    dash_mm: Option<Q>,
}

/// How one layer is drawn. The opacity applies to the layer as a whole, so
/// overlapping shapes of one layer do not darken each other.
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
#[serde(try_from = "Spec")]
pub struct Style {
    fill: Option<Color>,
    stroke: Option<Stroke>,
    opacity: Q,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StyleError {
    Color(String),
    /// Opacity outside `(0, 1]`.
    Opacity(Q),
    /// A stroke width that is not positive.
    Width(Q),
    /// A dash length that is not positive.
    Dash(Q),
    /// Neither a fill nor a stroke: the layer would be invisible.
    Invisible,
}

impl fmt::Display for StyleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StyleError::Color(text) => write!(f, "colour {text:?} is not #rrggbb"),
            StyleError::Opacity(o) => write!(f, "opacity {o} is not in (0, 1]"),
            StyleError::Width(w) => write!(f, "stroke width {w} mm is not positive"),
            StyleError::Dash(d) => write!(f, "dash length {d} mm is not positive"),
            StyleError::Invisible => write!(f, "a style needs a fill or a stroke"),
        }
    }
}

impl std::error::Error for StyleError {}

impl Color {
    pub fn parse(text: &str) -> Result<Self, StyleError> {
        let error = || StyleError::Color(text.to_owned());
        let digits = text.strip_prefix('#').ok_or_else(error)?;
        if digits.len() != 6 || !digits.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(error());
        }
        let channel =
            |k: usize| u8::from_str_radix(&digits[2 * k..2 * k + 2], 16).map_err(|_| error());
        Ok(Self([channel(0)?, channel(1)?, channel(2)?]))
    }

    /// `#rrggbb` in lower case.
    pub fn hex(&self) -> String {
        let [r, g, b] = self.0;
        format!("#{r:02x}{g:02x}{b:02x}")
    }
}

impl<'de> Deserialize<'de> for Color {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        Color::parse(&String::deserialize(deserializer)?).map_err(D::Error::custom)
    }
}

impl Stroke {
    pub fn new(color: Color, width_mm: Q) -> Result<Self, StyleError> {
        if width_mm <= Q::zero() {
            return Err(StyleError::Width(width_mm));
        }
        Ok(Self { color, width_mm, dash_mm: None })
    }
    /// The same stroke dashed, with dashes and gaps `dash_mm` long.
    pub fn dashed(self, dash_mm: Q) -> Result<Self, StyleError> {
        if dash_mm <= Q::zero() {
            return Err(StyleError::Dash(dash_mm));
        }
        Ok(Self { dash_mm: Some(dash_mm), ..self })
    }
    pub fn color(&self) -> Color {
        self.color
    }
    pub fn width_mm(&self) -> &Q {
        &self.width_mm
    }
    pub fn dash_mm(&self) -> Option<&Q> {
        self.dash_mm.as_ref()
    }
}

impl Style {
    pub fn new(
        fill: Option<Color>,
        stroke: Option<Stroke>,
        opacity: Q,
    ) -> Result<Self, StyleError> {
        if opacity <= Q::zero() || opacity > Q::one() {
            return Err(StyleError::Opacity(opacity));
        }
        if fill.is_none() && stroke.is_none() {
            return Err(StyleError::Invisible);
        }
        Ok(Self {
            fill,
            stroke,
            opacity,
        })
    }
    pub fn fill(&self) -> Option<Color> {
        self.fill
    }
    pub fn stroke(&self) -> Option<&Stroke> {
        self.stroke.as_ref()
    }
    pub fn opacity(&self) -> &Q {
        &self.opacity
    }
}

/// The JSON form of a style; every field is required, `null` meaning none.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Spec {
    #[serde(deserialize_with = "required")]
    fill: Option<Color>,
    #[serde(deserialize_with = "required")]
    stroke: Option<StrokeSpec>,
    opacity: Rational,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StrokeSpec {
    color: Color,
    width_mm: Rational,
    /// Optional: a solid line when absent.
    #[serde(default)]
    dash_mm: Option<Rational>,
}

impl TryFrom<Spec> for Style {
    type Error = StyleError;
    fn try_from(spec: Spec) -> Result<Self, StyleError> {
        let stroke = match spec.stroke {
            Some(s) => {
                let stroke = Stroke::new(s.color, s.width_mm.0)?;
                Some(match s.dash_mm {
                    Some(d) => stroke.dashed(d.0)?,
                    None => stroke,
                })
            }
            None => None,
        };
        Style::new(spec.fill, stroke, spec.opacity.0)
    }
}

/// One value for each outcome of a slice cell: every component, every
/// Exotic cover, unresolved and outside cells. All are required.
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PerOutcome<T> {
    #[serde(rename = "Domain")]
    pub domain: T,
    #[serde(rename = "Global")]
    pub global: T,
    #[serde(rename = "Local")]
    pub local: T,
    #[serde(rename = "Exotic")]
    pub exotic: PerCover<T>,
    pub unresolved: T,
    pub outside: T,
}

/// One value for every Exotic cover, keyed by the cover's name; every
/// cover is required, unknown names and repeated names are refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PerCover<T>(BTreeMap<ExoticCover, T>);

impl<T> PerCover<T> {
    pub fn get(&self, cover: ExoticCover) -> &T {
        &self.0[&cover]
    }
}

impl<'de, T: Deserialize<'de>> Deserialize<'de> for PerCover<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Entries<T>(std::marker::PhantomData<T>);
        impl<'de, T: Deserialize<'de>> serde::de::Visitor<'de> for Entries<T> {
            type Value = BTreeMap<ExoticCover, T>;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "a map from every Exotic cover")
            }
            fn visit_map<A: serde::de::MapAccess<'de>>(self, mut access: A) -> Result<Self::Value, A::Error> {
                let mut map = BTreeMap::new();
                while let Some((cover, value)) = access.next_entry::<ExoticCover, T>()? {
                    if map.insert(cover, value).is_some() {
                        return Err(A::Error::custom(format!("duplicate Exotic cover {cover}")));
                    }
                }
                Ok(map)
            }
        }
        let map = deserializer.deserialize_map(Entries(std::marker::PhantomData))?;
        if let Some(missing) = ExoticCover::ALL.into_iter().find(|c| !map.contains_key(c)) {
            return Err(D::Error::custom(format!("missing Exotic cover {missing}")));
        }
        Ok(Self(map))
    }
}

impl<T> PerOutcome<T> {
    pub fn get(&self, outcome: &Outcome) -> &T {
        match outcome {
            Outcome::Eliminated(Category::Domain) => &self.domain,
            Outcome::Eliminated(Category::Global) => &self.global,
            Outcome::Eliminated(Category::Local) => &self.local,
            Outcome::Eliminated(Category::Exotic { cover }) => self.exotic.get(*cover),
            Outcome::Unresolved => &self.unresolved,
            Outcome::Outside => &self.outside,
        }
    }
}

/// Every outcome, in the order in which cell layers draw them: the
/// components, the Exotic covers in their order, unresolved, outside.
pub fn outcomes() -> Vec<Outcome> {
    let components = [Category::Domain, Category::Global, Category::Local];
    let covers = ExoticCover::ALL.map(|cover| Category::Exotic { cover });
    components
        .into_iter()
        .chain(covers)
        .map(Outcome::Eliminated)
        .chain([Outcome::Unresolved, Outcome::Outside])
        .collect()
}

#[cfg(test)]
mod tests;
