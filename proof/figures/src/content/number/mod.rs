//! Exact numbers in figure files. A rational is a JSON string in the crate's
//! canonical spelling (`"3"`, `"-7/12"`); JSON numbers are refused, so no
//! value ever passes through floating point. An element `a + b√5` of Q(√5)
//! is the pair `["a", "b"]`; a closed interval is the pair `["lo", "hi"]`.
use rid::arithmetic::exact::{parse_rational, ExactError, Interval, QSqrt5, Q};
use serde::de::Error as _;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// A rational read from and written to its canonical string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rational(pub Q);

impl Serialize for Rational {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.0.to_string())
    }
}

impl<'de> Deserialize<'de> for Rational {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let text = String::deserialize(deserializer)?;
        parse_rational(&text).map(Rational).map_err(D::Error::custom)
    }
}

/// `[a, b]` as `a + b√5`.
pub fn field(pair: &[Rational; 2]) -> QSqrt5 {
    QSqrt5::new(pair[0].0.clone(), pair[1].0.clone())
}

/// `[lo, hi]` as a closed interval; refused when `lo > hi`.
pub fn interval(pair: &[Rational; 2]) -> Result<Interval, ExactError> {
    Interval::new(pair[0].0.clone(), pair[1].0.clone())
}

/// Makes an `Option` field required: serde would otherwise read a missing
/// field as `None`. Only an explicit `null` gives `None`.
pub fn required<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(deserializer)
}

#[cfg(test)]
mod tests;
