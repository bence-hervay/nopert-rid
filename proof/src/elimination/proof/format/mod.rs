//! The cover file: typed JSON with every field required, unknown fields
//! refused and one canonical spelling.
use crate::arithmetic::polynomial::Exponents;
use crate::elimination::zoom::cover::Parameters;
use crate::elimination::witness::Witness;
use serde::{Deserialize, Serialize};
use std::fmt;

/// The format tag, `"rid-cover/1"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Format {
    #[serde(rename = "rid-cover/1")]
    Version1,
}

/// For which configurations of the cover the statement holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Scope {
    /// No configuration of the cover is a fit: only support gaps are used.
    All,
    /// No configuration of the cover lying in D is a fit: domain
    /// inequalities (and `beyond`) may be used.
    Domain,
}

/// One zoom cover with its proof: its parameters, its witnesses and, per
/// zoom, the cell tree with one leaf entry per cell.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoverFile {
    pub format: Format,
    pub name: String,
    pub scope: Scope,
    pub parameters: Parameters,
    pub witnesses: Vec<Witness>,
    /// One entry per zoom of the cover, in the cover's zoom order.
    pub zooms: Vec<ZoomEntry>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ZoomEntry {
    pub name: String,
    /// The cell tree in preorder (see `zoom::tree`).
    pub tree: String,
    /// One entry per leaf, in preorder.
    pub leaves: Vec<Leaf>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase", deny_unknown_fields)]
pub enum Leaf {
    /// The lemma with witness `index` of the cover and this factor
    /// exponent vector over the zoom variables.
    Witness { index: usize, factor: Exponents },
    /// A cell handed over to the sheared family through the window.
    Delegated,
}

#[derive(Debug)]
pub enum FormatError {
    Json(serde_json::Error),
    /// The bytes are not the canonical spelling of the parsed value.
    NonCanonical,
}

impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FormatError::Json(error) => write!(f, "{error}"),
            FormatError::NonCanonical => write!(f, "the cover file is not in canonical form"),
        }
    }
}

impl std::error::Error for FormatError {}

impl CoverFile {
    /// The canonical bytes: compact JSON and one final newline.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = serde_json::to_vec(self).expect("cover files serialise");
        bytes.push(b'\n');
        bytes
    }

    /// Parses a cover file and accepts it only in canonical form.
    pub fn parse(bytes: &[u8]) -> Result<Self, FormatError> {
        let file: Self = serde_json::from_slice(bytes).map_err(FormatError::Json)?;
        if file.to_bytes() != bytes {
            return Err(FormatError::NonCanonical);
        }
        Ok(file)
    }
}

#[cfg(test)]
mod tests;
