//! Fixtures shared by the components' tests: pseudo-random numbers and
//! paths, boxes as the search meets them, and stand-in covers.
use crate::arithmetic::exact::Q;
use crate::components::collection::Collection;
use crate::components::covered::{Exotic, Local};
use crate::components::record::RecordData;
use crate::elimination::proof::ProvedSet;
use crate::problem::configuration::ConfigurationBox;

/// Deterministic pseudo-random numbers (SplitMix64).
pub struct Random(pub u64);

impl Random {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    /// A path of `depth` random bits.
    pub fn path(&mut self, depth: usize) -> String {
        (0..depth).map(|_| if self.below(2) == 0 { '0' } else { '1' }).collect()
    }
}

/// `2^-k` exactly.
pub fn two_to_minus(k: u32) -> Q {
    Q::new(num_bigint::BigInt::from(1), num_bigint::BigInt::from(1) << k)
}

pub fn from_path(path: &str) -> ConfigurationBox {
    ConfigurationBox::from_path(path).unwrap()
}

/// The collection with Domain and Global only.
pub fn bare() -> Collection {
    Collection::new(Exotic::new(vec![]).unwrap(), Local::new(vec![]).unwrap())
}

/// Subtree roots whose search by Domain and Global alone finishes after a few
/// dozen boxes, with records of both.
pub const SMALL_ROOTS: [&str; 3] = ["01110101000", "0100110011011", "10010001110"];

/// The search of `root` by `collection` in search order (depth, then path),
/// as `(path, record)` pairs, with the refused paths; `None` if more than
/// `budget` boxes are evaluated.
#[allow(clippy::type_complexity)]
pub fn search(
    collection: &Collection,
    root: &str,
    budget: usize,
) -> Option<(Vec<(String, RecordData)>, Vec<String>)> {
    let (mut records, mut refused) = (Vec::new(), Vec::new());
    let mut layer = vec![root.to_owned()];
    while !layer.is_empty() {
        let mut next = Vec::new();
        for path in layer {
            if records.len() + refused.len() >= budget {
                return None;
            }
            match collection.check(&from_path(&path)) {
                Some(data) => records.push((path, data)),
                None => {
                    next.push(format!("{path}0"));
                    next.push(format!("{path}1"));
                    refused.push(path);
                }
            }
        }
        layer = next;
    }
    Some((records, refused))
}

/// A cover given as a box, with no proof at all: for tests of the
/// components' logic, never for a certificate.
pub struct BoxCover {
    pub name: String,
    pub cover: ConfigurationBox,
}

impl BoxCover {
    pub fn boxed(name: &str, path: &str) -> Box<dyn ProvedSet> {
        Box::new(BoxCover {
            name: name.to_owned(),
            cover: from_path(path),
        })
    }
}

impl ProvedSet for BoxCover {
    fn name(&self) -> &str {
        &self.name
    }
    fn contains(&self, b: &ConfigurationBox) -> bool {
        self.cover.contains(b)
    }
}
