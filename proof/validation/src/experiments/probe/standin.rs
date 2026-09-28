//! A stand-in for the components, for tests only. Its decisions follow
//! simple exact rules, so tests can predict them: Domain is the crate's own
//! Domain component, Exotic and Local are boxes of configurations, and
//! Global accepts a box on which some rotation coordinate stays at least a
//! threshold away from zero. It makes no claim about fits. Faults can be
//! injected to test the experiments' refusals.
use super::Probes;
use rid::arithmetic::exact::{frac, Interval, Q};
use rid::components::domain::Domain;
use rid::components::Component;
use rid::elimination::witness::DomainInequality;
use rid::problem::configuration::{ConfigurationBox, AXES, R};
use rid::search::certificate::RecordData;
use rid::search::BoxError;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "component", deny_unknown_fields)]
pub(crate) enum Data {
    Domain { inequality: usize },
    Exotic { cover: String },
    Local { cover: u32 },
    Global { witness: usize },
}

impl RecordData for Data {}

/// The checkers' form of the stand-in's records.
pub(crate) const FORM: super::RecordForm = super::canonical::<Data>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Fault {
    None,
    /// This component's attempt returns a record for every box, which its
    /// verification then refuses where it does not hold.
    LyingAttempt(usize),
    /// The collection tries the components in reverse order.
    ReversedDecision,
    /// Global also accepts every box whose rotation part contains zero.
    UnsoundGlobal,
}

pub(crate) struct StandIn {
    pub exotic: Vec<(String, ConfigurationBox)>,
    pub local: Vec<(u32, ConfigurationBox)>,
    /// Global accepts a box on which some rotation coordinate keeps at least
    /// this distance from zero.
    pub global: Q,
    pub fault: Fault,
}

/// The box with the given `[lo, hi]` pairs.
pub(crate) fn boxed(ends: [(Q, Q); AXES]) -> ConfigurationBox {
    ConfigurationBox::new(ends.map(|(lo, hi)| Interval::new(lo, hi).unwrap()))
}

impl StandIn {
    /// Covers near the square view and the pentagon event, two overlapping
    /// aligned covers, and Global beyond rotation 1/8.
    pub(crate) fn typical() -> Self {
        // A view rectangle times a rotation cube of half-width `rotation`.
        let cover = |s: (Q, Q), t: (Q, Q), rotation: Q| {
            let r = || (-rotation.clone(), rotation.clone());
            boxed([s, t, r(), r(), r()])
        };
        Self {
            exotic: vec![
                (
                    "square".into(),
                    cover((frac(0, 1), frac(1, 16)), (frac(0, 1), frac(1, 16)), frac(1, 64)),
                ),
                (
                    "pentagon".into(),
                    cover((frac(1, 8), frac(7, 32)), (frac(1, 4), frac(5, 16)), frac(1, 64)),
                ),
            ],
            local: vec![
                (0, cover((frac(1, 32), frac(1, 3)), (frac(0, 1), frac(1, 4)), frac(1, 128))),
                (1, cover((frac(1, 4), frac(2, 3)), (frac(1, 16), frac(2, 5)), frac(1, 128))),
            ],
            global: frac(1, 8),
            fault: Fault::None,
        }
    }

    fn global_axis(&self, b: &ConfigurationBox) -> Option<usize> {
        let far = |axis: usize| {
            let a = &b.axes()[axis];
            *a.lo() >= self.global || *a.hi() <= -&self.global
        };
        let found = R.into_iter().find(|&axis| far(axis));
        if found.is_none() && self.fault == Fault::UnsoundGlobal {
            let zero = Q::zero();
            if R.into_iter().all(|axis| b.axes()[axis].contains(&zero)) {
                return Some(R[0]);
            }
        }
        found
    }

    fn honest(&self, index: usize, b: &ConfigurationBox) -> Option<Data> {
        match index {
            0 => Domain.check(b).map(|inequality| Data::Domain { inequality: inequality.index() }),
            1 => self
                .exotic
                .iter()
                .find(|(_, cover)| cover.contains(b))
                .map(|(name, _)| Data::Exotic { cover: name.clone() }),
            2 => self
                .local
                .iter()
                .find(|(_, cover)| cover.contains(b))
                .map(|(number, _)| Data::Local { cover: *number }),
            3 => self.global_axis(b).map(|witness| Data::Global { witness }),
            _ => panic!("component {index} out of range"),
        }
    }
}

impl Probes for StandIn {
    type Data = Data;

    fn components(&self) -> Vec<String> {
        ["Domain", "Exotic", "Local", "Global"].map(String::from).to_vec()
    }

    fn decide(&self, b: &ConfigurationBox) -> Result<Option<Data>, BoxError> {
        let mut order = vec![0, 1, 2, 3];
        if self.fault == Fault::ReversedDecision {
            order.reverse();
        }
        Ok(order.into_iter().find_map(|index| self.attempt(index, b)))
    }

    fn attempt(&self, index: usize, b: &ConfigurationBox) -> Option<Data> {
        if self.fault == Fault::LyingAttempt(index) {
            return Some(match index {
                0 => Data::Domain { inequality: 0 },
                1 => Data::Exotic { cover: self.exotic[0].0.clone() },
                2 => Data::Local { cover: self.local[0].0 },
                _ => Data::Global { witness: R[0] },
            });
        }
        self.honest(index, b)
    }

    fn covers(&self, index: usize) -> Vec<Data> {
        match index {
            1 => self
                .exotic
                .iter()
                .map(|(name, _)| Data::Exotic { cover: name.clone() })
                .collect(),
            2 => self
                .local
                .iter()
                .map(|(number, _)| Data::Local { cover: *number })
                .collect(),
            _ => Vec::new(),
        }
    }

    fn verify(&self, index: usize, b: &ConfigurationBox, data: &Data) -> Result<(), BoxError> {
        let holds = match (index, data) {
            (0, Data::Domain { inequality }) => {
                Domain.holds(b, &DomainInequality::new(*inequality)?).is_ok()
            }
            (1, Data::Exotic { cover }) => self
                .exotic
                .iter()
                .any(|(name, r)| name == cover && r.contains(b)),
            (2, Data::Local { cover }) => self
                .local
                .iter()
                .any(|(number, r)| number == cover && r.contains(b)),
            (3, Data::Global { witness }) => self.global_axis(b) == Some(*witness),
            _ => return Err(format!("component {index} cannot verify {data:?}").into()),
        };
        if holds {
            Ok(())
        } else {
            Err(format!("{data:?} does not hold for the box").into())
        }
    }
}
