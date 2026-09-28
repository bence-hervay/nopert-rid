//! The zoom covers of the crate: the parameters pinned in code for every
//! cover and its data file, and loading them all.
use super::format::Scope;
use super::{Pin, ProofError, ProvedCover};
use crate::arithmetic::exact::{frac, parse_rational, q, QSqrt5, Q};
use crate::elimination::zoom::cover::{Beyond, Face, Number, Parameters, Rational, Shape, UnitBox, Window};
use crate::elimination::zoom::{Coordinates, Sign};
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::AtomicBool;

/// The Exotic covers, in the order they are tried. A cover along the
/// arcs is written in the arc-plane coordinates of one sign `σ = ±1`; its
/// name ends in that sign.
pub const EXOTIC: [&str; 8] =
    ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"];

/// The number of Local covers; they are numbered from 0.
pub const LOCAL_COUNT: usize = 30;

/// How far a Local cover's view rectangle reaches beyond its row of the
/// table on every side (cut to the viewing rectangle `[0, 2/3] × [0, 2/5]`),
/// so that neighbouring covers overlap.
pub const MARGIN: (i64, i64) = (1, 1024);

/// The Local covers: rows `[s_lo, s_hi, t_lo, t_hi, radius]` with disjoint
/// interiors. Together with the square and pentagon covers, the view
/// rectangles (the rows enlarged by [`MARGIN`]) cover D's view triangle.
/// Cover `i` is the tube `[s_lo − m, s_hi + m] × [t_lo − m, t_hi + m]
/// × [−radius, radius]³` over the aligned set, with `m` = [`MARGIN`] and the
/// view rectangle cut to `[0, 2/3] × [0, 2/5]`.
const LOCAL: [[&str; 5]; LOCAL_COUNT] = [
    ["1/12", "1/6", "0", "1/10", "1/96"],
    ["0", "1/6", "1/10", "1/5", "1/80"],
    ["1/6", "1/3", "0", "1/10", "1/50"],
    ["1/6", "1/4", "1/10", "1/5", "1/50"],
    ["1/4", "1/3", "1/10", "3/20", "1/50"],
    ["1/4", "1/3", "3/20", "1/5", "1/50"],
    ["0", "1/12", "1/5", "3/10", "1/50"],
    ["1/12", "1/8", "1/5", "9/40", "1/50"],
    ["1/12", "1/8", "9/40", "1/4", "1/50"],
    ["1/8", "1/6", "1/5", "1/4", "1/200"],
    ["1/12", "1/8", "1/4", "3/10", "1/50"],
    ["1/8", "7/48", "1/4", "11/40", "1/50"],
    ["7/48", "1/6", "1/4", "21/80", "1/200"],
    ["7/48", "5/32", "21/80", "11/40", "1/50"],
    ["5/32", "31/192", "21/80", "43/160", "1/200"],
    ["1/8", "1/6", "11/40", "3/10", "1/50"],
    ["0", "1/6", "3/10", "2/5", "1/50"],
    ["1/6", "5/24", "1/5", "9/40", "1/50"],
    ["1/6", "3/16", "9/40", "1/4", "1/200"],
    ["3/16", "5/24", "9/40", "1/4", "1/50"],
    ["5/24", "1/4", "1/5", "9/40", "1/50"],
    ["5/24", "1/4", "9/40", "1/4", "1/200"],
    ["1/6", "3/16", "1/4", "21/80", "1/200"],
    ["3/16", "5/24", "1/4", "21/80", "1/200"],
    ["5/24", "11/48", "1/4", "11/40", "1/200"],
    ["1/4", "1/3", "1/5", "1/4", "1/50"],
    ["1/3", "1/2", "0", "1/10", "1/50"],
    ["1/3", "1/2", "1/10", "1/5", "1/50"],
    ["1/2", "7/12", "0", "1/10", "1/50"],
    ["7/12", "5/8", "0", "1/20", "1/20"],
];

macro_rules! data {
    ($($path:literal),* $(,)?) => {
        [$(include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/", $path)) as &[u8]),*]
    };
}

const EXOTIC_DATA: [&[u8]; 8] = data!(
    "exotic/square.json",
    "exotic/pentagon.json",
    "exotic/arc+.json",
    "exotic/arc-.json",
    "exotic/endpoint+.json",
    "exotic/endpoint-.json",
    "exotic/crossing+.json",
    "exotic/crossing-.json",
);

const LOCAL_DATA: [&[u8]; LOCAL_COUNT] = data!(
    "local/0.json",
    "local/1.json",
    "local/2.json",
    "local/3.json",
    "local/4.json",
    "local/5.json",
    "local/6.json",
    "local/7.json",
    "local/8.json",
    "local/9.json",
    "local/10.json",
    "local/11.json",
    "local/12.json",
    "local/13.json",
    "local/14.json",
    "local/15.json",
    "local/16.json",
    "local/17.json",
    "local/18.json",
    "local/19.json",
    "local/20.json",
    "local/21.json",
    "local/22.json",
    "local/23.json",
    "local/24.json",
    "local/25.json",
    "local/26.json",
    "local/27.json",
    "local/28.json",
    "local/29.json",
);

fn number(a: i64, b: i64, d: i64) -> Number {
    Number(QSqrt5::new(frac(a, d), frac(b, d)))
}

fn rational(n: i64, d: i64) -> Rational {
    Rational(frac(n, d))
}

/// The matrix `diag(1, w₁, …, w₄)`.
fn diagonal(w: [Q; 4]) -> [[Number; 5]; 5] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| match (i == j, i) {
            (true, 0) => number(1, 0, 1),
            (true, _) => Number(QSqrt5::from_rational(w[i - 1].clone())),
            _ => number(0, 0, 1),
        })
    })
}

fn zero_centre() -> [Number; 5] {
    std::array::from_fn(|_| number(0, 0, 1))
}

fn whole() -> UnitBox {
    vec![[-1, 1]; 3]
}

/// The offset unit box `(v, θ, ζ₁, ζ₃) ∈ [0, 1] × [−1, 1]³` of the arc covers.
fn arc_offset() -> UnitBox {
    vec![[0, 1], [-1, 1], [-1, 1], [-1, 1]]
}

fn interval(lo: Q, hi: Q) -> [Rational; 2] {
    [Rational(lo), Rational(hi)]
}

/// The square view: a point blow-up at `(s, t, r) = 0` over the quadrant
/// `s, t ≥ 0`, radii `1/8`, ratio 1.
fn square() -> Parameters {
    Parameters {
        coordinates: Coordinates::Configuration,
        centre: zero_centre(),
        map: diagonal([q(1), q(1), q(1), q(1)]),
        shape: Shape::Point {
            base: vec![[0, 1], [0, 1]],
            offset: whole(),
            radii: vec![rational(1, 8), rational(1, 8)],
            radius: rational(1, 8),
            ratio: rational(1, 1),
            window: None,
        },
        beyond: None,
    }
}

/// The pentagon event: a point blow-up at the pentagon view in the event
/// coordinates `(η, ξ) = (s − t/φ, s + φt − 1/φ)`, whose inverse is
/// `(s, t) = (s_P, t_P) + [[φ, 1/φ], [−1, 1]] (η, ξ) / √5`; base unit box
/// `[−1, 1] × [−1, 0]`, radii `1/32`, offset radius `1/16`, ratio 2. Beyond
/// `ξ = 0`, D's inequality 0 (`φξ ≤ 0`) is violated.
fn pentagon() -> Parameters {
    let mut map = diagonal([q(1), q(1), q(1), q(1)]);
    map[0][0] = number(5, 1, 10);
    map[0][1] = number(5, -1, 10);
    map[1][0] = number(0, -1, 5);
    map[1][1] = number(0, 1, 5);
    let mut centre = zero_centre();
    centre[0] = number(-5, 3, 10);
    centre[1] = number(5, -1, 10);
    Parameters {
        coordinates: Coordinates::Configuration,
        centre,
        map,
        shape: Shape::Point {
            base: vec![[-1, 1], [-1, 0]],
            offset: whole(),
            radii: vec![rational(1, 32); 3],
            radius: rational(1, 16),
            ratio: rational(2, 1),
            window: None,
        },
        beyond: Some(Beyond { axis: 1, inequality: 0 }),
    }
}

/// The tube along the arc `θ = 0`, `e ∈ [3/64, 3/16]`, offset
/// `(v, θ, ζ₁, ζ₃) ∈ (1/32)·(2, 1, 1, 1)·([0, 1] × [−1, 1]³)`.
fn arc(sign: Sign) -> Parameters {
    Parameters {
        coordinates: Coordinates::ArcPlane(sign),
        centre: zero_centre(),
        map: diagonal([q(2), q(1), q(1), q(1)]),
        shape: Shape::Tube {
            base: vec![interval(frac(3, 64), frac(3, 16))],
            offset: arc_offset(),
            radius: rational(1, 32),
        },
        beyond: None,
    }
}

/// The point blow-up at the arc's endpoint `e* = √5 − 2`, base radii `1/16`
/// below and `1/32` above, offset as the arc's tube, ratio 1.
fn endpoint(sign: Sign) -> Parameters {
    let mut centre = zero_centre();
    centre[0] = number(-2, 1, 1);
    Parameters {
        coordinates: Coordinates::ArcPlane(sign),
        centre,
        map: diagonal([q(2), q(1), q(1), q(1)]),
        shape: Shape::Point {
            base: vec![[-1, 1]],
            offset: arc_offset(),
            radii: vec![rational(1, 16), rational(1, 32)],
            radius: rational(1, 32),
            ratio: rational(1, 1),
            window: None,
        },
        beyond: None,
    }
}

/// The point blow-up at the crossing `e = 0`, radii `1/16`, offset
/// `(1/16)·(5/4, 1, 1, 1)·([0, 1] × [−1, 1]³)`, ratio `5/4`, with the window
/// hand-over on the base face `e > 0` to the family sheared by
/// `θ = θ' − e` with ratio at most `1/4`.
fn crossing(sign: Sign) -> Parameters {
    Parameters {
        coordinates: Coordinates::ArcPlane(sign),
        centre: zero_centre(),
        map: diagonal([frac(5, 4), q(1), q(1), q(1)]),
        shape: Shape::Point {
            base: vec![[-1, 1]],
            offset: arc_offset(),
            radii: vec![rational(1, 16), rational(1, 16)],
            radius: rational(1, 16),
            ratio: rational(5, 4),
            window: Some(Window {
                faces: vec![Face { axis: 0, side: 1 }],
                shear: vec![vec![rational(0, 1)], vec![rational(1, 1)], vec![rational(0, 1)], vec![rational(0, 1)]],
                radius: rational(1, 4),
            }),
        },
        beyond: None,
    }
}

/// The pin of an Exotic cover.
pub fn exotic_pin(name: &str) -> Option<Pin> {
    let parameters = match name {
        "square" => square(),
        "pentagon" => pentagon(),
        "arc+" => arc(Sign::Plus),
        "arc-" => arc(Sign::Minus),
        "endpoint+" => endpoint(Sign::Plus),
        "endpoint-" => endpoint(Sign::Minus),
        "crossing+" => crossing(Sign::Plus),
        "crossing-" => crossing(Sign::Minus),
        _ => return None,
    };
    Some(Pin {
        name: name.to_string(),
        scope: Scope::Domain,
        parameters,
    })
}

/// The view rectangle `[s_lo, s_hi] × [t_lo, t_hi]` of Local cover `number`:
/// its row enlarged by [`MARGIN`] and cut to `[0, 2/3] × [0, 2/5]`.
pub fn local_rectangle(number: usize) -> Option<[[Q; 2]; 2]> {
    let row = LOCAL.get(number)?;
    let value = |i: usize| parse_rational(row[i]).expect("canonical pinned rational");
    let m = frac(MARGIN.0, MARGIN.1);
    let (zero, s_max, t_max) = (Q::zero(), frac(2, 3), frac(2, 5));
    Some([
        [(value(0) - &m).max(zero.clone()), (value(1) + &m).min(s_max)],
        [(value(2) - &m).max(zero), (value(3) + &m).min(t_max)],
    ])
}

/// The pin of Local cover `number`: the tube over the aligned set above its
/// view rectangle ([`local_rectangle`]).
pub fn local_pin(number: usize) -> Option<Pin> {
    let [s, t] = local_rectangle(number)?;
    let radius = parse_rational(LOCAL[number][4]).expect("canonical pinned rational");
    let [s_lo, s_hi] = s;
    let [t_lo, t_hi] = t;
    Some(Pin {
        name: number.to_string(),
        scope: Scope::All,
        parameters: Parameters {
            coordinates: Coordinates::Configuration,
            centre: zero_centre(),
            map: diagonal([q(1), q(1), q(1), q(1)]),
            shape: Shape::Tube {
                base: vec![interval(s_lo, s_hi), interval(t_lo, t_hi)],
                offset: whole(),
                radius: Rational(radius),
            },
            beyond: None,
        },
    })
}

/// The data file of an Exotic cover.
pub fn exotic_data(name: &str) -> Option<&'static [u8]> {
    EXOTIC.iter().position(|n| *n == name).map(|i| EXOTIC_DATA[i])
}

/// The data file of Local cover `number`.
pub fn local_data(number: usize) -> Option<&'static [u8]> {
    LOCAL_DATA.get(number).copied()
}

/// A cover of the catalogue that failed to load.
#[derive(Debug)]
pub struct CatalogueError {
    pub cover: String,
    pub error: ProofError,
}

impl fmt::Display for CatalogueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "cover {}: {}", self.cover, self.error)
    }
}

impl std::error::Error for CatalogueError {}

fn load(pin: Pin, data: &[u8], threads: NonZeroUsize, stop: &AtomicBool) -> Result<ProvedCover, CatalogueError> {
    ProvedCover::load(&pin, data, threads, stop).map_err(|error| CatalogueError { cover: pin.name, error })
}

/// Loads and verifies the Exotic covers, in the order of [`EXOTIC`].
pub fn exotic(threads: NonZeroUsize, stop: &AtomicBool) -> Result<Vec<ProvedCover>, CatalogueError> {
    EXOTIC
        .iter()
        .map(|name| {
            let pin = exotic_pin(name).expect("pinned");
            load(pin, exotic_data(name).expect("data"), threads, stop)
        })
        .collect()
}

/// Loads and verifies the Local covers; cover `i` is at index `i`.
pub fn local(threads: NonZeroUsize, stop: &AtomicBool) -> Result<Vec<ProvedCover>, CatalogueError> {
    (0..LOCAL_COUNT)
        .map(|number| load(local_pin(number).expect("pinned"), LOCAL_DATA[number], threads, stop))
        .collect()
}

#[cfg(test)]
mod tests;
