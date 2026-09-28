//! Necessary inequalities (witnesses): support gaps and domain inequalities,
//! as data, in their canonical JSON form and as polynomials in any view and
//! rotation. A witness `g` satisfies `g ≥ 0` at every weakly contained
//! configuration of `D` whose view lies in a view set on which it is valid.
//!
use crate::arithmetic::exact::{parse_rational, ExactError, QSqrt5};
use crate::arithmetic::polynomial::{Polynomial, PolynomialError};
use crate::problem::domain::{self, DomainError, CONSTRAINT_COUNT};
use crate::problem::geometry::{self, Point, View, VERTEX_COUNT};
use serde::{Deserialize, Deserializer, Serialize};
use std::cmp::Ordering;
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WitnessError {
    VertexOutOfRange { index: usize },
    /// The ordered pair is not an edge of the RID.
    NotAnEdge { from: usize, to: usize },
    /// A planar direction `(0, 0)`.
    ZeroDirection,
    InequalityOutOfRange { index: usize },
    /// Validity was asked on an empty view set.
    NoViews,
    /// A view vector whose third coordinate is not positive.
    ViewNotPositive { view: usize },
    /// The support normal vanishes at a view.
    ZeroNormal { view: usize },
    /// `n(u)·(h - w) < 0` at view `view` for the vertex `vertex`.
    Invalid { view: usize, vertex: usize },
    /// Building the polynomial overflowed an exponent.
    Polynomial(PolynomialError),
    /// A number of the JSON form is not a canonical rational.
    Number(ExactError),
    /// A JSON object whose fields (listed) are not those of one form of the
    /// expected kind (`expected`).
    Form { fields: Vec<&'static str>, expected: &'static str },
}

impl fmt::Display for WitnessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WitnessError::VertexOutOfRange { index } => {
                write!(f, "vertex {index} is not below {VERTEX_COUNT}")
            }
            WitnessError::NotAnEdge { from, to } => write!(f, "[{from}, {to}] is not an edge"),
            WitnessError::ZeroDirection => write!(f, "the planar direction is zero"),
            WitnessError::InequalityOutOfRange { index } => {
                write!(f, "domain inequality {index} is not below {CONSTRAINT_COUNT}")
            }
            WitnessError::NoViews => write!(f, "validity needs at least one view"),
            WitnessError::ViewNotPositive { view } => {
                write!(f, "view {view} does not have a positive third coordinate")
            }
            WitnessError::ZeroNormal { view } => {
                write!(f, "the support normal vanishes at view {view}")
            }
            WitnessError::Invalid { view, vertex } => {
                write!(f, "vertex {vertex} lies beyond the support at view {view}")
            }
            WitnessError::Polynomial(error) => write!(f, "{error}"),
            WitnessError::Number(error) => write!(f, "{error}"),
            WitnessError::Form { fields, expected } => {
                write!(f, "an object with the fields {fields:?} is not {expected}")
            }
        }
    }
}

impl std::error::Error for WitnessError {}

impl From<PolynomialError> for WitnessError {
    fn from(error: PolynomialError) -> Self {
        WitnessError::Polynomial(error)
    }
}

fn check_vertex(index: usize) -> Result<(), WitnessError> {
    if index < VERTEX_COUNT {
        Ok(())
    } else {
        Err(WitnessError::VertexOutOfRange { index })
    }
}

/// An oriented edge of the RID, `from → to`. Its screen normal is
/// `n(u) = u × (v_to - v_from)` and its hole contact vertex is `from`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edge {
    from: usize,
    to: usize,
}

impl Edge {
    pub fn new(from: usize, to: usize) -> Result<Self, WitnessError> {
        check_vertex(from)?;
        check_vertex(to)?;
        if !geometry::is_edge(from, to) {
            return Err(WitnessError::NotAnEdge { from, to });
        }
        Ok(Self { from, to })
    }
    pub fn from(&self) -> usize {
        self.from
    }
    pub fn to(&self) -> usize {
        self.to
    }
    /// The edge vector `v_to - v_from`.
    pub fn direction(&self) -> Point {
        let v = geometry::vertices();
        geometry::difference(&v[self.to], &v[self.from])
    }
}

/// A fixed nonzero planar direction `(n₁, n₂)` in screen coordinates with a
/// hole contact vertex. Its support normal in space is
/// `n(u) = (n₁u₃, n₂u₃, -n₁u₁ - n₂u₂)` ([`Support::normal`]), so that
/// `n(u)·X = (n₁, n₂)·T_u(X)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Direction {
    screen_normal: [QSqrt5; 2],
    contact: usize,
}

impl Direction {
    pub fn new(screen_normal: [QSqrt5; 2], contact: usize) -> Result<Self, WitnessError> {
        check_vertex(contact)?;
        if screen_normal.iter().all(QSqrt5::is_zero) {
            return Err(WitnessError::ZeroDirection);
        }
        Ok(Self {
            screen_normal,
            contact,
        })
    }
    /// The planar direction `(n₁, n₂)` in screen coordinates (not the
    /// support normal `n(u)` in space).
    pub fn screen_normal(&self) -> &[QSqrt5; 2] {
        &self.screen_normal
    }
    pub fn contact(&self) -> usize {
        self.contact
    }
}

/// A support direction `n(u) ⟂ u`, linear in the view vector `u`, with its
/// hole contact vertex `h`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Support {
    Edge(Edge),
    Direction(Direction),
}

impl Support {
    /// The hole contact vertex `h`.
    pub fn contact(&self) -> usize {
        match self {
            Support::Edge(edge) => edge.from,
            Support::Direction(direction) => direction.contact,
        }
    }

    /// The support normal `n(u)` in space, for a polynomial view vector.
    pub fn normal(&self, u: &[Polynomial; 3]) -> [Polynomial; 3] {
        match self {
            Support::Edge(edge) => geometry::polynomial_cross(u, &edge.direction()),
            Support::Direction(direction) => {
                let [n1, n2] = &direction.screen_normal;
                [
                    u[2].scale(n1),
                    u[2].scale(n2),
                    -&(&u[0].scale(n1) + &u[1].scale(n2)),
                ]
            }
        }
    }

    /// The support normal `n(u)` in space, for an exact view vector.
    pub fn normal_at(&self, u: &Point) -> Point {
        match self {
            Support::Edge(edge) => geometry::cross(u, &edge.direction()),
            Support::Direction(direction) => {
                let [n1, n2] = &direction.screen_normal;
                [&u[2] * n1, &u[2] * n2, -&(&(&u[0] * n1) + &(&u[1] * n2))]
            }
        }
    }

    /// Validity on the view set generated by `views`: at every given view,
    /// `n(u) ≠ 0` and `n(u)·(h - w) ≥ 0` for all 60 vertices `w`. Each view
    /// must have `u₃ > 0`. Because `n(u)` is linear in `u`, validity then holds
    /// at every nonzero nonnegative combination of the views.
    pub fn check_valid(&self, views: &[Point]) -> Result<(), WitnessError> {
        check_views(views)?;
        let vertices = geometry::vertices();
        let h = &vertices[self.contact()];
        for (view, u) in views.iter().enumerate() {
            let n = self.normal_at(u);
            if n.iter().all(QSqrt5::is_zero) {
                return Err(WitnessError::ZeroNormal { view });
            }
            let support = geometry::dot(&n, h);
            for (vertex, w) in vertices.iter().enumerate() {
                if support < geometry::dot(&n, w) {
                    return Err(WitnessError::Invalid { view, vertex });
                }
            }
        }
        Ok(())
    }
}

fn check_views(views: &[Point]) -> Result<(), WitnessError> {
    if views.is_empty() {
        return Err(WitnessError::NoViews);
    }
    if let Some(view) = views.iter().position(|u| u[2].sign() != Ordering::Greater) {
        return Err(WitnessError::ViewNotPositive { view });
    }
    Ok(())
}

/// A support gap: a support with hole contact vertex `h` and a plug vertex `p`,
/// `g = (1 + |r|²) n(u)·h - n(u)·R̂(r) p`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "Fields", into = "GapForm")]
pub struct Gap {
    support: Support,
    plug: usize,
}

impl Gap {
    pub fn new(support: Support, plug: usize) -> Result<Self, WitnessError> {
        check_vertex(plug)?;
        Ok(Self { support, plug })
    }
    pub fn support(&self) -> &Support {
        &self.support
    }
    pub fn plug(&self) -> usize {
        self.plug
    }

    /// The gap polynomial at the view `view` and rotation vector `rotation`.
    /// It is linear in `u`: at `u = λ(s, t, 1)` it is `λ` times the affine gap.
    pub fn polynomial(
        &self,
        view: &View,
        rotation: &[Polynomial; 3],
    ) -> Result<Polynomial, WitnessError> {
        let vertices = geometry::vertices();
        let n = self.support.normal(&view.vector());
        let cayley = geometry::cayley_matrix(rotation)?;
        let turned = geometry::polynomial_matrix_point(&cayley, &vertices[self.plug]);
        let scale = &Polynomial::constant(QSqrt5::one()) + &geometry::squared_norm(rotation)?;
        let hole = geometry::polynomial_dot(&n, &vertices[self.support.contact()]);
        let mut plug = Polynomial::zero();
        for j in 0..3 {
            plug = &plug + &n[j].mul(&turned[j])?;
        }
        Ok(&scale.mul(&hole)? - &plug)
    }
}

/// The maximum form of a support gap: an ordered edge, whose normal
/// `n(u) = u × (v_to - v_from)` fixes a direction perpendicular to the view,
/// and a plug vertex `p`:
/// `G = max over the 60 hole vertices w of (1 + |r|²) n(u)·w - n(u)·R̂(r) p`.
/// `max_w n(u)·w` is the hole's own support value in the direction `n(u)`,
/// so `G ≥ 0` at every weakly contained configuration with no validity
/// hypothesis. It is not a polynomial; its members, one per hole
/// vertex `w`, are. Its JSON form is `{"edge":[from,to],"vertex":p}`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "MaximumGapForm", into = "MaximumGapForm")]
pub struct MaximumGap {
    edge: Edge,
    plug: usize,
}

impl MaximumGap {
    pub fn new(edge: Edge, plug: usize) -> Result<Self, WitnessError> {
        check_vertex(plug)?;
        Ok(Self { edge, plug })
    }
    pub fn edge(&self) -> &Edge {
        &self.edge
    }
    pub fn plug(&self) -> usize {
        self.plug
    }

    /// The 60 members `(1 + |r|²) n(u)·w - n(u)·R̂(r) p`, in vertex order, at
    /// the view `view` and rotation vector `rotation`. Member `w` is the gap
    /// polynomial of the support `n(u)` with contact vertex `w`; each is
    /// linear in `u`, so `G(λu) = λ G(u)` for `λ > 0`.
    pub fn members(
        &self,
        view: &View,
        rotation: &[Polynomial; 3],
    ) -> Result<Vec<Polynomial>, WitnessError> {
        let vertices = geometry::vertices();
        let n = Support::Edge(self.edge.clone()).normal(&view.vector());
        let cayley = geometry::cayley_matrix(rotation)?;
        let turned = geometry::polynomial_matrix_point(&cayley, &vertices[self.plug]);
        let scale = &Polynomial::constant(QSqrt5::one()) + &geometry::squared_norm(rotation)?;
        let mut plug = Polynomial::zero();
        for j in 0..3 {
            plug = &plug + &n[j].mul(&turned[j])?;
        }
        vertices
            .iter()
            .map(|w| Ok(&scale.mul(&geometry::polynomial_dot(&n, w))? - &plug))
            .collect()
    }
}

/// The JSON form of a maximum form, `{"edge":[from,to],"vertex":p}`.
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MaximumGapForm {
    edge: [usize; 2],
    vertex: usize,
}

impl From<MaximumGap> for MaximumGapForm {
    fn from(maximum: MaximumGap) -> Self {
        Self {
            edge: [maximum.edge.from, maximum.edge.to],
            vertex: maximum.plug,
        }
    }
}

impl TryFrom<MaximumGapForm> for MaximumGap {
    type Error = WitnessError;
    fn try_from(form: MaximumGapForm) -> Result<Self, WitnessError> {
        let [from, to] = form.edge;
        MaximumGap::new(Edge::new(from, to)?, form.vertex)
    }
}

/// Minus one of the inequalities of `D`: its sign conclusions hold only on `D`.
/// Its JSON form is the bare inequality index.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "usize", into = "usize")]
pub struct DomainInequality {
    index: usize,
}

impl DomainInequality {
    pub fn new(index: usize) -> Result<Self, WitnessError> {
        if index >= CONSTRAINT_COUNT {
            return Err(WitnessError::InequalityOutOfRange { index });
        }
        Ok(Self { index })
    }
    pub fn index(&self) -> usize {
        self.index
    }
}

/// A necessary inequality `g ≥ 0`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "Fields", into = "WitnessForm")]
pub enum Witness {
    Gap(Gap),
    Domain(DomainInequality),
}

impl Witness {
    /// The witness polynomial `g`. For a homogeneous view `u = λ(s, t, 1)` it is
    /// `λ^k` times the affine value, with `k` = 1 for a gap and
    /// `domain::view_degree` for a domain inequality.
    pub fn polynomial(
        &self,
        view: &View,
        rotation: &[Polynomial; 3],
    ) -> Result<Polynomial, WitnessError> {
        match self {
            Witness::Gap(gap) => gap.polynomial(view, rotation),
            Witness::Domain(inequality) => {
                match domain::polynomial(inequality.index, view, rotation) {
                    Ok(c) => Ok(-&c),
                    Err(DomainError::Polynomial(error)) => Err(WitnessError::Polynomial(error)),
                    Err(error) => unreachable!("index checked at construction: {error}"),
                }
            }
        }
    }

    /// Validity on the view set generated by `views` (see
    /// [`Support::check_valid`]). A domain inequality needs no support, so
    /// only the views themselves are checked.
    pub fn check_valid(&self, views: &[Point]) -> Result<(), WitnessError> {
        match self {
            Witness::Gap(gap) => gap.support.check_valid(views),
            Witness::Domain(_) => check_views(views),
        }
    }
}

// ---- The canonical JSON form ---------------------------------------------
//
// The only serialised form of witnesses, used by certificate records and
// cover files alike. Every field is required, unknown and repeated fields
// are refused, numbers of Q(√5) are pairs of canonical rational strings, and a
// parsed value is validated exactly like one built by the constructors.

/// `{"edge":[from,to],"vertex":plug}` or
/// `{"direction":[[a1,b1],[a2,b2]],"contact":h,"vertex":plug}` for the
/// planar direction `(a1 + b1√5, a2 + b2√5)`.
#[derive(Serialize)]
#[serde(untagged)]
enum GapForm {
    Edge {
        edge: [usize; 2],
        vertex: usize,
    },
    Direction {
        direction: [[String; 2]; 2],
        contact: usize,
        vertex: usize,
    },
}

/// A gap's own form, or `{"inequality":i}` for a domain inequality.
#[derive(Serialize)]
#[serde(untagged)]
enum WitnessForm {
    Gap(Gap),
    Domain { inequality: DomainInequality },
}

/// Every field of every form, each optional, so that parsing names the
/// reason of a refusal: the constructor's error, or the fields that form no
/// witness. (An untagged enum would only say that no form matched.) A field
/// that is present must hold a value of its type; `null` is refused, not
/// read as absent. Serde would also read a JSON array positionally, but a
/// present element is never absent, so no array has the fields of one form.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fields {
    #[serde(default, deserialize_with = "present")]
    edge: Option<[usize; 2]>,
    #[serde(default, deserialize_with = "present")]
    direction: Option<[[String; 2]; 2]>,
    #[serde(default, deserialize_with = "present")]
    contact: Option<usize>,
    #[serde(default, deserialize_with = "present")]
    vertex: Option<usize>,
    #[serde(default, deserialize_with = "present")]
    inequality: Option<usize>,
}

/// A present field: its value, never `null`.
fn present<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(deserializer).map(Some)
}

const GAP_FORMS: &str = "a support gap ({edge, vertex} or {direction, contact, vertex})";
const WITNESS_FORMS: &str =
    "a witness ({edge, vertex}, {direction, contact, vertex} or {inequality})";

impl Fields {
    fn names(&self) -> Vec<&'static str> {
        [
            ("edge", self.edge.is_some()),
            ("direction", self.direction.is_some()),
            ("contact", self.contact.is_some()),
            ("vertex", self.vertex.is_some()),
            ("inequality", self.inequality.is_some()),
        ]
        .into_iter()
        .filter_map(|(name, present)| present.then_some(name))
        .collect()
    }

    /// The witness these fields spell, through the constructors; `expected`
    /// names the accepted forms in the refusal of any other set of fields.
    fn parse(self, expected: &'static str) -> Result<Witness, WitnessError> {
        let fields = self.names();
        let (support, plug) = match self {
            Fields {
                edge: Some([from, to]),
                direction: None,
                contact: None,
                vertex: Some(plug),
                inequality: None,
            } => (Support::Edge(Edge::new(from, to)?), plug),
            Fields {
                edge: None,
                direction: Some([a, b]),
                contact: Some(contact),
                vertex: Some(plug),
                inequality: None,
            } => {
                let screen_normal = [parse_number(&a)?, parse_number(&b)?];
                (Support::Direction(Direction::new(screen_normal, contact)?), plug)
            }
            Fields {
                edge: None,
                direction: None,
                contact: None,
                vertex: None,
                inequality: Some(index),
            } => return Ok(Witness::Domain(DomainInequality::new(index)?)),
            _ => return Err(WitnessError::Form { fields, expected }),
        };
        Ok(Witness::Gap(Gap::new(support, plug)?))
    }
}

fn number_form(x: &QSqrt5) -> [String; 2] {
    [x.rational_part().to_string(), x.sqrt5_part().to_string()]
}

fn parse_number([a, b]: &[String; 2]) -> Result<QSqrt5, WitnessError> {
    let part = |text: &str| parse_rational(text).map_err(WitnessError::Number);
    Ok(QSqrt5::new(part(a)?, part(b)?))
}

impl From<Gap> for GapForm {
    fn from(gap: Gap) -> Self {
        match gap.support {
            Support::Edge(edge) => GapForm::Edge {
                edge: [edge.from, edge.to],
                vertex: gap.plug,
            },
            Support::Direction(direction) => GapForm::Direction {
                direction: [
                    number_form(&direction.screen_normal[0]),
                    number_form(&direction.screen_normal[1]),
                ],
                contact: direction.contact,
                vertex: gap.plug,
            },
        }
    }
}

impl TryFrom<Fields> for Gap {
    type Error = WitnessError;
    fn try_from(fields: Fields) -> Result<Self, WitnessError> {
        let names = fields.names();
        match fields.parse(GAP_FORMS)? {
            Witness::Gap(gap) => Ok(gap),
            Witness::Domain(_) => Err(WitnessError::Form {
                fields: names,
                expected: GAP_FORMS,
            }),
        }
    }
}

impl From<DomainInequality> for usize {
    fn from(inequality: DomainInequality) -> usize {
        inequality.index
    }
}

impl TryFrom<usize> for DomainInequality {
    type Error = WitnessError;
    fn try_from(index: usize) -> Result<Self, WitnessError> {
        DomainInequality::new(index)
    }
}

impl From<Witness> for WitnessForm {
    fn from(witness: Witness) -> Self {
        match witness {
            Witness::Gap(gap) => WitnessForm::Gap(gap),
            Witness::Domain(inequality) => WitnessForm::Domain { inequality },
        }
    }
}

impl TryFrom<Fields> for Witness {
    type Error = WitnessError;
    fn try_from(fields: Fields) -> Result<Self, WitnessError> {
        fields.parse(WITNESS_FORMS)
    }
}

#[cfg(test)]
mod tests;
