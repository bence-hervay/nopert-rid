//! The lemma, applied to one cell of one zoom: the only place in the crate
//! where it is applied.
use crate::arithmetic::exact::Q;
use crate::arithmetic::polynomial::{Polynomial, PolynomialError, Refusal, VARIABLES};
use crate::elimination::zoom::{Zoom, ZoomCell, ZoomFactor};
use crate::elimination::witness::{MaximumGap, Witness, WitnessError};
use crate::problem::geometry::Point;
use std::collections::HashSet;
use std::fmt;

/// Why a cell was not eliminated.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CellRefusal {
    /// The factor has a positive exponent on a variable that is not one of
    /// the zoom's scale variables.
    FactorOutsideScales { variable: usize },
    /// A divided variable takes negative values on the cell.
    NegativeScale { variable: usize },
    /// Building the pulled-back witness overflowed an exponent.
    Pullback(WitnessError),
    /// The zoom factor does not divide the pulled-back witness.
    NotDivisible(PolynomialError),
    /// `u₃ > 0` is not certified on the cell.
    ViewNotPositive(Refusal),
    /// The support is not valid at an exact corner view of the cell.
    InvalidSupport(WitnessError),
    /// The quotient is not certified strictly negative on the cell.
    NotNegative(Refusal),
    /// Maximum form: the zoom factor does not divide the member of hole
    /// vertex `vertex`.
    MemberNotDivisible { vertex: usize, error: PolynomialError },
    /// Maximum form: the quotient of the member of hole vertex `vertex` is
    /// not certified strictly negative on the cell.
    MemberNotNegative { vertex: usize, refusal: Refusal },
}

impl fmt::Display for CellRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CellRefusal::FactorOutsideScales { variable } => {
                write!(f, "variable {variable} is divided but is not a scale variable of the zoom")
            }
            CellRefusal::NegativeScale { variable } => {
                write!(f, "divided variable {variable} is negative on the cell")
            }
            CellRefusal::Pullback(error) => write!(f, "pull-back: {error}"),
            CellRefusal::NotDivisible(error) => write!(f, "{error}"),
            CellRefusal::ViewNotPositive(refusal) => write!(f, "u3 > 0 not certified: {refusal}"),
            CellRefusal::InvalidSupport(error) => write!(f, "support: {error}"),
            CellRefusal::NotNegative(refusal) => write!(f, "{refusal}"),
            CellRefusal::MemberNotDivisible { vertex, error } => write!(f, "member {vertex}: {error}"),
            CellRefusal::MemberNotNegative { vertex, refusal } => write!(f, "member {vertex}: {refusal}"),
        }
    }
}

impl std::error::Error for CellRefusal {}

/// A witness pulled back along a zoom and divided by the zoom factor:
/// `N = (g ∘ zoom) / y^factor`, exact. It can check any cell of the zoom.
#[derive(Clone, Debug)]
pub struct Quotient<'a> {
    zoom: &'a Zoom,
    witness: &'a Witness,
    factor: ZoomFactor,
    quotient: Polynomial,
}

impl<'a> Quotient<'a> {
    /// Refuses a factor outside the zoom's scale variables, and a
    /// monomial that does not divide the pulled-back witness exactly.
    pub fn new(zoom: &'a Zoom, witness: &'a Witness, factor: &ZoomFactor) -> Result<Self, CellRefusal> {
        check_factor(zoom, factor)?;
        let pulled_back = witness
            .polynomial(zoom.view(), zoom.rotation())
            .map_err(CellRefusal::Pullback)?;
        let quotient = pulled_back
            .divide_by_monomial(factor.exponents())
            .map_err(CellRefusal::NotDivisible)?;
        Ok(Self {
            zoom,
            witness,
            factor: *factor,
            quotient,
        })
    }

    /// The quotient `N` as a polynomial in the zoom variables.
    pub fn polynomial(&self) -> &Polynomial {
        &self.quotient
    }

    /// The lemma's remaining hypotheses on the closed `cell`: divided
    /// variables nonnegative, `u₃ > 0`, the witness valid at the exact corner
    /// views, and every tensor Bernstein coefficient of `N` strictly negative.
    pub fn check(&self, cell: &ZoomCell) -> Result<(), CellRefusal> {
        self.check_remembering(cell, &mut HashSet::new())
    }

    /// [`Quotient::check`] for several cells in order, stopping at the first
    /// refusal, which is returned with the cell's position. Neighbouring
    /// cells share corners, so the exact views at which the witness was found
    /// valid are remembered for the rest of this call.
    pub fn check_all<'c>(&self, cells: impl IntoIterator<Item = &'c ZoomCell>) -> Result<(), (usize, CellRefusal)> {
        let mut valid = HashSet::new();
        for (position, cell) in cells.into_iter().enumerate() {
            self.check_remembering(cell, &mut valid).map_err(|refusal| (position, refusal))?;
        }
        Ok(())
    }

    /// `valid` holds exact views at which the witness is known to be valid.
    fn check_remembering(&self, cell: &ZoomCell, valid: &mut HashSet<Point>) -> Result<(), CellRefusal> {
        check_signs(self.zoom, &self.factor, cell)?;
        for view in valid_views(self.witness, self.zoom.view_corners(cell)) {
            if !valid.contains(&view) {
                self.witness
                    .check_valid(std::slice::from_ref(&view))
                    .map_err(CellRefusal::InvalidSupport)?;
                valid.insert(view);
            }
        }
        self.quotient.certify_negative(cell).map_err(CellRefusal::NotNegative)
    }
}

/// Hypothesis 1, first half: only scale variables of the zoom are divided.
fn check_factor(zoom: &Zoom, factor: &ZoomFactor) -> Result<(), CellRefusal> {
    match (0..VARIABLES).find(|&j| factor.exponents()[j] > 0 && !zoom.is_scale(j)) {
        Some(variable) => Err(CellRefusal::FactorOutsideScales { variable }),
        None => Ok(()),
    }
}

/// Hypotheses 1, second half, and 2 on the cell: divided variables are
/// nonnegative, and `u₃ > 0`.
fn check_signs(zoom: &Zoom, factor: &ZoomFactor, cell: &ZoomCell) -> Result<(), CellRefusal> {
    if let Some(variable) = (0..VARIABLES).find(|&j| factor.exponents()[j] > 0 && cell[j].lo() < &Q::zero()) {
        return Err(CellRefusal::NegativeScale { variable });
    }
    let u3 = &zoom.view().vector()[2];
    (-u3).certify_negative(cell).map_err(CellRefusal::ViewNotPositive)
}

/// The corner views at which validity has to be checked. At a view where
/// the support normal vanishes, every validity condition `n(u)·(h − w) ≥ 0`
/// holds with equality, so such a view is dropped; the witness
/// module refuses zero normals because they carry no information. Domain
/// inequalities keep every view.
fn valid_views(witness: &Witness, views: Vec<Point>) -> Vec<Point> {
    match witness {
        Witness::Gap(gap) => views
            .into_iter()
            .filter(|u| !gap.support().normal_at(u).iter().all(|x| x.is_zero()))
            .collect(),
        Witness::Domain(_) => views,
    }
}

/// The lemma on one cell: if this returns `Ok`, no configuration of the
/// zoom's image of the closed `cell` is a fit lying in D.
pub fn check_cell(
    zoom: &Zoom,
    cell: &ZoomCell,
    witness: &Witness,
    factor: &ZoomFactor,
) -> Result<(), CellRefusal> {
    Quotient::new(zoom, witness, factor)?.check(cell)
}

/// The lemma for the maximum form on one cell:
/// the factor rule, divided variables nonnegative and `u₃ > 0` as for a
/// witness, then each of the 60 members, pulled back and divided by the same
/// zoom factor exactly, certified strictly negative. No validity is needed.
/// If this returns `Ok`, no configuration of the zoom's image of the closed
/// `cell` is a fit (in D or not); where the zoom factor is positive none is
/// even weakly contained (so with `ZoomFactor::ONE`, none at all), and where
/// it vanishes the configuration lies in a no-fit set.
pub fn check_maximum(
    zoom: &Zoom,
    cell: &ZoomCell,
    maximum: &MaximumGap,
    factor: &ZoomFactor,
) -> Result<(), CellRefusal> {
    let members = maximum.members(zoom.view(), zoom.rotation()).map_err(CellRefusal::Pullback)?;
    check_members(zoom, cell, &members, factor)
}

/// [`check_maximum`] with the members already built, for a caller that
/// reuses them across cells: `members` must be
/// `maximum.members(zoom.view(), zoom.rotation())` for the maximum form meant.
pub fn check_members(zoom: &Zoom, cell: &ZoomCell, members: &[Polynomial], factor: &ZoomFactor) -> Result<(), CellRefusal> {
    check_factor(zoom, factor)?;
    check_signs(zoom, factor, cell)?;
    for (vertex, member) in members.iter().enumerate() {
        let quotient = member
            .divide_by_monomial(factor.exponents())
            .map_err(|error| CellRefusal::MemberNotDivisible { vertex, error })?;
        quotient
            .certify_negative(cell)
            .map_err(|refusal| CellRefusal::MemberNotNegative { vertex, refusal })?;
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests;
