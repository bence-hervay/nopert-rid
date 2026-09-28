//! Cell trees: complete midpoint-bisection trees of a zoom's root box,
//! written in preorder, their leaf cells, and an independent exact check
//! that given cells tile a box.
use super::ZoomCell;
use crate::arithmetic::exact::Q;
use crate::arithmetic::polynomial::VARIABLES;
use std::fmt;

/// The token of a leaf; a split along variable `j` is the digit `j`.
pub const LEAF: u8 = b'.';

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TreeError {
    /// A byte that is neither a variable digit below 5 nor the leaf token.
    Token { position: usize, byte: u8 },
    /// The tokens end before every node has been given.
    Incomplete,
    /// Tokens remain after the tree is complete.
    Trailing { position: usize },
}

impl fmt::Display for TreeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TreeError::Token { position, byte } => write!(f, "tree byte {position} is {byte:#04x}"),
            TreeError::Incomplete => write!(f, "the tree ends early"),
            TreeError::Trailing { position } => write!(f, "tree bytes remain from {position}"),
        }
    }
}

impl std::error::Error for TreeError {}

/// The leaf cells, in preorder, of the tree `tokens` on the box `root`.
/// Each node is one byte: the digit `j` splits the node's box along
/// variable `j` at the exact midpoint into two closed halves, whose subtrees
/// follow (lower half first); [`LEAF`] ends a branch.
pub fn leaves(root: &ZoomCell, tokens: &str) -> Result<Vec<ZoomCell>, TreeError> {
    let mut pending = vec![root.clone()];
    let mut out = Vec::new();
    for (position, byte) in tokens.bytes().enumerate() {
        let cell = pending.pop().ok_or(TreeError::Trailing { position })?;
        match byte {
            LEAF => out.push(cell),
            b'0'..=b'4' => {
                let [lower, upper] = cell.bisect(usize::from(byte - b'0'));
                pending.push(upper);
                pending.push(lower);
            }
            _ => return Err(TreeError::Token { position, byte }),
        }
    }
    if pending.is_empty() {
        Ok(out)
    } else {
        Err(TreeError::Incomplete)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TilingError {
    /// The root box has zero width in some variable.
    FlatRoot,
    /// This cell is not inside the root box.
    Outside { cell: usize },
    /// These two cells share interior points.
    Overlap { first: usize, second: usize },
    /// The cells' volumes do not add up to the root's.
    Volume,
}

impl fmt::Display for TilingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TilingError::FlatRoot => write!(f, "the root box is flat"),
            TilingError::Outside { cell } => write!(f, "cell {cell} leaves the root box"),
            TilingError::Overlap { first, second } => write!(f, "cells {first} and {second} overlap"),
            TilingError::Volume => write!(f, "the cells do not fill the root box"),
        }
    }
}

impl std::error::Error for TilingError {}

fn volume(cell: &ZoomCell) -> Q {
    cell.iter().fold(Q::one(), |v, axis| v * axis.width())
}

/// Exactly whether the closed `cells` tile the closed `root` (positive
/// widths): every cell lies in the root, no two cells share an interior
/// point, and their volumes add up to the root's. Then their union is the
/// root, independently of how the cells were produced.
pub fn check_tiling(root: &ZoomCell, cells: &[ZoomCell]) -> Result<(), TilingError> {
    if root.iter().any(|axis| axis.width() == Q::zero()) {
        return Err(TilingError::FlatRoot);
    }
    for (index, cell) in cells.iter().enumerate() {
        if (0..VARIABLES).any(|j| cell[j].lo() < root[j].lo() || cell[j].hi() > root[j].hi()) {
            return Err(TilingError::Outside { cell: index });
        }
    }
    // Sweep along the variable with the most distinct lower ends (so that
    // few cells share a slab): cells whose intervals in that variable overlap
    // in more than a point are compared in all variables.
    let distinct = |j: usize| {
        let mut ends: Vec<&Q> = cells.iter().map(|cell| cell[j].lo()).collect();
        ends.sort();
        ends.dedup();
        ends.len()
    };
    let sweep = (0..VARIABLES).max_by_key(|&j| (distinct(j), std::cmp::Reverse(j))).expect("five variables");
    let mut order: Vec<usize> = (0..cells.len()).collect();
    order.sort_by(|&a, &b| cells[a][sweep].lo().cmp(cells[b][sweep].lo()));
    for (rank, &a) in order.iter().enumerate() {
        for &b in &order[rank + 1..] {
            if cells[b][sweep].lo() >= cells[a][sweep].hi() {
                break;
            }
            let interiors_meet = (0..VARIABLES)
                .all(|j| cells[a][j].lo() < cells[b][j].hi() && cells[b][j].lo() < cells[a][j].hi());
            if interiors_meet {
                return Err(TilingError::Overlap {
                    first: a.min(b),
                    second: a.max(b),
                });
            }
        }
    }
    let total = cells.iter().fold(Q::zero(), |sum, cell| sum + volume(cell));
    if total != volume(root) {
        return Err(TilingError::Volume);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
