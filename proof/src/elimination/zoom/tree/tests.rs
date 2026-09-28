//! Tests of cell trees against an independent recursive reading, and of the
//! tiling check against constructed failures.
use super::*;
use crate::arithmetic::exact::{frac, q, Interval};
use crate::elimination::zoom::tests::Random;

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn root() -> ZoomCell {
    ZoomCell([
        interval(q(0), frac(2, 3)),
        interval(frac(-1, 7), frac(3, 5)),
        interval(q(-1), q(1)),
        interval(q(0), frac(1, 1024)),
        interval(frac(5, 4), q(2)),
    ])
}

/// A random complete tree of at most `depth` levels, in preorder.
fn random_tree(random: &mut Random, depth: u32, out: &mut String) {
    if depth == 0 || random.below(3) == 0 {
        out.push('.');
        return;
    }
    out.push(char::from(b'0' + random.below(5) as u8));
    random_tree(random, depth - 1, out);
    random_tree(random, depth - 1, out);
}

/// Independent reading: recursion over the string, halves computed from
/// the parent's ends.
fn reference(cell: ZoomCell, tokens: &[u8], at: &mut usize, out: &mut Vec<ZoomCell>) {
    let token = tokens[*at];
    *at += 1;
    if token == b'.' {
        out.push(cell);
        return;
    }
    let axis = usize::from(token - b'0');
    let (lo, hi) = (cell[axis].lo().clone(), cell[axis].hi().clone());
    let middle = (&lo + &hi) * frac(1, 2);
    let mut lower = cell.clone();
    lower[axis] = interval(lo, middle.clone());
    let mut upper = cell;
    upper[axis] = interval(middle, hi);
    reference(lower, tokens, at, out);
    reference(upper, tokens, at, out);
}

#[test]
fn leaves_match_an_independent_reading_and_tile() {
    let mut random = Random::new(21);
    for _ in 0..300 {
        let mut tree = String::new();
        random_tree(&mut random, 9, &mut tree);
        let cells = leaves(&root(), &tree).unwrap();
        let mut expected = Vec::new();
        let mut at = 0;
        reference(root(), tree.as_bytes(), &mut at, &mut expected);
        assert_eq!(at, tree.len());
        assert_eq!(cells, expected);
        assert_eq!(cells.len(), tree.bytes().filter(|&b| b == LEAF).count());
        check_tiling(&root(), &cells).unwrap();
        // Every random point of the root lies in a leaf.
        for _ in 0..20 {
            let x = random.point(&root(), 12);
            assert!(cells.iter().any(|c| (0..5).all(|j| c[j].contains(&x[j]))));
        }
    }
}

#[test]
fn malformed_trees_are_refused() {
    let r = root();
    assert_eq!(leaves(&r, "."), Ok(vec![r.clone()]));
    assert_eq!(leaves(&r, ""), Err(TreeError::Incomplete));
    assert_eq!(leaves(&r, "0."), Err(TreeError::Incomplete));
    assert_eq!(leaves(&r, "0..."), Err(TreeError::Trailing { position: 3 }));
    assert_eq!(leaves(&r, "..").unwrap_err(), TreeError::Trailing { position: 1 });
    assert_eq!(leaves(&r, "5.."), Err(TreeError::Token { position: 0, byte: b'5' }));
    assert_eq!(leaves(&r, "0.x"), Err(TreeError::Token { position: 2, byte: b'x' }));
    assert_eq!(leaves(&r, "0. ."), Err(TreeError::Token { position: 2, byte: b' ' }));
    // Dropping or duplicating a subtree of a valid tree is always refused.
    let mut random = Random::new(22);
    for _ in 0..200 {
        let mut tree = String::new();
        random_tree(&mut random, 8, &mut tree);
        let cut = random.below(tree.len() as u64) as usize;
        let mut dropped = tree.clone();
        dropped.remove(cut);
        assert!(leaves(&r, &dropped).is_err(), "{tree} without byte {cut}");
        let mut doubled = tree.clone();
        doubled.push_str(&tree);
        assert!(leaves(&r, &doubled).is_err());
    }
}

#[test]
fn tiling_failures_are_found() {
    let r = root();
    let cells = leaves(&r, "01..2..").unwrap();
    assert_eq!(cells.len(), 4);
    check_tiling(&r, &cells).unwrap();
    // A dropped cell leaves a hole.
    assert_eq!(check_tiling(&r, &cells[1..]), Err(TilingError::Volume));
    // A repeated cell overlaps itself.
    let mut repeated = cells.clone();
    repeated.push(cells[2].clone());
    assert_eq!(check_tiling(&r, &repeated), Err(TilingError::Overlap { first: 2, second: 4 }));
    // An enlarged cell overlaps its neighbour.
    let mut enlarged = cells.clone();
    enlarged[0][1] = interval(r[1].lo().clone(), r[1].hi().clone());
    assert!(matches!(check_tiling(&r, &enlarged), Err(TilingError::Overlap { .. })));
    // A moved cell leaves the root.
    let mut moved = cells.clone();
    moved[3][4] = interval(frac(3, 2), frac(5, 2));
    assert_eq!(check_tiling(&r, &moved), Err(TilingError::Outside { cell: 3 }));
    // A flat root cannot be tiled by boxes with interiors.
    let mut flat = r.clone();
    flat[2] = Interval::point(q(0));
    assert_eq!(check_tiling(&flat, &[flat.clone()]), Err(TilingError::FlatRoot));
    // Cells meeting only in a face are fine; cells overlapping in a thin
    // slab are not.
    let mut slab = cells.clone();
    let h = slab[0][0].hi().clone();
    slab[0][0] = interval(slab[0][0].lo().clone(), &h + &frac(1, 1 << 30));
    assert!(check_tiling(&r, &slab).is_err());
}

#[test]
fn random_overlaps_and_gaps_are_found() {
    let mut random = Random::new(23);
    for _ in 0..200 {
        let mut tree = String::new();
        random_tree(&mut random, 7, &mut tree);
        let cells = leaves(&root(), &tree).unwrap();
        if cells.len() < 2 {
            continue;
        }
        // Replace one cell by a copy of another: an overlap and a gap.
        let (a, b) = (random.below(cells.len() as u64) as usize, random.below(cells.len() as u64) as usize);
        if a == b {
            continue;
        }
        let mut changed = cells.clone();
        changed[a] = cells[b].clone();
        assert!(check_tiling(&root(), &changed).is_err());
    }
}
