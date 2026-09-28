use super::*;
use rid::arithmetic::exact::frac;

fn worked() -> Vec<Q> {
    vec![q(-1), q(3), q(-3)]
}

#[test]
fn the_worked_example_has_the_documented_control_points() {
    let g = graph(&worked(), [q(0), q(1)], 8).unwrap();
    assert_eq!(g.points, vec![[q(0), q(-1)], [frac(1, 2), frac(1, 2)], [q(1), q(-1)]]);
    let left = graph(&worked(), [q(0), frac(1, 2)], 8).unwrap();
    assert_eq!(left.points, vec![[q(0), q(-1)], [frac(1, 4), frac(-1, 4)], [frac(1, 2), frac(-1, 4)]]);
    let right = graph(&worked(), [frac(1, 2), q(1)], 8).unwrap();
    assert_eq!(right.points, vec![[frac(1, 2), frac(-1, 4)], [frac(3, 4), frac(-1, 4)], [q(1), q(-1)]]);
}

#[test]
fn the_curve_passes_through_exact_values_and_the_polygon_joins_the_points() {
    let g = graph(&worked(), [q(0), q(1)], 4).unwrap();
    assert_eq!(g.curve.len(), 4);
    assert_eq!(g.curve[0], Shape::Segment([[q(0), q(-1)], [frac(1, 4), frac(-7, 16)]]));
    assert_eq!(g.curve[2], Shape::Segment([[frac(1, 2), frac(-1, 4)], [frac(3, 4), frac(-7, 16)]]));
    assert_eq!(g.polygon.len(), 2);
    assert_eq!(g.axis, vec![Shape::Segment([[q(0), q(0)], [q(1), q(0)]])]);
}

#[test]
fn every_value_on_the_interval_lies_between_the_extreme_control_points() {
    // The convex hull property, on random cubics and intervals.
    let mut seed = 7u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) % 21) as i64 - 10
    };
    for _ in 0..30 {
        let c: Vec<Q> = (0..4).map(|_| frac(next(), 3)).collect();
        let l = frac(next(), 4);
        let h = &l + &frac(next().abs() + 1, 5);
        let g = graph(&c, [l, h], 16).unwrap();
        let lo = g.points.iter().map(|p| p[1].clone()).min().unwrap();
        let hi = g.points.iter().map(|p| p[1].clone()).max().unwrap();
        for shape in &g.curve {
            let Shape::Segment([a, _]) = shape else { unreachable!() };
            assert!(lo <= a[1] && a[1] <= hi);
        }
    }
}

#[test]
fn degenerate_inputs_are_refused() {
    assert_eq!(graph(&[], [q(0), q(1)], 4), Err(GraphError::Degree));
    assert_eq!(graph(&worked(), [q(1), q(1)], 4), Err(GraphError::Interval));
    assert_eq!(graph(&worked(), [q(1), q(0)], 4), Err(GraphError::Interval));
}

#[test]
fn markers_are_diamonds_around_their_points() {
    let m = markers(&[[q(1), q(2)]], &frac(1, 10));
    assert_eq!(
        m,
        vec![Shape::Polygon(vec![
            [frac(11, 10), q(2)],
            [q(1), frac(21, 10)],
            [frac(9, 10), q(2)],
            [q(1), frac(19, 10)],
        ])]
    );
}
