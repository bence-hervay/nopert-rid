//! Adversarial tests of the layout: every point that `place`
//! creates (as opposed to maps) lies on the padded boundary, at least the
//! padding away from the viewport; a polygon around the viewport becomes
//! the padded rectangle.
use super::*;
use crate::random::Random;
use rid::arithmetic::exact::frac;

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn on_boundary(r: &Rectangle, p: &Point) -> bool {
    let inside = r[0].contains(&p[0]) && r[1].contains(&p[1]);
    inside && (&p[0] == r[0].lo() || &p[0] == r[0].hi() || &p[1] == r[1].lo() || &p[1] == r[1].hi())
}

#[test]
fn created_points_lie_on_the_padded_boundary() {
    let canvas = Canvas::new(frac(200, 1), frac(200, 1)).unwrap();
    let panel = Panel::new(
        &canvas,
        [interval(frac(20, 1), frac(120, 1)), interval(frac(30, 1), frac(80, 1))],
        [interval(frac(-2, 1), frac(2, 1)), interval(frac(-1, 1), frac(1, 1))],
    )
    .unwrap();
    let padding = frac(3, 10);
    let padded = panel.viewport().clone().map(|i| interval(i.lo() - &padding, i.hi() + &padding));
    let mut random = Random::new(904);
    let coordinate = |random: &mut Random| frac(random.integer(-600, 600), random.integer(1, 97));
    for case in 0..300 {
        let n = random.integer(1, 7) as usize;
        let points: Vec<Point> = (0..n).map(|_| [coordinate(&mut random), coordinate(&mut random)]).collect();
        let shape = if case % 2 == 0 && n >= 2 {
            Shape::Segment([points[0].clone(), points[1].clone()])
        } else {
            Shape::Polygon(points.clone())
        };
        let mapped: Vec<Point> = points.iter().map(|p| panel.to_page(p)).collect();
        for placed in panel.place(&[shape], &padding) {
            let out = match placed {
                Shape::Segment(p) => p.to_vec(),
                Shape::Polygon(p) => p,
                Shape::Circle { .. } => unreachable!("no circles were placed"),
            };
            for p in &out {
                assert!(mapped.contains(p) || on_boundary(&padded, p), "case {case}: {p:?} created off the padded boundary");
                assert!(padded[0].contains(&p[0]) && padded[1].contains(&p[1]));
            }
        }
    }
}

#[test]
fn a_polygon_around_the_viewport_becomes_the_padded_rectangle() {
    let canvas = Canvas::new(frac(200, 1), frac(200, 1)).unwrap();
    let panel = Panel::new(
        &canvas,
        [interval(frac(20, 1), frac(120, 1)), interval(frac(30, 1), frac(80, 1))],
        [interval(frac(-2, 1), frac(2, 1)), interval(frac(-1, 1), frac(1, 1))],
    )
    .unwrap();
    let big = frac(1000, 1);
    let around = Shape::Polygon(vec![[-&big, -&big], [big.clone(), -&big], [big.clone(), big.clone()], [-&big, big.clone()]]);
    let placed = panel.place(&[around], &frac(1, 1));
    let Shape::Polygon(points) = &placed[0] else { panic!() };
    let mut sorted = points.clone();
    sorted.sort();
    assert_eq!(sorted, vec![[frac(19, 1), frac(29, 1)], [frac(19, 1), frac(81, 1)], [frac(121, 1), frac(29, 1)], [frac(121, 1), frac(81, 1)]]);
}
