use super::*;
use crate::random::Random;
use rid::arithmetic::exact::frac;

fn interval(lo: Q, hi: Q) -> Interval {
    Interval::new(lo, hi).unwrap()
}

fn rect(x: (i64, i64), y: (i64, i64)) -> Rectangle {
    [interval(frac(x.0, 1), frac(x.1, 1)), interval(frac(y.0, 1), frac(y.1, 1))]
}

fn canvas() -> Canvas {
    Canvas::new(frac(180, 1), frac(100, 1)).unwrap()
}

fn point(x: i64, y: i64) -> Point {
    [frac(x, 1), frac(y, 1)]
}

#[test]
fn panels_have_one_scale_and_lie_on_the_canvas() {
    let c = canvas();
    assert!(Panel::new(&c, rect((10, 90), (10, 90)), rect((-5, 5), (-5, 5))).is_ok());
    assert_eq!(Panel::new(&c, rect((10, 90), (10, 80)), rect((-5, 5), (-5, 5))), Err(LayoutError::UnequalScale));
    assert_eq!(Panel::new(&c, rect((10, 190), (10, 90)), rect((-5, 5), (-5, 5))), Err(LayoutError::OutsideCanvas));
    assert_eq!(Panel::new(&c, rect((10, 10), (10, 90)), rect((-5, 5), (-5, 5))), Err(LayoutError::EmptyViewport));
    assert_eq!(Panel::new(&c, rect((10, 90), (10, 90)), rect((5, 5), (-5, 5))), Err(LayoutError::EmptyWorld));
    assert_eq!(Canvas::new(Q::zero(), Q::one()), Err(LayoutError::EmptyCanvas));
    let tiny = |bits: u32| Q::new(BigInt::from(1), BigInt::from(1) << bits);
    let world = |w: Q| [interval(Q::zero(), w.clone()), interval(Q::zero(), w)];
    assert!(Panel::new(&c, rect((0, 50), (0, 50)), world(tiny(MIN_WORLD_BITS))).is_ok());
    assert_eq!(Panel::new(&c, rect((0, 50), (0, 50)), world(tiny(MIN_WORLD_BITS + 1))), Err(LayoutError::TooNarrow));
}

#[test]
fn world_corners_map_to_viewport_corners_with_y_reversed() {
    let panel = Panel::new(&canvas(), rect((10, 90), (20, 60)), rect((-5, 5), (0, 5))).unwrap();
    assert_eq!(panel.to_page(&point(-5, 5)), point(10, 20));
    assert_eq!(panel.to_page(&point(5, 0)), point(90, 60));
    assert_eq!(panel.to_page(&point(0, 1)), point(50, 52));
}

#[test]
fn deep_magnifications_keep_exact_differences_at_large_offsets() {
    // A world window of width 2^-50 around (10^30, -10^30).
    let big = Q::new(BigInt::from(10).pow(30), BigInt::from(1));
    let width = Q::new(BigInt::from(1), BigInt::from(1) << 50);
    let world = [interval(big.clone(), &big + &width), interval(-&big, &width - &big)];
    let panel = Panel::new(&canvas(), rect((0, 64), (0, 64)), world).unwrap();
    let step = Q::new(BigInt::from(1), BigInt::from(1) << 55);
    let a = [&big + &step, &step - &big];
    let b = [&big + &(&step + &step), &step - &big];
    let (pa, pb) = (panel.to_page(&a), panel.to_page(&b));
    assert_eq!(&pb[0] - &pa[0], frac(2, 1));
    assert_eq!(pa[1], frac(62, 1));
}

/// `p + t(q - p)` for rational `t`.
fn along(p: &Point, q: &Point, t: &Q) -> Point {
    [&p[0] + &(&(&q[0] - &p[0]) * t), &p[1] + &(&(&q[1] - &p[1]) * t)]
}

fn inside(r: &Rectangle, p: &Point) -> bool {
    r[0].contains(&p[0]) && r[1].contains(&p[1])
}

#[test]
fn segments_are_cut_to_exactly_their_part_inside() {
    let mut random = Random::new(21);
    let r = rect((-3, 4), (-2, 5));
    let coordinate = |random: &mut Random| frac(random.integer(-12, 12), random.integer(1, 3));
    for case in 0..500 {
        let p = [coordinate(&mut random), coordinate(&mut random)];
        let q = if case % 9 == 0 { p.clone() } else { [coordinate(&mut random), coordinate(&mut random)] };
        let cut = clip_segment([p.clone(), q.clone()], &r);
        // Sample the segment at 65 parameters: a sample is inside the
        // rectangle exactly when it lies on the cut part.
        for k in 0..=64 {
            let t = frac(k, 64);
            let x = along(&p, &q, &t);
            let on_cut = match &cut {
                None => false,
                Some([a, b]) => {
                    let onto = |z: &Point| -> bool {
                        // x between a and b on the same line.
                        let dx = [&b[0] - &a[0], &b[1] - &a[1]];
                        let dz = [&z[0] - &a[0], &z[1] - &a[1]];
                        &dx[0] * &dz[1] == &dx[1] * &dz[0]
                            && (&dz[0] * &dx[0]) + (&dz[1] * &dx[1]) >= Q::zero()
                            && (&dz[0] * &dz[0]) + (&dz[1] * &dz[1]) <= (&dx[0] * &dx[0]) + (&dx[1] * &dx[1])
                    };
                    onto(&x)
                }
            };
            assert_eq!(inside(&r, &x), on_cut, "case {case} t {t}");
        }
        if let Some([a, b]) = &cut {
            assert!(inside(&r, a) && inside(&r, b));
        }
    }
}

/// Nonzero winding number of a closed polygon around `p`, decided exactly
/// (points on the boundary count as inside).
fn covers(polygon: &[Point], p: &Point) -> bool {
    let n = polygon.len();
    if n == 0 {
        return false;
    }
    let mut winding = 0i64;
    for k in 0..n {
        let (a, b) = (&polygon[k], &polygon[(k + 1) % n]);
        let side = (&(&b[0] - &a[0]) * &(&p[1] - &a[1])) - (&(&b[1] - &a[1]) * &(&p[0] - &a[0]));
        let between = |u: &Q, v: &Q, w: &Q| (u <= w && w <= v) || (v <= w && w <= u);
        if side == Q::zero() && between(&a[0], &b[0], &p[0]) && between(&a[1], &b[1], &p[1]) {
            return true;
        }
        if a[1] <= p[1] {
            if b[1] > p[1] && side > Q::zero() {
                winding += 1;
            }
        } else if b[1] <= p[1] && side < Q::zero() {
            winding -= 1;
        }
    }
    winding != 0
}

#[test]
fn polygons_are_cut_to_their_intersection_with_the_rectangle() {
    let mut random = Random::new(22);
    let r = rect((-3, 4), (-2, 5));
    for case in 0..200 {
        let count = random.integer(3, 7) as usize;
        let polygon: Vec<Point> = (0..count)
            .map(|_| [frac(random.integer(-10, 10), 1), frac(random.integer(-10, 10), 1)])
            .collect();
        let cut = clip_polygon(polygon.clone(), &r);
        // The cut covers exactly the covered points of the rectangle, for
        // concave and self-intersecting polygons alike.
        for x in -8..=8 {
            for y in -8..=8 {
                let p = [frac(2 * x + 1, 2), frac(2 * y + 1, 2)];
                assert_eq!(covers(&cut, &p), inside(&r, &p) && covers(&polygon, &p), "case {case} at {p:?}");
            }
        }
        assert!(cut.iter().all(|p| inside(&r, p)));
    }
    // Entirely outside, entirely inside, and degenerate polygons.
    assert!(clip_polygon(vec![point(10, 10), point(11, 10), point(10, 11)], &r).is_empty());
    let within = vec![point(0, 0), point(1, 0), point(0, 1)];
    assert_eq!(clip_polygon(within.clone(), &r), within);
    assert_eq!(clip_polygon(vec![point(0, 0)], &r), vec![point(0, 0)]);
}

#[test]
fn placing_cuts_far_geometry_to_the_padded_viewport() {
    let panel = Panel::new(&canvas(), rect((10, 90), (10, 90)), rect((-5, 5), (-5, 5))).unwrap();
    let far = Q::new(BigInt::from(10).pow(40), BigInt::from(1));
    let shapes = vec![
        Shape::Segment([[-&far, Q::zero()], [far.clone(), Q::zero()]]),
        Shape::Polygon(vec![[-&far, -&far], [far.clone(), -&far], [Q::zero(), far.clone()]]),
        Shape::Segment([point(100, 100), point(101, 100)]),
    ];
    let placed = panel.place(&shapes, &frac(1, 1));
    assert_eq!(placed.len(), 2, "the shape outside is dropped");
    let bound = rect((9, 91), (9, 91));
    for shape in &placed {
        let points = match shape {
            Shape::Segment(p) => p.to_vec(),
            Shape::Polygon(p) => p.clone(),
            Shape::Circle { .. } => unreachable!("no circles were placed"),
        };
        assert!(points.iter().all(|p| inside(&bound, p)));
    }
    assert_eq!(placed[0], Shape::Segment([point(9, 50), point(91, 50)]));
}

#[test]
fn inset_windows_frames_and_connectors_are_at_true_size() {
    let c = canvas();
    let overview = Panel::new(&c, rect((10, 90), (10, 90)), rect((-5, 5), (-5, 5))).unwrap();
    let detail = Panel::new(&c, rect((110, 170), (20, 80)), rect((1, 3), (2, 4))).unwrap();
    let pairs = [(Corner::TopRight, Corner::TopLeft), (Corner::BottomRight, Corner::BottomLeft)];
    let z = inset(&overview, &detail, &pairs, &frac(16, 1)).unwrap();
    // World x ∈ [1, 3] is page x ∈ [58, 74]; world y ∈ [2, 4] is page y ∈ [18, 34].
    assert_eq!(z.window, Shape::Polygon(vec![point(58, 18), point(74, 18), point(74, 34), point(58, 34)]));
    assert_eq!(z.frame, Shape::Polygon(vec![point(110, 20), point(170, 20), point(170, 80), point(110, 80)]));
    assert_eq!(
        z.connectors,
        vec![Shape::Segment([point(74, 18), point(110, 20)]), Shape::Segment([point(74, 34), point(110, 80)])]
    );
    let outside = Panel::new(&c, rect((110, 170), (20, 80)), rect((4, 6), (2, 4))).unwrap();
    assert_eq!(inset(&overview, &outside, &[], &Q::zero()), Err(LayoutError::WindowOutsideOverview));
    // A window on the overview's edge is inside.
    let edge = Panel::new(&c, rect((110, 170), (20, 80)), rect((3, 5), (3, 5))).unwrap();
    assert!(inset(&overview, &edge, &[], &Q::zero()).is_ok());
    // A window narrower than the mark in both directions is marked by a
    // square of the mark's side around its centre (page (66, 26)); one side
    // as long as the mark keeps the true window.
    let marked = inset(&overview, &detail, &pairs, &frac(20, 1)).unwrap();
    assert_eq!(marked.window, Shape::Polygon(vec![point(56, 16), point(76, 16), point(76, 36), point(56, 36)]));
    assert_eq!(marked.connectors[0], Shape::Segment([point(76, 16), point(110, 20)]));
    assert_eq!(inset(&overview, &detail, &pairs, &frac(16, 1)).unwrap().window, z.window);
}

mod review;

#[test]
fn circles_are_scaled_and_dropped_only_when_they_miss_the_window() {
    let canvas = Canvas::new(frac(100, 1), frac(100, 1)).unwrap();
    let panel = Panel::new(&canvas, rect((10, 90), (10, 90)), rect((0, 8), (0, 8))).unwrap();
    let circle = |x: i64, y: i64| Shape::Circle { centre: [frac(x, 1), frac(y, 1)], radius: frac(1, 2) };
    let placed = panel.place(&[circle(4, 4), circle(-2, 4)], &frac(1, 1));
    assert_eq!(placed, vec![Shape::Circle { centre: point(50, 50), radius: frac(5, 1) }]);
}
