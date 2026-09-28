use super::*;
use rid::arithmetic::exact::{frac, QSqrt5};
use num_bigint::BigInt;
use rid::problem::configuration::AXES;

/// The midpoint of the Global record `000000010000` of the final certificate,
/// `{"edge":[0,9],"vertex":1}`: the plug vertex 1 lies beyond the hole's
/// extent in the direction of the edge from 0 to 9.
fn global_example() -> ConfigurationScene {
    let x: [Q; AXES] = [frac(1, 24), frac(1, 40), frac(-3, 10), frac(-1, 10), frac(-3, 10)];
    ConfigurationScene::new(&x.map(QSqrt5::from_rational)).unwrap()
}

fn point_of(shape: &Shape, end: usize) -> Point {
    match shape {
        Shape::Segment(ends) => ends[end].clone(),
        Shape::Polygon(points) => points[end].clone(),
        Shape::Circle { centre, .. } => centre.clone(),
    }
}

#[test]
fn the_line_lies_at_the_holes_extent_and_the_normal_points_outward() {
    let scene = global_example();
    let s = support(&scene, [0, 9], 1, &frac(1, 2), &frac(1, 20)).unwrap();
    let (a, b) = (scene.position(Body::Hole, 0), scene.position(Body::Hole, 9));
    let n: Point = [-(&b[1] - &a[1]), &b[0] - &a[0]];
    let extent = (0..VERTEX_COUNT).map(|w| dot(&n, scene.position(Body::Hole, w))).max().unwrap();
    let base = point_of(&s.normal[0], 0);
    let tip = point_of(&s.normal[0], 1);
    assert_eq!(dot(&n, &base), extent, "the arrow starts on the line");
    // It starts at the foot of the perpendicular from the plug vertex.
    let p = scene.position(Body::Plug, 1);
    let offset: Point = [&p[0] - &base[0], &p[1] - &base[1]];
    assert_eq!(&offset[0] * &(&b[0] - &a[0]) + &offset[1] * &(&b[1] - &a[1]), q(0));
    // The silhouette segment, when drawn, lies on the line.
    for shape in &s.edge {
        for end in [point_of(shape, 0), point_of(shape, 1)] {
            assert!(&extent - &dot(&n, &end) <= Q::new(1.into(), BigInt::from(1) << 60));
        }
    }
    assert!(dot(&n, &tip) > extent, "the arrow points away from the hole");
    // Both ends of the line lie exactly at the extent: the line is
    // perpendicular to the normal, since n·d = 0.
    let (l0, l1) = (point_of(&s.line[0], 0), point_of(&s.line[0], 1));
    assert_eq!((dot(&n, &l0), dot(&n, &l1)), (extent.clone(), extent.clone()));
    // The witness: the plug vertex lies strictly beyond the line.
    assert!(dot(&n, scene.position(Body::Plug, 1)) > extent);
}

#[test]
fn the_parts_have_their_documented_shapes() {
    let s = support(&global_example(), [0, 9], 1, &frac(1, 2), &frac(1, 20)).unwrap();
    assert!(s.edge.len() <= 1);
    assert_eq!((s.line.len(), s.normal.len(), s.vertex.len()), (1, 3, 1));
    let plug = global_example().position(Body::Plug, 1).clone();
    assert!(matches!(&s.vertex[0], Shape::Circle { centre, radius } if *centre == plug && *radius == frac(1, 20)));
}

#[test]
fn invalid_supports_are_refused() {
    let scene = global_example();
    let (one, small) = (frac(1, 2), frac(1, 20));
    assert_eq!(support(&scene, [0, 60], 1, &one, &small), Err(SupportError::Vertex(60)));
    assert_eq!(support(&scene, [0, 9], 60, &one, &small), Err(SupportError::Vertex(60)));
    assert_eq!(support(&scene, [3, 3], 1, &one, &small), Err(SupportError::SameEnds));
    assert_eq!(support(&scene, [0, 9], 1, &q(0), &small), Err(SupportError::Length));
    assert_eq!(support(&scene, [0, 9], 1, &one, &frac(-1, 2)), Err(SupportError::Length));
}
