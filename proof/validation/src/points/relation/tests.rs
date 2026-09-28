use super::*;
use crate::points::catalogue::Catalogue;
use crate::testing::{shipped, Rng};
use rid::components::domain::Domain;
use rid::components::Component;
use rid::arithmetic::exact::{frac, q, Q};
use rid::problem::configuration::{ConfigurationBox, AXES};
use rid::problem::geometry::symmetries;

fn centre(values: [QSqrt5; AXES]) -> Centre {
    Centre::new(values).unwrap()
}

fn rational(values: [Q; AXES]) -> Centre {
    centre(values.map(QSqrt5::from_rational))
}

fn planar(x: i64, y: i64) -> Planar {
    [QSqrt5::integer(x), QSqrt5::integer(y)]
}

fn scaled(points: &[Planar], factor: &Q) -> Vec<Planar> {
    points.iter().map(|p| [p[0].scale(factor), p[1].scale(factor)]).collect()
}

/// Whether `p` lies in the triangle `abc` (closed), by orientations.
fn in_triangle(a: &Planar, b: &Planar, c: &Planar, p: &Planar) -> bool {
    let signs = [orientation(a, b, p), orientation(b, c, p), orientation(c, a, p)]
        .map(|o| o.sign());
    !(signs.contains(&Ordering::Greater) && signs.contains(&Ordering::Less))
}

/// Closed containment in the hull by Carathéodory: a point of the convex
/// hull of plane points lies in a triangle of three of them. Independent of
/// the edge-line argument of `classify`.
fn in_hull(hole: &[Planar], p: &Planar) -> bool {
    let n = hole.len();
    (0..n).any(|a| {
        (a + 1..n).any(|b| (b + 1..n).any(|c| in_triangle(&hole[a], &hole[b], &hole[c], p)))
    })
}

#[test]
fn synthetic_polygons_are_classified_exactly() {
    // A square with a repeated vertex and collinear edge points.
    let square = vec![
        planar(-2, -2),
        planar(2, -2),
        planar(2, 2),
        planar(-2, 2),
        planar(0, -2),
        planar(2, 0),
        planar(2, 2),
    ];
    assert_eq!(classify(&square, &square), Relation::Touch);
    assert_eq!(classify(&square, &scaled(&square, &frac(1, 2))), Relation::Fit);
    assert_eq!(classify(&square, &scaled(&square, &frac(1001, 1000))), Relation::Poke);
    assert_eq!(classify(&square, &[planar(0, 0), planar(1, 2)]), Relation::Touch);
    assert_eq!(classify(&square, &[planar(0, 0), planar(3, 0)]), Relation::Poke);
    // Irrational coordinates just inside and just outside an edge.
    let tiny = QSqrt5::new(q(161), q(-72)); // 1/(161 + 72√5) > 0
    let inside = [&QSqrt5::integer(2) - &tiny, QSqrt5::zero()];
    let outside = [&QSqrt5::integer(2) + &tiny, QSqrt5::zero()];
    assert_eq!(classify(&square, &[inside]), Relation::Fit);
    assert_eq!(classify(&square, &[outside]), Relation::Poke);
}

#[test]
#[should_panic(expected = "collinear")]
fn a_flat_hole_is_refused() {
    let flat = vec![planar(0, 0), planar(1, 1), planar(2, 2)];
    classify(&flat, &[planar(0, 0)]);
}

#[test]
fn random_polygons_agree_with_the_triangle_test() {
    let mut rng = Rng::new(7);
    for _ in 0..300 {
        let hole: Vec<Planar> = (0..3 + rng.below(6))
            .map(|_| [QSqrt5::from_rational(rng.small(4)), QSqrt5::from_rational(rng.small(4))])
            .collect();
        let n = hole.len();
        let flat = (0..n).all(|a| {
            (0..n).all(|b| (0..n).all(|c| orientation(&hole[a], &hole[b], &hole[c]).is_zero()))
        });
        if flat {
            continue;
        }
        let plug: Vec<Planar> = (0..1 + rng.below(4))
            .map(|_| [QSqrt5::from_rational(rng.small(4)), QSqrt5::from_rational(rng.small(4))])
            .collect();
        let weakly = plug.iter().all(|p| in_hull(&hole, p));
        let relation = classify(&hole, &plug);
        assert_eq!(relation != Relation::Poke, weakly);
        // Shrinking the hole around an interior point never gains containment,
        // and the hole contains its own shrunk copy strictly.
        let centroid: Planar = [0, 1].map(|j| {
            let sum = hole.iter().fold(QSqrt5::zero(), |s, h| &s + &h[j]);
            sum.scale(&Q::new((1).into(), (hole.len() as i64).into()))
        });
        let shrunk: Vec<Planar> = hole
            .iter()
            .map(|h| [0, 1].map(|j| &centroid[j] + &(&h[j] - &centroid[j]).scale(&frac(1, 2))))
            .collect();
        assert_eq!(classify(&hole, &shrunk), Relation::Fit);
        assert_eq!(classify(&hole, &hole), Relation::Touch);
    }
}

#[test]
fn every_symmetry_of_the_solid_gives_identical_shadows() {
    // A rotation g = (w, v) of the group with w ≠ 0 has the Cayley vector
    // v / w and maps K onto itself, so at every view the plug shadow is the
    // hole shadow: a touch, never a fit or a poke. Irrational rotations,
    // most of them outside the root box.
    let mut rng = Rng::new(8);
    let mut tested = 0;
    for g in symmetries() {
        let w = &g[0];
        if w.is_zero() {
            continue;
        }
        let inverse = QSqrt5::new(w.rational_part() / &w.norm(), -(w.sqrt5_part() / &w.norm()));
        assert_eq!(w * &inverse, QSqrt5::one());
        let r: [QSqrt5; 3] = std::array::from_fn(|j| &g[j + 1] * &inverse);
        for _ in 0..2 {
            let view = [rng.number(5), rng.number(5)];
            let (hole, plug) = shadows_at(view, r.clone());
            assert_eq!(classify(&hole, &plug), Relation::Touch);
        }
        tested += 1;
    }
    assert_eq!(tested, 45, "the 60 rotations minus the 15 half-turns");
    // The identity at an irrational view of the root box.
    let s = QSqrt5::new(frac(-2, 1), q(1)); // √5 - 2
    let t = QSqrt5::new(frac(1, 4), frac(-1, 20));
    let c = centre([s, t, QSqrt5::zero(), QSqrt5::zero(), QSqrt5::zero()]);
    assert_eq!(relation(&c), Relation::Touch);
}

#[test]
fn large_rotations_poke() {
    for values in [
        [frac(1, 4), frac(1, 8), frac(1, 4), q(0), q(0)],
        [frac(1, 5), frac(1, 10), frac(1, 8), frac(1, 16), frac(-1, 32)],
        [q(0), q(0), q(0), q(0), frac(2, 5)],
    ] {
        assert_eq!(relation(&rational(values)), Relation::Poke);
    }
}

#[test]
fn domain_membership_agrees_with_the_crate_box_test_at_rational_points() {
    // For a single configuration, `misses` finds an inequality violated there
    // exactly when the configuration is outside D; `positive_on` decides it
    // by its own interval rules, independently of evaluation.
    let mut rng = Rng::new(9);
    let root = ConfigurationBox::root();
    let mut outside = 0;
    for _ in 0..400 {
        let values: [Q; AXES] = std::array::from_fn(|j| {
            let a = &root.axes()[j];
            rng.rational(a.lo(), a.hi(), 6)
        });
        let point = ConfigurationBox::point(&values);
        let inside = in_domain(&rational(values));
        assert_eq!(inside, Domain.check(&point).is_none());
        outside += usize::from(!inside);
    }
    assert!(outside > 50 && outside < 390, "both sides are sampled: {outside}");
    assert!(in_domain(&rational([q(0), q(0), q(0), q(0), q(0)])));
    assert!(!in_domain(&rational([frac(2, 3), frac(2, 5), q(0), q(0), q(0)])));
}

#[test]
fn a_sample_of_the_catalogue_has_the_expected_properties() {
    let catalogue = Catalogue::load(&shipped("points.json")).unwrap();
    let find = |id: &str| &catalogue.points.iter().find(|p| p.id == id).unwrap().centre;
    for (id, domain, expected) in [
        ("square-center", true, Relation::Touch),
        ("pentagon-centre", true, Relation::Touch),
        ("crossing-r0", true, Relation::Touch),
        ("endpoint-r0", true, Relation::Touch),
        ("midpoint-r0", true, Relation::Touch),
        ("global-generic-mixed", true, Relation::Poke),
        ("global-quarter-x", false, Relation::Poke),
        ("square-second-touch-sheet", false, Relation::Touch),
        ("outside-root-corner", false, Relation::Touch),
    ] {
        assert_eq!(in_domain(find(id)), domain, "{id}");
        assert_eq!(relation(find(id)), expected, "{id}");
    }
}

/// Prints the exact properties of every catalogue point, one JSON line each:
/// `cargo test --release properties_of_the_catalogue -- --ignored --nocapture`.
#[test]
#[ignore]
fn properties_of_the_catalogue() {
    let catalogue = Catalogue::load(&shipped("points.json")).unwrap();
    for p in &catalogue.points {
        println!(
            "{}",
            serde_json::json!({"id": p.id, "in_domain": in_domain(&p.centre), "relation": relation(&p.centre)})
        );
    }
}
