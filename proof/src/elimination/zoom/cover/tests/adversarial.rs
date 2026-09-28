//! Adversarial tests of covers: every pinned
//! cover against an independent construction from the definitions (names,
//! order, roots, factors, maps, covered sets, zero sets), the inclusion
//! test against the exact 32-corner decision on adversarial boxes (cover
//! boxes themselves, faces, corners, zero widths, excesses of 2⁻²⁰⁰, views
//! beyond `s = −2`, far beyond D), and parameters whose families are empty.
use crate::arithmetic::exact::{frac, q, Interval, QSqrt5, Q};
use crate::arithmetic::polynomial::Polynomial;
use crate::elimination::zoom::cover::{ZoomCover, Face, Kind, Number, Parameters, Rational, Shape};
use crate::elimination::zoom::tests::independent::{self as model, rat, Family, ModelZoomCover, P5};
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::Coordinates;
use crate::elimination::proof::catalogue::{exotic_pin, local_pin, EXOTIC, LOCAL_COUNT};
use crate::problem::configuration::ConfigurationBox;
use crate::problem::geometry;

fn all_pins() -> Vec<(String, Parameters)> {
    let mut out = Vec::new();
    for name in EXOTIC {
        out.push((name.to_string(), exotic_pin(name).unwrap().parameters));
    }
    for n in 0..LOCAL_COUNT {
        out.push((format!("local {n}"), local_pin(n).unwrap().parameters));
    }
    out
}

fn identity_map() -> [[Number; 5]; 5] {
    std::array::from_fn(|i| std::array::from_fn(|j| Number(QSqrt5::integer(i64::from(i == j)))))
}

fn zero_centre() -> [Number; 5] {
    std::array::from_fn(|_| Number(QSqrt5::zero()))
}

/// DEFECT (soundness, latent): with an empty offset the family has no zoom,
/// so nothing is verified, yet the covered set is the whole base box and
/// `contains` accepts boxes in it. The coverage lemmas need an offset face.
#[test]
fn an_cover_without_offset_is_refused() {
    let unit = || [Rational(q(0)), Rational(q(1))];
    let tube = Parameters {
        coordinates: Coordinates::Configuration,
        centre: zero_centre(),
        map: identity_map(),
        shape: Shape::Tube { base: vec![unit(), unit(), unit(), unit(), unit()], offset: vec![], radius: Rational(q(1)) },
        beyond: None,
    };
    let point = Parameters {
        shape: Shape::Point {
            base: vec![[0, 1]; 5],
            offset: vec![],
            radii: vec![Rational(q(1)); 5],
            radius: Rational(q(1)),
            ratio: Rational(q(1)),
            window: None,
        },
        ..tube.clone()
    };
    for p in [tube, point] {
        if let Ok(cover) = ZoomCover::new(p.clone()) {
            // What the defect allows: no zooms, and a box of B₀ with every
            // kind of configuration (fits are not excluded) is accepted.
            let accepted = cover.contains(&ConfigurationBox::new(std::array::from_fn(|_| {
                Interval::new(frac(1, 8), frac(1, 4)).unwrap()
            })));
            panic!("accepted {:?} with {} zooms; contains a box: {accepted}", p.shape, cover.zooms().len());
        }
    }
}

fn kind_matches(kind: Kind, zoom: &model::ModelZoom) -> bool {
    let face = |f: Face| (f.axis, f.side);
    match kind {
        Kind::Tube { offset } => zoom.family == Family::Tube && face(offset) == zoom.offset_face,
        Kind::A { base, offset } => zoom.family == Family::A && Some(face(base)) == zoom.base_face && face(offset) == zoom.offset_face,
        Kind::B { offset } => zoom.family == Family::B && face(offset) == zoom.offset_face,
        Kind::Sheared { base, offset } => {
            zoom.family == Family::Sheared && Some(face(base)) == zoom.base_face && face(offset) == zoom.offset_face
        }
    }
}

#[test]
fn pinned_covers_equal_the_independent_construction() {
    for (name, p) in all_pins() {
        let cover = ZoomCover::new(p.clone()).unwrap();
        let expected = ModelZoomCover::new(&p);
        assert_eq!(cover.zooms().len(), expected.zooms.len(), "{name}");
        for (c, (zoom, e)) in cover.zooms().iter().zip(&expected.zooms).enumerate() {
            assert_eq!(zoom.name(), e.name, "{name}");
            assert!(kind_matches(cover.kind(c), e), "{name} {}", e.name);
            for j in 0..5 {
                assert_eq!((zoom.root()[j].lo(), zoom.root()[j].hi()), (&e.root[j].0, &e.root[j].1), "{name} {} root {j}", e.name);
                assert_eq!(zoom.is_scale(j), e.factors.contains(&j), "{name} {} factor {j}", e.name);
            }
            assert_eq!(zoom.map(), &e.x, "{name} {}", e.name);
            assert_eq!(zoom.coordinates(), p.coordinates);
            // Zero sets, independently: substituting y_d = 0 kills the no-fit
            // coordinates for every scale variable d.
            for &d in &e.factors {
                let images: [Polynomial; 5] =
                    std::array::from_fn(|j| if j == d { Polynomial::zero() } else { Polynomial::variables()[j].clone() });
                for i in model::no_fit_coordinates(p.coordinates) {
                    assert!(e.x[i].substitute(&images).unwrap().is_zero(), "{name} {} coordinate {i}", e.name);
                }
                assert!(e.root[d].0 >= Q::zero());
            }
            // Views multilinear (validity at corners); coordinates too.
            let (u, _) = model::view_rotation(p.coordinates, &e.x);
            assert!(u.iter().all(model::is_multilinear) && e.x.iter().all(model::is_multilinear), "{name}");
        }
        for i in 0..5 {
            assert_eq!(cover.bounds()[i], [expected.bounds[i].0.clone(), expected.bounds[i].1.clone()], "{name} bound {i}");
        }
    }
}

/// The arc-plane lemma's finite fact, with an independent rotation
/// formula: `max_w e₂·R̂(σr_A)w = (1 + |r_A|²) max_w e₂·w = (1 + |r_A|²)(2 + √5)`.
#[test]
fn arc_plane_support_value_is_shared() {
    let (a1, a3) = model::arc_constants();
    let top = geometry::vertices().iter().map(|w| w[1].clone()).max().unwrap();
    assert_eq!(top, QSqrt5::new(q(2), q(1)));
    for s in [1i64, -1] {
        let r = [&QSqrt5::integer(s) * &a1, QSqrt5::zero(), &QSqrt5::integer(s) * &a3];
        let rr = &(&r[0] * &r[0]) + &(&r[2] * &r[2]);
        let turned_top = geometry::vertices()
            .iter()
            .map(|p| {
                // Second coordinate of (1 − |r|²)p + 2(r·p)r + 2 r×p.
                let rp = &(&r[0] * &p[0]) + &(&r[2] * &p[2]);
                let cross1 = &(&r[2] * &p[0]) - &(&r[0] * &p[2]);
                &(&(&(&QSqrt5::one() - &rr) * &p[1]) + &(&(&QSqrt5::integer(2) * &rp) * &r[1])) + &(&QSqrt5::integer(2) * &cross1)
            })
            .max()
            .unwrap();
        assert_eq!(turned_top, &(&QSqrt5::one() + &rr) * &top);
    }
}

/// A rational close to `x`.
fn near(x: &QSqrt5) -> Q {
    x.rational_part() + &(x.sqrt5_part() * frac(2_236_067_977, 1_000_000_000))
}

fn cover_z(random: &mut Random, m: &ModelZoomCover) -> P5 {
    std::array::from_fn(|i| {
        let (lo, hi) = &m.bounds[i];
        let t = match random.below(3) {
            0 => frac(random.below(3) as i64, 2),
            _ => random.between(&q(0), &q(1), 8),
        };
        lo + &(&(hi - lo) * &rat(&t))
    })
}

fn tiny() -> Q {
    Q::one() / Q::from_integer(num_bigint::BigInt::from(2).pow(200))
}

fn check(name: &str, cover: &ZoomCover, m: &ModelZoomCover, b: &ConfigurationBox) -> bool {
    let decided = cover.contains(b);
    assert_eq!(decided, m.box_inside(b), "{name}: {:?}", b.axes());
    decided
}

#[test]
fn inclusion_agrees_with_the_corner_decision_on_adversarial_boxes() {
    let mut random = Random::new(0xad5e);
    for (name, p) in all_pins() {
        let cover = ZoomCover::new(p.clone()).unwrap();
        let m = ModelZoomCover::new(&p);
        let identity = p.coordinates == Coordinates::Configuration && p.map == identity_map();
        if identity {
            // The covered set itself, its faces, corners and zero-width slabs
            // are accepted; one face moved out by 2⁻²⁰⁰ is refused.
            let c = m.centre();
            let bounds: Vec<(Q, Q)> = (0..5)
                .map(|i| {
                    let lo = &m.bounds[i].0 + &c[i];
                    let hi = &m.bounds[i].1 + &c[i];
                    assert!(lo.sqrt5_part().is_zero() && hi.sqrt5_part().is_zero());
                    (lo.rational_part().clone(), hi.rational_part().clone())
                })
                .collect();
            let make = |f: &dyn Fn(usize) -> (Q, Q)| {
                ConfigurationBox::new(std::array::from_fn(|i| {
                    let (lo, hi) = f(i);
                    Interval::new(lo, hi).unwrap()
                }))
            };
            assert!(check(&name, &cover, &m, &make(&|i| bounds[i].clone())));
            for mask in 0..32u32 {
                let corner = |i: usize| if mask >> i & 1 == 1 { bounds[i].1.clone() } else { bounds[i].0.clone() };
                assert!(check(&name, &cover, &m, &make(&|i| (corner(i), corner(i)))), "{name} corner {mask}");
                // A slab of zero width on the faces of the mask, full elsewhere.
                assert!(check(&name, &cover, &m, &make(&|i| if mask >> i & 1 == 1 { (bounds[i].1.clone(), bounds[i].1.clone()) } else { bounds[i].clone() })));
            }
            for i in 0..5 {
                for side in 0..2 {
                    let b = make(&|j| {
                        let (lo, hi) = bounds[j].clone();
                        match (j == i, side) {
                            (true, 0) => (lo - tiny(), hi),
                            (true, _) => (lo, hi + tiny()),
                            _ => (lo, hi),
                        }
                    });
                    assert!(!check(&name, &cover, &m, &b), "{name} face {i}/{side} moved out");
                }
            }
        }
        // Random boxes around cover points, with widths that are zero, tiny
        // or a fraction of the cover, often straddling faces.
        let size = (0..5).map(|i| near(&(&m.bounds[i].1 - &m.bounds[i].0))).min().unwrap();
        let (mut yes, mut no) = (0, 0);
        for _ in 0..150 {
            let z = cover_z(&mut random, &m);
            let c = model::configuration_of(p.coordinates, &m.image(&z));
            let widths = |random: &mut Random| match random.below(4) {
                0 => Q::zero(),
                1 => tiny(),
                _ => &size * random.between(&q(0), &frac(1, 4), 6),
            };
            let axes: [Interval; 5] = std::array::from_fn(|j| {
                let centre = near(&c[j]);
                let (w1, w2) = (widths(&mut random), widths(&mut random));
                Interval::new(&centre - &w1, &centre + &w2).unwrap()
            });
            if check(&name, &cover, &m, &ConfigurationBox::new(axes)) {
                yes += 1;
            } else {
                no += 1;
            }
        }
        assert!(yes > 10 && no > 10, "{name}: {yes} accepted, {no} refused");
        // Views at or beyond s = −2 are refused by arc-plane covers without
        // panicking, and never accepted by any cover of the catalogue.
        for (lo, hi) in [(frac(-3, 1), frac(-2, 1)), (frac(-2, 1), frac(-2, 1)), (frac(-5, 2), frac(1, 1)), (frac(-201, 100), frac(-199, 100))] {
            let b = ConfigurationBox::new(std::array::from_fn(|j| {
                if j == 0 { Interval::new(lo.clone(), hi.clone()).unwrap() } else { Interval::point(q(0)) }
            }));
            assert!(!check(&name, &cover, &m, &b));
        }
    }
}

#[test]
fn boxes_far_beyond_the_pentagon_face_are_accepted_only_with_the_rest_inside() {
    let p = exotic_pin("pentagon").unwrap().parameters;
    let cover = ZoomCover::new(p.clone()).unwrap();
    let m = ModelZoomCover::new(&p);
    let phi = QSqrt5::new(frac(1, 2), frac(1, 2));
    let mut random = Random::new(7);
    let (mut yes, mut no) = (0, 0);
    for _ in 0..300 {
        // η inside (or slightly outside), ξ from inside the cover to far
        // beyond the face ξ = 0, r small.
        let mut z = cover_z(&mut random, &m);
        z[1] = rat(&random.between(&frac(-1, 64), &frac(1, 2), 8));
        if random.below(4) == 0 {
            z[0] = rat(&random.between(&frac(-1, 16), &frac(1, 16), 8));
        }
        let c = model::configuration_of(p.coordinates, &m.image(&z));
        let axes: [Interval; 5] = std::array::from_fn(|j| {
            let centre = near(&c[j]);
            let w = random.between(&q(0), &frac(1, 64), 5);
            Interval::new(&centre - &w, &centre + &w).unwrap()
        });
        let b = ConfigurationBox::new(axes);
        if check("pentagon", &cover, &m, &b) {
            yes += 1;
            // Every accepted corner is in the covered set or violates D's
            // inequality 0 (φs + φ²t > 1).
            for mask in 0..32u8 {
                let x = b.corner(mask).map(|v| rat(&v));
                let zc = m.cover_coordinates(&x);
                if zc[1] > m.bounds[1].1 {
                    let q0 = &(&(&phi * &x[0]) + &(&(&phi * &phi) * &x[1])) - &QSqrt5::one();
                    assert_eq!(q0.sign(), std::cmp::Ordering::Greater);
                }
            }
        } else {
            no += 1;
        }
    }
    assert!(yes > 20 && no > 20, "{yes} {no}");
}
