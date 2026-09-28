//! Tests of covers: the pinned families against the zoom formulas,
//! coverage by the constructive map, the exact inclusion test, the window
//! hand-over and the refusal of invalid parameters.
use super::*;
use crate::arithmetic::exact::{frac, q};
use crate::elimination::zoom::tests::Random;
use crate::elimination::zoom::Coordinates;
use crate::elimination::proof::catalogue::{exotic_pin, local_pin, EXOTIC, LOCAL_COUNT};

fn pinned() -> Vec<(String, Parameters)> {
    let mut out = Vec::new();
    for name in EXOTIC {
        out.push((name.to_string(), exotic_pin(name).unwrap().parameters));
    }
    for n in [0, 22, LOCAL_COUNT - 1] {
        out.push((format!("local {n}"), local_pin(n).unwrap().parameters));
    }
    out
}

/// The pinned Exotic cover `name`; for the arc covers, `i` = 0, 1
/// picks the sign `+`, `-`.
fn cover(name: &str, i: usize) -> ZoomCover {
    let name = if matches!(name, "square" | "pentagon") { name.to_string() } else { format!("{name}{}", ["+", "-"][i]) };
    ZoomCover::new(exotic_pin(&name).unwrap().parameters).unwrap()
}

#[test]
fn pinned_covers_have_the_expected_zooms() {
    let count = |a: &ZoomCover, f: fn(&Kind) -> bool| (0..a.zooms().len()).filter(|&c| f(&a.kind(c))).count();
    let expect = |name: &str, i: usize, a: usize, b: usize, t: usize, s: usize| {
        let cover = cover(name, i);
        assert_eq!(count(&cover, |k| matches!(k, Kind::A { .. })), a, "{name}");
        assert_eq!(count(&cover, |k| matches!(k, Kind::B { .. })), b, "{name}");
        assert_eq!(count(&cover, |k| matches!(k, Kind::Tube { .. })), t, "{name}");
        assert_eq!(count(&cover, |k| matches!(k, Kind::Sheared { .. })), s, "{name}");
    };
    expect("square", 0, 12, 6, 0, 0);
    expect("pentagon", 0, 18, 6, 0, 0);
    for i in 0..2 {
        expect("arc", i, 0, 0, 7, 0);
        expect("endpoint", i, 14, 7, 0, 0);
        expect("crossing", i, 14, 7, 0, 7);
    }
    let names: Vec<String> = cover("crossing", 0).zooms().iter().map(|c| c.name().to_string()).collect();
    assert_eq!(names[0], "A b0- o0+");
    assert_eq!(names[14], "B o0+");
    assert_eq!(names[27], "A' b0+ o3+");
    let local = ZoomCover::new(local_pin(7).unwrap().parameters).unwrap();
    let names: Vec<&str> = local.zooms().iter().map(|c| c.name()).collect();
    assert_eq!(names, ["T o0-", "T o0+", "T o1-", "T o1+", "T o2-", "T o2+"]);
}

fn number(x: &Number) -> QSqrt5 {
    x.0.clone()
}

/// `centre + map · z`, computed pointwise.
fn image(p: &Parameters, z: &[QSqrt5]) -> Five<QSqrt5> {
    std::array::from_fn(|i| {
        (0..5).fold(number(&p.centre[i]), |sum, j| &sum + &(&number(&p.map[i][j]) * &z[j]))
    })
}

/// The unit vector on `face` whose free coordinates are `free`.
fn unit(len: usize, face: Face, free: &mut impl Iterator<Item = QSqrt5>) -> Vec<QSqrt5> {
    (0..len)
        .map(|j| if j == face.axis { QSqrt5::integer(i64::from(face.side)) } else { free.next().unwrap() })
        .collect()
}

/// The cover coordinates of zoom variables `y`, from the zoom formulas.
fn table_coordinates(p: &Parameters, kind: Kind, y: &Five<QSqrt5>) -> Vec<QSqrt5> {
    let (k, offset) = match &p.shape {
        Shape::Tube { base, offset, .. } => (base.len(), offset.len()),
        Shape::Point { base, offset, .. } => (base.len(), offset.len()),
    };
    let scale = |c: &QSqrt5, v: &[QSqrt5]| v.iter().map(|x| c * x).collect::<Vec<_>>();
    match kind {
        Kind::Tube { offset: f } => {
            let mut rest = y[k + 1..].iter().cloned();
            let mut z = y[..k].to_vec();
            z.extend(scale(&y[k], &unit(offset, f, &mut rest)));
            z
        }
        Kind::A { base: g, offset: f } | Kind::Sheared { base: g, offset: f } => {
            let mut rest = y[2..].iter().cloned();
            let b = unit(k, g, &mut rest);
            let fh = unit(offset, f, &mut rest);
            let mut z = scale(&y[0], &b);
            let mut zf = scale(&(&y[0] * &y[1]), &fh);
            if let (Kind::Sheared { .. }, Shape::Point { window: Some(w), .. }) = (kind, &p.shape) {
                // z_offset = z'_offset − shear · z_base.
                for (l, zl) in zf.iter_mut().enumerate() {
                    for j in 0..k {
                        *zl = &*zl - &(&QSqrt5::from_rational(w.shear[l][j].0.clone()) * &z[j]);
                    }
                }
            }
            z.append(&mut zf);
            z
        }
        Kind::B { offset: f } => {
            let mut rest = y[1 + k..].iter().cloned();
            let mut z = scale(&y[0], &y[1..1 + k]);
            z.extend(scale(&y[0], &unit(offset, f, &mut rest)));
            z
        }
    }
}

fn exact_point(y: &Five<Q>) -> Five<QSqrt5> {
    std::array::from_fn(|j| QSqrt5::from_rational(y[j].clone()))
}

#[test]
fn zooms_follow_the_readme_table() {
    let mut random = Random::new(11);
    for (name, p) in pinned() {
        let cover = ZoomCover::new(p.clone()).unwrap();
        for (c, zoom) in cover.zooms().iter().enumerate() {
            let kind = cover.kind(c);
            // Scale variables: μ, ρ for A zooms, δ otherwise.
            let factors: Vec<usize> = (0..5).filter(|&j| zoom.is_scale(j)).collect();
            match kind {
                Kind::A { .. } | Kind::Sheared { .. } => assert_eq!(factors, [0, 1], "{name}"),
                Kind::B { .. } => assert_eq!(factors, [0], "{name}"),
                Kind::Tube { .. } => assert_eq!(factors.len(), 1, "{name}"),
            }
            for _ in 0..8 {
                let y = exact_point(&random.point(zoom.root(), 6));
                let z = table_coordinates(&p, kind, &y);
                assert_eq!(zoom.evaluate(&y), image(&p, &z), "{name} {}", zoom.name());
            }
        }
    }
}

/// A random point of the covered set, often on faces and ties.
fn cover_point(random: &mut Random, cover: &ZoomCover) -> Five<QSqrt5> {
    std::array::from_fn(|i| {
        let [lo, hi] = &cover.bounds()[i];
        let t = match random.below(4) {
            0 => frac(random.below(5) as i64, 4),
            _ => random.between(&q(0), &q(1), 10),
        };
        lo + &(&(hi - lo) * &QSqrt5::from_rational(t))
    })
}

#[test]
fn every_cover_point_is_covered() {
    let mut random = Random::new(12);
    for (name, p) in pinned() {
        let cover = ZoomCover::new(p.clone()).unwrap();
        let k = match &p.shape {
            Shape::Tube { base, .. } => base.len(),
            Shape::Point { base, .. } => base.len(),
        };
        let (mut zooms, mut zero) = (0, 0);
        for round in 0..400 {
            let mut z = cover_point(&mut random, &cover);
            if round % 10 == 0 {
                for v in &mut z[k..] {
                    *v = QSqrt5::zero();
                }
            }
            let x = image(&p, &z);
            assert_eq!(cover.cover_coordinates(&x), z);
            match cover.locate(&x).unwrap_or_else(|| panic!("{name}: a cover point is not located")) {
                Location::NoFit => {
                    assert!(z[k..].iter().all(|v| v.is_zero()), "{name}");
                    zero += 1;
                }
                Location::Zoom { zoom, y } => {
                    let root = cover.zooms()[zoom].root();
                    for j in 0..5 {
                        let v = &y[j];
                        assert!(QSqrt5::from_rational(root[j].lo().clone()) <= *v, "{name}");
                        assert!(*v <= QSqrt5::from_rational(root[j].hi().clone()), "{name}");
                    }
                    assert_eq!(cover.zooms()[zoom].evaluate(&y), x, "{name}");
                    assert!(!matches!(cover.kind(zoom), Kind::Sheared { .. }));
                    zooms += 1;
                }
                Location::OutsideDomain => panic!("{name}: inside the box is not beyond"),
            }
        }
        assert!(zooms > 200 && zero > 0, "{name}: {zooms} {zero}");
        // Just outside the box on every side: not located (except beyond).
        for i in 0..5 {
            for side in 0..2 {
                let mut z = cover_point(&mut random, &cover);
                let step = QSqrt5::from_rational(frac(1, 1 << 20));
                z[i] = if side == 0 { &cover.bounds()[i][0] - &step } else { &cover.bounds()[i][1] + &step };
                let located = cover.locate(&image(&p, &z));
                match p.beyond {
                    Some(b) if b.axis == i && side == 1 => assert_eq!(located, Some(Location::OutsideDomain)),
                    _ => assert_eq!(located, None, "{name} axis {i}"),
                }
            }
        }
    }
}

/// A rational close to `x`, for building boxes near irrational points.
fn near(x: &QSqrt5) -> Q {
    x.rational_part() + &(x.sqrt5_part() * frac(2_236_068, 1_000_000))
}

/// The configuration of coordinates `x` (exact, for `u₃ > 0`).
fn configuration(coordinates: Coordinates, x: &Five<QSqrt5>) -> Five<QSqrt5> {
    let constants: Five<Polynomial> = std::array::from_fn(|j| Polynomial::constant(x[j].clone()));
    let (view, rotation) = coordinates.configuration(&constants);
    let origin: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::zero());
    let u = view.vector().map(|p| p.evaluate(&origin));
    let r = rotation.map(|p| p.evaluate(&origin));
    let inverse = reciprocal(&u[2]);
    [&u[0] * &inverse, &u[1] * &inverse, r[0].clone(), r[1].clone(), r[2].clone()]
}

fn inside(cover: &ZoomCover, c: &Five<QSqrt5>) -> bool {
    let Some(x) = cover.parameters().coordinates.of_configuration(c) else { return false };
    let z = cover.cover_coordinates(&x);
    (0..5).all(|i| {
        let above = cover.parameters().beyond.is_some_and(|b| b.axis == i);
        cover.bounds()[i][0] <= z[i] && (above || z[i] <= cover.bounds()[i][1])
    })
}

#[test]
fn inclusion_test_is_exact() {
    let mut random = Random::new(13);
    for (name, p) in pinned() {
        let cover = ZoomCover::new(p.clone()).unwrap();
        let (mut accepted, mut refused) = (0, 0);
        let size = (0..5).map(|i| near(&(&cover.bounds()[i][1] - &cover.bounds()[i][0]))).min().unwrap();
        for _ in 0..300 {
            let centre = configuration(p.coordinates, &image(&p, &cover_point(&mut random, &cover)));
            let scale = &size * frac(1, 1 << (1 + random.below(7)));
            let axes: Five<Interval> = std::array::from_fn(|j| {
                let c = near(&centre[j]);
                let lo = &c - &(&scale * random.between(&q(0), &q(1), 4));
                let hi = &c + &(&scale * random.between(&q(0), &q(1), 4));
                Interval::new(lo, hi).unwrap()
            });
            let b = ConfigurationBox::new(axes.clone());
            let corners: Vec<Five<QSqrt5>> = (0..32u8).map(|m| exact_point(&b.corner(m))).collect();
            let all_corners = corners.iter().all(|c| inside(&cover, c));
            assert_eq!(cover.contains(&b), all_corners, "{name}");
            if all_corners {
                accepted += 1;
                for _ in 0..8 {
                    assert!(inside(&cover, &exact_point(&random.point(&axes, 6))), "{name}");
                }
            } else {
                refused += 1;
            }
        }
        assert!(accepted > 20 && refused > 20, "{name}: {accepted} accepted, {refused} refused");
    }
}

#[test]
fn window_is_exact_and_hands_over_the_same_configuration() {
    let mut random = Random::new(14);
    for i in 0..2 {
        let cover = cover("crossing", i);
        let shape = cover.parameters().shape.clone();
        let Shape::Point { window: Some(window), offset, .. } = &shape else { panic!() };
        let (mut held, mut failed, mut handed) = (0, 0, 0);
        for (c, zoom) in cover.zooms().iter().enumerate() {
            let kind = cover.kind(c);
            let on_window = matches!(kind, Kind::A { base, .. } if window.faces.contains(&base));
            for round in 0..40 {
                // Random dyadic sub-boxes of the root, half of them around the
                // second arc's factor point ρ = 1, free offset coordinates 0.
                let depth = 1 + random.below(4) as u32;
                let near = round % 2 == 0;
                let cell = ZoomCell(std::array::from_fn(|j| {
                    let root = &zoom.root()[j];
                    if near && j > 0 {
                        let centre = if j == 1 { q(1) } else { q(0) };
                        let centre = if root.contains(&centre) { centre } else { root.lo().clone() };
                        let w = frac(1, 1 << (3 + random.below(4)));
                        let lo = std::cmp::max(root.lo().clone(), &centre - &(&w * random.between(&q(0), &q(1), 3)));
                        let hi = std::cmp::min(root.hi().clone(), &lo + &w);
                        return Interval::new(lo, hi).unwrap();
                    }
                    let k = random.below(1 << depth);
                    let w = root.width() * frac(1, 1 << depth);
                    let lo = root.lo() + &(&w * q(k as i64));
                    Interval::new(lo.clone(), &lo + &w).unwrap()
                }));
                let holds = cover.window_holds(c, &cell);
                if !on_window {
                    assert!(!holds);
                    continue;
                }
                // The window vector is multilinear in the cell's variables: it
                // stays in the window iff it does at every corner.
                let at_corners = (0..32u32).all(|mask| {
                    let y: Five<QSqrt5> = std::array::from_fn(|j| {
                        QSqrt5::from_rational(if mask >> j & 1 == 1 { cell[j].hi() } else { cell[j].lo() }.clone())
                    });
                    let (_, w) = cover.window_vector(c, &y).unwrap();
                    w.iter().zip(offset).all(|(x, &[lo, hi])| {
                        let r = &window.radius.0;
                        QSqrt5::from_rational(r * q(i64::from(lo))) <= *x && *x <= QSqrt5::from_rational(r * q(i64::from(hi)))
                    })
                });
                assert_eq!(holds, at_corners);
                if holds {
                    held += 1;
                    for _ in 0..5 {
                        let y = exact_point(&random.point(&cell, 5));
                        match cover.hand_over(c, &y).unwrap() {
                            Location::Zoom { zoom: target, y: image } => {
                                assert!(matches!(cover.kind(target), Kind::Sheared { .. }));
                                let root = cover.zooms()[target].root();
                                assert!((0..5).all(|j| QSqrt5::from_rational(root[j].lo().clone()) <= image[j]
                                    && image[j] <= QSqrt5::from_rational(root[j].hi().clone())));
                                assert_eq!(cover.zooms()[target].evaluate(&image), zoom.evaluate(&y));
                                handed += 1;
                            }
                            Location::NoFit => {}
                            Location::OutsideDomain => panic!(),
                        }
                    }
                } else {
                    failed += 1;
                }
            }
        }
        assert!(held > 10 && failed > 10 && handed > 20, "{held} {failed} {handed}");
    }
}

#[test]
fn pentagon_cover_coordinates_are_the_event_coordinates() {
    let cover = cover("pentagon", 0);
    let phi = QSqrt5::new(frac(1, 2), frac(1, 2));
    let inverse_phi = &phi - &QSqrt5::one();
    let mut random = Random::new(15);
    for _ in 0..100 {
        let c: Five<QSqrt5> = std::array::from_fn(|_| QSqrt5::from_rational(random.small(40)));
        let z = cover.cover_coordinates(&c);
        assert_eq!(z[0], &c[0] - &(&c[1] * &inverse_phi));
        assert_eq!(z[1], &(&c[0] + &(&phi * &c[1])) - &inverse_phi);
        assert_eq!(z[2..], c[2..]);
        // Beyond ξ = 0 the triangle inequality of D is violated: q₀ = φξ.
        let q0 = &(&(&phi * &c[0]) + &(&(&phi * &phi) * &c[1])) - &QSqrt5::one();
        assert_eq!(q0, &phi * &z[1]);
    }
}

fn square_parameters() -> Parameters {
    exotic_pin("square").unwrap().parameters
}

#[test]
fn invalid_parameters_are_refused() {
    let refused = |change: &dyn Fn(&mut Parameters), expected: CoverError| {
        let mut p = square_parameters();
        change(&mut p);
        assert_eq!(ZoomCover::new(p).unwrap_err(), expected);
    };
    refused(&|p| p.map[3] = p.map[2].clone(), CoverError::Singular);
    refused(&|p| if let Shape::Point { base, .. } = &mut p.shape { base[0] = [0, 2] }, CoverError::UnitBox);
    refused(&|p| if let Shape::Point { offset, .. } = &mut p.shape { offset[1] = [1, 1] }, CoverError::UnitBox);
    refused(&|p| if let Shape::Point { offset, .. } = &mut p.shape { offset.push([-1, 1]) }, CoverError::Dimensions);
    refused(&|p| if let Shape::Point { radii, .. } = &mut p.shape { radii.pop(); }, CoverError::Dimensions);
    refused(&|p| if let Shape::Point { radius, .. } = &mut p.shape { radius.0 = q(0) }, CoverError::NotPositive);
    refused(&|p| if let Shape::Point { ratio, .. } = &mut p.shape { ratio.0 = frac(-1, 2) }, CoverError::NotPositive);
    refused(&|p| if let Shape::Point { radii, .. } = &mut p.shape { radii[1].0 = q(0) }, CoverError::NotPositive);
    let window = |faces: Vec<Face>, rows: usize| Window {
        faces,
        shear: vec![vec![Rational(q(0)); 2]; rows],
        radius: Rational(frac(1, 4)),
    };
    refused(&|p| if let Shape::Point { window: w, .. } = &mut p.shape { *w = Some(window(vec![Face { axis: 0, side: -1 }], 3)) }, CoverError::WindowFace);
    refused(&|p| if let Shape::Point { window: w, .. } = &mut p.shape { *w = Some(window(vec![Face { axis: 0, side: 1 }], 2)) }, CoverError::Dimensions);
    // Repeated window faces (duplicate sheared zooms and names) and a window
    // without faces.
    let twice = vec![Face { axis: 0, side: 1 }, Face { axis: 0, side: 1 }];
    refused(&|p| if let Shape::Point { window: w, .. } = &mut p.shape { *w = Some(window(twice.clone(), 3)) }, CoverError::WindowFace);
    refused(&|p| if let Shape::Point { window: w, .. } = &mut p.shape { *w = Some(window(Vec::new(), 3)) }, CoverError::WindowFace);
    assert!(ZoomCover::new({
        let mut p = square_parameters();
        if let Shape::Point { window: w, .. } = &mut p.shape {
            *w = Some(window(vec![Face { axis: 0, side: 1 }, Face { axis: 1, side: 1 }], 3));
        }
        p
    })
    .is_ok());
    // The square's centre moved off the aligned set: the zero sets leave it.
    refused(&|p| p.centre[2] = Number(QSqrt5::from_rational(frac(1, 1000))), CoverError::Zoom {
        zoom: "A b0+ o0-".into(),
        error: ZoomError::ZeroSet { variable: 0, coordinate: 2 },
    });
    // Beyond: only a positive multiple of the distance beyond the end.
    refused(&|p| p.beyond = Some(Beyond { axis: 1, inequality: 0 }), CoverError::Beyond);
    let mut pentagon = exotic_pin("pentagon").unwrap().parameters;
    assert!(ZoomCover::new(pentagon.clone()).is_ok());
    for beyond in [Beyond { axis: 0, inequality: 0 }, Beyond { axis: 1, inequality: 1 }, Beyond { axis: 5, inequality: 0 }, Beyond { axis: 1, inequality: 73 }] {
        pentagon.beyond = Some(beyond);
        assert_eq!(ZoomCover::new(pentagon.clone()).unwrap_err(), CoverError::Beyond, "{beyond:?}");
    }
    // A tube with a flat or reversed base interval.
    let mut arc = exotic_pin("arc+").unwrap().parameters;
    if let Shape::Tube { base, .. } = &mut arc.shape {
        base[0] = [Rational(frac(1, 8)), Rational(frac(1, 8))];
    }
    assert_eq!(ZoomCover::new(arc).unwrap_err(), CoverError::NotPositive);
}

/// Adversarial tests.
mod adversarial;

/// The permission to delegate, in isolation: every delegated cell of the
/// crossing's data (where the window holds on an unsheared A zoom of a
/// window face) is refused on every other zoom of its cover (A zooms of the
/// other base face, B zooms, and the sheared zooms, which receive
/// delegations and never delegate).
#[test]
fn only_a_zooms_on_window_faces_may_delegate() {
    use crate::elimination::zoom::tree;
    use crate::elimination::proof::catalogue::exotic_data;
    use crate::elimination::proof::format::{Leaf, CoverFile};
    let mut delegated = 0;
    for name in ["crossing+", "crossing-"] {
        let entry = CoverFile::parse(exotic_data(name).unwrap()).unwrap();
        let cover = ZoomCover::new(entry.parameters.clone()).unwrap();
        let Shape::Point { window: Some(window), .. } = &cover.parameters().shape else { panic!() };
        let permitted = |c: usize| matches!(cover.kind(c), Kind::A { base, .. } if window.faces.contains(&base));
        assert!((0..cover.zooms().len()).any(|c| matches!(cover.kind(c), Kind::Sheared { .. })));
        for (c, zoom) in entry.zooms.iter().enumerate() {
            let cells = tree::leaves(cover.zooms()[c].root(), &zoom.tree).unwrap();
            for (cell, leaf) in cells.iter().zip(&zoom.leaves) {
                if *leaf != Leaf::Delegated {
                    continue;
                }
                assert!(permitted(c) && cover.window_holds(c, cell));
                delegated += 1;
                for other in (0..cover.zooms().len()).filter(|&o| !permitted(o)) {
                    assert!(!cover.window_holds(other, cell), "{} delegates", cover.zooms()[other].name());
                }
            }
        }
    }
    assert!(delegated > 0);
}
