use super::*;
use crate::arithmetic::exact::QSqrt5;
use crate::problem::geometry::View;
use num_bigint::BigInt;

struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn path(&mut self, length: usize) -> String {
        (0..length)
            .map(|_| if self.next() % 2 == 0 { '0' } else { '1' })
            .collect()
    }
}

fn volume(b: &ConfigurationBox) -> Q {
    b.axes().iter().fold(q(1), |v, a| v * a.width())
}

#[test]
fn root_is_the_stated_box() {
    let root = ConfigurationBox::root();
    let expected = [
        (q(0), frac(2, 3)),
        (q(0), frac(2, 5)),
        (frac(-2, 5), frac(2, 5)),
        (frac(-2, 5), frac(2, 5)),
        (frac(-2, 5), frac(2, 5)),
    ];
    for (axis, (lo, hi)) in root.axes().iter().zip(expected) {
        assert_eq!((axis.lo().clone(), axis.hi().clone()), (lo, hi));
    }
    assert_eq!(ConfigurationBox::from_path("").unwrap(), root);
    assert_eq!([S, T, R[0], R[1], R[2]], [0, 1, 2, 3, 4]);
}

#[test]
fn cyclic_paths_are_closed_partitions_with_exact_volumes() {
    let root = ConfigurationBox::root();
    let mut layer = vec![(String::new(), root.clone())];
    for depth in 0..10 {
        let axis = split_axis(depth);
        assert_eq!(axis, depth % 5);
        let mut next = Vec::new();
        for (path, cover) in layer {
            assert_eq!(ConfigurationBox::from_path(&path).unwrap(), cover);
            assert!(root.contains(&cover));
            let [lower, upper] = cover.split(axis);
            assert_eq!(lower.axes()[axis].lo(), cover.axes()[axis].lo());
            assert_eq!(lower.axes()[axis].hi(), upper.axes()[axis].lo());
            assert_eq!(upper.axes()[axis].hi(), cover.axes()[axis].hi());
            assert_eq!(lower.axes()[axis].hi(), &cover.axes()[axis].midpoint());
            assert_eq!(volume(&lower) + volume(&upper), volume(&cover));
            assert!(cover.contains(&lower) && cover.contains(&upper));
            for other in (0..AXES).filter(|j| *j != axis) {
                assert_eq!(lower.axes()[other], cover.axes()[other]);
                assert_eq!(upper.axes()[other], cover.axes()[other]);
            }
            next.push((format!("{path}0"), lower));
            next.push((format!("{path}1"), upper));
        }
        assert_eq!(next.iter().fold(q(0), |v, (_, r)| v + volume(r)), volume(&root));
        layer = next;
    }
}

#[test]
fn random_paths_extend_by_one_split_and_contain_their_descendants() {
    let mut random = Random(5);
    for _ in 0..200 {
        let length = (random.next() % 60) as usize;
        let path = random.path(length);
        let cover = ConfigurationBox::from_path(&path).unwrap();
        for bit in 0..2 {
            let child = ConfigurationBox::from_path(&format!("{path}{bit}")).unwrap();
            assert_eq!(child, cover.split(split_axis(length))[bit].clone());
            assert!(cover.contains(&child));
        }
        let descendant = ConfigurationBox::from_path(&format!("{path}{}", random.path(7))).unwrap();
        assert!(cover.contains(&descendant));
        assert!(cover.contains_point(&descendant.midpoint()));
        assert!(cover.contains_point(&cover.midpoint()));
    }
}

#[test]
fn deep_paths_reconstruct_exact_widths() {
    let root = ConfigurationBox::root();
    let deep = ConfigurationBox::from_path(&"01".repeat(1500)).unwrap();
    let factor = Q::from_integer(BigInt::from(1) << 600);
    for j in 0..AXES {
        assert_eq!(deep.axes()[j].width(), root.axes()[j].width() / &factor);
        assert!(root.axes()[j].contains(deep.axes()[j].lo()));
    }
}

#[test]
fn invalid_paths_are_refused_with_their_first_bad_byte() {
    for (bad, position, byte) in [
        ("2", 0, b'2'),
        ("0 1", 1, b' '),
        ("01\n", 2, b'\n'),
        ("01x0", 2, b'x'),
        ("𝟎", 0, 0xf0),
    ] {
        assert_eq!(
            ConfigurationBox::from_path(bad),
            Err(PathError::InvalidByte { position, byte })
        );
    }
}

#[test]
fn paths_up_to_the_maximal_depth_are_reconstructed_and_longer_ones_refused_first() {
    let mut random = Random(8);
    let path = random.path(MAX_DEPTH);
    let deepest = ConfigurationBox::from_path(&path).unwrap();
    assert!(ConfigurationBox::root().contains(&deepest));
    assert_eq!(
        deepest,
        ConfigurationBox::from_path(&path[..MAX_DEPTH - 1])
            .unwrap()
            .split(split_axis(MAX_DEPTH - 1))[usize::from(path.as_bytes()[MAX_DEPTH - 1] - b'0')]
    );
    // A bad last byte is reported at its position.
    let mut bad = path[..MAX_DEPTH - 1].to_owned();
    bad.push('2');
    assert_eq!(
        ConfigurationBox::from_path(&bad),
        Err(PathError::InvalidByte { position: MAX_DEPTH - 1, byte: b'2' })
    );
    // One byte more is refused by its length alone, whatever the bytes.
    for longer in [format!("{path}0"), format!("{path}x"), "2".repeat(20_000)] {
        assert_eq!(
            ConfigurationBox::from_path(&longer),
            Err(PathError::TooDeep { length: longer.len() })
        );
    }
    let error = PathError::TooDeep { length: MAX_DEPTH + 1 }.to_string();
    assert!(error.contains("4097") && error.contains("4096"), "{error}");
    assert_eq!(MAX_DEPTH, crate::search::certificate::MAX_DEPTH);
}

/// Cost of reconstructing paths of the maximal depth, printed:
/// `cargo test --release --lib -- --ignored
/// problem::configuration::tests::path_cost --nocapture --test-threads 1`.
#[test]
#[ignore]
fn path_cost() {
    let mut random = Random(12);
    for length in [MAX_DEPTH / 4, MAX_DEPTH / 2, MAX_DEPTH] {
        let path = random.path(length);
        let start = std::time::Instant::now();
        ConfigurationBox::from_path(&path).unwrap();
        println!("from_path of {length} bytes: {:?}", start.elapsed());
    }
}

#[test]
fn the_path_test_agrees_with_path_parsing_and_the_order_is_depth_then_lexical() {
    // Up to the maximal depth; `is_path` does not check the length.
    for text in ["", "0", "1", "0110", "2", "01a", "0 1", "１", "0\n"] {
        assert_eq!(is_path(text), ConfigurationBox::from_path(text).is_ok(), "{text:?}");
    }
    assert!(is_path(&"0".repeat(MAX_DEPTH + 1)));
    let mut paths: Vec<String> = (0..64usize)
        .map(|n| format!("{:b}", n + 1)[1..].to_owned())
        .collect();
    paths.reverse();
    paths.sort_by(|a, b| search_order(a).cmp(&search_order(b)));
    // Breadth first: the root, then each depth in lexical order.
    assert_eq!(&paths[..7], ["", "0", "1", "00", "01", "10", "11"]);
    for pair in paths.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        assert!(a.len() < b.len() || (a.len() == b.len() && a < b), "{a} {b}");
    }
}

#[test]
fn degenerate_boxes_hold_one_value_on_their_flat_axes() {
    for axis in 0..AXES {
        let mut axes = ConfigurationBox::root().axes().clone();
        axes[axis] = Interval::point(frac(1, 7));
        let degenerate = ConfigurationBox::new(axes);
        assert_eq!(volume(&degenerate), q(0));
        assert_eq!(degenerate.split(axis)[0], degenerate);
        assert_eq!(degenerate.split(axis)[1], degenerate);
    }
    let point = [q(0), frac(1, 3), frac(-1, 5), q(0), frac(2, 5)];
    let single = ConfigurationBox::point(&point);
    assert_eq!(single.midpoint(), point);
    assert!((0..32).all(|mask| single.corner(mask) == point));
    assert!(ConfigurationBox::root().contains(&single));
    assert!(!single.contains(&ConfigurationBox::root()));
}

#[test]
fn containment_is_closed_and_exact_at_the_boundary() {
    let root = ConfigurationBox::root();
    for axis in 0..AXES {
        for (delta_lo, delta_hi, inside) in [
            (q(0), q(0), true),
            (-frac(1, 1_000_000_000), q(0), false),
            (q(0), frac(1, 1_000_000_000), false),
        ] {
            let mut axes = root.axes().clone();
            axes[axis] = Interval::new(axes[axis].lo() + delta_lo, axes[axis].hi() + delta_hi)
                .unwrap();
            let other = ConfigurationBox::new(axes);
            assert_eq!(root.contains(&other), inside);
            let corners_inside = (0..32).all(|mask| root.contains_point(&other.corner(mask)));
            assert_eq!(corners_inside, inside);
        }
    }
}

#[test]
fn corners_follow_the_mask_bits() {
    let root = ConfigurationBox::root();
    for mask in 0..32u8 {
        let corner = root.corner(mask);
        for j in 0..AXES {
            let expected = if mask & (1 << j) == 0 {
                root.axes()[j].lo()
            } else {
                root.axes()[j].hi()
            };
            assert_eq!(&corner[j], expected);
        }
    }
}

#[test]
#[should_panic(expected = "axis 5 out of range")]
fn splitting_a_missing_axis_panics() {
    ConfigurationBox::root().split(5);
}

#[test]
#[should_panic(expected = "corner mask 32 out of range")]
fn a_corner_mask_beyond_the_axes_panics() {
    ConfigurationBox::root().corner(32);
}

#[test]
fn coordinates_are_the_affine_view_and_the_rotation_variables() {
    let (view, rotation) = coordinates();
    let x = Polynomial::variables();
    assert_eq!(
        view,
        View::Affine {
            s: x[0].clone(),
            t: x[1].clone()
        }
    );
    assert_eq!(rotation, [x[2].clone(), x[3].clone(), x[4].clone()]);
    let point = [frac(1, 2), frac(1, 3), frac(1, 5), frac(1, 7), frac(1, 11)];
    let exact = point.clone().map(QSqrt5::from_rational);
    let u = view.vector();
    assert_eq!(u[0].evaluate(&exact), QSqrt5::from_rational(frac(1, 2)));
    assert_eq!(u[1].evaluate(&exact), QSqrt5::from_rational(frac(1, 3)));
    assert_eq!(u[2].evaluate(&exact), QSqrt5::one());
    for j in 0..3 {
        assert_eq!(rotation[j].evaluate(&exact), exact[j + 2]);
    }
}

#[test]
fn view_corners_are_the_four_affine_views_of_the_view_rectangle() {
    let b = ConfigurationBox::from_path("0110100111").unwrap();
    let corners = b.view_corners();
    for (k, u) in corners.iter().enumerate() {
        let corner = b.corner(k as u8);
        assert_eq!(u[0], QSqrt5::from_rational(corner[S].clone()));
        assert_eq!(u[1], QSqrt5::from_rational(corner[T].clone()));
        assert_eq!(u[2], QSqrt5::one());
    }
    let distinct: std::collections::BTreeSet<_> = corners.iter().collect();
    assert_eq!(distinct.len(), 4);
}
