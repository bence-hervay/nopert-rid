use super::*;
use crate::drawing::style::Stroke;
use crate::random::Random;
use rid::arithmetic::exact::{frac, Interval};

#[test]
fn decimals_round_halves_away_from_zero() {
    for (n, d, text) in [
        (0, 1, "0"),
        (1, 2, "0.5"),
        (-1, 2, "-0.5"),
        (1, 3, "0.3333"),
        (2, 3, "0.6667"),
        (-2, 3, "-0.6667"),
        (1, 20000, "0.0001"),
        (-1, 20000, "-0.0001"),
        (1, 20001, "0"),
        (-1, 20001, "0"),
        (123456789, 1000, "123456.789"),
        (7, 1, "7"),
        (-7, 1, "-7"),
        (199999, 200000, "1"),
    ] {
        assert_eq!(decimal(&frac(n, d)), text, "{n}/{d}");
    }
}

#[test]
fn decimals_are_the_nearest_multiples_of_the_unit() {
    let mut random = Random::new(31);
    let unit = frac(1, 10_000);
    let half = frac(1, 20_000);
    for _ in 0..2000 {
        let x = random.rational(&frac(-1000, 1), &frac(1000, 1), 40);
        let text = decimal(&x);
        // Parse the decimal back exactly.
        let (sign, digits) = text.strip_prefix('-').map_or((1, text.as_str()), |t| (-1, t));
        let (int, frac_part) = digits.split_once('.').unwrap_or((digits, ""));
        assert!(frac_part.len() <= DECIMALS && !frac_part.ends_with('0'));
        let scaled: BigInt = format!("{int}{frac_part:0<4}").parse().unwrap();
        let value = Q::from_integer(scaled * sign) * &unit;
        let error = &value - &x;
        assert!(error <= half && -&error <= half, "{x} as {text}");
    }
}

fn style(fill: Option<&str>, stroke: Option<(&str, Q)>, opacity: Q) -> Style {
    Style::new(
        fill.map(|c| Color::parse(c).unwrap()),
        stroke.map(|(c, w)| Stroke::new(Color::parse(c).unwrap(), w).unwrap()),
        opacity,
    )
    .unwrap()
}

#[test]
fn documents_have_millimetre_pages_and_one_group_per_layer() {
    let canvas = Canvas::new(frac(180, 1), frac(201, 2)).unwrap();
    let mut document = Document::new(&canvas, Some(Color::parse("#FFFFFF").unwrap()));
    let p = |x: i64, y: i64| [frac(x, 1), frac(y, 3)];
    let clip = [Interval::new(frac(10, 1), frac(90, 1)).unwrap(), Interval::new(frac(5, 1), frac(95, 1)).unwrap()];
    document.layer(
        &[Shape::Segment([p(1, 1), p(2, 2)]), Shape::Polygon(vec![p(0, 0), p(3, 0), p(0, 3)])],
        &style(None, Some(("#ff8000", frac(7, 50))), Q::one()),
        Some(&clip),
    );
    document.layer(&[Shape::Polygon(vec![p(5, 5)])], &style(Some("#0080FF"), None, frac(1, 6)), None);
    let svg = document.finish();
    let expected = "\
<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"180mm\" height=\"100.5mm\" viewBox=\"0 0 180 100.5\">
<rect x=\"0\" y=\"0\" width=\"180\" height=\"100.5\" fill=\"#ffffff\"/>
<clipPath id=\"clip-1\"><rect x=\"10\" y=\"5\" width=\"80\" height=\"90\"/></clipPath>
<g clip-path=\"url(#clip-1)\" opacity=\"1\" fill=\"none\" fill-rule=\"nonzero\" stroke=\"#ff8000\" stroke-width=\"0.14\" stroke-linecap=\"round\" stroke-linejoin=\"round\">
<path d=\"M1 0.3333L2 0.6667\"/>
<path d=\"M0 0L3 0L0 1Z\"/>
</g>
<g opacity=\"0.1667\" fill=\"#0080ff\" fill-rule=\"nonzero\" stroke=\"none\">
<path d=\"M5 1.6667Z\"/>
</g>
</svg>
";
    assert_eq!(svg, expected);
}

#[test]
fn transparent_pages_have_no_background() {
    let canvas = Canvas::new(frac(10, 1), frac(10, 1)).unwrap();
    let svg = Document::new(&canvas, None).finish();
    assert!(!svg.contains("<rect"));
    assert_eq!(svg.matches("<svg").count(), 1);
}
