//! `solids <output-directory> [figure ...]`: draws the article's
//! illustrations of other solids and of the rhombicosidodecahedron's
//! parameters, writing `<figure>.svg` and, with `rsvg-convert` on the path,
//! `<figure>.pdf` and a preview `<figure>.png`. Floating point throughout:
//! these figures illustrate and prove nothing.
mod figures;
mod geometry;
mod optimise;
mod passage;
mod svg;

use std::path::Path;
use std::process::{Command, ExitCode};

const FIGURES: [(&str, fn() -> String); 8] = [
    ("rid", figures::rid_emblem),
    ("nieuwland-cube", figures::nieuwland_cube),
    ("platonic", figures::platonic_figure),
    ("nopert214", figures::nopert214),
    ("triakis", figures::triakis_figure),
    ("rid-parametrisation", figures::rid_parametrisation),
    ("margin", figures::margin),
    ("reduced-ranges", figures::reduced_ranges),
];

fn convert(svg: &Path, format: &str, extra: &[&str]) -> bool {
    Command::new("rsvg-convert")
        .arg("-f")
        .arg(format)
        .args(extra)
        .arg("-o")
        .arg(svg.with_extension(format))
        .arg(svg)
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((output, chosen)) = args.split_first() else {
        eprintln!("usage: solids <output-directory> [figure ...]");
        return ExitCode::from(2);
    };
    for (name, draw) in FIGURES {
        if !chosen.is_empty() && !chosen.iter().any(|c| c == name) {
            continue;
        }
        let path = Path::new(output).join(format!("{name}.svg"));
        if let Err(e) = std::fs::write(&path, draw()) {
            eprintln!("{}: {e}", path.display());
            return ExitCode::FAILURE;
        }
        let converted = convert(&path, "pdf", &[]) && convert(&path, "png", &["--dpi-x=150", "--dpi-y=150"]);
        println!("wrote {}{}", path.display(), if converted { " (.pdf, .png)" } else { "" });
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::figures::*;
    use super::geometry::*;
    use rid::arithmetic::exact::{frac, QSqrt5};
    use std::cmp::Ordering;

    #[test]
    fn the_drawn_passages_fit() {
        let clearance = |p: &super::passage::Passage| p.clearances()[0].1;
        // Renshaw's passage for #214 clears by about 1.8·10⁻⁸.
        let c = clearance(&nopert214_passage());
        assert!(c > 1.7e-8 && c < 1.9e-8, "{c}");
        assert!(clearance(&triakis_passage()) > 3e-6);
        for (name, passage, scale) in platonic_passages() {
            assert!(scale < 0.99 && clearance(&passage) > 0.0, "{name}");
        }
    }

    #[test]
    fn the_drawn_configuration_lies_in_the_reduced_domain() {
        let point = [frac(1, 5), frac(1, 10), frac(1, 4), frac(-1, 12), frac(5, 24)].map(QSqrt5::from_rational);
        assert_eq!(VIEW, [0.2, 0.1]);
        assert_eq!(ROTATION, [0.25, -1.0 / 12.0, 5.0 / 24.0]);
        for p in rid::problem::domain::affine_polynomials() {
            assert_ne!(p.evaluate(&point).sign(), Ordering::Greater);
        }
    }

    #[test]
    fn hulls_and_faces() {
        assert_eq!(rid().faces.len(), 62);
        assert_eq!(cube(0.5).edges().len(), 12);
        let square = hull(&[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5], [1.0, 1.0], [0.0, 1.0]]);
        assert_eq!(square.len(), 4);
    }
}

