//! `rid-figures <description.json> <output-directory>`: renders one figure
//! description into a new directory holding `figure.svg`, its PNG preview
//! `figure.png` made by the description's converter, and
//! `<source>.partition.json` for every slice subdivided on the way. The
//! directory is assembled under `<output-directory>.incomplete` and renamed
//! only when everything succeeded.
use rid::search::BoxError;
use rid_figures::content::slice::Classifier;
use rid_figures::description::{is_name, Description, DescriptionError};
use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

/// The crate's components, loaded and verified the first time a slice is
/// subdivided.
fn components() -> Result<Box<dyn Classifier>, BoxError> {
    let threads = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    rid_figures::classifier::production(threads)
}

#[derive(Debug)]
enum Error {
    Usage,
    OutputExists(PathBuf),
    Io { path: PathBuf, error: std::io::Error },
    Description(DescriptionError),
    Converter { converter: PathBuf, detail: String },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Usage => write!(f, "usage: rid-figures <description.json> <output-directory>"),
            Error::OutputExists(path) => write!(f, "{} already exists", path.display()),
            Error::Io { path, error } => write!(f, "{}: {error}", path.display()),
            Error::Description(e) => write!(f, "{e}"),
            Error::Converter { converter, detail } => {
                write!(f, "converter {}: {detail}", converter.display())
            }
        }
    }
}

fn io(path: &Path) -> impl FnOnce(std::io::Error) -> Error + '_ {
    move |error| Error::Io {
        path: path.to_owned(),
        error,
    }
}

fn run(arguments: &[OsString]) -> Result<PathBuf, Error> {
    let [description, output] = arguments else {
        return Err(Error::Usage);
    };
    // Rebuilding the path from its components drops a trailing separator,
    // so the staging directory is a sibling of the output directory.
    let description = Path::new(description);
    let output: PathBuf = Path::new(output).components().collect();
    let mut staging = output.clone().into_os_string();
    staging.push(".incomplete");
    let staging = PathBuf::from(staging);
    for path in [&output, &staging] {
        if path.exists() {
            return Err(Error::OutputExists(path.clone()));
        }
    }
    let text = fs::read(description).map_err(io(description))?;
    let parsed = Description::parse(&text).map_err(Error::Description)?;
    // Everything that needs no computation is refused first: the
    // description's structure, then a missing or broken converter.
    parsed.check().map_err(Error::Description)?;
    probe_converter(&parsed.png.converter)?;
    let base = description.parent().unwrap_or(Path::new(""));
    let rendered = parsed
        .render(base, rid::POLICY, &components)
        .map_err(Error::Description)?;
    fs::create_dir(&staging).map_err(io(&staging))?;
    let written = (|| {
        let svg = staging.join("figure.svg");
        fs::write(&svg, &rendered.svg).map_err(io(&svg))?;
        for (name, partition) in &rendered.partitions {
            // Names are file-name safe by construction; a second guard.
            assert!(is_name(name), "source name {name:?} is not a file name");
            let file = staging.join(format!("{name}.partition.json"));
            fs::write(&file, partition.to_json(rid::POLICY)).map_err(io(&file))?;
        }
        convert(&parsed.png.converter, parsed.png.dpi.get(), &svg, &staging.join("figure.png"))?;
        fs::rename(&staging, &output).map_err(io(&output))
    })();
    if written.is_err() {
        // Only this run's own files are in the staging directory.
        let _ = fs::remove_dir_all(&staging);
    }
    written.map(|()| output)
}

/// Runs the converter with `--version`, which must succeed.
fn probe_converter(converter: &Path) -> Result<(), Error> {
    let failure = |detail: String| Error::Converter { converter: converter.to_owned(), detail };
    let output = Command::new(converter).arg("--version").output().map_err(|e| failure(e.to_string()))?;
    if !output.status.success() {
        return Err(failure(format!("--version failed: {}", output.status)));
    }
    Ok(())
}

fn convert(converter: &Path, dpi: u32, svg: &Path, png: &Path) -> Result<(), Error> {
    let failure = |detail: String| Error::Converter {
        converter: converter.to_owned(),
        detail,
    };
    let result = Command::new(converter)
        .arg("--format=png")
        .arg(format!("--dpi-x={dpi}"))
        .arg(format!("--dpi-y={dpi}"))
        .arg("--output")
        .arg(png)
        .arg(svg)
        .output()
        .map_err(|e| failure(e.to_string()))?;
    if !result.status.success() {
        return Err(failure(format!(
            "{}: {}",
            result.status,
            String::from_utf8_lossy(&result.stderr).trim()
        )));
    }
    if !png.is_file() {
        return Err(failure("no PNG was written".into()));
    }
    Ok(())
}

fn main() -> ExitCode {
    let arguments: Vec<OsString> = std::env::args_os().skip(1).collect();
    match run(&arguments) {
        Ok(output) => {
            println!("wrote {}", output.display());
            ExitCode::SUCCESS
        }
        Err(e @ Error::Usage) => {
            eprintln!("{e}");
            ExitCode::from(2)
        }
        Err(e) => {
            eprintln!("rid-figures: {e}");
            ExitCode::from(1)
        }
    }
}
