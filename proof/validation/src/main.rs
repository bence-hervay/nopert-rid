//! `rid-validation <completeness|irredundancy|pilot|targets|check|compare> <config.json>`:
//! runs one command and prints its summary as one JSON object.
use rid::components::collection::Collection;
use rid::components::record::{ComponentName, RecordData};
use rid::components::Component;
use rid::problem::configuration::ConfigurationBox;
use rid::search::command::wiring;
use rid::search::BoxError;
use rid_validation::config::{Command, Error as ConfigError};
use rid_validation::experiments::probe::{canonical, Probes};
use rid_validation::experiments::{self, completeness, irredundancy, pilot};
use serde::Serialize;
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::ExitCode;
use std::sync::atomic::AtomicBool;

/// Exit status: the summary was printed.
const SUCCESS: u8 = 0;
/// Exit status: a refusal or failure; the reason is on standard error.
const FAILURE: u8 = 1;
/// Exit status: a malformed command line or an unknown command.
const USAGE: u8 = 2;
/// Exit status: `compare` found a regression and was asked to fail on one.
const REGRESSION: u8 = 3;

/// The crate's collection as the experiments see it: its check, and each
/// component's check, `holds` and cover records in the collection's order.
struct Components(Collection);

impl Probes for Components {
    type Data = RecordData;
    fn components(&self) -> Vec<String> {
        ComponentName::ORDER.iter().map(ComponentName::to_string).collect()
    }
    fn decide(&self, b: &ConfigurationBox) -> Result<Option<RecordData>, BoxError> {
        Ok(self.0.check(b))
    }
    fn attempt(&self, index: usize, b: &ConfigurationBox) -> Option<RecordData> {
        let c = &self.0;
        match ComponentName::ORDER[index] {
            ComponentName::Domain => c.domain.check(b).map(|inequality| RecordData::Domain { inequality }),
            ComponentName::Exotic => c.exotic.check(b).map(|cover| RecordData::Exotic { cover }),
            ComponentName::Local => c.local.check(b).map(|cover| RecordData::Local { cover }),
            ComponentName::Global => c.global.check(b).map(RecordData::Global),
        }
    }
    fn covers(&self, index: usize) -> Vec<RecordData> {
        match ComponentName::ORDER[index] {
            ComponentName::Exotic => self.0.exotic.keys().into_iter().map(|cover| RecordData::Exotic { cover }).collect(),
            ComponentName::Local => self.0.local.keys().into_iter().map(|cover| RecordData::Local { cover }).collect(),
            ComponentName::Domain | ComponentName::Global => Vec::new(),
        }
    }
    /// The record holds for the component it names, which must be `index`.
    fn verify(&self, index: usize, b: &ConfigurationBox, data: &RecordData) -> Result<(), BoxError> {
        if data.component() != ComponentName::ORDER[index] {
            return Err(format!("a {} record given to {}", data.component(), ComponentName::ORDER[index]).into());
        }
        Ok(self.0.holds(b, data)?)
    }
}

/// Loads and verifies the proved covers of the crate's catalogue on
/// `threads` threads and builds the collection.
fn components(threads: NonZeroUsize) -> Result<Components, BoxError> {
    Ok(Components(wiring::prepare(threads, &AtomicBool::new(false))?))
}

fn print(value: &impl Serialize) -> Result<(), BoxError> {
    println!("{}", serde_json::to_string_pretty(value)?);
    Ok(())
}

fn execute(command: Command) -> Result<u8, BoxError> {
    match command {
        Command::Completeness(config) => print(&completeness::run(&components(config.threads)?, &config)?)?,
        Command::Irredundancy(config) => print(&irredundancy::run(&components(config.threads)?, &config)?)?,
        Command::Pilot(config) => print(&pilot::run(&components(config.threads)?, &config)?)?,
        Command::Targets(config) => print(&irredundancy::targets::write(&config)?)?,
        Command::Check(config) => print(&experiments::check(&config, canonical::<RecordData>)?)?,
        Command::Compare(config) => {
            let comparison = experiments::compare(&config, canonical::<RecordData>)?;
            print(&comparison)?;
            if config.fail_on_regression && comparison.regressions > 0 {
                return Ok(REGRESSION);
            }
        }
    }
    Ok(SUCCESS)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [name, file] = args.as_slice() else {
        eprintln!(
            "usage: rid-validation <completeness|irredundancy|pilot|targets|check|compare> <config.json>"
        );
        return ExitCode::from(USAGE);
    };
    let command = match Command::load(name, Path::new(file)) {
        Ok(command) => command,
        Err(error @ ConfigError::UnknownCommand(_)) => {
            eprintln!("rid-validation: {error}");
            return ExitCode::from(USAGE);
        }
        Err(error) => {
            eprintln!("rid-validation: {error}");
            return ExitCode::from(FAILURE);
        }
    };
    match execute(command) {
        Ok(status) => ExitCode::from(status),
        Err(error) => {
            eprintln!("rid-validation: {error}");
            ExitCode::from(FAILURE)
        }
    }
}
