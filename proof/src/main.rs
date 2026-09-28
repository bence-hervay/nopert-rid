//! `rid <prepare|search|check|generate> <config.json>`: installs the signal
//! handlers and runs one command. `prepare`, `search` and `check` run with the
//! production components; `generate`
//! writes cover files.
use rid::elimination::generate::command as generate;
use rid::search::command::{self, wiring, FAILURE};
use std::ffi::OsString;
use std::io;
use std::process::ExitCode;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;

fn main() -> ExitCode {
    let stop = Arc::new(AtomicBool::new(false));
    if let Err(e) = command::install_signal_handlers(&stop) {
        eprintln!("rid: cannot install signal handlers: {e}");
        return ExitCode::from(FAILURE);
    }
    let args: Vec<OsString> = std::env::args_os().skip(1).collect();
    if args.first().is_some_and(|name| name == "generate") {
        return ExitCode::from(generate::entry(&args[1..], &stop, &mut io::stdout().lock()));
    }
    ExitCode::from(command::entry(&args, &stop, wiring::prepare, &mut io::stdout().lock()))
}
