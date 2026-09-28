//! Interruption tests with real processes: the test binary starts itself as
//! a child that installs the real signal handlers and runs one command with
//! the stand-in components.

mod support;

use rid::search::certificate::{self, Header};
use rid::search::command::{entry, install_signal_handlers, FAILURE};
use rid::POLICY;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command as Process, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};
use support::*;

const CHILD: &str = "RID_INTERRUPTION_TEST_CHILD";

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ChildSpec {
    args: Vec<PathBuf>,
    stand_in: StandIn,
    /// Created by the child once its stop flag is set, so that the parent
    /// knows that a signal was handled.
    stop_marker: Option<PathBuf>,
}

/// Runs only as a child process of the tests below: installs the real signal
/// handlers and runs one invocation with the stand-in components.
#[test]
fn child() {
    let Some(spec) = std::env::var_os(CHILD) else {
        return;
    };
    let spec: ChildSpec = serde_json::from_str(spec.to_str().unwrap()).unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    install_signal_handlers(&stop).unwrap();
    if let Some(marker) = spec.stop_marker {
        let stop = Arc::clone(&stop);
        thread::spawn(move || {
            while !stop.load(Ordering::Acquire) {
                thread::sleep(Duration::from_millis(1));
            }
            fs::write(marker, b"").unwrap();
        });
    }
    let args: Vec<OsString> = spec.args.into_iter().map(OsString::from).collect();
    let status = entry(&args, &stop, |_, _| Ok(spec.stand_in), &mut io::stdout().lock());
    std::process::exit(i32::from(status));
}

fn spawn(args: &[OsString], stand_in: &StandIn) -> std::process::Child {
    spawn_with(args, stand_in, None)
}

fn spawn_with(
    args: &[OsString],
    stand_in: &StandIn,
    stop_marker: Option<&Path>,
) -> std::process::Child {
    let spec = ChildSpec {
        args: args.iter().map(PathBuf::from).collect(),
        stand_in: stand_in.clone_config(),
        stop_marker: stop_marker.map(Path::to_path_buf),
    };
    Process::new(std::env::current_exe().unwrap())
        .args(["--exact", "child", "--nocapture", "--quiet", "--test-threads=1"])
        .env(CHILD, serde_json::to_string(&spec).unwrap())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap()
}

/// The exit status, the summary line (if any) and standard error.
fn finish(child: std::process::Child) -> (std::process::ExitStatus, Option<Value>, String) {
    let output = child.wait_with_output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    let summary = stdout
        .lines()
        .filter(|line| line.starts_with('{'))
        .last()
        .map(|line| serde_json::from_str(line).unwrap());
    (output.status, summary, String::from_utf8_lossy(&output.stderr).into_owned())
}

#[test]
fn sigint_and_sigterm_stop_a_search_cooperatively() {
    let _processes = processes();
    let (expected, full) = reference_with_summary(11, 64);
    let total = full["decisions"].as_u64().unwrap();
    assert!(total > 200);
    for signal in [signal_hook::consts::SIGINT, signal_hook::consts::SIGTERM] {
        let dir = Folder::new();
        let mut stand_in = StandIn::new(11);
        stand_in.delay_us = 200;
        stand_in.raise_in_check = Some((60, signal));
        let (status, summary, stderr) = finish(spawn(&dir.search(3, 6, 64, None), &stand_in));
        assert_eq!(status.code(), Some(0), "{stderr}");
        let summary = summary.unwrap();
        assert_eq!(summary["halt"], "stopped");
        // The stop came at the 60th evaluation, long before the end.
        assert!(summary["evaluated"].as_u64().unwrap() >= 60);
        assert!(summary["decisions"].as_u64().unwrap() < 60 + 6);
        assert!(summary["unresolved"].as_u64().unwrap() > 0);
        // Everything committed was flushed, and it is a prefix of the full search.
        let bytes = fs::read(dir.certificate()).unwrap();
        assert!(expected.starts_with(&bytes) && bytes.ends_with(b"\n"));
        assert_eq!(succeed(&dir.check(64), 11)["records"], summary["records"]);
        let resumed = succeed(&dir.search(2, 5, 64, None), 11);
        assert_eq!(resumed["halt"], "complete");
        assert_eq!(fs::read(dir.certificate()).unwrap(), expected);
    }
}

/// Waits until `file` exists.
fn wait_for(file: &Path, failure: &str) {
    let started = Instant::now();
    while !file.exists() {
        assert!(started.elapsed() < Duration::from_secs(30), "{failure}");
        thread::sleep(Duration::from_millis(2));
    }
}

/// A second signal ends the process at once (status 1) even when the first
/// one cannot be honoured, here because the running evaluation takes a long
/// time. Each signal is sent only once the child is in the state it tests,
/// as the child reports through marker files: the first while an evaluation
/// is running (sent earlier, it would stop the search cleanly with status 0),
/// the second after the first was handled (sent earlier, the two could merge
/// into one pending signal).
#[test]
fn a_second_signal_ends_a_stuck_search() {
    let _processes = processes();
    let dir = Folder::new();
    let (deciding, stopping) = (dir.0.join("deciding"), dir.0.join("stopping"));
    let mut stand_in = StandIn::new(13);
    // At least 30 s: far longer than the steps below take.
    stand_in.delay_us = 30_000_000;
    stand_in.check_marker = Some(deciding.clone());
    let child = spawn_with(&dir.search(1, 1, 64, None), &stand_in, Some(&stopping));
    let pid = child.id().to_string();
    let terminate = || {
        assert!(Process::new("kill").args(["-TERM", &pid]).status().unwrap().success());
    };
    wait_for(&deciding, "no evaluation started");
    terminate();
    wait_for(&stopping, "the first signal was not handled");
    let second = Instant::now();
    terminate();
    let (status, summary, stderr) = finish(child);
    assert!(second.elapsed() < Duration::from_secs(10), "the second signal did not end it");
    assert_eq!(status.code(), Some(i32::from(FAILURE)), "{stderr}");
    assert!(summary.is_none());
    // The file left behind is a valid (here empty) certificate.
    assert_eq!(succeed(&dir.check(64), 13)["records"], 0);
}

#[test]
fn a_signal_during_verification_leaves_the_certificate_unchanged() {
    let _processes = processes();
    let dir = Folder::new();
    succeed(&dir.search(3, 4, 64, Some(120)), 12);
    let mut bytes = fs::read(dir.certificate()).unwrap();
    bytes.extend_from_slice(b"0123 {\"pa");
    fs::write(dir.certificate(), &bytes).unwrap();
    let before = snapshot(&dir.certificate());
    for args in [dir.search(3, 4, 64, None), dir.check(64)] {
        let mut stand_in = StandIn::new(12);
        stand_in.raise_in_holds = Some((5, signal_hook::consts::SIGTERM));
        // The signal sets the stop flag within microseconds; slow
        // verifications make sure that the other verification threads have
        // records left to check after that (with instant ones they could
        // finish all of them first, and the stop would come too late).
        stand_in.holds_delay_us = 2_000;
        let (status, summary, stderr) = finish(spawn(&args, &stand_in));
        assert_eq!(status.code(), Some(i32::from(FAILURE)), "{stderr}");
        assert!(summary.is_none());
        assert!(stderr.contains("stopped by a signal during certificate verification"), "{stderr}");
        assert_eq!(snapshot(&dir.certificate()), before);
    }
}

#[test]
fn a_long_layer_reaches_the_disk_before_its_barrier() {
    let _processes = processes();
    let dir = Folder::new();
    // Depths 0 to 5 refuse 63 boxes; depth 6 holds 64 records and takes
    // about 1.3 s on one thread, longer than the flush interval.
    let mut stand_in = StandIn::new(0);
    stand_in.accept_from_depth = Some(6);
    stand_in.delay_us = 10_000;
    let mut child = spawn(&dir.search(1, 1, 64, None), &stand_in);
    let header = Header::new(POLICY, "", 64).unwrap();
    let started = Instant::now();
    let partial = loop {
        assert!(started.elapsed() < Duration::from_secs(30), "no partial layer seen");
        let records = fs::read(dir.certificate())
            .map(|bytes| bytes.iter().filter(|&&b| b == b'\n').count().saturating_sub(1))
            .unwrap_or(0);
        assert!(records < 64, "the layer reached the disk only at its barrier");
        if records > 0 {
            break records;
        }
        thread::sleep(Duration::from_millis(5));
    };
    child.kill().unwrap();
    let (status, _, _) = finish(child);
    assert_eq!(status.signal(), Some(9));
    let checked = certificate::check(&dir.certificate(), &header, NonZeroUsize::MIN, |p, d| {
        rid::search::command::Components::holds(&stand_in, p, d)
    })
    .unwrap();
    assert!(checked.records as usize >= partial && checked.records < 64);
}

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// Kills searches with SIGKILL at random times and restarts them until one
/// completes, in episodes that each start from no file, until `kills` kills
/// have happened. After every kill the file is a prefix of the uninterrupted
/// certificate and checks as valid; each episode ends with exactly that
/// certificate. Returns the number of kills.
fn kill_campaign(seed: u64, kills: usize) -> usize {
    let expected = reference(seed, 64);
    let mut state = 0x9e37_79b9_7f4a_7c15 ^ seed;
    let mut killed = 0;
    while killed < kills {
        let dir = Folder::new();
        for attempt in 0.. {
            assert!(attempt < 200, "the search did not complete");
            let threads = 1 + attempt % 3;
            let mut stand_in = StandIn::new(seed);
            stand_in.delay_us = 400;
            let mut child = spawn(&dir.search(threads, 2 * threads, 64, None), &stand_in);
            // Some kills land in start-up: creation, verification, recovery.
            let delay = match random(&mut state) % 4 {
                0 => random(&mut state) % 20_000,
                _ => random(&mut state) % 250_000,
            };
            thread::sleep(Duration::from_micros(delay));
            child.kill().unwrap();
            let (status, summary, stderr) = finish(child);
            if status.signal() == Some(9) {
                killed += 1;
                let bytes = fs::read(dir.certificate()).unwrap_or_default();
                assert!(expected.starts_with(&bytes), "attempt {attempt}: not a prefix");
                eprintln!(
                    "seed {seed} attempt {attempt}: killed after {delay} us with {} of {} bytes",
                    bytes.len(),
                    expected.len()
                );
                if dir.certificate().exists() {
                    let checked = succeed(&dir.check(64), seed);
                    assert_eq!(checked["complete"], json!(bytes == expected));
                }
                // Now and then imitate a write torn inside a line.
                let rest = &expected[bytes.len()..];
                if let Some(end) = rest.iter().position(|&b| b == b'\n') {
                    if !bytes.is_empty() && random(&mut state) % 3 == 0 {
                        let torn = &rest[..random(&mut state) as usize % (end + 1)];
                        fs::OpenOptions::new()
                            .append(true)
                            .open(dir.certificate())
                            .and_then(|mut file| file.write_all(torn))
                            .unwrap();
                    }
                }
                continue;
            }
            // The child finished before the kill.
            assert_eq!(status.code(), Some(0), "{stderr}");
            if summary.unwrap()["halt"] == "complete" {
                break;
            }
        }
        assert_eq!(fs::read(dir.certificate()).unwrap(), expected);
        assert_eq!(succeed(&dir.check(64), seed)["complete"], true);
    }
    killed
}

#[test]
fn sigkill_at_random_times_then_restart_reproduces_the_certificate() {
    let _processes = processes();
    let killed = kill_campaign(21, 12);
    assert!(killed >= 12);
}

/// The heavier campaign: 200 kills over several trees. Run with `cargo test
/// --release --test interruption sigkill_campaign_heavy -- --ignored`.
#[test]
#[ignore]
fn sigkill_campaign_heavy() {
    let _processes = processes();
    for seed in 100..110 {
        assert!(kill_campaign(seed, 20) >= 20);
    }
}
