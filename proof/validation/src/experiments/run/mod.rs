//! Parallel evaluation of independent items with results committed in item
//! order, and the JSON-lines transcripts the experiments write and read.
//!
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc;
use std::thread;

/// Sets the flag if its thread unwinds, so that the other workers stop.
struct PanicGuard<'a>(&'a AtomicBool);

impl Drop for PanicGuard<'_> {
    fn drop(&mut self) {
        if thread::panicking() {
            self.0.store(true, Ordering::Release);
        }
    }
}

/// Evaluates every item on at most `threads` workers and passes each result
/// to `commit` on the calling thread, strictly in item order.
///
/// Workers take the next unevaluated item, so a slow item does not hold the
/// others back; results that finish early wait in memory for their turn. The
/// first error (of `evaluate` or `commit`, in item order) stops the workers
/// after their current item and is returned; nothing after it is committed.
/// A panic in `evaluate` stops the workers likewise and is propagated.
pub fn ordered<T, R, E>(
    items: &[T],
    threads: NonZeroUsize,
    evaluate: impl Fn(usize, &T) -> Result<R, E> + Sync,
    mut commit: impl FnMut(usize, R) -> Result<(), E>,
) -> Result<(), E>
where
    T: Sync,
    R: Send,
    E: Send,
{
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let (sender, receiver) = mpsc::channel::<(usize, Result<R, E>)>();
    thread::scope(|scope| {
        for _ in 0..threads.get().min(items.len()) {
            let (sender, next, failed, evaluate) = (sender.clone(), &next, &failed, &evaluate);
            scope.spawn(move || {
                let _guard = PanicGuard(failed);
                while !failed.load(Ordering::Acquire) {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    if index >= items.len() {
                        return;
                    }
                    let result = evaluate(index, &items[index]);
                    if result.is_err() {
                        failed.store(true, Ordering::Release);
                    }
                    if sender.send((index, result)).is_err() {
                        return;
                    }
                }
            });
        }
        drop(sender);
        let mut waiting = BTreeMap::new();
        let mut expected = 0;
        for (index, result) in receiver {
            waiting.insert(index, result);
            while let Some(result) = waiting.remove(&expected) {
                if let Err(error) = result.and_then(|value| commit(expected, value)) {
                    failed.store(true, Ordering::Release);
                    return Err(error);
                }
                expected += 1;
            }
        }
        // Fewer commits than items only after a worker panicked; the scope
        // propagates that panic when it ends.
        Ok(())
    })
}

/// A transcript being written: JSON lines, each synchronised to disk when
/// written, so an interruption leaves a prefix of complete lines.
pub struct Writer {
    file: File,
    path: PathBuf,
}

impl Writer {
    /// Creates the file; an existing file is refused, never overwritten.
    pub fn create(path: &Path) -> Result<Self, TranscriptError> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|source| TranscriptError::Io {
                path: path.to_owned(),
                source,
            })?;
        Ok(Self {
            file,
            path: path.to_owned(),
        })
    }

    /// Appends `value` as one line and synchronises the file.
    pub fn write<T: Serialize>(&mut self, value: &T) -> Result<(), TranscriptError> {
        let mut bytes = serde_json::to_vec(value).expect("transcript values serialise");
        bytes.push(b'\n');
        self.file
            .write_all(&bytes)
            .and_then(|()| self.file.sync_data())
            .map_err(|source| TranscriptError::Io {
                path: self.path.clone(),
                source,
            })
    }
}

#[derive(Debug)]
pub enum TranscriptError {
    Io { path: PathBuf, source: io::Error },
    /// The file does not end with a complete line, or has an empty line.
    Incomplete,
    /// Line `line` (counted from 1) does not parse into its type.
    Json { line: usize, source: serde_json::Error },
    /// Line `line` parses but is not in the canonical spelling.
    NotCanonical { line: usize },
    /// The lines do not make a header, entries and a summary.
    Shape(String),
}

impl fmt::Display for TranscriptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TranscriptError::Io { path, source } => write!(f, "{}: {source}", path.display()),
            TranscriptError::Incomplete => {
                write!(f, "the transcript is empty or ends inside a line")
            }
            TranscriptError::Json { line, source } => write!(f, "line {line}: {source}"),
            TranscriptError::NotCanonical { line } => {
                write!(f, "line {line} is not in canonical form")
            }
            TranscriptError::Shape(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for TranscriptError {}

/// Reads a whole file.
pub fn read(path: &Path) -> Result<Vec<u8>, TranscriptError> {
    std::fs::read(path).map_err(|source| TranscriptError::Io {
        path: path.to_owned(),
        source,
    })
}

/// The complete lines of a transcript, without their newlines. A transcript
/// ends with a newline and has no empty line.
pub fn lines(bytes: &[u8]) -> Result<Vec<&[u8]>, TranscriptError> {
    let body = bytes
        .strip_suffix(b"\n")
        .ok_or(TranscriptError::Incomplete)?;
    let lines: Vec<&[u8]> = body.split(|&b| b == b'\n').collect();
    if lines.iter().any(|line| line.is_empty()) {
        return Err(TranscriptError::Incomplete);
    }
    Ok(lines)
}

/// Parses line `number` (counted from 1) into `T` and refuses any spelling
/// other than the one `T` serialises to.
pub fn parse<T>(line: &[u8], number: usize) -> Result<T, TranscriptError>
where
    T: Serialize + DeserializeOwned,
{
    let value: T = serde_json::from_slice(line).map_err(|source| TranscriptError::Json {
        line: number,
        source,
    })?;
    if serde_json::to_vec(&value).expect("parsed values serialise") != line {
        return Err(TranscriptError::NotCanonical { line: number });
    }
    Ok(value)
}

/// The fields every transcript header has, checked alike by every
/// experiment.
pub struct Binding<'a> {
    pub format: &'a str,
    pub policy: &'a str,
    pub executable_sha256: &'a str,
    pub catalogue_sha256: &'a str,
    pub components: &'a [String],
}

impl Binding<'_> {
    /// The format tag is `format`, the header is bound to the catalogue
    /// with digest `catalogue_sha256`, the policy and executable digests are
    /// well formed, and the components are the decided four in order.
    /// Returns the header problem as a message.
    pub fn check(&self, format: &str, catalogue_sha256: &str) -> Result<(), String> {
        if self.format != format {
            return Err("wrong format".into());
        }
        if self.catalogue_sha256 != catalogue_sha256 {
            return Err("bound to another catalogue".into());
        }
        if !is_digest(self.policy) || !is_digest(self.executable_sha256) {
            return Err("the policy and the executable digest must be 64 hexadecimal digits".into());
        }
        if !crate::experiments::probe::decided(self.components) {
            return Err("the components must be Domain, Exotic, Local, Global in this order".into());
        }
        Ok(())
    }
}

/// Checks the frame of a transcript: a header line, one line per selected
/// catalogue entry in catalogue order, and a summary line. `header` checks
/// the parsed header and returns the selected entries and the empty
/// summary; `entry` checks one parsed entry line against its catalogue entry
/// (on `threads` workers, results in order); `add` accumulates the summary,
/// which must equal the last line (`mismatch` otherwise). `noun` names the
/// entries in messages.
pub fn check_frame<'c, H, C, R, S, X>(
    bytes: &[u8],
    threads: NonZeroUsize,
    noun: &str,
    header: impl FnOnce(&H) -> Result<(Vec<&'c C>, S), X>,
    entry: impl Fn(&H, &C, &R) -> Result<(), X> + Sync,
    mut add: impl FnMut(&mut S, &H, &R),
    mismatch: X,
) -> Result<(H, Vec<R>, S), X>
where
    H: Serialize + DeserializeOwned + Sync,
    C: Sync + 'c,
    R: Serialize + DeserializeOwned + Send,
    S: Serialize + DeserializeOwned + PartialEq,
    X: From<TranscriptError> + Send,
{
    let lines = lines(bytes)?;
    if lines.len() < 2 {
        return Err(TranscriptError::Shape("a transcript needs a header and a summary".into()).into());
    }
    let head: H = parse(lines[0], 1)?;
    let (selected, mut summary) = header(&head)?;
    if lines.len() != selected.len() + 2 {
        let message = format!("{} lines for {} {noun}", lines.len(), selected.len());
        return Err(TranscriptError::Shape(message).into());
    }
    let entries: Vec<(usize, &[u8], &C)> = lines[1..lines.len() - 1]
        .iter()
        .zip(&selected)
        .enumerate()
        .map(|(i, (line, c))| (i + 2, *line, *c))
        .collect();
    let mut results = Vec::with_capacity(selected.len());
    ordered(
        &entries,
        threads,
        |_, &(number, line, c)| {
            let result: R = parse(line, number)?;
            entry(&head, c, &result)?;
            Ok::<_, X>(result)
        },
        |_, result| {
            add(&mut summary, &head, &result);
            results.push(result);
            Ok(())
        },
    )?;
    let found: S = parse(lines[lines.len() - 1], lines.len())?;
    if found != summary {
        return Err(mismatch);
    }
    Ok((head, results, summary))
}

/// Whether `text` is 64 lowercase hexadecimal digits, the spelling of the
/// policy fingerprint and of SHA-256 digests in transcript headers.
pub fn is_digest(text: &str) -> bool {
    text.len() == 64 && text.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

/// SHA-256 of the running executable, recorded in transcript headers.
pub fn executable_sha256() -> Result<String, TranscriptError> {
    let path = std::env::current_exe().map_err(|source| TranscriptError::Io {
        path: PathBuf::from("<current executable>"),
        source,
    })?;
    Ok(crate::points::catalogue::sha256_hex(&read(&path)?))
}

#[cfg(test)]
mod tests;
