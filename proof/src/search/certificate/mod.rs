//! The certificate file: framing, canonical records, ordering, the append-only
//! store and recovery of the unresolved frontier. No geometry: the record data
//! is a caller type and its meaning is established only by the caller's
//! verifier.

use super::panic_message;
pub use super::BoxError;
use crate::problem::configuration::{self, is_path, search_order};
use fs2::FileExt;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::marker::PhantomData;
use std::num::NonZeroUsize;
use std::ops::Bound;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;

/// Largest depth a header may permit: the deepest box a path can name,
/// [`configuration::MAX_DEPTH`]. The root box has depth zero.
pub const MAX_DEPTH: usize = configuration::MAX_DEPTH;
/// Largest line, counting checksum, space, payload and newline.
pub const MAX_LINE: usize = 1 << 20;
/// Checksum characters at the start of each line.
const CHECKSUM_LEN: usize = 8;
/// Buffer of the append writer; `Store::flush` empties it.
const WRITE_BUFFER: usize = 1 << 16;

/// The typed data a component stores with the path of a box it eliminated.
///
/// Contract, checked by the store on every append and every load:
/// - `serde_json` serialises a value to a JSON object without a top-level
///   field named `path`;
/// - parsing that object and serialising the result reproduces its bytes
///   exactly, so every accepted line is in canonical form. This round trip is
///   what refuses unknown, reordered and repeated fields, even where serde
///   alone would ignore them; `#[serde(deny_unknown_fields)]` is advisable
///   for clearer messages but not relied on.
///
/// A record's meaning comes only from the caller's verifier.
pub trait RecordData: Serialize + DeserializeOwned + Send + Sync {}

/// Identifies this line format in the header.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Format {
    #[serde(rename = "rid-certificate/1")]
    V1,
}

/// The first line of every certificate. Everything a record's validity
/// depends on besides its own line: the format, the policy fingerprint of the
/// implementation and its fixed data, the root box and the depth limit.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    format: Format,
    policy: String,
    root: String,
    max_depth: usize,
}

impl Header {
    /// `policy` is 64 lowercase hexadecimal digits, `max_depth` is at most
    /// [`MAX_DEPTH`] and `root` is a binary path no deeper than `max_depth`.
    pub fn new(policy: &str, root: &str, max_depth: usize) -> Result<Self, Error> {
        let header = Self {
            format: Format::V1,
            policy: policy.to_owned(),
            root: root.to_owned(),
            max_depth,
        };
        header.validate()?;
        Ok(header)
    }

    /// Also run on every expected header, since a deserialised one skipped `new`.
    fn validate(&self) -> Result<(), Error> {
        let hex = |b: u8| matches!(b, b'0'..=b'9' | b'a'..=b'f');
        if self.policy.len() != 64 || !self.policy.bytes().all(hex) {
            return Err(Error::InvalidHeader(
                "policy must be 64 lowercase hexadecimal digits",
            ));
        }
        if self.max_depth > MAX_DEPTH {
            return Err(Error::InvalidHeader("max_depth exceeds the format limit"));
        }
        if !is_path(&self.root) || self.root.len() > self.max_depth {
            return Err(Error::InvalidHeader(
                "root must be a binary path no deeper than max_depth",
            ));
        }
        Ok(())
    }

    pub fn format(&self) -> Format {
        self.format
    }

    pub fn policy(&self) -> &str {
        &self.policy
    }

    pub fn root(&self) -> &str {
        &self.root
    }

    pub fn max_depth(&self) -> usize {
        self.max_depth
    }

    /// Whether `path` names a box of this certificate's tree.
    fn contains(&self, path: &str) -> bool {
        path.len() <= self.max_depth && path.starts_with(&self.root) && is_path(path)
    }

    fn line(&self) -> Vec<u8> {
        frame(&serde_json::to_vec(self).expect("a header always serialises"))
    }
}

/// What was found in a certificate file, before any repair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Inspection {
    /// False for an empty file or an interrupted header write.
    pub header_complete: bool,
    /// Complete records, all verified.
    pub records: u64,
    /// The path of the last complete record.
    pub last: Option<String>,
    /// Length of the complete lines.
    pub valid_bytes: u64,
    /// Length of an unterminated final line (an interrupted write).
    pub unterminated_bytes: u64,
    /// The number of minimal unresolved boxes: the maximal boxes below the
    /// root that neither are recorded nor contain a record. Counted without
    /// listing them; [`Store::frontier`] lists them.
    pub frontier: u64,
}

impl Inspection {
    /// A complete header, no unterminated line and no unresolved box.
    pub fn complete(&self) -> bool {
        self.header_complete && self.unterminated_bytes == 0 && self.frontier == 0
    }
}

/// The current state of an open store.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StoreReport {
    pub records: u64,
    /// File length including appends still buffered.
    pub bytes: u64,
    pub last: Option<String>,
}

/// Why a complete line is not a well-formed canonical line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Defect {
    /// Not `<8 lowercase hex> <payload>\n`.
    Frame,
    /// The checksum does not match the payload.
    Checksum,
    /// The payload does not parse as the expected type.
    Json(String),
    /// The payload parses, but serialising the result gives different bytes.
    NotCanonical,
}

#[derive(Debug)]
pub enum Error {
    /// The caller asked for an invalid header.
    InvalidHeader(&'static str),
    /// An operating-system operation failed.
    Io { action: &'static str, source: io::Error },
    /// Another handle holds a conflicting lock (or locking failed).
    Locked(io::Error),
    /// The file length changed while it was read under the lock.
    Changed,
    /// A line is longer than [`MAX_LINE`].
    TooLong { line: u64 },
    /// A complete line is malformed.
    Malformed { line: u64, defect: Defect },
    /// The complete header is well formed but is not the expected header.
    WrongHeader { found: Box<Header> },
    /// An unterminated first line is not a prefix of the expected header.
    WrongPartialHeader,
    /// A path is not binary, lies outside the root or is deeper than allowed.
    InvalidPath { line: u64, path: String },
    /// A path does not follow its predecessor in (depth, path) order.
    OutOfOrder { line: u64, path: String, previous: String },
    /// A path equals, contains or lies inside an earlier record's path.
    Overlap { line: u64, path: String, other: String },
    /// Record data has no canonical single-object form.
    Unrepresentable { path: String, reason: String },
    /// The caller's verifier refused a record (or panicked on it).
    Refuted { line: u64, path: String, source: BoxError },
    /// The certificate path names something other than a regular file (a
    /// directory or a FIFO, for example), which is refused before reading.
    NotAFile,
    /// An earlier write or synchronisation failed on this handle.
    Failed,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidHeader(reason) => write!(f, "invalid certificate header: {reason}"),
            Error::Io { action, source } => write!(f, "cannot {action}: {source}"),
            Error::Locked(source) => write!(f, "certificate is locked by another handle: {source}"),
            Error::Changed => write!(f, "certificate changed while it was read under a lock"),
            Error::TooLong { line } => write!(f, "line {line} exceeds {MAX_LINE} bytes"),
            Error::Malformed { line, defect } => match defect {
                Defect::Frame => write!(f, "line {line} is not a checksum, a space and a payload"),
                Defect::Checksum => write!(f, "line {line} has a wrong checksum"),
                Defect::Json(e) => write!(f, "line {line} has an invalid payload: {e}"),
                Defect::NotCanonical => write!(f, "line {line} is not in canonical form"),
            },
            Error::WrongHeader { found } => write!(
                f,
                "certificate header (policy {}, root {:?}, max_depth {}) is not the expected one",
                found.policy, found.root, found.max_depth
            ),
            Error::WrongPartialHeader => {
                write!(f, "unterminated first line is not a prefix of the expected header")
            }
            Error::InvalidPath { line, path } => {
                write!(f, "line {line}: path {path:?} is not a box of this certificate")
            }
            Error::OutOfOrder { line, path, previous } => write!(
                f,
                "line {line}: path {path:?} does not follow {previous:?} in (depth, path) order"
            ),
            Error::Overlap { line, path, other } => {
                write!(f, "line {line}: path {path:?} overlaps the recorded path {other:?}")
            }
            Error::Unrepresentable { path, reason } => {
                write!(f, "record data for {path:?} cannot be written: {reason}")
            }
            Error::Refuted { line, path, source } => {
                write!(f, "line {line}: record for {path:?} failed verification: {source}")
            }
            Error::NotAFile => write!(f, "the certificate path is not a regular file"),
            Error::Failed => write!(f, "an earlier write failed on this certificate handle"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io { source, .. } | Error::Locked(source) => Some(source),
            Error::Refuted { source, .. } => Some(source.as_ref()),
            _ => None,
        }
    }
}

fn io_error(action: &'static str) -> impl FnOnce(io::Error) -> Error {
    move |source| Error::Io { action, source }
}

/// The first 32 bits of SHA-256 as 8 lowercase hexadecimal characters.
fn checksum(payload: &[u8]) -> [u8; CHECKSUM_LEN] {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let digest = Sha256::digest(payload);
    let mut out = [0; CHECKSUM_LEN];
    for (i, byte) in digest[..CHECKSUM_LEN / 2].iter().enumerate() {
        out[2 * i] = HEX[usize::from(byte >> 4)];
        out[2 * i + 1] = HEX[usize::from(byte & 15)];
    }
    out
}

fn frame(payload: &[u8]) -> Vec<u8> {
    let mut line = Vec::with_capacity(payload.len() + CHECKSUM_LEN + 2);
    line.extend_from_slice(&checksum(payload));
    line.push(b' ');
    line.extend_from_slice(payload);
    line.push(b'\n');
    line
}

/// The payload of a complete line (one ending in its only newline).
fn unframe(line: &[u8]) -> Result<&[u8], Defect> {
    let body = line.strip_suffix(b"\n").ok_or(Defect::Frame)?;
    if body.len() <= CHECKSUM_LEN
        || body[CHECKSUM_LEN] != b' '
        || !body[..CHECKSUM_LEN]
            .iter()
            .all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
    {
        return Err(Defect::Frame);
    }
    let payload = &body[CHECKSUM_LEN + 1..];
    if body[..CHECKSUM_LEN] != checksum(payload) {
        return Err(Defect::Checksum);
    }
    Ok(payload)
}

/// `{"path":"<path>",<fields of data>}`: the path first, then the data's own
/// fields in their serialisation order. `path` must be binary.
fn encode_record<D: RecordData>(path: &str, data: &D) -> Result<Vec<u8>, String> {
    debug_assert!(is_path(path));
    let body = serde_json::to_vec(data).map_err(|e| format!("serialisation failed: {e}"))?;
    let object: serde_json::Map<String, serde_json::Value> = serde_json::from_slice(&body)
        .map_err(|_| "data does not serialise to a JSON object".to_owned())?;
    if object.contains_key("path") {
        return Err("data has a top-level field named \"path\"".into());
    }
    let mut payload = Vec::with_capacity(body.len() + path.len() + 11);
    payload.extend_from_slice(b"{\"path\":\"");
    payload.extend_from_slice(path.as_bytes());
    payload.push(b'"');
    if object.is_empty() {
        payload.push(b'}');
    } else {
        payload.push(b',');
        payload.extend_from_slice(&body[1..]);
    }
    Ok(payload)
}

/// Inverse of `encode_record`, accepting only its exact output.
fn decode_record<D: RecordData>(payload: &[u8], line: u64) -> Result<(String, D), Error> {
    let malformed = |defect| Error::Malformed { line, defect };
    let rest = payload
        .strip_prefix(b"{\"path\":\"")
        .ok_or(malformed(Defect::NotCanonical))?;
    let end = rest
        .iter()
        .position(|&b| b == b'"')
        .ok_or(malformed(Defect::NotCanonical))?;
    let (raw, tail) = (&rest[..end], &rest[end + 1..]);
    if !raw.iter().all(|&b| b == b'0' || b == b'1') {
        return Err(Error::InvalidPath {
            line,
            path: String::from_utf8_lossy(raw).into_owned(),
        });
    }
    let body = match tail {
        [b'}'] => b"{}".to_vec(),
        [b',', fields @ ..] => [b"{".as_slice(), fields].concat(),
        _ => return Err(malformed(Defect::NotCanonical)),
    };
    let data: D =
        serde_json::from_slice(&body).map_err(|e| malformed(Defect::Json(e.to_string())))?;
    let path = String::from_utf8(raw.to_vec()).expect("binary digits are ASCII");
    match encode_record(&path, &data) {
        Ok(encoded) if encoded == payload => Ok((path, data)),
        _ => Err(malformed(Defect::NotCanonical)),
    }
}

/// Explains why a complete first line is not the expected header line.
fn diagnose_header(line: &[u8], expected: &Header) -> Error {
    let malformed = |defect| Error::Malformed { line: 1, defect };
    let payload = match unframe(line) {
        Ok(payload) => payload,
        Err(defect) => return malformed(defect),
    };
    match serde_json::from_slice::<Header>(payload) {
        Err(e) => malformed(Defect::Json(e.to_string())),
        Ok(found) if found == *expected => malformed(Defect::NotCanonical),
        Ok(found) => Error::WrongHeader {
            found: Box::new(found),
        },
    }
}

/// Checks a new path against the header, its predecessor and the recorded
/// paths, without changing anything.
fn admissible(
    paths: &BTreeSet<String>,
    previous: Option<&str>,
    header: &Header,
    path: &str,
    line: u64,
) -> Result<(), Error> {
    if !header.contains(path) {
        return Err(Error::InvalidPath {
            line,
            path: path.to_owned(),
        });
    }
    if let Some(previous) = previous.filter(|p| search_order(path) <= search_order(p)) {
        return Err(Error::OutOfOrder {
            line,
            path: path.to_owned(),
            previous: previous.to_owned(),
        });
    }
    // In a prefix-free set, an ancestor of `path` (or `path` itself) can only
    // be its lexical predecessor and a descendant only its lexical successor.
    let before = paths
        .range::<str, _>((Bound::Unbounded, Bound::Included(path)))
        .next_back();
    let after = paths
        .range::<str, _>((Bound::Included(path), Bound::Unbounded))
        .next();
    let other = before
        .filter(|p| path.starts_with(p.as_str()))
        .or(after.filter(|p| p.starts_with(path)));
    if let Some(other) = other {
        return Err(Error::Overlap {
            line,
            path: path.to_owned(),
            other: other.clone(),
        });
    }
    Ok(())
}

/// Visits the maximal boxes of the tree below `root` that neither are
/// recorded nor contain a record, in preorder. Iterative: the stack holds at
/// most one pending sibling per depth, so memory is bounded by the square of
/// the depth whatever the number of boxes visited.
fn complement(root: &str, paths: &BTreeSet<String>, mut visit: impl FnMut(String)) {
    let mut saved = paths.iter().peekable();
    let mut stack = vec![root.to_owned()];
    while let Some(path) = stack.pop() {
        match saved.peek() {
            Some(next) if **next == path => {
                saved.next();
            }
            Some(next) if next.starts_with(path.as_str()) => {
                stack.push(format!("{path}1"));
                stack.push(format!("{path}0"));
            }
            _ => visit(path),
        }
    }
    debug_assert!(saved.next().is_none(), "every record lies below the root");
}

/// The frontier as a list in search order.
fn frontier(root: &str, paths: &BTreeSet<String>) -> Vec<String> {
    let mut frontier = Vec::new();
    complement(root, paths, |path| frontier.push(path));
    frontier.sort_unstable_by(|a, b| search_order(a).cmp(&search_order(b)));
    frontier
}

/// The line `append` writes for `data` at `path`; refused unless the line
/// round-trips exactly (never write what a later load would refuse) and
/// fits in [`MAX_LINE`].
fn saved_line<D: RecordData>(path: &str, data: &D, line: u64) -> Result<Vec<u8>, Error> {
    let unrepresentable = |reason: String| Error::Unrepresentable {
        path: path.to_owned(),
        reason,
    };
    if !is_path(path) {
        return Err(unrepresentable("the path is not binary".into()));
    }
    let payload = encode_record(path, data).map_err(unrepresentable)?;
    decode_record::<D>(&payload, line).map_err(|e| unrepresentable(format!("does not round-trip: {e}")))?;
    let bytes = frame(&payload);
    if bytes.len() > MAX_LINE {
        return Err(unrepresentable(format!("line exceeds {MAX_LINE} bytes")));
    }
    Ok(bytes)
}

struct Parsed<D> {
    header_complete: bool,
    /// In file order, so strictly increasing (depth, path).
    records: Vec<(String, D)>,
    paths: BTreeSet<String>,
    valid_bytes: u64,
    unterminated_bytes: u64,
}

fn read_line(reader: &mut impl BufRead, line: u64) -> Result<Vec<u8>, Error> {
    let mut bytes = Vec::new();
    reader
        .by_ref()
        .take(MAX_LINE as u64 + 1)
        .read_until(b'\n', &mut bytes)
        .map_err(io_error("read certificate"))?;
    if bytes.len() > MAX_LINE {
        return Err(Error::TooLong { line });
    }
    Ok(bytes)
}

/// Every structural check, with no verification and no side effect. Only an
/// unterminated final line is set aside; any complete defect refuses.
fn parse<D: RecordData>(reader: &mut impl BufRead, expected: &Header) -> Result<Parsed<D>, Error> {
    let header = expected.line();
    let first = read_line(reader, 1)?;
    let mut parsed = Parsed {
        header_complete: false,
        records: Vec::new(),
        paths: BTreeSet::new(),
        valid_bytes: 0,
        unterminated_bytes: 0,
    };
    if first.last() != Some(&b'\n') {
        // Nothing can follow an unterminated line; an interrupted header write
        // is recognised only as a prefix of the exact expected header.
        if !header.starts_with(&first) {
            return Err(Error::WrongPartialHeader);
        }
        parsed.unterminated_bytes = first.len() as u64;
        return Ok(parsed);
    }
    if first != header {
        return Err(diagnose_header(&first, expected));
    }
    parsed.header_complete = true;
    parsed.valid_bytes = first.len() as u64;
    for line in 2.. {
        let bytes = read_line(reader, line)?;
        if bytes.is_empty() {
            break;
        }
        if bytes.last() != Some(&b'\n') {
            parsed.unterminated_bytes = bytes.len() as u64;
            break;
        }
        let payload = unframe(&bytes).map_err(|defect| Error::Malformed { line, defect })?;
        let (path, data) = decode_record::<D>(payload, line)?;
        let previous = parsed.records.last().map(|(p, _)| p.as_str());
        admissible(&parsed.paths, previous, expected, &path, line)?;
        parsed.paths.insert(path.clone());
        parsed.valid_bytes += bytes.len() as u64;
        parsed.records.push((path, data));
    }
    Ok(parsed)
}

/// Runs `verify` on every record with `threads` workers. On failure, reports
/// the earliest failing record in file order, whatever the thread count: a
/// record is skipped only when an earlier one has already failed.
fn verify_all<D, V>(records: &[(String, D)], threads: NonZeroUsize, verify: &V) -> Result<(), Error>
where
    D: Sync,
    V: Fn(&str, &D) -> Result<(), BoxError> + Sync,
{
    let next = AtomicUsize::new(0);
    let earliest = AtomicUsize::new(usize::MAX);
    let failures = Mutex::new(Vec::new());
    let spawned = thread::scope(|scope| {
        for n in 0..threads.get().min(records.len()) {
            thread::Builder::new()
                .name(format!("verify-{n}"))
                .spawn_scoped(scope, || loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    if index >= records.len() || index > earliest.load(Ordering::Relaxed) {
                        return;
                    }
                    let (path, data) = &records[index];
                    let outcome = catch_unwind(AssertUnwindSafe(|| verify(path, data)))
                        .unwrap_or_else(|panic| {
                            let message = panic_message(panic.as_ref());
                            Err(format!("verifier panicked: {message}").into())
                        });
                    if let Err(source) = outcome {
                        earliest.fetch_min(index, Ordering::Relaxed);
                        failures
                            .lock()
                            .unwrap_or_else(|e| e.into_inner())
                            .push((index, source));
                    }
                })
                .map_err(io_error("start a verification thread"))?;
        }
        Ok(())
    });
    let failure = failures
        .into_inner()
        .unwrap_or_else(|e| e.into_inner())
        .into_iter()
        .min_by_key(|(index, _)| *index);
    if let Some((index, source)) = failure {
        return Err(Error::Refuted {
            line: index as u64 + 2,
            path: records[index].0.clone(),
            source,
        });
    }
    spawned
}

/// Parses and verifies an open, locked file without changing it.
fn load<D, V>(
    file: &File,
    expected: &Header,
    threads: NonZeroUsize,
    verify: &V,
) -> Result<Parsed<D>, Error>
where
    D: RecordData,
    V: Fn(&str, &D) -> Result<(), BoxError> + Sync,
{
    let length = file
        .metadata()
        .map_err(io_error("read certificate metadata"))?
        .len();
    let parsed = parse::<D>(&mut BufReader::with_capacity(1 << 16, file), expected)?;
    if parsed.valid_bytes + parsed.unterminated_bytes != length {
        return Err(Error::Changed);
    }
    verify_all(&parsed.records, threads, verify)?;
    Ok(parsed)
}

fn inspection<D>(expected: &Header, parsed: &Parsed<D>) -> Inspection {
    let mut boxes = 0;
    complement(expected.root(), &parsed.paths, |_| boxes += 1);
    Inspection {
        header_complete: parsed.header_complete,
        records: parsed.records.len() as u64,
        last: parsed.records.last().map(|(path, _)| path.clone()),
        valid_bytes: parsed.valid_bytes,
        unterminated_bytes: parsed.unterminated_bytes,
        frontier: boxes,
    }
}

/// Refuses a path that exists and is not a regular file. Opening a FIFO for
/// reading would block, so this runs before the open (and again on the open
/// handle, which closes the gap for the writer).
fn regular(file: &Path) -> Result<(), Error> {
    match std::fs::metadata(file) {
        Ok(metadata) if !metadata.is_file() => Err(Error::NotAFile),
        _ => Ok(()),
    }
}

/// Checks a certificate without creating, repairing or changing it. Takes a
/// shared lock, so it refuses while a writer holds the file. Every record
/// passes `verify`, run on `threads` workers.
pub fn check<D, V>(
    file: &Path,
    expected: &Header,
    threads: NonZeroUsize,
    verify: V,
) -> Result<Inspection, Error>
where
    D: RecordData,
    V: Fn(&str, &D) -> Result<(), BoxError> + Sync,
{
    expected.validate()?;
    regular(file)?;
    let handle = File::open(file).map_err(io_error("open certificate"))?;
    FileExt::try_lock_shared(&handle).map_err(Error::Locked)?;
    let parsed = load::<D, V>(&handle, expected, threads, &verify)?;
    Ok(inspection(expected, &parsed))
}

/// The only writer of one certificate file, holding an exclusive lock for its
/// lifetime. Appends are buffered; `flush` makes them durable.
pub struct Store<D> {
    file: BufWriter<File>,
    header: Header,
    paths: BTreeSet<String>,
    records: u64,
    bytes: u64,
    last: Option<String>,
    failed: bool,
    data: PhantomData<fn(&D)>,
}

impl<D: RecordData> Store<D> {
    /// Opens or creates a certificate; [`Store::frontier`] then lists its
    /// unresolved frontier.
    ///
    /// Every complete line is checked and every record verified before
    /// anything is written. Then an interrupted header is rewritten, or an
    /// unterminated final line removed; nothing else is ever changed. The
    /// returned inspection describes the file as it was found.
    pub fn open<V>(
        file: &Path,
        expected: &Header,
        threads: NonZeroUsize,
        verify: V,
    ) -> Result<(Self, Inspection), Error>
    where
        V: Fn(&str, &D) -> Result<(), BoxError> + Sync,
    {
        expected.validate()?;
        regular(file)?;
        let mut handle = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(file)
            .map_err(io_error("open certificate"))?;
        let metadata = handle.metadata().map_err(io_error("read certificate metadata"))?;
        if !metadata.is_file() {
            return Err(Error::NotAFile);
        }
        FileExt::try_lock_exclusive(&handle).map_err(Error::Locked)?;
        let parsed = load::<D, V>(&handle, expected, threads, &verify)?;
        let found = inspection(expected, &parsed);
        let mut bytes = parsed.valid_bytes;
        if !parsed.header_complete {
            let header = expected.line();
            handle.set_len(0).map_err(io_error("reset interrupted header"))?;
            handle
                .seek(SeekFrom::Start(0))
                .map_err(io_error("seek certificate"))?;
            handle
                .write_all(&header)
                .map_err(io_error("write certificate header"))?;
            handle.sync_all().map_err(io_error("sync certificate header"))?;
            bytes = header.len() as u64;
        } else if parsed.unterminated_bytes > 0 {
            handle
                .set_len(parsed.valid_bytes)
                .map_err(io_error("remove unterminated line"))?;
            handle.sync_all().map_err(io_error("sync repaired certificate"))?;
        }
        // A creator may have died before its directory entry became durable;
        // closing that window on every writable open costs one sync.
        let parent = file
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        File::open(parent)
            .and_then(|directory| directory.sync_all())
            .map_err(io_error("sync certificate directory"))?;
        handle
            .seek(SeekFrom::End(0))
            .map_err(io_error("seek certificate end"))?;
        let store = Self {
            file: BufWriter::with_capacity(WRITE_BUFFER, handle),
            header: expected.clone(),
            paths: parsed.paths,
            records: parsed.records.len() as u64,
            bytes,
            last: found.last.clone(),
            failed: false,
            data: PhantomData,
        };
        Ok((store, found))
    }

    /// Buffers one record. The store does not verify it; every later
    /// reading does. Refuses
    /// a path that is invalid, out of order or overlapping, and data without
    /// a canonical form; a refusal changes nothing.
    pub fn append(&mut self, path: &str, data: &D) -> Result<(), Error> {
        if self.failed {
            return Err(Error::Failed);
        }
        let line = self.records + 2;
        admissible(&self.paths, self.last.as_deref(), &self.header, path, line)?;
        let bytes = saved_line(path, data, line)?;
        if let Err(source) = self.file.write_all(&bytes) {
            self.failed = true;
            return Err(Error::Io {
                action: "append certificate record",
                source,
            });
        }
        self.paths.insert(path.to_owned());
        self.records += 1;
        self.bytes += bytes.len() as u64;
        self.last = Some(path.to_owned());
        Ok(())
    }

    /// Writes buffered records and synchronises the file: every record
    /// appended before this call is durable when it returns.
    pub fn flush(&mut self) -> Result<(), Error> {
        if self.failed {
            return Err(Error::Failed);
        }
        if let Err(source) = self.file.flush().and_then(|()| self.file.get_ref().sync_all()) {
            self.failed = true;
            return Err(Error::Io {
                action: "synchronise certificate",
                source,
            });
        }
        Ok(())
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The minimal unresolved boxes of the records saved so far, ordered by
    /// depth and then path: the maximal boxes below the root that neither
    /// are recorded nor contain a record.
    pub fn frontier(&self) -> Vec<String> {
        frontier(self.header.root(), &self.paths)
    }

    pub fn report(&self) -> StoreReport {
        StoreReport {
            records: self.records,
            bytes: self.bytes,
            last: self.last.clone(),
        }
    }
}

#[cfg(test)]
mod tests;
