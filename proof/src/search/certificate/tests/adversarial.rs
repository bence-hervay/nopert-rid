//! Adversarial tests of the certificate layer:
//! a measure-based frontier oracle, resource limits, fault injection.

use super::*;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::AtomicU64;

static SERIAL: AtomicU64 = AtomicU64::new(0);

struct Folder(PathBuf);

impl Folder {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "rid-certificate-review-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    fn file(&self) -> PathBuf {
        self.0.join("certificate")
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Tag {
    tag: u32,
}

impl RecordData for Tag {}

fn fnv(text: &str) -> u64 {
    text.bytes()
        .fold(0xcbf29ce484222325, |h, b| (h ^ u64::from(b)).wrapping_mul(0x100000001b3))
}

fn tag(path: &str) -> Tag {
    Tag {
        tag: (fnv(path) % 1000) as u32,
    }
}

fn verify_tag(path: &str, data: &Tag) -> Result<(), BoxError> {
    if *data == tag(path) {
        Ok(())
    } else {
        Err("false tag".into())
    }
}

fn one() -> NonZeroUsize {
    NonZeroUsize::new(1).unwrap()
}

fn policy() -> String {
    "c".repeat(64)
}

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn certificate_bytes(header: &Header, paths: &[String]) -> Vec<u8> {
    let mut bytes = header.line();
    for path in paths {
        bytes.extend(frame(&encode_record(path, &tag(path)).unwrap()));
    }
    bytes
}

/// Measure-based oracle, independent of any tree traversal: records and
/// frontier together are prefix-free, lie below the root, and their dyadic
/// measures add up to the root's measure (so they tile it); every non-root
/// frontier box has a record inside its parent (maximality); search order.
fn assert_tiling(root: &str, records: &[String], frontier: &[String]) {
    const D: usize = 127;
    let all: Vec<&String> = records.iter().chain(frontier).collect();
    for (i, a) in all.iter().enumerate() {
        assert!(a.starts_with(root) && a.len() <= D, "{a} outside {root}");
        for b in &all[i + 1..] {
            assert!(!a.starts_with(b.as_str()) && !b.starts_with(a.as_str()), "{a} and {b}");
        }
    }
    let total: u128 = all.iter().map(|p| 1u128 << (D - p.len())).sum();
    assert_eq!(total, 1u128 << (D - root.len()), "measure");
    for f in frontier {
        if f.as_str() != root {
            let parent = &f[..f.len() - 1];
            assert!(records.iter().any(|r| r.starts_with(parent)), "{f} is not maximal");
        }
    }
    assert!(frontier.windows(2).all(|w| search_order(&w[0]) < search_order(&w[1])));
}

/// Random prefix-free leaves of a random deep subdivision below `root`.
fn deep_leaves(state: &mut u64, root: &str, max_depth: usize) -> Vec<String> {
    let mut leaves = vec![root.to_owned()];
    // Split mostly along one random branch to reach large depths.
    for _ in 0..(20 + random(state) % 60) {
        let at = if random(state) % 4 == 0 {
            random(state) as usize % leaves.len()
        } else {
            leaves
                .iter()
                .enumerate()
                .max_by_key(|(_, p)| p.len())
                .map(|(i, _)| i)
                .unwrap()
        };
        let path = leaves.swap_remove(at);
        if path.len() < max_depth {
            let first = random(state) % 2;
            leaves.push(format!("{path}{first}"));
            leaves.push(format!("{path}{}", 1 - first));
        } else {
            leaves.push(path);
        }
    }
    leaves
}

#[test]
fn frontiers_of_deep_random_certificates_tile_the_root_by_measure() {
    let _files = super::files();
    let dir = Folder::new();
    let mut state = 0x0bad_5eed_1234_5678;
    for round in 0..400 {
        let root: String = (0..random(&mut state) % 5)
            .map(|_| if random(&mut state) % 2 == 0 { '0' } else { '1' })
            .collect();
        let max_depth = 20 + (random(&mut state) as usize % 100);
        let header = Header::new(&policy(), &root, max_depth).unwrap();
        let mut leaves = deep_leaves(&mut state, &root, max_depth);
        leaves.retain(|_| random(&mut state) % 3 != 0);
        leaves.sort_unstable_by(|a, b| search_order(a).cmp(&search_order(b)));
        let bytes = certificate_bytes(&header, &leaves);
        let cut = if round % 2 == 0 {
            bytes.len()
        } else {
            random(&mut state) as usize % (bytes.len() + 1)
        };
        fs::write(dir.file(), &bytes[..cut]).unwrap();
        let checked = check(&dir.file(), &header, one(), verify_tag).unwrap();
        let kept = &leaves[..checked.records as usize];
        // The frontier of the file after `open`, counted by `check`.
        let (store, found) = Store::<Tag>::open(&dir.file(), &header, one(), verify_tag).unwrap();
        let frontier = store.frontier();
        assert_tiling(&root, kept, &frontier);
        assert_eq!((found.frontier, checked.frontier), (frontier.len() as u64, found.frontier));
        drop(store);
    }
}

/// `check` counts the frontier without listing it, so a small certificate of
/// deep records (whose frontier as strings would take about depth²/2 bytes per
/// record, here over 100 MB) is checked in memory proportional to its size.
/// The count is checked against an independent formula: the frontier of a
/// prefix-free set below the root "" has `I - L + 1` boxes, with `L` records
/// and `I` distinct proper prefixes of records (the inner nodes of the tree
/// they span; each has two children, and all but the root have a parent).
#[test]
fn deep_frontiers_are_counted_without_being_listed() {
    let _files = super::files();
    let dir = Folder::new();
    let header = Header::new(&policy(), "", MAX_DEPTH).unwrap();
    let mut state = 0x1357_9bdf_2468_ace0;
    let mut paths: Vec<String> = (0..16)
        .map(|_| {
            (0..MAX_DEPTH)
                .map(|_| if random(&mut state) % 2 == 0 { '0' } else { '1' })
                .collect()
        })
        .collect();
    paths.sort_unstable_by(|a, b| search_order(a).cmp(&search_order(b)));
    paths.dedup();
    fs::write(dir.file(), certificate_bytes(&header, &paths)).unwrap();
    let checked = check(&dir.file(), &header, one(), verify_tag).unwrap();
    // Distinct nonempty prefixes: in lexical order, each path adds those
    // longer than its common prefix with the previous one. The proper ones
    // are these without the L records themselves, plus the empty root.
    let mut lexical = paths.clone();
    lexical.sort();
    let mut prefixes = 0;
    let mut previous = "";
    for path in &lexical {
        let common = path.bytes().zip(previous.bytes()).take_while(|(a, b)| a == b).count();
        prefixes += path.len() - common;
        previous = path;
    }
    let inner = prefixes - paths.len() + 1;
    assert_eq!(checked.frontier, (inner - paths.len() + 1) as u64);
    assert_eq!(checked.records, paths.len() as u64);
}

#[test]
fn oversized_record_data_is_refused_on_append_without_change() {
    let _files = super::files();
    #[derive(Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Big {
        bytes: Vec<u8>,
    }
    impl RecordData for Big {}
    let dir = Folder::new();
    let header = Header::new(&policy(), "", 64).unwrap();
    let (mut store, _) = Store::<Big>::open(&dir.file(), &header, one(), |_, _| Ok(())).unwrap();
    let before = store.report();
    let big = Big {
        bytes: vec![255; MAX_LINE / 3],
    };
    assert!(matches!(store.append("0", &big), Err(Error::Unrepresentable { .. })));
    assert_eq!(store.report(), before);
    // Just below the limit is written and read back.
    let fits = Big {
        bytes: vec![1; (MAX_LINE - 64) / 2],
    };
    store.append("0", &fits).unwrap();
    store.flush().unwrap();
    drop(store);
    let checked = check::<Big, _>(&dir.file(), &header, one(), |_, _| Ok(())).unwrap();
    assert_eq!(checked.records, 1);
}

#[test]
fn flattened_inner_path_fields_are_refused() {
    let _files = super::files();
    #[derive(Serialize, Deserialize)]
    struct Inner {
        path: String,
    }
    #[derive(Serialize, Deserialize)]
    struct Outer {
        #[serde(flatten)]
        inner: Inner,
        x: u8,
    }
    impl RecordData for Outer {}
    let dir = Folder::new();
    let header = Header::new(&policy(), "", 64).unwrap();
    let (mut store, _) =
        Store::<Outer>::open(&dir.file(), &header, one(), |_, _| Ok(())).unwrap();
    let value = Outer {
        inner: Inner { path: "1".into() },
        x: 1,
    };
    assert!(matches!(store.append("0", &value), Err(Error::Unrepresentable { .. })));
    // A hand-written line whose data carries a second "path" is refused.
    drop(store);
    let mut bytes = header.line();
    bytes.extend(frame(br#"{"path":"0","path":"1","x":1}"#));
    fs::write(dir.file(), &bytes).unwrap();
    assert!(check::<Outer, _>(&dir.file(), &header, one(), |_, _| Ok(())).is_err());
}

#[test]
fn invalid_utf8_inside_a_payload_with_a_correct_checksum_is_refused() {
    let header = Header::new(&policy(), "", 64).unwrap();
    let mut bytes = header.line();
    bytes.extend(frame(b"{\"path\":\"0\",\"tag\":1\xff}"));
    assert!(matches!(
        parse::<Tag>(&mut &bytes[..], &header),
        Err(Error::Malformed { line: 2, .. })
    ));
}

#[test]
fn a_missing_parent_directory_creates_nothing() {
    let dir = Folder::new();
    let file = dir.0.join("missing").join("certificate");
    let header = Header::new(&policy(), "", 64).unwrap();
    assert!(matches!(
        Store::<Tag>::open(&file, &header, one(), verify_tag),
        Err(Error::Io { .. })
    ));
    assert!(!dir.0.join("missing").exists());
}

const FSIZE_CHILD: &str = "RID_REVIEW_FSIZE_CHILD";

/// A write failure (the file-size limit, EFBIG) poisons the handle, and the
/// file left behind recovers to a valid prefix of the appended records.
#[test]
fn a_failed_write_poisons_the_store_and_leaves_a_recoverable_file() {
    let _processes = super::processes();
    let dir = Folder::new();
    let limit = 4096;
    // The child's output goes to pipes: the size limit would also apply to
    // an inherited standard output that is a file.
    let child = "search::certificate::tests::adversarial::fsize_child";
    let output = Command::new("sh")
        .arg("-c")
        .arg("trap '' XFSZ; exec prlimit --fsize=$1 -- \"$0\" --exact $2 --test-threads=1")
        .arg(std::env::current_exe().unwrap())
        .arg(limit.to_string())
        .arg(child)
        .env(FSIZE_CHILD, dir.file())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let report = fs::read_to_string(dir.0.join("report")).unwrap();
    assert_eq!(report, "flush:io append:failed flush:failed");
    let bytes = fs::read(dir.file()).unwrap();
    assert_eq!(bytes.len(), limit);
    let header = Header::new(&policy(), "", 64).unwrap();
    let (store, found) = Store::<Tag>::open(&dir.file(), &header, one(), verify_tag).unwrap();
    assert!(found.unterminated_bytes > 0 || bytes.ends_with(b"\n"));
    let expected: Vec<String> = (0..found.records).map(|n| format!("{n:012b}")).collect();
    assert_eq!(found.last, expected.last().cloned());
    drop(store);
    // After recovery the file is exactly the complete prefix.
    assert_eq!(fs::read(dir.file()).unwrap(), certificate_bytes(&header, &expected));
}

/// Runs only as the child of the test above, under a file-size limit.
#[test]
fn fsize_child() {
    let Some(file) = std::env::var_os(FSIZE_CHILD) else {
        return;
    };
    let file = PathBuf::from(file);
    let header = Header::new(&policy(), "", 64).unwrap();
    let (mut store, _) = Store::<Tag>::open(&file, &header, one(), verify_tag).unwrap();
    for n in 0..400u32 {
        let path = format!("{n:012b}");
        store.append(&path, &tag(&path)).unwrap();
    }
    let mut report = Vec::new();
    report.push(match store.flush() {
        Err(Error::Io { .. }) => "flush:io",
        Ok(()) => "flush:ok",
        Err(_) => "flush:other",
    });
    let path = format!("{:012b}", 400);
    report.push(match store.append(&path, &tag(&path)) {
        Err(Error::Failed) => "append:failed",
        _ => "append:other",
    });
    report.push(match store.flush() {
        Err(Error::Failed) => "flush:failed",
        _ => "flush:other",
    });
    drop(store);
    fs::write(file.parent().unwrap().join("report"), report.join(" ")).unwrap();
}

/// A certificate path naming a FIFO or a directory is refused at once by
/// `check` and `Store::open` (a FIFO would otherwise block them forever).
#[test]
fn non_regular_files_are_refused_promptly() {
    let processes = super::processes();
    let dir = Folder::new();
    assert!(Command::new("mkfifo").arg(dir.file()).status().unwrap().success());
    drop(processes);
    let header = Header::new(&policy(), "", 64).unwrap();
    let (fifo, directory) = (dir.file(), dir.0.clone());
    let cases = [(fifo.clone(), false), (fifo, true)]
        .into_iter()
        .chain([(directory.clone(), false), (directory, true)]);
    for (file, writer) in cases {
        let (send, receive) = std::sync::mpsc::channel();
        let header = header.clone();
        std::thread::spawn(move || {
            let refused = if writer {
                let opened = Store::<Tag>::open(&file, &header, one(), verify_tag);
                matches!(opened, Err(Error::NotAFile))
            } else {
                matches!(check::<Tag, _>(&file, &header, one(), verify_tag), Err(Error::NotAFile))
            };
            let _ = send.send(refused);
        });
        let refused = receive.recv_timeout(std::time::Duration::from_secs(20));
        assert_eq!(refused, Ok(true), "writer {writer}");
    }
}
