use super::*;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::AtomicU64;
use std::sync::{Arc, Barrier, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::time::Duration;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Spawning a process (fork, then exec) lends every open file of this process
/// to the child until its exec, including a certificate's lock: a lock can
/// seem held just after its handle was dropped. Tests that spawn processes
/// hold this guard exclusively; tests that use certificate files, shared.
static PROCESSES: RwLock<()> = RwLock::new(());

pub(super) fn files() -> RwLockReadGuard<'static, ()> {
    PROCESSES.read().unwrap_or_else(|e| e.into_inner())
}

pub(super) fn processes() -> RwLockWriteGuard<'static, ()> {
    PROCESSES.write().unwrap_or_else(|e| e.into_inner())
}

static SERIAL: AtomicU64 = AtomicU64::new(0);

/// A fresh directory under the system temporary directory, removed on drop.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "rid-certificate-tests-{}-{}",
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

/// Shaped like component record data: an internally tagged enum with nested
/// structs, arrays, integers, names with signs and a unit variant.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "component", deny_unknown_fields)]
enum Data {
    Domain { inequality: u8 },
    Global { witness: Witness },
    Local { cover: u32 },
    Exotic { cover: CoverName },
    Marker,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Witness {
    edge: [u8; 2],
    vertex: u8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
enum CoverName {
    #[serde(rename = "square")]
    Square,
    #[serde(rename = "pentagon")]
    Pentagon,
    #[serde(rename = "arc+")]
    ArcPlus,
    #[serde(rename = "arc-")]
    ArcMinus,
    #[serde(rename = "endpoint+")]
    EndpointPlus,
    #[serde(rename = "crossing-")]
    CrossingMinus,
}

impl RecordData for Data {}

const COVERS: [CoverName; 6] = [
    CoverName::Square,
    CoverName::Pentagon,
    CoverName::ArcPlus,
    CoverName::ArcMinus,
    CoverName::EndpointPlus,
    CoverName::CrossingMinus,
];

fn hash(text: &str) -> u64 {
    text.bytes().fold(0xcbf29ce484222325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
    })
}

/// The only data `verify` accepts for `path`.
fn data_for(path: &str) -> Data {
    let h = hash(path);
    let pick = |shift: u32, modulus: u64| ((h >> shift) % modulus) as u8;
    match h % 5 {
        0 => Data::Domain {
            inequality: pick(8, 73),
        },
        1 => Data::Global {
            witness: Witness {
                edge: [pick(8, 120), pick(16, 120)],
                vertex: pick(24, 60),
            },
        },
        2 => Data::Local {
            cover: u32::from(pick(8, 40)),
        },
        3 => Data::Exotic {
            cover: COVERS[usize::from(pick(8, 6))],
        },
        _ => Data::Marker,
    }
}

fn verify(path: &str, data: &Data) -> Result<(), BoxError> {
    if *data == data_for(path) {
        Ok(())
    } else {
        Err("false record".into())
    }
}

fn one() -> NonZeroUsize {
    NonZeroUsize::new(1).unwrap()
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn header() -> Header {
    Header::new(&"a".repeat(64), "", MAX_DEPTH).unwrap()
}

fn record_line(path: &str, data: &Data) -> Vec<u8> {
    frame(&encode_record(path, data).unwrap())
}

fn certificate_with(header: &Header, paths: &[&str]) -> Vec<u8> {
    let mut bytes = header.line();
    for path in paths {
        bytes.extend(record_line(path, &data_for(path)));
    }
    bytes
}

fn certificate(paths: &[&str]) -> Vec<u8> {
    certificate_with(&header(), paths)
}

/// Lines framed with a correct checksum around an arbitrary payload.
fn with_payload(payload: &str) -> Vec<u8> {
    let mut bytes = certificate(&["0"]);
    bytes.extend(frame(payload.as_bytes()));
    bytes
}

fn parse_bytes(bytes: &[u8]) -> Result<Parsed<Data>, Error> {
    parse::<Data>(&mut &bytes[..], &header())
}

fn open(file: &Path) -> Result<(Store<Data>, Inspection), Error> {
    Store::open(file, &header(), one(), verify)
}

fn inspect(file: &Path) -> Result<Inspection, Error> {
    check(file, &header(), one(), verify)
}

/// The file's bytes and modification time, to prove it was not touched.
fn snapshot(file: &Path) -> (Vec<u8>, std::time::SystemTime) {
    (
        fs::read(file).unwrap(),
        fs::metadata(file).unwrap().modified().unwrap(),
    )
}

/// Independent oracle, not calling the production complement: every cell at
/// `depth` below `root` lies in exactly one recorded or frontier box, and no
/// cell outside `root` lies in any; every frontier box other than the root has
/// a record in its parent (so it is maximal); the frontier is in search order.
fn assert_minimal_cover(records: &[String], frontier: &[String], depth: usize, root: &str) {
    let all: Vec<&String> = records.iter().chain(frontier).collect();
    for n in 0..(1usize << depth) {
        let cell = format!("{n:0depth$b}");
        let expected = usize::from(cell.starts_with(root));
        let count = all.iter().filter(|p| cell.starts_with(p.as_str())).count();
        assert_eq!(count, expected, "cell {cell}");
    }
    for path in frontier {
        if path.as_str() != root {
            let parent = &path[..path.len() - 1];
            assert!(
                records.iter().any(|r| r.starts_with(parent)),
                "frontier box {path} could be merged with its sibling"
            );
        }
    }
    assert!(frontier.windows(2).all(|p| search_order(&p[0]) < search_order(&p[1])));
}

/// The frontier of a certificate file as a list, read without changing it:
/// the production complement over the parsed records. Its length must be the
/// count `check` reports.
fn listed<D: RecordData>(file: &Path, header: &Header) -> Vec<String> {
    let mut reader = std::io::BufReader::new(File::open(file).unwrap());
    let parsed = parse::<D>(&mut reader, header).unwrap();
    frontier(header.root(), &parsed.paths)
}

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// A random prefix-free leaf set of a random subdivision below `root`.
fn random_leaves(state: &mut u64, root: &str, depth: usize) -> Vec<String> {
    let mut leaves = vec![root.to_owned()];
    for _ in 0..24 {
        let at = random(state) as usize % leaves.len();
        let path = leaves.swap_remove(at);
        if path.len() < depth {
            leaves.push(format!("{path}0"));
            leaves.push(format!("{path}1"));
        } else {
            leaves.push(path);
        }
    }
    leaves.sort_unstable_by(|a, b| search_order(a).cmp(&search_order(b)));
    leaves
}

// ---------------------------------------------------------------------------
// Framing and canonical form
// ---------------------------------------------------------------------------

#[test]
fn checksum_is_the_first_32_bits_of_sha256_in_lowercase_hex() {
    // Published SHA-256 test vectors.
    assert_eq!(&checksum(b"abc"), b"ba7816bf");
    assert_eq!(&checksum(b""), b"e3b0c442");
    assert_eq!(
        &checksum(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        b"248d6a61"
    );
}

#[test]
fn header_line_is_exact() {
    let policy = "0123456789abcdef".repeat(4);
    let header = Header::new(&policy, "0110", 64).unwrap();
    let payload = format!(
        r#"{{"format":"rid-certificate/1","policy":"{policy}","root":"0110","max_depth":64}}"#
    );
    let expected = format!(
        "{} {payload}\n",
        std::str::from_utf8(&checksum(payload.as_bytes())).unwrap()
    );
    assert_eq!(header.line(), expected.as_bytes());
    assert_eq!(
        (header.format(), header.policy(), header.root(), header.max_depth()),
        (Format::V1, policy.as_str(), "0110", 64)
    );
}

#[test]
fn record_payloads_put_the_path_first_then_the_data_fields() {
    let cases = [
        (
            "00011010110",
            Data::Domain { inequality: 13 },
            r#"{"path":"00011010110","component":"Domain","inequality":13}"#,
        ),
        (
            "0001101011010",
            Data::Global {
                witness: Witness {
                    edge: [16, 17],
                    vertex: 3,
                },
            },
            concat!(
                r#"{"path":"0001101011010","component":"Global","#,
                r#""witness":{"edge":[16,17],"vertex":3}}"#
            ),
        ),
        (
            "01101001000110",
            Data::Local { cover: 9 },
            r#"{"path":"01101001000110","component":"Local","cover":9}"#,
        ),
        (
            "0110100100011101",
            Data::Exotic {
                cover: CoverName::ArcMinus,
            },
            r#"{"path":"0110100100011101","component":"Exotic","cover":"arc-"}"#,
        ),
        ("", Data::Marker, r#"{"path":"","component":"Marker"}"#),
    ];
    for (path, data, payload) in cases {
        assert_eq!(encode_record(path, &data).unwrap(), payload.as_bytes());
        let (decoded, value) = decode_record::<Data>(payload.as_bytes(), 2).unwrap();
        assert_eq!((decoded.as_str(), value), (path, data));
    }
    // Data serialising to an empty object gives a path-only record.
    #[derive(Debug, PartialEq, Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Nothing {}
    impl RecordData for Nothing {}
    assert_eq!(encode_record("01", &Nothing {}).unwrap(), br#"{"path":"01"}"#);
    let (path, value) = decode_record::<Nothing>(br#"{"path":"01"}"#, 2).unwrap();
    assert_eq!((path.as_str(), value), ("01", Nothing {}));
    assert!(decode_record::<Nothing>(br#"{"path":"01",}"#, 2).is_err());
}

#[test]
fn random_records_round_trip() {
    let mut state = 0x5eed_1234_abcd_0001;
    for _ in 0..2000 {
        let depth = random(&mut state) as usize % 64;
        let path: String = (0..depth)
            .map(|_| if random(&mut state) % 2 == 0 { '0' } else { '1' })
            .collect();
        let data = data_for(&format!("{path}{}", random(&mut state)));
        let line = record_line(&path, &data);
        let payload = unframe(&line).unwrap();
        assert_eq!(decode_record::<Data>(payload, 2).unwrap(), (path, data));
    }
}

#[test]
fn header_constructor_refuses_invalid_fields() {
    let good = "a".repeat(64);
    for (policy, root, max_depth) in [
        ("a".repeat(63), "", 10),
        ("a".repeat(65), "", 10),
        ("A".repeat(64), "", 10),
        ("g".repeat(64), "", 10),
        (format!("{} ", "a".repeat(63)), "", 10),
        (good.clone(), "", MAX_DEPTH + 1),
        (good.clone(), "2", 10),
        (good.clone(), "0 ", 10),
        (good.clone(), "0101", 3),
    ] {
        assert!(
            matches!(Header::new(&policy, root, max_depth), Err(Error::InvalidHeader(_))),
            "{policy:?} {root:?} {max_depth}"
        );
    }
    assert!(Header::new(&good, &"1".repeat(MAX_DEPTH), MAX_DEPTH).is_ok());
    assert!(Header::new(&good, "", 0).is_ok());
}

#[test]
fn a_deserialised_invalid_header_is_refused_before_any_file_is_touched() {
    let _files = files();
    let dir = Folder::new();
    let payload = serde_json::to_string(&header()).unwrap();
    for text in [
        payload.replace("aaaa", "AAAA"),
        payload.replace("4096", "5000"),
        payload.replace(r#""root":"""#, r#""root":"2""#),
    ] {
        let invalid: Header = serde_json::from_str(&text).unwrap();
        assert!(matches!(
            Store::<Data>::open(&dir.file(), &invalid, one(), verify),
            Err(Error::InvalidHeader(_))
        ));
        assert!(!dir.file().exists());
        fs::write(dir.file(), invalid.line()).unwrap();
        assert!(matches!(
            check::<Data, _>(&dir.file(), &invalid, one(), verify),
            Err(Error::InvalidHeader(_))
        ));
        fs::remove_file(dir.file()).unwrap();
    }
}

#[test]
fn wrong_checksums_and_frames_are_refused() {
    // A record whose checksum has a letter, so that its uppercase differs.
    let (payload, good) = (0..)
        .map(|n| {
            let path = format!("{n:b}");
            let payload = encode_record(&path, &data_for(&path)).unwrap();
            let payload = String::from_utf8(payload).unwrap();
            let good = std::str::from_utf8(&checksum(payload.as_bytes())).unwrap().to_owned();
            (payload, good)
        })
        .find(|(_, good)| good.bytes().any(|b| b.is_ascii_lowercase()))
        .unwrap();
    let flipped = format!("{}{}", &good[..7], if &good[7..] == "0" { "1" } else { "0" });
    let cases = [
        (format!("{flipped} {payload}\n"), Defect::Checksum),
        (format!("{} {payload}\n", good.to_uppercase()), Defect::Frame),
        (format!("{} {payload}\n", &good[..7]), Defect::Frame),
        (format!("{good}0 {payload}\n"), Defect::Frame),
        (format!("{good}  {payload}\n"), Defect::Checksum),
        (format!("{good}\t{payload}\n"), Defect::Frame),
        (format!("{good}{payload}\n"), Defect::Frame),
        (format!("{good} \n"), Defect::Checksum),
        ("\n".to_owned(), Defect::Frame),
        (format!(" {good} {payload}\n"), Defect::Frame),
        (format!("{good} {payload}\r\n"), Defect::Checksum),
    ];
    for (line, defect) in cases {
        let mut bytes = certificate(&[]);
        bytes.extend_from_slice(line.as_bytes());
        match parse_bytes(&bytes) {
            Err(Error::Malformed { line: 2, defect: found }) => {
                assert_eq!(found, defect, "{line:?}")
            }
            other => panic!("{line:?} gave {:?}", other.map(|p| p.records.len())),
        }
    }
}

#[test]
fn non_canonical_payloads_are_refused_even_with_correct_checksums() {
    let canonical = r#"{"path":"1","component":"Domain","inequality":13}"#;
    assert!(parse_bytes(&with_payload(canonical)).is_ok());
    for payload in [
        r#"{"path":"1","component":"Domain","inequality":13} "#,
        r#" {"path":"1","component":"Domain","inequality":13}"#,
        r#"{"path": "1","component":"Domain","inequality":13}"#,
        r#"{"path":"1", "component":"Domain","inequality":13}"#,
        r#"{"path":"1","component":"Domain","inequality":13 }"#,
        r#"{"path":"1","inequality":13,"component":"Domain"}"#,
        r#"{"component":"Domain","path":"1","inequality":13}"#,
        r#"{"component":"Domain","inequality":13,"path":"1"}"#,
        r#"{"path":"1","component":"Domain","inequality":13.0}"#,
        r#"{"path":"1","component":"Domain","inequality":1.3e1}"#,
        r#"{"path":"1","component":"Domain","inequality":013}"#,
        r#"{"path":"1","component":"Domain","inequality":+13}"#,
        r#"{"path":"1","component":"Domain","inequality":"13"}"#,
        r#"{"path":"1","component":"Domain","inequality":13,"extra":0}"#,
        r#"{"path":"1","component":"Domain","inequality":13,"inequality":13}"#,
        r#"{"path":"1","component":"Domain","component":"Domain","inequality":13}"#,
        r#"{"path":"1","path":"1","component":"Domain","inequality":13}"#,
        r#"{"path":"1","component":"\u0044omain","inequality":13}"#,
        r#"{"path":"1","component":"domain","inequality":13}"#,
        r#"{"path":"1","component":"Domain","inequality":13}{}"#,
        r#"{"path":"1","component":"Domain","inequality":13},"#,
        r#"{"path":"1","component":"Domain"}"#,
        r#"{"path":"1","component":"Domain","inequality":null}"#,
        r#"{"path":"1","component":"Domain","inequality":-13}"#,
        r#"{"path":"1","component":"Domain","inequality":256}"#,
        r#"{"path":"1"}"#,
        r#"{"path":"1",}"#,
        r#"{"path":"1""#,
        r#"{"path":1,"component":"Domain","inequality":13}"#,
        r#"{"path":null,"component":"Domain","inequality":13}"#,
        r#"{"component":"Domain","inequality":13}"#,
        r#"[{"path":"1","component":"Domain","inequality":13}]"#,
        r#""{\"path\":\"1\"}""#,
        "",
        "{}",
        "null",
        // Unit variants of internally tagged enums ignore unknown fields in
        // serde; only the canonical comparison refuses this one.
        r#"{"path":"1","component":"Marker","x":1}"#,
        r#"{"path":"1","component":"Exotic","cover":"crossing-","cover":"crossing-"}"#,
        r#"{"path":"1","component":"Exotic","cover":"cro\u0073sing-"}"#,
        r#"{"path":"1","component":"Global","witness":{"vertex":3,"edge":[16,17]}}"#,
        r#"{"path":"1","component":"Global","witness":{"edge":[16, 58],"vertex":3}}"#,
        r#"{"path":"1","component":"Global","witness":{"edge":[16,17,1],"vertex":3}}"#,
    ] {
        match parse_bytes(&with_payload(payload)) {
            Err(Error::Malformed { line: 3, .. }) => {}
            other => panic!("{payload} gave {:?}", other.map(|p| p.records.len())),
        }
    }
    // Escaped or non-binary path spellings are not paths.
    for payload in [
        r#"{"path":"\u0031","component":"Domain","inequality":13}"#,
        r#"{"path":"2","component":"Domain","inequality":13}"#,
        r#"{"path":"1 ","component":"Domain","inequality":13}"#,
        r#"{"path":"1\n","component":"Domain","inequality":13}"#,
    ] {
        match parse_bytes(&with_payload(payload)) {
            Err(Error::InvalidPath { line: 3, .. }) => {}
            other => panic!("{payload} gave {:?}", other.map(|p| p.records.len())),
        }
    }
}

#[test]
fn data_without_a_canonical_object_form_is_never_written() {
    let _files = files();
    fn refused<D: RecordData>(value: D) {
        let dir = Folder::new();
        let (mut store, _) =
            Store::<D>::open(&dir.file(), &header(), one(), |_, _| Ok(())).unwrap();
        let before = store.report();
        assert!(
            matches!(store.append("0", &value), Err(Error::Unrepresentable { .. })),
            "{}",
            std::any::type_name::<D>()
        );
        assert_eq!(store.report(), before);
        store.flush().unwrap();
        drop(store);
        assert_eq!(fs::read(dir.file()).unwrap(), header().line());
    }
    #[derive(Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct HasPath {
        path: String,
    }
    impl RecordData for HasPath {}
    refused(HasPath { path: "1".into() });

    #[derive(Serialize, Deserialize)]
    struct Number(u32);
    impl RecordData for Number {}
    refused(Number(3));

    #[derive(Serialize, Deserialize)]
    struct Pair(u8, u8);
    impl RecordData for Pair {}
    refused(Pair(1, 2));

    // Serialises a field its deserialiser does not accept.
    #[derive(Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Asymmetric {
        #[serde(rename(serialize = "a", deserialize = "b"))]
        value: u8,
    }
    impl RecordData for Asymmetric {}
    refused(Asymmetric { value: 1 });

    // A NaN serialises as null, which does not read back as a float.
    #[derive(Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Float {
        x: f64,
    }
    impl RecordData for Float {}
    refused(Float { x: f64::NAN });

    // Reads back, but serialises differently the second time.
    #[derive(Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Lossy {
        x: f64,
    }
    impl RecordData for Lossy {}
    let dir = Folder::new();
    let (mut store, _) =
        Store::<Lossy>::open(&dir.file(), &header(), one(), |_, _| Ok(())).unwrap();
    store.append("0", &Lossy { x: 0.1 }).unwrap();
    store.append("1", &Lossy { x: 1e300 }).unwrap();
    store.flush().unwrap();
    drop(store);
    let checked = check::<Lossy, _>(&dir.file(), &header(), one(), |_, _| Ok(())).unwrap();
    assert!(checked.complete());
}

// ---------------------------------------------------------------------------
// Interruption, corruption and recovery
// ---------------------------------------------------------------------------

const MIXED: [&str; 6] = ["00", "10", "011", "110", "0101", "1111"];

fn line_ends(bytes: &[u8]) -> Vec<usize> {
    bytes
        .iter()
        .enumerate()
        .filter_map(|(n, &b)| (b == b'\n').then_some(n + 1))
        .collect()
}

#[test]
fn every_byte_cut_keeps_exactly_the_complete_prefix() {
    let _files = files();
    let dir = Folder::new();
    for root in ["", "1"] {
        let header = Header::new(&"a".repeat(64), root, 16).unwrap();
        let paths: Vec<&str> = MIXED.iter().copied().filter(|p| p.starts_with(root)).collect();
        let original = certificate_with(&header, &paths);
        let ends = line_ends(&original);
        for cut in 0..=original.len() {
            fs::write(dir.file(), &original[..cut]).unwrap();
            let checked = check(&dir.file(), &header, one(), verify).unwrap();
            assert_eq!(fs::read(dir.file()).unwrap(), &original[..cut]);
            let retained = ends.iter().skip(1).filter(|&&end| end <= cut).count();
            let complete_lines = ends.iter().filter(|&&end| end <= cut).count();
            let valid = if complete_lines == 0 { 0 } else { ends[complete_lines - 1] };
            assert_eq!(
                (
                    checked.header_complete,
                    checked.records as usize,
                    checked.valid_bytes as usize,
                    checked.unterminated_bytes as usize,
                    checked.last.as_deref(),
                ),
                (
                    cut >= ends[0],
                    retained,
                    valid,
                    cut - valid,
                    retained.checked_sub(1).map(|i| paths[i]),
                ),
                "cut {cut}"
            );
            let records: Vec<String> = paths[..retained].iter().map(|p| p.to_string()).collect();
            let frontier = listed::<Data>(&dir.file(), &header);
            assert_eq!(frontier.len() as u64, checked.frontier);
            assert_minimal_cover(&records, &frontier, 5, root);
            assert_eq!(checked.complete(), cut == original.len() && frontier.is_empty());
            let (store, found) = Store::<Data>::open(&dir.file(), &header, one(), verify).unwrap();
            assert_eq!(found, checked);
            assert_eq!(store.report().records as usize, retained);
            drop(store);
            let kept = if cut < ends[0] { ends[0] } else { valid };
            assert_eq!(fs::read(dir.file()).unwrap(), &original[..kept], "cut {cut}");
        }
    }
}

#[test]
fn every_single_byte_mutation_is_refused() {
    let original = certificate(&MIXED);
    let last = original.len() - 1;
    let last_line = original.len() - line_ends(&original)[MIXED.len() - 1];
    let complete = parse_bytes(&original).unwrap();
    assert_eq!(complete.records.len(), MIXED.len());
    for at in 0..original.len() {
        for value in 0..=255u8 {
            if value == original[at] {
                continue;
            }
            let mut damaged = original.clone();
            damaged[at] = value;
            match parse_bytes(&damaged) {
                // Only a lost final newline reads as an interrupted write, and
                // then the whole last record is set aside, never reinterpreted.
                Ok(parsed) if at == last => {
                    assert_eq!(parsed.records.len(), MIXED.len() - 1);
                    assert_eq!(parsed.unterminated_bytes as usize, last_line);
                }
                Err(e) if at == last => panic!("a lost final newline was refused: {e}"),
                Ok(_) => panic!("mutation at {at} to {value} was accepted"),
                Err(_) => {}
            }
        }
    }
}

#[test]
fn every_bit_flip_on_disk_is_refused_without_modification() {
    let _files = files();
    let dir = Folder::new();
    let original = certificate(&MIXED);
    for at in 0..original.len() - 1 {
        let mut damaged = original.clone();
        damaged[at] ^= 1;
        fs::write(dir.file(), &damaged).unwrap();
        let before = snapshot(&dir.file());
        assert!(inspect(&dir.file()).is_err(), "byte {at}");
        assert!(open(&dir.file()).is_err(), "byte {at}");
        assert_eq!(snapshot(&dir.file()), before, "byte {at}");
    }
}

#[test]
fn complete_invalid_lines_are_never_treated_as_an_interrupted_tail() {
    let _files = files();
    let dir = Folder::new();
    let damaged = {
        let mut line = record_line("10", &data_for("10"));
        line[20] ^= 4;
        line
    };
    let cases: [(Vec<u8>, fn(&Error) -> bool); 6] = [
        (record_line("10", &Data::Local { cover: 1000 }), |e| {
            matches!(e, Error::Refuted { line: 3, .. })
        }),
        (damaged, |e| matches!(e, Error::Malformed { line: 3, defect: Defect::Checksum })),
        (b"not a record\n".to_vec(), |e| {
            matches!(e, Error::Malformed { line: 3, defect: Defect::Frame })
        }),
        (frame(b"{broken}"), |e| matches!(e, Error::Malformed { line: 3, .. })),
        (record_line("01", &data_for("01")), |e| matches!(e, Error::Overlap { line: 3, .. })),
        (record_line("", &data_for("")), |e| matches!(e, Error::OutOfOrder { line: 3, .. })),
    ];
    for (bad, expected) in cases {
        let mut bytes = certificate(&["0", "11"]);
        bytes.truncate(bytes.len() - record_line("11", &data_for("11")).len());
        bytes.extend(bad);
        bytes.extend(record_line("11", &data_for("11")));
        bytes.extend_from_slice(b"partial trailing wr");
        fs::write(dir.file(), &bytes).unwrap();
        let before = snapshot(&dir.file());
        let error = inspect(&dir.file()).unwrap_err();
        assert!(expected(&error), "{error}");
        let error = open(&dir.file()).err().unwrap();
        assert!(expected(&error), "{error}");
        assert_eq!(snapshot(&dir.file()), before);
    }
}

#[test]
fn paths_out_of_order_duplicate_overlapping_or_invalid_refuse_to_load() {
    let _files = files();
    let dir = Folder::new();
    let rooted = Header::new(&"a".repeat(64), "1", 4).unwrap();
    let cases: Vec<(Header, Vec<&str>, fn(&Error) -> bool)> = vec![
        (header(), vec!["0", "0"], |e| matches!(e, Error::OutOfOrder { .. })),
        (header(), vec!["1", "0"], |e| matches!(e, Error::OutOfOrder { .. })),
        (header(), vec!["01", "1"], |e| matches!(e, Error::OutOfOrder { .. })),
        (header(), vec!["00", "0"], |e| matches!(e, Error::OutOfOrder { .. })),
        (header(), vec!["0", "01"], |e| matches!(e, Error::Overlap { .. })),
        (header(), vec!["", "0"], |e| matches!(e, Error::Overlap { .. })),
        (header(), vec!["10", "11", "100"], |e| matches!(e, Error::Overlap { .. })),
        (header(), vec!["0", "10", "011"], |e| matches!(e, Error::Overlap { .. })),
        (rooted.clone(), vec!["0"], |e| matches!(e, Error::InvalidPath { .. })),
        (rooted.clone(), vec![""], |e| matches!(e, Error::InvalidPath { .. })),
        (rooted.clone(), vec!["10", "0110"], |e| matches!(e, Error::InvalidPath { .. })),
        (rooted.clone(), vec!["11", "10000"], |e| matches!(e, Error::InvalidPath { .. })),
    ];
    for (header, paths, expected) in cases {
        let bytes = certificate_with(&header, &paths);
        fs::write(dir.file(), &bytes).unwrap();
        let before = snapshot(&dir.file());
        let error = check(&dir.file(), &header, one(), verify).unwrap_err();
        assert!(expected(&error), "{paths:?}: {error}");
        let error = Store::<Data>::open(&dir.file(), &header, one(), verify).err().unwrap();
        assert!(expected(&error), "{paths:?}: {error}");
        assert_eq!(snapshot(&dir.file()), before);
    }
}

#[test]
fn appends_refuse_bad_paths_without_any_change() {
    let _files = files();
    let dir = Folder::new();
    let header = Header::new(&"a".repeat(64), "0", 6).unwrap();
    let (mut store, found) = Store::<Data>::open(&dir.file(), &header, one(), verify).unwrap();
    assert_eq!((found.frontier, store.frontier()), (1, vec!["0".to_owned()]));
    assert!(!found.header_complete);
    store.append("000", &data_for("000")).unwrap();
    let before = store.report();
    let cases: [(&str, fn(&Error) -> bool); 10] = [
        ("000", |e| matches!(e, Error::OutOfOrder { .. })),
        ("00", |e| matches!(e, Error::OutOfOrder { .. })),
        ("0", |e| matches!(e, Error::OutOfOrder { .. })),
        ("0000", |e| matches!(e, Error::Overlap { .. })),
        ("0001", |e| matches!(e, Error::Overlap { .. })),
        ("1000", |e| matches!(e, Error::InvalidPath { .. })),
        ("0002", |e| matches!(e, Error::InvalidPath { .. })),
        ("0x00", |e| matches!(e, Error::InvalidPath { .. })),
        ("0100000", |e| matches!(e, Error::InvalidPath { .. })),
        ("", |e| matches!(e, Error::InvalidPath { .. })),
    ];
    for (path, expected) in cases {
        let error = store.append(path, &data_for(path)).unwrap_err();
        assert!(expected(&error), "{path:?}: {error}");
        assert_eq!(store.report(), before);
    }
    let added = ["001", "0100", "010100"];
    for path in added {
        store.append(path, &data_for(path)).unwrap();
    }
    assert_eq!(store.report().records, 4);
    store.flush().unwrap();
    drop(store);
    let checked = check(&dir.file(), &header, one(), verify).unwrap();
    assert_eq!(checked.records, 4);
    let records: Vec<String> = ["000"].iter().chain(&added).map(|p| p.to_string()).collect();
    let frontier = listed::<Data>(&dir.file(), &header);
    assert_eq!(frontier.len() as u64, checked.frontier);
    assert_minimal_cover(&records, &frontier, 6, "0");
}

#[test]
fn randomized_recovery_matches_an_independent_minimal_cover_oracle() {
    let _files = files();
    let dir = Folder::new();
    let mut state = 0x721a_c823_bcdf_3579;
    for round in 0..300 {
        let root = ["", "0", "10", "011"][round % 4];
        let header = Header::new(&"b".repeat(64), root, 7).unwrap();
        let mut leaves = random_leaves(&mut state, root, 7);
        leaves.retain(|_| random(&mut state) % 3 != 0);
        let paths: Vec<&str> = leaves.iter().map(String::as_str).collect();
        let complete = certificate_with(&header, &paths);
        let cut = random(&mut state) as usize % (complete.len() + 1);
        fs::write(dir.file(), &complete[..cut]).unwrap();
        let checked = check(&dir.file(), &header, one(), verify).unwrap();
        let retained = checked.records as usize;
        let (mut store, found) =
            Store::<Data>::open(&dir.file(), &header, threads(3), verify).unwrap();
        assert_eq!(found, checked);
        let frontier = store.frontier();
        assert_eq!(frontier.len() as u64, checked.frontier);
        assert_minimal_cover(&leaves[..retained], &frontier, 7, root);
        for path in &leaves[retained..] {
            store.append(path, &data_for(path)).unwrap();
        }
        store.flush().unwrap();
        drop(store);
        assert_eq!(fs::read(dir.file()).unwrap(), complete);
        let full = check(&dir.file(), &header, one(), verify).unwrap();
        let frontier = listed::<Data>(&dir.file(), &header);
        assert_eq!(frontier.len() as u64, full.frontier);
        assert_minimal_cover(&leaves, &frontier, 7, root);
    }
}

#[test]
fn complement_examples() {
    let set = |paths: &[&str]| paths.iter().map(|p| p.to_string()).collect::<BTreeSet<_>>();
    assert_eq!(frontier("", &set(&[])), [""]);
    assert!(frontier("", &set(&[""])).is_empty());
    assert_eq!(frontier("", &set(&["00"])), ["1", "01"]);
    assert_eq!(frontier("", &set(&["10"])), ["0", "11"]);
    assert_eq!(frontier("", &set(&["0", "11"])), ["10"]);
    assert_eq!(frontier("01", &set(&["0110"])), ["010", "0111"]);
    assert_eq!(frontier("01", &set(&[])), ["01"]);
    // The traversal itself visits in preorder: lower halves first.
    let mut visited = Vec::new();
    complement("", &set(&["0110", "10"]), |path| visited.push(path));
    assert_eq!(visited, ["00", "010", "0111", "11"]);
}

#[test]
fn deep_records_recover_without_recursion() {
    let _files = files();
    let dir = Folder::new();
    let path = "0".repeat(MAX_DEPTH);
    fs::write(dir.file(), certificate(&[&path])).unwrap();
    let checked = inspect(&dir.file()).unwrap();
    assert_eq!(checked.frontier, MAX_DEPTH as u64);
    let (mut store, _) = open(&dir.file()).unwrap();
    let frontier = store.frontier();
    assert_eq!(frontier.len(), MAX_DEPTH);
    assert_eq!(frontier[0], "1");
    assert_eq!(frontier[MAX_DEPTH - 1], format!("{}1", "0".repeat(MAX_DEPTH - 1)));
    let deeper = format!("1{}", "0".repeat(MAX_DEPTH));
    assert!(matches!(
        store.append(&deeper, &data_for(&deeper)),
        Err(Error::InvalidPath { .. })
    ));
}

#[test]
fn oversized_lines_are_refused_without_modification() {
    let _files = files();
    let dir = Folder::new();
    let mut bytes = certificate(&["0"]);
    bytes.extend(std::iter::repeat(b'x').take(MAX_LINE + 1));
    fs::write(dir.file(), &bytes).unwrap();
    let before = snapshot(&dir.file());
    assert!(matches!(inspect(&dir.file()), Err(Error::TooLong { line: 3 })));
    assert!(matches!(open(&dir.file()), Err(Error::TooLong { line: 3 })));
    assert_eq!(snapshot(&dir.file()), before);
    // A line of exactly the limit is read (and refused only for its content).
    let mut bytes = certificate(&["0"]);
    bytes.extend(std::iter::repeat(b'x').take(MAX_LINE - 1));
    bytes.push(b'\n');
    fs::write(dir.file(), &bytes).unwrap();
    assert!(matches!(
        inspect(&dir.file()),
        Err(Error::Malformed { line: 3, defect: Defect::Frame })
    ));
}

// ---------------------------------------------------------------------------
// Headers
// ---------------------------------------------------------------------------

#[test]
fn headers_must_match_exactly() {
    let _files = files();
    let dir = Folder::new();
    let good = certificate(&["0"]);
    fs::write(dir.file(), &good).unwrap();
    let before = snapshot(&dir.file());
    for other in [
        Header::new(&"b".repeat(64), "", MAX_DEPTH).unwrap(),
        Header::new(&"a".repeat(64), "0", MAX_DEPTH).unwrap(),
        Header::new(&"a".repeat(64), "", 64).unwrap(),
    ] {
        assert!(matches!(
            Store::<Data>::open(&dir.file(), &other, one(), verify),
            Err(Error::WrongHeader { .. })
        ));
        assert!(matches!(
            check::<Data, _>(&dir.file(), &other, one(), verify),
            Err(Error::WrongHeader { .. })
        ));
    }
    assert_eq!(snapshot(&dir.file()), before);
    let payload = serde_json::to_string(&header()).unwrap();
    for (text, expected) in [
        (payload.replace(",", ", "), Defect::NotCanonical),
        (payload.replace("rid-certificate/1", "rid-certificate/2"), Defect::Json(String::new())),
        (payload.replace("}", r#","extra":1}"#), Defect::Json(String::new())),
        (payload.replace(r#""root":"","#, ""), Defect::Json(String::new())),
        (payload.replace("4096", "4096.0"), Defect::Json(String::new())),
        (payload.replace("4096", "\"4096\""), Defect::Json(String::new())),
    ] {
        let mut bytes = frame(text.as_bytes());
        bytes.extend(record_line("0", &data_for("0")));
        fs::write(dir.file(), &bytes).unwrap();
        let before = snapshot(&dir.file());
        let error = open(&dir.file()).err().unwrap();
        match (&error, &expected) {
            (Error::Malformed { line: 1, defect: Defect::Json(_) }, Defect::Json(_)) => {}
            (Error::Malformed { line: 1, defect }, expected) => {
                assert_eq!(defect, expected, "{text}")
            }
            _ => panic!("{text}: {error}"),
        }
        assert_eq!(snapshot(&dir.file()), before);
    }
}

#[test]
fn interrupted_headers_are_rewritten_only_when_they_are_exact_prefixes() {
    let _files = files();
    let dir = Folder::new();
    let line = header().line();
    for cut in 0..line.len() {
        fs::write(dir.file(), &line[..cut]).unwrap();
        let checked = inspect(&dir.file()).unwrap();
        assert_eq!(
            (
                checked.header_complete,
                checked.records,
                checked.frontier,
                checked.complete()
            ),
            (false, 0, 1, false)
        );
        assert_eq!(checked.unterminated_bytes as usize, cut);
        assert_eq!(fs::read(dir.file()).unwrap(), &line[..cut]);
        let (store, found) = open(&dir.file()).unwrap();
        assert_eq!(found, checked);
        assert_eq!(store.report().bytes as usize, line.len());
        drop(store);
        assert_eq!(fs::read(dir.file()).unwrap(), line);
    }
    let other = Header::new(&"b".repeat(64), "", MAX_DEPTH).unwrap().line();
    for partial in [&other[..other.len() - 1], b"unrelated".as_slice(), &line[1..40]] {
        fs::write(dir.file(), partial).unwrap();
        let before = snapshot(&dir.file());
        assert!(matches!(inspect(&dir.file()), Err(Error::WrongPartialHeader)));
        assert!(matches!(open(&dir.file()), Err(Error::WrongPartialHeader)));
        assert_eq!(snapshot(&dir.file()), before);
    }
    // A missing file is created by the writer, never by the checker.
    fs::remove_file(dir.file()).unwrap();
    assert!(matches!(inspect(&dir.file()), Err(Error::Io { .. })));
    assert!(!dir.file().exists());
    let (store, found) = open(&dir.file()).unwrap();
    assert_eq!((found.header_complete, found.frontier), (false, 1));
    assert_eq!(store.frontier(), [""]);
    drop(store);
    assert_eq!(fs::read(dir.file()).unwrap(), line);
}

// ---------------------------------------------------------------------------
// Verification
// ---------------------------------------------------------------------------

fn many_paths(depth: usize) -> Vec<String> {
    (0..1usize << depth).map(|n| format!("{n:0depth$b}")).collect()
}

#[test]
fn every_record_is_verified_exactly_once_on_any_thread_count() {
    let _files = files();
    let dir = Folder::new();
    let paths = many_paths(6);
    let refs: Vec<&str> = paths.iter().map(String::as_str).collect();
    fs::write(dir.file(), certificate(&refs)).unwrap();
    for n in 1..=5 {
        let seen = Mutex::new(Vec::new());
        let checked = check(&dir.file(), &header(), threads(n), |path, data: &Data| {
            seen.lock().unwrap().push(path.to_owned());
            verify(path, data)
        })
        .unwrap();
        assert!(checked.complete());
        let mut seen = seen.into_inner().unwrap();
        seen.sort();
        assert_eq!(seen, paths, "threads {n}");
    }
}

#[test]
fn the_earliest_false_record_is_reported_whatever_the_thread_count() {
    let _files = files();
    let dir = Folder::new();
    let paths = many_paths(5);
    let mut bytes = header().line();
    for (n, path) in paths.iter().enumerate() {
        let data = if n == 7 || n == 19 || n == 20 {
            Data::Local { cover: 999 }
        } else {
            data_for(path)
        };
        bytes.extend(record_line(path, &data));
    }
    bytes.extend_from_slice(b"interrupted");
    fs::write(dir.file(), &bytes).unwrap();
    let before = snapshot(&dir.file());
    for n in 1..=8 {
        let slow = |path: &str, data: &Data| {
            // Later records finish first, to tempt a wrong report.
            std::thread::sleep(Duration::from_micros(200 * (32 - hash(path) % 32)));
            verify(path, data)
        };
        match check(&dir.file(), &header(), threads(n), slow) {
            Err(Error::Refuted { line: 9, path, .. }) => assert_eq!(path, paths[7]),
            other => panic!("threads {n}: {:?}", other.map(|i| i.records)),
        }
        match Store::open(&dir.file(), &header(), threads(n), slow) {
            Err(Error::Refuted { line: 9, .. }) => {}
            other => panic!("threads {n}: {:?}", other.map(|s| s.1.records)),
        }
        assert_eq!(snapshot(&dir.file()), before);
    }
}

#[test]
fn verifier_panics_are_refusals() {
    let _files = files();
    let dir = Folder::new();
    let mut bytes = certificate(&["0", "10", "11"]);
    bytes.extend_from_slice(b"tail");
    fs::write(dir.file(), &bytes).unwrap();
    let before = snapshot(&dir.file());
    for n in 1..=3 {
        let result = Store::open(&dir.file(), &header(), threads(n), |path: &str, _: &Data| {
            if path == "10" {
                panic!("verifier bug");
            }
            Ok(())
        });
        match result {
            Err(Error::Refuted { line: 3, source, .. }) => {
                assert!(source.to_string().contains("verifier bug"))
            }
            other => panic!("{:?}", other.map(|s| s.1.records)),
        }
        assert_eq!(snapshot(&dir.file()), before);
    }
}

// ---------------------------------------------------------------------------
// The store
// ---------------------------------------------------------------------------

#[test]
fn a_store_extends_and_completes_a_certificate() {
    let _files = files();
    let dir = Folder::new();
    let (mut store, found) = open(&dir.file()).unwrap();
    assert_eq!((found.records, found.frontier, store.frontier()), (0, 1, vec![String::new()]));
    for path in ["0", "10"] {
        store.append(path, &data_for(path)).unwrap();
    }
    let report = store.report();
    assert_eq!(report.records, 2);
    assert_eq!(report.last.as_deref(), Some("10"));
    store.flush().unwrap();
    assert_eq!(report.bytes, fs::metadata(dir.file()).unwrap().len());
    drop(store);
    let (mut store, found) = open(&dir.file()).unwrap();
    assert_eq!((found.frontier, store.frontier()), (1, vec!["11".to_owned()]));
    store.append("11", &data_for("11")).unwrap();
    // Dropping without flush still hands the buffer to the operating system.
    drop(store);
    let checked = inspect(&dir.file()).unwrap();
    assert!(checked.complete());
    let bytes = fs::read(dir.file()).unwrap();
    assert_eq!(bytes, certificate(&["0", "10", "11"]));
    let (store, found) = open(&dir.file()).unwrap();
    assert!(found.frontier == 0 && store.frontier().is_empty() && store.header() == &header());
    drop(store);
    assert_eq!(fs::read(dir.file()).unwrap(), bytes);
}

#[test]
fn read_only_check_changes_nothing_in_any_state() {
    let _files = files();
    let dir = Folder::new();
    let mut tail = certificate(&["0", "1"]);
    tail.extend_from_slice(b"unfinished");
    let mut false_record = certificate(&["0"]);
    false_record.extend(record_line("1", &Data::Marker));
    let states: Vec<Vec<u8>> = vec![
        Vec::new(),
        header().line()[..10].to_vec(),
        certificate(&[]),
        certificate(&["0", "1"]),
        tail,
        false_record,
        b"garbage\n".to_vec(),
    ];
    for bytes in states {
        fs::write(dir.file(), &bytes).unwrap();
        let before = snapshot(&dir.file());
        std::thread::sleep(Duration::from_millis(5));
        let _ = inspect(&dir.file());
        let _ = check(&dir.file(), &header(), threads(3), |_: &str, _: &Data| {
            Err("refuse everything".into())
        });
        assert_eq!(snapshot(&dir.file()), before);
    }
}

#[test]
fn an_active_writer_excludes_writers_and_readers_in_another_process() {
    let _processes = processes();
    let dir = Folder::new();
    let (mut store, _) = open(&dir.file()).unwrap();
    store.append("0", &data_for("0")).unwrap();
    store.flush().unwrap();
    assert!(matches!(open(&dir.file()), Err(Error::Locked(_))));
    assert!(matches!(inspect(&dir.file()), Err(Error::Locked(_))));
    let status = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "search::certificate::tests::lock_child", "--nocapture"])
        .arg("--test-threads=1")
        .env("RID_CERTIFICATE_LOCK_CHILD", dir.file())
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(fs::read(dir.0.join("child")).unwrap(), b"both refused");
    drop(store);
    // Readers share their lock with each other but not with a writer.
    let reader = File::open(dir.file()).unwrap();
    FileExt::try_lock_shared(&reader).unwrap();
    assert!(inspect(&dir.file()).is_ok());
    assert!(matches!(open(&dir.file()), Err(Error::Locked(_))));
    drop(reader);
    assert!(open(&dir.file()).is_ok());
}

/// Runs only as the child of the test above.
#[test]
fn lock_child() {
    let Some(file) = std::env::var_os("RID_CERTIFICATE_LOCK_CHILD") else {
        return;
    };
    let file = Path::new(&file);
    assert!(matches!(open(file), Err(Error::Locked(_))));
    assert!(matches!(inspect(file), Err(Error::Locked(_))));
    fs::write(file.parent().unwrap().join("child"), b"both refused").unwrap();
}

#[test]
fn simultaneous_creators_have_exactly_one_owner() {
    let _files = files();
    let dir = Folder::new();
    for _ in 0..32 {
        let start = Arc::new(Barrier::new(5));
        let holding = Arc::new(Barrier::new(5));
        let winners = Arc::new(AtomicU64::new(0));
        std::thread::scope(|scope| {
            for _ in 0..4 {
                let (start, holding, winners) = (start.clone(), holding.clone(), winners.clone());
                let file = dir.file();
                scope.spawn(move || {
                    start.wait();
                    let owner = open(&file);
                    if owner.is_ok() {
                        winners.fetch_add(1, Ordering::SeqCst);
                    } else {
                        assert!(matches!(owner, Err(Error::Locked(_))));
                    }
                    holding.wait();
                    holding.wait();
                    drop(owner);
                });
            }
            start.wait();
            holding.wait();
            assert_eq!(winners.load(Ordering::SeqCst), 1);
            holding.wait();
        });
        assert_eq!(fs::read(dir.file()).unwrap(), header().line());
        fs::remove_file(dir.file()).unwrap();
    }
}

/// Measures the hot paths. Run with
/// `cargo test --release --lib search::certificate::tests::measure -- --ignored
/// --nocapture`.
#[test]
#[ignore]
fn measure() {
    let _files = files();
    let dir = Folder::new();
    let count = 1_000_000usize;
    let depth = 20;
    let paths: Vec<String> = (0..count).map(|n| format!("{n:0depth$b}")).collect();
    let started = std::time::Instant::now();
    let (mut store, _) = open(&dir.file()).unwrap();
    for path in &paths {
        store.append(path, &data_for(path)).unwrap();
    }
    store.flush().unwrap();
    drop(store);
    let append = started.elapsed();
    let bytes = fs::metadata(dir.file()).unwrap().len();
    for n in [1, 3] {
        let started = std::time::Instant::now();
        let checked = check(&dir.file(), &header(), threads(n), verify).unwrap();
        println!(
            "check with {n} threads: {:.2} s for {} records",
            started.elapsed().as_secs_f64(),
            checked.records
        );
    }
    // A verifier costing about 50 us of CPU, on the first 100000 records.
    let busy = |path: &str, data: &Data| {
        let mut digest = Sha256::digest(path.as_bytes());
        for _ in 0..800 {
            digest = Sha256::digest(digest);
        }
        std::hint::black_box(digest);
        verify(path, data)
    };
    let small = dir.0.join("small");
    let head: Vec<&str> = paths[..100_000].iter().map(String::as_str).collect();
    fs::write(&small, certificate(&head)).unwrap();
    for n in [1, 3] {
        let started = std::time::Instant::now();
        check(&small, &header(), threads(n), busy).unwrap();
        println!(
            "check with a busy verifier on {n} threads: {:.2} s for 100000 records",
            started.elapsed().as_secs_f64()
        );
    }
    let started = std::time::Instant::now();
    let (_store, found) = Store::<Data>::open(&dir.file(), &header(), threads(3), verify).unwrap();
    println!(
        "open with 3 threads: {:.2} s, frontier {}",
        started.elapsed().as_secs_f64(),
        found.frontier
    );
    println!(
        "append+flush: {:.2} s for {count} records ({:.2} us each), {bytes} bytes",
        append.as_secs_f64(),
        append.as_secs_f64() * 1e6 / count as f64
    );
}

mod adversarial;
