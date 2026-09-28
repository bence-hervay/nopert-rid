//! Tests of the cover file format: round trips of every file of the
//! crate, a hand-written file, and refusals of every other spelling.
use super::*;
use crate::elimination::proof::catalogue::{exotic_data, local_data, EXOTIC, LOCAL_COUNT};

const SMALL: &str = concat!(
    r#"{"format":"rid-cover/1","name":"7","scope":"all","#,
    r#""parameters":{"coordinates":"configuration","#,
    r#""centre":[["0","0"],["0","0"],["0","0"],["0","0"],["0","0"]],"#,
    r#""map":[[["1","0"],["0","0"],["0","0"],["0","0"],["0","0"]],[["0","0"],["1","0"],["0","0"],["0","0"],["0","0"]],[["0","0"],["0","0"],["1","0"],["0","0"],["0","0"]],[["0","0"],["0","0"],["0","0"],["1","0"],["0","0"]],[["0","0"],["0","0"],["0","0"],["0","0"],["1","0"]]],"#,
    r#""shape":{"tube":{"base":[["0","1/6"],["1/10","1/5"]],"offset":[[-1,1],[-1,1],[-1,1]],"radius":"1/80"}},"beyond":null},"#,
    r#""witnesses":[{"edge":[3,5],"vertex":3},{"direction":[["-1","0"],["2","1"]],"contact":22,"vertex":22}],"#,
    r#""zooms":[{"name":"T o0-","tree":"3..","leaves":[{"witness":{"index":0,"factor":[0,0,1,0,0]}},{"witness":{"index":1,"factor":[0,0,1,0,0]}}]}]}"#,
    "\n"
);

#[test]
fn every_cover_file_round_trips() {
    let files = EXOTIC.iter().map(|n| exotic_data(n).unwrap()).chain((0..LOCAL_COUNT).map(|n| local_data(n).unwrap()));
    for bytes in files {
        let file = CoverFile::parse(bytes).unwrap();
        assert_eq!(file.to_bytes(), bytes);
        assert_eq!(file.format, Format::Version1);
    }
}

#[test]
fn a_small_file_parses_to_its_values() {
    let file = CoverFile::parse(SMALL.as_bytes()).unwrap();
    assert_eq!(file.name, "7");
    assert_eq!(file.scope, Scope::All);
    assert_eq!(file.witnesses.len(), 2);
    assert_eq!(file.parameters.coordinates, crate::elimination::zoom::Coordinates::Configuration);
    assert_eq!(file.parameters.beyond, None);
    let zoom = &file.zooms[0];
    assert_eq!((zoom.name.as_str(), zoom.tree.as_str()), ("T o0-", "3.."));
    assert_eq!(zoom.leaves[1], Leaf::Witness { index: 1, factor: [0, 0, 1, 0, 0] });
    assert_eq!(file.to_bytes(), SMALL.as_bytes());
    let delegated = SMALL.replacen(r#"{"witness":{"index":1,"factor":[0,0,1,0,0]}}"#, r#""delegated""#, 1);
    assert_eq!(CoverFile::parse(delegated.as_bytes()).unwrap().zooms[0].leaves[1], Leaf::Delegated);
    let arc = SMALL.replacen(r#""coordinates":"configuration""#, r#""coordinates":{"arc-plane":"minus"}"#, 1);
    assert_eq!(
        CoverFile::parse(arc.as_bytes()).unwrap().parameters.coordinates,
        crate::elimination::zoom::Coordinates::ArcPlane(crate::elimination::zoom::Sign::Minus)
    );
}

#[test]
fn other_spellings_and_malformed_files_are_refused() {
    let noncanonical = [
        SMALL.replacen("{", "{ ", 1),
        SMALL.trim_end().to_string(),
        format!("{SMALL}\n"),
        SMALL.replacen("\n", "\r\n", 1),
        SMALL.replacen(r#""format":"rid-cover/1","name":"7""#, r#""name":"7","format":"rid-cover/1""#, 1),
        SMALL.replacen(r#""vertex":3},{"direction""#, r#""vertex":3 },{"direction""#, 1),
    ];
    for text in &noncanonical {
        assert!(matches!(CoverFile::parse(text.as_bytes()), Err(FormatError::NonCanonical)), "{text}");
    }
    let malformed = [
        SMALL.replacen("rid-cover/1", "rid-cover/0", 1),
        SMALL.replacen(r#""scope":"all""#, r#""scope":"everywhere""#, 1),
        SMALL.replacen(r#""scope":"all","#, "", 1),
        SMALL.replacen(r#""scope":"all","#, r#""scope":"all","scope":"all","#, 1),
        SMALL.replacen(r#""beyond":null"#, r#""beyond":null,"window":null"#, 1),
        SMALL.replacen(r#","beyond":null"#, "", 1),
        SMALL.replacen(r#""radius":"1/80""#, r#""radius":"2/160""#, 1),
        SMALL.replacen(r#""radius":"1/80""#, r#""radius":"+1/80""#, 1),
        SMALL.replacen(r#""radius":"1/80""#, r#""radius":0.0125"#, 1),
        SMALL.replacen(r#"["1","0"]"#, r#"["1"]"#, 1),
        SMALL.replacen(r#""coordinates":"configuration""#, r#""coordinates":"physical""#, 1),
        SMALL.replacen(r#""coordinates":"configuration""#, r#""coordinates":{"arc-plane":"zero"}"#, 1),
        SMALL.replacen(r#"{"witness":{"index":1,"factor":[0,0,1,0,0]}}"#, r#""outside""#, 1),
        SMALL.replacen(r#""factor":[0,0,1,0,0]"#, r#""factor":[0,0,1,0]"#, 1),
        SMALL.replacen(r#""factor":[0,0,1,0,0]"#, r#""factor":[0,0,256,0,0]"#, 1),
        SMALL.replacen(r#""index":0,"#, r#""index":-1,"#, 1),
        SMALL.replacen(r#"{"edge":[3,5],"vertex":3}"#, r#"{"edge":[3,4],"vertex":3}"#, 1),
        SMALL.replacen(r#""offset":[[-1,1],"#, r#""offset":[[-1.0,1],"#, 1),
        SMALL.replacen(r#""tree":"3..""#, r#""tree":["3",".","."]"#, 1),
        SMALL[..SMALL.len() / 3].to_string(),
        "null\n".to_string(),
    ];
    for text in &malformed {
        assert!(matches!(CoverFile::parse(text.as_bytes()), Err(FormatError::Json(_))), "{text}");
    }
}
