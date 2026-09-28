//! Floating point in the crate's production code only where it proposes:
//! the float view of the data (`src/elimination/float/`), the cover
//! generator (`src/elimination/generate/`) and the components' proposals
//! (`src/components/proposal/`). Every other Rust file below `src/`, except
//! test code (files named `*tests.rs` and anything below a `tests`
//! directory, the same rule as the policy fingerprint's), is scanned for
//! floating-point types, conversions and literals, with comments, strings
//! and character literals skipped. The scanner is itself tested on
//! adversarial snippets.
use std::fs;
use std::path::{Path, PathBuf};

/// The code of `text` with comments, string literals (raw ones included)
/// and character literals replaced by spaces, so that only code is scanned.
fn code(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    let at = |i: usize| chars.get(i).copied().unwrap_or('\0');
    while i < chars.len() {
        let c = chars[i];
        let identifier_before = i > 0 && (chars[i - 1].is_alphanumeric() || chars[i - 1] == '_');
        if c == '/' && at(i + 1) == '/' {
            while i < chars.len() && chars[i] != '\n' {
                i += 1;
            }
        } else if c == '/' && at(i + 1) == '*' {
            let mut depth = 0;
            loop {
                if at(i) == '/' && at(i + 1) == '*' {
                    depth += 1;
                    i += 2;
                } else if at(i) == '*' && at(i + 1) == '/' {
                    depth -= 1;
                    i += 2;
                    if depth == 0 {
                        break;
                    }
                } else if i >= chars.len() {
                    break;
                } else {
                    i += 1;
                }
            }
            out.push(' ');
        } else if !identifier_before && (c == 'r' || (c == 'b' && at(i + 1) == 'r')) && {
            let start = if c == 'b' { i + 2 } else { i + 1 };
            let hashes = chars[start..].iter().take_while(|&&h| h == '#').count();
            at(start + hashes) == '"'
        } {
            // A raw string: r"…", r#"…"#, br"…".
            let start = if c == 'b' { i + 2 } else { i + 1 };
            let hashes = chars[start..].iter().take_while(|&&h| h == '#').count();
            i = start + hashes + 1;
            while i < chars.len() && !(chars[i] == '"' && chars[i + 1..].iter().take(hashes).filter(|&&h| h == '#').count() == hashes) {
                i += 1;
            }
            i += 1 + hashes;
            out.push(' ');
        } else if c == '"' {
            i += 1;
            while i < chars.len() && chars[i] != '"' {
                i += if chars[i] == '\\' { 2 } else { 1 };
            }
            i += 1;
            out.push(' ');
        } else if c == '\'' && (at(i + 1) == '\\' || at(i + 2) == '\'') {
            // A character literal ('x' or an escape); a lifetime has no
            // closing quote after one character.
            i += 1;
            while i < chars.len() && chars[i] != '\'' {
                i += if chars[i] == '\\' { 2 } else { 1 };
            }
            i += 1;
            out.push(' ');
        } else {
            out.push(c);
            i += 1;
        }
    }
    out
}

/// The floating-point tokens of `code` (already stripped by [`code`]):
/// identifiers with a part `f32` or `f64` (`f64`, `to_f64`, `as_secs_f64`)
/// and number literals with a fractional part, an exponent or a float
/// suffix. Tuple indices (`x.0.1`), ranges (`0..3`) and methods on integers
/// (`1.max(2)`) are not literals of a float.
fn floating(code: &str) -> Vec<String> {
    let chars: Vec<char> = code.chars().collect();
    let mut found = Vec::new();
    let mut i = 0;
    let word = |c: char| c.is_alphanumeric() || c == '_';
    while i < chars.len() {
        let c = chars[i];
        let before = if i > 0 { chars[i - 1] } else { ' ' };
        if c.is_alphabetic() || c == '_' {
            let start = i;
            while i < chars.len() && word(chars[i]) {
                i += 1;
            }
            let identifier: String = chars[start..i].iter().collect();
            if identifier.split('_').any(|part| part == "f32" || part == "f64") {
                found.push(identifier);
            }
        } else if c.is_ascii_digit() && !word(before) && before != '.' {
            let start = i;
            while i < chars.len() && word(chars[i]) {
                i += 1;
            }
            let mut literal: String = chars[start..i].iter().collect();
            let radix = literal.starts_with("0x") || literal.starts_with("0o") || literal.starts_with("0b");
            let mut float = !radix && (literal.ends_with("f32") || literal.ends_with("f64"));
            // An exponent: digits, then e or E, then an optional sign and a digit.
            if !radix && literal.contains(['e', 'E']) {
                let (mantissa, _) = literal.split_once(['e', 'E']).unwrap();
                float |= mantissa.chars().all(|d| d.is_ascii_digit() || d == '_');
            }
            if !radix && chars.get(i) == Some(&'.') {
                let next = chars.get(i + 1).copied().unwrap_or(' ');
                if next.is_ascii_digit() || !(next == '.' || next.is_alphabetic() || next == '_') {
                    float = true;
                    literal.push('.');
                    i += 1;
                    while i < chars.len() && (word(chars[i]) || ((chars[i] == '+' || chars[i] == '-') && matches!(chars[i - 1], 'e' | 'E'))) {
                        literal.push(chars[i]);
                        i += 1;
                    }
                }
            }
            if float {
                found.push(literal);
            }
        } else {
            i += 1;
        }
    }
    found
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

#[test]
fn production_code_has_no_floating_point_outside_the_proposals() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut files = Vec::new();
    rust_files(&root.join("src"), &mut files);
    files.sort();
    let mut scanned = 0;
    let mut offending = Vec::new();
    for file in &files {
        let relative = file.strip_prefix(root).unwrap();
        let test_code = relative.components().any(|c| c.as_os_str() == "tests")
            || relative.file_name().unwrap().to_str().unwrap().ends_with("tests.rs");
        let proposes = ["src/elimination/float", "src/elimination/generate", "src/components/proposal"]
            .iter()
            .any(|allowed| relative.starts_with(allowed));
        if test_code || proposes {
            continue;
        }
        scanned += 1;
        let tokens = floating(&code(&fs::read_to_string(file).unwrap()));
        if !tokens.is_empty() {
            offending.push(format!("{}: {tokens:?}", relative.display()));
        }
    }
    assert!(scanned >= 30, "only {scanned} files scanned");
    assert!(offending.is_empty(), "floating point in production code:\n{}", offending.join("\n"));
}

#[test]
fn the_scanner_finds_floating_point_and_nothing_else() {
    let flagged = [
        ("let x = 1.0;", "1.0"),
        ("let x = 0.5 * y;", "0.5"),
        ("const T: f64 = 1e-9;", "f64"),
        ("let x = 2E5;", "2E5"),
        ("let x = 1e9 as u64;", "1e9"),
        ("let x = 1_000.25;", "1_000.25"),
        ("let x = 1f64;", "1f64"),
        ("let x = 3f32;", "3f32"),
        ("let x = y as f64;", "f64"),
        ("let x = q.to_f64();", "to_f64"),
        ("let s = t.elapsed().as_secs_f64();", "as_secs_f64"),
        ("fn f(x: f32) {}", "f32"),
        ("let v: Vec<f64> = vec![];", "f64"),
        ("let x = 1.;", "1."),
        ("let x = 2.5e-3;", "2.5e-3"),
        ("x.min(1.0)", "1.0"),
        ("let p = std::f64::consts::PI;", "f64"),
    ];
    for (snippet, token) in flagged {
        let found = floating(&code(snippet));
        assert!(found.iter().any(|f| f == token), "{snippet}: {found:?}");
    }
    let clean = [
        "for i in 0..3 {}",
        "let y = x.0.1;",
        "let m = 1.max(2);",
        "let n = 10u64 + 5i64 as u64;",
        "// a comment with 1.0 and f64\nlet x = 1;",
        "/* block 2.5 /* nested f32 */ still */ let x = 2;",
        r#"let s = "1.0 and f64 in a string";"#,
        r####"let s = r#"raw 3.5 "quoted" f64"#;"####,
        "let h = 0x1e9;",
        "let b = 0b1010;",
        "let c = 'e'; let d = '.'; let e = '\\'';",
        "fn f<'a>(x: &'a str) -> &'a str { x }",
        "let byte = b'e';",
        "let big = 1_000_000;",
        "let identifier_e5 = 1;",
        "let f64ish = 2;",
        "let x = a.e1;",
        "let r = 1..=9;",
    ];
    for snippet in clean {
        let found = floating(&code(snippet));
        assert!(found.is_empty(), "{snippet}: {found:?}");
    }
}
