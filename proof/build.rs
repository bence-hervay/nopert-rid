use sha2::{Digest, Sha256};
use std::{
    env, fs,
    path::{Path, PathBuf},
};

fn files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).expect("read source directory") {
        let path = entry.expect("source entry").path();
        assert!(!path.is_symlink(), "source symlinks are not allowed");
        if path.is_dir() {
            files(&path, out);
        } else if path.is_file() {
            out.push(path);
        }
    }
}

fn main() {
    let root = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let mut input = vec![
        root.join("Cargo.toml"),
        root.join("Cargo.lock"),
        root.join("build.rs"),
    ];
    files(&root.join("src"), &mut input);
    files(&root.join("data"), &mut input);
    // Directory watches also notice additions, not only changes to known files.
    // The dependencies are fixed by Cargo.lock, which pins their exact
    // versions and checksums.
    for directory in ["src", "data"] {
        println!("cargo:rerun-if-changed={directory}");
    }
    // Of the sources, only library and binary code counts: documentation and
    // test code (files named `*tests.rs` and anything below a `tests`
    // directory) cannot change a decision.
    input.retain(|p| {
        let relative = p.strip_prefix(&root).unwrap();
        let test_code = relative.components().any(|c| c.as_os_str() == "tests")
            || p.file_name().unwrap().to_str().unwrap().ends_with("tests.rs");
        !relative.starts_with("src") || (p.extension().is_some_and(|e| e == "rs") && !test_code)
    });
    input.sort();
    let mut hash = Sha256::new();
    hash.update(b"rid-exact-certificate-policy\0");
    for path in input {
        let name = path.strip_prefix(&root).unwrap().to_str().unwrap();
        println!("cargo:rerun-if-changed={name}");
        let bytes = fs::read(&path).expect("read policy source");
        hash.update((name.len() as u64).to_be_bytes());
        hash.update(name.as_bytes());
        hash.update((bytes.len() as u64).to_be_bytes());
        hash.update(bytes);
    }
    println!("cargo:rustc-env=RID_POLICY={:x}", hash.finalize());
}
