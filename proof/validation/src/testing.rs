//! Helpers shared by the tests: a seeded generator of exact numbers and
//! scratch directories.
use rid::arithmetic::exact::{Q, QSqrt5};
use num_bigint::BigInt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// SplitMix64: deterministic test data for a seed.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n`.
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// A rational in `[lo, hi]` with denominator `2^bits`.
    pub fn rational(&mut self, lo: &Q, hi: &Q, bits: u32) -> Q {
        let scale = Q::from_integer(BigInt::from(1) << bits);
        let fraction = Q::from_integer(BigInt::from(self.below((1 << bits) + 1))) / scale;
        lo + &((hi - lo) * fraction)
    }

    /// A rational with numerator and denominator below `2^bits`, of either sign.
    pub fn small(&mut self, bits: u32) -> Q {
        let n = self.below(1 << bits) as i64 - (1 << (bits - 1));
        let d = self.below(1 << bits) as i64 + 1;
        Q::new(BigInt::from(n), BigInt::from(d))
    }

    /// `a + b√5` with small random `a` and `b`.
    pub fn number(&mut self, bits: u32) -> QSqrt5 {
        QSqrt5::new(self.small(bits), self.small(bits))
    }
}

/// A fresh empty directory, removed when dropped.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(name: &str) -> Self {
        static COUNT: AtomicUsize = AtomicUsize::new(0);
        let n = COUNT.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rid-validation-{}-{name}-{n}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    pub fn join(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The shipped catalogue files.
pub fn shipped(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("catalogue").join(name)
}
