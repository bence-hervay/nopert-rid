//! A deterministic random source for the tests (SplitMix64), so that every
//! randomised test is reproducible from its seed.
use rid::arithmetic::exact::Q;
use num_bigint::BigInt;

pub struct Random(u64);

impl Random {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// An integer in `[lo, hi]`.
    pub fn integer(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }

    /// A rational `lo + (hi - lo)·k/2^bits` with a random `k ∈ [0, 2^bits]`.
    pub fn rational(&mut self, lo: &Q, hi: &Q, bits: u32) -> Q {
        let k = self.next() % ((1u64 << bits) + 1);
        lo + &((hi - lo) * Q::new(BigInt::from(k), BigInt::from(1u64 << bits)))
    }
}
