//! SplitMix64.
//!
//! Used two ways, and the reasons differ:
//!
//! - In the cross-check, both languages implement it identically so the two sides see
//!   bit-identical random configurations. Each side's native RNG would make the comparison
//!   meaningless.
//! - In the ensemble, it mixes `(i, j, seed)` into a per-pixel stream, per BRIEF §7 —
//!   "never from a global RNG", so any pixel is reproducible in isolation.
//!
//! The constants are written out rather than imported, so a dependency bump cannot move the
//! initial conditions underneath us.

pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)` from the top 53 bits. Exact, and reproduced identically in
    /// `tools/xcheck/cases.py`.
    #[inline]
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    #[inline]
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}
