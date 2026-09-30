//! The minimal discrete surface, in the two shapes that matter.
//!
//! Every function here exists twice: a `_frozen`/`_sq` form obeying the comparison-only rule, and
//! a `_trans`/`_sqrt` twin that deliberately breaks it. The twin is the NEGATIVE CONTROL. A test
//! that cannot fail is indistinguishable from a test that passes -- if the violating arm also
//! agrees across backends, the harness is not sensitive enough and the clean arm's agreement means
//! nothing.

pub mod gpu;

/// Substep-bucket parameters, from the spike brief's `N_sub` definition.
pub const R_SUB: f64 = 0.05;
pub const GAMMA: f64 = 1.5;
pub const N_MAX: u32 = 64;

/// The frozen threshold table. `THR[n]` is the squared distance at which the bucket count becomes
/// `n + 1`, i.e. `(r_sub / (n+1)^(2/3))^2`, computed **once in f64 libm** and rounded to f32.
///
/// This is the whole mechanism: the `pow` lives here, at table-generation time, and every backend
/// then compares against the identical bit pattern. Nothing at runtime computes a transcendental.
pub fn thresholds() -> [f32; N_MAX as usize] {
    let mut t = [0f32; N_MAX as usize];
    for n in 0..N_MAX as usize {
        let d = R_SUB / ((n + 1) as f64).powf(2.0 / 3.0);
        t[n] = (d * d) as f32;
    }
    t
}

/// Comparison-only bucket: `d2` against frozen constants, descending. No sqrt, no div, no pow.
/// Returns `N_sub` in `1..=N_MAX`.
#[inline]
pub fn bucket_frozen(d2: f32, thr: &[f32; N_MAX as usize]) -> u32 {
    // Walk from the largest threshold down; the first one d2 is NOT below fixes the bucket.
    // Written as a loop for legibility; the contract's generated form is an unrolled comparison
    // tree, which lowers to the same comparisons without the runtime index.
    let mut n: u32 = 1;
    let mut i: usize = 0;
    while i < N_MAX as usize {
        if d2 < thr[i] {
            n = (i + 1) as u32 + 1;
        }
        i += 1;
    }
    if n > N_MAX { N_MAX } else { n }
}

/// NEGATIVE CONTROL: the spike-brief formula, evaluated at runtime.
/// `N_sub = min(N_max, max(1, ceil((r_sub / r_min)^gamma)))` -- a sqrt, a divide and a pow, all
/// feeding a `ceil` whose integer output is a branch decision. This is the arm that forked 5/272.
#[inline]
pub fn bucket_trans(d2: f32) -> u32 {
    let r_min = d2.sqrt();
    let q = (R_SUB as f32) / r_min;
    let p = q.powf(GAMMA as f32);
    let c = p.ceil();
    let c = if c < 1.0 { 1.0 } else { c };
    let c = if c > N_MAX as f32 { N_MAX as f32 } else { c };
    c as u32
}

/// Comparison-only collision: squared distance against a squared constant.
#[inline]
pub fn collision_sq(d2: f32, r_coll_sq: f32) -> bool { d2 < r_coll_sq }

/// NEGATIVE CONTROL: prin-rs's actual shape -- `d < r_coll` with `d` from a sqrt
/// (`src/vec2.rs:42`, `src/integrate/heggie/driver.rs:431`).
#[inline]
pub fn collision_sqrt(d2: f32, r_coll: f32) -> bool { d2.sqrt() < r_coll }

/// The crate's only bit-packing, copied verbatim from `src/outcome.rs:157` -- not imported, so the
/// spike stays independent of prin-rs.
#[inline]
pub fn pack(state: u32, detail: u32) -> u32 { ((state & 7) << 2) | (detail & 3) }

/// The six-variant code table (`src/outcome.rs:104-112`). Returns 0xFF for an unknown code, which
/// is the spike's stand-in for `None`.
#[inline]
pub fn from_bits(code: u32) -> u32 {
    let s = (code >> 2) & 7;
    if s <= 5 { s } else { 0xFF }
}

/// One input's full discrete output. Every field is an integer: this is the Tier-B surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Discrete {
    pub bucket_frozen: u32,
    pub bucket_trans: u32,
    pub coll_sq: u32,
    pub coll_sqrt: u32,
    pub packed: u32,
    pub roundtrip: u32,
    pub pairsum_flag: u32,
    pub pairsum_brk: u32,
    /// The same packing fed by the R1-violating bucket. Not a control in itself: the PROPAGATION
    /// arm, showing that a descriptor is only as deterministic as the branches feeding it.
    pub packed_ctl: u32,
    pub roundtrip_ctl: u32,
}

/// The triangular pair loop, flag-tested, zero `break` -- the contract's canonical wrapper.
#[inline]
pub fn pairsum_flag(seed: u32) -> u32 {
    let mut acc = 0u32;
    let mut i = 0u32;
    let mut di = false;
    while !di {
        let mut j = i + 1;
        let mut dj = false;
        while !dj {
            acc = acc.wrapping_add(seed + i * 16 + j);
            j += 1;
            dj = j >= 3;
        }
        i += 1;
        di = i >= 2;
    }
    acc
}

/// NEGATIVE CONTROL: counted `for` with mid-body `break`, two levels.
#[inline]
pub fn pairsum_brk(seed: u32) -> u32 {
    let mut acc = 0u32;
    for i in 0u32..3 {
        if i >= 2 { break; }
        for j in 0u32..3 {
            if j <= i { continue; }
            if j >= 3 { break; }
            acc = acc.wrapping_add(seed + i * 16 + j);
        }
    }
    acc
}

pub fn eval_cpu(d2: f32, thr: &[f32; N_MAX as usize], r_coll: f32) -> Discrete {
    let bf = bucket_frozen(d2, thr);
    let bt = bucket_trans(d2);
    // The packed descriptor is fed by a branch decision, so a forked bucket shows up here too --
    // which is the point: the descriptor is only as deterministic as the branches feeding it.
    let state = bf % 6;
    let packed = pack(state, bf & 3);
    let packed_ctl = pack(state, bt & 3);
    Discrete {
        bucket_frozen: bf,
        bucket_trans: bt,
        coll_sq: collision_sq(d2, r_coll * r_coll) as u32,
        coll_sqrt: collision_sqrt(d2, r_coll) as u32,
        packed,
        roundtrip: from_bits(packed),
        pairsum_flag: pairsum_flag(bf),
        pairsum_brk: pairsum_brk(bf),
        packed_ctl,
        roundtrip_ctl: from_bits(packed_ctl),
    }
}

/// The collision threshold the sweep must straddle, or `coll_sqrt` is a control with no subject.
pub const R_COLL: f32 = 0.02;

/// The boundary sweep: the states where `ceil()` flips. A fork hides exactly here and nowhere else,
/// so a sweep of ordinary states would report a clean pass while learning nothing.
pub fn boundary_sweep() -> Vec<f32> {
    let mut v = Vec::new();
    for n in 1..=N_MAX {
        let d_exact = (R_SUB / (n as f64).powf(2.0 / 3.0)) as f32;
        // Step +/- k ulp around the exact boundary.
        for k in -4i32..=4 {
            let bits = d_exact.to_bits() as i32;
            let d = f32::from_bits((bits + k) as u32);
            v.push(d * d);
        }
    }
    // And the COLLISION boundary. Without this the sweep never approaches `r_coll` and `coll_sqrt`
    // is a control with nothing to decide -- it read 0 forks for that reason and not because sqrt
    // agrees. WGSL specifies sqrt to 1 ULP, not correctly-rounded, so it CAN fork; the first cut of
    // this sweep simply never asked it to.
    for k in -64i32..=64 {
        let bits = R_COLL.to_bits() as i32;
        let d = f32::from_bits((bits + k) as u32);
        v.push(d * d);
    }
    v
}

/// The correct triangular pair sum, computed in closed form: pairs (0,1), (0,2), (1,2).
/// A loop that agrees across backends while being WRONG on both is not a passing control --
/// it is the miscompile happening identically everywhere.
#[inline]
pub fn pairsum_expected(seed: u32) -> u32 { seed.wrapping_mul(3).wrapping_add(21) }
