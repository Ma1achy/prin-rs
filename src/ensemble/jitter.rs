//! Per-pixel ensemble construction.
//!
//! BRIEF §7: seed the jitter per pixel from `(i, j, seed)`, **never from a global RNG**, so
//! any pixel is reproducible in isolation. This is incompatible with the reference, which
//! draws one PCG64 stream across the whole array — which is why the Step 4 cross-check runs
//! on nominal copies only, where no RNG participates on either side.
//!
//! Initial conditions are always built in `f64` and cast down. Generating them separately
//! per precision would make an IC difference indistinguishable from a genuine f32 arithmetic
//! effect, which is the decomposition the f32 question depends on.

use crate::grid::Slice;
use crate::physics::Cart;
use crate::rng::SplitMix64;
use crate::{Real, Vec2};

/// Mix `(i, j, seed)` into a stream seed. Written out rather than delegated, so a dependency
/// bump cannot move the initial conditions underneath us.
pub fn pixel_seed(i: usize, j: usize, seed: u64) -> u64 {
    let mut z = seed;
    z ^= (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = z.rotate_left(31);
    z ^= (j as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    z = z.rotate_left(17);
    z ^= 0x1656_67B1_9E37_79F9;
    SplitMix64::new(z).next_u64()
}

/// The `n_extra + 1` copies of one pixel.
///
/// **Copy 0 is always the un-jittered nominal.** Load-bearing, and asserted in a test: it is
/// what makes a nominal-only cross-check possible, and it is the copy whose reference-body
/// choices the shared-reference policy propagates.
///
/// Jitter is `jitter_frac * cell width`, **per axis**. The reference computes only `hx` and
/// uses it for both axes — latent on square grids, wrong on any other. Scaling with cell
/// width is required: a fixed perturbation would make measured spreads drift with resolution
/// for a purely trivial reason (BRIEF §3).
pub fn copies<T: Real>(
    slice: &Slice,
    idx: usize,
    n_extra: usize,
    jitter_frac: f64,
    seed: u64,
) -> Vec<Cart<T>> {
    let (hx, hy) = slice.cell_widths();
    let jx = jitter_frac * hx;
    let jy = jitter_frac * hy;
    let (i, j) = (idx % slice.nx, idx / slice.nx);
    let mut rng = SplitMix64::new(pixel_seed(i, j, seed));

    let base = slice.nominal::<f64>(idx);
    let mut out = Vec::with_capacity(n_extra + 1);
    out.push(base.cast::<T>());
    for _ in 0..n_extra {
        let mut c = base;
        c.r[slice.body] += Vec2::new(rng.range(-jx, jx), rng.range(-jy, jy));
        out.push(c.cast::<T>());
    }
    out
}
