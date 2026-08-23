//! Per-pixel ensemble construction.
//!
//! BRIEF §7 required per-pixel seeding from `(i, j, seed)` rather than a global stream, so any
//! pixel is reproducible in isolation. The spec's actual scheme is stronger and simpler: a
//! **fixed** offset prefix, which is reproducible in isolation because it does not depend on
//! the pixel at all. Per-pixel seeding survives only as [`Scheme::Pcg`].
//!
//! Either way the Step 4 cross-check runs on nominal copies only, where no offset scheme
//! participates on either side.
//!
//! Initial conditions are always built in `f64` and cast down. Generating them separately
//! per precision would make an IC difference indistinguishable from a genuine f32 arithmetic
//! effect, which is the decomposition the f32 question depends on.
//!
//! # Two schemes
//!
//! The spec calls for a **fixed low-discrepancy Halton (2,3) prefix indexed by `copy_index`**.
//! Two properties, both load-bearing and neither provided by the reference's global PCG64
//! stream:
//!
//! - **Quasi-random, not pseudo-random.** Halton fills the footprint evenly by construction.
//!   Independent uniform draws clump and leave gaps, which at `E+1 = 4` or `8` is most of the
//!   variance in a spread estimate.
//! - **Fixed, not per-pixel.** Copy `k` sits at the same offset in *every* footprint at *every*
//!   refinement level, so a parent ensemble and its children share the perturbation pattern.
//!   That is common random numbers by construction, and the sampling noise largely cancels in
//!   the parent/child ratio that the refinement exponent is computed from.
//!
//! [`Scheme::Pcg`] keeps the reference's behaviour available. Prior results were produced with
//! it and reproducing them needs it.

use crate::decode::{self, Path};
use crate::grid::Slice;
use crate::physics::Cart;
use crate::rng::SplitMix64;
use crate::Real;

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

/// How copy offsets are chosen.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Scheme {
    /// Fixed Halton (2,3) prefix indexed by copy index. The spec, and the default.
    #[default]
    Halton,
    /// The reference's per-pixel PCG64 stream. Kept so prior results reproduce.
    Pcg,
}

/// Radical inverse of `n` in `base` — the `b`-ary digits of `n`, reflected about the point.
///
/// This is the whole of Halton: `phi_2(n)` and `phi_3(n)` are each equidistributed in `[0,1)`
/// and, being built on coprime bases, jointly fill the square without the lattice alignment
/// that a single base would give.
pub fn radical_inverse(mut n: usize, base: usize) -> f64 {
    let inv = 1.0 / base as f64;
    let mut f = inv;
    let mut r = 0.0;
    while n > 0 {
        r += (n % base) as f64 * f;
        n /= base;
        f *= inv;
    }
    r
}

/// The `k`-th Halton (2,3) point, mapped to `[-1, 1)^2`.
///
/// Indexed from 1, because `radical_inverse(0, _)` is 0 and would place a copy exactly on the
/// nominal — the one offset the ensemble already has.
pub fn halton_offset(k: usize) -> (f64, f64) {
    let n = k + 1;
    (
        2.0 * radical_inverse(n, 2) - 1.0,
        2.0 * radical_inverse(n, 3) - 1.0,
    )
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
    copies_with(slice, idx, n_extra, jitter_frac, seed, Scheme::default())
}

/// As [`copies`], with an explicit scheme.
pub fn copies_with<T: Real>(
    slice: &Slice,
    idx: usize,
    n_extra: usize,
    jitter_frac: f64,
    seed: u64,
    scheme: Scheme,
) -> Vec<Cart<T>> {
    copies_with_path(slice, idx, n_extra, jitter_frac, seed, scheme, Path::DirectF64)
}

/// The same, with the arithmetic that forms each copy's initial condition named explicitly.
///
/// **The jitter is in CHART coordinates, not in one body's Cartesian position.** The earlier
/// form added the offset to `c.r[slice.body]` directly, which is the same thing only because
/// `Chart::BodyPlane`'s decode writes `(u, v)` straight into `r[body]`. On an oblique or shape
/// chart it would perturb a body's position rather than a sub-cell sample of the chart, and the
/// ensemble would no longer be sampling the cell it is supposed to be sampling. For
/// `BodyPlane` the two are **bitwise identical** — `x + du*jx` either way — which is what
/// `tests/seeding_golden.rs` holds.
pub fn copies_with_path<T: Real>(
    slice: &Slice,
    idx: usize,
    n_extra: usize,
    jitter_frac: f64,
    seed: u64,
    scheme: Scheme,
    path: Path,
) -> Vec<Cart<T>> {
    let (hx, hy) = slice.cell_widths();
    let jx = jitter_frac * hx;
    let jy = jitter_frac * hy;
    let (x, y) = slice.decode_pos(idx);

    // Only built when a non-default path asks for it: five extra f64 decodes per footprint.
    let lin = (path != Path::DirectF64)
        .then(|| decode::linearise(&slice.chart, slice.body, slice.cx, slice.cy, slice.half));
    let make = |u: f64, v: f64| -> Cart<f64> {
        match &lin {
            None => slice.decode_state(u, v),
            Some(l) => decode::sample(
                path, &slice.chart, slice.body, slice.cx, slice.cy, slice.half,
                (u - slice.cx) / slice.half, (v - slice.cy) / slice.half, l,
            ),
        }
    };

    let mut out = Vec::with_capacity(n_extra + 1);
    out.push(make(x, y).cast::<T>());

    match scheme {
        // No RNG, no per-pixel seed, no ordering. Copy k is at the same place in every
        // footprint at every level, which is what makes parent and child ensembles share a
        // perturbation pattern.
        Scheme::Halton => {
            for k in 0..n_extra {
                let (u, v) = halton_offset(k);
                out.push(make(x + u * jx, y + v * jy).cast::<T>());
            }
        }
        Scheme::Pcg => {
            let (i, j) = (idx % slice.nx, idx / slice.nx);
            let mut rng = SplitMix64::new(pixel_seed(i, j, seed));
            for _ in 0..n_extra {
                out.push(make(x + rng.range(-jx, jx), y + rng.range(-jy, jy)).cast::<T>());
            }
        }
    }
    out
}
