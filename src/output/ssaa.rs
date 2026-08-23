//! SSAA resolve: many sub-pixel samples, one pixel colour.
//!
//! The ensemble has been used for exactly one thing so far — `ensemble_spread`, a
//! **disagreement** statistic that drives scheduling. Its other job is **resolve**, an
//! *average* that drives display, and that path had never run.
//!
//! Keeping the two apart is not pedantry. A footprint whose copies split 4/4 between two
//! outcomes has a large spread *and* a blended colour; the spread says "refine here", the
//! blend says "this is what the pixel looks like". Substituting either for the other loses a
//! real distinction — and a resolved image that happened to equal the nominal-copy image would
//! mean the ensemble is doing no anti-aliasing at all, which is a finding, not a null.
//!
//! The copies are jittered by `jitter_frac * cell_width`, so their sub-pixel footprint scales
//! with the texel size **by construction** — a level-3 leaf's copies spread over 4x the world
//! distance of a level-5 leaf's, matching its 4x texel. Asserted rather than assumed.

use crate::ensemble::pixel::PixelOut;
use crate::outcome::State;
use crate::output::png::outcome_rgb;

/// Colour for a packed `(state, detail)` byte — the same palette the nominal-copy image uses,
/// so a resolved render and a nominal render differ only by the averaging.
pub fn packed_rgb(packed: u8) -> [u8; 3] {
    let (state, detail) = (packed >> 2, packed & 0b11);
    let base = match State::from_bits(state) {
        Some(State::Escape) => [220, 80, 60],
        Some(State::Collision) => [110, 190, 110],
        Some(State::Bounded) => [70, 150, 220],
        Some(State::Running) => [200, 190, 90],
        Some(State::SimFailed) => return [255, 0, 255],
        _ => [40, 40, 48],
    };
    let k = 0.55 + 0.15 * detail as f64;
    [
        (base[0] as f64 * k).min(255.0) as u8,
        (base[1] as f64 * k).min(255.0) as u8,
        (base[2] as f64 * k).min(255.0) as u8,
    ]
}

/// The resolve: the mean colour over the `E+1` copies.
///
/// Falls back to the nominal-copy colour when the copies were not retained, so a caller that
/// forgot [`crate::ensemble::pixel::EnsembleCfg::keep_copy_outcomes`] gets the nominal image
/// rather than a black one — and [`resolve_available`] says which happened.
pub fn resolve_rgb(p: &PixelOut) -> [u8; 3] {
    if p.copy_outcomes.is_empty() {
        return outcome_rgb(p);
    }
    let n = p.copy_outcomes.len() as f64;
    let mut acc = [0.0f64; 3];
    for &c in &p.copy_outcomes {
        let rgb = packed_rgb(c);
        for k in 0..3 {
            acc[k] += rgb[k] as f64;
        }
    }
    [
        (acc[0] / n).round() as u8,
        (acc[1] / n).round() as u8,
        (acc[2] / n).round() as u8,
    ]
}

pub fn resolve_available(px: &[PixelOut]) -> bool {
    px.iter().any(|p| !p.copy_outcomes.is_empty())
}

/// How far the resolved image moves from the nominal-copy image.
///
/// Reported per pixel, not only as an aggregate: **an aggregate can only say the distribution
/// did not move, never that the pixels did not.** Returns
/// `(mean |dRGB|, max |dRGB|, fraction of pixels that moved at all)`.
pub fn resolve_effect(px: &[PixelOut]) -> (f64, f64, f64) {
    let (mut sum, mut worst, mut moved) = (0.0f64, 0.0f64, 0usize);
    for p in px {
        let (a, b) = (resolve_rgb(p), outcome_rgb(p));
        let d = (0..3).map(|k| (a[k] as f64 - b[k] as f64).abs()).fold(0.0, f64::max);
        sum += d;
        worst = worst.max(d);
        if d > 0.0 {
            moved += 1;
        }
    }
    let n = px.len().max(1) as f64;
    (sum / n, worst, moved as f64 / n)
}
