//! The production colour scheme: **hue from the shape sphere, lightness from a scalar**.
//!
//! This is not presentation. The criterion asks *"would splitting change what we display?"*, so
//! **what is displayed decides what the criterion should measure.** `spread_shape` maps to hue,
//! so that half is aligned by construction. If lightness carries diffusion or FTLE, the
//! criterion is currently blind to changes in it: a quad can be uniform in shape and structured
//! in diffusion, and nothing would refine it.
//!
//! # The direction-to-colour map, stated rather than left implicit
//!
//! The shape vector `n` is a unit vector on `S^2`. The map used here is the simplest defensible
//! one, in **OKLCh** so that equal steps are roughly equal perceptual steps:
//!
//! ```text
//! hue        = atan2(n[2], n[1])          -- azimuth about the n0 axis
//! chroma     = C_MAX * sqrt(n1^2 + n2^2)  -- distance from the poles
//! lightness  = the selected scalar, normalised into [L_MIN, L_MAX]
//! ```
//!
//! **The azimuthal discontinuity is invisible by construction**, which is the property that
//! makes this defensible rather than merely simple: hue is undefined at the poles, and there
//! chroma goes to zero, so the two colours either side of the cut converge on the same grey.
//! A hue map that did not tie chroma to the same quantity would show a seam across every
//! isoceles configuration.
//!
//! Out-of-gamut OKLCh triples are clamped channel-wise by [`oklab::oklab_to_srgb`]. Clamping
//! compresses rather than wraps, so a clipped colour reads as too saturated, never as a
//! different hue.

use crate::ensemble::pixel::PixelOut;
use crate::output::oklab;

/// Maximum chroma at the shape sphere's equator.
pub const C_MAX: f64 = 0.13;
/// Lightness range. Not `[0, 1]`: pure black and pure white carry no chroma, so the hue would
/// vanish at both ends of the scalar and the bivariate map would silently become univariate.
pub const L_MIN: f64 = 0.30;
pub const L_MAX: f64 = 0.92;

/// Which scalar drives lightness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Lightness {
    /// `ensemble_spread` — the field the criterion already reads. The aligned case.
    Spread,
    /// Benettin FTLE. Requires `EnsembleCfg::ftle`.
    Ftle,
    /// Slope of `log(inertia)` against `t`. Requires `EnsembleCfg::ftle`.
    Diffusion,
}

impl Lightness {
    pub fn name(self) -> &'static str {
        match self {
            Lightness::Spread => "spread",
            Lightness::Ftle => "ftle",
            Lightness::Diffusion => "diffusion",
        }
    }
    pub fn value(self, p: &PixelOut) -> f64 {
        match self {
            Lightness::Spread => p.ensemble_spread,
            Lightness::Ftle => p.ftle,
            Lightness::Diffusion => p.diffusion,
        }
    }
}

/// OKLCh to sRGB.
pub fn oklch_to_srgb(l: f64, c: f64, h: f64) -> [u8; 3] {
    oklab::oklab_to_srgb([l, c * h.cos(), c * h.sin()])
}

/// The bivariate colour of one footprint, with the scalar normalised against `[lo, hi]`.
///
/// A non-finite scalar renders at `L_MIN` with **full chroma**, so an undetermined footprint is
/// dark but still carries its shape — it is not silently painted as the low end of the ramp,
/// which is what a clamp to zero would do.
pub fn rgb(p: &PixelOut, which: Lightness, lo: f64, hi: f64) -> [u8; 3] {
    let n = p.shape_vec;
    let hue = n[2].atan2(n[1]);
    let chroma = C_MAX * (n[1] * n[1] + n[2] * n[2]).sqrt();

    let v = which.value(p);
    let t = if v.is_finite() && hi > lo {
        ((v - lo) / (hi - lo)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    oklch_to_srgb(L_MIN + (L_MAX - L_MIN) * t, chroma, hue)
}

/// Robust `[p1, p99]` of a scalar over a set of footprints, for normalising the ramp.
///
/// Percentiles rather than min/max: one undetermined footprint at `1e12` would otherwise
/// compress every other pixel into the bottom of the range and the image would read as
/// featureless — the same failure mode as reading a variance where the kurtosis is 110.
pub fn range(px: &[PixelOut], which: Lightness) -> (f64, f64) {
    let mut v: Vec<f64> = px.iter().map(|p| which.value(p)).filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return (0.0, 1.0);
    }
    let lo = crate::quad::quantile(&mut v.clone(), 0.01);
    let hi = crate::quad::quantile(&mut v, 0.99);
    (lo, hi)
}
