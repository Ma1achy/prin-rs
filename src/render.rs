//! Rayon over pixels, and the single point where precision is dispatched.
//!
//! Both arms compile into one binary, so BRIEF §8 experiment 2 runs f32 and f64 on identical
//! initial conditions in a single process — no rebuild, and no chance of a config drifting
//! between runs. A feature flag would make the two comparable only across separate builds,
//! where an IC-generation difference would contaminate the comparison silently and look
//! exactly like a real f32 effect.

use rayon::prelude::*;

use crate::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use crate::grid::Slice;
use crate::Real;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    F32,
    F64,
}

impl Precision {
    pub fn name(&self) -> &'static str {
        match self {
            Precision::F32 => "f32",
            Precision::F64 => "f64",
        }
    }
}

fn render_typed<T: Real>(slice: &Slice, cfg: &EnsembleCfg) -> Vec<PixelOut> {
    (0..slice.npix())
        .into_par_iter()
        .map(|i| evaluate::<T>(slice, i, cfg))
        .collect()
}

pub fn render(slice: &Slice, cfg: &EnsembleCfg, precision: Precision) -> Vec<PixelOut> {
    match precision {
        Precision::F32 => render_typed::<f32>(slice, cfg),
        Precision::F64 => render_typed::<f64>(slice, cfg),
    }
}

/// Summary statistics over a rendered grid.
///
/// `error_ratio` aggregates **across pixels by max**, per BRIEF §4 — max tracks damage at
/// Spearman +0.956 against +0.599 for median. That result is about this aggregation and
/// nothing else; the separate choice of maximum deviation as the estimator *within* a
/// footprint rests on its own measurement (see [`crate::ensemble::stats`]). Two decisions,
/// both spelled "max", independently arrived at.
///
/// The magnitude is unstable, so it is a boolean flag; the median and p99 are reported
/// alongside so the distribution is visible rather than reduced to one untrustworthy number.
#[derive(Clone, Debug, Default)]
pub struct Summary {
    pub n: usize,
    pub error_ratio_max: f64,
    pub error_ratio_median: f64,
    pub error_ratio_p99: f64,
    pub error_ratio_argmax: usize,
    pub drift_median: f64,
    pub drift_max: f64,
    pub d_min_gap_median: f64,
    pub d_min_gap_max: f64,
    pub ref_disagree_total: u64,
    /// Fraction of pixels in each [`crate::outcome::State`], by the nominal copy's label.
    pub state_fracs: [f64; 6],
    pub t_end_median: f64,
    pub n_refined: usize,
    /// Pixels still above the refinement threshold after the last pass. Reported, not hidden.
    pub n_still_flagged: usize,
    pub drift_max_coarse: f64,
    pub n_pixels_with_nonfinite: usize,
    pub sigma_e_0_median: f64,
}

fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * q).round() as usize]
}

pub fn summarise(pixels: &[PixelOut]) -> Summary {
    let finite = |x: &f64| x.is_finite();
    let mut er: Vec<f64> = pixels.iter().map(|p| p.error_ratio).filter(finite).collect();
    let mut dr: Vec<f64> = pixels.iter().map(|p| p.energy_drift_max).filter(finite).collect();
    let mut gap: Vec<f64> = pixels.iter().map(|p| p.d_min_gap).filter(finite).collect();
    let mut s0: Vec<f64> = pixels.iter().map(|p| p.sigma_e_0).filter(finite).collect();
    let mut te: Vec<f64> = pixels.iter().map(|p| p.t_end).filter(finite).collect();

    let (mut argmax, mut best) = (0usize, f64::NEG_INFINITY);
    for (i, p) in pixels.iter().enumerate() {
        if p.error_ratio.is_finite() && p.error_ratio > best {
            best = p.error_ratio;
            argmax = i;
        }
    }

    Summary {
        n: pixels.len(),
        error_ratio_max: quantile(&mut er.clone(), 1.0),
        error_ratio_median: quantile(&mut er, 0.5),
        error_ratio_p99: quantile(&mut er.clone(), 0.99),
        error_ratio_argmax: argmax,
        drift_median: quantile(&mut dr.clone(), 0.5),
        drift_max: quantile(&mut dr, 1.0),
        d_min_gap_median: quantile(&mut gap.clone(), 0.5),
        d_min_gap_max: quantile(&mut gap, 1.0),
        ref_disagree_total: pixels.iter().map(|p| p.ref_disagree as u64).sum(),
        state_fracs: {
            let mut f = [0.0f64; 6];
            for p in pixels {
                if (p.state as usize) < 6 {
                    f[p.state as usize] += 1.0;
                }
            }
            let n = pixels.len().max(1) as f64;
            f.map(|x| x / n)
        },
        t_end_median: quantile(&mut te, 0.5),
        n_refined: pixels.iter().filter(|p| p.refined).count(),
        n_still_flagged: pixels.iter().filter(|p| p.error_ratio > 10.0).count(),
        drift_max_coarse: quantile(
            &mut pixels.iter().map(|p| p.energy_drift_max_coarse).filter(finite).collect(),
            1.0,
        ),
        n_pixels_with_nonfinite: pixels.iter().filter(|p| p.n_nonfinite > 0).count(),
        sigma_e_0_median: quantile(&mut s0, 0.5),
    }
}
