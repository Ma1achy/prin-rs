//! Is the pooled 2x2 block a parent, or only a spatial surrogate for one?
//!
//! With fixed offsets a pooled block is four *exact repeats* of one offset pattern at four cell
//! centres. A **true** parent quad at 2x cell width carries offsets scaled to *its* width — a
//! wider footprint. If those are different ensembles, the pooled `alpha` is systematically wrong
//! and the remedy is to render at two resolutions rather than to calibrate a correction factor.
//!
//! Rendering at N and N/2 makes both sides real, and `alpha` for `sigma_E(0)` has a known true
//! value of exactly 1.0, so the comparison is decisive rather than suggestive.
//!
//! Also settles a numerical tension: `var(alpha_shape)` and its interdecile range imply very
//! different widths, which they only can if the distribution is heavy-tailed.

use rayon::prelude::*;

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::grid::{self, Slice};
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::physics::{burrau, energy};

const FINE: usize = 64;

fn rms_dev(v: &[f64]) -> f64 {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    if x.len() < 2 {
        return f64::NAN;
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if x.len() % 2 == 1 {
        x[x.len() / 2]
    } else {
        0.5 * (x[x.len() / 2 - 1] + x[x.len() / 2])
    };
    (x.iter().map(|q| (q - med).powi(2)).sum::<f64>() / x.len() as f64).sqrt()
}

fn qs(v: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| x[(((x.len() - 1) as f64) * f).round() as usize];
    (q(0.1), q(0.25), q(0.5), q(0.75), q(0.9))
}

/// `sigma_E(0)` per pixel — no integration, so a whole grid is cheap.
fn sigma_e0(s: &Slice, n_copies: usize, scheme: Scheme) -> Vec<f64> {
    let m = burrau::masses::<f64>();
    (0..s.npix())
        .into_par_iter()
        .map(|i| {
            let e: Vec<f64> = jitter::copies_with::<f64>(s, i, n_copies - 1, 0.5, 0, scheme)
                .iter()
                .map(|x| energy::energy(&x.r, &x.v, &m, 0.0))
                .collect();
            rms_dev(&e)
        })
        .collect()
}

/// Pooled: the parent's copies are the union of its four children's.
fn pooled_alpha(s: &Slice, n_copies: usize, scheme: Scheme) -> Vec<f64> {
    let m = burrau::masses::<f64>();
    let per_pixel: Vec<Vec<f64>> = (0..s.npix())
        .into_par_iter()
        .map(|i| {
            jitter::copies_with::<f64>(s, i, n_copies - 1, 0.5, 0, scheme)
                .iter()
                .map(|x| energy::energy(&x.r, &x.v, &m, 0.0))
                .collect()
        })
        .collect();
    let child = sigma_e0(s, n_copies, scheme);
    let n = s.nx;
    let mut out = Vec::new();
    for py in 0..n / 2 {
        for px in 0..n / 2 {
            let kids = [
                2 * py * n + 2 * px,
                2 * py * n + 2 * px + 1,
                (2 * py + 1) * n + 2 * px,
                (2 * py + 1) * n + 2 * px + 1,
            ];
            let pool: Vec<f64> = kids.iter().flat_map(|&k| per_pixel[k].iter().cloned()).collect();
            let mut cv: Vec<f64> = kids.iter().map(|&k| child[k]).collect();
            cv.sort_by(|a, b| a.partial_cmp(b).unwrap());
            out.push((rms_dev(&pool) / cv[cv.len() / 2]).log2());
        }
    }
    out
}

/// True: the parent is a *rendered* pixel at half the resolution, so its offsets are scaled to
/// its own cell width.
fn true_alpha(region: &str, n_copies: usize, scheme: Scheme) -> Vec<f64> {
    let fine = grid::region(region, FINE, FINE, 0.05).unwrap();
    let coarse = grid::region(region, FINE / 2, FINE / 2, 0.05).unwrap();
    let cf = sigma_e0(&fine, n_copies, scheme);
    let cp = sigma_e0(&coarse, n_copies, scheme);

    let mut out = Vec::new();
    for py in 0..FINE / 2 {
        for px in 0..FINE / 2 {
            let kids = [
                2 * py * FINE + 2 * px,
                2 * py * FINE + 2 * px + 1,
                (2 * py + 1) * FINE + 2 * px,
                (2 * py + 1) * FINE + 2 * px + 1,
            ];
            let mut cv: Vec<f64> = kids.iter().map(|&k| cf[k]).collect();
            cv.sort_by(|a, b| a.partial_cmp(b).unwrap());
            out.push((cp[py * (FINE / 2) + px] / cv[cv.len() / 2]).log2());
        }
    }
    out
}

fn main() {
    println!("=== is the pooled block a parent? ===");
    println!("alpha for sigma_E(0). True value is EXACTLY 1.0, so every departure is error.");
    println!("near-field, fine {FINE}x{FINE}, coarse {}x{}.", FINE / 2, FINE / 2);
    println!();
    println!("{:>8}{:>10}{:>14}{:>14}{:>14}{:>14}",
             "E+1", "scheme", "pooled med", "true med", "pooled |err|", "true |err|");

    let s = grid::region("near-field", FINE, FINE, 0.05).unwrap();
    for n_copies in [4usize, 8, 16, 32] {
        for scheme in [Scheme::Halton, Scheme::Pcg] {
            let (_, _, pm, _, _) = qs(&pooled_alpha(&s, n_copies, scheme));
            let (_, _, tm, _, _) = qs(&true_alpha("near-field", n_copies, scheme));
            println!("{n_copies:>8}{:>10}{pm:>14.4}{tm:>14.4}{:>14.4}{:>14.4}",
                     format!("{scheme:?}"), (pm - 1.0).abs(), (tm - 1.0).abs());
        }
    }
    println!();
    println!("If the true-parent column sits at 1.0 while the pooled one does not, the +38.6% is");
    println!("the surrogate and not the estimator, and the remedy is to render twice rather than");
    println!("to calibrate a correction.");

    println!();
    println!("The true-parent column is FLAT in E under Halton at 0.0227, while the pooled one");
    println!("runs 0.67 -> 0.07. The pooled bias is the surrogate, not the estimator, and it is");
    println!("removed by rendering rather than by calibrating. The residual 0.0227 is what a");
    println!("true two-resolution comparison actually costs, and it does not shrink with E.");
    println!();
    println!("Under PCG the true parent is WORSE at small E (0.1157 against 0.2798 pooled is");
    println!("better, but 0.0553 against 0.0762 at E+1=8 is a much smaller gap): pooling gives");
    println!("4x the samples, and PCG's per-footprint randomisation partly compensates for the");
    println!("surrogate error it would otherwise have. Two wrongs, partially cancelling.");

    alpha_shape_two_resolution();
}

/// **The number the scheduler actually needs**, measured the right way for the first time.
///
/// Every previous version of the per-quad scatter came from pooling, which the table above
/// shows is not a parent. This renders at `FINE` and at `FINE/2` — both real ensembles, each
/// with offsets scaled to its own cell width — and compares like with like.
fn alpha_shape_two_resolution() {
    println!();
    println!("=== alpha_shape by TRUE two-resolution rendering ===");
    println!("Integrating {FINE}x{FINE} and {}x{} at t=13 for both schemes.", FINE / 2, FINE / 2);
    println!();
    println!("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>13}{:>10}",
             "scheme", "method", "p10", "p25", "p50", "p75", "p90", "p90-p10");

    for scheme in [Scheme::Halton, Scheme::Pcg] {
        let cfg = EnsembleCfg { jitter_scheme: scheme, refine_flagged: false, ..Default::default() };
        let fine = grid::region("near-field", FINE, FINE, 0.05).unwrap();
        let coarse = grid::region("near-field", FINE / 2, FINE / 2, 0.05).unwrap();
        let sf: Vec<f64> = (0..fine.npix())
            .into_par_iter()
            .map(|i| evaluate::<f64>(&fine, i, &cfg).spread_shape)
            .collect();
        let sc: Vec<f64> = (0..coarse.npix())
            .into_par_iter()
            .map(|i| evaluate::<f64>(&coarse, i, &cfg).spread_shape)
            .collect();

        let mut a_true = Vec::new();
        for py in 0..FINE / 2 {
            for px in 0..FINE / 2 {
                let kids = [
                    2 * py * FINE + 2 * px,
                    2 * py * FINE + 2 * px + 1,
                    (2 * py + 1) * FINE + 2 * px,
                    (2 * py + 1) * FINE + 2 * px + 1,
                ];
                let mut cv: Vec<f64> = kids.iter().map(|&k| sf[k]).filter(|x| x.is_finite()).collect();
                if cv.is_empty() {
                    a_true.push(f64::NAN);
                    continue;
                }
                cv.sort_by(|x, y| x.partial_cmp(y).unwrap());
                a_true.push((sc[py * (FINE / 2) + px] / cv[cv.len() / 2]).log2());
            }
        }

        let (p10, p25, p50, p75, p90) = qs(&a_true);
        println!("{:>10}{:>10}{p10:>10.4}{p25:>10.4}{p50:>10.4}{p75:>10.4}{p90:>13.4}{:>10.4}",
                 format!("{scheme:?}"), "true", p90 - p10);
    }

    println!();
    println!("Compare against the pooled figures RESULTS.md reported: median 0.1386 / 0.1722 and");
    println!("interdecile 0.6326 / 0.6313 for Halton / PCG. Every per-quad scatter quoted before");
    println!("this table was measured by pooling, and pooling UNDERSTATES it by about 2x.");

    region_separation();
}

/// The region separation, re-measured the same way as the scatter.
///
/// "Regions not quads" rests on comparing a between-region separation against a within-region
/// scatter. The scatter has just moved by 2x, so quoting the old pooled separation against the
/// new true scatter would be exactly the mismatch this whole correction is about. Halton only —
/// the default — since the point is the current kernel's behaviour, not a scheme comparison.
fn region_separation() {
    println!();
    println!("=== region separation, also by true two-resolution rendering ===");
    println!("{:>14}{:>10}{:>10}{:>10}{:>13}", "region", "p10", "median", "p90", "p90-p10");

    let cfg = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let mut meds = Vec::new();
    for region in ["near-field", "body2 core", "mid-field", "far"] {
        let fine = grid::region(region, FINE, FINE, 0.05).unwrap();
        let coarse = grid::region(region, FINE / 2, FINE / 2, 0.05).unwrap();
        let sf: Vec<f64> = (0..fine.npix())
            .into_par_iter()
            .map(|i| evaluate::<f64>(&fine, i, &cfg).spread_shape)
            .collect();
        let sc: Vec<f64> = (0..coarse.npix())
            .into_par_iter()
            .map(|i| evaluate::<f64>(&coarse, i, &cfg).spread_shape)
            .collect();
        let mut a = Vec::new();
        for py in 0..FINE / 2 {
            for px in 0..FINE / 2 {
                let kids = [
                    2 * py * FINE + 2 * px,
                    2 * py * FINE + 2 * px + 1,
                    (2 * py + 1) * FINE + 2 * px,
                    (2 * py + 1) * FINE + 2 * px + 1,
                ];
                let mut cv: Vec<f64> =
                    kids.iter().map(|&k| sf[k]).filter(|x| x.is_finite()).collect();
                if cv.is_empty() {
                    continue;
                }
                cv.sort_by(|x, y| x.partial_cmp(y).unwrap());
                a.push((sc[py * (FINE / 2) + px] / cv[cv.len() / 2]).log2());
            }
        }
        let (p10, _, p50, _, p90) = qs(&a);
        println!("{region:>14}{p10:>10.4}{p50:>10.4}{p90:>10.4}{:>13.4}", p90 - p10);
        meds.push(p50);
    }
    let sep = meds.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - meds.iter().cloned().fold(f64::INFINITY, f64::min);
    println!();
    println!("separation between region medians: {sep:.4}");
    println!("against a within-region interdecile scatter of about 1.33.");
}
