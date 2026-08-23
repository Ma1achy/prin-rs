//! Does `spread_shape` inherit the LC branch-cut sensitivity?
//!
//! The concern: copies of one pixel differ in configuration, so they can straddle the cut
//! differently, and a spread *across copies* would then partly measure registration error
//! rather than dynamics.
//!
//! Answer: at f64 the leak is negligible. **At f32 the unstable branch destroys the
//! statistic** — inflating it by an order of magnitude or returning NaN — while the
//! conditioned branch tracks the f64 answer to about 1%.
//!
//! That is the exact shape of the unresolved f32 dispute in the prior numpy work: single
//! trajectory drift looks acceptable, the ensemble diagnostic breaks early. The brief's
//! working hypothesis (BRIEF §3, §8 experiment 2) was reference-body switching across
//! copies. This is a demonstrated alternative that has nothing to do with switching, and
//! nothing to do with arithmetic.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

/// Per-pixel spreads, not a mean over them.
///
/// The mean was the original statistic and a single pixel controlled it: at 5x5 near-field under
/// the fixed Halton prefix, pixel 24 has an f64 spread of `2.131e-4` against an f32 spread of
/// `1.251e-1` — a chaotic divergence over a near-zero denominator — and it moved the mean by
/// 30% while the median relative error was 1.6%. This project's own rule about aggregates,
/// applied to its own test.
fn spreads<T: prin_rs::Real>(region: &str, lc_stable: bool) -> (Vec<f64>, usize) {
    let s = grid::region(region, 5, 5, 0.05).unwrap();
    let cfg = EnsembleCfg {
        t_max: 13.0,
        lc_stable,
        refine_flagged: false,
        ..Default::default()
    };
    let mut v = Vec::new();
    let mut n_bad = 0usize;
    for i in 0..s.npix() {
        let p = evaluate::<T>(&s, i, &cfg);
        if p.spread_shape.is_finite() {
            v.push(p.spread_shape);
        } else {
            n_bad += 1;
            v.push(f64::NAN);
        }
    }
    (v, n_bad)
}

/// Median over pixels of `|f32 - f64| / f64`.
fn median_rel(a: &[f64], b: &[f64]) -> f64 {
    let mut r: Vec<f64> = a
        .iter()
        .zip(b)
        .filter(|(x, y)| x.is_finite() && y.is_finite() && **x > 0.0)
        .map(|(x, y)| (y - x).abs() / x)
        .collect();
    r.sort_by(|p, q| p.partial_cmp(q).unwrap());
    r[r.len() / 2]
}

fn mean_spread<T: prin_rs::Real>(region: &str, lc_stable: bool) -> (f64, usize) {
    let s = grid::region(region, 5, 5, 0.05).unwrap();
    // Refinement off. It is triggered by a threshold on `error_ratio`, and the two precisions
    // flag different pixel sets — so with it on, this would compare which pixels happened to
    // get a second pass rather than comparing the arithmetic. Any precision comparison has to
    // hold the pipeline fixed.
    let cfg = EnsembleCfg { t_max: 13.0, lc_stable, refine_flagged: false, ..Default::default() };
    let mut acc = 0.0;
    let mut n_bad = 0usize;
    for i in 0..s.npix() {
        let p = evaluate::<T>(&s, i, &cfg);
        if p.spread_shape.is_finite() {
            acc += p.spread_shape;
        } else {
            n_bad += 1;
        }
    }
    (acc / s.npix() as f64, n_bad)
}

#[test]
fn at_f64_the_branch_cut_barely_leaks_into_spread_shape() {
    for region in ["near-field", "body2 core", "body1 slice", "mid-field"] {
        let (u, _) = mean_spread::<f64>(region, false);
        let (s, _) = mean_spread::<f64>(region, true);
        let rel = (u - s).abs() / s;
        println!("{region:>14}  unstable {u:.6e}  stable {s:.6e}  rel {rel:.2e}");
        assert!(rel < 1e-3, "{region}: f64 leak of {rel:e} is larger than expected");
    }
}

#[test]
fn at_f32_the_unstable_branch_destroys_spread_shape() {
    println!("{:>14}{:>14}{:>14}{:>14}{:>10}", "region", "f64 truth", "f32 unstable", "f32 stable", "NaN pix");
    let mut any_broken = false;
    for region in ["near-field", "body2 core", "body1 slice", "mid-field"] {
        let (truth, _) = mean_spread::<f64>(region, true);
        let (u32, bad_u) = mean_spread::<f32>(region, false);
        let (s32, bad_s) = mean_spread::<f32>(region, true);
        println!("{region:>14}{truth:>14.4e}{u32:>14.4e}{s32:>14.4e}{:>10}", bad_u);

        // The conditioned branch must track the f64 answer and produce no NaN pixels.
        assert_eq!(bad_s, 0, "{region}: stable f32 produced {bad_s} non-finite spreads");
        // Gated on the MEDIAN over pixels, not the mean: see `spreads`. The mean is controlled
        // by whichever single pixel has the smallest f64 denominator.
        let (f64s, _) = spreads::<f64>(region, true);
        let (f32s, _) = spreads::<f32>(region, true);
        let rel_s = median_rel(&f64s, &f32s);
        println!("{:>14}  median per-pixel rel err, conditioned f32 vs f64: {rel_s:.4}", "");
        assert!(rel_s < 0.05, "{region}: stable f32 median is {rel_s:e} from the f64 answer");

        // The unstable branch is expected to fail — either NaN, or badly wrong.
        let rel_u = (u32 - truth).abs() / truth;
        if bad_u > 0 || rel_u > 0.5 {
            any_broken = true;
        }
    }
    println!();
    println!("The conditioned branch holds at f32. The unstable one does not: it returns NaN");
    println!("pixels or inflates the spread by an order of magnitude, while single-trajectory");
    println!("energy drift stays superficially reasonable. That is the reported f32 symptom.");
    assert!(any_broken, "the unstable branch was expected to break spread_shape at f32");
}
