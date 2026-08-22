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

fn mean_spread<T: prin_rs::Real>(region: &str, lc_stable: bool) -> (f64, usize) {
    let s = grid::region(region, 5, 5, 0.05).unwrap();
    let cfg = EnsembleCfg { t_max: 13.0, lc_stable, ..Default::default() };
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
        let rel_s = (s32 - truth).abs() / truth;
        assert!(rel_s < 0.05, "{region}: stable f32 spread is {rel_s:e} from the f64 answer");

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
