//! **§4.2 — q7's threshold sweep, re-run under the screen floor.**
//!
//! PR #11's answer was *"`alpha_hi` from 0.20 to 0.50 collapses the tree 80x, while `tau` is inert
//! over four orders."* That was measured over levels 0–12. Under the veto the descent has
//! `bootstrap_levels = 2` to level 6 — **four discretionary levels** — and the prediction on record
//! is that `alpha_hi`'s effect shrinks by more than an order.
//!
//! If `alpha_hi` still dominates over four levels, that is the *stronger* result, not a
//! disappointment. Both are reported; neither is tuned toward.
//!
//! Run: `cargo run --release --example sweep_screen [budget]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn main() {
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("q7 under the screen floor. budget {budget} quads, N=8, E+1=8, t=13, agg=median,");
    println!("viewport 512x512, camera framing the root box. Compare each block against the");
    println!("same sweep without a veto in results/output/sched_sweep.txt.\n");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        println!("=== {region} ===");
        println!("{:>10}{:>10}{:>9}{:>9}{:>8}{:>8}{:>8}{:>8}{:>7}{:>9}",
                 "tau", "alpha_hi", "quads", "leaves", "split*", "floor", "keep", "screen",
                 "depth", "wall s");
        let mut leaf_counts: Vec<(f64, f64, usize)> = Vec::new();
        // **The ladder was measured, and the first cut of it was wrong.** Pooled over the 89,088
        // committed leaves the spread median is 6.6e-4 and `tau = 1e-4` sits at the 0.4th
        // percentile, so the top of the old ladder measured nothing: 1e-6, 1e-4 and 3e-4 give a
        // BITWISE IDENTICAL tree in near-field. It now runs up to 1e-1, past the point where the
        // predicate goes false everywhere.
        //
        // **The bottom rungs are NOT redundant, and dropping 1e-8 broke `far`.** The regional
        // spread medians span six orders -- 4.26e-8 in `far`, 9.45e-5 in `deep interior`,
        // 9.75e-4 in near-field -- so which rung is degenerate is a fact about the REGION, not
        // about the ladder. `1e-8` is the only rung below `far`'s bulk; without it `far` reads 16
        // leaves at every cell and the sweep says "tau is inert here", which is a statement about
        // the ladder. Both low rungs stay, as labelled degenerate controls for different regions.
        // See `examples/threshold_diagnosis.rs` for the percentiles.
        for tau in [1e-8f64, 1e-6, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1] {
            for alpha_hi in [0.2f64, 0.5, 0.8, 1.0] {
                let cfg = SchedCfg {
                    budget,
                    tau_display: tau,
                    alpha_hi,
                    alpha_lo: alpha_hi * 0.4,
                    camera: Some(Camera::framing(root.cx, root.cy, 0.05, 512)),
                    ..Default::default()
                };
                let (t, st) = scheduler::descend(
                    root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
                let leaves: Vec<usize> = t.leaves().collect();
                let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
                let splits = t.nodes.iter().filter(|q| q.decision == D::Split).count();
                let depth = t.depth_histogram().len().saturating_sub(1);
                leaf_counts.push((tau, alpha_hi, leaves.len()));
                println!("{tau:>10.0e}{alpha_hi:>10.2}{:>9}{:>9}{splits:>8}{:>8}{:>8}{:>8}{depth:>7}{:>9.1}",
                         st.quads_computed, leaves.len(), c(D::Floor), c(D::Keep),
                         c(D::ScreenFloor), st.wall_seconds);
            }
        }
        // The two ratios PR #11 quoted, recomputed here so the comparison is arithmetic.
        let at = |t: f64, a: f64| {
            leaf_counts.iter().find(|x| x.0 == t && x.1 == a).map(|x| x.2).unwrap_or(0) as f64
        };
        let alpha_span = at(1e-4, 0.2) / at(1e-4, 0.5).max(1.0);
        // The tau span is taken over the WHOLE ladder at the one alpha_hi where tau is live,
        // not between two adjacent rungs. Two adjacent rungs can both sit below the region's bulk
        // and read "identical", which is a fact about the rungs; the max/min over the ladder is
        // the honest span. Regions differ by six orders in spread, so a fixed pair cannot serve
        // all three -- that is the mistake this line used to make.
        let live: Vec<f64> =
            leaf_counts.iter().filter(|x| x.1 == 0.2).map(|x| x.2 as f64).collect();
        let lo = live.iter().cloned().fold(f64::INFINITY, f64::min).max(1.0);
        let hi = live.iter().cloned().fold(0.0f64, f64::max);
        println!("  alpha_hi 0.20 -> 0.50 at tau=1e-4: leaf count x{alpha_span:.2}  (PR #11: x80)");
        println!("  tau over the whole ladder at alpha_hi=0.20: leaf count x{:.2}  ({lo} .. {hi})",
                 hi / lo);
        println!();
    }

    println!("split* counts SPLIT over all quads; a leaf can never carry Split.");
    println!();
    println!("With only four discretionary levels below the bootstrap, a threshold has far less");
    println!("room to express itself. A small alpha ratio here does NOT mean alpha stopped");
    println!("mattering — it means the veto stopped the descent before alpha could.");
}
