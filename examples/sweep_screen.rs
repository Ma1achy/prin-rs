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
        for tau in [1e-8f64, 1e-6, 1e-4, 1e-3, 1e-2] {
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
        let tau_span = at(1e-8, 0.2) / at(1e-6, 0.2).max(1.0);
        println!("  alpha_hi 0.20 -> 0.50 at tau=1e-4: leaf count x{alpha_span:.2}  (PR #11: x80)");
        println!("  tau 1e-8 -> 1e-6 at alpha_hi=0.20: leaf count x{tau_span:.2}  (PR #11: identical)");
        println!();
    }

    println!("split* counts SPLIT over all quads; a leaf can never carry Split.");
    println!();
    println!("With only four discretionary levels below the bootstrap, a threshold has far less");
    println!("room to express itself. A small alpha ratio here does NOT mean alpha stopped");
    println!("mattering — it means the veto stopped the descent before alpha could.");
}
