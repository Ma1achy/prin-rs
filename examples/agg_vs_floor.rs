//! **Does the aggregation fix collide with the screen floor?**
//!
//! Raised in the PR #12 review. `deep interior` under mean or p90 reaches **depth 7**, and the
//! screen floor at `N = 8` on a 512² viewport sits at **level 6**. If the fix wants a level the
//! view will not grant, the improvement is capped in exactly the configuration production runs.
//!
//! The floor is arithmetic: one sample one tile, so a fully-refined tree at level `L` holds
//! `4^L * N²` samples, and `4^L * 64 <= V²` gives `L = 6` at 512, **7 at 1024**, 8 at 2048. So the
//! question has an answer that does not depend on the physics — the viewport at which level 7
//! becomes displayable is 1024² — and this run confirms it end to end and measures what the
//! truncation actually costs each aggregation.
//!
//! **Read the depth histograms, not the leaf counts.** A capped tree and an uncapped one can have
//! similar leaf totals while differing entirely in where those leaves sit.
//!
//! **And read BOTH cap columns.** A first version of this reported only `screen`, and at 1024²
//! that column fell to zero while the tree did not grow by a single quad — which reads as "the
//! viewport is inert". It is not: `MAX_REL_DEPTH = 6` had taken over as the binding cap, and it is
//! a *policy default*, not arithmetic. The two coincide at 512² and diverge above it, so a table
//! that shows one of them describes the wrong system. The `rel` rows below raise `MAX_REL_DEPTH`
//! to the viewport'"'"'s own floor, which is what the contract'"'"'s `MAX_REL_DEPTH <= screen floor`
//! actually permits.
//!
//! Run: `cargo run --release --example agg_vs_floor [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::{Agg, Decision as D};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let budget: usize = arg(1, 6000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    const N: usize = 8;

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N={N}, E+1=8, t=13, f64.\n");
    println!("Screen floor by viewport, from 4^L * N^2 <= V^2:");
    for v in [512usize, 1024, 2048, 4096] {
        let l = (0..20).take_while(|&l| 4usize.pow(l) * N * N <= v * v).count().saturating_sub(1);
        println!("  viewport {v:>5}^2 -> deepest displayable level {l}");
    }
    println!();
    println!("{:>14} {:>9} {:>5} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>7} {:>7} {:>28}",
             "region", "viewport", "relD", "agg", "quads", "leaves", "depth", "floor", "keep",
             "screen", "relcap", "depth histogram");

    for region in ["deep interior", "near-field"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        // (viewport, max_rel_depth). `None` viewport is the no-camera control; the paired rows
        // at each viewport hold MAX_REL_DEPTH at its default 6 and then at the viewport's own
        // screen-floor level, so the two caps can be told apart.
        for (vp, rel) in [(512usize, 6u32), (1024, 6), (1024, 7), (2048, 6), (2048, 8), (0, 0)] {
            for agg in [Agg::Median, Agg::Mean, Agg::P90] {
                let cam = (vp > 0).then(|| {
                    let mut c = Camera::framing(root.cx, root.cy, 0.05, vp);
                    c.max_rel_depth = Some(rel);
                    c
                });
                let cfg = SchedCfg {
                    n: N, budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4, agg,
                    camera: cam,
                    ..Default::default()
                };
                let (t, st) = scheduler::descend(
                    root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
                let leaves: Vec<usize> = t.leaves().collect();
                let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
                let hist = t.depth_histogram();
                let hs: String = hist.iter().enumerate().filter(|(_, &n)| n > 0)
                    .map(|(l, n)| format!("{l}:{n} ")).collect();
                println!("{:>14} {:>9} {:>5} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>7} {:>7} {:>28}",
                         region,
                         if vp == 0 { "none".into() } else { format!("{vp}^2") },
                         if vp == 0 { "-".into() } else { rel.to_string() },
                         agg.name(), st.quads_computed, leaves.len(),
                         hist.len().saturating_sub(1), c(D::Floor), c(D::Keep),
                         c(D::ScreenFloor), c(D::MaxRelDepth), hs.trim());
            }
            println!();
        }
    }

    println!("If a row's histogram is identical to the no-viewport row up to the cap level and");
    println!("then piles up at it, the veto TRUNCATED that aggregation rather than changing it.");
    println!("The share of leaves in that pile is what the truncation cost.");
    println!();
    println!("`screen` and `relcap` are different caps and a row can be bound by either. Where");
    println!("relcap is nonzero and screen is zero, the SCREEN FLOOR was not what stopped the");
    println!("descent — MAX_REL_DEPTH was, and that is a policy default rather than arithmetic.");
}
