//! **§4 question 7: threshold sensitivity.** And it runs first, because it is what sets `tau`.
//!
//! The measured quad spread spans six orders across regions — `~4e-8` in `far`, `~2e-3` in
//! near-field, median `7.5e-5` but p90 `1.4e-1` in `deep interior`. No single `tau_display` can
//! serve all three, so any `tau` picked before this sweep would be an arbitrary constant, which is
//! the defect that has already disqualified two candidate designs on this project.
//!
//! **The sweep is the result. The picture is a diagnostic.**
//!
//! A criterion whose output is dominated by an arbitrary threshold is not a criterion — that
//! finding would matter, so it is measured rather than assumed away.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

const BUDGET: usize = 2000;

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        println!("=== {region} ===  budget {BUDGET} quads, N=8, E+1=8, t=13, agg=median");
        println!("{:>10}{:>10}{:>9}{:>9}{:>8}{:>8}{:>8}{:>8}{:>7}{:>9}",
                 "tau", "alpha_hi", "quads", "leaves", "split*", "floor", "keep", "budget",
                 "depth", "wall s");

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
                    budget: BUDGET,
                    tau_display: tau,
                    alpha_hi,
                    alpha_lo: alpha_hi * 0.4,
                    ..Default::default()
                };
                let (t, st) =
                    scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
                let leaves: Vec<usize> = t.leaves().collect();
                let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
                // A leaf can never be Split — deciding to split makes it an internal node — so
                // the split count is taken over ALL quads, not over leaves.
                let splits = t.nodes.iter().filter(|q| q.decision == D::Split).count();
                let depth = t.depth_histogram().len().saturating_sub(1);
                println!("{tau:>10.0e}{alpha_hi:>10.2}{:>9}{:>9}{splits:>8}{:>8}{:>8}{:>8}{depth:>7}{:>9.1}",
                         st.quads_computed, leaves.len(), c(D::Floor), c(D::Keep),
                         c(D::BudgetExhausted), st.wall_seconds);
            }
        }
        println!();
    }

    println!("split* counts SPLIT decisions over all quads: a leaf can never be Split, since");
    println!("deciding to split turns it into an internal node.");
    println!();
    println!("Read the columns, not a single row. If leaf count and depth move by orders across");
    println!("the sweep, the threshold is doing the work rather than the criterion. If they");
    println!("plateau over a range, that plateau is where tau can be set with a reason.");
}
