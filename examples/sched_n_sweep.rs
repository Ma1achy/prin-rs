//! **The N sweep: does a coarse quad misclassify itself as coherent by undersampling its area?**
//!
//! `N` is the quality/compute driver and it interacts with the criterion directly, so the question
//! is whether the tree changes **shape**, not just cost. If a coarse `N` systematically
//! under-splits, that matters before any tier defaults are set.
//!
//! **`N = 7` is in the sweep for a reason that is not the obvious one.** The child sample grid is a
//! strict refinement of the parent's with a shared origin, so *every* parent sample inside a child
//! is also a child sample — for every `N >= 2`, verified. Copy 0 is un-jittered, so at those
//! footprints the parent and child nominal trajectories are **identical**: common random numbers
//! between the two scales `alpha` is a ratio of, arriving by accident of the grid.
//!
//! It is **not** parity-dependent, and odd `N` does not remove it — it *strengthens* it. Measured
//! overlap: exactly **25.00%** at every even `N`, **32.65% at N=7**, 30.86% at 9. So an odd `N` is
//! the only available lever that *varies* CRN strength, which is what separates "coarse `N`
//! under-splits" from "parent-child CRN is doing the work". Two effects that would otherwise be
//! confounded across an all-even sweep.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

/// Fraction of a child's footprints that coincide exactly with a parent sample.
fn crn_fraction(n: usize) -> f64 {
    let per_axis = ((n - 1) / 2) + 1;
    (per_axis * per_axis) as f64 / (n * n) as f64
}

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(4000);
    let tau: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1e-4);

    println!("fixed budget {budget} QUADS (not trajectories), tau={tau:.0e}, t=13, f64");
    println!("cost per quad is N^2*(E+1) trajectories, so equal budget is NOT equal compute");
    println!();
    println!("{:>14}{:>4}{:>8}{:>10}{:>9}{:>7}{:>9}{:>9}{:>10}{:>10}",
             "region", "N", "traj/q", "crn frac", "leaves", "depth", "floor", "keep",
             "med alpha", "wall s");

    for region in ["far", "near-field", "deep interior"] {
        for n in [4usize, 7, 8, 16] {
            let root = grid::region(region, 2, 2, 0.05).unwrap();
            let cfg = SchedCfg { n, budget, tau_display: tau, ..Default::default() };
            let (t, st) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
            let mut a: Vec<f64> = t
                .nodes
                .iter()
                .filter_map(|q| q.alpha)
                .filter(|x| x.is_finite())
                .collect();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let med = if a.is_empty() { f64::NAN } else { a[a.len() / 2] };
            println!("{region:>14}{n:>4}{:>8}{:>9.2}%{:>9}{:>7}{:>9}{:>9}{med:>10.4}{:>10.1}",
                     n * n * (ens.n_extra + 1), 100.0 * crn_fraction(n), leaves.len(),
                     t.depth_histogram().len().saturating_sub(1), c(D::Floor), c(D::Keep),
                     st.wall_seconds);
        }
        println!();
    }

    println!("If leaf count and depth track N monotonically, coarse N under-splits and the tier");
    println!("default matters. If N=7 breaks the trend against N=4/8/16, the parent-child CRN is");
    println!("doing work that the even-N sweep alone could not have separated.");
}
