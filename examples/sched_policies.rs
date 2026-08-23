//! **§4 question 5: does the sibling-spread policy beat the alpha policy?** Equal budget, both.
//!
//! Separation in `alpha`'s **value** is 0.9862 against a chaotic scatter of 1.1–1.3 — marginal.
//! Separation in `alpha`'s **reliability** is 0.001 against 1.2 — three orders. A split produces
//! four children and therefore four exponents; their spread is a per-quad reliability estimate at
//! no extra cost, and it is currently discarded.
//!
//! Where the four scatter, the unreliability *is* the answer: floor, without needing a trustworthy
//! `alpha` at all. That removes the awkwardness of thresholding a quantity whose noise is
//! comparable to its range.
//!
//! **Caveat carried forward:** the range of four samples is itself a noisy statistic. If this
//! policy looks promising, that noise is the next thing to characterise rather than the first
//! thing to trust.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::{Decision as D, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Policy, SchedCfg};

/// How much of the budget landed where the spread is high — the thing a scheduler is for.
fn spend_profile(t: &QuadTree) -> (f64, f64, f64) {
    let mut v: Vec<(f64, f64)> = t
        .nodes
        .iter()
        .filter(|q| q.red.n_footprints > 0)
        .map(|q| (q.red.spread_median, 1.0))
        .collect();
    if v.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let total: f64 = v.iter().map(|x| x.1).sum();
    let top = |f: f64| -> f64 {
        let k = ((v.len() as f64) * f).ceil() as usize;
        v.iter().take(k).map(|x| x.1).sum::<f64>() / total
    };
    // Fraction of quads spent in the top decile of spread, and the spread at the median quad.
    let med = v[v.len() / 2].0;
    (top(0.1), med, v[0].0)
}

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(10_000);
    let tau: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1e-4);
    let alpha_hi: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.2);

    println!("equal budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, t=13, f64");
    println!("A comparison at a threshold where neither policy descends is not a comparison, so");
    println!("alpha_hi comes from the sweep.");
    println!();
    println!("{:>14}{:>10}{:>9}{:>9}{:>8}{:>8}{:>8}{:>7}{:>12}{:>12}",
             "region", "policy", "quads", "leaves", "split", "floor", "keep", "depth",
             "med spread", "max spread");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        let mut trees = Vec::new();
        for policy in [Policy::Alpha, Policy::Sibling] {
            let cfg = SchedCfg { budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4, policy, ..Default::default() };
            let (t, st) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
            let splits = t.nodes.iter().filter(|q| q.decision == D::Split).count();
            let (_, med, max) = spend_profile(&t);
            println!("{region:>14}{:>10}{:>9}{:>9}{splits:>8}{:>8}{:>8}{:>7}{med:>12.4e}{max:>12.4e}",
                     policy.name(), st.quads_computed, leaves.len(), c(D::Floor), c(D::Keep),
                     t.depth_histogram().len().saturating_sub(1));
            trees.push(t);
        }
        // Where the two trees actually differ: leaves present in one and not the other.
        let leaf_key = |t: &QuadTree, i: usize| {
            (t.nodes[i].level, format!("{:.12e},{:.12e}", t.nodes[i].cx, t.nodes[i].cy))
        };
        let a: std::collections::HashSet<_> =
            trees[0].leaves().map(|i| leaf_key(&trees[0], i)).collect();
        let b: std::collections::HashSet<_> =
            trees[1].leaves().map(|i| leaf_key(&trees[1], i)).collect();
        println!("{:>14}  shared leaves {}, alpha-only {}, sibling-only {}",
                 "", a.intersection(&b).count(), a.difference(&b).count(), b.difference(&a).count());
    }

    println!();
    println!("Equal budget, so the question is where each SPENT it, not how much it used. A policy");
    println!("that floors the chaotic sea earlier has budget left for structure; one that does not");
    println!("starves the structure to feed the sea.");
}
