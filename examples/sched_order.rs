//! **§4 question 6: does the queue order matter?**
//!
//! Priority-ordered against shuffled, same budget. If the trees come out the same, ordering is
//! not doing any work and the priority function is a free parameter that could be dropped; if
//! they differ, the ordering is part of the criterion's behaviour and has to be reported as such.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::QuadTree;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Order, SchedCfg};

fn leaf_set(t: &QuadTree) -> std::collections::HashSet<(u32, String)> {
    t.leaves()
        .map(|i| (t.nodes[i].level, format!("{:.12e},{:.12e}", t.nodes[i].cx, t.nodes[i].cy)))
        .collect()
}

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(10_000);
    let tau: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1e-4);
    let alpha_hi: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.2);

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, t=13, f64");
    println!();
    println!("**Order can only matter when the budget BINDS.** If the descent terminates before the");
    println!("cap, every ordering computes the same set and a jaccard of 1.0000 is structural, not");
    println!("a result. The budget-exhausted column says whether this run could see anything.");
    println!();
    println!("{:>14}{:>14}{:>9}{:>9}{:>7}{:>10}{:>10}{:>10}",
             "region", "order", "quads", "leaves", "depth", "cap hit", "vs spread", "jaccard");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        let mut base: Option<std::collections::HashSet<(u32, String)>> = None;
        for (k, order) in [Order::Spread, Order::SpreadArea, Order::Shuffled].into_iter().enumerate()
        {
            let cfg = SchedCfg {
                budget,
                tau_display: tau,
                alpha_hi,
                alpha_lo: alpha_hi * 0.4,
                order,
                seed: 12345,
                ..Default::default()
            };
            let (t, st) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let s = leaf_set(&t);
            let (diff, jac) = match &base {
                None => (0usize, 1.0f64),
                Some(b) => (
                    b.symmetric_difference(&s).count(),
                    b.intersection(&s).count() as f64 / b.union(&s).count().max(1) as f64,
                ),
            };
            println!("{region:>14}{:>14}{:>9}{:>9}{:>7}{:>10}{diff:>10}{jac:>10.4}",
                     order.name(), st.quads_computed, s.len(),
                     t.depth_histogram().len().saturating_sub(1), st.budget_exhausted);
            if k == 0 {
                base = Some(s);
            }
        }
    }

    println!();
    println!("'vs spread' is the symmetric difference of the leaf sets against the spread-ordered");
    println!("run; jaccard is their overlap. 1.0000 means the order made no difference at all.");
}
