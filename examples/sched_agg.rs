//! **§3.4: which aggregation, and how much does the choice change the decisions?**
//!
//! A quad holds `N²` footprints, each with its own `ensemble_spread`; `QuadReduction` needs one
//! number. With excess kurtosis 110 a mean is dominated by a single footprint — but a **median is
//! blind to a thin filament crossing a quad**, because most of that quad's footprints are still in
//! the smooth sea.
//!
//! That is not hypothetical here. The spread overlay shows the descent refining coherent bands
//! while leaving the brightest, thinnest spread filaments in coarse quads. This measures it.
//!
//! The brief asks for the **decision-level** disagreement, not just the three numbers.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::{Agg, Decision as D, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn leaf_set(t: &QuadTree) -> std::collections::HashSet<(u32, String)> {
    t.leaves()
        .map(|i| (t.nodes[i].level, format!("{:.12e},{:.12e}", t.nodes[i].cx, t.nodes[i].cy)))
        .collect()
}

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(6000);
    let tau: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1e-4);
    let alpha_hi: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.2);

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, t=13, f64");
    println!();
    println!("{:>14}{:>9}{:>9}{:>9}{:>7}{:>8}{:>8}{:>10}{:>12}",
             "region", "agg", "quads", "leaves", "depth", "floor", "keep", "cap hit", "vs median");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        let mut base: Option<std::collections::HashSet<(u32, String)>> = None;
        let mut rows = Vec::new();

        for agg in [Agg::Median, Agg::Mean, Agg::P90] {
            let cfg = SchedCfg {
                budget,
                tau_display: tau,
                alpha_hi,
                alpha_lo: alpha_hi * 0.4,
                agg,
                ..Default::default()
            };
            let (t, st) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
            let s = leaf_set(&t);
            let jac = match &base {
                None => 1.0,
                Some(b) => b.intersection(&s).count() as f64 / b.union(&s).count().max(1) as f64,
            };
            println!("{region:>14}{:>9}{:>9}{:>9}{:>7}{:>8}{:>8}{:>10}{jac:>12.4}",
                     agg.name(), st.quads_computed, leaves.len(),
                     t.depth_histogram().len().saturating_sub(1),
                     c(D::Floor), c(D::Keep), st.budget_exhausted);
            if base.is_none() {
                base = Some(s);
            }
            rows.push((agg, t));
        }

        // Decision-level disagreement: over quads present in BOTH trees, how many decided
        // differently? A leaf-set jaccard says the trees differ; this says the DECISIONS do.
        for k in 1..rows.len() {
            let (a, b) = (&rows[0].1, &rows[k].1);
            let map_a: std::collections::HashMap<(u32, String), D> = a
                .nodes
                .iter()
                .map(|q| ((q.level, format!("{:.12e},{:.12e}", q.cx, q.cy)), q.decision))
                .collect();
            let (mut common, mut differ) = (0usize, 0usize);
            for q in &b.nodes {
                let key = (q.level, format!("{:.12e},{:.12e}", q.cx, q.cy));
                if let Some(&d) = map_a.get(&key) {
                    common += 1;
                    if d != q.decision {
                        differ += 1;
                    }
                }
            }
            println!("{:>14}  {} vs median: {} of {} shared quads decided differently ({:.1}%)",
                     "", rows[k].0.name(), differ, common,
                     100.0 * differ as f64 / common.max(1) as f64);
        }
        println!();
    }

    println!("A median is blind to a thin filament crossing a quad; a p90 is not, and a mean sits");
    println!("between while being hostage to one footprint. If the decisions move, the aggregation");
    println!("is a real parameter and cannot be picked silently.");
}
