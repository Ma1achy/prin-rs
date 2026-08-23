//! **§4 questions 1 and 2: does the descent terminate, and does the floor engage?**
//!
//! Question 1 is the one that could invalidate the design. With Wada-like boundaries dense at
//! every scale, the spread may stay high however far you refine, and the descent may never stop.
//! It was flagged at the outset of this work and never tested.
//!
//! **If it does not terminate, that is the headline result** — where it kept splitting and what
//! `alpha` was doing there, not a cap and a shrug.
//!
//! Run with **no `max_level`**, so only the budget or the criterion can stop it. A descent that
//! stops at the precision floor (level ~36) has hit a *numerical* limit, not a physical one, and
//! must be reported as such.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn main() {
    let budget: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(50_000);
    let tau: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1e-4);

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    println!("budget {budget} quads, NO max_level, tau={tau:.0e}, alpha_hi=0.5, N=8, t=13, f64");
    println!();

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        let cfg = SchedCfg { budget, tau_display: tau, max_level: None, ..Default::default() };
        let (t, st) =
            scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);

        let leaves: Vec<usize> = t.leaves().collect();
        let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
        // Over all quads: a leaf can never carry Split, since splitting makes it internal.
        let splits = t.nodes.iter().filter(|q| q.decision == D::Split).count();
        let hist = t.depth_histogram();

        println!("=== {region} ===");
        println!("  {} quads computed of {budget}, {} leaves, {} iterations, {:.1} s",
                 st.quads_computed, leaves.len(), st.iterations, st.wall_seconds);
        println!("  budget exhausted: {}", st.budget_exhausted);
        println!("  {splits} quads decided to split (over all quads; a leaf never carries Split)");
        println!("  leaf decisions: floor {} keep {} precision_floor {} budget {}",
                 c(D::Floor), c(D::Keep), c(D::PrecisionFloor), c(D::BudgetExhausted));

        let terminal = c(D::Floor) + c(D::Keep);
        let frac_floor = c(D::Floor) as f64 / leaves.len().max(1) as f64;
        println!("  TERMINATED of its own accord: {terminal}/{} leaves ({:.1}%)",
                 leaves.len(), 100.0 * terminal as f64 / leaves.len().max(1) as f64);
        println!("  floored specifically:         {:.1}%  <- if ~0, the floor branch is not working",
                 100.0 * frac_floor);

        print!("  depth histogram:");
        for (l, n) in hist.iter().enumerate() {
            if *n > 0 {
                print!(" {l}:{n}");
            }
        }
        println!();
        print!("  leaves per iteration:");
        for n in st.leaves_per_iteration.iter().take(40) {
            print!(" {n}");
        }
        if st.leaves_per_iteration.len() > 40 {
            print!(" ...");
        }
        println!();

        // Where did it keep splitting, and what was alpha doing there?
        let deepest = hist.len().saturating_sub(1);
        let deep: Vec<usize> = leaves
            .iter()
            .cloned()
            .filter(|&i| t.nodes[i].level as usize + 2 >= deepest)
            .collect();
        if !deep.is_empty() && deepest > 3 {
            let mut a: Vec<f64> =
                deep.iter().filter_map(|&i| t.nodes[i].alpha).filter(|x| x.is_finite()).collect();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let mut s: Vec<f64> = deep.iter().map(|&i| t.nodes[i].red.spread_median).collect();
            s.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let q = |v: &Vec<f64>, f: f64| v[(((v.len() - 1) as f64) * f).round() as usize];
            println!("  in the deepest two levels ({} leaves):", deep.len());
            if !a.is_empty() {
                println!("    alpha  p10 {:+.3}  median {:+.3}  p90 {:+.3}",
                         q(&a, 0.1), q(&a, 0.5), q(&a, 0.9));
            }
            println!("    spread p10 {:.3e}  median {:.3e}  p90 {:.3e}",
                     q(&s, 0.1), q(&s, 0.5), q(&s, 0.9));
        }
        if c(D::PrecisionFloor) > 0 {
            println!("  ** {} leaves hit the PRECISION floor (level ~36). That is a numerical",
                     c(D::PrecisionFloor));
            println!("     limit, not a physical one, and is not evidence about termination. **");
        }
        println!();
    }

    println!("A descent that exhausts the budget with most leaves still wanting to split has NOT");
    println!("terminated, and no cap makes that a different answer.");
}
