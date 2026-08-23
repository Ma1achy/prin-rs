//! **§4 question 4: does per-quad noise cause thrash?**
//!
//! Chaotic-region `alpha` scatters over 1.1–1.3, so neighbouring quads can get opposite decisions
//! from the same underlying physics. Quantified as: of adjacent leaf pairs whose spreads are
//! **similar**, what fraction sit at **different levels**?
//!
//! **This figure is under-reported, and the amount is known.** `Slice::axis` is endpoint-inclusive,
//! so adjacent siblings share a whole column of footprints — same level, same cell width, same
//! Halton offsets, therefore *identical copies*. That is `1/N` of a quad's data: **12.5% at N=8**,
//! 25% at N=4, 6.25% at N=16. Shared footprints make neighbours more alike than the physics makes
//! them, so measured thrash is a **lower bound**. The overlap fraction is printed beside every
//! figure rather than left implicit.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::QuadTree;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

/// Do two leaf boxes share an edge (touch, not overlap)?
fn adjacent(t: &QuadTree, i: usize, j: usize) -> bool {
    let (a, b) = (&t.nodes[i], &t.nodes[j]);
    let (dx, dy) = ((a.cx - b.cx).abs(), (a.cy - b.cy).abs());
    let (tx, ty) = (a.half + b.half, a.half + b.half);
    let e = 1e-12 * tx;
    // Touching in x and overlapping in y, or the transpose.
    ((dx - tx).abs() < e && dy < ty - e) || ((dy - ty).abs() < e && dx < tx - e)
}

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let tau: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(1e-4);
    let alpha_hi: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.2);

    println!("adjacent leaf pairs, tau={tau:.0e}, alpha_hi={alpha_hi}, budget 4000 quads, t=13, f64");
    println!();
    println!("**A uniform tree cannot thrash.** Where the descent stops at the bootstrap every leaf");
    println!("is the same level, so `diff level` is structurally 0 and a thrash of 0.0000 means the");
    println!("tree never descended - not that neighbours agreed. Read the depth column with it.");
    println!("'similar spread' = within a factor of 1.5 of each other");
    println!();
    println!("{:>14}{:>5}{:>8}{:>7}{:>10}{:>9}{:>10}{:>9}{:>12}",
             "region", "N", "leaves", "depth", "adj pairs", "similar", "diff lvl", "thrash", "edge share");

    for region in ["far", "near-field", "deep interior"] {
        for n in [4usize, 8, 16] {
            let root = grid::region(region, 2, 2, 0.05).unwrap();
            let cfg = SchedCfg { n, budget: 4000, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4, ..Default::default() };
            let (t, _) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();

            let (mut adj, mut similar, mut diff) = (0usize, 0usize, 0usize);
            for (u, &i) in leaves.iter().enumerate() {
                for &j in leaves.iter().skip(u + 1) {
                    if !adjacent(&t, i, j) {
                        continue;
                    }
                    adj += 1;
                    let (si, sj) = (t.nodes[i].red.spread_median, t.nodes[j].red.spread_median);
                    if si > 0.0 && sj > 0.0 && (si / sj).max(sj / si) < 1.5 {
                        similar += 1;
                        if t.nodes[i].level != t.nodes[j].level {
                            diff += 1;
                        }
                    }
                }
            }
            let thrash = if similar > 0 { diff as f64 / similar as f64 } else { f64::NAN };
            println!("{region:>14}{n:>5}{:>8}{:>7}{adj:>10}{similar:>9}{diff:>10}{thrash:>9.4}{:>11.1}%",
                     leaves.len(), t.depth_histogram().len().saturating_sub(1), 100.0 / n as f64);
        }
    }

    println!();
    println!("'thrash' is diff-level over similar-spread adjacent pairs. The last column is the");
    println!("fraction of a quad's footprints it shares identically with each edge neighbour, so");
    println!("every thrash figure is a LOWER BOUND and the bound is loosest at small N.");
}
