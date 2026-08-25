//! §8 — two costings, and neither is an implementation.
//!
//! # Cost-aware priority
//!
//! Some quads are far more expensive than others: close encounters take more substeps. Ranking
//! by `spread / compute_cost` spends a budget better **if costs vary widely**. So the cost
//! distribution is reported first: if it is narrow, the idea is moot and stops here. `Rank`
//! carries a `GreedyOraclePerCost` variant so the ceiling can be read per unit cost too.
//!
//! # Anisotropic splitting
//!
//! A boundary running diagonally through a quad gets four children, three of them mostly
//! wasted. Splitting two ways along the disagreement direction is strictly cheaper. **Costing
//! only**: how many splits produce children that would immediately `keep`? A large fraction
//! makes it worth building; a small one does not.
//!
//! # How to misread this
//!
//! **`p99/p50` is the number that decides cost-aware priority, not the mean.** A cost
//! distribution with a long tail and a tight bulk buys nothing from cost-weighting for the
//! typical quad, which is what a scheduler is deciding about.
//!
//! **"Children that would immediately keep" is measured against a threshold**, so it inherits
//! `tau`'s region-dependence. It is reported at three `tau` values for that reason: a single
//! value would make the anisotropy case look settled when it is a function of a knob the
//! vertical slice already showed moves trees 64x.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Colouring, Rank};
use prin_rs::quad::{Agg, Criterion};
use prin_rs::stats;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let levels: u32 = arg(1, 5);
    let n: usize = arg(2, 8);
    let res = (1usize << levels) * n;
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("level {levels}, N={n}, res {res}^2, {full} quads per region, t=13\n");

    for &(region, cx, cy, body) in grid::REGIONS
        .iter()
        .filter(|r| matches!(r.0, "near-field" | "deep interior" | "far"))
    {
        let cache = metric::build(
            region, cx, cy, 0.05, body, Chart::BodyPlane, levels, n, res, 1e-4, &ens,
            Colouring::Outcome,
        );

        // ---- cost distribution ----
        let cost: Vec<f64> = cache.quads.values().map(|q| q.red.total_substeps as f64).collect();
        let (p10, p50, p90, _) = stats::interdecile(&cost);
        let mut sorted = cost.clone();
        let p99 = prin_rs::quad::quantile(&mut sorted, 0.99);
        let mx = cost.iter().cloned().fold(0.0f64, f64::max);
        println!(
            "--- {region} --- substeps per quad: p10 {p10:.0}  p50 {p50:.0}  p90 {p90:.0}  \
             p99 {p99:.0}  max {mx:.0}"
        );
        println!(
            "    p90/p50 = {:.2}   p99/p50 = {:.2}   max/p50 = {:.2}",
            p90 / p50,
            p99 / p50,
            mx / p50
        );

        // Does cost-weighting change the ceiling? If the cost distribution is narrow this is
        // the same curve twice.
        let budgets = [21usize, 85, 341, full];
        for r in [Rank::GreedyOracle, Rank::GreedyOraclePerCost] {
            let pts = metric::replay(&cache, r, full);
            print!("{:>24}", r.name());
            for e in metric::curve_at(&pts, &budgets) {
                print!(" {e:>9.5}");
            }
            println!();
        }

        // ---- anisotropy costing ----
        //
        // A split is "wasted" in proportion to how many of its four children would immediately
        // keep. Measured against tau, and therefore reported across tau rather than at one.
        for tau in [1e-5, 1e-4, 1e-3] {
            let (mut splits, mut wasted3, mut wasted4) = (0usize, 0usize, 0usize);
            for q in cache.quads.values() {
                if q.key.0 >= cache.levels {
                    continue;
                }
                let keeps = metric::Cache::children(q.key)
                    .iter()
                    .filter(|&&c| {
                        let v = cache.get(c).red.signal(Criterion::Within, Agg::Median);
                        !(v > tau)
                    })
                    .count();
                splits += 1;
                if keeps >= 3 {
                    wasted3 += 1;
                }
                if keeps == 4 {
                    wasted4 += 1;
                }
            }
            println!(
                "    tau={tau:e}: {:.1}% of splits give >=3 children that immediately keep, \
                 {:.1}% give all four",
                100.0 * wasted3 as f64 / splits.max(1) as f64,
                100.0 * wasted4 as f64 / splits.max(1) as f64
            );
        }
        println!();
    }

    println!(
        "Cost-aware priority is worth building only if the cost distribution is WIDE. Read\n\
         p99/p50, not the mean: a long tail with a tight bulk buys nothing for the typical quad,\n\
         which is what a scheduler decides about. If greedy_oracle and greedy_oracle/cost give\n\
         the same curve, cost-weighting has nothing to move.\n\
         \n\
         Anisotropy: a high `all four keep` fraction means an isotropic split is mostly wasted\n\
         and a 2-way split along the disagreement direction -- which §3.2's layout already\n\
         measures -- would be strictly cheaper. This is a COSTING. Nothing anisotropic is\n\
         implemented, and the tau sweep is there because a single tau would make the case look\n\
         settled when it is a function of a knob that moves trees 64x."
    );
}
