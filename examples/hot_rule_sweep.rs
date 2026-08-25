//! **The hot rule, swept — and what it can and cannot move.**
//!
//! The mask threshold feeds the *shape* statistics (`n_components`, `largest_component`,
//! `perimeter_ratio`) and the criteria built on them. Under the shipped absolute rule with
//! `tau = 1e-4` those are saturated: pooled over the 89,088 committed leaves, `n_hot == N^2` in
//! 90.2% and `n_components <= 1` in 95.0%. One blob covering the whole quad, everywhere.
//!
//! # What this sweep does NOT show, stated first
//!
//! **The tree is identical across every row**, and that is by construction, not a result. The
//! split decision reads `cfg.criterion`, which is `Within` here; the hot rule only reaches the
//! decision when the criterion is one of the mask-derived ones. A leaf count that did not move
//! would otherwise read as "the hot rule is inert", which is a measurement of the wiring.
//!
//! So the tree column is printed as the **control** — it must be constant — and the measurement
//! is the mask, plus the distinct-value count of each criterion the mask feeds.
//!
//! # Two things this measured that were not expected
//!
//! **The saturation is not uniform across regions, and in `far` the absolute mask is EMPTY, not
//! full.** `far`'s leaf-spread median is `4.26e-8` against `tau = 1e-4`, so nothing clears the
//! cut: `n_hot == 0`, `perimeter_ratio` is `NaN` by the empty-set convention, and every criterion
//! built on it takes **one** distinct value over all 16 leaves. `deep interior` already resolves
//! a median of 2 components under the absolute rule. It is near-field and the latent charts where
//! the mask is full. "Saturated everywhere" is the pooled number, not the regional one, and the
//! two failure modes — full mask and empty mask — are the same threshold landing on either side
//! of the distribution.
//!
//! **The relative rule desaturates the mask and COARSENS the ordering.** In near-field the median
//! component count runs 1 -> 5 as the rule goes absolute -> `q[0.50]`, which is the mask finally
//! describing something; but `Criterion::LayoutRel`'s distinct-value count falls 78 -> 26 -> 17 ->
//! 9 across `abs, q[0.50], q[0.75], q[0.90]` while `Criterion::Layout` holds at 58. With `n_hot`
//! pinned by the rule, `largest/n_hot` can only take as many values as there are component sizes.
//!
//! That is reported, not hidden, and it does not settle anything by itself: the standing result
//! is that **signal resolution is not what makes a ranking good** — `frac_hot_between` is the
//! best criterion measured here on 65 distinct values, beating a 4994-valued one. `error(B)`
//! decides. But a criterion whose ordering coarsens as its input improves is worth watching.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::{Agg, Criterion};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::spatial::HotRule;
use prin_rs::stats;

const BUDGET: usize = 600;
const TAU: f64 = 1e-4;

fn main() {
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let rules = [
        HotRule::AbsTau(TAU),
        HotRule::Quantile(0.5),
        HotRule::Quantile(0.75),
        HotRule::Quantile(0.9),
    ];

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        println!("=== {region} ===  budget {BUDGET}, N=8, E+1=8, t=13, tau={TAU:.0e}, criterion=within");
        println!("{:>12}{:>8}{:>10}{:>10}{:>10}{:>11}{:>11}{:>11}{:>10}",
                 "hot rule", "leaves", "sat%", "1comp%", "med comp", "perim p10", "perim p90",
                 "d(layout)", "d(l_rel)");

        let mut base_leaves = None;
        for rule in rules {
            let cfg = SchedCfg {
                budget: BUDGET,
                tau_display: TAU,
                hot_rule: rule,
                alpha_hi: 0.2,
                alpha_lo: 0.2,
                ..Default::default()
            };
            let (t, _) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let m = leaves.len() as f64;

            // Which layout the rule actually moved: absolute rows read `layout_within`,
            // quantile rows read `layout_rel_within`. Both are always computed.
            let lay = |i: usize| match rule {
                HotRule::AbsTau(_) => t.nodes[i].red.layout_within,
                HotRule::Quantile(_) => t.nodes[i].red.layout_rel_within,
            };
            let full = (cfg.n * cfg.n) as u32;
            let sat = leaves.iter().filter(|&&i| lay(i).n_hot == full).count() as f64 / m;
            let one = leaves.iter().filter(|&&i| lay(i).n_components <= 1).count() as f64 / m;
            let mut comps: Vec<f64> = leaves.iter().map(|&i| lay(i).n_components as f64).collect();
            let mut per: Vec<f64> =
                leaves.iter().map(|&i| lay(i).perimeter_ratio).filter(|x| x.is_finite()).collect();
            let (p10, _, p90, _) = stats::interdecile(&per);
            let med = stats::quantile(&mut comps, 0.5);
            let _ = &mut per;

            let distinct = |c: Criterion| {
                let mut v: Vec<u64> = leaves
                    .iter()
                    .map(|&i| t.nodes[i].red.signal(c, Agg::Median).to_bits())
                    .collect();
                v.sort_unstable();
                v.dedup();
                v.len()
            };

            println!("{:>12}{:>8}{:>10.1}{:>10.1}{:>10.1}{:>11.3}{:>11.3}{:>11}{:>10}",
                     rule.name(), leaves.len(), 100.0 * sat, 100.0 * one, med, p10, p90,
                     distinct(Criterion::Layout), distinct(Criterion::LayoutRel));

            match base_leaves {
                None => base_leaves = Some(leaves.len()),
                Some(b) => assert_eq!(
                    b,
                    leaves.len(),
                    "the tree moved with the hot rule while criterion=within -- the control failed, \
                     so the hot rule has reached the decision path and every row above is confounded"
                ),
            }
        }
        println!();
    }

    println!("sat%     : leaves whose mask covers every footprint. The saturation being corrected.");
    println!("1comp%   : leaves the mask reads as a single blob.");
    println!("d(layout): distinct values of Criterion::Layout over the leaves -- the ordering's");
    println!("           resolution. A criterion with few distinct values has no ordering to offer,");
    println!("           whatever its error(B) curve does; that is a property to read, not to hide.");
    println!();
    println!("The leaf count is CONSTANT down each block, asserted. criterion=within, so the hot");
    println!("rule cannot reach the split decision; a moving tree here would mean the wiring leaked.");
}
