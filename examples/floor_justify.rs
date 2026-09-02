//! **Does `Decision::Floor` still fire on a clean substrate, and if not, was it detecting a bug?**
//!
//! `corpus_scope` measured `Floor` collapsing **17 -> 0, 16 -> 1, 9 -> 1** between the corpus
//! kernel and production. `quad.rs` calls `Floor` *"the branch that must be shown to engage"*, so
//! a branch that stops engaging needs re-justifying rather than quietly staying in the descent.
//!
//! # The mechanism, and why it predicts exactly this
//!
//! `Floor` fires at `alpha < alpha_lo`, where `alpha = log2(spread_parent / spread_child)`: refining
//! did not reduce the spread. A spread inflated by **integration failure** is exactly what looks
//! irreducible, because the failure belongs to the trajectory rather than to the cell, so halving
//! the cell does not halve it. Under the old kernel that inflation was everywhere; under the new
//! one it is gone.
//!
//! So the prediction is not merely "`Floor` fires less". It is that the **`alpha` distribution
//! moves up** — away from zero — and that is the discriminating measurement, because it separates
//! two very different readings of the same collapse:
//!
//! - **`alpha` rose**: the old `Floor` firings were the bug. Refinement really does reduce the
//!   spread on a clean substrate, and a predicate that only fired on broken data has no subject.
//! - **`alpha` is unchanged and the threshold merely sits elsewhere**: `Floor` is still measuring
//!   something and the collapse is a scaling artefact.
//!
//! # The guards
//!
//! - **`alpha` must exist to be thresholded.** A quad with no parent, or whose children were never
//!   computed, has `alpha = None` and can never reach the branch. Count those separately: a
//!   `Floor` rate over quads that could not have fired is a rate with the wrong denominator.
//! - **`alpha_lo` is swept.** Reading one threshold cannot tell "the distribution moved" from
//!   "the threshold is in the wrong place", and *`alpha_hi` from 0.20 to 0.50 collapses the tree
//!   80x* is already on record. The sweep is what makes the claim about the distribution.
//! - **Both kernels, one code path.** The only thing that differs between the arms is
//!   `EnsembleCfg`.
use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::ensemble::provenance::Override;
use prin_rs::integrate::az::driver::{DtauMode, StepLimit};
use prin_rs::integrate::Integrator;
use prin_rs::outcome::EscapeRule;
use prin_rs::quad::Decision;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::spatial::HotRule;
use prin_rs::stats;

fn pre_kernel() -> EnsembleCfg {
    let mut c = EnsembleCfg::production().with_overrides(&[
        Override::Integrator(Integrator::Az),
        Override::DtauMode(DtauMode::FixedPerInterval),
        Override::ClampFinalStep(false),
        Override::StepLimit(StepLimit::None),
        Override::EscapeRule(EscapeRule::Reference),
        Override::EscapeConfirm(false),
        Override::RefineFlagged(false),
    ]);
    c.land_iterate = false;
    c
}

fn q(v: &[f64], p: f64) -> f64 {
    let mut f: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if f.is_empty() {
        return f64::NAN;
    }
    f.sort_by(|a, b| a.partial_cmp(b).unwrap());
    f[(((f.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let regions: Vec<String> = {
        let r: Vec<String> = std::env::args().skip(2).collect();
        if r.is_empty() {
            ["near-field", "body2 core", "deep interior", "mid-field", "far"]
                .iter().map(|s| s.to_string()).collect()
        } else { r }
    };
    let pre = pre_kernel();
    let now = EnsembleCfg::production().with_overrides(&[Override::RefineFlagged(false)]);
    println!("DOES `Decision::Floor` STILL FIRE ON A CLEAN SUBSTRATE?  budget {budget} quads\n");
    println!("  pre config: {}", pre.provenance());
    println!("  now config: {}\n", now.provenance());

    // ---- the alpha distribution, which is what decides the reading -------------------------
    println!("THE `alpha` DISTRIBUTION over quads that HAVE one. `alpha = log2(spread_parent /");
    println!("spread_child)`, so `Floor` at `alpha < alpha_lo` fires when refining did not pay.\n");
    println!("{:<14} {:>4} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
             "region", "arm", "n_alpha", "p10", "p25", "p50", "p75", "p90", "frac<0.2");
    let mut trees = Vec::new();
    for r in &regions {
        let Some((_, cx, cy, body)) = prin_rs::grid::REGIONS.iter().find(|x| x.0 == *r).cloned()
        else { continue };
        let cfg = SchedCfg {
            n: 8, bootstrap_levels: 2, budget, tau_display: 1e-4,
            hot_rule: HotRule::Quantile(0.5), alpha_hi: 0.2,
            camera: Some(Camera::framing(cx, cy, 0.05, 1024)),
            ..Default::default()
        };
        let mut per_arm = Vec::new();
        for (tag, ens) in [("pre", &pre), ("now", &now)] {
            let (t, _) = scheduler::descend(cx, cy, 0.05, body, &cfg, ens, Precision::F64);
            let al: Vec<f64> =
                t.nodes.iter().filter_map(|n| n.alpha).filter(|a| a.is_finite()).collect();
            let no_alpha = t.nodes.iter().filter(|n| n.red.n_footprints > 0 && n.alpha.is_none()).count();
            let below = al.iter().filter(|&&a| a < 0.2).count() as f64 / al.len().max(1) as f64;
            println!("{:<14} {tag:>4} {:>7} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {below:>8.3}",
                     if tag == "pre" { r.as_str() } else { "" }, al.len(),
                     q(&al, 0.10), q(&al, 0.25), q(&al, 0.50), q(&al, 0.75), q(&al, 0.90));
            per_arm.push((tag, t, al, no_alpha));
        }
        trees.push((r.clone(), per_arm));
    }

    println!("\n  QUADS WITH NO `alpha` AT ALL -- no parent, or children never computed. They can");
    println!("  never reach the branch, so a Floor rate including them has the wrong denominator:");
    for (r, arms) in &trees {
        let s: Vec<String> =
            arms.iter().map(|(t, tr, al, na)| {
                format!("{t} {na}/{} computed, {} with alpha",
                        tr.nodes.iter().filter(|n| n.red.n_footprints > 0).count(), al.len())
            }).collect();
        println!("    {r:<14} {}", s.join("   |   "));
    }

    // ---- the sweep: distribution moved, or threshold in the wrong place? --------------------
    println!("\nTHE `alpha_lo` SWEEP. Reading one threshold cannot tell \"the distribution moved\"");
    println!("from \"the threshold is in the wrong place\". Leaf counts stopped by `Floor`:\n");
    let ladder = [0.05f64, 0.1, 0.2, 0.35, 0.5, 0.8];
    print!("{:<14} {:>4}", "region", "arm");
    for a in ladder { print!("{:>10}", format!("lo={a}")); }
    println!("{:>10}", "leaves");
    for r in &regions {
        let Some((_, cx, cy, body)) = prin_rs::grid::REGIONS.iter().find(|x| x.0 == *r).cloned()
        else { continue };
        for (tag, ens) in [("pre", &pre), ("now", &now)] {
            print!("{:<14} {tag:>4}", if tag == "pre" { r.as_str() } else { "" });
            let mut leaves = 0;
            for a_lo in ladder {
                let cfg = SchedCfg {
                    n: 8, bootstrap_levels: 2, budget, tau_display: 1e-4,
                    hot_rule: HotRule::Quantile(0.5), alpha_hi: 0.2, alpha_lo: a_lo,
                    camera: Some(Camera::framing(cx, cy, 0.05, 1024)),
                    ..Default::default()
                };
                let (t, _) = scheduler::descend(cx, cy, 0.05, body, &cfg, ens, Precision::F64);
                let f = t.nodes.iter().filter(|n| n.is_leaf() && n.decision == Decision::Floor).count();
                leaves = t.leaves().count();
                print!("{f:>10}");
            }
            println!("{leaves:>10}");
        }
    }

    println!("\nHOW TO READ IT");
    println!("  If `alpha` ROSE between the arms, the old Floor firings were the integration bug:");
    println!("  refinement really does reduce the spread on a clean substrate, and a predicate that");
    println!("  only ever fired on broken data has no subject. If `alpha` is UNCHANGED and only the");
    println!("  count moved, Floor still measures something and `alpha_lo` is merely mis-placed.");
    println!("  The sweep is what separates them: a distribution that moved shows Floor firing again");
    println!("  at a HIGHER `alpha_lo` under `now`; one that did not shows the same curve shifted.");
    let _ = stats::quantile(&mut vec![0.0], 0.5);
}
