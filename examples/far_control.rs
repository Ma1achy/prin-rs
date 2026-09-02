//! **Is `far` still the control, and what is the tail?**
//!
//! `far` is on record as *what a featureless field looks like*: thirteen non-greedy criteria
//! produce the **identical leaf set**, `err_sum` is flat through level 3, and ranking by spread on
//! a smooth field *is* breadth-first — so its degenerating is correct behaviour rather than a
//! region where the criteria mysteriously tie.
//!
//! `corpus_scope` then found that under the fixed integrator `far` gains a **two-decade tail above
//! p90** — p99 `1.917e-8 -> 1.926e-6`, max `1.921e-8 -> 4.403e-6` — while its bulk stays flat to
//! 1.4%. A tail is exactly what a criterion ranks on. So the description was written against a
//! field that no longer exists, and there are two possible endings: the tail is characterised and
//! `far` stays the control, or it does not and something else has to be.
//!
//! # Three questions, in the order that makes the third answerable
//!
//! **1. Where is the tail, and is it real?** *Amplitude cannot tell a small real signal from
//! noise; coherence can.* A scatter of isolated pixels is round-off; a coherent region is physics.
//! Component sizes and lag-1 neighbour correlation decide it — and *read a mask's component sizes
//! before measuring anything about it*, which caught a dust mask once already.
//!
//! **2. What is it made of?** Cross-tabulated against `energy_drift_max`, `t_end`, terminal state
//! and `d_min`, and split by **which arm** dominates: `ensemble_spread` is
//! `max(spread_shape, spread_event)`, and the event arm has five distinct values and dominates a
//! small tail by construction, so a tail that is entirely event arm is a staircase and not a
//! field.
//!
//! **3. Does the criterion still degenerate there?** The operative question, and the only one that
//! decides whether `far` is still the control. Run the descent under every `Criterion` and compare
//! **leaf sets**, not error numbers.
//!
//! **`tau` must sit below `far`'s bulk or this measures nothing.** At `tau = 1e-4` `far` reads 16
//! leaves, `keep:16`, the standing below-the-bulk degeneracy — every criterion "agrees" because
//! none of them ever runs. *Which rung of a sweep is degenerate is a fact about the region*, and
//! `far`'s median is `4.26e-8`, so the ladder has to reach `1e-8`.
use rayon::prelude::*;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::integrate::az::driver::{DtauMode, StepLimit};
use prin_rs::integrate::Integrator;
use prin_rs::outcome::EscapeRule;
use prin_rs::quad::Criterion;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::spatial::HotRule;

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
    if f.is_empty() { return f64::NAN; }
    f.sort_by(|a, b| a.partial_cmp(b).unwrap());
    f[(((f.len() - 1) as f64) * p).round() as usize]
}

/// Connected components of a boolean mask, 4-connected, largest first.
fn components(m: &[bool], n: usize) -> Vec<usize> {
    let mut seen = vec![false; m.len()];
    let mut out = Vec::new();
    for s in 0..m.len() {
        if !m[s] || seen[s] { continue; }
        let (mut stack, mut sz) = (vec![s], 0usize);
        seen[s] = true;
        while let Some(k) = stack.pop() {
            sz += 1;
            let (x, y) = (k % n, k / n);
            for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (nx, ny) = (x as i64 + dx, y as i64 + dy);
                if nx < 0 || ny < 0 || nx >= n as i64 || ny >= n as i64 { continue; }
                let j = ny as usize * n + nx as usize;
                if m[j] && !seen[j] { seen[j] = true; stack.push(j); }
            }
        }
        out.push(sz);
    }
    out.sort_unstable_by(|a, b| b.cmp(a));
    out
}

/// Lag-1 neighbour correlation of `log10(v)` over the finite pairs, both axes pooled.
fn coherence(v: &[f64], n: usize) -> f64 {
    let l: Vec<f64> = v.iter().map(|x| if *x > 0.0 { x.log10() } else { f64::NAN }).collect();
    let mut p = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            if x + 1 < n { p.push((l[k], l[k + 1])); }
            if y + 1 < n { p.push((l[k], l[k + n])); }
        }
    }
    p.retain(|(a, b)| a.is_finite() && b.is_finite());
    if p.len() < 3 { return f64::NAN; }
    let m = |f: &dyn Fn(&(f64, f64)) -> f64| p.iter().map(f).sum::<f64>() / p.len() as f64;
    let (ma, mb) = (m(&|t: &(f64, f64)| t.0), m(&|t: &(f64, f64)| t.1));
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for (a, b) in &p {
        num += (a - ma) * (b - mb);
        da += (a - ma) * (a - ma);
        db += (b - mb) * (b - mb);
    }
    num / (da * db).sqrt()
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(128);
    let budget: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let pre = pre_kernel();
    let now = EnsembleCfg::production().with_overrides(&[Override::RefineFlagged(false)]);
    println!("IS `far` STILL THE CONTROL?  {n}^2 field, budget {budget} quads\n");
    println!("  pre config: {}", pre.provenance());
    println!("  now config: {}\n", now.provenance());

    let sl = prin_rs::grid::region("far", n, n, 0.05).unwrap();
    let ev = |c: &EnsembleCfg| -> Vec<PixelOut> {
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, c)).collect()
    };
    let (a, b) = (ev(&pre), ev(&now));
    let sa: Vec<f64> = a.iter().map(|p| p.ensemble_spread).collect();
    let sb: Vec<f64> = b.iter().map(|p| p.ensemble_spread).collect();

    // ---- 1. the ladder, and where the tail begins ------------------------------------------
    println!("1. THE LADDER. The bulk is the standing description; the tail is what is new.\n");
    println!("   {:>4} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
             "arm", "p01", "p50", "p90", "p99", "p999", "max", "p99/p01");
    for (t, v) in [("pre", &sa), ("now", &sb)] {
        println!("   {t:>4} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e}",
                 q(v, 0.01), q(v, 0.50), q(v, 0.90), q(v, 0.99), q(v, 0.999), q(v, 1.0),
                 q(v, 0.99) / q(v, 0.01));
    }

    // The tail, defined against the OLD field's own maximum: everything the pre-fix kernel could
    // not produce at all. A threshold picked from the new field would be fitted to it.
    let cut = q(&sa, 1.0);
    let mask: Vec<bool> = sb.iter().map(|&x| x > cut).collect();
    let n_tail = mask.iter().filter(|&&m| m).count();
    println!("\n2. WHERE IT IS. The cut is the PRE arm's own maximum ({cut:.3e}) -- everything the");
    println!("   old kernel could not produce at all, so the threshold is not fitted to the new field.\n");
    let comp = components(&mask, n);
    println!("   tail pixels {n_tail}/{} = {:.4}", sb.len(), n_tail as f64 / sb.len() as f64);
    println!("   components {}  largest {:?}", comp.len(), &comp[..comp.len().min(6)]);
    println!("   lag-1 coherence of log10(spread):  pre {:.4}   now {:.4}",
             coherence(&sa, n), coherence(&sb, n));
    let verdict = if comp.is_empty() {
        "NO TAIL at this resolution -- the corpus_scope finding does not reproduce here"
    } else if comp[0] as f64 / n_tail.max(1) as f64 > 0.5 {
        "COHERENT -- one component holds over half the tail, so this is structure not round-off"
    } else if comp.len() as f64 / n_tail.max(1) as f64 > 0.5 {
        "DUST -- mostly isolated pixels; component sizes say this is not a region"
    } else {
        "MIXED -- read the component list rather than this line"
    };
    println!("   -> {verdict}");

    // ---- 3. what is it made of --------------------------------------------------------------
    println!("\n3. WHAT IT IS MADE OF. `ensemble_spread = max(spread_shape, spread_event)`, and the");
    println!("   event arm has five distinct values -- a tail that is all event arm is a staircase.\n");
    let idx: Vec<usize> = (0..sb.len()).filter(|&k| mask[k]).collect();
    let rest: Vec<usize> = (0..sb.len()).filter(|&k| !mask[k]).collect();
    println!("   {:>10} {:>7} {:>11} {:>11} {:>11} {:>11} {:>9} {:>9}",
             "population", "n", "spread p50", "drift p50", "t_end p50", "d_min p50", "ev doms", "term");
    for (t, ix) in [("tail", &idx), ("bulk", &rest)] {
        if ix.is_empty() { println!("   {t:>10}       0   -- empty --"); continue; }
        let g = |f: &dyn Fn(&PixelOut) -> f64| -> Vec<f64> { ix.iter().map(|&k| f(&b[k])).collect() };
        let ev_dom = ix.iter().filter(|&&k| b[k].spread_event >= b[k].spread_shape).count();
        let term = ix.iter().filter(|&&k| !b[k].censored).count();
        println!("   {t:>10} {:>7} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>9.3} {:>9.3}",
                 ix.len(),
                 q(&g(&|p| p.ensemble_spread), 0.5), q(&g(&|p| p.energy_drift_max), 0.5),
                 q(&g(&|p| p.t_end), 0.5), q(&g(&|p| p.d_min_true), 0.5),
                 ev_dom as f64 / ix.len() as f64, term as f64 / ix.len() as f64);
    }

    // ---- 4. does the criterion still degenerate ---------------------------------------------
    println!("\n4. THE OPERATIVE QUESTION: do the criteria still produce the IDENTICAL leaf set?");
    println!("   `tau` must sit below the bulk or every criterion ties for want of ever running --");
    println!("   `far`'s median is ~4e-8, so the ladder reaches 1e-8.\n");
    let crits = [
        ("within", Criterion::Within), ("between", Criterion::Between),
        ("maxboth", Criterion::MaxOfBoth), ("frac_hot_w", Criterion::FracHotWithin),
        ("frac_hot_b", Criterion::FracHotBetween), ("layout", Criterion::Layout),
        ("layout_rel", Criterion::LayoutRel), ("grad_rms", Criterion::GradRms),
        ("running_max", Criterion::RunningMax), ("first_div", Criterion::FirstDivergence),
        ("term_grad", Criterion::TerminationGradient),
    ];
    let (_, cx, cy, body) = *prin_rs::grid::REGIONS.iter().find(|r| r.0 == "far").unwrap();
    for tau in [1e-4f64, 1e-6, 1e-8, 1e-9] {
        let mut sets: std::collections::HashMap<String, Vec<&str>> = Default::default();
        let mut leaves0 = 0usize;
        for (name, c) in crits {
            let cfg = SchedCfg {
                n: 8, bootstrap_levels: 2, budget, tau_display: tau, criterion: c,
                hot_rule: HotRule::Quantile(0.5), alpha_hi: 0.2,
                camera: Some(Camera::framing(cx, cy, 0.05, 1024)),
                ..Default::default()
            };
            let (t, _) = scheduler::descend(cx, cy, 0.05, body, &cfg, &now, Precision::F64);
            let mut key: Vec<String> = t
                .leaves()
                .map(|i| {
                    let q = &t.nodes[i];
                    format!("{}:{}:{}", q.level, (q.cx / q.half).round() as i64,
                            (q.cy / q.half).round() as i64)
                })
                .collect();
            key.sort();
            leaves0 = key.len();
            sets.entry(key.join(",")).or_default().push(name);
        }
        let mut groups: Vec<&Vec<&str>> = sets.values().collect();
        groups.sort_by_key(|g| std::cmp::Reverse(g.len()));
        println!("   tau={tau:.0e}  leaves {leaves0:>5}  DISTINCT LEAF SETS: {}", groups.len());
        for g in groups {
            println!("      [{}] {}", g.len(), g.join(" "));
        }
    }

    println!("\nHOW TO READ IT");
    println!("  ONE distinct leaf set at a tau below the bulk = the standing result survives, `far`");
    println!("  is still the control, and the tail does not reach the descent. MORE than one = the");
    println!("  tail is doing work, `far` orders criteria differently now, and every result leaning");
    println!("  on \"thirteen criteria give the identical leaf set\" needs re-reading.");
    println!("  At tau=1e-4 the answer is one set FOR THE WRONG REASON -- 16 leaves, all `keep`,");
    println!("  nothing ever ranked. That row is the degeneracy control, not a result.");
}
