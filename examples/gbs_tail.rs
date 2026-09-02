//! **What is in GBS's `err>10` tail: the method, or a truncated trajectory?**
//!
//! `logh_arms` reports `config_stability` under `logh_gbs` at `err>10 = 131` with `budget = 35`,
//! and under `logh_gbs_nolim` at `err>10 = 718` with `budget = 0` and `over = 234`. Two different
//! failure shapes, and the obvious question — how much of the tail is budget artefact — cannot be
//! asked of that table.
//!
//! # The column it cannot be asked of, and why
//!
//! `logh_arms` counts `budget` over `px.iter().chain(dpx.iter())` — **both passes**, 2 x `npix`
//! pixel-outs — while `err>10`, `drift`, `nonfin`, `steps`, `evals` and `over` are all over the
//! diagnostic pass alone. So the budget column has twice the denominator of every other column in
//! its own row and mixes in a population the error columns never see. A truncated run in the
//! science pass cannot contribute to `err>10` at all, so "35 of 131" was never a subset claim.
//!
//! **This harness runs one pass and reports every column over it.**
//!
//! # And the counter that would answer it did not exist above the driver
//!
//! `LhOut::gbs_unconverged` — the macro-step took the extrapolated state without meeting
//! tolerance — was computed on every march since GBS landed and read by nothing: it stopped at
//! `LhOut` and never reached `MarchOut` or `PixelOut`, the same way `ab_floored` and `ab_min` did.
//! It is the *advance-anyway* site for this stepper, and unlike `budget_exhausted` it is not
//! terminal.
//!
//! # Base rates are printed above the lifts
//!
//! A predictor that fires everywhere has a lift of exactly 1.000 by arithmetic, and this project
//! has read one of those as a null once already. Every row carries its own base rate, and a
//! saturated or empty predictor is labelled rather than scored.
//!
//! Args: `res case max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::Integrator;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::ensemble::provenance::Override;
use prin_rs::integrate::az::StepLimit;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// `P(hot | pred) / P(hot)`, with the guards that keep it from being read as a finding.
fn lift(px: &[PixelOut], pred: &dyn Fn(&PixelOut) -> bool) -> (usize, f64, String) {
    let n = px.len() as f64;
    let hot = px.iter().filter(|p| p.error_ratio > 10.0).count() as f64;
    let k = px.iter().filter(|p| pred(p)).count();
    let both = px.iter().filter(|p| pred(p) && p.error_ratio > 10.0).count() as f64;
    if k == 0 {
        return (0, f64::NAN, "EMPTY -- fires nowhere, cannot score".into());
    }
    if k as f64 == n {
        return (k, 1.0, "SATURATED -- lift is 1.000 by arithmetic".into());
    }
    if hot == 0.0 {
        return (k, f64::NAN, "no hot pixels -- nothing to predict".into());
    }
    let l = (both / k as f64) / (hot / n);
    (k, l, format!("covers {:.4} of the hot set", both / hot))
}

fn main() {
    let res: usize = arg(1, 256);
    let case: String = arg(2, "config_stability".to_string());
    let max_steps: usize = arg(3, 400_000);

    let (chart, cx, cy, half, body, t_max) = if case == "config_stability" {
        let (c, x, y, h) = Chart::config_stability();
        (c, x, y, h, 0usize, 50.0f64)
    } else {
        let s = grid::region(&case.replace('_', " "), 4, 4, 0.05).expect("unknown case");
        (s.chart, s.cx, s.cy, s.half, s.body, 13.0f64)
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
    let n_sync = (t_max / 0.4).round().max(4.0) as usize;

    println!(
        "{case} at {res}^2, DIAGNOSTIC PASS ONLY (termination off, r_coll = 0), \
         refine_flagged = false, max_steps = {max_steps}.\n"
    );
    println!(
        "  Every column below is over the SAME {} pixels. `logh_arms`' `budget` column is not:\n  \
         it is counted over both passes and therefore over twice this denominator.\n",
        sl.npix()
    );

    for (label, limit) in [("logh_gbs", true), ("logh_gbs_nolim", false)] {
        let ens = EnsembleCfg::production().with_overrides(&[
            Override::TMax(t_max),
            Override::NSync(n_sync),
            Override::RCollFrac(0.0),
            Override::EscapeRule(EscapeRule::Closure(CLOSURE_TAU)),
            Override::ClosureK(1),
            Override::Integrator(Integrator::LogHGbs),
            Override::Eta(EnsembleCfg::production().eta * 0.25),
            Override::MaxSteps(max_steps),
            Override::RefineFlagged(false),
            Override::StopOnEvent(false),
            Override::StepLimit(if limit { StepLimit::Predictive } else { StepLimit::None }),
        ]);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();

        let n = px.len();
        let hot = px.iter().filter(|p| p.error_ratio > 10.0).count();
        let bud = px.iter().filter(|p| p.budget_exhausted).count();
        let ovr = px.iter().filter(|p| p.n_overshoot > 0).count();
        let unc = px.iter().filter(|p| p.gbs_unconverged > 0).count();
        let unc_tot: u64 = px.iter().map(|p| p.gbs_unconverged).sum();
        let steps_tot: u64 = px.iter().map(|p| p.total_substeps as u64).sum();

        println!("  == {label} ==");
        println!(
            "    n {n}   err>10 {hot} ({:.6})   budget {bud}   overshoot>0 {ovr}   unconverged>0 {unc}",
            hot as f64 / n as f64
        );
        println!(
            "    unconverged macro-steps {unc_tot} of {steps_tot} total steps ({:.3e} of steps)",
            unc_tot as f64 / steps_tot.max(1) as f64
        );
        // **THE RATE, NOT THE BOOLEAN.** `gbs_unconverged > 0` fires on 83% of pixels and
        // therefore cannot discriminate: it scores a lift of 1.203 while covering 1.0000 of the
        // hot set, which is close to what chance alone gives. The reason is arithmetic and is the
        // standing lesson about a count capped below its own decision threshold -- a trajectory of
        // ~5e5 macro-steps almost always accumulates at least one unconverged step, at a rate of
        // only 6.3e-4 per step.
        //
        // The rate is unbounded and continuous, so it can rank. Reported as the fraction of hot
        // pixels in each rate quintile against the frame base rate: a signal that carries nothing
        // gives a flat column, which is a result rather than an absence.
        let rate: Vec<f64> = px
            .iter()
            .map(|p| {
                if p.total_substeps > 0 {
                    p.gbs_unconverged as f64 / p.total_substeps as f64
                } else {
                    f64::NAN
                }
            })
            .collect();
        let nd = {
            let mut w: Vec<u64> =
                rate.iter().filter(|x| x.is_finite()).map(|x| x.to_bits()).collect();
            w.sort_unstable();
            w.dedup();
            w.len()
        };
        let mut order: Vec<usize> = (0..px.len()).filter(|&i| rate[i].is_finite()).collect();
        order.sort_by(|&a, &b| rate[a].partial_cmp(&rate[b]).unwrap());
        println!("    unconverged RATE: {nd} distinct values over {} pixels", order.len());
        let base = hot as f64 / n as f64;
        for k in 0..5 {
            let lo = order.len() * k / 5;
            let hi = order.len() * (k + 1) / 5;
            let bucket = &order[lo..hi];
            let h = bucket.iter().filter(|&&i| px[i].error_ratio > 10.0).count();
            println!(
                "      quintile {} rate<={:.3e}  hot {h:>5} of {:>5}  ({:.5}, lift {:>6.2})",
                k + 1,
                rate[*bucket.last().unwrap_or(&0)],
                bucket.len(),
                h as f64 / bucket.len().max(1) as f64,
                (h as f64 / bucket.len().max(1) as f64) / base.max(1e-30)
            );
        }

        for (nm, f) in [
            ("budget_exhausted", &(|p: &PixelOut| p.budget_exhausted) as &dyn Fn(&PixelOut) -> bool),
            ("n_overshoot > 0", &|p: &PixelOut| p.n_overshoot > 0),
            ("gbs_unconverged > 0", &|p: &PixelOut| p.gbs_unconverged > 0),
        ] {
            let (k, l, note) = lift(&px, f);
            println!("    {nm:>22}  fires {k:>6}  lift {l:>8.3}  {note}");
        }
        // The residual: hot pixels no advance-anyway counter explains.
        let unexplained = px
            .iter()
            .filter(|p| {
                p.error_ratio > 10.0 && !p.budget_exhausted && p.n_overshoot == 0
                    && p.gbs_unconverged == 0
            })
            .count();
        println!(
            "    hot with NO counter set: {unexplained} of {hot} -- the method's own share\n"
        );
    }

    println!(
        "HOW TO READ THIS\n\n\
         **A lift needs its base rate above it.** A predictor firing on every pixel scores 1.000 \n\
         by arithmetic and one firing on none scores NaN; both are labelled rather than quoted.\n\n\
         **`hot with NO counter set` is the number the comparison turns on.** Budget exhaustion \n\
         truncates a trajectory and an overshoot corrupts one; neither is GBS being inaccurate. \n\
         What is left over is."
    );
}
