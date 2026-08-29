//! **Can a mechanism be REMOVED?** AZ's `dt = A*B*dtau` is already adaptive; `PerStepInterval`'s
//! cap is a second controller on top, and a per-step limit would be a third. Three step
//! controllers fighting is worse than one that works, and removing a mechanism is worth more than
//! tuning one.
//!
//! # What is asked, and what is deliberately not
//!
//! Only the **cap** is a candidate for removal. `clamp_final_step` is not: it is a *correctness*
//! property, not a step-size one — it lands the final step on the boundary and takes the measured
//! convergence order from 1.06 to 2.08 on the figure-eight. Removing it would reintroduce a
//! first-order error at every boundary, which is a different defect from the one under study.
//!
//! # Already known, and it is why this exists
//!
//! `StepLimit::AbGrowth` is **bitwise inert** under the shipped `PerStepInterval` and bites only
//! under `FixedPerInterval`, because `dtau = eta*dt_left/(A*B)` recomputed per step already holds
//! `dt ~ eta*dt_left` however much `A*B` grows. **`PerStepInterval` IS an `A*B` growth clamp at
//! `C = 1`.** So candidate C was already shipped under another name. The open question is the
//! other direction: with the winning per-step limit on, is the cap still doing anything?
//!
//! Args: `res cap root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::az::{DtauMode, StepLimit};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let res: usize = arg(1, 128);
    let cap: usize = arg(2, 64);
    let root: String = std::env::args().nth(3).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    let (chart, cx, cy, half) = Chart::config_stability();
    let mut cases: Vec<(&str, Slice, EnsembleCfg)> = vec![(
        "config_stability",
        grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart),
        EnsembleCfg {
            refine_flagged: false,
            t_max: 50.0,
            n_sync: (50.0f64 / WINDOW).round() as usize,
            r_coll_frac: 0.005,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU),
            closure_k: 1,
            stop_on_escape: false,
            ..Default::default()
        },
    )];
    cases.push((
        "deep interior",
        grid::region("deep interior", res, res, 0.05).unwrap().with_chart(Chart::BodyPlane),
        EnsembleCfg { refine_flagged: false, ..Default::default() },
    ));

    println!(
        "Two step controllers, or one? `dtau_mode` x `step_limit`, on the flagged-plus-control\n\
         sampled set. **`PerStepInterval` is an `A*B` growth clamp at C = 1** -- so if the\n\
         per-step limit alone reproduces the cap's quality, the cap is redundant and can go.\n"
    );

    for (name, sl, base) in cases {
        println!("\n================ {name} ================");
        let full: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &base))
            .collect();
        let hot: Vec<usize> =
            (0..full.len()).filter(|&i| full[i].error_ratio > base.refine_threshold).collect();
        let cool: Vec<usize> = (0..full.len())
            .filter(|&i| full[i].error_ratio <= 2.0 && full[i].error_ratio.is_finite())
            .collect();
        let take = |v: &[usize], k: usize| -> Vec<usize> {
            if v.len() <= k { v.to_vec() } else { (0..k).map(|j| v[j * v.len() / k]).collect() }
        };
        let idx: Vec<usize> = take(&hot, cap).into_iter().chain(take(&cool, cap)).collect();
        println!("  sampled {} ({} flagged of {} in frame)\n", idx.len(), take(&hot, cap).len(), hot.len());

        println!(
            "  {:>18} {:>14} {:>8} {:>11} {:>11} {:>9} {:>11} {:>10}",
            "dtau_mode", "step_limit", "secs", "err p90", "err p99", "err>10", "steps p50",
            "overshoot"
        );
        for dm in [DtauMode::PerStepInterval, DtauMode::FixedPerInterval] {
            for (lim, f) in [(StepLimit::None, 0.0), (StepLimit::Predictive, 0.02)] {
                let cfg = EnsembleCfg {
                    dtau_mode: dm,
                    step_limit: lim,
                    step_limit_f: f,
                    ..base
                };
                let t = std::time::Instant::now();
                let px: Vec<PixelOut> =
                    idx.par_iter().map(|&k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
                let secs = t.elapsed().as_secs_f64();
                let mut e: Vec<f64> =
                    px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
                let mut st: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
                println!(
                    "  {:>18} {:>14} {secs:>8.1} {:>11.3e} {:>11.3e} {:>9.4} {:>11.3e} {:>10}",
                    format!("{dm:?}"),
                    format!("{lim:?}"),
                    q(&mut e.clone(), 0.9),
                    q(&mut e, 0.99),
                    px.iter().filter(|p| p.error_ratio > 10.0).count() as f64 / px.len() as f64,
                    q(&mut st, 0.5),
                    px.iter().map(|p| p.n_overshoot).sum::<u64>(),
                );
            }
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         The row that decides is **`FixedPerInterval` + `Predictive`** against\n\
         **`PerStepInterval` + `Predictive`**. If they agree, the cap is redundant once the\n\
         per-step limit is in force and one of the three controllers can be deleted. If the\n\
         `Fixed` row is worse, the cap is doing work the limit does not, and both stay --\n\
         which is a real answer and not a failure of the experiment.\n\n\
         `FixedPerInterval` + `None` is the historical behaviour and is here as the fourth cell,\n\
         because a 2x2 read from three cells is an inference."
    );
}
