//! **Does the secant landing help on a real field, as it does on the figure-eight?**
//!
//! `examples/landing_correction.rs` shows `land_iterate` taking the landing residual from `5.4e-6`
//! to `1.8e-15` and RK4's figure-eight closure from `5.43e-5` to `2.77e-8` — a factor of 1960 for
//! 7% more force evaluations — with KDK unchanged as the control. That is one periodic orbit with
//! no close encounters. **Whether it transfers to a chaotic field with collisions in it is a
//! different question**, and it is the one that decides whether porting the correction to AZ and
//! Heggie is worth invalidating the corpus for.
//!
//! # Why this calls the driver directly
//!
//! `land_iterate` is deliberately not reachable from `EnsembleCfg`: it is off in `pixel.rs`
//! because AZ and Heggie have no landing correction at all, and exposing it as a config knob
//! would invite exactly the arm-asymmetry that whole comparison exists to avoid. So this marches
//! decoded copies directly, the way `heggie_machinery.rs` does, and stays out of the ensemble
//! layer entirely.
//!
//! Copies are decoded **once** and every arm marches the same ones, so nothing between the chart
//! and the march can differ between arms.
//!
//! # What would make this a null, stated first
//!
//! Two things, and both are printed rather than assumed:
//!
//!   - **`corr` must be non-zero.** If the landing step is never the binding one — because the
//!     interval always divides evenly, say — there is nothing to correct and every row is the
//!     same run twice.
//!   - **`land resid` must fall.** That is the quantity being fixed. If it does not move, the
//!     correction is not working and any drift difference is something else.
//!
//! And **KDK is carried as the control**, for the reason it is in the figure-eight version: it is
//! second order, an `O(h^2)` landing was never its constraint, and if the correction improves it
//! too then it is not removing an order-two cap.
//!
//! Termination is **off** and `r_coll = 0`: a trajectory stopped by an event is parked at a close
//! approach where the Cartesian energy is a cancellation of enormous terms, and that flag alone
//! produced five wrong conclusions in a row once.
//!
//! Args: `res max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::logh::{integrate_lh, LhOpts, Stepper};
use prin_rs::physics::Ic;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &[f64], p: f64) -> f64 {
    let mut w: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if w.is_empty() {
        return f64::NAN;
    }
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    w[(((w.len() - 1) as f64) * p).round() as usize]
}

struct Case {
    name: &'static str,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    t_max: f64,
}

fn cases() -> Vec<Case> {
    let (chart, cx, cy, half) = Chart::config_stability();
    let mut v = vec![Case { name: "config_stability", chart, cx, cy, half, body: 0, t_max: 50.0 }];
    for name in ["deep interior", "near-field"] {
        if let Some(s) = grid::region(name, 4, 4, 0.05) {
            v.push(Case {
                name: if name == "deep interior" { "deep_interior" } else { "near-field" },
                chart: s.chart, cx: s.cx, cy: s.cy, half: s.half, body: s.body,
                t_max: 13.0,
            });
        }
    }
    v
}

fn main() {
    let res: usize = arg(1, 96);
    let max_steps: usize = arg(2, 400_000);
    let cfg = EnsembleCfg::production();

    println!(
        "{res}^2, termination OFF and r_coll = 0, predictive limit off, eta = {}, copies = {}\n",
        cfg.eta,
        cfg.n_extra + 1
    );
    println!(
        "  **`corr` and `land resid` are the guards.** If `corr` is zero the landing step is never\n  \
         the binding one and every pair of rows is one run twice; if `land resid` does not fall,\n  \
         the correction is not working and a drift difference is something else.\n"
    );
    println!(
        "  {:>18} {:>6} {:>6} {:>11} {:>11} {:>11} {:>11} {:>7} {:>8}",
        "case", "step", "land", "drift p50", "drift p99", "land resid", "evals p50", "corr", "nonfin"
    );

    for c in cases() {
        let sl = grid::Slice::body_plane(res, res, c.cx, c.cy, c.half, c.body).with_chart(c.chart);
        let n_sync = (c.t_max / 0.4).round().max(4.0) as usize;
        // Decoded once; every arm marches the same copies.
        let copies: Vec<Vec<Ic<f64>>> = (0..sl.npix())
            .into_par_iter()
            .map(|k| {
                jitter::copies_with_path::<f64>(
                    &sl, k, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme,
                    cfg.decode_path,
                )
            })
            .collect();

        for stepper in [Stepper::Rk4, Stepper::Gbs, Stepper::Kdk] {
            for land in [false, true] {
                let eta = if stepper == Stepper::Kdk { cfg.eta * 0.25 } else { cfg.eta };
                let out: Vec<(f64, f64, f64, u64, bool)> = copies
                    .par_iter()
                    .map(|cs| {
                        // `energy_drift_max` over the copies, matching `PixelOut`'s own reduction.
                        // A non-finite copy contributes `inf` rather than being dropped: it is a
                        // measurement outcome, not missing data.
                        let (mut d, mut resid, mut ev, mut corr, mut fin) =
                            (0.0f64, 0.0f64, 0u64, 0u64, true);
                        for ic in cs {
                            let o = integrate_lh(
                                ic.s, &ic.m, c.t_max, n_sync, eta, max_steps,
                                &LhOpts {
                                    stepper,
                                    land_iterate: land,
                                    r_coll_frac: 0.0,
                                    stop_on_event: false,
                                    step_limit_f: 0.0,
                                    ..Default::default()
                                },
                            );
                            d = d.max(if o.finite { o.drift } else { f64::INFINITY });
                            resid = resid.max(o.land_residual_max);
                            ev += o.force_evals as u64;
                            corr += o.land_iters;
                            fin &= o.finite;
                        }
                        (d, resid, ev as f64, corr, fin)
                    })
                    .collect();

                let dr: Vec<f64> = out.iter().map(|x| x.0).collect();
                let rs: Vec<f64> = out.iter().map(|x| x.1).collect();
                let ev: Vec<f64> = out.iter().map(|x| x.2).collect();
                let corr: u64 = out.iter().map(|x| x.3).sum();
                let nonfin = out.iter().filter(|x| !x.4).count();
                println!(
                    "  {:>18} {:>6} {:>6} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {corr:>7} {nonfin:>8}",
                    c.name, format!("{stepper:?}"), land,
                    q(&dr, 0.50), q(&dr, 0.99), q(&rs, 0.50), q(&ev, 0.50)
                );
            }
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         **`land resid` is the quantity and `drift` is the consequence.** The figure-eight says\n\
         the residual falls nine orders and a better-than-second-order stepper gains ~1960x. This\n\
         asks whether that survives a chaotic field with close approaches in it.\n\n\
         **KDK is the control.** Second order, so an O(h^2) landing was never its binding\n\
         constraint. If it improves here too, the correction is doing something other than\n\
         removing an order-two cap and the figure-eight account is wrong.\n\n\
         **If the RK4 and GBS rows improve materially, that is the case for porting the secant\n\
         landing to AZ and Heggie** -- whose measured orders of 2.08 and 2.40 are exactly what an\n\
         O(h^2) landing allows. That port would invalidate every committed number in `results/`,\n\
         so it is a decision and not a consequence."
    );
}
