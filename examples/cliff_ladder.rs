//! **Is the failure a slope or a floor?** `eta` over decades, on the affected pixels only.
//!
//! `BUG_REPORT.md` §8 records the number this exists to explain: `error_ratio` p99 is **35.6
//! after three passes at `eta/4`** — 64x finer stepping, and the tail does not clear. Either the
//! refinement ladder simply does not go far enough (a **slope**, and the answer is a finer step),
//! or something in the march does not care about the step size at all (a **floor**, and the
//! answer is elsewhere). *A quantity that does not converge under refinement is measuring the
//! sampling rather than the system.*
//!
//! # Read the slope, not the error
//!
//! An error falls for many reasons. Only the slope of `log(error_ratio)` against `log(eta)` says
//! whether the leading term is truncation.
//!
//! # The decoy, and why three columns are printed instead of one ratio
//!
//! `error_ratio = sigma_E(t) / sigma_E(0)`. The denominator is taken at `t = 0`, so it is
//! **`eta`-independent by construction** and every bit of the ratio's `eta`-dependence lives in
//! the numerator — which is what makes it a legitimate convergence probe. But `sigma_E(0)` is
//! proportional to the jitter and so to the cell width, and it is *small*: at fine `eta` the
//! numerator reaches the round-off floor `~eps*|E|` and the ratio then floors at
//! `eps*|E|/sigma_E(0)` for a purely arithmetic reason.
//!
//! So the expected round-off floor is computed per pixel and printed beside the measured value.
//! **A floor above it is a mechanism. A floor at it is arithmetic.** Reading the ratio alone
//! would score the second as the first.
//!
//! **MEASURED, AND THE GUARD WAS AIMED AT THE WRONG SIDE.** The round-off floor comes out at
//! `4.7e-13`, and the ratio bottoms out at **exactly 1.000** — which is not a floor at all but
//! `error_ratio`'s *converged* value, since `sigma_E(t) -> sigma_E(0)` under exact dynamics. The
//! statistic is normalised to 1 by construction, so it can never reach an arithmetic floor
//! beneath it. The decoy is real for an unnormalised quantity and inert for this one; kept, with
//! the reason, because a guard that cannot fire has to be labelled rather than deleted.
//!
//! # This is also the refinement-convergence question, at no extra cost
//!
//! `refine_max_passes = N` **is** this ladder at `eta/4^k` with an early exit. So the run reports
//! at which rung each flagged pixel first clears `refine_threshold` and **what fraction never
//! does** — without a separate experiment, and without `refine_max_passes = 20`, which is
//! `eta/4^20 ~ 1e-12 eta` and absurd as a remedy. If a fraction never clears, the repair pass
//! does not repair either; it reduces the count below visual threshold, which is a different
//! claim.
//!
//! # Costs and caps, stated
//!
//! Each rung is ~4x the steps of the last, so the ladder is ~85x one pass per pixel at four
//! rungs. The pixel set is therefore **capped and the cap is printed** — no silent truncation.
//! Selection runs at 256^2: `error_ratio`, `drift` and `steps` are per-trajectory statistics and
//! are not biased by a coarse grid the way a chord ratio is. **No chord is quoted here.**
//!
//! Args: `res rungs cap root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::energy;

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
    let res: usize = arg(1, 256);
    let rungs: usize = arg(2, 4);
    let cap: usize = arg(3, 64);
    let root: String = std::env::args().nth(4).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    let (chart, cx, cy, half) = Chart::config_stability();
    let (t_max, r_coll) = (50.0, 0.005);
    let n_sync = (t_max / WINDOW).round().max(4.0) as usize;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        r_coll_frac: r_coll,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    println!("config_stability {res}^2 selection, then an eta ladder on the selected pixels.");
    println!("config: {}\n", ens.provenance());

    // --- selection ---------------------------------------------------------------------
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let t0 = std::time::Instant::now();
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();
    println!("selection pass {:.1}s", t0.elapsed().as_secs_f64());

    let flagged_all: Vec<usize> =
        (0..px.len()).filter(|&i| px[i].error_ratio > ens.refine_threshold).collect();
    let control_all: Vec<usize> =
        (0..px.len()).filter(|&i| px[i].error_ratio <= 2.0 && px[i].error_ratio.is_finite()).collect();

    // Evenly spaced, not random: the same pixels every run, and no seed to report.
    let take = |v: &[usize], n: usize| -> Vec<usize> {
        if v.len() <= n {
            return v.to_vec();
        }
        (0..n).map(|k| v[k * v.len() / n]).collect()
    };
    let flagged = take(&flagged_all, cap);
    let control = take(&control_all, cap);
    println!(
        "flagged (error_ratio > {}): {} of {} pixels; laddering {} (cap {cap})\n\
         control (error_ratio <= 2): {} of {}; laddering {}\n\
         **The cap is printed because a silently truncated set reads as full coverage.**\n",
        ens.refine_threshold, flagged_all.len(), px.len(), flagged.len(),
        control_all.len(), px.len(), control.len()
    );
    if flagged.is_empty() {
        println!("NOTHING FLAGGED at this resolution. The ladder has no subject; stopping.");
        return;
    }

    // The round-off floor, per pixel: eps*|E| / sigma_E(0). Needs no integration.
    let e_abs = |i: usize| {
        let (x, y) = sl.decode_pos(i);
        let s = grid::decode_state(&chart, 0, x, y);
        energy::energy(&s.s.r, &s.s.v, &s.m, 0.0).abs()
    };
    let floor_of = |i: usize| f64::EPSILON * e_abs(i) / px[i].sigma_e_0;

    let mut fl: Vec<f64> = flagged.iter().map(|&i| floor_of(i)).collect();
    println!(
        "expected ROUND-OFF floor on the flagged set, eps*|E|/sigma_E(0): \
         p50 {:.3e}  p10 {:.3e}  p90 {:.3e}",
        q(&mut fl.clone(), 0.5), q(&mut fl.clone(), 0.1), q(&mut fl, 0.9)
    );
    println!("**A measured floor ABOVE that is a mechanism; a floor AT it is arithmetic.**\n");

    // --- the ladder --------------------------------------------------------------------
    // (eta, per-pixel error_ratio, sigma_E(t), drift, steps, dt_max, cap hits, floored)
    let mut table: Vec<(f64, Vec<[f64; 6]>, Vec<[f64; 6]>)> = Vec::new();
    for r in 0..rungs {
        let eta = ens.eta * 0.25f64.powi(r as i32);
        let t1 = std::time::Instant::now();
        let run = |set: &[usize]| -> Vec<[f64; 6]> {
            set.par_iter()
                .map(|&i| {
                    let p = pixel::evaluate_at::<f64>(&sl, i, &ens, eta);
                    [
                        p.error_ratio, p.sigma_e_t, p.energy_drift_max,
                        p.total_substeps as f64, p.dt_max, p.ab_floored as u8 as f64,
                    ]
                })
                .collect()
        };
        let f = run(&flagged);
        let c = run(&control);
        println!("  rung eta={eta:.3e} in {:.1}s", t1.elapsed().as_secs_f64());
        table.push((eta, f, c));
    }

    println!("\n== THE LADDER ==");
    println!(
        "  `slope` is d log10(err p50) / d log10(eta) against the previous rung. **Negative and\n\
         near-constant is a slope; ~0 is a floor.** RK4 truncation in the energy would give\n\
         about +4 in `drift` as eta falls, i.e. a slope of +4 here by sign convention.\n"
    );
    for (tag, sel) in [("FLAGGED", 0usize), ("CONTROL", 1usize)] {
        println!("  -- {tag} --");
        println!(
            "  {:>10} {:>11} {:>11} {:>11} {:>11} {:>8} {:>11} {:>10} {:>7}",
            "eta", "err p50", "err p90", "sigE(t) p50", "drift p50", "slope", "steps p50",
            "dt_max p50", "floored"
        );
        let mut prev: Option<(f64, f64)> = None;
        for (eta, f, c) in &table {
            let d = if sel == 0 { f } else { c };
            let col = |k: usize| -> Vec<f64> {
                d.iter().map(|r| r[k]).filter(|x| x.is_finite()).collect()
            };
            let e50 = q(&mut col(0), 0.5);
            let slope = match prev {
                Some((pe, pv)) if pv > 0.0 && e50 > 0.0 => {
                    format!("{:+.2}", (e50.log10() - pv.log10()) / (eta.log10() - pe.log10()))
                }
                _ => "--".into(),
            };
            println!(
                "  {eta:>10.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {slope:>8} {:>11.3e} \
                 {:>10.3e} {:>7.3}",
                e50, q(&mut col(0), 0.9), q(&mut col(1), 0.5), q(&mut col(2), 0.5),
                q(&mut col(3), 0.5), q(&mut col(4), 0.5),
                d.iter().map(|r| r[5]).sum::<f64>() / d.len() as f64
            );
            prev = Some((*eta, e50));
        }
        println!();
    }

    // --- does the refinement pass ever converge? ---------------------------------------
    println!("== DOES THE REPAIR PASS CONVERGE? ==");
    println!(
        "  `refine_max_passes = N` IS the ladder above at eta/4^k with an early exit, so this\n\
         needs no separate run. `cleared at rung k` counts flagged pixels whose error_ratio first\n\
         falls to or below {} at rung k. Rung 3 is the shipped `refine_max_passes = 3`.\n",
        ens.refine_threshold
    );
    let mut first: Vec<Option<usize>> = vec![None; flagged.len()];
    for (r, (_, f, _)) in table.iter().enumerate() {
        for (j, row) in f.iter().enumerate() {
            if first[j].is_none() && row[0] <= ens.refine_threshold {
                first[j] = Some(r);
            }
        }
    }
    for r in 0..rungs {
        let n = first.iter().filter(|x| **x == Some(r)).count();
        println!(
            "  cleared at rung {r} (eta = {:.3e}): {n:>4}  ({:.4} cumulative)",
            ens.eta * 0.25f64.powi(r as i32),
            first.iter().filter(|x| x.map(|k| k <= r).unwrap_or(false)).count() as f64
                / flagged.len() as f64
        );
    }
    let never = first.iter().filter(|x| x.is_none()).count();
    println!(
        "  NEVER cleared: {never} of {} ({:.4})",
        flagged.len(), never as f64 / flagged.len() as f64
    );
    println!(
        "\n  If that last figure is materially above zero the repair pass does not repair --\n\
         it reduces the count below visual threshold, which is a different claim, and the\n\
         cap is not on `eta`."
    );
}
