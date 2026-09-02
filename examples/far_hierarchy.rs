//! **Why does `far` beat Heggie on every one of 65536 pixels?**
//!
//! It is the only AZ win in 32 cases, and it is total: `frac better 0.000`, a flat 0.7-0.9
//! decades. A result that only ever wins is weaker evidence than one that loses where a mechanism
//! says it should, so this is worth more than its size.
//!
//! # The hypothesis that is already dead
//!
//! **Not scale.** `far` spans body positions to 13 units where the latent charts sit at `R = 1`,
//! and `Gamma*` is degree six in the coordinates where AZ's `Gamma` is linear in `A` and `B`, so
//! conditioning-with-scale was the obvious guess. It is refuted by this project's own test:
//! `the_march_respects_the_scale_gauge` reads **exactly `0.000e0`** at `alpha = 0.25` and `4`.
//! Heggie's error is invariant under a global rescaling, bitwise. "far is wide" cannot be it.
//!
//! # The hypothesis under test: HIERARCHY, which is scale-invariant
//!
//! In `far` one body sits distant and **stays** distant. The configuration is hierarchical
//! throughout, so AZ's reference-body choice is ideal and never has to change: the two
//! regularised pairs are the two short sides and the unregularised side is genuinely long for the
//! whole run. **AZ is at its best exactly when it never re-registers.** Heggie treats all three
//! symmetrically and carries `R_3 ~ 13` through every term of `Gamma*` for nothing.
//!
//! That closes the loop with the finding this port was built on -- re-registration costs 0.444
//! decades at fixed step size -- because the same mechanism would explain the win *and* the loss.
//!
//! # Four measurements, and the third is the one with teeth
//!
//! 1. **Per-pixel gain against AZ's `switches`, BLOCKED WITHIN EACH CASE.** Pooled across cases it
//!    would be confounded by case identity: `switches` and difficulty are both properties of a
//!    region, so a pooled correlation would report that hard regions are hard.
//! 2. **The `switches` distribution per case.** If `far` is not near zero, or `config_stability`
//!    is not high, the hypothesis dies here for the price of one render.
//! 3. **The control: reduce AZ's switching on a case where it switches a lot, and see whether AZ
//!    improves.** `ref_hysteresis` is the graded knob and `forced_refs` the extreme one. A
//!    correlation cannot separate "switching hurts AZ" from "AZ switches where it is struggling
//!    anyway"; **making it switch less can.**
//! 4. **A hierarchy measure that is not an AZ internal** -- `d[2]/d[0]` at `t = 0`, pure geometry
//!    -- so the answer does not rest on one quantity that might be measuring something else.
//!
//! # The way this fails to fire
//!
//! `far` may have **zero variance** in `switches` -- every pixel at 0 -- which is exactly the case
//! where a within-case correlation cannot be computed. That is not a null result, it is an
//! undefined one, and the distinct-value count is printed so it cannot be read as a null. If it
//! happens, the correlation is measured on the cases that do vary and `far` is only the endpoint.
//!
//! Args: `n root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::physics::{newton, Ic};

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// Spearman rank correlation. `NaN` when either side has fewer than two distinct values —
/// **undefined, not zero**: a constant `x` has no correlation and saying so is the point.
fn spearman(x: &[f64], y: &[f64]) -> (f64, usize, usize) {
    let n = x.len();
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
        let mut r = vec![0.0; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0;
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    };
    let dx = {
        let mut u: Vec<f64> = x.to_vec();
        u.sort_by(|a, b| a.partial_cmp(b).unwrap());
        u.dedup();
        u.len()
    };
    let dy = {
        let mut u: Vec<f64> = y.to_vec();
        u.sort_by(|a, b| a.partial_cmp(b).unwrap());
        u.dedup();
        u.len()
    };
    if dx < 2 || dy < 2 || n < 3 {
        return (f64::NAN, dx, dy);
    }
    let (rx, ry) = (rank(x), rank(y));
    let (mx, my) = (rx.iter().sum::<f64>() / n as f64, ry.iter().sum::<f64>() / n as f64);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (a, b) = (rx[i] - mx, ry[i] - my);
        sxy += a * b;
        sxx += a * a;
        syy += b * b;
    }
    (sxy / (sxx * syy).sqrt(), dx, dy)
}

fn az_opts<'a>(hyst: f64, forced: Option<&'a [u8]>) -> AzOpts<'a, f64> {
    AzOpts {
        stop_on_event: false,
        r_coll_frac: 0.0,
        step_limit: az::StepLimit::Predictive,
        step_limit_f: 0.02,
        ref_hysteresis: hyst,
        forced_refs: forced,
        ..Default::default()
    }
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
    let mut v = Vec::new();
    let (chart, cx, cy, half) = Chart::config_stability();
    v.push(Case { name: "config_stability", chart, cx, cy, half, body: 0, t_max: 50.0 });
    for name in ["far", "mid-field", "near-field", "deep interior", "body2 core"] {
        if let Some(sl) = grid::region(name, 4, 4, 0.05) {
            v.push(Case {
                name: Box::leak(name.to_string().into_boxed_str()),
                chart: sl.chart, cx: sl.cx, cy: sl.cy, half: sl.half, body: sl.body,
                t_max: 13.0,
            });
        }
    }
    v
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(128);
    let cfg = EnsembleCfg::production();

    println!("{n}^2 per case, nominal copy only, no termination, predictive limit on both.\n");
    println!("== §2  SWITCHES, and the hierarchy at t = 0 ==");
    println!(
        "  {:>18} {:>9} {:>9} {:>9} {:>10} {:>12} {:>12}",
        "case", "sw p10", "sw p50", "sw p90", "sw distinct", "hier p50", "gain p50"
    );

    let mut keep: Vec<(&'static str, Vec<f64>, Vec<f64>, Vec<f64>)> = Vec::new();

    for c in cases() {
        let n_sync = (c.t_max / 0.4).round().max(4.0) as usize;
        let sl = grid::Slice::body_plane(n, n, c.cx, c.cy, c.half, c.body).with_chart(c.chart);
        let ics: Vec<Ic<f64>> = (0..sl.npix()).map(|k| sl.nominal_ic::<f64>(k)).collect();

        let out: Vec<(f64, f64, f64, f64)> = ics
            .par_iter()
            .map(|ic| {
                let a = az::integrate_az_opts(
                    ic.s, &ic.m, c.t_max, n_sync, cfg.eta, 400_000, &az_opts(0.0, None),
                );
                let h = integrate_hg(
                    ic.s, &ic.m, c.t_max, n_sync, cfg.eta, 400_000,
                    &HgOpts { r_coll_frac: 0.0, stop_on_event: false, ..Default::default() },
                );
                // Hierarchy at t = 0: longest side over shortest. Pure geometry, scale-invariant,
                // and not an AZ internal — so §1's answer does not rest on one quantity.
                let mut d = newton::pair_dists(&ic.s.r);
                d.sort_by(|x, y| x.partial_cmp(y).unwrap());
                let hier = d[2] / d[0].max(f64::MIN_POSITIVE);
                let gain = if a.drift > 0.0 && h.drift > 0.0 && a.finite && h.finite {
                    (a.drift / h.drift).log10()
                } else {
                    f64::NAN
                };
                (a.switches as f64, hier, gain, a.drift)
            })
            .collect();

        let sw: Vec<f64> = out.iter().map(|x| x.0).collect();
        let hier: Vec<f64> = out.iter().map(|x| x.1).collect();
        let gain: Vec<f64> = out.iter().map(|x| x.2).collect();
        let mut swd = sw.clone();
        swd.sort_by(|a, b| a.partial_cmp(b).unwrap());
        swd.dedup();
        println!(
            "  {:>18} {:>9.0} {:>9.0} {:>9.0} {:>10} {:>12.3e} {:>12.2}",
            c.name,
            q(&mut sw.clone(), 0.10),
            q(&mut sw.clone(), 0.50),
            q(&mut sw.clone(), 0.90),
            swd.len(),
            q(&mut hier.clone(), 0.5),
            q(&mut gain.clone().into_iter().filter(|x| x.is_finite()).collect(), 0.5),
        );
        keep.push((c.name, sw, hier, gain));
    }

    println!("\n== §1  GAIN AGAINST SWITCHES, BLOCKED WITHIN EACH CASE ==");
    println!("  A pooled correlation across cases would report that hard regions are hard:");
    println!("  `switches` and difficulty are both properties of a region.\n");
    println!(
        "  {:>18} {:>10} {:>10} {:>12} {:>10} {:>10}",
        "case", "rho(sw)", "rho(hier)", "n", "sw distinct", "note"
    );
    for (name, sw, hier, gain) in &keep {
        let idx: Vec<usize> = (0..gain.len()).filter(|&i| gain[i].is_finite()).collect();
        let g: Vec<f64> = idx.iter().map(|&i| gain[i]).collect();
        let s: Vec<f64> = idx.iter().map(|&i| sw[i]).collect();
        let h: Vec<f64> = idx.iter().map(|&i| hier[i]).collect();
        let (rs, ds, _) = spearman(&s, &g);
        let (rh, _, _) = spearman(&h, &g);
        let note = if ds < 2 { "UNDEFINED: switches constant" } else { "" };
        println!("  {name:>18} {rs:>10.3} {rh:>10.3} {:>12} {ds:>10} {note:>10}", g.len());
    }
    println!(
        "\n  **A constant `switches` gives UNDEFINED, not zero.** `far` is expected to sit at one\n\
         value, which is exactly the case a correlation cannot speak about — and printing the\n\
         distinct count is what stops that being read as a null."
    );

    println!("\n== §3  THE CONTROL: make AZ switch LESS, and see whether AZ improves ==");
    println!("  A correlation cannot separate `switching hurts AZ` from `AZ switches where it is");
    println!("  struggling anyway`. Reducing the switching can. `ref_hysteresis` is the graded");
    println!("  knob; a frozen reference is the extreme, and it is expected to be WORSE — a bad");
    println!("  fixed choice costs more than re-registration, which is why AZ re-chooses at all.\n");
    println!(
        "  {:>18} {:>10} {:>11} {:>11} {:>11}",
        "case", "hysteresis", "sw p50", "AZ drift p50", "vs baseline"
    );
    for c in cases() {
        if c.name != "config_stability" && c.name != "deep interior" {
            continue;
        }
        let n_sync = (c.t_max / 0.4).round().max(4.0) as usize;
        let sl = grid::Slice::body_plane(n / 2, n / 2, c.cx, c.cy, c.half, c.body)
            .with_chart(c.chart);
        let ics: Vec<Ic<f64>> = (0..sl.npix()).map(|k| sl.nominal_ic::<f64>(k)).collect();
        let mut base = f64::NAN;
        for hyst in [0.0f64, 0.02, 0.10, 0.50] {
            let out: Vec<(f64, f64)> = ics
                .par_iter()
                .map(|ic| {
                    let a = az::integrate_az_opts(
                        ic.s, &ic.m, c.t_max, n_sync, cfg.eta, 400_000, &az_opts(hyst, None),
                    );
                    (a.switches as f64, if a.finite { a.drift } else { f64::NAN })
                })
                .collect();
            let mut sw: Vec<f64> = out.iter().map(|x| x.0).collect();
            let mut dr: Vec<f64> =
                out.iter().map(|x| x.1).filter(|x| x.is_finite()).collect();
            let d50 = q(&mut dr, 0.5);
            if hyst == 0.0 {
                base = d50;
            }
            println!(
                "  {:>18} {hyst:>10.2} {:>11.1} {d50:>11.3e} {:>11.2}",
                c.name,
                q(&mut sw, 0.5),
                d50 / base
            );
        }
        println!();
    }

    println!(
        "HOW TO READ IT.\n\n\
         **The hypothesis predicts three things at once, and all three must hold.** `far` sits at\n\
         ~0 switches and high hierarchy; the cases Heggie wins sit at high switches; and within a\n\
         case that varies, gain rises with switches. Any one of them alone is a coincidence\n\
         waiting to be found out.\n\n\
         **§3 is the part that can refute it.** If cutting AZ's switching leaves its drift alone,\n\
         then switching is a symptom of difficulty rather than a cause of error, the `far` story\n\
         is wrong, and the 0.444-decade re-registration result does not extend to explaining WHERE\n\
         each integrator wins -- only that the cadence moves the field."
    );
}
