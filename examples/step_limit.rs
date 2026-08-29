//! **Four step-control candidates, measured against each other. The numbers decide.**
//!
//! The defect: `dt = A*B*dtau` is emergent — `dt/dtau = A*B` is integrated *by* the stepper — so
//! the step taken is not the step predicted. One step advanced the physical clock by `2.209e128`
//! against a sync interval of `0.4` and the march recorded a clean landing. **Nothing asked
//! whether the step it just took was one it could afford.**
//!
//! The two batch remedies are **characterisation, not fixes**: `refine_flagged` re-integrates
//! from `t = 0`, which a live playhead cannot do, and a global `eta/256` pays 256x everywhere for
//! a local failure. That `eta/256` brings every flagged pixel to `error_ratio` 1.000 is what
//! proves this is ordinary under-resolution — and it is why `eta/256` is used here as the
//! **ground truth**, not as a candidate.
//!
//! # The comparison that decides
//!
//! **`error_ratio` p99 against wall clock.** The winner reaches `error_ratio ~ 1` for the least
//! compute, *not* the lowest error. B and C should be near-free, A pays only where it fires, D
//! pays everywhere. If a mode plateaus above 1.0 it does not address the defect and the write-up
//! says so.
//!
//! # Controls, because a mode tuned to one slice is not a fix
//!
//! `preset_shape` (clean — all four must leave it alone) and `deep interior` (known-hard)
//! alongside the affected `config_stability`.
//!
//! # `step_limit_f` means four different things
//!
//! A fraction of `d_min` (A), a crossing-time fraction (B), an `A*B` growth factor (C), an `eta`
//! multiplier (D). One number with four meanings cannot be swept across modes, so each mode
//! carries its own ladder and the table prints the meaning in force.
//!
//! Args: `stage res cap root`, `stage` in `sample | frame`. Both write to `<root>/output/`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::az::StepLimit;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
/// The ground truth's `eta` divisor. From the ladder: at `eta/256` every flagged pixel reads
/// `error_ratio` 1.000, so this is the reference answer rather than merely a finer one.
const TRUTH_DIV: f64 = 256.0;

/// `(mode, meaning of f, ladder)`. Each mode's ladder is its own; see the module docs.
const LADDERS: [(StepLimit, &str, [f64; 3]); 4] = [
    (StepLimit::Reject, "fraction of d_min a step may move", [0.5, 0.1, 0.02]),
    (StepLimit::Predictive, "fraction of a crossing time", [0.5, 0.1, 0.02]),
    (StepLimit::AbGrowth, "A*B growth factor", [4.0, 2.0, 1.1]),
    (StepLimit::Global, "eta multiplier", [0.5, 0.25, 0.125]),
];

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

/// `(name, slice-builder, base config)`.
fn targets(n: usize) -> Vec<(&'static str, Slice, EnsembleCfg)> {
    let mut out: Vec<(&'static str, Slice, EnsembleCfg)> = Vec::new();

    let (chart, cx, cy, half) = Chart::config_stability();
    let n_sync = (50.0f64 / WINDOW).round() as usize;
    out.push((
        "config_stability",
        grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart),
        EnsembleCfg {
            refine_flagged: false,
            t_max: 50.0,
            n_sync,
            r_coll_frac: 0.005,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU),
            closure_k: 1,
            stop_on_escape: false,
            ..Default::default()
        },
    ));

    if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == "preset_shape") {
        out.push((
            "preset_shape",
            grid::Slice::body_plane(n, n, c.2, c.3, c.4, 0).with_chart(c.1),
            EnsembleCfg { refine_flagged: false, ..Default::default() },
        ));
    }

    out.push((
        "deep interior",
        grid::region("deep interior", n, n, 0.05).unwrap().with_chart(Chart::BodyPlane),
        EnsembleCfg { refine_flagged: false, ..Default::default() },
    ));
    out
}

struct Row {
    label: String,
    secs: f64,
    err: [f64; 4],
    err_hot: f64,
    drift: [f64; 3],
    steps: [f64; 3],
    esc: f64,
    col: f64,
    overshoot: u64,
    retry: u64,
    retry_exh: usize,
    chord: [f64; 2],
    /// Chord against ground truth over the **control** half only. The pooled figure is saturated
    /// by chaotic divergence -- over `t = 50` any change of step size gives a different
    /// trajectory, so `flips` reads 1.0 for a correct mode and a broken one alike. The tame half
    /// is the arm that can still discriminate.
    chord_tame: f64,
    flips: f64,
}

fn summarise(
    label: String,
    secs: f64,
    px: &[PixelOut],
    truth: Option<&[PixelOut]>,
    n_flagged: usize,
) -> Row {
    let mut e: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
    let mut d: Vec<f64> =
        px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
    let mut st: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
    let n = px.len() as f64;
    let (mut chord, mut flips, mut chord_tame) = ([f64::NAN; 2], f64::NAN, f64::NAN);
    if let Some(t) = truth {
        let mut c: Vec<f64> = px
            .iter()
            .zip(t)
            .map(|(a, b)| {
                (0..3).map(|k| (a.shape_vec[k] - b.shape_vec[k]).powi(2)).sum::<f64>().sqrt()
            })
            .filter(|x| x.is_finite())
            .collect();
        let mut tame: Vec<f64> = c.iter().skip(n_flagged).cloned().collect();
        chord_tame = q(&mut tame, 0.5);
        chord = [q(&mut c.clone(), 0.5), q(&mut c, 1.0)];
        flips = px.iter().zip(t).filter(|(a, b)| a.outcome != b.outcome).count() as f64 / n;
    }
    Row {
        label,
        secs,
        err: [
            q(&mut e.clone(), 0.5),
            q(&mut e.clone(), 0.9),
            q(&mut e.clone(), 0.99),
            q(&mut e, 1.0),
        ],
        err_hot: px.iter().filter(|p| p.error_ratio > 10.0).count() as f64 / n,
        drift: [q(&mut d.clone(), 0.5), q(&mut d.clone(), 0.99), q(&mut d, 1.0)],
        steps: [q(&mut st.clone(), 0.5), q(&mut st.clone(), 0.99), q(&mut st, 1.0)],
        esc: px.iter().filter(|p| p.state == 0).count() as f64 / n,
        col: px.iter().filter(|p| p.state == 2).count() as f64 / n,
        overshoot: px.iter().map(|p| p.n_overshoot).sum(),
        retry: px.iter().map(|p| p.n_retry).sum(),
        retry_exh: px.iter().filter(|p| p.retry_exhausted).count(),
        chord,
        chord_tame,
        flips,
    }
}

fn header() {
    println!(
        "  {:>26} {:>8} {:>10} {:>10} {:>10} {:>9} {:>10} {:>10} {:>11} {:>8} {:>8} {:>9} {:>10} {:>8}",
        "mode / f", "secs", "err p50", "err p90", "err p99", "err>10", "drift p50", "drift p99",
        "steps p50", "escape", "chord tame", "flips", "overshoot", "retries"
    );
}

fn print_row(r: &Row) {
    println!(
        "  {:>26} {:>8.1} {:>10.3e} {:>10.3e} {:>10.3e} {:>9.4} {:>10.3e} {:>10.3e} {:>11.3e} \
         {:>8.4} {:>8.3e} {:>9.4} {:>10} {:>8}",
        r.label, r.secs, r.err[0], r.err[1], r.err[2], r.err_hot, r.drift[0], r.drift[1],
        r.steps[0], r.esc, r.chord_tame, r.flips, r.overshoot, r.retry
    );
}

fn main() {
    let stage: String = std::env::args().nth(1).unwrap_or_else(|| "sample".into());
    let res: usize = arg(2, 192);
    let cap: usize = arg(3, 96);
    let root: String = std::env::args().nth(4).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    println!(
        "STEP-CONTROL CANDIDATES -- stage `{stage}`, {res}^2.\n\n\
         `err` is `error_ratio`; its healthy value is exactly 1.0 and `err>10` is the project's\n\
         own flag for *this pixel is not data*. `chord50` and `flips` are against the eta/{:.0}\n\
         GROUND TRUTH, which is the reference answer because that run brings every flagged pixel\n\
         to 1.000. **The winner reaches err ~ 1 for the least wall clock, not the lowest error.**\n",
        TRUTH_DIV
    );

    for (name, sl, base) in targets(res) {
        println!("\n================ {name} ================");
        println!("base config: {}\n", base.provenance());

        // --- the pixel set --------------------------------------------------------------
        let t0 = std::time::Instant::now();
        let full: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &base))
            .collect();
        let base_secs = t0.elapsed().as_secs_f64();

        let hot: Vec<usize> =
            (0..full.len()).filter(|&i| full[i].error_ratio > base.refine_threshold).collect();
        let cool: Vec<usize> = (0..full.len())
            .filter(|&i| full[i].error_ratio <= 2.0 && full[i].error_ratio.is_finite())
            .collect();
        // Evenly spaced, never random: the same pixels every run and no seed to report.
        let take = |v: &[usize], k: usize| -> Vec<usize> {
            if v.len() <= k { v.to_vec() } else { (0..k).map(|j| v[j * v.len() / k]).collect() }
        };
        let n_flagged = take(&hot, cap).len();
        let idx: Vec<usize> =
            take(&hot, cap).into_iter().chain(take(&cool, cap)).collect();
        println!(
            "  flagged {} of {} pixels; sampled set is {} ({} flagged + {} control). \
             **The cap is printed because a silently truncated set reads as full coverage.**\n\
             FULL-FRAME err>10 is {:.4} in {base_secs:.1}s. **The sampled `err>10` starts at\n\
             {:.4} BY CONSTRUCTION** -- the set is half flagged on purpose -- so read the sampled\n\
             column as a fall from that, never as a frame statistic.",
            hot.len(), full.len(), idx.len(), n_flagged, take(&cool, cap).len(),
            hot.len() as f64 / full.len() as f64,
            n_flagged as f64 / idx.len() as f64
        );

        let run_on = |cfg: &EnsembleCfg, set: &[usize]| -> (Vec<PixelOut>, f64) {
            let t = std::time::Instant::now();
            let v: Vec<PixelOut> =
                set.par_iter().map(|&k| pixel::evaluate::<f64>(&sl, k, cfg)).collect();
            (v, t.elapsed().as_secs_f64())
        };

        if stage == "sample" {
            // --- ground truth ------------------------------------------------------------
            let t = std::time::Instant::now();
            let truth: Vec<PixelOut> = idx
                .par_iter()
                .map(|&k| pixel::evaluate_at::<f64>(&sl, k, &base, base.eta / TRUTH_DIV))
                .collect();
            let truth_secs = t.elapsed().as_secs_f64();
            let mut te: Vec<f64> =
                truth.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
            println!(
                "  ground truth at eta/{:.0} in {truth_secs:.1}s -- err p50 {:.4} p99 {:.4}, \
                 which is what makes it the reference.\n",
                TRUTH_DIV,
                q(&mut te.clone(), 0.5),
                q(&mut te, 0.99)
            );

            header();
            let (b, s) = run_on(&base, &idx);
            print_row(&summarise("None (baseline)".into(), s, &b, Some(&truth), n_flagged));
            for (mode, meaning, ladder) in LADDERS {
                println!("  -- {mode:?}: f is a {meaning}");
                for f in ladder {
                    let cfg = EnsembleCfg { step_limit: mode, step_limit_f: f, ..base };
                    let (v, s) = run_on(&cfg, &idx);
                    let mut r = summarise(format!("{mode:?} f={f}"), s, &v, Some(&truth), n_flagged);
                    r.label = format!("{mode:?} f={f}");
                    print_row(&r);
                    if r.retry_exh > 0 {
                        println!(
                            "  {:>26}   undetermined (retry budget exhausted): {} pixels",
                            "", r.retry_exh
                        );
                    }
                }
            }
        } else {
            // --- full frame, for wall clock and distributions ----------------------------
            println!("  full frame, baseline in {base_secs:.1}s\n");
            header();
            print_row(&summarise("None (baseline)".into(), base_secs, &full, None, 0));
            for (mode, _, ladder) in LADDERS {
                // **Which rungs run at frame scale is a cost decision and is stated.** The two
                // live candidates get their middle and finest rungs, because the sample stage
                // showed the middle one understates B. `AbGrowth` gets one -- it is bitwise
                // inert under the shipped `DtauMode` and a second rung would print the baseline
                // again. `Global` is the control and its whole point is the cost curve, so it
                // keeps two; its finest rung is 8x the baseline and is not run at frame scale.
                let rungs: &[f64] = match mode {
                    StepLimit::AbGrowth => &ladder[1..2],
                    StepLimit::Global => &ladder[0..2],
                    _ => &ladder[1..3],
                };
                for &f in rungs {
                let cfg = EnsembleCfg { step_limit: mode, step_limit_f: f, ..base };
                let t = std::time::Instant::now();
                let v: Vec<PixelOut> = (0..sl.npix())
                    .into_par_iter()
                    .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
                    .collect();
                print_row(&summarise(
                    format!("{mode:?} f={f}"),
                    t.elapsed().as_secs_f64(),
                    &v,
                    None,
                    0,
                ));
                }
            }
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **`overshoot` must be zero.** It counts steps that carried the interval clock past twice\n\
         its own interval -- `dt > dt_left` is a bug, not a condition to handle. A nonzero\n\
         baseline column is the defect being measured; a nonzero column under a candidate means\n\
         that candidate does not address it.\n\n\
         **`err>10` against `secs` is the decision.** A mode that reaches the baseline's error at\n\
         the baseline's cost has done nothing; one that reaches ~1.0 for a few percent has won.\n\
         `Global` is the control and is expected to pay in proportion to its `eta` multiplier --\n\
         if it does not beat the others on cost it should not ship, and if it does, it should.\n\n\
         **The controls decide whether a win generalises.** A mode that only helps\n\
         `config_stability` is tuned to it. `preset_shape` must be untouched by all four."
    );
}
