//! **STAGE 1 — reversibility as an instrument, normalised by the amplification.**
//!
//! March to `t_max`, negate the velocities, march the same span, compare COM-centred positions
//! against the initial condition. One scalar per pixel, no reference trajectory, 2x cost.
//!
//! # Why the raw number must not be reported alone
//!
//! On a chaotic field it measures **round-off integrated through the tangent flow** — Lyapunov
//! amplification times machine epsilon. That is what "irreversible to the Planck length" *is*, and
//! quoting it as integrator error would be reading the field's chaos as the integrator's fault.
//! So both forms are reported:
//!
//!   - `rev`        -- the raw residual, normalised by the initial hyperradius (scale-invariant).
//!   - `rev/amp`    -- divided by `exp(lambda * t_max)` with `lambda` the pixel's own FTLE. This
//!                     is round-off **per unit amplification**, and it is the integrator-quality
//!                     arm.
//!
//! **And the normalisation carries its own falsification test.** If `rev/amp` still tracks FTLE as
//! strongly as `rev` does, the division did not remove the amplification and the second column is
//! not integrator quality — it is the first column wearing a different name. That comparison is
//! printed, not assumed.
//!
//! # Three things enter and none can be removed
//!
//! 1. Round-off amplified by the flow (above).
//! 2. **Stepper time-symmetry.** KDK retraces exactly in exact arithmetic; classical RK4 does not.
//!    So reversibility compared **across steppers scores symmetry, not accuracy**. `az`, `heggie`
//!    and `logh_rk4` are all RK4 and are mutually comparable. Read down a column.
//! 3. **Step-control asymmetry.** The reverse leg picks its own step sequence from its own states
//!    and does not mirror the forward one, so even a symmetric stepper does not retrace under
//!    adaptive `eta`. Measured: GBS saturates rather than converging, and *worsens* under
//!    refinement.
//!
//! # The confound in the FTLE correlation, stated rather than fixed
//!
//! `src/physics/ftle.rs` integrates the tangent vector on the **plain leapfrog**. So an arm using
//! a leapfrog shares a stepper with the field it is correlated against and an RK4 arm does not.
//! A **half-frame shifted control** is on every row: it destroys any real spatial relationship
//! while preserving both marginal distributions, so a correlation that survives the shift was
//! never spatial. AZ's published drift-vs-FTLE reads `-0.0820` against a shifted `-0.1022` — that
//! is a **null**, not a negative correlation, and it is why the control is printed beside every
//! figure rather than quoted once.
//!
//! Args: `res case max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts, StepLimit};
use prin_rs::integrate::heggie::{self, HgOpts};
use prin_rs::integrate::logh::{integrate_lh, LhOpts, LhTime, Stepper};
use prin_rs::integrate::Integrator;
use prin_rs::physics::{energy, Cart, Ic};
use prin_rs::vec2::Vec2;

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

/// Spearman rank correlation over the pairs where **both** are finite, with the count returned:
/// a coefficient over an unstated population is not a measurement.
fn spearman(a: &[f64], b: &[f64]) -> (f64, usize) {
    let pairs: Vec<(f64, f64)> =
        a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| (*x, *y)).collect();
    let n = pairs.len();
    if n < 8 {
        return (f64::NAN, n);
    }
    let rank = |get: &dyn Fn(&(f64, f64)) -> f64| {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&i, &j| get(&pairs[i]).partial_cmp(&get(&pairs[j])).unwrap());
        let mut r = vec![0.0f64; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && get(&pairs[idx[j + 1]]) == get(&pairs[idx[i]]) {
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
    let (ra, rb) = (rank(&|p| p.0), rank(&|p| p.1));
    let m = (n as f64 - 1.0) / 2.0;
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (x, y) = (ra[i] - m, rb[i] - m);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    (num / (da.sqrt() * db.sqrt()), n)
}

/// The half-frame shifted control: same marginals, spatial relationship destroyed.
fn shift_half(v: &[f64], res: usize) -> Vec<f64> {
    let mut o = vec![0.0; v.len()];
    for y in 0..res {
        for x in 0..res {
            o[y * res + x] = v[((y + res / 2) % res) * res + x];
        }
    }
    o
}

fn flip(s: Cart<f64>) -> Cart<f64> {
    let mut o = s;
    for i in 0..3 {
        o.v[i] = -o.v[i];
    }
    o
}

/// COM-removed positions. **Not tidiness.** AZ and Heggie reconstruct from relative coordinates
/// and so place the COM at the origin; logH integrates absolute positions and leaves it where the
/// decode put it. On a `body_plane` slice the decode is not COM-centred, so the two families sit a
/// constant translation apart -- measured `(-0.0125, +2.4875)` on `far`, identical for all three
/// bodies. Every other comparison in this project is translation-invariant, so this is the first
/// quantity here that has to say which frame it is in.
fn centred(s: &Cart<f64>, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let mt = m[0] + m[1] + m[2];
    let com = (s.r[0] * m[0] + s.r[1] * m[1] + s.r[2] * m[2]) / mt;
    [s.r[0] - com, s.r[1] - com, s.r[2] - com]
}

fn sep(a: &Cart<f64>, b: &Cart<f64>, m: &[f64; 3]) -> f64 {
    let (x, y) = (centred(a, m), centred(b, m));
    (0..3).map(|i| (x[i] - y[i]).norm()).fold(0.0f64, f64::max)
}

#[derive(Clone, Copy)]
enum Arm {
    Az,
    Heggie,
    Logh(LhTime, Stepper),
}

/// Returns the final state, whether it stayed finite, and the **force evaluations** spent.
/// Steps are not comparable across these -- RK4 spends four per step and KDK one.
fn march(
    arm: Arm, s: Cart<f64>, m: &[f64; 3], t_max: f64, n_sync: usize, eta: f64, ms: usize,
) -> (Cart<f64>, bool, u64) {
    match arm {
        Arm::Az => {
            let o = az::integrate_az_opts(
                s, m, t_max, n_sync, eta, ms,
                &AzOpts {
                    r_coll_frac: 0.0,
                    stop_on_event: false,
                    step_limit: StepLimit::Predictive,
                    step_limit_f: 0.02,
                    ..Default::default()
                },
            );
            (o.state, o.finite && !o.budget_exhausted, o.steps as u64 * 4)
        }
        Arm::Heggie => {
            let o = heggie::integrate_hg(
                s, m, t_max, n_sync, eta, ms,
                &HgOpts { r_coll_frac: 0.0, stop_on_event: false, ..Default::default() },
            );
            (o.state, o.finite && !o.budget_exhausted, o.steps as u64 * 4)
        }
        Arm::Logh(time, stepper) => {
            let o = integrate_lh(
                s, m, t_max, n_sync, eta, ms,
                &LhOpts {
                    time,
                    stepper,
                    r_coll_frac: 0.0,
                    stop_on_event: false,
                    step_limit_f: 0.02,
                    ..Default::default()
                },
            );
            (o.state, o.finite && !o.budget_exhausted, o.force_evals as u64)
        }
    }
}

fn main() {
    let res: usize = arg(1, 64);
    let case: String = arg(2, "near-field".to_string());
    let max_steps: usize = arg(3, 400_000);
    let cfg0 = EnsembleCfg::production();

    let (chart, cx, cy, half, body, t_max) = if case == "config_stability" {
        let (c, x, y, h) = Chart::config_stability();
        (c, x, y, h, 0usize, 50.0f64)
    } else {
        let s = grid::region(&case.replace('_', " "), 4, 4, 0.05).expect("unknown case");
        (s.chart, s.cx, s.cy, s.half, s.body, 13.0f64)
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
    let n_sync = (t_max / 0.4).round().max(4.0) as usize;

    println!("STAGE 1: reversibility on {case} at {res}^2, t_max = {t_max}, n_sync = {n_sync}.");
    println!("  Termination OFF, r_coll = 0, nominal copy only, COM-centred residual.\n");

    // --- the FTLE field, once, shared by every arm --------------------------------------------
    let ftle_cfg = EnsembleCfg::production().with_overrides(&[
        Override::TMax(t_max),
        Override::NSync(n_sync),
        Override::RCollFrac(0.0),
        Override::StopOnEvent(false),
        Override::RefineFlagged(false),
        Override::Integrator(Integrator::Az),
        Override::MaxSteps(max_steps),
        Override::Ftle(Some(Default::default())),
    ]);
    let fpx: Vec<PixelOut> =
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ftle_cfg)).collect();
    let ftle: Vec<f64> = fpx.iter().map(|p| p.ftle).collect();
    let renorm_ok = fpx.iter().filter(|p| p.ftle_renorm > 0).count();
    let drift_ref: Vec<f64> = fpx.iter().map(|p| p.energy_drift_max).collect();
    println!(
        "  FTLE field: {renorm_ok} of {} pixels completed at least one renormalisation. \
         **Read this first** -- `ftle` is meaningless without it.",
        fpx.len()
    );
    println!("  FTLE p10 {:.4}  p50 {:.4}  p90 {:.4}\n", q(&ftle, 0.10), q(&ftle, 0.50), q(&ftle, 0.90));

    let ics: Vec<Ic<f64>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, 0, cfg0.jitter_frac, cfg0.seed, cfg0.jitter_scheme, cfg0.decode_path,
            )[0]
        })
        .collect();

    println!(
        "  {:>10} {:>11} {:>11} {:>11} {:>11} {:>10} {:>7} {:>7}",
        "arm", "rev p50", "rev p90", "rev/amp p50", "disp p50", "evals p50", "nonfin", "secs"
    );

    let arms: [(&str, Arm, f64); 5] = [
        ("az", Arm::Az, 1.0),
        ("heggie", Arm::Heggie, 1.0),
        ("logh_rk4", Arm::Logh(LhTime::LogH, Stepper::Rk4), 1.0),
        // The leapfrog arms spend one force evaluation per step where RK4 spends four, so they run
        // at eta/4 for a nominal evaluation match. Read each row at the evaluations it spent.
        ("logh_lf", Arm::Logh(LhTime::LogH, Stepper::Kdk), 0.25),
        ("logh_gbs", Arm::Logh(LhTime::LogH, Stepper::Gbs), 0.25),
    ];

    let mut corr_rows: Vec<(String, f64, usize, f64, usize, f64, usize)> = Vec::new();

    for (label, arm, sc) in arms {
        let eta = cfg0.eta * sc;
        let t0 = std::time::Instant::now();
        let out: Vec<(f64, f64, f64, bool)> = ics
            .par_iter()
            .map(|ic| {
                let r0 = energy::hyperradius(&ic.s.r, &ic.m);
                let (s1, f1, e1) = march(arm, ic.s, &ic.m, t_max, n_sync, eta, max_steps);
                let (s2, f2, e2) = march(arm, flip(s1), &ic.m, t_max, n_sync, eta, max_steps);
                let rev = sep(&s2, &ic.s, &ic.m) / r0;
                let disp = sep(&s1, &ic.s, &ic.m) / r0;
                (rev, disp, (e1 + e2) as f64, f1 && f2)
            })
            .collect();
        let secs = t0.elapsed().as_secs_f64();

        let rev: Vec<f64> = out.iter().map(|x| x.0).collect();
        let disp: Vec<f64> = out.iter().map(|x| x.1).collect();
        let ev: Vec<f64> = out.iter().map(|x| x.2).collect();
        let nonfin = out.iter().filter(|x| !x.3).count();
        // Round-off per unit amplification. `exp(lambda t)` overflows for a large FTLE, so this is
        // formed in logs and exponentiated once -- a ratio of two huge numbers is not the quantity.
        let rev_amp: Vec<f64> = rev
            .iter()
            .zip(ftle.iter())
            .map(|(r, l)| if l.is_finite() { (r.ln() - l * t_max).exp() } else { f64::NAN })
            .collect();

        println!(
            "  {label:>10} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>10.3e} {nonfin:>7} {secs:>7.1}",
            q(&rev, 0.50), q(&rev, 0.90), q(&rev_amp, 0.50), q(&disp, 0.50), q(&ev, 0.50)
        );

        let (c_raw, n_raw) = spearman(&rev, &ftle);
        let (c_amp, n_amp) = spearman(&rev_amp, &ftle);
        let (c_ctl, n_ctl) = spearman(&rev, &shift_half(&ftle, res));
        corr_rows.push((label.into(), c_raw, n_raw, c_amp, n_amp, c_ctl, n_ctl));
    }

    // --- the correlations, with the drift row for comparison ----------------------------------
    let (d_raw, d_n) = spearman(&drift_ref, &ftle);
    let (d_ctl, _) = spearman(&drift_ref, &shift_half(&ftle, res));
    println!("\n  SPEARMAN vs FTLE. `shifted` is the control: same marginals, no spatial relation.");
    println!("  A coefficient no larger than its own shifted control is a NULL, not a correlation.");
    println!(
        "  {:>10} {:>10} {:>10} {:>10} {:>8}",
        "arm", "rev", "rev/amp", "shifted", "n"
    );
    for (l, cr, nr, ca, _na, cc, _nc) in &corr_rows {
        println!("  {l:>10} {cr:>10.4} {ca:>10.4} {cc:>10.4} {nr:>8}");
    }
    println!("  {:>10} {d_raw:>10.4} {:>10} {d_ctl:>10.4} {d_n:>8}   <- energy drift, AZ", "drift", "--");

    println!(
        "\nHOW TO READ THIS\n\n\
         **PREDICTION 1, recorded before the run: reversibility tracks FTLE more strongly than\n\
         drift does.** Drift-vs-FTLE is a published null for AZ (-0.0820 against a shifted\n\
         -0.1022).\n\n\
         **PREDICTION 1b: `rev/amp` should track FTLE MUCH LESS than `rev` does.** That is the\n\
         test of the normalisation itself. If the two columns correlate alike, the division did\n\
         not remove the amplification and `rev/amp` is not integrator quality -- it is `rev` under\n\
         another name, and must not be reported as the integrator arm.\n\n\
         **Compare arms DOWN a column only.** RK4 is not time-symmetric and KDK is; across\n\
         steppers this statistic scores symmetry rather than accuracy."
    );
}
