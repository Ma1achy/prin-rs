//! **A per-pixel error measure that needs no reference trajectory.**
//!
//! March to `t_max`, negate the velocities, march the same span again, negate back, and compare
//! against the initial condition. In exact arithmetic with a time-symmetric method the return is
//! the identity, so the residual is a pure error — and unlike every instrument this project
//! currently has, it needs no second trajectory to compare against.
//!
//! The gap it is aimed at is on record. `eta/256` as a ground truth came back **saturated**:
//! chord 2.000, antipodal, for a correct mode and a broken one alike. Energy drift is
//! **documented blind** to a displacement in time — the overshoot clamp bought 24,000x on the
//! figure-eight while moving `near-field`'s median drift 37x the wrong way. `error_ratio` works
//! but needs the whole ensemble. This needs one trajectory and one extra march.
//!
//! # What it actually measures, stated before any number
//!
//! Not integrator error alone. Three things enter, and none of them can be removed:
//!
//!   1. **Round-off, amplified by the Lyapunov exponent.** This is the dominant term on a chaotic
//!      field and it is the subject of Portegies Zwart & Boekholt's irreversibility result. It is
//!      a property of the *problem*, not of the integrator.
//!   2. **Stepper asymmetry.** KDK is time-symmetric and retraces exactly in exact arithmetic;
//!      **classical RK4 is not**, and neither is GBS over an adaptive base. So a reversibility
//!      number compared **across steppers scores time-symmetry, not accuracy**. Read it down a
//!      column — AZ, Heggie and `logh_rk4` are all RK4 and are directly comparable; `logh_lf` and
//!      `logh_gbs` are not comparable to them on this statistic.
//!   3. **Step-control asymmetry.** The reverse leg chooses its own step sequence from its own
//!      states and does not mirror the forward one, so even a symmetric stepper does not retrace
//!      exactly under adaptive `eta`. This is real and is not a defect of the harness.
//!
//! Same discipline as `rho/gamma` in `logh_arms`: one column, one meaning, never across.
//!
//! # The two guards, and what makes each of them fire
//!
//! **`disp` is the control.** A trajectory that barely moved returns to its start trivially, and
//! a small residual would read as a triumph. `disp` is the forward displacement over the same
//! normalisation, so `rev/disp` is the honest ratio and a row with `disp` near zero is a row with
//! no subject. *A difference can be small because both sides are right or because both are dead.*
//!
//! **The `eta` ladder separates the two terms that matter.** If the residual **falls** with `eta`
//! it is truncation and the integrator owns it. If it **saturates** it is round-off amplified by
//! the flow, and no step size buys it off — which is the answer for a chaotic pixel and is a
//! measurement rather than a failure. A single `eta` cannot tell those apart, so this never runs
//! at one.
//!
//! Termination is **off** and `r_coll = 0` throughout: a run stopped by an event has legs of
//! different lengths and the comparison is meaningless. That flag alone produced five wrong
//! conclusions in a row once.
//!
//! Args: `res case max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts, StepLimit};
use prin_rs::integrate::heggie::{self, HgOpts};
use prin_rs::integrate::logh::{integrate_lh, LhOpts, LhTime, Stepper};
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

/// Velocities negated. The reverse leg is a *forward* march of this state: every time
/// transformation in this codebase gives `dt > 0`, so running time backwards is done in the
/// state, never in the step.
fn flip(s: Cart<f64>) -> Cart<f64> {
    let mut o = s;
    for i in 0..3 {
        o.v[i] = -o.v[i];
    }
    o
}

/// Positions with the centre of mass removed.
///
/// **This is not tidiness; without it the diagnostic measures the frame.** AZ and Heggie
/// reconstruct Cartesian positions from relative coordinates, which places the COM at the origin.
/// logH integrates absolute positions and leaves the COM where the decode put it. On a slice that
/// moves one body — every `body_plane` region does — the decoded configuration is not COM-centred,
/// so the two families sit a constant translation apart: measured `(-0.0125, +2.4875)` on `far`,
/// identical for all three bodies.
///
/// Every other comparison in this project is translation-invariant — energy drift, `shape_vec`,
/// outcome labels — so this has never mattered before. A reversibility residual on absolute
/// positions is the first quantity here that is not, and read raw it reports the COM offset as a
/// flat, `eta`-independent error: `4.016e-1` on `far`, unchanged over a 16x refinement.
fn centred(s: &Cart<f64>, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let mt = m[0] + m[1] + m[2];
    let com = (s.r[0] * m[0] + s.r[1] * m[1] + s.r[2] * m[2]) / mt;
    [s.r[0] - com, s.r[1] - com, s.r[2] - com]
}

/// `max_i |a_i - b_i|` over COM-centred positions, the sup norm over bodies. Not an RMS: one body
/// returning to the wrong place is the failure, and averaging it against two that did not is how
/// it would be hidden.
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

fn march(arm: Arm, s: Cart<f64>, m: &[f64; 3], t_max: f64, n_sync: usize, eta: f64, ms: usize)
    -> (Cart<f64>, bool)
{
    match arm {
        Arm::Az => {
            let o = az::integrate_az_opts(
                s, m, t_max, n_sync, eta, ms,
                // **Named, not inherited.** `AzOpts::default()` carries `StepLimit::None` while
                // `HgOpts::default()` carries `step_limit_f = 0.02`; taking both defaults would
                // run the two arms under different step control and call it a comparison.
                &AzOpts {
                    r_coll_frac: 0.0,
                    stop_on_event: false,
                    step_limit: StepLimit::Predictive,
                    step_limit_f: 0.02,
                    ..Default::default()
                },
            );
            (o.state, o.finite)
        }
        Arm::Heggie => {
            let o = heggie::integrate_hg(
                s, m, t_max, n_sync, eta, ms,
                &HgOpts { r_coll_frac: 0.0, stop_on_event: false, ..Default::default() },
            );
            (o.state, o.finite)
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
            (o.state, o.finite)
        }
    }
}

fn main() {
    let res: usize = arg(1, 64);
    let case: String = arg(2, "near-field".to_string());
    let max_steps: usize = arg(3, 400_000);
    let cfg = EnsembleCfg::production();

    let (chart, cx, cy, half, body, t_max) = if case == "config_stability" {
        let (c, x, y, h) = Chart::config_stability();
        (c, x, y, h, 0usize, 50.0f64)
    } else {
        let s = grid::region(&case.replace('_', " "), 4, 4, 0.05).expect("unknown case");
        (s.chart, s.cx, s.cy, s.half, s.body, 13.0)
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
    let n_sync = (t_max / 0.4).round().max(4.0) as usize;

    // Nominal copy only. Reversibility is a per-trajectory property; the ensemble is what
    // `error_ratio` is for, and spending `E+1` marches here would measure the same thing eight
    // times.
    let ics: Vec<Ic<f64>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, 0, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )[0]
        })
        .collect();

    println!(
        "{case} at {res}^2, t_max = {t_max}, n_sync = {n_sync}, termination OFF, r_coll = 0, \
         nominal copy only.\n"
    );
    println!(
        "  **Read this DOWN a column, never across.** RK4 is not time-symmetric and KDK is, so a\n  \
         reversibility number compared between steppers scores symmetry rather than accuracy.\n  \
         `az`, `heggie` and `logh_rk4` are all RK4 and are mutually comparable; the leapfrog and\n  \
         GBS arms are not comparable to them on this statistic.\n"
    );
    println!(
        "  **`disp` is the guard.** It is the forward displacement in the same units. A row whose\n  \
         `disp` is near zero has no subject, and `rev/disp` is the honest ratio.\n"
    );
    println!(
        "  **The `eta` ladder is the discriminator.** Falling with `eta` is truncation and the\n  \
         integrator owns it; saturating is round-off amplified by the flow and no step size buys\n  \
         it off. One `eta` cannot tell those apart.\n"
    );

    println!(
        "  {:>10} {:>9} {:>11} {:>11} {:>11} {:>11} {:>8} {:>8}",
        "arm", "eta", "rev p50", "rev p90", "disp p50", "rev/disp", "nonfin", "secs"
    );

    let arms: [(&str, Arm); 5] = [
        ("az", Arm::Az),
        ("heggie", Arm::Heggie),
        ("logh_rk4", Arm::Logh(LhTime::LogH, Stepper::Rk4)),
        ("logh_lf", Arm::Logh(LhTime::LogH, Stepper::Kdk)),
        ("logh_gbs", Arm::Logh(LhTime::LogH, Stepper::Gbs)),
    ];

    for (label, arm) in arms {
        for mult in [1.0, 0.25, 0.0625] {
            // The leapfrog arms spend one force evaluation per step where RK4 spends four, so
            // they run at eta/4 for a nominal evaluation match -- the same convention as
            // `logh_arms`, carried here so a row is not cheaper than the row above it.
            let sc = if matches!(arm, Arm::Logh(_, Stepper::Kdk) | Arm::Logh(_, Stepper::Gbs)) {
                0.25
            } else {
                1.0
            };
            let eta = cfg.eta * sc * mult;
            let t0 = std::time::Instant::now();
            let out: Vec<(f64, f64, bool)> = ics
                .par_iter()
                .map(|ic| {
                    let r0 = energy::hyperradius(&ic.s.r, &ic.m);
                    let (s1, f1) = march(arm, ic.s, &ic.m, t_max, n_sync, eta, max_steps);
                    let (s2, f2) = march(arm, flip(s1), &ic.m, t_max, n_sync, eta, max_steps);
                    // `s2` is the returned state with velocities still reversed; positions are
                    // what is compared, so the second flip is unnecessary and is not done.
                    let rev = sep(&s2, &ic.s, &ic.m) / r0;
                    let disp = sep(&s1, &ic.s, &ic.m) / r0;
                    (rev, disp, f1 && f2)
                })
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let rev: Vec<f64> = out.iter().map(|x| x.0).collect();
            let disp: Vec<f64> = out.iter().map(|x| x.1).collect();
            let ratio: Vec<f64> = out.iter().map(|x| x.0 / x.1).collect();
            let nonfin = out.iter().filter(|x| !x.2).count();
            println!(
                "  {label:>10} {eta:>9.2e} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {nonfin:>8} {secs:>8.1}",
                q(&rev, 0.50), q(&rev, 0.90), q(&disp, 0.50), q(&ratio, 0.50)
            );
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         **A saturating column is the expected result on a chaotic region, not a failure.** It\n\
         says the residual is round-off amplified by the flow, which is a property of the field.\n\
         A column that falls with `eta` is one where truncation still dominates and the\n\
         integrator has room.\n\n\
         **`rev/disp` near 1 means the trajectory has lost its way entirely** -- the return is as\n\
         far from the start as the forward leg went. Near 0 means the march retraced.\n\n\
         **What this cannot do is rank steppers.** Time-symmetry is worth orders here and is not\n\
         accuracy. The comparison this supports is regularisation at fixed stepper: `az` against\n\
         `heggie` against `logh_rk4`, all three RK4."
    );
}
