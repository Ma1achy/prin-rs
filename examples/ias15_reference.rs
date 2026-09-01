//! **Score the integrators against a machine-precision reference.**
//!
//! This project has never had a reference. `eta/256` was tried and came back **saturated** — chord
//! 2.000, antipodal, for a correct mode and a broken one alike — so it could not separate the arms
//! it was asked to rank. IAS15 (Rein & Spiegel 2015) conserves energy to `5.2e-15` over `t = 200`
//! with adaptive stepping, which is what makes it usable as ground truth.
//!
//! # The reference must be shown to be converged, or this repeats the failure it exists to fix
//!
//! A reference that is not itself converged scores its own error, not the arms'. So the reference
//! is computed at **two tolerances** and the disagreement between them is printed **above** the
//! table. If the two references differ by anything approaching the arm errors, the row says so and
//! the numbers below it are not readable.
//!
//! # And the horizon is swept, because saturation is the expected failure
//!
//! Chaotic divergence eventually carries every arm to the diameter of the shape sphere, at which
//! point the comparison is between two numbers that mean nothing. Rather than pick a horizon and
//! hope, several are run and each is labelled **LIVE** or **SATURATED** by the fraction of pixels
//! whose error has reached the scale of the configuration itself.
//!
//! Termination is off and `r_coll = 0` throughout: a run stopped by an event is parked at a close
//! approach, and every arm must integrate the same span for the comparison to mean anything.
//!
//! Args: `res case`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts, StepLimit};
use prin_rs::integrate::heggie::{self, HgOpts};
use prin_rs::integrate::ias15;
use prin_rs::integrate::logh::chain::{rk4 as chain_rk4, ChainOrder, ChainState};
use prin_rs::integrate::logh::hamiltonian::LhTime;
use prin_rs::integrate::logh::{integrate_lh, LhOpts, Stepper};
use prin_rs::physics::{energy, Cart, Ic};
use prin_rs::Vec2;

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

fn com_centre(c: &Cart<f64>, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let mt = m[0] + m[1] + m[2];
    let com = (c.r[0] * m[0] + c.r[1] * m[1] + c.r[2] * m[2]) / mt;
    [c.r[0] - com, c.r[1] - com, c.r[2] - com]
}

/// Sup-norm position difference over bodies, COM-removed and normalised by the initial
/// hyperradius so it is scale-invariant. **Not an RMS**: one body in the wrong place is the
/// failure, and averaging it against two that are right is how it would be hidden.
fn err(a: &Cart<f64>, b: &Cart<f64>, m: &[f64; 3], r0: f64) -> f64 {
    let (x, y) = (com_centre(a, m), com_centre(b, m));
    (0..3).map(|i| (x[i] - y[i]).norm()).fold(0.0f64, f64::max) / r0
}

/// The reference: IAS15, adaptive, to `t_max`. Returns the state, the force evaluations spent,
/// and **whether it actually reached `t_max`**.
///
/// That third value is not bookkeeping. The step budget used to be a silent `return`, so a
/// reference that ran out returned wherever it had got to and the harness compared a finished
/// trajectory against an unfinished one -- measured self-gap `1.510e-1` at `t = 1`, thirteen
/// orders worse than the same run at a looser tolerance, which reads as a tolerance catastrophe
/// and is a truncation. *An unfinished run reported as a value* is the failure this file has now
/// produced three times in different places; the only fix that holds is to return the fact.
fn reference(c: &Cart<f64>, m: &[f64; 3], t_max: f64, eps: f64) -> (Cart<f64>, u64, bool) {
    let (mut r, mut v) = (c.r, c.v);
    let mut t = 0.0f64;
    let mut dt = 1e-3f64;
    let mut evals = 0u64;
    let mut guard = 0usize;
    while t < t_max && guard < 4_000_000 {
        guard += 1;
        let step = dt.min(t_max - t);
        // 24 corrector iterations, not 12: at a close encounter the predictor-corrector needs
        // more, and capping it low is an *advance-anyway* site -- the step is taken with an
        // unconverged `b`, silently, which is what a reference must never do.
        let (out, e, _, b_last) = ias15::step(&r, &v, m, step, 24, 1e-16);
        if !out.r.iter().all(|x| x.is_finite()) {
            return (Cart { r, v }, evals, false);
        }
        r = out.r;
        v = out.v;
        t += step;
        evals += e as u64;
        // Floor lowered to 1e-13. At 1e-9 a close encounter that needs a smaller step cannot get
        // one, so the reference stops resolving the encounter and its two tolerances part company
        // -- measured, the self-gap went 4.2e-16 at t = 1 to 1.9e-10 at t = 2.
        dt = ias15::next_dt(&r, m, dt, b_last, eps).clamp(1e-13, 0.5);
    }
    (Cart { r, v }, evals, t >= t_max - 1e-12 * t_max.max(1.0))
}

fn main() {
    let res: usize = arg(1, 32);
    let case: String = arg(2, "near-field".to_string());
    let cfg = EnsembleCfg::production();

    let (chart, cx, cy, half, body) = if case == "config_stability" {
        let (c, x, y, h) = Chart::config_stability();
        (c, x, y, h, 0usize)
    } else {
        let s = grid::region(&case.replace('_', " "), 4, 4, 0.05).expect("unknown case");
        (s.chart, s.cx, s.cy, s.half, s.body)
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);

    println!("Scoring against an IAS15 reference. {case} at {res}^2, nominal copy, termination OFF.\n");
    println!(
        "  **The reference is checked against itself first.** Two tolerances, and the row prints\n  \
         their disagreement. A reference that is not converged scores its own error -- which is\n  \
         exactly how `eta/256` failed, coming back saturated at chord 2.000 for a correct mode and\n  \
         a broken one alike.\n"
    );

    let ics: Vec<Ic<f64>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, 0, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )[0]
        })
        .collect();

    for t_max in [1.0f64, 2.0, 4.0] {
        let n_sync = (t_max / 0.4).round().max(4.0) as usize;
        let rows: Vec<(f64, f64, [f64; 4], [u64; 4], u64)> = ics
            .par_iter()
            .map(|ic| {
                let m = ic.m;
                let r0 = energy::hyperradius(&ic.s.r, &m);
                // **(1e-11, 1e-12), and the upper bound is measured rather than chosen.**
                //
                // Tightened from (1e-9, 1e-11), which failed its own guard at t = 2. The first
                // retightening went to (1e-12, 1e-14) and made things far worse -- self-gap
                // `1.510e-1`, which reads as a tolerance catastrophe and is nothing of the kind:
                // at `eps <= 1e-13` the step control drives `dt` onto its `1e-13` floor and burns
                // all 4,000,000 steps covering **0.4% of the span**. Measured, at t_max = 1:
                //
                //     eps 1e-11   36 steps      completed
                //     eps 1e-12   96 steps      completed
                //     eps 1e-13   4,000,000     BUDGET, reached t = 0.0044
                //     eps 1e-14   4,000,000     BUDGET, reached t = 0.0023
                //
                // So this is the tightest pair that both completes and still differs in work.
                let (ref_a, ev_a, ok_a) = reference(&ic.s, &m, t_max, 1e-11);
                let (ref_b, ev_b, ok_b) = reference(&ic.s, &m, t_max, 1e-12);
                // A reference that did not finish is not a reference. `NaN` rather than a number,
                // so it is dropped from the quantiles instead of poisoning them, and counted.
                let self_gap =
                    if ok_a && ok_b { err(&ref_a, &ref_b, &m, r0) } else { f64::NAN };
                // **Both arms' costs are carried.** The previous cut recorded only `ev_a`, and I
                // then used it to argue the *other* arm had not run out of budget -- diagnosing one
                // arm by reading the other's cost. `ev_b` exists so that cannot recur.
                if std::env::var("IAS_DEBUG").is_ok() {
                    eprintln!("ok_a={ok_a} ok_b={ok_b} ev_a={ev_a} ev_b={ev_b}");
                }

                let base_az = AzOpts {
                    r_coll_frac: 0.0,
                    stop_on_event: false,
                    step_limit: StepLimit::Predictive,
                    step_limit_f: 0.02,
                    ..Default::default()
                };
                let a = az::integrate_az_opts(ic.s, &m, t_max, n_sync, cfg.eta, 4_000_000, &base_az);
                let h = heggie::integrate_hg(
                    ic.s, &m, t_max, n_sync, cfg.eta, 4_000_000,
                    &HgOpts { r_coll_frac: 0.0, stop_on_event: false, ..Default::default() },
                );
                let l = integrate_lh(
                    ic.s, &m, t_max, n_sync, cfg.eta, 4_000_000,
                    &LhOpts {
                        time: LhTime::LogH,
                        stepper: Stepper::Rk4,
                        r_coll_frac: 0.0,
                        stop_on_event: false,
                        step_limit_f: 0.02,
                        ..Default::default()
                    },
                );
                // Chain, fixed fictitious step, matched in force evaluations to the logH arm.
                let o = ChainOrder::select(&ic.s.r);
                let mut cs = ChainState::from_cart(&ic.s, o);
                let bc = energy::potential_pos(&ic.s.r, &m, 0.0) - energy::kinetic(&ic.s.v, &m);
                // **March to the CLOCK, not to a step count.** `ChainState` carries its own
                // physical time, and under LogH the physical time per fictitious step varies along
                // the trajectory -- so a fixed step count lands wherever it lands. Two earlier
                // cuts of this harness matched a step count instead and reported the shortfall as
                // an *error*: `1.509e-1`, then `1.363e-2` after one repair, both of them a
                // trajectory that had not finished rather than one that had gone wrong. The
                // clock removes the guess.
                let u0 = energy::potential_pos(&ic.s.r, &m, 0.0);
                let hstep = cfg.eta * t_max * u0 / (t_max / cfg.eta).max(1.0);
                let mut ok = true;
                let mut nsteps = 0usize;
                while cs.t < t_max && nsteps < 4_000_000 {
                    // Land exactly on the horizon: the last step is shortened in fictitious time
                    // by the ratio the clock still needs, which is the same first-order landing the
                    // other occupants use.
                    let dtds = 1.0 / (energy::kinetic(&cs.to_cart(&m, o).v, &m) + bc).max(1e-300);
                    let want = (t_max - cs.t) / dtds;
                    let h = hstep.min(want.max(0.0));
                    if !(h > 0.0) {
                        break;
                    }
                    let (n, _) = chain_rk4(&m, &cs, o, bc, LhTime::LogH, h);
                    cs = n;
                    nsteps += 1;
                    if !cs.is_finite() {
                        ok = false;
                        break;
                    }
                }
                // The guard that says the span actually matched. A chain arm that stopped short is
                // not a chain result, and this is what the two earlier cuts lacked.
                if (cs.t - t_max).abs() > 1e-6 * t_max {
                    ok = false;
                }
                let cc = cs.to_cart(&m, o);

                let e = [
                    err(&a.state, &ref_a, &m, r0),
                    err(&h.state, &ref_a, &m, r0),
                    err(&l.state, &ref_a, &m, r0),
                    if ok { err(&cc, &ref_a, &m, r0) } else { f64::INFINITY },
                ];
                let ev = [
                    a.steps as u64 * 4,
                    h.steps as u64 * 4,
                    l.force_evals as u64,
                    nsteps as u64 * 4,
                ];
                (self_gap, r0, e, ev, ev_a.max(ev_b))
            })
            .collect();

        let gaps: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let incomplete = gaps.iter().filter(|x| !x.is_finite()).count();
        let names = ["az", "heggie", "logh_rk4", "chain"];
        let errs: Vec<Vec<f64>> =
            (0..4).map(|k| rows.iter().map(|r| r.2[k]).collect()).collect();
        // Saturated when the error has reached the scale of the configuration itself.
        let sat = errs
            .iter()
            .map(|e| e.iter().filter(|x| **x > 0.5).count() as f64 / e.len() as f64)
            .fold(0.0f64, f64::max);
        let ref_ev = q(&rows.iter().map(|r| r.4 as f64).collect::<Vec<_>>(), 0.50);

        let verdict = if incomplete * 2 > gaps.len() {
            "REFERENCE DID NOT FINISH on over half the pixels -- it ran out of step budget rather \
             than disagreeing, and a truncated reference reads as a tolerance catastrophe"
        } else if q(&gaps, 0.90) > q(&errs[1], 0.10) {
            "REFERENCE NOT CONVERGED -- its own two tolerances disagree by more than the best arm's\n     \
             error, so nothing below is readable"
        } else if sat > 0.5 {
            "SATURATED -- over half the pixels have diverged to the configuration scale"
        } else {
            "LIVE"
        };
        println!(
            "  t_max {t_max:>4.1}   ref self-gap p50 {:.3e} p90 {:.3e}   incomplete {incomplete}/{}   ref evals p50 {ref_ev:.2e}   {verdict}",
            q(&gaps, 0.50), q(&gaps, 0.90), gaps.len()
        );
        for (k, n) in names.iter().enumerate() {
            println!(
                "      {n:>9}  err p50 {:>10.3e}  p90 {:>10.3e}  evals p50 {:>9.2e}",
                q(&errs[k], 0.50),
                q(&errs[k], 0.90),
                q(&rows.iter().map(|r| r.3[k] as f64).collect::<Vec<_>>(), 0.50)
            );
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         **Read the verdict on each horizon before the arms.** A reference whose two tolerances\n\
         disagree is scoring itself; a saturated horizon is comparing two meaningless numbers.\n\
         Both failures have happened in this project already.\n\n\
         **Errors are against TRUTH, not against each other.** That is what this reference buys and\n\
         what a difference between two arms could never say: `smaller` here means closer to the\n\
         real trajectory, not merely closer to some other integrator's answer.\n\n\
         **Read the evaluation counts beside the errors.** The arms are matched at the same `eta`,\n\
         not at matched cost, and an arm spending twice the evaluations for a better number has\n\
         made a different claim."
    );
}
