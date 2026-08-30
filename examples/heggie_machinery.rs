//! **The falsifiable test. Does removing re-registration remove the sensitivity?**
//!
//! `examples/az_machinery.rs` established the mechanism behind the `config_stability` wedges:
//! doubling the sync cadence at **fixed step size** moves AZ's drift field by 0.444 decades,
//! 6000x the effect of changing *which* reference body is chosen. Not which chart is chosen — how
//! often the state is passed through one.
//!
//! Heggie's global regularisation has no reference body and therefore no re-registration at all:
//! `HgSystem` depends on nothing but the masses, and `to_reg` runs once at `t = 0`. So the
//! prediction is specific: **Heggie's field must be insensitive to `n_sync` at fixed step size.**
//!
//! **The null is the point.** If Heggie also moves ~0.444 decades the mechanism was misidentified,
//! the conclusion on record is wrong, and this port was built on it. That is a result worth
//! having, and it is why this runs before anything is plumbed into the ensemble layer.
//!
//! # Both integrators in one harness, on the same copies
//!
//! The AZ arms are re-run here rather than quoted from the committed table. `pixel::evaluate`'s
//! `energy_drift_max` and `HgOut::drift` are the same definition, but "the same definition" is an
//! argument and this makes it an identity: one slice, one set of jittered copies, one statistic,
//! and the only thing that differs between the two blocks is which integrator marched them.
//!
//! # The confounds, both of them, from the AZ harness's own record
//!
//! `dt ~ eta * t_max / n_sync`, so raising `n_sync` at fixed `eta` also shrinks every step, and
//! "more boundaries" becomes inseparable from "finer stepping". **`eta` is scaled with `n_sync`**
//! and `steps p50` is printed as the check that it worked. A **deliberately confounded arm** is
//! included so the controlled rows are a demonstration rather than a set of numbers.
//!
//! The escape-window confound does not apply here: termination is **off in both integrators**, so
//! every trajectory runs to `t_max` and the comparison is of the marching kernel alone. That is
//! also why `r_coll` is zero — a collision stop would make `t_end` differ between arms and the
//! drift would be reporting when each run stopped.
//!
//! # The control arm
//!
//! `refresh_h_at_boundary` re-derives the frozen energy at each boundary. It is the one thing the
//! Heggie driver can do that reintroduces a boundary-dependent quantity into an otherwise
//! uninterrupted march. **A control that changes nothing proves nothing**, so its effect is
//! reported: if it moves the field, the harness can see boundary sensitivity and a Heggie null is
//! informative; if it moves nothing, the null might be the harness.
//!
//! Args: `res root max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::physics::Ic;

const WINDOW: f64 = 0.4;
const T_MAX: f64 = 50.0;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;
/// **Not the production 30,000.** At that budget Heggie exhausts it on 27% of the frame — every
/// one of them a budget exhaustion and **not a divergence**, since it needs ~22% more steps than
/// AZ for the same trajectory. Those pixels would then be pinned at the ramp cap in *every* arm
/// and contribute exactly zero to the chord, so the headline number would be small because a
/// quarter of the field was dead rather than because the field did not move. Raised until the
/// `budget` column reads zero, and the column is printed so the reader can see that it does.
const MAX_STEPS_DEFAULT: usize = 400_000;

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

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.015], [0.34, 0.06, 0.43], [0.72, 0.21, 0.33],
        [0.98, 0.55, 0.04], [0.99, 1.0, 0.64],
    ];
    let t = x.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    std::array::from_fn(|k| (255.0 * (S[i][k] * (1.0 - f) + S[i + 1][k] * f)).clamp(0.0, 255.0) as u8)
}

/// Which integrator, and what varies. `refresh` is Heggie's control arm.
#[derive(Clone, Copy)]
enum Arm {
    Az,
    Hg { refresh: bool },
}

/// Max drift over the `E+1` copies of one pixel, and the summed step count.
///
/// **Max, not median.** That is `PixelOut::energy_drift_max`'s own reduction, and it is the one
/// that tracks damage; a median over eight copies hides the wild one, which is the whole reason
/// `error_ratio` aggregates by max.
fn pixel(cs: &[Ic<f64>], arm: Arm, n_sync: usize, eta: f64, max_steps: usize) -> (f64, f64, bool) {
    let (mut d, mut st) = (0.0f64, 0.0f64);
    // **Budget exhaustion is separated from divergence.** "Heggie fails on a quarter of the
    // frame" and "Heggie wants a larger step budget on a quarter of the frame" are different
    // claims, and a single `nonfin` column cannot tell them apart.
    let mut budget = false;
    for c in cs {
        let (drift, steps, finite) = match arm {
            Arm::Az => {
                let o = az::integrate_az_opts(
                    c.s, &c.m, T_MAX, n_sync, eta, max_steps,
                    &AzOpts {
                        // Termination off: the kernel is what is being compared, and a stop would
                        // make `t_end` differ between arms so the drift would report when each
                        // run stopped rather than how it marched.
                        stop_on_event: false,
                        r_coll_frac: 0.0,
                        // **Step control held constant across the two integrators.** This is
                        // production for AZ and the default for Heggie; `AzOpts::default()` is
                        // `StepLimit::None`, so taking the default here would have Heggie paying
                        // for a limit AZ was not paying for, and the step counts would be
                        // comparing step control rather than regularisation.
                        step_limit: az::StepLimit::Predictive,
                        step_limit_f: 0.02,
                        ..Default::default()
                    },
                );
                budget |= o.budget_exhausted;
                (o.drift, o.steps as f64, o.finite)
            }
            Arm::Hg { refresh } => {
                let o = integrate_hg(
                    c.s, &c.m, T_MAX, n_sync, eta, max_steps,
                    &HgOpts { refresh_h_at_boundary: refresh, ..Default::default() },
                );
                budget |= o.budget_exhausted;
                (o.drift, o.steps as f64, o.finite)
            }
        };
        st += steps;
        // A non-finite copy is a measurement outcome, never discarded: it carries the top of the
        // ramp so it cannot quietly lower a max.
        d = d.max(if finite && drift.is_finite() { drift } else { f64::INFINITY });
    }
    (d, st, budget)
}

fn main() {
    let res: usize = arg(1, 256);
    let max_steps: usize = arg(3, MAX_STEPS_DEFAULT);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/heggie_machinery");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let n0 = (T_MAX / WINDOW).round() as usize; // 125, the production cadence
    let cfg = EnsembleCfg::production();
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let eta0 = cfg.eta;

    println!("config_stability {res}^2, termination OFF in both integrators.");
    println!("t_max={T_MAX} eta={eta0} max_steps={max_steps} copies={}\n", cfg.n_extra + 1);

    // Decode once. Both integrators march the SAME copies, so nothing between the chart and the
    // march can differ between the two blocks.
    let copies: Vec<Vec<Ic<f64>>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )
        })
        .collect();

    let arms: Vec<(String, Arm, usize, f64)> = vec![
        ("AZ  baseline n=125".into(), Arm::Az, n0, eta0),
        ("AZ  n=250 CONTROLLED".into(), Arm::Az, n0 * 2, eta0 * 2.0),
        ("AZ  n=500 CONTROLLED".into(), Arm::Az, n0 * 4, eta0 * 4.0),
        ("AZ  n=250 confounded".into(), Arm::Az, n0 * 2, eta0),
        ("HG  baseline n=125".into(), Arm::Hg { refresh: false }, n0, eta0),
        ("HG  n=250 CONTROLLED".into(), Arm::Hg { refresh: false }, n0 * 2, eta0 * 2.0),
        ("HG  n=500 CONTROLLED".into(), Arm::Hg { refresh: false }, n0 * 4, eta0 * 4.0),
        ("HG  n=250 confounded".into(), Arm::Hg { refresh: false }, n0 * 2, eta0),
        ("HG  refresh_h n=125".into(), Arm::Hg { refresh: true }, n0, eta0),
        ("HG  refresh_h n=250".into(), Arm::Hg { refresh: true }, n0 * 2, eta0 * 2.0),
    ];

    println!(
        "  {:>22} {:>7} {:>8} {:>11} {:>11} {:>7} {:>7} {:>9} {:>10}",
        "arm", "n_sync", "eta", "steps p50", "drift p50", "nonfin", "budget", "hot lift", "chord p50"
    );

    // Each integrator is compared against ITS OWN baseline. Comparing Heggie against AZ's would
    // be measuring the difference between two integrators, not the sensitivity of either.
    let mut base_az: Option<(Vec<bool>, Vec<f64>)> = None;
    let mut base_hg: Option<(Vec<bool>, Vec<f64>)> = None;

    for (label, arm, ns, eta) in arms {
        let out: Vec<(f64, f64, bool)> =
            copies.par_iter().map(|cs| pixel(cs, arm, ns, eta, max_steps)).collect();
        let n = out.len();

        let lg = |x: f64| if x.is_finite() && x > 0.0 { x.log10() } else { DHI.log10() };
        let d: Vec<f64> = out.iter().map(|(x, _, _)| lg(*x)).collect();
        let mut dr: Vec<f64> = out.iter().map(|(x, _, _)| *x).filter(|x| x.is_finite()).collect();
        let nonfin = n - dr.len();
        let budget = out.iter().filter(|(_, _, b)| *b).count();
        let cut = q(&mut dr.clone(), 0.75);
        let hot: Vec<bool> = out.iter().map(|(x, _, _)| !x.is_finite() || *x > cut).collect();
        let mut st: Vec<f64> = out.iter().map(|(_, s, _)| *s).collect();

        let slot = match arm {
            Arm::Az => &mut base_az,
            Arm::Hg { .. } => &mut base_hg,
        };
        let (mut lift, mut chord) = (f64::NAN, f64::NAN);
        if let Some((h0, d0)) = slot.as_ref() {
            let n_h0 = h0.iter().filter(|x| **x).count().max(1) as f64;
            let bse = hot.iter().filter(|x| **x).count() as f64 / n as f64;
            let p = (0..n).filter(|&i| h0[i] && hot[i]).count() as f64 / n_h0;
            lift = p / bse.max(f64::MIN_POSITIVE);
            let mut c: Vec<f64> =
                (0..n).map(|i| (d[i] - d0[i]).abs()).filter(|x| x.is_finite()).collect();
            chord = q(&mut c, 0.5);
        }

        println!(
            "  {label:>22} {ns:>7} {eta:>8.4} {:>11.3e} {:>11.3e} {nonfin:>7} {budget:>7} {lift:>9.3} {chord:>10.3e}",
            q(&mut st, 0.5),
            q(&mut dr, 0.5),
        );

        let buf: Vec<u8> = out
            .iter()
            .flat_map(|(x, _, _)| {
                if x.is_finite() && *x > 0.0 {
                    ramp((lg(*x) - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    [255, 0, 255]
                }
            })
            .collect();
        let slug = label.replace(' ', "_").replace('=', "");
        let _ = prin_rs::output::adaptive::save_rect(
            &format!("{dir}/drift_{slug}.png"), res, res, &buf,
        );

        if slot.is_none() {
            *slot = Some((hot, d));
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **`steps p50` is the control and must be roughly flat across each block's CONTROLLED\n\
         rows.** If it moves, `eta` did not hold the step size and this is a step-size comparison\n\
         wearing a re-registration label. The `confounded` rows show what the uncontrolled\n\
         experiment would have claimed.\n\n\
         `hot lift` is against each integrator's OWN baseline hot set: base rate 0.25 by\n\
         construction, so 1.0 is chance and 4.0 is perfect agreement. `chord p50` is the median\n\
         |delta log10 drift|.\n\n\
         **THE RESULT IS THE HG CONTROLLED ROWS.** AZ moved 0.444 decades under this change. If\n\
         Heggie's chord is orders smaller, removing re-registration removed the sensitivity. If\n\
         it is comparable, the mechanism was misidentified and the conclusion on record is wrong.\n\n\
         **Read `refresh_h` before believing a Heggie null.** It is the arm that reintroduces a\n\
         boundary-dependent quantity, and if it also reads zero then this harness cannot see\n\
         boundary sensitivity in Heggie at all and the null is the instrument, not the physics."
    );
}
