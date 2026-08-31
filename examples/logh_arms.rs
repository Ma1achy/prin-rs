//! **The falsification test for the re-registration mechanism.**
//!
//! The standing claim, from `38bd303` and `INVESTIGATION.md` §2.3: doubling the sync cadence at
//! fixed step size moves AZ's drift field by 0.44-0.52 decades and Heggie's by 0.048, so it is
//! not *which* chart is chosen but *how often the state is passed through one*. Heggie then wins
//! 31 of 32 gallery cases. That is consistent with the mechanism and does not establish it:
//! Heggie changes the chart **and** removes the re-registration in one move.
//!
//! logH separates them. Algorithmic regularisation has **no coordinate transformation at all** —
//! a time transformation and a good integrator — which is a strictly stronger form of the
//! property the Heggie win is attributed to.
//!
//! ```text
//!   If the mechanism is real, logH matches or beats Heggie.
//!   If logH loses to Heggie, the mechanism is wrong and the win comes from something else.
//! ```
//!
//! Both outcomes are informative. **Do not tune toward the first.** The predictions are recorded
//! in the plan and in `NOTES.md` before any number here came back.
//!
//! # Why six arms and not two
//!
//! **There is no common stepper, and pretending otherwise would repeat the confound that already
//! cost this project a factor of 6000.** `Az + leapfrog` and `Heggie + leapfrog` do not exist:
//! their `Gamma` couples position and momentum, so neither Hamiltonian is separable. And logH's
//! regularisation **is** the leapfrog — Mikkola & Merritt, *"in these new methods the
//! regularization is achieved by using the leapfrog"* — already confirmed in
//! `tests/logh_march.rs`, where KDK traverses the radial collision and RK4 does not at any step
//! size tried.
//!
//! So logH is run under both steppers and the arms are matched on **force evaluations, not
//! steps**: RK4 spends four per step and KDK one, so the leapfrog arms run at `eta/4`.
//! Accuracy per unit compute is the operationally meaningful axis anyway, and it does not
//! require crippling any method to obtain.
//!
//! **`eta/4` is a NOMINAL match and the measured column is the truth.** Steps per interval go as
//! `1/eta` only while nothing else is resizing the step, and the predictive limit is not an `eta`
//! multiple — its bound `ds <= f d_min (K+B)/|v_rel|` is absolute, so it binds *less* at a
//! smaller `eta` and the step count comes out below `4x`. Measured on `far`: `logh_rk4` spends
//! `2.345e5` evaluations and `logh_lf` at `eta/4` spends `1.477e5`, a shortfall of 1.6x rather
//! than parity. So `evals p50` is printed on every row and **a row is read at the evaluations it
//! actually spent**, not at the ones its `eta` was supposed to buy. Quoting the nominal match
//! and ignoring the column would be the `steps p50` control of `az_machinery.rs` printed and
//! then not read.
//!
//! And the **stepper-only control** is what makes the rest readable: `plain_rk4` and `plain_lf`
//! are the same code path with the time transformation switched off. If they differ substantially
//! at equal evaluations, the stepper is contributing to the comparison and the amount is now
//! measured. If they agree, the regularisation comparison is clean. Either way it is quantified
//! rather than assumed.
//!
//! # The predictive step limit is swept, not held
//!
//! It is production for AZ and Heggie at `f = 0.02` and **it defeats logH at an exact collision**
//! — the bound is `ds <= f d_min U/|v_rel| ~ f sqrt(d)` while the unbounded step wants `1/d`, so
//! it is *do not shrink the fictitious step at close approach* at a third site. On Burrau it is a
//! large improvement instead. A knob held fixed for fairness is a knob whose effect is
//! unattributed, so both logH arms appear twice, with the limit and without.
//!
//! # Two passes per arm, for the reason `integrator_gallery` has two
//!
//! A trajectory stopped by `stop_on_event` is parked at a close approach, where the Cartesian
//! energy is a cancellation of two enormous terms. That flag alone produced five wrong
//! conclusions in a row once. The science pass has termination on and feeds `_uniform` and
//! `_outcome`; the diagnostic pass runs at `r_coll = 0`, `stop_on_event = false` and feeds
//! `_drift` and the gain maps.
//!
//! `refine_flagged` is **off** in both. It re-integrates from `t = 0` at finer `eta` and has no
//! live-playhead analogue, and it removes the very population a comparison of integrators is
//! about — the first `integrator_gallery` render had it on and showed no wedges under *either*
//! integrator, which said nothing.
//!
//! Args: `res root only max_steps`.

use std::collections::HashMap;

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::StepLimit;
use prin_rs::integrate::Integrator;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png, provenance_sidecar};

/// The drift ramp, FIXED across every case and every arm — the same constants
/// `integrator_gallery`, `az_machinery` and `heggie_machinery` use, so panels from all four are
/// directly comparable. A per-arm window on a question about which arm is cleaner would
/// manufacture or hide exactly the thing being measured.
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

/// The gain map's half-range, in decades. Symmetric and FIXED, never auto-ranged: an auto-ranged
/// diverging map centres itself on whatever it is given and paints a null and a rout alike.
const GAIN_DECADES: f64 = 4.0;

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

fn qs(v: &[f64], p: f64) -> f64 {
    q(&mut v.to_vec(), p)
}

/// `log10(reference/arm)` on a fixed symmetric ramp: blue where `arm` is lower (better), red where
/// the reference is, near-white where they agree, magenta where either is undetermined.
fn gain_rgb(reference: f64, arm: f64) -> [u8; 3] {
    if !reference.is_finite() || !arm.is_finite() || reference <= 0.0 || arm <= 0.0 {
        return [255, 0, 255];
    }
    let g = ((reference / arm).log10() / GAIN_DECADES).clamp(-1.0, 1.0);
    let t = g.abs();
    let base = 0.82;
    let c = |x: f64| (255.0 * x.clamp(0.0, 1.0)) as u8;
    if g >= 0.0 {
        [c(base - 0.62 * t), c(base - 0.42 * t), c(base + 0.18 * t)]
    } else {
        [c(base + 0.18 * t), c(base - 0.55 * t), c(base - 0.55 * t)]
    }
}

struct Case {
    name: String,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    t_max: f64,
    r_coll: f64,
}

/// **Six cases, and they are the six the 1024^2 gallery finished.**
///
/// Not a sample of convenience: `far` is the one case AZ wins outright (all 65536 pixels, a flat
/// 0.7-0.9 decades) and is prediction 2's subject; `preset_shape_h1` is the case Heggie helps
/// most (+4.42 decades, AZ 2222 not-data pixels against zero) and is also the only tree in the
/// corpus the refinement criterion actually controls; `deep_interior` and `config_stability` are
/// where the wedges live. Their AZ and Heggie panels already exist at both resolutions, so every
/// row here has a reference.
fn cases() -> Vec<Case> {
    let mut v = Vec::new();
    let (chart, cx, cy, half) = Chart::config_stability();
    v.push(Case {
        name: "config_stability".into(),
        chart, cx, cy, half, body: 0,
        t_max: 50.0,
        r_coll: 0.005,
    });
    for name in ["near-field", "far", "deep interior"] {
        if let Some(sl) = grid::region(name, 4, 4, 0.05) {
            v.push(Case {
                name: name.replace(' ', "_"),
                chart: sl.chart, cx: sl.cx, cy: sl.cy, half: sl.half, body: sl.body,
                t_max: 13.0,
                r_coll: 0.001,
            });
        }
    }
    for (name, chart, cx, cy, half) in grid::gallery_cases() {
        if name == "preset_shape" || name == "preset_shape_h1" {
            v.push(Case {
                name: name.into(),
                chart, cx, cy, half, body: 0,
                t_max: 13.0,
                r_coll: 0.001,
            });
        }
    }
    v
}

/// One arm: what to run, and how to make it comparable.
struct Arm {
    /// Panel stem and table label.
    label: &'static str,
    integ: Integrator,
    /// Multiplier on `eta`. **`0.25` for the leapfrog arms**, so force evaluations match: RK4
    /// spends four per step and KDK one, so a quarter of the step size is the same work.
    eta_scale: f64,
    /// The predictive step limit. Swept rather than held, because it is production for AZ and
    /// Heggie and is fatal to logH at an exact collision.
    limit: bool,
}

fn arms() -> Vec<Arm> {
    vec![
        Arm { label: "az", integ: Integrator::Az, eta_scale: 1.0, limit: true },
        Arm { label: "heggie", integ: Integrator::Heggie, eta_scale: 1.0, limit: true },
        Arm { label: "logh_rk4", integ: Integrator::LogHRk4, eta_scale: 1.0, limit: true },
        Arm { label: "logh_rk4_nolim", integ: Integrator::LogHRk4, eta_scale: 1.0, limit: false },
        Arm { label: "logh_lf", integ: Integrator::LogHLeapfrog, eta_scale: 0.25, limit: true },
        Arm { label: "logh_lf_nolim", integ: Integrator::LogHLeapfrog, eta_scale: 0.25, limit: false },
        // The stepper-only control, read FIRST. If these two differ at equal evaluations, the
        // stepper is in the comparison and by how much is now measured.
        Arm { label: "plain_rk4", integ: Integrator::PlainRk4, eta_scale: 1.0, limit: false },
        Arm { label: "plain_lf", integ: Integrator::PlainLeapfrog, eta_scale: 0.25, limit: false },
    ]
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    let max_steps: usize = arg(4, 400_000);
    let dir = format!("{root}/logh_arms");
    let _ = std::fs::create_dir_all(&dir);

    println!(
        "{res}^2, six occupants, shared colour windows, max_steps={max_steps}, \
         refine_flagged=false.\n"
    );
    println!(
        "  **The UNMASKED kernel.** The repair pass re-integrates from t = 0, has no live\n  \
         analogue, and removes the population this comparison is about.\n"
    );
    println!(
        "  **Arms are matched on FORCE EVALUATIONS, not steps.** RK4 spends 4 per step and KDK 1,\n  \
         so the leapfrog arms run at eta/4. That is a NOMINAL match: the predictive limit is an\n  \
         absolute bound, not an `eta` multiple, so it binds less at smaller `eta` and the step\n  \
         count lands below 4x. **Read each row at the evaluations it actually spent.**\n"
    );
    println!(
        "  `rho/gamma` is each occupant's own regularised-Hamiltonian residual and is NOT one\n  \
         quantity across the table: AZ's is |Gamma|/largest term, Heggie's the same for Gamma*,\n  \
         logH's |K+B-U|/U -- which is the energy defect normalised by U, not an independent\n  \
         constraint. Compare it down a column, never across one.\n"
    );

    let hdr = format!(
        "  {:>18} {:>15} {:>10} {:>10} {:>10} {:>6} {:>10} {:>10} {:>6} {:>6} {:>6} {:>6} {:>7}",
        "case", "arm", "drift p50", "drift p99", "err p99", "err>10", "steps p50", "evals p50",
        "nonfin", "budget", "over", "escape", "secs"
    );

    for c in cases() {
        if only != "all" && only != c.name {
            continue;
        }
        // The last panel written for the last arm is the resume sentinel, exactly as in
        // `integrator_gallery`. Nothing is resumed across a code change, because nothing here
        // can detect one.
        let done = format!("{dir}/{}_{}_gain_vs_heggie.png", c.name, arms().last().unwrap().label);
        if std::path::Path::new(&done).exists() {
            println!("  {:>18}  (already rendered, skipping)", c.name);
            continue;
        }

        let n_sync = (c.t_max / 0.4).round().max(4.0) as usize;
        let sl = grid::Slice::body_plane(res, res, c.cx, c.cy, c.half, c.body).with_chart(c.chart);
        let m_here = grid::decode_state(&c.chart, c.body, c.cx, c.cy).m;
        let sites = colour::landmarks(&m_here);

        // Windows taken ONCE, from the `az` arm, and shared. Held in Options so a later arm
        // cannot silently fall back to its own range.
        let mut window: Option<(f64, f64)> = None;
        let mut extra_win: HashMap<&str, (f64, f64)> = HashMap::new();
        // The two references every gain map is drawn against.
        let mut az_drift: Option<Vec<f64>> = None;
        let mut hg_drift: Option<Vec<f64>> = None;

        println!("\n{hdr}");

        for a in arms() {
            let t0 = std::time::Instant::now();
            let eta = EnsembleCfg::production().eta * a.eta_scale;
            let ov = |r_coll: f64, stop: bool| {
                vec![
                    Override::TMax(c.t_max),
                    Override::NSync(n_sync),
                    Override::RCollFrac(r_coll),
                    Override::EscapeRule(EscapeRule::Closure(CLOSURE_TAU)),
                    Override::ClosureK(1),
                    Override::Integrator(a.integ),
                    Override::Eta(eta),
                    Override::MaxSteps(max_steps),
                    Override::RefineFlagged(false),
                    Override::StopOnEvent(stop),
                    Override::StepLimit(if a.limit { StepLimit::Predictive } else { StepLimit::None }),
                ]
            };

            let ens = EnsembleCfg::production().with_overrides(&ov(c.r_coll, true));
            let px: Vec<PixelOut> =
                (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();
            let dens = EnsembleCfg::production().with_overrides(&ov(0.0, false));
            let dpx: Vec<PixelOut> =
                (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &dens)).collect();
            let secs = t0.elapsed().as_secs_f64();

            // --- the row ---
            let dr: Vec<f64> =
                dpx.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
            let nonfin = dpx.len() - dr.len();
            let er: Vec<f64> =
                dpx.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
            let hot = dpx.iter().filter(|p| p.error_ratio > 10.0).count();
            let budget = px.iter().chain(dpx.iter()).filter(|p| p.budget_exhausted).count();
            let over: u64 = dpx.iter().map(|p| p.n_overshoot as u64).sum();
            let steps: Vec<f64> = dpx.iter().map(|p| p.total_substeps as f64).collect();
            let evals: Vec<f64> = dpx.iter().map(|p| p.total_force_evals as f64).collect();
            let esc = px.iter().filter(|p| p.state == 0).count() as f64 / px.len() as f64;
            println!(
                "  {:>18} {:>15} {:>10.3e} {:>10.3e} {:>10.3e} {:>6} {:>10.3e} {:>10.3e} {:>6} \
                 {:>6} {:>6} {:>6.4} {:>7.1}",
                c.name, a.label, qs(&dr, 0.50), qs(&dr, 0.99), qs(&er, 0.99), hot,
                qs(&steps, 0.50), qs(&evals, 0.50), nonfin, budget, over, esc, secs
            );

            // --- panels ---
            let (lo, hi) = *window.get_or_insert_with(|| colour::range(&px, Scalar::ShapeSpread));
            let mut sbuf = Vec::with_capacity(px.len() * 3);
            let mut obuf = Vec::with_capacity(px.len() * 3);
            let mut dbuf = Vec::with_capacity(dpx.len() * 3);
            for p in &px {
                sbuf.extend(colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
                obuf.extend(png::outcome_rgb(p));
            }
            for p in &dpx {
                dbuf.extend(colour::drift_rgb(p, DLO, DHI));
            }
            let ew = *extra_win
                .entry("errratio")
                .or_insert_with(|| colour::range(&dpx, Scalar::ErrorRatio));
            let mut ebuf = Vec::with_capacity(dpx.len() * 3);
            for p in &dpx {
                ebuf.extend(colour::rgb(p, Scalar::ErrorRatio, &sites, ew.0, ew.1));
            }

            let stem = format!("{dir}/{}_{}", c.name, a.label);
            let mut panels: Vec<(&str, Vec<u8>)> = vec![
                ("uniform", sbuf),
                ("outcome", obuf),
                ("drift", dbuf),
                ("errratio", ebuf),
            ];

            // The gain maps, against BOTH references. A quantile ladder says the distribution
            // moved and cannot say **where** — the standing rule, and the one I broke by reading
            // a count as though it located anything.
            let this: Vec<f64> = dpx.iter().map(|p| p.energy_drift_max).collect();
            for (name, reference) in
                [("gain_vs_az", &az_drift), ("gain_vs_heggie", &hg_drift)]
            {
                if let Some(r) = reference {
                    let buf: Vec<u8> =
                        (0..this.len()).flat_map(|i| gain_rgb(r[i], this[i])).collect();
                    let g: Vec<f64> = (0..this.len())
                        .filter(|&i| r[i] > 0.0 && this[i] > 0.0 && r[i].is_finite() && this[i].is_finite())
                        .map(|i| (r[i] / this[i]).log10())
                        .collect();
                    let better = g.iter().filter(|x| **x > 0.0).count() as f64 / g.len().max(1) as f64;
                    println!(
                        "  {:>18} {:>15}   {name}: p10 {:+.2}  p50 {:+.2}  p90 {:+.2}  \
                         frac better {better:.4}  n {}",
                        "", a.label, qs(&g, 0.10), qs(&g, 0.50), qs(&g, 0.90), g.len()
                    );
                    panels.push((name, buf));
                }
            }

            for (suffix, buf) in panels {
                let path = format!("{stem}_{suffix}.png");
                let _ = adaptive::save_rect(&path, res, res, &buf);
                let _ = provenance_sidecar(
                    &path,
                    if suffix == "uniform" || suffix == "outcome" { &ens } else { &dens },
                    &format!(
                        "res={res}x{res}\ncase={}\narm={}\nintegrator={}\nfield={suffix}\n\
                         eta={eta:e} (scale {} on production, to match FORCE EVALUATIONS)\n\
                         predictive step limit={}\n\
                         shape window=({lo:e},{hi:e}) taken from the `az` arm and SHARED\n\
                         drift ramp=({DLO:e},{DHI:e}) FIXED across every case and arm\n\
                         gain ramp=+/-{GAIN_DECADES} decades, FIXED and symmetric\n\
                         pass: uniform/outcome from the SCIENCE run (termination on);\n\
                         drift/errratio/gain from the DIAGNOSTIC run (r_coll=0, stop_on_event=false)\n",
                        c.name, a.label, a.integ.name(), a.eta_scale, a.limit
                    ),
                );
            }

            if a.label == "az" {
                az_drift = Some(this);
            } else if a.label == "heggie" {
                hg_drift = Some(this);
            }
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **Read the stepper-only control first.** `plain_rk4` against `plain_lf` at equal\n\
         evaluations is the stepper's own contribution. If they differ substantially, that\n\
         amount is inside every other comparison in the table; if they agree, the\n\
         regularisation comparison is clean.\n\n\
         **Then `logh_lf` against `heggie`.** That is the falsification test. logH has no chart\n\
         and therefore no re-registration at all, which is a stronger form of the property\n\
         Heggie's win is attributed to. If it loses, the attribution is wrong or incomplete and\n\
         the next candidates are the KS square-root's own round-off and the per-boundary energy\n\
         re-freeze -- a real fork that deserves its own run, not a paragraph.\n\n\
         **`far` is prediction 2.** AZ wins it outright today, on all 65536 pixels by a flat\n\
         0.7-0.9 decades, and that win is attributed to AZ never re-registering there. A method\n\
         that never re-registers anywhere should not lose it.\n\n\
         **`logh_rk4` against `logh_lf` is prediction 4.** On shell `K + B == U`, so an\n\
         integrator that evaluates both denominators at the same point sees only a Sundman\n\
         transformation; the leapfrog sees them at different points and that asymmetry is the\n\
         method. If the two arms come out alike, that prediction is wrong and so is Mikkola &\n\
         Merritt's sentence it rests on.\n\n\
         **`err>10` is the project's own flag for *this pixel is not data*.** A count locates\n\
         nothing -- read the gain map for where, and the decile ladder for how much."
    );
}
