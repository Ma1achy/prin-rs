//! **Does the secant landing improve the TRAJECTORY on a field?** The previous attempt asked with
//! the wrong instrument.
//!
//! `examples/landing_field.rs` measured `energy_drift_max` and found the correction worth 0.03%,
//! with both its guards passing — the residual fell nine orders and the correction fired 7.5
//! million times. That is not a null, it is a **blind diagnostic**, and `CLAUDE.md:1292` says so
//! about the same defect one revision earlier: *"energy drift is blind to this one... the
//! overshoot displaces the state in TIME and the energy is nearly stationary along the flow."*
//!
//! A landing residual is a displacement **in time**. Energy barely notices; position does. The
//! figure-eight, where the correction is worth 1960x, measures *closure* — a position and
//! velocity error — which is exactly the quantity a time displacement moves.
//!
//! So this measures the shape chord against a converged reference.
//!
//! # Two things that would make it meaningless, handled rather than hoped
//!
//! **Chaos saturates a chord.** Over a long horizon in a chaotic region every arm diverges from
//! the reference and the chord goes to 2.000 — antipodal — which this project already records as
//! *not evidence of anything*. So the horizon is swept **short**, `t = 1, 2, 4, 8`, and
//! `chord p50` is printed at each: where it saturates the row is dead and says so on its face.
//!
//! **A reference computed under one arm's settings favours that arm.** The reference here is
//! **GBS with the landing on at `eta/16`** — the most accurate configuration available, chosen
//! so that neither `land` arm is its own yardstick. Its own landing residual is negligible, so
//! it is not "the land-on answer" but the converged one.
//!
//! Args: `res`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::logh::{integrate_lh, LhOpts, Stepper};
use prin_rs::physics::{shape, Cart, Ic};

fn q(v: &[f64], p: f64) -> f64 {
    let mut w: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if w.is_empty() {
        return f64::NAN;
    }
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    w[(((w.len() - 1) as f64) * p).round() as usize]
}

fn chord(a: &Cart<f64>, b: &Cart<f64>, m: &[f64; 3]) -> f64 {
    let (x, y) = (shape::shape_vec(&a.r, m), shape::shape_vec(&b.r, m));
    ((0..3).map(|i| (x[i] - y[i]).powi(2)).sum::<f64>()).sqrt()
}

fn run(ic: &Ic<f64>, t_max: f64, n_sync: usize, eta: f64, stepper: Stepper, land: bool) -> Cart<f64> {
    let o = integrate_lh(
        ic.s, &ic.m, t_max, n_sync, eta, 20_000_000,
        &LhOpts {
            stepper,
            land_iterate: land,
            r_coll_frac: 0.0,
            stop_on_event: false,
            step_limit_f: 0.0,
            gbs_tol: 1e-13,
            gbs_k_max: 8,
            ..Default::default()
        },
    );
    if o.finite { o.state } else { Cart::new([Default::default(); 3], [Default::default(); 3]) }
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(48);
    let cfg = EnsembleCfg::production();
    let (chart, cx, cy, half) = Chart::config_stability();
    let mut cases: Vec<(&str, Chart, f64, f64, f64, usize)> =
        vec![("config_stability", chart, cx, cy, half, 0)];
    for n in ["near-field", "deep interior"] {
        if let Some(s) = grid::region(n, 4, 4, 0.05) {
            let nm = if n == "deep interior" { "deep_interior" } else { "near-field" };
            cases.push((nm, s.chart, s.cx, s.cy, s.half, s.body));
        }
    }

    println!("{res}^2, shape chord against a converged reference (GBS, land ON, eta/16)\n");
    println!(
        "  **`chord` saturates at 2.000 (antipodal) once chaos dominates, and a saturated row is\n  \
         dead.** The horizon is swept short for that reason: where `off` and `on` are both near\n  \
         2.0 the comparison says nothing, whatever the ratio between them.\n"
    );
    println!(
        "  {:>17} {:>6} {:>5} {:>11} {:>11} {:>9} {:>11}",
        "case", "step", "t", "chord off", "chord on", "gain", "ref chord"
    );

    for (name, ch, cx, cy, half, body) in cases {
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(ch);
        let ics: Vec<Ic<f64>> = (0..sl.npix())
            .into_par_iter()
            .map(|k| {
                jitter::copies_with_path::<f64>(
                    &sl, k, 0, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
                )[0]
            })
            .collect();

        for t_max in [1.0f64, 2.0, 4.0, 8.0] {
            let n_sync = (t_max / 0.4).round().max(4.0) as usize;
            for stepper in [Stepper::Rk4, Stepper::Kdk] {
                let rows: Vec<(f64, f64, f64)> = ics
                    .par_iter()
                    .map(|ic| {
                        let re = run(ic, t_max, n_sync, cfg.eta / 16.0, Stepper::Gbs, true);
                        let off = run(ic, t_max, n_sync, cfg.eta, stepper, false);
                        let on = run(ic, t_max, n_sync, cfg.eta, stepper, true);
                        // A second reference at half the step: if the two references disagree as
                        // much as the arms do, the reference is not converged and no row means
                        // anything.
                        let re2 = run(ic, t_max, n_sync, cfg.eta / 32.0, Stepper::Gbs, true);
                        (
                            chord(&off, &re, &ic.m),
                            chord(&on, &re, &ic.m),
                            chord(&re2, &re, &ic.m),
                        )
                    })
                    .collect();
                let (o, n, r): (Vec<f64>, Vec<f64>, Vec<f64>) = (
                    rows.iter().map(|x| x.0).collect(),
                    rows.iter().map(|x| x.1).collect(),
                    rows.iter().map(|x| x.2).collect(),
                );
                let (mo, mn) = (q(&o, 0.5), q(&n, 0.5));
                println!(
                    "  {name:>17} {:>6} {t_max:>5} {mo:>11.4e} {mn:>11.4e} {:>9.3} {:>11.4e}",
                    format!("{stepper:?}"),
                    if mn > 0.0 { mo / mn } else { f64::NAN },
                    q(&r, 0.5)
                );
            }
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         **`ref chord` is the guard.** It is the distance between the reference and a second\n\
         reference at half its step. If it is not far below `chord off`, the reference is not\n\
         converged and the row is measuring reference error rather than the arms.\n\n\
         **`gain` is `off / on`.** Above 1 the correction helps the trajectory. The figure-eight\n\
         says 1960x for RK4; if a field says ~1, the correction is a periodic-orbit result and\n\
         does not transfer, and the case for porting it to AZ and Heggie fails on the evidence\n\
         rather than on a blind instrument.\n\n\
         **KDK remains the control**: second order, never limited by an O(h^2) landing, so its\n\
         gain must stay near 1 whatever the others do."
    );
}
