//! **Why is Heggie 15,000x worse than AZ on `far`, and is it the control term?**
//!
//! The gallery shows Heggie sitting at a floor of `1.5e-7` on `mid-field` and `5.6e-8` on `far`
//! where AZ reaches `2e-11` and `4e-12`. A floor, not a scaling — so it is not resolution.
//!
//! The suspect is Heggie's Eq. (24) control term `3 Gamma* Q_i / S`, which Eq. (20)/(21) does not
//! have. `Gamma*` is zero on the solution path, so the term is *formally* zero and is retained
//! only for its stabilising effect near triple collision (Heggie §3: it is the difference between
//! having `R_i^{-1}`-growing modes and having none). But numerically `Gamma*` is round-off times
//! the magnitude of its largest term — and on a wide configuration those terms go like `R_j R_k`
//! and `h R1 R2 R3`, so on `far` they are enormous and the "zero" is not small.
//!
//! **The discriminator is direct**: run the same pixels under `Product` (no control term) and
//! under `SumPow32 { keep_gamma_term }` both ways. If the floor is the control term, `Product`
//! and `keep_gamma_term: false` clear it and `true` does not.
//!
//! `eta` is swept as the control. **A floor and a slope are what separate a round-off mechanism
//! from a truncation one**, and a single `eta` cannot tell them apart — this project's own
//! headline diagnostic, applied to a candidate rather than to a bug.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts, HgTime};
use prin_rs::physics::Ic;

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(48);
    let cfg = EnsembleCfg::production();

    println!("Burrau regions, {n}^2, nominal copy only, t_max = 13, no termination.\n");
    println!(
        "  {:>14} {:>22} {:>9} {:>11} {:>11}",
        "region", "arm", "eta", "drift p50", "reg p50"
    );

    // **Termination is the variable.** The gallery runs at production `r_coll_frac` with
    // `stop_on_event`, so a colliding trajectory is parked AT a close approach and its Cartesian
    // energy is a cancellation of two enormous terms — the ill-conditioning already measured on
    // the two-body test, where `drift` ran 280x `drift_reg`. This example ran with `r_coll = 0`
    // and disagreed with the gallery by orders, which is what pointed here.
    for r_coll in [0.0f64, 0.001] {
    println!("=== r_coll_frac = {r_coll}, stop_on_event = {} ===\n", r_coll > 0.0);
    for region in ["far", "mid-field", "deep interior", "near-field"] {
        let Some(sl) = grid::region(region, n, n, 0.05) else { continue };
        let ics: Vec<Ic<f64>> = (0..sl.npix()).map(|k| sl.nominal_ic::<f64>(k)).collect();

        for e in 0..3 {
            let eta = cfg.eta / 4f64.powi(e);
            let mut rows: Vec<(String, Vec<f64>, Vec<f64>)> = Vec::new();

            // AZ, the reference.
            let r: Vec<(f64, f64)> = ics
                .par_iter()
                .map(|ic| {
                    let o = az::integrate_az_opts(
                        ic.s, &ic.m, 13.0, 32, eta, 4_000_000,
                        &AzOpts {
                            stop_on_event: r_coll > 0.0,
                            r_coll_frac: r_coll,
                            step_limit: az::StepLimit::Predictive,
                            step_limit_f: 0.02,
                            ..Default::default()
                        },
                    );
                    // AZ has no regularised-drift analogue exposed; the column repeats its
                    // Cartesian drift so the row cannot be misread as a second measurement.
                    (o.drift, f64::NAN)
                })
                .collect();
            rows.push((
                "AZ".into(),
                r.iter().map(|x| x.0).filter(|x| x.is_finite()).collect(),
                r.iter().map(|x| x.1).filter(|x| x.is_finite()).collect(),
            ));

            for (label, time) in [
                ("HG Eq.(22)-(24)", HgTime::SumPow32 { keep_gamma_term: true }),
                ("HG Eq.(22),(23),(25)", HgTime::SumPow32 { keep_gamma_term: false }),
                ("HG Eq.(20)/(21)", HgTime::Product),
            ] {
                let r: Vec<(f64, f64)> = ics
                    .par_iter()
                    .map(|ic| {
                        let o = integrate_hg(
                            ic.s, &ic.m, 13.0, 32, eta, 4_000_000,
                            &HgOpts {
                                time,
                                r_coll_frac: r_coll,
                                stop_on_event: r_coll > 0.0,
                                ..Default::default()
                            },
                        );
                        // **Both drift measures.** `drift` is the returned Cartesian state's and
                        // is what the gallery plots; `drift_reg` is the integration's own. Their
                        // ratio is the readout's conditioning, not the integrator's error.
                        (o.drift, o.drift_reg)
                    })
                    .collect();
                rows.push((
                    label.into(),
                    r.iter().map(|x| x.0).filter(|x| x.is_finite()).collect(),
                    r.iter().map(|x| x.1).filter(|x| x.is_finite()).collect(),
                ));
            }

            for (label, mut d, mut st) in rows {
                println!(
                    "  {:>14} {:>22} {:>9.2e} {:>11.3e} {:>11.3e}",
                    region,
                    label,
                    eta,
                    q(&mut d, 0.5),
                    q(&mut st, 0.5)
                );
            }
            println!();
        }
    }
    }
    println!(
        "HOW TO READ IT. **A floor is round-off, a slope is truncation.** If the Eq.(22)-(24) row\n\
         is flat across the `eta` ladder while the other two fall, the control term is injecting\n\
         `Gamma*`'s round-off and the fix is to drop it -- at the cost Heggie §3 names, which is\n\
         the R_i inverse modes near triple collision. If ALL Heggie rows are flat, the control term\n\
         is not the mechanism and the floor is somewhere else."
    );
}
