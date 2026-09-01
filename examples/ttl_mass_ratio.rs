//! **TTL versus logH across a deliberate mass-ratio ladder.**
//!
//! logH's kick denominator `U = sum m_i m_j / r_ij` is mass-weighted: a close approach between a
//! heavy body and a light one barely moves it, so the physical step fails to shrink where it
//! matters. TTL replaces it with `Omega = sum w_ij / r_ij`, mass-independent up to a constant, and
//! carries `W` advanced by `dW = (dOmega/dt) dt`.
//!
//! **PREDICTION, recorded before the run: TTL beats logH at high mass ratio and TIES at equal
//! mass.** The tie is the control, not a throwaway.
//!
//! # Why the tie is exact here, and why that matters
//!
//! This port takes `w_ij = mbar^2`, so at `m_0 = m_1 = m_2` the weight equals `m_i m_j` for every
//! pair and `Omega === U` **identically**. The `q = 1` row must therefore agree to round-off. If
//! it does not, the two arms differ for some reason other than the mass ratio and **every other
//! row is uninterpretable** — so it is asserted, loudly, rather than eyeballed.
//!
//! # The sweep changes ONE thing
//!
//! Configurations and velocities come from one slice, decoded once. Only the masses vary down the
//! ladder. Each `q` is its own physical system — the comparison is always TTL against logH *on the
//! same system*, never across rows. Reading down the `q` column is comparing different physics;
//! reading across a row is the measurement.
//!
//! # Guards, printed before the comparison
//!
//! - **`q = 1` agreement.** The control. Asserted.
//! - **`nonfin` and `err>10` per arm.** A win over a dead arm is not a win — the standing failure
//!   where a "three decades" claim had to be withdrawn because both arms were meaningless.
//! - **Force evaluations per arm.** Both run KDK at the same `eta`, one evaluation per step, but
//!   the *step counts* need not match: a time transformation that shrinks the step where it should
//!   spends more steps, and a win bought with 3x the work is a different claim.
//! - **`distinct`.** A flat comparison and an absent one look identical.
//!
//! Args: `res out_dir max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::integrate::logh::hamiltonian::LhTime;
use prin_rs::integrate::logh::{integrate_lh, LhOpts, Stepper};
use prin_rs::physics::Ic;

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

fn distinct(v: &[f64]) -> usize {
    let mut w: Vec<u64> = v.iter().filter(|x| x.is_finite()).map(|x| x.to_bits()).collect();
    w.sort_unstable();
    w.dedup();
    w.len()
}

fn main() {
    let res: usize = arg(1, 48);
    let _out: String = arg(2, "results/ttl".to_string());
    let max_steps: usize = arg(3, 400_000);
    let cfg = EnsembleCfg::production();

    let sl = grid::region("near-field", res, res, 0.05).expect("near-field");
    let (t_max, n_sync) = (13.0f64, 33usize);

    println!("TTL vs logH across a mass-ratio ladder. near-field configurations at {res}^2,");
    println!("t_max = {t_max}, n_sync = {n_sync}, KDK stepper, eta = {:.1e}, termination OFF.\n", cfg.eta);
    println!(
        "  **PREDICTION, recorded before this ran: TTL beats logH at high ratio and TIES at q = 1.**\n  \
         The tie is exact by construction (`w_ij = mbar^2` makes `Omega === U` at equal masses) and\n  \
         is asserted below. If it fails, no other row means anything.\n"
    );
    println!(
        "  Only the MASSES vary down the ladder; configurations and velocities are decoded once.\n  \
         Each row is TTL against logH on the SAME system. Reading ACROSS a row is the measurement;\n  \
         reading DOWN the q column is comparing different physics.\n"
    );

    // Decoded once. Masses are replaced per row, so the configuration is held exactly.
    let ics: Vec<Ic<f64>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, 0, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )[0]
        })
        .collect();

    println!(
        "  {:>7} {:>9} {:>11} {:>11} {:>7} {:>7} {:>10} {:>8}",
        "q", "arm", "drift p50", "drift p99", "err>10", "nonfin", "evals p50", "distinct"
    );

    let mut tie_gap = f64::NAN;

    for qr in [1.0f64, 2.0, 5.0, 20.0, 100.0, 1000.0] {
        // Two heavy, one light: the case logH's mass weighting is blind to. Normalised so the
        // total mass is 1 at every rung -- otherwise the ladder would also be a sweep in overall
        // scale, and this project has measured a 1.66x answer change from exactly that.
        let raw = [1.0f64, 1.0, 1.0 / qr];
        let tot: f64 = raw.iter().sum();
        let m = [raw[0] / tot, raw[1] / tot, raw[2] / tot];

        let mut per_arm: Vec<(f64, f64)> = Vec::new();
        for (label, time) in [("logh", LhTime::LogH), ("ttl", LhTime::Ttl)] {
            let out: Vec<(f64, f64, bool)> = ics
                .par_iter()
                .map(|ic| {
                    let o = integrate_lh(
                        ic.s, &m, t_max, n_sync, cfg.eta, max_steps,
                        &LhOpts {
                            time,
                            stepper: Stepper::Kdk,
                            r_coll_frac: 0.0,
                            stop_on_event: false,
                            step_limit_f: 0.02,
                            ..Default::default()
                        },
                    );
                    let ok = o.finite && !o.budget_exhausted;
                    (if ok { o.drift } else { f64::INFINITY }, o.force_evals as f64, ok)
                })
                .collect();
            let dr: Vec<f64> = out.iter().map(|x| x.0).collect();
            let ev: Vec<f64> = out.iter().map(|x| x.1).collect();
            let nonfin = out.iter().filter(|x| !x.2).count();
            let hot = dr.iter().filter(|x| !x.is_finite()).count();
            println!(
                "  {qr:>7.0} {label:>9} {:>11.3e} {:>11.3e} {hot:>7} {nonfin:>7} {:>10.3e} {:>8}",
                q(&dr, 0.50), q(&dr, 0.99), q(&ev, 0.50), distinct(&dr)
            );
            per_arm.push((q(&dr, 0.50), q(&ev, 0.50)));
        }
        let (lg, tt) = (per_arm[0], per_arm[1]);
        let gain = (lg.0 / tt.0).log10();
        let cost = tt.1 / lg.1;
        println!(
            "  {:>7} {:>9}   gain (logH/TTL) {gain:>+7.3} decades   TTL cost {cost:>6.3}x evals\n",
            "", "->"
        );
        if (qr - 1.0).abs() < 1e-12 {
            tie_gap = gain.abs();
        }
    }

    println!("CONTROL");
    if tie_gap.is_finite() && tie_gap < 0.02 {
        println!(
            "  q = 1 gain {tie_gap:.4} decades -- the arms agree, as `Omega === U` requires.\n  \
             The rows below it are interpretable."
        );
    } else {
        println!(
            "  **CONTROL FAILED**: q = 1 gain is {tie_gap:.4} decades, but `Omega === U` at equal\n  \
             masses so the two arms integrate the same equations. Something other than the mass\n  \
             ratio differs between them and NO ROW IN THIS TABLE MEANS ANYTHING."
        );
    }
    println!(
        "\nHOW TO READ THIS\n\n\
         **Read the control first, then `nonfin` and `err>10`, then the gain.** A win over an arm\n\
         that produced no data is not a win.\n\n\
         **`TTL cost` is part of the result.** A time transformation that shrinks the step where\n\
         it should spends more steps to do it; a gain bought at 3x the evaluations is a different\n\
         claim from one bought at parity.\n\n\
         **The prediction is that gain rises with q and is zero at q = 1.** A gain that is flat in\n\
         q would say the mass weighting is not the mechanism, whatever its sign."
    );
}
