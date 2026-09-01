//! **STAGE 0 — is AZ's `far` win a CONDITIONING effect? Ask f32.**
//!
//! AZ beats Heggie on `far` on all 65536 pixels by a flat 0.7-0.9 decades. The surviving
//! explanation is conditioning: `far` spans body positions to ~13 units where the latent charts
//! sit at `R = 1` by algebraic identity, and Heggie's `Gamma*` is **degree six** in the
//! coordinates where AZ's `Gamma` is linear in `A` and `B`. That explanation is labelled a guess
//! in `CLAUDE.md`, and it is precision-sensitive by construction — which makes it cheap to test.
//!
//! # The decision rule, and why I think the brief has it backwards
//!
//! `Gamma*` forms intermediates of order `13^6 ~ 4.8e6` and cancels them to `O(1)`, costing ~6.7
//! decimal digits. f64 has ~15.95 and survives with ~9.2; **f32 has ~7.22 and has nothing left**.
//! So if conditioning is the mechanism, AZ's advantage should **WIDEN** at f32. A collapse would
//! mean the penalty does not care how many digits are available, which is nearly the definition of
//! *not* a precision effect.
//!
//! # The guard that decides whether this run says anything at all
//!
//! `far`'s f64 drifts are `2.8e-13` and `2.2e-12`. f32's documented median drift on this project
//! is `9.3e-6` — **seven orders above both**. If f32 round-off dominates for both arms they land
//! on a common floor and the ratio is a difference between two meaningless numbers. That is the
//! failure already on record, where a "three decades on `far`" claim had to be withdrawn because
//! both arms carried `err>10` on every pixel.
//!
//! So this prints, **before** the comparison: distinct drift values per arm, the overlap of the
//! two distributions, and `err>10` per arm. *A difference can be small because both sides are
//! right or because both are dead.*
//!
//! Args: `res out_dir max_steps`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::integrate::Integrator;

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

/// How many distinct values a signal takes. **Read this before any curve**: a flat comparison and
/// an absent one look identical, and only the count separates them.
fn distinct(v: &[f64]) -> usize {
    let mut w: Vec<u64> = v.iter().filter(|x| x.is_finite()).map(|x| x.to_bits()).collect();
    w.sort_unstable();
    w.dedup();
    w.len()
}

fn main() {
    let res: usize = arg(1, 128);
    let out_dir: String = arg(2, "results/overnight".to_string());
    let max_steps: usize = arg(3, 400_000);

    let sl = grid::region("far", res, res, 0.05).expect("far region");
    let t_max = 13.0f64;
    let n_sync = (t_max / 0.4).round().max(4.0) as usize;

    println!("STAGE 0: `far` at {res}^2, AZ vs Heggie, f32 vs f64. Nothing else changed.\n");
    println!(
        "  Diagnostic pass: termination OFF, r_coll = 0, refine_flagged OFF. A run stopped by an\n  \
         event is parked at a close approach where the Cartesian energy is a cancellation of two\n  \
         enormous terms -- that flag alone produced five wrong conclusions in a row once.\n"
    );
    println!(
        "  PREDICTION (recorded in PREDICTIONS.md before this ran): if conditioning is the\n  \
         mechanism, AZ's advantage WIDENS at f32. The brief expects a collapse; I think that\n  \
         implication runs the wrong way and this run adjudicates it.\n"
    );

    let cfg_for = |i: Integrator| {
        EnsembleCfg::production().with_overrides(&[
            Override::TMax(t_max),
            Override::NSync(n_sync),
            Override::RCollFrac(0.0),
            Override::StopOnEvent(false),
            Override::RefineFlagged(false),
            Override::Integrator(i),
            Override::MaxSteps(max_steps),
        ])
    };

    let mut rows: Vec<(String, String, Vec<f64>, Vec<f64>, usize, usize)> = Vec::new();

    for (iname, integ) in [("az", Integrator::Az), ("heggie", Integrator::Heggie)] {
        let cfg = cfg_for(integ);
        for prec in ["f64", "f32"] {
            let t0 = std::time::Instant::now();
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| {
                    if prec == "f64" {
                        pixel::evaluate::<f64>(&sl, k, &cfg)
                    } else {
                        pixel::evaluate::<f32>(&sl, k, &cfg)
                    }
                })
                .collect();
            let secs = t0.elapsed().as_secs_f64();
            let dr: Vec<f64> = px.iter().map(|p| p.energy_drift_max).collect();
            let er: Vec<f64> = px.iter().map(|p| p.error_ratio).collect();
            let hot = px.iter().filter(|p| p.error_ratio > 10.0).count();
            let ev: Vec<f64> = px.iter().map(|p| p.total_force_evals as f64).collect();
            println!(
                "  {iname:>7} {prec:>4}  drift p50 {:>10.3e}  p99 {:>10.3e}  err>10 {hot:>6}  \
                 distinct {:>6}  evals p50 {:>9.3e}  {secs:>7.1}s",
                q(&dr, 0.50), q(&dr, 0.99), distinct(&dr), q(&ev, 0.50)
            );
            rows.push((iname.into(), prec.into(), dr, er, hot, distinct(&ev)));
        }
    }

    // --- the saturation guard, printed BEFORE the comparison ---------------------------------
    println!("\n  SATURATION GUARD -- does this comparison have a subject?");
    for prec in ["f64", "f32"] {
        let a = rows.iter().find(|r| r.0 == "az" && r.1 == prec).unwrap();
        let h = rows.iter().find(|r| r.0 == "heggie" && r.1 == prec).unwrap();
        let (da, dh) = (distinct(&a.2), distinct(&h.2));
        // Overlap of the two drift distributions: the fraction of AZ pixels whose drift lies
        // inside Heggie's p1-p99 range. Near 1 means the arms are indistinguishable.
        let (lo, hi) = (q(&h.2, 0.01), q(&h.2, 0.99));
        let ov = a.2.iter().filter(|x| x.is_finite() && **x >= lo && **x <= hi).count() as f64
            / a.2.len() as f64;
        let verdict = if da < 8 || dh < 8 {
            "DEAD -- a drift column with <8 distinct values is a floor, not a measurement"
        } else if ov > 0.9 {
            "SATURATED -- the two arms overlap almost completely; any ratio is noise"
        } else {
            "LIVE -- the arms are separable at this precision"
        };
        println!(
            "    {prec}: distinct az {da:>6} hg {dh:>6}   az-inside-hg p1..p99 {ov:>6.4}   {verdict}"
        );
    }

    // --- the comparison ------------------------------------------------------------------------
    println!("\n  AZ ADVANTAGE  gain = log10(drift_heggie / drift_az), positive = AZ better");
    for prec in ["f64", "f32"] {
        let a = rows.iter().find(|r| r.0 == "az" && r.1 == prec).unwrap();
        let h = rows.iter().find(|r| r.0 == "heggie" && r.1 == prec).unwrap();
        let g: Vec<f64> = a
            .2
            .iter()
            .zip(h.2.iter())
            .map(|(x, y)| (y / x).log10())
            .filter(|x| x.is_finite())
            .collect();
        let better = g.iter().filter(|x| **x > 0.0).count() as f64 / g.len().max(1) as f64;
        println!(
            "    {prec}: gain p10 {:>7.3}  p50 {:>7.3}  p90 {:>7.3}  frac az better {better:>6.4}  n {}",
            q(&g, 0.10), q(&g, 0.50), q(&g, 0.90), g.len()
        );
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **Read the saturation guard first.** If either precision reads DEAD or SATURATED, the\n\
         gain line beneath it for that precision is a difference between two numbers that are not\n\
         measuring anything, and it does not adjudicate the conditioning question.\n\n\
         **The decision, if the f32 row is LIVE:** an AZ advantage that GROWS at f32 supports\n\
         conditioning and justifies Stage 2 (chain coordinates). One that SHRINKS says the\n\
         mechanism is not precision, and chain is unlikely to fix `far` either.\n\n\
         **`err>10` is the project's own flag for *this pixel is not data*.** An f32 arm carrying\n\
         it on a large share of pixels is not a weaker measurement, it is a different one."
    );
    let _ = out_dir;
}
