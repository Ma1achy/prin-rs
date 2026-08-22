//! Step 6: the question the whole build exists to settle — is Aarseth–Zare usable at f32?
//!
//! Four combinations, {f32, f64} x shared-reference {on, off}, plus the three specific
//! questions: whether the conditioned LC branch removes the f32 `spread_shape` inflation,
//! whether the shared-reference flag matters, and whether the branch cut reaches the outcome
//! encoding at f32 as it does not at f64.
//!
//! **Initial conditions are generated once in f64 and cast down** (`ensemble/jitter.rs`), so
//! nothing here is an IC difference wearing an f32 costume.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Slice};
use prin_rs::integrate::az::RefPolicy;
use prin_rs::render::Precision;

const SIZE: usize = 32;

fn render(p: Precision, cfg: &EnsembleCfg, s: &Slice) -> Vec<PixelOut> {
    (0..s.npix())
        .into_par_iter()
        .map(|i| match p {
            Precision::F32 => evaluate::<f32>(s, i, cfg),
            Precision::F64 => evaluate::<f64>(s, i, cfg),
        })
        .collect()
}

fn q(v: impl Iterator<Item = f64>, f: f64) -> f64 {
    let mut x: Vec<f64> = v.filter(|q| q.is_finite()).collect();
    if x.is_empty() {
        return f64::NAN;
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x[(((x.len() - 1) as f64) * f).round() as usize]
}

fn main() {
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();

    println!("near-field {SIZE}x{SIZE}, t=13, E+1=8, eta=0.01, r_coll=1e-3 R");
    println!("ICs built once in f64 and cast down, so no row differs by its initial condition.");
    println!();
    println!("=== the four combinations ===");
    println!("{:>6}{:>10}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}{:>9}",
             "prec", "shared", "drift med", "drift max", "er med", "er max",
             "ens_sp med", "sp_shape med", "refdis");

    let mut store: Vec<(Precision, RefPolicy, Vec<PixelOut>)> = Vec::new();
    for prec in [Precision::F64, Precision::F32] {
        for policy in [RefPolicy::PerCopy, RefPolicy::Shared] {
            let cfg = EnsembleCfg { ref_policy: policy, ..Default::default() };
            let px = render(prec, &cfg, &s);
            println!("{:>6}{:>10}{:>12.3e}{:>12.3e}{:>12.4e}{:>12.4e}{:>12.4e}{:>12.4e}{:>9}",
                     prec.name(),
                     if matches!(policy, RefPolicy::Shared) { "on" } else { "off" },
                     q(px.iter().map(|p| p.energy_drift_max), 0.5),
                     q(px.iter().map(|p| p.energy_drift_max), 1.0),
                     q(px.iter().map(|p| p.error_ratio), 0.5),
                     q(px.iter().map(|p| p.error_ratio), 1.0),
                     q(px.iter().map(|p| p.ensemble_spread), 0.5),
                     q(px.iter().map(|p| p.spread_shape), 0.5),
                     px.iter().map(|p| p.ref_disagree as u64).sum::<u64>());
            store.push((prec, policy, px));
        }
    }

    println!();
    println!("{:>6}{:>10}{:>14}{:>14}{:>12}{:>14}",
             "prec", "shared", "sp_event nz", "sp_event max", "nonfinite", "outcome flips");
    let f64_unshared: Vec<u8> = store[0].2.iter().map(|p| p.outcome).collect();
    for (prec, policy, px) in &store {
        let flips = px
            .iter()
            .zip(f64_unshared.iter())
            .filter(|(a, b)| a.outcome != **b)
            .count();
        println!("{:>6}{:>10}{:>14}{:>14.4}{:>12}{:>14}",
                 prec.name(),
                 if matches!(policy, RefPolicy::Shared) { "on" } else { "off" },
                 format!("{}/{}", px.iter().filter(|p| p.spread_event > 0.0).count(), px.len()),
                 q(px.iter().map(|p| p.spread_event), 1.0),
                 px.iter().filter(|p| p.n_nonfinite > 0).count(),
                 flips);
    }
    println!();
    println!("'outcome flips' is against the f64 unshared row - the reference configuration.");

    println!();
    println!("=== the drift tail, which the medians hide ===");
    println!("{:>6}{:>10}{:>12}{:>12}{:>12}{:>12}{:>14}{:>14}",
             "prec", "shared", "p50", "p90", "p99", "max", "> 1e-3", "> 1");
    for (prec, policy, px) in &store {
        println!("{:>6}{:>10}{:>12.3e}{:>12.3e}{:>12.3e}{:>12.3e}{:>14}{:>14}",
                 prec.name(),
                 if matches!(policy, RefPolicy::Shared) { "on" } else { "off" },
                 q(px.iter().map(|p| p.energy_drift_max), 0.5),
                 q(px.iter().map(|p| p.energy_drift_max), 0.9),
                 q(px.iter().map(|p| p.energy_drift_max), 0.99),
                 q(px.iter().map(|p| p.energy_drift_max), 1.0),
                 px.iter().filter(|p| p.energy_drift_max > 1e-3).count(),
                 px.iter().filter(|p| p.energy_drift_max > 1.0).count());
    }
    println!();
    println!("The last two columns are the ones that decide usability: a pixel with |dE/E| > 1");
    println!("has lost more than the total energy of the system and its trajectory means");
    println!("nothing. error_ratio flags them, which is what it is for.");

    // --- question 1 -------------------------------------------------------------------
    println!();
    println!("=== 1. does the conditioned LC branch fix spread_shape at f32? ===");
    println!("{:>6}{:>12}{:>16}{:>16}{:>16}{:>12}",
             "prec", "lc", "sp_shape med", "sp_shape p99", "sp_shape max", "nonfinite");
    let mut truth = f64::NAN;
    for (prec, lc_stable) in [
        (Precision::F64, true),
        (Precision::F64, false),
        (Precision::F32, true),
        (Precision::F32, false),
    ] {
        let cfg = EnsembleCfg { lc_stable, ..Default::default() };
        let px = render(prec, &cfg, &s);
        let med = q(px.iter().map(|p| p.spread_shape), 0.5);
        if matches!(prec, Precision::F64) && lc_stable {
            truth = med;
        }
        println!("{:>6}{:>12}{:>16.4e}{:>16.4e}{:>16.4e}{:>12}",
                 prec.name(),
                 if lc_stable { "conditioned" } else { "reference" },
                 med,
                 q(px.iter().map(|p| p.spread_shape), 0.99),
                 q(px.iter().map(|p| p.spread_shape), 1.0),
                 px.iter().filter(|p| !p.spread_shape.is_finite() || p.n_nonfinite > 0).count());
    }
    println!();
    println!("f64 conditioned is the truth column: {truth:.6e}");

    // --- question 2 -------------------------------------------------------------------
    println!();
    println!("=== 2. does the shared-reference flag matter? ===");
    println!("Read the four-combination table above. Summarised, shared against unshared:");
    for prec in [Precision::F64, Precision::F32] {
        let un = &store.iter().find(|(p, q, _)| *p == prec && matches!(q, RefPolicy::PerCopy)).unwrap().2;
        let sh = &store.iter().find(|(p, q, _)| *p == prec && matches!(q, RefPolicy::Shared)).unwrap().2;
        let r = |a: f64, b: f64| b / a;
        println!("  {:>3}: drift med x{:.3}   drift max x{:.3}   er max x{:.3}   sp_shape med x{:.4}",
                 prec.name(),
                 r(q(un.iter().map(|p| p.energy_drift_max), 0.5), q(sh.iter().map(|p| p.energy_drift_max), 0.5)),
                 r(q(un.iter().map(|p| p.energy_drift_max), 1.0), q(sh.iter().map(|p| p.energy_drift_max), 1.0)),
                 r(q(un.iter().map(|p| p.error_ratio), 1.0), q(sh.iter().map(|p| p.error_ratio), 1.0)),
                 r(q(un.iter().map(|p| p.spread_shape), 0.5), q(sh.iter().map(|p| p.spread_shape), 0.5)));
    }
    println!("Ratios are shared/unshared, so >1 means sharing made it worse.");

    // --- question 3 -------------------------------------------------------------------
    println!();
    println!("=== 3. does the branch cut reach the outcome encoding at f32? ===");
    println!("(at f64 it did not: 0 label flips of 4096 at every r_coll)");
    println!("{:>6}{:>12}{:>16}{:>16}", "prec", "r_coll/R", "label flips", "of");
    for prec in [Precision::F64, Precision::F32] {
        for r_coll_frac in [1e-4f64, 1e-3, 1e-2] {
            let stable = render(prec, &EnsembleCfg { r_coll_frac, lc_stable: true, ..Default::default() }, &s);
            let unstable = render(prec, &EnsembleCfg { r_coll_frac, lc_stable: false, ..Default::default() }, &s);
            let flips = stable
                .iter()
                .zip(unstable.iter())
                .filter(|(a, b)| a.outcome != b.outcome)
                .count();
            println!("{:>6}{r_coll_frac:>12.0e}{flips:>16}{:>16}", prec.name(), stable.len());
        }
    }
}
