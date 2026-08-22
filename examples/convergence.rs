//! **Experiment B (BRIEF §8): which conclusions survive at large `n`?**
//!
//! Every prior measurement in this project used 16–64 trajectories across 4–8 regions, and
//! conclusions repeatedly flipped when re-tested at larger `n` — a field excluded on a 1.2x
//! effect turned out to be 18.8x. This measures the instability directly.
//!
//! Two things, and the second is the one that matters:
//!
//! 1. **The same slice at increasing resolution.** Reports where each statistic settles. But a
//!    resolution sweep confounds two effects: more samples, *and* a different physical
//!    quantity, since the jitter scales with the cell width and finer cells sample the region
//!    differently. It cannot separate "the estimate converged" from "the thing being estimated
//!    changed".
//!
//! 2. **Subsampling one fixed grid.** Draw `n` pixels from the finest grid and report how far a
//!    statistic computed on `n` pixels scatters from its full-grid value. Same physical
//!    quantity throughout, so this isolates sampling error — and it answers the question
//!    directly: *would a conclusion drawn from `n` pixels have held?*

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::rng::SplitMix64;

const FINEST: usize = 128;

/// **Refinement off throughout.** This experiment characterises the kernel whose behaviour
/// motivated the second pass; measuring the repaired kernel would hide the thing being
/// measured. `examples/refine_pass.rs` reports the repaired numbers.
fn render(size: usize) -> Vec<PixelOut> {
    let s = grid::region("near-field", size, size, 0.05).unwrap();
    let cfg = EnsembleCfg { refine_flagged: false, ..Default::default() };
    (0..s.npix())
        .into_par_iter()
        .map(|i| evaluate::<f64>(&s, i, &cfg))
        .collect()
}

fn quant(v: &mut Vec<f64>, f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * f).round() as usize]
}

/// The statistics a conclusion would actually be drawn from.
fn summarise(px: &[PixelOut]) -> Vec<(&'static str, f64)> {
    let g = |f: fn(&PixelOut) -> f64| -> Vec<f64> {
        px.iter().map(f).filter(|x| x.is_finite()).collect()
    };
    let mut drift = g(|p| p.energy_drift_max);
    let mut er = g(|p| p.error_ratio);
    let mut ss = g(|p| p.spread_shape);
    let mut dm = g(|p| p.d_min_true);
    let n = px.len() as f64;
    vec![
        ("drift median", quant(&mut drift.clone(), 0.5)),
        ("drift p99", quant(&mut drift.clone(), 0.99)),
        ("drift max", quant(&mut drift, 1.0)),
        ("error_ratio median", quant(&mut er.clone(), 0.5)),
        ("error_ratio p99", quant(&mut er, 0.99)),
        ("spread_shape median", quant(&mut ss.clone(), 0.5)),
        ("spread_shape p99", quant(&mut ss, 0.99)),
        ("d_min_true median", quant(&mut dm, 0.5)),
        ("frac collision", px.iter().filter(|p| p.state == 2).count() as f64 / n),
        ("frac er > 10", px.iter().filter(|p| p.error_ratio > 10.0).count() as f64 / n),
        ("frac drift > 1e-3", px.iter().filter(|p| p.energy_drift_max > 1e-3).count() as f64 / n),
        ("frac spread_event > 0", px.iter().filter(|p| p.spread_event > 0.0).count() as f64 / n),
    ]
}

fn main() {
    println!("Experiment B: statistical convergence. near-field, t=13, E+1=8, eta=0.01, f64.");
    println!();
    println!("=== 1. the same slice at increasing resolution ===");
    let sizes = [8usize, 16, 32, 64, FINEST];
    let runs: Vec<Vec<(&str, f64)>> = sizes.iter().map(|&n| summarise(&render(n))).collect();

    print!("{:>24}", "quantity");
    for n in sizes {
        print!("{:>14}", format!("{n}x{n}"));
    }
    println!("{:>12}", "|last/64-1|");
    for (k, (name, _)) in runs[0].iter().enumerate() {
        print!("{name:>24}");
        for r in &runs {
            print!("{:>14.4e}", r[k].1);
        }
        let a = runs[runs.len() - 2][k].1;
        let b = runs[runs.len() - 1][k].1;
        println!("{:>12.3}", (b / a - 1.0).abs());
    }
    println!();
    println!("The last column is how much the 64x64 answer moved on going to {FINEST}x{FINEST}.");
    println!("It is NOT a convergence error: the jitter scales with cell width, so the finer");
    println!("grid measures a different physical ensemble. Read column 2 for that.");

    // --- 2 ----------------------------------------------------------------------------
    println!();
    println!("=== 2. subsampling one fixed grid — the same physical quantity throughout ===");
    let full = render(FINEST);
    let truth = summarise(&full);
    println!("truth = the full {}x{} grid ({} pixels).", FINEST, FINEST, full.len());
    println!("Each cell is the interdecile spread of the statistic over 200 random draws of n");
    println!("pixels, as a fraction of the truth. Below ~0.1 the conclusion is stable.");
    println!();

    print!("{:>24}{:>14}", "quantity", "truth");
    let ns = [16usize, 64, 256, 1024, 4096];
    for n in ns {
        print!("{:>12}", format!("n={n}"));
    }
    println!();

    for (k, (name, t)) in truth.iter().enumerate() {
        print!("{name:>24}{t:>14.4e}");
        for n in ns {
            let mut vals = Vec::with_capacity(200);
            let mut rng = SplitMix64::new(0xC0FFEE ^ (n as u64) ^ ((k as u64) << 32));
            for _ in 0..200 {
                let sample: Vec<PixelOut> = (0..n)
                    .map(|_| full[(rng.next_u64() as usize) % full.len()].clone())
                    .collect();
                vals.push(summarise(&sample)[k].1);
            }
            let lo = quant(&mut vals.clone(), 0.1);
            let hi = quant(&mut vals, 0.9);
            let rel = if t.abs() > 0.0 { (hi - lo) / t.abs() } else { f64::NAN };
            print!("{rel:>12.3}");
        }
        println!();
    }
    println!();
    println!("A value of 1.0 means the interdecile scatter is as wide as the quantity itself:");
    println!("two independent studies at that n would routinely disagree by a factor of two or");
    println!("more, which is how a 1.2x effect turns out to be 18.8x.");
}
