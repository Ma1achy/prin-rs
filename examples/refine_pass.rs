//! Does flag-then-re-integrate actually fix the pixels BRIEF §2.5 records?
//!
//! The 128x128 near-field grid carries 7 pixels of 16384 with `|dE/E| > 1`, worst `1.49e4`.
//! `error_ratio` flags all 7. This runs the grid with the second pass off and on, and reports
//! what it cost and what it bought.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;

const SIZE: usize = 128;

fn go(cfg: &EnsembleCfg) -> (Vec<PixelOut>, f64) {
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let t = std::time::Instant::now();
    let px = (0..s.npix())
        .into_par_iter()
        .map(|i| evaluate::<f64>(&s, i, cfg))
        .collect();
    (px, t.elapsed().as_secs_f64())
}

fn q(v: &[PixelOut], f: fn(&PixelOut) -> f64, p: f64) -> f64 {
    let mut x: Vec<f64> = v.iter().map(f).filter(|q| q.is_finite()).collect();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x[(((x.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let off = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let on = EnsembleCfg::default();
    let (a, ta) = go(&off);
    let (b, tb) = go(&on);

    println!("near-field {SIZE}x{SIZE}, t=13, E+1=8, eta=1e-2, f64");
    println!("second pass: error_ratio > {} re-integrated at eta x {}",
             on.refine_threshold, on.refine_eta_factor);
    println!();
    println!("{:>28}{:>16}{:>16}", "", "no refinement", "with refinement");
    println!("{:>28}{:>16.4e}{:>16.4e}", "drift max",
             q(&a, |p| p.energy_drift_max, 1.0), q(&b, |p| p.energy_drift_max, 1.0));
    println!("{:>28}{:>16.4e}{:>16.4e}", "drift p99",
             q(&a, |p| p.energy_drift_max, 0.99), q(&b, |p| p.energy_drift_max, 0.99));
    println!("{:>28}{:>16.4e}{:>16.4e}", "drift median",
             q(&a, |p| p.energy_drift_max, 0.5), q(&b, |p| p.energy_drift_max, 0.5));
    println!("{:>28}{:>16.4e}{:>16.4e}", "error_ratio max",
             q(&a, |p| p.error_ratio, 1.0), q(&b, |p| p.error_ratio, 1.0));
    println!("{:>28}{:>16}{:>16}", "pixels |dE/E| > 1",
             a.iter().filter(|p| p.energy_drift_max > 1.0).count(),
             b.iter().filter(|p| p.energy_drift_max > 1.0).count());
    println!("{:>28}{:>16}{:>16}", "pixels |dE/E| > 1e-3",
             a.iter().filter(|p| p.energy_drift_max > 1e-3).count(),
             b.iter().filter(|p| p.energy_drift_max > 1e-3).count());
    println!("{:>28}{:>16}{:>16}", "pixels re-integrated",
             0, b.iter().filter(|p| p.refined).count());
    println!("{:>28}{:>16.2}{:>16.2}", "wall clock (s)", ta, tb);
    println!();

    let refined: Vec<&PixelOut> = b.iter().filter(|p| p.refined).collect();
    println!("{} of {} pixels re-integrated ({:.2}% of the grid), costing {:.0}% wall clock.",
             refined.len(), b.len(),
             100.0 * refined.len() as f64 / b.len() as f64,
             100.0 * (tb / ta - 1.0));
    println!();
    println!("The worst 8 refined pixels — coarse against refined, both recorded:");
    let mut by: Vec<&PixelOut> = refined.clone();
    by.sort_by(|x, y| y.energy_drift_max_coarse.partial_cmp(&x.energy_drift_max_coarse).unwrap());
    println!("{:>14}{:>14}{:>16}{:>16}{:>12}",
             "drift coarse", "drift fine", "er coarse", "er fine", "eta used");
    for p in by.iter().take(8) {
        println!("{:>14.4e}{:>14.4e}{:>16.4e}{:>16.4e}{:>12.1e}",
                 p.energy_drift_max_coarse, p.energy_drift_max,
                 p.error_ratio_coarse, p.error_ratio, p.eta_used);
    }
    println!();
    let still = refined.iter().filter(|p| p.energy_drift_max > 1e-3).count();
    println!("refined pixels still above 1e-3 drift: {still} of {}", refined.len());
    println!("Both values are dumped for every pixel, so a refinement that did not help is");
    println!("visible rather than hidden behind a replaced number.");
}
