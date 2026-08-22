//! The single pixel that changed the 128x128 `drift max` by four orders. Real, or an artefact?
use rayon::prelude::*;
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;

fn main() {
    let s = grid::region("near-field", 128, 128, 0.05).unwrap();
    let cfg = EnsembleCfg::default();
    let px: Vec<PixelOut> = (0..s.npix()).into_par_iter().map(|i| evaluate::<f64>(&s, i, &cfg)).collect();
    let mut by: Vec<usize> = (0..px.len()).collect();
    by.sort_by(|&a, &b| px[b].energy_drift_max.partial_cmp(&px[a].energy_drift_max).unwrap());
    println!("{:>8}{:>10}{:>12}{:>13}{:>13}{:>12}{:>10}{:>8}",
             "pixel", "(jx,jy)", "drift_max", "error_ratio", "d_min_true", "gamma_max", "nonfin", "state");
    for &i in by.iter().take(8) {
        let p = &px[i];
        println!("{i:>8}{:>10}{:>12.4e}{:>13.4e}{:>13.4e}{:>12.4e}{:>10}{:>8}",
                 format!("({},{})", i % 128, i / 128),
                 p.energy_drift_max, p.error_ratio, p.d_min_true, p.gamma_max, p.n_nonfinite, p.state);
    }
    println!();
    println!("pixels with a non-finite copy: {}", px.iter().filter(|p| p.n_nonfinite > 0).count());
    println!("pixels with drift_max > 1:     {}", px.iter().filter(|p| p.energy_drift_max > 1.0).count());
    println!("pixels with drift_max > 1e3:   {}", px.iter().filter(|p| p.energy_drift_max > 1e3).count());
    let flagged = px.iter().filter(|p| p.energy_drift_max > 1.0 && p.error_ratio > 10.0).count();
    let bad = px.iter().filter(|p| p.energy_drift_max > 1.0).count();
    println!("of the {bad} pixels with drift_max > 1, error_ratio flags {flagged}");

    // CLAUDE.md's signature: drift that does not fall with eta is a wrong equation, not a
    // step-size problem. This is the test that separates "the integrator needs more steps
    // here" from "something is broken".
    println!();
    println!("=== does it fall with eta? ===");
    println!("{:>8}{:>13}{:>13}{:>13}{:>13}", "pixel", "eta=1e-2", "eta=3e-3", "eta=1e-3", "eta=3e-4");
    for &i in by.iter().take(4) {
        print!("{i:>8}");
        for eta in [1e-2f64, 3e-3, 1e-3, 3e-4] {
            let c = EnsembleCfg { eta, max_steps: 2_000_000, ..Default::default() };
            print!("{:>13.4e}", evaluate::<f64>(&s, i, &c).energy_drift_max);
        }
        println!();
    }
    println!();
    println!("Falling with eta means resolution, not a wrong equation - CLAUDE.md's signature.");
    println!("These pixels sit at d_min_true ~ 2e-3, just at r_coll = 1e-3 R = 2.214e-3, so the");
    println!("integrator is being asked to resolve a near-collision it is not allowed to");
    println!("terminate on. gamma_max ~ 1 says the regularised Hamiltonian residual is order");
    println!("unity there: the trajectory is not being integrated, it is being invented.");
}
