//! Two things BRIEF §8 needs measured before its experiments can be interpreted.
//!
//! 1. `error_ratio = sigma_E(t)/sigma_E(0)` and `sigma_E(0)` is proportional to the jitter,
//!    hence to cell width. Refine the grid and the ratio inflates for a purely trivial
//!    reason. Experiment 1 (synthesising parents by 2x2 aggregation) compares across
//!    different effective cell widths, so it is directly exposed.
//!
//! 2. `ensemble_spread` should be free of this: `spread_shape` normalises by the chord bound
//!    2, a geometric constant, and carries no dependence on sigma_E at all. Confirmed with a
//!    number rather than an argument.
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

fn q(v: &mut Vec<f64>, f: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * f).round() as usize]
}

fn main() {
    println!("Refining the grid at fixed box size: cell width halves at each row.\n");
    println!(
        "{:>6}{:>12}{:>14}{:>14}{:>14}{:>14}",
        "size", "cell width", "sigma_E(0)", "sigma_E(t)", "error_ratio", "ens_spread"
    );
    let mut first: Option<(f64, f64, f64)> = None;
    for size in [8usize, 16, 32, 64] {
        let s = grid::region("near-field", size, size, 0.05).unwrap();
        let cfg = EnsembleCfg { t_max: 4.0, n_sync: 10, refine_flagged: false, ..Default::default() };
        let px: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect();
        let (hx, _) = s.cell_widths();

        let mut s0: Vec<f64> = px.iter().map(|p| p.sigma_e_0).filter(|x| x.is_finite()).collect();
        let mut st: Vec<f64> = px.iter().map(|p| p.sigma_e_t).filter(|x| x.is_finite()).collect();
        let mut er: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
        let mut es: Vec<f64> = px.iter().map(|p| p.ensemble_spread).filter(|x| x.is_finite()).collect();
        let (m0, mt, me, mes) = (q(&mut s0, 0.5), q(&mut st, 0.5), q(&mut er, 0.5), q(&mut es, 0.5));
        println!("{size:>6}{hx:>12.4e}{m0:>14.4e}{mt:>14.4e}{me:>14.6}{mes:>14.4e}");
        if first.is_none() {
            first = Some((m0, me, mes));
        }
    }
    let (f0, fe, fes) = first.unwrap();
    println!();
    println!("sigma_E(0) tracks cell width, as it must - it is the jitter. The ratio is the");
    println!("quantity that must not be compared across resolutions without correcting for");
    println!("it; the dump carries sigma_E(0) and sigma_E(t) separately so it can be.");
    println!();
    println!("first-row references: sigma_E(0) {f0:.4e}, error_ratio {fe:.6}, ens_spread {fes:.4e}");

    println!();
    println!("=== ref_disagree vs error_ratio (NOTES §1, unshared run) ===");
    let s = grid::region("near-field", 24, 24, 0.05).unwrap();
    let cfg = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let px: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect();
    let with: Vec<&_> = px.iter().filter(|p| p.ref_disagree > 0).collect();
    let without: Vec<&_> = px.iter().filter(|p| p.ref_disagree == 0).collect();
    let med = |v: &[&prin_rs::ensemble::pixel::PixelOut], f: fn(&prin_rs::ensemble::pixel::PixelOut) -> f64| {
        let mut x: Vec<f64> = v.iter().map(|p| f(p)).filter(|y| y.is_finite()).collect();
        if x.is_empty() { return f64::NAN; }
        x.sort_by(|a, b| a.partial_cmp(b).unwrap());
        x[x.len() / 2]
    };
    println!("pixels where copies disagreed on the reference body: {} of {}", with.len(), px.len());
    println!("  median error_ratio   with disagreement {:.6}   without {:.6}",
             med(&with, |p| p.error_ratio), med(&without, |p| p.error_ratio));
    println!("  median drift_max     with disagreement {:.3e}   without {:.3e}",
             med(&with, |p| p.energy_drift_max), med(&without, |p| p.energy_drift_max));
    println!("  median ens_spread    with disagreement {:.4e}   without {:.4e}",
             med(&with, |p| p.ensemble_spread), med(&without, |p| p.ensemble_spread));

    println!();
    println!("=== d_min_gap across regions (NOTES §2.1) ===");
    println!("{:>14}{:>16}{:>16}", "region", "median gap", "max gap");
    for region in ["near-field", "mid-field", "body2 core", "body1 slice", "far"] {
        let s = grid::region(region, 16, 16, 0.05).unwrap();
        let px: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &EnsembleCfg { refine_flagged: false, ..Default::default() })).collect();
        let mut g: Vec<f64> = px.iter().map(|p| p.d_min_gap).filter(|x| x.is_finite()).collect();
        println!("{region:>14}{:>16.3e}{:>16.3e}", q(&mut g.clone(), 0.5), q(&mut g, 1.0));
    }
    println!();
    println!("d_min_ref tracks only the two regularised pairs; d_min_true includes the");
    println!("unregularised side. A gap of zero means the reference's blind spot never bit.");
}
