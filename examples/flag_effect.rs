//! Do `lc_stable` and `ref_policy` actually change anything at f64?
//!
//! Several rows of `examples/f32_report.rs` are identical to every printed digit, which is
//! either a real result or a flag that is not reaching the kernel. This checks per pixel
//! rather than reading a quantile — the answer is that both flags reach it and the
//! *distribution* genuinely does not move, which is a different and more interesting claim.
//!
//! It also surfaces what an aggregate hides: shared references change `spread_shape` by up to
//! 1.9x on an individual pixel while shifting the median by 1%. That is exactly the failure
//! NOTES §1 warned about — a difference of aggregates cannot see it, which is why
//! `ref_disagree` is dumped per pixel.
use rayon::prelude::*;
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::integrate::az::RefPolicy;

fn main() {
    let s = grid::region("near-field", 32, 32, 0.05).unwrap();
    let go = |cfg: EnsembleCfg| -> Vec<PixelOut> {
        (0..s.npix()).into_par_iter().map(|i| evaluate::<f64>(&s, i, &cfg)).collect()
    };
    let base = go(EnsembleCfg::default());

    for (name, cfg) in [
        ("lc_stable=false", EnsembleCfg { lc_stable: false, ..Default::default() }),
        ("shared refs", EnsembleCfg { ref_policy: RefPolicy::Shared, ..Default::default() }),
    ] {
        let v = go(cfg);
        let n_drift = base.iter().zip(&v).filter(|(a, b)| a.energy_drift_max != b.energy_drift_max).count();
        let n_shape = base.iter().zip(&v).filter(|(a, b)| a.spread_shape != b.spread_shape).count();
        let n_er = base.iter().zip(&v).filter(|(a, b)| a.error_ratio != b.error_ratio).count();
        let worst_shape = base
            .iter()
            .zip(&v)
            .map(|(a, b)| (a.spread_shape - b.spread_shape).abs() / a.spread_shape.abs().max(1e-300))
            .fold(0.0f64, f64::max);
        println!("f64, {name}:");
        println!("  pixels whose drift changed at all: {n_drift}/{}", base.len());
        println!("  pixels whose spread_shape changed: {n_shape}/{}", base.len());
        println!("  pixels whose error_ratio changed:  {n_er}/{}", base.len());
        println!("  worst relative change in spread_shape: {worst_shape:.3e}");
        if worst_shape > 0.5 {
            println!("  ^ an individual pixel moved by {:.0}% while the median moved by ~1%.",
                     100.0 * worst_shape);
        }
    }
}
