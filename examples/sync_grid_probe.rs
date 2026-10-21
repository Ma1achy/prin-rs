//! Can the sync-boundary grid be made as fine as a display frame under Heggie, and what does it
//! cost?
//!
//! The Principia corpus marches a lockstep frame loop in macro-steps, while the regularised
//! integrators land on a separate sync grid `t_max / n_sync`. Under Heggie a boundary is "a
//! sampling and event cadence, and nothing more" (no re-registration), so a fine fixed grid is a
//! candidate for the playhead's quantum. This holds the PHYSICAL step size fixed (`eta` scaled
//! with `n_sync`, as the controlled re-registration measurement did) and the closure window fixed
//! in time (`closure_k` scaled with `n_sync`), so the only change is how often the march stops to
//! land and sample.
//!
//! Per arm against n_sync = 32 (production): substeps (cost), outcome-label flips, and the
//! shape-vector chord (p50, max) of the nominal copy.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use rayon::prelude::*;

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(32);
    let mut targets: Vec<(String, Chart, f64, f64, f64, usize)> = Vec::new();
    let (c, cx, cy, h) = Chart::config_stability();
    targets.push(("config_stability".into(), c, cx, cy, h, 0));
    for &(n, cx, cy, body) in grid::REGIONS.iter() {
        if n == "near-field" || n == "far" {
            targets.push((n.into(), Chart::BodyPlane, cx, cy, 0.05, body));
        }
    }
    let base = EnsembleCfg::production();
    println!("res {res}^2, t_max {}, physical step held fixed (eta x n/32), closure window held fixed in time", base.t_max);
    println!("{:>16} {:>6} {:>9} {:>12} {:>9} {:>10} {:>10} {:>8}",
        "target", "n_sync", "dt_sync", "substeps", "x cost", "flips", "chord p50", "max");
    for (name, chart, cx, cy, half, body) in &targets {
        let sl = grid::Slice::body_plane(res, res, *cx, *cy, *half, *body).with_chart(*chart);
        let mut reference: Option<(Vec<PixelOut>, u64)> = None;
        for &n in &[32usize, 128, 512, 1040] {
            let ens = EnsembleCfg {
                n_sync: n,
                eta: base.eta * n as f64 / 32.0,
                closure_k: base.closure_k * n / 32,
                ..base.clone()
            };
            let px: Vec<PixelOut> =
                (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();
            let steps: u64 = px.iter().map(|p| p.total_substeps).sum();
            match &reference {
                None => {
                    println!("{:>16} {:>6} {:>9.5} {:>12} {:>9} {:>10} {:>10} {:>8}", name, n, base.t_max / n as f64, steps, "1.000", "-", "-", "-");
                    reference = Some((px, steps));
                }
                Some((r, s0)) => {
                    let flips = px.iter().zip(r).filter(|(a, b)| a.outcome != b.outcome).count();
                    let mut ch: Vec<f64> = px.iter().zip(r).map(|(a, b)| {
                        let d = (0..3).map(|i| (a.shape_vec[i] - b.shape_vec[i]).powi(2)).sum::<f64>().sqrt();
                        if d.is_finite() { d } else { 2.0 }
                    }).collect();
                    let p50 = q(&mut ch, 0.5);
                    let mx = *ch.last().unwrap();
                    println!("{:>16} {:>6} {:>9.5} {:>12} {:>9.3} {:>10} {:>10.2e} {:>8.3}",
                        name, n, base.t_max / n as f64, steps, steps as f64 / *s0 as f64, format!("{flips}/{}", px.len()), p50, mx);
                }
            }
        }
    }
}
