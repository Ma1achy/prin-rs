//! Does the step budget, sized at `t_max = 13`, still hold at the Principia corpus's default
//! horizon of 50 -- and what is the step RATE a live playhead would need to budget instead?
//!
//! The production budget is a whole-trajectory count (`max_steps = 480_000`, measured at
//! `t_max = 13`, `results/step_budget/`). A live playhead has no fixed `t_max`, so a whole-run
//! count is not even defined there; a per-unit-time rate is. This measures, per pixel, the
//! substeps per copy per unit physical time, on the step budget's two worst targets and the
//! clean control, at 13 and at 50, with `n_sync` scaled with `t_max` (the standing rule).
//!
//! Arms per target, one fixed grid:
//!   A  t = 13, 480k     the production cell (control: must reproduce the sweep's 0 flagged)
//!   B  t = 50, 480k     the same budget at the corpus horizon
//!   C  t = 50, 8M       effectively unbounded: the true step need at 50

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use rayon::prelude::*;

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(24);
    let mut targets: Vec<(String, Chart, f64, f64, f64)> = Vec::new();
    let (c, cx, cy, h) = Chart::config_stability();
    targets.push(("config_stability".into(), c, cx, cy, h));
    for (n, c, cx, cy, h) in grid::gallery_cases() {
        if n == "burrau_nu_k" || n == "latent_mixed_h3" || n == "latent_shape" {
            targets.push((n.into(), c, cx, cy, h));
        }
    }
    println!("res {res}^2; per pixel, steps are summed over the E+1 copies and divided by E+1");
    println!("{:>18} {:>6} {:>9} {:>8} {:>8} {:>11} {:>11} {:>11} {:>11}",
        "target", "t_max", "budget", "flagged", "budget!", "steps/copy p50", "p99", "max", "rate p99 /t");
    for (name, chart, cx, cy, half) in &targets {
        for &(t_max, ms) in &[(13.0, 480_000usize), (50.0, 480_000), (50.0, 8_000_000)] {
            let base = EnsembleCfg::production();
            let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round() as usize;
            let ens = EnsembleCfg { t_max, n_sync, max_steps: ms, ..base };
            let sl = grid::Slice::body_plane(res, res, *cx, *cy, *half, 0).with_chart(*chart);
            let px: Vec<PixelOut> =
                (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();
            let e1 = (ens.n_extra + 1) as f64;
            let flagged = px.iter().filter(|p| p.n_nonfinite > 0).count();
            let budget = px.iter().filter(|p| p.budget_exhausted).count();
            let mut per: Vec<f64> = px.iter().map(|p| p.total_substeps as f64 / e1).collect();
            let p50 = q(&mut per, 0.5);
            let p99 = q(&mut per, 0.99);
            let mx = *per.last().unwrap();
            println!("{:>18} {:>6} {:>9} {:>8} {:>8} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e}",
                name, t_max, ms, flagged, budget, p50, p99, mx, p99 / t_max);
        }
    }
}
