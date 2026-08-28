//! Does `hot` fall when `eta` is halved on THIS slice? At 86% of the frame above `1e-6` it
//! should, if that is truncation error. If it does not, the field is not converged at
//! horizon 50 and no code change reaches it -- a different finding, and one that needs saying.
//!
//! `hot`, drift and the state census are PER-TRAJECTORY statistics, so a coarse grid does not
//! bias them the way it biases a chord ratio. The chord columns here are between `eta` rungs at
//! ONE resolution, and are labelled with it.
use rayon::prelude::*;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn qt(v: &mut Vec<f64>, p: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0; q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx, cy, half) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let st = grid::decode_state(&chart, 0, cx, cy);
    println!("masses {:.5?}  -- the gate; equal masses would mean the decode is overridden\n", st.m);

    println!("{res}^2, HEAD, horizon 50, r_coll 0.005, n_sync 32 held FIXED across the ladder\n\
              (scaling it with eta would compare different discretisations).\n");
    println!("{:>10} {:>8} {:>8} {:>10} {:>10} {:>8} {:>8} {:>8} {:>9} {:>10} {:>10} {:>11} {:>11}",
        "eta", "nonfin", "hot", "drift p50", "drift p99", "escape", "bounded", "collis",
        "|n0|>0.9", "steps/copy", "dt p50", "chord p50", "chord max");
    let mut prev: Option<Vec<[f64; 3]>> = None;
    for eta in [1e-2f64, 5e-3, 2.5e-3, 1.25e-3] {
        let ens = EnsembleCfg {
            t_max: 50.0, n_sync: 32, r_coll_frac: 0.005, eta,
            refine_flagged: false, ..Default::default()
        };
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let px: Vec<PixelOut> = (0..sl.npix()).into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens)).collect();
        let n = px.len() as f64;
        let mut dv: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
        let (p50, p99) = (qt(&mut dv, 0.5), qt(&mut dv, 0.99));
        let sv: Vec<[f64; 3]> = px.iter().map(|p| p.shape_vec).collect();
        let (mut cp50, mut cmax) = (f64::NAN, f64::NAN);
        if let Some(pv) = &prev {
            let mut c: Vec<f64> = pv.iter().zip(sv.iter())
                .map(|(a, b)| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt())
                .filter(|x| x.is_finite()).collect();
            cmax = c.iter().cloned().fold(0.0, f64::max);
            cp50 = qt(&mut c, 0.5);
        }
        println!("{eta:>10.2e} {:>8} {:>8.4} {p50:>10.3e} {p99:>10.3e} {:>8.4} {:>8.4} {:>8.4} {:>9.4} {:>10} {:>10.3e} {cp50:>11.3e} {cmax:>11.3e}",
            px.iter().filter(|p| p.n_nonfinite > 0).count(),
            px.iter().filter(|p| !(p.energy_drift_max <= 1e-6)).count() as f64 / n,
            px.iter().filter(|p| p.state == 0).count() as f64 / n,
            px.iter().filter(|p| p.state == 1).count() as f64 / n,
            px.iter().filter(|p| p.state == 2).count() as f64 / n,
            px.iter().filter(|p| p.shape_vec[0].abs() > 0.9).count() as f64 / n,
            { let ncopy = (ens.n_extra + 1) as f64;
              let mut v: Vec<f64> = px.iter().map(|p| p.total_substeps as f64 / ncopy).collect();
              qt(&mut v, 0.5) },
            // `total_substeps` is summed over ALL E+1 copies, so the per-copy step count is
            // `total_substeps / (n_extra + 1)`. Dividing by the raw sum understates the physical
            // step EIGHT-fold -- which matters, because the number is being compared against the
            // reference's `dtMacro = 0.002` to settle the collision-cadence question.
            { let ncopy = (ens.n_extra + 1) as f64;
              let mut v: Vec<f64> = px.iter().filter(|p| p.total_substeps > 0)
                .map(|p| p.t_end * ncopy / p.total_substeps as f64).filter(|x| x.is_finite()).collect();
              if v.is_empty() { f64::NAN } else { qt(&mut v, 0.5) } });
        prev = Some(sv);
    }
}
