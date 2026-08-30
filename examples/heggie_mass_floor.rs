//! **How small does a decoded mass get, and how large does `1/m` get in `Gamma*`?**
//!
//! Heggie's regularised Hamiltonian carries `1/m_1, 1/m_2, 1/m_3` in its coupling term, and his
//! §4 states plainly that Eq. (21) is inapplicable if any mass vanishes — that is why he needs a
//! separate treatment for the restricted problem. AZ has no such term, so this is a risk the
//! Heggie port introduces and the AZ port never had.
//!
//! `Chart::Latent` decodes masses through a softmax over logits saturated at
//! `MU_MAX * (2 sigmoid(z) - 1)`, so nothing is exactly zero by construction. **That is an
//! argument, and this is the measurement.** The question is not whether `1/m` is finite; it is
//! whether it is large enough to cost precision in a sum against terms of order one.
//!
//! Read-only: it decodes and reports, and integrates nothing. Cheap enough to run at the
//! shipping resolution, which matters — this project has both understated a maximum eightfold and
//! overstated a median twenty-sixfold by reading a coarse grid.
//!
//! The presets are the **control** and they carry nothing: on an equal-mass slice `m = 1/3` and
//! `1/m = 3` exactly, so a preset reading clean says nothing on its own. `config_stability` is
//! the case with a live answer.
use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::physics::Ic;

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(1024);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1],
        z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]],
        z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let cs = Chart::Latent { z0, q1, q2 };
    let (cx, cy, half) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);

    let mut cases: Vec<(String, Chart, f64, f64, f64)> =
        vec![("config_stability".into(), cs, cx, cy, half)];
    for w in ["preset_shape", "preset_prho", "preset_plambda", "preset_shape_pl", "shape_sphere"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            cases.push((w.into(), c.1, c.2, c.3, c.4));
        }
    }

    let cfg = EnsembleCfg::production();
    println!(
        "{res}^2, every pixel, all {} copies -- {} systems per case.\n\
         `mu` is Heggie's reduced mass mu_jk; `1/m` is the coefficient of the coupling term\n\
         `- p_j . p_k / m_i` in his Eq. (6) and of `-(1/4)(R_i/m_i) W_j . W_k` in Gamma*.\n",
        cfg.n_extra + 1,
        res * res * (cfg.n_extra + 1),
    );
    println!(
        "{:>18} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "case", "min m", "max 1/m", "min mu", "max ratio", "nonfinite"
    );

    for (name, chart, cx, cy, half) in cases {
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let acc = (0..sl.npix())
            .into_par_iter()
            .map(|k| {
                let cs: Vec<Ic<f64>> = jitter::copies_with_path::<f64>(
                    &sl, k, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme,
                    cfg.decode_path,
                );
                // [min m, min mu, max m_max/m_min, nonfinite]
                let mut o = [f64::INFINITY, f64::INFINITY, 0.0, 0.0];
                for c in &cs {
                    if !c.m.iter().all(|x| x.is_finite() && *x > 0.0) {
                        o[3] += 1.0;
                        continue;
                    }
                    let lo = c.m.iter().cloned().fold(f64::INFINITY, f64::min);
                    let hi = c.m.iter().cloned().fold(0.0f64, f64::max);
                    o[0] = o[0].min(lo);
                    o[2] = o[2].max(hi / lo);
                    for i in 0..3 {
                        let (j, k2) = ((i + 1) % 3, (i + 2) % 3);
                        o[1] = o[1].min(c.m[j] * c.m[k2] / (c.m[j] + c.m[k2]));
                    }
                }
                o
            })
            .reduce(
                || [f64::INFINITY, f64::INFINITY, 0.0, 0.0],
                |a, b| [a[0].min(b[0]), a[1].min(b[1]), a[2].max(b[2]), a[3] + b[3]],
            );

        println!(
            "{name:>18} {:>12.4e} {:>12.4e} {:>12.4e} {:>12.3} {:>10}",
            acc[0],
            1.0 / acc[0],
            acc[1],
            acc[2],
            acc[3] as u64
        );
    }

    println!(
        "\nHOW TO READ IT. `1/m` is a coefficient, not a singularity: it multiplies a term that\n\
         sits in a sum with terms of order `m_j m_k R_j R_k`. What costs precision is the RATIO\n\
         between the largest and smallest term in that sum, which is what `max ratio` bounds. A\n\
         `1/m` of ~10 alongside `min mu` of the same order is unremarkable; a `1/m` of 1e5 would\n\
         mean the coupling term dominates Gamma* by five orders on those pixels and the\n\
         `gamma_residual` normalisation would be reading that term and nothing else.\n\
         The presets sit at exactly `m = 1/3, 1/m = 3, mu = 1/6` and are the control."
    );
}
