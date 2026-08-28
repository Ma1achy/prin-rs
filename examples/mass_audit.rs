//! **The mass path, audited over every pixel and every copy the integrator receives.**
//!
//! If any path still reads a global mass instead of the decoded one, it integrates a *different*
//! three-body problem from the one the initial conditions were built for: `sum p != 0`, the
//! centre of mass drifts, and AZ reconstructs assuming it at the origin. That failure is
//! **invisible on an equal-mass slice** and severe on an unequal one — `config_stability`
//! decodes to `(0.32735, 0.42763, 0.24502)` while every preset sits at `(1/3, 1/3, 1/3)`.
//!
//! The audit is on `jitter::copies_with_path` called with **exactly** `evaluate_at`'s arguments,
//! because `evaluate_at`'s next lines are `integrate_az_opts(c.s, &c.m, ..)` and
//! `energy(&c.s.r, &c.s.v, &c.m, ..)` — the copy's own mass, straight through. So this is the
//! integrator's input and not a parallel reconstruction of it.
//!
//! Four quantities per copy, reported as the max over the frame:
//!
//! - `|m - m_expect|` — the decode gives the masses the chart says it does;
//! - `|sum m|` — they are normalised;
//! - `|sum m_i v_i|` — total momentum is zero, so the COM does not drift;
//! - `|sum m_i r_i|` — the COM is *at* the origin, which is what AZ's reconstruction assumes.
//!
//! The last two are separate assertions on purpose. A construction that assumes a COM-centred
//! input returns a drifting system without one, and zero momentum does not imply zero first
//! moment — that is a standing result in this project at a different site.
use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::physics::Ic;

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
const M_EXPECT: [f64; 3] = [0.327_35, 0.427_63, 0.245_02];

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

    // `config_stability` first, then the presets as the control: the failure being looked for is
    // one that CANNOT show on an equal-mass slice, so a preset reading clean says nothing on its
    // own and reading dirty would be a second, larger finding.
    let mut cases: Vec<(String, Chart, f64, f64, f64, Option<[f64; 3]>)> =
        vec![("config_stability".into(), cs, cx, cy, half, Some(M_EXPECT))];
    for w in ["preset_shape", "preset_prho", "preset_plambda", "preset_shape_pl"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            cases.push((w.into(), c.1, c.2, c.3, c.4, Some([1.0 / 3.0; 3])));
        }
    }

    let cfg = EnsembleCfg::default();
    println!(
        "{res}^2, every pixel, all {} copies -- {} systems per case.\n\
         Built by `jitter::copies_with_path` with `evaluate_at`'s own arguments.\n",
        cfg.n_extra + 1,
        res * res * (cfg.n_extra + 1),
    );
    println!(
        "{:>18} {:>14} {:>11} {:>11} {:>11} {:>11} {:>11}",
        "case", "masses (px 0)", "max|dm|", "max|sum m-1|", "max|sum p|", "max|sum m r|", "m spread"
    );

    let mut bad = 0usize;
    for (name, chart, cx, cy, half, expect) in cases {
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let acc = (0..sl.npix())
            .into_par_iter()
            .map(|k| {
                let cs: Vec<Ic<f64>> = jitter::copies_with_path::<f64>(
                    &sl, k, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme,
                    cfg.decode_path,
                );
                let mut o = [0.0f64; 5];
                for c in &cs {
                    let ms: f64 = c.m.iter().sum();
                    let mut p = [0.0f64; 2];
                    let mut r = [0.0f64; 2];
                    for i in 0..3 {
                        p[0] += c.m[i] * c.s.v[i].x;
                        p[1] += c.m[i] * c.s.v[i].y;
                        r[0] += c.m[i] * c.s.r[i].x;
                        r[1] += c.m[i] * c.s.r[i].y;
                    }
                    if let Some(e) = expect {
                        for i in 0..3 {
                            o[0] = o[0].max((c.m[i] - e[i]).abs());
                        }
                    }
                    o[1] = o[1].max((ms - 1.0).abs());
                    o[2] = o[2].max(p[0].hypot(p[1]));
                    o[3] = o[3].max(r[0].hypot(r[1]));
                    // How much the masses move ACROSS the footprint. On a configuration chart
                    // this is exactly zero, which is what makes `copies[0].m` harmless there;
                    // on a mass chart it is not, and the shortcut would be a real approximation.
                    for i in 0..3 {
                        o[4] = o[4].max((c.m[i] - cs[0].m[i]).abs());
                    }
                }
                o
            })
            .reduce(
                || [0.0f64; 5],
                |a, b| [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2]), a[3].max(b[3]), a[4].max(b[4])],
            );

        let m0 = sl.decode_state(cx, cy).m;
        println!(
            "{name:>18} {:>14} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e}",
            format!("{:.5},{:.5}", m0[0], m0[1]),
            acc[0], acc[1], acc[2], acc[3], acc[4],
        );
        // 1e-4 on the masses matches the precision the reference UI quotes them to; the three
        // conservation quantities are exact identities of the construction and get round-off
        // tolerances, not physical ones.
        if acc[0] > 1e-4 || acc[1] > 1e-12 || acc[2] > 1e-12 || acc[3] > 1e-12 {
            println!("{:>18}   ^^ FAILS", "");
            bad += 1;
        }
    }
    println!();
    if bad == 0 {
        println!(
            "PASS. The masses the integrator receives are the decoded ones on every pixel and\n\
             every copy, they sum to 1, total momentum is zero and the COM is at the origin.\n\
             **The mass path is not the bug**, and an equal-mass control cannot have shown that."
        );
    } else {
        println!("{bad} case(s) FAIL -- the mass path is wrong and no bisect is needed.");
        std::process::exit(1);
    }
}
