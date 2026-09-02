//! **What sets the band spacing? The candidate is the bound pair's orbital phase.**
//!
//! Three-body scattering: a bound pair plus a third body. Every outcome -- which body leaves,
//! when, how close the approach -- is set by the pair's **phase at the encounter**. Move through
//! IC space and that phase winds; each 2*pi of winding is one band. The mechanism predicts one
//! orientation, every observable banded together, continuous `t_end`, and a self-similar
//! cascade, all of which is measured.
//!
//! **Two things it must survive.**
//!
//! 1. `T_pair`, measured from the trajectory, must match the per-band increment of `t_end`
//!    (1.487e-2 at the `z1` window). An earlier attempt used the period at `d_min` -- the
//!    CLOSEST APPROACH, separation 9.8e-4 -- and missed by 160x. That is the wrong orbit: the
//!    pair between encounters is O(0.03), not O(1e-3).
//! 2. Across exactly one band, the pair's phase must advance by **2*pi**. A period that merely
//!    matches in magnitude is suggestive; a phase that winds once per band is the mechanism.
//!
//! Termination is OFF and `r_coll = 0` so the series runs the full horizon -- a run stopped by
//! an event is parked at a close approach and its late samples do not exist. `n_sync` is set far
//! above the candidate frequency; that changes the step size (`dtau = eta*dt_left/(A*B)`), so
//! these trajectories are NOT bitwise the production ones, and `t_end` is printed to show the
//! discretisation did not move the physics.
use std::io::Write;

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn main() {
    let dir: String = std::env::args().nth(1).unwrap_or_else(|| "results/osc/phase".into());
    let zf: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.060);
    let nsy: usize = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(20_000);
    // Samples spanning ONE band. lambda = 12.64 px of a 256 panel at zf = 0.060.
    let nsamp: usize = std::env::args().nth(4).and_then(|s| s.parse().ok()).unwrap_or(24);
    let lam_px: f64 = std::env::args().nth(5).and_then(|s| s.parse().ok()).unwrap_or(12.64);
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0], z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (fcx, fcy, fh) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (cx, cy) = (fcx + fh * (2.0 * 0.10 - 1.0), fcy + fh * (2.0 * 0.45 - 1.0));
    let half = fh * zf;

    // Step ALONG THE BAND NORMAL. The bands run at 69.8 deg, so the normal is -20.2 deg; a cut
    // along the bands crosses no band at all and would show nothing, which is how two earlier
    // row profiles came back flat.
    let ang = 69.8f64.to_radians() - std::f64::consts::FRAC_PI_2;
    let (ux, uy) = (ang.cos(), ang.sin());
    let px_chart = 2.0 * half / 256.0;

    println!("BAND PHASE at zf = {zf}, n_sync = {nsy}, {nsamp} samples across ONE band");
    println!("  band = {lam_px} px = {:.4e} chart units; step = {:.4e}", lam_px * px_chart,
             lam_px * px_chart / nsamp as f64);
    println!();
    println!("{:>6} {:>12} {:>12} {:>10} {:>12}", "i", "t_end", "d_min", "steps", "state");

    let opts = HgOpts::<f64> {
        r_coll_frac: 0.0,
        stop_on_event: false,
        keep_boundary_shapes: true,
        ..Default::default()
    };
    let mut meta = Vec::new();
    for i in 0..nsamp {
        let f = (i as f64 / nsamp as f64 - 0.5) * lam_px * px_chart;
        let (x, y) = (cx + f * ux, cy + f * uy);
        let c = grid::decode_state(&chart, 0, x, y);
        let o = integrate_hg(c.s, &c.m, 50.0, nsy, 1e-2, 4_000_000, &opts);
        let mut fbuf = std::io::BufWriter::new(
            std::fs::File::create(format!("{dir}/shape_{i:03}.f64")).unwrap());
        for s in &o.boundary_shapes {
            for k in 0..3 {
                fbuf.write_all(&s[k].to_le_bytes()).unwrap();
            }
        }
        println!("{i:>6} {:>12.6} {:>12.4e} {:>10} {:>12}",
                 o.t, o.d_min, o.steps, if o.finite { "ok" } else { "NONFINITE" });
        meta.push(o.t);
    }
    let (lo, hi) = (meta.iter().cloned().fold(f64::INFINITY, f64::min),
                    meta.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    println!();
    println!("  t_end across the band: {lo:.6} .. {hi:.6}   span {:.4e}", hi - lo);
    println!("  If the mechanism holds, that span is ONE pair period and the shape series");
    println!("  advances by exactly 2*pi of phase from sample 0 to sample {}.", nsamp - 1);
}
