//! **Does the pale region resolve under magnification, or stay a one-pixel mixture?**
//!
//! Measured on the 1024^2 panel, the pale class has **no interior**: median inscribed radius
//! `1.00 px`, max `2.24`, and **0.0%** survives a 5x5 opening -- against the green ribbon's
//! median `7.00` / 83.7% and the red band's `32.25` / 93.6%. Every pale pixel is within one pixel
//! of a coloured one. So what reads as a solid wedge with a sharp edge is a **pixel-scale
//! mixture** whose *density* changes over a few pixels, not a region with an inside.
//!
//! That measurement cannot say whether the mixture is physics or arithmetic, because an
//! unresolved fractal boundary and a numerical failure scattered along a sensitive set look
//! identical at one resolution. **Resolution is the discriminator**, and it is this project's
//! standing form: *stable across resolution is the tell that it is the chart, not the grid*.
//!
//! - If the mixture is a genuine basin boundary, magnifying the window keeps the inscribed
//!   radius at ~1 px: the structure refines exactly as fast as the grid.
//! - If the features have a fixed size in chart space, the radius grows with the zoom and the
//!   1024^2 render was simply under-resolving them.
//!
//! Renders one sub-window at the same pixel count, so the cost is one render and the linear
//! resolution rises by the zoom factor. Args: `res centre_u centre_v zoom out`, where the centre
//! is in FRACTIONS of the full panel, origin top-left, so it can be read straight off the mask.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
/// `closure_render`'s own window, so the sub-window is a magnification of the committed panel and
/// not a different experiment.
const WINDOW: f64 = 0.4;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 1024);
    let fu: f64 = arg(2, 0.87);
    let fv: f64 = arg(3, 0.25);
    let zoom: f64 = arg(4, 8.0);
    let out: String = std::env::args().nth(5).unwrap_or_else(|| "wedge_zoom".into());
    let _ = std::fs::create_dir_all(&out);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx0, cy0, half0) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);

    // `Slice::axis` runs low-to-high with index, and `decode_pos` uses `idx / nx` for the row, so
    // fractional v measured from the TOP of the saved image maps to `1 - v` on the axis.
    let cx = cx0 + (2.0 * fu - 1.0) * half0;
    let cy = cy0 + (2.0 * (1.0 - fv) - 1.0) * half0;
    let half = half0 / zoom;
    let (t_max, r_coll) = (50.0, 0.005);
    let n_sync = (t_max / WINDOW).round().max(4.0) as usize;

    println!(
        "{res}^2 over a window {zoom}x smaller than the committed panel.\n\
         centre ({cx:.6}, {cy:.6})  half {half:.6}  n_sync {n_sync}  t_max {t_max}  \
         r_coll {r_coll}\n\
         cell width {:.3e} against the panel's {:.3e} -- {:.1}x finer.\n",
        2.0 * half / res as f64,
        2.0 * half0 / 1024.0,
        (2.0 * half0 / 1024.0) / (2.0 * half / res as f64),
    );

    let m = grid::decode_state(&chart, 0, cx, cy).m;
    println!("masses {m:.5?}  -- the gate\n");
    assert!((m[0] - 0.32735).abs() < 1e-4, "decode overridden");

    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        r_coll_frac: r_coll,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();

    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
    let sites = colour::landmarks(&m);
    let mut buf = Vec::with_capacity(px.len() * 3);
    let mut obuf = Vec::with_capacity(px.len() * 3);
    for p in &px {
        buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
        obuf.extend_from_slice(&png::outcome_rgb(p));
    }
    let stem = format!("{out}/wedge_z{}", zoom as i64);
    let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
    let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);
    println!(
        "{:.1}s  nonfin {}  escape {:.4} collision {:.4} bounded {:.4}  ramp ({lo:.3e}, {hi:.3e})",
        t0.elapsed().as_secs_f64(),
        px.iter().filter(|p| p.n_nonfinite > 0).count(),
        px.iter().filter(|p| p.state == 0).count() as f64 / px.len() as f64,
        px.iter().filter(|p| p.state == 2).count() as f64 / px.len() as f64,
        px.iter().filter(|p| p.state == 1).count() as f64 / px.len() as f64,
    );
    println!("wrote {stem}_uniform.png and {stem}_outcome.png");
}
