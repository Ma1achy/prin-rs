//! **Every marked box, four panels each, one consistent pass.**
//!
//! `uniform` (the production bivariate map), `outcome` (the terminal class), `drift` (the
//! diagnostic field -- *when a numerical defect is suspected, render the diagnostic field, not
//! the science field*), and `tend`.
//!
//! Each box is rendered at ITS OWN window, at `closure_render`'s settings, at one resolution for
//! all of them so the panels are comparable box to box. The zoom factor differs per box because
//! the boxes differ in size; the cell width is printed with each so no comparison is made across
//! them without it.
//!
//! **The drift and `t_end` ramps are FIXED constants shared by every box.** Auto-ranging them per
//! box would stretch each one's own p1-p99 to full scale and make sixteen incomparable pictures --
//! the standing result about auto-ranged ramps, at sixteen sites at once. The `uniform` panel keeps
//! its own auto-range because that is what the committed renders do, and its window is printed.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
const WINDOW: f64 = 0.4;
const DRIFT_LO: f64 = 1e-12;
const DRIFT_HI: f64 = 1e2;
const NAN_RGB: [u8; 3] = [255, 0, 255];

const BOXES: [(&str, f64, f64, f64); 17] = [
    ("B1", 0.890, 0.245, 0.0517), ("B2", 0.862, 0.332, 0.0571),
    ("B3", 0.379, 0.446, 0.0326), ("B4", 0.476, 0.497, 0.0294),
    ("B5", 0.590, 0.437, 0.0294), ("B6", 0.608, 0.528, 0.0337),
    ("B7", 0.419, 0.664, 0.0381), ("B8", 0.383, 0.748, 0.0403),
    ("B9", 0.807, 0.838, 0.0566), ("B10", 0.942, 0.789, 0.0522),
    ("P1", 0.539, 0.199, 0.0354), ("P2", 0.447, 0.426, 0.0354),
    ("P3", 0.533, 0.457, 0.0305), ("P4", 0.335, 0.682, 0.0408),
    ("P5", 0.428, 0.718, 0.0381), ("P6", 0.510, 0.742, 0.0397),
    ("FRAME", 0.5, 0.5, 0.5),
];

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [[0.0, 0.0, 0.015], [0.34, 0.06, 0.43], [0.72, 0.21, 0.33],
                              [0.98, 0.55, 0.04], [0.99, 1.0, 0.64]];
    let t = x.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    let mut o = [0u8; 3];
    for k in 0..3 { o[k] = (255.0 * (S[i][k] * (1.0 - f) + S[i + 1][k] * f)).clamp(0.0, 255.0) as u8; }
    o
}
/// Cyclic, so `t_end` bands are visible as bands rather than as a smooth wash -- the quantisation
/// this field is known to carry is the point of looking at it.
fn cyc(x: f64) -> [u8; 3] {
    let t = (x.clamp(0.0, 1.0) * 6.0).fract();
    let h = t * std::f64::consts::TAU;
    [(128.0 + 110.0 * h.cos()) as u8,
     (128.0 + 110.0 * (h + 2.094).cos()) as u8,
     (128.0 + 110.0 * (h + 4.189).cos()) as u8]
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(320);
    let out: String = std::env::args().nth(2).unwrap_or_else(|| "box_panels".into());
    let _ = std::fs::create_dir_all(&out);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0; q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx0, cy0, half0) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (t_max, r_coll) = (50.0f64, 0.005f64);
    let n_sync = (t_max / WINDOW).round().max(4.0) as usize;

    println!("{res}^2 per box. drift ramp FIXED at ({DRIFT_LO:e}, {DRIFT_HI:e}) log10, \
              t_end cyclic over [0, {t_max}].\n");
    for (name, u, v, h) in BOXES {
        let cx = cx0 + (2.0 * u - 1.0) * half0;
        let cy = cy0 + (2.0 * v - 1.0) * half0;
        let half = 2.0 * h * half0;
        let ens = EnsembleCfg {
            refine_flagged: false, t_max, n_sync, r_coll_frac: r_coll,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU), closure_k: 1,
            stop_on_escape: false, ..Default::default()
        };
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let t0 = std::time::Instant::now();
        let px: Vec<PixelOut> = (0..sl.npix()).into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens)).collect();
        let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
        let m = grid::decode_state(&chart, 0, cx, cy).m;
        let sites = colour::landmarks(&m);
        let (mut a, mut b, mut c, mut d) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let l10 = (DRIFT_LO.log10(), DRIFT_HI.log10());
        for p in &px {
            a.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
            b.extend_from_slice(&png::outcome_rgb(p));
            c.extend_from_slice(&if p.n_nonfinite > 0 || !p.energy_drift_max.is_finite() { NAN_RGB }
                else { ramp((p.energy_drift_max.max(1e-300).log10() - l10.0) / (l10.1 - l10.0)) });
            d.extend_from_slice(&if p.t_end.is_finite() { cyc(p.t_end / t_max) } else { NAN_RGB });
        }
        for (k, buf) in [("uniform", &a), ("outcome", &b), ("drift", &c), ("tend", &d)] {
            let _ = adaptive::save_rect(&format!("{out}/{name}_{k}.png"), res, res, buf);
        }
        println!(
            "{name:>6} {:.1}s  cell {:.3e}  nonfin {:>5}  esc {:.4} col {:.4} bnd {:.4}  \
             drift p50 {:.2e}  spread ramp ({lo:.3e}, {hi:.3e})",
            t0.elapsed().as_secs_f64(), 2.0 * half / res as f64,
            px.iter().filter(|p| p.n_nonfinite > 0).count(),
            px.iter().filter(|p| p.state == 0).count() as f64 / px.len() as f64,
            px.iter().filter(|p| p.state == 2).count() as f64 / px.len() as f64,
            px.iter().filter(|p| p.state == 1).count() as f64 / px.len() as f64,
            { let mut v: Vec<f64> = px.iter().map(|p| p.energy_drift_max)
                .filter(|x| x.is_finite()).collect();
              v.sort_by(|x, y| x.partial_cmp(y).unwrap()); v[v.len() / 2] },
        );
    }
}
