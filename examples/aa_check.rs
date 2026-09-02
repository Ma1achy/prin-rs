//! **The ensemble samples the renderer was throwing away.** `rgb` against `rgb_aa`, one march.
//!
//! Every footprint carries `E+1` copies jittered across the whole cell. `spread_shape` reduces
//! over all of them and sets the lightness; the hue reads `shapes[0]` alone. This renders the
//! same march both ways -- **one integration, two colourings** -- so any difference is the
//! sampling and nothing else.
//!
//! Banding is measured the same way `examples/moire.rs` does: per row, detrend, power at each
//! integer wavelength 2..64 px, averaged over rows, peak reported with its **prominence**. A peak
//! at prominence ~1 is the largest bin of a flat spectrum and is not a band.
//!
//! ```text
//! cargo run --release --example aa_check -- [res] [out_dir] [zoom_frac]
//! ```
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn banding(v: &[f64], n: usize) -> (usize, f64) {
    let mut power = vec![0.0f64; 65];
    let mut rows = 0usize;
    for y in 0..n {
        let row: Vec<f64> = (0..n).map(|x| v[y * n + x]).collect();
        if row.iter().any(|x| !x.is_finite()) {
            continue;
        }
        let m = row.iter().sum::<f64>() / n as f64;
        let mid = (n as f64 - 1.0) / 2.0;
        let (mut sxy, mut sxx) = (0.0f64, 0.0f64);
        for (i, &r) in row.iter().enumerate() {
            let dx = i as f64 - mid;
            sxy += dx * (r - m);
            sxx += dx * dx;
        }
        let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
        let d: Vec<f64> =
            row.iter().enumerate().map(|(i, &r)| r - m - slope * (i as f64 - mid)).collect();
        for lam in 2..=64usize {
            let w = std::f64::consts::TAU / lam as f64;
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (i, &x) in d.iter().enumerate() {
                re += x * (w * i as f64).cos();
                im += x * (w * i as f64).sin();
            }
            power[lam] += (re * re + im * im) / (n * n) as f64;
        }
        rows += 1;
    }
    if rows == 0 {
        return (0, f64::NAN);
    }
    for p in power.iter_mut() {
        *p /= rows as f64;
    }
    let mut srt: Vec<f64> = power[2..=64].to_vec();
    srt.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = srt[srt.len() / 2];
    let (mut best, mut bp) = (0usize, 0.0f64);
    for lam in 2..=64usize {
        if power[lam] > bp {
            bp = power[lam];
            best = lam;
        }
    }
    (best, if med > 0.0 { bp / med } else { f64::NAN })
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(512);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/aa".into());
    let zf: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (fcx, fcy, fhalf) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let half = fhalf * zf;
    let (cx, cy) = if zf < 1.0 { (fcx - fhalf * 0.45, fcy + fhalf * 0.42) } else { (fcx, fcy) };

    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);

    // ONE march. `keep_copy_shapes` is the only setting that differs from a normal render, and it
    // changes no trajectory -- it retains a per-copy quantity that was already computed.
    let cfg = EnsembleCfg::production().with_overrides(&[
        Override::TMax(50.0),
        Override::NSync(125),
        Override::RefineFlagged(false),
        Override::MaxSteps(2_000_000),
        Override::KeepCopyShapes(true),
    ]);
    let px: Vec<PixelOut> =
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
    let copies = px.iter().map(|p| p.copy_shapes.len()).max().unwrap_or(0);

    println!("AA CHECK on config_stability, {res}^2, zoom_frac {zf}, t_max = 50.");
    println!("ONE march, two colourings -- {copies} samples per pixel, already computed.");
    println!();
    println!("{:>10} {:>10} {:>10}", "render", "lambda", "promin");

    let mut out: Vec<(&str, Vec<u8>)> = Vec::new();
    for (name, aa) in [("point", false), ("resolved", true)] {
        let mut rgb = Vec::with_capacity(px.len() * 3);
        for p in &px {
            let c = if aa {
                colour::rgb_resolved(p, Scalar::ShapeSpread, &sites, lo, hi)
            } else {
                colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi)
            };
            rgb.extend_from_slice(&c);
        }
        let lum: Vec<f64> = (0..px.len())
            .map(|i| {
                0.2126 * rgb[3 * i] as f64 + 0.7152 * rgb[3 * i + 1] as f64
                    + 0.0722 * rgb[3 * i + 2] as f64
            })
            .collect();
        let (lam, prom) = banding(&lum, res);
        println!("{name:>10} {lam:>10} {prom:>10.2}");
        let _ = adaptive::save_rect(&format!("{dir}/{name}.png"), res, res, &rgb);
        out.push((name, rgb));
    }

    // How many pixels actually changed, and by how much. A render that moves nothing is inert.
    let (a, b) = (&out[0].1, &out[1].1);
    let mut moved = 0usize;
    let mut worst = 0i32;
    for i in 0..res * res {
        let d = (0..3).fold(0i32, |w, k| {
            w.max((a[3 * i + k] as i32 - b[3 * i + k] as i32).abs())
        });
        if d > 0 {
            moved += 1;
        }
        worst = worst.max(d);
    }
    println!();
    println!(
        "pixels changed {moved} of {} ({:.4}), worst channel delta {worst} of 255",
        res * res,
        moved as f64 / (res * res) as f64
    );
    println!();
    println!("**A render that moves nothing is inert and the banding is not aliasing.** If the");
    println!("bands are `t_end` quantised to sync boundaries, every copy in a footprint snaps to");
    println!("the same boundary and averaging them cannot help -- that is the discriminating");
    println!("prediction, and it is why this is a toggle rather than a replacement.");
}
