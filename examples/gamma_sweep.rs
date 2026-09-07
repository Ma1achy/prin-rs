//! **`tilt_plambda` at every gamma from −10 to +10 degrees, one degree apart.**
//!
//! A visual sweep, not a measurement: the port's structure matches the reference GLSL and its
//! orientation does not, and the candidate causes — a sign convention on gamma, a different
//! handedness in the UI's rotation, a transposed basis — all show up as *which angle looks
//! right*. So the instrument is a contact sheet and the reading is by eye.
//!
//! Run: `cargo run --release --example gamma_sweep <out_root> [res] [t_max] [lo] [hi] [step]`
//!
//! **The output root is argument one and has no default.** This writes 21 panels of one slice at
//! a small raster; a run like that landing in `results/` is the failure this project has paid for
//! twice, and a default is what makes it happen.
//!
//! # One window over the whole sweep
//!
//! Every panel takes the p1–p99 of **all 21 passes pooled**. A ramp re-ranged per panel would
//! stretch each rotation's own quantiles to full scale, and the question here is which rotation
//! looks right — the one comparison an auto-ranged ramp is guaranteed to corrupt.
//!
//! # The sheet
//!
//! Tiles are nearest-neighbour upscaled (the panels are 64² and a viewer would smooth them), laid
//! out row-major from `lo` to `hi`. The tile nearest the **shipped** gamma carries a white frame,
//! so the sheet has an anchor without needing glyphs drawn into it.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelSlim};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::adaptive;
use rayon::prelude::*;

const Z10: [f64; 10] = [2.23, -0.56, -0.05, 0.0, 0.05, -0.04, -0.02, 0.12, -0.1, 0.02];
const ZOOM: f64 = 0.167_798_447_231_782_42;
const PAN: (f64, f64) = (0.516_965_993_956_604_7, 0.538_322_670_350_332_3);
const SHIPPED: f64 = 2.0;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let root: String = std::env::args().nth(1).expect(
        "argument one is the output root and has no default: 21 panels of one slice at a small \
         raster must not land in results/",
    );
    let res: usize = arg(2, 64);
    let t_max: f64 = arg(3, EnsembleCfg::production().t_max);
    let (lo, hi, step): (f64, f64, f64) = (arg(4, -10.0), arg(5, 10.0), arg(6, 1.0));

    let dir = format!("{root}/gamma_sweep");
    let _ = std::fs::create_dir_all(&dir);

    let base = EnsembleCfg::production();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { t_max, n_sync, ..base };

    let mut angles: Vec<f64> = Vec::new();
    let mut g = lo;
    while g <= hi + 1e-9 {
        angles.push((g * 1e6).round() / 1e6);
        g += step;
    }
    println!("tilt_plambda gamma sweep: {} angles from {lo} to {hi} by {step}, {res}^2, t={t_max}",
             angles.len());
    println!("config: {}", ens.provenance());
    println!("reproduce: cargo run --release --example gamma_sweep {root} {res} {t_max} {lo} {hi} {step}");

    // Every pass first, then the pooled window, then the paint. The panels cannot be written as
    // they are computed without giving each its own ramp.
    let mut fields: Vec<(f64, Vec<PixelSlim>)> = Vec::with_capacity(angles.len());
    let t0 = std::time::Instant::now();
    for &deg in &angles {
        let (chart, cx, cy, half) =
            Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], deg, ZOOM, PAN);
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let px: Vec<PixelSlim> = (0..sl.npix())
            .into_par_iter()
            .map(|k| PixelSlim::thin(&pixel::evaluate::<f64>(&sl, k, &ens)))
            .collect();
        fields.push((deg, px));
    }
    println!("{} passes in {:.1}s", fields.len(), t0.elapsed().as_secs_f64());

    let (wlo, whi) = colour::range_q_of(
        fields.iter().flat_map(|(_, px)| px.iter().map(|p| p.spread_shape)),
        0.01,
        0.99,
    );
    println!("pooled ramp window over all {} passes: ({wlo:.4e}, {whi:.4e})", fields.len());

    let (chart0, cx0, cy0, _) = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], 0.0, ZOOM, PAN);
    let sites = colour::landmarks(&grid::decode_state(&chart0, 0, cx0, cy0).m);

    // Panels, and the tiles for the sheet in the same pass.
    let scale = (256 / res).max(1);
    let tile = res * scale;
    let mut tiles: Vec<Vec<u8>> = Vec::with_capacity(fields.len());
    for (deg, px) in &fields {
        let mut buf = Vec::with_capacity(px.len() * 3);
        for p in px.iter().map(|p| p.fat()) {
            buf.extend_from_slice(&colour::rgb_veto(
                &p, Scalar::ShapeSpread, &sites, wlo, whi, colour::Veto::None,
            ));
        }
        let name = format!("{dir}/gamma_{}{:04.1}.png", if *deg < 0.0 { "m" } else { "p" }, deg.abs());
        let _ = adaptive::save_rect(&name, res, res, &buf);

        let mut up = vec![0u8; tile * tile * 3];
        for y in 0..tile {
            for x in 0..tile {
                let s = ((y / scale) * res + (x / scale)) * 3;
                up[(y * tile + x) * 3..(y * tile + x) * 3 + 3].copy_from_slice(&buf[s..s + 3]);
            }
        }
        tiles.push(up);
    }

    // The anchor: a white frame on the tile nearest the shipped gamma. Without it the sheet has
    // no landmark and a reader has to count tiles to say which angle they are looking at.
    let anchor = angles
        .iter()
        .enumerate()
        .min_by(|a, b| (a.1 - SHIPPED).abs().partial_cmp(&(b.1 - SHIPPED).abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    for y in 0..tile {
        for x in 0..tile {
            if x < 2 || y < 2 || x >= tile - 2 || y >= tile - 2 {
                let o = (y * tile + x) * 3;
                tiles[anchor][o..o + 3].copy_from_slice(&[255, 255, 255]);
            }
        }
    }

    let cols = 5usize;
    let rows = fields.len().div_ceil(cols);
    let gap = 8usize;
    let (sw, sh) = (cols * tile + (cols + 1) * gap, rows * tile + (rows + 1) * gap);
    let mut sheet = vec![24u8; sw * sh * 3];
    for (i, t) in tiles.iter().enumerate() {
        let (c, r) = (i % cols, i / cols);
        let (ox, oy) = (gap + c * (tile + gap), gap + r * (tile + gap));
        for y in 0..tile {
            let dst = ((oy + y) * sw + ox) * 3;
            sheet[dst..dst + tile * 3].copy_from_slice(&t[y * tile * 3..(y + 1) * tile * 3]);
        }
    }
    let spath = format!("{dir}/contact_sheet.png");
    let _ = adaptive::save_rect(&spath, sw, sh, &sheet);
    println!("sheet {sw}x{sh}, {cols} cols, row-major {lo} -> {hi} by {step}; \
              white frame on gamma = {:+.1} (shipped)", angles[anchor]);
    println!("wrote {spath}");
}
