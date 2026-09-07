//! **`tilt_plambda` over a grid of zoom against pan, at the shipped gamma.**
//!
//! The gamma sweep said the mismatch against the reference GLSL is not the rotation, so this is
//! the other two window controls. Rows tighten the zoom; columns walk the pan toward the
//! **bottom-right** of the saved panel.
//!
//! **Which way is "bottom right" is a fact about the writer, not a guess.** `Slice::axis` runs
//! low-to-high with index and `save_rect` writes rows in buffer order, so PNG row 0 is the
//! **minimum** `v` — the top of the image is low `v`. So bottom-right is `cx` **and** `cy`
//! increasing, and pushing the view that way takes the structure at the top of the frame off the
//! top. A harness that assumed the opposite would walk the window away from the target and read
//! as "panning does not help"; this project has already paid for one flipped vertical axis.
//!
//! Run: `cargo run --release --example window_sweep <out_root> [res] [t_max]`
//!
//! **A one-row mode**, for walking a chosen cell along a single axis once the grid has picked it:
//! `... <out_root> <res> <t_max> <zoom_mult> <u_lo> <u_hi> <u_step> [v_shift]`. Shifts stay in
//! units of that zoom's own half-window, so a number means the same fraction of a frame here as
//! it does in the grid.
//!
//! The output root is argument one and has **no default**, for the same reason as `gamma_sweep`.
//!
//! # One window over the whole grid, and what that costs here
//!
//! Every panel takes the p1–p99 of all passes pooled, so the comparison is not about the ramp.
//! Note the cost this time: a tighter zoom sees less of the field, so its own quantiles are
//! narrower and it will read flatter under a pooled ramp than it would alone. That is the honest
//! direction — it understates the tight rows rather than flattering them.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelSlim};
use prin_rs::grid::{self, Chart};
use prin_rs::output::adaptive;
use prin_rs::output::colour::{self, Scalar};
use rayon::prelude::*;

const Z10: [f64; 10] = [2.23, -0.56, -0.05, 0.0, 0.05, -0.04, -0.02, 0.12, -0.1, 0.02];
const ZOOM: f64 = 0.167_798_447_231_782_42;
const PAN: (f64, f64) = (0.516_965_993_956_604_7, 0.538_322_670_350_332_3);
const GAMMA: f64 = 2.0;

/// Zoom multipliers, coarsest first, so row 0 is the shipped window.
const ZOOMS: [f64; 4] = [1.0, 0.75, 0.55, 0.40];
/// Pan shifts toward bottom-right, in units of the **row's own** half-window, so a column means
/// the same fraction of a frame at every zoom.
const SHIFTS: [f64; 4] = [0.0, 0.5, 1.0, 1.5];

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let root: String = std::env::args().nth(1).expect("argument one is the output root, no default");
    let res: usize = arg(2, 64);
    let t_max: f64 = arg(3, EnsembleCfg::production().t_max);

    let dir = format!("{root}/window_sweep");
    let _ = std::fs::create_dir_all(&dir);

    let base = EnsembleCfg::production();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { t_max, n_sync, ..base };

    println!("tilt_plambda window sweep at gamma {GAMMA}: {} zooms x {} pans, {res}^2, t={t_max}",
             ZOOMS.len(), SHIFTS.len());
    println!("config: {}", ens.provenance());
    println!("reproduce: cargo run --release --example window_sweep {root} {res} {t_max}");
    println!("shipped window: cx {:.6} cy {:.6} half {:.6}", 2.0 * PAN.0 - 1.0, 2.0 * PAN.1 - 1.0, ZOOM);

    // The one-row mode: a zoom multiplier and a range along +u, with `v` held. Present iff a
    // fourth argument is given, so the default invocation is the 4x4 grid exactly as before.
    let row: Option<(f64, Vec<f64>, f64)> = std::env::args().nth(4).and_then(|s| s.parse().ok()).map(
        |zm: f64| {
            let (ulo, uhi, ustep): (f64, f64, f64) = (arg(5, 0.0), arg(6, 1.0), arg(7, 0.25));
            let vsh: f64 = arg(8, 0.0);
            let mut us = Vec::new();
            let mut u = ulo;
            while u <= uhi + 1e-9 {
                us.push((u * 1e6).round() / 1e6);
                u += ustep;
            }
            (zm, us, vsh)
        },
    );

    let plan: Vec<(f64, f64, f64)> = match &row {
        Some((zm, us, vsh)) => us.iter().map(|&u| (*zm, u, *vsh)).collect(),
        None => ZOOMS.iter().flat_map(|&z| SHIFTS.iter().map(move |&f| (z, f, f))).collect(),
    };
    if let Some((zm, us, vsh)) = &row {
        println!("one-row mode: zoom x{zm}, u shifts {us:?}, v held at {vsh}");
    }

    let mut cells: Vec<(f64, f64, f64, f64, Vec<PixelSlim>)> = Vec::new();
    let t0 = std::time::Instant::now();
    {
        for &(z, f, fv) in &plan {
            let half = ZOOM * z;
            // Bottom-right of the saved panel is +u and +v; see the module note. In the grid
            // the two shifts are equal (a diagonal walk); the one-row mode drives them apart.
            let (cx, cy) = (2.0 * PAN.0 - 1.0 + f * half, 2.0 * PAN.1 - 1.0 + fv * half);
            let (chart, ..) = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], GAMMA, ZOOM, PAN);
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
            let px: Vec<PixelSlim> = (0..sl.npix())
                .into_par_iter()
                .map(|k| PixelSlim::thin(&pixel::evaluate::<f64>(&sl, k, &ens)))
                .collect();
            cells.push((z, f, cx, cy, px));
        }
    }
    println!("{} passes in {:.1}s", cells.len(), t0.elapsed().as_secs_f64());

    let (wlo, whi) = colour::range_q_of(
        cells.iter().flat_map(|(.., px)| px.iter().map(|p| p.spread_shape)),
        0.01,
        0.99,
    );
    println!("pooled ramp window: ({wlo:.4e}, {whi:.4e})");

    let (chart0, cx0, cy0, _) = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], GAMMA, ZOOM, PAN);
    let sites = colour::landmarks(&grid::decode_state(&chart0, 0, cx0, cy0).m);

    let scale = (256 / res).max(1);
    let tile = res * scale;
    let mut tiles: Vec<Vec<u8>> = Vec::with_capacity(cells.len());
    for (z, f, cx, cy, px) in &cells {
        let mut buf = Vec::with_capacity(px.len() * 3);
        for p in px.iter().map(|p| p.fat()) {
            buf.extend_from_slice(&colour::rgb_veto(
                &p, Scalar::ShapeSpread, &sites, wlo, whi, colour::Veto::None,
            ));
        }
        let name = format!("{dir}/z{:03.0}_pan{:03.0}.png", z * 100.0, f * 100.0);
        let _ = adaptive::save_rect(&name, res, res, &buf);
        println!("  zoom x{z:.2} pan +{f:.1} half -> cx {cx:+.6} cy {cy:+.6} half {:.6}  {name}",
                 ZOOM * z);

        let mut up = vec![0u8; tile * tile * 3];
        for y in 0..tile {
            for x in 0..tile {
                let s = ((y / scale) * res + (x / scale)) * 3;
                up[(y * tile + x) * 3..(y * tile + x) * 3 + 3].copy_from_slice(&buf[s..s + 3]);
            }
        }
        tiles.push(up);
    }

    // The anchor is the shipped window: zoom x1, no shift, which is cell 0.
    for y in 0..tile {
        for x in 0..tile {
            if x < 2 || y < 2 || x >= tile - 2 || y >= tile - 2 {
                let o = (y * tile + x) * 3;
                tiles[0][o..o + 3].copy_from_slice(&[255, 255, 255]);
            }
        }
    }

    let (cols, rows, gap) = match &row {
        // Capped at four across: a single row of eight is a 2120x272 strip, which is a shape no
        // one can read on a phone. The order is unchanged and stated with the sheet.
        Some((_, us, _)) => (us.len().min(4), us.len().div_ceil(us.len().min(4)), 8usize),
        None => (SHIFTS.len(), ZOOMS.len(), 8usize),
    };
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
    match &row {
        Some((zm, us, vsh)) => println!(
            "sheet {sw}x{sh}: {cols} across, row-major, zoom x{zm}, u shifts {us:?} of a \
             half-window to the RIGHT, v held at {vsh}; white frame is the first cell"
        ),
        None => println!(
            "sheet {sw}x{sh}: rows are zoom {ZOOMS:?}, columns are pan {SHIFTS:?} of a \
             half-window toward bottom-right; white frame is the shipped window"
        ),
    }
    println!("wrote {spath}");
}
