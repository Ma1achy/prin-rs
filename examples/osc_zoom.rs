//! **Is the ribbon oscillation chart-locked (physics) or pixel-locked (sampling)?**
//!
//! One centre, three magnifications, **one shared colour window**. The test needs no statistic:
//!
//! - **physics** — a fixed period in chart coordinates, so 2x magnification makes the bands
//!   **2x wider in pixels**;
//! - **pixel-locked** (a beat against the fixed Halton offsets, or against the pixel grid) —
//!   the period stays the same **in pixels** and the bands do not move;
//! - **step-count level sets** — the spacing tracks `eta`, tested by the second row of arms.
//!
//! The window is taken **once** at the coarsest zoom and shared. Auto-ranging per panel is what
//! manufactures or hides the difference a comparison exists to show, and doing it here would
//! rescale each magnification independently and make the bands incomparable by construction.
//!
//! **But a shared window STARVES the deepest zoom, and a starved panel is a dead arm.** Zooming in
//! sees less of the field's variation, so at `z4` the shared window leaves a range of 0.056 --
//! about 14 of 255 levels -- and the panel's content is quantisation noise. Its registration
//! against `z2` then reads 0.0256 against a shifted control of 0.0693: a NULL, and one that would
//! be misread as "the structure is not chart-locked at this depth" when what it means is *this
//! side no longer resolves anything*.
//!
//! Argument five selects the window. `shared` is the default and is the arm to quote for
//! AMPLITUDE. `own` gives each panel its own p1-p99 and is the arm to quote for PERIOD and
//! REGISTRATION -- both are invariant under a monotone rescale, and `range_norm` is affine, so
//! correlation is preserved exactly up to clipping and the 8-bit floor. The two arms answer
//! different questions and the header names which one ran.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, colour as c2};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/osc".into());
    let fx: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.10);
    let fy: f64 = std::env::args().nth(4).and_then(|s| s.parse().ok()).unwrap_or(0.45);
    let own = std::env::args().nth(5).map(|s| s == "own").unwrap_or(false);
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0], z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (fcx, fcy, fh) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (cx, cy) = (fcx + fh * (2.0 * fx - 1.0), fcy + fh * (2.0 * fy - 1.0));
    let m = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = c2::landmarks(&m);

    println!("OSCILLATION ZOOM at panel fraction ({fx}, {fy}), {res}^2, {}.",
             if own { "PER-PANEL WINDOW (period/registration arm)" } else { "ONE SHARED WINDOW (amplitude arm)" });
    println!();
    println!("{:>8} {:>9} {:>12} {:>12} {:>10} {:>9}", "half", "eta", "steps p50", "t range", "turns/row", "px/half");

    let mut window: Option<(f64, f64)> = None;
    for (tag, zf, eta) in [
        ("z1", 0.060f64, 1e-2f64),
        ("z2", 0.030, 1e-2),
        ("z4", 0.015, 1e-2),
        ("z1_eta4", 0.060, 2.5e-3),
        ("z2_eta4", 0.030, 2.5e-3),
    ] {
        let half = fh * zf;
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let cfg = EnsembleCfg::production().with_overrides(&[
            Override::TMax(50.0), Override::NSync(125), Override::Eta(eta),
            Override::RefineFlagged(false), Override::MaxSteps(4_000_000),
        ]);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        // Shared: the coarsest zoom sets it and every other panel reuses it.
        let (lo, hi) = if own {
            colour::range(&px, Scalar::ShapeSpread)
        } else {
            *window.get_or_insert_with(|| colour::range(&px, Scalar::ShapeSpread))
        };

        let mut rgb = Vec::with_capacity(px.len() * 3);
        for p in &px {
            rgb.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
        }
        let stem = if own { format!("{tag}_own") } else { tag.to_string() };
        let _ = adaptive::save_rect(&format!("{dir}/{stem}.png"), res, res, &rgb);

        // Turning points per row of the normalised field, averaged. A monotone gradient gives 0.
        let (mut turns, mut rows) = (0usize, 0usize);
        for y in 0..res {
            let t: Vec<f64> = (0..res)
                .map(|x| colour::range_norm(Scalar::ShapeSpread, px[y * res + x].spread_shape, lo, hi)
                    .unwrap_or(f64::NAN))
                .collect();
            if t.iter().any(|v| !v.is_finite()) {
                continue;
            }
            let d: Vec<f64> = t.windows(2).map(|w| w[1] - w[0]).collect();
            turns += d.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
            rows += 1;
        }
        let tpr = turns as f64 / rows.max(1) as f64;
        let mut st: Vec<u64> = px.iter().map(|p| p.total_substeps).collect();
        st.sort_unstable();
        let mut fin: Vec<f64> = px.iter()
            .map(|p| colour::range_norm(Scalar::ShapeSpread, p.spread_shape, lo, hi).unwrap_or(f64::NAN))
            .filter(|x| x.is_finite()).collect();
        fin.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let rng = if fin.is_empty() { f64::NAN } else { fin[fin.len() - 1] - fin[0] };
        println!("{:>8} {eta:>9.2e} {:>12.3e} {rng:>12.4} {tpr:>10.1} {:>9.1}",
                 format!("{zf:.3}"), st[st.len() / 2] as f64,
                 if tpr > 0.0 { res as f64 / tpr } else { f64::NAN });
    }
    println!();
    println!("`px/half` is pixels per half-oscillation. z1 -> z2 -> z4 is 2x magnification each");
    println!("step, so PHYSICS DOUBLES it each row and PIXEL-LOCKED leaves it flat. The `_eta4`");
    println!("rows repeat z1 and z2 at a quarter step: step-count level sets would move there.");
    println!();
    println!("**Read `steps p50` first.** If it differs by more than ~2x between z1 and z2 the");
    println!("two panels are not the same structure magnified and nothing below it is comparable.");
}
