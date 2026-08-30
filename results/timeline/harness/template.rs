//! **Bisect harness.** One slice, one set of parameters, rendered at each commit in a walk.
//!
//! Everything that can be pinned is pinned *here*, in a file that is regenerated per commit
//! from one template rather than checked out with the code. Defaults are exactly what changed
//! between these commits, so a harness that read them would measure the defaults and not the
//! code. The `EnsembleCfg` block below is emitted field-by-field by `run.sh`, which greps the
//! checked-out `src/ensemble/pixel.rs` for each name: a field absent at a commit is *listed in
//! the run log as absent*, which is the bisect signal rather than a gap.
//!
//! The slice is built from the reference UI's own ten-slot `z0`, `zoom`, `pan` and `mag` as
//! literals, so no change in `grid.rs` can move it. `Chart::config_stability()` is deliberately
//! NOT called -- it did not exist before `e53223d`, and a window that is nearly right reads as
//! a physics disagreement.
//!
//! The two colour windows are **fixed constants shared by the whole strip**. An auto-ranged ramp
//! per commit would stretch each panel's own p1-p99 to full scale, which on a question about
//! *bleaching* would manufacture or hide the very thing being looked for. Each panel's own
//! auto-range is printed beside it instead.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

// ---- THE SLICE. The user's saved config, as literals. -------------------------------------
const Z_GLSL: [f64; 10] =
    [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
const MAG: f64 = 1.0;
/// The gate. If a run reports equal masses the decode is being overridden and the render is void.
const M_EXPECT: [f64; 3] = [0.327_35, 0.427_63, 0.245_02];

// ---- FIXED colour windows, shared by every commit. ----------------------------------------
const SPREAD_LO: f64 = 6.85e-5;
const SPREAD_HI: f64 = 4.955e-1;
const DRIFT_LO: f64 = 1.0e-8;
const DRIFT_HI: f64 = 4.0e7;

const NAN_RGB: [u8; 3] = [255, 0, 255];

/// Inferno-like ramp on `log10(drift)`, implemented here rather than called from `colour.rs`,
/// because `Scalar::Drift` and `drift_rgb` did not exist before `5cc8dec` and a panel that
/// changes its colouring mid-strip is not a bisect.
fn drift_px(v: f64, nonfinite: bool) -> [u8; 3] {
    if nonfinite || !v.is_finite() {
        return NAN_RGB;
    }
    let l = (DRIFT_LO.log10(), DRIFT_HI.log10());
    let x = ((v.max(1e-300).log10() - l.0) / (l.1 - l.0)).clamp(0.0, 1.0);
    // five-stop inferno
    const STOPS: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.015],
        [0.34, 0.06, 0.43],
        [0.72, 0.21, 0.33],
        [0.98, 0.55, 0.04],
        [0.99, 1.0, 0.64],
    ];
    let t = x * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    let mut o = [0u8; 3];
    for k in 0..3 {
        o[k] = (255.0 * (STOPS[i][k] * (1.0 - f) + STOPS[i + 1][k] * f)).clamp(0.0, 255.0) as u8;
    }
    o
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let tag: String = std::env::args().nth(1).unwrap_or_else(|| "unknown".into());
    let res: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1024);
    let dir: String = std::env::args().nth(3).unwrap_or_else(|| "bisect_out".into());
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z_GLSL[1],
        z_beta: Z_GLSL[0],
        z_q: [Z_GLSL[4], Z_GLSL[5], Z_GLSL[6], Z_GLSL[7]],
        z_mu: [Z_GLSL[8], Z_GLSL[9]],
    };
    let mut q1 = [0.0f64; 8];
    let mut q2 = [0.0f64; 8];
    q1[1] = MAG; // GLSL dimH = 0 = beta = spec index 1
    q2[0] = MAG; // GLSL dimV = 1 = alpha = spec index 0
    let chart = Chart::Latent { z0, q1, q2 };
    let cx = 2.0 * PAN.0 - 1.0 + ZOOM;
    let cy = 2.0 * PAN.1 - 1.0 + ZOOM;
    let half = ZOOM;

    // THE GATE, before anything is integrated.
    let st = grid::decode_state(&chart, 0, cx, cy);
    let dm = (0..3).map(|k| (st.m[k] - M_EXPECT[k]).abs()).fold(0.0, f64::max);
    println!(
        "[{tag}] masses {:.5?}  max|dm| {dm:.2e}  window z_beta [{:+.4}, {:+.4}] z_alpha [{:+.4}, {:+.4}]",
        st.m,
        Z_GLSL[0] + cx - half,
        Z_GLSL[0] + cx + half,
        Z_GLSL[1] + cy - half,
        Z_GLSL[1] + cy + half,
    );
    assert!(dm < 1e-4, "DECODE OVERRIDDEN -- render void");

    let ens = EnsembleCfg {
// @@CFG@@
    };

    let t0 = std::time::Instant::now();
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();
    let secs = t0.elapsed().as_secs_f64();

    let n = px.len() as f64;
    let nonfin = px.iter().filter(|p| p.n_nonfinite > 0).count();
    let simfail = px.iter().filter(|p| p.state == 6).count();
    let hot = px.iter().filter(|p| !(p.energy_drift_max <= 1e-6)).count();
    let esc = px.iter().filter(|p| p.state == 0).count();
    let col = px.iter().filter(|p| p.state == 2).count();
    let bnd = px.iter().filter(|p| p.state == 1).count();

    // Each panel's OWN auto-range, printed but not used, so the strip stays comparable and the
    // window is on record beside it.
    let (alo, ahi) = colour::range(&px, Scalar::ShapeSpread);
    let mut dv: Vec<f64> =
        px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
    let (dlo, dhi) = (q(&mut dv, 0.02), q(&mut dv, 0.98));
    let dmed = q(&mut dv, 0.5);

    let sites = colour::landmarks(&st.m);
    let mut ubuf = Vec::with_capacity(px.len() * 3);
    let mut obuf = Vec::with_capacity(px.len() * 3);
    let mut dbuf = Vec::with_capacity(px.len() * 3);
    let mut dump: Vec<u8> = Vec::with_capacity(px.len() * 24);
    for p in &px {
        ubuf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, SPREAD_LO, SPREAD_HI));
        obuf.extend_from_slice(&png::outcome_rgb(p));
        dbuf.extend_from_slice(&drift_px(p.energy_drift_max, p.n_nonfinite > 0));
        dump.push(p.outcome);
        dump.push(p.state);
        dump.push(p.n_nonfinite);
        dump.push(0);
        dump.extend_from_slice(&(p.energy_drift_max as f32).to_le_bytes());
        for k in 0..3 {
            dump.extend_from_slice(&(p.shape_vec[k] as f32).to_le_bytes());
        }
        dump.extend_from_slice(&(p.t_end as f32).to_le_bytes());
    }
    let stem = format!("{dir}/{tag}");
    let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &ubuf);
    let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);
    let _ = adaptive::save_rect(&format!("{stem}_drift.png"), res, res, &dbuf);
    let _ = std::fs::write(format!("{stem}.bin"), &dump);

    println!(
        "STAT\t{tag}\t{res}\t{secs:.1}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}",
        nonfin,
        simfail,
        hot,
        hot as f64 / n,
        esc as f64 / n,
        col as f64 / n,
        bnd as f64 / n,
        dlo,
        dhi,
        dmed,
        alo,
        ahi,
    );
}
