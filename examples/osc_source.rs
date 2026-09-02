//! **Is the ripple in the FIELD, or in the ESTIMATOR that samples it?**
//!
//! `osc_zoom` settled two of three candidates: the ripples magnify with the chart (high-passed
//! registration 0.66-0.89 against shifted controls of -0.08 to -0.25) and they do not move when
//! the step size falls 3.07x (0.96-0.99, on an arm proven not inert). So they are neither a
//! pixel-grid beat nor integer step-count level sets.
//!
//! **It cannot separate the remaining one.** `spread_shape` is a spread over `E+1` copies jittered
//! within the cell, and `halton_offset` scales with the cell width -- so the sample pattern is
//! *self-similar under zoom*. A beat between that fixed 8-point pattern and the field's own
//! structure would magnify exactly like real structure and register at 0.9 all the same.
//!
//! Two arms separate them, and they fail in opposite directions:
//!
//! - **Change the offsets, keep the field.** `Scheme::Pcg` is a different set entirely; `E+1` of
//!   4 and 16 are different counts *and* different sets. Real structure survives all three. A
//!   sampling beat cannot -- its period is set by the pattern that just changed.
//! - **Remove the ensemble.** `shape_vec`, `t_end`, `d_min_true` and `energy_drift_nominal` are
//!   the NOMINAL copy: one trajectory, no jitter, no spread. A ripple there cannot be a sampling
//!   artefact of any kind, because nothing was sampled.
//!
//! A third arm varies the sampling EXTENT at fixed zoom. `jitter_frac` sets how far the copies
//! spread inside the cell, so 0.25 and 1.0 sample the same field at a different scale without
//! moving the window. Real structure keeps its period and only changes amplitude; a beat's period
//! is set by the sampling geometry and must move with it. This is orthogonal to changing the
//! offset SET -- one varies which points, the other how far apart they are.
//!
//! `jitter_frac = 0` is the guard that keeps the second arm from being vacuous: it collapses every
//! copy onto the cell centre, so `spread_shape` must read identically zero. If it does not, the
//! field is not the ensemble quantity this whole argument assumes it is.
//!
//! Fields are written as raw `f64` so the spectrum is measured on the float and the 8-bit
//! question never arises -- it has already cost this project one wrong conclusion.
use rayon::prelude::*;
use std::io::Write;

use prin_rs::ensemble::jitter::Scheme;
use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn dump(path: &str, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&x.to_le_bytes()).unwrap();
    }
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/osc/fields".into());
    let (fx, fy) = (0.10f64, 0.45f64);
    // The zoom is an ARGUMENT because the subject is not uniform across it. The cross-hatch's
    // high-frequency share of the lightness channel runs 0.0364 / 0.0424 / 0.2908 at zf =
    // 0.060 / 0.030 / 0.015, so a run pinned to the coarsest window measures where the effect is
    // weakest by a factor of eight. Same defect as `pan_sequence`'s hardcoded viewport.
    let zf: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.060);
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
    let half = fh * zf;

    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let base = |ov: &[Override]| {
        let mut all = vec![
            Override::TMax(50.0), Override::NSync(125), Override::Eta(1e-2),
            Override::RefineFlagged(false), Override::MaxSteps(4_000_000),
        ];
        all.extend_from_slice(ov);
        EnsembleCfg::production().with_overrides(&all)
    };

    println!("OSCILLATION SOURCE at zf = {zf}, {res}^2, fields dumped as raw f64.");
    println!();
    println!("{:>14} {:>10} {:>14} {:>14}", "arm", "E+1", "spread p50", "spread sd");

    for (tag, ov) in [
        ("halton_E8", vec![]),
        ("pcg_E8", vec![Override::JitterScheme(Scheme::Pcg)]),
        ("halton_E4", vec![Override::NExtra(3)]),
        ("halton_E16", vec![Override::NExtra(15)]),
        ("jitter_q", vec![Override::JitterFrac(0.25)]),
        ("jitter_f", vec![Override::JitterFrac(1.0)]),
        ("jitter0", vec![Override::JitterFrac(0.0)]),
    ] {
        let cfg = base(&ov);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let s: Vec<f64> = px.iter().map(|p| p.spread_shape).collect();
        dump(&format!("{dir}/{tag}_spread.f64"), &s);
        if tag == "halton_E8" {
            // The nominal copy alone. No jitter enters any of these.
            dump(&format!("{dir}/nominal_shape0.f64"), &px.iter().map(|p| p.shape_vec[0]).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_shape1.f64"), &px.iter().map(|p| p.shape_vec[1]).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_shape2.f64"), &px.iter().map(|p| p.shape_vec[2]).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_tend.f64"), &px.iter().map(|p| p.t_end).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_dmin.f64"), &px.iter().map(|p| p.d_min_true).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_drift.f64"), &px.iter().map(|p| p.energy_drift_nominal).collect::<Vec<_>>());
            dump(&format!("{dir}/nominal_steps.f64"), &px.iter().map(|p| p.total_substeps as f64).collect::<Vec<_>>());
        }
        let mut f: Vec<f64> = s.iter().copied().filter(|x| x.is_finite()).collect();
        f.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = f.len().max(1);
        let mean = f.iter().sum::<f64>() / n as f64;
        let sd = (f.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
        let ecount = cfg.n_extra + 1;
        println!("{tag:>14} {ecount:>10} {:>14.6e} {sd:>14.6e}",
                 if f.is_empty() { f64::NAN } else { f[n / 2] });
    }
    println!();
    println!("`jitter0` must read spread p50 and sd EXACTLY 0 -- every copy is the cell centre.");
    println!("If it does not, `spread_shape` is not the ensemble quantity the argument assumes.");
}
