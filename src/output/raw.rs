//! The raw dump: a header plus one fixed-size record per pixel.
//!
//! Every numeric field is `f64` regardless of kernel precision, so an f32 run and an f64 run
//! produce directly comparable dumps — which is what BRIEF §8 experiment 2 needs.
//!
//! The image is a diagnostic; **this file is the product**.

use std::io::{self, Write};

use crate::ensemble::pixel::{EnsembleCfg, PixelOut};
use crate::grid::Slice;

pub const MAGIC: &[u8; 4] = b"PRIN";
pub const VERSION: u32 = 1;

/// Field order in each record. Written into the header so a reader never has to guess.
pub const FIELDS: &[&str] = &[
    "legacy_class", "binary_id", "n_nonfinite", "censored",
    "t_end", "d_min_ref", "d_min_true", "d_min_gap",
    "energy_drift_nominal", "energy_drift_max", "gamma_max",
    "shape_x", "shape_y", "shape_z",
    "spread_shape", "svar", "spread_event", "ensemble_spread",
    "sigma_e_0", "sigma_e_t", "error_ratio", "error_ratio_mad",
    "switches", "ref_disagree",
];

pub fn write<W: Write>(
    w: &mut W,
    slice: &Slice,
    cfg: &EnsembleCfg,
    precision: &str,
    pixels: &[PixelOut],
) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    // A text header, newline-terminated and length-prefixed. Self-describing beats compact:
    // a dump that outlives the code that wrote it has to carry its own provenance.
    let header = format!(
        "nx={} ny={} cx={} cy={} half={} body={}\n\
         t_max={} n_sync={} eta={} max_steps={} n_copies={} jitter_frac={} seed={}\n\
         ref_policy={:?} lc_stable={} precision={} eps=0\n\
         fields={}\n",
        slice.nx, slice.ny, slice.cx, slice.cy, slice.half, slice.body,
        cfg.t_max, cfg.n_sync, cfg.eta, cfg.max_steps, cfg.n_extra + 1,
        cfg.jitter_frac, cfg.seed,
        cfg.ref_policy, cfg.lc_stable, precision,
        FIELDS.join(","),
    );
    w.write_all(&(header.len() as u32).to_le_bytes())?;
    w.write_all(header.as_bytes())?;
    w.write_all(&(pixels.len() as u64).to_le_bytes())?;
    w.write_all(&(FIELDS.len() as u32).to_le_bytes())?;

    for p in pixels {
        for v in record(p) {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    Ok(())
}

pub fn record(p: &PixelOut) -> [f64; 24] {
    [
        p.legacy_class as f64,
        p.binary_id as f64,
        p.n_nonfinite as f64,
        if p.censored { 1.0 } else { 0.0 },
        p.t_end,
        p.d_min_ref,
        p.d_min_true,
        p.d_min_gap,
        p.energy_drift_nominal,
        p.energy_drift_max,
        p.gamma_max,
        p.shape_vec[0],
        p.shape_vec[1],
        p.shape_vec[2],
        p.spread_shape,
        p.svar,
        p.spread_event,
        p.ensemble_spread,
        p.sigma_e_0,
        p.sigma_e_t,
        p.error_ratio,
        p.error_ratio_mad,
        p.switches as f64,
        p.ref_disagree as f64,
    ]
}
