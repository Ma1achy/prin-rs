//! `PRQC` — the complete-tree cache, dumped so every §2 curve can be recomputed offline.
//!
//! The metric integrates one complete tree per region and then replays every criterion, both
//! controls and the whole `error(B)` curve over it. **Without this dump that tree lives only in
//! RAM for the length of one process**, and reproducing any table in §10 means paying the
//! 2.8-million-trajectory integration again. Every criterion is dumped whatever the run was
//! ranking on, which is what makes offline comparison real rather than aspirational.
//!
//! Same self-describing shape as `PRIN` and `PRNQ`: magic, version, a length-prefixed text
//! header naming every parameter, then a record count, a field count, and one `f64` per field
//! per quad. A reader never guesses.
//!
//! `err_sum` is the field that makes the replay possible: quads are disjoint, so a quad's
//! contribution to the image error is a **constant** independent of what the rest of the tree
//! does. That is why the greedy replay is a static priority queue and why the whole curve is a
//! traversal rather than a re-render.

use std::io::{self, Write};

use crate::metric::Cache;
use crate::quad::{Agg, Criterion};

pub const MAGIC: &[u8; 4] = b"PRQC";
pub const VERSION: u32 = 1;

pub const FIELDS: &[&str] = &[
    "level", "ix", "iy", "cx", "cy", "half",
    "err_sum", "gain",
    "spread_mean", "spread_median", "spread_p90",
    "between_shape", "between_event", "between_spread", "between_matched", "within_pooled",
    "n_hot_within", "n_components_within", "largest_component_within", "perimeter_ratio_within",
    "n_hot_between", "n_components_between", "largest_component_between", "perimeter_ratio_between",
    "frac_above_tau_within", "frac_above_tau_between",
    "terminated_fraction", "escape_fraction", "t_end_gradient",
    "running_max_divergence", "divergence_trend", "frac_diverged", "first_divergence_median",
    "error_ratio_max", "worst_energy_drift", "total_substeps", "n_distinct_ic", "n_nonfinite",
    // Every criterion's scalar, so a reader can rank offline without reimplementing `signal`.
    "sig_within_median", "sig_within_mean", "sig_within_p90",
    "sig_between", "sig_max_of_both",
    "sig_frac_hot_within", "sig_frac_hot_between", "sig_layout",
    "sig_running_max", "sig_first_div", "sig_term_grad",
    "contrast_within", "contrast_between",
];

fn record(c: &Cache, k: crate::metric::Key) -> Vec<f64> {
    let q = c.get(k);
    let r = &q.red;
    let (l, ix, iy) = k;
    let h = c.half / (1u64 << l) as f64;
    let cx = c.cx - c.half + (2 * ix + 1) as f64 * h;
    let cy = c.cy - c.half + (2 * iy + 1) as f64 * h;
    let sig = |cr: Criterion, a: Agg| r.signal(cr, a);
    vec![
        l as f64, ix as f64, iy as f64, cx, cy, h,
        q.err_sum, c.gain(k),
        r.spread_mean, r.spread_median, r.spread_p90,
        r.between_shape, r.between_event, r.between_spread, r.between_matched, r.within_pooled,
        r.layout_within.n_hot as f64,
        r.layout_within.n_components as f64,
        r.layout_within.largest_component as f64,
        r.layout_within.perimeter_ratio,
        r.layout_between.n_hot as f64,
        r.layout_between.n_components as f64,
        r.layout_between.largest_component as f64,
        r.layout_between.perimeter_ratio,
        r.frac_above_tau_within, r.frac_above_tau_between,
        r.terminated_fraction, r.escape_fraction, r.t_end_gradient,
        r.running_max_divergence_median, r.divergence_trend_median,
        r.frac_diverged, r.first_divergence_median,
        r.error_ratio_max, r.worst_energy_drift,
        r.total_substeps as f64, r.n_distinct_ic as f64, r.n_nonfinite as f64,
        sig(Criterion::Within, Agg::Median),
        sig(Criterion::Within, Agg::Mean),
        sig(Criterion::Within, Agg::P90),
        sig(Criterion::Between, Agg::Median),
        sig(Criterion::MaxOfBoth, Agg::Median),
        sig(Criterion::FracHotWithin, Agg::Median),
        sig(Criterion::FracHotBetween, Agg::Median),
        sig(Criterion::Layout, Agg::Median),
        sig(Criterion::RunningMax, Agg::Median),
        sig(Criterion::FirstDivergence, Agg::Median),
        sig(Criterion::TerminationGradient, Agg::Median),
        c.contrast(k, Criterion::Within, Agg::Median),
        c.contrast(k, Criterion::Between, Agg::Median),
    ]
}

pub fn write<W: Write>(w: &mut W, c: &Cache, ens: &crate::ensemble::pixel::EnsembleCfg, tau: f64) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    let header = format!(
        "region={} body={} cx={:?} cy={:?} half={:?} levels={} n={} res={}\n\
         colouring={} ramp_lo={:?} ramp_hi={:?} tau={:?}\n\
         t_max={} n_sync={} eta={} n_copies={} jitter_frac={} r_coll_frac={} \
         jitter_scheme={:?} precision=f64\n\
         quads={} trajectories={}\n\
         note=err_sum is this quad's SUMMED OKLab distance to the reference were it drawn as a \
leaf; it is a constant of the quad because quads are disjoint, which is what makes the replay \
exact. error(tree) = sum(err_sum over leaves) / res^2.\n\
         note=error=0 means MATCHES THIS SAMPLING, not correct: the reference is the \
fully-refined tree at one sample per pixel, and at the screen floor sub-pixel structure is \
sampled arbitrarily.\n\
         fields={}\n",
        c.region, c.body, c.cx, c.cy, c.half, c.levels, c.n, c.res,
        c.colouring.name(), c.ramp.0, c.ramp.1, tau,
        ens.t_max, ens.n_sync, ens.eta, ens.n_extra + 1, ens.jitter_frac, ens.r_coll_frac,
        ens.jitter_scheme,
        c.quads.len(), c.trajectories,
        FIELDS.join(","),
    );
    let hb = header.as_bytes();
    w.write_all(&(hb.len() as u32).to_le_bytes())?;
    w.write_all(hb)?;

    // Sorted by (level, iy, ix) so the dump is stable across runs and diffable.
    let mut keys: Vec<crate::metric::Key> = c.quads.keys().cloned().collect();
    keys.sort_by_key(|&(l, ix, iy)| (l, iy, ix));

    w.write_all(&(keys.len() as u64).to_le_bytes())?;
    w.write_all(&(FIELDS.len() as u32).to_le_bytes())?;
    for k in keys {
        for v in record(c, k) {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    Ok(())
}
