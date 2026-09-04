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
pub const VERSION: u32 = 2;

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
    // --- v2: the relative mask, the threshold-free gradient, and their contrasts ---
    "n_hot_rel_within", "n_components_rel_within", "largest_component_rel_within",
    "perimeter_ratio_rel_within",
    "n_hot_rel_between", "n_components_rel_between", "largest_component_rel_between",
    "perimeter_ratio_rel_between",
    "grad_rms_within", "grad_rms_between",
    "sig_layout_rel", "sig_grad_rms",
    "contrast_within", "contrast_between", "contrast_layout_rel", "contrast_grad_rms",
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
        r.layout_rel_within.n_hot as f64,
        r.layout_rel_within.n_components as f64,
        r.layout_rel_within.largest_component as f64,
        r.layout_rel_within.perimeter_ratio,
        r.layout_rel_between.n_hot as f64,
        r.layout_rel_between.n_components as f64,
        r.layout_rel_between.largest_component as f64,
        r.layout_rel_between.perimeter_ratio,
        r.grad_rms_within, r.grad_rms_between,
        sig(Criterion::LayoutRel, Agg::Median),
        sig(Criterion::GradRms, Agg::Median),
        c.contrast(k, Criterion::Within, Agg::Median),
        c.contrast(k, Criterion::Between, Agg::Median),
        c.contrast(k, Criterion::LayoutRel, Agg::Median),
        c.contrast(k, Criterion::GradRms, Agg::Median),
    ]
}

pub fn write<W: Write>(w: &mut W, c: &Cache, ens: &crate::ensemble::pixel::EnsembleCfg, tau: f64) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    let header = format!(
        "region={} chart={} body={} cx={:?} cy={:?} half={:?} levels={} n={} res={}\n\
         chart_params={}\n\
         colouring={} metric={} ramp_lo={:?} ramp_hi={:?} tau={:?}\n\
         t_max={} n_sync={} eta={} n_copies={} jitter_frac={} r_coll_frac={} escape_rule={:?} closure_k={} stop_on_escape={} dtau_mode={:?} clamp_final={} \
         jitter_scheme={:?} precision=f64\n\
         quads={} trajectories={}\n\
         note=err_sum is this quad's SUMMED OKLab distance to the reference were it drawn as a \
leaf; it is a constant of the quad because quads are disjoint, which is what makes the replay \
exact. error(tree) = sum(err_sum over leaves) / res^2.\n\
         note=error=0 means MATCHES THIS SAMPLING, not correct: the reference is the \
fully-refined tree at one sample per pixel, and at the screen floor sub-pixel structure is \
sampled arbitrarily.\n\
         fields={}\n",
        c.region, c.chart.name(), c.body, c.cx, c.cy, c.half, c.levels, c.n, c.res, c.chart.params(),
        c.colouring.name(), c.metric.name(), c.ramp.0, c.ramp.1, tau,
        ens.t_max, ens.n_sync, ens.eta, ens.n_extra + 1, ens.jitter_frac, ens.r_coll_frac, ens.escape_rule, ens.closure_k, ens.stop_on_escape, ens.dtau_mode, ens.clamp_final_step,
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

// ---------------------------------------------------------------------------------------------
// The reader.
//
// The module doc has said "a reader never guesses" since the writer was written, and until now
// there was no reader at all: the committed `.qcache` files could be produced and not consumed,
// so `total_substeps` -- the one machine-independent cost column in this project -- sat on disk
// for a fortnight unreadable. That is the gap this closes, and the cost ledger is what needed it.
//
// It does NOT reconstruct a `Cache`: `err_sum`, `rgb`, `payload`, `reference` and
// `reference_spread` are not in the file, so a returned `Cache` would be a `Cache` with holes and
// every consumer would have to know which. What comes back is the table as written -- the header
// verbatim, the field names, and one `f64` row per quad -- and the caller reads columns by name.

/// The `PRQC` table as written: the header verbatim, the field names, one row of `f64` per quad.
///
/// Deliberately **not** a [`Cache`]. Half of `Cache`'s fields are not in the file, and a struct
/// with silent holes in it is how a consumer ends up reading a default as a measurement.
pub struct QuadRows {
    /// The header, verbatim. Everything the run recorded about itself; parse what you need.
    pub header: String,
    pub region: String,
    pub levels: u32,
    /// Footprint grid edge, so a quad holds `n * n` footprints.
    pub n: usize,
    /// `E + 1`, so a quad holds `n * n * n_copies` trajectories.
    pub n_copies: usize,
    pub fields: Vec<String>,
    pub rows: Vec<Vec<f64>>,
}

impl QuadRows {
    /// Column index by name, or `None` if this file's version does not carry it.
    ///
    /// By **name**, never by position: `FIELDS` has been appended to once already (v2 added the
    /// relative mask and the gradient), and a positional read of a v1 file under v2's indices
    /// would return a different column silently.
    pub fn col(&self, name: &str) -> Option<usize> {
        self.fields.iter().position(|f| f == name)
    }

    /// One cell, `NaN` if the column is absent — an absent field is not a zero, and a sum over
    /// `NaN` is loud where a sum over zero is a plausible wrong answer.
    pub fn get(&self, row: usize, name: &str) -> f64 {
        self.col(name).map_or(f64::NAN, |c| self.rows[row][c])
    }

    /// `(level, ix, iy)` — the same [`crate::metric::Key`] the cache is indexed by.
    pub fn key(&self, row: usize) -> crate::metric::Key {
        (
            self.get(row, "level") as u32,
            self.get(row, "ix") as u32,
            self.get(row, "iy") as u32,
        )
    }

    /// Trajectories per quad: `n * n * n_copies`. The memory unit the plan's `mem_all` counts.
    pub fn trajectories_per_quad(&self) -> u64 {
        (self.n * self.n * self.n_copies) as u64
    }
}

pub fn read<R: std::io::Read>(r: &mut R) -> io::Result<QuadRows> {
    let bad = |m: String| io::Error::new(io::ErrorKind::InvalidData, m);

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(bad("not a PRQC file".into()));
    }
    let mut u32b = [0u8; 4];
    r.read_exact(&mut u32b)?;
    let v = u32::from_le_bytes(u32b);
    if v > VERSION {
        return Err(bad(format!("PRQC version {v}, this build reads up to {VERSION}")));
    }
    r.read_exact(&mut u32b)?;
    let hlen = u32::from_le_bytes(u32b) as usize;
    let mut hb = vec![0u8; hlen];
    r.read_exact(&mut hb)?;
    let header = String::from_utf8_lossy(&hb).into_owned();

    let tok = |name: &str| -> Option<String> {
        header
            .split_whitespace()
            .find_map(|t| t.strip_prefix(&format!("{name}=")).map(str::to_string))
    };
    // `fields=` is the last line and its value contains no spaces, but `chart_params` and the
    // notes do, so it is taken line-wise rather than by whitespace token.
    let fields: Vec<String> = header
        .lines()
        .find_map(|l| l.strip_prefix("fields="))
        .ok_or_else(|| bad("PRQC header carries no fields= line".into()))?
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();

    let mut u64b = [0u8; 8];
    r.read_exact(&mut u64b)?;
    let n_rows = u64::from_le_bytes(u64b) as usize;
    r.read_exact(&mut u32b)?;
    let n_fields = u32::from_le_bytes(u32b) as usize;
    // The header's `fields=` and the record's field count are written from the same `FIELDS`, so
    // a disagreement means the file is truncated or from a build whose header and body diverged.
    // Refuse rather than read the shorter of the two, which would shift every column.
    if n_fields != fields.len() {
        return Err(bad(format!(
            "PRQC field count {n_fields} but the header names {} fields",
            fields.len()
        )));
    }

    let mut buf = vec![0u8; n_rows * n_fields * 8];
    r.read_exact(&mut buf)?;
    let rows: Vec<Vec<f64>> = buf
        .chunks_exact(n_fields * 8)
        .map(|c| c.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect())
        .collect();

    // **`region` is not a whitespace token.** The header writes `region={} chart={} ...` on one
    // shared line and three region names carry a space -- `deep interior`, `body2 core`,
    // `mid-field` does not -- so a whitespace read returns `deep` and the header stops being
    // self-describing. Exactly the truncation `fcache`'s `line_field` was written to avoid, at a
    // site that had no reader to notice. Taken as the text between `region=` and the ` chart=`
    // that follows it, which is exact for every file already written.
    let region = header
        .split_once("region=")
        .and_then(|(_, rest)| rest.split_once(" chart=").map(|(r, _)| r.to_string()))
        .or_else(|| tok("region"))
        .unwrap_or_default();

    Ok(QuadRows {
        region,
        levels: tok("levels").and_then(|s| s.parse().ok()).unwrap_or(0),
        n: tok("n").and_then(|s| s.parse().ok()).unwrap_or(0),
        n_copies: tok("n_copies").and_then(|s| s.parse().ok()).unwrap_or(0),
        header,
        fields,
        rows,
    })
}
