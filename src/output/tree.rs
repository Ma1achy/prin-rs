//! The tree dump, and the leaf-boundary overlay.
//!
//! **The tree and the decisions are the output** (SCHEDULER_BRIEF §4). The overlay is a diagnostic
//! — the honest check is that it is dense at fractal boundaries and sparse in smooth regions — but
//! **the threshold sweep is the result**, and a threshold chosen because the picture looked right
//! is an arbitrary constant.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::ensemble::pixel::{EnsembleCfg, PixelOut};
use crate::quad::QuadTree;
use crate::scheduler::{SchedCfg, SchedStats};

pub const MAGIC: &[u8; 4] = b"PRNQ";
pub const VERSION: u32 = 1;

/// One record per quad. Same self-describing shape as the pixel dump: a reader never guesses.
pub const FIELDS: &[&str] = &[
    "index", "level", "parent", "sib_index", "iteration",
    "cx", "cy", "half", "cell_width", "is_leaf",
    "spread_mean", "spread_median", "spread_p90",
    "spread_shape_median", "spread_event_median",
    "alpha", "alpha_mean", "alpha_p90", "alpha_sibling_spread",
    "error_ratio_max", "worst_energy_drift", "n_nonfinite", "n_footprints",
    "decision",
];

pub fn record(t: &QuadTree, i: usize) -> [f64; 24] {
    let q = &t.nodes[i];
    let nan = f64::NAN;
    [
        i as f64,
        q.level as f64,
        q.parent.map(|p| p as f64).unwrap_or(-1.0),
        q.sib_index as f64,
        q.iteration as f64,
        q.cx,
        q.cy,
        q.half,
        q.cell_width(t.n),
        if q.is_leaf() { 1.0 } else { 0.0 },
        q.red.spread_mean,
        q.red.spread_median,
        q.red.spread_p90,
        q.red.spread_shape_median,
        q.red.spread_event_median,
        q.alpha.unwrap_or(nan),
        q.alpha_mean.unwrap_or(nan),
        q.alpha_p90.unwrap_or(nan),
        q.alpha_sibling_spread.unwrap_or(nan),
        q.red.error_ratio_max,
        q.red.worst_energy_drift,
        q.red.n_nonfinite as f64,
        q.red.n_footprints as f64,
        q.decision.code() as f64,
    ]
}

#[allow(clippy::too_many_arguments)]
pub fn write<W: Write>(
    w: &mut W,
    tree: &QuadTree,
    cfg: &SchedCfg,
    ens: &EnsembleCfg,
    st: &SchedStats,
    region: &str,
    precision: &str,
) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    let header = format!(
        "region={} body={} n_samples_per_axis={} n_copies={} budget={} bootstrap_levels={}\n\
         tau_display={} alpha_hi={} alpha_lo={} sib_tau={} policy={} order={} agg={} max_level={:?}\n\
         t_max={} eta={} n_sync={} r_coll_frac={} lc_stable={} jitter_scheme={:?} precision={}\n\
         quads_computed={} footprints={} iterations={} budget_exhausted={} wall_seconds={:.3}\n\
         trajectories_per_quad={} sibling_edge_overlap_frac={:.6}\n\
         fields={}\n",
        region, tree.body, tree.n, ens.n_extra + 1, cfg.budget, cfg.bootstrap_levels,
        cfg.tau_display, cfg.alpha_hi, cfg.alpha_lo, cfg.sib_tau,
        cfg.policy.name(), cfg.order.name(), cfg.agg.name(), cfg.max_level,
        ens.t_max, ens.eta, ens.n_sync, ens.r_coll_frac, ens.lc_stable, ens.jitter_scheme,
        precision,
        st.quads_computed, st.footprints, st.iterations, st.budget_exhausted, st.wall_seconds,
        tree.n * tree.n * (ens.n_extra + 1),
        1.0 / tree.n as f64,
        FIELDS.join(","),
    );
    w.write_all(&(header.len() as u32).to_le_bytes())?;
    w.write_all(header.as_bytes())?;
    w.write_all(&(tree.nodes.len() as u64).to_le_bytes())?;
    w.write_all(&(FIELDS.len() as u32).to_le_bytes())?;

    for i in 0..tree.nodes.len() {
        for v in record(tree, i) {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    Ok(())
}

fn save(path: &Path, w: u32, h: u32, data: &[u8]) -> io::Result<()> {
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), w, h);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header()?.write_image_data(data)?;
    Ok(())
}

/// Leaf boundaries drawn over a base image of the same box.
///
/// `base` is a uniform render of the *same* region at `res × res`, so the tree can be checked
/// against the structure it is supposed to be tracking. §4 question 3.
///
/// **Which base matters.** The outcome image is nearly uniform in near-field (97.7% bounded), and
/// the tree does not track outcome labels — it tracks `ensemble_spread`. Overlaying on the outcome
/// image shows the tree refining where nothing appears to be happening, which is correct behaviour
/// read against the wrong picture. The spread base is the direct check.
pub fn overlay(
    stem: &str,
    suffix: &str,
    tree: &QuadTree,
    base: &[PixelOut],
    res: usize,
    base_rgb: impl Fn(&PixelOut) -> [u8; 3],
) -> io::Result<()> {
    let root = &tree.nodes[0];
    let (x0, y0) = (root.cx - root.half, root.cy - root.half);
    let span = 2.0 * root.half;

    let mut img = vec![0u8; res * res * 3];
    for (k, p) in base.iter().enumerate().take(res * res) {
        // The dump is row-major with y increasing upward; PNG rows go downward, so flip.
        let (jx, jy) = (k % res, k / res);
        let row = res - 1 - jy;
        let o = (row * res + jx) * 3;
        img[o..o + 3].copy_from_slice(&base_rgb(p));
    }

    // Dim the base so the boundaries read, then draw every leaf's edges.
    for b in img.iter_mut() {
        *b = (*b as f64 * 0.55) as u8;
    }
    let to_px = |x: f64, y: f64| -> (i64, i64) {
        let fx = (x - x0) / span * res as f64;
        let fy = (y - y0) / span * res as f64;
        (fx.round() as i64, (res as f64 - fy).round() as i64)
    };
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as usize) < res && (y as usize) < res {
            let o = (y as usize * res + x as usize) * 3;
            img[o..o + 3].copy_from_slice(&c);
        }
    };

    for i in tree.leaves() {
        let q = &tree.nodes[i];
        let (ax, ay) = to_px(q.cx - q.half, q.cy - q.half);
        let (bx, by) = to_px(q.cx + q.half, q.cy + q.half);
        // Deeper leaves brighter, so depth is visible without a legend.
        let t = (q.level as f64 / 12.0).min(1.0);
        let c = [
            (90.0 + 165.0 * t) as u8,
            (255.0 - 120.0 * t) as u8,
            (90.0 + 60.0 * (1.0 - t)) as u8,
        ];
        for x in ax.min(bx)..=ax.max(bx) {
            put(x, ay, c);
            put(x, by, c);
        }
        for y in ay.min(by)..=ay.max(by) {
            put(ax, y, c);
            put(bx, y, c);
        }
    }

    save(Path::new(&format!("{stem}_{suffix}.png")), res as u32, res as u32, &img)
}
