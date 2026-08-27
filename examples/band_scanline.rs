//! Which field is the banding actually in? One column, every quantity, no guessing.
//!
//! # Why this exists
//!
//! §21 established that `t_end` is quantised to the sync cadence wherever escape terminates, and
//! that decoupling the escape test removes it: `preset_plambda` at 1024^2 goes from **2497
//! distinct `t_end` values to 692456**, with the fraction landing exactly on a sync boundary
//! falling **99.74% -> 0.19%**.
//!
//! **The rendered image did not change by one bit.** `md5` of the `_uniform` panel is identical
//! at `escape_every` 0 and 4, on all four latent charts. So a 277x change in `t_end` resolution
//! moves the picture not at all, and the concentric arcs have a different cause entirely.
//!
//! Rather than guess which quantity carries them, this walks one image column and prints every
//! field the colouring could depend on, plus the rendered lightness itself. The banded quantity
//! is the one whose distinct-value count is small and whose steps line up with the arcs.
//!
//! # The controls that make it readable
//!
//! - **`lightness` is recomputed from the same ramp the render used**, so a step in the printed
//!   column is a step in the image and not in a proxy for it.
//! - **`t_end` stays in the table** even though it is now known not to be the cause. Dropping it
//!   would remove the negative control, and the whole point is that a field can be quantised and
//!   irrelevant at the same time.
//! - Distinct-value counts per column, because a staircase is a small count and an eye is not a
//!   measurement.
//!
//! # Writes
//!
//! stdout only.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::colour::{self, Scalar};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn distinct(v: &[f64]) -> usize {
    let mut b: Vec<u64> = v.iter().map(|x| x.to_bits()).collect();
    b.sort_unstable();
    b.dedup();
    b.len()
}

fn main() {
    let res: usize = arg(1, 1024);
    let col: usize = arg(2, 200);
    let ev: usize = arg(3, 0);
    let case: String = std::env::args().nth(4).unwrap_or_else(|| "preset_plambda".into());
    let rows: usize = arg(5, 60);

    let (name, chart, cx, cy, half) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == case)
        .expect("case is in the gallery");

    let ens = EnsembleCfg { refine_flagged: false, escape_every: ev, ..Default::default() };
    let dt_sync = ens.t_max / ens.n_sync as f64;
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);

    // The whole slice, because the ramp is a slice-wide p1-p99 and a column cannot set it.
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();
    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);

    // The column, in image order. `Slice` is row-major over `res x res`.
    let idx: Vec<usize> = (0..res).map(|y| y * res + col).collect();
    let take = |f: &dyn Fn(&PixelOut) -> f64| -> Vec<f64> {
        idx.iter().map(|k| f(&px[*k])).collect()
    };

    let t_end = take(&|p| p.t_end);
    let sshape = take(&|p| p.spread_shape);
    let sev = take(&|p| p.spread_event);
    let sv0 = take(&|p| p.shape_vec[0]);
    let sv1 = take(&|p| p.shape_vec[1]);
    let sv2 = take(&|p| p.shape_vec[2]);
    // **The actual rendered bytes, from `colour::rgb` itself.** The first cut of this
    // reconstructed the lightness from the ramp by hand, which puts a second implementation of
    // the colouring between the measurement and the thing measured -- and the whole question is
    // whether the steps are in the field or in the output. These ARE the bytes in the PNG.
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);
    let rgb: Vec<[u8; 3]> =
        idx.iter().map(|k| colour::rgb(&px[*k], Scalar::ShapeSpread, &sites, lo, hi)).collect();
    // A luma summary, only so there is one number per pixel in the printed column. Edges are
    // detected on the full triple, so this is a display column and not the instrument.
    let _light: Vec<f64> = rgb
        .iter()
        .map(|c| (0.2126 * c[0] as f64 + 0.7152 * c[1] as f64 + 0.0722 * c[2] as f64) / 255.0)
        .collect();

    println!(
        "{name} at {res}^2, column x={col}, escape_every={ev}, n_sync={} (interval \
         {dt_sync:.6}), ramp ({lo:.4e}, {hi:.4e})\n",
        ens.n_sync
    );
    println!(
        "distinct over the column ({} px): t_end {}  spread_shape {}  spread_event {}  \
         shape_vec[0] {}  RENDERED RGB {}",
        res,
        distinct(&t_end),
        distinct(&sshape),
        distinct(&sev),
        distinct(&sv0),
        {
            let mut b = rgb.clone();
            b.sort_unstable();
            b.dedup();
            b.len()
        }
    );

    // Step positions in the rendered lightness: where the quantised 8-bit value changes. That is
    // literally what an arc edge is, so it is detected rather than eyeballed.
    // On the FULL triple: an arc edge is any change in the emitted colour, not only a change in
    // a luma summary that could hide two channels moving against each other.
    let edges: Vec<usize> = (1..res).filter(|&y| rgb[y] != rgb[y - 1]).collect();
    println!(
        "\n8-bit lightness edges down this column: {} of {} rows\n  first 24 at y = {:?}",
        edges.len(),
        res,
        &edges[..edges.len().min(24)]
    );
    if edges.len() > 2 {
        let gaps: Vec<usize> = edges.windows(2).map(|w| w[1] - w[0]).collect();
        let mut g = gaps.clone();
        g.sort_unstable();
        println!(
            "  edge spacing: min {} median {} max {}",
            g[0],
            g[g.len() / 2],
            g[g.len() - 1]
        );
    }

    println!(
        "\n{:>5} {:>12} {:>6} {:>13} {:>13} {:>13} {:>10} {:>5}",
        "y", "t_end", "on_b", "spread_shape", "spread_event", "shape_vec[0]", "lightness", "8bit"
    );
    let step = (res / rows).max(1);
    for y in (0..res).step_by(step) {
        let k = (t_end[y] / dt_sync).round();
        let on_b = (t_end[y] - k * dt_sync).abs() <= 1e-9 * ens.t_max;
        println!(
            "{y:>5} {:>12.6} {:>6} {:>13.6e} {:>13.6e} {:>13.6e}   {:>3} {:>3} {:>3}",
            t_end[y],
            if on_b { "yes" } else { "no" },
            sshape[y],
            sev[y],
            sv0[y],
            rgb[y][0],
            rgb[y][1],
            rgb[y][2]
        );
    }
    let _ = (sv1, sv2);

    println!(
        "\nHOW TO READ THIS\n\n\
         The banded quantity is the one with a SMALL distinct count whose steps line up with the\n\
         8-bit lightness edges. `t_end` is kept in the table as the NEGATIVE CONTROL: it is\n\
         genuinely quantised and the render is bitwise identical with and without that\n\
         quantisation, so a field can be stepped and irrelevant at the same time.\n\n\
         If `spread_shape` itself has a small distinct count, the arcs are in the ensemble\n\
         statistic. If it is smooth and the 8-BIT column steps anyway, the arcs are the display\n\
         quantisation of a smooth field -- contours of a gentle gradient at 1/255 intervals,\n\
         which is a rendering property and not a physics one."
    );
}
