//! Confirm the banding visually: one chart, two escape-test cadences, side by side.
//!
//! §21 established the *mechanism* -- `t_end` is quantised to `n_sync` values wherever escape
//! terminates, and decoupling the escape test from the sync boundary takes `preset_plambda` from
//! **16 distinct `t_end` values to 2623** with the fraction landing exactly on a boundary falling
//! **99.52% -> 0.26%**. It did **not** establish that the arcs are gone from the image, and a
//! 56x finer step is not the same claim as a clean render. This is the confirmation.
//!
//! # Why the lightness field carries a `t_end` artefact at all
//!
//! The committed colouring is hue from the shape sphere, lightness from `spread_shape` -- neither
//! of which is `t_end`. The coupling is termination: a copy's shape vector is read **where that
//! copy ended**, so quantising `t_end` quantises where each copy is sampled, and the cross-copy
//! spread inherits it. That is why the arcs appear in a field that never mentions `t_end`.
//!
//! # What each panel is for
//!
//! `_uniform` is the continuous field and is where the arcs live. `_uniform_outcome` is the
//! categorical control: the arcs must be **absent** from it at both cadences (they are a
//! continuous-field artefact) while the crisp polygonal edges must **survive** (they are real
//! regime boundaries). Rendering only the first would leave a change in the second unnoticed.
//!
//! One sample per pixel, no interpolation, and the ramp window is printed -- a false-colour
//! image without its scale is decoration.
//!
//! # Writes
//!
//! `<root>/postfix/`. **Never** `results/charts`, which holds the committed pre-fix renders: they
//! are the "before" for this comparison and for every claim the fix will make. A reduced-argument
//! run has already destroyed committed 1024^2 artefacts once on this project, and the same run
//! would destroy the only copy of the evidence here.
//!
//! Both cadences are re-rendered rather than only the new one, so the pair differs in **exactly
//! one thing**. Reusing the committed image as the "before" would also be reusing whatever else
//! has moved since it was made.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::{adaptive, colour, png};
use prin_rs::output::colour::Scalar;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 1024);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    // Two cadences, and stride 4 rather than 1 because §21.7 measured the labels converged
    // there -- `deep interior` reads 0.1564 at both 4 and 1, so 1 buys nothing and costs a
    // `to_cartesian` per RK4 step.
    //
    // **Stride 4 is not the shipped default.** `EnsembleCfg::escape_every` is still 0; making it
    // nonzero changes every future render against the committed corpus, which is a decision to
    // take deliberately rather than a side effect of confirming a fix. The panels are named for
    // the stride so nothing here reads as "the default".
    let strides: Vec<usize> = vec![0, 4];

    let dir = format!("{root}/postfix");
    let _ = std::fs::create_dir_all(&dir);

    let base = EnsembleCfg::default();
    let dt_sync = base.t_max / base.n_sync as f64;

    println!(
        "{res}^2, one sample per pixel, E+1={}, t={}, n_sync={} (interval {dt_sync:.6})\n\
         strides {strides:?}; escape_confirm ON at both (inert at 0)\n\
         hue = shape sphere, lightness = spread_shape over each panel's own p1-p99\n",
        base.n_extra + 1,
        base.t_max,
        base.n_sync
    );

    // In priority order: the worst case first, then its family, then the one carrying BOTH
    // artefacts, then its window twin, then the region whose terminal class actually moved.
    const WANT: [&str; 4] =
        ["preset_plambda", "preset_prho", "preset_shape_pl_h1", "preset_shape_pl"];
    let mut targets: Vec<(String, grid::Chart, f64, f64, f64, usize)> = Vec::new();
    for w in WANT {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            targets.push((c.0.into(), c.1, c.2, c.3, c.4, 0));
        }
    }
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if region == "deep interior" {
            targets.push((region.replace(' ', "_"), grid::Chart::BodyPlane, cx, cy, 0.05, body));
        }
    }

    for (name, chart, cx, cy, half, body) in targets {
        let name = name.as_str();
        if only != "all" && only != name {
            continue;
        }
        for &ev in &strides {
            let ens = EnsembleCfg { refine_flagged: false, escape_every: ev, ..Default::default() };
            let t0 = std::time::Instant::now();
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            // The number beside the image, so the panel is never read on its own.
            let te: Vec<f64> = px.iter().map(|p| p.t_end).collect();
            let (distinct, on_b) = {
                let mut b: Vec<u64> = te.iter().map(|x| x.to_bits()).collect();
                b.sort_unstable();
                b.dedup();
                let ob = te
                    .iter()
                    .filter(|x| x.is_finite())
                    .filter(|x| {
                        let k = (*x / dt_sync).round();
                        (*x - k * dt_sync).abs() <= 1e-9 * base.t_max
                    })
                    .count();
                (b.len(), ob)
            };
            let esc = px.iter().filter(|p| p.state == 0).count();

            let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
            let m_here = grid::decode_state(&chart, body, cx, cy).m;
            let sites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(px.len() * 3);
            let mut obuf = Vec::with_capacity(px.len() * 3);
            for p in &px {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
                obuf.extend_from_slice(&png::outcome_rgb(p));
            }
            let stem = format!("{dir}/{name}_e{ev}");
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
            let _ = adaptive::save_rect(&format!("{stem}_uniform_outcome.png"), res, res, &obuf);

            println!(
                "{name:>20} escape_every={ev:<3} {secs:>7.1}s   escape {:.4}   \
                 t_end distinct {distinct:>6}  on a sync boundary {:>6.2}%   \
                 ramp ({lo:.3e}, {hi:.3e})",
                esc as f64 / px.len() as f64,
                on_b as f64 / px.len() as f64 * 100.0
            );
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         The `_uniform` pair is the test: the arcs must be present at stride 0 and absent at\n\
         stride 4. The `_uniform_outcome` pair is the CONTROL: the arcs must be absent from BOTH\n\
         (they are a continuous-field artefact) while the crisp polygonal edges survive in both\n\
         (they are outcome-class boundaries, so real regime structure). A change in the control\n\
         would mean the cadence moved the physics rather than the sampling.\n\n\
         `t_end distinct` and `on a sync boundary` are printed beside each panel because the\n\
         image is the confirmation and the counts are the measurement -- neither substitutes for\n\
         the other, and §21 was careful not to claim the first from the second."
    );
}
