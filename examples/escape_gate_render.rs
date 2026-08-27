//! Measurement 3: render before and after the escape **distance gate**, one variable changed.
//!
//! The prediction on record is explicit -- the domes, wedges and pale seam-outlines soften or
//! vanish and the ribbons become continuous. It is a prediction about an *image*, so an image is
//! what settles it, and a finding read off a wireframe or an image is a finding about an
//! appearance until the mechanism is tested. The counts printed beside each panel are the
//! measurement; neither substitutes for the other.
//!
//! # The chain being tested
//!
//! Under `stop_on_event` a terminated trajectory's `shape_vec` is the state at **its own**
//! `t_end`, and termination happens at a sync boundary -- so the rendered field is a patchwork of
//! `n_sync` time strata stitched at hard seams. A ribbon is a level set at a *common* time; where
//! neighbouring pixels froze at different times it breaks and resumes. The GLSL's ribbon modes
//! terminate on **collision only** and freeze ~0.15% of `preset_plambda`; ungated, this port
//! freezes essentially all of it. The distance gate attacks the exposure, not the mechanism:
//! fewer spurious escapes means fewer frozen pixels.
//!
//! So the honest reading of a *partial* improvement is that the gate reduced exposure without
//! removing the patchwork. That is stated up front so it is not discovered as a disappointment.
//!
//! # Panels
//!
//! `_uniform` is the continuous field (hue = shape sphere, lightness = `spread_shape` over the
//! panel's own p1-p99) and is where the arcs live. `_outcome` is the categorical control: for the
//! **basin** config slice the colour IS the terminal outcome, so freezing cannot corrupt it and
//! the two implementations should agree exactly. A disagreement there is a second bug.
//!
//! # Writes
//!
//! `<root>/escgate/`. **Never** `results/charts` -- those are the committed pre-fix renders and
//! the "before" for every claim made here. A reduced-argument validation run has destroyed
//! committed artefacts on this project once already.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 1024);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    let dir = format!("{root}/escgate");
    let _ = std::fs::create_dir_all(&dir);

    // (name, chart, cx, cy, half, body, t_max, r_coll)
    let mut targets: Vec<(String, Chart, f64, f64, f64, usize, f64, f64)> = Vec::new();
    for w in ["preset_plambda", "preset_shape_pl_h1"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            targets.push((c.0.into(), c.1, c.2, c.3, c.4, 0, 13.0, 1e-3));
        }
    }
    // The user's saved Config slices, at their own settings. `horizon = 50` is past the f64
    // measurement horizon (~52 at lambda = 0.7), so the deepest structure is expected to be
    // unresolvable -- a horizon statement, not a disagreement with the reference.
    let (cb, bx, by, bh) = Chart::config_basin();
    targets.push(("config_basin".into(), cb, bx, by, bh, 0, 50.0, 0.02));
    let (cs, sx, sy, sh) = Chart::config_stability();
    targets.push(("config_stability".into(), cs, sx, sy, sh, 0, 50.0, 0.005));

    println!(
        "{res}^2, one sample per pixel, stride 0 (the reference cadence) at both r_esc.\n\
         r_esc = 0 is the shipped-until-now ungated test; r_esc = 5 is the GLSL's gate in\n\
         canonical units (R = 1 on the latent charts, so its literal transfers unconverted).\n"
    );

    for (name, chart, cx, cy, half, body, t_max, r_coll) in targets {
        let name = name.as_str();
        if only != "all" && only != name {
            continue;
        }
        for &r_esc in &[0.0f64, 5.0] {
            let ens = EnsembleCfg {
                refine_flagged: false,
                t_max,
                r_coll_frac: r_coll,
                r_esc_frac: r_esc,
                escape_all_bodies: true,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let dt_sync = t_max / ens.n_sync as f64;
            let mut bits: Vec<u64> = px.iter().map(|p| p.t_end.to_bits()).collect();
            bits.sort_unstable();
            bits.dedup();
            let on_b = px.iter().map(|p| p.t_end).filter(|x| x.is_finite())
                .filter(|x| {
                    let k = (*x / dt_sync).round();
                    (*x - k * dt_sync).abs() <= 1e-9 * t_max
                })
                .count();
            let esc = px.iter().filter(|p| p.state == 0).count();
            let col = px.iter().filter(|p| p.state == 2).count();
            // Freezing exposure: the fraction whose shape_vec is read anywhere but the horizon.
            let frozen = px.iter().filter(|p| p.t_end < t_max - 1e-9).count();

            let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
            let m_here = grid::decode_state(&chart, body, cx, cy).m;
            let sites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(px.len() * 3);
            let mut obuf = Vec::with_capacity(px.len() * 3);
            for p in &px {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
                obuf.extend_from_slice(&png::outcome_rgb(p));
            }
            let stem = format!("{dir}/{name}_resc{}", r_esc as i32);
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
            let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);

            println!(
                "{name:>20} r_esc={r_esc:<5.1} {secs:>7.1}s  escape {:.4} collision {:.4}  \
                 frozen {:.4}  t_end distinct {:>6}  on boundary {:>6.2}%  ramp ({lo:.3e}, {hi:.3e})",
                esc as f64 / px.len() as f64,
                col as f64 / px.len() as f64,
                frozen as f64 / px.len() as f64,
                bits.len(),
                on_b as f64 / px.len() as f64 * 100.0
            );
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         `frozen` is the exposure the patchwork mechanism acts through -- the fraction of pixels\n\
         whose shape_vec is read somewhere other than the horizon. If the gate works, it falls,\n\
         and the ribbons in `_uniform` become continuous. If `frozen` falls and the arcs stay,\n\
         the gate reduced exposure without removing the mechanism, and the next candidates are\n\
         the flatness-based escape criterion or not terminating on escape at all.\n\n\
         `config_basin_*_outcome.png` is the CONTROL. In basin mode the colour is the terminal\n\
         outcome, so freezing cannot corrupt it; the two implementations should agree exactly and\n\
         a disagreement is a SECOND, independent bug."
    );
}
