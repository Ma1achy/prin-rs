//! Before/after renders for the **closure-and-energy** criterion, both toggle states.
//!
//! # The prediction, stated before the run
//!
//! The patchwork is `stop_on_escape` freezing `shape_vec` at each pixel's own `t_end` while the
//! displayed quantity was **still moving**. A ribbon is a level set of `theta` at a common time;
//! where neighbouring pixels froze at different times it breaks and resumes, and since escape is
//! detected at sync boundaries there are only `n_sync` stopping times -- a mosaic of time strata
//! stitched at hard seams.
//!
//! Under this criterion escape fires on a trajectory whose shape has already **stopped moving**
//! (`|dn/dt| ~ 1/t^3`), so the state it freezes is the state it would have had anyway. So:
//!
//! - the domes, tents, wedges and pale seams should be gone, and
//! - **the two toggle rows should look nearly identical**.
//!
//! `shape d` is the per-pixel median chord between the two rows' shape vectors, printed rather
//! than eyeballed -- *an aggregate can only say the distribution did not move; it cannot say the
//! pixels did not*, and that has been read wrong twice on this project.
//!
//! If they are not identical the criterion is still firing before convergence, and this example
//! says so rather than a second fix being stacked on top.
//!
//! # Two colour modes, deliberately
//!
//! `_uniform` is the shipping continuous colouring (hue from the shape sphere, lightness from
//! `spread_shape` over the region's own p1-p99); `_outcome` is the discrete event-class map.
//! Comparing a continuous field against a categorical one is how a rendering choice gets mistaken
//! for a physics bug, so both are written and the mode is in the filename.
//!
//! # Writes
//!
//! `<root>/closure/` -- a **new** directory. The committed renders are the "before" and are not
//! touched; no validation pass writes into `results/charts`.

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};
use prin_rs::outcome::{closure, EscapeRule, CLOSURE_TAU};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The reference's window, in time units; `n_sync` is derived so every case realises it.
const WINDOW: f64 = 0.4;

fn main() {
    let res: usize = arg(1, 1024);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    let sub: String = std::env::args().nth(4).unwrap_or_else(|| "closure".into());
    let tau: f64 = arg(5, CLOSURE_TAU);
    // The subdirectory is an argument, not a constant -- so a regeneration lands beside the
    // committed set rather than over it.
    let dir = format!("{root}/{sub}");
    let _ = std::fs::create_dir_all(&dir);

    // (name, chart, cx, cy, half, body, t_max, r_coll)
    let mut targets: Vec<(String, Chart, f64, f64, f64, usize, f64, f64)> = Vec::new();
    for w in ["preset_plambda", "preset_shape_pl_h1"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            targets.push((c.0.into(), c.1, c.2, c.3, c.4, 0, 13.0, 1e-3));
        }
    }
    // The saved Config slices at their own settings. `horizon = 50` is past the f64 measurement
    // horizon (~52 at lambda = 0.7), so the deepest structure is expected to be unresolvable --
    // a horizon statement, not a disagreement with the reference.
    let (cb, bx, by, bh) = Chart::config_basin();
    targets.push(("config_basin".into(), cb, bx, by, bh, 0, 50.0, 0.02));
    let (cs, sx, sy, sh) = Chart::config_stability();
    targets.push(("config_stability".into(), cs, sx, sy, sh, 0, 50.0, 0.005));

    println!("{res}^2, one sample per pixel, tau={tau:e}, closure_k=1.");
    println!("`n_sync` is derived from `t_max` so every case realises the same ~{WINDOW} window --");
    println!("holding it fixed would compare different discretisations and different criteria.\n");

    for (name, chart, cx, cy, half, body, t_max, r_coll) in targets {
        let name = name.as_str();
        if only != "all" && only != name {
            continue;
        }
        let n_sync = (t_max / WINDOW).round().max(4.0) as usize;
        let mut prev: Option<Vec<[f64; 3]>> = None;
        for &stop in &[false, true] {
            let ens = EnsembleCfg {
                refine_flagged: false,
                t_max,
                n_sync,
                r_coll_frac: r_coll,
                escape_rule: EscapeRule::Closure(tau),
                closure_k: 1,
                stop_on_escape: stop,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, body).with_chart(chart);
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let dt_sync = t_max / n_sync as f64;
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
            // The exposure the patchwork mechanism acts through: pixels whose `shape_vec` is read
            // anywhere other than the horizon.
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
            let stem = format!("{dir}/{name}_stop{}", if stop { 1 } else { 0 });
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
            let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);

            let shapes: Vec<[f64; 3]> = px.iter().map(|p| p.shape_vec).collect();
            let (dmed, dmax, moved) = prev.as_ref().map_or((f64::NAN, f64::NAN, 0), |p| {
                let mut v: Vec<f64> = p.iter().zip(shapes.iter())
                    .map(|(a, b)| closure(a, b))
                    .filter(|x| x.is_finite())
                    .collect();
                let moved = v.iter().filter(|x| **x > 1e-12).count();
                let med = if v.is_empty() { f64::NAN } else { prin_rs::stats::quantile(&mut v, 0.5) };
                let mx = v.iter().cloned().fold(0.0f64, f64::max);
                (med, mx, moved)
            });
            prev = Some(shapes);

            println!(
                "{name:>18} stop={stop:<5} {secs:>7.1}s  escape {:.4} collision {:.4}  \
                 frozen {:.4}  t_end distinct {:>6}  on bdry {:>6.2}%  ramp ({lo:.3e}, {hi:.3e})",
                esc as f64 / px.len() as f64,
                col as f64 / px.len() as f64,
                frozen as f64 / px.len() as f64,
                bits.len(),
                on_b as f64 / px.len() as f64 * 100.0
            );
            if stop {
                println!("{:>18}   shape d: median {dmed:.3e}  max {dmax:.3e}  pixels moved {moved}",
                         "");
            }
        }
        println!();
    }

    println!("HOW TO READ THIS\n");
    println!("`frozen` is the exposure the patchwork acts through. It does NOT have to be small:");
    println!("under this criterion freezing a converged trajectory is close to a no-op, so a high");
    println!("`frozen` with a `shape d` near zero is the *expected* result, not a warning.");
    println!();
    println!("`shape d` decides. Median AND max AND the moved count, because an aggregate can only");
    println!("say the distribution did not move -- twice on this project a row identical to five");
    println!("digits hid every pixel moving, worst 6.7% and 1.86x.");
    println!();
    println!("`config_basin_*_outcome.png` is the CONTROL. In basin mode the colour IS the terminal");
    println!("outcome, so freezing cannot corrupt it and the two toggle rows must agree exactly.");
    println!("A disagreement there is a SECOND fault, in the physics or the event detection, and");
    println!("it needs its own investigation rather than being folded into this one.");
}
