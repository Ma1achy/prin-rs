//! Before/after renders for the `dtau` step-control fix, plus the diagnostic drift panel.
//!
//! # What is being compared
//!
//! **One variable.** The option set is `closure_render`'s -- the settings the committed "before"
//! images were made under -- and only `dtau_mode` changes between the two rows. Anything else
//! moving would make the pair a comparison of two changes.
//!
//! # The three panels
//!
//! `_uniform` is the continuous science field (hue = shape sphere, lightness = `spread_shape`)
//! and `_outcome` is the categorical one. `_drift` is the **diagnostic**: `energy_drift_max` on
//! an inferno ramp, magenta where there is no value, auto-ranged over the field's own p2-p98.
//!
//! **The diagnostic panel is the point.** A science field only shows a numerical defect once it
//! has propagated into an outcome or a spread; the drift map shows it at source, as coherent
//! arcs with the non-finite pixels sitting inside them. That is how this bug was found, and it
//! is why the panel is now standard rather than ad hoc.
//!
//! # Writes
//!
//! `<root>/<sub>/`, both arguments, defaulting to `results/dtau_fix`. **Never** into a committed
//! directory: the existing renders are the "before" for every claim made here, and a reduced
//! validation run has destroyed committed artefacts on this project once already.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::DtauMode;
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// (name, chart, cx, cy, half, t_max, r_coll, both_modes)
type Target = (String, Chart, f64, f64, f64, f64, f64, bool);

fn main() {
    let res: usize = arg(1, 1024);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    let sub: String = std::env::args().nth(4).unwrap_or_else(|| "dtau_fix".into());
    let dir = format!("{root}/{sub}");
    let _ = std::fs::create_dir_all(&dir);

    let mut targets: Vec<Target> = Vec::new();
    // The user's slice, and the two images the artefact was seen on. Both modes: this is the
    // four-panel confirmation the whole change is judged on.
    let (cs, sx, sy, sh) = Chart::config_stability();
    targets.push(("config_stability".into(), cs, sx, sy, sh, 50.0, 0.005, true));
    // The standard preset gallery at the corrected window, fixed mode ON, so it stops carrying
    // the defect. `preset_shape` is the one tree in the corpus the criterion actually controls.
    for w in ["preset_shape", "preset_prho", "preset_plambda", "preset_shape_pl"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            targets.push((c.0.into(), c.1, c.2, c.3, c.4, 13.0, 1e-3, false));
        }
    }

    println!(
        "{res}^2, one sample per pixel. Only `dtau_mode` differs between rows; every other\n\
         setting is `closure_render`'s, which is what the committed \"before\" images used.\n\n\
         `nonfin` is the magenta count -- pixels with a non-finite copy. `simfail` is the\n\
         budget-exhausted count, SEPARATE on purpose: a fix that traded one for the other would\n\
         not have fixed anything. `hot` is the fraction above 1e-6 energy drift -- the speckle.\n\
         `ramp` is each panel's own auto-range window, and it belongs with the image: a clean\n\
         field and a blown-up one both fill the ramp, so the window is where the magnitude is.\n"
    );

    for (name, chart, cx, cy, half, t_max, r_coll, both) in targets {
        let name = name.as_str();
        if only != "all" && only != name {
            continue;
        }
        let modes: Vec<(DtauMode, &str)> = if both {
            vec![(DtauMode::FixedPerInterval, "fixoff"), (DtauMode::PerStepInterval, "fixon")]
        } else {
            vec![(DtauMode::PerStepInterval, "fixon")]
        };
        let mut prev: Option<Vec<PixelOut>> = None;
        for (mode, tag) in modes {
            let ens = EnsembleCfg {
                refine_flagged: false,
                t_max,
                r_coll_frac: r_coll,
                dtau_mode: mode,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let nonfin = px.iter().filter(|p| p.n_nonfinite > 0).count();
            // Separate from `nonfin` **on purpose**: a fix that swapped a diverged copy for a
            // budget-exhausted one would have moved the magenta count and fixed nothing. They are
            // two different failures and they get two columns.
            let simfail = px.iter().filter(|p| p.state == 6).count();
            let hot = px.iter().filter(|p| !(p.energy_drift_max <= 1e-6)).count();
            let esc = px.iter().filter(|p| p.state == 0).count();
            let col = px.iter().filter(|p| p.state == 2).count();

            let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
            // p2-p98 for the diagnostic, stated rather than inherited from the science default.
            let (dlo, dhi) = colour::range_q(&px, Scalar::Drift, 0.02, 0.98);
            let m_here = grid::decode_state(&chart, 0, cx, cy).m;
            let sites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(px.len() * 3);
            let mut obuf = Vec::with_capacity(px.len() * 3);
            let mut dbuf = Vec::with_capacity(px.len() * 3);
            for p in &px {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
                obuf.extend_from_slice(&png::outcome_rgb(p));
                dbuf.extend_from_slice(&colour::drift_rgb(p, dlo, dhi));
            }
            let stem = format!("{dir}/{name}_{tag}");
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
            let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);
            let _ = adaptive::save_rect(&format!("{stem}_drift.png"), res, res, &dbuf);

            println!(
                "{name:>20} {tag:<6} {secs:>7.1}s  escape {:.4} collision {:.4}  \
                 nonfin {:>7}  simfail {:>7}  hot {:.4}  spread ramp ({lo:.3e}, {hi:.3e})  \
                 drift ramp ({dlo:.3e}, {dhi:.3e})",
                esc as f64 / px.len() as f64,
                col as f64 / px.len() as f64,
                nonfin,
                simfail,
                hot as f64 / px.len() as f64,
            );

            // Per-pixel, never the aggregate alone: twice on this project a row identical to five
            // digits hid every pixel moving.
            if let Some(a) = &prev {
                let mut d: Vec<f64> = Vec::with_capacity(px.len());
                let mut moved = 0usize;
                let mut lab = 0usize;
                for (x, y) in a.iter().zip(px.iter()) {
                    let s: f64 = (0..3).map(|k| (x.shape_vec[k] - y.shape_vec[k]).powi(2)).sum();
                    let s = s.sqrt();
                    if s > 0.0 {
                        moved += 1;
                    }
                    if x.outcome != y.outcome {
                        lab += 1;
                    }
                    d.push(s);
                }
                let mut f: Vec<f64> = d.iter().cloned().filter(|x| x.is_finite()).collect();
                let med = prin_rs::stats::quantile(&mut f.clone(), 0.5);
                println!(
                    "{:>20}        shape d: median {med:.3e}  max {:.3e}  pixels moved {moved}  \
                     outcome labels flipped {lab}",
                    "",
                    f.iter().cloned().fold(0.0, f64::max)
                );
            }
            prev = Some(px);
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         THE PREDICTION, stated so it can fail: the magenta clusters shrink sharply or vanish,\n\
         the speckle halos around them go with them and the outcome bands become SOLID as in the\n\
         reference GLSL, the continuous ribbons stop being interrupted, and nothing new appears.\n\n\
         `nonfin` and `hot` are the numbers behind the first two clauses. `outcome labels\n\
         flipped` is the number behind the third: the speckle IS label flipping, so if the\n\
         bands are solidifying that count is large and concentrated where the halos were.\n\n\
         If the magenta survives, that is a SECOND cause and it is worth knowing before\n\
         anything is built on top of this one."
    );
}
