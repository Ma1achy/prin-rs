//! **Why is THIS slice fragile?** -- reference-body switching under asymmetric masses.
//!
//! `config_stability` carries ~50x the density of sharp discontinuities in its drift map that
//! `preset_prho` does. A chaotic field gives smooth gradients and fractal boundaries; a *discrete*
//! decision boundary gives a step. AZ chooses its reference body from the longest side, and every
//! change re-derives the Levi-Civita registration mid-trajectory -- a discrete choice, made on a
//! slice whose masses are `(0.32735, 0.42763, 0.24502)` rather than the presets' equal thirds.
//!
//! **The reference is chosen once per SYNC BOUNDARY, not per RK4 step** (`driver.rs`, the
//! `for kk in 0..n_sync` head). So one of the obvious remedies is already the shipped behaviour,
//! and switch times are quantised to `t_max / n_sync` by construction.
//!
//! Four measurements, and the fourth is the only one that can catch the mechanism in the act:
//!
//! 1. switch **count** per pixel, as its own panel;
//! 2. **first-switch time** per pixel, as a field -- if it carries the drift map's structure they
//!    are the same phenomenon;
//! 3. the **alignment test**: is a large drift step between neighbours enriched where the two
//!    pixels have different switch histories? Reported against the matched-history population,
//!    and against a **shuffled** control, because on a field where both quantities are spatially
//!    smooth any two maps agree somewhat;
//! 4. **drift accumulated per switch**: the paired increment `|drift[k] - drift[k-1]|` at the
//!    boundaries where the reference changed, against the boundaries where it held. A correlation
//!    between two maps can be produced by both tracking a third thing; a paired increment across
//!    the switch cannot.
//!
//! `preset_prho` is the control on every one of them. A metric that does not separate the two
//! slices is not measuring what differs.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::{adaptive, png};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn qt(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.015], [0.34, 0.06, 0.43], [0.72, 0.21, 0.33],
        [0.98, 0.55, 0.04], [0.99, 1.0, 0.64],
    ];
    let t = x.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    let mut o = [0u8; 3];
    for k in 0..3 {
        o[k] = (255.0 * (S[i][k] * (1.0 - f) + S[i + 1][k] * f)).clamp(0.0, 255.0) as u8;
    }
    o
}

/// Luminance of the committed drift panel, so the neighbour-gradient census reproduces the
/// measurement that prompted this rather than inventing a new one.
fn drift_lum(v: f64, nonfin: bool) -> f64 {
    if nonfin || !v.is_finite() {
        return f64::NAN;
    }
    let (lo, hi) = (1.0e-8f64.log10(), 4.0e7f64.log10());
    let x = ((v.max(1e-300).log10() - lo) / (hi - lo)).clamp(0.0, 1.0);
    let c = ramp(x);
    0.2126 * c[0] as f64 + 0.7152 * c[1] as f64 + 0.0722 * c[2] as f64
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(512);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "switch_out".into());
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let cs = Chart::Latent { z0, q1, q2 };
    let mut cases: Vec<(&'static str, Chart, f64, f64, f64, f64, f64)> = vec![(
        "config_stability", cs,
        2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM, 50.0, 0.005,
    )];
    // `preset_prho` is the smooth control the hypothesis is aimed at. `preset_shape` and
    // `preset_plambda` are the ones that DISCRIMINATE: measured on the committed drift panels,
    // `preset_shape` carries MORE sharp gradients than `config_stability` (0.0032 against 0.0026
    // above 200), so the contrast is not "this slice against the presets" -- it is `preset_prho`
    // against everything. If switch rate is the mechanism, `preset_shape` must switch too.
    for w in ["preset_prho", "preset_shape", "preset_plambda"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            cases.push((w, c.1, c.2, c.3, c.4, 13.0, 1e-3));
        }
    }

    println!("{res}^2. `keep_drift_hist` on, so the drift series and `refs` share one cadence.\n");

    for (name, chart, cx, cy, half, t_max, r_coll) in cases {
        let n_sync = 32usize;
        let ens = EnsembleCfg {
            t_max, n_sync, r_coll_frac: r_coll, refine_flagged: false,
            keep_drift_hist: true, ..Default::default()
        };
        let m = grid::decode_state(&chart, 0, cx, cy).m;
        let t0 = std::time::Instant::now();
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let px: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
            .collect();
        let secs = t0.elapsed().as_secs_f64();
        let n = px.len();

        // ---- panels -----------------------------------------------------------------------
        let sw_max = px.iter().map(|p| p.switches).max().unwrap_or(1).max(1) as f64;
        let mut sbuf = Vec::with_capacity(n * 3);
        let mut tbuf = Vec::with_capacity(n * 3);
        let mut dbuf = Vec::with_capacity(n * 3);
        for p in &px {
            sbuf.extend_from_slice(&ramp(p.switches as f64 / sw_max));
            tbuf.extend_from_slice(&if p.t_first_switch.is_finite() {
                ramp(p.t_first_switch / t_max)
            } else {
                [255, 0, 255]
            });
            let l = drift_lum(p.energy_drift_max, p.n_nonfinite > 0);
            dbuf.extend_from_slice(&if l.is_finite() {
                ramp(((p.energy_drift_max.max(1e-300).log10() - (-8.0)) / (7.602 + 8.0)).clamp(0.0, 1.0))
            } else {
                [255, 0, 255]
            });
        }
        let stem = format!("{dir}/{name}");
        let _ = adaptive::save_rect(&format!("{stem}_switches.png"), res, res, &sbuf);
        let _ = adaptive::save_rect(&format!("{stem}_tfirstswitch.png"), res, res, &tbuf);
        let _ = adaptive::save_rect(&format!("{stem}_drift.png"), res, res, &dbuf);
        let mut obuf = Vec::with_capacity(n * 3);
        for p in &px {
            obuf.extend_from_slice(&png::outcome_rgb(p));
        }
        let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);

        // ---- 1. switch census -------------------------------------------------------------
        let mut sw: Vec<f64> = px.iter().map(|p| p.switches as f64).collect();
        let hist: Vec<usize> = (0..6)
            .map(|k| px.iter().filter(|p| p.switches as usize == k).count())
            .collect();
        println!(
            "{name:>18}  masses {m:.5?}  {secs:.1}s\n\
             {:>18}  switches: p50 {:.0} p90 {:.0} max {:.0}  mean {:.3}  \
             hist[0..5] {hist:?} (of {n})",
            "", qt(&mut sw.clone(), 0.5), qt(&mut sw.clone(), 0.9),
            qt(&mut sw, 1.0),
            px.iter().map(|p| p.switches as f64).sum::<f64>() / n as f64,
        );
        let never = px.iter().filter(|p| p.switches == 0).count();
        let mut tf: Vec<f64> = px.iter().map(|p| p.t_first_switch).filter(|x| x.is_finite()).collect();
        println!(
            "{:>18}  never switched {never} ({:.4})   t_first p10 {:.3} p50 {:.3} p90 {:.3}  \
             (quantised to t_max/n_sync = {:.4})",
            "", never as f64 / n as f64,
            qt(&mut tf.clone(), 0.1), qt(&mut tf.clone(), 0.5), qt(&mut tf, 0.9),
            t_max / n_sync as f64,
        );

        // ---- 4. the paired increment ------------------------------------------------------
        let mut sj: Vec<f64> = px.iter().map(|p| p.switch_jump_med).filter(|x| x.is_finite()).collect();
        let mut hj: Vec<f64> = px.iter().map(|p| p.hold_jump_med).filter(|x| x.is_finite()).collect();
        let mut sx: Vec<f64> = px.iter().map(|p| p.switch_jump_max).filter(|x| x.is_finite()).collect();
        let mut hx: Vec<f64> = px.iter().map(|p| p.hold_jump_max).filter(|x| x.is_finite()).collect();
        // Paired, over pixels that HAVE both arms -- the unpaired medians above compare
        // different populations, which is the defect already on record for a selection-
        // conditioned median.
        let mut pr: Vec<f64> = px
            .iter()
            .filter(|p| p.switch_jump_med.is_finite() && p.hold_jump_med.is_finite() && p.hold_jump_med > 0.0)
            .map(|p| p.switch_jump_med / p.hold_jump_med)
            .filter(|x| x.is_finite())
            .collect();
        let npair = pr.len();
        println!(
            "{:>18}  jump med: switch {:.3e} (n={})  hold {:.3e} (n={})   \
             max: switch {:.3e} hold {:.3e}",
            "", qt(&mut sj.clone(), 0.5), sj.len(), qt(&mut hj.clone(), 0.5), hj.len(),
            qt(&mut sx, 0.5), qt(&mut hx, 0.5),
        );
        println!(
            "{:>18}  PAIRED switch/hold ratio over {npair} pixels with both arms: \
             p10 {:.3} p50 {:.3} p90 {:.3}   frac>1 {:.4}",
            "", qt(&mut pr.clone(), 0.1), qt(&mut pr.clone(), 0.5), qt(&mut pr.clone(), 0.9),
            pr.iter().filter(|x| **x > 1.0).count() as f64 / npair.max(1) as f64,
        );

        // ---- 3. the alignment test --------------------------------------------------------
        // Neighbour pairs, both directions. A "step" is a large jump in the RENDERED drift
        // luminance, which is the quantity the 6.19%/0.12% census was taken on.
        let lum: Vec<f64> = px.iter().map(|p| drift_lum(p.energy_drift_max, p.n_nonfinite > 0)).collect();
        // The SHIFTED control. Both maps are spatially smooth, so any two of them agree
        // somewhat; shifting the switch map by a prime stride preserves its own spatial
        // statistics and destroys only the alignment. An enrichment that survives the shift is
        // measuring smoothness, not coincidence.
        for thr in [100.0f64, 200.0] {
            let (mut pairs, mut steps) = (0usize, 0usize);
            let (mut cd, mut csd, mut cs, mut css) = (0usize, 0usize, 0usize, 0usize);
            let (mut diff_hist, mut step_and_diff) = (0usize, 0usize);
            let (mut same_hist, mut step_and_same) = (0usize, 0usize);
            let mut push = |i: usize, j: usize,
                            pairs: &mut usize, steps: &mut usize,
                            dh: &mut usize, sd: &mut usize, sh: &mut usize, ss: &mut usize| {
                let (a, b) = (&px[i], &px[j]);
                if !lum[i].is_finite() || !lum[j].is_finite() {
                    return;
                }
                *pairs += 1;
                let step = (lum[i] - lum[j]).abs() > thr;
                if step {
                    *steps += 1;
                }
                // "Different switch history" = different count, or a first switch at a
                // different boundary. NaN-vs-finite counts as different.
                let ta = a.t_first_switch;
                let tb = b.t_first_switch;
                let tdiff = if ta.is_nan() != tb.is_nan() {
                    true
                } else if ta.is_nan() {
                    false
                } else {
                    (ta - tb).abs() > 1e-9
                };
                if a.switches != b.switches || tdiff {
                    *dh += 1;
                    if step {
                        *sd += 1;
                    }
                } else {
                    *sh += 1;
                    if step {
                        *ss += 1;
                    }
                }
            };
            for y in 0..res {
                for x in 0..res {
                    let i = y * res + x;
                    if x + 1 < res {
                        push(i, i + 1, &mut pairs, &mut steps, &mut diff_hist,
                             &mut step_and_diff, &mut same_hist, &mut step_and_same);
                    }
                    if y + 1 < res {
                        push(i, i + res, &mut pairs, &mut steps, &mut diff_hist,
                             &mut step_and_diff, &mut same_hist, &mut step_and_same);
                    }
                }
            }
            // Shifted control: same census, but each pixel's switch history is taken from a
            // pixel 37 rows and 53 columns away.
            let sh = |i: usize| -> usize {
                let (x, y) = (i % res, i / res);
                ((y + 37) % res) * res + (x + 53) % res
            };
            for y in 0..res {
                for x in 0..res {
                    let i = y * res + x;
                    for j in [
                        if x + 1 < res { Some(i + 1) } else { None },
                        if y + 1 < res { Some(i + res) } else { None },
                    ]
                    .into_iter()
                    .flatten()
                    {
                        if !lum[i].is_finite() || !lum[j].is_finite() {
                            continue;
                        }
                        let step = (lum[i] - lum[j]).abs() > thr;
                        let (a, b) = (&px[sh(i)], &px[sh(j)]);
                        let (ta, tb) = (a.t_first_switch, b.t_first_switch);
                        let tdiff = if ta.is_nan() != tb.is_nan() {
                            true
                        } else if ta.is_nan() {
                            false
                        } else {
                            (ta - tb).abs() > 1e-9
                        };
                        if a.switches != b.switches || tdiff {
                            cd += 1;
                            if step {
                                csd += 1;
                            }
                        } else {
                            cs += 1;
                            if step {
                                css += 1;
                            }
                        }
                    }
                }
            }
            let p_step = steps as f64 / pairs.max(1) as f64;
            let p_diff = step_and_diff as f64 / diff_hist.max(1) as f64;
            let p_same = step_and_same as f64 / same_hist.max(1) as f64;
            println!(
                "{:>18}  |grad lum|>{thr:.0}: {:.4} of {pairs} pairs.  \
                 P(step | history DIFFERS) {p_diff:.4} over {diff_hist}   \
                 P(step | history MATCHES) {p_same:.4} over {same_hist}   \
                 enrichment {:.2}x   lift-vs-all {:.2}x",
                "", p_step,
                if p_same > 0.0 { p_diff / p_same } else { f64::NAN },
                if p_step > 0.0 { p_diff / p_step } else { f64::NAN },
            );
            let (qd, qs) = (csd as f64 / cd.max(1) as f64, css as f64 / cs.max(1) as f64);
            println!(
                "{:>18}    shifted control (history taken 37 rows / 53 cols away): \
                 P(step | differs) {qd:.4}  P(step | matches) {qs:.4}  enrichment {:.2}x",
                "",
                if qs > 0.0 { qd / qs } else { f64::NAN },
            );
        }
        println!();
    }
}
