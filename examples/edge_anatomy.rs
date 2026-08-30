//! **What draws the hard edges in the drift field, now that the step limit is in?**
//!
//! The panel under `StepLimit::Predictive` is mostly smooth swirls with sharp-edged wedges and
//! speckled filaments cutting through them. Four candidates, and each gets a test that can come
//! back negative:
//!
//! ```text
//!   1  the ORDER STATISTIC   `energy_drift_max` is a max over 8 copies -- eight overlapping
//!                            sheets, and which copy wins changes pixel to pixel. Compare
//!                            against `energy_drift_nominal`, which carries no order statistic.
//!   2  the REFERENCE FLIP    `choose_reference` is a bare argmax, re-evaluated at all 125 sync
//!                            boundaries, no hysteresis. Neighbours that take different paths
//!                            integrate in different coordinates. `ref_path_hash` makes this
//!                            exact rather than a proxy on the switch COUNT -- two neighbours can
//!                            switch equally often at different boundaries.
//!   3  the STEP COUNT        an integer. If its jumps do not line up with the drift edges, the
//!                            "landing sawtooth" is a non-problem and should not be built for.
//!   4  chaos                 genuine divergence of neighbouring ICs. Not removable, and the
//!                            null hypothesis the other three have to beat.
//! ```
//!
//! # The statistic
//!
//! For each candidate field, the **edge set** is the top decile of `|grad log10|` on the drift
//! field, and the test is whether the candidate's own discontinuity set coincides with it.
//! Reported as a lift — `P(candidate edge | drift edge) / P(candidate edge)` — because a
//! candidate that fires on half the frame has a lift of ~1 by arithmetic whatever it explains.
//! The base rate is printed first for that reason.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::StepBlend;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
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

/// `|grad|` of a scalar field by max forward difference, in raster order.
fn grad(f: &[f64], res: usize) -> Vec<f64> {
    (0..f.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let mut g: f64 = 0.0;
            if x + 1 < res {
                g = g.max((f[i + 1] - f[i]).abs());
            }
            if y + 1 < res {
                g = g.max((f[i + res] - f[i]).abs());
            }
            g
        })
        .collect()
}

/// True where the pixel differs from its right or lower neighbour.
fn differs<T: PartialEq + Copy>(f: &[T], res: usize) -> Vec<bool> {
    (0..f.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            (x + 1 < res && f[i + 1] != f[i]) || (y + 1 < res && f[i + res] != f[i])
        })
        .collect()
}

fn lg0(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 {
        x.log10()
    } else {
        DLO.log10()
    }
}

fn main() {
    let res: usize = arg(1, 384);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/edges");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let cfg = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        keep_ref_path: true,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", cfg.provenance());

    // Three arms in one run: the shipped hard `min`, and the soft minimum at two exponents.
    // `p -> inf` is the `min`, `p = 1` is the harmonic form, `p = 4` is the compromise.
    let arms: Vec<(String, EnsembleCfg)> = vec![
        ("Min (shipped)".into(), cfg),
        ("SoftMin p=4".into(), EnsembleCfg { step_blend: StepBlend::SoftMin, blend_p: 4.0, ..cfg }),
        ("SoftMin p=1".into(), EnsembleCfg { step_blend: StepBlend::SoftMin, blend_p: 1.0, ..cfg }),
    ];
    println!(
        "  {:>16} {:>9} {:>12} {:>12} {:>12} {:>12} {:>11}",
        "arm", "secs", "rough p50", "rough p90", "steps p50", "err p99", "overshoot"
    );
    let mut first: Option<Vec<PixelOut>> = None;
    for (label, acfg) in &arms {
        let t = std::time::Instant::now();
        let v: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, acfg))
            .collect();
        let secs = t.elapsed().as_secs_f64();
        let d: Vec<f64> = v.iter().map(|p| lg0(p.energy_drift_max)).collect();
        let mut g = grad(&d, res);
        let mut st: Vec<f64> = v.iter().map(|p| p.total_substeps as f64).collect();
        let mut er: Vec<f64> =
            v.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
        println!(
            "  {label:>16} {secs:>9.1} {:>12.4} {:>12.4} {:>12.3e} {:>12.3e} {:>11}",
            q(&mut g.clone(), 0.5),
            q(&mut g, 0.9),
            q(&mut st, 0.5),
            q(&mut er, 0.99),
            v.iter().map(|p| p.n_overshoot).sum::<u64>()
        );
        // The blend panel, on the same fixed ramp as everything else.
        let name = if label.contains("p=4") {
            "drift_softmin_p4"
        } else if label.contains("p=1") {
            "drift_softmin_p1"
        } else {
            "drift_min"
        };
        let buf: Vec<u8> = d
            .iter()
            .flat_map(|&x| ramp((x - DLO.log10()) / (DHI.log10() - DLO.log10())).into_iter())
            .collect();
        let path = format!("{dir}/{name}.png");
        let _ = prin_rs::output::adaptive::save_rect(&path, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &path, acfg, &format!("res={res}x{res}\ncase=config_stability\narm={label}\n"),
        );
        if first.is_none() {
            first = Some(v);
        }
    }
    println!(
        "\n  **`rough` is the median |grad log10 drift| -- a per-pixel number, not an image\n\
         impression.** If the soft minimum removes the constraint-switching creases, roughness\n\
         falls. If it does not, the creases were not what was drawing the edges and the cost is\n\
         being paid for nothing.\n"
    );
    let px = first.unwrap();

    let dmax: Vec<f64> = px.iter().map(|p| lg0(p.energy_drift_max)).collect();
    let dnom: Vec<f64> = px.iter().map(|p| lg0(p.energy_drift_nominal)).collect();

    // --- 1. the order statistic ----------------------------------------------------------
    println!("== 1. IS IT THE ORDER STATISTIC? ==");
    println!(
        "  `energy_drift_max` is a max over {} copies; `energy_drift_nominal` is one copy and\n\
         carries no order statistic. **If the edges are the max, they soften here.** Roughness is\n\
         the median of |grad log10| -- a per-pixel number, not an image impression.\n",
        cfg.n_extra + 1
    );
    let gmax = grad(&dmax, res);
    let gnom = grad(&dnom, res);
    let (mut a, mut b) = (gmax.clone(), gnom.clone());
    let (rm, rn) = (q(&mut a, 0.5), q(&mut b, 0.5));
    let (mut a9, mut b9) = (gmax.clone(), gnom.clone());
    println!(
        "  roughness (median |grad log10 drift|):  max {rm:.4}   nominal {rn:.4}   ratio {:.3}",
        rm / rn.max(f64::MIN_POSITIVE)
    );
    println!(
        "  p90 |grad|:                             max {:.4}   nominal {:.4}",
        q(&mut a9, 0.9),
        q(&mut b9, 0.9)
    );

    // --- 2. the reference flip -------------------------------------------------------------
    let hashes: Vec<u64> = px.iter().map(|p| p.ref_path_hash).collect();
    let flip = differs(&hashes, res);
    let mut steps: Vec<u64> = px.iter().map(|p| p.total_substeps).collect();
    let stepdiff = differs(&steps, res);

    // The drift edge set: top decile of the gradient.
    let mut gs = gmax.clone();
    let cut = q(&mut gs, 0.90);
    let edge: Vec<bool> = gmax.iter().map(|&g| g >= cut).collect();

    let frac = |m: &[bool]| m.iter().filter(|x| **x).count() as f64 / m.len() as f64;
    let lift = |m: &[bool]| {
        let n_e = edge.iter().filter(|x| **x).count().max(1) as f64;
        let p_e = (0..m.len()).filter(|&i| edge[i] && m[i]).count() as f64 / n_e;
        p_e / frac(m).max(f64::MIN_POSITIVE)
    };

    println!("\n== 2 & 3. WHAT COINCIDES WITH THE DRIFT EDGES? ==");
    println!(
        "  Base rates first: **a candidate covering half the frame has a lift of ~1 by\n\
         arithmetic**, whatever it explains.\n"
    );
    println!("  {:>34} {:>10} {:>12} {:>8}", "candidate", "base rate", "P(cand|edge)", "lift");
    let n_e = edge.iter().filter(|x| **x).count().max(1) as f64;
    for (name, m) in [
        ("reference path differs (flip)", &flip),
        ("step count differs", &stepdiff),
    ] {
        let p_e = (0..m.len()).filter(|&i| edge[i] && m[i]).count() as f64 / n_e;
        println!("  {name:>34} {:>10.4} {:>12.4} {:>8.3}", frac(m), p_e, lift(m));
    }
    println!(
        "  {:>34} {:>10.4}",
        "drift edge set (top decile)",
        frac(&edge)
    );

    // --- the graded arms, because both binary masks are near-saturated -------------------
    // Hamming distance between neighbouring reference SEQUENCES: how many of the `n_sync`
    // boundaries two adjacent pixels disagree at. A binary differs-or-not mask reads 0.72 here
    // and cannot discriminate; this can.
    let ham: Vec<f64> = (0..px.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let d = |j: usize| -> f64 {
                px[i]
                    .ref_path
                    .iter()
                    .zip(px[j].ref_path.iter())
                    .filter(|(a, b)| a != b)
                    .count() as f64
            };
            let mut m: f64 = 0.0;
            if x + 1 < res {
                m = m.max(d(i + 1));
            }
            if y + 1 < res {
                m = m.max(d(i + res));
            }
            m
        })
        .collect();
    let dstep: Vec<f64> = (0..px.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let d = |j: usize| (px[i].total_substeps as f64 - px[j].total_substeps as f64).abs();
            let mut m: f64 = 0.0;
            if x + 1 < res {
                m = m.max(d(i + 1));
            }
            if y + 1 < res {
                m = m.max(d(i + res));
            }
            m
        })
        .collect();

    println!("\n== THE GRADED FORM, and it is the one that decides ==");
    println!(
        "  Both masks above are near-saturated, so their lifts cannot discriminate -- that is a\n\
         defect in the test, stated rather than worked around. These are magnitudes.\n"
    );
    let mut hs = ham.clone();
    let mut ds = dstep.clone();
    println!(
        "  neighbour reference-path Hamming distance (of {} boundaries): p50 {:.1} p90 {:.1} max {:.0}",
        cfg.n_sync,
        q(&mut hs.clone(), 0.5),
        q(&mut hs.clone(), 0.9),
        q(&mut hs, 1.0)
    );
    println!(
        "  neighbour |delta step count|:                             p50 {:.0} p90 {:.0} max {:.0}",
        q(&mut ds.clone(), 0.5),
        q(&mut ds.clone(), 0.9),
        q(&mut ds, 1.0)
    );

    // Conditional roughness: inside a coherent cell against across its boundary. **This is the
    // test.** If the reference partition draws the edges, the drift gradient is much larger where
    // the paths differ, and a step-size change cannot touch it.
    let same: Vec<usize> = (0..px.len()).filter(|&i| ham[i] == 0.0).collect();
    let diff: Vec<usize> = (0..px.len()).filter(|&i| ham[i] > 0.0).collect();
    let mut gsame: Vec<f64> = same.iter().map(|&i| gmax[i]).collect();
    let mut gdiff: Vec<f64> = diff.iter().map(|&i| gmax[i]).collect();
    println!(
        "\n  median |grad log10 drift| INSIDE a coherent cell (same reference path): {:.4}  (n = {})",
        q(&mut gsame, 0.5),
        same.len()
    );
    println!(
        "  median |grad log10 drift| ACROSS a path difference:                     {:.4}  (n = {})",
        q(&mut gdiff, 0.5),
        diff.len()
    );
    let ratio = q(&mut gdiff.clone(), 0.5) / q(&mut gsame.clone(), 0.5).max(f64::MIN_POSITIVE);
    println!(
        "  ratio: {ratio:.3}\n\n  **A large ratio means the reference-body argmax draws the edges.**\n\
         It is a coordinate choice the exact solution does not care about, so the jump is purely\n\
         the DIFFERENCE in integration error between two coordinate systems -- shrinkable by\n\
         making both errors smaller, and not removable by any step-size rule."
    );

    println!(
        "\n  **If `step count differs` has a lift near 1, the landing sawtooth is a non-problem**\n\
         and the uniform-N change should not be built on that premise. If `reference path\n\
         differs` carries the lift, the edges are the argmax and no step-size change removes them."
    );

    // --- panels ----------------------------------------------------------------------------
    let mut save = |name: &str, buf: &[u8]| {
        let path = format!("{dir}/{name}.png");
        let _ = prin_rs::output::adaptive::save_rect(&path, res, res, buf);
        let _ = prin_rs::output::provenance_sidecar(
            &path,
            &cfg,
            &format!("res={res}x{res}\ncase=config_stability\nfield={name}\n"),
        );
    };
    let ramped = |v: &[f64]| -> Vec<u8> {
        v.iter()
            .flat_map(|&x| {
                ramp((x - DLO.log10()) / (DHI.log10() - DLO.log10())).into_iter()
            })
            .collect()
    };
    save("drift_max", &ramped(&dmax));
    save("drift_nominal", &ramped(&dnom));
    let mask = |m: &[bool]| -> Vec<u8> {
        m.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect()
    };
    save("edge_drift", &mask(&edge));
    save("flip_reference_path", &mask(&flip));
    save("flip_step_count", &mask(&stepdiff));
    // Step count as a field, so the integer terraces are visible rather than inferred.
    let mut sv: Vec<f64> = steps.iter().map(|&s| s as f64).collect();
    let (slo, shi) = (q(&mut sv.clone(), 0.02), q(&mut sv, 0.98));
    save(
        "field_step_count",
        &px.iter()
            .flat_map(|p| {
                ramp((p.total_substeps as f64 - slo) / (shi - slo).max(1.0)).into_iter()
            })
            .collect::<Vec<u8>>(),
    );
    steps.sort_unstable();

    println!("\nWrote {dir}/ -- drift_max, drift_nominal, edge_drift, flip_*, field_step_count.");
}
