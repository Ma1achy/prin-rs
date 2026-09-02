//! **Does Heggie actually remove the wedges?** The claim the port was FOR, measured.
//!
//! Every AZ-vs-Heggie number in this corpus is an **error magnitude** -- drift quantiles, `err>10`,
//! gain maps. But the motivating defect was never the error: it was the *spatial artefact* that
//! AZ's reference-body choice writes into the field. `config_stability`'s wedges are a **discrete
//! decision boundary** printed across a continuous slice, and the case that Heggie removes them
//! has rested on looking at the panels.
//!
//! `switch_study.rs` has the machinery -- neighbour-gradient census, alignment test, shuffled
//! control -- and it is **AZ-only**. This is that census with an integrator arm.
//!
//! # The confound that decides the design: the metric MUST be scale-free
//!
//! `switch_study` measures neighbour steps in the **luminance of the committed drift panel**, a
//! fixed `1e-8 .. 4e7` log ramp. That is right for comparing two slices under one integrator and
//! **wrong for comparing integrators**: Heggie's drift is 1-4 orders lower, so it occupies a
//! different part of the ramp and would show fewer sharp luminance steps *for that reason alone*.
//! It would score artefact-free without being it -- the mirror of this project's auto-ranged-ramp
//! trap.
//!
//! So the field is `log10(drift)` and a step is measured in **decades**. A uniform k-fold change
//! of the whole field shifts `log10` by a constant and leaves every neighbour difference
//! unchanged. That invariance is the point, and it is asserted below: the census is re-run on a
//! deliberately rescaled copy of the field and must return the identical count.
//!
//! # Density alone cannot tell a decision boundary from chaos
//!
//! A fractal boundary is *also* full of sharp neighbour steps. What separates them is **shape**:
//! a discrete decision boundary is a long connected curve, chaotic mixing is scattered dust.
//! So the census reports the sharp set's **largest connected component** beside its density --
//! this project's own lesson from the wedge mask, where a threshold selecting 11640 components of
//! largest 211 px would have "confirmed" a geometric edge that was not there. Component size
//! caught it; the dimension number did not.
//!
//! # Controls, and a metric that fails them is not measuring the artefact
//!
//! - **`preset_prho` is the slice control.** Equal masses, so the reference-body choice is not
//!   ambiguous and there are no wedges. Any integrator must read LOW there. A metric that does not
//!   separate `config_stability` from `preset_prho` **under AZ** is not measuring this.
//! - **`plain_rk4` is the integrator control**: no regularisation at all, so no chart to
//!   re-register into and no wedge mechanism -- but a wrecked field, which is what says the metric
//!   is reading structure and not damage.
//! - The threshold is **swept**, not chosen. A single constant would be a knob picked to make a
//!   picture look right.
//!
//! ```text
//! cargo run --release --example artefact_census -- [res] [out_dir]
//! ```
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::Integrator;

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

/// Neighbour steps above `thr` decades, as a fraction of 4-adjacent pairs where both are finite.
///
/// `NaN` pixels are excluded from the DENOMINATOR as well as the numerator: an undetermined pixel
/// has no gradient, and counting it as smooth would reward a field for being broken.
fn sharp_pairs(f: &[f64], n: usize, thr: f64) -> (f64, usize, usize) {
    let (mut sharp, mut pairs) = (0usize, 0usize);
    let at = |x: usize, y: usize| f[y * n + x];
    for y in 0..n {
        for x in 0..n {
            let a = at(x, y);
            if !a.is_finite() {
                continue;
            }
            for (dx, dy) in [(1usize, 0usize), (0, 1)] {
                let (u, v) = (x + dx, y + dy);
                if u >= n || v >= n {
                    continue;
                }
                let b = at(u, v);
                if !b.is_finite() {
                    continue;
                }
                pairs += 1;
                if (a - b).abs() > thr {
                    sharp += 1;
                }
            }
        }
    }
    (sharp as f64 / pairs.max(1) as f64, sharp, pairs)
}

/// The sharp SET -- pixels touching a sharp pair -- as a mask.
fn sharp_mask(f: &[f64], n: usize, thr: f64) -> Vec<bool> {
    let mut hot = vec![false; n * n];
    let at = |x: usize, y: usize| f[y * n + x];
    for y in 0..n {
        for x in 0..n {
            let a = at(x, y);
            if !a.is_finite() {
                continue;
            }
            for (dx, dy) in [(1i32, 0i32), (0, 1), (-1, 0), (0, -1)] {
                let (u, v) = (x as i32 + dx, y as i32 + dy);
                if u < 0 || v < 0 || u >= n as i32 || v >= n as i32 {
                    continue;
                }
                let b = at(u as usize, v as usize);
                if b.is_finite() && (a - b).abs() > thr {
                    hot[y * n + x] = true;
                }
            }
        }
    }
    hot
}

/// **Straightness of the sharp set, which is what separates a decision boundary from chaos.**
///
/// Density and connectedness both rise for a fractal boundary *and* for a wrecked field --
/// measured: `plain_rk4`, which has no reference body and therefore no wedge mechanism at all,
/// scores a large connected component just as AZ does. Only a *decision* boundary is straight.
///
/// Returns `(straightness of the largest component, median over components >= MIN, count)`.
/// Components under `MIN` px are excluded and the count is printed, because a short component is
/// straight by construction and a median over dust would read as a boundary.
fn straight_stats(mask: &[bool], n: usize) -> (f64, f64, usize) {
    const MIN: usize = 20;
    let comps = prin_rs::spatial::components(mask, n);
    let largest = comps.first().map(|c| prin_rs::spatial::straightness(c)).unwrap_or(f64::NAN);
    let mut ss: Vec<f64> = comps
        .iter()
        .filter(|c| c.len() >= MIN)
        .map(|c| prin_rs::spatial::straightness(c))
        .filter(|x| x.is_finite())
        .collect();
    let k = ss.len();
    ss.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if ss.is_empty() { f64::NAN } else { ss[ss.len() / 2] };
    (largest, med, k)
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/artefact".into());
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
    for w in ["preset_prho"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            cases.push((w, c.1, c.2, c.3, c.4, 13.0, 1e-3));
        }
    }

    let arms = [
        (Integrator::Az, "az"),
        (Integrator::Heggie, "heggie"),
        (Integrator::LogHRk4, "logh_rk4"),
        (Integrator::PlainRk4, "plain_rk4"),
    ];
    let thrs = [0.5f64, 1.0, 1.5, 2.0];

    println!("SPATIAL ARTEFACT CENSUS at {res}^2. Diagnostic pass: termination OFF, r_coll = 0,");
    println!("refine_flagged OFF -- the repair pass removes the population this is about.");
    println!();
    println!("  Field is `log10(energy_drift_max)` and a step is in DECADES, so the census is");
    println!("  invariant to a uniform rescaling of the field. That is asserted per row: the");
    println!("  census is re-run on the field times 1e6 and must return the SAME count.");
    println!();
    println!("  `preset_prho` is the slice control (equal masses, no wedges) and `plain_rk4` the");
    println!("  integrator control (no chart at all). A metric that fails to separate");
    println!("  config_stability from preset_prho UNDER AZ is not measuring the artefact.");
    println!();

    for (name, chart, cx, cy, half, t_max, r_coll) in cases {
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        println!("== {name}  (t_max = {t_max})");
        println!(
            "{:>10} {:>11} {:>6} {:>9} {:>8} {:>10} {:>9} {:>8} {:>6} {:>10} {:>3}",
            "arm", "nonfin/bud", "thr", "sharp", "largest", "lgst/shrp", "str_lgst", "str_med", "ncomp", "drift p50", "sf"
        );
        for (integ, label) in arms {
            let cfg = EnsembleCfg::production().with_overrides(&[
                Override::Integrator(integ),
                Override::TMax(t_max),
                Override::RCollFrac(r_coll),
                Override::StopOnEvent(false),
                Override::RefineFlagged(false),
                // **The budget must not exclude anyone, or each arm is scored on a different
                // population.** At production's 30_000 and `t_max = 50` the census lost 1569,
                // 1565 and 1566 pixels of 9216 on Heggie, logH and plain while AZ lost ZERO --
                // and `nonfin == budget` exactly on all three, so it was truncation and not
                // divergence. Three unrelated integrators losing near-identical counts is the
                // tell that it is the cap and not the method. If that excluded 17% is the wedge
                // region then Heggie's artefact would read as removed *because it was dropped*.
                Override::MaxSteps(2_000_000),
            ]);
            let px: Vec<PixelOut> =
                (0..sl.npix()).into_par_iter().map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
            let lg: Vec<f64> = px
                .iter()
                .map(|p| {
                    let d = p.energy_drift_max;
                    if d.is_finite() && d > 0.0 { d.log10() } else { f64::NAN }
                })
                .collect();
            let nonfin = lg.iter().filter(|x| !x.is_finite()).count();
            // **The exclusion has to be accounted for or the census is scored on a different
            // population per arm.** Non-finite pixels leave both numerator and denominator, so an
            // arm that loses 17% of the frame is being measured on the 83% that survived -- and if
            // the lost 17% IS the wedge region, the artefact would read as removed because it was
            // dropped. `budget` separates "the integrator diverged" from "the run was truncated",
            // which after the no-discard fix both read as non-finite.
            let budget = px.iter().filter(|p| p.budget_exhausted).count();
            let mut fin: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
            fin.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p50 = if fin.is_empty() { f64::NAN } else { fin[fin.len() / 2] };
            // The scale-free assertion, per row: shift every finite value by +6 decades.
            let shifted: Vec<f64> = lg.iter().map(|x| x + 6.0).collect();
            for thr in thrs {
                let (frac, sharp, _pairs) = sharp_pairs(&lg, res, thr);
                let (_f2, s2, _) = sharp_pairs(&shifted, res, thr);
                let hot = sharp_mask(&lg, res, thr);
                let total_hot = hot.iter().filter(|&&h| h).count();
                let largest = prin_rs::spatial::components(&hot, res)
                    .first().map(|c| c.len()).unwrap_or(0);
                let (str_lg, str_med, n_comp) = straight_stats(&hot, res);
                println!(
                    "{:>10} {:>11} {:>6.1} {:>9.5} {:>8} {:>10.4} {:>9.4} {:>8.4} {:>6} {:>10.3e} {:>3}",
                    if thr == thrs[0] { label } else { "" },
                    if thr == thrs[0] { format!("{nonfin}/{budget}") } else { String::new() },
                    thr,
                    frac,
                    largest,
                    largest as f64 / total_hot.max(1) as f64,
                    str_lg,
                    str_med,
                    n_comp,
                    if thr == thrs[0] { p50 } else { f64::NAN },
                    if s2 == sharp { "ok" } else { "BRK" }
                );
            }
        }
        println!();
    }

    println!("HOW TO READ THIS");
    println!();
    println!("**`scalefree` first.** If any row reads BROKEN the census is measuring magnitude");
    println!("and every number above it is void.");
    println!();
    println!("**Then the AZ row on `config_stability` against `preset_prho`.** If they do not");
    println!("separate, the metric has no subject and the integrator columns say nothing.");
    println!();
    println!("**`lgst/shrp` is not enough either, and `plain_rk4` is why.** It has NO reference");
    println!("body and therefore no wedge mechanism, and it still scores a large connected");
    println!("component -- a wrecked field is connected too. Density and connectedness cannot");
    println!("tell a decision boundary from damage.");
    println!();
    println!("**`str_lgst` and `str_med` are the discriminator.** 0.0 is a perfect line, 1.0 is");
    println!("isotropic; `spatial::straightness` is total least squares in closed form, and");
    println!("`tests/straightness.rs` shows it scores a line 0.000000 against a wandering curve");
    println!("0.367927 at MATCHED extent and count, and is rotation-invariant. A decision");
    println!("boundary is STRAIGHT; a fractal boundary and a damage blob are not.");
    println!();
    println!("**`ncomp` guards the median.** It counts components of >= 20 px; a short component");
    println!("is straight by construction, so a median over dust would read as a boundary.");
}
