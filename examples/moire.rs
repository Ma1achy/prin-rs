//! **Are the striations in the ribbons a stepping artefact, a cadence artefact, or physics?**
//!
//! The `_uniform` panels carry fine parallel banding inside the otherwise solid ribbons. Three
//! candidates, and they are separable by a causal test rather than by inspection:
//!
//! 1. **Sync-cadence quantisation.** This project's own standing result: escape is sampled only
//!    at sync boundaries, so `t_end` takes `n_sync` values wherever escape terminates and *every
//!    derived field draws those steps*. Measured once at `escape_every` 0 -> 1: `preset_plambda`
//!    went **16 -> 2623** distinct `t_end` with 99.52% -> 0.26% landing exactly on a boundary.
//! 2. **Step-size resonance.** Rauch & Holman 1999 (AJ 117, 1087) on the Wisdom-Holman mapping:
//!    artificial chaos arises from the **overlap of step-size resonances**, cured when the step
//!    resolves pericentre -- they quote ~`T_p/20` with `T_p = 2*pi/f_dot` at pericentre. This
//!    port's `StepLimit::Predictive` at `f = 0.02` is that condition at 1/50, so if the banding
//!    is a step resonance it should already be suppressed and should move with `eta`.
//! 3. **Physics.** Real fine structure in initial-condition space.
//!
//! # The separation, and the confound it has to avoid
//!
//! `dtau = eta*dt_left/(A0*B0)`, so **changing `n_sync` at fixed `eta` also changes the step
//! size** -- the standing rule that `n_sync` fixed while `t_max` varies compares different
//! discretisations, in its other direction. To vary the *cadence* alone, `eta` is scaled with
//! `n_sync`, exactly as `az_machinery`'s CONTROLLED rows do. `steps p50` is printed as the check
//! that it worked.
//!
//! ```text
//!   cadence  : n_sync x2, eta x2   -> step size HELD. Bands move => the cadence draws them.
//!   step     : eta /2 at fixed n_sync -> cadence HELD. Bands move => the stepper draws them.
//!   neither moves => physics.
//! ```
//!
//! # Reading the banding
//!
//! Per row of the window, the lightness is detrended and its power computed at each integer
//! wavelength from 2 to 64 px; the rows are averaged and the peak reported with its **prominence**
//! (peak over median power). **A peak with prominence near 1 is not a band**, it is the largest
//! bin of a flat spectrum -- printed so a wavelength is never quoted off a featureless spectrum.
//!
//! ```text
//! cargo run --release --example moire -- [res] [out_dir]
//! ```
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

/// Peak spatial wavelength in px and its prominence over the median bin.
fn banding(v: &[f64], n: usize) -> (usize, f64) {
    let mut power = vec![0.0f64; 65];
    let mut rows = 0usize;
    for y in 0..n {
        let row: Vec<f64> = (0..n).map(|x| v[y * n + x]).filter(|x| x.is_finite()).collect();
        if row.len() < n {
            continue; // a row with a hole would alias; skip it rather than pad
        }
        // Detrend: remove the linear fit, so a smooth gradient is not read as a long wave.
        let m = row.iter().sum::<f64>() / n as f64;
        let mid = (n as f64 - 1.0) / 2.0;
        let (mut sxy, mut sxx) = (0.0f64, 0.0f64);
        for (i, &r) in row.iter().enumerate() {
            let dx = i as f64 - mid;
            sxy += dx * (r - m);
            sxx += dx * dx;
        }
        let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
        let d: Vec<f64> = row.iter().enumerate().map(|(i, &r)| r - m - slope * (i as f64 - mid)).collect();
        for lam in 2..=64usize {
            let w = std::f64::consts::TAU / lam as f64;
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (i, &x) in d.iter().enumerate() {
                re += x * (w * i as f64).cos();
                im += x * (w * i as f64).sin();
            }
            power[lam] += (re * re + im * im) / (n * n) as f64;
        }
        rows += 1;
    }
    if rows == 0 {
        return (0, f64::NAN);
    }
    for p in power.iter_mut() {
        *p /= rows as f64;
    }
    let mut sorted: Vec<f64> = power[2..=64].to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = sorted[sorted.len() / 2];
    let (mut best, mut bp) = (0usize, 0.0f64);
    for lam in 2..=64usize {
        if power[lam] > bp {
            bp = power[lam];
            best = lam;
        }
    }
    (best, if med > 0.0 { bp / med } else { f64::NAN })
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/moire".into());
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    // A window inside the striated red band, upper-left of the full view: a quarter of the
    // half-width, offset so it sits in the ribbon rather than across its edge.
    let (fcx, fcy, fhalf) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let half = fhalf * 0.18;
    let (cx, cy) = (fcx - fhalf * 0.45, fcy + fhalf * 0.42);
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);

    const N0: usize = 125;
    const ETA0: f64 = 1e-2;
    let arms: Vec<(&str, usize, f64)> = vec![
        ("baseline", N0, ETA0),
        ("cadence x2", 2 * N0, 2.0 * ETA0),   // step size HELD
        ("cadence x4", 4 * N0, 4.0 * ETA0),   // step size HELD
        ("step /2", N0, ETA0 / 2.0),          // cadence HELD
        ("step /4", N0, ETA0 / 4.0),          // cadence HELD
    ];

    println!("MOIRE / BANDING on a ribbon window of config_stability, {res}^2, t_max = 50.");
    println!("half = {half:.5} (0.18 of the full view), centre ({cx:.5}, {cy:.5}).");
    println!();
    println!("  `cadence` scales eta WITH n_sync so the step size is held -- `steps p50` is the");
    println!("  check. `step` holds n_sync and shrinks eta. Bands moving under the first means");
    println!("  the sync cadence draws them; under the second, the stepper; under neither,");
    println!("  physics.");
    println!();
    println!(
        "{:>12} {:>7} {:>9} {:>10} {:>8} {:>9} {:>9} {:>9} {:>8}",
        "arm", "n_sync", "eta", "steps p50", "lambda", "promin", "t_end dst", "on bnd", "nonfin"
    );

    let mut window: Option<(f64, f64)> = None;
    for (label, n_sync, eta) in arms {
        let cfg = EnsembleCfg::production().with_overrides(&[
            Override::TMax(50.0),
            Override::NSync(n_sync),
            Override::Eta(eta),
            Override::RefineFlagged(false),
            Override::MaxSteps(4_000_000),
        ]);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let (lo, hi) = *window.get_or_insert_with(|| colour::range(&px, Scalar::ShapeSpread));

        let mut rgb = Vec::with_capacity(px.len() * 3);
        for p in &px {
            rgb.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
        }
        let stem = label.replace([' ', '/'], "_");
        let _ = adaptive::save_rect(&format!("{dir}/{stem}.png"), res, res, &rgb);

        // Band the LIGHTNESS of what is displayed, which is what the eye reads.
        let lum: Vec<f64> = (0..px.len())
            .map(|i| {
                0.2126 * rgb[3 * i] as f64 + 0.7152 * rgb[3 * i + 1] as f64
                    + 0.0722 * rgb[3 * i + 2] as f64
            })
            .collect();
        let (lam, prom) = banding(&lum, res);

        let dt = 50.0 / n_sync as f64;
        let te: Vec<f64> = px.iter().map(|p| p.t_end).filter(|x| x.is_finite()).collect();
        let mut d: Vec<u64> = te.iter().map(|t| (t / dt * 1e6).round() as u64).collect();
        d.sort_unstable();
        d.dedup();
        let on_b = te.iter().filter(|t| ((*t / dt) - (*t / dt).round()).abs() < 1e-9).count();
        let mut steps: Vec<u64> = px.iter().map(|p| p.total_substeps).collect();
        steps.sort_unstable();

        println!(
            "{:>12} {:>7} {:>9.2e} {:>10.3e} {:>8} {:>9.2} {:>9} {:>9.4} {:>8}",
            label,
            n_sync,
            eta,
            steps[steps.len() / 2] as f64,
            lam,
            prom,
            d.len(),
            on_b as f64 / te.len().max(1) as f64,
            px.iter().filter(|p| !p.spread_shape.is_finite()).count()
        );
    }

    println!();
    println!("HOW TO READ THIS");
    println!();
    println!("**`promin` first.** A peak with prominence near 1 is the largest bin of a FLAT");
    println!("spectrum, not a band. `lambda` means nothing without it.");
    println!();
    println!("**Then `steps p50` across the two `cadence` rows.** If it is not roughly flat the");
    println!("step size was not held and those rows measure the stepper, not the cadence --");
    println!("`dtau = eta*dt_left/(A0*B0)` is why scaling eta with n_sync is required at all.");
    println!();
    println!("**`t_end dst` and `on bnd` are the standing cadence diagnostic.** Escape is sampled");
    println!("only at sync boundaries, so where it terminates `t_end` takes n_sync values and");
    println!("every derived field draws those steps. Read the DELTA across cadences, never the");
    println!("level: `near-field` reads 97.85% on-boundary while being completely clean, because");
    println!("its footprints sit at t_end = t_max, and the horizon is a boundary time.");
}
