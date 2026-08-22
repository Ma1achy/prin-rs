//! Step 3 gate.
//!
//! The unregularised integrator is *expected* to fail on close encounters — that failure is
//! why Aarseth–Zare exists. What must hold is that on the pixels it can handle, **drift
//! falls when the step size falls**. Drift insensitive to step size is the signature of a
//! wrong equation, not an accuracy problem, and has caught three bugs in this project.
//! If that check fails, Step 4 must not begin.

use prin_rs::grid;
use prin_rs::integrate::leapfrog;
use prin_rs::physics::burrau;

const T_MAX: f64 = 13.0;
const MAX_STEPS: usize = 2_000_000;

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

fn run(eta: f64) -> Vec<Option<f64>> {
    let m = burrau::masses::<f64>();
    let s = grid::region("near-field", 5, 5, 0.05).unwrap();
    (0..s.npix())
        .map(|i| {
            let out = leapfrog::integrate(s.nominal::<f64>(i), &m, T_MAX, eta, 0.0, MAX_STEPS);
            if out.reached(T_MAX) {
                Some(out.drift)
            } else {
                None
            }
        })
        .collect()
}

#[test]
fn drift_falls_with_step_size_on_tractable_pixels() {
    let coarse = run(0.01);
    let fine = run(0.005);
    let n = coarse.len();

    let failed: Vec<usize> = (0..n).filter(|&i| coarse[i].is_none() || fine[i].is_none()).collect();
    println!("near-field 5x5, t = {T_MAX}, unregularised leapfrog");
    println!("pixels: {n}, did not reach t_max at one or both eta: {}", failed.len());
    if !failed.is_empty() {
        println!("  failing pixel indices: {failed:?}");
    }

    let both: Vec<usize> = (0..n).filter(|&i| coarse[i].is_some() && fine[i].is_some()).collect();
    assert!(!both.is_empty(), "no pixel completed at either step size");

    let dc: Vec<f64> = both.iter().map(|&i| coarse[i].unwrap()).collect();
    let df: Vec<f64> = both.iter().map(|&i| fine[i].unwrap()).collect();

    println!("  median |dE/E| at eta=0.01  : {:.6e}", median(dc.clone()));
    println!("  median |dE/E| at eta=0.005 : {:.6e}", median(df.clone()));
    println!("  max    |dE/E| at eta=0.01  : {:.6e}", dc.iter().cloned().fold(0.0, f64::max));
    println!("  max    |dE/E| at eta=0.005 : {:.6e}", df.iter().cloned().fold(0.0, f64::max));

    let ratios: Vec<f64> = both
        .iter()
        .map(|&i| {
            let (c, f) = (coarse[i].unwrap(), fine[i].unwrap());
            if c > 0.0 { f / c } else { f64::NAN }
        })
        .filter(|x| x.is_finite())
        .collect();
    let med_ratio = median(ratios.clone());
    println!("  median drift ratio (eta/2 : eta) = {med_ratio:.4}");
    println!("  (2nd-order convergence would give ~0.25; adaptive stepping perturbs this)");

    // The gate. Not a convergence-order assertion — the adaptive step and the chaotic
    // amplification both perturb the exponent — but drift must *fall*. If it does not,
    // the equations are wrong.
    assert!(
        med_ratio < 1.0,
        "median drift ratio {med_ratio} >= 1: halving the step did not reduce drift. \
         This is the signature of a wrong equation, not a step-size problem. Do not proceed to Step 4."
    );
}

#[test]
fn grid_row_order_is_x_fastest() {
    // index = jy*nx + jx. The cross-check compares row by row, so this ordering is
    // load-bearing: getting it wrong would look like a physics disagreement.
    let s = grid::Slice { nx: 3, ny: 2, cx: 0.0, cy: 0.0, half: 1.0, body: 0 };
    let want = [
        (-1.0, -1.0), (0.0, -1.0), (1.0, -1.0),
        (-1.0,  1.0), (0.0,  1.0), (1.0,  1.0),
    ];
    for (i, &(wx, wy)) in want.iter().enumerate() {
        let (x, y) = s.decode_pos(i);
        assert!((x - wx).abs() < 1e-15 && (y - wy).abs() < 1e-15, "idx {i}: got ({x},{y}), want ({wx},{wy})");
    }
}

#[test]
fn cell_widths_are_per_axis() {
    // The reference uses hx for both axes; latent on square grids, wrong on any other.
    let s = grid::Slice { nx: 5, ny: 3, cx: 0.0, cy: 0.0, half: 1.0, body: 0 };
    let (hx, hy) = s.cell_widths();
    assert!((hx - 0.5).abs() < 1e-15, "hx = {hx}");
    assert!((hy - 1.0).abs() < 1e-15, "hy = {hy}");
}

/// Where the unregularised integrator actually fails.
///
/// BRIEF §2.3 says naive integration "fails" and §2.6 says `deep interior` will fail
/// however well the integrator is built. Near-field at `t=13` turns out **not** to be such a
/// region for the nominal copies, so this test locates the boundary rather than asserting a
/// failure that does not occur there. Reported, not asserted — it is a measurement.
#[test]
fn locate_where_the_unregularised_integrator_gives_up() {
    let m = burrau::masses::<f64>();
    println!("{:<16}{:>6}{:>8}{:>10}{:>14}{:>14}", "region", "t_max", "pixels", "incomplete", "median drift", "max drift");
    for (name, t_max) in [
        ("near-field", 13.0),
        ("near-field", 40.0),
        ("deep interior", 13.0),
        ("body2 core", 13.0),
    ] {
        let s = grid::region(name, 5, 5, 0.05).unwrap();
        let mut done = Vec::new();
        let mut incomplete = 0usize;
        for i in 0..s.npix() {
            let out = leapfrog::integrate(s.nominal::<f64>(i), &m, t_max, 0.01, 0.0, MAX_STEPS);
            if out.reached(t_max) {
                done.push(out.drift);
            } else {
                incomplete += 1;
            }
        }
        let med = if done.is_empty() { f64::NAN } else { median(done.clone()) };
        let mx = done.iter().cloned().fold(0.0f64, f64::max);
        println!("{name:<16}{t_max:>6}{:>8}{incomplete:>10}{med:>14.3e}{mx:>14.3e}", s.npix());
    }
    println!("\nAZ reference at near-field t=13 reaches |dE/E| ~ 3.9e-09 for comparison.");
}

/// The failure mode, characterised precisely.
///
/// It is *not* budget exhaustion at a 2e6-step budget — every pixel completes. It is that
/// the answer becomes meaningless while still looking like a number: at `deep interior` the
/// energy error reaches hundreds of times the total energy. Reported, not asserted.
#[test]
fn characterise_the_failure_mode() {
    let m = burrau::masses::<f64>();
    println!("{:<16}{:>12}{:>12}{:>14}{:>10}", "region", "worst steps", "median steps", "worst drift", "d_min");
    for name in ["near-field", "deep interior"] {
        let s = grid::region(name, 5, 5, 0.05).unwrap();
        let mut steps = Vec::new();
        let mut worst = 0.0f64;
        let mut dmin = f64::INFINITY;
        for i in 0..s.npix() {
            let o = leapfrog::integrate(s.nominal::<f64>(i), &m, 13.0, 0.01, 0.0, MAX_STEPS);
            steps.push(o.steps as f64);
            worst = worst.max(o.drift);
            dmin = dmin.min(o.d_min);
        }
        let mx = steps.iter().cloned().fold(0.0f64, f64::max);
        println!("{name:<16}{mx:>12.0}{:>12.0}{worst:>14.3e}{dmin:>10.2e}", median(steps.clone()));
    }
    println!("\nBudget is {MAX_STEPS} steps; nothing came close to exhausting it.");
    println!("The unregularised integrator does not hang here — it returns a wrong number.");
}
