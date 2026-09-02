//! **IAS15 — is the port right, and does it earn the name "reference arm"?**
//!
//! Rein & Spiegel, MNRAS **446** (2015) 1424. The claim is machine precision with an error that
//! random-walks rather than drifting, so the tests are about *accuracy*, not about order: a
//! fifteenth-order method reaches round-off within one or two step halvings and an order fit would
//! be measuring the floor.
//!
//! It is a **reference arm and not a production candidate** — its predictor-corrector iterates a
//! variable number of times per step, which is the per-lane variable work already measured as
//! fatal on GPU here. The iteration count is asserted to be *variable*, because that is the
//! property, not a defect to hide.

use prin_rs::integrate::ias15::{conversion, next_dt, step, H, N};
use prin_rs::physics::{energy, Cart};
use prin_rs::{Real, Vec2};

/// **The conversion matrix must satisfy its DEFINING IDENTITY, not merely be non-empty.**
///
/// `sum_j c[k][j] t^j == prod_{i<k} (t - h_i)`. Checked at random points rather than against a
/// transcribed table: a table would test that two copies of the same digits agree, which is not
/// the property that matters.
#[test]
fn the_conversion_matrix_satisfies_its_defining_identity() {
    let c = conversion();
    let mut worst = 0.0f64;
    for s in 0..64 {
        let t = -1.5 + 3.0 * (s as f64) / 63.0;
        for k in 0..N {
            let poly: f64 = (0..N).map(|j| c[k][j] * t.powi(j as i32 + 1)).sum();
            let prod: f64 = (0..=k).map(|i| t - H[i]).product();
            worst = worst.max((poly - prod).abs());
        }
    }
    println!("conversion identity worst |poly - prod| = {worst:.3e}");
    assert!(worst < 1e-12, "the g->b matrix does not expand the Newton basis: {worst:.3e}");

    // Negative control: a perturbed matrix must NOT satisfy it, or the check is vacuous.
    let mut bad = c;
    bad[3][2] += 1e-3;
    let t = 0.7f64;
    let poly: f64 = (0..N).map(|j| bad[3][j] * t.powi(j as i32 + 1)).sum();
    let prod: f64 = (0..=3).map(|i| t - H[i]).product();
    assert!((poly - prod).abs() > 1e-6, "a perturbed matrix still satisfies the identity");
}

/// The nodes must be strictly increasing on `(0, 1]` with `h_0 = 0`. A transposed or duplicated
/// node would make a divided difference divide by zero, which is worth catching by construction.
#[test]
fn the_gauss_radau_nodes_are_ordered_and_distinct() {
    assert_eq!(H[0], 0.0);
    for k in 1..8 {
        assert!(H[k] > H[k - 1], "node {k} is not greater than its predecessor");
        assert!(H[k] <= 1.0, "node {k} lies outside (0, 1]");
    }
}

fn two_body<T: Real>() -> (Cart<T>, [T; 3]) {
    // A bound eccentric pair plus a distant near-massless third body, so `accel` is exercised in
    // full while the dynamics stay close to an exactly-known two-body problem.
    let f = |x: f64| T::lit(x);
    (
        Cart {
            r: [
                Vec2::new(f(-0.5), f(0.0)),
                Vec2::new(f(0.5), f(0.0)),
                Vec2::new(f(40.0), f(0.0)),
            ],
            v: [
                Vec2::new(f(0.0), f(-0.6)),
                Vec2::new(f(0.0), f(0.6)),
                Vec2::new(f(0.0), f(0.02)),
            ],
        },
        [f(1.0), f(1.0), f(1e-8)],
    )
}

/// **The headline property: energy conserved at machine precision over many orbits.** This is what
/// makes it usable as ground truth, and it is the reason `eta/256` could not be — that reference
/// came back saturated, scoring a correct mode and a broken one alike.
#[test]
fn energy_is_conserved_at_machine_precision_over_many_orbits() {
    let (c, m) = two_body::<f64>();
    let e0 = energy::energy(&c.r, &c.v, &m, 0.0);
    let (mut r, mut v) = (c.r, c.v);
    let mut dt = 0.05f64;
    let mut evals = 0usize;
    let mut iters_seen: Vec<usize> = Vec::new();
    let mut t = 0.0f64;
    while t < 200.0 {
        let (out, e, it, b_last) = step(&r, &v, &m, dt, 12, 1e-14);
        assert!(out.r.iter().all(|q| q.is_finite()), "went non-finite at t = {t}");
        r = out.r;
        v = out.v;
        t += dt;
        evals += e;
        iters_seen.push(it);
        dt = next_dt(&r, &m, dt, b_last, 1e-9).min(0.5);
    }
    let e1 = energy::energy(&r, &v, &m, 0.0);
    let drift = ((e1 - e0) / e0).abs();
    let mean_it = iters_seen.iter().sum::<usize>() as f64 / iters_seen.len() as f64;
    println!(
        "IAS15 over t = 200: drift {drift:.3e}, {} steps, {evals} evals, mean {mean_it:.2} iters",
        iters_seen.len()
    );
    assert!(drift < 1e-12, "IAS15 energy drift {drift:.3e} is not machine-precision");
}

/// **The iteration count must VARY.** That is the property that makes this a reference arm rather
/// than a production candidate — per-lane variable work is what gave reject-and-retry
/// `warps hit 1.0000` here. A constant count would mean the corrector is not converging on its own
/// terms, and the GPU objection would not apply as stated.
#[test]
fn the_corrector_iteration_count_is_variable() {
    let (c, m) = two_body::<f64>();
    let (mut r, mut v) = (c.r, c.v);
    let mut counts: Vec<usize> = Vec::new();
    for _ in 0..200 {
        let (out, _, it, _) = step(&r, &v, &m, 0.05, 12, 1e-14);
        r = out.r;
        v = out.v;
        counts.push(it);
    }
    let lo = *counts.iter().min().unwrap();
    let hi = *counts.iter().max().unwrap();
    println!("corrector iterations over 200 steps: min {lo}, max {hi}");
    assert!(hi > lo, "the corrector always ran {lo} iterations -- the count is not state-dependent");
}

/// A halved step must not make the answer worse. Not an order fit: at fifteenth order the method
/// reaches round-off within a rung or two and a fitted slope would read the floor, which is the
/// failure this project has now hit on the figure-eight three separate times.
#[test]
fn halving_the_step_does_not_degrade_the_answer() {
    let (c, m) = two_body::<f64>();
    let e0 = energy::energy(&c.r, &c.v, &m, 0.0);
    let run = |dt: f64, n: usize| {
        let (mut r, mut v) = (c.r, c.v);
        for _ in 0..n {
            let (out, _, _, _) = step(&r, &v, &m, dt, 12, 1e-14);
            r = out.r;
            v = out.v;
        }
        let e1 = energy::energy(&r, &v, &m, 0.0);
        ((e1 - e0) / e0).abs()
    };
    let coarse = run(0.1, 200);
    let fine = run(0.05, 400);
    println!("drift at dt = 0.1: {coarse:.3e}   at dt = 0.05: {fine:.3e}");
    assert!(coarse.is_finite() && fine.is_finite());
    assert!(
        fine <= coarse * 10.0,
        "halving the step made the drift WORSE by more than a decade ({coarse:.3e} -> {fine:.3e})"
    );
}
