//! Link 4: finite-difference `Gamma` against the analytic `deriv`.
//!
//! This is the test CLAUDE.md mandates — the one that caught two sign errors in the
//! reference that were otherwise invisible. It only means anything because `gamma` is
//! independently anchored by `tests/az_identities.rs`; without that, a sign error present in
//! both functions would pass silently.
//!
//! Two assertions, not one:
//!   - agreement at a small step, and
//!   - **the error falls as h^2**. An FD error insensitive to `h` is a wrong derivative, not
//!     a truncation problem. Same diagnostic logic as the rest of this project, applied one
//!     level up.

use prin_rs::integrate::az::hamiltonian::{deriv, gamma};
use prin_rs::integrate::az::lc;
use prin_rs::integrate::az::reference_body::triple;
use prin_rs::integrate::az::{AzState, AzSystem};
use prin_rs::physics::burrau;
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

const N: usize = 64;

/// `dGamma/ds[k]` for the eight non-time components, from the analytic derivative.
///
/// The state array is `[u1x,u1y,p1x,p1y,u2x,u2y,p2x,p2y,t]`, and `deriv` returns
/// `du = +dGamma/dp`, `dp = -dGamma/du`. So the signs alternate by block.
fn analytic_grad(sys: &AzSystem<f64>, s: &AzState<f64>, e: f64) -> [f64; 8] {
    let d = deriv(sys, s, e);
    [
        -d.p1.x, -d.p1.y, // dGamma/du1
        d.u1.x, d.u1.y,   // dGamma/dp1
        -d.p2.x, -d.p2.y, // dGamma/du2
        d.u2.x, d.u2.y,   // dGamma/dp2
    ]
}

/// The step is scaled per component, `h_k = h * max(|s_k|, 1)`. A fixed absolute step is
/// wrong here: the states span orders of magnitude in `|u|` and `|p|`, so one absolute `h`
/// is simultaneously too coarse for the small components and too fine for the large ones.
fn fd_grad<F>(g: F, s: &AzState<f64>, h: f64) -> [f64; 8]
where
    F: Fn(&AzState<f64>) -> f64,
{
    let base = s.to_array9();
    let mut out = [0.0; 8];
    for k in 0..8 {
        let hk = h * base[k].abs().max(1.0);
        let mut hi = base;
        let mut lo = base;
        hi[k] += hk;
        lo[k] -= hk;
        out[k] = (g(&AzState::from_array9(hi)) - g(&AzState::from_array9(lo))) / (2.0 * hk);
    }
    out
}

/// States spanning several orders of magnitude in `|u|` and `|p|`, kept away from the LC
/// branch cut so the measurement is of the derivative and not of `u_of_rho`'s conditioning.
fn random_states(n: usize, seed: u64) -> Vec<(AzState<f64>, f64)> {
    let mut rng = SplitMix64::new(seed);
    (0..n)
        .map(|_| {
            let su = 10f64.powf(rng.range(-1.0, 1.0));
            let sp = 10f64.powf(rng.range(-1.0, 1.0));
            (
                AzState {
                    u1: Vec2::new(rng.range(0.2, 2.0) * su, rng.range(-2.0, 2.0) * su),
                    p1: Vec2::new(rng.range(-2.0, 2.0) * sp, rng.range(-2.0, 2.0) * sp),
                    u2: Vec2::new(rng.range(0.2, 2.0) * su, rng.range(-2.0, 2.0) * su),
                    p2: Vec2::new(rng.range(-2.0, 2.0) * sp, rng.range(-2.0, 2.0) * sp),
                    t: 0.0,
                },
                rng.range(-30.0, 5.0),
            )
        })
        .collect()
}

fn systems() -> Vec<AzSystem<f64>> {
    let m = burrau::masses::<f64>();
    (0..3).map(|a| { let (x, y, z) = triple(a); AzSystem::new(x, y, z, m) }).collect()
}

/// Error normalised by the **magnitude of the gradient**, not by the individual component.
///
/// Per-component relative error is the wrong measure here: a gradient component passing
/// through zero makes it blow up while the absolute discrepancy is negligible. The object
/// under test is the gradient as a vector, so that is what sets the scale.
fn worst_rel(sys: &AzSystem<f64>, s: &AzState<f64>, e: f64, h: f64) -> f64 {
    let an = analytic_grad(sys, s, e);
    let fd = fd_grad(|st| gamma(sys, st, e), s, h);
    let scale = (0..8)
        .map(|k| an[k].abs().max(fd[k].abs()))
        .fold(0.0f64, f64::max)
        .max(1e-300);
    (0..8).map(|k| (an[k] - fd[k]).abs() / scale).fold(0.0, f64::max)
}

#[test]
fn deriv_matches_finite_differenced_gamma() {
    let mut worst = 0.0f64;
    let mut worst_at = String::new();
    let mut all = Vec::new();
    for sys in systems() {
        for (s, e) in random_states(N, 0x5160) {
            let w = worst_rel(&sys, &s, e, 1e-5);
            all.push(w);
            if w > worst {
                worst = w;
                let r3 = (lc::rho_of_u(s.u2) - lc::rho_of_u(s.u1)).norm();
                worst_at = format!("a={} A={:.2e} B={:.2e} |R3|={:.2e}", sys.a, s.a(), s.b(), r3);
            }
        }
    }
    all.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let med = all[all.len() / 2];
    let p99 = all[(all.len() * 99) / 100];
    println!("FD vs analytic, h = 1e-5 * max(|s_k|,1), {} states:", all.len());
    println!("  median = {med:.3e}   p99 = {p99:.3e}   worst = {worst:.3e}");
    println!("  worst at {worst_at}");
    println!("  worst cases sit at small |R3|, where the unregularised 1/|R3| term has a");
    println!("  large third derivative and central differencing truncates harder. That is a");
    println!("  property of the geometry, not of the derivative — the h^2 test below is the");
    println!("  correctness statement; this one bounds the conditioning.");
    assert!(med < 1e-9, "median = {med:e} — the typical case should be far tighter than this");
    assert!(worst < 1e-6, "worst = {worst:e} at {worst_at}");
}

#[test]
fn finite_difference_error_falls_as_h_squared() {
    // Deliberately in the truncation-dominated regime. At h = 1e-5 roundoff and truncation
    // are comparable and the exponent is muddied; at 1e-3 truncation dominates cleanly.
    let (h1, h2) = (1e-3, 5e-4);
    let mut ratios = Vec::new();
    for sys in systems() {
        for (s, e) in random_states(N, 0x1234) {
            let (a, b) = (worst_rel(&sys, &s, e, h1), worst_rel(&sys, &s, e, h2));
            if a > 1e-12 {
                ratios.push(b / a);
            }
        }
    }
    ratios.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let med = ratios[ratios.len() / 2];
    println!("FD error ratio (h/2 : h) over {} states: median = {med:.4}", ratios.len());
    println!("  second-order central differencing predicts 0.25");
    assert!(
        med < 0.4,
        "FD error ratio {med} does not fall as h^2. An error insensitive to h means a wrong \
         derivative, not truncation."
    );
}

/// The FD test must have teeth. Flip the sign of the `R3` sub-term in `g2` — the exact
/// hazard flagged in CLAUDE.md, and one with no visible symptom in a trajectory — and the
/// test has to catch it.
#[test]
fn the_fd_test_detects_the_r3_sign_error() {
    let sys = systems()[0];
    let mut worst_correct = 0.0f64;
    let mut worst_flipped = 0.0f64;

    for (s, e) in random_states(N, 0x9999) {
        let fd = fd_grad(|st| gamma(&sys, st, e), &s, 1e-5);
        let an = analytic_grad(&sys, &s, e);

        // Reproduce deriv's g2, then flip the second R3 sub-piece from - to +.
        let r1 = lc::rho_of_u(s.u1);
        let r2 = lc::rho_of_u(s.u2);
        let r3v = r2 - r1;
        let r3 = r3v.norm().max(f64::MIN_POSITIVE);
        let corr = lc::lt_apply(s.u2, r3v) * (2.0 * s.a() * s.b() / (r3 * r3 * r3))
            * (sys.mb * sys.mc)
            * 2.0;
        // an[4], an[5] are dGamma/du2 = -deriv.p2
        let flipped = [an[4] + corr.x, an[5] + corr.y];

        let scale = (0..8)
            .map(|k| an[k].abs().max(fd[k].abs()))
            .fold(0.0f64, f64::max)
            .max(1e-300);
        for (i, k) in [4usize, 5].iter().enumerate() {
            worst_correct = worst_correct.max((an[*k] - fd[*k]).abs() / scale);
            worst_flipped = worst_flipped.max((flipped[i] - fd[*k]).abs() / scale);
        }
    }
    println!("dGamma/du2 — correct sign: {worst_correct:.3e}, flipped sign: {worst_flipped:.3e}");
    assert!(worst_correct < 1e-7, "the correct derivative did not agree: {worst_correct:e}");
    assert!(
        worst_flipped > 1e-3,
        "a flipped R3 sign still agreed with FD at {worst_flipped:e} — the test has no teeth"
    );
}
