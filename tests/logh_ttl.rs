//! **TTL — the time-transformed leapfrog, and the two things that must hold before any comparison.**
//!
//! Mikkola & Aarseth, CMDA **84** (2002) 343. logH's kick denominator `U = sum m_i m_j / r_ij` is
//! mass-weighted, so a close approach between a heavy body and a light one barely moves it and the
//! physical step fails to shrink. TTL swaps in `Omega = sum w_ij / r_ij` with free weights and
//! carries an auxiliary `W` advanced by `dW = (dOmega/dt) dt`.
//!
//! This port takes `w_ij = mbar^2`, which makes `Omega === U` **identically** at equal masses. So
//! the equal-mass control is an algebraic identity and not an approximation, and it is the first
//! test here: a TTL arm that differs on an equal-mass slice is measuring something other than the
//! mass ratio, and every number after it would be uninterpretable.

use prin_rs::integrate::logh::hamiltonian::{omega, omega_dot, ttl_weight, LhTime};
use prin_rs::integrate::logh::{integrate_lh, LhOpts, LhState, LhSystem, Stepper};
use prin_rs::physics::{energy, Cart};
use prin_rs::Vec2;

fn cart(r: [[f64; 2]; 3], v: [[f64; 2]; 3]) -> Cart<f64> {
    Cart {
        r: [Vec2::new(r[0][0], r[0][1]), Vec2::new(r[1][0], r[1][1]), Vec2::new(r[2][0], r[2][1])],
        v: [Vec2::new(v[0][0], v[0][1]), Vec2::new(v[1][0], v[1][1]), Vec2::new(v[2][0], v[2][1])],
    }
}

/// A COM-centred, non-degenerate triangle with unequal separations.
fn fixture() -> Cart<f64> {
    cart(
        [[1.0, 0.3], [-0.7, 0.9], [-0.3, -1.2]],
        [[0.05, -0.11], [-0.02, 0.07], [-0.03, 0.04]],
    )
}

/// **THE CONTROL, AND IT IS AN IDENTITY.** At `m_0 = m_1 = m_2` the weight `mbar^2` equals
/// `m_i m_j` for every pair, so `Omega` and `U` are the same expression. If this ever drifts, the
/// weight convention has changed and the equal-mass tie in the sweep stops meaning anything.
#[test]
fn omega_equals_the_potential_exactly_at_equal_masses() {
    let m = [0.4f64, 0.4, 0.4];
    let sys = LhSystem::new(m);
    let c = fixture();
    let u = energy::potential_pos(&c.r, &m, 0.0);
    let om = omega(&sys, &c.r);
    assert!(
        (om - u).abs() <= 1e-15 * u.abs(),
        "Omega {om:.17e} against U {u:.17e} at equal masses -- these must be the same expression"
    );
    assert!((ttl_weight(&sys) - m[0] * m[1]).abs() <= 1e-16);
}

/// **And the control has an arm that says it is not vacuous.** At UNEQUAL masses the two must
/// differ, or `Omega === U` always and TTL is logH under another name.
#[test]
fn omega_differs_from_the_potential_at_unequal_masses() {
    let m = [0.9f64, 0.09, 0.01];
    let sys = LhSystem::new(m);
    let c = fixture();
    let u = energy::potential_pos(&c.r, &m, 0.0);
    let om = omega(&sys, &c.r);
    let rel = (om - u).abs() / u.abs();
    assert!(rel > 0.1, "Omega and U differ by only {rel:.3e} at a 90:1 mass ratio");
}

/// `omega_dot` is the exact time derivative and it **advances the time transformation itself**,
/// so an error here is not a diagnostic error. Central-differenced against `omega` along the flow.
///
/// `h = 1e-3` and not `1e-8`: `Omega` is smooth and `O(1)` here, so the truncation term is tiny
/// while a smaller step would surrender digits to cancellation for nothing. *Ask what order the
/// function is in the variable before choosing the step* — a small step is not the safe choice.
#[test]
fn omega_dot_is_the_derivative_of_omega() {
    let m = [0.5f64, 0.3, 0.2];
    let sys = LhSystem::new(m);
    let c = fixture();
    let h = 1e-3;
    let adv = |dt: f64| {
        let mut r = c.r;
        for i in 0..3 {
            r[i] = r[i] + c.v[i] * dt;
        }
        omega(&sys, &r)
    };
    let fd = (adv(h) - adv(-h)) / (2.0 * h);
    let exact = omega_dot(&sys, &c.r, &c.v);
    let rel = (fd - exact).abs() / exact.abs().max(1e-30);
    assert!(rel < 1e-6, "omega_dot {exact:.9e} against FD {fd:.9e}, rel {rel:.3e}");

    // The mutation arm: a sign error is the plausible transcription failure and must be caught.
    let wrong = -exact;
    assert!(
        (fd - wrong).abs() / exact.abs() > 1.0,
        "a sign flip in omega_dot is not detected by this test"
    );
}

/// **`W` must track `Omega` along the march, and the test is CONVERGENCE, not a threshold.**
///
/// They are equal at registration and separate only by integration error. Because `W` is the
/// *drift denominator*, a `W` that wandered would still march and would quietly produce wrong
/// physics — nothing would fail, the trajectory would just be someone else's.
///
/// An absolute tolerance cannot distinguish "second-order error, working as designed" from "the
/// update is wrong": the first attempt at this asserted `< 1e-6`, measured `3.7e-6`, and would
/// have been passed by halving the step or loosened to keep green. The leapfrog is second order,
/// so the honest statement is that **halving the step quarters the gap**. A wrong `dW` term does
/// not converge at all, and a first-order one converges at 2x rather than 4x.
#[test]
fn the_gap_between_w_and_omega_converges_at_second_order() {
    let m = [0.5f64, 0.3, 0.2];
    let sys = LhSystem::new(m);
    let c = fixture();
    let b = sys.b_of(&c);
    // Same span in fictitious time at every rung, so the rungs are the same trajectory at
    // different resolutions and not different trajectories.
    let span = 2.0f64;
    let gap = |h: f64| {
        let n = (span / h).round() as usize;
        let mut s = LhState::from_cart(&c);
        s.w = omega(&sys, &s.r);
        for _ in 0..n {
            let (nx, _) = prin_rs::integrate::logh::step::kdk(&sys, &s, b, LhTime::Ttl, h);
            s = nx;
        }
        assert!(s.t > 0.0 && s.is_finite(), "NO SUBJECT: the march at h = {h:.1e} did not advance");
        let om = omega(&sys, &s.r);
        ((s.w - om).abs() / om.abs(), s.w)
    };
    let (g1, _) = gap(2e-3);
    let (g2, _) = gap(1e-3);
    let (g3, _) = gap(5e-4);
    let (r1, r2) = (g1 / g2, g2 / g3);
    println!("gap {g1:.3e} -> {g2:.3e} -> {g3:.3e}   ratios {r1:.2} {r2:.2}");
    assert!(g3 > 0.0, "the gap is exactly zero -- W is not being advanced at all");
    for (r, lbl) in [(r1, "first"), (r2, "second")] {
        assert!(
            r > 3.0 && r < 5.5,
            "{lbl} halving moved the W-Omega gap by {r:.2}x, not the ~4x a second-order \
             leapfrog requires -- 2x would mean the dW term is first order and ~1x that it is wrong"
        );
    }
}

/// **The equal-mass tie, end to end.** With `Omega === U`, TTL and logH integrate the same
/// equations, so a full march must agree to round-off. This is the control the mass-ratio sweep
/// rests on: if it fails, a TTL-vs-logH difference at unequal masses cannot be attributed to the
/// mass ratio.
#[test]
fn ttl_and_logh_agree_at_equal_masses_and_differ_at_a_large_ratio() {
    let c = fixture();
    let run = |m: [f64; 3], time: LhTime| {
        integrate_lh(
            c, &m, 4.0, 10, 1e-3, 200_000,
            &LhOpts {
                time,
                stepper: Stepper::Kdk,
                r_coll_frac: 0.0,
                stop_on_event: false,
                ..Default::default()
            },
        )
    };
    let sep = |a: &Cart<f64>, b: &Cart<f64>| {
        (0..3).map(|i| (a.r[i] - b.r[i]).norm()).fold(0.0f64, f64::max)
    };

    let eq = [0.4f64, 0.4, 0.4];
    let (a, b) = (run(eq, LhTime::LogH), run(eq, LhTime::Ttl));
    assert!(a.finite && b.finite, "equal-mass control did not complete");
    let d_eq = sep(&a.state, &b.state);
    assert!(d_eq < 1e-9, "equal-mass: TTL and logH differ by {d_eq:.3e}, must agree to round-off");

    // The arm that keeps the control from being vacuous.
    let un = [0.9f64, 0.09, 0.01];
    let (c1, c2) = (run(un, LhTime::LogH), run(un, LhTime::Ttl));
    assert!(c1.finite && c2.finite, "unequal-mass arms did not complete");
    let d_un = sep(&c1.state, &c2.state);
    assert!(
        d_un > 1e-6,
        "at a 90:1 mass ratio TTL and logH agree to {d_un:.3e} -- the mode is inert and the \
         equal-mass agreement above proves nothing"
    );
    println!("equal-mass |dr| {d_eq:.3e}   90:1 ratio |dr| {d_un:.3e}");
}
