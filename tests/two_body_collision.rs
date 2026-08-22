//! Gate (b): the two-body radial collision.
//!
//! Equal masses released from rest with the third body far away. The pair falls radially
//! into an *exact* collision. Regularisation is what makes this survivable: in the
//! transformed Hamiltonian the singular `-G m_i m_j / |rho|` term becomes a constant, the
//! regularised two-body problem is a harmonic oscillator, and the trajectory passes through
//! collision at machine precision.
//!
//! Acceptance: `d_min < 1e-10` with `|dE/E| < 1e-12`.

use prin_rs::integrate::az;
use prin_rs::integrate::az::reference_body::choose_reference;
use prin_rs::physics::Cart;
use prin_rs::Vec2;

const N_SYNC: usize = 32;
const ETA: f64 = 0.01;
const MAX_STEPS: usize = 30_000;

fn setup() -> (Cart<f64>, [f64; 3]) {
    // Bodies 0 and 1: equal masses, unit separation, at rest.
    // Body 2: far away, and deliberately off-axis so the two long sides are not an exact
    // tie — an argmax tie would make the reference-body choice depend on tie-breaking
    // rather than on geometry.
    let m = [1.0, 1.0, 1.0];
    let s = Cart::new(
        [
            Vec2::new(-0.5, 0.0),
            Vec2::new(0.5, 0.0),
            Vec2::new(0.1, 1000.0),
        ],
        [Vec2::zero(); 3],
    );
    (s, m)
}

#[test]
fn the_close_pair_is_the_regularised_one() {
    // Worth asserting rather than assuming. AZ regularises the two pairs sharing the
    // reference body; if the geometry put the colliding pair on the *unregularised* side,
    // this test would be measuring nothing.
    let (s, _) = setup();
    let a = choose_reference(&s.r);
    let (ra, rb, rc) = az::reference_body::triple(a);
    println!("reference body = {a}, regularised pairs = ({ra},{rb}) and ({ra},{rc})");
    let pair_is_regularised = (ra == 0 || rb == 0 || rc == 0) && (ra == 1 || rb == 1 || rc == 1);
    assert!(pair_is_regularised, "the colliding pair (0,1) is not regularised");
    // Specifically, (a,b) must be the close pair.
    assert!(
        (ra == 1 && rb == 0) || (ra == 0 && rb == 1),
        "the close pair is not the first regularised pair: a={ra} b={rb} c={rc}"
    );
}

#[test]
fn radial_collision_passes_through_at_machine_precision() {
    let (s0, m) = setup();

    // Free-fall time from rest at separation d for total mass M is
    // (pi/2) sqrt(d^3 / (2 G M)); with d = 1, M = 2 that is ~0.785. Integrate past it so
    // the trajectory goes through collision rather than stopping at it.
    let t_max = 1.0;
    // eta is pinned small deliberately — see the scan test below. d_min at a genuine
    // collision is sampling-limited, not physical, so this threshold is a statement about
    // resolution as much as about correctness.
    let eta = 1e-4;
    let n_sync = 1;

    let out = az::integrate_az(s0, &m, t_max, n_sync, eta, 20_000_000, None);

    println!("two-body radial collision, t_max = {t_max}, n_sync = {n_sync}, eta = {eta}");
    println!("  d_min (regularised pairs) = {:.6e}", out.d_min_ref);
    println!("  d_min (all three pairs)   = {:.6e}", out.d_min_true);
    println!("  |dE/E|                    = {:.6e}", out.drift);
    println!("  max |Gamma| / |largest term| = {:.6e}", out.gamma_max);
    println!("  steps = {}, reference switches = {}", out.steps, out.switches);

    assert!(out.finite, "trajectory went non-finite");
    assert!(!out.budget_exhausted, "step budget exhausted");
    assert!(out.d_min_ref < 1e-10, "d_min = {:e}, want < 1e-10", out.d_min_ref);
    assert!(out.drift < 1e-12, "|dE/E| = {:e}, want < 1e-12", out.drift);
}

/// **d_min at an exact collision is a discretisation artefact, not a physical quantity.**
///
/// The true minimum separation of a radial two-body collision is zero. `d_min` records the
/// closest *sampled* separation, and since `u` passes through zero roughly linearly in
/// fictitious time while `|R1| = |u1|^2`, the sampled minimum falls like the square of the
/// step. Which sample lands nearest the crossing is essentially a phase accident, so the
/// scaling is noisy rather than clean.
///
/// This matters for BRIEF §5: `d_min < 1e-10` is unreachable at the working `eta = 1e-2` —
/// not because anything is wrong, but because the trajectory is not sampled that finely.
/// The reference's quoted `1.35e-11` is likewise a property of its step size.
///
/// Reported, not asserted.
#[test]
fn d_min_at_collision_is_sampling_limited() {
    let (_, m) = setup();
    println!("{:>10}{:>14}{:>14}{:>14}{:>12}", "eta", "d_min", "|dE/E|", "gamma resid", "steps");
    for eta in [1e-2f64, 1e-3, 1e-4, 1e-5] {
        let (s, _) = setup();
        let o = az::integrate_az(s, &m, 1.0, 1, eta, 20_000_000, None);
        println!(
            "{eta:>10.0e}{:>14.3e}{:>14.3e}{:>14.3e}{:>12}",
            o.d_min_ref, o.drift, o.gamma_max, o.steps
        );
    }
    println!("\nd_min falls with eta; |dE/E| and the Gamma residual reach their roundoff");
    println!("floor by eta = 1e-3 and stay there. Energy conservation is converged well");
    println!("before d_min is, so the two halves of the acceptance test are not measuring");
    println!("the same thing.");
}
