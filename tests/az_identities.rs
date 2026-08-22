//! Links 1-3 of the Step 4 validation chain.
//!
//! An FD test alone is not sufficient: a sign error present in **both** `gamma` and `deriv`
//! passes it silently. So `gamma` is anchored here to `energy_phys`, which is anchored to
//! the Cartesian `energy`, which is anchored to the Burrau constants. A sign error in
//! `gamma` cannot survive this chain. Only then does the FD test say anything about `deriv`.

use prin_rs::integrate::az::hamiltonian::gamma;
use prin_rs::integrate::az::reference_body::triple;
use prin_rs::integrate::az::{AzState, AzSystem};
use prin_rs::physics::{burrau, energy, Cart};
use prin_rs::rng::SplitMix64;
use prin_rs::{Real, Vec2};

const N: usize = 512;

/// Random configurations **in the COM frame**.
///
/// The zero-total-momentum condition is not cosmetic: `energy_phys` is the energy in
/// relative coordinates, so it equals the Cartesian energy only when the centre of mass is
/// at rest. Comparing them on a drifting configuration would fail for a reason that has
/// nothing to do with the algebra.
fn random_com_configs(n: usize, seed: u64) -> Vec<Cart<f64>> {
    let m = burrau::masses::<f64>();
    let mtot: f64 = m.iter().sum();
    let mut rng = SplitMix64::new(seed);
    (0..n)
        .map(|_| {
            let mut r = [Vec2::zero(); 3];
            let mut v = [Vec2::zero(); 3];
            for k in 0..3 {
                r[k] = Vec2::new(rng.range(-4.0, 4.0), rng.range(-4.0, 4.0));
                v[k] = Vec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
            }
            let mut pc = Vec2::zero();
            for k in 0..3 {
                pc += v[k] * m[k];
            }
            pc = pc / mtot;
            for k in 0..3 {
                v[k] -= pc;
            }
            Cart::new(r, v)
        })
        .collect()
}

fn systems() -> Vec<AzSystem<f64>> {
    let m = burrau::masses::<f64>();
    (0..3)
        .map(|a| {
            let (x, y, z) = triple(a);
            AzSystem::new(x, y, z, m)
        })
        .collect()
}

/// Link 2a: the registered energy equals the Cartesian energy.
///
/// This validates the velocity-space mass matrix `K`, `mu1`/`mu2`, the `P1.P2/ma` kinetic
/// cross term, and the LC map — all against a quantity already anchored by the Burrau
/// constants.
#[test]
fn energy_phys_matches_cartesian_energy() {
    let m = burrau::masses::<f64>();
    let mut worst = 0.0f64;
    for sys in systems() {
        for c in random_com_configs(N, 0xA21) {
            let (_, e) = sys.to_reg(&c);
            let ec = energy::energy(&c.r, &c.v, &m, 0.0);
            let rel = (e - ec).abs() / ec.abs().max(1e-30);
            worst = worst.max(rel);
        }
    }
    println!("energy_phys vs Cartesian energy: worst relative deviation = {worst:.3e}");
    assert!(worst < 1e-13, "worst = {worst:e}");
}

/// Link 2b: registration round-trips.
///
/// The tolerance here is set from measurement, not by eye. `u_of_rho` computes
/// `u0 = sqrt((|rho| + rho.x)/2)`, which cancels catastrophically when `rho` points along
/// negative x; `u1 = rho.y/(2 u0)` then amplifies the damage. `tests/lc_conditioning.rs`
/// measures the loss as roughly `eps / (u0/|u|)^2` — 1.9e-13 at 179 degrees, 3.5e-9 at
/// 179.99. Over random configurations the worst case lands near 1e-10.
///
/// This is a property of the LC branch and is present identically in the reference, so it
/// is transcribed rather than fixed: match first, change one thing at a time. The stable
/// variant (compute whichever of `u0`, `u1` is larger, derive the other) is a candidate for
/// after the Step 4 gate, and is likely to matter a great deal at f32.
#[test]
fn to_reg_to_cartesian_round_trips() {
    let m = burrau::masses::<f64>();
    let mtot: f64 = m.iter().sum();
    let mut worst = 0.0f64;
    let mut worst_cond = 1.0f64;

    for sys in systems() {
        for c in random_com_configs(N, 0xB33) {
            let (s, _) = sys.to_reg(&c);
            let back = sys.to_cartesian(&s);

            // to_cartesian returns the COM frame, so compare against the COM-shifted input.
            let mut rc = Vec2::zero();
            for k in 0..3 {
                rc += c.r[k] * m[k];
            }
            rc = rc / mtot;

            let mut here = 0.0f64;
            for k in 0..3 {
                let want = c.r[k] - rc;
                here = here.max((back.r[k] - want).norm() / want.norm().max(1.0));
                here = here.max((back.v[k] - c.v[k]).norm() / c.v[k].norm().max(1.0));
            }
            if here > worst {
                worst = here;
                // How close the registration sat to the LC branch cut.
                let c1 = s.u1.x / s.u1.norm().max(1e-300);
                let c2 = s.u2.x / s.u2.norm().max(1e-300);
                worst_cond = c1.abs().min(c2.abs());
            }
        }
    }
    println!("to_cartesian(to_reg(.)) round trip: worst relative deviation = {worst:.3e}");
    println!("  u0/|u| at the worst case = {worst_cond:.3e}  (small means near the LC branch cut)");
    println!("  see tests/lc_conditioning.rs — this is the LC branch, not the AZ algebra");
    assert!(worst < 1e-9, "worst = {worst:e}, beyond the measured LC conditioning floor");
}

/// Link 3, the one that breaks the circularity: `Gamma == A*B*(H - E)`.
///
/// Random `s` **and** random `E` — not `E` taken from registration — so the identity is
/// tested as an algebraic statement about `gamma`, off-shell as well as on. A sign error in
/// `gamma` cannot survive this.
#[test]
fn gamma_equals_ab_times_energy_defect() {
    let mut rng = SplitMix64::new(0xC0FFEE);
    let mut worst = 0.0f64;
    let mut worst_at = String::new();

    for sys in systems() {
        for _ in 0..N {
            let s = AzState {
                u1: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
                p1: Vec2::new(rng.range(-3.0, 3.0), rng.range(-3.0, 3.0)),
                u2: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
                p2: Vec2::new(rng.range(-3.0, 3.0), rng.range(-3.0, 3.0)),
                t: 0.0,
            };
            let e = rng.range(-30.0, 10.0);

            let lhs = gamma(&sys, &s, e);
            let rhs = s.a() * s.b() * (sys.energy_of(&s) - e);
            let scale = lhs.abs().max(rhs.abs()).max(1e-300);
            let rel = (lhs - rhs).abs() / scale;
            if rel > worst {
                worst = rel;
                worst_at = format!("a={} A={:.3e} B={:.3e} E={:.3}", sys.a, s.a(), s.b(), e);
            }
        }
    }
    println!("Gamma == A*B*(H - E): worst relative deviation = {worst:.3e}   at {worst_at}");
    assert!(worst < 1e-12, "worst = {worst:e} at {worst_at}");
}

/// The identity must have teeth: a deliberately wrong Gamma has to fail it.
#[test]
fn the_gamma_identity_detects_a_sign_error() {
    let sys = systems()[0];
    let mut rng = SplitMix64::new(0xDEAD);
    let mut worst = 0.0f64;
    for _ in 0..N {
        let s = AzState {
            u1: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
            p1: Vec2::new(rng.range(-3.0, 3.0), rng.range(-3.0, 3.0)),
            u2: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
            p2: Vec2::new(rng.range(-3.0, 3.0), rng.range(-3.0, 3.0)),
            t: 0.0,
        };
        let e = rng.range(-30.0, 10.0);
        // Flip the sign of the LC cross term — the single most plausible transcription slip,
        // and one with no visible symptom in a trajectory.
        let bad = gamma(&sys, &s, e)
            - <f64 as Real>::lit(2.0)
                * prin_rs::integrate::az::lc::l_apply(s.u1, s.p1)
                    .dot(prin_rs::integrate::az::lc::l_apply(s.u2, s.p2))
                / (<f64 as Real>::lit(4.0) * sys.ma);
        let rhs = s.a() * s.b() * (sys.energy_of(&s) - e);
        let rel = (bad - rhs).abs() / bad.abs().max(rhs.abs()).max(1e-300);
        worst = worst.max(rel);
    }
    println!("sign-flipped Gamma: worst relative deviation = {worst:.3e} (must be large)");
    assert!(worst > 1e-3, "a sign-flipped Gamma still passed at {worst:e} — the identity has no teeth");
}
