//! The anchor chain for Heggie's global regularisation.
//!
//! `tests/heggie_hamiltonian_fd.rs` finite-differences `gamma` against `deriv`. That test is
//! **worthless on its own**: a sign error present in both functions passes it silently. It means
//! something only because `gamma` is independently anchored, and this file is that anchor —
//! exactly the role `az_identities.rs` plays for AZ.
//!
//! Heggie's chain is one link longer than AZ's, because of the enlarged phase space:
//!
//! ```text
//!   L0  Eq. (18) is the Jacobian of Eq. (17), and both reduce to Levi-Civita in the plane
//!   L1  q_i = (1/2) A_i^T Q_i, and the LC round trip
//!   L2  Eqs. (10)+(12) recover the Cartesian state  <- the crossed-mass hazard lives here
//!   L3  Eq. (6) equals the Cartesian energy Eq. (1)
//!   L4  Gamma* == (energy_enlarged - h) * R1 R2 R3, off shell
//!   L5  sum q_i == 0 at registration
//! ```
//!
//! Every link carries a **mutation arm**: a deliberately wrong version that the same assertion
//! must reject. A test that cannot fail is indistinguishable from a test that passes, and three
//! of these links are checking identities that hold under a wide class of wrong constants.

use prin_rs::integrate::heggie::system::cyc;
use prin_rs::integrate::heggie::{hamiltonian, HgState, HgSystem};
use prin_rs::integrate::az::lc;
use prin_rs::physics::{energy, Cart};
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

const N: usize = 256;

/// Burrau's masses. **Unequal, and that is load-bearing**: Eq. (10) is
/// `(m_j q_k - m_k q_j)/M`, precisely the crossed-mass shape this project has on record as
/// invisible to a `sum p_i = 0` check. On an equal-mass slice a mass error is invisible by
/// construction, so an equal-mass control could not produce these results.
fn masses() -> [f64; 3] {
    [3.0, 4.0, 5.0]
}

/// Random Cartesian states, **centred**. Heggie's Eq. (8b) reduces to its simple form only in
/// the centre-of-mass frame, and `cart_from_enlarged` reconstructs into that frame, so an
/// uncentred input would fail the round trip for a reason that is not an algebra error.
fn random_carts(n: usize, seed: u64, m: &[f64; 3]) -> Vec<Cart<f64>> {
    let mut rng = SplitMix64::new(seed);
    let mtot: f64 = m.iter().sum();
    (0..n)
        .map(|_| {
            let sr = 10f64.powf(rng.range(-1.0, 1.0));
            let sv = 10f64.powf(rng.range(-1.0, 1.0));
            let mut r = [Vec2::zero(); 3];
            let mut v = [Vec2::zero(); 3];
            for i in 0..3 {
                r[i] = Vec2::new(rng.range(-2.0, 2.0) * sr, rng.range(-2.0, 2.0) * sr);
                v[i] = Vec2::new(rng.range(-2.0, 2.0) * sv, rng.range(-2.0, 2.0) * sv);
            }
            let (mut rc, mut vc) = (Vec2::<f64>::zero(), Vec2::<f64>::zero());
            for i in 0..3 {
                rc += r[i] * m[i];
                vc += v[i] * m[i];
            }
            rc = rc / mtot;
            vc = vc / mtot;
            for i in 0..3 {
                r[i] -= rc;
                v[i] -= vc;
            }
            Cart::new(r, v)
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// L0 — the planar reduction, which is the one step the paper does not state.

/// Heggie Eq. (17), written literally in three dimensions from a 4-vector `Q`.
fn eq17(q: [f64; 4]) -> [f64; 3] {
    [
        q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3],
        2.0 * (q[0] * q[1] - q[2] * q[3]),
        2.0 * (q[0] * q[2] + q[1] * q[3]),
    ]
}

/// Heggie Eq. (18), written literally as the 4x3 matrix `A_i`, rows indexed by `Q` component.
fn eq18(q: [f64; 4]) -> [[f64; 3]; 4] {
    [
        [2.0 * q[0], 2.0 * q[1], 2.0 * q[2]],
        [-2.0 * q[1], 2.0 * q[0], 2.0 * q[3]],
        [-2.0 * q[2], -2.0 * q[3], 2.0 * q[0]],
        [2.0 * q[3], -2.0 * q[2], 2.0 * q[1]],
    ]
}

/// **Eq. (18) is the Jacobian of Eq. (17)** — checked in full 3D by finite differences, so a
/// transcription error in either equation is caught without reference to the other.
///
/// This runs in four dimensions on purpose. Checking it only on the planar slice would leave the
/// two out-of-plane columns untested, and those columns are what make the *reduction* a claim
/// rather than a definition.
#[test]
fn eq18_is_the_jacobian_of_eq17() {
    let mut rng = SplitMix64::new(0x18_17);
    let mut worst = 0.0f64;
    for _ in 0..N {
        let q: [f64; 4] = std::array::from_fn(|_| rng.range(-2.0, 2.0));
        let a = eq18(q);
        for b in 0..4 {
            let h = 1e-6 * q[b].abs().max(1.0);
            let (mut hi, mut lo) = (q, q);
            hi[b] += h;
            lo[b] -= h;
            let (fh, fl) = (eq17(hi), eq17(lo));
            for c in 0..3 {
                let fd = (fh[c] - fl[c]) / (2.0 * h);
                let scale = a[b][c].abs().max(fd.abs()).max(1.0);
                worst = worst.max((a[b][c] - fd).abs() / scale);
            }
        }
    }
    println!("Eq. (18) against the finite-differenced Eq. (17), 4D: worst = {worst:.3e}");
    assert!(worst < 1e-8, "Eq. (18) is not the Jacobian of Eq. (17): {worst:e}");
}

/// The planar reduction: at `Q_3 = Q_4 = 0`, Eq. (17) is `lc::rho_of_u` and the leading 2x2 block
/// of Eq. (18) is `2 L(Q)^T`, which is what `HgSystem::a_transpose_apply` implements.
///
/// The mutation arm transposes the block. `L(Q)` and `L(Q)^T` agree whenever `Q.y = 0`, so the
/// arm is only meaningful off that line and the test says so by asserting a large discrepancy.
#[test]
fn the_planar_reduction_is_levi_civita() {
    let mut rng = SplitMix64::new(0x2d);
    let mut worst_q = 0.0f64;
    let mut worst_a = 0.0f64;
    let mut worst_transposed = 0.0f64;
    let mut worst_identity = 0.0f64;

    for _ in 0..N {
        let u = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0));
        let q4 = [u.x, u.y, 0.0, 0.0];

        // Eq. (17) restricted is rho_of_u, and the third component vanishes.
        let e17 = eq17(q4);
        let rho = lc::rho_of_u(u);
        worst_q = worst_q.max((e17[0] - rho.x).abs()).max((e17[1] - rho.y).abs()).max(e17[2].abs());

        // Eq. (18)'s leading 2x2 block is 2 L(u)^T.
        let a = eq18(q4);
        let lt = [
            [2.0 * u.x, 2.0 * u.y],
            [-2.0 * u.y, 2.0 * u.x],
        ];
        for b in 0..2 {
            for c in 0..2 {
                worst_a = worst_a.max((a[b][c] - lt[b][c]).abs());
                worst_transposed = worst_transposed.max((a[b][c] - lt[c][b]).abs());
            }
        }

        // Heggie's identity q_i = (1/2) A_i^T Q_i, through the shipped helper.
        let via_a = HgSystem::<f64>::a_transpose_apply(u, u) * 0.5;
        worst_identity = worst_identity.max((via_a - rho).norm());
    }

    println!("planar reduction: Eq.(17) vs rho_of_u {worst_q:.3e}, Eq.(18) vs 2L^T {worst_a:.3e}");
    println!("  q = (1/2) A^T Q via the shipped helper: {worst_identity:.3e}");
    println!("  transposed block (the mutation arm): {worst_transposed:.3e}");
    assert!(worst_q < 1e-14, "Eq. (17) restricted is not rho_of_u: {worst_q:e}");
    assert!(worst_a < 1e-14, "Eq. (18)'s planar block is not 2 L(Q)^T: {worst_a:e}");
    assert!(worst_identity < 1e-14, "q = (1/2) A^T Q failed: {worst_identity:e}");
    assert!(
        worst_transposed > 1e-2,
        "the transposed block agreed at {worst_transposed:e} — this arm has no teeth"
    );
}

// ---------------------------------------------------------------------------------------------
// L2 — the reconstruction, and the crossed-mass hazard.

/// `Cart -> (q_i, p_i) -> (q_i*, p_i*)` recovers the original state, Heggie Eqs. (8b), (10), (12).
///
/// The mutation arm swaps `m_j` and `m_k` inside Eq. (10) only. That is the exact defect this
/// project has recorded as invisible to `sum p_i = 0`, and it is invisible on equal masses too —
/// which is why this runs at Burrau's `(3, 4, 5)`.
#[test]
fn the_enlarged_reconstruction_recovers_the_cartesian_state() {
    let m = masses();
    let mtot: f64 = m.iter().sum();
    let sys = HgSystem::new(m);
    let mut worst_r = 0.0f64;
    let mut worst_v = 0.0f64;
    let mut worst_swapped = 0.0f64;

    for c in random_carts(N, 0xA10, &m) {
        let (q, p) = sys.enlarged_from_cart(&c);
        let back = sys.cart_from_enlarged(&q, &p);
        for i in 0..3 {
            let sr = c.r[i].norm().max(1.0);
            let sv = c.v[i].norm().max(1.0);
            worst_r = worst_r.max((back.r[i] - c.r[i]).norm() / sr);
            worst_v = worst_v.max((back.v[i] - c.v[i]).norm() / sv);

            // Eq. (10) with the two masses crossed.
            let (j, k) = cyc(i);
            let swapped = (q[k] * m[k] - q[j] * m[j]) / mtot;
            worst_swapped = worst_swapped.max((swapped - c.r[i]).norm() / sr);
        }
    }
    println!("Eqs. (10)/(12) round trip at m = (3,4,5): worst dr = {worst_r:.3e}, dv = {worst_v:.3e}");
    println!("  crossed masses in Eq. (10), the mutation arm: {worst_swapped:.3e}");
    assert!(worst_r < 1e-13, "position round trip failed: {worst_r:e}");
    assert!(worst_v < 1e-13, "velocity round trip failed: {worst_v:e}");
    assert!(
        worst_swapped > 1e-2,
        "crossed masses still round-tripped at {worst_swapped:e} — this arm has no teeth"
    );
}

/// `sum q_i == 0` — Heggie's Eq. (9), at registration. The free invariant the march will watch.
#[test]
fn sum_q_vanishes_at_registration() {
    let m = masses();
    let sys = HgSystem::new(m);
    let mut worst = 0.0f64;
    let mut worst_state = 0.0f64;
    for c in random_carts(N, 0x59, &m) {
        let (q, _) = sys.enlarged_from_cart(&c);
        let scale = q.iter().fold(0.0f64, |a, v| a.max(v.norm())).max(1e-300);
        worst = worst.max((q[0] + q[1] + q[2]).norm() / scale);
        let (s, _) = sys.to_reg(&c);
        worst_state = worst_state.max(sys.sum_q_residual(&s));
    }
    println!("sum q_i / max|q_i| at registration: {worst:.3e} direct, {worst_state:.3e} via LC");
    assert!(worst < 1e-14, "sum q_i does not vanish: {worst:e}");
    assert!(worst_state < 1e-13, "sum q_i does not vanish through the LC map: {worst_state:e}");
}

// ---------------------------------------------------------------------------------------------
// L3 — the enlarged Hamiltonian is the physical energy.

/// Heggie Eq. (6) on `(q_i, p_i)` equals Eq. (1) on the Cartesian state.
///
/// Two mutation arms, each aimed at a pairing that fails silently:
///   - **rotate `mu`**, putting `mu_{31}` on `|p_1|^2`. The cyclic form in `hamiltonian.rs` exists
///     to prevent exactly this, and the test has to prove it would be caught.
///   - **drop the coupling term** `- p_j . p_k / m_i`. Heggie's `q_i` are relative-to-each-other,
///     not Jacobi, so the cross terms are mandatory; dropping them is the AZ hazard documented at
///     `system.rs:107` transplanted to three vectors.
#[test]
fn the_enlarged_hamiltonian_is_the_cartesian_energy() {
    let m = masses();
    let sys = HgSystem::new(m);
    let mut worst = 0.0f64;
    let mut worst_rotated = 0.0f64;
    let mut worst_nocoupling = 0.0f64;

    for c in random_carts(N, 0xE6, &m) {
        let e_cart = sys.energy_cartesian(&c);
        let (q, p) = sys.enlarged_from_cart(&c);
        let e_enl = sys.energy_enlarged(&q, &p);
        let scale = e_cart.abs().max(1.0);
        worst = worst.max((e_enl - e_cart).abs() / scale);

        // Arm 1: mu rotated by one.
        let mut rot = 0.0;
        for i in 0..3 {
            let (j, k) = cyc(i);
            rot += p[i].norm_sq() / (2.0 * sys.mu[j]) - p[j].dot(p[k]) / m[i]
                - sys.mm[i] / q[i].norm();
        }
        worst_rotated = worst_rotated.max((rot - e_cart).abs() / scale);

        // Arm 2: the coupling dropped.
        let mut nc = 0.0;
        for i in 0..3 {
            nc += p[i].norm_sq() / (2.0 * sys.mu[i]) - sys.mm[i] / q[i].norm();
        }
        worst_nocoupling = worst_nocoupling.max((nc - e_cart).abs() / scale);
    }
    println!("Eq. (6) against Eq. (1): worst relative = {worst:.3e}");
    println!("  mu rotated by one: {worst_rotated:.3e}");
    println!("  coupling dropped:  {worst_nocoupling:.3e}");
    assert!(worst < 1e-13, "Eq. (6) is not the Cartesian energy: {worst:e}");
    assert!(worst_rotated > 1e-2, "a rotated mu still agreed at {worst_rotated:e}");
    assert!(worst_nocoupling > 1e-2, "dropping the coupling still agreed at {worst_nocoupling:e}");
}

/// The energy survives the round trip through the LC map, so `to_reg`'s frozen `h` is the
/// physical energy and not merely something self-consistent.
#[test]
fn the_frozen_energy_is_the_cartesian_energy() {
    let m = masses();
    let sys = HgSystem::new(m);
    let mut worst_h = 0.0f64;
    let mut worst_back = 0.0f64;
    for c in random_carts(N, 0x7C, &m) {
        let (s, h) = sys.to_reg(&c);
        let e_cart = energy::energy(&c.r, &c.v, &m, 0.0);
        let scale = e_cart.abs().max(1.0);
        worst_h = worst_h.max((h - e_cart).abs() / scale);
        let back = sys.to_cartesian(&s);
        for i in 0..3 {
            worst_back = worst_back
                .max((back.r[i] - c.r[i]).norm() / c.r[i].norm().max(1.0))
                .max((back.v[i] - c.v[i]).norm() / c.v[i].norm().max(1.0));
        }
    }
    println!("to_reg's frozen h against the Cartesian energy: {worst_h:.3e}");
    println!("Cart -> regularised -> Cart round trip:         {worst_back:.3e}");
    assert!(worst_h < 1e-13, "the frozen h is not the physical energy: {worst_h:e}");
    assert!(worst_back < 1e-12, "the full round trip failed: {worst_back:e}");
}

// ---------------------------------------------------------------------------------------------
// L4 — Gamma* is the energy defect, times the time factor.

/// Random regularised states, spanning orders of magnitude, with **each of the three pairs in
/// turn made the closest**.
///
/// AZ's random states are implicitly conditioned by the reference-body choice: its `u1`, `u2` are
/// the two regularised sides and `R3` is always the longest. Heggie has no such ordering, and
/// sampling the same way would leave the globality — the property the whole method exists for —
/// untested.
fn random_states(n: usize, seed: u64) -> Vec<(HgState<f64>, f64)> {
    let mut rng = SplitMix64::new(seed);
    (0..n)
        .map(|i| {
            let tight = i % 4; // 3 means "no pair singled out"
            let su = 10f64.powf(rng.range(-1.0, 1.0));
            let sp = 10f64.powf(rng.range(-1.0, 1.0));
            let mut s = HgState { u: [Vec2::zero(); 3], p: [Vec2::zero(); 3], t: 0.0 };
            for a in 0..3 {
                let shrink = if a == tight { 1e-2 } else { 1.0 };
                s.u[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * (su * shrink);
                s.p[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * sp;
            }
            (s, rng.range(-30.0, 5.0))
        })
        .collect()
}

/// `Gamma* == (energy_enlarged - h) * R1 R2 R3`, with **random off-shell `h`**.
///
/// Off-shell is the point: on the solution path `Gamma* = 0` identically, so an on-shell check
/// would compare zero against zero and pass for any `gamma` whatsoever.
#[test]
fn gamma_equals_the_energy_defect_times_the_time_factor() {
    let sys = HgSystem::new(masses());
    let mut worst = 0.0f64;
    let mut worst_sign = 0.0f64;
    for (s, h) in random_states(N, 0x4A11) {
        let g = hamiltonian::gamma(&sys, &s, h);
        let rp = s.r(0) * s.r(1) * s.r(2);
        let want = (sys.energy_of(&s) - h) * rp;
        let scale = g.abs().max(want.abs()).max(1e-300);
        worst = worst.max((g - want).abs() / scale);
        // Mutation arm: the energy term entering with the wrong sign.
        let flipped = g + 2.0 * h * rp;
        worst_sign = worst_sign.max((flipped - want).abs() / scale);
    }
    println!("Gamma* against (H - h) R1R2R3, off shell: worst = {worst:.3e}");
    println!("  energy term sign flipped: {worst_sign:.3e}");
    assert!(worst < 1e-12, "Gamma* is not the energy defect: {worst:e}");
    assert!(worst_sign > 1e-2, "a flipped energy-term sign still agreed at {worst_sign:e}");
}

/// `Gamma*` vanishes on shell — for **every** state, since `h` is then defined from the state.
/// This is the residual the march will watch, so it must start at zero.
#[test]
fn gamma_vanishes_on_shell() {
    let sys = HgSystem::new(masses());
    let mut worst = 0.0f64;
    for (s, _) in random_states(N, 0x0F5) {
        let h = sys.energy_of(&s);
        worst = worst.max(hamiltonian::gamma_residual(&sys, &s, h));
    }
    println!("gamma_residual on shell: worst = {worst:.3e}");
    assert!(worst < 1e-13, "Gamma* does not vanish on shell: {worst:e}");
}
