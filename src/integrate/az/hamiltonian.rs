//! The regularised Hamiltonian and its derivatives, kept **adjacent on purpose**.
//!
//! `Gamma` is not implemented in the reference — it exists only in `tb_az.py`'s docstring.
//! Writing it here is what makes the finite-difference test possible at all. But an FD test
//! alone is not sufficient: a sign error present in *both* `gamma` and `deriv` passes it
//! silently. So `gamma` is independently anchored by the identity
//!
//!     Gamma(s, E) == A * B * (energy_phys(s) - E)
//!
//! where `energy_phys` is itself anchored to the Cartesian energy, and thence to the Burrau
//! constants. A sign error in `Gamma` cannot survive that chain. Only then does the FD test
//! mean anything about `deriv`.
//!
//! Transcribed from `reference/tb_az.py:deriv` (lines 112-142) and the docstring (27-33).

use crate::Real;

use super::lc;
use super::state::AzState;
use super::system::AzSystem;

/// `Gamma = A B (H - E)`, from the docstring:
///
/// ```text
/// Gamma = B|p1|^2/(8 mu1) + A|p2|^2/(8 mu2) + (L1 p1).(L2 p2)/(4 m_a)
///         - m_a m_b B - m_a m_c A - A B m_b m_c/|R3| - E A B
/// ```
///
/// Both `-m_a m_b/|R1|` and `-m_a m_c/|R2|` have become constant-order terms, linear in `B`
/// and `A`. Nothing is singular unless `|R3| -> 0`, which is a genuine triple collision.
pub fn gamma<T: Real>(sys: &AzSystem<T>, s: &AzState<T>, e: T) -> T {
    let a = s.a();
    let b = s.b();
    let r1 = lc::rho_of_u(s.u1);
    let r2 = lc::rho_of_u(s.u2);
    let r3 = (r2 - r1).norm().max(T::TINY);

    let l1p1 = lc::l_apply(s.u1, s.p1);
    let l2p2 = lc::l_apply(s.u2, s.p2);

    let eight = T::lit(8.0);
    let four = T::lit(4.0);

    b * s.p1.norm_sq() / (eight * sys.mu1)
        + a * s.p2.norm_sq() / (eight * sys.mu2)
        + l1p1.dot(l2p2) / (four * sys.ma)
        - sys.ma * sys.mb * b
        - sys.ma * sys.mc * a
        - a * b * sys.mb * sys.mc / r3
        - e * a * b
}

/// `dGamma/dp1`, `-dGamma/du1`, `dGamma/dp2`, `-dGamma/du2`, `dt/dtau = A*B`.
///
/// `A` and `B` are **not** floored here, though `phys_from_state` floors them. The asymmetry
/// is the reference's and is transcribed as-is.
pub fn deriv<T: Real>(sys: &AzSystem<T>, s: &AzState<T>, e: T) -> AzState<T> {
    let a = s.a();
    let b = s.b();
    let r1 = lc::rho_of_u(s.u1);
    let r2 = lc::rho_of_u(s.u2);
    let r3v = r2 - r1;
    let r3 = r3v.norm().max(T::TINY);

    let l1p1 = lc::l_apply(s.u1, s.p1);
    let l2p2 = lc::l_apply(s.u2, s.p2);
    let n1 = s.p1.norm_sq();
    let n2 = s.p2.norm_sq();
    let mbc = sys.mb * sys.mc;

    let four = T::lit(4.0);
    let two = T::lit(2.0);

    // dGamma/dp1 and dGamma/dp2. `lt_apply(u1, l2p2)` is L(u1)^T applied to L(u2)p2.
    let du1 = s.p1 * (b / (four * sys.mu1)) + lc::lt_apply(s.u1, l2p2) / (four * sys.ma);
    let du2 = s.p2 * (a / (four * sys.mu2)) + lc::lt_apply(s.u2, l1p1) / (four * sys.ma);

    // dGamma/du1. Five terms.
    //   1. the *other* pair's |p|^2 and mu — A appears in the p2 term, so d/du1 hits it
    //   2. the LC cross term, with p1 in the MATRIX slot: L(u)w = L(w)u, so this is exact
    //   3. the constant-order binary term, carrying the CROSS mass pair ma*mc
    //   4. the unregularised R3 term, two sub-pieces
    //   5. the energy term
    let g1 = s.u1 * (n2 / (four * sys.mu2))
        + lc::lt_apply(s.p1, l2p2) / (four * sys.ma)
        - s.u1 * (two * sys.ma * sys.mc)
        - (s.u1 * (two * b / r3)
            + lc::lt_apply(s.u1, r3v) * (two * a * b / (r3 * r3 * r3)))
            * mbc
        - s.u1 * (two * e * b);

    // dGamma/du2. Same five, with two differences that are easy to miss:
    //   - term 3 carries ma*mb here, not ma*mc
    //   - the second R3 sub-piece flips SIGN, because dR3/dR1 = -I while dR3/dR2 = +I
    let g2 = s.u2 * (n1 / (four * sys.mu1))
        + lc::lt_apply(s.p2, l1p1) / (four * sys.ma)
        - s.u2 * (two * sys.ma * sys.mb)
        - (s.u2 * (two * a / r3)
            - lc::lt_apply(s.u2, r3v) * (two * a * b / (r3 * r3 * r3)))
            * mbc
        - s.u2 * (two * e * a);

    AzState {
        u1: du1,
        p1: -g1,
        u2: du2,
        p2: -g2,
        t: a * b,
    }
}

/// A free integration-quality residual.
///
/// `E` is the true energy at registration, so `Gamma == 0` along the exact trajectory. The
/// residual is normalised to `|H - E| / |E|`, which puts it on the same footing as
/// `energy_drift` and makes it directly comparable.
///
/// This is *not* the withdrawn claim that `Gamma` vanishes identically as a function — it
/// does not, and `dGamma/du` is exactly what drives the motion. It vanishes only *along the
/// trajectory*, which is what makes it a residual worth watching.
pub fn gamma_residual<T: Real>(sys: &AzSystem<T>, s: &AzState<T>, e: T) -> T {
    let ab = s.a() * s.b();
    if !ab.is_finite() || ab <= T::zero() {
        return T::infinity();
    }
    (gamma(sys, s, e) / ab).abs() / e.abs().max(T::DRIFT_FLOOR)
}
