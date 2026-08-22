//! Classic RK4 in fictitious time.
//!
//! **Not symplectic and not time-symmetric.** The reference chose it to prove the physics,
//! not to ship. Match it at f64 first, then change one thing at a time.

use crate::Real;

use super::hamiltonian::deriv;
use super::state::AzState;
use super::system::AzSystem;

/// One step: `s + (h/6)(k1 + 2 k2 + 2 k3 + k4)`.
///
/// The accumulation order matches the reference's `(k1 + 2*k2 + 2*k3 + k4)` evaluated
/// left to right, and `h/6` is formed once before scaling — both matter at the ulp level,
/// which is where the cross-check lives.
#[inline]
pub fn step<T: Real>(sys: &AzSystem<T>, s: &AzState<T>, e: T, h: T) -> AzState<T> {
    let half = T::lit(0.5);
    let two = T::lit(2.0);
    let six = T::lit(6.0);

    let k1 = deriv(sys, s, e);
    let k2 = deriv(sys, &s.axpy(half * h, &k1), e);
    let k3 = deriv(sys, &s.axpy(half * h, &k2), e);
    let k4 = deriv(sys, &s.axpy(h, &k3), e);

    let acc = AzState {
        u1: k1.u1 + k2.u1 * two + k3.u1 * two + k4.u1,
        p1: k1.p1 + k2.p1 * two + k3.p1 * two + k4.p1,
        u2: k1.u2 + k2.u2 * two + k3.u2 * two + k4.u2,
        p2: k1.p2 + k2.p2 * two + k3.p2 * two + k4.p2,
        t: k1.t + k2.t * two + k3.t * two + k4.t,
    };
    s.axpy(h / six, &acc)
}
