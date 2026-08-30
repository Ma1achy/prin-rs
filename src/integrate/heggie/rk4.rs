//! Classic RK4 in fictitious time, over the thirteen-component Heggie state.
//!
//! **Not symplectic and not time-symmetric**, exactly as AZ's is. The two integrators must be
//! compared under the same stepper or the comparison is scoring the stepper.
//!
//! The accumulation order matches `az::rk4::step` — `(k1 + 2 k2 + 2 k3 + k4)` left to right, with
//! `h/6` formed once before scaling. That is ulp-load-bearing there because of the NumPy
//! cross-check; here it is kept so that an AZ/Heggie difference is a difference of *equations*
//! and not of floating-point accumulation order.

use crate::Real;

use super::hamiltonian::{deriv_time, HgTime};
use super::state::HgState;
use super::system::HgSystem;

/// One step: `s + (h/6)(k1 + 2 k2 + 2 k3 + k4)`.
#[inline]
pub fn step<T: Real>(
    sys: &HgSystem<T>,
    s: &HgState<T>,
    h_energy: T,
    time: HgTime,
    h: T,
) -> HgState<T> {
    let half = T::lit(0.5);
    let two = T::lit(2.0);
    let six = T::lit(6.0);

    let k1 = deriv_time(sys, s, h_energy, time);
    let k2 = deriv_time(sys, &s.axpy(half * h, &k1), h_energy, time);
    let k3 = deriv_time(sys, &s.axpy(half * h, &k2), h_energy, time);
    let k4 = deriv_time(sys, &s.axpy(h, &k3), h_energy, time);

    let mut acc = HgState {
        u: [Default::default(); 3],
        p: [Default::default(); 3],
        t: k1.t + k2.t * two + k3.t * two + k4.t,
    };
    for i in 0..3 {
        acc.u[i] = k1.u[i] + k2.u[i] * two + k3.u[i] * two + k4.u[i];
        acc.p[i] = k1.p[i] + k2.p[i] * two + k3.p[i] * two + k4.p[i];
    }
    s.axpy(h / six, &acc)
}
