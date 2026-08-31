//! The two steppers, and the reason there are two.
//!
//! Mikkola & Merritt (2008) are explicit that *"in these new methods the regularization is
//! achieved by using the leapfrog, hence the name algorithmic regularization"*. The
//! transformation alone is not the method — a leapfrog on `Lambda` is **exact** for the two-body
//! problem, and that exactness is the regularisation.
//!
//! So logH is run under both:
//!
//! - [`kdk`] — the method as designed. One force evaluation per step.
//! - [`rk4`] — the stepper AZ and Heggie both use, so this arm is the one directly comparable
//!   to them. Four force evaluations per step.
//!
//! **The RK4 arm is expected to lose most of the method, and the prediction is specific.** On
//! shell `K + B == U`, so an integrator that evaluates both denominators at the same point sees
//! only `dt/ds = 1/U` — a plain Sundman transformation. The leapfrog sees them at *different*
//! points, which is the entire asymmetry. If the two arms come out alike, that prediction is
//! wrong and so is the sentence from Mikkola & Merritt it rests on.
//!
//! # There is no common stepper, and pretending otherwise would be a confound
//!
//! `Az + kdk` and `Heggie + kdk` do not exist: their `Gamma` couples position and momentum, so
//! neither Hamiltonian is separable and KDK does not apply. Matching on *steps* across steppers
//! is therefore meaningless — RK4 does four force evaluations where KDK does one. Every entry
//! point here returns its evaluation count, and the harnesses match on that.

use crate::Real;

use super::hamiltonian::{denominators, deriv, deriv_with, Dens, LhTime};
use super::state::LhState;
use super::system::LhSystem;

/// One RK4 step in fictitious time. **Four force evaluations**, returned.
///
/// The accumulation order matches `az::rk4::step` and `heggie::rk4::step` — `(k1 + 2 k2 + 2 k3 +
/// k4)` left to right, with `h/6` formed once before scaling. Kept identical so that a
/// difference between integrators is a difference of *equations* and not of floating-point
/// accumulation order.
#[inline]
pub fn rk4<T: Real>(
    sys: &LhSystem<T>,
    s: &LhState<T>,
    b: T,
    time: LhTime,
    h: T,
) -> (LhState<T>, usize) {
    let half = T::lit(0.5);
    let two = T::lit(2.0);
    let six = T::lit(6.0);

    let k1 = deriv(sys, s, b, time);
    let k2 = deriv(sys, &s.axpy(half * h, &k1), b, time);
    let k3 = deriv(sys, &s.axpy(half * h, &k2), b, time);
    let k4 = deriv(sys, &s.axpy(h, &k3), b, time);

    let mut acc = LhState {
        r: [Default::default(); 3],
        v: [Default::default(); 3],
        t: k1.t + k2.t * two + k3.t * two + k4.t,
    };
    for i in 0..3 {
        acc.r[i] = k1.r[i] + k2.r[i] * two + k3.r[i] * two + k4.r[i];
        acc.v[i] = k1.v[i] + k2.v[i] * two + k3.v[i] * two + k4.v[i];
    }
    (s.axpy(h / six, &acc), 4)
}

/// One drift of fictitious length `h`: `dt = h / (K + B)`, `r += v dt`, `t += dt`.
///
/// **No force evaluation.** `K` is a function of the velocities, which a drift does not change,
/// so the denominator is a read.
#[inline]
pub fn drift<T: Real>(sys: &LhSystem<T>, s: &mut LhState<T>, b: T, time: LhTime, h: T) {
    let d = denominators(sys, s, b, time).drift.max(T::TINY);
    let dt = h / d;
    for i in 0..3 {
        s.r[i] = s.r[i] + s.v[i] * dt;
    }
    s.t = s.t + dt;
}

/// One kick of fictitious length `h`: `dt = h / U`, `v += a dt`.
///
/// **One force evaluation**, returned. `U` is a function of the positions, which a kick does not
/// change, so this is the only half that touches `accel`.
#[inline]
pub fn kick<T: Real>(sys: &LhSystem<T>, s: &mut LhState<T>, b: T, time: LhTime, h: T) -> usize {
    let d = denominators(sys, s, b, time).kick.max(T::TINY);
    // Reuse the equations of motion rather than re-deriving the kick, so a sign convention lives
    // in exactly one place. `deriv_with`'s `.v` is `a_i / kick`; a unit drift denominator makes
    // the `.r` and `.t` components it also computes irrelevant here, and they are discarded.
    let k = deriv_with(sys, s, Dens { drift: T::one(), kick: d });
    for i in 0..3 {
        s.v[i] = s.v[i] + k.v[i] * h;
    }
    1
}

/// One drift-kick-drift step. **One force evaluation**, returned.
///
/// `X(h/2) V(h) X(h/2)`, Mikkola's own composition. Time-symmetric, and the symmetry is what
/// makes the two-body case exact — not the transformation on its own.
#[inline]
pub fn kdk<T: Real>(
    sys: &LhSystem<T>,
    s: &LhState<T>,
    b: T,
    time: LhTime,
    h: T,
) -> (LhState<T>, usize) {
    let half = T::lit(0.5);
    let mut out = *s;
    drift(sys, &mut out, b, time, half * h);
    let n = kick(sys, &mut out, b, time, h);
    drift(sys, &mut out, b, time, half * h);
    (out, n)
}

/// Which stepper a march uses.
///
/// `Rk4` is a fair comparison against AZ and Heggie and an unfair one against the method;
/// `Kdk` is the reverse. Both are run, and the gap between them is a measurement of how much of
/// logH's behaviour is the stepper — which no single arm can report.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Stepper {
    /// Classic RK4, four force evaluations per step. What AZ and Heggie run.
    #[default]
    Rk4,
    /// Drift-kick-drift leapfrog, one force evaluation per step. What logH is designed for.
    Kdk,
}

impl Stepper {
    /// Force evaluations per step. Used to match arms on evaluations rather than steps.
    pub fn evals_per_step(self) -> usize {
        match self {
            Stepper::Rk4 => 4,
            Stepper::Kdk => 1,
        }
    }
}

/// One step under the chosen stepper, returning the new state and the force evaluations spent.
#[inline]
pub fn step<T: Real>(
    sys: &LhSystem<T>,
    s: &LhState<T>,
    b: T,
    time: LhTime,
    stepper: Stepper,
    h: T,
) -> (LhState<T>, usize) {
    match stepper {
        Stepper::Rk4 => rk4(sys, s, b, time, h),
        Stepper::Kdk => kdk(sys, s, b, time, h),
    }
}
