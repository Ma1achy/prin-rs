//! Adaptive-timestep KDK leapfrog, **unregularised**.
//!
//! This is BRIEF §2.3(a) without (b). It is expected to fail on close encounters: gravity
//! diverges as two bodies approach, the adaptive step is driven toward zero, and the step
//! budget is exhausted. That failure is the entire reason Aarseth–Zare exists, so it is
//! recorded here rather than tuned around.
//!
//! Note the adaptive step breaks the symplecticity of plain KDK. That is accepted — this
//! integrator is a stepping stone, not the kernel.

use crate::physics::{energy, newton, Cart};
use crate::Real;

use super::TrajOut;

pub fn integrate<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    eta: T,
    eps2: T,
    max_steps: usize,
) -> TrajOut<T> {
    let mut r = s0.r;
    let mut v = s0.v;
    let e0 = energy::energy(&r, &v, m, eps2);
    let mut t = T::zero();
    let mut d_min = T::infinity();
    let mut steps = 0usize;
    let mut finite = true;

    let half = T::lit(0.5);
    let mut a = newton::accel(&r, m, eps2);

    while t < t_max && steps < max_steps {
        let mut dt = newton::adaptive_dt(&r, m, eta);
        // NaN >= x is false, so a non-finite dt would never satisfy any loop guard and the
        // trajectory would burn its whole budget: measured 354 s against 3 s nominal
        // (BRIEF §5.3). Test is_finite explicitly.
        if !dt.is_finite() || dt <= T::zero() {
            finite = false;
            break;
        }
        if t + dt > t_max {
            dt = t_max - t;
        }

        for k in 0..3 {
            v[k] += a[k] * (half * dt);
            r[k] += v[k] * dt;
        }
        a = newton::accel(&r, m, eps2);
        for k in 0..3 {
            v[k] += a[k] * (half * dt);
        }
        t += dt;
        steps += 1;

        let d = newton::pair_dists(&r);
        let step_min = d[0].min(d[1]).min(d[2]);
        // Only fold a finite separation into d_min. min(x, NaN) is NaN in the reference,
        // which silently poisons the field for every consumer downstream.
        if step_min.is_finite() && step_min < d_min {
            d_min = step_min;
        }

        if !r.iter().chain(v.iter()).all(|p| p.is_finite()) {
            finite = false;
            break;
        }
    }

    let e1 = energy::energy(&r, &v, m, eps2);
    let drift = ((e1 - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs();

    TrajOut {
        state: Cart::new(r, v),
        t,
        drift,
        d_min,
        steps,
        finite,
        budget_exhausted: steps >= max_steps && t < t_max,
    }
}
