//! The Heggie march.
//!
//! **What this file does not contain is the point.** `az::driver::integrate_az_opts` re-chooses a
//! reference body at every sync boundary, rebuilds `AzSystem`, round-trips the state through
//! Cartesian and re-registers. Measured, doubling that count at fixed step size moves the drift
//! field by **0.444 decades** — 6000x the effect of changing *which* reference is chosen.
//!
//! Here there is no reference body. `HgSystem` depends on nothing but the masses, `to_reg` is
//! called **once at `t = 0`**, and the sync loop is a *sampling* cadence: it records the residuals
//! and the closest approach, and it does not touch the state. The regularised march runs
//! uninterrupted from start to finish.
//!
//! Two free residuals come out of it, and one of them has no AZ analogue:
//!   - `|Gamma*|` against its largest term — zero along the exact trajectory, as AZ's is;
//!   - `|sum q_i| / max|q_i|` — Heggie's Eq. (9), an integral of the enlarged motion. AZ has
//!     nothing like it because AZ has no enlarged phase space.

use crate::physics::Cart;
use crate::Real;

use super::hamiltonian::{gamma_residual, HgTime};
use super::rk4;
use super::state::HgState;
use super::system::{cyc, HgSystem};

/// How `dtau` is sized within a sync interval. The same three arms as
/// [`DtauMode`](crate::integrate::az::DtauMode), named identically so a comparison between the
/// two integrators can hold this constant.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HgDtauMode {
    /// Sized once at interval entry and held.
    FixedPerInterval,
    /// `eta * (dt_left - s.t) / (dt/dtau)`. **Zeno by arithmetic** — the remaining time in the
    /// numerator gives `rem_{n+1} = rem_n (1 - eta)`, so the interval is approached geometrically
    /// and completed only through the relative landing tolerance. Named and kept as an axis,
    /// never a candidate; carried here so the AZ result reproduces on this integrator too.
    PerStepRemaining,
    /// Recomputed per step from the current `dt/dtau` with `dt_left` **held fixed**, capped at
    /// the interval's entry value. AZ's default and this one's.
    #[default]
    PerStepInterval,
}

#[derive(Clone, Copy, Debug)]
pub struct HgOpts<T> {
    /// Which of Heggie's time transformations. Default is Eqs. (22)-(24).
    pub time: HgTime,
    pub dtau_mode: HgDtauMode,
    /// Land the final step of each interval **on** the boundary rather than past it. A
    /// correctness property, not a preference: on AZ it is worth 1.06 -> 2.08 in measured
    /// convergence order.
    pub clamp_final_step: bool,
    /// `f` for the predictive per-step limit, `dtau <= f d_min / (|dq/dt|_max dt/dtau)`. Zero
    /// disables it. AZ ships `0.02`.
    pub step_limit_f: T,
    /// Collision radius, absolute. Zero disables the test.
    pub r_coll: T,
    pub stop_on_collision: bool,
    /// Re-derive the frozen `h` from the state at each sync boundary.
    ///
    /// **Default off, and it exists to be measured rather than used.** Refreshing `h` is the one
    /// thing this driver could do that would reintroduce a boundary-dependent quantity into an
    /// otherwise uninterrupted march — so it is the control arm for the claim that removing
    /// re-registration is what matters. A control that changes nothing proves nothing, so its
    /// effect is reported, not assumed.
    pub refresh_h_at_boundary: bool,
    pub lc_stable: bool,
}

impl<T: Real> Default for HgOpts<T> {
    fn default() -> Self {
        Self {
            time: HgTime::default(),
            dtau_mode: HgDtauMode::default(),
            clamp_final_step: true,
            step_limit_f: T::lit(0.02),
            r_coll: T::zero(),
            stop_on_collision: false,
            refresh_h_at_boundary: false,
            lc_stable: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct HgOut<T> {
    pub state: Cart<T>,
    /// Physical time actually reached.
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|`.
    pub drift: T,
    /// Closest approach over all three pairs, sampled at every RK4 step.
    ///
    /// **All three, symmetrically.** AZ's `d_min_ref` is over its two regularised pairs and needs
    /// `d_min_true` beside it because the third side is treated differently. Here the three `q_i`
    /// are the three pairs and there is only one number to report.
    pub d_min: T,
    pub steps: usize,
    pub finite: bool,
    pub budget_exhausted: bool,
    /// Running max of `|Gamma*|` against its largest term.
    pub gamma_max: T,
    /// Running max of `|sum q_i| / max|q_i|`, Heggie's Eq. (9). **No AZ analogue.**
    pub sum_q_max: T,
    /// Largest physical step actually taken, as an `s.t` difference across one step.
    ///
    /// The tripwire. On AZ a single step advanced the clock by `2.209e128` against an interval of
    /// `0.4` and the march recorded a clean landing, because the clamp corrects the *clock* and
    /// cannot un-take the step. Nothing recorded it until this field existed.
    pub dt_max: T,
    /// Steps whose `s.t` exceeded twice their own interval. Must be zero.
    pub n_overshoot: u32,
    /// `R1 R2 R3` hit the `T::TINY` floor, so `dtau` divided by a fabricated denominator and the
    /// step advanced anyway. An advance-anyway site with no telemetry is indistinguishable from
    /// one that never fires.
    pub r_floored: bool,
    /// Set if a collision was detected and `stop_on_collision` was on.
    pub collided: bool,
    /// Physical `d|q_i|/dt` at the end, all three pairs. Kept because it is free — the march
    /// needs it for the step limit anyway.
    pub q_dot_max: T,
}

/// `dt/dtau` under the time transformation in force: `R1 R2 R3` for Eq. (20), divided by
/// `S^{3/2}` for Eq. (22).
#[inline]
fn dt_dtau<T: Real>(s: &HgState<T>, time: HgTime) -> T {
    let rp = s.r_prod();
    match time {
        HgTime::Product => rp,
        HgTime::SumPow32 { .. } => {
            let ss = s.s().max(T::TINY);
            rp / (ss * ss.sqrt())
        }
    }
}

/// `dq_i/dt` from Heggie's Eq. (7): `q_i' = p_i/mu_i - p_k/m_j - p_j/m_k`.
///
/// The index pattern is the one place this is easy to get wrong — the momentum with the
/// **further** cyclic index carries the reciprocal of the **nearer** mass, and vice versa.
/// `tests/heggie_march.rs` finite-differences it against Eq. (6) rather than trusting it.
#[inline]
pub fn q_dot<T: Real>(sys: &HgSystem<T>, p: &[crate::Vec2<T>; 3]) -> [crate::Vec2<T>; 3] {
    std::array::from_fn(|i| {
        let (j, k) = cyc(i);
        p[i] / sys.mu[i] - p[k] * sys.inv_m[j] - p[j] * sys.inv_m[k]
    })
}

/// The predictive step limit, in `dtau` units.
///
/// `dt = (dt/dtau) dtau`, so a physical bound `dt <= f d_min / |dq/dt|_max` is
/// `dtau <= f d_min / (|dq/dt|_max (dt/dtau))`. Every input is already in hand: the three `q_i`
/// **are** the three pair separations, with no reference body making one of them special.
///
/// `+inf` when nothing is in force — the absence of a bound, not a step of zero.
#[inline]
fn predictive_dtau<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, f: T, time: HgTime) -> T {
    if f <= T::zero() {
        return T::infinity();
    }
    let (q, p) = sys.phys_from_state(s);
    let v = q_dot(sys, &p);
    let d_min = q.iter().fold(T::infinity(), |a, x| a.min(x.norm()));
    let v_max = v.iter().fold(T::zero(), |a, x| a.max(x.norm()));
    let denom = v_max * dt_dtau(s, time);
    if denom > T::zero() && d_min.is_finite() && denom.is_finite() {
        f * d_min / denom
    } else {
        T::infinity()
    }
}

/// The full entry point.
///
/// `n_sync` is a **sampling** cadence here and nothing more: the loop records residuals and the
/// closest approach at each boundary and does not re-register. That is the whole structural
/// difference from `integrate_az_opts`, and it is why this function is a fifth the length.
pub fn integrate_hg<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &HgOpts<T>,
) -> HgOut<T> {
    let mut sys = HgSystem::new(*m);
    sys.lc_stable = opts.lc_stable;

    // Registration. Once, at t = 0, and never again.
    let (mut s, h0) = sys.to_reg(&s0);
    let mut h = h0;

    let n = n_sync.max(1);
    let dt_sync = t_max / T::lit(n as f64);

    let mut out = HgOut {
        state: s0,
        t: T::zero(),
        drift: T::zero(),
        d_min: T::infinity(),
        steps: 0,
        finite: true,
        budget_exhausted: false,
        gamma_max: T::zero(),
        sum_q_max: T::zero(),
        dt_max: T::zero(),
        n_overshoot: 0,
        r_floored: false,
        collided: false,
        q_dot_max: T::zero(),
    };

    let mut t_now = T::zero();
    'outer: for _ in 0..n {
        let dt_left = dt_sync;
        let tol = if opts.clamp_final_step { dt_left * T::LAND_EPS_REL } else { T::zero() };
        s.t = T::zero();

        let entry = dt_dtau(&s, opts.time);
        if entry <= T::zero() || !entry.is_finite() {
            out.finite = false;
            break;
        }
        let dtau_entry = eta * dt_left / entry;

        loop {
            // `is_finite` explicitly: `NaN >= x` is false, so a diverged trajectory never
            // satisfies the loop exit and burns its entire step budget.
            if !s.is_finite() {
                out.finite = false;
                break 'outer;
            }
            if s.t >= dt_left - tol {
                break;
            }
            if out.steps >= max_steps {
                out.budget_exhausted = true;
                out.finite = false;
                break 'outer;
            }

            let rp_raw = s.r_prod();
            if s.r(0) < T::TINY || s.r(1) < T::TINY || s.r(2) < T::TINY {
                out.r_floored = true;
            }
            let d = dt_dtau(&s, opts.time).max(T::TINY);

            let mut dtau = match opts.dtau_mode {
                HgDtauMode::FixedPerInterval => dtau_entry,
                HgDtauMode::PerStepRemaining => eta * (dt_left - s.t).max(T::zero()) / d,
                HgDtauMode::PerStepInterval => (eta * dt_left / d).min(dtau_entry),
            };
            let lim = predictive_dtau(&sys, &s, opts.step_limit_f, opts.time);
            if lim < dtau {
                dtau = lim;
            }
            // Land ON the boundary. A one-sided reduction of the final step only, using the same
            // `d` the sizing used — recomputing it here would let the clamp and the step disagree.
            if opts.clamp_final_step {
                dtau = dtau.min((dt_left - s.t).max(T::zero()) / d);
            }
            let _ = rp_raw;

            let before = s.t;
            s = rk4::step(&sys, &s, h, opts.time, dtau);
            out.steps += 1;
            let took = s.t - before;
            if took > out.dt_max {
                out.dt_max = took;
            }
            if opts.clamp_final_step && s.t > dt_left * T::lit(2.0) {
                out.n_overshoot += 1;
                debug_assert!(
                    false,
                    "step overshot its interval: s.t = {} against dt_left = {}",
                    s.t, dt_left
                );
            }

            // Closest approach, sampled every step rather than at boundaries: a collision that
            // happens between two boundaries is not one that did not happen.
            let (q, p) = sys.phys_from_state(&s);
            let dm = q.iter().fold(T::infinity(), |a, x| a.min(x.norm()));
            if dm < out.d_min {
                out.d_min = dm;
            }
            let v = q_dot(&sys, &p);
            out.q_dot_max = out.q_dot_max.max(v.iter().fold(T::zero(), |a, x| a.max(x.norm())));
            if opts.stop_on_collision && opts.r_coll > T::zero() && dm < opts.r_coll {
                out.collided = true;
                t_now += s.t;
                break 'outer;
            }
        }

        t_now += s.t;

        // The boundary: sample, do not re-register.
        out.gamma_max = out.gamma_max.max(gamma_residual(&sys, &s, h));
        out.sum_q_max = out.sum_q_max.max(sys.sum_q_residual(&s));
        if opts.refresh_h_at_boundary {
            h = sys.energy_of(&s);
        }
    }

    out.t = t_now;
    if out.finite {
        out.state = sys.to_cartesian(&s);
        let e = sys.energy_of(&s);
        out.drift = ((e - h0) / h0.abs().max(T::DRIFT_FLOOR)).abs();
        out.finite = out.state.is_finite() && out.drift.is_finite();
    }
    out
}
