//! The logarithmic Hamiltonian and its equations of motion, kept adjacent as in the other two
//! integrator modules.
//!
//! Mikkola & Tanikawa, *Algorithmic regularization of the few-body problem*, MNRAS **310** (1999)
//! 745; and independently Preto & Tremaine, AJ **118** (1999) 2532. Implementation form follows
//! Mikkola & Merritt, AJ **135** (2008) 2398 (AR-CHAIN).
//!
//! ```text
//!   U   = +sum_pairs G m_i m_j / |r_i - r_j|        > 0
//!   K   = 0.5 sum_i m_i |v_i|^2
//!   E   = K - U                                     the physical energy
//!   B   = U - K = -E                                frozen at registration, CONSTANT
//!
//!   Lambda = ln(K + B) - ln(U)
//! ```
//!
//! `K` is written for the kinetic energy throughout this module rather than `T`, because `T` is
//! the generic scalar parameter everywhere else in this crate and one letter cannot be both.
//!
//! # Why there is no coordinate transformation here
//!
//! This is *algorithmic* regularisation: a time transformation and a good integrator, and
//! nothing else. There is no chart, no reference body, no chain, and therefore **nothing to
//! re-select** — which is the property this module exists to test. AZ re-registers at every sync
//! boundary and Heggie never does; logH has no registration to perform even once.
//!
//! # On shell the two denominators coincide, and that is load-bearing twice
//!
//! `K - U + B = 0` on the solution path, so `K + B` and `U` are **the same number** there. Off
//! shell they are not, and the difference is the whole mechanism:
//!
//! - It is why a leapfrog on `Lambda` is exact for the two-body problem. The drift and the kick
//!   are evaluated at different points, so they see the two denominators differing, and the
//!   asymmetry is what corrects the step.
//! - It is why **the most plausible transcription error is invisible on shell.** Swapping which
//!   denominator each half uses changes nothing along an exact solution. So the finite-difference
//!   test draws `B` *independently of the state* — the same reason the `Gamma*` test uses a
//!   random off-shell `h` — and [`Dens`] exists as a named pair so the swap can be constructed
//!   and asserted to fire.
//!
//! # The residual is the energy error, and saying otherwise would be overselling it
//!
//! `rho = |K + B - U| / U` is available every step for free. It is **not** an independent
//! constraint: substituting `B = U_0 - K_0` gives `K + B - U = E(t) - E(0)` exactly, so `rho` is
//! the absolute energy defect normalised by `U` instead of by `|E_0|`. That normalisation is the
//! natural one here — `U` is what sits under the kick — and it is finite where `E_0` is near
//! zero, which is worth having. But it measures energy conservation and nothing else.
//!
//! **logH has no analogue of Heggie's `sum q_i = 0`.** That one is a genuine second instrument
//! because Heggie's enlarged phase space carries a constraint; logH's phase space is the
//! physical one and has nothing to constrain. Stated here so the comparison table does not read
//! as though a check went missing.

use crate::physics::{energy, newton};
use crate::{Real, Vec2};

use super::state::LhState;
use super::system::LhSystem;

/// Whether the time transformation is in force.
///
/// `None` is the **control**, and it is deliberately the same code path rather than a separate
/// integrator: with both denominators set to one, `Lambda` degenerates to the physical
/// Hamiltonian `K - U` and `dt/ds = 1`, so the march becomes ordinary unregularised Cartesian
/// dynamics under whichever stepper is in force. The control is then literally the
/// regularisation removed and nothing else — it inherits the same event sampling, the same
/// closure escape rule, the same boundary cadence and the same `t_end`, so its labels are
/// comparable with every other occupant's.
///
/// `src/integrate/leapfrog.rs` is the older, standalone unregularised KDK. It returns a bare
/// `TrajOut` with no outcome machinery at all, which is why it cannot serve as this control.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LhTime {
    /// `Lambda = ln(K + B) - ln(U)`; `dt = ds/(K+B)` on the drift, `ds/U` on the kick.
    #[default]
    LogH,
    /// No transformation: both denominators are `1`, `dt = ds`, and the generator is `K - U`.
    None,
    /// **Time-transformed leapfrog** (Mikkola & Aarseth, CMDA **84** (2002) 343).
    ///
    /// `dt = ds/W` on the drift and `ds/Omega(r)` on the kick, with `W` carried in the state and
    /// advanced by `dW = (dOmega/dt) dt`.
    ///
    /// # Why it exists, in one sentence
    ///
    /// logH's kick denominator is `U = sum m_i m_j / r_ij`, which is **mass-weighted**: a close
    /// approach between a heavy body and a light one barely moves `U`, so the physical step fails
    /// to shrink when it should. TTL replaces it with `Omega = sum w_ij / r_ij` for freely chosen
    /// weights, and this port takes `w_ij = mbar^2` with `mbar` the mean mass.
    ///
    /// # That weight choice makes the equal-mass control EXACT, which is the point
    ///
    /// At `m_0 = m_1 = m_2 = m` we have `mbar^2 = m^2 = m_i m_j` for every pair, so
    /// `Omega === U` identically and TTL degenerates to logH. The control is then an algebraic
    /// identity rather than a near-identity, and a TTL arm that differs on an equal-mass slice
    /// is measuring something other than the mass ratio. `w_ij = 1` would have been simpler and
    /// would have left `Omega` on a different scale from `U`, so a comparison at fixed `eta`
    /// would have scored the step size instead of the transformation.
    Ttl,
}

/// The two time-transformation denominators, as a named pair.
///
/// They are separated rather than folded into [`deriv`] for one reason: on shell they are equal,
/// so a test that cannot construct them independently cannot detect a swap. `deriv_with` takes
/// this by value precisely so `tests/logh_hamiltonian_fd.rs` can hand it a swapped pair and
/// assert the finite differences disagree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dens<T> {
    /// `K + B` under [`LhTime::LogH`], `1` under [`LhTime::None`]. Divides the drift and `dt`.
    pub drift: T,
    /// `U` under [`LhTime::LogH`], `1` under [`LhTime::None`]. Divides the kick.
    pub kick: T,
}

impl<T: Real> Dens<T> {
    /// True if either denominator was non-positive or non-finite before flooring.
    ///
    /// `K + B` is `U > 0` on shell, so a non-positive value means the trajectory has wandered far
    /// enough off shell that the time transformation itself has failed — an *advance-anyway*
    /// site in the sense of `results/saturation/README.md`, not a terminal one, and therefore one
    /// that has to be recorded or it is invisible.
    pub fn degenerate(&self) -> bool {
        !(self.drift > T::zero()) || !(self.kick > T::zero())
    }
}

/// TTL's pair weight, `mbar^2` with `mbar` the mean mass. See [`LhTime::Ttl`] for why this and
/// not `1`: it makes `Omega === U` at equal masses, so the control is an identity.
#[inline]
pub fn ttl_weight<T: Real>(sys: &LhSystem<T>) -> T {
    let mbar = (sys.masses[0] + sys.masses[1] + sys.masses[2]) / T::lit(3.0);
    mbar * mbar
}

/// `Omega = w sum_pairs 1/|r_i - r_j|`, the **mass-independent** time-transformation function.
pub fn omega<T: Real>(sys: &LhSystem<T>, r: &[Vec2<T>; 3]) -> T {
    let w = ttl_weight(sys);
    let d = newton::pair_dists(r);
    let mut acc = T::zero();
    for x in d.iter() {
        acc = acc + w / x.max(T::TINY);
    }
    acc
}

/// `dOmega/dt = -w sum_pairs (r_ij . v_ij) / |r_ij|^3`, exact rather than differenced.
///
/// This is what advances `W`, so an error here is an error in the time transformation itself and
/// not merely in a diagnostic. `tests/logh_ttl.rs` finite-differences it against [`omega`].
pub fn omega_dot<T: Real>(sys: &LhSystem<T>, r: &[Vec2<T>; 3], v: &[Vec2<T>; 3]) -> T {
    let w = ttl_weight(sys);
    let mut acc = T::zero();
    for &(i, j) in crate::physics::PAIRS.iter() {
        let dr = r[j] - r[i];
        let dv = v[j] - v[i];
        let d = dr.norm().max(T::TINY);
        acc = acc - w * (dr.x * dv.x + dr.y * dv.y) / (d * d * d);
    }
    acc
}

/// The denominators for a state, before flooring.
pub fn denominators<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T, time: LhTime) -> Dens<T> {
    match time {
        LhTime::None => Dens { drift: T::one(), kick: T::one() },
        LhTime::LogH => Dens {
            drift: energy::kinetic(&s.v, &sys.masses) + b,
            kick: energy::potential_pos(&s.r, &sys.masses, T::zero()),
        },
        // The carried `W` drifts and the coordinate function `Omega` kicks. They are equal at
        // registration and diverge off shell -- the same asymmetry that makes the logH leapfrog
        // exact for two bodies, with `W` playing the part `K + B` plays there.
        LhTime::Ttl => Dens { drift: s.w, kick: omega(sys, &s.r) },
    }
}

/// `Lambda` itself. Only the finite-difference test calls this; the march never needs it.
///
/// Returns `T - U` under [`LhTime::None`] so the FD harness exercises both arms through one
/// path — the control arm is what says the harness is right before the LogH number is read.
pub fn lambda<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T, time: LhTime) -> T {
    let k = energy::kinetic(&s.v, &sys.masses);
    let u = energy::potential_pos(&s.r, &sys.masses, T::zero());
    match time {
        LhTime::None => k - u,
        LhTime::LogH => (k + b).ln() - u.ln(),
        // Written for completeness of the FD harness. TTL's generator is not a function of the
        // state alone -- it carries `W` -- so this is `ln W - ln Omega`, the direct analogue.
        LhTime::Ttl => s.w.max(T::TINY).ln() - omega(sys, &s.r).max(T::TINY).ln(),
    }
}

/// `rho = |K + B - U| / U`, the energy defect normalised by the potential.
///
/// See the module doc: this **is** `|E(t) - E(0)| / U`, not an independent residual.
pub fn residual<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T) -> T {
    let k = energy::kinetic(&s.v, &sys.masses);
    let u = energy::potential_pos(&s.r, &sys.masses, T::zero());
    if !(u > T::zero()) || !u.is_finite() {
        return T::infinity();
    }
    ((k + b - u) / u).abs()
}

/// The equations of motion, given denominators.
///
/// ```text
///   dr_i/ds = v_i / drift        dv_i/ds = a_i / kick        dt/ds = 1 / drift
/// ```
///
/// Each is the corresponding Hamiltonian derivative: `dr_i/ds = dLambda/dp_i` is
/// `(1/(K+B)) m_i v_i / m_i`, and `dv_i/ds = -(1/m_i) dLambda/dr_i` is `(1/(m_i U)) dU/dr_i`,
/// which is `a_i / U` because `dU/dr_i = m_i a_i`. That last identity is where the sign
/// convention in [`energy::potential_pos`] earns its name.
///
/// **One force evaluation**, and the caller counts it. Nothing here infers a count from a step
/// number: RK4 spends four of these per step and KDK one, so a derived `steps * k` stops being
/// true the moment two steppers share a table.
pub fn deriv_with<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, d: Dens<T>) -> LhState<T> {
    let a = newton::accel(&s.r, &sys.masses, T::zero());
    let inv_drift = T::one() / d.drift.max(T::TINY);
    let inv_kick = T::one() / d.kick.max(T::TINY);
    let mut out =
        LhState { r: [Vec2::zero(); 3], v: [Vec2::zero(); 3], t: inv_drift, w: T::zero() };
    for i in 0..3 {
        out.r[i] = s.v[i] * inv_drift;
        out.v[i] = a[i] * inv_kick;
    }
    out
}

/// [`deriv_with`] at the denominators implied by `time`. One force evaluation.
///
/// **`dW/ds` is set here and not in [`deriv_with`]**, because it is the one component whose rate
/// depends on the *mode* rather than on the denominators. `deriv_with` takes a bare [`Dens`] so
/// the finite-difference test can hand it a swapped pair; giving it a mode as well would give the
/// swap somewhere to hide. Under `LogH` and `None` the rate is exactly zero, so `W` is inert and
/// carried rather than branched on.
pub fn deriv<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T, time: LhTime) -> LhState<T> {
    let d = denominators(sys, s, b, time);
    let mut out = deriv_with(sys, s, d);
    if time == LhTime::Ttl {
        // `dW = (dOmega/dt) dt` and `dt/ds = 1/Omega` on the kick, so `dW/ds = omega_dot/Omega`.
        out.w = omega_dot(sys, &s.r, &s.v) / d.kick.max(T::TINY);
    }
    out
}
