//! The logH state: `[r0(2), v0(2), r1(2), v1(2), r2(2), v2(2), t(1)]`, thirteen numbers.
//!
//! **These are ordinary Cartesian positions and velocities.** That is the whole point of the
//! method and the reason this file is short: there is no chart, no regularised variable, and
//! nothing to transform into or out of. `to_reg` has no analogue here, and neither does
//! `phys_from_state` — the state *is* the physical state, at every step, and a boundary sample
//! costs a read rather than a round trip.
//!
//! Thirteen components, deliberately the same count as [`HgState`](crate::integrate::heggie::
//! HgState), so the finite-difference harness perturbs component `k` the same generic way in
//! both modules and the two tests are the same test over different algebra.
//!
//! `B` is **not** part of the state, for the same reason `h` is not part of `HgState` and `E` is
//! not part of `AzState`: it is frozen at registration and enters `Lambda` as a constant. Under
//! autonomous Newtonian gravity `dB/ds` is identically zero — the `(s/U) dU/dt` term in
//! Mikkola's general form is non-zero only for velocity-dependent or explicitly time-dependent
//! forces — so evolving it here would be evolving a quantity whose derivative is zero, which is
//! a slower way of holding it constant and a place for a bug to hide.

use crate::{Real, Vec2};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LhState<T> {
    /// Cartesian positions.
    pub r: [Vec2<T>; 3],
    /// Cartesian velocities. Momenta are `m_i v_i`; the mass factor lives in the system.
    pub v: [Vec2<T>; 3],
    /// Accumulated **physical** time since the start of the current interval. `dt/ds` is
    /// `1/(T + B)` under [`LhTime::LogH`](super::hamiltonian::LhTime::LogH) and exactly `1`
    /// under [`LhTime::None`](super::hamiltonian::LhTime::None).
    pub t: T,
    /// **TTL's auxiliary variable.** Under [`LhTime::Ttl`](super::hamiltonian::LhTime::Ttl) this
    /// is integrated alongside the state and used as the *drift* denominator in place of `K + B`.
    ///
    /// It exists because TTL's time transformation is driven by a function `Omega(r)` of the
    /// coordinates alone, which has no conserved partner to read off the state the way logH reads
    /// `K + B`. So `W` is carried, initialised to `Omega(r_0)`, and advanced by
    /// `dW = dOmega/dt * dt` — it tracks `Omega` along the flow and diverges from it exactly to
    /// the extent the integration has erred, which is the same off-shell asymmetry that makes the
    /// logH leapfrog work.
    ///
    /// **Unused and held constant under `LogH` and `None`.** Not `Option`: a branch in the inner
    /// loop to skip one addition would cost more than the addition, and a state whose component
    /// count depends on a mode is a state two code paths can disagree about.
    pub w: T,
}

impl<T: Real> LhState<T> {
    pub fn is_finite(&self) -> bool {
        self.r.iter().chain(self.v.iter()).all(|q| q.is_finite())
            && self.t.is_finite()
            && self.w.is_finite()
    }

    /// Flat form, so the finite-difference test can perturb component `k` generically while the
    /// physics stays written as vector algebra. Layout is `r0 v0 r1 v1 r2 v2 t w`.
    ///
    /// **Fourteen, not thirteen.** It matched `HgState`'s count until TTL needed `W`; that
    /// coincidence was never load-bearing and the name now says the length so a mismatch is a
    /// compile error rather than a silent index shift.
    pub fn to_array14(&self) -> [T; 14] {
        [
            self.r[0].x, self.r[0].y, self.v[0].x, self.v[0].y,
            self.r[1].x, self.r[1].y, self.v[1].x, self.v[1].y,
            self.r[2].x, self.r[2].y, self.v[2].x, self.v[2].y,
            self.t, self.w,
        ]
    }

    pub fn from_array14(a: [T; 14]) -> Self {
        Self {
            r: [Vec2::new(a[0], a[1]), Vec2::new(a[4], a[5]), Vec2::new(a[8], a[9])],
            v: [Vec2::new(a[2], a[3]), Vec2::new(a[6], a[7]), Vec2::new(a[10], a[11])],
            t: a[12],
            w: a[13],
        }
    }

    /// `self + h * d`, componentwise. The RK4 stages and nothing else.
    pub fn axpy(&self, h: T, d: &Self) -> Self {
        let mut out = *self;
        for i in 0..3 {
            out.r[i] = self.r[i] + d.r[i] * h;
            out.v[i] = self.v[i] + d.v[i] * h;
        }
        out.t = self.t + d.t * h;
        out.w = self.w + d.w * h;
        out
    }

    pub fn from_cart(c: &crate::physics::Cart<T>) -> Self {
        Self { r: c.r, v: c.v, t: T::zero(), w: T::zero() }
    }

    pub fn to_cart(&self) -> crate::physics::Cart<T> {
        crate::physics::Cart::new(self.r, self.v)
    }
}
