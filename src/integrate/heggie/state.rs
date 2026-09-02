//! The regularised state: `[u0(2), p0(2), u1(2), p1(2), u2(2), p2(2), t(1)]`, thirteen numbers.
//!
//! `u[i]` is Heggie's 4-vector `Q_i` reduced to the plane (see the module doc); `p[i]` is his
//! `P_i`. The physical relative vectors `q_i` are never stored — they are `lc::rho_of_u(u[i])`,
//! which is Heggie's own identity `q_i = (1/2) A_i^T Q_i`.
//!
//! Indices are **0-based here and 1-based in the paper**. Heggie's `q_1 = q_2' - q_3'` is this
//! module's `u[0]`, the vector between bodies 1 and 2. Every cyclic formula is written with
//! `i`, `j = (i+1)%3`, `k = (i+2)%3`, which is his `i, i+1, i+2`.
//!
//! `h` is deliberately **not** part of the state, for the same reason `E` is not part of
//! `AzState`: it is the physical energy frozen at registration and enters `Gamma*` as a
//! constant. Evolving it would be a different system.

use crate::{Real, Vec2};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HgState<T> {
    /// Heggie's `Q_i`, planar.
    pub u: [Vec2<T>; 3],
    /// Heggie's `P_i`, planar.
    pub p: [Vec2<T>; 3],
    /// Accumulated **physical** time since registration. `dt/dtau` depends on the time
    /// transformation in force — `R1 R2 R3` under Eq. (20), divided by `S^{3/2}` under Eq. (22).
    pub t: T,
}

impl<T: Real> HgState<T> {
    /// `R_i = Q_i^T Q_i = |q_i|`. Deliberately **unfloored**, matching `AzState::a`/`b`: the
    /// floor belongs at the sites that divide, not at the site that reports.
    #[inline(always)]
    pub fn r(&self, i: usize) -> T {
        self.u[i].norm_sq()
    }

    /// `S = R1 + R2 + R3`, the sum that appears in Heggie's Eq. (22) time transformation.
    #[inline(always)]
    pub fn s(&self) -> T {
        self.r(0) + self.r(1) + self.r(2)
    }

    /// `R1 R2 R3`, `dt/dtau` under Eq. (20).
    #[inline(always)]
    pub fn r_prod(&self) -> T {
        self.r(0) * self.r(1) * self.r(2)
    }

    pub fn is_finite(&self) -> bool {
        self.u.iter().chain(self.p.iter()).all(|v| v.is_finite()) && self.t.is_finite()
    }

    /// Flat form, so the finite-difference test can perturb component `k` generically while the
    /// physics stays written as vector algebra. Layout is `u0 p0 u1 p1 u2 p2 t`.
    pub fn to_array13(&self) -> [T; 13] {
        [
            self.u[0].x, self.u[0].y, self.p[0].x, self.p[0].y,
            self.u[1].x, self.u[1].y, self.p[1].x, self.p[1].y,
            self.u[2].x, self.u[2].y, self.p[2].x, self.p[2].y,
            self.t,
        ]
    }

    pub fn from_array13(a: [T; 13]) -> Self {
        Self {
            u: [Vec2::new(a[0], a[1]), Vec2::new(a[4], a[5]), Vec2::new(a[8], a[9])],
            p: [Vec2::new(a[2], a[3]), Vec2::new(a[6], a[7]), Vec2::new(a[10], a[11])],
            t: a[12],
        }
    }

    /// `self + h * d`, componentwise. The RK4 stages and nothing else.
    pub fn axpy(&self, h: T, d: &Self) -> Self {
        let mut out = *self;
        for i in 0..3 {
            out.u[i] = self.u[i] + d.u[i] * h;
            out.p[i] = self.p[i] + d.p[i] * h;
        }
        out.t = self.t + d.t * h;
        out
    }
}
