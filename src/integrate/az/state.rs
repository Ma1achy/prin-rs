//! The regularised state: `[u1(2), p1(2), u2(2), p2(2), t(1)]`, nine numbers.
//!
//! `E` is deliberately **not** part of the state. It is the physical energy at the moment of
//! registration, frozen for the whole sub-interval, and it enters the regularised
//! Hamiltonian as a constant. Evolving it would be a different (and wrong) system.

use crate::{Real, Vec2};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AzState<T> {
    pub u1: Vec2<T>,
    pub p1: Vec2<T>,
    pub u2: Vec2<T>,
    pub p2: Vec2<T>,
    /// Accumulated **physical** time since registration, since `dt/dtau = A*B`.
    pub t: T,
}

impl<T: Real> AzState<T> {
    /// `A = |u1|^2 = |R1|`. Deliberately unfloored — the reference floors `A` and `B` in
    /// `phys_from_state` but **not** in `deriv`. Asymmetric, and transcribed as-is.
    #[inline(always)]
    pub fn a(&self) -> T {
        self.u1.norm_sq()
    }

    /// `B = |u2|^2 = |R2|`. Unfloored, as above.
    #[inline(always)]
    pub fn b(&self) -> T {
        self.u2.norm_sq()
    }

    pub fn is_finite(&self) -> bool {
        self.u1.is_finite() && self.p1.is_finite() && self.u2.is_finite() && self.p2.is_finite()
            && self.t.is_finite()
    }

    /// Flat form, so the finite-difference test can perturb component `k` generically while
    /// the physics stays written as vector algebra.
    pub fn to_array9(&self) -> [T; 9] {
        [
            self.u1.x, self.u1.y, self.p1.x, self.p1.y,
            self.u2.x, self.u2.y, self.p2.x, self.p2.y,
            self.t,
        ]
    }

    pub fn from_array9(a: [T; 9]) -> Self {
        Self {
            u1: Vec2::new(a[0], a[1]),
            p1: Vec2::new(a[2], a[3]),
            u2: Vec2::new(a[4], a[5]),
            p2: Vec2::new(a[6], a[7]),
            t: a[8],
        }
    }

    /// `self + h * d`, componentwise. The RK4 stages and nothing else.
    pub fn axpy(&self, h: T, d: &Self) -> Self {
        Self {
            u1: self.u1 + d.u1 * h,
            p1: self.p1 + d.p1 * h,
            u2: self.u2 + d.u2 * h,
            p2: self.p2 + d.p2 * h,
            t: self.t + d.t * h,
        }
    }
}
