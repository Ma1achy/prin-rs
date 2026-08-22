//! The AZ system for one choice of reference body.
//!
//! Reference body `a`; regularised pairs `(a,b)` and `(a,c)`; unregularised side `(b,c)`.
//! Because `(b,c)` is chosen to be the longest side, `|R3| >= max(|R1|,|R2|)`, so `R3 -> 0`
//! only in a genuine triple collision — which is provably non-regularisable anyway.
//!
//! Transcribed from `reference/tb_az.py:AZSystem`.

use crate::physics::Cart;
use crate::{Real, Vec2};

use super::lc;
use super::state::AzState;

#[derive(Clone, Copy, Debug)]
pub struct AzSystem<T> {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub ma: T,
    pub mb: T,
    pub mc: T,
    pub masses: [T; 3],
    pub mtot: T,
    pub mu1: T,
    pub mu2: T,
    /// Velocity-space mass matrix for `(Rdot1, Rdot2)`; its determinant is `ma*mb*mc/M`.
    pub k11: T,
    pub k12: T,
    pub k22: T,
}

impl<T: Real> AzSystem<T> {
    pub fn new(a: usize, b: usize, c: usize, masses: [T; 3]) -> Self {
        let (ma, mb, mc) = (masses[a], masses[b], masses[c]);
        let mtot = masses[0] + masses[1] + masses[2];
        Self {
            a,
            b,
            c,
            ma,
            mb,
            mc,
            masses,
            mtot,
            mu1: ma * mb / (ma + mb),
            mu2: ma * mc / (ma + mc),
            k11: mb * (ma + mc) / mtot,
            k22: mc * (ma + mb) / mtot,
            k12: -(mb * mc) / mtot,
        }
    }

    /// Cartesian -> regularised. Returns the state and the frozen energy `E`.
    pub fn to_reg(&self, s: &Cart<T>) -> (AzState<T>, T) {
        let (a, b, c) = (self.a, self.b, self.c);
        let r1 = s.r[b] - s.r[a];
        let r2 = s.r[c] - s.r[a];
        let v1 = s.v[b] - s.v[a];
        let v2 = s.v[c] - s.v[a];
        let p1v = v1 * self.k11 + v2 * self.k12;
        let p2v = v1 * self.k12 + v2 * self.k22;

        let u1 = lc::u_of_rho(r1);
        let u2 = lc::u_of_rho(r2);
        let two = T::lit(2.0);
        let st = AzState {
            u1,
            p1: lc::lt_apply(u1, p1v) * two,
            u2,
            p2: lc::lt_apply(u2, p2v) * two,
            t: T::zero(),
        };
        (st, self.energy_phys(r1, r2, p1v, p2v))
    }

    /// Regularised -> physical relative coordinates and momenta.
    ///
    /// Here `A` and `B` **are** floored, unlike in `deriv`. The asymmetry is the
    /// reference's; keeping it is what makes the cross-check an equality.
    pub fn phys_from_state(&self, s: &AzState<T>) -> (Vec2<T>, Vec2<T>, Vec2<T>, Vec2<T>) {
        let a = s.a().max(T::TINY);
        let b = s.b().max(T::TINY);
        let two = T::lit(2.0);
        (
            lc::rho_of_u(s.u1),
            lc::rho_of_u(s.u2),
            lc::l_apply(s.u1, s.p1) / (two * a),
            lc::l_apply(s.u2, s.p2) / (two * b),
        )
    }

    /// Physical energy in the relative coordinates.
    ///
    /// `R1, R2` are relative to a common body, **not** Jacobi, so the kinetic energy carries
    /// the cross term `P1.P2/ma`. Dropping it is an easy and silent error.
    pub fn energy_phys(&self, r1: Vec2<T>, r2: Vec2<T>, p1: Vec2<T>, p2: Vec2<T>) -> T {
        let d1 = r1.norm().max(T::TINY);
        let d2 = r2.norm().max(T::TINY);
        let d3 = (r2 - r1).norm().max(T::TINY);
        let two = T::lit(2.0);
        let kin = p1.norm_sq() / (two * self.mu1) + p1.dot(p2) / self.ma
            + p2.norm_sq() / (two * self.mu2);
        let pot = self.ma * self.mb / d1 + self.ma * self.mc / d2 + self.mb * self.mc / d3;
        kin - pot
    }

    /// Energy directly from a regularised state.
    pub fn energy_of(&self, s: &AzState<T>) -> T {
        let (r1, r2, p1, p2) = self.phys_from_state(s);
        self.energy_phys(r1, r2, p1, p2)
    }

    /// Regularised -> Cartesian, in the COM frame.
    pub fn to_cartesian(&self, s: &AzState<T>) -> Cart<T> {
        let (r1, r2, p1, p2) = self.phys_from_state(s);
        let det = self.k11 * self.k22 - self.k12 * self.k12;
        let v1 = (p1 * self.k22 - p2 * self.k12) / det;
        let v2 = (p2 * self.k11 - p1 * self.k12) / det;

        let mut r = [Vec2::zero(); 3];
        let mut v = [Vec2::zero(); 3];
        r[self.b] = r1;
        r[self.c] = r2;
        v[self.b] = v1;
        v[self.c] = v2;

        let mut rc = Vec2::zero();
        let mut vc = Vec2::zero();
        for k in 0..3 {
            rc += r[k] * self.masses[k];
            vc += v[k] * self.masses[k];
        }
        rc = rc / self.mtot;
        vc = vc / self.mtot;
        for k in 0..3 {
            r[k] -= rc;
            v[k] -= vc;
        }
        Cart::new(r, v)
    }
}
