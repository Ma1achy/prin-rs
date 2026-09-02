//! **Chain coordinates** — Mikkola & Aarseth, CMDA **57** (1993) 439.
//!
//! Hold the system as inter-particle vectors along a chain rather than as three positions, so a
//! small separation is a **sum of small quantities** instead of a **difference of large ones**.
//! For three bodies the chain is two vectors:
//!
//! ```text
//!   X1 = r_b - r_a        X2 = r_c - r_b        r_c - r_a = X1 + X2
//! ```
//!
//! where `(a, b, c)` is the chain ordering. The third separation is never formed by subtracting
//! two positions, which is where the digits go on a wide configuration.
//!
//! # Why this exists here, and what it is aimed at
//!
//! `far` spans body positions to ~13 units where the latent charts sit at `R = 1`. AZ wins it at
//! f64 by 0.890 decades on every pixel and **loses it at f32 by 0.768 decades on every pixel** —
//! measured, guard LIVE at both precisions. That inversion says the `far` result is
//! precision-limited, and chain is the published repair for precisely that.
//!
//! It is a **round-off** fix, so it may be invisible at f64 and must be run at both precisions or
//! the result is uninterpretable.
//!
//! # The ordering is FROZEN, and that is not a detail
//!
//! Choosing a chain ordering is a **re-registration** — the same class of act as AZ picking a
//! reference body, which is the mechanism this whole logH investigation exists to isolate. A chain
//! that re-selects its ordering at sync boundaries reintroduces exactly what logH was built to
//! have none of, and the comparison stops being clean.
//!
//! So [`ChainOrder::select`] is called **once, at registration**, and the ordering is carried. The
//! re-selecting variant is a separate named arm and never the default — kept so the question can
//! be asked, not so it can be assumed.
//!
//! # This is a DIAGNOSTIC integrator, not an ensemble occupant
//!
//! It carries no events, no closure escape rule, no `t_end` and no outcome classification, so it
//! **cannot produce labels and its numbers are not comparable with the committed corpus**. It
//! exists to answer one question — does holding the coordinates as differences reduce round-off on
//! a wide configuration — and it is deliberately too small to be mistaken for an arm of the
//! gallery.

use crate::physics::{energy, Cart};
use crate::{Real, Vec2};

use super::hamiltonian::LhTime;

/// The chain ordering, as body indices `(a, b, c)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChainOrder(pub [usize; 3]);

impl ChainOrder {
    /// Pick the ordering: the **tightest pair adjacent**, with the remaining body on the end
    /// nearer to it. That is the 3-body case of Mikkola & Aarseth's rule — build the chain from
    /// the shortest inter-particle vector outward — and it is what puts the badly-conditioned
    /// separation on `X1` or `X2` rather than on the sum.
    ///
    /// **Called once, at registration.** See the module doc: re-selecting is a re-registration.
    pub fn select<T: Real>(r: &[Vec2<T>; 3]) -> Self {
        let d = [
            (r[1] - r[0]).norm(), // pair (0,1)
            (r[2] - r[0]).norm(), // pair (0,2)
            (r[2] - r[1]).norm(), // pair (1,2)
        ];
        let mut k = 0usize;
        for i in 1..3 {
            if d[i] < d[k] {
                k = i;
            }
        }
        let (a, b, c) = match k {
            0 => (0usize, 1usize, 2usize),
            1 => (0, 2, 1),
            _ => (1, 2, 0),
        };
        // Attach the third body to whichever end it is closer to, so the LONGER chain vector is
        // the one that spans the wide gap and the sum `X1 + X2` is never the tight pair.
        let da = (r[c] - r[a]).norm();
        let db = (r[c] - r[b]).norm();
        if da < db {
            ChainOrder([c, a, b])
        } else {
            ChainOrder([a, b, c])
        }
    }
}

/// Chain state: two inter-particle vectors, their rates, the interval clock and TTL's `W`.
///
/// **Ten components, not fourteen.** The centre of mass is absent by construction rather than
/// carried and ignored — a chain holds only differences, so there is no COM degree of freedom to
/// integrate, and three body positions cannot be recovered from it without one. `to_cart`
/// therefore returns a **COM-centred** configuration, which is the frame every comparison in this
/// project already uses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChainState<T> {
    /// `X1 = r_b - r_a`, `X2 = r_c - r_b`.
    pub x: [Vec2<T>; 2],
    /// `U1 = v_b - v_a`, `U2 = v_c - v_b`.
    pub u: [Vec2<T>; 2],
    pub t: T,
    pub w: T,
}

impl<T: Real> ChainState<T> {
    pub fn is_finite(&self) -> bool {
        self.x.iter().chain(self.u.iter()).all(|q| q.is_finite())
            && self.t.is_finite()
            && self.w.is_finite()
    }

    pub fn from_cart(c: &Cart<T>, o: ChainOrder) -> Self {
        let [a, b, k] = o.0;
        Self {
            x: [c.r[b] - c.r[a], c.r[k] - c.r[b]],
            u: [c.v[b] - c.v[a], c.v[k] - c.v[b]],
            t: T::zero(),
            w: T::zero(),
        }
    }

    /// Back to COM-centred Cartesian, given the masses in **body** order.
    pub fn to_cart(&self, m: &[T; 3], o: ChainOrder) -> Cart<T> {
        let [a, b, k] = o.0;
        let mut r = [Vec2::zero(); 3];
        let mut v = [Vec2::zero(); 3];
        // Place body `a` at the origin, then walk the chain. The COM is removed afterwards, so
        // the arbitrary anchor cancels exactly.
        r[a] = Vec2::zero();
        r[b] = self.x[0];
        r[k] = self.x[0] + self.x[1];
        v[a] = Vec2::zero();
        v[b] = self.u[0];
        v[k] = self.u[0] + self.u[1];
        let mt = m[0] + m[1] + m[2];
        let rc = (r[0] * m[0] + r[1] * m[1] + r[2] * m[2]) / mt;
        let vc = (v[0] * m[0] + v[1] * m[1] + v[2] * m[2]) / mt;
        for i in 0..3 {
            r[i] = r[i] - rc;
            v[i] = v[i] - vc;
        }
        Cart { r, v }
    }

    /// The three separations, **in body-pair order** `[(0,1), (0,2), (1,2)]` to match
    /// [`crate::physics::PAIRS`], computed from the chain vectors alone.
    ///
    /// This is the whole mechanism: the pair spanned by `X1 + X2` is formed by **adding** two
    /// chain vectors, never by subtracting two positions.
    pub fn seps(&self, o: ChainOrder) -> [T; 3] {
        let [a, b, k] = o.0;
        let d_ab = self.x[0].norm();
        let d_bk = self.x[1].norm();
        let d_ak = (self.x[0] + self.x[1]).norm();
        let mut out = [T::zero(); 3];
        let put = |out: &mut [T; 3], i: usize, j: usize, d: T| {
            let idx = match (i.min(j), i.max(j)) {
                (0, 1) => 0,
                (0, 2) => 1,
                _ => 2,
            };
            out[idx] = d;
        };
        put(&mut out, a, b, d_ab);
        put(&mut out, b, k, d_bk);
        put(&mut out, a, k, d_ak);
        out
    }
}

/// `d(X)/ds` and `d(U)/ds` under the time transformation.
///
/// The relative accelerations are formed from the chain vectors directly. Writing
/// `A1 = a_b - a_a` in terms of the separations avoids ever building an absolute position.
///
/// **One force evaluation**, and the caller counts it.
pub fn deriv<T: Real>(
    m: &[T; 3],
    s: &ChainState<T>,
    o: ChainOrder,
    b_const: T,
    time: LhTime,
) -> ChainState<T> {
    let [ia, ib, ik] = o.0;
    let x1 = s.x[0];
    let x2 = s.x[1];
    let x3 = x1 + x2; // r_k - r_a, by construction a SUM

    let g = |d: Vec2<T>| {
        let n = d.norm().max(T::TINY);
        d / (n * n * n)
    };
    let (g1, g2, g3) = (g(x1), g(x2), g(x3));

    // `a_b - a_a` and `a_k - a_b`, with `G = 1` as everywhere in this crate. Written out rather
    // than derived in place, because AZ-family algebra **fails silently** — a sign error here
    // produces trajectories that look like physics. `tests/chain_coords.rs` differences these
    // against `newton::accel` and carries a swapped-pair negative control, and the first cut of
    // this expression had the `g1` sign inverted in BOTH lines and was caught by it.
    //
    //   a_a = +m_b g1 + m_k g3
    //   a_b = -m_a g1 + m_k g2
    //   a_k = -m_a g3 - m_b g2
    let a1 = g1 * -(m[ia] + m[ib]) + g2 * m[ik] - g3 * m[ik];
    let a2 = g1 * m[ia] - g2 * (m[ib] + m[ik]) - g3 * m[ia];

    // Denominators. `U` and `K` are formed from the COM-centred reconstruction, which is exact
    // arithmetic on the chain vectors; the point of the chain is the SEPARATIONS, and those are
    // never reconstructed.
    let (drift, kick) = match time {
        LhTime::None => (T::one(), T::one()),
        LhTime::LogH => {
            let c = s.to_cart(m, o);
            let k = energy::kinetic(&c.v, m);
            let d = s.seps(o);
            let mut u = T::zero();
            for (idx, &(i, j)) in crate::physics::PAIRS.iter().enumerate() {
                u = u + m[i] * m[j] / d[idx].max(T::TINY);
            }
            (k + b_const, u)
        }
        // TTL's `Omega` is mass-independent and reads the chain separations directly.
        LhTime::Ttl => {
            let mbar = (m[0] + m[1] + m[2]) / T::lit(3.0);
            let wgt = mbar * mbar;
            let d = s.seps(o);
            let mut om = T::zero();
            for x in d.iter() {
                om = om + wgt / x.max(T::TINY);
            }
            (s.w, om)
        }
    };
    let inv_d = T::one() / drift.max(T::TINY);
    let inv_k = T::one() / kick.max(T::TINY);

    ChainState {
        x: [s.u[0] * inv_d, s.u[1] * inv_d],
        u: [a1 * inv_k, a2 * inv_k],
        t: inv_d,
        w: T::zero(),
    }
}

/// `self + h * d`, componentwise.
pub fn axpy<T: Real>(s: &ChainState<T>, h: T, d: &ChainState<T>) -> ChainState<T> {
    ChainState {
        x: [s.x[0] + d.x[0] * h, s.x[1] + d.x[1] * h],
        u: [s.u[0] + d.u[0] * h, s.u[1] + d.u[1] * h],
        t: s.t + d.t * h,
        w: s.w + d.w * h,
    }
}

/// One RK4 step. **Four force evaluations**, returned, and the accumulation order matches the
/// other three integrators so a difference is a difference of equations.
pub fn rk4<T: Real>(
    m: &[T; 3],
    s: &ChainState<T>,
    o: ChainOrder,
    b: T,
    time: LhTime,
    h: T,
) -> (ChainState<T>, usize) {
    let half = T::lit(0.5);
    let two = T::lit(2.0);
    let six = T::lit(6.0);
    let k1 = deriv(m, s, o, b, time);
    let k2 = deriv(m, &axpy(s, half * h, &k1), o, b, time);
    let k3 = deriv(m, &axpy(s, half * h, &k2), o, b, time);
    let k4 = deriv(m, &axpy(s, h, &k3), o, b, time);
    let acc = ChainState {
        x: [
            k1.x[0] + k2.x[0] * two + k3.x[0] * two + k4.x[0],
            k1.x[1] + k2.x[1] * two + k3.x[1] * two + k4.x[1],
        ],
        u: [
            k1.u[0] + k2.u[0] * two + k3.u[0] * two + k4.u[0],
            k1.u[1] + k2.u[1] * two + k3.u[1] * two + k4.u[1],
        ],
        t: k1.t + k2.t * two + k3.t * two + k4.t,
        w: k1.w + k2.w * two + k3.w * two + k4.w,
    };
    (axpy(s, h / six, &acc), 4)
}
