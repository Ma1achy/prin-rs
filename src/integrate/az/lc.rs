//! Levi-Civita primitives, transcribed from `reference/tb_lc.py`.
//!
//! `L(u)` is the LC matrix `[[u.x, -u.y], [u.y, u.x]]`. Two properties are used throughout
//! and both matter for the derivative algebra:
//!
//! - `L(u)^T L(u) = |u|^2 I`, which is what inverts `p = 2 L(u)^T P`.
//! - `L(u) w = L(w) u` — bilinear and symmetric in its two arguments. This is why
//!   `dGamma/du1` of the cross term `(L(u1)p1).(L(u2)p2)` is `L(p1)^T (L(u2)p2)`, with `p1`
//!   in the *matrix* slot. It looks like a typo in the reference. It is not.

use crate::{Real, Vec2};

/// `L(u) @ w`
#[inline(always)]
pub fn l_apply<T: Real>(u: Vec2<T>, w: Vec2<T>) -> Vec2<T> {
    Vec2::new(u.x * w.x - u.y * w.y, u.y * w.x + u.x * w.y)
}

/// `L(u)^T @ w`
#[inline(always)]
pub fn lt_apply<T: Real>(u: Vec2<T>, w: Vec2<T>) -> Vec2<T> {
    Vec2::new(u.x * w.x + u.y * w.y, -u.y * w.x + u.x * w.y)
}

/// `rho = u^2` as a complex square. `|rho| = |u|^2`.
#[inline(always)]
pub fn rho_of_u<T: Real>(u: Vec2<T>) -> Vec2<T> {
    Vec2::new(u.x * u.x - u.y * u.y, T::lit(2.0) * u.x * u.y)
}

/// Inverse LC map. The branch choice is irrelevant: `rho(u) = rho(-u)`.
#[inline(always)]
pub fn u_of_rho<T: Real>(rho: Vec2<T>) -> Vec2<T> {
    let half = T::lit(0.5);
    let r = rho.norm();
    let u0 = (half * (r + rho.x)).max(T::zero()).sqrt();
    let u1 = if u0 > T::TINY {
        rho.y / (T::lit(2.0) * u0.max(T::TINY))
    } else {
        (half * (r - rho.x)).max(T::zero()).sqrt()
    };
    Vec2::new(u0, u1)
}
