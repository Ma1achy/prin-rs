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

/// Inverse LC map, **numerically stable**.
///
/// The reference always computes `u0 = sqrt((|rho| + rho.x)/2)` first and derives `u1` from
/// it. When `rho` points along negative x that sum cancels catastrophically, and the
/// division `u1 = rho.y/(2 u0)` then amplifies the damage — measured 1.9e-13 relative loss
/// at 179 degrees, 3.5e-9 at 179.99.
///
/// The fix is standard: compute whichever component is **larger** directly, and derive the
/// other from it. `u0` is larger when `rho.x >= 0`, `u1` when `rho.x < 0`. Both branches
/// satisfy `rho_of_u(u) == rho`; the branch choice is otherwise irrelevant since
/// `rho(u) = rho(-u)`.
///
/// This matters beyond precision. The reference's branch cut is fixed along negative x, so
/// its accuracy depends on the **absolute orientation** of a configuration in the coordinate
/// frame — the physics is rotationally invariant, that implementation of it is not. And
/// registration happens at every sync boundary, so the loss is injected `n_sync` times per
/// trajectory rather than once.
#[inline(always)]
pub fn u_of_rho<T: Real>(rho: Vec2<T>) -> Vec2<T> {
    let half = T::lit(0.5);
    let two = T::lit(2.0);
    let r = rho.norm();
    if rho.x >= T::zero() {
        let u0 = (half * (r + rho.x)).max(T::zero()).sqrt();
        if u0 > T::TINY {
            Vec2::new(u0, rho.y / (two * u0))
        } else {
            Vec2::new(u0, (half * (r - rho.x)).max(T::zero()).sqrt())
        }
    } else {
        let mag = (half * (r - rho.x)).max(T::zero()).sqrt();
        let u1 = if rho.y < T::zero() { -mag } else { mag };
        if u1.abs() > T::TINY {
            Vec2::new(rho.y / (two * u1), u1)
        } else {
            Vec2::new((half * (r + rho.x)).max(T::zero()).sqrt(), u1)
        }
    }
}

/// The reference's inverse LC map, transcribed exactly. Retained so the change above can be
/// measured as a single controlled variable against the cross-check.
#[inline(always)]
pub fn u_of_rho_reference<T: Real>(rho: Vec2<T>) -> Vec2<T> {
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
