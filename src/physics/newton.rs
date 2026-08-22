//! Softened Newtonian acceleration and pair separations.
//!
//! Note the softening argument is **eps squared**, matching `tb.accel(r, eps2)`. Softening
//! enters as `|d|^2 + eps2` inside both the `-3/2` power here and the `sqrt` in the
//! potential, so the force stays the exact gradient of the softened potential.
//!
//! Softening is a *different force law*: BRIEF §2.5 requires `eps > 0` data to be tagged
//! and never mixed with `eps = 0` data. The AZ path always passes `eps2 = 0` — removing the
//! singularity by coordinate change is the whole point of regularising.

use crate::physics::PAIRS;
use crate::{Real, Vec2};

/// `a_i = sum_{j != i} G m_j (r_j - r_i) / (|r_j - r_i|^2 + eps2)^{3/2}`
pub fn accel<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3], eps2: T) -> [Vec2<T>; 3] {
    let mut a = [Vec2::zero(); 3];
    for &(i, j) in PAIRS.iter() {
        let d = r[j] - r[i];
        let d2 = d.norm_sq() + eps2;
        let inv3 = d2.powf(T::lit(-1.5));
        let f = d * inv3;
        a[i] += f * m[j];
        a[j] -= f * m[i];
    }
    a
}

/// Unsoftened pair separations, in `PAIRS` order: `|r1-r0|, |r2-r0|, |r2-r1|`.
pub fn pair_dists<T: Real>(r: &[Vec2<T>; 3]) -> [T; 3] {
    let mut d = [T::zero(); 3];
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        d[k] = (r[j] - r[i]).norm();
    }
    d
}

/// Adaptive timestep: the shortest local two-body free-fall time over the three pairs.
///
/// `dt = eta * min_pairs( r_ij^{3/2} / sqrt(G (m_i + m_j)) )`
///
/// Scale-covariant by construction: under `r -> alpha r` this scales as `alpha^{3/2}`,
/// exactly as time does, so it introduces **no fixed length or time scale** (BRIEF §2.3).
/// That is what lets the project quotient out Newtonian scale invariance.
pub fn adaptive_dt<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3], eta: T) -> T {
    let d = pair_dists(r);
    let mut best = T::infinity();
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        let rij = d[k].max(T::TINY);
        let t = rij.powf(T::lit(1.5)) / (T::lit(crate::physics::G) * (m[i] + m[j])).sqrt();
        if t < best {
            best = t;
        }
    }
    eta * best
}
