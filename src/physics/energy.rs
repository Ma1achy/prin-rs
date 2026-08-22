//! Energy, moment of inertia, and the hyperradius that sets canonical units.

use crate::physics::{PAIRS, G};
use crate::{Real, Vec2};

/// Total energy `T + U`, with the same softening convention as [`super::newton::accel`].
///
/// `ke = 0.5 sum_k m_k |v_k|^2`, `pe = -sum_pairs G m_i m_j / sqrt(|d|^2 + eps2)`.
pub fn energy<T: Real>(r: &[Vec2<T>; 3], v: &[Vec2<T>; 3], m: &[T; 3], eps2: T) -> T {
    let mut ke = T::zero();
    for k in 0..3 {
        ke += m[k] * v[k].norm_sq();
    }
    ke = ke * T::lit(0.5);

    let mut pe = T::zero();
    for &(i, j) in PAIRS.iter() {
        let d2 = (r[j] - r[i]).norm_sq() + eps2;
        pe -= T::lit(G) * m[i] * m[j] / d2.sqrt();
    }
    ke + pe
}

/// Centre of mass.
pub fn com<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3]) -> Vec2<T> {
    let mtot = m[0] + m[1] + m[2];
    (r[0] * m[0] + r[1] * m[1] + r[2] * m[2]) / mtot
}

/// Moment of inertia about the centre of mass: `I = sum_k m_k |r_k - R_com|^2`.
///
/// The reference's `tb_ftle.inertia` reads `tb.M` from module scope and silently ignores
/// any mass override passed to its caller. Here the masses are an argument.
pub fn inertia<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3]) -> T {
    let c = com(r, m);
    let mut i = T::zero();
    for k in 0..3 {
        i += m[k] * (r[k] - c).norm_sq();
    }
    i
}

/// Mass-weighted hyperradius `R = sqrt(I / M)`.
///
/// **Evaluated once at `t = 0` and never updated.** `r_coll` and `epsilon` are expressed as
/// fractions of this. A co-moving length makes the Hamiltonian time-dependent and destroys
/// energy conservation — measured `|dE/E| = 3.06e-02`, *identical* at `dt = 1e-4` and
/// `dt = 2e-5`, which is the signature of a wrong equation rather than an accuracy problem
/// (BRIEF §2.5).
pub fn hyperradius<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3]) -> T {
    let mtot = m[0] + m[1] + m[2];
    (inertia(r, m) / mtot).sqrt()
}

/// Crossing time `sqrt(R^3 / M)`. One time unit is roughly one crossing time.
pub fn crossing_time<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3]) -> T {
    let mtot = m[0] + m[1] + m[2];
    let rg = hyperradius(r, m);
    (rg.powf(T::lit(3.0)) / mtot).sqrt()
}
