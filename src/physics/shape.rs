//! The shape sphere: `shape_vec` maps a configuration to a point on the unit 2-sphere,
//! quotienting out translation, rotation and scale so only the *shape* of the triangle
//! remains. Mass-weighted Jacobi coordinates, then the Hopf map (BRIEF §4).
//!
//! Transcribed from `reference/refine_test.py:shape_vec`. The Jacobi pair is hard-wired as
//! (0,1) with body 2 outer, matching the reference.

use crate::{Real, Vec2};

/// Unit vector on the shape sphere.
///
/// `I = a + b` is the mass-weighted moment of inertia and is **not floored** in the
/// reference, so a triple collision gives 0/0. Transcribed as-is: a NaN here is a genuine
/// measurement outcome ("this could not be determined"), not missing data, and the caller
/// records it rather than discarding the copy.
pub fn shape_vec<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3]) -> [T; 3] {
    let (m0, m1, m2) = (m[0], m[1], m[2]);
    let mtot = m0 + m1 + m2;

    let rho = r[1] - r[0];
    let com01 = (r[0] * m0 + r[1] * m1) / (m0 + m1);
    let lam = r[2] - com01;

    let mu_rho = m0 * m1 / (m0 + m1);
    let mu_lam = m2 * (m0 + m1) / mtot;

    let rt = rho * mu_rho.sqrt();
    let lt = lam * mu_lam.sqrt();

    let a = rt.norm_sq();
    let b = lt.norm_sq();
    let i = a + b;

    let p = rt.x * lt.x + rt.y * lt.y;
    let q = rt.y * lt.x - rt.x * lt.y;

    let two = T::lit(2.0);
    let n = [(a - b) / i, two * p / i, two * q / i];

    // Algebraically |n| == 1 already; the renormalisation only damps roundoff. Kept
    // because the reference keeps it, and dropping it would move the last bits.
    let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    [n[0] / norm, n[1] / norm, n[2] / norm]
}

/// Spherical variance of a set of unit vectors, `1 - |mean|`, in `[0, 1]`.
///
/// This is `refine_test.svar` — the statistic the reference actually uses. It is **not**
/// what BRIEF §4 defines `spread_shape` to be (mean distance from the centroid, halved).
/// Both are computed and dumped: the brief's is the spec, this one is the one with a
/// reference, and reporting both makes the discrepancy measurable rather than a silent
/// choice.
pub fn svar<T: Real>(n: &[[T; 3]]) -> T {
    if n.len() < 2 {
        return T::zero();
    }
    let cnt = T::lit(n.len() as f64);
    let mut mean = [T::zero(); 3];
    for v in n {
        for k in 0..3 {
            mean[k] += v[k];
        }
    }
    for k in 0..3 {
        mean[k] = mean[k] / cnt;
    }
    let norm = (mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]).sqrt();
    T::one() - norm
}

/// BRIEF §4's `spread_shape`: mean distance of the copies' `shape_vec` from their centroid,
/// divided by 2.
///
/// The divisor is the chord bound — two antipodal points on the unit sphere are 2 apart —
/// so the result is normalised into `[0, 1]` by a *geometric* constant. It carries no
/// dependence on `sigma_E(0)` and therefore none on cell width, which is what keeps
/// `ensemble_spread` free of the resolution confound that afflicts `error_ratio`.
pub fn spread_shape<T: Real>(n: &[[T; 3]]) -> T {
    if n.len() < 2 {
        return T::zero();
    }
    let cnt = T::lit(n.len() as f64);
    let mut c = [T::zero(); 3];
    for v in n {
        for k in 0..3 {
            c[k] += v[k];
        }
    }
    for k in 0..3 {
        c[k] = c[k] / cnt;
    }
    let mut acc = T::zero();
    for v in n {
        let d = [v[0] - c[0], v[1] - c[1], v[2] - c[2]];
        acc += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    }
    acc / cnt / T::lit(2.0)
}
