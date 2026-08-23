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

// ---------------------------------------------------------------------------------------
// The inverse: shape sphere -> configuration.
//
// `shape_vec` above is the forward Hopf map from mass-weighted Jacobi coordinates. This is
// its closed-form inverse, and it is an **inverse, not an invention** — the round-trip
// `shape_vec(from_shape(n, I, phi)) == n` is a test that can fail, and does if any sign is
// wrong.
//
// It exists because the vertical-slice build needs a chart whose decode is genuinely
// NONLINEAR. `Slice::decode_pos` is a linspace and `Slice::nominal` writes (x, y) into
// `r[body]`, so that decode is affine: `J_D` is constant and the linearised path
// `x = x0 + J_D.delta` is exact rather than approximate. Measuring "where does the
// linearisation start to matter" on an affine chart returns "never", at every depth, in
// exact arithmetic — a measurement that cannot fail. This chart gives the curvature term
// something to be.
// ---------------------------------------------------------------------------------------

/// The two Jacobi reduced masses, `(mu_rho, mu_lambda)`, for the hard-wired (0,1)+2 pairing.
pub fn reduced_masses(m: &[f64; 3]) -> (f64, f64) {
    let (m0, m1, m2) = (m[0], m[1], m[2]);
    let mtot = m0 + m1 + m2;
    (m0 * m1 / (m0 + m1), m2 * (m0 + m1) / mtot)
}

/// Mass-weighted moment of inertia `I = |r~|^2 + |l~|^2` of a configuration. The scale that
/// `shape_vec` quotients out, recovered so a shape chart can put it back.
pub fn inertia(r: &[Vec2<f64>; 3], m: &[f64; 3]) -> f64 {
    let (m0, m1, m2) = (m[0], m[1], m[2]);
    let (mu_rho, mu_lam) = reduced_masses(m);
    let _ = m2;
    let rho = r[1] - r[0];
    let com01 = (r[0] * m0 + r[1] * m1) / (m0 + m1);
    let lam = r[2] - com01;
    rho.norm_sq() * mu_rho + lam.norm_sq() * mu_lam
}

/// Configuration from a point on the shape sphere, a scale, and the fibre phase.
///
/// ```text
/// a = I(1+n0)/2   b = I(1-n0)/2   p = I*n1/2   q = I*n2/2      (p^2 + q^2 = ab identically)
/// r~ = sqrt(a) (cos phi, sin phi)
/// l~ = sqrt(b) (cos psi, sin psi),   psi = phi + atan2(-q, p)
/// ```
///
/// **The sign in `atan2(-q, p)` is load-bearing.** `shape_vec` computes
/// `q = rt.y*lt.x - rt.x*lt.y`, which is the *negative* of the standard 2D cross product
/// `rt x lt`. Writing `atan2(q, p)` reflects the configuration and the round-trip test
/// fires on it.
///
/// The centre of mass is placed at the origin, and the bodies are released from rest —
/// every configuration in this project is.
pub fn from_shape(n: [f64; 3], inertia: f64, phase: f64, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let (m0, m1, m2) = (m[0], m[1], m[2]);
    let mtot = m0 + m1 + m2;
    let (mu_rho, mu_lam) = reduced_masses(m);

    let a = inertia * (1.0 + n[0]) / 2.0;
    let b = inertia * (1.0 - n[0]) / 2.0;
    let p = inertia * n[1] / 2.0;
    let q = inertia * n[2] / 2.0;

    let psi = phase + (-q).atan2(p);
    let rt = Vec2::new(a.max(0.0).sqrt() * phase.cos(), a.max(0.0).sqrt() * phase.sin());
    let lt = Vec2::new(b.max(0.0).sqrt() * psi.cos(), b.max(0.0).sqrt() * psi.sin());

    let rho = rt / mu_rho.sqrt();
    let lam = lt / mu_lam.sqrt();

    // COM at the origin: c01 = -m2*lam/M, then r2 = c01 + lam.
    let c01 = lam * (-m2 / mtot);
    [
        c01 - rho * (m1 / (m0 + m1)),
        c01 + rho * (m0 / (m0 + m1)),
        c01 + lam,
    ]
}

/// Exponential map on the unit 2-sphere: step from `n0` along the tangent vector `t`.
///
/// This is where a shape chart's **curvature** lives. `n0 cos|t| + (t/|t|) sin|t|` is not
/// affine in `t`, so `J_D` genuinely varies across a quad and the linearised decode path
/// genuinely approximates.
pub fn exp_map(n0: [f64; 3], t: [f64; 3]) -> [f64; 3] {
    let th = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
    if th == 0.0 {
        return n0;
    }
    let (c, s) = (th.cos(), th.sin() / th);
    let mut n = [
        n0[0] * c + t[0] * s,
        n0[1] * c + t[1] * s,
        n0[2] * c + t[2] * s,
    ];
    let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    for v in n.iter_mut() {
        *v /= norm;
    }
    n
}

/// An orthonormal tangent frame at `n0`, chosen deterministically so a chart is reproducible.
pub fn tangent_frame(n0: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Cross with whichever axis is least aligned with n0, so the frame never degenerates.
    let k = (0..3).min_by(|&i, &j| n0[i].abs().partial_cmp(&n0[j].abs()).unwrap()).unwrap();
    let mut ax = [0.0; 3];
    ax[k] = 1.0;
    let mut e1 = [
        n0[1] * ax[2] - n0[2] * ax[1],
        n0[2] * ax[0] - n0[0] * ax[2],
        n0[0] * ax[1] - n0[1] * ax[0],
    ];
    let n1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
    for v in e1.iter_mut() {
        *v /= n1;
    }
    let e2 = [
        n0[1] * e1[2] - n0[2] * e1[1],
        n0[2] * e1[0] - n0[0] * e1[2],
        n0[0] * e1[1] - n0[1] * e1[0],
    ];
    (e1, e2)
}
