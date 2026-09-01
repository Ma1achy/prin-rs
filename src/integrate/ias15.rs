//! **IAS15** — Rein & Spiegel, MNRAS **446** (2015) 1424. A 15th-order Gauss-Radau integrator with
//! adaptive step control, accurate to machine precision, whose error follows Brouwer's law: a
//! random walk rather than systematic drift.
//!
//! # What it is here for, and what it is NOT here for
//!
//! **A reference arm.** This project does not have one. `eta/256` was tried and came back
//! *saturated* — chord 2.000, antipodal, for a correct mode and a broken one alike — so it could
//! not separate arms it was asked to rank. IAS15 reaches machine precision at a cost that makes it
//! usable as ground truth on a small grid.
//!
//! **It is not a production candidate, and that verdict is already measured rather than assumed.**
//! Its predictor-corrector iterates a variable number of times per step, which is per-lane variable
//! work — the exact property that killed reject-and-retry here, where every warp contained a
//! retrying lane (`warps hit 1.0000`) and the worst lane retried 5.2 million times. Do not
//! re-derive that; it is on record.
//!
//! # The conversion matrix is COMPUTED, not transcribed
//!
//! The `g -> b` conversion is 21 constants in the reference implementations. Two sign errors in
//! this project's AZ algebra were invisible until someone finite-differenced the Hamiltonian —
//! wrong algebra of this kind produces trajectories that look like physics. So the matrix is built
//! here by expanding the Newton basis into the monomial basis:
//!
//! ```text
//!   sum_k b_k t^k  ==  sum_k g_k prod_{j<k} (t - h_j)
//! ```
//!
//! which is a short polynomial multiplication with no magic numbers in it, and
//! `tests/ias15.rs::the_conversion_matrix_satisfies_its_defining_identity` checks it against that
//! identity directly at random points. A transcribed table would need a test to catch a typo; a
//! computed one needs a test to catch a wrong *definition*, which is the more useful thing to be
//! forced to write down.
//!
//! # Substep formulas
//!
//! With `h` a Gauss-Radau node in `[0, 1]` and `dt` the step:
//!
//! ```text
//!   x(h) = x0 + h dt v0 + (h dt)^2/2 [ a0 + sum_k b_k h^(k+1) * 2/((k+2)(k+3)) ]
//!   v(h) = v0 + h dt    [ a0 + sum_k b_k h^(k+1) / (k+2) ]
//! ```

use crate::physics::{newton, Cart};
use crate::{Real, Vec2};

/// Gauss-Radau spacings, `h_0 = 0` and seven interior nodes. These ARE transcribed — they are the
/// roots of a Legendre-family polynomial and there is no short expansion to compute them from — so
/// `tests/ias15.rs` checks the order they deliver rather than the digits themselves. A wrong node
/// costs order, which is measurable; a wrong digit in the fifteenth place is not.
pub const H: [f64; 8] = [
    0.0,
    0.056_262_560_536_922_146,
    0.180_240_691_736_892_36,
    0.352_624_717_113_169_64,
    0.547_153_626_330_555_4,
    0.734_210_177_215_410_5,
    0.885_320_946_839_095_8,
    0.977_520_613_561_287_5,
];

/// Stages, i.e. the number of `b` coefficients. Seven gives fifteenth order.
pub const N: usize = 7;

/// The `g -> b` matrix, computed from [`H`] by expanding the Newton basis.
///
/// `c[k][j]` is the coefficient of `t^(j+1)` in `prod_{i=0..=k} (t - h_i)`, so
/// `b_j = sum_k c[k][j] g_k`.
///
/// # The `j + 1` is not cosmetic and the first cut got it wrong
///
/// `g_k` multiplies `prod_{i=0..=k}(t - h_i)`, and `h_0 = 0`, so **every basis polynomial carries
/// a factor of `t`**: the interpolant of `a(h) - a0` has no constant term, and IAS15's `b_k` is
/// the coefficient of `h^(k+1)` rather than of `h^k`. Recording the product *before* multiplying
/// by `(t - h_k)`, and indexing from `t^0`, shifts the whole matrix one place. It compiles, it
/// runs, and it integrates: measured energy drift **7.8e-1** at `t = 200` where the method should
/// give `1e-15`, and a halving that bought 3.6x rather than the ~32000x of fifteenth order.
///
/// `tests/ias15.rs::the_conversion_matrix_satisfies_its_defining_identity` checks this against the
/// product form directly, with a perturbed-matrix negative control.
pub fn conversion() -> [[f64; N]; N] {
    let mut c = [[0.0f64; N]; N];
    // poly holds the coefficients of prod_{i=0..k}(t - h_i), lowest order first.
    let mut poly = [0.0f64; N + 2];
    poly[0] = 1.0;
    let mut deg = 0usize;
    for k in 0..N {
        // Multiply FIRST: g_k's basis includes the (t - h_k) factor.
        let hk = H[k];
        let mut next = [0.0f64; N + 2];
        for j in 0..=deg {
            next[j + 1] += poly[j];
            next[j] -= hk * poly[j];
        }
        poly = next;
        deg += 1;
        // Then record, from t^1 upward -- there is no constant term to record.
        for j in 0..N {
            c[k][j] = poly[j + 1];
        }
    }
    c
}

/// One IAS15 step of size `dt` from `(r, v)`.
///
/// Returns the new state, the **force evaluations** spent, the predictor-corrector iterations
/// used, and `max |b_6|` — the highest coefficient, which is what [`next_dt`] reads.
///
/// The counts are returned rather than derived: the iteration count is variable, which is the
/// whole reason this integrator is a reference arm and not a candidate. **`b_6` is returned rather
/// than recomputed** because a step controller handed a constant in its place silently degrades to
/// "scale `dt` by a fixed factor every step" — measured, that turned a machine-precision method
/// into `1.282e-7` drift and looked like an integrator fault.
#[allow(clippy::needless_range_loop)]
pub fn step<T: Real>(
    r0: &[Vec2<T>; 3],
    v0: &[Vec2<T>; 3],
    m: &[T; 3],
    dt: T,
    max_iter: usize,
    tol: T,
) -> (Cart<T>, usize, usize, T) {
    let cmat = conversion();
    let a0 = newton::accel(r0, m, T::zero());
    let mut evals = 1usize;

    // b and g coefficients, [stage][body][component-as-Vec2].
    let mut b = [[Vec2::<T>::zero(); 3]; N];
    let mut g = [[Vec2::<T>::zero(); 3]; N];

    let mut iters = 0usize;
    for _ in 0..max_iter {
        iters += 1;
        let b_prev = b;
        // Accelerations at the seven interior nodes.
        let mut at = [[Vec2::<T>::zero(); 3]; N];
        for n in 0..N {
            let h = T::lit(H[n + 1]);
            let mut rn = [Vec2::<T>::zero(); 3];
            for i in 0..3 {
                // x(h) = x0 + h dt v0 + (h dt)^2/2 [ a0 + sum_k b_k h^(k+1) 2/((k+2)(k+3)) ]
                let mut acc = a0[i];
                let mut hp = h;
                for k in 0..N {
                    let w = T::lit(2.0 / (((k + 2) * (k + 3)) as f64));
                    acc = acc + b[k][i] * (hp * w);
                    hp = hp * h;
                }
                let hdt = h * dt;
                rn[i] = r0[i] + v0[i] * hdt + acc * (hdt * hdt * T::lit(0.5));
            }
            at[n] = newton::accel(&rn, m, T::zero());
            evals += 1;

            // **Newton divided differences**, which is what `g` is by definition:
            //
            //   g_n = ( ... ((a_{n+1} - a_0)/(h_{n+1} - h_0) - g_0)/(h_{n+1} - h_1) ... - g_{n-1})
            //         / (h_{n+1} - h_n)
            //
            // Written as the recurrence rather than as a table, for the same reason `conversion`
            // is computed: a transcribed coefficient that is wrong produces a trajectory that
            // still looks like physics.
            for i in 0..3 {
                let mut acc = (at[n][i] - a0[i]) / T::lit(H[n + 1] - H[0]);
                for k in 0..n {
                    acc = (acc - g[k][i]) / T::lit(H[n + 1] - H[k + 1]);
                }
                g[n][i] = acc;
            }
        }
        // b from g.
        for j in 0..N {
            for i in 0..3 {
                let mut s = Vec2::<T>::zero();
                for k in 0..N {
                    s = s + g[k][i] * T::lit(cmat[k][j]);
                }
                b[j][i] = s;
            }
        }
        // Convergence on the highest coefficient, which is what the step control also reads.
        let mut delta = T::zero();
        let mut scale = T::zero();
        for i in 0..3 {
            delta = delta.max((b[N - 1][i] - b_prev[N - 1][i]).norm());
            scale = scale.max(a0[i].norm());
        }
        if scale > T::zero() && delta / scale < tol {
            break;
        }
    }

    // Final update at h = 1.
    let mut r = [Vec2::<T>::zero(); 3];
    let mut v = [Vec2::<T>::zero(); 3];
    for i in 0..3 {
        let mut ax = a0[i];
        let mut av = a0[i];
        for k in 0..N {
            ax = ax + b[k][i] * T::lit(2.0 / (((k + 2) * (k + 3)) as f64));
            av = av + b[k][i] * T::lit(1.0 / ((k + 2) as f64));
        }
        r[i] = r0[i] + v0[i] * dt + ax * (dt * dt * T::lit(0.5));
        v[i] = v0[i] + av * dt;
    }
    let b_last = (0..3).fold(T::zero(), |w, i| w.max(b[N - 1][i].norm()));
    (Cart { r, v }, evals, iters, b_last)
}

/// Rein & Spiegel's step-size estimator: the ratio of the highest coefficient to the acceleration
/// scale, taken to the `1/7` power because `b_6` enters at order `dt^7` relative to `a0`.
pub fn next_dt<T: Real>(r: &[Vec2<T>; 3], m: &[T; 3], dt: T, b_last: T, eps: T) -> T {
    let a = newton::accel(r, m, T::zero());
    let scale = a.iter().fold(T::zero(), |w, x| w.max(x.norm()));
    if !(scale > T::zero()) || !b_last.is_finite() || !(b_last > T::zero()) {
        return dt;
    }
    let ratio = b_last / scale;
    let f = (eps / ratio).to_f64().unwrap().powf(1.0 / 7.0);
    // Clamped: an unbounded growth factor is how a step control produces the 2.209e128 step this
    // project has already recorded once.
    dt * T::lit(f.clamp(0.2, 4.0))
}
