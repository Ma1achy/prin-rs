//! Gragg-Bulirsch-Stoer extrapolation over the logH leapfrog.
//!
//! Mikkola & Merritt, AJ **135** (2008) 2398, run the algorithmic-regularisation leapfrog under
//! GBS rather than alone, and that is the configuration they recommend. Everything measured in
//! `results/output/logh_arms.txt` is the leapfrog *bare*, which is not how the method is meant to
//! be used — so this is the rematch, and the plan said it had to be a separate experiment rather
//! than something quietly folded into the RK4 comparison.
//!
//! # Why the leapfrog specifically, and why this is not free
//!
//! Extrapolation in `h^2` is valid only if the base method's error expansion contains **even
//! powers of `h` alone**. That holds for a time-symmetric method and fails for a general one. The
//! drift-kick-drift leapfrog is symmetric, and `n` equal DKD steps compose to a symmetric map, so
//! each extrapolation level buys **two** orders rather than one. RK4 is not symmetric; putting it
//! under the same extrapolation would gain one order per level and cost four evaluations a step
//! to do it.
//!
//! That symmetry is a claim about this implementation and not only about the literature, so
//! `tests/logh_gbs.rs` finite-differences it: the macro-step must reverse, and the observed order
//! must rise as `2k` with the number of levels. **If the order does not rise, the extrapolation
//! is not working and every number from it is a slower leapfrog.**
//!
//! # The sequence and the extrapolation
//!
//! Substep counts `n_k = 2, 4, 6, 8, ...` — even, so `h_k = H/n_k` and the expansion is in
//! `h_k^2`. Aitken-Neville to `h = 0`, in the form that needs only the ratio of substep counts:
//!
//! ```text
//!   T[k][j] = T[k][j-1] + (T[k][j-1] - T[k-1][j-1]) / ((n_k/n_{k-j})^2 - 1)
//! ```
//!
//! The whole 13-component state is extrapolated, **`t` included**. That is deliberate: `t` is a
//! dependent variable of the fictitious-time march like any other, and extrapolating the state
//! while leaving the clock at its unextrapolated value would return a state and a time that
//! belong to different trajectories.
//!
//! # Cost is counted, never derived
//!
//! A macro-step accepted at level `k` has spent `sum_{j<=k} n_j = k(k+1)` force evaluations —
//! 42 at `k = 6`. So `steps` is meaningless here in a way it is merely incomparable elsewhere,
//! and [`GbsOut::evals`] is returned rather than inferred. `MarchOut::force_evals` is what any
//! table must use.
//!
//! # And the advance-anyway site is recorded
//!
//! Standard GBS reduces the macro-step when the tolerance cannot be met at `k_max`. This one does
//! not — the macro-step is set by `eta` so the arms stay comparable — so it can reach `k_max`
//! still above tolerance and **advance anyway**. That is exactly the class of site
//! `results/saturation/README.md` is about: it is not terminal, the march continues, and if it is
//! not counted it is invisible. [`GbsOut::converged`] is false there and the driver totals it.

use crate::Real;

use super::hamiltonian::LhTime;
use super::state::LhState;
use super::step;
use super::system::LhSystem;

/// Substep counts. Even throughout, so every `h_k` has an `h^2` expansion to extrapolate in.
pub const SEQ: [usize; 12] = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24];

/// The largest level the sequence supports.
pub const K_MAX: usize = SEQ.len();

#[derive(Clone, Copy, Debug)]
pub struct GbsOut<T> {
    pub state: LhState<T>,
    /// Force evaluations actually spent, summed over every level attempted.
    pub evals: usize,
    /// Levels used, `1..=k_max`.
    pub k_used: usize,
    /// The last error estimate, `|T[k][k] - T[k][k-1]|` in the scaled norm.
    pub err: T,
    /// False if `k_max` was reached with the estimate still above tolerance, and the step was
    /// taken anyway. **An advance-anyway site**, counted by the driver.
    pub converged: bool,
}

/// Relative-plus-absolute scaled max norm over the thirteen components.
///
/// A plain max would be dominated by whichever coordinate happens to be largest, which on these
/// charts spans orders of magnitude; a plain relative one blows up wherever a component passes
/// through zero. The mixed form is the standard ODE-solver scaling and is what makes one
/// tolerance mean the same thing across `far` and `deep interior`.
#[inline]
fn scaled_diff<T: Real>(a: &[T; 13], b: &[T; 13], atol: T) -> T {
    let mut w = T::zero();
    for i in 0..13 {
        let scale = a[i].abs().max(b[i].abs()).max(atol);
        let d = (a[i] - b[i]).abs() / scale;
        if d > w {
            w = d;
        }
    }
    w
}

/// `n` leapfrog substeps of size `h/n`. Returns the state and the evaluations spent (`= n`).
///
/// Consecutive half-drifts are **not** merged. `drift` advances `r += v*(ds/(K+B))` with `K`
/// unchanged, so `D(a)` then `D(b)` is `D(a+b)` up to rounding, and the merged form would buy
/// only arithmetic. Keeping the plain composition means this shares `step::kdk` with the bare
/// leapfrog arm exactly, so a GBS-against-leapfrog difference is the extrapolation and not a
/// second implementation of the stepper.
#[inline]
fn substeps<T: Real>(
    sys: &LhSystem<T>,
    s: &LhState<T>,
    b: T,
    time: LhTime,
    h: T,
    n: usize,
) -> (LhState<T>, usize) {
    let hs = h / T::lit(n as f64);
    let mut cur = *s;
    let mut evals = 0usize;
    for _ in 0..n {
        let (next, e) = step::kdk(sys, &cur, b, time, hs);
        cur = next;
        evals += e;
    }
    (cur, evals)
}

/// One extrapolated macro-step of fictitious length `h`.
///
/// `tol` is the relative tolerance on the extrapolation error estimate; `k_max` caps the levels.
/// `atol` floors the scaling in [`scaled_diff`].
pub fn macro_step<T: Real>(
    sys: &LhSystem<T>,
    s: &LhState<T>,
    b: T,
    time: LhTime,
    h: T,
    tol: T,
    k_max: usize,
) -> GbsOut<T> {
    let kmax = k_max.clamp(1, K_MAX);
    let atol = T::lit(1e-30);
    // `tab[j]` holds row `k`'s entries as they are built; `prev[j]` is row `k-1`.
    let mut prev: Vec<[T; 13]> = Vec::with_capacity(kmax);
    let mut evals = 0usize;
    let mut best = *s;
    let mut err = T::infinity();

    for k in 0..kmax {
        let (raw, e) = substeps(sys, s, b, time, h, SEQ[k]);
        evals += e;
        let mut row: Vec<[T; 13]> = Vec::with_capacity(k + 1);
        row.push(raw.to_array13());

        for j in 1..=k {
            // `(n_k / n_{k-j})^2 - 1`, the Aitken-Neville denominator for an `h^2` expansion.
            let r = T::lit(SEQ[k] as f64) / T::lit(SEQ[k - j] as f64);
            let den = r * r - T::one();
            let (cur, up) = (row[j - 1], prev[j - 1]);
            let mut next = [T::zero(); 13];
            for i in 0..13 {
                next[i] = cur[i] + (cur[i] - up[i]) / den;
            }
            row.push(next);
        }

        let top = row[k];
        best = LhState::from_array13(top);
        if k > 0 {
            err = scaled_diff(&top, &row[k - 1], atol);
            if err <= tol && best.is_finite() {
                return GbsOut { state: best, evals, k_used: k + 1, err, converged: true };
            }
        }
        if !best.is_finite() {
            // Extrapolating a diverged row produces a diverged answer; stop rather than spend
            // the remaining levels on it. `NaN >= x` is false, so this must be tested and not
            // inferred from a comparison.
            return GbsOut { state: best, evals, k_used: k + 1, err, converged: false };
        }
        prev = row;
    }
    GbsOut { state: best, evals, k_used: kmax, err, converged: false }
}
