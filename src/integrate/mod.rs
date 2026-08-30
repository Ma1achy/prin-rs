//! Trajectory integration.

pub mod az;
pub mod heggie;
pub mod leapfrog;

use crate::physics::Cart;
use crate::Real;

/// What one trajectory reports back.
///
/// `finite` is carried explicitly rather than inferred from NaNs downstream. A copy that
/// diverged is a *measurement outcome* — "this could not be determined" — and BRIEF §3
/// forbids discarding it. Making the flag explicit keeps that information without letting a
/// NaN contaminate every `min`/`max` it touches (NOTES §2.3).
#[derive(Clone, Copy, Debug)]
pub struct TrajOut<T> {
    pub state: Cart<T>,
    /// Physical time actually reached.
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|`, the integration-quality measure.
    pub drift: T,
    /// Closest approach over all three pairs.
    pub d_min: T,
    pub steps: usize,
    /// False if the trajectory went non-finite, or ran out of step budget.
    pub finite: bool,
    /// True if the step budget was exhausted before `t_max`.
    pub budget_exhausted: bool,
}

impl<T: Real> TrajOut<T> {
    pub fn reached(&self, t_max: T) -> bool {
        self.finite && !self.budget_exhausted && self.t >= t_max
    }
}

/// Which integrator marches a trajectory.
///
/// Both regularise; they differ in **what they regularise around**. AZ picks a reference body —
/// the one not in the longest side — regularises the two pairs sharing it, and re-chooses at
/// every sync boundary. Heggie regularises all three relative vectors symmetrically and has no
/// reference body to choose, so its march runs uninterrupted from `t = 0`.
///
/// Measured on `config_stability` at 256^2, doubling the sync cadence at fixed step size moves
/// AZ's drift field by **0.516 decades** and Heggie's by **0.048** — the reason this enum exists.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Integrator {
    /// Aarseth-Zare. The default, and what every committed number in `results/` was taken under.
    #[default]
    Az,
    /// Heggie 1974 global regularisation.
    Heggie,
}

impl Integrator {
    pub fn name(self) -> &'static str {
        match self {
            Integrator::Az => "az",
            Integrator::Heggie => "heggie",
        }
    }
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "az" | "AZ" => Some(Integrator::Az),
            "heggie" | "hg" | "HG" => Some(Integrator::Heggie),
            _ => None,
        }
    }
}

/// What the ensemble layer needs from a march, whichever integrator produced it.
///
/// # Why a struct and not a trait
///
/// `pixel::evaluate_at` reads nineteen fields off `AzOut` and then runs two hundred lines of
/// statistics over them. A trait would need nineteen methods and would buy nothing: there is no
/// third implementation in prospect and no behaviour to dispatch, only data to carry.
///
/// # Absent is not zero
///
/// Six of the nineteen have **no Heggie analogue at all** — `refs`, `switches`, `ref_tie` are
/// about a reference body it does not have, and `n_retry`, `retry_exhausted`, `n_cap_hits` are
/// about step-control machinery it does not run. They are empty vectors, `0`, or **non-finite**,
/// never a plausible-looking zero: `evaluate_at` already filters non-finite out of its
/// reductions, so an absence is dropped rather than folded in as a value. A `0.0` for `ab_min`
/// would read as "the product hit its floor on every step", which is the opposite of the truth.
#[derive(Clone, Debug)]
pub struct MarchOut<T> {
    pub state: Cart<T>,
    pub drift: T,
    pub d_min_ref: T,
    pub d_min_true: T,
    pub gamma_max: T,
    pub steps: usize,
    pub finite: bool,
    pub budget_exhausted: bool,
    pub events: crate::outcome::Events<T>,
    pub t_end: T,
    pub tight: Vec<u8>,
    pub boundary_shapes: Vec<[T; 3]>,
    pub drift_hist: Vec<T>,
    pub tie_ratio: Vec<T>,
    pub dt_max: T,
    pub n_overshoot: u32,
    pub ab_min: T,
    pub ab_floored: bool,
    // --- AZ-only, and absent rather than zero under Heggie ---
    pub refs: Vec<u8>,
    pub switches: u32,
    pub ref_tie: Vec<T>,
    pub n_cap_hits: u32,
    pub n_retry: u32,
    pub retry_exhausted: bool,
}

impl<T: Real> From<az::AzOut<T>> for MarchOut<T> {
    fn from(o: az::AzOut<T>) -> Self {
        Self {
            state: o.state,
            drift: o.drift,
            d_min_ref: o.d_min_ref,
            d_min_true: o.d_min_true,
            gamma_max: o.gamma_max,
            steps: o.steps,
            finite: o.finite,
            budget_exhausted: o.budget_exhausted,
            events: o.events,
            t_end: o.t_end,
            tight: o.tight,
            boundary_shapes: o.boundary_shapes,
            drift_hist: o.drift_hist,
            tie_ratio: o.tie_ratio,
            dt_max: o.dt_max,
            n_overshoot: o.n_overshoot,
            ab_min: o.ab_min,
            ab_floored: o.ab_floored,
            refs: o.refs,
            switches: o.switches,
            ref_tie: o.ref_tie,
            n_cap_hits: o.n_cap_hits,
            n_retry: o.n_retry,
            retry_exhausted: o.retry_exhausted,
        }
    }
}

impl<T: Real> From<heggie::HgOut<T>> for MarchOut<T> {
    fn from(o: heggie::HgOut<T>) -> Self {
        Self {
            state: o.state,
            // **The Cartesian drift, not `drift_reg`.** AZ reports the energy of the state it
            // returns and so must this, or the two integrators are compared on two different
            // measurements — which they were, for the whole of the Phase 4 table, flattering
            // Heggie by up to 280x on a deep collision.
            drift: o.drift,
            // No unregularised side, so the two `d_min` measures coincide **by construction**.
            // In AZ their gap measures how well the reference-switching cadence tracks
            // encounters; here it is identically zero because there is no cadence to track with.
            d_min_ref: o.d_min,
            d_min_true: o.d_min,
            gamma_max: o.gamma_max,
            steps: o.steps,
            finite: o.finite,
            budget_exhausted: o.budget_exhausted,
            events: o.events,
            t_end: o.t_end,
            tight: o.tight,
            boundary_shapes: o.boundary_shapes,
            drift_hist: o.drift_hist,
            tie_ratio: o.tie_ratio,
            dt_max: o.dt_max,
            n_overshoot: o.n_overshoot,
            // `A*B` does not exist here; `R1 R2 R3` is the analogue and is not the same quantity,
            // so this is an absence. Non-finite, because `evaluate_at` filters those out of its
            // reduction and a zero would read as "floored on every step".
            ab_min: T::nan(),
            ab_floored: o.r_floored,
            refs: Vec::new(),
            switches: 0,
            ref_tie: Vec::new(),
            n_cap_hits: 0,
            n_retry: 0,
            retry_exhausted: false,
        }
    }
}
