//! Trajectory integration.

pub mod az;
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
