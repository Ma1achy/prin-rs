//! Trajectory integration.

pub mod az;
pub mod heggie;
pub mod ias15;
pub mod leapfrog;
pub mod logh;

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

/// Which regularisation and which stepper march a trajectory.
///
/// Flat rather than a product of `(regularisation, stepper)`, because **the product has invalid
/// cells and saying so is part of the point**: there is no `Az + Leapfrog` and no
/// `Heggie + Leapfrog`. Their `Gamma` couples position and momentum, so neither Hamiltonian is
/// separable and a drift-kick-drift composition does not apply. There is no common stepper that
/// leaves all three methods intact, and matching arms on *steps* across steppers would be a
/// confound rather than a fairness measure — see [`Profile::evals_per_step`].
///
/// # The registry
///
/// ```text
///   variant          chart          re-registrations   owns dt/ds   stepper
///   Az               2 KS pairs     every boundary     yes          RK4
///   Heggie           3 KS vectors   none               yes          RK4
///   LogHLeapfrog     none           none               yes          KDK
///   LogHRk4          none           none               yes          RK4
///   LogHGbs          none           none               yes          GBS over KDK
///   PlainLeapfrog    none           none               NO           KDK
///   PlainRk4         none           none               NO           RK4
///   PlainGbs         none           none               NO           GBS over KDK
/// ```
///
/// `principia_integrator_contract.md` is the **GLSL app's** contract and is not in this repo;
/// its `substep_bucket`/`N_sub`/`N_max`/descriptor bit 5 appear nowhere in `src/`. So the
/// registry it describes is built here in this codebase's own terms, as [`Profile`], and
/// asserted in `tests/logh_seam.rs` rather than borrowed from a document this port does not
/// implement.
///
/// # Why the two logH arms and the two controls exist
///
/// Measured on `config_stability` at 256^2 with the step size held fixed by scaling `eta` with
/// `n_sync`, doubling the sync cadence moves AZ's drift field by **0.516 decades** and Heggie's
/// by **0.048** — the re-registration mechanism, and the reason this enum existed at all.
///
/// logH is the **falsification test** for it: no coordinate transformation of any kind, which is
/// a strictly stronger form of the property Heggie's win is attributed to. Two steppers, because
/// Mikkola & Merritt are explicit that in these methods *"the regularization is achieved by using
/// the leapfrog"* — already confirmed here on the radial collision, which KDK traverses and RK4
/// does not. And two `Plain*` controls, so the stepper's own contribution is measured rather than
/// assumed: they are the **same code path** with the time transformation switched off, which is a
/// tighter control than a separate integrator could be.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Integrator {
    /// Aarseth-Zare. **No longer the default**, and what every number committed to `results/`
    /// before 2026-09-02 was taken under. Kept, and not only for history: it is the **only
    /// occupant with an independent implementation to check against** — `xcheck` runs it against
    /// `reference/tb_az.py` and reads `0.000e+00` over all 32 reference points. `reference_opts`
    /// pins it, so that gate is unaffected by this default. AZ is also genuinely better on
    /// **sustained hierarchy**, where it never re-registers: `far`, on all 65536 pixels, by a
    /// flat 0.7–0.9 decades.
    Az,
    /// Heggie 1974 global regularisation. **The default from 2026-09-02.**
    ///
    /// # Why, and the reason is not the scoreboard
    ///
    /// The win is large — 31 of 32 gallery cases, `err>10` 3915 against 74, and AZ's worst drift
    /// decile fixed on **100% of its pixels by a median 1900x** — and it survived a full corpus
    /// regeneration under two later fixes with `drift p50` ratios at 0.987–1.004. But a default
    /// justified by a win count does not survive the next investigation. The mechanism is:
    ///
    /// **Heggie has no reference body, so the state is never re-registered into a chart mid-run,
    /// and re-registration is what the drift field's boundary-cadence sensitivity costs.**
    ///
    /// **It is NOT what removed the visible wedges, and an earlier version of this comment said
    /// it was.** Measured directly (`examples/wedge_census.rs`, `results/wedge/`): on
    /// `config_stability` at `t_max = 50`, AZ with the step-control fixes *off* carries a
    /// dense-pale set of 0.0026 of the frame -- the record's own 0.0018 at 1024^2 -- and turning
    /// them on takes it to **0.0001** with non-finite `1153 -> 0`. `heggie` then reads
    /// **identical to post-fix AZ to three digits**, as do `logh_rk4` and `plain_rk4`. Post-fix
    /// the slice is integrator-insensitive at the median: every pixel moves between arms, but
    /// `chord p50` is `1.2e-6` against a `chord max` of 1.907. The wedges were a numerical
    /// artefact and the `dtau` fix, the landing clamp and the predictive step limit are what
    /// removed them. Measured
    /// unconfounded at 256², termination off, step control held constant, every arm `nonfin <= 1`:
    ///
    /// ```text
    ///     AZ  n=250 CONTROLLED   lift 2.061   chord 5.162e-1
    ///     HG  n=250 CONTROLLED   lift 3.824   chord 4.832e-2   <- 10.7x smaller
    ///     HG  n=250 confounded   lift 2.327   chord 4.185e-1   <- the harness CAN see a large chord
    ///     HG  refresh_h n=125    lift 2.909   chord 5.570e-1   <- and sees BOUNDARY change
    /// ```
    ///
    /// Three guards, and the last is what makes the null mean anything: `steps p50` flat to 1.3%
    /// says `eta` held the step size, the confounded arm says the harness is not blind, and
    /// `refresh_h` — the one arm reintroducing a boundary-dependent quantity — says it sees
    /// **boundary** sensitivity specifically.
    ///
    /// # What does NOT support it, and must not be cited as if it did
    ///
    /// **`LogHRk4` and `LogHLeapfrog` losing does not establish the mechanism**, though they were
    /// built as its falsification test. They lose decisively — 29x to 1494x (RK4) and 16,000x to
    /// 215,000x (KDK) against Heggie at matched force evaluations, better on 0.0–3.6% of pixels,
    /// and third against the IAS15 reference. But logH differs from Heggie in **two** ways, not
    /// one: no re-registration *and* no coordinate transformation. The second is fatal here —
    /// logH reaches `d_min < 1e-10` and **never** `|dE/E| < 1e-12`, its drift *rising* with
    /// penetration depth where Heggie's `drift_reg` is flat at `4.4e-15`, because a KS map removes
    /// the `1/r` from the Hamiltonian and a time transformation only slows the clock. Every region
    /// tested is collision-dominated (`far` collides on 1048576 of 1048576 pixels), so logH was
    /// graded almost entirely on the one thing it is known not to do.
    /// **A chartless method is not a coordinate regularisation with the coordinates left out.**
    #[default]
    Heggie,
    /// logH under drift-kick-drift. **The method as designed**: one force evaluation per step.
    LogHLeapfrog,
    /// logH under RK4. Not how the method is meant to be used, and the arm directly comparable
    /// to AZ and Heggie. Four force evaluations per step.
    LogHRk4,
    /// Control: no regularisation, KDK. `dt/ds = 1`, the same code path as `LogHLeapfrog`.
    PlainLeapfrog,
    /// Control: no regularisation, RK4. `dt/ds = 1`, the same code path as `LogHRk4`.
    PlainRk4,
    /// logH under Gragg-Bulirsch-Stoer extrapolation over the leapfrog. **The configuration
    /// Mikkola & Merritt actually recommend**, and the one every earlier logH number was not.
    LogHGbs,
    /// Control: no regularisation, GBS. Says how much of the GBS arm's result is the
    /// extrapolation rather than the time transformation.
    PlainGbs,
}

/// What an occupant of the [`Integrator`] seam is, in fields a table can print.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Profile {
    /// How many coordinate charts the state passes through. `0` is algorithmic regularisation.
    pub charts: usize,
    /// Whether the state is re-registered into a chart during the march, and how often. `None`
    /// is "never"; `Some(n)` is "at every sync boundary" with `n` charts rebuilt.
    ///
    /// **This is the quantity the whole logH experiment is about.**
    pub re_registers: bool,
    /// Whether the integrator supplies its own time transformation, as opposed to marching in
    /// physical time. The `Plain*` controls are the only occupants that do not.
    pub owns_time_mapping: bool,
    /// Force evaluations per step. `steps` is only comparable between occupants sharing this
    /// number; `MarchOut::force_evals` is what a cross-stepper table must use.
    pub evals_per_step: usize,
}

impl Integrator {
    pub fn name(self) -> &'static str {
        match self {
            Integrator::Az => "az",
            Integrator::Heggie => "heggie",
            Integrator::LogHLeapfrog => "logh_lf",
            Integrator::LogHRk4 => "logh_rk4",
            Integrator::PlainLeapfrog => "plain_lf",
            Integrator::PlainRk4 => "plain_rk4",
            Integrator::LogHGbs => "logh_gbs",
            Integrator::PlainGbs => "plain_gbs",
        }
    }
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "az" | "AZ" => Some(Integrator::Az),
            "heggie" | "hg" | "HG" => Some(Integrator::Heggie),
            "logh_lf" | "logh-lf" | "loghlf" => Some(Integrator::LogHLeapfrog),
            "logh_rk4" | "logh-rk4" | "loghrk4" => Some(Integrator::LogHRk4),
            "plain_lf" | "plain-lf" | "none_lf" => Some(Integrator::PlainLeapfrog),
            "plain_rk4" | "plain-rk4" | "none_rk4" => Some(Integrator::PlainRk4),
            "logh_gbs" | "logh-gbs" | "gbs" => Some(Integrator::LogHGbs),
            "plain_gbs" | "plain-gbs" | "none_gbs" => Some(Integrator::PlainGbs),
            _ => None,
        }
    }
    pub fn profile(self) -> Profile {
        use Integrator::*;
        match self {
            Az => Profile {
                charts: 2,
                re_registers: true,
                owns_time_mapping: true,
                evals_per_step: az::rk4::EVALS_PER_STEP,
            },
            Heggie => Profile {
                charts: 3,
                re_registers: false,
                owns_time_mapping: true,
                evals_per_step: heggie::rk4::EVALS_PER_STEP,
            },
            LogHLeapfrog => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: true,
                evals_per_step: logh::Stepper::Kdk.evals_per_step(),
            },
            LogHRk4 => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: true,
                evals_per_step: logh::Stepper::Rk4.evals_per_step(),
            },
            PlainLeapfrog => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: false,
                evals_per_step: logh::Stepper::Kdk.evals_per_step(),
            },
            PlainRk4 => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: false,
                evals_per_step: logh::Stepper::Rk4.evals_per_step(),
            },
            // `evals_per_step` is **0** for both GBS occupants, and that is the honest value
            // rather than a missing one: a macro-step accepted at level `k` costs `k(k+1)`
            // evaluations with `k` chosen adaptively, so any `steps * evals_per_step` a caller
            // derives comes out zero and obviously wrong instead of plausibly wrong.
            LogHGbs => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: true,
                evals_per_step: logh::Stepper::Gbs.evals_per_step(),
            },
            PlainGbs => Profile {
                charts: 0,
                re_registers: false,
                owns_time_mapping: false,
                evals_per_step: logh::Stepper::Gbs.evals_per_step(),
            },
        }
    }

    /// The `(time transformation, stepper)` pair for the four occupants that run through
    /// `logh::integrate_lh`, and `None` for AZ and Heggie.
    ///
    /// Exhaustive with no `_` arm, so a seventh variant breaks the build here rather than
    /// silently falling through to a default.
    pub fn logh_arms(self) -> Option<(logh::LhTime, logh::Stepper)> {
        use logh::{LhTime, Stepper};
        match self {
            Integrator::Az | Integrator::Heggie => None,
            Integrator::LogHLeapfrog => Some((LhTime::LogH, Stepper::Kdk)),
            Integrator::LogHRk4 => Some((LhTime::LogH, Stepper::Rk4)),
            Integrator::PlainLeapfrog => Some((LhTime::None, Stepper::Kdk)),
            Integrator::PlainRk4 => Some((LhTime::None, Stepper::Rk4)),
            Integrator::LogHGbs => Some((LhTime::LogH, Stepper::Gbs)),
            Integrator::PlainGbs => Some((LhTime::None, Stepper::Gbs)),
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
    /// Force evaluations. **`steps` is not comparable across steppers** — RK4 spends four per
    /// step and a drift-kick-drift leapfrog one — so any table with more than one stepper in it
    /// has to match and report on this instead.
    pub force_evals: usize,
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
    /// **An advance-anyway site, and its two zeros mean different things.** The GBS macro-step
    /// hit `gbs_k_max` without meeting `gbs_tol`; the extrapolated state is taken regardless.
    ///
    /// `0` under `Stepper::Gbs` says every macro-step converged. `0` under AZ, Heggie, RK4 or
    /// KDK says there is no extrapolation to converge — the same number for the absence and for
    /// the clean result. A harness must gate this column on the arm, and the denominator is
    /// `steps`, not pixels: this counts macro-steps, summed over copies.
    pub gbs_unconverged: u32,
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
            // Exact rather than an estimate: every `rk4::step` in the march is paired with
            // `steps += 1`, retries included, and the driver calls `deriv` nowhere else.
            force_evals: o.steps * az::rk4::EVALS_PER_STEP,
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
            gbs_unconverged: 0,
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
            force_evals: o.steps * heggie::rk4::EVALS_PER_STEP,
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
            // No extrapolation here, so this is an absence rather than a clean result.
            gbs_unconverged: 0,
        }
    }
}

impl<T: Real> From<logh::LhOut<T>> for MarchOut<T> {
    fn from(o: logh::LhOut<T>) -> Self {
        Self {
            state: o.state,
            drift: o.drift,
            // No chart, so there is no reference-relative separation to be distinct from the true
            // one. As in Heggie the two coincide, and here they coincide because there was never
            // a second definition rather than because two definitions agree.
            d_min_ref: o.d_min,
            d_min_true: o.d_min,
            // `rho = |K + B - U|/U`. Occupying `gamma_max` because it is the same *kind* of
            // quantity the other two put there — the regularised Hamiltonian's running residual.
            // It is the energy defect normalised by `U`, not an independent constraint; logH has
            // no analogue of Heggie's `sum q_i = 0` because its phase space is the physical one.
            gamma_max: o.gamma_max,
            steps: o.steps,
            // **Counted, not derived.** Its two steppers spend different numbers per step, so
            // there is no single multiplier to derive one from the other.
            force_evals: o.force_evals,
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
            // No `A*B` and no `R1 R2 R3`: **non-finite, never zero**, so `evaluate_at` drops it
            // from its reductions instead of folding in a value that would read as "the product
            // hit its floor on every step".
            ab_min: T::nan(),
            // The denominator degeneracy is the analogous *advance-anyway* site: `K + B` is `U`
            // on shell, so a non-positive value means the transformation itself has failed.
            ab_floored: o.den_degenerate,
            refs: Vec::new(),
            switches: 0,
            ref_tie: Vec::new(),
            n_cap_hits: 0,
            n_retry: 0,
            retry_exhausted: false,
            gbs_unconverged: o.gbs_unconverged,
        }
    }
}
