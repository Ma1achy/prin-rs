//! The sync loop: choose a reference body, register, RK4-march to the boundary, return to
//! Cartesian, repeat.
//!
//! Transcribed from `reference/tb_az.py:integrate_az`.

use crate::outcome::{collision_pairs_from, Events};
use crate::physics::{energy, Cart};
use crate::Real;

use super::hamiltonian::gamma_residual;
use super::reference_body::{choose_reference, triple, RefPolicy};
use super::rk4;
use super::system::AzSystem;

/// How the fictitious step `dtau` is sized within a sync interval.
///
/// `dt = A*B*dtau`, so the physical step is `eta*dt_left` **only while `A*B` stays near its
/// entry value**. A trajectory sitting at a close encounter *at a sync boundary* has a tiny
/// `A0*B0`, so `dtau` is enormous; as the bodies separate through the interval `A*B` grows by
/// orders and `dt` grows with it. Giant physical steps immediately after an encounter, on a thin
/// set -- encounters coinciding with a boundary -- which is why the damage clusters spatially
/// rather than tracking `d_min`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DtauMode {
    /// `dtau = eta*dt_left/(A0*B0)`, computed **once** at the interval's entry and never
    /// updated. The reference's behaviour (`reference/tb_az.py:integrate_az`) and the shipped
    /// behaviour until this change, so **every committed Rust and NumPy number was taken with
    /// it on**. Kept, named for what it is, and never deleted.
    FixedPerInterval,
    /// `dtau = eta*(dt_left - s.t)/(A*B)` recomputed per step -- the obvious repair, and it is
    /// **Zeno by arithmetic**: `dt ~ eta*rem` gives `rem_{n+1} = rem_n (1 - eta)`, so the
    /// interval is approached geometrically and the loop can never satisfy `s.t >= dt_left`.
    /// At `eta = 0.01` it runs to `max_steps` on *every* interval. Kept as a measurement axis
    /// precisely because it looks right; it is not a candidate default.
    PerStepRemaining,
    /// `dtau = min(eta*dt_left/(A*B), dtau_entry)` -- `dt_left` held fixed, only `A*B`
    /// recomputed. `dt ~ eta*dt_left` throughout, so the step count stays at `~1/eta`.
    ///
    /// **The cap is one-sided in the right direction.** When `A*B` grows the recomputed value
    /// falls and the blow-up is removed; when `A*B` *falls* at a close approach the cap holds
    /// `dtau` at nominal, so `dt = A*B*dtau` shrinks with the separation -- which is what
    /// regularisation buys and what the original comment wanted. It removes the over-correction
    /// without reintroducing the thing that was over-corrected *for* (shrinking `dtau` with
    /// separation drove `dt -> 1e-13` and produced a false "intractable region").
    PerStepInterval,
}

impl Default for DtauMode {
    fn default() -> Self {
        DtauMode::PerStepInterval
    }
}

#[derive(Clone, Debug)]
pub struct AzOut<T> {
    pub state: Cart<T>,
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|`, Cartesian and unsoftened, exactly as the reference.
    pub drift: T,
    /// `min(|R1|, |R2|)` — the two **regularised** pairs only. This is what the reference
    /// tracks, and it is what the cross-check compares.
    pub d_min_ref: T,
    /// `min` over all three pairs, including the unregularised side. BRIEF §4 defines
    /// `d_min` this way; the reference does not compute it. The gap between the two
    /// measures how well the reference-switching cadence tracks encounters (NOTES §2.1).
    pub d_min_true: T,
    pub switches: u32,
    /// `max |H - E| / |E|` over the run — the free residual from `Gamma == 0` along the
    /// trajectory.
    pub gamma_max: T,
    /// Reference body chosen at each sync boundary. Compared as a column in the
    /// cross-check: a near-tie in `argmax` broken differently by the two implementations
    /// fails the comparison while looking exactly like a transcription error.
    pub refs: Vec<u8>,
    pub steps: usize,
    pub finite: bool,
    pub budget_exhausted: bool,
    /// Termination events, per BRIEF §2.4. Collisions are sampled **inside** the RK4 loop;
    /// escape at sync boundaries, which is where the reference samples it.
    pub events: Events<T>,
    /// Identity of the **currently tightest pair** at each completed sync boundary, as an
    /// index into [`crate::physics::PAIRS`].
    ///
    /// This is what `spread_event` is built on, and it is deliberately *not* the terminal
    /// outcome. A terminal label is terminal-grain: early in the march nothing has
    /// terminated, so every copy agrees, so a spread built on it reports maximum confidence
    /// at exactly the playhead where least is known. The tightest-pair identity is defined at
    /// every playhead, needs no gate, and moves earlier. See NOTES §2.8.
    pub tight: Vec<u8>,
    /// Second-tightest over tightest pair separation, at each completed sync boundary.
    ///
    /// A value near 1 is a **near-tie**: two pairs are nearly equally close, so copies can
    /// disagree about which is *tightest* without their trajectories having diverged at all.
    /// That distinction does not exist for a continuous divergence measure and it decides
    /// whether `spread_event` may latch — a running max would make a near-tie permanent and
    /// it would never clear. See NOTES §5.
    pub tie_ratio: Vec<T>,
    /// The **shape vector at each completed sync boundary**, when
    /// `AzOpts::keep_boundary_shapes` is set; empty otherwise.
    ///
    /// Needed because the temporal accumulators are a *cross-copy* statistic: the spread at
    /// boundary `k` is over the copies, so each copy's own history must survive until they can
    /// be compared. `AzOut::state` is the final state only — `cart` is overwritten in place
    /// every boundary — so there was no extension point and this is one.
    ///
    /// **Ragged by construction.** Copies terminate at different boundaries under
    /// `stop_on_event`, so these vectors have different lengths and any reader must handle that
    /// or it silently reads short. `stats::event_class_at` already does, with `tight.get(k)`
    /// and a terminal fallback.
    pub boundary_shapes: Vec<[T; 3]>,
    /// Energy drift at each completed sync boundary, when `AzOpts::keep_drift_hist` is set.
    ///
    /// Same cadence and same index as [`Self::refs`], so `refs[k] != refs[k-1]` selects exactly
    /// the boundaries at which the Levi-Civita registration was re-derived, and the paired
    /// increment `drift_hist[k] - drift_hist[k-1]` is what that switch cost.
    pub drift_hist: Vec<T>,
    /// Whether the **escape condition holds** at each completed sync boundary.
    ///
    /// Not "has escaped": the instantaneous candidacy, sampled on the same cadence as
    /// [`Self::tight`]. It exists because the escape condition turned out **not to be
    /// absorbing**: `escape_candidate` is relative energy `> 0` and receding, and in a
    /// collision-rich region a pair can be transiently unbound and then re-bind. Measured in
    /// `deep interior`, **885 of 895** trajectories that escape only under an in-loop test are
    /// re-bound one sync interval later.
    ///
    /// A persistence guard therefore cannot be written against a single later instant, and this
    /// is the record it has to be written against instead -- one integration, one
    /// discretisation, the whole history. Reading persistence by re-running to `t_e + w` with
    /// `n_sync` rescaled makes every window a different discretisation, which is the same defect
    /// as holding `n_sync` fixed while `t_max` varies.
    pub escape_flags: Vec<bool>,
    /// `|n(t_k) - n(t_{k-closure_k})|` at each completed boundary, `NaN` until the window is
    /// full. **`NaN`, never `0`** — an unfilled window is "not yet determined", and `NaN < tau`
    /// is false so it cannot fire, where a `0` would read as perfectly settled and fire on
    /// everything at `t ~ 0`.
    ///
    /// Recorded whatever the rule, so the closure distribution can be measured on trajectories
    /// the criterion did not label — which is what setting `tau` from a measured gap needs.
    pub closure_hist: Vec<T>,
    /// Per body, is it unbound from the other two, at each completed boundary?
    ///
    /// **The persistence question is "does it RE-BIND", and that is the energy arm alone.**
    /// [`Self::escape_flags`] records full criterion candidacy, which under `Closure` also
    /// carries the closure gate -- and closure is a difference of neighbouring samples, so it
    /// fluctuates above `tau` on a perfectly settled escape. Reading persistence off it would
    /// score ordinary jitter as a re-binding and report a correct criterion as broken.
    pub unbound_flags: Vec<[bool; 3]>,
    /// The escape's energy condition already held on **entry** to the interval in which the
    /// conjunction fired, so the replay found no crossing and `t_end` is the entry boundary.
    ///
    /// The number that says whether the `t_end` refinement is doing anything. See the firing
    /// block in `integrate_az_opts`.
    pub t_end_at_entry: bool,
    /// Smallest `A*B` seen at any step, **before** the `T::TINY` floor is applied.
    ///
    /// `dtau` divides by this, so it is the quantity the blow-up is a function of. Recorded
    /// rather than argued about: the floor is `1e-300` at f64 and `1e-37` at f32, and
    /// `TINY*TINY` **underflows at f32** -- so a doubly-degenerate state gives `dtau = inf`,
    /// caught by the explicit `is_finite` test and *not* by the floor the guard is named for.
    pub ab_min: T,
    /// Whether either factor was ever clamped to `T::TINY`. Says which guard did the work.
    pub ab_floored: bool,
    /// Time of the first terminating event, or the time reached if none fired. Distinct from
    /// `t` only when `stop_on_event` is off, where the run continues past the event.
    pub t_end: T,
}

/// The options that grew past what a positional argument list can carry legibly.
#[derive(Clone, Copy, Debug)]
pub struct AzOpts<'a, T> {
    pub forced_refs: Option<&'a [u8]>,
    /// Conditioned inverse LC branch. `false` reproduces the reference's original.
    pub lc_stable: bool,
    /// **Canonical**: a fraction of the initial hyperradius `R`, evaluated once at `t = 0`
    /// from this trajectory's own initial condition and never updated. Zero disables
    /// collision detection.
    ///
    /// A fraction, not a length. An absolute `r_coll` would break the scale invariance the
    /// project quotients out — measured, the same physical system gave answers differing by
    /// 1.66x purely from an arbitrary overall size. A *co-moving* one would make the
    /// Hamiltonian time-dependent and destroy energy conservation.
    pub r_coll_frac: T,
    /// Stop at the first terminating event, per §2.4. With it off the run continues to
    /// `t_max` and the event is recorded but not acted on — which is what the reference does,
    /// and what keeps every copy's continuous fields evaluated at a common playhead.
    pub stop_on_event: bool,
    /// Terminate on **escape**. Default **off**, and separately from [`Self::stop_on_event`],
    /// which continues to govern collision.
    ///
    /// **Collision stays terminal unconditionally** — it is a singularity, and there is genuinely
    /// nothing to integrate past. Escape is a *heuristic*, and freezing on one is what built the
    /// patchwork: a terminated trajectory's `shape_vec` is the state at its own `t_end`, escape
    /// is detected at sync boundaries, so the rendered field became a mosaic of `n_sync` time
    /// strata stitched at hard seams. A ribbon is a level set of `theta` **at a common time**;
    /// where neighbouring pixels froze at different times it breaks and resumes.
    ///
    /// Under [`crate::outcome::EscapeRule::Closure`] freezing should barely matter -- it fires on
    /// a trajectory whose shape has already stopped moving, so the state it freezes is the state
    /// it would have had anyway. That is the prediction, and it is measured before this flips.
    pub stop_on_escape: bool,
    /// Which escape condition to use. See [`crate::outcome::EscapeRule`] for the three and
    /// where each comes from.
    ///
    /// **[`Distance`](crate::outcome::EscapeRule::Distance) carries a canonical fraction**, not a
    /// length: a multiple of the initial hyperradius `R` evaluated once at `t = 0`, which this
    /// driver multiplies out at `r0`. The reference's saved configs use `rEsc = 5` and `12` as
    /// absolute lengths in the latent decode's own units, and every latent decode is normalised
    /// to `M = 1` with `I = 1` an algebraic identity — so `R = 1` there and the literal transfers.
    /// It does not transfer to Burrau's near-field, whose `R = 2.2361`. Expressing it as a
    /// fraction is the scale-gauge rule (BRIEF §2.5) at a third constant.
    ///
    /// **[`Closure`](crate::outcome::EscapeRule::Closure)'s `tau` is dimensionless** — a chord on
    /// the unit sphere — so it needs no conversion. Its *window* does, and that is [`Self::closure_k`].
    ///
    /// The body set is a property of the rule, each matching its own reference exactly, rather
    /// than a separate knob: `Reference` labels only the body outside the tightest pair, the
    /// other two test all three. That arm was measured on its own axis in the previous round and
    /// it is **not free** — it alone moved near-field's ungated escape fraction 0.0000 -> 1.0000.
    pub escape_rule: crate::outcome::EscapeRule<T>,
    /// The closure window, in **sync boundaries**. `|n(t_k) - n(t_{k-closure_k})|`.
    ///
    /// The reference measures a 0.4 time-unit window and reads only the two **ends** of it
    /// (`reference/escape_criterion.py` buffers `nbuf` samples and uses `buf[-1]`, `buf[0]`), so
    /// boundary sampling is a transcription rather than an approximation — at `t_max = 13,
    /// n_sync = 32` the realised window is 0.406 against 0.400.
    ///
    /// **It is a time, and `n_sync` is not scaled with `t_max`.** At `t_max = 50, n_sync = 32`
    /// the interval is 1.5625, so `closure_k = 1` there is a 3.9x wider window and a different
    /// criterion. The standing `n_sync`/`t_max` trap, landing on the new constant: scale `n_sync`
    /// with `t_max` or state the realised window with every figure.
    pub closure_k: usize,
    /// How often the **escape** test runs inside the RK4 loop, in steps. `0` is the
    /// reference's cadence: boundaries only.
    ///
    /// **This is the one place `t_end` is quantised, and it is a rendering defect as well as a
    /// measurement one.** Collision is sampled inside the loop (`tc = t + s.t`) and carries
    /// RK4-step resolution; escape is sampled only where the state is already Cartesian and
    /// every trajectory shares a playhead, which is what the reference does. So with
    /// `n_sync = 32` at `t_max = 13`, an escape-terminated `t_end` takes **32 possible values
    /// across a whole chart**, and any field derived from it renders those steps as concentric
    /// contour bands.
    ///
    /// That matters exactly where escape terminates. On Burrau's near-field at `t = 13` the
    /// escape arm is silent and every termination is a collision, so `t_end` there is already
    /// continuous. On the latent charts `escape_fraction` runs 0.9894-1.0000, so essentially
    /// every `t_end` is a boundary time. **The prediction is that the banding appears on the
    /// second set and not the first**, and it is measured rather than argued.
    ///
    /// **Default `0`, and it must stay there.** Turning this on changes results: the
    /// cross-check against `reference/tb_az.py` and the horizon table were both measured at
    /// the coarse cadence, and the reference has no finer one to compare against. This is a
    /// spec change behind a flag, not a tidy-up.
    pub escape_every: usize,
    /// Require an **in-loop** escape detection to still hold at the next sync boundary before
    /// it is accepted as terminal.
    ///
    /// **Measured, and it is not a precaution — without it the in-loop test is simply wrong in a
    /// collision-rich region.** `escape_candidate` is relative energy `> 0` and receding, and
    /// during a close encounter that is transiently true. In `deep interior`, of the 895
    /// trajectories that escape under `escape_every = 1` and not at the reference cadence,
    /// **0.000 are still unbound one boundary later** — and 0.000 at +2, +3, +4 and +8. All 895
    /// are transients, and latching them took the escape fraction from 0.0945 to 0.5494.
    ///
    /// It applies to in-loop detections **only**. Boundary detections are the reference's own
    /// arm and are left exactly as they are; on the latent charts, where escape genuinely
    /// terminates, the finer stride adds **zero** new escapes, so there is nothing there for a
    /// guard to catch and nothing for it to break.
    ///
    /// The committed time is the **first crossing**, not the confirmation — the guard decides
    /// whether the event is real, not when it happened.
    ///
    /// Same shape as `spread_event_latched`'s `LATCH_RUN`, which exists for this reason on a
    /// different field; the convention is reused rather than a new one invented.
    pub escape_confirm: bool,
    /// How `dtau` is sized within an interval. See [`DtauMode`]; the default is the fix and
    /// [`DtauMode::FixedPerInterval`] is the behaviour every committed number was taken under.
    pub dtau_mode: DtauMode,
    /// Clamp the **final** step of each sync interval so the state lands *on* the boundary
    /// instead of past it. Default **on**; `false` is the behaviour every committed number
    /// was taken under, kept for the same reason [`DtauMode::FixedPerInterval`] is.
    ///
    /// Without it the march overshoots — the loop exits on `s.t >= dt_left` — and only the
    /// *clock* is corrected (`t += s.t.min(dt_left)`), while the Cartesian state written back
    /// is the overshot one. That is a **first-order** error injected at every boundary. On the
    /// Chenciner-Montgomery figure-eight it dominates the whole integration: measured over one
    /// period at `n_sync = 32`, closure error `2.976e-3` against `3.595e-9` at `eta = 1e-3`, and
    /// the convergence order across `eta in [0.02, 0.001]` **1.13 against 3.06**. Read the order,
    /// not the error -- an error falls for many reasons; only the order says the leading term
    /// changed.
    ///
    /// **It is the partner of [`DtauMode::PerStepInterval`] and must not ship without it.**
    /// Under `FixedPerInterval` the overshoot is a fixed slice of fictitious time, so
    /// neighbouring trajectories overshoot alike and the error is large but spatially *smooth*.
    /// Under `PerStepInterval` the last step's size is a function of local `A*B`, so the
    /// overshoot varies pixel to pixel — a spatially-varying error injected at every boundary.
    /// Fixing the step control alone trades a smooth large error for a structured one.
    ///
    /// **The nested-arc banding this was first proposed to explain is NOT caused by it.** All
    /// four arms carry it, including the one predating both changes (`RESULTS §24.8`), and under
    /// outcome-class colouring it vanishes — a colouring artefact, per `RESULTS §21`. The defect
    /// here is real and independently measured; it is not the cause of that appearance.
    pub clamp_final_step: bool,
    /// Record the shape vector at every sync boundary, for the temporal accumulators (§5).
    ///
    /// Off by default: `n_sync` triples per copy is ~70x the size of a `PixelOut`, and it is
    /// reduced and dropped inside one footprint's evaluation, so the peak cost is one
    /// footprint's worth rather than the tree's.
    pub keep_boundary_shapes: bool,
    /// Record `|E(t) - E(0)| / |E(0)|` at every sync boundary.
    ///
    /// Off by default and gated for the same reason as [`Self::keep_boundary_shapes`]: it is a
    /// diagnostic, and a production run should not pay for it. It exists to answer one question
    /// the aggregate `drift` cannot -- whether drift arrives **at** the reference-body switches
    /// or accumulates smoothly between them. `AzOut::refs` already carries the switch record at
    /// this same cadence, so the two series line up index for index.
    pub keep_drift_hist: bool,
}

impl<T: Real> Default for AzOpts<'_, T> {
    fn default() -> Self {
        Self {
            forced_refs: None,
            lc_stable: true,
            r_coll_frac: T::zero(),
            escape_rule: crate::outcome::EscapeRule::Reference,
            closure_k: 1,
            stop_on_event: true,
            stop_on_escape: false,
            escape_every: 0,
            escape_confirm: true,
            dtau_mode: DtauMode::default(),
            clamp_final_step: true,
            keep_boundary_shapes: false,
            keep_drift_hist: false,
        }
    }
}

/// `forced_refs`, when supplied, overrides the per-sync choice — this is how the shared
/// reference policy is applied, by handing every copy the nominal copy's choices.
pub fn integrate_az<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    forced_refs: Option<&[u8]>,
) -> AzOut<T> {
    integrate_az_lc(s0, m, t_max, n_sync, eta, max_steps, forced_refs, true)
}

/// As [`integrate_az`], with an explicit choice of inverse LC branch.
///
/// `lc_stable = false` reproduces the reference's branch exactly and is what the cross-check
/// uses; `true` is the production kernel. The two produce genuinely different trajectories —
/// the stable branch is more accurate, so it necessarily stops agreeing bit-for-bit with the
/// reference. Keeping both makes that a measured difference rather than a lost gate.
#[allow(clippy::too_many_arguments)]
pub fn integrate_az_lc<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    forced_refs: Option<&[u8]>,
    lc_stable: bool,
) -> AzOut<T> {
    // No r_coll and no early stop: these two entry points are the reference-matching path,
    // and the reference has no termination logic at all. Events are still *recorded* — that
    // costs nothing and changes no arithmetic — but nothing acts on them.
    integrate_az_opts(
        s0, m, t_max, n_sync, eta, max_steps,
        &reference_opts(forced_refs, lc_stable, DtauMode::default(), true),
    )
}

/// The reference-matching option set, in one place so the cross-check binary can vary the one
/// axis it needs to without restating the other nine.
///
/// The reference has no distance gate, no closure arm and no persistence guard: the cross-check
/// measures `EscapeRule::Reference` and this hardcodes it, so nothing added since can reach that
/// path. `dtau_mode` is the exception and is a parameter, because **both sides carry the same
/// defect and both are being fixed** -- the comparison that means anything is old-against-old and
/// new-against-new, and either alone would pass while the two transcriptions diverged.
pub fn reference_opts<T: Real>(
    forced_refs: Option<&[u8]>,
    lc_stable: bool,
    dtau_mode: DtauMode,
    clamp_final_step: bool,
) -> AzOpts<'_, T> {
    AzOpts {
        forced_refs,
        lc_stable,
        r_coll_frac: T::zero(),
        escape_rule: crate::outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_event: false,
        stop_on_escape: false,
        escape_every: 0,
        escape_confirm: true,
        dtau_mode,
        clamp_final_step,
        keep_boundary_shapes: false,
        keep_drift_hist: false,
    }
}

/// The step size for the *next* RK4 step, given the interval's entry sizing and the current
/// regularised state.
///
/// Returns `(dtau_now, ab_raw, floored)` -- `ab_raw` is `A*B` **before** the `T::TINY` floor,
/// which is what the blow-up is a function of and what the `TINY` report reads.
#[inline]
fn dtau_for_step<T: Real>(
    mode: DtauMode,
    clamp_final: bool,
    eta: T,
    dt_left: T,
    dtau_entry: T,
    s: &super::state::AzState<T>,
) -> (T, T, bool) {
    let a = s.a();
    let b = s.b();
    let ab_raw = a * b;
    let floored = a < T::TINY || b < T::TINY;
    let ab = a.max(T::TINY) * b.max(T::TINY);
    let dtau = match mode {
        DtauMode::FixedPerInterval => dtau_entry,
        // `rem`, not `dt_left`: geometric decay, kept so the measurement can show it.
        DtauMode::PerStepRemaining => {
            let rem = (dt_left - s.t).max(T::zero());
            eta * rem / ab
        }
        DtauMode::PerStepInterval => (eta * dt_left / ab).min(dtau_entry),
    };
    // Land ON the boundary rather than past it. `dt = A*B*dtau` to leading order, so
    // `(dt_left - s.t)/ab` is the fictitious time still owed; `.min` makes this a one-sided
    // reduction of the last step only and leaves every interior step untouched. **`ab` is the
    // same floored product the mode used** -- recomputing it here would let the clamp and the
    // step disagree. The landing is exact only to the accuracy with which `A*B` predicts the time
    // increment over the step -- first order -- so the residual is `O(h^2)` per boundary and the
    // measured global order is **3.06 under `FixedPerInterval` and 2.08 under `PerStepInterval`**,
    // against 1.13 and 1.06 without. Not the stepper's own four, and stated rather than assumed.
    let dtau = if clamp_final {
        dtau.min((dt_left - s.t).max(T::zero()) / ab)
    } else {
        dtau
    };
    (dtau, ab_raw, floored)
}

/// Refine an escape's `t_end` by **replaying** the sync interval it fired in.
///
/// The closure criterion is evaluated at boundaries, so a firing at `t_k` says only that the
/// conjunction holds by then. Closure itself has no finer resolution — it is defined from the
/// boundary series — but the *energy* arm does, and if it crossed inside `(t_prev, t_k]` that
/// crossing is the honest escape time. Re-run the interval with the same reference body, the same
/// `dtau` rule and the same stepper, and take the first sub-step at which body `b` is unbound.
///
/// Returns the boundary time unchanged when there is no saved entry state (the first interval),
/// and the **entry** time with `at_entry` set when the body was already unbound on entry — energy
/// flickers, closure is what just settled, and there is then no crossing inside the interval.
#[allow(clippy::too_many_arguments)]
fn refine_escape_time<T: Real>(
    m: &[T; 3],
    prev: &Option<(Cart<T>, T, usize, T)>,
    b: usize,
    t_boundary: T,
    lc_stable: bool,
    mode: DtauMode,
    clamp_final: bool,
    eta: T,
    at_entry: &mut bool,
) -> T {
    let Some((c0, t0, a, dtau)) = *prev else {
        return t_boundary;
    };
    if crate::outcome::unbound(&c0, m, b) {
        // **Keep the boundary time, not the entry time.** The criterion is the CONJUNCTION, and
        // it first held here; the energy arm crossed at some earlier, unknown boundary. Reporting
        // `t0` would claim an escape at a time the criterion had not yet concluded one.
        *at_entry = true;
        return t_boundary;
    }
    let (ab, bb, cb) = triple(a);
    let sys = if lc_stable {
        AzSystem::new(ab, bb, cb, *m)
    } else {
        AzSystem::new(ab, bb, cb, *m).with_reference_lc()
    };
    let (mut s, e) = sys.to_reg(&c0);
    let dt_left = t_boundary - t0;
    let mut steps = 0usize;
    // The same budget shape as the main loop: bounded, and `is_finite` tested explicitly so a
    // diverged replay cannot spin (NaN >= x is false).
    // **The replay must step exactly as the march it is replaying.** A different step rule here
    // returns a `t_end` that is wrong for a reason that looks like a criterion result -- and the
    // escape criterion is precisely what this quantity is read as evidence about.
    while s.is_finite() && s.t < dt_left - land_tol(clamp_final, dt_left) && steps < REPLAY_MAX_STEPS {
        let (h, _, _) = dtau_for_step(mode, clamp_final, eta, dt_left, dtau, &s);
        s = rk4::step(&sys, &s, e, h);
        steps += 1;
        if crate::outcome::unbound(&sys.to_cartesian(&s), m, b) {
            return t0 + s.t;
        }
    }
    t_boundary
}

/// Step ceiling for [`refine_escape_time`]. One sync interval at the driver's own `dtau` is
/// `~1/eta` steps; this is far above that and exists only so a pathological replay terminates.
const REPLAY_MAX_STEPS: usize = 100_000;

/// The interval is complete once the state is within this of the boundary.
///
/// Zero without the clamp — the loop exits by *overshooting*, so any positive tolerance would
/// change which step is the last one. With it the final step lands to the stepper's own order,
/// and without a tolerance the residual would be paid for by a cascade of ever-tinier steps
/// that cannot reach equality in floating point.
///
/// **Relative to `dt_left`, never absolute.** All times rescale by `alpha^{3/2}` under the
/// project's scale gauge, so an absolute slack is a different tolerance at every scale — measured,
/// it broke the bitwise scale-invariance test at `4.24e-15`.
#[inline]
fn land_tol<T: Real>(clamp_final: bool, dt_left: T) -> T {
    if clamp_final {
        dt_left * T::LAND_EPS_REL
    } else {
        T::zero()
    }
}

/// The full entry point.
#[allow(clippy::too_many_arguments)]
pub fn integrate_az_opts<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &AzOpts<T>,
) -> AzOut<T> {
    let (forced_refs, lc_stable) = (opts.forced_refs, opts.lc_stable);
    let mut cart = s0;
    let e0 = energy::energy(&s0.r, &s0.v, m, T::zero());
    let mut t = T::zero();
    let mut d_min_ref = T::infinity();
    let mut d_min_true = T::infinity();
    let mut gamma_max = T::zero();
    let mut switches = 0u32;
    let mut prev_ref: Option<usize> = None;
    let mut refs = Vec::with_capacity(n_sync);
    let mut tight = Vec::with_capacity(n_sync);
    let mut escape_flags = Vec::with_capacity(n_sync);
    let mut tie_ratio = Vec::with_capacity(n_sync);
    let mut boundary_shapes: Vec<[T; 3]> =
        Vec::with_capacity(if opts.keep_boundary_shapes { n_sync } else { 0 });
    let mut drift_hist: Vec<T> = Vec::with_capacity(if opts.keep_drift_hist { n_sync } else { 0 });
    let mut total_steps = 0usize;
    let mut finite = true;
    let mut budget_exhausted = false;
    let mut events = Events::default();
    let mut t_end: Option<T> = None;
    // An in-loop escape awaiting confirmation at the next boundary. See `AzOpts::escape_confirm`.
    let mut pending_escape: Option<(u8, T)> = None;

    // Canonical and fixed at t=0: a fraction of *this* trajectory's initial hyperradius,
    // evaluated once, before anything moves. Never recomputed from the instantaneous
    // configuration — a co-moving length makes the Hamiltonian time-dependent.
    let r0 = energy::hyperradius(&s0.r, m);
    let r_coll = opts.r_coll_frac * r0;
    // `Distance`'s gate is canonical too — a fraction of R fixed at t = 0, multiplied out here
    // and never recomputed from the instantaneous configuration.
    let rule = opts.escape_rule;
    let esc = |c: &Cart<T>, cl: Option<T>| {
        crate::outcome::escape_candidate_rule(c, m, rule, r0, cl)
    };
    // The closure window's ring buffer: the last `closure_k + 1` boundary shape vectors. Pushed
    // unconditionally rather than under `keep_boundary_shapes`, because the criterion reads it
    // and a criterion that depends on a diagnostic flag is a criterion with two behaviours.
    let kw = opts.closure_k.max(1);
    let mut nbuf: std::collections::VecDeque<[T; 3]> = std::collections::VecDeque::with_capacity(kw + 1);
    // The closure value at each boundary, NaN until the window is full. `NaN < tau` is false, so
    // an unfilled window cannot fire — correct, since nothing has settled at t ~ 0.
    let mut closure_hist: Vec<T> = Vec::with_capacity(n_sync);
    let mut unbound_flags: Vec<[bool; 3]> = Vec::with_capacity(n_sync);
    // The Cartesian state and time at the PREVIOUS boundary, for the replay refinement of
    // `t_end`. See the firing block below.
    // Always written before the boundary block that reads it -- the `None` is unreachable and
    // exists so `refine_escape_time` has a total signature rather than an unwrap.
    #[allow(unused_assignments)]
    let mut prev_boundary: Option<(Cart<T>, T, usize, T)> = None;
    let mut t_end_at_entry = false;
    let mut ab_min = T::infinity();
    let mut ab_floored = false;

    for kk in 0..n_sync {
        let t_target = T::lit((kk + 1) as f64) * t_max / T::lit(n_sync as f64);

        // `f.get(kk)`, not `f[kk]`. Since Step 5b the nominal copy can terminate early, so its
        // `refs` record is shorter than `n_sync` and the shared policy has no opinion past
        // that point. Falling back to the per-copy choice is the only defensible reading:
        // sharing applies where the nominal has a choice to share. (This indexed out of
        // bounds the first time the shared policy met a terminating run.)
        let a = match forced_refs.and_then(|f| f.get(kk)) {
            Some(&f) => f as usize,
            None => choose_reference(&cart.r),
        };
        refs.push(a as u8);
        if let Some(p) = prev_ref {
            if p != a {
                switches += 1;
            }
        }
        prev_ref = Some(a);

        // Matches the reference's row selection `t < t_target - 1e-15`.
        if t >= t_target - T::SYNC_EPS {
            continue;
        }

        let (ab, bb, cb) = triple(a);
        let sys = if lc_stable {
            AzSystem::new(ab, bb, cb, *m)
        } else {
            AzSystem::new(ab, bb, cb, *m).with_reference_lc()
        };
        let (mut s, e) = sys.to_reg(&cart);
        let dt_left = t_target - t;

        // dt = A*B*dtau, so the A*B factor ALREADY shrinks the physical step at close
        // approach — that is what regularisation buys. dtau must therefore NOT also shrink
        // with the separation: doing so drove dt -> 1e-13 and exhausted the step budget,
        // producing a false "this region is intractable". Size dtau so the FIRST physical
        // step is a fixed fraction of the interval; the rest follows.
        // The interval's **entry** sizing. Under `FixedPerInterval` it is the whole step
        // control; under `PerStepInterval` it is the one-sided cap that keeps the physical step
        // shrinking at a close approach instead of being held at `eta*dt_left`.
        let a0 = s.a().max(T::TINY);
        let b0 = s.b().max(T::TINY);
        let dtau = eta * dt_left / (a0 * b0);

        // The interval's entry state, kept for the replay refinement of `t_end`. `cart` is
        // overwritten at the boundary below, so it has to be taken now.
        prev_boundary = Some((cart, t, a, dtau));

        let mut steps = 0usize;
        // Whether the march reached the boundary, as against breaking out non-finite or out of
        // budget. Only then may the clock be set to the boundary exactly.
        let mut landed = false;
        loop {
            // NaN >= x is false, so a non-finite trajectory never satisfies `done` and the
            // loop would burn the whole budget (measured 354 s against 3 s nominal).
            // Test is_finite explicitly.
            let bad = !s.is_finite();
            if bad {
                finite = false;
            }
            if s.t >= dt_left - land_tol::<T>(opts.clamp_final_step, dt_left) || bad {
                landed = !bad;
                break;
            }
            if steps >= max_steps {
                budget_exhausted = true;
                break;
            }

            let (h, ab_raw, floored) =
                dtau_for_step(opts.dtau_mode, opts.clamp_final_step, eta, dt_left, dtau, &s);
            if ab_raw.is_finite() && ab_raw < ab_min {
                ab_min = ab_raw;
            }
            ab_floored |= floored;
            s = rk4::step(&sys, &s, e, h);
            steps += 1;

            let (r1, r2, _, _) = sys.phys_from_state(&s);
            let d1 = r1.norm();
            let d2 = r2.norm();
            let d3 = (r2 - r1).norm();
            // Only fold finite separations in. `min(x, NaN)` is NaN in the reference, which
            // silently poisons d_min for every consumer downstream (NOTES §2.3). The copy
            // is still kept — `finite` records the outcome.
            let m_ref = d1.min(d2);
            if m_ref.is_finite() && m_ref < d_min_ref {
                d_min_ref = m_ref;
            }
            let m_true = m_ref.min(d3);
            if m_true.is_finite() && m_true < d_min_true {
                d_min_true = m_true;
            }

            // Sampled here, inside the RK4 loop, not at sync boundaries: with n_sync = 32 and
            // t_max = 13 the boundaries are 0.4 apart, and a close encounter passes between
            // two of them entirely unseen.
            if events.collision.is_none() && r_coll > T::zero() {
                let mask = collision_pairs_from(ab, bb, cb, d1, d2, d3, r_coll);
                if mask != 0 {
                    let tc = t + s.t;
                    events.collision = Some((mask, tc));
                    t_end.get_or_insert(tc);
                    if opts.stop_on_event {
                        break;
                    }
                }
            }

            // The escape test at RK4-step resolution, when asked for. Off by default: this is
            // the reference's cadence and changing it changes results. `to_cartesian` per
            // tested step is the cost, which is why it is strided rather than unconditional.
            // Under `Closure` this arm is **off** regardless of `escape_every`: closure is only
            // defined where the state is Cartesian, so there is no value to test between
            // boundaries, and the window is already the persistence guard `escape_confirm`
            // exists to be. Both flags become vacuous under that rule; said, not silently dropped.
            if opts.escape_every > 0
                && !matches!(rule, crate::outcome::EscapeRule::Closure(_))
                && events.escape.is_none()
                && pending_escape.is_none()
                && steps % opts.escape_every == 0
            {
                let c = sys.to_cartesian(&s);
                if let Some(b) = esc(&c, None) {
                    let te = t + s.t;
                    if opts.escape_confirm {
                        // Provisional. Do NOT break: the trajectory has to reach the next
                        // boundary for the condition to be re-tested, and breaking here is
                        // what turned 895 transients into terminal escapes.
                        pending_escape = Some((b, te));
                    } else {
                        events.escape = Some((b, te));
                        t_end.get_or_insert(te);
                        if opts.stop_on_escape {
                            break;
                        }
                    }
                }
            }

            let g = gamma_residual(&sys, &s, e);
            if g.is_finite() && g > gamma_max {
                gamma_max = g;
            }
        }
        total_steps += steps;

        cart = sys.to_cartesian(&s);
        // Under `clamp_final_step` the state IS the boundary state to within `LAND_EPS_REL`, so the
        // clock is set to the boundary exactly and the two cannot disagree. Without it, the
        // overshoot past the boundary is clipped in the time bookkeeping *only* and the state
        // written back is the overshot one -- the reference's behaviour, and a first-order error.
        t += if landed && opts.clamp_final_step {
            dt_left
        } else {
            s.t.min(dt_left)
        };

        tight.push(crate::outcome::binary_id(&cart));
        // Computed unconditionally: it is ~20 flops against a whole RK4 interval, the closure
        // criterion reads it, and gating a criterion on a diagnostic flag gives it two behaviours.
        let n_now = crate::physics::shape::shape_vec(&cart.r, m);
        if opts.keep_boundary_shapes {
            boundary_shapes.push(n_now);
        }
        if opts.keep_drift_hist {
            let ek = energy::energy(&cart.r, &cart.v, m, T::zero());
            drift_hist.push(((ek - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs());
        }
        // The closure window reads only the two ENDS — `buf[-1]` and `buf[0]` in the reference,
        // never the interior — so a ring buffer of `kw + 1` is the whole state it needs.
        nbuf.push_back(n_now);
        if nbuf.len() > kw + 1 {
            nbuf.pop_front();
        }
        let cl = if nbuf.len() == kw + 1 {
            Some(crate::outcome::closure(&n_now, &nbuf[0]))
        } else {
            None
        };
        closure_hist.push(cl.unwrap_or(T::nan()));
        unbound_flags.push([
            crate::outcome::unbound(&cart, m, 0),
            crate::outcome::unbound(&cart, m, 1),
            crate::outcome::unbound(&cart, m, 2),
        ]);
        {
            let mut d = crate::physics::newton::pair_dists(&cart.r);
            d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            tie_ratio.push(d[1] / d[0].max(T::TINY));
        }

        // Instantaneous candidacy at this boundary, recorded whether or not it has already
        // fired — this is the history a persistence guard reads, and it must not stop being
        // written once `events.escape` is set.
        let candidate_now = esc(&cart, cl);
        escape_flags.push(candidate_now.is_some());

        // Confirm or discard a provisional in-loop escape. The committed time is the FIRST
        // crossing, not this boundary: the guard decides whether the event was real, not when
        // it happened.
        if let Some((b, te)) = pending_escape.take() {
            if candidate_now.is_some() {
                events.escape = Some((b, te));
                t_end.get_or_insert(te);
            }
        }

        // The escape test runs at the sync boundary, where the state is Cartesian and every
        // trajectory shares a playhead — the reference's cadence, transcribed.
        if events.escape.is_none() {
            if let Some(b) = candidate_now {
                // **Refine `t_end` by replay, not by sampling finer.** The conjunction first
                // holds at this boundary, so the energy crossing that completes it lies in
                // `(t_prev, t]`. Replay that one interval from the saved entry state with the
                // same stepper and the same `dtau`, and take the first sub-step at which body
                // `b` is unbound. One extra interval per escaping trajectory, no second
                // sampling path in the main loop, and no new numerics.
                //
                // **The honest limit, and it is counted rather than assumed:** closure's own
                // resolution is still the boundary cadence, and `spec > 0` may already hold on
                // ENTRY to the interval — energy flickers, closure is what just settled. Then
                // there is no crossing to find, `t_end` is the entry time, and
                // `t_end_at_entry` records it. If that fires on nearly everything the
                // refinement is decoration, and the measurement says so.
                //
                // **Only under `Closure`.** There the boundary cadence is intrinsic -- closure is
                // defined from the boundary series and has no finer resolution -- and energy is
                // the one continuous arm. Under `Reference`/`Distance` both arms are continuous,
                // so refining on energy alone would return a time at which the full condition may
                // not yet hold, and `escape_every` is the knob that already samples those finer,
                // with committed semantics and a cross-check measuring them.
                let te = if matches!(rule, crate::outcome::EscapeRule::Closure(_)) {
                    refine_escape_time(
                        m, &prev_boundary, b as usize, t, lc_stable, opts.dtau_mode,
                        opts.clamp_final_step, eta, &mut t_end_at_entry,
                    )
                } else {
                    t
                };
                events.escape = Some((b, te));
                t_end.get_or_insert(te);
            }
        }

        if budget_exhausted
            || (opts.stop_on_event && events.collision.is_some())
            || (opts.stop_on_escape && events.escape.is_some())
        {
            break;
        }
    }

    let e1 = energy::energy(&cart.r, &cart.v, m, T::zero());
    let drift = ((e1 - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs();

    AzOut {
        state: cart,
        t,
        drift,
        d_min_ref,
        d_min_true,
        switches,
        gamma_max,
        refs,
        tight,
        tie_ratio,
        escape_flags,
        closure_hist,
        unbound_flags,
        t_end_at_entry,
        ab_min,
        ab_floored,
        boundary_shapes,
        drift_hist,
        steps: total_steps,
        finite,
        budget_exhausted,
        events,
        // **The EARLIEST recorded event, not the first one inserted.** `get_or_insert` gave
        // whichever arm happened to be sampled first, which is the same ordering error
        // `classify` carried: escape is tested only where the state is Cartesian, so a
        // collision detected mid-interval could take `t_end` from an escape that preceded it.
        // Taking the min makes `t_end` agree with the state `classify` returns by construction.
        t_end: match (events.collision, events.escape) {
            (Some((_, a)), Some((_, b))) => a.min(b),
            (Some((_, a)), None) => a,
            (None, Some((_, b))) => b,
            (None, None) => t_end.unwrap_or(t),
        },
    }
}

/// Convenience: run with the shared-reference policy by first computing the nominal copy's
/// choices, then forcing them on the copy.
pub fn integrate_with_policy<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    policy: RefPolicy,
    nominal_refs: Option<&[u8]>,
) -> AzOut<T> {
    let forced = match (policy, nominal_refs) {
        (RefPolicy::Shared, Some(f)) => Some(f),
        _ => None,
    };
    integrate_az(s0, m, t_max, n_sync, eta, max_steps, forced)
}
