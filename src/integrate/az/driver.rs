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
/// How the competing step-size constraints are combined.
///
/// # Why a `min` is not neutral
///
/// The step size is `min(mode, entry cap, per-step limit)`. A `min` is C^0 but **not C^1**: where
/// the active constraint switches, the step size has a *crease*, and its derivative with respect
/// to the initial condition jumps. Accumulated over ~10^5 steps those creases can print as edges
/// in a field that is otherwise smooth — the constraint-switching surface is a codimension-1 sheet
/// in IC space, which is what an edge in a rendered field is.
///
/// # The construction, and why not the plain harmonic mean
///
/// `SoftMin` is the p-norm soft minimum `(sum x_i^-p)^(-1/p)`, which is C^infinity in every
/// argument and **always <= min**, so it can only ever make the step more conservative.
///
/// At `p = 1` it is exactly the reciprocal sum `1/(1/a + 1/b + 1/c)` — the harmonic form — and
/// that is where the cost lives: with `n` constraints *tied* it returns `min/n`, so three
/// comparable constraints give a **third** of the step and three times the work. The plain
/// harmonic *mean* `n/sum(1/x)` fixes the tie case and is **not conservative** — with one large
/// argument it exceeds the min, which is the one thing a step limit may never do.
///
/// `p` recovers both ends: `p -> inf` is the hard `min`, `p = 1` is the harmonic form, and at
/// `p = 4` three tied constraints cost `3^(1/4) = 1.32x` rather than `3x` while the map stays
/// smooth. **`p` is a measured parameter, not a chosen one** — see `results/step_control/`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum StepBlend {
    /// Today's behaviour, and what every committed number was taken under.
    #[default]
    Min,
    /// The p-norm soft minimum. `p` is [`AzOpts::blend_p`].
    SoftMin,
}

/// Which per-step limit, if any, bounds the step the march is about to take.
///
/// # Why this exists
///
/// `dt = A*B*dtau` is emergent: `dt/dtau = A*B` is integrated *by the RK4 stepper*, so the step
/// actually taken is not the one predicted. Measured on `config_stability`, **one step advanced
/// the physical clock by `2.209e128` against a sync interval of `0.4`** and the march recorded a
/// clean landing -- `1e128` is finite so the divergence guard passes, `s.t >= dt_left` is
/// satisfied by 128 orders, and `t += dt_left` then corrects the *clock* while keeping the
/// *state*. **Nothing asked whether the step it just took was one it could afford.**
///
/// The two batch remedies -- `refine_flagged`, and a global `eta/256` -- are characterisation,
/// not fixes: the first re-integrates from `t = 0`, which a live playhead cannot do, and the
/// second pays 256x everywhere for a local failure. Only a per-step mechanism survives contact
/// with marching. These are the four candidates, decided by measurement rather than by argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum StepLimit {
    /// Today's behaviour. **Every committed number in the corpus was taken under this**, so it
    /// stays the default and a test asserts it is bitwise unchanged.
    #[default]
    None,
    /// **A.** Take the step, test it, and if the largest separation moved by more than
    /// `f * d_min` restore and halve. Local and live-compatible; on a GPU a rejection is a
    /// divergent branch, which `examples/warp_divergence.rs` measures rather than assumes.
    Reject,
    /// **B.** Branch-free, never rejects: bound `dtau` by the closest pair's crossing time,
    /// `f * d_min / |v_rel|_max`, computed from values already in registers before stepping.
    /// `f` is a fraction of a crossing time and has physical meaning.
    Predictive,
    /// **C.** The narrowest patch: the blow-up is `A*B` growing mid-interval, so cap the growth
    /// at `f` times its interval-entry value. **Distinct from `PerStepInterval`'s cap**, which
    /// bounds `dtau` at its entry value rather than the product.
    AbGrowth,
    /// **D.** The dumb control: scale `eta` by `f` and pay everywhere. **This port has no
    /// substep-bucket table** -- `substep_bucket`, `N_sub`, `N_max` and descriptor bit 5 are the
    /// GLSL app's contract and appear nowhere in `src/`. The faithful stand-in for "widen the
    /// table" is the knob that buys resolution uniformly, and in an AZ march that is `eta`, since
    /// steps per interval go as `1/eta`. A stand-in, and labelled as one.
    Global,
}

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
    /// Second-**longest** over **longest** pair separation, at each completed sync boundary.
    ///
    /// Distinct from [`AzOut::tie_ratio`], which is about the two *tightest* pairs and decides
    /// which binary is which. This one is about the two *longest*, and it is the quantity
    /// `choose_reference` turns on: the reference is `THIRD[argmax d]`, so it flips exactly where
    /// the longest side changes identity — where this ratio reaches 1. It is the coordinate of
    /// the chart-switching surface, and the two must not be confused: reading `tie_ratio` for
    /// this question would measure the wrong pair entirely.
    pub ref_tie: Vec<T>,
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
    ///
    /// **This is an advance-anyway site**: the step is taken with a denominator the code knows is
    /// fabricated, and nothing terminates. Read it as saturation, not as a warning.
    pub ab_floored: bool,
    /// Largest **physical** step taken over the whole run, as an actual `s.t` difference.
    ///
    /// `ab_min` records the worst denominator; this records the step it produced. A step-control
    /// cliff shows here and in no other recorded quantity.
    pub dt_max: T,
    /// Steps taken at `PerStepInterval`'s `.min(dtau_entry)` cap -- the mode asked for a larger
    /// step than the interval's entry sizing and was refused. The landing clamp is excluded.
    pub n_cap_hits: u32,
    /// **THE TRIPWIRE.** Steps after which the interval-local clock exceeded `2 * dt_left`.
    ///
    /// `dt > dt_left` is a **bug**, not a condition to handle: a legitimate overshoot under fixed
    /// `dtau` is at most one nominal step, ~1% of `dt_left`, so `2x` is unambiguous. It went
    /// undetected for six days because nothing asserted it. `debug_assert`ed in debug and counted
    /// here in release -- and **conditioned on `clamp_final_step`**, because with the clamp off
    /// overshoot is the expected behaviour of a named measurement axis, and an assert that fires
    /// on a deliberate mode is a broken assert.
    pub n_overshoot: u32,
    /// Steps retried under [`StepLimit::Reject`]. Zero under every other mode.
    pub n_retry: u32,
    /// A step exhausted [`MAX_RETRIES`] under `Reject`. The trajectory is **undetermined**, not
    /// discarded -- a measurement outcome, and counted separately from `budget_exhausted` so a
    /// failure swapped for a differently-named failure is visible rather than absorbed.
    pub retry_exhausted: bool,
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
    /// Which per-step limit bounds the step. See [`StepLimit`].
    pub step_limit: StepLimit,
    /// **Hysteresis on the reference-body choice.** `0` is the plain `argmax`.
    ///
    /// The selector normally switches the instant another pair becomes the longest. With
    /// `eps > 0` it keeps the current reference until a rival exceeds the current reference's
    /// own opposite side by a factor `1 + eps`, so the switching surface **moves** and small
    /// perturbations no longer flip the chart back and forth.
    ///
    /// It is an intervention, not a proposed default: every trajectory remains a legitimate
    /// integration, but the itinerary is now path-dependent and the NumPy cross-check pins
    /// `0`. Its purpose is a falsifiable picture — **if the rendered wedges are chart-selection
    /// artefacts, displacing the surfaces displaces them; if they are dynamical, the field is
    /// invariant and only the itinerary changes.**
    pub ref_hysteresis: f64,
    /// How the competing constraints are combined. See [`StepBlend`].
    pub step_blend: StepBlend,
    /// The soft-minimum exponent. `1.0` is the harmonic form; large is the hard `min`.
    pub blend_p: f64,
    /// The limit's single parameter. Its meaning is **per mode** and deliberately not shared:
    /// a fraction of `d_min` for `Reject`, a fraction of a crossing time for `Predictive`, an
    /// `A*B` growth factor for `AbGrowth`, and an `eta` multiplier for `Global`. One number with
    /// four meanings is a knob that cannot be swept across modes, so the harness sweeps it
    /// per mode and prints which meaning is in force.
    pub step_limit_f: f64,
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
            step_limit: StepLimit::None,
            step_limit_f: 0.0,
            ref_hysteresis: 0.0,
            step_blend: StepBlend::Min,
            blend_p: 4.0,
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
        step_limit: StepLimit::None,
        step_limit_f: 0.0,
        ref_hysteresis: 0.0,
        step_blend: StepBlend::Min,
        blend_p: 4.0,
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

/// What [`dtau_for_step`] decided, and the two things it decided *despite*.
///
/// Returned as a struct rather than a tuple because two of the four fields exist only to be
/// recorded: **an advance-anyway site with no telemetry is indistinguishable from one that never
/// fires**, and both of these advance with a step the code knows is not the one it asked for.
#[derive(Clone, Copy)]
pub(crate) struct StepSize<T> {
    /// The fictitious step actually taken.
    pub dtau: T,
    /// `A*B` **before** the `T::TINY` floor.
    pub ab_raw: T,
    /// Either factor was clamped to `T::TINY`, so `dtau` divides by a fabricated denominator
    /// and the step advances anyway.
    pub floored: bool,
    /// `PerStepInterval`'s `.min(dtau_entry)` bound, and only that one -- **not** the landing
    /// clamp, which is a deliberate one-sided reduction of the final step. When this binds the
    /// mode wanted a *larger* step than the interval's entry sizing and was refused it.
    pub capped: bool,
}

/// Retries a single step may take under [`StepLimit::Reject`] before the trajectory is called
/// **undetermined**.
///
/// Halving eight times is a 256x reduction -- the same factor the `eta` ladder needed to bring
/// every flagged pixel to `error_ratio` 1.000, so a step that still fails here is not failing for
/// want of resolution.
pub const MAX_RETRIES: u32 = 8;

/// **B.** The closest pair's crossing time, as a `dtau` bound.
///
/// `dt = A*B*dtau`, so a physical bound `dt <= f*d_min/|v_rel|_max` is
/// `dtau <= f*d_min/(|v_rel|_max*A*B)`. `phys_from_state` returns `(R1, R2, V1, V2)` and the
/// third pair is the difference of each -- everything is already in registers, and the whole
/// limit is one divide with no trial step, no retry and no branch.
///
/// Returns `+inf` when nothing is moving, which is the correct absence of a bound rather than a
/// zero step.
fn predictive_dtau<T: Real>(
    sys: &AzSystem<T>,
    s: &super::state::AzState<T>,
    f: T,
) -> T {
    let (r1, r2, v1, v2) = sys.phys_from_state(s);
    let d_min = r1.norm().min(r2.norm()).min((r2 - r1).norm());
    let v_max = v1.norm().max(v2.norm()).max((v2 - v1).norm());
    let ab = s.a().max(T::TINY) * s.b().max(T::TINY);
    let denom = v_max * ab;
    // `f <= 0` is not "the tightest possible bound", it is a step of zero and a march that never
    // advances. A mode selected without its parameter would otherwise burn its whole budget and
    // report `budget_exhausted`, which reads as a physics failure. Treated as no bound.
    if f > T::zero() && denom > T::zero() && d_min.is_finite() && denom.is_finite() {
        f * d_min / denom
    } else {
        T::infinity()
    }
}

/// Combine the competing `dtau` constraints. See [`StepBlend`].
///
/// Non-finite and non-positive entries are skipped — `+inf` is how "this constraint is not in
/// force" is spelled, and a zero or negative bound is not a tighter limit but a broken one.
/// Returns `+inf` when nothing is in force, which the caller reads as no bound.
fn blend_dtau<T: Real>(blend: StepBlend, p: T, c: &[T]) -> T {
    let live = || c.iter().copied().filter(|x| x.is_finite() && *x > T::zero());
    match blend {
        StepBlend::Min => live().fold(T::infinity(), |a, b| if b < a { b } else { a }),
        StepBlend::SoftMin => {
            // (sum x^-p)^(-1/p). One `powf` per constraint; the constraints number two or three.
            let s: T = live().map(|x| x.powf(-p)).sum();
            if s > T::zero() && s.is_finite() {
                s.powf(-T::one() / p)
            } else {
                live().fold(T::infinity(), |a, b| if b < a { b } else { a })
            }
        }
    }
}

/// **C.** The `A*B` growth clamp, as a `dtau` bound.
///
/// `dt = min(A*B, ab_entry*C) * dtau_entry`, realised as a *reduction of `dtau`* because
/// `dt/dtau = A*B` is integrated by the stepper and cannot be multiplied after the fact. Relative
/// to the interval's **entry** sizing, so it composes with whichever [`DtauMode`] is in force
/// instead of replacing it — and `+inf` while `A*B` has not grown, which is the absence of a
/// bound rather than a step of zero.
fn ab_growth_dtau<T: Real>(ab_entry: T, ab: T, dtau_entry: T, c: T) -> T {
    let cap = ab_entry * c;
    if ab > cap && ab > T::zero() {
        dtau_entry * (cap / ab)
    } else {
        T::infinity()
    }
}

/// **A.** Was the step one the geometry could afford?
///
/// Rejects when either relative position moved further than `f * d_min` **as it was before the
/// step** -- the pre-step `d_min` is the one that bounds what a step may do, and using the
/// post-step value would let a step that destroyed the configuration justify itself. A non-finite
/// candidate is a rejection, not an acceptance: that is the case the retry exists for.
fn step_accepted<T: Real>(
    sys: &AzSystem<T>,
    before: &super::state::AzState<T>,
    after: &super::state::AzState<T>,
    f: T,
) -> bool {
    let (a1, a2, _, _) = sys.phys_from_state(before);
    let (b1, b2, _, _) = sys.phys_from_state(after);
    if !b1.is_finite() || !b2.is_finite() {
        return false;
    }
    let d_min = a1.norm().min(a2.norm()).min((a2 - a1).norm());
    let moved = (b1 - a1).norm().max((b2 - a2).norm());
    moved <= f * d_min
}

/// The step size for the *next* RK4 step, given the interval's entry sizing and the current
/// regularised state.
///
/// `ab_raw` is `A*B` **before** the `T::TINY` floor, which is what the blow-up is a function of
/// and what the `TINY` report reads.
#[inline]
fn dtau_for_step<T: Real>(
    mode: DtauMode,
    blend: StepBlend,
    blend_p: T,
    clamp_final: bool,
    eta: T,
    dt_left: T,
    dtau_entry: T,
    // The per-step limit, already in `dtau` units, or `+inf` when no limit is in force.
    // **Applied before the landing clamp**, so the clamp stays the last word and the final step
    // still lands on the boundary rather than stopping short of it.
    extra_limit: T,
    s: &super::state::AzState<T>,
) -> StepSize<T> {
    let a = s.a();
    let b = s.b();
    let ab_raw = a * b;
    let floored = a < T::TINY || b < T::TINY;
    let ab = a.max(T::TINY) * b.max(T::TINY);
    let mut capped = false;
    // **Every competing constraint is gathered, then combined once.** Written this way rather
    // than as a chain of `.min()` calls so the blend has something to blend: a chain hard-codes
    // the hard minimum into the control flow, which is the crease this exists to remove.
    let mut cons = [T::infinity(); 3];
    match mode {
        DtauMode::FixedPerInterval => cons[0] = dtau_entry,
        // `rem`, not `dt_left`: geometric decay, kept so the measurement can show it.
        DtauMode::PerStepRemaining => {
            let rem = (dt_left - s.t).max(T::zero());
            cons[0] = eta * rem / ab;
        }
        DtauMode::PerStepInterval => {
            cons[0] = eta * dt_left / ab;
            cons[1] = dtau_entry;
            capped = cons[0] > cons[1];
        }
    }
    cons[2] = extra_limit;
    let dtau = blend_dtau(blend, blend_p, &cons);
    // Land ON the boundary rather than past it. `dt = A*B*dtau` to leading order, so
    // `(dt_left - s.t)/ab` is the fictitious time still owed; `.min` makes this a one-sided
    // reduction of the last step only and leaves every interior step untouched. **`ab` is the
    // same floored product the mode used** -- recomputing it here would let the clamp and the
    // step disagree. The landing is exact only to the accuracy with which `A*B` predicts the time
    // increment over the step -- first order -- so the residual is `O(h^2)` per boundary and the
    // measured global order is **3.06 under `FixedPerInterval` and 2.08 under `PerStepInterval`**,
    // against 1.13 and 1.06 without. Not the stepper's own four, and stated rather than assumed.
    // **The landing clamp stays a HARD `min`, deliberately.** It is not a step-size preference
    // to be traded against the others -- it is an exact endpoint condition, and it is worth
    // 1.06 -> 2.08 in measured convergence order. Softening it would leave every final step
    // short of the boundary, costing the exactness for a crease that occurs once per interval at
    // essentially the same place for neighbouring pixels, rather than on the constraint-switching
    // surfaces that sweep across the frame.
    let dtau = if clamp_final {
        dtau.min((dt_left - s.t).max(T::zero()) / ab)
    } else {
        dtau
    };
    StepSize { dtau, ab_raw, floored, capped }
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
        // No step limit in the replay: it re-walks an interval the march already took, and a
        // different step rule would return a `t_end` from a different trajectory.
        let h = dtau_for_step(mode, StepBlend::Min, T::lit(4.0), clamp_final, eta, dt_left, dtau, T::infinity(), &s).dtau;
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
    // **D**, and the only limit that acts here rather than per step. Scaling `eta` buys
    // resolution uniformly and pays for it uniformly -- which is the property that makes it the
    // control the other three have to beat, not a defect of it.
    let eta = if opts.step_limit == StepLimit::Global {
        eta * T::lit(opts.step_limit_f)
    } else {
        eta
    };
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
    let mut ref_tie = Vec::with_capacity(n_sync);
    let mut boundary_shapes: Vec<[T; 3]> =
        Vec::with_capacity(if opts.keep_boundary_shapes { n_sync } else { 0 });
    let mut drift_hist: Vec<T> = Vec::with_capacity(if opts.keep_drift_hist { n_sync } else { 0 });
    let mut total_steps = 0usize;
    let mut finite = true;
    let mut budget_exhausted = false;
    // Recording only. `dt_max` is the largest physical step any interval took; `n_cap_hits`
    // counts steps at `PerStepInterval`'s entry-sizing cap. Neither changes a trajectory.
    let mut dt_max = T::zero();
    let mut n_cap_hits = 0u32;
    let mut n_overshoot = 0u32;
    let mut n_retry = 0u32;
    let mut retry_exhausted = false;
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
            None => {
                let want = choose_reference(&cart.r);
                // Hysteresis: hold the current reference until a rival beats **its** opposite
                // side by `1 + eps`. `THIRD[k] = 2 - k`, so the pair opposite reference `p` is
                // index `2 - p`. At `eps = 0` this is the plain argmax, bitwise.
                match prev_ref {
                    Some(p) if opts.ref_hysteresis > 0.0 && p != want => {
                        let d = crate::physics::newton::pair_dists(&cart.r);
                        let cur = d[2 - p];
                        let rival = d[0].max(d[1]).max(d[2]);
                        if rival > cur * (T::one() + T::lit(opts.ref_hysteresis)) {
                            want
                        } else {
                            p
                        }
                    }
                    _ => want,
                }
            }
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

            // The per-step limit, in `dtau` units and `+inf` when none is in force. `Global`
            // does not appear here: it scales `eta` at entry, which is what makes it the control
            // that pays everywhere rather than where the geometry asks.
            let extra = match opts.step_limit {
                StepLimit::Predictive => predictive_dtau(&sys, &s, T::lit(opts.step_limit_f)),
                StepLimit::AbGrowth => ab_growth_dtau(
                    a0 * b0,
                    s.a().max(T::TINY) * s.b().max(T::TINY),
                    dtau,
                    T::lit(opts.step_limit_f),
                ),
                StepLimit::None | StepLimit::Reject | StepLimit::Global => T::infinity(),
            };
            let ss =
                dtau_for_step(
                opts.dtau_mode, opts.step_blend, T::lit(opts.blend_p),
                opts.clamp_final_step, eta, dt_left, dtau, extra, &s,
            );
            if ss.ab_raw.is_finite() && ss.ab_raw < ab_min {
                ab_min = ss.ab_raw;
            }
            ab_floored |= ss.floored;
            n_cap_hits += ss.capped as u32;
            // The PHYSICAL increment, taken as a difference rather than as `A*B*dtau`: the
            // latter is the first-order predictor, and the whole question is how far the step
            // actually went. `s.t` is the interval-local clock, so this is a clean difference.
            let t_before = s.t;
            // **A** takes the step, tests it, and restores on failure. `AzState` is `Copy` and
            // nine numbers, so the save is free. Every attempt counts against the step budget --
            // a retry is real work and hiding it would make the mode look cheaper than it is.
            let s_save = s;
            let mut h = ss.dtau;
            let mut tries = 0u32;
            s = loop {
                let cand = rk4::step(&sys, &s_save, e, h);
                steps += 1;
                if opts.step_limit != StepLimit::Reject
                    || step_accepted(&sys, &s_save, &cand, T::lit(opts.step_limit_f))
                {
                    break cand;
                }
                tries += 1;
                n_retry += 1;
                if tries >= MAX_RETRIES {
                    // **Undetermined, not discarded.** Counted apart from `budget_exhausted` so
                    // one failure swapped for another is visible rather than absorbed.
                    retry_exhausted = true;
                    break cand;
                }
                h = h * T::lit(0.5);
            };
            let dt_took = s.t - t_before;
            if dt_took.is_finite() && dt_took > dt_max {
                dt_max = dt_took;
            }
            // **THE TRIPWIRE.** Not a remedy and not a condition to handle: a step that carries
            // the interval-local clock past twice the interval is a bug. A legitimate overshoot
            // under fixed `dtau` is at most one nominal step, ~1% of `dt_left`. Conditioned on
            // `clamp_final_step`, because with the clamp off overshoot is the *expected*
            // behaviour of a named measurement axis and an assert that fires on a deliberate
            // mode is a broken assert.
            if opts.clamp_final_step && s.t > dt_left * T::lit(2.0) {
                n_overshoot += 1;
                debug_assert!(
                    false,
                    "step overshot its interval: s.t = {} against dt_left = {}",
                    s.t, dt_left
                );
            }

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
            // Ascending, so `d[2]` is the longest and `d[1]` the second: 1 is a tie for the
            // LONGEST side, which is where the reference flips.
            ref_tie.push(d[1] / d[2].max(T::TINY));
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
        ref_tie,
        escape_flags,
        closure_hist,
        unbound_flags,
        t_end_at_entry,
        ab_min,
        ab_floored,
        dt_max,
        n_cap_hits,
        n_overshoot,
        n_retry,
        retry_exhausted,
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

#[cfg(test)]
mod step_control_tests {
    use super::*;
    use crate::Vec2;

    fn state(u1: f64, u2: f64) -> super::super::state::AzState<f64> {
        super::super::state::AzState {
            u1: Vec2::new(u1, 0.0),
            p1: Vec2::new(0.0, 0.0),
            u2: Vec2::new(u2, 0.0),
            p2: Vec2::new(0.0, 0.0),
            t: 0.0,
        }
    }

    /// **The floor fires, and the step is taken anyway.**
    ///
    /// This is the advance-anyway site the saturation question is about, and the test has to show
    /// both halves: that `floored` is set, *and* that the returned step is a usable finite number
    /// rather than a terminal. `ab_raw` is exactly zero here — `(1e-200)^2` underflows at f64, so
    /// the doubly-degenerate hole is open at **both** precisions and not only f32 — which is what
    /// makes the floor load-bearing rather than decorative: without it `dtau` would be `inf`.
    #[test]
    fn the_tiny_floor_fires_and_the_march_advances_anyway() {
        let s = state(1e-200, 1e-200);
        assert_eq!(s.a(), 0.0, "underflow is the premise of this test");
        let ss = dtau_for_step(DtauMode::PerStepInterval, StepBlend::Min, 4.0, false, 0.01, 1.0, 1e-3, f64::INFINITY, &s);
        assert!(ss.floored, "the floor must report itself");
        assert_eq!(ss.ab_raw, 0.0);
        assert!(ss.dtau.is_finite(), "the step is taken anyway -- that is the whole point");
    }

    /// The negative control. A guard that always fires passes as easily as one that never does.
    #[test]
    fn a_healthy_state_does_not_report_the_floor() {
        let ss = dtau_for_step(DtauMode::PerStepInterval, StepBlend::Min, 4.0, false, 0.01, 1.0, 1e-3, f64::INFINITY, &state(1.0, 1.0));
        assert!(!ss.floored);
        assert!(ss.ab_raw > 0.0);
    }

    /// `capped` is `PerStepInterval`'s entry-sizing bound and **only** that one.
    ///
    /// Under `FixedPerInterval` there is no cap to hit, so the same state must report `false` —
    /// otherwise the counter would be reading the mode rather than the step.
    #[test]
    fn the_cap_is_reported_only_where_it_exists() {
        // `A*B` has fallen below its entry value, so the mode wants a larger step than entry.
        let s = state(0.1, 0.1);
        let (eta, dt_left, dtau_entry) = (0.01, 1.0, 0.01f64);
        let want = eta * dt_left / (s.a() * s.b());
        assert!(want > dtau_entry, "the premise: the mode asks for more than entry sizing");

        let per = dtau_for_step(DtauMode::PerStepInterval, StepBlend::Min, 4.0, false, eta, dt_left, dtau_entry, f64::INFINITY, &s);
        assert!(per.capped);
        assert_eq!(per.dtau, dtau_entry, "capped means held AT entry, not near it");

        let fixed = dtau_for_step(DtauMode::FixedPerInterval, StepBlend::Min, 4.0, false, eta, dt_left, dtau_entry, f64::INFINITY, &s);
        assert!(!fixed.capped, "there is no cap under FixedPerInterval to hit");
    }

    /// And it does **not** fire when the mode gets the step it asked for.
    #[test]
    fn the_cap_does_not_fire_when_the_step_is_granted() {
        let s = state(3.0, 3.0);
        let ss = dtau_for_step(DtauMode::PerStepInterval, StepBlend::Min, 4.0, false, 0.01, 1.0, 1.0, f64::INFINITY, &s);
        assert!(!ss.capped);
    }

    /// A three-body system in the regularised coordinates, built from a Cartesian one, so the
    /// limit tests act on a state the integrator could actually be in.
    fn sys_and_state(sep: f64, speed: f64) -> (AzSystem<f64>, super::super::state::AzState<f64>) {
        let m = [3.0, 4.0, 5.0];
        let sys = AzSystem::new(0, 1, 2, m);
        let cart = Cart {
            r: [Vec2::new(0.0, 0.0), Vec2::new(sep, 0.0), Vec2::new(0.0, 2.0)],
            v: [Vec2::new(0.0, 0.0), Vec2::new(-speed, 0.0), Vec2::new(0.0, 0.0)],
        };
        let (st, _) = sys.to_reg(&cart);
        (sys, st)
    }

    /// **B's bound is on the PHYSICAL step, and in `dtau` units it is nearly separation-blind.**
    ///
    /// This is the regularisation doing its job, not a defect, and the first cut of this test got
    /// it wrong: `dt = A*B*dtau` with `A = |R1|`, so `f*d_min/(|v|*A*B)` has the separation cancel
    /// out of the numerator and denominator together, and a tight pair and a wide one return the
    /// **same** `dtau` to sixteen digits. Asserting "tighter is bounded harder" in `dtau` failed
    /// for that reason. It is worth knowing before the measurement: in the one unit AZ does not
    /// already adapt, B may be adding little — which is what the `error_ratio`-against-cost curve
    /// is for.
    ///
    /// What the bound must do is scale as a crossing time in `|v_rel|`: halving the speed at fixed
    /// geometry doubles it. That arm is not circular and is the one kept.
    #[test]
    fn the_predictive_limit_is_a_crossing_time_in_the_speed() {
        let f = 0.1;
        let (sys, tight) = sys_and_state(0.01, 1.0);
        let (_, wide) = sys_and_state(1.0, 1.0);
        let (lt, lw) = (predictive_dtau(&sys, &tight, f), predictive_dtau(&sys, &wide, f));
        assert!(lt.is_finite() && lw.is_finite());
        assert!(
            (lt - lw).abs() <= 1e-12 * lw,
            "the separation is expected to cancel in dtau units: {lt:e} against {lw:e}"
        );

        let (sys2, half) = sys_and_state(0.01, 0.5);
        let ratio = predictive_dtau(&sys2, &half, f) / lt;
        assert!((ratio - 2.0).abs() < 0.15, "expected ~2x on halving |v|, got {ratio:.4}");

        // And it is linear in `f`, so the knob means what its name says.
        let double_f = predictive_dtau(&sys, &tight, 2.0 * f);
        assert!((double_f / lt - 2.0).abs() < 1e-12);
    }

    /// And a state at rest has no crossing time, so it takes **no** bound rather than a zero step.
    #[test]
    fn the_predictive_limit_is_absent_when_nothing_moves() {
        let (sys, still) = sys_and_state(1.0, 0.0);
        assert!(predictive_dtau(&sys, &still, 0.1).is_infinite());
    }

    /// **C engages once `A*B` has grown past the factor, and by exactly the right amount.**
    ///
    /// The exact-value arm is the one with teeth: a clamp that merely *reduced* the step would
    /// satisfy an inequality while holding `dt` at the wrong place.
    #[test]
    fn the_ab_growth_clamp_engages_at_the_factor_and_not_before() {
        let (entry, dtau) = (1.0f64, 0.5f64);
        assert!(ab_growth_dtau(entry, 1.5, dtau, 2.0).is_infinite());
        assert!(ab_growth_dtau(entry, 2.0, dtau, 2.0).is_infinite(), "at the factor, not past it");
        let l = ab_growth_dtau(entry, 8.0, dtau, 2.0);
        assert!((l - dtau * 0.25).abs() < 1e-15, "got {l:e}");
        assert!((8.0 * l - 2.0 * entry * dtau).abs() < 1e-15, "dt held at cap * dtau_entry");
    }

    /// **A rejects a step that moves a body further than its own closest approach, and accepts
    /// one that does not.** The accept arm is what says it is not rejecting everything.
    #[test]
    fn the_acceptance_test_separates_a_large_step_from_a_small_one() {
        let (sys, st) = sys_and_state(0.05, 1.0);
        let small = rk4::step(&sys, &st, 0.0, 1e-6);
        let huge = rk4::step(&sys, &st, 0.0, 1e3);
        assert!(step_accepted(&sys, &st, &small, 0.25), "a tiny step must be accepted");
        assert!(!step_accepted(&sys, &st, &huge, 0.25), "a step of 1e3 must not be");
    }

    /// A non-finite candidate is a **rejection**. `NaN <= x` is `false`, so this could hold by
    /// accident; asserted so it cannot silently become an acceptance if the test is rewritten.
    #[test]
    fn a_non_finite_candidate_is_rejected() {
        let (sys, st) = sys_and_state(0.05, 1.0);
        let mut bad = st;
        bad.u1 = Vec2::new(f64::NAN, 0.0);
        assert!(!step_accepted(&sys, &st, &bad, 1e30));
    }

    /// The extra limit reaches the step, and `+inf` leaves it **bitwise** alone.
    ///
    /// The second arm is the one that matters: `StepLimit::None` passes `+inf`, and every
    /// committed number in the corpus was taken under it.
    #[test]
    fn the_extra_limit_binds_and_infinity_is_bitwise_inert() {
        let s = state(1.0, 1.0);
        let free =
            dtau_for_step(DtauMode::FixedPerInterval, StepBlend::Min, 4.0, false, 0.01, 1.0, 1e-3, f64::INFINITY, &s);
        let bound = dtau_for_step(DtauMode::FixedPerInterval, StepBlend::Min, 4.0, false, 0.01, 1.0, 1e-3, 1e-5, &s);
        assert_eq!(free.dtau.to_bits(), 1e-3f64.to_bits());
        assert_eq!(bound.dtau, 1e-5);
    }
}

/// What [`integrate_softref`] returns. Deliberately thin: this is a **diagnostic** integrator for
/// one question — does smoothing the reference-body choice smooth the drift field — and it does
/// not carry events, escape or outcome classification.
#[derive(Clone, Copy, Debug)]
pub struct SoftRefOut<T> {
    pub state: Cart<T>,
    pub drift: T,
    /// Arms integrated, summed over boundaries. `n_sync` means it never blended (pure argmax);
    /// `2*n_sync` means it blended two arms at every boundary. **This is the cost column.**
    pub arms: u64,
    pub finite: bool,
}

/// One sync interval under one reference body, returning the endpoint **and the interval-local
/// physical time it actually reached**.
///
/// The achieved time is not `dt_left`. The landing clamp sizes the final step from the
/// *instantaneous* `A*B`, a first-order predictor of the time increment, so the landing residual
/// is `O(h^2)` — which is exactly why the measured convergence order is 2.08 and not 3.06. Two
/// charts have different `A*B` and therefore land at *different* physical times, and comparing
/// their endpoints without accounting for that would manufacture a difference out of the time
/// transformations themselves. Callers must have the number to correct or to report.
///
/// `None` on a non-finite state or an exhausted budget — those are outcomes the caller must see,
/// not values to substitute for.
///
/// `StepLimit::Reject` is **not** supported here and is treated as no limit: a retry loop makes
/// the arm's cost depend on its own rejections, which would confound the arm-count column that
/// this experiment is measured on. Stated rather than silently mapped.
fn march_interval<T: Real>(
    cart: &Cart<T>,
    m: &[T; 3],
    a: usize,
    dt_left: T,
    eta: T,
    max_steps: usize,
    opts: &AzOpts<T>,
) -> Option<(Cart<T>, T)> {
    let (ab, bb, cb) = triple(a);
    let sys = if opts.lc_stable {
        AzSystem::new(ab, bb, cb, *m)
    } else {
        AzSystem::new(ab, bb, cb, *m).with_reference_lc()
    };
    let (mut s, e) = sys.to_reg(cart);
    let a0 = s.a().max(T::TINY);
    let b0 = s.b().max(T::TINY);
    let dtau = eta * dt_left / (a0 * b0);
    let mut steps = 0usize;
    loop {
        if !s.is_finite() {
            return None;
        }
        if s.t >= dt_left - land_tol::<T>(opts.clamp_final_step, dt_left) {
            break;
        }
        if steps >= max_steps {
            return None;
        }
        let extra = match opts.step_limit {
            StepLimit::Predictive => predictive_dtau(&sys, &s, T::lit(opts.step_limit_f)),
            StepLimit::AbGrowth => ab_growth_dtau(
                a0 * b0,
                s.a().max(T::TINY) * s.b().max(T::TINY),
                dtau,
                T::lit(opts.step_limit_f),
            ),
            StepLimit::None | StepLimit::Reject | StepLimit::Global => T::infinity(),
        };
        let ss = dtau_for_step(
            opts.dtau_mode, opts.step_blend, T::lit(opts.blend_p), opts.clamp_final_step,
            eta, dt_left, dtau, extra, &s,
        );
        s = rk4::step(&sys, &s, e, ss.dtau);
        steps += 1;
    }
    Some((sys.to_cartesian(&s), s.t))
}

/// **A softmax over the reference-body choice, instead of an argmax.**
///
/// # Is this even well posed?
///
/// You cannot be in two regularised charts at once, so the blend cannot happen inside a step.
/// But at every sync boundary the state is **Cartesian and chart-free**, and the reference choice
/// governs only how the *next* interval is integrated. All three choices approximate the same
/// true trajectory, so a convex combination of their endpoints is another approximation to it,
/// and it converges to the same limit as `dtau -> 0`. The blend is smooth in the geometry because
/// the weights are — which is the whole point: `choose_reference` is a bare `argmax` with no
/// hysteresis, and its cell boundaries in initial-condition space are measurably where the drift
/// field's edges are.
///
/// # The weights
///
/// `w_k ∝ exp((d_k - d_max) / (temp * d_max))` over the three pair separations, with the
/// reference body `THIRD[k]` as usual. **Relative, not absolute**: `d` carries units and this
/// project quotients out overall scale, so an absolute temperature would mean different things at
/// different hyperradii. `temp <= 0` is exact `argmax`, including its first-maximum tie-break, so
/// the zero-temperature limit reproduces the shipped path bitwise rather than approximately.
///
/// Arms below `W_TOL` are pruned and the rest renormalised, so away from a tie this costs exactly
/// one arm and the expense is confined to a thin shell around the tie surface.
///
/// # What it does not do
///
/// No events, no escape, no outcome. A blended state has no single terminal classification, and
/// inventing one would put a discontinuity straight back in via a different door. This measures
/// the **drift field** only, which is the field the edges were seen in.
pub fn integrate_softref<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &AzOpts<T>,
    temp: T,
) -> SoftRefOut<T> {
    /// Weight below which an arm is not worth integrating.
    const W_TOL: f64 = 1e-3;

    let e0 = energy::energy(&s0.r, &s0.v, m, T::zero());
    let mut cart = s0;
    let mut arms = 0u64;
    let mut finite = true;

    for k in 0..n_sync {
        let t = t_max * T::lit(k as f64 / n_sync as f64);
        let t_target = t_max * T::lit((k + 1) as f64 / n_sync as f64);
        let dt_left = t_target - t;

        let d = crate::physics::newton::pair_dists(&cart.r);
        let d_max = d[0].max(d[1]).max(d[2]);

        let mut w = [T::zero(); 3];
        if temp > T::zero() && d_max > T::zero() && d_max.is_finite() {
            for j in 0..3 {
                w[j] = ((d[j] - d_max) / (temp * d_max)).exp();
            }
        } else {
            // Exact argmax, first maximum — `choose_reference`'s own convention, so `temp = 0`
            // reproduces the shipped path rather than merely resembling it.
            w[crate::integrate::az::reference_body::longest_index(&d)] = T::one();
        }
        let total: T = w.iter().copied().sum();
        if !(total > T::zero()) || !total.is_finite() {
            finite = false;
            break;
        }
        for x in w.iter_mut() {
            *x = *x / total;
        }
        // Prune, then renormalise over the survivors. Away from a tie exactly one arm survives.
        let live: T = w.iter().copied().filter(|x| *x >= T::lit(W_TOL)).sum();

        let mut blended = Cart { r: [crate::Vec2::new(T::zero(), T::zero()); 3], v: [crate::Vec2::new(T::zero(), T::zero()); 3] };
        let mut any = false;
        for j in 0..3 {
            if w[j] < T::lit(W_TOL) {
                continue;
            }
            let a = crate::physics::THIRD[j];
            match march_interval(&cart, m, a, dt_left, eta, max_steps, opts) {
                Some((c, _)) => {
                    arms += 1;
                    let wj = w[j] / live;
                    for i in 0..3 {
                        blended.r[i] = blended.r[i] + c.r[i] * wj;
                        blended.v[i] = blended.v[i] + c.v[i] * wj;
                    }
                    any = true;
                }
                None => {
                    // An arm that failed is an outcome, not a value to substitute for. The
                    // blend is abandoned rather than silently re-weighted onto the survivors:
                    // re-weighting would hide a failed integration inside a plausible state.
                    finite = false;
                    any = false;
                    break;
                }
            }
        }
        if !any {
            finite = false;
            break;
        }
        cart = blended;
    }

    let e1 = energy::energy(&cart.r, &cart.v, m, T::zero());
    SoftRefOut {
        state: cart,
        drift: ((e1 - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs(),
        arms,
        finite,
    }
}

/// What [`branch_jump`] measured at one chart-switching event.
///
/// **Position and velocity are reported separately, never as one norm.** They are dimensionally
/// different, so a combined Euclidean norm is arbitrary unless phase space has been explicitly
/// non-dimensionalised, which it has not been here.
#[derive(Clone, Copy, Debug)]
pub struct BranchJump<T> {
    /// The boundary was reached and every arm integrated. `false` means no measurement, not a
    /// measurement of zero.
    pub ok: bool,
    /// `|| r_win - r_alt ||` at a **common physical time** — the jump crossing the argmax surface.
    pub dr_chart: T,
    /// `|| v_win - v_alt ||` at the same common time.
    pub dv_chart: T,
    /// The same interval, same chart, `eta` against `eta/2`: the ORDINARY local step error, and
    /// the normaliser. Without it an absolute jump is unreadable — a large number on a violent
    /// interval and a small one on a quiet interval say the same thing.
    pub dr_step: T,
    pub dv_step: T,
    /// **The confound, measured rather than assumed.** The two charts land at different physical
    /// times because the landing residual is `O(h^2)` and `A*B` differs between them. This is
    /// `|t_win - t_alt|`.
    pub dt_mismatch: T,
    /// `|| r_win - r_alt ||` **without** the common-time correction. If this differs materially
    /// from `dr_chart`, the time mismatch was doing the work and the raw comparison would have
    /// manufactured a branch discrepancy out of the time transformations themselves.
    pub dr_raw: T,
    /// Largest speed at the endpoint — with `dt_mismatch`, the displacement the mismatch alone
    /// could explain.
    pub speed: T,
    /// Second-longest over longest at the entry state: distance to the selector's surface.
    pub ref_tie: T,
    /// The argmax winner and the runner-up.
    pub ref_win: usize,
    pub ref_alt: usize,
}

/// **How large a perturbation does crossing the argmax surface actually inject?**
///
/// March normally to sync boundary `n`, then from that one Cartesian state integrate the *same*
/// interval under both the chosen reference and the runner-up, and measure the difference. The
/// normaliser is the ordinary local step error — the same interval, same chart, at `eta` against
/// `eta/2` — so the answer is a ratio: *the chart jump is N times the step error*.
///
/// This is the amplitude that `ref_tie` cannot give. Together they are the mechanism:
/// `d_i - d_j -> 0` locates the surface, `delta_chart / delta_step` says the jump is real and
/// how big, and chaotic amplification over the remaining `t_max - t_n` does the rest.
///
/// A ratio near 1 would mean crossing the surface costs no more than an ordinary step, and the
/// selector is then a symptom rather than a seed.
#[allow(clippy::too_many_arguments)]
pub fn branch_jump<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &AzOpts<T>,
    n: usize,
) -> BranchJump<T> {
    let none = BranchJump {
        ok: false,
        dr_chart: T::nan(),
        dv_chart: T::nan(),
        dr_step: T::nan(),
        dv_step: T::nan(),
        dt_mismatch: T::nan(),
        dr_raw: T::nan(),
        speed: T::nan(),
        ref_tie: T::nan(),
        ref_win: 0,
        ref_alt: 0,
    };
    if n >= n_sync {
        return none;
    }
    let mut cart = s0;
    for k in 0..n {
        let t = t_max * T::lit(k as f64 / n_sync as f64);
        let dt_left = t_max * T::lit((k + 1) as f64 / n_sync as f64) - t;
        let a = super::reference_body::choose_reference(&cart.r);
        match march_interval(&cart, m, a, dt_left, eta, max_steps, opts) {
            Some((c, _)) => cart = c,
            None => return none,
        }
    }

    // The two candidate charts at the entry state: the longest side and the runner-up.
    let d = crate::physics::newton::pair_dists(&cart.r);
    let mut idx = [0usize, 1, 2];
    idx.sort_by(|&a, &b| d[b].partial_cmp(&d[a]).unwrap_or(std::cmp::Ordering::Equal));
    let (win, alt) = (crate::physics::THIRD[idx[0]], crate::physics::THIRD[idx[1]]);
    let ref_tie = d[idx[1]] / d[idx[0]].max(T::TINY);

    let t = t_max * T::lit(n as f64 / n_sync as f64);
    let dt_left = t_max * T::lit((n + 1) as f64 / n_sync as f64) - t;

    // Separate norms, per the note on `BranchJump`.
    let dr = |a: &Cart<T>, b: &Cart<T>| -> T {
        (0..3).map(|i| (a.r[i] - b.r[i]).norm_sq()).fold(T::zero(), |x, y| x + y).sqrt()
    };
    let dv = |a: &Cart<T>, b: &Cart<T>| -> T {
        (0..3).map(|i| (a.v[i] - b.v[i]).norm_sq()).fold(T::zero(), |x, y| x + y).sqrt()
    };
    // **Bring an arm to the common target time.** The landing residual is `O(h^2)` and differs
    // between charts, so the arms stop at different physical times; a first-order drift by the
    // endpoint velocity removes that, and its own error is `O(dt^2 * accel)` with `dt ~ h^2`,
    // i.e. `O(h^4)` — far below what is being measured. The uncorrected value is returned too,
    // so the size of the correction is visible rather than trusted.
    let to_time = |c: &Cart<T>, reached: T, target: T| -> Cart<T> {
        let d = target - reached;
        let mut o = *c;
        for i in 0..3 {
            o.r[i] = o.r[i] + o.v[i] * d;
        }
        o
    };

    let two = T::lit(2.0);
    match (
        march_interval(&cart, m, win, dt_left, eta, max_steps, opts),
        march_interval(&cart, m, alt, dt_left, eta, max_steps, opts),
        march_interval(&cart, m, win, dt_left, eta / two, max_steps * 2, opts),
    ) {
        (Some((a, ta)), Some((b, tb)), Some((a2, ta2))) => {
            let (ca, cb) = (to_time(&a, ta, dt_left), to_time(&b, tb, dt_left));
            let ca2 = to_time(&a2, ta2, dt_left);
            let speed = (0..3)
                .map(|i| a.v[i].norm_sq().sqrt())
                .fold(T::zero(), |x, y| if y > x { y } else { x });
            BranchJump {
                ok: true,
                dr_chart: dr(&ca, &cb),
                dv_chart: dv(&ca, &cb),
                // The step-error normaliser gets the SAME common-time treatment, or the two
                // would be measured under different conventions and the ratio would be a
                // comparison of methods rather than of magnitudes.
                dr_step: dr(&ca, &ca2),
                dv_step: dv(&ca, &ca2),
                dt_mismatch: (ta - tb).abs(),
                dr_raw: dr(&a, &b),
                speed,
                ref_tie,
                ref_win: win,
                ref_alt: alt,
            }
        }
        _ => none,
    }
}
