//! Per-pixel evaluation: build the copies, integrate them, reduce.

use crate::decode::Path;
use crate::grid::Slice;
use crate::integrate::az::{self, AzOpts, RefPolicy};
use crate::outcome::{self, Outcome, State};
use crate::physics::{energy, shape};
use crate::Real;

use super::jitter::Scheme;
use super::{jitter, stats};

#[derive(Clone, Copy, Debug)]
pub struct EnsembleCfg {
    /// `E` extra copies; the pixel always carries `E + 1`.
    pub n_extra: usize,
    pub jitter_frac: f64,
    /// How copy offsets are chosen. Default is the spec's fixed Halton (2,3) prefix; `Pcg`
    /// reproduces the reference's per-pixel stream and every result measured before the switch.
    pub jitter_scheme: Scheme,
    pub seed: u64,
    pub t_max: f64,
    pub n_sync: usize,
    /// Escape distance gate as a **fraction of the initial hyperradius** `R`, fixed at `t = 0`.
    ///
    /// **This is the condition the port did not have.** The numpy reference tests unbound and
    /// receding; the GLSL also tests `dist > r_esc`, and without it an escape is declared at
    /// any distance including mid-encounter. See
    /// [`crate::integrate::az::AzOpts::r_esc_frac` ] for the derivation of the value and
    /// `examples/escape_gate.rs` for the sweep.
    ///
    /// **Default `5.0`, and every committed result predates it.** The GLSL's saved configs use
    /// 5 and 12; on the latent charts `R = 1` identically so those are already fractions.
    /// `--r-esc 0` restores the numpy behaviour every dump in `results/` was made under.
    pub r_esc_frac: f64,
    /// Test all three bodies for escape, as the GLSL does, rather than only the body outside
    /// the tightest pair, as the numpy reference does.
    pub escape_all_bodies: bool,
    /// Escape-test stride inside the RK4 loop; `0` is the reference's boundary-only cadence.
    /// See [`crate::integrate::az::AzOpts::escape_every`] -- this is the knob that decides
    /// whether `t_end` carries 32 distinct values across a chart or RK4-step resolution.
    pub escape_every: usize,
    /// Require an in-loop escape to still hold at the next sync boundary. See
    /// [`crate::integrate::az::AzOpts::escape_confirm`] — without it the in-loop test latches
    /// transients, measured at **895 of 895** in `deep interior`.
    pub escape_confirm: bool,
    pub eta: f64,
    pub max_steps: usize,
    pub ref_policy: RefPolicy,
    /// Conditioned inverse LC branch. Default true; false reproduces the reference's
    /// original branch, for measuring what the conditioning is worth.
    pub lc_stable: bool,
    /// `r_coll` as a **fraction of the initial hyperradius** `R`, fixed at `t = 0`. Never an
    /// absolute length and never co-moving (BRIEF §2.5).
    ///
    /// The default has no reference — `r_coll` appears nowhere in the numpy tree — so it is
    /// reported as a measurement instead: see `examples/r_coll_sweep.rs` for how the outcome
    /// fractions move across `{1e-4, 1e-3, 1e-2}`.
    pub r_coll_frac: f64,
    /// Stop at the first terminating event (BRIEF §2.4). Off keeps every copy integrating to
    /// `t_max`, which is the reference's behaviour.
    pub stop_on_event: bool,

    /// Re-integrate flagged pixels at finer `eta` after the first pass. **BRIEF §2.5.**
    ///
    /// `eta = 1e-2` is not sufficient above ~64x64: at 128x128 near-field, 7 pixels of 16384
    /// carry `|dE/E| > 1`, worst `1.49e4`. The failure is a **cliff, not a slope** — the same
    /// pixels give `1.49e4 -> 5.96e-9` for a 3.3x change in `eta` — so an adaptive-`eta`
    /// policy that interpolates between regions has nothing to interpolate along. Flagging and
    /// re-integrating does work, because `error_ratio` flagged 7 of 7 and one refinement step
    /// was enough.
    ///
    /// This is one extra pass over a small flagged subset, not a scheduler. Measured cost on
    /// the 128x128 near-field grid: ~1.6% of pixels flagged.
    pub refine_flagged: bool,
    /// `error_ratio` above which a pixel is re-integrated. Set from the measurement: at
    /// `> 10` it caught 7 of the 7 pixels with `|dE/E| > 1` and flagged ~1.6% of the grid.
    pub refine_threshold: f64,
    /// Factor applied to `eta` on each refinement pass. 1/4 sits past the measured cliff,
    /// which closed between `1e-2` and `3e-3`.
    pub refine_eta_factor: f64,
    /// Which arithmetic forms each copy's initial condition. Default [`Path::DirectF64`] —
    /// every result before the vertical slice. The other paths exist for the deep-zoom
    /// measurement, where a **collapsed** decode is the failure to watch for: identical
    /// footprints give a spread of exactly zero, which the criterion reads as
    /// "perfectly resolved" rather than as "no data".
    pub decode_path: Path,
    /// Keep each copy's packed outcome, for the SSAA resolve. Off by default: it is only
    /// wanted at render time, and it makes `PixelOut` allocate.
    pub keep_copy_outcomes: bool,
    /// Keep each copy's `shape_vec`, for the criterion's **matched-count** comparison.
    ///
    /// Off by default, for the same reason as [`Self::keep_copy_outcomes`]: it allocates
    /// `E+1` triples per footprint. It exists because comparing an `E+1`-sample within-cell
    /// spread against an `N^2`-sample between-quad spread conflates *scale* with *sample
    /// count*, and a spread estimator's expectation depends on the count — measured, `E+1 = 2`
    /// reports 0.539 of `E+1 = 32`'s spread in near-field and 0.131 in `far`. `within_pooled`
    /// needs every copy, not just the nominal, to hold the count fixed while the extent moves.
    pub keep_copy_shapes: bool,
    /// Record each copy's `shape_vec` at every sync boundary, for the §5 temporal
    /// accumulators. Off by default; reduced and dropped inside one footprint's evaluation.
    pub keep_boundary_shapes: bool,
    /// Also run the nominal copy through the Benettin/diffusion accumulators, for the
    /// production colouring's lightness field. `None` skips it entirely.
    ///
    /// **Off by default and expensive**: it is a second, fixed-step, unregularised march at
    /// `dt` alongside the AZ one, so at `t = 13, dt = 1e-4` it is 130,000 steps against AZ's
    /// few thousand. Only the nominal copy is run — FTLE is a per-point scalar, not an
    /// ensemble statistic — so the cost is per footprint, not per copy.
    ///
    /// It is the **unregularised** integrator, so the result is only trustworthy away from a
    /// close approach. `d_min` from the AZ march is the column to read beside it.
    pub ftle: Option<crate::physics::ftle::FtleOpts>,
    /// Fixed step for the FTLE march. The reference's `1e-4`.
    pub ftle_dt: f64,
    /// Maximum refinement passes. **Bounded on purpose, and not a scheduler**: each pass is
    /// one extra evaluation of a shrinking flagged subset, with no tree and no state carried
    /// between pixels.
    ///
    /// One pass suffices in the ordinary regions — near-field 256x256 goes from `1.38e7` to
    /// `3.12e-4` — but **not everywhere**. `deep interior` needs more: one pass takes it from
    /// `1.10e12` to `1.99e1`, still far above any usable bound, because 14% of that region is
    /// flagged and its close approaches are much deeper. A pixel still flagged after the last
    /// pass keeps its `error_ratio`, so an unrepaired pixel is reported, never silently
    /// accepted.
    pub refine_max_passes: u8,
}

impl Default for EnsembleCfg {
    fn default() -> Self {
        Self {
            n_extra: 7, // E + 1 = 8, per BRIEF §3
            jitter_frac: 0.5,
            jitter_scheme: Scheme::Halton,
            seed: 0,
            t_max: 13.0,
            n_sync: 32,
            escape_every: 0,
            escape_confirm: true,
            r_esc_frac: 5.0,
            escape_all_bodies: true,
            eta: 0.01,
            max_steps: 30_000,
            ref_policy: RefPolicy::PerCopy,
            lc_stable: true,
            r_coll_frac: 1e-3,
            stop_on_event: true,
            refine_flagged: true,
            refine_threshold: 10.0,
            refine_eta_factor: 0.25,
            refine_max_passes: 3,
            keep_copy_outcomes: false,
            keep_copy_shapes: false,
            keep_boundary_shapes: false,
            ftle: None,
            ftle_dt: 1e-4,
            decode_path: Path::DirectF64,
        }
    }
}

/// Everything BRIEF §4 asks for, plus the fields that make its confounds visible.
///
/// All stored as `f64` regardless of kernel precision, so an f32 run and an f64 run produce
/// directly comparable dumps.
#[derive(Clone, Debug, Default)]
pub struct PixelOut {
    /// BRIEF §2.4's packed encoding: `state` in the high 3 bits, `detail` in the low 2.
    pub outcome: u8,
    /// Unpacked for readers who would rather not shift. Same information.
    pub state: u8,
    pub detail: u8,
    /// `tb.classify`'s labelling — kept because it is the one classification with a
    /// reference, and `spread_event_legacy` is checkable against it.
    pub legacy_class: u8,
    pub binary_id: u8,
    pub t_end: f64,
    pub censored: bool,

    /// `min(|R1|,|R2|)` — the reference's blind spot, kept for cross-check comparability.
    pub d_min_ref: f64,
    /// All three pairs, including the unregularised side. BRIEF §4's actual definition.
    pub d_min_true: f64,
    /// `d_min_ref - d_min_true`. Measures how well the reference-switching cadence tracks
    /// encounters, instead of arguing about it (NOTES §2.1).
    pub d_min_gap: f64,

    pub energy_drift_nominal: f64,
    pub energy_drift_max: f64,
    pub gamma_max: f64,

    pub shape_vec: [f64; 3],
    /// BRIEF §4: mean distance from the centroid, halved. The spec.
    pub spread_shape: f64,
    /// `refine_test.svar`: `1 - |mean|`. A *different* statistic, and the one with a
    /// reference. Both are dumped so the discrepancy is measurable rather than a silent
    /// choice.
    pub svar: f64,
    /// Disagreement over the **event class** — the currently tightest pair at the final sync
    /// boundary, joined with the terminal outcome for copies that have terminated.
    ///
    /// **Not the terminal `(state, detail)`.** That quantity is terminal-grain and inverts
    /// under lockstep: early in the march nothing has terminated, so every copy agrees and
    /// the field reports maximum confidence where least is known. Measured on near-field
    /// 32x32: at `t_max = 8` the corrected field fires on 110 of 1024 pixels and the terminal
    /// one on **none**; at `t_max = 13`, 35 against 22. See NOTES §2.8.
    pub spread_event: f64,
    /// Running **max** of the event-class spread over every boundary up to the playhead.
    ///
    /// Dumped because the playhead value is a snapshot and can *un*-fire: the tightest-pair
    /// identity fluctuates, so copies that disagreed at one boundary can agree again at the
    /// next. Measured, the playhead value fires on 110 of 1024 pixels at `t_max = 8` and only
    /// 35 at `t_max = 13` — non-monotone in the horizon, which is not what a confidence flag
    /// should do. This one is monotone. Which of the two `ensemble_spread` should use is a
    /// judgement, so both are dumped and the spec one is the default.
    pub spread_event_max: f64,
    /// The latching field: the running max over boundaries, **guarded by persistence**, joined
    /// with the playhead value.
    ///
    /// An unguarded running max is wrong for a *discrete* label. Measured, near-field 32x32 at
    /// `t = 13`: of 165 pixels that ever disagree, 130 disagree and then re-agree, and **129 of
    /// those 130 were at a near-tie** (second-tightest/tightest below 1.1, median 1.0030) at
    /// the boundary where they first disagreed. The copies disagreed about which pair is
    /// *tightest* without having diverged. An unguarded latch would light 79% of the firing
    /// pixels permanently for a labelling artefact.
    ///
    /// The tie ratio does not separate the populations cleanly enough to threshold on —
    /// genuine disagreements also sit near 1 (median 1.0797). **Persistence does**: artefacts
    /// last one boundary (median run 1, max 2), genuine divergence persists (median run 10).
    /// A run of [`LATCH_RUN`] admits 0 of 130 artefacts; joining with the playhead value picks
    /// up the genuine disagreements too recent to have persisted yet, which is censoring at
    /// the horizon rather than a miss.
    pub spread_event_latched: f64,
    /// The terminal-outcome version, kept so the correction stays a measured difference.
    pub spread_event_terminal: f64,
    /// Over `classify_legacy`. Dumped because it is the reference-checkable one.
    pub spread_event_legacy: f64,
    /// Minimum over copies of (second-tightest / tightest) pair separation, at the boundary
    /// where the copies' event classes **first disagree**. NaN if they never do.
    ///
    /// This is the measurement that decides whether `spread_event` may latch. Near 1 means a
    /// **near-tie**: the copies disagree about which pair is tightest without having diverged,
    /// so a running max would latch an artefact that never clears. Well above 1 means the
    /// disagreement is genuine divergence and latching is correct. See NOTES §5.
    pub tie_ratio_at_disagree: f64,
    /// Number of sync boundaries at which the copies' event classes disagree, and the longest
    /// **consecutive** run of them.
    ///
    /// A near-tie produces isolated single-boundary disagreements that clear immediately; a
    /// genuine divergence persists. This is the lever a latching field would need if it is not
    /// to latch an artefact, and it is dumped so the guard can be chosen from data.
    pub n_disagree: u16,
    pub longest_disagree_run: u16,
    /// First sync-boundary time at which the copies' event classes disagree; **NaN** if they
    /// never do — not `t_max`, which would be indistinguishable from disagreeing at the last
    /// boundary. This is the property the event class was chosen for — it fires while the
    /// march is still running rather than only at a terminal label.
    pub t_spread_event: f64,
    pub ensemble_spread: f64,

    /// Dumped separately from the ratio. `sigma_E(0)` is proportional to the jitter and so
    /// to cell width; as resolution rises it shrinks while integration error does not, and
    /// `error_ratio` inflates for a purely trivial reason. Burying that in a ratio would
    /// hide it (threatens BRIEF §8 experiments 1 and 3).
    pub sigma_e_0: f64,
    pub sigma_e_t: f64,
    /// Built on the **maximum deviation** from the median, not the MAD. See
    /// [`crate::ensemble::stats`] for why: MAD is robust to exactly the single wild copy this
    /// field exists to flag, at a measured damaged/healthy separation of 1.06 against 59.51.
    pub error_ratio: f64,
    /// The MAD-based ratio, BRIEF §4's original wording. Dumped so the change stays a
    /// measured difference rather than a silent replacement.
    pub error_ratio_mad: f64,

    pub switches: u32,
    /// Per-sync count of copies whose chosen reference body differs from the nominal copy's.
    /// Lets `error_ratio` be conditioned on reference disagreement, so §8 experiment 2 gets
    /// an attributed answer from one run rather than a difference of aggregates (NOTES §1).
    pub ref_disagree: u32,
    /// Copies whose `(state, detail)` differs from the nominal copy's. `spread_event`
    /// normalised away; this is the raw count.
    pub n_outcome_disagree: u8,
    /// Copies that went non-finite. **Never discarded** — a copy that could not be
    /// determined is a measurement outcome, and this records it explicitly rather than
    /// letting a NaN contaminate every aggregate it touches.
    pub n_nonfinite: u8,

    /// The `eta` this pixel's reported values were computed at. Equal to `cfg.eta` unless the
    /// pixel was flagged and re-integrated, in which case it is the refined one.
    pub eta_used: f64,
    /// `error_ratio` and `energy_drift_max` from the **first** pass, kept whether or not a
    /// refinement happened. A refinement that silently replaced the coarse value would hide
    /// exactly the pixels this mechanism exists to find.
    pub error_ratio_coarse: f64,
    pub energy_drift_max_coarse: f64,
    /// Whether the second pass ran on this pixel.
    pub refined: bool,

    /// Every copy's packed `(state, detail)`, in copy order. Empty unless
    /// [`EnsembleCfg::keep_copy_outcomes`] is set.
    ///
    /// **This is the SSAA input, and it is not `spread_event`.** The `E+1` copies serve two
    /// jobs that must not be confused: `spread_*` is a *disagreement* statistic and drives
    /// scheduling; resolve is an *average* and drives display. A pixel where the copies split
    /// 4/4 has a large spread and a blended colour, and neither number substitutes for the
    /// other.
    pub copy_outcomes: Vec<u8>,

    /// The **nominal copy's event class** at the final sync boundary: the identity of its
    /// currently-tightest pair, or `TERMINAL_TAG + terminal` once it has stopped.
    ///
    /// One byte, always retained, and it exists so a *between*-footprint event statistic can
    /// be built without reinstating a quantity the project has already rejected. The obvious
    /// implementation of `between_event` is `spread_event` over the `N^2` footprints'
    /// [`Self::outcome`] — and that is **the terminal outcome**, which is terminal-grain and
    /// inverts under lockstep: early in the march nothing has terminated, every footprint
    /// agrees, and the field reports maximum confidence at exactly the playhead where least is
    /// known. The event class is defined at every playhead. Building the between-footprint arm
    /// on `outcome` would be that regression at a new level, so the class is carried instead.
    pub event_class: u8,

    /// Integrator substeps summed over the `E+1` copies — the cost side of a cost-aware
    /// priority. `AzOut::steps` was computed on every march and never read; ranking by
    /// `spread / cost` only pays if the cost distribution is wide, and this is what says
    /// whether it is.
    pub total_substeps: u64,

    /// Every copy's `shape_vec`, in copy order. Empty unless
    /// [`EnsembleCfg::keep_copy_shapes`] is set.
    ///
    /// The input to `within_pooled` — the within-footprint arm evaluated at the *same* sample
    /// count as the between-footprint arm, so the two can be differenced without the
    /// small-sample bias standing in for a scale effect.
    pub copy_shapes: Vec<[f64; 3]>,

    // -----------------------------------------------------------------------------------
    // §5 — the temporal accumulators, shape arm.
    //
    // **The event arm already has all three**, contrary to the brief: `spread_event_max` is a
    // running max over boundaries, `t_spread_event` is a first-divergence time that is NaN
    // rather than `t_max` when it never fires, and `spread_event_latched` is the
    // persistence-guarded latch. What was missing is the CONTINUOUS arm, below.
    //
    // All three are `NaN` unless `EnsembleCfg::keep_boundary_shapes` is set — never 0, which
    // would read as "no divergence" on a quantity that was not measured.
    // -----------------------------------------------------------------------------------
    /// Running max of `spread_shape` over every boundary up to the playhead.
    ///
    /// Empirically justified rather than assumed: the shape spread was observed to **fall 6x**
    /// between `t = 6` and `t = 8` in one region, so an instantaneous read genuinely misses
    /// divergence that has already happened. Max-updated, never decayed.
    pub running_max_divergence: f64,
    /// OLS slope of `spread_shape` against boundary time — is divergence still growing at the
    /// playhead, or has it saturated? NaN below two boundaries.
    pub divergence_trend: f64,
    /// First boundary time at which `spread_shape` crosses [`DIVERGENCE_TRIGGER`] of its
    /// achievable maximum; **NaN if it never does** — not `t_max`, which would be
    /// indistinguishable from crossing at the last boundary.
    ///
    /// **The one signal that cannot saturate.** Every instantaneous spread saturates once the
    /// copies fill the accessible space, and then reports `lambda ~ 0` for the *most* chaotic
    /// regions — the inversion this project has now met three times. A crossing time cannot do
    /// that.
    pub first_divergence_t: f64,

    /// Benettin FTLE of the **nominal** copy, `NaN` unless `EnsembleCfg::ftle` is set.
    pub ftle: f64,
    /// Slope of `log(inertia)` against `t` for the nominal copy, `NaN` unless enabled.
    pub diffusion: f64,
    /// Renormalisations completed. **Assert this is nonzero before reading `ftle`**: without
    /// renormalisation the shadow saturates and the estimator reports `lambda ~ 0` for the most
    /// chaotic regions, which is the inversion rather than a small error.
    pub ftle_renorm: u64,
}

/// Fraction of `spread_shape`'s achievable maximum that counts as diverged.
///
/// `spread_shape` is a mean chord distance halved, so it is bounded by 1 and this is an
/// absolute fraction of a known ceiling rather than a tuned constant.
pub const DIVERGENCE_TRIGGER: f64 = 0.1;

/// Consecutive boundaries a disagreement must survive before the latch counts it.
///
/// Chosen from the measurement in `examples/latching_decision.rs`, not by eye: at 2 it admits
/// 1 of 130 artefacts, at 3 it admits none, and raising it further changes nothing. On this
/// slice, at this `n_sync`. It is a named constant so a future change to it shows in a diff.
pub const LATCH_RUN: u16 = 3;

/// One pass at the configured `eta`, then — if the pixel is flagged — one more at finer `eta`.
///
/// The refinement is deliberately **flag-driven, not region-driven**. See
/// [`EnsembleCfg::refine_flagged`]: the failure it addresses is a cliff, and a cliff gives an
/// adaptive-`eta`-by-region policy nothing to interpolate along.
pub fn evaluate<T: Real>(slice: &Slice, idx: usize, cfg: &EnsembleCfg) -> PixelOut {
    let coarse = evaluate_at::<T>(slice, idx, cfg, cfg.eta);
    if !cfg.refine_flagged || !(coarse.error_ratio > cfg.refine_threshold) {
        return coarse;
    }

    let mut out = coarse.clone();
    let mut eta = cfg.eta;
    for _ in 0..cfg.refine_max_passes {
        eta *= cfg.refine_eta_factor;
        out = evaluate_at::<T>(slice, idx, cfg, eta);
        out.eta_used = eta;
        out.refined = true;
        if !(out.error_ratio > cfg.refine_threshold) {
            break;
        }
    }
    // Carried forward, not overwritten: the coarse pair is what a run without refinement would
    // have reported, and keeping it makes the refinement measurable rather than silent. A pixel
    // still flagged after the last pass keeps its error_ratio and is reported as such.
    out.error_ratio_coarse = coarse.error_ratio;
    out.energy_drift_max_coarse = coarse.energy_drift_max;
    out
}

/// The single pass. `eta` is explicit so the refinement pass can differ from `cfg.eta`.
pub fn evaluate_at<T: Real>(slice: &Slice, idx: usize, cfg: &EnsembleCfg, eta_v: f64) -> PixelOut {
    let copies = jitter::copies_with_path::<T>(
        slice, idx, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
    );
    let n = copies.len();

    // **Masses come from the decode, not from a global.** They used to be `burrau::masses()`
    // read here; four of the five chart families vary them. `m` is the NOMINAL copy's, used for
    // every whole-footprint statistic; per-copy masses are used where a copy is integrated or
    // its energy taken, because on a mass chart two jittered copies are two different systems.
    let m = copies[0].m;

    let t_max = T::lit(cfg.t_max);
    let eta = T::lit(eta_v);

    let base = AzOpts::<T> {
        keep_boundary_shapes: cfg.keep_boundary_shapes,
        forced_refs: None,
        lc_stable: cfg.lc_stable,
        r_coll_frac: T::lit(cfg.r_coll_frac),
        stop_on_event: cfg.stop_on_event,
        escape_every: cfg.escape_every,
        escape_confirm: cfg.escape_confirm,
        r_esc_frac: T::lit(cfg.r_esc_frac),
        escape_all_bodies: cfg.escape_all_bodies,
    };

    // The nominal copy first: its reference-body choices are what the shared policy hands to
    // the others.
    let nominal =
        az::integrate_az_opts(copies[0].s, &copies[0].m, t_max, cfg.n_sync, eta, cfg.max_steps, &base);
    let nominal_refs = nominal.refs.clone();

    let mut outs = Vec::with_capacity(n);
    outs.push(nominal);
    for c in copies.iter().skip(1) {
        let forced = match cfg.ref_policy {
            RefPolicy::Shared => Some(nominal_refs.as_slice()),
            RefPolicy::PerCopy => None,
        };
        outs.push(az::integrate_az_opts(
            c.s, &c.m, t_max, cfg.n_sync, eta, cfg.max_steps,
            &AzOpts { forced_refs: forced, ..base },
        ));
    }

    let e0: Vec<T> = copies
        .iter()
        .map(|c| energy::energy(&c.s.r, &c.s.v, &c.m, T::zero()))
        .collect();
    let et: Vec<T> = outs
        .iter()
        .map(|o| energy::energy(&o.state.r, &o.state.v, &m, T::zero()))
        .collect();
    let (ratio, s0, st) = stats::error_ratio(&e0, &et);
    let ratio_mad = stats::error_ratio_mad(&e0, &et);

    let shapes: Vec<[T; 3]> = outs
        .iter()
        .map(|o| shape::shape_vec(&o.state.r, &m))
        .collect();
    let classes: Vec<u8> = outs
        .iter()
        .map(|o| outcome::classify_legacy(&o.state, &m))
        .collect();
    let outcomes: Vec<Outcome> = outs
        .iter()
        .map(|o| outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted))
        .collect();
    let packed: Vec<u8> = outcomes.iter().map(|o| o.pack()).collect();

    let sp_shape = shape::spread_shape(&shapes).to_f64().unwrap();
    let sv = shape::svar(&shapes).to_f64().unwrap();
    // The event class at each sync boundary, per copy. Evaluated at every boundary rather
    // than only at the end: `t_spread_event` is the whole reason this quantity was chosen
    // over the terminal one.
    let tights: Vec<&[u8]> = outs.iter().map(|o| o.tight.as_slice()).collect();
    let ev_at = |k: usize| -> Vec<u8> {
        tights
            .iter()
            .zip(packed.iter())
            .map(|(t, &term)| stats::event_class_at(t, term, k))
            .collect()
    };
    let per_boundary: Vec<f64> = (0..cfg.n_sync)
        .map(|k| stats::spread_event::<T>(&ev_at(k)).to_f64().unwrap())
        .collect();
    // ---- §5, shape arm: per-boundary spread over the copies ----
    //
    // Ragged by construction: copies terminate at different boundaries under `stop_on_event`,
    // so a copy's record is short and its LAST recorded shape is carried forward. Truncating to
    // the shortest instead would silently discard the boundaries where the surviving copies are
    // doing the diverging, which is the whole quantity.
    let (t_run_max, t_trend, t_first_div) = if cfg.keep_boundary_shapes {
        let bs: Vec<&[[T; 3]]> = outs.iter().map(|o| o.boundary_shapes.as_slice()).collect();
        let n_b = bs.iter().map(|b| b.len()).max().unwrap_or(0);
        let mut series: Vec<(f64, f64)> = Vec::with_capacity(n_b);
        for k in 0..n_b {
            let at: Vec<[T; 3]> = bs
                .iter()
                .filter_map(|b| if b.is_empty() { None } else { Some(b[k.min(b.len() - 1)]) })
                .collect();
            if at.len() < 2 {
                continue;
            }
            let t = cfg.t_max * (k + 1) as f64 / n_b as f64;
            series.push((t, shape::spread_shape(&at).to_f64().unwrap()));
        }
        if series.is_empty() {
            (f64::NAN, f64::NAN, f64::NAN)
        } else {
            let run_max = series.iter().map(|p| p.1).fold(0.0f64, f64::max);
            let first = series
                .iter()
                .find(|p| p.1 >= DIVERGENCE_TRIGGER)
                .map(|p| p.0)
                .unwrap_or(f64::NAN);
            let trend = if series.len() >= 2 {
                let n = series.len() as f64;
                let (st, sy): (f64, f64) =
                    series.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1));
                let (stt, sty): (f64, f64) = series
                    .iter()
                    .fold((0.0, 0.0), |a, p| (a.0 + p.0 * p.0, a.1 + p.0 * p.1));
                let den = n * stt - st * st;
                if den.abs() > 1e-12 { (n * sty - st * sy) / den } else { f64::NAN }
            } else {
                f64::NAN
            };
            (run_max, trend, first)
        }
    } else {
        (f64::NAN, f64::NAN, f64::NAN)
    };

    // The lightness field, nominal copy only. A separate integrator with a separate reference
    // (`tb_ftle.py` sits on `tb.py`), cross-checked at 8.88e-16 -- see `tests/xcheck.rs`.
    let ftle_out = cfg.ftle.as_ref().map(|o| {
        crate::physics::ftle::integrate_full::<T>(
            copies[0].s,
            &m,
            T::lit(cfg.t_max),
            T::lit(cfg.ftle_dt),
            o,
            &crate::physics::ftle::unit_perturbation::<T>(cfg.seed),
        )
    });

    let sp_event = per_boundary[cfg.n_sync - 1];
    let sp_event_max = per_boundary.iter().cloned().fold(0.0f64, f64::max);
    let n_disagree = per_boundary.iter().filter(|&&x| x > 0.0).count() as u16;
    let longest_disagree_run = {
        let (mut best, mut run) = (0u16, 0u16);
        for &x in &per_boundary {
            run = if x > 0.0 { run + 1 } else { 0 };
            best = best.max(run);
        }
        best
    };
    // The latch: the largest value inside any run of at least LATCH_RUN consecutive
    // disagreeing boundaries, joined with the playhead value.
    let spread_event_latched = {
        let mut best = sp_event;
        let (mut run_start, mut run) = (0usize, 0usize);
        for (k, &x) in per_boundary.iter().enumerate() {
            if x > 0.0 {
                if run == 0 {
                    run_start = k;
                }
                run += 1;
                if run >= LATCH_RUN as usize {
                    for &y in &per_boundary[run_start..=k] {
                        best = best.max(y);
                    }
                }
            } else {
                run = 0;
            }
        }
        best
    };

    let k_first = per_boundary.iter().position(|&x| x > 0.0);
    let tie_ratio_at_disagree = k_first
        .map(|k| {
            outs.iter()
                .filter_map(|o| o.tie_ratio.get(k).map(|x| x.to_f64().unwrap()))
                .fold(f64::INFINITY, f64::min)
        })
        .filter(|x| x.is_finite())
        .unwrap_or(f64::NAN);
    let t_spread_event = per_boundary
        .iter()
        .position(|&x| x > 0.0)
        .map(|k| (k + 1) as f64 * cfg.t_max / cfg.n_sync as f64)
        .unwrap_or(f64::NAN);
    let sp_event_terminal = stats::spread_event::<T>(&packed).to_f64().unwrap();
    let sp_event_legacy = stats::spread_event::<T>(&classes).to_f64().unwrap();
    let n_outcome_disagree = packed.iter().filter(|&&c| c != packed[0]).count() as u8;

    // d_min over the ensemble, from finite states only.
    let mut d_ref = f64::INFINITY;
    let mut d_true = f64::INFINITY;
    let mut drift_max: f64 = 0.0;
    let mut gamma_max: f64 = 0.0;
    let mut n_nonfinite = 0u8;
    let mut ref_disagree = 0u32;

    for o in &outs {
        if !o.finite {
            n_nonfinite += 1;
        } else {
            let a = o.d_min_ref.to_f64().unwrap();
            let b = o.d_min_true.to_f64().unwrap();
            if a.is_finite() {
                d_ref = d_ref.min(a);
            }
            if b.is_finite() {
                d_true = d_true.min(b);
            }
        }
        let dr = o.drift.to_f64().unwrap();
        if dr.is_finite() {
            drift_max = drift_max.max(dr);
        }
        let g = o.gamma_max.to_f64().unwrap();
        if g.is_finite() {
            gamma_max = gamma_max.max(g);
        }
        for (k, r) in o.refs.iter().enumerate() {
            if nominal_refs.get(k).is_some_and(|nr| nr != r) {
                ref_disagree += 1;
            }
        }
    }

    let nom = &outs[0];
    PixelOut {
        outcome: packed[0],
        state: outcomes[0].state as u8,
        detail: outcomes[0].detail,
        legacy_class: classes[0],
        binary_id: outcome::binary_id(&nom.state),
        t_end: nom.t_end.to_f64().unwrap(),
        censored: outcomes[0].state != State::Collision && outcomes[0].state != State::Escape,
        d_min_ref: d_ref,
        d_min_true: d_true,
        d_min_gap: d_ref - d_true,
        energy_drift_nominal: nom.drift.to_f64().unwrap(),
        energy_drift_max: drift_max,
        gamma_max,
        shape_vec: [
            shapes[0][0].to_f64().unwrap(),
            shapes[0][1].to_f64().unwrap(),
            shapes[0][2].to_f64().unwrap(),
        ],
        spread_shape: sp_shape,
        svar: sv,
        spread_event: sp_event,
        spread_event_max: sp_event_max,
        spread_event_latched,
        spread_event_terminal: sp_event_terminal,
        spread_event_legacy: sp_event_legacy,
        n_disagree,
        longest_disagree_run,
        tie_ratio_at_disagree,
        t_spread_event,
        ensemble_spread: sp_shape.max(sp_event),
        sigma_e_0: s0.to_f64().unwrap(),
        sigma_e_t: st.to_f64().unwrap(),
        error_ratio: ratio.to_f64().unwrap(),
        error_ratio_mad: ratio_mad.to_f64().unwrap(),
        switches: nom.switches,
        ref_disagree,
        n_outcome_disagree,
        n_nonfinite,
        eta_used: eta_v,
        error_ratio_coarse: ratio.to_f64().unwrap(),
        energy_drift_max_coarse: drift_max,
        refined: false,
        copy_outcomes: if cfg.keep_copy_outcomes { packed.clone() } else { Vec::new() },
        event_class: stats::event_class_at(&outs[0].tight, packed[0], cfg.n_sync - 1),
        total_substeps: outs.iter().map(|o| o.steps as u64).sum(),
        copy_shapes: if cfg.keep_copy_shapes {
            shapes
                .iter()
                .map(|s| {
                    [
                        s[0].to_f64().unwrap(),
                        s[1].to_f64().unwrap(),
                        s[2].to_f64().unwrap(),
                    ]
                })
                .collect()
        } else {
            Vec::new()
        },
        running_max_divergence: t_run_max,
        divergence_trend: t_trend,
        first_divergence_t: t_first_div,
        ftle: ftle_out.map(|o| o.ftle.to_f64().unwrap()).unwrap_or(f64::NAN),
        diffusion: ftle_out.map(|o| o.diffusion.to_f64().unwrap()).unwrap_or(f64::NAN),
        ftle_renorm: ftle_out.map(|o| o.n_renorm).unwrap_or(0),
    }
}
