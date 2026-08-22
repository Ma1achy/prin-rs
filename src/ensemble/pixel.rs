//! Per-pixel evaluation: build the copies, integrate them, reduce.

use crate::grid::Slice;
use crate::integrate::az::{self, AzOpts, RefPolicy};
use crate::outcome::{self, Outcome, State};
use crate::physics::{burrau, energy, shape};
use crate::Real;

use super::{jitter, stats};

#[derive(Clone, Copy, Debug)]
pub struct EnsembleCfg {
    /// `E` extra copies; the pixel always carries `E + 1`.
    pub n_extra: usize,
    pub jitter_frac: f64,
    pub seed: u64,
    pub t_max: f64,
    pub n_sync: usize,
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
}

impl Default for EnsembleCfg {
    fn default() -> Self {
        Self {
            n_extra: 7, // E + 1 = 8, per BRIEF §3
            jitter_frac: 0.5,
            seed: 0,
            t_max: 13.0,
            n_sync: 32,
            eta: 0.01,
            max_steps: 30_000,
            ref_policy: RefPolicy::PerCopy,
            lc_stable: true,
            r_coll_frac: 1e-3,
            stop_on_event: true,
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
}

/// Consecutive boundaries a disagreement must survive before the latch counts it.
///
/// Chosen from the measurement in `examples/latching_decision.rs`, not by eye: at 2 it admits
/// 1 of 130 artefacts, at 3 it admits none, and raising it further changes nothing. On this
/// slice, at this `n_sync`. It is a named constant so a future change to it shows in a diff.
pub const LATCH_RUN: u16 = 3;

pub fn evaluate<T: Real>(slice: &Slice, idx: usize, cfg: &EnsembleCfg) -> PixelOut {
    let m = burrau::masses::<T>();
    let copies = jitter::copies::<T>(slice, idx, cfg.n_extra, cfg.jitter_frac, cfg.seed);
    let n = copies.len();

    let t_max = T::lit(cfg.t_max);
    let eta = T::lit(cfg.eta);

    let base = AzOpts::<T> {
        forced_refs: None,
        lc_stable: cfg.lc_stable,
        r_coll_frac: T::lit(cfg.r_coll_frac),
        stop_on_event: cfg.stop_on_event,
    };

    // The nominal copy first: its reference-body choices are what the shared policy hands to
    // the others.
    let nominal = az::integrate_az_opts(copies[0], &m, t_max, cfg.n_sync, eta, cfg.max_steps, &base);
    let nominal_refs = nominal.refs.clone();

    let mut outs = Vec::with_capacity(n);
    outs.push(nominal);
    for c in copies.iter().skip(1) {
        let forced = match cfg.ref_policy {
            RefPolicy::Shared => Some(nominal_refs.as_slice()),
            RefPolicy::PerCopy => None,
        };
        outs.push(az::integrate_az_opts(
            *c, &m, t_max, cfg.n_sync, eta, cfg.max_steps,
            &AzOpts { forced_refs: forced, ..base },
        ));
    }

    let e0: Vec<T> = copies
        .iter()
        .map(|c| energy::energy(&c.r, &c.v, &m, T::zero()))
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
    }
}
