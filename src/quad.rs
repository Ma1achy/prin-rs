//! The quadtree: quads, their reductions, and the geometry that makes `alpha` mean something.
//!
//! **A quad is a square patch of IC space holding `N × N` footprints** (SCHEDULER_BRIEF §2.1). It
//! is the unit of scheduling and the thing that splits. A footprint is one nominal initial
//! condition carrying `E+1` copies, and it is the unit `ensemble_spread` is defined on.
//!
//! The geometry is load-bearing. A quad at half-width `h` sampled `N × N` has
//! `cell_width = 2h/(N-1)`; its children at `h/2` have `h/(N-1)` — **exactly half**. That factor of
//! two is what `alpha = log2(spread_parent / spread_child)` is measured against, and it holds by
//! construction rather than by convention.
//!
//! **A quad is never synthesised by pooling its children.** With fixed offsets a pooled block is
//! four exact repeats of one pattern at four cell centres, not a wider-footprint ensemble; the
//! measured surrogate error is +38.6%, flat in `E`. The uniform kernel had to pool because it has
//! no tree. This one computes every quad as a real quad at its own cell width.

use crate::grid::{Chart, Slice};

/// Minimum samples per quad axis.
///
/// **Not stylistic.** `Slice::axis` takes an `n <= 1` branch returning `c - half` — the lower-left
/// *corner*, not the centre — and `Slice::cell_widths` clamps to `nx.max(2) - 1`, giving `2*half`
/// rather than a cell width. A one-sample quad would therefore be corner-anchored and jittered by
/// the whole box. Rejected at construction rather than left as a trap.
pub const MIN_SAMPLES_PER_AXIS: usize = 2;

/// Relative margin above `f64::EPSILON` at which the descent is stopped for precision.
///
/// At `half0 = 0.05`, `N = 8` the cell width crosses `f64::EPSILON` at level **45.87**; below that
/// the copies are no longer distinct initial conditions and the spread is pure noise. `1e3` puts
/// the trigger at level **35.90** — comfortably above any physically meaningful descent, and
/// derived from the arithmetic rather than picked as a round number.
pub const PRECISION_MARGIN: f64 = 1e3;

/// What one quad reduces to. One number per field, from `N²` footprints.
#[derive(Clone, Copy, Debug, Default)]
pub struct QuadReduction {
    /// `ensemble_spread` aggregated three ways. **All three are dumped and none is silently
    /// picked** (§3.4): with excess kurtosis 110 a mean is dominated by a single footprint, so the
    /// choice changes decisions and the report has to say by how much.
    pub spread_mean: f64,
    pub spread_median: f64,
    pub spread_p90: f64,
    /// The two contributors, at the median, so a spread can be attributed.
    pub spread_shape_median: f64,
    pub spread_event_median: f64,
    /// **Trust, measured alongside and never enforced** (§2.1). `error_ratio` sees spread, so a
    /// correlated drift is invisible to it — hence the worst absolute drift beside it.
    pub error_ratio_max: f64,
    pub worst_energy_drift: f64,
    pub n_nonfinite: u32,
    pub n_footprints: u32,

    // ---------------------------------------------------------------------------------
    // The between-footprint arm.
    //
    // `spread_*` above are statistics over the `E+1` copies of one footprint, aggregated over
    // the quad. The brief's §1 reads that as a category error — that refinement can only
    // reduce *between*-footprint variation, so a *within*-footprint statistic is the wrong
    // quantity. **The premise does not describe this implementation**, and the reason is
    // measurable: `jitter_frac` is 0.5 and `halton_offset` returns `[-1, 1)^2` scaled by cell
    // width, so the copies span the **whole cell, edge to edge**. They are a quasi-random
    // sample of exactly the area the footprint stands for, not a cloud around a point. The
    // corroboration is already on record: the Halton control's true `alpha` is exactly 1.0,
    // and an irreducible within-point statistic would have `alpha == 0` by construction,
    // because splitting would not shrink it.
    //
    // What genuinely differs is **scale** (cell against quad), **sample count** (`E+1` against
    // `N^2`), and the aggregation. All four numbers below are carried so those can be
    // separated rather than argued about.
    // ---------------------------------------------------------------------------------
    /// `spread_shape` over the `N^2` **nominal** (copy 0, un-jittered) shape vectors. §1.4 as
    /// briefed. Copy 0 only, so between-footprint variation is not contaminated by the
    /// within-footprint jitter.
    pub between_shape: f64,
    /// `spread_event` over the `N^2` nominals' **event class** — not their terminal outcome.
    /// See `PixelOut::event_class` for why that distinction is load-bearing here.
    pub between_event: f64,
    /// `max` of the two, the between-arm analogue of `ensemble_spread`.
    pub between_spread: f64,
    /// [`Self::between_shape`] over only the first `E+1` nominals, **holding the sample count
    /// fixed** so the comparison against `spread_*` isolates scale rather than count.
    ///
    /// Required, not decorative: a spread estimator's expectation depends on its sample size,
    /// and `E+1 = 2` reports 0.539 of `E+1 = 32`'s value in near-field and 0.131 in `far`.
    /// Differencing an 8-sample statistic against a 64-sample one would read that bias as a
    /// scale effect.
    pub between_matched: f64,
    /// `spread_shape` over **all** `N^2 * (E+1)` copies pooled — the within arm at the
    /// between arm's sample count, holding the count fixed while the extent moves. `NaN`
    /// unless `EnsembleCfg::keep_copy_shapes` is set.
    pub within_pooled: f64,

    /// Where the hot footprints sit, on the within-footprint field (`ensemble_spread > tau`).
    pub layout_within: crate::spatial::Layout,
    /// The same, on the per-footprint contribution to [`Self::between_shape`] — each nominal's
    /// distance from the quad's nominal centroid, halved. That is the one between-arm quantity
    /// that is defined *per footprint* and so can carry a mask at all.
    pub layout_between: crate::spatial::Layout,
    /// Fraction of footprints above `tau`, both arms. §3.1: the direct form of "does this quad
    /// contain a boundary?" is a **count in the tail**, not a quantile of the distribution. A
    /// quad with 5% hot footprints has a filament; one with a high median is uniformly
    /// blurred, and every quantile conflates them.
    pub frac_above_tau_within: f64,
    pub frac_above_tau_between: f64,

    /// The same two layouts under the **relative** hot rule — above the quad's own quantile
    /// rather than above `tau`. See [`crate::spatial::HotRule`] for why both are carried.
    ///
    /// **`n_hot` here is a constant**, `N^2/2` at the median, on every quad by construction. Do
    /// not read `frac_hot` off these: the signal is entirely `n_components`,
    /// `largest_component` and `perimeter_ratio`, which is what desaturating the mask buys.
    pub layout_rel_within: crate::spatial::Layout,
    pub layout_rel_between: crate::spatial::Layout,
    /// RMS forward-difference gradient of the two per-footprint fields. The magnitude companion
    /// to the layouts, and the only structure measure here that needs **no threshold at all** —
    /// which is why it is worth having beside two that do. `NaN`, never 0, when no adjacent
    /// pair is finite.
    pub grad_rms_within: f64,
    pub grad_rms_between: f64,

    /// Footprints whose nominal **terminated** before the horizon — collision *or* escape.
    ///
    /// **Not "escaped", which is what an earlier draft of this called it and got wrong.** §3.5
    /// asks for a `t_end` gradient, and `t_end` is set by whichever terminating event came
    /// first. In `deep interior` this reads **0.99** while the escape arm is silent: those are
    /// collisions. Quoting it as an escape fraction would have contradicted a standing result
    /// ("zero of 1024 near-field pixels escape at `t_max = 13`") while agreeing with it, which
    /// is worse than either.
    ///
    /// [`Self::escape_fraction`] is carried separately so the two cannot be confused again.
    pub terminated_fraction: f64,
    /// Footprints whose nominal terminated **by escape** specifically.
    pub escape_fraction: f64,
    /// Mean absolute spatial gradient of nominal `t_end`, over the **terminated** subset only.
    ///
    /// **`NaN` when nothing terminated**, never 0. Censoring is the failure mode: `t_end` pinned
    /// at the horizon carries no gradient information, so a gradient computed over censored
    /// footprints would be a gradient of the horizon constant — exactly zero, everywhere, and
    /// indistinguishable from a smooth field.
    pub t_end_gradient: f64,
    /// Total integrator substeps over the quad — the cost side of §8's cost-aware priority.
    pub total_substeps: u64,
    /// Distinct **initial conditions** among the `N^2` footprints, by exact bitwise comparison
    /// of the full state. Equal to `n_footprints` on any healthy quad.
    ///
    /// Read before any spread is read. `decode::distinct` is the existing guard, and the rule
    /// it enforces is that a difference can be small because both sides are right or because
    /// both are dead — count distinct ICs first, read divergence second.
    pub n_distinct_ic: u32,

    /// §5, shape arm, aggregated over footprints. `NaN` unless
    /// `EnsembleCfg::keep_boundary_shapes` is set — never 0.
    pub running_max_divergence_median: f64,
    pub divergence_trend_median: f64,
    /// Fraction of footprints whose copies ever crossed the divergence trigger.
    ///
    /// The *time* is the quantity that cannot saturate; a fraction can, at 1.0. Both are
    /// carried: the fraction is what a ranking can read without a NaN convention, and
    /// [`Self::first_divergence_median`] is the one to read when asking whether the signal is
    /// still informative in a saturated region.
    pub frac_diverged: f64,
    /// Median first-divergence time over the footprints that crossed. `NaN` if none did.
    pub first_divergence_median: f64,
}

impl QuadReduction {
    /// The aggregate a decision reads, by policy.
    pub fn spread(&self, agg: Agg) -> f64 {
        match agg {
            Agg::Mean => self.spread_mean,
            Agg::Median => self.spread_median,
            Agg::P90 => self.spread_p90,
        }
    }

    /// The scalar a decision reads, by criterion and aggregation.
    ///
    /// `agg` is ignored by every criterion but [`Criterion::Within`]: the between-footprint and
    /// layout signals are already one number per quad, with no `N^2` distribution left to
    /// aggregate. That asymmetry is the point — half the current parameter surface exists only
    /// because the within arm keeps a distribution it then throws away.
    pub fn signal(&self, criterion: Criterion, agg: Agg) -> f64 {
        match criterion {
            Criterion::Within => self.spread(agg),
            Criterion::Between => self.between_spread,
            Criterion::MaxOfBoth => self.spread(agg).max(self.between_spread),
            Criterion::FracHotWithin => self.frac_above_tau_within,
            Criterion::FracHotBetween => self.frac_above_tau_between,
            Criterion::RunningMax => self.running_max_divergence_median,
            Criterion::FirstDivergence => self.frac_diverged,
            Criterion::TerminationGradient => {
                // NaN where nothing escaped. Deliberately NOT mapped to 0 here: a caller
                // ranking on this must decide what an undetermined quad means, and silently
                // calling it "no structure" is the failure this signal is most prone to.
                self.t_end_gradient
            }
            Criterion::Layout => {
                // Thin and connected reads as a boundary and must outrank scatter of the same
                // count; scattered hot footprints are chaos, which no refinement resolves.
                let l = self.layout_within;
                if l.n_hot == 0 {
                    0.0
                } else {
                    l.frac_hot(self.n_side()) * (l.largest_component as f64 / l.n_hot as f64)
                }
            }
            Criterion::LayoutRel => {
                // The same reading on the relative mask -- but `frac_hot` is a constant there,
                // so it is dropped rather than carried as a scale factor that varies with
                // nothing. What is left is connectedness alone, which is the whole point: a
                // relative mask turns a magnitude statistic into a shape one.
                let l = self.layout_rel_within;
                if l.n_hot == 0 {
                    0.0
                } else {
                    l.largest_component as f64 / l.n_hot as f64
                }
            }
            Criterion::GradRms => self.grad_rms_within,
        }
    }

    /// **How much this quad looks like a BOUNDARY rather than a uniform sea**, in `[0, 1]`.
    ///
    /// The §1 diagnosis: `ensemble_spread` measures *uncertainty*, and a uniformly chaotic quad
    /// and a filament are both uncertain. Only one of them repays refining. This is the term
    /// that separates them, and it is two factors because either alone is fooled:
    ///
    /// - **connectedness**, `largest_component / n_hot`. Scattered hot footprints are chaos;
    ///   one run of them is a structure. Without this a checkerboard scores maximum thinness.
    /// - **thinness**, `perimeter_ratio / 2`, clamped at 1. A one-cell-wide filament reads
    ///   exactly `2.0` under the internal-edges convention, a compact blob `~4/sqrt(A)`, and a
    ///   featureless fully-hot quad exactly **0**. Without this a fully-hot quad scores maximum
    ///   connectedness.
    /// - **extent**, `largest_component / N`. **This third factor was not anticipated; the test
    ///   found it.** A single isolated hot cell is trivially connected (it *is* the largest
    ///   component) and maximally thin (`perimeter_ratio == 4`), so the first two factors scored
    ///   it **1.0** — maximum structure, for one cell. A boundary crossing a quad spans it;
    ///   `Layout::looks_like_boundary` already encodes that as `largest_component >= N/2`, and
    ///   this is the graded form. An isolated cell now reads `1/N`.
    ///
    /// Read it on the **relative** mask: on the absolute one `n_hot == N^2` in 98.8% of committed
    /// leaves, so `perimeter_ratio` is 0 and this is identically zero — a term that cannot fire.
    ///
    /// `NaN` when nothing is hot, following `perimeter_ratio`'s own convention. **Not 0**: an
    /// empty mask is "not determined", and `far` is the case that makes the difference — its
    /// absolute mask is empty on every leaf, and a 0 there would read as "no structure found"
    /// rather than "not measured".
    pub fn structure(&self, relative: bool) -> f64 {
        let l = if relative { self.layout_rel_within } else { self.layout_within };
        if l.n_hot == 0 {
            return f64::NAN;
        }
        let connected = l.largest_component as f64 / l.n_hot as f64;
        let thin = (l.perimeter_ratio / 2.0).clamp(0.0, 1.0);
        let extent = (l.largest_component as f64 / self.n_side().max(1) as f64).clamp(0.0, 1.0);
        connected * thin * extent
    }

    /// The scalar a decision reads once [`StructureMode`] is applied.
    ///
    /// **`Multiply` is "uncertain AND structured"**, which floors the uniform sea by
    /// construction; `Replace` says structure is the whole answer. The recommendation on record
    /// is multiply, and the recommendation does not decide — `error(B)` does.
    ///
    /// A `NaN` structure term propagates rather than being coerced to 0 or 1. Under `Multiply`
    /// that demotes an undetermined quad to the bottom of any ranking, which is the wrong
    /// standing for "could not be measured" and is why the mode is a measured choice rather than
    /// an obvious one. It is stated here so it is read off the code rather than discovered.
    pub fn signal_with(&self, criterion: Criterion, agg: Agg, mode: StructureMode) -> f64 {
        let base = self.signal(criterion, agg);
        match mode {
            StructureMode::Off => base,
            StructureMode::Replace => self.structure(true),
            StructureMode::Multiply => base * self.structure(true),
        }
    }

    /// `N`, recovered from the footprint count. The reduction does not carry the grid width
    /// separately, and every quad is square by construction (`Quad::slice` builds `n x n`).
    pub fn n_side(&self) -> usize {
        (self.n_footprints as f64).sqrt().round() as usize
    }

    /// Is this quad's between-footprint arm **collapsed** — every nominal bitwise identical?
    ///
    /// A collapsed decode gives a spread of exactly zero, which reads as "perfectly resolved"
    /// and stops the descent with a small tidy tree built from nothing. Treated as
    /// **undetermined**, the same way a non-finite copy is a measurement outcome rather than
    /// missing data.
    pub fn between_collapsed(&self) -> bool {
        self.n_distinct_ic < self.n_footprints
    }
}

/// Whether the spatial-structure term enters the decision, and how.
///
/// §2.2's open question, implemented as both variants because it is the sort of thing this
/// project has settled by measurement rather than argument.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StructureMode {
    /// The signal alone. Today's behaviour, and the control.
    #[default]
    Off,
    /// Structure is the whole answer.
    ///
    /// **This mode has no criterion axis**: `signal_with(_, _, Replace)` discards both arguments,
    /// so `replace x within` and `replace x between` are the same ranking, and both are
    /// identically `Rank::StructureOnly`. Measured and confirmed — their `error(B)` curves match
    /// to five digits on every target. That is a structural identity, not a finding, and it is
    /// stated here so a table carrying both does not read as two independent rows agreeing.
    Replace,
    /// *Uncertain **and** structured.* Keeps the determinacy question `ensemble_spread` answers
    /// while adding the structure question it cannot.
    Multiply,
}

impl StructureMode {
    pub fn name(self) -> &'static str {
        match self {
            StructureMode::Off => "off",
            StructureMode::Replace => "replace",
            StructureMode::Multiply => "multiply",
        }
    }
    pub fn parse(s: &str) -> Option<StructureMode> {
        Some(match s {
            "off" => StructureMode::Off,
            "replace" => StructureMode::Replace,
            "multiply" => StructureMode::Multiply,
            _ => return None,
        })
    }
    pub const ALL: [StructureMode; 3] =
        [StructureMode::Off, StructureMode::Replace, StructureMode::Multiply];
}

/// Which signal a split decision reads.
///
/// All are computed and dumped on every run whatever this is set to — the point is to compare
/// criteria offline without re-integrating, and the marginal cost of every one of them is
/// `O(N^2)` against 512 trajectories per quad.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Criterion {
    /// `ensemble_spread` aggregated over footprints — the criterion as it stands. The default
    /// until the §2 metric says otherwise; a change of default is a measured decision, not a
    /// tidy-up.
    #[default]
    Within,
    /// The between-footprint arm alone.
    Between,
    /// `max` of the two. Both quantities are wanted and they answer different questions —
    /// within is the trust / display-honesty signal, between is the refinement signal — so
    /// this is the conservative join rather than a replacement.
    MaxOfBoth,
    /// Fraction of footprints above `tau`, within arm (§3.1).
    FracHotWithin,
    /// The same on the between arm.
    FracHotBetween,
    /// Hot-set layout (§3.2): connectedness weighted by count.
    Layout,
    /// Running max of the shape spread over boundaries (§5) — catches divergence that has
    /// already happened and then subsided. The shape spread was measured falling **6x**
    /// between `t = 6` and `t = 8`, so an instantaneous read genuinely misses it.
    RunningMax,
    /// Fraction of footprints whose copies ever crossed the divergence trigger (§5).
    ///
    /// **The sign of its usefulness is not obvious and is left to the measurement.** A quad
    /// that diverges early is uncertain, but chaotic uncertainty is exactly what refinement
    /// cannot reduce. Ranking it high and ranking it low are both defensible before the fact,
    /// so the §2 curve decides rather than an argument.
    FirstDivergence,
    /// Spatial gradient of nominal `t_end` (§3.5) — a boundary detector needing no ensemble.
    ///
    /// Admissible only where `terminated_fraction` is nonzero, and that column is reported
    /// beside it. In `near-field` at `t = 13` it is `NaN` on **97.1%** of quads, which is a
    /// property to read, not a defect to hide: NaN never wins a comparison, so the ranking is
    /// decided entirely by the 2.9% it does score.
    TerminationGradient,
    /// Hot-set layout on the **relative** mask. The desaturated twin of [`Criterion::Layout`];
    /// connectedness alone, since `frac_hot` is constant under a quantile rule.
    LayoutRel,
    /// RMS spatial gradient of `ensemble_spread` across the footprint grid.
    ///
    /// **The only candidate here with no threshold in it.** That makes it the control on the
    /// whole hot-mask family: if a masked signal cannot beat it, the mask is not earning its
    /// parameter.
    GradRms,
}

impl Criterion {
    pub fn name(self) -> &'static str {
        match self {
            Criterion::Within => "within",
            Criterion::Between => "between",
            Criterion::MaxOfBoth => "max_of_both",
            Criterion::FracHotWithin => "frac_hot_within",
            Criterion::FracHotBetween => "frac_hot_between",
            Criterion::Layout => "layout",
            Criterion::RunningMax => "running_max",
            Criterion::FirstDivergence => "first_div",
            Criterion::TerminationGradient => "term_grad",
            Criterion::LayoutRel => "layout_rel",
            Criterion::GradRms => "grad_rms",
        }
    }
    pub fn parse(s: &str) -> Option<Criterion> {
        Some(match s {
            "within" => Criterion::Within,
            "between" => Criterion::Between,
            "max_of_both" => Criterion::MaxOfBoth,
            "frac_hot_within" => Criterion::FracHotWithin,
            "frac_hot_between" => Criterion::FracHotBetween,
            "layout" => Criterion::Layout,
            "running_max" => Criterion::RunningMax,
            "first_div" => Criterion::FirstDivergence,
            "term_grad" => Criterion::TerminationGradient,
            "layout_rel" => Criterion::LayoutRel,
            "grad_rms" => Criterion::GradRms,
            _ => return None,
        })
    }
    /// Every variant, for sweeps that must not silently omit one.
    pub const ALL: [Criterion; 11] = [
        Criterion::Within,
        Criterion::Between,
        Criterion::MaxOfBoth,
        Criterion::FracHotWithin,
        Criterion::FracHotBetween,
        Criterion::Layout,
        Criterion::RunningMax,
        Criterion::FirstDivergence,
        Criterion::TerminationGradient,
        Criterion::LayoutRel,
        Criterion::GradRms,
    ];
}

/// Which aggregation of the `N²` footprint spreads a decision uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Agg {
    Mean,
    /// Default: robust to the heavy tail (excess kurtosis 110).
    #[default]
    Median,
    P90,
}

impl Agg {
    pub fn name(self) -> &'static str {
        match self {
            Agg::Mean => "mean",
            Agg::Median => "median",
            Agg::P90 => "p90",
        }
    }
    pub fn parse(s: &str) -> Option<Agg> {
        Some(match s {
            "mean" => Agg::Mean,
            "median" => Agg::Median,
            "p90" => Agg::P90,
            _ => return None,
        })
    }
}

/// What the scheduler decided about a quad, and why.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Decision {
    /// Computed but not yet decided.
    #[default]
    Pending,
    /// Refining reduced the spread; go deeper.
    Split,
    /// Spread is high and refining will not reduce it. **Flagged and not re-queued** — this is the
    /// branch that must be shown to engage.
    Floor,
    /// Spread is already low at both scales.
    Keep,
    /// The cell width is approaching `f64::EPSILON`. A numerical stop, not a physical one, and it
    /// must never be reported as "the descent did not terminate".
    PrecisionFloor,
    MaxLevel,
    /// Wanted to split; the budget ran out first. Reported, never silently dropped.
    BudgetExhausted,
    /// The quad's tiles have shrunk to pixel size. **The everyday stop** — in normal
    /// exploration this fires far shallower than any precision floor, and PR #11 had no
    /// camera to apply it. View-relative: the same quad refines again when zoomed into.
    ScreenFloor,
    /// `level >= camera_depth + MAX_REL_DEPTH`. Replaces absolute `MaxLevel`, which caps
    /// infinite zoom at ~14. Scheduler state, never on the sim key.
    MaxRelDepth,
    /// **The decode has collapsed**: fewer distinct initial conditions than footprints, so the
    /// quad's `N^2` samples are repeats and every spread computed over them is spread over
    /// nothing.
    ///
    /// A separate decision rather than a `Floor`, because the failure it names is invisible
    /// otherwise: identical footprints give a spread of exactly **zero**, which reads as
    /// "perfectly resolved" and terminates the descent with a small tidy tree built from no
    /// information. Treated as *undetermined* — the same standing that a non-finite copy has as
    /// a measurement outcome rather than missing data — and it must be countable in the dump,
    /// which a `Floor` would not be.
    ///
    /// Tested on **initial conditions**, never on the spread being zero: a genuinely uniform
    /// region has a zero spread over perfectly distinct ICs, and conflating the two would flag
    /// the physics as a numerical failure.
    Collapsed,
    /// **Split to satisfy the 2:1 balance constraint, not because the criterion asked.**
    ///
    /// No two adjacent leaves may differ by more than one level or the adaptive render has
    /// cracks. That forces splits the criterion declined, and they are spent from the same
    /// budget — so a run where most of the budget went on geometry rather than physics must be
    /// *countable*, not inferred. A `Split` here would be indistinguishable from a
    /// criterion-driven one, which is the same failure the stop-reason column exists to prevent.
    BalanceForced,
}

impl Decision {
    pub fn name(self) -> &'static str {
        match self {
            Decision::Pending => "pending",
            Decision::Split => "split",
            Decision::Floor => "floor",
            Decision::Keep => "keep",
            Decision::PrecisionFloor => "precision_floor",
            Decision::MaxLevel => "max_level",
            Decision::BudgetExhausted => "budget_exhausted",
            Decision::ScreenFloor => "screen_floor",
            Decision::MaxRelDepth => "max_rel_depth",
            Decision::Collapsed => "collapsed",
            Decision::BalanceForced => "balance",
        }
    }
    pub fn code(self) -> u8 {
        match self {
            Decision::Pending => 0,
            Decision::Split => 1,
            Decision::Floor => 2,
            Decision::Keep => 3,
            Decision::PrecisionFloor => 4,
            Decision::MaxLevel => 5,
            Decision::BudgetExhausted => 6,
            Decision::ScreenFloor => 7,
            Decision::MaxRelDepth => 8,
            Decision::Collapsed => 9,
            Decision::BalanceForced => 10,
        }
    }
}

/// One node. A leaf is a quad with no children.
#[derive(Clone, Debug)]
pub struct Quad {
    pub level: u32,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    pub parent: Option<usize>,
    pub children: Option<[usize; 4]>,
    /// Which child of its parent this is, `0..4`; `0` for the root.
    pub sib_index: u8,
    /// Descent iteration at which this quad was computed — gives leaf-count against iteration.
    pub iteration: u32,
    pub red: QuadReduction,
    /// `log2(spread_parent / spread_self)` at the decision aggregation. `None` at the root, which
    /// has no parent and therefore no exponent — that is why the first level or two split
    /// unconditionally.
    pub alpha: Option<f64>,
    /// The same exponent under the other two aggregations, so §3.4's sensitivity is measurable
    /// without a second run.
    pub alpha_mean: Option<f64>,
    pub alpha_p90: Option<f64>,
    /// **The reliability signal** (§3.3). Set on a *parent* once all four children are computed:
    /// the **range** (max − min) of their four `alpha` values.
    ///
    /// Range, stated as such rather than dressed as a robust estimator — with four samples an
    /// interdecile is meaningless. It is itself a noisy statistic, and if this policy looks
    /// promising that noise is the next thing to characterise rather than the first thing to
    /// trust.
    pub alpha_sibling_spread: Option<f64>,
    pub decision: Decision,
}

impl Quad {
    /// The `Slice` that samples this quad's `N × N` footprints.
    ///
    /// Panics below [`MIN_SAMPLES_PER_AXIS`] rather than silently corner-anchoring.
    pub fn slice(&self, n: usize, body: usize, chart: Chart) -> Slice {
        assert!(
            n >= MIN_SAMPLES_PER_AXIS,
            "samples per quad axis must be >= {MIN_SAMPLES_PER_AXIS}; \
             Slice::axis corner-anchors at n <= 1 and cell_widths clamps to 2*half"
        );
        Slice::body_plane(n, n, self.cx, self.cy, self.half, body).with_chart(chart)
    }

    /// Footprint spacing: `2*half/(n-1)`. Exactly halves per level.
    pub fn cell_width(&self, n: usize) -> f64 {
        2.0 * self.half / (n - 1) as f64
    }

    /// Has the cell width fallen far enough that the copies are no longer distinct ICs?
    ///
    /// Scaled by the coordinate magnitude, since epsilon is relative.
    pub fn below_precision_floor(&self, n: usize) -> bool {
        let scale = self.cx.abs().max(self.cy.abs()).max(1.0);
        self.cell_width(n) < PRECISION_MARGIN * f64::EPSILON * scale
    }

    /// The four child boxes, in `(jy, jx)` order: lower-left, lower-right, upper-left, upper-right.
    pub fn child_boxes(&self) -> [(f64, f64, f64); 4] {
        let q = self.half / 2.0;
        [
            (self.cx - q, self.cy - q, q),
            (self.cx + q, self.cy - q, q),
            (self.cx - q, self.cy + q, q),
            (self.cx + q, self.cy + q, q),
        ]
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }
}

/// The tree. A flat arena; indices are stable and children are contiguous.
#[derive(Clone, Debug, Default)]
pub struct QuadTree {
    pub nodes: Vec<Quad>,
    /// Samples per quad axis, `N`.
    pub n: usize,
    pub body: usize,
    /// The chart every quad in this tree decodes through. One tree, one chart.
    pub chart: Chart,
}

impl QuadTree {
    pub fn new(cx: f64, cy: f64, half: f64, n: usize, body: usize) -> Self {
        Self::with_chart(cx, cy, half, n, body, Chart::BodyPlane)
    }

    pub fn with_chart(cx: f64, cy: f64, half: f64, n: usize, body: usize, chart: Chart) -> Self {
        assert!(n >= MIN_SAMPLES_PER_AXIS, "samples per quad axis must be >= {MIN_SAMPLES_PER_AXIS}");
        let root = Quad {
            level: 0,
            cx,
            cy,
            half,
            parent: None,
            children: None,
            sib_index: 0,
            iteration: 0,
            red: QuadReduction::default(),
            alpha: None,
            alpha_mean: None,
            alpha_p90: None,
            alpha_sibling_spread: None,
            decision: Decision::Pending,
        };
        QuadTree { nodes: vec![root], n, body, chart }
    }

    /// Create four children of `i`. Returns their indices. Does **not** compute them.
    pub fn split(&mut self, i: usize, iteration: u32) -> [usize; 4] {
        assert!(self.nodes[i].children.is_none(), "quad {i} already split");
        let (level, boxes) = (self.nodes[i].level + 1, self.nodes[i].child_boxes());
        let base = self.nodes.len();
        for (k, (cx, cy, half)) in boxes.into_iter().enumerate() {
            self.nodes.push(Quad {
                level,
                cx,
                cy,
                half,
                parent: Some(i),
                children: None,
                sib_index: k as u8,
                iteration,
                red: QuadReduction::default(),
                alpha: None,
                alpha_mean: None,
                alpha_p90: None,
                alpha_sibling_spread: None,
                decision: Decision::Pending,
            });
        }
        let kids = [base, base + 1, base + 2, base + 3];
        self.nodes[i].children = Some(kids);
        kids
    }

    pub fn leaves(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.nodes.len()).filter(|&i| self.nodes[i].is_leaf())
    }

    /// The same-or-coarser neighbour across one edge, or `None` at the root box's border.
    ///
    /// Descends from the root toward a point just outside the edge, stopping at the deepest
    /// existing node that is **no finer than the query quad**. Same-or-coarser is the right
    /// answer rather than a limitation: a finer neighbour is several quads, and picking one of
    /// them would make the contrast depend on which.
    ///
    /// `O(level)`, no neighbour pointers, nothing cached. The only adjacency code that existed
    /// before this was an `O(n^2)` geometric box-touch test inside `examples/sched_thrash.rs`;
    /// `tests/criterion.rs` checks this against that predicate over a whole tree, which is a
    /// comparison against an independent implementation rather than against itself.
    pub fn neighbour(&self, i: usize, dir: Dir) -> Option<usize> {
        let q = &self.nodes[i];
        // Just outside the edge. A fraction of the quad's own half-width, so the probe scales
        // with the box and cannot fall through a neighbour at any depth.
        let e = q.half * 1e-6;
        let (px, py) = match dir {
            Dir::NegX => (q.cx - q.half - e, q.cy),
            Dir::PosX => (q.cx + q.half + e, q.cy),
            Dir::NegY => (q.cx, q.cy - q.half - e),
            Dir::PosY => (q.cx, q.cy + q.half + e),
        };

        let root = &self.nodes[0];
        if (px - root.cx).abs() > root.half || (py - root.cy).abs() > root.half {
            return None;
        }

        let mut cur = 0usize;
        while self.nodes[cur].level < q.level {
            let Some(kids) = self.nodes[cur].children else { break };
            let node = &self.nodes[cur];
            let jx = usize::from(px >= node.cx);
            let jy = usize::from(py >= node.cy);
            cur = kids[jy * 2 + jx];
        }
        Some(cur)
    }

    /// Neighbour contrast: `max` over the four edges of `|signal_self - signal_neighbour|`.
    ///
    /// **Computed at decision time and never stored on a `Quad`.** It is a relative quantity,
    /// and freezing one onto a node is the mistake the screen floor's "never cached as a quad
    /// fact" rule exists to prevent — a neighbour can be split after this is read.
    ///
    /// Returns the contrast and how many edges contributed. The count matters: a quad on the
    /// root border has fewer neighbours, so its contrast is a max over a smaller set and is
    /// biased low by construction. Reported, never silently absorbed.
    pub fn contrast(&self, i: usize, criterion: Criterion, agg: Agg) -> (f64, u8) {
        let me = self.nodes[i].red.signal(criterion, agg);
        let (mut best, mut count) = (0.0f64, 0u8);
        for d in Dir::ALL {
            if let Some(j) = self.neighbour(i, d) {
                if self.nodes[j].red.n_footprints == 0 {
                    continue; // not yet computed; not a zero contrast
                }
                let v = self.nodes[j].red.signal(criterion, agg);
                if v.is_finite() && me.is_finite() {
                    best = best.max((me - v).abs());
                    count += 1;
                }
            }
        }
        (if count == 0 { f64::NAN } else { best }, count)
    }

    pub fn depth_histogram(&self) -> Vec<usize> {
        let mut h = Vec::new();
        for i in self.leaves() {
            let l = self.nodes[i].level as usize;
            if h.len() <= l {
                h.resize(l + 1, 0);
            }
            h[l] += 1;
        }
        h
    }
}

/// One of the four edges of a quad.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dir {
    NegX,
    PosX,
    NegY,
    PosY,
}

impl Dir {
    pub const ALL: [Dir; 4] = [Dir::NegX, Dir::PosX, Dir::NegY, Dir::PosY];
    pub fn opposite(self) -> Dir {
        match self {
            Dir::NegX => Dir::PosX,
            Dir::PosX => Dir::NegX,
            Dir::NegY => Dir::PosY,
            Dir::PosY => Dir::NegY,
        }
    }
}

/// Aggregate `N²` footprint spreads into the three numbers a `QuadReduction` carries.
pub fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[(((v.len() - 1) as f64) * q).round() as usize]
}
