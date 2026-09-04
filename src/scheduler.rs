//! The descent loop — the one thing the uniform kernel was deliberately built without.
//!
//! Every measurement before this ran the criterion on **one split in isolation**. Everything here
//! exists because the remaining questions are dynamic: does the descent terminate, does the floor
//! engage, does a budget get spent well, does per-quad noise cause thrash.
//!
//! **Scope discipline**: no eviction, no caching, no async, no promotion, no interaction. A quad is
//! computed once and the tree keeps it — that is the tree holding its own data, not a cache. If any
//! of the others appears here, it is a bug.
//!
//! The **camera is now in scope**, and only in scope as a veto: SCHEDULER_BRIEF §6 excluded it, and
//! that exclusion is exactly what made PR #11's q1/q2/q3/q7 describe a regime the real system never
//! enters. See [`crate::camera`].

use rayon::prelude::*;

use crate::camera::Camera;
use crate::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use crate::grid::Chart;
use crate::ensemble::stats;
use crate::physics::shape;
use crate::quad::{quantile, Agg, Criterion, Decision, Dir, QuadReduction, QuadTree, StructureMode};
use crate::spatial::{self, HotRule, Layout};
use crate::render::Precision;
use crate::rng::SplitMix64;

/// How a quad that passed every guard is decided.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Policy {
    /// **The tolerance policy, the default from Phase 2 of the refinement rebuild.** A quad is
    /// `Keep` when every footprint is resolved at `tau_display` — its copies, which span the
    /// whole cell, agree on the shape sphere to within the tolerance and on the event class
    /// exactly — and `Split` otherwise. No exponent, no aggregate: a count in the tail, so a
    /// filament crossing one footprint column of a quad refines it where a median reads it
    /// resolved. A quad outranked by `k_frac` is `Deferred` and re-decided next round, never
    /// dropped as `Keep`.
    ///
    /// The legacy policies split only where `alpha >= alpha_hi`, i.e. where halving the cell had
    /// *already* halved the spread — smooth, converging regions — and floored or kept exactly the
    /// quads whose spread does not fall: discontinuities and fractal mixing at every scale coarser
    /// than their filaments. Measured on `preset_shape_h1`: the refinement went into the smooth
    /// regular island and floored the fractal core at level 2, Spearman(depth, terminated) −0.68.
    /// That is the inverse of "refine the filaments", and it was the design.
    #[default]
    Tolerance,
    /// Threshold on `alpha`'s **value**. Separation between region types is 0.9862 against a
    /// chaotic scatter of 1.1–1.3 — marginal. **Legacy**: every `.prnq` committed before Phase 2
    /// was cut with it, and `tests/live_decision.rs` pins that it still reproduces them.
    Alpha,
    /// Threshold on `alpha_sibling_spread`, the range of the four children's exponents.
    /// Separation in `alpha`'s **reliability** is 0.001 against 1.2 — three orders. Where the four
    /// scatter, the unreliability *is* the answer, and no trustworthy `alpha` is needed.
    Sibling,
}

impl Policy {
    pub fn name(self) -> &'static str {
        match self {
            Policy::Tolerance => "tolerance",
            Policy::Alpha => "alpha",
            Policy::Sibling => "sibling",
        }
    }
    pub fn parse(s: &str) -> Option<Policy> {
        Some(match s {
            "tolerance" => Policy::Tolerance,
            "alpha" => Policy::Alpha,
            "sibling" => Policy::Sibling,
            _ => return None,
        })
    }
}

/// Queue order, and the control that says whether order mattered at all.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Order {
    /// By quad spread.
    #[default]
    Spread,
    /// By spread × area — pays for what is visible rather than what is merely uncertain.
    SpreadArea,
    /// The control (§3.6): same budget, no priority.
    Shuffled,
}

impl Order {
    pub fn name(self) -> &'static str {
        match self {
            Order::Spread => "spread",
            Order::SpreadArea => "spread_area",
            Order::Shuffled => "shuffled",
        }
    }
    pub fn parse(s: &str) -> Option<Order> {
        Some(match s {
            "spread" => Order::Spread,
            "spread-area" | "spread_area" => Order::SpreadArea,
            "shuffled" => Order::Shuffled,
            _ => return None,
        })
    }
}

/// The two modes, which are **one mechanism at different budgets** (§3.1).
///
/// The reframe that matters: *the screen floor stops things; the criterion decides what gets
/// attention first.* The criterion was never a stop condition.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mode {
    /// Every quad that passes the guards splits, to the veto. The criterion is **off**, and this
    /// is the control that must degenerate to uniform depth — without it, a balanced mode that
    /// was merely frozen would look the same in a depth-variance plot.
    Uniform,
    /// The same descent, frontier **ranked**, top `k` per round gets budget. The criterion is a
    /// **priority ordering**, never a threshold, so it cannot land above or below the
    /// distribution the way a fixed `tau` does.
    #[default]
    Balanced,
}

impl Mode {
    pub fn name(self) -> &'static str {
        match self {
            Mode::Uniform => "uniform",
            Mode::Balanced => "balanced",
        }
    }
    pub fn parse(s: &str) -> Option<Mode> {
        Some(match s {
            "uniform" => Mode::Uniform,
            "balanced" => Mode::Balanced,
            _ => return None,
        })
    }
}

/// **The uniform-mode control.** `k_frac = 1.0` takes the top 100% of the split-eligible
/// frontier, so [`Mode::Balanced`] ranks the queue and then refines all of it. The ranking runs
/// and changes nothing.
///
/// It is not *identical* to [`Mode::Uniform`] — uniform returns `Split` unconditionally and
/// bypasses the `tau` and `alpha` gates, where balanced still applies both — but the **rank
/// truncation**, which is the mechanism the criterion work exists to add, never engages. Keep it
/// available and keep it labelled; it was the silent default through PR #21 and every committed
/// dump outside `results/sweep` carries it.
pub const K_FRAC_UNRANKED: f64 = 1.0;

/// **The default.** Chosen from the widened `tau x k_frac` sweep, not from a picture.
///
/// Measured over `near-field`, `deep interior` and `preset_shape`: `tau` moves leaf-level depth
/// variance by essentially nothing (1.900 -> 1.866 across a decade at `k = 0.5`), while `k_frac`
/// from 1.0 to 0.25 moves it 0.577 -> 2.109. `tau` is a threshold and cannot land inside a
/// distribution whose median moves six orders between regions; `k_frac` is a rank and always
/// cuts through it. See `RESULTS.md` §19 and `examples/criterion_sweep.rs`.
pub const K_FRAC_RANKED: f64 = 0.25;

/// **Refuse to write a headline artefact from the degenerate configuration.**
///
/// `Mode::Balanced` with [`K_FRAC_UNRANKED`] is uniform mode wearing balanced mode's name: the
/// priority is computed, the queue is sorted, and every element of it is refined anyway. A render
/// made that way looks exactly like a render of the ranked frontier and is not one — which has
/// now happened twice in this project, first in PR #18 and again in every committed chart dump.
///
/// This is the `preset_control.rs` pattern at a second site: a configuration that silently
/// reproduces the old behaviour needs a guard, not a convention. Call it in any example that
/// writes under `results/`; pass `--allow-unranked` (or set `allow` yourself) only when the
/// unranked run **is** the control being measured, as it is in `criterion_sweep`.
pub fn assert_not_uniform_in_disguise(cfg: &SchedCfg, path: &str, allow: bool) {
    if allow || cfg.mode != Mode::Balanced || cfg.k_frac < K_FRAC_UNRANKED {
        return;
    }
    if !path.replace('\\', "/").split('/').any(|c| c == "results") {
        return;
    }
    panic!(
        "refusing to write `{path}`: mode=balanced with k_frac={} takes the top 100% of the \
         frontier, so the ranking runs and changes nothing. That is the uniform-mode control, \
         not a render of the ranked frontier. Set k_frac < 1 (default {K_FRAC_RANKED}), or pass \
         the allow flag if the unranked run is deliberately the control.",
        cfg.k_frac
    );
}

/// **Refuse to write a tree under `results/` from any integration kernel but the production one.**
///
/// The structural answer to "is the refinement mechanism running on the fixed integrator". The
/// scheduler integrates through whatever `EnsembleCfg` it is handed, and twelve harnesses took
/// the integrator from the default while two pinned `Az` deliberately; a setting correct where
/// it was born and copied into a tree-writing harness would produce a corpus of the wrong
/// physics wearing the right filenames, which is what `refine_flagged: false` did for six days.
/// This checks the four kernel knobs that moved between the superseded corpus and the current
/// one — `integrator`, `step_limit`, `dtau_mode`, `clamp_final_step` — and nothing else:
/// `refine_flagged` is a legitimate named argument, `t_max` and `eta` are experiment axes.
///
/// Same shape as [`assert_not_uniform_in_disguise`]: a configuration that silently reproduces
/// the old behaviour needs a guard, not a convention. Call it from every harness that writes a
/// `.prnq` or a tree render under `results/`.
pub fn assert_production_kernel(ens: &EnsembleCfg, path: &str) {
    if !path.replace('\\', "/").split('/').any(|c| c == "results") {
        return;
    }
    let p = EnsembleCfg::production();
    let mut bad: Vec<String> = Vec::new();
    if ens.integrator != p.integrator {
        bad.push(format!("integrator={:?} (production {:?})", ens.integrator, p.integrator));
    }
    if format!("{:?}", ens.step_limit) != format!("{:?}", p.step_limit)
        || ens.step_limit_f != p.step_limit_f
    {
        bad.push(format!(
            "step_limit={:?} f={} (production {:?} f={})",
            ens.step_limit, ens.step_limit_f, p.step_limit, p.step_limit_f
        ));
    }
    if format!("{:?}", ens.dtau_mode) != format!("{:?}", p.dtau_mode) {
        bad.push(format!("dtau_mode={:?} (production {:?})", ens.dtau_mode, p.dtau_mode));
    }
    if ens.clamp_final_step != p.clamp_final_step {
        bad.push(format!(
            "clamp_final_step={} (production {})",
            ens.clamp_final_step, p.clamp_final_step
        ));
    }
    if !bad.is_empty() {
        panic!(
            "refusing to write `{path}`: the integration kernel is not production's -- {}. A \
             tree built on another kernel is the superseded corpus over again; write it to a \
             scratch root, or fix the config.",
            bad.join("; ")
        );
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SchedCfg {
    /// `N`, samples per quad axis. The quality/compute driver: `N²(E+1)` trajectories per quad.
    pub n: usize,
    /// Levels split unconditionally before any decision. Level 0 has no parent and therefore no
    /// `alpha`, so a bootstrap is unavoidable — but *how many* levels is a policy choice, not a
    /// physical one, hence a parameter.
    pub bootstrap_levels: u32,
    /// Cap on **quads computed**, not trajectories. At `N=8`, `E+1=8` one quad is 512 trajectories
    /// and ~47 ms.
    pub budget: usize,
    /// Absolute depth cap. **Superseded by the camera's relative-depth predicate** and kept
    /// only so PR #11's runs reproduce exactly. `None` runs to the budget, which is what
    /// §4 question 1 required *before there was a screen floor*.
    pub max_level: Option<u32>,
    /// The view. `None` reproduces PR #11 — the criterion **minus its principal stop
    /// condition**, which is the regime the real system never enters.
    ///
    /// The camera is read in exactly one place, [`Camera::veto`], which cannot return
    /// `Decision::Split`. Complexity stays the sole trigger.
    pub camera: Option<Camera>,
    pub tau_display: f64,
    /// **The stationarity stop** (Phase 2b), `Policy::Tolerance` only. An unresolved quad that is
    /// white at the footprint scale, whose class mixture matches its parent's and whose quadrants
    /// agree, and whose shape spread did not fall from the parent, is a homogeneous sea: going
    /// deeper re-samples it and resolves nothing. Such a quad is `Decision::Stationary` -- a keep
    /// for now, re-tested at every playhead. Off, it refines to the floor.
    ///
    /// **Off by default, by measurement.** It passes its analytic tests (a white sea stops at the
    /// bootstrap, a filament through it refines), and on the real fields at `eps = 0.01` it fires
    /// 0, 1 and 34 times on `near-field`, `deep interior` and `preset_shape_h1` -- and on the sea
    /// chart the 34 saved 5.5% of the quads while leaving 0.94% of *resolvable* pixels
    /// unresolved (resolvable-form error 0.0094 against 0.0002 with the stop off, 1.60x the
    /// optimum against 1.11x). The real mixing region is a coherent sponge at the footprint
    /// scale, not white noise, and coherence reads it as structure, which it is. The arms are
    /// computed and dumped on every quad regardless, for the `c_stat`/`delta_mix`/`eps` sweep.
    pub stationary: bool,
    /// Coherence below which a quad reads as white: lag-1 neighbour correlation of the nominal
    /// shape and class fields, `max` of the two arms. The white-noise floor at `N = 8` has sd
    /// about 0.09 over 112 neighbour pairs, so 0.3 is about three sigma.
    pub c_stat: f64,
    /// Total-variation distance below which two class mixtures read as the same mixture -- the
    /// quad against its parent, and each quadrant against the quad. Sixteen footprints per
    /// quadrant at `N = 8` put the sampling noise near 0.15 for a three-class mixture.
    pub delta_mix: f64,
    /// The ranked-frontier fraction used by the **post-horizon rounds** of a live descent. The
    /// `k_frac` throttle exists to spend a budget well while the playhead moves; once it has
    /// stopped there is nothing to defer for, and at `k_frac = 0.25` near-field took seventeen
    /// rounds at the horizon where four would do. `1.0` takes the whole want-list per round.
    pub k_frac_post: f64,
    /// **Merging, in the live descent.** A parent whose four children are all leaves that did not
    /// split this boundary is merged back when it has become resolved or its split shows no gain
    /// (`alpha_area < alpha_lo`): the children are released and the parent is the leaf again.
    /// Off is the running-union control -- the tree only grows.
    pub merge: bool,
    /// **The agreement arm**: the area exponent and the noise stop read the unresolved area that
    /// is structure (`QuadReduction::structured_weight`) rather than all of it, so a sea is
    /// uninteresting and a filament through it is not. Off, every unresolved footprint counts.
    pub agreement: bool,
    /// **The dimension floor.** Whether the no-gain test (`no_gain`: `alpha_area` and
    /// `alpha_spread_set` both under `alpha_lo`) floors a quad and merges its parent. Off, only
    /// the noise stop floors -- an unresolved quad with no structured footprint at all -- so a
    /// sea is stopped by its lack of neighbour agreement and a fat fractal, whose footprints
    /// agree with their neighbours at every scale, is refined like any other structure.
    /// Measured on `config_stability`, whose mixing region has box dimension about 1.94 over
    /// the measurable levels: the dimension floor at any `alpha_lo` in 0.05-0.3 stops 200-260
    /// boxes there and the tree lands above uniform. `alpha_lo = 0` disables both floors.
    pub dim_floor: bool,
    /// How a footprint is called hot for the **shape** statistics.
    ///
    /// Separate from `tau_display`, which still drives the split gate and the absolute mask.
    /// The absolute mask is not replaced: `frac_hot` is identically constant under any quantile
    /// rule, and `frac_hot_between` is the best criterion measured here. See
    /// [`crate::spatial::HotRule`].
    pub hot_rule: HotRule,
    /// Split above this exponent, floor below `alpha_lo`. Between them: keep.
    pub alpha_hi: f64,
    /// Under `Policy::Tolerance`: **the floor on the area exponent**. A split whose children hold
    /// no less structured unresolved area than the coarse end did, to within `2^-alpha_lo`, did
    /// not pay, and its children are `Floor` (re-tested every boundary). `alpha_area = 2 - d` for
    /// an unresolved set of box dimension `d`, so `0.2` refines where the set is thinner than
    /// `d = 1.8` and floors where it is fatter. **The default is `0.005`, a noise margin**: a split
    /// is floored only where it resolved nothing. Measured on six charts, the dimension rung 0.2
    /// floors a fat fractal (`config_stability`, `d ~ 1.94`, 16% of whose area resolves at level
    /// 6) exactly as it floors a sea, costs 12% of that chart's resolvable pixels and lands the
    /// tree above uniform; 0.005 keeps the sea chart's saving (39% against 44%), halves the
    /// other costs and puts every chart at or under uniform. No rung separates a sea from a
    /// sponge that thins only below the sampled levels -- that is a bet on depth not bought, by
    /// construction. `0.0` allows full depth everywhere, the configuration the user must opt
    /// into, because the alternative is uniform depth on every sea; `dim_floor = false` keeps
    /// the noise stop alone. Under the legacy policies it thresholds the spread exponent, as
    /// before, and their pins carry `0.2` explicitly.
    pub alpha_lo: f64,
    /// Floor above this sibling range, under [`Policy::Sibling`].
    pub sib_tau: f64,
    pub policy: Policy,
    pub order: Order,
    pub agg: Agg,
    /// Which signal the split decision reads. Every criterion is computed and dumped whatever
    /// this is set to, so criteria can be compared offline without re-integrating.
    pub criterion: Criterion,
    /// Whether the spatial-structure term enters the signal, and how. §2.2.
    pub structure: StructureMode,
    /// Uniform or balanced. See [`Mode`].
    pub mode: Mode,
    /// **Balanced mode's budget**: the fraction of the split-eligible frontier refined per round.
    ///
    /// Defaults to [`K_FRAC_RANKED`]. [`K_FRAC_UNRANKED`] (`1.0`) takes the top 100% of the
    /// frontier, so the ranking runs and changes nothing — it is the **uniform-mode control**,
    /// and it was the silent default through PR #21. Every dump in `results/charts`,
    /// `results/criterion` and `results/vertical` carries it, which is why the criterion work
    /// has no committed *after* outside `results/sweep`.
    ///
    /// It has **no recommended value in the sense of a tuned one** — the default is the value the
    /// sweep supports, and the sweep is the result. Picking one because a tree looked right is
    /// the constant-tuning defect in its most tempting form, and this is the knob most exposed to
    /// it. See [`assert_not_uniform_in_disguise`].
    pub k_frac: f64,
    /// **Camera bias in the PRIORITY** (§4.3): rank on the visible part of a quad, not the whole
    /// quad. A quad half off-screen has structure the viewer cannot see, and ranking on it spends
    /// budget on nothing.
    ///
    /// `None` disables it, which is every run before this. The value is the viewport margin in
    /// quad-widths — §4.3's honest baseline, against which any prediction model must justify
    /// itself. **The camera enters `priority` here and `veto` never; a `Quad` gains no camera
    /// field.**
    pub camera_bias: Option<f64>,
    /// Enforce the **2:1 balance constraint** — no two adjacent leaves more than one level
    /// apart, or the adaptive render has cracks.
    ///
    /// Off by default so every prior run reproduces byte for byte. It is a *rendering*
    /// requirement, separate from the neighbour-contrast idea that shares the same lookup, and
    /// the splits it forces are marked [`Decision::BalanceForced`] so the share of the budget
    /// spent on geometry rather than physics is countable.
    pub balance: bool,
    /// The chart every quad decodes through. One tree, one chart.
    pub chart: Chart,
    /// Retain each quad's `N²` footprints for the adaptive render. Not a cache — the run is
    /// over when `descend` returns, and nothing is reused across runs.
    pub keep_pixels: bool,
    pub seed: u64,
}

impl Default for SchedCfg {
    fn default() -> Self {
        Self {
            n: 8,
            bootstrap_levels: 2,
            budget: 2000,
            max_level: None,
            camera: None,
            tau_display: 1e-2,
            stationary: false,
            k_frac_post: 1.0,
            merge: true,
            agreement: true,
            dim_floor: true,
            c_stat: 0.3,
            delta_mix: 0.25,
            hot_rule: HotRule::Quantile(0.5),
            alpha_hi: 0.5,
            alpha_lo: 0.005,
            sib_tau: 0.5,
            // The enum's `#[default]`, so the struct and the enum cannot disagree on it again:
            // they did, and every harness built on `..Default::default()` ran the legacy policy
            // while the enum said `Tolerance`.
            policy: Policy::default(),
            order: Order::Spread,
            agg: Agg::Median,
            criterion: Criterion::Within,
            structure: StructureMode::Off,
            mode: Mode::Balanced,
            k_frac: K_FRAC_RANKED,
            camera_bias: None,
            balance: false,
            chart: Chart::BodyPlane,
            keep_pixels: false,
            seed: 0,
        }
    }
}

/// What the descent did, beyond the tree itself.
#[derive(Clone, Debug, Default)]
pub struct SchedStats {
    pub iterations: u32,
    pub quads_computed: usize,
    pub leaves_per_iteration: Vec<usize>,
    pub budget_exhausted: bool,
    pub wall_seconds: f64,
    /// Per-node footprints, kept only when [`SchedCfg::keep_pixels`] is set. **The adaptive
    /// render needs the samples, not the reductions** — a level-3 leaf's `N²` samples are what
    /// it rasterises across its own screen footprint. Indexed by node; empty for nodes not
    /// computed. Off by default, because at 4096 leaves this is ~100 MB.
    pub pixels: Vec<Vec<PixelOut>>,
    /// Footprints integrated, and the share of them duplicated at shared sibling edges. The
    /// duplication is `1/N` of a quad and is a *known cost*, reported rather than fixed: keeping
    /// `Slice` shared with the uniform kernel is worth more than the saving.
    pub footprints: usize,
    /// Quads split to satisfy 2:1 rather than because the criterion asked. **Reported**: if this
    /// is a large share of `quads_computed`, the budget went on geometry rather than physics,
    /// and that is a fact about the run rather than a detail of it.
    pub balance_forced: usize,
    /// The growth curve of a **live** descent, one point per recorded boundary. Empty for a
    /// static descent.
    pub live: Vec<LivePoint>,
    /// The leaf set at each recorded boundary of a live descent, so a test can assert the tree
    /// only grows. Empty for a static descent.
    pub live_leaves: Vec<Vec<usize>>,
    /// **What a late split costs.** Substeps the children of every split requested at boundary
    /// `j > 0` would have spent bringing themselves from `t = 0` to `t_j` — estimated as their
    /// total substeps scaled by `t_j / t_max`, because the march does not record per-boundary
    /// step counts. Zero for a static descent, whose every quad is requested at the start.
    pub catchup_substeps: u64,
    /// Children released by merges over a live descent.
    pub merged: usize,
    /// The most quads resident at once over a live descent, and the count at the end.
    pub resident_peak: usize,
    pub resident_final: usize,
}

/// One boundary of a live descent.
#[derive(Clone, Copy, Debug, Default)]
pub struct LivePoint {
    pub j: usize,
    pub t: f64,
    pub computed: usize,
    pub leaves: usize,
    pub split: usize,
    pub keep: usize,
    pub stationary: usize,
    pub deferred: usize,
    pub merged: usize,
    pub resident: usize,
    /// Splits the 2:1 pass forced this boundary, which the criterion did not ask for.
    ///
    /// Separate from `split` on purpose: §4.4 wants the geometry share of the budget countable,
    /// and it cannot be read off the stop-reason breakdown — [`Decision::BalanceForced`] is set on
    /// the quad being *split*, which immediately gains children, so it is never a leaf and
    /// `QuadTree::stop_breakdown` (which walks leaves) can never report it.
    pub balance_forced: usize,
}

/// A footprint sampler: what fills one footprint of a slice. The integrator in production;
/// an analytic field under test (`crate::testing`), so a policy can be checked against a tree
/// whose right shape is known in advance.
pub type Sampler<'a> = &'a (dyn Fn(&crate::grid::Slice, usize) -> PixelOut + Sync);

/// [`compute_quad`] with the sampler injected. **The one place footprints are produced for the
/// tree**; `descend` and `descend_with` both come through here.
fn compute_quad_with(
    tree: &QuadTree,
    i: usize,
    n: usize,
    tau: f64,
    hot_rule: HotRule,
    t_max: f64,
    sampler: Sampler<'_>,
) -> (QuadReduction, Vec<PixelOut>) {
    let slice = tree.nodes[i].slice(n, tree.body, tree.chart);
    let px: Vec<PixelOut> = (0..slice.npix()).into_par_iter().map(|k| sampler(&slice, k)).collect();
    // Distinctness before divergence: N^2 decodes, no integration, and it is the only test
    // that separates a collapsed decode from a genuinely uniform region.
    let ics: Vec<crate::physics::Cart<f64>> = (0..slice.npix()).map(|k| slice.nominal::<f64>(k)).collect();
    let mut red = reduce(&px, n, tau, hot_rule, t_max);
    red.n_distinct_ic = crate::decode::distinct(&ics) as u32;
    (red, px)
}

/// Max over footprints, with a non-finite footprint treated as a **measurement outcome**.
///
/// §4.3's no-discard rule at the aggregation layer. The earlier form was
/// `.filter(|x| x.is_finite()).fold(0.0, f64::max)`, which discarded exactly the footprints the
/// two statistics using it exist to flag — and the `pixel.rs` no-discard fix made that strictly
/// **worse**, not better. Before it, a budget-truncated pixel contributed a finite,
/// healthy-looking drift and at least reached the `max`; after it that pixel is `+inf`, the
/// filter dropped it entirely, and a quad whose footprints were *all* undetermined folded to
/// `0.0` — reading as perfectly clean. *A statistic can report maximum confidence precisely when
/// it is least informed*, at a second site and caused by the repair of the first.
///
/// **`NaN` and `+inf` are not the same and are not collapsed**, following `band_of`'s
/// convention — `NaN` to the bottom, `+inf` to the top:
///
/// - `+inf` — a copy went non-finite, so [`crate::ensemble::stats::max_dev`] returned an
///   infinite deviation *by design*, or `pixel.rs` marked the pixel unusable. The quad is
///   undetermined and propagates that.
/// - `NaN` — `error_ratio` is `0/0` because `sigma_E(0) == 0`: a collapsed decode, or a
///   configuration family where the statistic is structurally undefined (released from rest,
///   every `L_z` is zero — the reason there is no `L_z` version of `error_ratio` at all). Not
///   evidence of damage. It neither sets the result nor is silently counted as clean.
/// - nothing finite and no `+inf` — `NaN`, **never `0.0`**. A quad nothing is known about must
///   not report the best value the scale admits.
///
/// Note the asymmetry is real and not tidiness: `f64::max` already ignores `NaN`, so dropping
/// the filter alone would have looked correct while still folding an all-`NaN` quad to `0.0`.
pub fn max_no_discard(it: impl Iterator<Item = f64>) -> f64 {
    let mut worst = f64::NEG_INFINITY;
    let mut any = false;
    for x in it {
        if x.is_nan() {
            continue;
        }
        any = true;
        if x > worst {
            worst = x;
        }
    }
    if any {
        worst
    } else {
        f64::NAN
    }
}

/// Can the within-arm signal be read from this footprint at all?
///
/// **Two causes, both counted, and the second is the one that fires.**
///
/// 1. `ensemble_spread` is not a number. `spread_shape` over a copy with a non-finite shape
///    vector — a triple collision, or a diverged trajectory. Real, and at the shipped step
///    control it is rare: measured across `deep interior`, `near-field` and `far`, a quad with
///    every footprint budget-exhausted has **zero** footprints failing this test.
/// 2. **`ensemble_spread` swallowed a `NaN`.** It is `sp_shape.max(sp_event)`, and Rust's
///    `f64::max` **ignores `NaN`** — so a footprint whose shape spread is undetermined reports its
///    *event* spread as an ordinary number. Measured on `deep interior` under the pre-fix kernel:
///    **11 footprints carry a `NaN` `spread_shape` and all 11 report a finite `ensemble_spread`.**
///    A triple collision reaches this with every copy still flagged usable, since `shape_vec` is
///    `NaN` at `I = 0` while the state stays finite.
/// 3. **A copy was flagged unusable and the spread was computed anyway.** `PixelOut::n_nonfinite`
///    counts copies the driver marked `!finite` — budget exhaustion included — and by the
///    standing *never discard an ensemble copy* rule the shape spread is taken over all `E+1`
///    regardless. When the copy stopped early its shape vector is a perfectly finite number that
///    is simply not the number the statistic claims, so the quad reports an ordinary spread over
///    a sample it does not have.
///
/// The third is what makes `Decision::Undetermined` reachable at all; a predicate written on the
/// first alone is dead code wearing a guard's name. **The second is tested explicitly rather than
/// left to the third to cover**: on this corpus all 11 of those footprints also carried an
/// unusable copy, so the guard caught them by coincidence — and coincidence is not coverage.
///
/// The `f64::max` swallowing is a defect in `pixel.rs` and is **not repaired there**: propagating
/// the `NaN` would change `ensemble_spread` itself, which moves every tree and every render, and
/// it wants its own attribution rather than riding along with this one.
///
/// **Deliberately strict: one unusable copy of `E+1` marks the footprint.** The alternative is a
/// fraction with a threshold in it, and a fixed threshold picked without measurement is what this
/// project keeps having to withdraw. It costs nothing where the integration succeeds — all three
/// regions above read `n_nonfinite = 0` at production settings, so no production tree moves.
pub fn footprint_undetermined(p: &PixelOut) -> bool {
    !p.ensemble_spread.is_finite() || !p.spread_shape.is_finite() || p.n_nonfinite > 0
}

/// Reduce `N x N` footprints to one quad number per field.
///
/// `tau` is needed here and not only at decision time because the §3.1/§3.2 signals are
/// **counts and shapes of the hot set**, which have no meaning without a threshold. That makes
/// `tau` an input to the *measurement*, not only to the *decision* — a real widening of what
/// `tau` does, and worth saying out loud given the vertical slice promoted it to the dominant
/// knob under the screen floor.
pub fn reduce(px: &[PixelOut], n: usize, tau: f64, hot_rule: HotRule, t_max: f64) -> QuadReduction {
    let finite = |x: &f64| x.is_finite();
    let mut sp: Vec<f64> = px.iter().map(|p| p.ensemble_spread).filter(finite).collect();
    let mut sh: Vec<f64> = px.iter().map(|p| p.spread_shape).filter(finite).collect();
    let mut ev: Vec<f64> = px.iter().map(|p| p.spread_event).filter(finite).collect();
    // How much of the quad the quantiles below are actually speaking for. Counted here, beside
    // the filter, rather than recomputed downstream from a second pass that could drift out of
    // agreement with it.
    let n_undetermined = px.iter().filter(|p| footprint_undetermined(p)).count() as u32;
    // The tolerance arm. `ensemble_spread > tau` fires on any class disagreement as well as on
    // a shape spread past the tolerance (one dissenting copy of eight reads 0.143), and an
    // unreadable footprint is unresolved by definition. Counts, never quantiles.
    let unresolved = |p: &PixelOut| footprint_undetermined(p) || p.ensemble_spread > tau;
    let n_unresolved = px.iter().filter(|p| unresolved(p)).count() as u32;
    let unresolved_weight: f64 = px
        .iter()
        .enumerate()
        .filter(|(_, p)| unresolved(p))
        .map(|(i, _)| {
            let (ix, iy) = (i % n, i / n);
            let wx = if ix == 0 || ix + 1 == n { 0.5 } else { 1.0 };
            let wy = if iy == 0 || iy + 1 == n { 0.5 } else { 1.0 };
            wx * wy
        })
        .sum();
    let n_unresolved_event_only = px
        .iter()
        .filter(|p| !footprint_undetermined(p) && !(p.spread_shape > tau) && p.spread_event > 0.0)
        .count() as u32;
    let spread_max = max_no_discard(px.iter().map(|p| p.ensemble_spread));
    let (coh_shape, coh_class, class_hist, mix_tv_quadrants) = stationarity_arms(px, n);
    // **Structure per footprint**: unresolved, with a 4-neighbour of the same class whose
    // nominal shape agrees to within `STRUCTURE_AGREE`. Then the weights, whole and per
    // quadrant; a sample on the midline of an odd grid straddles two quadrants and gives each a
    // half.
    let idx = |jx: usize, jy: usize| jy * n + jx;
    let agrees = |a: &PixelOut, b: &PixelOut| -> bool {
        if a.event_class != b.event_class {
            return false;
        }
        let d = (0..3).map(|c| (a.shape_vec[c] - b.shape_vec[c]).powi(2)).sum::<f64>().sqrt() * 0.5;
        d.is_finite() && d <= STRUCTURE_AGREE
    };
    let mut structured_weight = 0.0;
    let mut unresolved_quadrant = [0f32; 4];
    let mut structured_quadrant = [0f32; 4];
    let mut unresolved_edge = [0f32; 4];
    let mut structured_edge = [0f32; 4];
    for (i, p) in px.iter().enumerate() {
        if !unresolved(p) {
            continue;
        }
        let (ix, iy) = (i % n, i / n);
        // Two agreeing neighbours of eight. One of four let a sea footprint through on the
        // chance cap of a random neighbour, a few percent per quad, and those diluted the
        // edge share below; two of eight is a few in ten thousand, and a filament's interior
        // footprints still have the two along it.
        let mut agreeing = 0usize;
        for dy in -1i64..=1 {
            for dx in -1i64..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let (jx, jy) = (ix as i64 + dx, iy as i64 + dy);
                if jx < 0 || jy < 0 || jx >= n as i64 || jy >= n as i64 {
                    continue;
                }
                if agrees(p, &px[idx(jx as usize, jy as usize)]) {
                    agreeing += 1;
                }
            }
        }
        let structured = agreeing >= 2;
        let wx = if ix == 0 || ix + 1 == n { 0.5 } else { 1.0 };
        let wy = if iy == 0 || iy + 1 == n { 0.5 } else { 1.0 };
        if structured {
            structured_weight += wx * wy;
        }
        for (e, on) in [(0, ix == 0), (1, ix + 1 == n), (2, iy == 0), (3, iy + 1 == n)] {
            if on {
                unresolved_edge[e] += (wx * wy) as f32;
                if structured {
                    structured_edge[e] += (wx * wy) as f32;
                }
            }
        }
        let side = |j: usize| -> Vec<(usize, f64)> {
            if n % 2 == 1 && j == n / 2 {
                vec![(0, 0.5), (1, 0.5)]
            } else if 2 * j < n {
                vec![(0, if j == 0 { 0.5 } else { 1.0 })]
            } else {
                vec![(1, if j + 1 == n { 0.5 } else { 1.0 })]
            }
        };
        for (qx, sx) in side(ix) {
            for (qy, sy) in side(iy) {
                unresolved_quadrant[qx + 2 * qy] += (sx * sy) as f32;
                if structured {
                    structured_quadrant[qx + 2 * qy] += (sx * sy) as f32;
                }
            }
        }
    }
    let nfin = sp.len().max(1) as f64;
    let mean = sp.iter().sum::<f64>() / nfin;
    let Between {
        shape: b_shape,
        event: b_event,
        matched: b_matched,
        pooled,
        lay_within,
        lay_between,
        lay_rel_within,
        lay_rel_between,
        grad_within,
        grad_between,
    } = between(px, n, tau, hot_rule);
    let (term_frac, esc_frac, grad) = termination_gradient(px, n, t_max);
    QuadReduction {
        spread_mean: mean,
        spread_median: quantile(&mut sp.clone(), 0.5),
        spread_p90: quantile(&mut sp, 0.9),
        spread_shape_median: quantile(&mut sh, 0.5),
        spread_event_median: quantile(&mut ev, 0.5),
        error_ratio_max: max_no_discard(px.iter().map(|p| p.error_ratio)),
        worst_energy_drift: max_no_discard(px.iter().map(|p| p.energy_drift_max)),
        n_nonfinite: px.iter().map(|p| p.n_nonfinite as u32).sum(),
        n_footprints: px.len() as u32,
        n_undetermined,

        between_shape: b_shape,
        between_event: b_event,
        between_spread: b_shape.max(b_event),
        between_matched: b_matched,
        within_pooled: pooled,

        layout_within: lay_within,
        layout_between: lay_between,
        frac_above_tau_within: lay_within.frac_hot(n),
        frac_above_tau_between: lay_between.frac_hot(n),

        layout_rel_within: lay_rel_within,
        layout_rel_between: lay_rel_between,
        grad_rms_within: grad_within,
        grad_rms_between: grad_between,

        terminated_fraction: term_frac,
        escape_fraction: esc_frac,
        t_end_gradient: grad,
        total_substeps: px.iter().map(|p| p.total_substeps as u64).sum(),
        // Overwritten by `compute_quad`, which has the slice. Defaulting to the footprint count
        // keeps a hand-built reduction from reading as collapsed.
        n_distinct_ic: px.len() as u32,

        running_max_divergence_median: quantile(
            &mut px.iter().map(|p| p.running_max_divergence).filter(finite).collect(),
            0.5,
        ),
        divergence_trend_median: quantile(
            &mut px.iter().map(|p| p.divergence_trend).filter(finite).collect(),
            0.5,
        ),
        // A footprint that never crossed is a measurement outcome, not missing data: it counts
        // in the denominator. Only footprints whose accumulators were never computed at all
        // are excluded, and then the fraction is NaN rather than 0.
        frac_diverged: if px.iter().all(|p| p.running_max_divergence.is_nan()) {
            f64::NAN
        } else {
            px.iter().filter(|p| p.first_divergence_t.is_finite()).count() as f64
                / px.len().max(1) as f64
        },
        first_divergence_median: quantile(
            &mut px.iter().map(|p| p.first_divergence_t).filter(finite).collect(),
            0.5,
        ),

        n_unresolved,
        n_unresolved_undetermined: n_undetermined,
        n_unresolved_event_only,
        unresolved_weight,
        spread_max,
        max_excess: spread_max - tau,

        coh_shape,
        coh_class,
        class_hist,
        mix_tv_quadrants,
        structured_weight,
        unresolved_quadrant,
        structured_quadrant,
        unresolved_edge,
        structured_edge,
        mix_tv_parent: f64::NAN,
    }
}

/// Pearson correlation of `(a, b)` pairs; `NaN` when either side does not vary.
fn pearson(pairs: &[(f64, f64)]) -> f64 {
    let n = pairs.len() as f64;
    if pairs.len() < 3 {
        return f64::NAN;
    }
    let (ma, mb) = (
        pairs.iter().map(|p| p.0).sum::<f64>() / n,
        pairs.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (a, b) in pairs {
        sab += (a - ma) * (b - mb);
        saa += (a - ma) * (a - ma);
        sbb += (b - mb) * (b - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        f64::NAN
    } else {
        sab / (saa * sbb).sqrt()
    }
}

/// The two coherence arms, the class histogram and the quadrant mixture distance of one quad.
///
/// Coherence is read on the **nominal** fields, copy 0, so it is a statement about the field at
/// the footprint scale and not about the jitter. The class arm is agreement between neighbours
/// above what a random permutation of the same histogram would give: `(agree - expect) /
/// (1 - expect)`, `NaN` on a single class, where `expect = sum p_i^2`.
fn stationarity_arms(px: &[PixelOut], n: usize) -> (f64, f64, [u16; 27], f64) {
    let idx = |jx: usize, jy: usize| jy * n + jx;
    // Shape: lag-1 Pearson per component over horizontal and vertical neighbour pairs, meaned
    // over the components that vary.
    let mut comp = Vec::new();
    for c in 0..3 {
        let mut pairs = Vec::with_capacity(2 * n * n);
        for jy in 0..n {
            for jx in 0..n {
                let a = px[idx(jx, jy)].shape_vec[c];
                if jx + 1 < n {
                    let b = px[idx(jx + 1, jy)].shape_vec[c];
                    if a.is_finite() && b.is_finite() {
                        pairs.push((a, b));
                    }
                }
                if jy + 1 < n {
                    let b = px[idx(jx, jy + 1)].shape_vec[c];
                    if a.is_finite() && b.is_finite() {
                        pairs.push((a, b));
                    }
                }
            }
        }
        let r = pearson(&pairs);
        if r.is_finite() {
            comp.push(r);
        }
    }
    let coh_shape = if comp.is_empty() { f64::NAN } else { comp.iter().sum::<f64>() / comp.len() as f64 };

    // Class: histogram, then neighbour agreement above chance.
    let mut hist = [0u16; 27];
    for p in px {
        let k = crate::output::png::event_class_ordinal(p.event_class).unwrap_or(26);
        hist[k.min(26)] = hist[k.min(26)].saturating_add(1);
    }
    let total: f64 = hist.iter().map(|&c| c as f64).sum::<f64>().max(1.0);
    // **Class-conditional**, and the max over classes: for each class present at least twice,
    // the fraction of its footprints' neighbours that are the same class, above the class's own
    // base rate, normalised to 1. A global agreement-above-chance statistic is nearly blind to
    // a thin filament -- one coherent column of eight moves agreement over a hundred pairs by a
    // few percent -- while the filament's own class clusters at 0.6 against a base rate of an
    // eighth. A sea reads near zero on every class; the max over classes is what a filament
    // violates and a sea does not.
    let mut coh_class = f64::NAN;
    for c in 0..27usize {
        if hist[c] < 2 {
            continue;
        }
        let p_c = hist[c] as f64 / total;
        if p_c >= 1.0 - 1e-12 {
            continue;
        }
        let (mut adj, mut same) = (0usize, 0usize);
        for jy in 0..n {
            for jx in 0..n {
                let a = crate::output::png::event_class_ordinal(px[idx(jx, jy)].event_class).unwrap_or(26).min(26);
                if a != c {
                    continue;
                }
                let mut nb = Vec::with_capacity(4);
                if jx > 0 { nb.push(idx(jx - 1, jy)); }
                if jx + 1 < n { nb.push(idx(jx + 1, jy)); }
                if jy > 0 { nb.push(idx(jx, jy - 1)); }
                if jy + 1 < n { nb.push(idx(jx, jy + 1)); }
                for k in nb {
                    adj += 1;
                    let b = crate::output::png::event_class_ordinal(px[k].event_class).unwrap_or(26).min(26);
                    same += (b == c) as usize;
                }
            }
        }
        if adj == 0 {
            continue;
        }
        let v = (same as f64 / adj as f64 - p_c) / (1.0 - p_c);
        coh_class = if coh_class.is_nan() { v } else { coh_class.max(v) };
    }

    // Quadrants: the MEAN TV distance between a quadrant's mixture and the whole quad's. Not
    // the max: at `N = 8` a quadrant holds sixteen footprints, and the largest of four
    // multinomial deviations on a three-class sea reaches 0.33 by sampling alone -- measured on
    // the synthetic sea, where a pure sea split on its own noise. The mean sits near 0.14
    // there, and a filament that touches one quadrant still lifts it.
    let whole = {
        let mut m = [0.0f64; 27];
        for (k, &c) in hist.iter().enumerate() {
            m[k] = c as f64 / total;
        }
        m
    };
    let half = n / 2;
    let (mut acc, mut nq) = (0.0f64, 0usize);
    if half >= 1 {
        for qy in 0..2 {
            for qx in 0..2 {
                let mut h = [0.0f64; 27];
                let mut cnt = 0.0;
                for jy in qy * half..((qy + 1) * half).min(n) {
                    for jx in qx * half..((qx + 1) * half).min(n) {
                        let k = crate::output::png::event_class_ordinal(px[idx(jx, jy)].event_class).unwrap_or(26);
                        h[k.min(26)] += 1.0;
                        cnt += 1.0;
                    }
                }
                if cnt > 0.0 {
                    for v in h.iter_mut() {
                        *v /= cnt;
                    }
                    acc += QuadReduction::mix_tv(&h, &whole);
                    nq += 1;
                }
            }
        }
    }
    let mean_tv = if nq == 0 { f64::NAN } else { acc / nq as f64 };
    (coh_shape, coh_class, hist, mean_tv)
}

/// **The stationarity test**, `Policy::Tolerance` only, on an unresolved quad above the bootstrap.
///
/// Three arms, all required: white at the footprint scale, the same class mixture at two scales
/// and across its quadrants, and a shape spread that did not fall from the parent. A filament
/// fails the first two (organised, and its quadrants differ); a regular island's ribbon fails the
/// third (`alpha ~ 1`); a fat mixing cross at level 2 fails the first (its arms are coherent at
/// that scale). Only a region white at the footprint scale with the same mixture at two scales
/// passes -- and it says "nothing resolvable at this sampling", not "nothing". A `NaN` arm never
/// passes: an unmeasured arm is not evidence of a sea.
pub fn stationary(q: &crate::quad::Quad, cfg: &SchedCfg) -> bool {
    let r = &q.red;
    let coh = r.coherence();
    let white = coh.is_finite() && coh < cfg.c_stat;
    let same_mix = r.mix_tv_parent.is_finite()
        && r.mix_tv_parent < cfg.delta_mix
        && r.mix_tv_quadrants.is_finite()
        && r.mix_tv_quadrants < cfg.delta_mix;
    let flat = matches!(q.alpha, Some(a) if a.is_finite() && a < cfg.alpha_lo);
    white && same_mix && flat
}

/// The between-footprint arm, the matched-count controls, and the layout fields.
///
/// Split out of [`reduce`] so the ordering is visible: the nominals are collected once, the
/// centroid distances are the per-footprint between-field, and the hot masks are built from
/// that same field rather than from a second pass with a different definition.
struct Between {
    shape: f64,
    event: f64,
    matched: f64,
    pooled: f64,
    lay_within: Layout,
    lay_between: Layout,
    lay_rel_within: Layout,
    lay_rel_between: Layout,
    grad_within: f64,
    grad_between: f64,
}

fn between(px: &[PixelOut], n: usize, tau: f64, hot_rule: HotRule) -> Between {
    // Copy 0 only. The nominal is un-jittered, so between-footprint variation is not
    // contaminated by the within-footprint jitter — which would otherwise put the same
    // perturbation into both arms and make their correlation partly an artefact of sharing an
    // input.
    let nominals: Vec<[f64; 3]> = px.iter().map(|p| p.shape_vec).collect();
    let classes: Vec<u8> = px.iter().map(|p| p.event_class).collect();

    let shape = shape::spread_shape(&nominals);
    let event: f64 = stats::spread_event(&classes);

    // Matched count: the first `E+1` nominals. `E+1` is read from the ensemble the footprints
    // actually carry rather than from cfg, so this cannot silently disagree with them.
    let e1 = px
        .first()
        .map(|p| p.copy_shapes.len().max(p.copy_outcomes.len()))
        .filter(|&k| k > 1)
        .unwrap_or(0);
    let matched = if e1 >= 2 && e1 <= nominals.len() {
        shape::spread_shape(&nominals[..e1])
    } else {
        f64::NAN
    };

    // Pooled: every copy of every footprint. NaN unless the copies were kept — reported as
    // "not measured", never as zero.
    let pooled = if px.iter().all(|p| !p.copy_shapes.is_empty()) && !px.is_empty() {
        let all: Vec<[f64; 3]> = px.iter().flat_map(|p| p.copy_shapes.iter().cloned()).collect();
        shape::spread_shape(&all)
    } else {
        f64::NAN
    };

    // The per-footprint between-field: each nominal's distance from the quad's nominal
    // centroid, halved to match `spread_shape`'s chord convention. This is the only between-arm
    // quantity defined per footprint, and therefore the only one that can carry a mask.
    let cnt = nominals.len().max(1) as f64;
    let mut c = [0.0f64; 3];
    for v in &nominals {
        for k in 0..3 {
            c[k] += v[k];
        }
    }
    for k in 0..3 {
        c[k] /= cnt;
    }
    let dev: Vec<f64> = nominals
        .iter()
        .map(|v| {
            let d = [v[0] - c[0], v[1] - c[1], v[2] - c[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / 2.0
        })
        .collect();

    // **Both hot rules, on both arms.** Non-finite is hot under either -- a footprint that could
    // not be determined is not evidence of calm, and treating it as cold would make the
    // pathological case invisible to exactly the statistic built to find structure.
    //
    // The absolute pair keeps `frac_above_tau_*` and the `frac_hot_*` criteria untouched; the
    // relative pair is what desaturates the shape statistics. Measured on the committed corpus,
    // the absolute mask reads `n_hot == N^2` in 98.8% of the 75,359 `charts/` leaves and 87.1%
    // over all 92,880: one blob covering the whole quad, nearly everywhere, which is no
    // measurement at all. (Two scopes, both stated -- the chart dumps are the saturated end and
    // the zoom ladders the unsaturated one.)
    let field_w: Vec<f64> = px.iter().map(|p| p.ensemble_spread).collect();
    let hot_w = spatial::hot_mask(&field_w, HotRule::AbsTau(tau));
    let hot_b = spatial::hot_mask(&dev, HotRule::AbsTau(tau));
    let rel_w = spatial::hot_mask(&field_w, hot_rule);
    let rel_b = spatial::hot_mask(&dev, hot_rule);

    Between {
        shape,
        event,
        matched,
        pooled,
        lay_within: spatial::layout(&hot_w, n),
        lay_between: spatial::layout(&hot_b, n),
        lay_rel_within: spatial::layout(&rel_w, n),
        lay_rel_between: spatial::layout(&rel_b, n),
        grad_within: spatial::grad_rms(&field_w, n),
        grad_between: spatial::grad_rms(&dev, n),
    }
}

/// Mean absolute spatial gradient of nominal `t_end`, over the **terminated** footprints only.
///
/// Returns `(terminated_fraction, escape_fraction, gradient)`; the gradient is `NaN` when fewer
/// than two adjacent terminated footprints exist. **Not 0** — a zero would be a null that could
/// not have failed, reported as though it were a measurement about the field.
///
/// Terminated means collision **or** escape, because `t_end` is set by whichever came first.
/// The two are counted separately because they are not interchangeable: `deep interior` reads
/// `terminated = 0.99` with the escape arm silent, and calling that an escape fraction would
/// contradict a standing result while appearing to agree with it.
fn termination_gradient(px: &[PixelOut], n: usize, t_max: f64) -> (f64, f64, f64) {
    use crate::outcome::State;
    // `t_end` pinned at the horizon is the censoring case and carries no gradient information.
    let esc: Vec<bool> = px
        .iter()
        .map(|p| !p.censored && p.t_end.is_finite() && p.t_end < t_max * (1.0 - 1e-12))
        .collect();
    let n_term = esc.iter().filter(|&&e| e).count();
    let frac = n_term as f64 / px.len().max(1) as f64;
    let esc_only = px
        .iter()
        .zip(&esc)
        .filter(|(p, &e)| e && State::from_bits(p.state) == Some(State::Escape))
        .count() as f64
        / px.len().max(1) as f64;

    let idx = |jx: usize, jy: usize| jy * n + jx;
    let mut acc = 0.0;
    let mut pairs = 0usize;
    for jy in 0..n {
        for jx in 0..n {
            let a = idx(jx, jy);
            if !esc[a] {
                continue;
            }
            if jx + 1 < n && esc[idx(jx + 1, jy)] {
                acc += (px[a].t_end - px[idx(jx + 1, jy)].t_end).abs();
                pairs += 1;
            }
            if jy + 1 < n && esc[idx(jx, jy + 1)] {
                acc += (px[a].t_end - px[idx(jx, jy + 1)].t_end).abs();
                pairs += 1;
            }
        }
    }
    (frac, esc_only, if pairs == 0 { f64::NAN } else { acc / pairs as f64 })
}

/// **The 2:1 balance pass.** Split any leaf more than one level coarser than a neighbour.
///
/// Uses the existing [`QuadTree::neighbour`], which returns the *same-or-coarser* neighbour by
/// root descent — so a deep leaf's probe lands on the coarse quad that needs splitting, which is
/// exactly the direction this needs.
///
/// Iterated to a fixed point, because splitting a coarse quad creates children that may
/// themselves be two levels under one of *their* neighbours. Bounded by `room` (in quads, not
/// splits) and by a hard iteration cap: an unbounded loop inside a scheduler is not a failure
/// mode worth leaving available, and if the cap is ever reached that is a bug rather than a
/// budget.
///
/// Returns the newly created nodes, which still need computing.
fn balance_pass(tree: &mut QuadTree, iteration: u32, room: usize) -> Vec<usize> {
    let mut made: Vec<usize> = Vec::new();
    for _ in 0..64 {
        let mut want: Vec<usize> = Vec::new();
        for i in tree.leaves() {
            let lv = tree.nodes[i].level;
            for d in Dir::ALL {
                if let Some(j) = tree.neighbour(i, d) {
                    // `j` is a leaf that is at least two levels coarser: it must split.
                    if tree.nodes[j].children.is_none() && tree.nodes[j].level + 1 < lv {
                        want.push(j);
                    }
                }
            }
        }
        want.sort_unstable();
        want.dedup();
        if want.is_empty() {
            return made;
        }
        for i in want {
            if made.len() + 4 > room {
                return made;
            }
            tree.nodes[i].decision = Decision::BalanceForced;
            made.extend_from_slice(&tree.split(i, iteration));
        }
    }
    debug_assert!(false, "balance pass did not reach a fixed point in 64 rounds");
    made
}

/// Run the descent with the integrator. Returns the tree and what it did.
pub fn descend(
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    cfg: &SchedCfg,
    ens: &EnsembleCfg,
    precision: Precision,
) -> (QuadTree, SchedStats) {
    match precision {
        Precision::F32 => {
            let sampler = |sl: &crate::grid::Slice, k: usize| evaluate::<f32>(sl, k, ens);
            descend_with(cx, cy, half, body, cfg, ens.t_max, &sampler)
        }
        Precision::F64 => {
            let sampler = |sl: &crate::grid::Slice, k: usize| evaluate::<f64>(sl, k, ens);
            descend_with(cx, cy, half, body, cfg, ens.t_max, &sampler)
        }
    }
}

/// Run the descent with an injected footprint sampler — the integrator, or an analytic field
/// from `crate::testing` whose right tree is known in advance.
pub fn descend_with(
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    cfg: &SchedCfg,
    t_max: f64,
    sampler: Sampler<'_>,
) -> (QuadTree, SchedStats) {
    let t0 = std::time::Instant::now();
    let mut tree = QuadTree::with_chart(cx, cy, half, cfg.n, body, cfg.chart);
    let mut st = SchedStats::default();
    let mut pending = vec![0usize];
    // **Quads that wanted to split and were outranked**, carried across rounds and re-decided
    // each time. `Policy::Tolerance` only: under it `Keep` means *resolved*, and a dropped
    // unresolved quad wearing that label was the conflation the stop-reason column exists to
    // prevent. The legacy policies drop them as `Keep`, bitwise as before.
    let mut deferred: Vec<usize> = Vec::new();
    let mut iteration = 0u32;

    while !pending.is_empty() {
        // ---- compute ------------------------------------------------------------------
        if st.quads_computed + pending.len() > cfg.budget {
            let room = cfg.budget.saturating_sub(st.quads_computed);
            for &i in pending.iter().skip(room) {
                tree.nodes[i].decision = Decision::BudgetExhausted;
            }
            pending.truncate(room);
            st.budget_exhausted = true;
        }
        if pending.is_empty() {
            break;
        }

        let reds: Vec<(QuadReduction, Vec<PixelOut>)> = pending
            .iter()
            .map(|&i| compute_quad_with(&tree, i, cfg.n, cfg.tau_display, cfg.hot_rule, t_max, sampler))
            .collect();
        for (&i, (r, px)) in pending.iter().zip(reds) {
            tree.nodes[i].red = r;
            tree.nodes[i].iteration = iteration;
            st.footprints += r.n_footprints as usize;
            if cfg.keep_pixels {
                if st.pixels.len() <= i {
                    st.pixels.resize(i + 1, Vec::new());
                }
                st.pixels[i] = px;
            }
        }
        st.quads_computed += pending.len();

        // ---- alpha, against the quad's OWN parent -------------------------------------
        for &i in &pending {
            if let Some(p) = tree.nodes[i].parent {
                let (pr, cr) = (tree.nodes[p].red, tree.nodes[i].red);
                tree.nodes[i].alpha = ratio_log2(pr.spread(cfg.agg), cr.spread(cfg.agg));
                tree.nodes[i].alpha_mean = ratio_log2(pr.spread_mean, cr.spread_mean);
                tree.nodes[i].alpha_p90 = ratio_log2(pr.spread_p90, cr.spread_p90);
                // The two-scale mixture arm: this quad's class mixture against its parent's, at
                // the same playhead -- the reason the parent is kept marching.
                tree.nodes[i].red.mix_tv_parent =
                    QuadReduction::mix_tv(&cr.class_mix(), &pr.class_mix());
            }
        }

        // ---- the reliability signal, once a parent's four children all exist ----------
        let parents: Vec<usize> = {
            let mut v: Vec<usize> = pending.iter().filter_map(|&i| tree.nodes[i].parent).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        for &p in &parents {
            if let Some(kids) = tree.nodes[p].children {
                let a: Vec<f64> = kids
                    .iter()
                    .filter_map(|&k| tree.nodes[k].alpha)
                    .filter(|x| x.is_finite())
                    .collect();
                if a.len() == 4 {
                    let lo = a.iter().cloned().fold(f64::INFINITY, f64::min);
                    let hi = a.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    tree.nodes[p].alpha_sibling_spread = Some(hi - lo);
                }
            }
        }

        // ---- decide -------------------------------------------------------------------
        //
        // The frontier is this round's computed quads plus everything deferred by the ranking
        // in earlier rounds. `decide` is pure on the reduction, so re-deciding costs nothing and
        // nothing is recomputed.
        // ---- the gain of each split, once its four children are computed ----------------
        for &p in &parents {
            tree.nodes[p].alpha_area = area_exponent(&tree, p, cfg);
            tree.nodes[p].alpha_spread_set = spread_exponent(&tree, p, cfg);
        }

        let mut want: Vec<usize> = Vec::new();
        // A deferred quad can have been split by the balance pass since it was deferred; it is
        // no longer a leaf and is not re-decided -- `split` on it would be a second split.
        let frontier: Vec<usize> = pending
            .iter()
            .cloned()
            .chain(deferred.drain(..))
            .filter(|&i| tree.nodes[i].is_leaf())
            .collect();
        for &i in &frontier {
            let d = decide(&tree, i, cfg);
            tree.nodes[i].decision = d;
            if d == Decision::Split {
                want.push(i);
            }
        }

        st.leaves_per_iteration.push(tree.leaves().count());
        iteration += 1;

        // ---- order, then split ---------------------------------------------------------
        order_queue(&mut want, &tree, cfg);

        // **The frontier is ranked and the top `k` gets budget.** This is where the criterion
        // stops being a threshold and becomes a priority: a quad that falls down the ranking is
        // simply not spent on, which is the demotion §3.1 asks for -- no merging, no eviction.
        //
        // `k_frac = 1.0` refines the whole eligible frontier, reproducing the unranked descent
        // exactly, so every prior run stays byte-identical. Deferred quads are **`Keep`, not
        // `BudgetExhausted`**: they were not refused for want of budget, they were outranked, and
        // conflating the two would hide the mechanism inside the stop-reason column that exists
        // to expose it.
        // **Uniform mode is exempt, and leaving it in was a control-destroying bug.**
        //
        // `decide` short-circuits uniform mode to `Split` -- the criterion is *off*, not set
        // permissive -- but the truncation below runs afterwards and would demote the outranked
        // quads to `Keep` anyway. That applies a ranking to the arm whose whole purpose is to
        // have none, and it read as the control mysteriously matching the treatment: measured on
        // `balanced_march`, near-field at `t = 4` gave 40 leaves and depth variance 0.6900 under
        // BOTH arms, to four digits, with the budget never exhausted. Two arms agreeing to the
        // digit is the same tell as three unrelated charts agreeing, one level up.
        if cfg.mode != Mode::Uniform && cfg.k_frac < 1.0 && !want.is_empty() {
            // **The bootstrap is never truncated, and leaving it truncatable was a bug.**
            //
            // Levels below `bootstrap_levels` split *unconditionally*, because level 0 has no
            // parent and therefore no `alpha`: there is no signal to rank them by. Ranking them
            // anyway ranks on a quantity that does not exist yet, and demoting one to `Keep`
            // stops the tree ever reaching the depth where the criterion can decide anything.
            //
            // Measured before the fix: at every `k_frac < 1`, `near-field`, `deep interior` and
            // `preset_shape` returned **byte-identical** leaf counts and depth variances --
            // 16/1/0.000, 10/2/0.160, 7/2/0.245. Three unrelated charts agreeing to the digit is
            // never physics; it was `k_frac` eating the bootstrap, which is chart-independent
            // arithmetic. The `split` column said so too: rows read `split = 2` and `split = 3`
            // where the bootstrap alone requires 5.
            let (boot, rest): (Vec<usize>, Vec<usize>) =
                want.iter().partition(|&&i| tree.nodes[i].level < cfg.bootstrap_levels);
            let k = ((rest.len() as f64 * cfg.k_frac).ceil() as usize).min(rest.len());
            // `min`, not `clamp(1, ..)`: with an empty `rest` there is nothing to take, and
            // forcing one would reintroduce a split the ranking declined.
            for &i in rest.iter().skip(k) {
                if cfg.policy == Policy::Tolerance {
                    tree.nodes[i].decision = Decision::Deferred;
                    deferred.push(i);
                } else {
                    tree.nodes[i].decision = Decision::Keep;
                }
            }
            want = boot.into_iter().chain(rest.into_iter().take(k)).collect();
        }

        let room = cfg.budget.saturating_sub(st.quads_computed) / 4;
        if want.len() > room {
            for &i in want.iter().skip(room) {
                tree.nodes[i].decision = Decision::BudgetExhausted;
            }
            want.truncate(room);
            st.budget_exhausted = true;
        }

        pending = Vec::new();
        for i in want {
            tree.nodes[i].alpha_area = None;
            tree.nodes[i].alpha_spread_set = None;
            tree.nodes[i].no_gain_weight = None;
            pending.extend_from_slice(&tree.split(i, iteration));
        }

        // **After the criterion's splits, not instead of them.** Balance is a rendering
        // requirement and never a reason to refine, so it runs last and can only add.
        if cfg.balance {
            let room = cfg.budget.saturating_sub(st.quads_computed + pending.len());
            let forced = balance_pass(&mut tree, iteration, room);
            st.balance_forced += forced.len();
            pending.extend_from_slice(&forced);
        }
    }

    st.iterations = iteration;
    st.wall_seconds = t0.elapsed().as_secs_f64();
    (tree, st)
}

fn ratio_log2(parent: f64, child: f64) -> Option<f64> {
    if parent > 0.0 && child > 0.0 && parent.is_finite() && child.is_finite() {
        Some((parent / child).log2())
    } else {
        None
    }
}

/// **The area exponent of a quad's split**: `log2(unresolved_area(parent) / sum
/// unresolved_area(children))`, from the reductions as they stand, the edge footprints weighed by
/// the share of their cell inside the quad (`QuadReduction::unresolved_weight`). `None` when the
/// parent has no children or nothing unresolved (there is no gain to measure); `+inf` when the
/// children resolved everything. Judged over two levels, from the grandparent's quadrant to the
/// children, per level. A line reads 1, a sea 0, a boundary of box dimension `d` reads `2 - d`
/// -- its unresolved area scales as `cell^(2-d)`.
///
/// With the whiteness arm on (`cfg.stationary`) the area is `structured_weight`, not the whole
/// unresolved weight: a sea is unresolved over its whole area, so by area alone a filament
/// through it buys no gain, and it is the whiteness arm that tells the sea from the structure
/// crossing it. Per footprint, never per child -- a child that holds a coherent column and sea
/// around it is one eighth structure, not all of it and not none.
pub fn area_exponent(tree: &QuadTree, parent: usize, cfg: &SchedCfg) -> Option<f64> {
    let p = &tree.nodes[parent];
    let kids = p.children?;
    if p.red.n_footprints == 0 || kids.iter().any(|&k| tree.nodes[k].red.n_footprints == 0) {
        return None;
    }
    // **Two levels, not one, and the coarse end is the grandparent's grid restricted to this
    // parent's quadrant.** A footprint cell on a quad boundary is shared by both quads at half
    // weight, and the finer grid below locates the same structure on one side at full weight:
    // per parent that reads as no gain on one side and infinite gain on the other, and only the
    // sum over both is right. The grandparent's cells split cleanly at its midlines, so its
    // quadrant is an honest coarse measurement of this parent's box that skips the parent's
    // own straddle. Measured: a shore in a level-1 edge cell floored its level-2 children.
    // `w * half^2` is the area up to one constant shared by every level.
    let (pa, levels, coarse_unresolved) = match p.parent {
        Some(g) if tree.nodes[g].red.n_footprints > 0 => {
            let gq = &tree.nodes[g];
            let quadrant = (p.cx > gq.cx) as usize + 2 * (p.cy > gq.cy) as usize;
            (
                structured_weight_quadrant(&gq.red, quadrant, false, cfg) * gq.half * gq.half,
                2.0,
                structured_weight_quadrant(&gq.red, quadrant, true, cfg),
            )
        }
        _ => (structured_weight(&p.red, cfg) * p.half * p.half, 1.0, p.red.unresolved_weight),
    };
    // **On the parent's boundary the exponent is not a measurement.** A structure within a
    // quarter of a cell of the boundary is caught by the overhang of the children on both sides
    // at half weight each, while the grandparent's cells, which split cleanly at their midline,
    // assign it wholly to one side: the far side reads a loss and the near side a windfall, and
    // only their sum is right. Where more than half the children's structure lies on the
    // parent's outer edges, decline -- no floor, no merge for no gain. A sea's edge share is its
    // perimeter over its area, well under a half.
    let (mut ca_w, mut edge_w) = (0.0, 0.0);
    for &k in &kids {
        let c = &tree.nodes[k];
        let (e, w) = if cfg.agreement {
            (&c.red.structured_edge, above_chance(c.red.structured_weight))
        } else {
            (&c.red.unresolved_edge, c.red.unresolved_weight)
        };
        if w <= 0.0 {
            continue;
        }
        ca_w += w;
        edge_w += (if c.cx > p.cx { e[1] } else { e[0] }) as f64 + (if c.cy > p.cy { e[3] } else { e[2] }) as f64;
    }
    if ca_w > 0.0 && edge_w > 0.5 * ca_w {
        return None;
    }
    if !(pa > 0.0) {
        // Nothing structured at the coarse end. With the whiteness arm off that means nothing
        // unresolved and there is no gain to judge; with it on, an unresolved region with no
        // coherent population is noise by the arm's own reading, and a split of noise shows no
        // gain by definition -- a sea quad that slips past the stationarity stop on the mixture
        // arm's sampling noise would otherwise run to the cap.
        return if coarse_unresolved > 0.0 { Some(0.0) } else { None };
    }
    let ca: f64 = kids
        .iter()
        .map(|&k| {
            let c = &tree.nodes[k];
            structured_weight(&c.red, cfg) * c.half * c.half
        })
        .sum();
    if !(ca > 0.0) {
        Some(f64::INFINITY)
    } else {
        Some((pa / ca).log2() / levels)
    }
}

/// **Agreement between neighbouring footprints, in chord/2 units of the nominal shape** -- what
/// makes an unresolved footprint structure rather than noise (`QuadReduction::structured_weight`,
/// two agreeing neighbours of eight). 0.1 is about eleven degrees on the shape sphere: the cap a
/// random neighbour lands in by chance is a percent of the sphere, so two of eight is a few in
/// ten thousand. A smooth field reads structured once its neighbours are within it; a steeper
/// one is carried by its spread arm.
pub const STRUCTURE_AGREE: f64 = 0.1;

/// **Below this, a structured weight is chance.** Under the one-of-four rule a 64-footprint sea
/// quad carried about one chance-structured footprint (Poisson, sd about one) and a ratio of two
/// such counts floored a sea at random; under two-of-eight the chance count is far below this
/// and the floor is a guard. A one-column filament weighs 6 interior (its two end rows have one
/// neighbour along it) and 3 on an edge.
pub const STRUCTURE_MIN_WEIGHT: f64 = 2.5;

fn above_chance(w: f64) -> f64 {
    if w >= STRUCTURE_MIN_WEIGHT { w } else { 0.0 }
}

/// The unresolved area that is structure: `structured_weight` under the agreement arm (zero
/// below `STRUCTURE_MIN_WEIGHT`), the whole unresolved weight without it.
pub fn structured_weight(r: &QuadReduction, cfg: &SchedCfg) -> f64 {
    if cfg.agreement { above_chance(r.structured_weight) } else { r.unresolved_weight }
}

/// The same over one quadrant; `all` counts every unresolved footprint, which is the
/// quadrant's unresolved weight.
pub fn structured_weight_quadrant(r: &QuadReduction, quadrant: usize, all: bool, cfg: &SchedCfg) -> f64 {
    if all || !cfg.agreement {
        r.unresolved_quadrant[quadrant] as f64
    } else {
        above_chance(r.structured_quadrant[quadrant] as f64)
    }
}

/// **The spread exponent of a quad's split**: `log2(spread(parent) / mean spread(children))`
/// under `cfg.agg`, one number for the sibling set where `Quad::alpha` is one per child. A
/// smooth field under a tolerance below its cell spread is unresolved over its whole area at
/// every level until the level at which it resolves everywhere at once -- no area is gained,
/// but the spread halves per level and the split is paying. A sea's spread does not fall.
/// `None` unless both spreads are finite and positive.
pub fn spread_exponent(tree: &QuadTree, parent: usize, cfg: &SchedCfg) -> Option<f64> {
    let p = &tree.nodes[parent];
    let kids = p.children?;
    if p.red.n_footprints == 0 || kids.iter().any(|&k| tree.nodes[k].red.n_footprints == 0) {
        return None;
    }
    let ps = p.red.spread(cfg.agg);
    let cs = kids.iter().map(|&k| tree.nodes[k].red.spread(cfg.agg)).sum::<f64>() / 4.0;
    if ps.is_finite() && cs.is_finite() && ps > 0.0 && cs > 0.0 {
        Some((ps / cs).log2())
    } else {
        None
    }
}

/// **A split showed no gain** when neither its unresolved area nor its spread fell by
/// `2^alpha_lo`: noise at this scale. A split whose parent had nothing unresolved is not judged,
/// and at `alpha_lo = 0` -- the opt-in that allows full depth -- nothing is: a structured area
/// can grow with resolution where neighbours agree only once the cells are small enough, and
/// a negative exponent read as no gain would floor exactly the emergence the opt-in exists to
/// follow. Measured: 131 floors on the sea chart at `alpha_lo = 0` before this guard.
pub fn no_gain(q: &crate::quad::Quad, cfg: &SchedCfg) -> bool {
    if cfg.alpha_lo <= 0.0 || !cfg.dim_floor {
        return false;
    }
    match q.alpha_area {
        Some(a) => a < cfg.alpha_lo && q.alpha_spread_set.map_or(true, |s| s < cfg.alpha_lo),
        None => false,
    }
}

/// §3.2. **Guards first, and the default is keep.**
pub fn decide(tree: &QuadTree, i: usize, cfg: &SchedCfg) -> Decision {
    let q = &tree.nodes[i];

    // A numerical stop, distinct from a physical one, and checked before anything else so it can
    // never be mistaken for "the descent did not terminate".
    if q.below_precision_floor(tree.n) {
        return Decision::PrecisionFloor;
    }
    // **Under the tolerance policy a quad that does not want to split is decided before any cap
    // or veto.** The camera floor and the depth cap are stops for a quad that *wanted* to split;
    // a quad whose every footprint is resolved did not, and neither did a stationary sea, and
    // reporting either as `MaxLevel` or `ScreenFloor` attributes the stop to the cap when the
    // criterion had already decided. The stop-reason breakdown exists to say which one fired.
    // The bootstrap, a collapsed decode and an undetermined quad are decided by their own rules
    // below, as before. Legacy policies keep their order bitwise.
    if cfg.policy == Policy::Tolerance
        && cfg.mode != Mode::Uniform
        && q.level >= cfg.bootstrap_levels
        && !q.red.between_collapsed()
        && !q.red.within_undetermined()
    {
        if q.red.n_unresolved == 0 {
            return Decision::Keep;
        }
        if cfg.stationary && stationary(q, cfg) {
            return Decision::Stationary;
        }
        // **Noise at this sampling.** Unresolved, under the agreement arm nothing in the quad is
        // structure, and its spread did not fall from its parent's: a split cannot pay, whatever
        // its parent's split did. A pure-sea quad carried past the stationarity stop by a
        // sibling's gain would otherwise run to the cap. `alpha_lo = 0` is the opt-in that
        // allows full depth, and it disables this stop too.
        if cfg.alpha_lo > 0.0
            && cfg.agreement
            && q.red.unresolved_weight > 0.0
            && structured_weight(&q.red, cfg) <= 0.0
            && !q.alpha.is_some_and(|a| a >= cfg.alpha_lo)
        {
            return Decision::Floor;
        }
        // **The area floor.** The split that made this quad bought no less unresolved area than
        // its parent held, to within `2^-alpha_lo`: the region is noise at this scale, and going
        // deeper re-samples it. A keep for now, re-tested every boundary.
        if let Some(p) = q.parent {
            if no_gain(&tree.nodes[p], cfg) {
                return Decision::Floor;
            }
        }
        // A quad merged back for no gain remembers its own exponents, so it is not re-split
        // into the same four children every boundary -- for as long as its unresolved area
        // stands within a factor of two of where it was when the split was judged. Past that the
        // region has changed (a band collapsing to a filament) and the memory expires.
        // The memory is keyed on the STRUCTURED weight, the quantity the exponent judged: on a
        // sea the unresolved weight never moves while structure can appear from nothing. A merge
        // judged at zero structure expires the moment any appears.
        if q.is_leaf() && no_gain(q, cfg) {
            let stands = q.no_gain_weight.map_or(true, |w0| {
                let w = structured_weight(&q.red, cfg);
                if w0 <= 0.0 { w <= 0.0 } else { w > 0.5 * w0 && w < 2.0 * w0 }
            });
            if stands {
                return Decision::Floor;
            }
        }
    }
    // **The veto.** Evaluated live from (quad, camera) and never stored on the quad: zoom in
    // and the same patch regrows above pixel size and refines with real new samples. It sits
    // ahead of the bootstrap too — an unconditional split past the screen floor would be the
    // same error one level up.
    if let Some(cam) = cfg.camera {
        if let Some(d) = cam.veto(q, tree.n, tree.nodes[0].half) {
            return d;
        }
    }
    if let Some(m) = cfg.max_level {
        if q.level >= m {
            return Decision::MaxLevel;
        }
    }
    // No parent, so no exponent to read. Split blind for the bootstrap levels.
    if q.level < cfg.bootstrap_levels {
        return Decision::Split;
    }

    // Undetermined, not resolved. Placed with the precision floor rather than among the policy
    // branches because it is a property of the samples, not of the signal read from them — a
    // collapsed quad is collapsed under every criterion at once.
    if q.red.between_collapsed() {
        return Decision::Collapsed;
    }

    // **The second way to be undetermined**, and until this landed it had no decision at all.
    // A quad none of whose footprints could be read fell through to the spread gate below and
    // came out `Keep` — *refinement does not pay*, about a patch where nothing integrated. Not
    // by the `NaN` route the defect was first written up as: measured, the spread is an ordinary
    // finite number computed over truncated copies, and in `near-field` it is 5.6x *smaller* than
    // the healthy quad's. See `QuadReduction::within_undetermined`.
    //
    // Tested after `Collapsed` because identical ICs are the cause and divergence is downstream
    // of them; both can hold and the more fundamental label wins. Stopping is right — the
    // trajectories are hard for physical reasons and four smaller quads are four harder ones —
    // but it must be **countable**, which is the whole difference from `Keep`.
    if q.red.within_undetermined() {
        return Decision::Undetermined;
    }

    // **Uniform mode turns the criterion off entirely.** Not "sets a permissive threshold" --
    // off. It is the control for §3.2's depth-variance test and must be able to reach the veto
    // on every branch, or the test is comparing two criteria rather than criterion against none.
    if cfg.mode == Mode::Uniform {
        return Decision::Split;
    }

    // **The tolerance policy.** Resolved means every footprint resolved; anything else splits.
    // No aggregate runs ahead of it: a quad with one hot footprint of 64 has a median below
    // `tau` and would read `Keep` under the gate below, which is *median under-refines thin
    // structure* at full strength.
    if cfg.policy == Policy::Tolerance {
        // Resolved and stationary were returned above, ahead of the caps; what is left wants
        // to split.
        return Decision::Split;
    }

    let spread = q.red.signal_with(cfg.criterion, cfg.agg, cfg.structure);
    if !(spread > cfg.tau_display) {
        return Decision::Keep;
    }

    match cfg.policy {
        Policy::Sibling => {
            // Read the parent's sibling range: the reliability of the exponent this quad was
            // handed. Scattered siblings mean no alpha here is worth acting on.
            let sib = q.parent.and_then(|p| tree.nodes[p].alpha_sibling_spread);
            match sib {
                Some(s) if s > cfg.sib_tau => Decision::Floor,
                _ => alpha_branch(q.alpha, cfg),
            }
        }
        Policy::Alpha => alpha_branch(q.alpha, cfg),
        // Returned above, before any aggregate or exponent is read.
        Policy::Tolerance => unreachable!("the tolerance policy decides before the legacy branch"),
    }
}

fn alpha_branch(alpha: Option<f64>, cfg: &SchedCfg) -> Decision {
    match alpha {
        Some(a) if a >= cfg.alpha_hi => Decision::Split,
        Some(a) if a < cfg.alpha_lo => Decision::Floor,
        // Between the thresholds, or no exponent at all: the default is keep.
        _ => Decision::Keep,
    }
}

/// The quantity the queue orders by.
///
/// **This used to read `red.spread(agg)` and ignore `cfg.criterion` entirely.** Every
/// `--order spread` run in the corpus therefore ordered by the within arm whatever its header
/// said, so those orderings were measured on a different quantity than they claimed. Fixed here,
/// and recorded rather than quietly corrected: it means a prior `order` result compared the
/// *budget-truncation point* under one signal while the header named another.
pub fn priority(tree: &QuadTree, i: usize, cfg: &SchedCfg) -> f64 {
    let q = &tree.nodes[i];
    let v = q.red.signal_with(cfg.criterion, cfg.agg, cfg.structure);
    let v = match cfg.order {
        Order::SpreadArea => v * q.half.powi(2),
        _ => v,
    };
    // **A product of two terms, never either alone** (§4.3). Structure changes only when a quad
    // is recomputed or the zoom changes; relevance changes on every frame the camera moves --
    // which is the split the persistent frontier is built around, and the reason this is
    // computed here rather than stored.
    match (cfg.camera_bias, cfg.camera) {
        (Some(margin), Some(cam)) => v * cam.relevance(q.cx, q.cy, q.half, margin),
        _ => v,
    }
}

fn order_queue(want: &mut [usize], tree: &QuadTree, cfg: &SchedCfg) {
    match cfg.order {
        Order::Shuffled => {
            let mut rng = SplitMix64::new(cfg.seed ^ 0x5EED_C0DE_5EED_C0DE);
            for j in (1..want.len()).rev() {
                let k = (rng.next_u64() % (j as u64 + 1)) as usize;
                want.swap(j, k);
            }
        }
        // Descending priority. NaN sorts LAST rather than blocking: a signal that declines to
        // score must not outrank one that did, and must not stop the ones that did from being
        // ordered among themselves.
        _ => want.sort_by(|&a, &b| {
            let (pa, pb) = (priority(tree, a, cfg), priority(tree, b, cfg));
            match (pa.is_nan(), pb.is_nan()) {
                (true, true) => a.cmp(&b),
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                _ => pb.partial_cmp(&pa).unwrap_or(std::cmp::Ordering::Equal),
            }
        }),
    }
}


// -------------------------------------------------------------------------------------------
// Phase 2c: the live descent.
// -------------------------------------------------------------------------------------------

/// **The footprint as it was known at recorded boundary `j`** — the live view.
///
/// Everything the tolerance policy reads is taken from the live series at `j`: the copies'
/// shape and event spreads, the nominal shape and class. A footprint that had terminated by
/// `t_j` keeps its terminal state; one that had not is `Bounded` at `t_j`, censored, because at
/// `t_j` that is all that was known. Every accumulator that reads past `t_j` is masked to `NaN`
/// rather than left at its terminal value: the running maxima, the first-divergence time, the
/// terminal energy statistics. `n_nonfinite` is left at its terminal value, which is
/// **conservative** — a copy that failed later than `t_j` reads as unusable earlier than it was,
/// so the view splits sooner than a perfectly live one would, never later. `total_substeps` is
/// scaled to `t_j / t_max` as an estimate.
///
/// A decision made on this view cannot read the future, which is the live-playhead contract as
/// a function rather than a thing to remember.
pub fn project_at(p: &PixelOut, j: usize) -> PixelOut {
    use crate::outcome::State;
    assert!(j < p.live_t.len(), "no live series entry {j}: keep_live_series was off, or j is past the end");
    let t_j = p.live_t[j];
    let mut q = p.clone();
    q.spread_shape = p.live_spread_shape[j];
    q.spread_event = p.live_spread_event[j];
    q.ensemble_spread = q.spread_shape.max(q.spread_event);
    q.shape_vec = p.live_shape[j];
    q.event_class = p.live_class[j];
    // `n_nonfinite` is a verdict on the whole march; the live view takes the count known **at**
    // this boundary. Without this a copy that diverges at `t = 12` paints its footprint
    // undetermined in the frame at `t = 0.8`, and both the render and `footprint_undetermined`
    // read it. Empty for a series written before the field existed: fall back to the run's count
    // rather than silently reporting zero, which would read as "nothing is wrong here".
    q.n_nonfinite = p.live_nonfinite.get(j).copied().unwrap_or(p.n_nonfinite);
    let terminated = !p.censored && p.t_end <= t_j * (1.0 + 1e-12);
    if !terminated {
        q.state = State::Bounded as u8;
        q.detail = 0;
        q.outcome = (State::Bounded as u8) << 2;
        q.t_end = t_j;
        q.censored = true;
    }
    q.running_max_divergence = f64::NAN;
    q.divergence_trend = f64::NAN;
    q.first_divergence_t = f64::NAN;
    q.spread_event_max = f64::NAN;
    q.spread_event_latched = f64::NAN;
    q.t_spread_event = f64::NAN;
    q.error_ratio = f64::NAN;
    q.error_ratio_mad = f64::NAN;
    q.energy_drift_max = f64::NAN;
    q.energy_drift_nominal = f64::NAN;
    let frac = if p.t_end > 0.0 { (t_j / p.live_t.last().cloned().unwrap_or(t_j)).clamp(0.0, 1.0) } else { 1.0 };
    q.total_substeps = (p.total_substeps as f64 * frac).round() as u64;
    q.total_force_evals = (p.total_force_evals as f64 * frac).round() as u64;
    q
}

/// Run the **live** descent with the integrator: one march per footprint, the tree grown
/// boundary by boundary from the live series. `ens.keep_live_series` must be on.
pub fn descend_live(
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    cfg: &SchedCfg,
    ens: &EnsembleCfg,
    precision: Precision,
) -> (QuadTree, SchedStats) {
    assert!(ens.keep_live_series, "descend_live needs EnsembleCfg::keep_live_series");
    match precision {
        Precision::F32 => {
            let sampler = |sl: &crate::grid::Slice, k: usize| evaluate::<f32>(sl, k, ens);
            descend_live_with(cx, cy, half, body, cfg, ens.t_max, &sampler)
        }
        Precision::F64 => {
            let sampler = |sl: &crate::grid::Slice, k: usize| evaluate::<f64>(sl, k, ens);
            descend_live_with(cx, cy, half, body, cfg, ens.t_max, &sampler)
        }
    }
}

/// **The live descent**: the tree as a playhead would have built it.
///
/// Every quad is marched to `t_max` once and keeps its live series. The tree is then grown
/// boundary by boundary: at each recorded boundary `j`, every live leaf is re-decided on its
/// footprints **projected to `j`** ([`project_at`]), with its parent's projection beside it for
/// `alpha` and the two-scale mixture arm; a `Split` computes the four children — which join the
/// frontier at the *same* boundary, since they too are marched to `t_max` — and is
/// **irreversible**; `Keep`, `Stationary` and `Deferred` are re-tested at the next boundary;
/// the caps, the precision floor, a collapsed decode and an undetermined quad are terminal for
/// the quad. So the tree can only grow, `Quad::iteration` records the boundary at which a quad
/// was requested, and `SchedStats::catchup_substeps` prices what a late split costs. At the last
/// boundary every projection equals the terminal footprint, so the finished tree's reductions
/// are the same numbers the static descent would carry.
///
/// **The bootstrap is at `j = 0`**: nothing is known before the first boundary, and the
/// unconditional splits below `bootstrap_levels` are by fiat, exactly as in the static descent.
pub fn descend_live_with(
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    cfg: &SchedCfg,
    t_max: f64,
    sampler: Sampler<'_>,
) -> (QuadTree, SchedStats) {
    assert!(cfg.policy == Policy::Tolerance, "the live descent is defined for Policy::Tolerance");
    let t0 = std::time::Instant::now();
    let mut tree = QuadTree::with_chart(cx, cy, half, cfg.n, body, cfg.chart);
    let mut st = SchedStats::default();
    let mut px_of: Vec<Vec<PixelOut>> = Vec::new();

    // Compute a quad at boundary `j`: march, store, record the request time.
    let compute = |tree: &mut QuadTree, st: &mut SchedStats, px_of: &mut Vec<Vec<PixelOut>>, i: usize, j: usize| {
        let (r, px) = compute_quad_with(tree, i, cfg.n, cfg.tau_display, cfg.hot_rule, t_max, sampler);
        assert!(!px.is_empty() && !px[0].live_t.is_empty(), "the sampler produced no live series");
        tree.nodes[i].red = r;
        tree.nodes[i].iteration = j as u32;
        st.footprints += r.n_footprints as usize;
        st.quads_computed += 1;
        if px_of.len() <= i {
            px_of.resize(i + 1, Vec::new());
        }
        px_of[i] = px;
    };

    compute(&mut tree, &mut st, &mut px_of, 0, 0);
    let n_b = px_of[0][0].live_t.len();

    // The bootstrap, by fiat, at j = 0.
    let mut frontier: Vec<usize> = vec![0];
    for _ in 0..cfg.bootstrap_levels {
        let mut next = Vec::new();
        for &i in &frontier {
            if st.quads_computed + 4 > cfg.budget {
                tree.nodes[i].decision = Decision::BudgetExhausted;
                st.budget_exhausted = true;
                continue;
            }
            tree.nodes[i].decision = Decision::Split;
            let kids = tree.split(i, 0);
            for &k in &kids {
                compute(&mut tree, &mut st, &mut px_of, k, 0);
            }
            next.extend_from_slice(&kids);
        }
        frontier = next;
    }

    // Terminal decisions leave the frontier; everything else is re-tested every boundary. **A
    // cap is not terminal.** A leaf stopped by `MaxLevel`, `ScreenFloor` or `MaxRelDepth` wanted
    // to split and could not; its region can still resolve at a later boundary, and since a
    // resolved quad is decided ahead of the caps, re-testing it then reads `Keep` -- which is
    // what lets its parent merge it. Measured on the pulse under `alpha_lo = 0.005`: with caps
    // terminal, 24 parents read zero unresolved footprints at the horizon and still held 96
    // capped children, so the live tree ended at 149 quads where the static tree at the same
    // playhead had 69; the no-gain merges at 0.2 had hidden it by merging those parents earlier.
    let terminal = |d: Decision| {
        matches!(
            d,
            Decision::PrecisionFloor | Decision::Collapsed | Decision::Undetermined | Decision::BudgetExhausted
        )
    };

    // Boundaries `0..n_b`, then **post-horizon rounds** at the last boundary until nothing
    // wants to split: the playhead stops at the horizon, the tree does not. Children requested
    // at boundary `j` can first be decided at `j + 1` (they have to catch up to the playhead),
    // so the tree gains at most one level per boundary during the march; at the horizon it
    // continues in the static regime, and a child requested there costs a full march.
    let mut j = 0usize;
    let mut post = 0usize;
    loop {
        let t_j = px_of[0][0].live_t[j];
        // Project every live leaf and its parent to `j`, so the decision reads the boundary.
        for &i in &frontier {
            let proj: Vec<PixelOut> = px_of[i].iter().map(|p| project_at(p, j)).collect();
            let mut r = reduce(&proj, cfg.n, cfg.tau_display, cfg.hot_rule, t_max);
            r.n_distinct_ic = tree.nodes[i].red.n_distinct_ic;
            tree.nodes[i].red = r;
        }
        for &i in &frontier {
            if let Some(pi) = tree.nodes[i].parent {
                let pproj: Vec<PixelOut> = px_of[pi].iter().map(|p| project_at(p, j)).collect();
                let mut pr = reduce(&pproj, cfg.n, cfg.tau_display, cfg.hot_rule, t_max);
                pr.n_distinct_ic = tree.nodes[pi].red.n_distinct_ic;
                let cr = tree.nodes[i].red;
                tree.nodes[i].alpha = ratio_log2(pr.spread(cfg.agg), cr.spread(cfg.agg));
                tree.nodes[i].alpha_mean = ratio_log2(pr.spread_mean, cr.spread_mean);
                tree.nodes[i].alpha_p90 = ratio_log2(pr.spread_p90, cr.spread_p90);
                tree.nodes[i].red.mix_tv_parent = QuadReduction::mix_tv(&cr.class_mix(), &pr.class_mix());
            }
        }

        // The area exponent of every frontier quad's parent, on the projections at `j`: the
        // parent was projected above, and so were its children -- every child of a frontier
        // quad's parent is a leaf (it was split together) and so is in the frontier or terminal.
        {
            let mut parents: Vec<usize> = frontier.iter().filter_map(|&i| tree.nodes[i].parent).collect();
            parents.sort_unstable();
            parents.dedup();
            for p in parents {
                // The parent as it stands at `j` -- the loop above read it and did not keep it,
                // and the exponent must compare parent and children at the same boundary. The
                // grandparent too: its quadrant is the exponent's coarse end.
                if let Some(g) = tree.nodes[p].parent {
                    let gproj: Vec<PixelOut> = px_of[g].iter().map(|q| project_at(q, j)).collect();
                    let mut gr = reduce(&gproj, cfg.n, cfg.tau_display, cfg.hot_rule, t_max);
                    gr.n_distinct_ic = tree.nodes[g].red.n_distinct_ic;
                    gr.mix_tv_parent = tree.nodes[g].red.mix_tv_parent;
                    tree.nodes[g].red = gr;
                }
                let pproj: Vec<PixelOut> = px_of[p].iter().map(|q| project_at(q, j)).collect();
                let mut pr = reduce(&pproj, cfg.n, cfg.tau_display, cfg.hot_rule, t_max);
                pr.n_distinct_ic = tree.nodes[p].red.n_distinct_ic;
                let pmix = pr.class_mix();
                tree.nodes[p].red = pr;
                if let Some(kids) = tree.nodes[p].children {
                    for k in kids {
                        if !frontier.contains(&k) {
                            let proj: Vec<PixelOut> = px_of[k].iter().map(|q| project_at(q, j)).collect();
                            let mut r = reduce(&proj, cfg.n, cfg.tau_display, cfg.hot_rule, t_max);
                            r.n_distinct_ic = tree.nodes[k].red.n_distinct_ic;
                            r.mix_tv_parent = QuadReduction::mix_tv(&r.class_mix(), &pmix);
                            tree.nodes[k].red = r;
                        }
                    }
                }
                tree.nodes[p].alpha_area = area_exponent(&tree, p, cfg);
                tree.nodes[p].alpha_spread_set = spread_exponent(&tree, p, cfg);
            }
        }

        let mut want: Vec<usize> = Vec::new();
        let mut point = LivePoint { j, t: t_j, ..Default::default() };
        for &i in &frontier {
            let d = decide(&tree, i, cfg);
            tree.nodes[i].decision = d;
            match d {
                Decision::Split => want.push(i),
                Decision::Keep => point.keep += 1,
                Decision::Stationary => point.stationary += 1,
                _ => {}
            }
        }
        order_queue(&mut want, &tree, cfg);
        // While the playhead moves the frontier is throttled by `k_frac`; at the horizon by
        // `k_frac_post`, which is 1.0 by default because there is nothing left to defer for.
        let kf = if post == 0 { cfg.k_frac } else { cfg.k_frac_post };
        if cfg.mode != Mode::Uniform && kf < 1.0 && !want.is_empty() {
            let k = ((want.len() as f64 * kf).ceil() as usize).min(want.len());
            for &i in want.iter().skip(k) {
                tree.nodes[i].decision = Decision::Deferred;
                point.deferred += 1;
            }
            want.truncate(k);
        }
        let room = cfg.budget.saturating_sub(st.quads_computed) / 4;
        if want.len() > room {
            for &i in want.iter().skip(room) {
                tree.nodes[i].decision = Decision::BudgetExhausted;
            }
            want.truncate(room);
            st.budget_exhausted = true;
        }
        point.split = want.len();

        // **The merge pass**, the reverse of the split rule at a later boundary. A parent whose
        // four children are all leaves that did not split this round is merged when it has become
        // resolved (uniform) or when its split has stopped paying (`alpha_area < alpha_lo`: noise).
        // The children are released -- their trajectories go, the parent's were kept marching --
        // and the parent is the leaf again, re-tested next boundary; it may split again later.
        {
            let mut parents: Vec<usize> = frontier.iter().filter_map(|&i| tree.nodes[i].parent).collect();
            parents.sort_unstable();
            parents.dedup();
            let mut rejoin: Vec<usize> = Vec::new();
            for p in parents {
                let Some(kids) = tree.nodes[p].children else { continue };
                // Settled: a leaf that did not split this round, whatever stopped it. A child at a
                // cap is as settled as one that kept -- the parent's own resolution or no-gain is
                // what the merge reads, and the static rule would never have split that parent.
                let all_settled = kids.iter().all(|&k| {
                    let d = tree.nodes[k].decision;
                    tree.nodes[k].is_leaf() && !tree.nodes[k].merged
                        && matches!(
                            d,
                            Decision::Keep | Decision::Floor | Decision::Stationary | Decision::Deferred
                                | Decision::MaxLevel | Decision::ScreenFloor | Decision::MaxRelDepth
                                | Decision::PrecisionFloor | Decision::Collapsed | Decision::Undetermined
                        )
                });
                if !cfg.merge || !all_settled || tree.nodes[p].level < cfg.bootstrap_levels {
                    continue;
                }
                let resolved = tree.nodes[p].red.n_unresolved == 0 && tree.nodes[p].red.n_footprints > 0;
                let stalled = no_gain(&tree.nodes[p], cfg);
                if resolved || stalled {
                    for k in kids {
                        tree.nodes[k].merged = true;
                        tree.nodes[k].decision = Decision::Merged;
                        st.merged += 1;
                        point.merged += 1;
                    }
                    tree.nodes[p].children = None;
                    tree.nodes[p].decision = if resolved { Decision::Keep } else { Decision::Floor };
                    if resolved {
                        tree.nodes[p].alpha_area = None;
                        tree.nodes[p].alpha_spread_set = None;
                        tree.nodes[p].no_gain_weight = None;
                    } else {
                        tree.nodes[p].no_gain_weight = Some(structured_weight(&tree.nodes[p].red, cfg));
                    }
                    rejoin.push(p);
                }
            }
            frontier.retain(|&i| !tree.nodes[i].merged);
            frontier.extend(rejoin);
        }

        let mut next: Vec<usize> = Vec::new();
        for &i in &frontier {
            if !terminal(tree.nodes[i].decision) && tree.nodes[i].decision != Decision::Split {
                next.push(i);
            }
        }
        for &i in &want {
            tree.nodes[i].alpha_area = None;
            tree.nodes[i].alpha_spread_set = None;
            tree.nodes[i].no_gain_weight = None;
            let kids = tree.split(i, j as u32);
            for &k in &kids {
                compute(&mut tree, &mut st, &mut px_of, k, j);
                if j > 0 || post > 0 {
                    let steps: u64 = px_of[k].iter().map(|p| p.total_substeps as u64).sum();
                    st.catchup_substeps += (steps as f64 * (t_j / t_max).clamp(0.0, 1.0)).round() as u64;
                }
            }
            next.extend_from_slice(&kids);
        }
        frontier = next;

        // **The 2:1 balance pass, live.** The static descent has run this since it was written and
        // this one never did, so every tree the live descent produced was unbalanced -- and the
        // live descent is the one closest to the target design. It goes *here*, after the split
        // loop rather than beside `decide`, because a forced child has to be computed **and caught
        // up to the playhead** exactly like a child the criterion asked for; putting it next to the
        // decision would have created quads that no boundary ever marched.
        //
        // It runs after the merge pass, which matters: a merge un-splits a parent and can itself
        // create a violation against a neighbour, and running balance afterwards repairs that in
        // the same round rather than leaving a cracked frame until the next one.
        if cfg.balance {
            let room = cfg.budget.saturating_sub(st.quads_computed);
            let forced = balance_pass(&mut tree, j as u32, room);
            st.balance_forced += forced.len();
            point.balance_forced += forced.len();
            for &k in &forced {
                compute(&mut tree, &mut st, &mut px_of, k, j);
                if j > 0 || post > 0 {
                    let steps: u64 = px_of[k].iter().map(|p| p.total_substeps as u64).sum();
                    st.catchup_substeps += (steps as f64 * (t_j / t_max).clamp(0.0, 1.0)).round() as u64;
                }
            }
            frontier.extend_from_slice(&forced);
        }

        point.j = j + post;
        point.computed = st.quads_computed;
        point.leaves = tree.leaves().count();
        point.resident = tree.resident();
        st.resident_peak = st.resident_peak.max(point.resident);
        st.resident_final = point.resident;
        st.live.push(point);
        st.live_leaves.push(tree.leaves().collect());
        st.leaves_per_iteration.push(point.leaves);
        st.iterations = (j + post + 1) as u32;
        if j + 1 < n_b {
            j += 1;
        } else {
            post += 1;
            let merged_now = st.live.last().map_or(0, |p| p.merged);
            if (want.is_empty() && merged_now == 0) || frontier.is_empty() || post > 64 {
                break;
            }
        }
    }

    if cfg.keep_pixels {
        st.pixels = px_of;
    }
    st.wall_seconds = t0.elapsed().as_secs_f64();
    (tree, st)
}
