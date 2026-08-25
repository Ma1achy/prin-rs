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
use crate::quad::{quantile, Agg, Criterion, Decision, QuadReduction, QuadTree};
use crate::spatial::{self, HotRule, Layout};
use crate::render::Precision;
use crate::rng::SplitMix64;

/// Which signal the floor decision reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Policy {
    /// Threshold on `alpha`'s **value**. Separation between region types is 0.9862 against a
    /// chaotic scatter of 1.1–1.3 — marginal.
    #[default]
    Alpha,
    /// Threshold on `alpha_sibling_spread`, the range of the four children's exponents.
    /// Separation in `alpha`'s **reliability** is 0.001 against 1.2 — three orders. Where the four
    /// scatter, the unreliability *is* the answer, and no trustworthy `alpha` is needed.
    Sibling,
}

impl Policy {
    pub fn name(self) -> &'static str {
        match self {
            Policy::Alpha => "alpha",
            Policy::Sibling => "sibling",
        }
    }
    pub fn parse(s: &str) -> Option<Policy> {
        Some(match s {
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
    /// How a footprint is called hot for the **shape** statistics.
    ///
    /// Separate from `tau_display`, which still drives the split gate and the absolute mask.
    /// The absolute mask is not replaced: `frac_hot` is identically constant under any quantile
    /// rule, and `frac_hot_between` is the best criterion measured here. See
    /// [`crate::spatial::HotRule`].
    pub hot_rule: HotRule,
    /// Split above this exponent, floor below `alpha_lo`. Between them: keep.
    pub alpha_hi: f64,
    pub alpha_lo: f64,
    /// Floor above this sibling range, under [`Policy::Sibling`].
    pub sib_tau: f64,
    pub policy: Policy,
    pub order: Order,
    pub agg: Agg,
    /// Which signal the split decision reads. Every criterion is computed and dumped whatever
    /// this is set to, so criteria can be compared offline without re-integrating.
    pub criterion: Criterion,
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
            hot_rule: HotRule::Quantile(0.5),
            alpha_hi: 0.5,
            alpha_lo: 0.2,
            sib_tau: 0.5,
            policy: Policy::Alpha,
            order: Order::Spread,
            agg: Agg::Median,
            criterion: Criterion::Within,
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
}

/// Compute one quad: `N²` footprints, each an `E+1` ensemble, reduced to one `QuadReduction`.
///
/// **Never pooled from children.** This is the quad's own ensemble at its own cell width.
fn compute_quad<T: crate::Real>(
    tree: &QuadTree,
    i: usize,
    ens: &EnsembleCfg,
    n: usize,
    tau: f64,
    hot_rule: HotRule,
) -> (QuadReduction, Vec<PixelOut>) {
    let slice = tree.nodes[i].slice(n, tree.body, tree.chart);
    let px: Vec<PixelOut> = (0..slice.npix())
        .into_par_iter()
        .map(|k| evaluate::<T>(&slice, k, ens))
        .collect();
    // Distinctness before divergence: N^2 decodes, no integration, and it is the only test
    // that separates a collapsed decode from a genuinely uniform region.
    let ics: Vec<crate::physics::Cart<f64>> = (0..slice.npix()).map(|k| slice.nominal::<f64>(k)).collect();
    let mut red = reduce(&px, n, tau, hot_rule, ens.t_max);
    red.n_distinct_ic = crate::decode::distinct(&ics) as u32;
    (red, px)
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
        error_ratio_max: px
            .iter()
            .map(|p| p.error_ratio)
            .filter(|x| x.is_finite())
            .fold(0.0f64, f64::max),
        worst_energy_drift: px
            .iter()
            .map(|p| p.energy_drift_max)
            .filter(|x| x.is_finite())
            .fold(0.0f64, f64::max),
        n_nonfinite: px.iter().map(|p| p.n_nonfinite as u32).sum(),
        n_footprints: px.len() as u32,

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
    }
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
    // the absolute mask reads `n_hot == N^2` in 98.8% of leaves and `n_components == 1` in
    // 99.5%: one blob covering the whole quad, everywhere, which is no measurement at all.
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

/// Run the descent. Returns the tree and what it did.
pub fn descend(
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    cfg: &SchedCfg,
    ens: &EnsembleCfg,
    precision: Precision,
) -> (QuadTree, SchedStats) {
    let t0 = std::time::Instant::now();
    let mut tree = QuadTree::with_chart(cx, cy, half, cfg.n, body, cfg.chart);
    let mut st = SchedStats::default();
    let mut pending = vec![0usize];
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
            .map(|&i| match precision {
                Precision::F32 => {
                    compute_quad::<f32>(&tree, i, ens, cfg.n, cfg.tau_display, cfg.hot_rule)
                }
                Precision::F64 => {
                    compute_quad::<f64>(&tree, i, ens, cfg.n, cfg.tau_display, cfg.hot_rule)
                }
            })
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
            }
        }

        // ---- the reliability signal, once a parent's four children all exist ----------
        let parents: Vec<usize> = {
            let mut v: Vec<usize> = pending.iter().filter_map(|&i| tree.nodes[i].parent).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        for p in parents {
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
        let mut want: Vec<usize> = Vec::new();
        for &i in &pending {
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
            pending.extend_from_slice(&tree.split(i, iteration));
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

/// §3.2. **Guards first, and the default is keep.**
pub fn decide(tree: &QuadTree, i: usize, cfg: &SchedCfg) -> Decision {
    let q = &tree.nodes[i];

    // A numerical stop, distinct from a physical one, and checked before anything else so it can
    // never be mistaken for "the descent did not terminate".
    if q.below_precision_floor(tree.n) {
        return Decision::PrecisionFloor;
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

    let spread = q.red.signal(cfg.criterion, cfg.agg);
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

fn order_queue(want: &mut [usize], tree: &QuadTree, cfg: &SchedCfg) {
    match cfg.order {
        Order::Spread => want.sort_by(|&a, &b| {
            tree.nodes[b]
                .red
                .spread(cfg.agg)
                .partial_cmp(&tree.nodes[a].red.spread(cfg.agg))
                .unwrap_or(std::cmp::Ordering::Equal)
        }),
        Order::SpreadArea => want.sort_by(|&a, &b| {
            let w = |i: usize| tree.nodes[i].red.spread(cfg.agg) * tree.nodes[i].half.powi(2);
            w(b).partial_cmp(&w(a)).unwrap_or(std::cmp::Ordering::Equal)
        }),
        Order::Shuffled => {
            let mut rng = SplitMix64::new(cfg.seed ^ 0x5EED_C0DE_5EED_C0DE);
            for j in (1..want.len()).rev() {
                let k = (rng.next_u64() % (j as u64 + 1)) as usize;
                want.swap(j, k);
            }
        }
    }
}
