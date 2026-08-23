//! The descent loop — the one thing the uniform kernel was deliberately built without.
//!
//! Every measurement before this ran the criterion on **one split in isolation**. Everything here
//! exists because the remaining questions are dynamic: does the descent terminate, does the floor
//! engage, does a budget get spent well, does per-quad noise cause thrash.
//!
//! **Scope discipline** (SCHEDULER_BRIEF §6): no eviction, no caching, no async, no promotion, no
//! camera, no interaction. A quad is computed once and the tree keeps it — that is the tree holding
//! its own data, not a cache. If any of the others appears here, it is a bug.

use rayon::prelude::*;

use crate::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use crate::quad::{quantile, Agg, Decision, QuadReduction, QuadTree};
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
    /// `None` runs to the budget, which is what §4 question 1 requires.
    pub max_level: Option<u32>,
    pub tau_display: f64,
    /// Split above this exponent, floor below `alpha_lo`. Between them: keep.
    pub alpha_hi: f64,
    pub alpha_lo: f64,
    /// Floor above this sibling range, under [`Policy::Sibling`].
    pub sib_tau: f64,
    pub policy: Policy,
    pub order: Order,
    pub agg: Agg,
    pub seed: u64,
}

impl Default for SchedCfg {
    fn default() -> Self {
        Self {
            n: 8,
            bootstrap_levels: 2,
            budget: 2000,
            max_level: None,
            tau_display: 1e-2,
            alpha_hi: 0.5,
            alpha_lo: 0.2,
            sib_tau: 0.5,
            policy: Policy::Alpha,
            order: Order::Spread,
            agg: Agg::Median,
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
) -> QuadReduction {
    let slice = tree.nodes[i].slice(n, tree.body);
    let px: Vec<PixelOut> = (0..slice.npix())
        .into_par_iter()
        .map(|k| evaluate::<T>(&slice, k, ens))
        .collect();
    reduce(&px)
}

pub fn reduce(px: &[PixelOut]) -> QuadReduction {
    let finite = |x: &f64| x.is_finite();
    let mut sp: Vec<f64> = px.iter().map(|p| p.ensemble_spread).filter(finite).collect();
    let mut sh: Vec<f64> = px.iter().map(|p| p.spread_shape).filter(finite).collect();
    let mut ev: Vec<f64> = px.iter().map(|p| p.spread_event).filter(finite).collect();
    let n = sp.len().max(1) as f64;
    let mean = sp.iter().sum::<f64>() / n;
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
    }
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
    let mut tree = QuadTree::new(cx, cy, half, cfg.n, body);
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

        let reds: Vec<QuadReduction> = pending
            .iter()
            .map(|&i| match precision {
                Precision::F32 => compute_quad::<f32>(&tree, i, ens, cfg.n),
                Precision::F64 => compute_quad::<f64>(&tree, i, ens, cfg.n),
            })
            .collect();
        for (&i, r) in pending.iter().zip(reds) {
            tree.nodes[i].red = r;
            tree.nodes[i].iteration = iteration;
            st.footprints += r.n_footprints as usize;
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
fn decide(tree: &QuadTree, i: usize, cfg: &SchedCfg) -> Decision {
    let q = &tree.nodes[i];

    // A numerical stop, distinct from a physical one, and checked before anything else so it can
    // never be mistaken for "the descent did not terminate".
    if q.below_precision_floor(tree.n) {
        return Decision::PrecisionFloor;
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

    let spread = q.red.spread(cfg.agg);
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
