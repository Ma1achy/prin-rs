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

/// Aggregate `N²` footprint spreads into the three numbers a `QuadReduction` carries.
pub fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[(((v.len() - 1) as f64) * q).round() as usize]
}
