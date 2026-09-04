//! **`PRNF` — the per-frame record. A frame time is meaningless without the work it did.**
//!
//! That sentence is the telemetry contract's own rule and it is the whole design: every record
//! carries what was asked of it, or the log is uninterpretable six weeks later.
//!
//! **Profiling is first-class, not a debug mode.** *"Instrumentation that can be compiled out will
//! be, and the numbers it produces will not match the build people actually run."* This project
//! has the precedent on both sides: `ab_floored` and `ab_min` were computed on every march and read
//! by nothing — a sticky bit nothing reads is indistinguishable from one that never fires — and the
//! inverse failure, a counter that exists only in a debug build, measures the debug build. So the
//! measurement is unconditional and only the *reporting* is a choice. A timestamp per stage against
//! a 16.7 ms budget is nanoseconds.
//!
//! **`stage_ms` is the load-bearing field.** Without it you know a frame was slow and not **which
//! resource ran out**, and different devices bottleneck differently: integrated graphics is
//! bandwidth-bound, a large discrete part at small workloads may never approach its throughput at
//! all. A tier derived on the assumption that everything is compute-bound is wrong on both.
//!
//! **A stage this build does not have is `NaN`, never `0.0`.** There is no GPU and no window here,
//! so `upload` and `present` are not measured — and a zero there reads as *instant* where the truth
//! is *absent*, which is the same conflation as an empty mask reading as "no structure found", or
//! `Absent` pooled with `Evicted` one module over.
//!
//! **One format for the batch path and the interactive one.** A `prin` render emits the same record
//! with `camera_delta = 0` and no present stage: one parser, one set of percentile code. And the
//! headline is `frac_frames_over_41.7ms` **measured during motion** — static frames may take longer
//! without anyone minding, a dropped frame mid-pan is immediately visible, and `camera_delta > 0`
//! is the free discriminator.
//!
//! Same self-describing shape as `PRNQ` and `PRQC`: magic, version, a length-prefixed text header
//! naming every parameter, then a record count, a field count, and one `f64` per field per frame.

use std::io::{self, Write};

pub const MAGIC: &[u8; 4] = b"PRNF";
pub const VERSION: u32 = 1;

/// 60 fps. The **goal**.
pub const BUDGET_60: f64 = 16.7;
/// 24 fps. The **hard floor, not a goal** — below it the thing is not interactive. A tier that
/// lands here has failed to find a good setting, not succeeded at finding an acceptable one.
pub const BUDGET_24: f64 = 41.7;

/// Per-stage wall time, in milliseconds.
#[derive(Clone, Copy, Debug)]
pub struct StageMs {
    pub integrate: f64,
    pub reduce: f64,
    pub decide: f64,
    pub evict: f64,
    pub colour: f64,
    /// `NaN` here: there is no GPU in this build, and `0.0` would read as *instant*.
    pub upload: f64,
    /// `NaN` here: there is no window in this build.
    pub present: f64,
}

impl Default for StageMs {
    fn default() -> Self {
        StageMs {
            integrate: 0.0,
            reduce: 0.0,
            decide: 0.0,
            evict: 0.0,
            colour: 0.0,
            upload: f64::NAN,
            present: f64::NAN,
        }
    }
}

/// One frame.
#[derive(Clone, Copy, Debug, Default)]
pub struct FrameRecord {
    pub frame: u64,
    /// **Measured, not budgeted.** The quota is a count; this is what that count cost on this
    /// machine. Nothing in the scheduler branches on it — `the_quota_is_not_a_deadline` pins that.
    pub frame_ms: f64,
    pub stage_ms: StageMs,
    pub quads_computed: usize,
    /// Leaves the frame used without computing: already decided, payload resident.
    pub quads_reused: usize,
    /// Computed **again** because the payload had been released. **Never pooled with
    /// `quads_computed`**: one is progress and the other is the price of the cap, and a session
    /// reporting their sum reports a budget being spent well while it thrashes.
    pub quads_recomputed: usize,
    pub quads_evicted: usize,
    pub resident_payloads: usize,
    /// The cap could not be met without evicting into the exempt class. Reported, never silently
    /// ignored — a cap that cannot be met and does not say so reads as a cap that is being met.
    pub cap_unsatisfiable: bool,
    /// Footprints integrated. `trajectories = samples * (E+1)` is carried too, because "samples"
    /// has meant three different things in this project.
    pub samples: usize,
    pub trajectories: u64,
    /// The honest cost measure. Steps are not comparable across steppers; substeps are.
    pub substeps_total: u64,
    /// How far the playhead moved. **Separate from `camera`**, because a measurement taken during
    /// a gesture has two causes and one that cannot attribute a churn spike to either is not
    /// evidence.
    pub playhead_dt: f64,
    /// Camera motion. `0` for a static frame, and the discriminator for the headline statistic.
    pub camera_delta: f64,
    pub camera_zoom_octaves: f64,
    /// Cursor motion — the **third** delta, so foveation churn is separable from the other two.
    pub cursor_delta: f64,
    pub tree_depth_max: u32,
    pub leaf_count: usize,
    /// Share of the raster painted by an **ancestor** rather than the leaf owning the pixel — the
    /// direct measure of the tree lagging the camera. If it stays high after motion stops, the
    /// scheduler is not converging, and nothing could see that before.
    pub ancestor_fill_fraction: f64,
    /// Painted by an ancestor because the leaf was **evicted**, split out from the above for the
    /// same reason `Evicted` is not `Absent`.
    pub evicted_fill_fraction: f64,
    pub background_fraction: f64,
    pub balance_forced_fraction: f64,
    pub rounds: u32,
    /// Which quota bound, as a code — never inferred from the counts.
    pub quota_hit: u8,
    /// `NaN` on frames where the frontier audit did not run. **Never `1.0` by default**: a check
    /// that did not run must not report a pass.
    pub frontier_agrees: f64,
}

impl FrameRecord {
    /// Was this frame taken during motion? The free discriminator for the headline statistic.
    pub fn moving(&self) -> bool {
        self.camera_delta != 0.0 || self.camera_zoom_octaves != 0.0
    }

    pub fn over_floor(&self) -> bool {
        self.frame_ms > BUDGET_24
    }

    pub fn over_goal(&self) -> bool {
        self.frame_ms > BUDGET_60
    }
}

pub const FIELDS: &[&str] = &[
    "frame", "frame_ms",
    "integrate_ms", "reduce_ms", "decide_ms", "evict_ms", "colour_ms", "upload_ms", "present_ms",
    "quads_computed", "quads_reused", "quads_recomputed", "quads_evicted", "resident_payloads",
    "cap_unsatisfiable", "samples", "trajectories", "substeps_total",
    "playhead_dt", "camera_delta", "camera_zoom_octaves", "cursor_delta",
    "tree_depth_max", "leaf_count",
    "ancestor_fill_fraction", "evicted_fill_fraction", "background_fraction",
    "balance_forced_fraction", "rounds", "quota_hit", "frontier_agrees",
];

/// Compile-time tie between the field names and the record, so adding one to either alone breaks
/// the build rather than shifting every column at read time. `tree.rs` already establishes the
/// precedent, after a hand-kept array length went stale.
pub const N_FIELDS: usize = FIELDS.len();

pub fn record(r: &FrameRecord) -> [f64; N_FIELDS] {
    [
        r.frame as f64,
        r.frame_ms,
        r.stage_ms.integrate,
        r.stage_ms.reduce,
        r.stage_ms.decide,
        r.stage_ms.evict,
        r.stage_ms.colour,
        r.stage_ms.upload,
        r.stage_ms.present,
        r.quads_computed as f64,
        r.quads_reused as f64,
        r.quads_recomputed as f64,
        r.quads_evicted as f64,
        r.resident_payloads as f64,
        if r.cap_unsatisfiable { 1.0 } else { 0.0 },
        r.samples as f64,
        r.trajectories as f64,
        r.substeps_total as f64,
        r.playhead_dt,
        r.camera_delta,
        r.camera_zoom_octaves,
        r.cursor_delta,
        r.tree_depth_max as f64,
        r.leaf_count as f64,
        r.ancestor_fill_fraction,
        r.evicted_fill_fraction,
        r.background_fraction,
        r.balance_forced_fraction,
        r.rounds as f64,
        r.quota_hit as f64,
        r.frontier_agrees,
    ]
}

/// The percentiles a frame log is read by. **Never a mean** — a mean hides stutter, and stutter is
/// what makes a thing feel broken.
pub struct Percentiles {
    pub p50: f64,
    pub p95: f64,
    pub p99: f64,
    pub max: f64,
    pub n: usize,
}

pub fn percentiles(v: &[f64]) -> Percentiles {
    if v.is_empty() {
        return Percentiles { p50: f64::NAN, p95: f64::NAN, p99: f64::NAN, max: f64::NAN, n: 0 };
    }
    let mut s: Vec<f64> = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q = |p: f64| {
        let i = ((s.len() - 1) as f64 * p).round() as usize;
        s[i]
    };
    Percentiles { p50: q(0.50), p95: q(0.95), p99: q(0.99), max: s[s.len() - 1], n: s.len() }
}

/// **The headline: the share of frames over the 41.7 ms floor, measured DURING MOTION.**
///
/// Static frames may take longer without anyone minding; a dropped frame mid-pan is immediately
/// visible. `NaN` when no frame moved — which is a statement that the question was not asked, not
/// a score of zero, and a harness reporting `0.0` there would be claiming a pass it never earned.
pub fn frac_over_floor_moving(log: &[FrameRecord]) -> f64 {
    let moving: Vec<&FrameRecord> = log.iter().filter(|r| r.moving()).collect();
    if moving.is_empty() {
        return f64::NAN;
    }
    moving.iter().filter(|r| r.over_floor()).count() as f64 / moving.len() as f64
}

pub fn write<W: Write>(w: &mut W, log: &[FrameRecord], header: &str) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;
    let h = format!("{header}\nfields={}\n", FIELDS.join(","));
    let hb = h.as_bytes();
    w.write_all(&(hb.len() as u32).to_le_bytes())?;
    w.write_all(hb)?;
    w.write_all(&(log.len() as u64).to_le_bytes())?;
    w.write_all(&(N_FIELDS as u32).to_le_bytes())?;
    for r in log {
        for v in record(r) {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    Ok(())
}
