//! **Phase 1 of the refinement rebuild: the physics-space ground truth, and the memory number.**
//!
//! Every `error(B)` curve before this scored a tree by OKLab distance under the shipping
//! colouring, whose lightness is auto-ranged to each region's own p1–p99 — so a smooth region's
//! `1e-8` residual was stretched to full contrast and counted as error at every depth, and
//! breadth-first came out near-optimal by construction of the metric. Here a tree is scored on
//! the **payload**: the nominal copy's point on the shape sphere and its event class, against a
//! fixed-scale tolerance `eps` (chord/2, the units of `spread_shape`). Same number under every
//! colouring anyone might choose, which is the user's constraint.
//!
//! Three forms of the same metric are built from one integration and printed side by side:
//! `indicator` (the fraction of the image unresolved — the headline), `hinge` (how far past the
//! tolerance), and `resolvable` (unresolved pixels a finer grid could resolve). **The gap
//! between `indicator` and `resolvable` is the sea cost**, and `sea_fraction(eps)` is printed
//! before any curve: the fraction of the frame that no depth resolves at this `eps`.
//!
//! Four modes:
//!
//! - `build <target> [levels] [n] [eps] [t_max] [root]` — integrate one complete tree for a
//!   Burrau region (`grid::REGIONS`), `config_stability`, or a gallery chart
//!   (`grid::gallery_cases`), write the v2 footprint file and the quad cache, and print every
//!   table for every form plus the shipping OKLab metric as the comparison arm. Signal rankings
//!   and a static descent under the shipped policy (mapped onto the cache) are scored here.
//! - `replay <file.fcache> [eps] [root]` — zero trajectories: the floor (`uniform`), the
//!   reference (`greedy_lookahead_1`), the random band and the ceiling (`dp_optimal`) from a
//!   committed footprint file. A v1 file is replayed under the **outcome** arm only and says so.
//! - `live <file.fcache> [policy] [eps] [k_frac] [root] [stationary] [tau] [alpha_lo] [agreement] [dim_floor]` — one static descent under
//!   the named policy, integrated fresh, mapped onto the cache and scored against the ceiling:
//!   the memory ratio at matched quality, `Policy::Tolerance` against `Policy::Alpha`.
//! - `march <file.fcache> [eps] [k_frac] [root] [stationary] [live_stride] [alpha_lo] [merge] [agreement] [dim_floor]` — the **live**
//!   descent (`scheduler::descend_live`): the tree grown boundary by boundary from one march per
//!   quad, with its growth curve and catch-up cost, scored the same way.
//!
//! The memory number is the `B_needed` table: the budget each strategy needs to bring the
//! unresolved fraction to a target, and `B_needed(row) / B_needed(dp)` — the memory ratio at
//! matched quality. `steps` (total substeps at that budget) is printed beside every `B`.

use std::collections::HashMap;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::logln;
use prin_rs::metric::{self, Cache, ClassArm, Colouring, Key, Metric, PayloadForm, Point, Rank};
use prin_rs::output::colour::Scalar;
use prin_rs::output::Log;
use prin_rs::quad::{Agg, Criterion, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg, SchedStats};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const SHIPPING: Colouring = Colouring::Bivariate(Scalar::ShapeSpread);
const TARGETS: [f64; 4] = [0.05, 0.02, 0.01, 0.005];
const EPS_LADDER: [f64; 5] = [2e-3, 5e-3, 1e-2, 2e-2, 5e-2];

struct Target {
    name: String,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
}

fn target(name: &str) -> Option<Target> {
    if let Some(&(n, cx, cy, body)) = grid::REGIONS.iter().find(|r| r.0 == name) {
        return Some(Target { name: n.into(), chart: Chart::BodyPlane, cx, cy, half: 0.05, body });
    }
    // Every named slice with its own window, from the one table in `grid`.
    if let Some((chart, cx, cy, half)) = grid::named_slice(name) {
        return Some(Target { name: name.into(), chart, cx, cy, half, body: 0 });
    }
    grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == name)
        .map(|(n, chart, cx, cy, half)| Target { name: n.into(), chart, cx, cy, half, body: 0 })
}

/// `4k+1` plus every complete-level count, so a rung lands ON a complete tree as well as past it.
fn budgets(full: usize) -> Vec<usize> {
    let mut b = vec![5usize];
    while *b.last().unwrap() * 2 < full {
        b.push(b.last().unwrap() * 2 + 1);
    }
    b.push(full);
    let (mut lvl, mut acc) = (1usize, 1usize);
    while acc < full {
        b.push(acc);
        lvl *= 4;
        acc += lvl;
    }
    b.sort_unstable();
    b.dedup();
    b
}

fn fmt_b(b: Option<Point>) -> String {
    match b {
        // A replayed cache has no reductions and so no substeps; `--` rather than a zero that
        // reads as "free".
        Some(p) if p.cost == 0 => format!("{:>7} {:>10}", p.budget, "(steps --)"),
        Some(p) => format!("{:>7} ({:.2e})", p.budget, p.cost as f64),
        None => format!("{:>7} {:>10}", "--", ""),
    }
}

/// The whole table set for one cache.
fn tables(log: &Log, c: &Cache, ranks: &[Rank], live: &[(String, Vec<Key>, usize, u64, String)]) {
    let full = c.quads.len();
    let max_splits = (full - 1) / 4;
    let budgets = budgets(full);
    let deepest: Vec<Key> = {
        let w = 1u32 << c.levels;
        (0..w).flat_map(|iy| (0..w).map(move |ix| (c.levels, ix, iy))).collect()
    };
    logln!(log, "=== metric {} ===", c.metric.name());
    logln!(log, "  error(root) = {:.5}   error(full) = {:.1e}   (full = matches this sampling, not correct)",
           c.error_of(&[(0, 0, 0)]), c.error_of(&deepest));
    if c.metric.is_payload() {
        let row: Vec<String> = EPS_LADDER.iter().map(|&e| format!("{e:.0e}:{:.4}", c.sea_fraction(e))).collect();
        logln!(log, "  sea_fraction(eps) -- the frame no depth resolves: {}", row.join("  "));
    }

    let dp = c.dp_optimal(max_splits);
    logln!(log, "  dp_optimal: {} splits in {:.2}s, prefix-min binds at {} budgets",
           dp.max_splits, dp.elapsed_s, dp.prefix_min_binds.len());

    let runs: Vec<(String, Vec<Point>)> =
        ranks.iter().map(|&r| (r.name(), metric::replay(c, r, full))).collect();
    // The DP bound, asserted at scale -- an assert that has never run is a gap that survives.
    let mut worst = 0.0f64;
    for (name, pts) in &runs {
        for p in pts {
            if dp.covers(p.budget) {
                let m = p.error - dp.at_budget(p.budget);
                if m < worst {
                    worst = m;
                    if m < -1e-9 {
                        panic!("{name} beats dp_optimal at B={} by {m:e}: the ceiling is not a ceiling", p.budget);
                    }
                }
            }
        }
    }
    logln!(log, "  dp bound holds over {} rankings; worst margin {worst:+.1e}", runs.len());

    // ---- error(B) ----
    let hdr: Vec<String> = budgets.iter().map(|b| format!("{b:>8}")).collect();
    logln!(log, "  {:<26} {}", "error(B)", hdr.join(" "));
    let dp_row: Vec<String> = budgets
        .iter()
        .map(|&b| if dp.covers(b) { format!("{:>8.5}", dp.at_budget(b)) } else { format!("{:>8}", "--") })
        .collect();
    logln!(log, "  {:<26} {}   CEILING", "dp_optimal", dp_row.join(" "));
    let uni = runs.iter().find(|r| r.0 == "uniform").map(|r| metric::curve_at(&r.1, &budgets));
    let mut rand: Vec<Vec<f64>> = Vec::new();
    for (name, pts) in &runs {
        let cv = metric::curve_at(pts, &budgets);
        if name.starts_with("random") {
            rand.push(cv);
            continue;
        }
        let row: Vec<String> = cv.iter().map(|e| format!("{e:>8.5}")).collect();
        logln!(log, "  {:<26} {}", name, row.join(" "));
    }
    if !rand.is_empty() {
        let lo: Vec<String> = (0..budgets.len())
            .map(|j| format!("{:>8.5}", rand.iter().map(|r| r[j]).fold(f64::INFINITY, f64::min)))
            .collect();
        let hi: Vec<String> = (0..budgets.len())
            .map(|j| format!("{:>8.5}", rand.iter().map(|r| r[j]).fold(f64::NEG_INFINITY, f64::max)))
            .collect();
        logln!(log, "  {:<26} {}   FLOOR band, {} seeds", "random lo", lo.join(" "), rand.len());
        logln!(log, "  {:<26} {}", "random hi", hi.join(" "));
    }

    // ---- captured = (uniform - row)/(uniform - dp), and the headroom it is normalised by ----
    if let Some(uni) = &uni {
        let denom: Vec<f64> = budgets
            .iter()
            .zip(uni.iter())
            .map(|(&b, &u)| if dp.covers(b) { u - dp.at_budget(b) } else { f64::NAN })
            .collect();
        let hr: Vec<String> = denom.iter().map(|d| format!("{d:>8.1e}")).collect();
        logln!(log, "  {:<26} {}   the room", "headroom uni-dp", hr.join(" "));
        let hre: Vec<String> = denom
            .iter()
            .zip(uni.iter())
            .map(|(d, u)| format!("{:>8.1e}", d / u.abs().max(1e-300)))
            .collect();
        logln!(log, "  {:<26} {}   the room as a share of uniform's error", "headroom/err", hre.join(" "));
        for (name, pts) in &runs {
            if name == "uniform" || name.starts_with("random") {
                continue;
            }
            let cv = metric::curve_at(pts, &budgets);
            let row: Vec<String> = (0..budgets.len())
                .map(|j| {
                    let d = denom[j];
                    if !d.is_finite() || d.abs() <= 1e-6 * uni[j].abs().max(1e-30) {
                        format!("{:>8}", "--")
                    } else {
                        format!("{:>8.4}", (uni[j] - cv[j]) / d)
                    }
                })
                .collect();
            logln!(log, "  {:<26} {}", format!("captured {name}"), row.join(" "));
        }
    }

    // ---- B_needed: the memory number ----
    let th: Vec<String> = TARGETS.iter().map(|t| format!("{:>18}", format!("<= {t}"))).collect();
    logln!(log, "  {:<26} {}", "B_needed (steps)", th.join(" "));
    let dp_b: Vec<String> = TARGETS
        .iter()
        .map(|&t| match dp.budget_needed(t) {
            Some(b) => format!("{:>18}", b),
            None => format!("{:>18}", "--"),
        })
        .collect();
    logln!(log, "  {:<26} {}   CEILING", "dp_optimal", dp_b.join(" "));
    for (name, pts) in &runs {
        if name.starts_with("random") {
            continue;
        }
        let row: Vec<String> = TARGETS.iter().map(|&t| fmt_b(Cache::budget_needed(pts, t))).collect();
        logln!(log, "  {:<26} {}", name, row.join(" "));
    }
    for (name, leaves, b, steps, stop) in live {
        let e = c.error_of(leaves);
        let need = dp.budget_needed(e);
        logln!(log, "  live tree `{name}`: B={b} steps={steps:.2e} error={e:.5} stop[{stop}]  dp needs B={} for that error -> memory ratio {}",
               need.map(|x| x.to_string()).unwrap_or("--".into()),
               need.map(|x| format!("{:.2}x", *b as f64 / x as f64)).unwrap_or("--".into()));
    }
    logln!(log);
}

/// A real descent on the cache's box, mapped onto its keys. Leaves that do not map are counted
/// and the row is refused rather than scored over a hole.
fn live_tree(
    c: &Cache,
    t: &Target,
    ens: &EnsembleCfg,
    res: usize,
    label: &str,
) -> Option<(String, Vec<Key>, usize, u64, String)> {
    let cfg = SchedCfg {
        n: c.n,
        camera: Some(Camera::framing(t.cx, t.cy, t.half, res)),
        max_level: Some(c.levels),
        chart: t.chart,
        keep_pixels: false,
        ..Default::default()
    };
    let (tree, st): (QuadTree, SchedStats) =
        scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, ens, Precision::F64);
    let leaves: Vec<usize> = tree.leaves().collect();
    let mut keys = Vec::with_capacity(leaves.len());
    let mut lost = 0usize;
    for &i in &leaves {
        let q = &tree.nodes[i];
        match c.key_of(q.cx, q.cy, q.level) {
            Some(k) => keys.push(k),
            None => lost += 1,
        }
    }
    if lost > 0 {
        eprintln!("  live tree `{label}`: {lost} of {} leaves do not map onto the cache; row refused", leaves.len());
        return None;
    }
    let steps: u64 = tree.nodes.iter().map(|q| q.red.total_substeps as u64).sum();
    Some((
        format!("{label} policy={:?} tau={:e} k_frac={}", cfg.policy, cfg.tau_display, cfg.k_frac),
        keys,
        st.quads_computed,
        steps,
        tree.stop_breakdown(),
    ))
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_default();
    match mode.as_str() {
        "build" => build(),
        "replay" => replay(),
        "live" => live(),
        "march" => march(),
        _ => {
            eprintln!("usage: payload_metric build <target> [levels=6] [n=8] [eps=0.01] [t_max=13] [root=results]");
            eprintln!("       payload_metric replay <file.fcache> [eps=0.01] [root=results]");
            eprintln!("       payload_metric live <file.fcache> [policy=tolerance|alpha] [eps=0.01] [k_frac=0.25] [root=results] [stationary=0] [tau=eps] [alpha_lo=0.2] [agreement=1] [dim_floor=1] [ttl=-1] [camera=1]");
            eprintln!("       payload_metric march <file.fcache> [eps=0.01] [k_frac=0.25] [root=results] [stationary=0] [live_stride=4] [alpha_lo=0.2] [merge=1] [agreement=1] [dim_floor=1] [no_gain_ttl=-1]");
            std::process::exit(2);
        }
    }
}

/// **The live descent, scored against a committed cache.** The tree as a playhead would have
/// built it: grown boundary by boundary from one march per quad, decisions that can only add.
/// Prints the growth curve, the catch-up cost of late splits, and the final tree's error against
/// the same ceiling the static descents are scored against — so `mem_all(live) / mem_all(static)`
/// is the price of monotonicity.
fn march() {
    let file: String = std::env::args().nth(2).expect("a .fcache path");
    let eps: f64 = arg(3, 0.01);
    let k_frac: f64 = arg(4, prin_rs::scheduler::K_FRAC_RANKED);
    let root: String = std::env::args().nth(5).unwrap_or_else(|| "results".into());
    let stationary: bool = std::env::args().nth(6).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().stationary);
    let live_stride: usize = arg(7, 4);
    // **The area floor and merging.** `alpha_lo` is the floor on the area exponent (0 allows
    // full depth -- the opt-in); `merge` releases the children of a parent that has become
    // resolved or whose split shows no gain.
    let alpha_lo: f64 = arg(8, SchedCfg::default().alpha_lo);
    let merge: bool = std::env::args().nth(9).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().merge);
    // The agreement arm: off is the floor-only control, where a sea and a filament through it
    // read alike.
    let agreement: bool = std::env::args().nth(10).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().agreement);
    // The dimension floor: off leaves only the noise stop, so a sea floors and a fat fractal
    // does not. Measured against on, it says which of the two is doing the flooring.
    let dim_floor: bool = std::env::args().nth(11).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().dim_floor);
    // **The no-gain memory's time-to-live**, the second of the two live-compatible expiries the
    // record names. The shipped one is a *state* rule keyed on the structured weight; this is a
    // *clock* rule. `-1` means the shipped behaviour (no clock), and they compose when both are
    // on, because a memory must satisfy both to stand.
    let ttl_arg: i64 = arg(12, -1);
    let no_gain_ttl: Option<u32> = (ttl_arg >= 0).then(|| ttl_arg as u32);
    let fp = {
        let f = std::fs::File::open(&file).expect("open fcache");
        prin_rs::output::fcache::read(&mut std::io::BufReader::new(f)).expect("read fcache")
    };
    let t = target(&fp.region).unwrap_or_else(|| panic!("the file's region `{}` is not a known target", fp.region));
    let stem = std::path::Path::new(&file).file_stem().unwrap().to_string_lossy().to_string();
    let tag = format!("{}{}{}{}{}",
        if stationary { "" } else { "_nostat" },
        if alpha_lo != SchedCfg::default().alpha_lo { format!("_alo{alpha_lo}") } else { String::new() },
        if merge { "" } else { "_nomerge" },
        if agreement { "" } else { "_noagree" },
        if dim_floor { "" } else { "_nodim" });
    // **A self-describing filename has to carry every setting that is swept**, or the last writer
    // wins over a stem that says nothing about the difference -- the `criterion_sweep` failure.
    let tag = format!("{tag}{}", no_gain_ttl.map_or(String::new(), |t| format!("_ttl{t}")));
    let log = Log::tee(&format!("{root}/output/payload_march_{stem}{tag}.txt"));
    let log = &log;
    let class = if fp.has_event_class() { ClassArm::EventClass } else { ClassArm::Outcome };
    logln!(log, "payload_metric march: {file} -- region {}, levels {} N={} res {}, t_max {}; tolerance policy, stationary {stationary}, eps {eps:e} k_frac {k_frac}, live_stride {live_stride}, alpha_lo {alpha_lo} merge {merge} agreement {agreement} dim_floor {dim_floor} no_gain_ttl {no_gain_ttl:?}; class arm {}",
           fp.region, fp.levels, fp.n, fp.res, fp.t_max, class.name());

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * fp.t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max: fp.t_max,
        n_sync,
        keep_live_series: true,
        live_stride,
        ..Default::default()
    };
    logln!(log, "  config: {}", ens.provenance());
    let metrics = [
        Metric::Payload { eps, class, form: PayloadForm::Indicator },
        Metric::Payload { eps, class, form: PayloadForm::Resolvable },
    ];
    let caches: Vec<Cache> = metrics.iter().map(|&m| Cache::from_footprints(&fp, m).expect("cache")).collect();

    let cfg = SchedCfg {
        n: fp.n,
        tau_display: eps,
        policy: prin_rs::scheduler::Policy::Tolerance,
        stationary,
        alpha_lo,
        merge,
        no_gain_ttl,
        agreement,
        dim_floor,
        k_frac,
        budget: fp.quads.len() * 2,
        camera: Some(Camera::framing(t.cx, t.cy, t.half, fp.res)),
        max_level: Some(fp.levels),
        chart: t.chart,
        keep_pixels: false,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let (tree, st) = scheduler::descend_live(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64);
    let leaves: Vec<usize> = tree.leaves().collect();
    let depth = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
    let steps: u64 = tree.nodes.iter().map(|q| q.red.total_substeps as u64).sum();
    logln!(log, "  live descent: {} quads ({} leaves, depth {depth}) in {:.1}s, {steps:.3e} substeps, catch-up {:.3e} substeps ({:.1}% of the total), stop [{}]",
           st.quads_computed, leaves.len(), t0.elapsed().as_secs_f64(), st.catchup_substeps as f64,
           100.0 * st.catchup_substeps as f64 / steps.max(1) as f64, tree.stop_breakdown());
    logln!(log, "  memory: {} quads ever computed (mem_all), resident peak {} and final {} (mem_resident), {} children merged",
           st.quads_computed, st.resident_peak, st.resident_final, st.merged);
    logln!(log, "  {:>3} {:>7} {:>9} {:>8} {:>7} {:>6} {:>6} {:>10} {:>8} {:>6}", "j", "t", "computed", "resident", "leaves", "split", "keep", "stationary", "deferred", "merged");
    for p in &st.live {
        logln!(log, "  {:>3} {:>7.3} {:>9} {:>8} {:>7} {:>6} {:>6} {:>10} {:>8} {:>6}", p.j, p.t, p.computed, p.resident, p.leaves, p.split, p.keep, p.stationary, p.deferred, p.merged);
    }
    let by_level = {
        let mut m = std::collections::BTreeMap::new();
        for &i in &leaves { *m.entry(tree.nodes[i].level).or_insert(0usize) += 1; }
        m.into_iter().map(|(l, n)| format!("{l}:{n}")).collect::<Vec<_>>().join(" ")
    };
    logln!(log, "  leaves by level: {by_level}");
    let mut keys = Vec::with_capacity(leaves.len());
    let mut lost = 0usize;
    for &i in &leaves {
        let q = &tree.nodes[i];
        match caches[0].key_of(q.cx, q.cy, q.level) {
            Some(k) => keys.push(k),
            None => lost += 1,
        }
    }
    assert_eq!(lost, 0, "{lost} leaves do not map onto the cache");
    for c in &caches {
        let e = c.error_of(&keys);
        let dp = c.dp_optimal((c.quads.len() - 1) / 4);
        let need = dp.budget_needed(e);
        let uni = metric::replay(c, Rank::Uniform, c.quads.len());
        let uni_need = Cache::budget_needed(&uni, e).map(|p| p.budget);
        logln!(log, "  {}: final error {e:.5}; dp needs B={} (ratio {}); uniform needs B={} (ratio {}); sea_fraction {:.4}",
               c.metric.name(),
               need.map(|x| x.to_string()).unwrap_or("--".into()),
               need.map(|x| format!("{:.2}x", st.quads_computed as f64 / x as f64)).unwrap_or("--".into()),
               uni_need.map(|x| x.to_string()).unwrap_or("--".into()),
               uni_need.map(|x| format!("{:.2}x", st.quads_computed as f64 / x as f64)).unwrap_or("--".into()),
               c.sea_fraction(eps));
    }
}

/// **A real descent, scored against a committed cache.** The cache supplies the reference and
/// the ceiling from its footprint file alone (no reductions needed for `error_of` or the DP);
/// the descent integrates its own quads — at most a few percent of the cache's cost — under the
/// named policy, and is mapped onto the cache's keys. Phase 3.1: `Policy::Tolerance` against
/// `Policy::Alpha` on the same box, the same eps, the same ceiling.
fn live() {
    let file: String = std::env::args().nth(2).expect("a .fcache path");
    let policy = std::env::args().nth(3).unwrap_or_else(|| "tolerance".into());
    let policy = prin_rs::scheduler::Policy::parse(&policy).expect("policy: tolerance | alpha | sibling");
    let eps: f64 = arg(4, 0.01);
    let k_frac: f64 = arg(5, prin_rs::scheduler::K_FRAC_RANKED);
    let root: String = std::env::args().nth(6).unwrap_or_else(|| "results".into());
    // The stationarity stop; the struct's default is the one default (off, by measurement).
    // `1` is the arm the sea-chart control was measured against.
    let stationary: bool = std::env::args().nth(7).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().stationary);
    // **The descent's tolerance, separately from the metric's.** `spread_shape` is a MEAN
    // deviation from the copies' centroid, halved, so a footprint reading `spread <= tau` can
    // still hold a pixel whose chord to the texel exceeds `tau` by a factor near three: measured
    // on near-field, the tree built at `tau = 0.003` leaves 5% of pixels unresolved at
    // `eps = 0.003` and none at `eps = 0.01`. The factor between the two is a calibration to
    // measure, not a constant to assume; this argument is how.
    let tau: f64 = arg(8, eps);
    // The floor on the area exponent; `0` is the opt-in that allows full depth on a sea.
    let alpha_lo: f64 = arg(9, SchedCfg::default().alpha_lo);
    let agreement: bool = std::env::args().nth(10).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().agreement);
    // The dimension floor: off leaves only the noise stop, so a sea floors and a fat fractal
    // does not. Measured against on, it says which of the two is doing the flooring.
    let dim_floor: bool = std::env::args().nth(11).map(|v| v != "0" && v != "false").unwrap_or(SchedCfg::default().dim_floor);
    // **The no-gain memory's time-to-live**, the second of the two live-compatible expiries the
    // record names. The shipped one is a *state* rule keyed on the structured weight; this is a
    // *clock* rule. `-1` means the shipped behaviour (no clock), and they compose when both are
    // on, because a memory must satisfy both to stand.
    let ttl_arg: i64 = arg(12, -1);
    let no_gain_ttl: Option<u32> = (ttl_arg >= 0).then(|| ttl_arg as u32);
    // **The camera, as an arm rather than a constant.** With it on, a quad whose texels fall below
    // display resolution is stopped by `ScreenFloor` -- a veto, not a criterion decision -- so a
    // saving quoted in quads can be a fact about the viewport. Off, the descent runs to `max_level`
    // and the stop breakdown is the criterion's alone. Neither is the honest number by itself: the
    // veto is what a real frame does and the criterion is what is being measured. Run both.
    let camera_on: bool = std::env::args().nth(13).map(|v| v != "0" && v != "false").unwrap_or(true);
    let fp = {
        let f = std::fs::File::open(&file).expect("open fcache");
        prin_rs::output::fcache::read(&mut std::io::BufReader::new(f)).expect("read fcache")
    };
    let t = target(&fp.region).unwrap_or_else(|| panic!("the file's region `{}` is not a known target", fp.region));
    let stem = std::path::Path::new(&file).file_stem().unwrap().to_string_lossy().to_string();
    let tag = format!("{}{}{}{}{}", if camera_on { "" } else { "nocam_" }, policy.name(),
        if policy == prin_rs::scheduler::Policy::Tolerance && !stationary { "_nostat" } else { "" },
        if tau != eps { format!("_tau{tau:e}") } else { String::new() },
        if alpha_lo != SchedCfg::default().alpha_lo { format!("_alo{alpha_lo}") } else { String::new() })
        + if agreement { "" } else { "_noagree" }
        + if dim_floor { "" } else { "_nodim" };
    let log = Log::tee(&format!("{root}/output/payload_live_{stem}_{tag}.txt"));
    let log = &log;
    let class = if fp.has_event_class() { ClassArm::EventClass } else { ClassArm::Outcome };
    logln!(log, "payload_metric live: {file} -- PRQF v{}, region {}, levels {} N={} res {}, t_max {}; policy {} camera {camera_on} stationary {stationary} eps {eps:e} tau {tau:e} alpha_lo {alpha_lo} agreement {agreement} dim_floor {dim_floor} k_frac {k_frac}; class arm {}",
           fp.version, fp.region, fp.levels, fp.n, fp.res, fp.t_max, policy.name(), class.name());

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * fp.t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { refine_flagged: false, t_max: fp.t_max, n_sync, ..Default::default() };
    logln!(log, "  config: {}", ens.provenance());

    let metrics = [
        Metric::Payload { eps, class, form: PayloadForm::Indicator },
        Metric::Payload { eps, class, form: PayloadForm::Resolvable },
    ];
    let caches: Vec<Cache> = metrics.iter().map(|&m| Cache::from_footprints(&fp, m).expect("cache")).collect();

    let cfg = SchedCfg {
        n: fp.n,
        tau_display: tau,
        policy,
        stationary,
        alpha_lo,
        agreement,
        dim_floor,
        k_frac,
        budget: fp.quads.len() * 2,
        camera: camera_on.then(|| Camera::framing(t.cx, t.cy, t.half, fp.res)),
        max_level: Some(fp.levels),
        chart: t.chart,
        keep_pixels: false,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let (tree, st) = scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64);
    let leaves: Vec<usize> = tree.leaves().collect();
    let mut keys = Vec::with_capacity(leaves.len());
    let mut lost = 0usize;
    for &i in &leaves {
        let q = &tree.nodes[i];
        match caches[0].key_of(q.cx, q.cy, q.level) {
            Some(k) => keys.push(k),
            None => lost += 1,
        }
    }
    assert_eq!(lost, 0, "{lost} leaves do not map onto the cache");
    let steps: u64 = tree.nodes.iter().map(|q| q.red.total_substeps as u64).sum();
    let depth = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
    logln!(log, "  descent: {} quads ({} leaves, depth {depth}) in {:.1}s, {steps:.3e} substeps, stop [{}]",
           st.quads_computed, leaves.len(), t0.elapsed().as_secs_f64(), tree.stop_breakdown());
    logln!(log, "  mem_all = {} quads x {} trajectories; mem_leaf = {} quads",
           st.quads_computed, fp.n * fp.n * (ens.n_extra + 1), leaves.len());
    let by_level = {
        let mut m = std::collections::BTreeMap::new();
        for &i in &leaves { *m.entry(tree.nodes[i].level).or_insert(0usize) += 1; }
        m.into_iter().map(|(l, n)| format!("{l}:{n}")).collect::<Vec<_>>().join(" ")
    };
    logln!(log, "  leaves by level: {by_level}");

    // **The second budget line.** An `Undetermined` leaf wants finer `eta`, not finer cells, so it
    // is not bought by subdivision at all -- a saving quoted in quads counts only the first kind of
    // work. Printed as a fraction of leaves rather than a count, because the leaf count is exactly
    // what the comparison varies.
    let n_undet = leaves.iter().filter(|&&i| tree.nodes[i].decision == prin_rs::quad::Decision::Undetermined).count();
    let n_collapsed = leaves.iter().filter(|&&i| tree.nodes[i].decision == prin_rs::quad::Decision::Collapsed).count();
    logln!(log, "  undetermined {n_undet}/{} ({:.4}); collapsed {n_collapsed} ({:.4})",
           leaves.len(), n_undet as f64 / leaves.len() as f64, n_collapsed as f64 / leaves.len() as f64);

    // `alpha_area` at the quads the floor actually stopped. It is `2 - d` for a boundary of box
    // dimension `d`, so this column is a **measured box dimension** and not only a threshold check.
    // Reported as a spread, never a mean: the record's standing rule for this exponent is that its
    // variance lives in the tails.
    // **Read it off the PARENT, not off the floored leaf.** `alpha_area` is recorded on the quad
    // whose split was judged; its children are what carry `Decision::Floor`. Two populations are
    // printed, because they answer different questions: every judged quad is the measured box
    // dimension of the field, and the subset with a floored child is what the threshold acted on.
    // Requiring all four children to be `Floor` was the first cut and returned an empty set on
    // every tree -- `near-field` has three floored leaves in total, so no parent qualifies, and an
    // empty column reads as "the floor never fired".
    {
        let q = |a: &Vec<f64>, f: f64| a[((a.len() - 1) as f64 * f).round() as usize];
        let mut all: Vec<f64> = Vec::new();
        let mut acted: Vec<f64> = Vec::new();
        for i in 0..tree.nodes.len() {
            let Some(x) = tree.nodes[i].alpha_area else { continue };
            if !x.is_finite() { continue }
            all.push(x);
            if tree.nodes[i].children.map_or(false, |ch| {
                ch.iter().any(|&k| tree.nodes[k].decision == prin_rs::quad::Decision::Floor)
            }) {
                acted.push(x);
            }
        }
        all.sort_by(|x, y| x.partial_cmp(y).unwrap());
        acted.sort_by(|x, y| x.partial_cmp(y).unwrap());
        if all.is_empty() {
            logln!(log, "  alpha_area: -- (no quad recorded one: the exponent was never judged)");
        } else {
            logln!(log, "  alpha_area, all judged quads: n {} p10 {:+.4} p50 {:+.4} p90 {:+.4}  (box dim d = 2 - alpha, p50 d = {:.4})",
                   all.len(), q(&all, 0.1), q(&all, 0.5), q(&all, 0.9), 2.0 - q(&all, 0.5));
            if acted.is_empty() {
                logln!(log, "  alpha_area, quads with a floored child: n 0 (the floor acted on nothing)");
            } else {
                logln!(log, "  alpha_area, quads with a floored child: n {} p10 {:+.4} p50 {:+.4} p90 {:+.4}  (p50 d = {:.4})",
                       acted.len(), q(&acted, 0.1), q(&acted, 0.5), q(&acted, 0.9), 2.0 - q(&acted, 0.5));
            }
        }
    }

    for c in &caches {
        let e = c.error_of(&keys);
        let dp = c.dp_optimal((c.quads.len() - 1) / 4);
        let need = dp.budget_needed(e);
        let uni = metric::replay(c, Rank::Uniform, c.quads.len());
        let uni_need = Cache::budget_needed(&uni, e).map(|p| p.budget);
        logln!(log, "  {}: error {e:.5}; dp needs B={} (ratio {}); uniform needs B={} (ratio {}); sea_fraction {:.4}",
               c.metric.name(),
               need.map(|x| x.to_string()).unwrap_or("--".into()),
               need.map(|x| format!("{:.2}x", st.quads_computed as f64 / x as f64)).unwrap_or("--".into()),
               uni_need.map(|x| x.to_string()).unwrap_or("--".into()),
               uni_need.map(|x| format!("{:.2}x", st.quads_computed as f64 / x as f64)).unwrap_or("--".into()),
               c.sea_fraction(eps));
    }
    // **The panels, from the cache alone -- no integration enters here.** Written only when a
    // panel root is given as argument 14, so every cell already computed stays valid and
    // nothing has to be re-run. Three, always: the tolerance tree's wire, `Rank::Uniform` at
    // the SAME quad count -- the equal-budget comparison the ratio columns are about -- and the
    // full-depth reference the two are approximating. The wire says where the tree cut; the
    // reference says what there was to cut around, and neither substitutes for the other.
    if let Some(proot) = std::env::args().nth(14).filter(|s| !s.is_empty()) {
        let c = &caches[0];
        let full = (c.quads.len() - 1) / 4;
        let deepest: Vec<metric::Key> = {
            let w = 1u32 << c.levels;
            (0..w).flat_map(|iy| (0..w).map(move |ix| (c.levels, ix, iy))).collect()
        };
        let uni_eq = c.leaves_at(Rank::Uniform, st.quads_computed.min(full));
        let _ = std::fs::create_dir_all(&proot);
        // A magenta count per panel rather than a silent render: `Cache::render` paints a
        // non-finite shape `DEBUG_NAN`, and a reader has to be told how much of a figure is the
        // instrument reporting rather than the field.
        let magenta = |img: &[u8]| -> usize {
            img.chunks_exact(3).filter(|p| *p == prin_rs::output::colour::DEBUG_NAN).count()
        };
        for (name, img) in [
            ("tolerance_wire", c.render_wire(&keys)),
            ("uniform_wire", c.render_wire(&uni_eq)),
            ("reference", c.render(&deepest)),
        ] {
            logln!(log, "  panel {name}: {} magenta of {}", magenta(&img), c.res * c.res);
            let _ = prin_rs::output::adaptive::save(&format!("{proot}/{stem}_eps{eps:e}_{name}.png"), c.res, &img);
        }
    }
}

fn build() {
    let name: String = arg(2, "near-field".to_string());
    let levels: u32 = arg(3, 6);
    let n: usize = arg(4, 8);
    let eps: f64 = arg(5, 0.01);
    let t_max: f64 = arg(6, 13.0);
    let root: String = std::env::args().nth(7).unwrap_or_else(|| "results".into());
    let t = target(&name).unwrap_or_else(|| panic!("unknown target `{name}`: a grid::REGIONS name or a gallery case"));
    let res = (1usize << levels) * n;
    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { refine_flagged: false, t_max, n_sync, keep_boundary_shapes: true, ..Default::default() };
    let dir = format!("{root}/payload");
    let _ = std::fs::create_dir_all(&dir);
    scheduler::assert_production_kernel(&ens, &dir);
    let stem = format!("{dir}/{}_t{t_max}_L{levels}", t.name.replace(' ', "_"));
    let log = Log::tee(&format!("{root}/output/payload_metric_{}_L{levels}.txt", t.name.replace(' ', "_")));
    let log = &log;

    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;
    logln!(log, "payload_metric build: target {} ({} at ({}, {}) half {}), levels {levels}, N={n}, E+1={}, res {res}^2, eps={eps:e}, t={t_max} n_sync={n_sync}",
           t.name, t.chart.name(), t.cx, t.cy, t.half, ens.n_extra + 1);
    logln!(log, "  config: {}", ens.provenance());
    logln!(log, "  {full} quads, {} trajectories", full * n * n * (ens.n_extra + 1));

    let metrics = [
        Metric::Payload { eps, class: ClassArm::EventClass, form: PayloadForm::Indicator },
        Metric::Payload { eps, class: ClassArm::EventClass, form: PayloadForm::Hinge },
        Metric::Payload { eps, class: ClassArm::EventClass, form: PayloadForm::Resolvable },
        Metric::Payload { eps, class: ClassArm::Outcome, form: PayloadForm::Indicator },
        Metric::Colour(SHIPPING),
    ];
    let t0 = std::time::Instant::now();
    let (caches, px_of) = metric::build_metrics_with_footprints(
        &t.name, t.cx, t.cy, t.half, t.body, t.chart, levels, n, res, SchedCfg::default().tau_display, &ens, &metrics,
    );
    logln!(log, "  built in {:.1}s, {} trajectories", t0.elapsed().as_secs_f64(), caches[0].trajectories);

    // **The truncation counter, added 2026-09-06.** This header printed no footprint health at
    // all, so a build at `max_steps = 30_000` and one at 480_000 were distinguishable only by the
    // trees they later produced -- and 15 of the 26 gallery trees moved on that field alone.
    //
    // `budget` is the count the step budget stopped; `nonfin` is the count with any unusable copy,
    // which is the superset. They are printed **separately** and not as a ratio: a budget stop is
    // curable by `max_steps`, a non-finite copy from a triple collision is the instrument
    // reporting, and pooling them would put the two in the column that exists to separate them.
    // A build reading `budget 0` predicts bitwise reproduction under a budget change; one reading
    // `budget > 0` predicts a move whose size tracks the count.
    let (mut n_fp, mut n_budget, mut n_nonfin) = (0usize, 0usize, 0usize);
    for v in px_of.values() {
        for p in v {
            n_fp += 1;
            if p.budget_exhausted {
                n_budget += 1;
            }
            if p.n_nonfinite > 0 {
                n_nonfin += 1;
            }
        }
    }
    logln!(log, "  footprints {n_fp}: budget-stopped {n_budget}, any unusable copy {n_nonfin}");

    // The footprints (v2, with the event class) and the quad cache of the headline.
    if let Ok(f) = std::fs::File::create(format!("{stem}.fcache")) {
        let mut w = std::io::BufWriter::new(f);
        let _ = prin_rs::output::fcache::write(&mut w, &caches[0].footprints_from(&px_of, t_max));
    }
    if let Ok(f) = std::fs::File::create(format!("{stem}.qcache")) {
        let mut w = std::io::BufWriter::new(f);
        let _ = prin_rs::output::qcache::write(&mut w, &caches[0], &ens, SchedCfg::default().tau_display);
    }
    // The replay path, checked live: remeasuring the headline from its own footprints must
    // reproduce it bitwise, or every replay number is suspect.
    {
        let fp = caches[0].footprints_from(&px_of, t_max);
        let r = caches[0].remeasure(&fp, metrics[0]).expect("remeasure");
        let same = r.quads.iter().all(|(k, q)| q.err_sum.to_bits() == caches[0].quads[k].err_sum.to_bits());
        logln!(log, "  replay check: remeasure from footprints reproduces err_sum {}", if same { "BITWISE" } else { "**DIFFERENTLY**" });
    }
    drop(px_of);

    let ranks = vec![
        Rank::Uniform,
        Rank::GreedyLookahead1,
        Rank::Signal(Criterion::Within, Agg::Median),
        Rank::Signal(Criterion::Between, Agg::Median),
        Rank::Signal(Criterion::FracHotBetween, Agg::Median),
        Rank::Signal(Criterion::PerimeterWithin, Agg::Median),
        Rank::Signal(Criterion::FirstDivergence, Agg::Median),
        Rank::Signal(Criterion::GradRms, Agg::Median),
        Rank::Random(1),
        Rank::Random(2),
        Rank::Random(3),
    ];
    // The static live tree under the shipped policy, once; scored under every metric.
    let live: Vec<_> = live_tree(&caches[0], &t, &ens, res, "descend").into_iter().collect();
    for c in &caches {
        tables(log, c, &ranks, &live);
    }
    logln!(log, "wrote {stem}.fcache (PRQF v2) and {stem}.qcache");
}

fn replay() {
    let file: String = std::env::args().nth(2).expect("a .fcache path");
    let eps: f64 = arg(3, 0.01);
    let root: String = std::env::args().nth(4).unwrap_or_else(|| "results".into());
    let fp = {
        let f = std::fs::File::open(&file).expect("open fcache");
        prin_rs::output::fcache::read(&mut std::io::BufReader::new(f)).expect("read fcache")
    };
    let stem = std::path::Path::new(&file).file_stem().unwrap().to_string_lossy().to_string();
    let log = Log::tee(&format!("{root}/output/payload_replay_{stem}.txt"));
    let log = &log;
    logln!(log, "payload_metric replay: {file} -- PRQF v{}, region {}, chart {}, levels {} N={} res {}, t_max {}",
           fp.version, fp.region, fp.chart, fp.levels, fp.n, fp.res, fp.t_max);
    logln!(log, "  ZERO trajectories. Reductions are absent from a footprint file, so only the floor \
                 (uniform), the reference (greedy), the random band and the ceiling (dp) are scored.");
    let class = if fp.has_event_class() { ClassArm::EventClass } else { ClassArm::Outcome };
    if !fp.has_event_class() {
        logln!(log, "  v1 file: no event class stored, replaying under ClassArm::Outcome, which is \
                     saturated at t = 13 -- a machinery check, not the measurement.");
    }
    let metrics = [
        Metric::Payload { eps, class, form: PayloadForm::Indicator },
        Metric::Payload { eps, class, form: PayloadForm::Hinge },
        Metric::Payload { eps, class, form: PayloadForm::Resolvable },
        Metric::Colour(SHIPPING),
    ];
    let ranks = vec![Rank::Uniform, Rank::GreedyLookahead1, Rank::Random(1), Rank::Random(2), Rank::Random(3)];
    let live: Vec<(String, Vec<Key>, usize, u64, String)> = Vec::new();
    for m in metrics {
        match Cache::from_footprints(&fp, m) {
            Ok(c) => tables(log, &c, &ranks, &live),
            Err(e) => logln!(log, "=== metric {} === REFUSED: {e}", m.name()),
        }
    }
    let _ = HashMap::<Key, ()>::new();
}
