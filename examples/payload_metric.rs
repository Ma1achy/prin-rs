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
//! Two modes:
//!
//! - `build <target> [levels] [n] [eps] [t_max] [root]` — integrate one complete tree for a
//!   Burrau region (`grid::REGIONS`) or a gallery chart (`grid::gallery_cases`), write the v2
//!   footprint file and the quad cache, and print every table for every form plus the shipping
//!   OKLab metric as the comparison arm. Signal rankings and the **static live tree** (a real
//!   `scheduler::descend` on the same box, mapped onto the cache) are scored here.
//! - `replay <file.fcache> [eps] [root]` — zero trajectories: the floor (`uniform`), the
//!   reference (`greedy_lookahead_1`), the random band and the ceiling (`dp_optimal`) from a
//!   committed footprint file. A v1 file is replayed under the **outcome** arm only and says so.
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
        _ => {
            eprintln!("usage: payload_metric build <target> [levels=6] [n=8] [eps=0.01] [t_max=13] [root=results]");
            eprintln!("       payload_metric replay <file.fcache> [eps=0.01] [root=results]");
            std::process::exit(2);
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
