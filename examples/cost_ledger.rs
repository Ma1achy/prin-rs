//! **What the refinement criterion costs against uniform-at-max-depth: quads, substeps, bytes.**
//!
//! The plan's objective is **memory**, and every table in `results/payload` reports it as a quad
//! count. A quad count is a memory model and it is **not** a compute model, because quads do not
//! cost the same: `total_substeps` per trajectory varies across a frame, and the criterion
//! refines exactly where the physics is hard, which is where trajectories are expensive. So the
//! two ratios can differ, and this harness is what says by how much.
//!
//! Three arms, and the middle one is the control the user named:
//!
//! - **uniform at max depth** — the `4^L` deepest quads and nothing else, one sample per pixel.
//!   "Just go to max depth everywhere." It has no parents, so it is the *cheapest* uniform arm.
//! - **uniform tree** — every quad at every level, `(4^(L+1)-1)/3`. What `Rank::Uniform` counts
//!   and what `B` counts, because the design keeps parents marching. A third more than the above.
//! - **the tolerance tree** — the descent under `Policy::Tolerance`, parents included.
//!
//! **The cost axis is `total_substeps`, never seconds** — the standing rule, and this table is
//! where it bites: a substep is machine-independent, and the same tree timed on a loaded machine
//! reads faster or slower than one doing more work. Force evaluations are `4 x substeps` for RK4
//! to within 0.03% (measured: Heggie `steps p50 1.536e5` against `evals p50 6.146e5` in
//! `results/output/logh_arms.txt`, the excess being the secant landing), so substeps and evals
//! are one column here; they are **not** for leapfrog, which reads 1:1, and any table with two
//! steppers in it has to carry both.
//!
//! **The core-seconds column is a floor, not a prediction.** It prices the `deriv` calls alone at
//! the measured 27.40 ns for Heggie's Eqs. (22)-(24), and the measured cost per substep is
//! several times that: the build of one levels-6 cache ran 1.520e10 substeps in 681.9 s on 11
//! cores, about 493 ns per substep against 110 ns of `deriv`. What transfers between machines is
//! the **ratio**; the constant does not, and it is printed so the gap is visible rather than
//! implied.
//!
//! **The guard.** The adaptive arm's substeps come from a fresh descent and the uniform arm's
//! from the committed cache, so they are two runs. Every quad the descent computes is compared
//! against the cache's own count for the same key: *a difference can be small because both sides
//! are right or because both are dead*, and two arms measured under different kernels would
//! produce a plausible ratio out of nothing. The per-quad agreement is printed before the table.
//!
//! **The error column is what makes it a comparison.** A cheaper tree that displays worse has
//! not saved anything, so the adaptive arm is scored against the same cache, on the same metric
//! the cache was written under (`payload/event_class/indicator/eps=1e-2`): `sum err_sum` over its
//! leaves, divided by `res^2`. Uniform at max depth is exactly `0` there by construction -- it
//! *is* the reference -- so a saving is only real where the adaptive column reads zero too.
//!
//! **Bytes, two of them, and they answer different questions.** `march MB` is `HgState<f64>` per
//! trajectory: what a live playhead must hold to keep every trajectory advancing, which is the
//! quantity the plan's memory objective is about. `record MB` is `PixelOut` per footprint: what
//! the finished frame costs to keep. Neither includes the live series, which is optional and
//! whose size is a stride choice.
//!
//! Run: `cargo run --release --example cost_ledger -- <dir=results/payload> [charts] [descend=1]`
//! With `descend=0` it reads only the caches and prints the uniform ladder — zero trajectories.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::metric::Key;
use prin_rs::output::qcache;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

/// Measured single-threaded cost of one Heggie `deriv` call, Eqs. (22)-(24). §31.
const DERIV_NS: f64 = 27.40;
/// RK4 stages. Force evaluations per substep, exact for this stepper to 0.03%.
const EVALS_PER_STEP: f64 = 4.0;
/// Bytes a marching trajectory must hold: `HgState<f64>` is `[u(2), p(2)] x 3 + t`, thirteen
/// numbers. Not `size_of::<PixelOut>()`, which is the *result* record and 6.3x larger.
const MARCH_BYTES: usize = std::mem::size_of::<prin_rs::integrate::heggie::state::HgState<f64>>();

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn quantile(v: &mut [f64], p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = (v.len() - 1) as f64 * p;
    let (lo, hi) = (i.floor() as usize, (i.ceil() as usize).min(v.len() - 1));
    v[lo] + (v[hi] - v[lo]) * (i - i.floor())
}

struct Target {
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
}

fn target(name: &str) -> Option<Target> {
    if let Some(&(_, cx, cy, body)) = grid::REGIONS.iter().find(|r| r.0 == name) {
        return Some(Target { chart: Chart::BodyPlane, cx, cy, half: 0.05, body });
    }
    // Every named slice with its own window, from the one table in `grid`.
    if let Some((chart, cx, cy, half)) = grid::named_slice(name) {
        return Some(Target { chart, cx, cy, half, body: 0 });
    }
    grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == name)
        .map(|(_, chart, cx, cy, half)| Target { chart, cx, cy, half, body: 0 })
}

fn main() {
    let dir: String = std::env::args().nth(1).unwrap_or_else(|| "results/payload".into());
    let charts: Vec<String> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "near-field,deep_interior,far,preset_shape_h1".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let descend: bool = std::env::args().nth(3).map(|v| v != "0" && v != "false").unwrap_or(true);
    let eps: f64 = arg(4, SchedCfg::default().tau_display);
    let k_frac: f64 = arg(5, prin_rs::scheduler::K_FRAC_RANKED);

    println!("cost_ledger: {dir}, charts {charts:?}, descend {descend}, eps {eps:e}, k_frac {k_frac}");
    println!("  cost axis = total_substeps (machine-independent); evals = 4 x substeps for RK4;");
    println!("  core-s floor = substeps x {EVALS_PER_STEP} x {DERIV_NS} ns -- deriv arithmetic only, ~22% of the measured substep cost.");
    println!("  resident bytes per trajectory: HgState<f64> {} B marching, PixelOut {} B as a result record.",
             std::mem::size_of::<prin_rs::integrate::heggie::state::HgState<f64>>(),
             std::mem::size_of::<PixelOut>());

    // ---- the uniform ladder, from the committed caches: zero trajectories -------------------
    println!("\n=== the uniform ladder, per level, read from the caches ===");
    println!("  {:>16} {:>3} {:>6} {:>13} {:>11} {:>10} {:>9} {:>9} {:>9} {:>8}",
             "region", "lvl", "quads", "trajectories", "substeps", "st/traj", "p10", "p90", "max", "p90/p10");
    let mut ladders: Vec<(String, qcache::QuadRows)> = Vec::new();
    for name in &charts {
        let path = format!("{dir}/{name}_t13_L6.qcache");
        let Ok(f) = std::fs::File::open(&path) else {
            println!("  {name:>16}  no cache at {path}, skipped");
            continue;
        };
        let qr = match qcache::read(&mut std::io::BufReader::new(f)) {
            Ok(q) => q,
            Err(e) => {
                println!("  {name:>16}  {e}");
                continue;
            }
        };
        let tr = qr.trajectories_per_quad();
        for lv in 0..=qr.levels {
            let rows: Vec<usize> = (0..qr.rows.len()).filter(|&i| qr.key(i).0 == lv).collect();
            let sub: f64 = rows.iter().map(|&i| qr.get(i, "total_substeps")).sum();
            let mut per: Vec<f64> = rows.iter().map(|&i| qr.get(i, "total_substeps") / tr as f64).collect();
            let (p10, p90, mx) = (quantile(&mut per, 0.10), quantile(&mut per, 0.90),
                                  per.last().copied().unwrap_or(f64::NAN));
            println!("  {:>16} {lv:>3} {:>6} {:>13} {:>11.4e} {:>10.0} {p10:>9.0} {p90:>9.0} {mx:>9.0} {:>8.2}",
                     if lv == 0 { name.as_str() } else { "" }, rows.len(), rows.len() as u64 * tr,
                     sub, sub / (rows.len() as u64 * tr) as f64, p90 / p10.max(1.0));
        }
        ladders.push((name.clone(), qr));
    }

    // ---- the ledger --------------------------------------------------------------------------
    println!("\n=== the ledger: the tolerance tree against uniform at max depth ===");
    println!("  {:>16} | {:>6} {:>11} {:>8} {:>8} {:>8} {:>8} | {:>6} {:>11} {:>8} {:>8} | {:>7} {:>7} {:>7}",
             "region", "quads", "substeps", "st/traj", "error", "march MB", "rec MB",
             "quads", "substeps", "march MB", "rec MB", "mem", "cpu", "cpu/mem");
    println!("  {:>16} | {:^54} | {:^38} | {:^23}",
             "", "----------------- tolerance tree -----------------",
             "------- uniform at max depth -------", "ratios to uniform");

    for (name, qr) in &ladders {
        let tr = qr.trajectories_per_quad();
        let deep: Vec<usize> = (0..qr.rows.len()).filter(|&i| qr.key(i).0 == qr.levels).collect();
        let u_quads = deep.len() as u64;
        let u_sub: f64 = deep.iter().map(|&i| qr.get(i, "total_substeps")).sum();
        let t_quads = qr.rows.len() as u64;
        let t_sub: f64 = (0..qr.rows.len()).map(|i| qr.get(i, "total_substeps")).sum();

        if !descend {
            println!("  {name:>16} |  (descend=0: uniform max depth {u_quads} quads {u_sub:.4e} substeps; \
                      uniform tree {t_quads} quads {t_sub:.4e} substeps)");
            continue;
        }

        let Some(t) = target(name.replace('_', " ").as_str()).or_else(|| target(name)) else {
            println!("  {name:>16}  not a known target, skipped");
            continue;
        };
        let base = EnsembleCfg::default();
        let t_max: f64 = qr.header.split_whitespace()
            .find_map(|s| s.strip_prefix("t_max=")).and_then(|s| s.parse().ok()).unwrap_or(base.t_max);
        let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
        let ens = EnsembleCfg { refine_flagged: false, t_max, n_sync, ..Default::default() };
        let res: usize = qr.header.split_whitespace()
            .find_map(|s| s.strip_prefix("res=")).and_then(|s| s.parse().ok()).unwrap_or(512);
        let cfg = SchedCfg {
            n: qr.n,
            tau_display: eps,
            policy: scheduler::Policy::Tolerance,
            k_frac,
            budget: qr.rows.len() * 2,
            camera: Some(Camera::framing(t.cx, t.cy, t.half, res)),
            max_level: Some(qr.levels),
            chart: t.chart,
            keep_pixels: false,
            ..Default::default()
        };
        let t0 = std::time::Instant::now();
        let (tree, st) = scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64);
        let secs = t0.elapsed().as_secs_f64();
        let a_quads = st.quads_computed as u64;
        let a_sub: f64 = tree.nodes.iter().map(|q| q.red.total_substeps as f64).sum();

        // **The guard.** Two runs, one ratio: assert the descent and the cache integrated the
        // same trajectories before differencing their costs. A quad the cache does not hold is
        // counted separately -- that is a mapping failure, not a disagreement.
        let mut cache_by_key: std::collections::HashMap<Key, f64> = std::collections::HashMap::new();
        let mut err_by_key: std::collections::HashMap<Key, f64> = std::collections::HashMap::new();
        for i in 0..qr.rows.len() {
            cache_by_key.insert(qr.key(i), qr.get(i, "total_substeps"));
            err_by_key.insert(qr.key(i), qr.get(i, "err_sum"));
        }
        let (mut matched, mut unmapped, mut worst) = (0usize, 0usize, 0.0f64);
        for q in tree.nodes.iter() {
            let h = t.half / (1u64 << q.level) as f64;
            let ix = ((q.cx - (t.cx - t.half)) / (2.0 * h) - 0.5).round();
            let iy = ((q.cy - (t.cy - t.half)) / (2.0 * h) - 0.5).round();
            match cache_by_key.get(&(q.level, ix as u32, iy as u32)) {
                Some(&c) if c > 0.0 => {
                    matched += 1;
                    worst = worst.max((q.red.total_substeps as f64 - c).abs() / c);
                }
                _ => unmapped += 1,
            }
        }
        println!("  {name:>16}   guard: {matched} of {} quads matched to the cache ({unmapped} unmapped), \
                  worst per-quad substep disagreement {:.3e}", tree.nodes.len(), worst);

        // **Scored against the same cache, on the metric the cache was written under.** A
        // cheaper tree that displays worse has saved nothing, and uniform at max depth reads
        // exactly 0 here because it *is* the reference -- so a saving counts only where the
        // adaptive column reads zero too.
        let mut err = 0.0f64;
        let mut unscored = 0usize;
        for i in tree.leaves() {
            let q = &tree.nodes[i];
            let h = t.half / (1u64 << q.level) as f64;
            let ix = ((q.cx - (t.cx - t.half)) / (2.0 * h) - 0.5).round();
            let iy = ((q.cy - (t.cy - t.half)) / (2.0 * h) - 0.5).round();
            match err_by_key.get(&(q.level, ix as u32, iy as u32)) {
                Some(&e) => err += e,
                None => unscored += 1,
            }
        }
        let err = if unscored > 0 { f64::NAN } else { err / (res * res) as f64 };

        let (mem, cpu) = (a_quads as f64 / u_quads as f64, a_sub / u_sub);
        let mb = |quads: u64, per: usize| (quads * tr) as f64 * per as f64 / 1.048576e6;
        let rec = |quads: u64| (quads * (qr.n * qr.n) as u64) as f64
            * std::mem::size_of::<PixelOut>() as f64 / 1.048576e6;
        println!("  {name:>16} | {a_quads:>6} {a_sub:>11.4e} {:>8.0} {err:>8.5} {:>8.1} {:>8.1} | \
                  {u_quads:>6} {u_sub:>11.4e} {:>8.1} {:>8.1} | {mem:>7.4} {cpu:>7.4} {:>7.2}",
                 a_sub / (a_quads * tr) as f64, mb(a_quads, MARCH_BYTES), rec(a_quads),
                 mb(u_quads, MARCH_BYTES), rec(u_quads), cpu / mem);
        if unscored > 0 {
            println!("  {:>16} |  {unscored} leaves do not map onto the cache; the error column is refused", "");
        }
        let core_s = (u_sub - a_sub) * EVALS_PER_STEP * DERIV_NS * 1e-9;
        println!("  {:>16} |  deriv-arithmetic floor on the compute saved: {core_s:.1} core-s", "");
        println!("  {:>16} |  uniform tree with parents: {t_quads} quads {t_sub:.4e} substeps \
                  -> mem {:.4} cpu {:.4}; descent wall {secs:.1}s on {} cores",
                 "", a_quads as f64 / t_quads as f64, a_sub / t_sub,
                 std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0));
    }
    println!("\n`cpu/mem` is the finding: 1.0 means the quad count is an honest compute model, and");
    println!("above 1.0 means the criterion is refining into the expensive trajectories, so the");
    println!("memory saving overstates the compute saving by that factor.");
}
