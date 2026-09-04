//! **The frame loop, marched: playhead AND camera moving together, with the record written out.**
//!
//! §9's acceptance test, and it is specific about what makes it able to fail. Depth variance alone
//! cannot tell a **balanced** tree from a **frozen** one — both give a flat curve — so per-quad
//! churn is plotted beside it, and `Mode::Uniform` runs in the same figure as the control that
//! proves the test discriminates. §4.5 adds the other half: every measurement taken during a
//! gesture has **two causes**, so the camera delta and the playhead delta are logged separately, or
//! a churn spike cannot be attributed to either.
//!
//! **Churn is read over SHARED quads and the count is printed.** A quad present at one frame and
//! not the other has not changed its decision; counting it folds the tree's growth into a statistic
//! about its stability. `near-field` at `t = 16` shares fourteen quads, so a churn of 0.4286 is
//! six of fourteen — thin, and labelled thin.
//!
//! **`n_sync` scales with `t_max`.** `dtau = eta*dt_left/(A0*B0)`, so holding `n_sync` fixed while
//! `t_max` varies changes the step size and the rows are different discretisations rather than one
//! trajectory at several playheads.
//!
//! The record is `PRNF` beside a `.cfg.txt`, the same division `.prnq` already has. `frame_ms` is
//! **measured, never budgeted** — the quota is a count, and `tests/session.rs` pins that a sleeping
//! sampler does not move a single decision.
//!
//! Run: `cargo run --release --example session_frames -- [root=results] [chart] [frames=24] [res=256]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::output::frame::{self, FrameRecord, StageMs};
use prin_rs::quad::{Decision, QuadTree};
use prin_rs::scheduler::{Mode, SchedCfg};
use prin_rs::session::{FrameQuota, Regrow, Session, SessionCfg};
use std::collections::HashMap;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn target(name: &str) -> Option<(String, Chart, f64, f64, f64, usize)> {
    if let Some(&(n, cx, cy, body)) = grid::REGIONS.iter().find(|r| r.0 == name) {
        return Some((n.into(), Chart::BodyPlane, cx, cy, 0.05, body));
    }
    if name == "config_stability" {
        let (chart, cx, cy, half) = Chart::config_stability();
        return Some((name.into(), chart, cx, cy, half, 0));
    }
    grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == name)
        .map(|(n, chart, cx, cy, half)| (n.into(), chart, cx, cy, half, 0))
}

/// Variance of leaf depth — the quantity §3.2 plots against `t`.
fn depth_variance(t: &QuadTree) -> f64 {
    let d: Vec<f64> = t.leaves().map(|i| t.nodes[i].level as f64).collect();
    if d.is_empty() {
        return f64::NAN;
    }
    let m = d.iter().sum::<f64>() / d.len() as f64;
    d.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / d.len() as f64
}

/// Decisions keyed by box, so two frames of a growing tree are comparable.
fn by_box(t: &QuadTree) -> HashMap<(u32, i64, i64), Decision> {
    t.nodes
        .iter()
        .filter(|q| q.is_leaf() && !q.merged)
        .map(|q| ((q.level, (q.cx / 1e-12).round() as i64, (q.cy / 1e-12).round() as i64), q.decision))
        .collect()
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let chart: String = std::env::args().nth(2).unwrap_or_else(|| "near-field".into());
    let frames: usize = arg(3, 24);
    let res: usize = arg(4, 256);

    let Some((name, ch, cx, cy, half, body)) = target(&chart) else {
        eprintln!("unknown target `{chart}`");
        std::process::exit(2);
    };
    let dir = format!("{root}/session");
    let _ = std::fs::create_dir_all(&dir);

    // **`n_sync` scales with `t_max`.** Holding it fixed while the horizon moves compares different
    // discretisations, which is a standing result on this project and not a precaution.
    let base = EnsembleCfg::default();
    let t_max = base.t_max;
    let n_sync = base.n_sync;
    // The live series is what makes the playhead real: without it the session can only re-read a
    // quad at the horizon, the field never changes between frames, and the march is camera-only.
    let live_stride = 4;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        keep_live_series: true,
        live_stride,
        ..Default::default()
    };
    let boundaries = n_sync / live_stride;

    println!("session_frames: {name} at {res}^2, {frames} frames, t={t_max} n_sync={n_sync}");
    println!("  config: {}", ens.provenance());
    println!("  the camera pans AND the playhead marches; both deltas are logged separately,");
    println!("  because a churn spike with one cause recorded cannot be attributed to either.\n");

    // Two arms. `Mode::Uniform` is the control that proves the depth-variance test discriminates:
    // the criterion is OFF there, so it must degenerate while the balanced arm does not.
    for mode in [Mode::Balanced, Mode::Uniform] {
        let cfg = SessionCfg {
            sched: SchedCfg {
                n: 8,
                budget: 20000,
                mode,
                camera: Some(Camera::framing(cx, cy, half, res)),
                camera_bias: Some(0.5),
                chart: ch,
                // Payloads must be kept for `set_playhead` to re-reduce from them; an evicted quad
                // keeps the reduction it had, which is the caching contract's resume point.
                keep_pixels: true,
                ..Default::default()
            },
            quota: FrameQuota { quads: 24, substeps: None, rounds: 4 },
            payload_cap: Some(4096),
            regrow: Regrow::Upto(2),
            ..Default::default()
        };
        let mut s = Session::new(cx, cy, half, body, Camera::framing(cx, cy, half, res), cfg, t_max);

        let mut log: Vec<FrameRecord> = Vec::with_capacity(frames);
        let mut prev = by_box(s.tree());
        println!("{:>10} {:>6} {:>8} {:>7} {:>7} {:>8} {:>8} {:>7} {:>7} {:>8} {:>9}",
                 "mode", "frame", "frame_ms", "quads", "leaves", "dvar", "churn", "shared",
                 "reproj", "cam_d", "stop");

        for f in 0..frames {
            // The camera pans across a quarter of the viewport over the run, so `camera_delta` is
            // non-zero on every frame and the headline statistic has a population.
            let pan = half * 0.25 * (f as f64) / frames.max(1) as f64;
            let delta = s.set_camera(Camera::framing(cx + pan, cy, half, res));
            // **The playhead marches too.** One boundary per frame under lockstep, so a churn
            // spike has two possible causes and both are logged.
            let j = (f * boundaries / frames.max(1)).min(boundaries.saturating_sub(1));
            let reprojected = s.set_playhead(j);

            let t0 = std::time::Instant::now();
            let (hit, spend) = s.step(&|sl, k| {
                prin_rs::ensemble::pixel::evaluate::<f64>(sl, k, &ens)
            });
            let integrate_ms = t0.elapsed().as_secs_f64() * 1e3;
            let t1 = std::time::Instant::now();
            let evicted = s.evict_after_frame(&|_| false).len();
            let evict_ms = t1.elapsed().as_secs_f64() * 1e3;

            let now = by_box(s.tree());
            let shared: Vec<_> = prev.keys().filter(|k| now.contains_key(*k)).collect();
            let moved = shared.iter().filter(|k| prev[**k] != now[**k]).count();
            let churn = if shared.is_empty() { f64::NAN } else { moved as f64 / shared.len() as f64 };

            let depth = s.tree().leaves().map(|i| s.tree().nodes[i].level).max().unwrap_or(0);
            let rec = FrameRecord {
                frame: f as u64,
                frame_ms: integrate_ms + evict_ms,
                stage_ms: StageMs { integrate: integrate_ms, evict: evict_ms, ..Default::default() },
                quads_computed: spend.quads,
                quads_evicted: evicted,
                resident_payloads: s.store().resident(),
                substeps_total: spend.substeps,
                // The playhead marches with the frame under lockstep: one fixed dt per frame.
                playhead_dt: t_max / frames.max(1) as f64,
                quads_reused: reprojected,
                camera_delta: delta.d_centre,
                camera_zoom_octaves: delta.d_zoom_octaves,
                tree_depth_max: depth,
                leaf_count: s.tree().leaves().count(),
                rounds: spend.rounds,
                quota_hit: hit as u8,
                frontier_agrees: f64::NAN,
                ..Default::default()
            };
            // `reproj` is the LIVENESS arm for the playhead: a churn of 0.0000 means a steady
            // state only if the quads were actually re-read at the new boundary. Zero here would
            // mean the playhead is not wired, and the two are indistinguishable in the churn
            // column alone.
            println!("{:>10} {f:>6} {:>8.2} {:>7} {:>7} {:>8.4} {:>8.4} {:>7} {reprojected:>7} {:>7.2e} {:>9}",
                     mode.name(), rec.frame_ms, rec.quads_computed, rec.leaf_count,
                     depth_variance(s.tree()), churn, shared.len(), delta.d_centre, hit.name());
            log.push(rec);
            prev = now;
        }

        let ms: Vec<f64> = log.iter().map(|r| r.frame_ms).collect();
        let p = frame::percentiles(&ms);
        let over = frame::frac_over_floor_moving(&log);
        println!("\n  {} : p50 {:.2} p95 {:.2} p99 {:.2} max {:.2} ms over {} frames",
                 mode.name(), p.p50, p.p95, p.p99, p.max, p.n);
        println!("  frac over the {:.1} ms floor, DURING MOTION: {over:.4}{}",
                 frame::BUDGET_24,
                 if over.is_nan() { "  (no frame moved -- the question was not asked)" } else { "" });
        println!("  final depth variance {:.4}\n", depth_variance(s.tree()));

        let stem = format!("{dir}/{}_{}", name.replace(' ', "_"), mode.name());
        if let Ok(mut w) = std::fs::File::create(format!("{stem}.prnf")) {
            let header = format!(
                "region={name} mode={} res={res} frames={frames} quota_quads=24 quota_rounds=4\n\
                 t_max={t_max} n_sync={n_sync} payload_cap=4096 regrow=Upto(2) camera_bias=0.5\n\
                 note=upload_ms and present_ms are NaN: there is no GPU and no window in this \
build, and 0.0 would read as instant where the truth is absent.\n\
                 note=frame_ms is MEASURED, never budgeted. The quota is a count; a sleeping \
sampler does not move a decision (tests/session.rs).\n\
                 config={}",
                mode.name(),
                ens.provenance()
            );
            let _ = frame::write(&mut w, &log, &header);
        }
        let _ = prin_rs::output::provenance_sidecar(
            &format!("{stem}.prnf"),
            &ens,
            &format!("frames={frames} res={res} mode={} pan_total={:.4}\n", mode.name(), half * 0.25),
        );
    }

    println!("{dir}/<region>_<mode>.prnf, with the balanced arm and the uniform control.");
    println!("Depth variance alone cannot tell balanced from FROZEN -- read churn beside it.");
}
