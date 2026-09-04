//! **Does ranking the frame's frontier change the tree, or is it `k_frac = 1.0` again?**
//!
//! The frame quota truncates `pending` — the quads about to be *computed* — and before this it
//! truncated by **position**, which is split order. `Session::rank_pending` puts the persistent
//! frontier there. Whether that moves anything is not obvious and must not be assumed:
//!
//! - A round whose `pending` comes from **one** split is a **total tie**. All four children
//!   inherit their parent's stored term, ties break by id ascending, and id order *is* split
//!   order. The ranking runs and changes nothing.
//! - It can only bite where `pending` spans parents whose signals differ, or where the camera
//!   term separates quads the physics ties.
//!
//! This project has shipped two mechanisms that computed, sorted and changed nothing — `k_frac`
//! at 1.0 for a whole corpus, and `Frontier::top_k` sorting the buckets it had just built. So the
//! measurement is the deliverable, and `rank_frame = false` is the named control.
//!
//! Three arms per case, all else held:
//!
//! | arm | `rank_frame` | `camera_bias` | what it isolates |
//! |---|---|---|---|
//! | `off` | false | — | the pre-wiring tree: quota truncates by split order |
//! | `rank` | true | none | the physics ordering alone |
//! | `rank+cam` | true | 0.5 | physics × camera relevance |
//!
//! **`scan/len` is printed for every arm**, because a ranking that walks its whole frontier has
//! bought nothing over a sort whatever it did to the tree.
//!
//! # The regime is set by the arithmetic, not by a knob
//!
//! The quota truncates `pending`, and `pending` is ~`4s` times the previous round's quota, where
//! `s` is the fraction that split. So `k/n ~ 1/(4s)` sits in `[0.25, 1]` **by construction** —
//! the large-`k` regime, which `results/frontier/` measured as saving 0–67%. The 83–94% headline
//! there was taken at `k/n = 0.01` on the **visible frontier**, a different population entirely.
//! Quoting it for this site would be quoting a number about another measurement.
//!
//! # Two regimes, and the second one is a correctness check
//!
//! **Truncated**: stop after `frames` frames with the quota still binding. The trees may differ,
//! and whether they do is the question. **Drained**: run until the frontier empties. There the
//! ranking has only reordered work that all happened, so the trees must be **identical** — a
//! ranking that changed the destination would be changing the criterion, not the schedule.
//!
//! Run: `cargo run --release --example frame_rank -- [root=results] [frames=24] [quota=8]`

use prin_rs::camera::Camera;
use prin_rs::grid::Slice;
use prin_rs::quad::{Decision, QuadTree};
use prin_rs::scheduler::{Mode, Policy, SchedCfg};
use prin_rs::session::{FrameQuota, Regrow, Session, SessionCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// Decisions keyed by box, so two trees of different sizes are comparable.
fn by_box(t: &QuadTree) -> std::collections::HashMap<(u32, i64, i64), Decision> {
    t.nodes
        .iter()
        .filter(|q| q.is_leaf() && !q.merged)
        .map(|q| {
            ((q.level, (q.cx / 1e-12).round() as i64, (q.cy / 1e-12).round() as i64), q.decision)
        })
        .collect()
}

fn depth_variance(t: &QuadTree) -> f64 {
    let d: Vec<f64> = t.leaves().map(|i| t.nodes[i].level as f64).collect();
    if d.is_empty() {
        return f64::NAN;
    }
    let m = d.iter().sum::<f64>() / d.len() as f64;
    d.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / d.len() as f64
}

type Boxes = std::collections::HashMap<(u32, i64, i64), Decision>;

struct Arm {
    name: &'static str,
    leaves: usize,
    quads: usize,
    dvar: f64,
    scan: f64,
    len: f64,
    /// Frames whose quota actually bound. Printed, because a truncated-regime comparison run at a
    /// non-binding quota is the same regime measured twice — the `camera_probe` failure exactly.
    bound: usize,
    boxes: Boxes,
    stops: std::collections::BTreeMap<&'static str, usize>,
    /// The tree after running to drain, and the frames it took.
    drained: Boxes,
    drain_leaves: usize,
    drain_frames: usize,
}

fn run<F>(name: &'static str, rank: bool, bias: Option<f64>, frames: usize, quota: usize, f: &F) -> Arm
where
    F: Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync,
{
    // Zoomed to a corner, so `relevance` actually varies. `Camera::framing` would set
    // `half_world` to the root half and make the camera arm an identity -- the fixture defect
    // this project has now met four times.
    let base = Camera::framing(0.0, 0.0, 1.0, 512);
    let cam = Camera { cx: 0.4, cy: 0.4, half_world: 0.35, ..base };
    let cfg = SessionCfg {
        sched: SchedCfg {
            n: 8,
            budget: 20000,
            tau_display: 1e-2,
            policy: Policy::Tolerance,
            mode: Mode::Balanced,
            max_level: Some(6),
            camera: Some(cam),
            camera_bias: bias,
            ..Default::default()
        },
        quota: FrameQuota { quads: quota, substeps: None, rounds: 2 },
        rank_frame: rank,
        regrow: Regrow::Off,
        audit_every: 4,
        ..Default::default()
    };
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg, 13.0);
    let (mut scan, mut len, mut n, mut bound) = (0.0, 0.0, 0.0, 0usize);
    for _ in 0..frames {
        let (hit, _) = s.step(f);
        bound += usize::from(hit != prin_rs::session::QuotaHit::Drained);
        let (sc, ln, ag) = s.frontier_telemetry();
        // The audit must never disagree, and it is asserted here rather than printed: a silent
        // disagreement is the invisible-staleness mode the rebuild path exists to catch.
        assert!(ag.is_nan() || ag == 1.0, "the frontier disagreed with its rebuild");
        if ln > 0 {
            scan += sc as f64;
            len += ln as f64;
            n += 1.0;
        }
    }
    let truncated = by_box(s.tree());
    let (leaves, quads, dvar) =
        (s.tree().leaves().count(), s.stats().quads_computed, depth_variance(s.tree()));
    // **Never a leaf count without its stop-reason breakdown.** A tree that a cap or a veto
    // completed is a fact about the cap, and a criterion comparison run on one is comparing two
    // schedules that both had nothing to decide.
    let mut stops: std::collections::BTreeMap<&'static str, usize> = Default::default();
    for i in s.tree().leaves() {
        *stops.entry(s.tree().nodes[i].decision.name()).or_default() += 1;
    }

    // Now run the SAME session on to drain. The destination must not depend on the schedule.
    let mut extra = 0usize;
    while extra < 4000 {
        let (hit, sp) = s.step(f);
        extra += 1;
        if hit == prin_rs::session::QuotaHit::Drained && sp.quads == 0 {
            break;
        }
    }
    Arm {
        name,
        leaves,
        quads,
        dvar,
        scan: if n > 0.0 { scan / n } else { f64::NAN },
        len: if n > 0.0 { len / n } else { f64::NAN },
        bound,
        boxes: truncated,
        stops,
        drained: by_box(s.tree()),
        drain_leaves: s.tree().leaves().count(),
        drain_frames: frames + extra,
    }
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let frames: usize = arg(2, 24);
    let quota: usize = arg(3, 8);
    let dir = format!("{root}/session");
    let _ = std::fs::create_dir_all(&dir);

    println!("frame_rank: does ranking the frame frontier move the tree?");
    println!("  frames={frames} quota={quota} quads/frame, rounds=2, budget non-binding\n");

    // Two analytic fields. `step` is a single discontinuity -- a filament through an otherwise
    // resolved plane, where the physics ordering has something to say. `filament_through_sea` is
    // the harder one: structure inside noise, where a mis-ranking is invisible in a leaf count.
    for (case, f) in [
        ("step", Box::new(prin_rs::testing::step(0.137, 13.0))
            as Box<dyn Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync>),
        ("filament_through_sea", Box::new(prin_rs::testing::filament_through_sea(0.137, 0xC0FFEE, 13.0))),
    ] {
        let arms = [
            run("off", false, None, frames, quota, &f),
            run("rank", true, None, frames, quota, &f),
            run("rank+cam", true, Some(0.5), frames, quota, &f),
        ];
        println!("== {case}");
        println!("{:>10} {:>7} {:>7} {:>8} {:>6} {:>6} {:>9} {:>8} {:>9} {:>8} {:>8}",
                 "arm", "leaves", "quads", "dvar", "scan", "len", "scan/len", "bound",
                 "only-off", "only-arm", "moved");
        for a in &arms {
            // Three separate counts, because they are three different statements. A box in one
            // tree and not the other is a SIZE change; a decision differing on a box both hold is
            // a RANKING change; and pooling them reports a growing tree as a moving one.
            let only_off = arms[0].boxes.keys().filter(|k| !a.boxes.contains_key(*k)).count();
            let only_arm = a.boxes.keys().filter(|k| !arms[0].boxes.contains_key(*k)).count();
            let moved = a
                .boxes
                .iter()
                .filter(|(k, d)| arms[0].boxes.get(*k).is_some_and(|o| o != *d))
                .count();
            let (a1, a2, a3) = if a.name == "off" {
                ("-".into(), "-".into(), "-".to_string())
            } else {
                (only_off.to_string(), only_arm.to_string(), moved.to_string())
            };
            println!("{:>10} {:>7} {:>7} {:>8.4} {:>6.1} {:>6.1} {:>9.4} {:>8} {:>9} {:>8} {:>8}",
                     a.name, a.leaves, a.quads, a.dvar, a.scan, a.len,
                     a.scan / a.len, a.bound, a1, a2, a3);
        }
        for a in &arms {
            let b: Vec<String> =
                a.stops.iter().map(|(k, v)| format!("{k}:{v}")).collect();
            println!("  stops {:>8}: {}", a.name, b.join(" "));
        }
        for a in &arms[1..] {
            println!("  truncated: {:>8} {} the tree", a.name,
                     if a.boxes == arms[0].boxes { "does NOT move" } else { "MOVES" });
        }
        // **The correctness arm.** Run to drain, the ranking has only reordered work that all
        // happened, so the destination must be identical. A difference here would mean the
        // ranking is changing the criterion rather than the schedule.
        for a in &arms[1..] {
            assert_eq!(a.drained, arms[0].drained,
                       "{}: the DRAINED tree differs from the control -- the ranking changed the \
                        destination, not the schedule", a.name);
        }
        println!("  drained:   all three arms reach the SAME {} leaves ({}-{} frames)",
                 arms[0].drain_leaves,
                 arms.iter().map(|a| a.drain_frames).min().unwrap(),
                 arms.iter().map(|a| a.drain_frames).max().unwrap());
        println!();
    }

    println!("Read `scan/len` beside `vs off`: a ranking that walks its whole frontier has bought");
    println!("nothing over a sort whatever it did to the tree, and one that moves no decision has");
    println!("bought nothing at all. Both failures are on this project's record.");
    println!("\n{dir}/output/frame_rank.txt");
}
