//! **§18: does foveating the frontier toward the pointer beat not doing it?**
//!
//! `Camera::foveation` **modulates camera relevance** — it is not a third factor, it lives in
//! priority and never in the veto, and no cursor field goes on a `Quad`. With no cursor, or a
//! cursor still in motion (`dwell = 0`), it returns exactly `1.0` and `priority` is §4.3's product
//! unchanged: **the fallback is the default path**, which every keyboard, touch, unfocused-window
//! and headless route takes.
//!
//! It ships **off**, and this is the measurement that decides whether it ever ships on. The plan is
//! explicit and so is this harness: *if it does not beat uniform on time-to-resolve-at-cursor by a
//! clear margin it is dropped.* It costs a factor on the hot path and a mode in the telemetry, and
//! "it might be cool" is not a measurement.
//!
//! # The obvious way for this to be vacuous, named before it is run
//!
//! **A cursor parked at the frame centre makes the fovea concentric with the viewport**, so the
//! bias is very nearly `relevance` itself and the foveated arm reproduces the unfoveated one — a
//! null that reads as "no effect" and is an *identity*. So `centre` is run as a named arm and is
//! expected to be near-inert; the arm that can show anything is `off-centre` with `dwell = 1`.
//! Same shape as `a_pan_is_an_identity_until_camera_bias_is_switched_on`, already in the suite for
//! exactly this reason.
//!
//! # What is measured, and the cost side is measured in the same run
//!
//! `resolved(region, frame)` is the share of leaves overlapping a disc whose decision is not
//! `Split`, `Pending` or `Deferred` — i.e. the descent there has settled. **Time-to-resolve** is
//! the first frame at which that reaches 1.0.
//!
//! Two discs: one at the **cursor** (the gain) and one at the **frame edge** (the cost). A
//! foveation that resolves the cursor sooner by starving the edge has not improved anything; it
//! has moved the same budget. Both columns are printed side by side, and `edge` is what says which
//! happened.
//!
//! Run: `cargo run --release --example cursor_bias -- [root=results] [frames=40] [quota=12]`

use prin_rs::camera::{Camera, Cursor};
use prin_rs::grid::Slice;
use prin_rs::quad::{Decision, QuadTree};
use prin_rs::scheduler::{Mode, Policy, SchedCfg};
use prin_rs::session::{FrameQuota, Session, SessionCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// Share of leaves overlapping the disc that have settled, and how many there were.
///
/// `NaN` on an empty disc, never `1.0`: a region with no leaves has not resolved, and reporting a
/// pass there is the *"a check that did not run must not report a pass"* conflation.
fn resolved(t: &QuadTree, cx: f64, cy: f64, r: f64) -> (f64, usize) {
    let mut n = 0usize;
    let mut done = 0usize;
    for i in t.leaves() {
        let q = &t.nodes[i];
        // Box-disc overlap, by the nearest point of the box to the centre.
        let dx = (cx - q.cx).abs() - q.half;
        let dy = (cy - q.cy).abs() - q.half;
        let d2 = dx.max(0.0).powi(2) + dy.max(0.0).powi(2);
        if d2 > r * r {
            continue;
        }
        n += 1;
        done += usize::from(!matches!(
            q.decision,
            Decision::Split | Decision::Pending | Decision::Deferred
        ));
    }
    if n == 0 { (f64::NAN, 0) } else { (done as f64 / n as f64, n) }
}

struct Arm {
    name: &'static str,
    t_cursor: Option<usize>,
    t_edge: Option<usize>,
    leaves: usize,
    quads: usize,
    depth_cursor: f64,
    depth_edge: f64,
    /// Leaves whose decision differs from the `off` arm's, over boxes both hold.
    moved: usize,
    shared: usize,
    boxes: std::collections::HashMap<(u32, i64, i64), Decision>,
}

fn mean_depth(t: &QuadTree, cx: f64, cy: f64, r: f64) -> f64 {
    let d: Vec<f64> = t
        .leaves()
        .filter(|&i| {
            let q = &t.nodes[i];
            let dx = (cx - q.cx).abs() - q.half;
            let dy = (cy - q.cy).abs() - q.half;
            dx.max(0.0).powi(2) + dy.max(0.0).powi(2) <= r * r
        })
        .map(|i| t.nodes[i].level as f64)
        .collect();
    if d.is_empty() { f64::NAN } else { d.iter().sum::<f64>() / d.len() as f64 }
}

#[allow(clippy::too_many_arguments)]
fn run<F>(
    name: &'static str,
    cursor: Option<Cursor>,
    cap: f64,
    frames: usize,
    quota: usize,
    f: &F,
    probe: (f64, f64),
    edge: (f64, f64),
    r: f64,
) -> Arm
where
    F: Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync,
{
    let base = Camera::framing(0.0, 0.0, 1.0, 512);
    // Zoomed, so `relevance` varies at all -- `Camera::framing` sets `half_world` to the root half
    // and makes every quad fully visible, which would make the whole comparison an identity.
    let cam = Camera { cx: 0.0, cy: 0.0, half_world: 0.6, ..base };
    let cfg = SessionCfg {
        sched: SchedCfg {
            n: 8,
            budget: 20000,
            tau_display: 1e-2,
            policy: Policy::Tolerance,
            mode: Mode::Balanced,
            max_level: Some(6),
            camera: Some(cam),
            camera_bias: Some(0.5),
            cursor,
            fovea_cap: cap,
            ..Default::default()
        },
        quota: FrameQuota { quads: quota, substeps: None, rounds: 2 },
        audit_every: 8,
        ..Default::default()
    };
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg, 13.0);
    let (mut t_cursor, mut t_edge) = (None, None);
    for fr in 0..frames {
        s.step(f);
        let (rc, nc) = resolved(s.tree(), probe.0, probe.1, r);
        let (re, ne) = resolved(s.tree(), edge.0, edge.1, r);
        if t_cursor.is_none() && nc > 0 && rc >= 1.0 {
            t_cursor = Some(fr);
        }
        if t_edge.is_none() && ne > 0 && re >= 1.0 {
            t_edge = Some(fr);
        }
    }
    let boxes = s
        .tree()
        .nodes
        .iter()
        .filter(|q| q.is_leaf() && !q.merged)
        .map(|q| {
            ((q.level, (q.cx / 1e-12).round() as i64, (q.cy / 1e-12).round() as i64), q.decision)
        })
        .collect();
    Arm {
        name,
        t_cursor,
        t_edge,
        leaves: s.tree().leaves().count(),
        quads: s.stats().quads_computed,
        depth_cursor: mean_depth(s.tree(), probe.0, probe.1, r),
        depth_edge: mean_depth(s.tree(), edge.0, edge.1, r),
        moved: 0,
        shared: 0,
        boxes,
    }
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let frames: usize = arg(2, 40);
    let quota: usize = arg(3, 12);
    let dir = format!("{root}/cursor");
    let _ = std::fs::create_dir_all(format!("{dir}/output"));

    // **BOTH probes sit ON the structure.** The first cut put the edge probe at an empty corner,
    // where `step`'s field is smooth: it resolved at frame 1 in every arm at mean depth 2.000, so
    // the cost side had no subject and a foveation that starved the edge would have looked free.
    // `step(0.137)` is a vertical discontinuity at `x = 0.137`, so both discs are placed on that
    // line, far apart in `y` -- the cursor near the top, the probe near the bottom.
    let (probe, edge, r) = ((0.137, 0.40), (0.137, -0.40), 0.12);

    println!("cursor_bias: §18, foveation against its own off-state");
    println!("  frames={frames} quota={quota}, camera half_world=0.6 at the origin");
    println!("  cursor probe at {probe:?}, edge probe at {edge:?}, disc r={r}\n");
    println!("  time-to-resolve is the FIRST frame the disc's leaves have all settled;");
    println!("  `-` means it never did within {frames} frames.\n");

    for (case, f) in [
        ("step", Box::new(prin_rs::testing::step(0.137, 13.0))
            as Box<dyn Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync>),
        ("filament_through_sea",
         Box::new(prin_rs::testing::filament_through_sea(0.137, 0xC0FFEE, 13.0))),
    ] {
        let arms = vec![
            run("off", None, 4.0, frames, quota, &f, probe, edge, r),
            // **The named vacuous cell.** A cursor at the frame centre is concentric with the
            // viewport, so the fovea is nearly `relevance` itself; this arm is expected to be
            // near-inert and is what says the off-centre arm is not measuring an identity.
            run("centre", Some(Cursor { cx: 0.0, cy: 0.0, dwell: 1.0 }), 4.0,
                frames, quota, &f, probe, edge, r),
            // **Moving, so `dwell = 0`.** `foveation` returns exactly 1.0 there, so this must
            // reproduce `off` BITWISE -- the fallback-is-the-default-path claim, as a test.
            run("moving", Some(Cursor { cx: probe.0, cy: probe.1, dwell: 0.0 }), 4.0,
                frames, quota, &f, probe, edge, r),
            run("fovea x4", Some(Cursor { cx: probe.0, cy: probe.1, dwell: 1.0 }), 4.0,
                frames, quota, &f, probe, edge, r),
            run("fovea x16", Some(Cursor { cx: probe.0, cy: probe.1, dwell: 1.0 }), 16.0,
                frames, quota, &f, probe, edge, r),
            // **Dwell is the low-pass and it is a CONTINUOUS knob**, so a half-settled pointer must
            // land between `moving` and `fovea x4` or the weighting is a switch wearing a float's
            // type. This is the arm that says so.
            run("dwell 0.5", Some(Cursor { cx: probe.0, cy: probe.1, dwell: 0.5 }), 4.0,
                frames, quota, &f, probe, edge, r),
            // **THE CONTROL THAT SAYS THE EFFECT IS THE CURSOR AND NOT THE PLACE.** The two probes
            // are geometrically symmetric and their baselines are NOT: the tie-break is
            // lexicographic on `(level, ix, iy)`, so the lower-`y` probe is reached first and the
            // `off` arm resolves it two frames sooner. So a gain measured at one probe could be
            // the scan order rather than the fovea. Putting the cursor on the OTHER probe must
            // move the advantage with it -- and the columns are reported unswapped, so `t_edge`
            // is the foveated one in this row.
            run("x16 @edge", Some(Cursor { cx: edge.0, cy: edge.1, dwell: 1.0 }), 16.0,
                frames, quota, &f, probe, edge, r),
        ];

        println!("== {case}");
        println!("{:>10} {:>9} {:>9} {:>8} {:>8} {:>9} {:>9} {:>12}",
                 "arm", "t_cursor", "t_edge", "leaves", "quads", "d_cursor", "d_edge", "vs off");
        for a in &arms {
            let shared: usize =
                a.boxes.keys().filter(|k| arms[0].boxes.contains_key(*k)).count();
            let moved = a
                .boxes
                .iter()
                .filter(|(k, d)| arms[0].boxes.get(*k).is_some_and(|o| o != *d))
                .count();
            let only = a.boxes.keys().filter(|k| !arms[0].boxes.contains_key(*k)).count();
            let vs = if a.name == "off" {
                "-".to_string()
            } else {
                format!("{moved}m {only}new/{shared}")
            };
            let tc = a.t_cursor.map_or("-".into(), |x| x.to_string());
            let te = a.t_edge.map_or("-".into(), |x| x.to_string());
            println!("{:>10} {tc:>9} {te:>9} {:>8} {:>8} {:>9.3} {:>9.3} {vs:>12}",
                     a.name, a.leaves, a.quads, a.depth_cursor, a.depth_edge);
        }

        // **The fallback claim, asserted rather than printed.** A cursor in motion has
        // `dwell = 0`, `foveation` returns exactly 1.0, and the tree must be the `off` tree.
        assert_eq!(arms[2].boxes, arms[0].boxes,
                   "a MOVING cursor changed the tree -- `foveation` is not returning 1.0 at \
                    dwell = 0, and the fallback is not the default path");
        println!("  moving cursor (dwell = 0) reproduces `off` bitwise: the fallback IS the \
default path");
        let centre_inert = arms[1].boxes == arms[0].boxes;
        println!("  centre cursor {} the tree -- the named vacuous cell",
                 if centre_inert { "does NOT move" } else { "MOVES" });
        // **`fovea_cap` saturates, and the reason is arithmetic.** At this separation the
        // Gaussian weight at the edge probe is ~1e-13, so `foveation` is essentially `1/cap` --
        // 0.25 at x4 and 0.0625 at x16. Both put the periphery LAST in the ranking, and once it is
        // last, demoting it further cannot change an order. So the cap is a rank knob with a
        // ceiling, not a dial, and a table showing x4 and x16 identical is that and not a bug.
        for a in &arms[3..5] {
            let gain = match (a.t_cursor, arms[0].t_cursor) {
                (Some(x), Some(y)) => format!("{}", y as i64 - x as i64),
                _ => "n/a".into(),
            };
            let cost = match (a.t_edge, arms[0].t_edge) {
                (Some(x), Some(y)) => format!("{}", x as i64 - y as i64),
                (None, Some(_)) => "edge NEVER resolved".into(),
                _ => "n/a".into(),
            };
            println!("  {:>9}: frames saved at the cursor = {gain}, frames LOST at the edge = {cost}",
                     a.name);
        }
        // The swapped control, read the other way round: the fovea is on `edge` here, so the
        // advantage must appear in the `t_edge` column and the disadvantage in `t_cursor`.
        let sw = arms.last().unwrap();
        match (sw.t_edge, arms[0].t_edge, sw.t_cursor, arms[0].t_cursor) {
            (Some(a), Some(b), Some(c), Some(d)) => println!(
                "  x16 @edge: with the fovea moved ONTO the edge probe, edge {} -> {} \
({:+}) and cursor {} -> {} ({:+}). The advantage follows the CURSOR if the first \
number improves.", b, a, b as i64 - a as i64, d, c, c as i64 - d as i64),
            _ => println!("  x16 @edge: one probe never resolved -- no comparison"),
        }
        println!();
    }

    println!("THE DECISION RULE, stated before the numbers were seen: foveation ships only if it");
    println!("beats its own off-state on time-to-resolve at the cursor BY A CLEAR MARGIN, and the");
    println!("edge column is what says whether a saving is a gain or the same budget moved.");
    println!("\n{dir}/output/cursor_bias.txt");
}
