//! §7 — the slippy map: what a pan does, and what it would cost.
//!
//! # The structural answer, before any number
//!
//! `Camera::veto` reads `tile_size_px`, which depends on the quad's width and the camera's
//! `half_world` and `viewport` — **and not on `cx`/`cy`**. So in the current design **panning
//! changes no scheduling decision at all.** There is no view culling, and there is no cache to
//! invalidate: `SchedStats` is dropped when `descend` returns.
//!
//! That makes three of §7's four questions identities rather than measurements, and reporting
//! "the tree persists perfectly across a pan" as a finding would be reporting an identity. So
//! this measures the two things that are actually open:
//!
//! 1. **What would be evictable**, via `Camera::covers` — a pure predicate that the scheduler
//!    does not consult. Nothing is evicted. Adding culling would be a caching decision and this
//!    build has no cache, by scope.
//! 2. **Whether a re-entering quad comes back the same.** The scheduler is deterministic and
//!    the ensemble offsets are a fixed Halton prefix indexed by copy, not by pixel or by time,
//!    so a quad recomputed after leaving view must be **bitwise** what it was. That is a claim
//!    that can fail, and it is what makes a future cache sound: a cached reduction and a
//!    recomputed one have to be interchangeable.
//!
//! # How to misread this
//!
//! **`recompute%` is 100% by construction and is not a criticism of the criterion.** Nothing is
//! cached, so everything in view is recomputed every frame. The column worth reading is
//! `would-cache%`: the fraction that a cache *would* have served, which is what a caching design
//! would be buying.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use std::collections::{HashMap, HashSet};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// A quad's identity across runs: level and box centre, quantised so float formatting cannot
/// make the same box look like two.
fn key(level: u32, cx: f64, cy: f64) -> (u32, i64, i64) {
    (level, (cx * 1e12).round() as i64, (cy * 1e12).round() as i64)
}

fn main() {
    let steps: usize = arg(1, 9);
    let budget: usize = arg(2, 2000);
    let viewport: usize = arg(3, 512);
    let region: String = std::env::args().nth(4).unwrap_or_else(|| "near-field".into());

    let s = grid::region(&region, 2, 2, 0.05).unwrap();
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let half_world = 0.05;

    println!(
        "pan across {region}, {steps} steps, budget {budget}, viewport {viewport}^2.\n\
         Camera::veto does not read cx/cy, so panning changes NO scheduling decision here.\n\
         `covers` is a predicate the scheduler never consults; nothing is evicted.\n"
    );
    println!(
        "{:>5} {:>9} {:>7} {:>7} {:>9} {:>12} {:>11} {:>9}",
        "step", "cam cx", "quads", "in view", "new", "would-cache%", "recompute%", "floored"
    );

    // Frames of the adaptive render at each camera position, so the pan is watchable and the
    // "nothing changes" result is visible rather than only tabulated.
    let mut frames: Vec<Vec<u8>> = Vec::new();
    let mut wire_frames: Vec<Vec<u8>> = Vec::new();
    // Follows the viewport. It was hardcoded at 384 while `viewport` set only the camera, so
    // asking for a larger pan rendered the same small raster -- and a small raster upscaled by a
    // viewer reads as a blurry render rather than as a stale size.
    let frame_res = viewport;

    let mut prev: HashSet<(u32, i64, i64)> = HashSet::new();
    // Every quad ever computed, with its reduction, so a re-entering quad can be compared with
    // what it was the first time it was seen.
    let mut seen: HashMap<(u32, i64, i64), (f64, f64, u8)> = HashMap::new();
    let mut mismatches = 0usize;
    let mut rechecks = 0usize;
    let mut floor_flips = 0usize;

    for k in 0..steps {
        // Pan across one full viewport width over the sequence.
        let frac = (k as f64 / (steps - 1).max(1) as f64) - 0.5;
        let cam_cx = s.cx + 2.0 * half_world * frac;
        let cam = Camera { cx: cam_cx, cy: s.cy, half_world, viewport, max_rel_depth: None };

        let cfg = SchedCfg {
            budget,
            tau_display: 1e-4,
            alpha_hi: 0.2,
            alpha_lo: 0.2,
            camera: Some(cam),
            keep_pixels: true,
            ..Default::default()
        };
        let (t, st) = scheduler::descend(s.cx, s.cy, half_world, s.body, &cfg, &ens, Precision::F64);

        // Adaptive, at true per-quad texel sizes. A uniform render would show where boundaries
        // fell rather than what the system displays, which is the instrument this project has
        // already been caught using once.
        let (img, _tex) = prin_rs::output::adaptive::render(
            &t,
            &st.pixels,
            &cam,
            frame_res,
            prin_rs::output::adaptive::TexelMode::Adaptive,
            prin_rs::output::png::outcome_rgb,
        );
        let stem = format!("results/criterion/pan_{}_{k:02}", region.replace(' ', "_"));
        let _ = prin_rs::output::adaptive::save(&format!("{stem}.png"), frame_res, &img);
        let mut wimg = img.clone();
        {
            let boxes = prin_rs::output::wire::boxes_from_tree(&t, &cam, frame_res);
            let deepest = boxes.iter().map(|b| b.level).max().unwrap_or(1);
            prin_rs::output::wire::draw(&mut wimg, frame_res, frame_res, &boxes, deepest.max(1));
        }
        let _ = prin_rs::output::adaptive::save(&format!("{stem}_wire.png"), frame_res, &wimg);
        wire_frames.push(wimg);
        if let Ok(f) = std::fs::File::create(format!("{stem}.prnq")) {
            let mut w = std::io::BufWriter::new(f);
            let _ = prin_rs::output::tree::write(&mut w, &t, &cfg, &ens, &st, &region, "f64");
        }
        frames.push(img);

        let computed: Vec<usize> =
            (0..t.nodes.len()).filter(|&i| t.nodes[i].red.n_footprints > 0).collect();
        let cur: HashSet<_> = computed.iter().map(|&i| {
            let q = &t.nodes[i];
            key(q.level, q.cx, q.cy)
        }).collect();

        let in_view = computed
            .iter()
            .filter(|&&i| {
                let q = &t.nodes[i];
                cam.covers(q.cx, q.cy, q.half)
            })
            .count();
        let new = cur.difference(&prev).count();
        let would_cache = cur.intersection(&prev).count();
        let floored = computed
            .iter()
            .filter(|&&i| t.nodes[i].decision == Decision::ScreenFloor)
            .count();

        // Determinism: a quad seen before must come back bitwise identical, or a cache would be
        // unsound. Compared on the reduction AND the decision, because a decision that moved
        // while the payload did not would be scheduler state leaking into the sim key.
        for &i in &computed {
            let q = &t.nodes[i];
            let kk = key(q.level, q.cx, q.cy);
            let now = (q.red.spread_median, q.red.between_spread, q.decision.code());
            match seen.get(&kk) {
                Some(&was) => {
                    rechecks += 1;
                    if was.0.to_bits() != now.0.to_bits() || was.1.to_bits() != now.1.to_bits() {
                        mismatches += 1;
                    }
                    if was.2 != now.2 {
                        floor_flips += 1;
                    }
                }
                None => {
                    seen.insert(kk, now);
                }
            }
        }

        println!(
            "{k:>5} {cam_cx:>9.5} {:>7} {in_view:>7} {new:>9} {:>11.1}% {:>10.1}% {floored:>9}",
            computed.len(),
            100.0 * would_cache as f64 / computed.len().max(1) as f64,
            100.0,
        );
        prev = cur;
    }

    let _ = prin_rs::output::apng::write(
        &format!("results/animated/pan_{}_animated.png", region.replace(' ', "_")),
        frame_res,
        frame_res,
        &frames,
        1,
        3,
    );
    let _ = prin_rs::output::apng::write(
        &format!("results/animated/pan_{}_wire_animated.png", region.replace(' ', "_")),
        frame_res,
        frame_res,
        &wire_frames,
        1,
        3,
    );

    println!(
        "\n{rechecks} quads recomputed after having been seen before.\n\
         {mismatches} came back with a DIFFERENT reduction. {floor_flips} with a different decision.\n"
    );
    if mismatches == 0 {
        println!(
            "Zero payload mismatches is the property a cache would need: a recomputed reduction\n\
             and a cached one are interchangeable. It holds because the ensemble offsets are a\n\
             fixed Halton prefix indexed by COPY -- not by pixel, not by camera, not by time --\n\
             so a quad's ensemble does not know the camera exists."
        );
    } else {
        println!("NON-DETERMINISM: a cache would be unsound. Investigate before caching anything.");
    }
    println!(
        "\n`recompute%` is 100 by construction: nothing is cached, so everything is recomputed\n\
         every frame. `would-cache%` is what a cache would have served, and is the number a\n\
         caching design would be buying. No eviction policy is implemented and none is implied.\n\
         \n\
         `in view` uses Camera::covers, which the scheduler never consults. If it ever does,\n\
         that is view culling and a scope change -- the screen floor is a veto on SCALE, and\n\
         adding a position term to it would make a quad's decision depend on where the camera\n\
         is pointing, which is the thing `never cached as a quad fact` exists to keep out."
    );
}
