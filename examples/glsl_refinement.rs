//! **The four GLSL preset charts, refining over time.** One animation each, nothing else.
//!
//! These are the reference's own default slices (`Ma1achy/principia-ii`, `src/state.ts:71-76`) at
//! `z0 = 0`, which decodes to the equilateral Lagrange configuration — the four pictures a person
//! can recognise rather than tabulate. The window is `half = 3.0`, the reference UI's
//! `Slice +/- 3.0e+0`.
//!
//! # Frames are quads, not levels
//!
//! `results/animated/<case>_levels.png` steps one level per frame, so a chart that reaches depth 6
//! is a six-frame animation — too few to read as motion. Here a frame is emitted every `n` quads
//! **in the order the scheduler computed them**, so the picture sharpens continuously and the
//! frame count is a parameter rather than an accident of the tree's depth.
//!
//! One descent per chart. Every `Quad` records the `iteration` it was computed in, so the whole
//! sequence is reconstructed from the finished tree — the animation costs a raster per frame and
//! no extra physics.

use std::collections::HashSet;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, wire};
use prin_rs::quad::{Criterion, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Mode, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The leaf set once the first `n` computed quads exist.
///
/// A revealed node is drawn as a leaf when none of its children have been revealed yet — which is
/// exactly what the scheduler was displaying at that moment.
fn leaves_after(t: &QuadTree, order: &[usize], n: usize) -> Vec<usize> {
    let live: HashSet<usize> = order.iter().take(n).cloned().collect();
    order
        .iter()
        .take(n)
        .filter(|&&i| match t.nodes[i].children {
            None => true,
            Some(k) => !live.contains(&k[0]),
        })
        .cloned()
        .collect()
}

fn main() {
    let budget: usize = arg(1, 40000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let res: usize = arg(4, 1024);
    let frames_wanted: usize = arg(5, 72);
    // **The configuration the sweep found, not the defaults.** `k_frac` is what makes the
    // criterion a priority ordering rather than a gate; at 1.0 the ranking takes the whole
    // frontier and changes nothing, which is the defect that made every dump in PR #18 a
    // pre-fix run. RESULTS.md §18: at `k = 0.25` near-field's depth variance doubles and its
    // veto share falls 61% -> 13%.
    let k_frac: f64 = arg(6, 0.25);
    // `grad_rms` is the criterion that unlocks `preset_shape` -- the only one of these four
    // that `within` cannot move at ANY tau, k_frac or alpha, including a gate-off control.
    // Measured: 16 leaves at one level under `within`, 31 at five under `grad_rms`.
    let crit = match std::env::args().nth(7).as_deref() {
        Some(c) => Criterion::parse(c).expect("criterion"),
        None => Criterion::GradRms,
    };

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let dir = "results/glsl";
    let _ = std::fs::create_dir_all(dir);

    let half = Chart::preset_shape().default_half();
    let cases: [(&str, Chart); 4] = [
        ("shape", Chart::preset_shape()),
        ("prho", Chart::preset_prho()),
        ("plambda", Chart::preset_plambda()),
        ("shape_pl", Chart::preset_shape_pl()),
    ];

    println!("the four GLSL preset slices, refining. {res}^2, half {half}, budget {budget}, \
              tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, E+1={}, t={}, f64",
             ens.n_extra + 1, ens.t_max);
    println!("criterion={}, k_frac={k_frac}, mode=balanced -- the ranked frontier, which is the \
              mechanism.", crit.name());
    if k_frac >= 1.0 {
        println!("WARNING: k_frac = 1 takes the top 100% of the frontier, so the ranking runs and");
        println!("changes nothing. That is the PRE-FIX configuration, not the new system.");
    }
    println!("{:>10} {:>7} {:>7} {:>6} {:>7} {:>8}", "case", "quads", "leaves", "depth", "frames",
             "wall s");

    for (name, chart) in cases {
        let t0 = std::time::Instant::now();
        let cam = Camera::framing(0.0, 0.0, half, res);
        let cfg = SchedCfg {
            budget,
            tau_display: tau,
            alpha_hi,
            alpha_lo: alpha_hi,
            // The criterion measured best in §16, and the whole frontier each round so the tree
            // fills out rather than being thinned by a partial budget -- this is a picture of
            // refinement, not of the demotion mechanism.
            criterion: crit,
            mode: Mode::Balanced,
            k_frac,
            camera: Some(cam),
            keep_pixels: true,
            chart,
            ..Default::default()
        };
        let (t, st) = scheduler::descend(0.0, 0.0, half, 0, &cfg, &ens, Precision::F64);

        // Computed order: the scheduler's own, round by round.
        let mut order: Vec<usize> =
            (0..t.nodes.len()).filter(|&i| t.nodes[i].red.n_footprints > 0).collect();
        order.sort_by_key(|&i| (t.nodes[i].iteration, i));

        // One ramp for the whole chart, from the finished tree. Re-ranging per frame would
        // animate the ramp rather than the refinement.
        let all_px: Vec<PixelOut> =
            t.leaves().flat_map(|i| st.pixels.get(i).cloned().unwrap_or_default()).collect();
        let (lo, hi) = colour::range(&all_px, Scalar::ShapeSpread);
        let sites = colour::landmarks(&grid::decode_state(&chart, 0, 0.0, 0.0).m);
        let rgb = |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);

        let step = (order.len() / frames_wanted).max(1);
        let mut frames = Vec::new();
        let mut wframes = Vec::new();
        let mut n = 1usize;
        loop {
            let lv = leaves_after(&t, &order, n.min(order.len()));
            let mut shadow = t.clone();
            let keep: HashSet<usize> = lv.iter().cloned().collect();
            for i in 0..shadow.nodes.len() {
                if keep.contains(&i) {
                    shadow.nodes[i].children = None;
                }
            }
            let f = adaptive::render(
                &shadow, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, &rgb,
            )
            .0;
            let mut wf = f.clone();
            wire::draw(&mut wf, res, res, &wire::boxes_from_tree(&shadow, &cam, res), 1);
            frames.push(f);
            wframes.push(wf);
            if n >= order.len() {
                break;
            }
            n = (n + step).min(order.len());
        }
        // Hold the finished frame, so the loop reads as an ending rather than a snap back.
        for _ in 0..8 {
            frames.push(frames.last().unwrap().clone());
            wframes.push(wframes.last().unwrap().clone());
        }

        let _ = apng::write(&format!("{dir}/{name}.png"), res, res, &frames, 1, 12);
        let _ = apng::write(&format!("{dir}/{name}_wire.png"), res, res, &wframes, 1, 12);

        let leaves: Vec<usize> = t.leaves().collect();
        let depth = leaves.iter().map(|&i| t.nodes[i].level).max().unwrap_or(0);
        println!("{name:>10} {:>7} {:>7} {depth:>6} {:>7} {:>8.1}",
                 st.quads_computed, leaves.len(), frames.len(), t0.elapsed().as_secs_f64());
    }

    println!();
    println!("{dir}/<case>.png and <case>_wire.png. A frame every {} quads or so, in the order",
             "few");
    println!("the scheduler computed them, so the picture sharpens continuously.");
}
