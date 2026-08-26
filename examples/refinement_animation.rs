//! **The new refinement mechanism, animated — three views, every chart.**
//!
//! `<case>_levels.png` in `results/animated/` truncates one descent **by depth**. That was the
//! right picture for a criterion that was a *stop condition*: the tree grows level by level and
//! the animation shows how deep it got.
//!
//! This mechanism is not that. The criterion is a **priority ordering** over a ranked frontier,
//! and a depth ladder cannot show it — a quad that is refined last and a quad that is never
//! refined sit at the same depth in every frame. What is new is *which* quads get the budget and
//! *when*, so the three animations here advance the three things that actually vary.
//!
//! # 1. `<case>_budget.png` — the frontier being spent
//!
//! One frame per descent round. The tree at round `k` is reconstructed from the single completed
//! descent, because every `Quad` carries the `iteration` it was computed in — so this costs one
//! descent, not one per frame. **This is the animation the old ladder could not be**: it shows
//! the budget landing somewhere, round by round, rather than the tree deepening.
//!
//! # 2. `<case>_oldnew.png` — the shipped criterion against the measured-best one
//!
//! Two panels at the same budget, `within/median` on the left and `frac_hot_between/median` on
//! the right, both under the ranked frontier. §16 measured the second beating the random band in
//! all three targets and reaching `0.07038` against greedy's `0.06881` on `preset_shape`, while
//! the first is beaten by random at every budget — on **31 distinct values against 5418**.
//!
//! **A caveat that decides how to read every frame of this one.** On 23 of 26 charts a *camera
//! veto* stops 95%+ of leaves, so the two panels are largely showing the same cap being reached
//! by two routes. The charts where the difference is a criterion difference are the ones whose
//! veto share is low, and the printed table names them.
//!
//! # 3. `<case>_kfrac.png` — the demotion mechanism itself
//!
//! Same chart, same budget, `k_frac` stepping `0.25 → 1.0`: the fraction of the eligible frontier
//! that gets spent on each round. `1.0` reproduces the unranked descent exactly, so the last
//! frame is the control. The frames before it are quads being **outranked rather than refused** —
//! marked `Keep`, never `BudgetExhausted`, which is what keeps the two distinguishable in a dump.
//!
//! # What these are and are not
//!
//! **Diagnostics, not measurements.** They run at a smaller viewport than the committed stills so
//! the descent is a few minutes rather than an hour per chart, which means the screen floor bites
//! one level shallower and the trees are *not* the committed ones. The measurements are
//! `results/output/structure_metric.txt` and the `error(B)` curves; these say what a tree looks
//! like while it is being built. **Do not read a leaf count off a frame.**

use std::collections::HashSet;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, wire};
use prin_rs::quad::{Agg, Criterion, Decision, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Mode, SchedCfg, SchedStats};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The leaf set **as it stood at the end of round `cap`**, from a completed descent.
///
/// A node belongs to that frontier if it had been computed by then (`iteration <= cap`) and had
/// no children computed by then. One descent, every frame — the same trick the depth ladder uses
/// on `level`, applied to the axis that actually varies here.
fn leaves_at_round(t: &QuadTree, cap: u32) -> Vec<usize> {
    (0..t.nodes.len())
        .filter(|&i| {
            let q = &t.nodes[i];
            if q.iteration > cap {
                return false;
            }
            match q.children {
                None => true,
                Some(k) => t.nodes[k[0]].iteration > cap,
            }
        })
        .collect()
}

/// A tree whose leaf set is exactly `leaves`.
///
/// Built rather than filtered because `wire::Box2` carries only geometry and a level, not the
/// node it came from — so a truncated wireframe has to come from a truncated *tree*. Which is
/// also the honest way round: the wire must describe the same tree the colour frame does, and
/// deriving both from one shadow makes that structural rather than a thing to remember.
fn shadow_of(t: &QuadTree, leaves: &[usize]) -> QuadTree {
    let mut shadow = t.clone();
    let keep: HashSet<usize> = leaves.iter().cloned().collect();
    for i in 0..shadow.nodes.len() {
        if keep.contains(&i) {
            shadow.nodes[i].children = None;
        }
    }
    shadow
}

/// Samples for the revealed set only.
///
/// **The shadow tree alone does not truncate the render.** `adaptive::render` draws every node
/// that has samples, coarsest first, so quads outside the set paint last and win -- which made
/// every frame of every animation here the finished image, byte-identical.
fn mask(pixels: &[Vec<PixelOut>], leaves: &[usize]) -> Vec<Vec<PixelOut>> {
    let keep: HashSet<usize> = leaves.iter().cloned().collect();
    (0..pixels.len())
        .map(|i| if keep.contains(&i) { pixels[i].clone() } else { Vec::new() })
        .collect()
}

fn render_leaves(
    t: &QuadTree,
    pixels: &[Vec<PixelOut>],
    cam: &Camera,
    res: usize,
    leaves: &[usize],
    rgb: &dyn Fn(&PixelOut) -> [u8; 3],
) -> Vec<u8> {
    let shadow = shadow_of(t, leaves);
    let masked = mask(pixels, leaves);
    adaptive::render(&shadow, &masked, cam, res, adaptive::TexelMode::Adaptive, |p| rgb(p)).0
}

/// Two panels side by side with a one-pixel divider, so a frame is one image.
fn side_by_side(a: &[u8], b: &[u8], res: usize) -> Vec<u8> {
    let w = res * 2 + 1;
    let mut out = vec![90u8; w * res * 3];
    for y in 0..res {
        let dst = y * w * 3;
        out[dst..dst + res * 3].copy_from_slice(&a[y * res * 3..(y + 1) * res * 3]);
        let dst2 = dst + (res + 1) * 3;
        out[dst2..dst2 + res * 3].copy_from_slice(&b[y * res * 3..(y + 1) * res * 3]);
    }
    out
}

struct Run {
    tree: QuadTree,
    st: SchedStats,
}

fn descend(
    chart: &Chart,
    cx: f64,
    cy: f64,
    half: f64,
    res: usize,
    budget: usize,
    tau: f64,
    alpha_hi: f64,
    criterion: Criterion,
    k_frac: f64,
    ens: &EnsembleCfg,
) -> Run {
    let cfg = SchedCfg {
        budget,
        tau_display: tau,
        alpha_hi,
        alpha_lo: alpha_hi,
        criterion,
        mode: Mode::Balanced,
        k_frac,
        camera: Some(Camera::framing(cx, cy, half, res)),
        keep_pixels: true,
        chart: *chart,
        ..Default::default()
    };
    let (tree, st) = scheduler::descend(cx, cy, half, 0, &cfg, ens, Precision::F64);
    Run { tree, st }
}

fn main() {
    let budget: usize = arg(1, 40000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    // **512, not the stills' 1024, and the module header says why.** The descent cost is set by
    // the quad count, which the screen floor sets, which the viewport sets: halving it is a 4x
    // saving per descent and this example runs five of them per chart.
    let res: usize = arg(4, 512);

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let dir = "results/refinement";
    let _ = std::fs::create_dir_all(dir);

    // The measured-best criterion (§16) and the shipped one, so the pair is the real comparison
    // rather than two arbitrary settings.
    const NEW: Criterion = Criterion::FracHotBetween;
    const OLD: Criterion = Criterion::Within;
    const K_STEPS: [f64; 4] = [0.25, 0.5, 0.75, 1.0];
    const K_MAIN: f64 = 0.5;

    println!("refinement animations. budget {budget}, tau={tau:.0e}, alpha_hi={alpha_hi}, \
              N=8, E+1={}, t={}, f64, {res}^2",
             ens.n_extra + 1, ens.t_max);
    println!("new = {}/median, old = {}/median, both Mode::Balanced. k_frac main {K_MAIN}, \
              sweep {K_STEPS:?}", NEW.name(), OLD.name());
    // **Print the raster, do not name it.** A hardcoded "512" in this line would survive a run
    // at any other size and describe the wrong thing -- the defect already on record from
    // `pan_sequence`'s hardcoded viewport and `between_vs_within`'s literal 512.
    let stills = 1024usize;
    println!("DIAGNOSTICS, not measurements: the viewport is {res} against the stills' {stills},");
    println!("so the screen floor bites {} level(s) shallower and these are NOT the committed trees.",
             ((stills as f64 / res as f64).log2().max(0.0)).round() as i32);
    println!();
    println!("{:>20} {:>7} {:>7} {:>6} {:>7} {:>8} {:>9} {:>9}",
             "case", "rounds", "leaves", "depth", "veto%", "old lvs", "new lvs", "wall s");

    for (name, chart, cx, cy, half) in grid::gallery_cases() {
        let t0 = std::time::Instant::now();
        let cam = Camera::framing(cx, cy, half, res);

        // The main run: the new mechanism at the main k_frac. Every frame of animation 1 comes
        // out of this one descent.
        let new = descend(&chart, cx, cy, half, res, budget, tau, alpha_hi, NEW, K_MAIN, &ens);

        // One ramp and one site set per chart, from the pixels this tree produced. Per-frame
        // normalisation would make a quad's colour depend on which round it is being drawn in,
        // and the animation would show the ramp moving rather than the tree.
        let all_px: Vec<PixelOut> = new
            .tree
            .leaves()
            .flat_map(|i| new.st.pixels.get(i).cloned().unwrap_or_default())
            .collect();
        let (lo, hi) = colour::range(&all_px, Scalar::ShapeSpread);
        let sites = colour::landmarks(&grid::decode_state(&chart, 0, cx, cy).m);
        let rgb = move |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);

        // ---- 1. the frontier being spent, one frame per round ---------------------------
        let rounds = new.st.iterations;
        let mut budget_frames = Vec::new();
        let mut budget_wire = Vec::new();
        for k in 0..=rounds {
            let lv = leaves_at_round(&new.tree, k);
            let shadow = shadow_of(&new.tree, &lv);
            let masked = mask(&new.st.pixels, &lv);
            let f = adaptive::render(
                &shadow, &masked, &cam, res, adaptive::TexelMode::Adaptive, |p| rgb(p),
            )
            .0;
            let mut wf = f.clone();
            wire::draw(&mut wf, res, res, &wire::boxes_from_leaves(&new.tree, &cam, res, &lv), 1);
            budget_frames.push(f);
            budget_wire.push(wf);
        }
        let _ = apng::write(&format!("{dir}/{name}_budget.png"), res, res, &budget_frames, 1, 2);
        let _ = apng::write(&format!("{dir}/{name}_budget_wire.png"), res, res, &budget_wire, 1, 2);

        // ---- 2. the shipped criterion against the measured-best one ----------------------
        let old = descend(&chart, cx, cy, half, res, budget, tau, alpha_hi, OLD, K_MAIN, &ens);
        let rounds2 = rounds.max(old.st.iterations);
        let mut pair_frames = Vec::new();
        for k in 0..=rounds2 {
            let la = leaves_at_round(&old.tree, k.min(old.st.iterations));
            let lb = leaves_at_round(&new.tree, k.min(rounds));
            let fa = render_leaves(&old.tree, &old.st.pixels, &cam, res, &la, &rgb);
            let fb = render_leaves(&new.tree, &new.st.pixels, &cam, res, &lb, &rgb);
            pair_frames.push(side_by_side(&fa, &fb, res));
        }
        let _ = apng::write(
            &format!("{dir}/{name}_oldnew.png"), res * 2 + 1, res, &pair_frames, 1, 2,
        );

        // ---- 3. the demotion mechanism: k_frac 0.25 -> 1.0 -------------------------------
        //
        // A separate descent per step, because k_frac changes what gets computed and cannot be
        // reconstructed from one run the way the round axis can.
        let mut k_frames = Vec::new();
        let mut k_leaves = Vec::new();
        for &k in &K_STEPS {
            let r = if (k - K_MAIN).abs() < 1e-12 {
                None
            } else {
                Some(descend(&chart, cx, cy, half, res, budget, tau, alpha_hi, NEW, k, &ens))
            };
            let run = r.as_ref().unwrap_or(&new);
            let lv: Vec<usize> = run.tree.leaves().collect();
            k_leaves.push(lv.len());
            k_frames.push(render_leaves(&run.tree, &run.st.pixels, &cam, res, &lv, &rgb));
        }
        let _ = apng::write(&format!("{dir}/{name}_kfrac.png"), res, res, &k_frames, 1, 1);

        // ---- the row ---------------------------------------------------------------------
        let lv: Vec<usize> = new.tree.leaves().collect();
        let depth = lv.iter().map(|&i| new.tree.nodes[i].level).max().unwrap_or(0);
        let veto = lv
            .iter()
            .filter(|&&i| {
                matches!(
                    new.tree.nodes[i].decision,
                    Decision::ScreenFloor | Decision::MaxRelDepth
                )
            })
            .count() as f64
            / lv.len().max(1) as f64;
        println!("{name:>20} {rounds:>7} {:>7} {depth:>6} {:>6.0}% {:>8} {:>9} {:>9.1}",
                 lv.len(), 100.0 * veto, old.tree.leaves().count(), lv.len(),
                 t0.elapsed().as_secs_f64());
        println!("{:>20}   k_frac {K_STEPS:?} -> leaves {k_leaves:?}", "");
        let _ = Agg::Median;
    }

    println!();
    println!("Three animations per chart, in {dir}/:");
    println!("  _budget.png / _budget_wire.png  one frame per descent round -- the frontier being");
    println!("                                  spent. What the depth ladder cannot show.");
    println!("  _oldnew.png                     two panels, {}/median | {}/median.",
             OLD.name(), NEW.name());
    println!("  _kfrac.png                      k_frac {K_STEPS:?}; the last frame is the control.");
    println!();
    println!("READ THE veto% COLUMN FIRST. Where it is high the two panels of _oldnew are largely");
    println!("the same cap reached by two routes, and the difference between them is not a");
    println!("criterion difference. The charts worth looking at are the ones with a low veto share.");
}
