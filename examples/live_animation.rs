//! **The live descent, as a picture: the tree growing while the simulation runs.**
//!
//! Every animation before this advanced a budget, a depth cap, a `k_frac` or a chart index; none
//! advanced *time*. This one does. One `descend_live` per chart records the leaf set at every
//! sync boundary and every post-horizon round; each frame is that leaf set rendered on the
//! footprints **as they stood at that boundary** (`scheduler::project_at`: the shape and spread
//! the copies had reached by then, a copy's terminal class only once it had terminated), with
//! the wire over the coarse-ancestor fill. The playhead is the frame axis. Post-horizon rounds
//! are frames at `t = t_max`, held on the last boundary's footprints.
//!
//! It is the closest CPU picture of the target design. With merging on (the default) the tree
//! grows and merges back: a leaf that appears in one frame is in every later frame, has been
//! split, or has been merged into its parent, which is then the leaf again.
//!
//! **The debug flag is off.** `colour::DEBUG_NAN` marks a footprint with no value so it cannot be
//! read as a dark one; that is a diagnostic device and it renders as a magenta speckle. These
//! panels use `colour::Veto::Quiet` instead -- the nominal copy's own hue at the floor of the
//! lightness ramp -- so an undetermined footprint is **indistinguishable from a resolved dark
//! one** here. The count it hides is printed per chart and in every sidecar; the diagnostic
//! renders keep the flag.
//!
//! **Diagnostic, at a small viewport.** The stills are 1024²; this runs at 256 by default so a
//! chart is a minute, and the sidecar says so. Do not read a leaf count off a frame.
//!
//! Run: `cargo run --release --example live_animation -- <root> [charts=near-field,...] [res=256]
//!       [budget=4000] [eps=0.01] [k_frac=0.25] [live_stride=2]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, wire};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

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
    if name == "config_stability" {
        let (chart, cx, cy, half) = Chart::config_stability();
        return Some(Target { name: name.into(), chart, cx, cy, half, body: 0 });
    }
    grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == name)
        .map(|(n, chart, cx, cy, half)| Target { name: n.into(), chart, cx, cy, half, body: 0 })
}

fn main() {
    let Some(root) = std::env::args().nth(1) else {
        eprintln!("usage: live_animation <root> [charts] [res] [budget] [eps] [k_frac] [live_stride]");
        std::process::exit(2);
    };
    let charts: Vec<String> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "near-field,deep interior,preset_shape_h1".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let res: usize = arg(3, 256);
    let budget: usize = arg(4, 4000);
    let eps: f64 = arg(5, SchedCfg::default().tau_display);
    let k_frac: f64 = arg(6, prin_rs::scheduler::K_FRAC_RANKED);
    let live_stride: usize = arg(7, 2);
    let dir = format!("{root}/live");
    let _ = std::fs::create_dir_all(&dir);

    let ens = EnsembleCfg { refine_flagged: false, keep_live_series: true, live_stride, ..Default::default() };
    scheduler::assert_production_kernel(&ens, &dir);
    println!("live animation: {res}^2, budget {budget}, eps {eps:e}, k_frac {k_frac}, live_stride {live_stride}, t={} n_sync={}",
             ens.t_max, ens.n_sync);
    println!("  config: {}", ens.provenance());
    println!("{:>18} {:>7} {:>7} {:>6} {:>7} {:>8} {:>7}  stop", "chart", "quads", "leaves", "depth", "frames", "catchup%", "wall s");

    for name in &charts {
        let Some(t) = target(name) else {
            println!("{name:>18}  unknown target, skipped");
            continue;
        };
        let t0 = std::time::Instant::now();
        let cam = Camera::framing(t.cx, t.cy, t.half, res);
        let cfg = SchedCfg {
            budget,
            tau_display: eps,
            k_frac,
            camera: Some(cam),
            keep_pixels: true,
            chart: t.chart,
            ..Default::default()
        };
        let (tree, st) = scheduler::descend_live(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64);
        let leaves: Vec<usize> = tree.leaves().collect();
        let depth = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
        let steps: u64 = tree.nodes.iter().map(|q| q.red.total_substeps as u64).sum();

        // One window for the whole animation, from the finished tree's terminal footprints, so
        // the ramp does not move between frames.
        let all_px: Vec<PixelOut> = leaves.iter().flat_map(|&i| st.pixels.get(i).cloned().unwrap_or_default()).collect();
        let (lo, hi) = colour::range(&all_px, Scalar::ShapeSpread);
        let sites = colour::landmarks(&grid::decode_state(&t.chart, t.body, t.cx, t.cy).m);
        // **The debug flag is not a presentation colour.** `colour::DEBUG_NAN` exists so an
        // undetermined footprint cannot be mistaken for a dark one, and every diagnostic render
        // wants it; in an animation it is a magenta speckle the eye reads before the field. These
        // panels take `Veto::Quiet` -- the nominal copy's own hue at the floor of the ramp -- and
        // the count it hides is printed here and in every sidecar. Measured: the vetoed footprints
        // are copies that exhausted `max_steps`, and their nominal `shape_vec` is finite, so the
        // hue is real.
        let rgb = |p: &PixelOut| colour::rgb_veto(p, Scalar::ShapeSpread, &sites, lo, hi, colour::Veto::Quiet);
        let vetoed = all_px.iter().filter(|p| colour::vetoed(p, Scalar::ShapeSpread, &sites, lo, hi)).count();

        let n_b = st.pixels.get(0).and_then(|p| p.first()).map(|p| p.live_t.len()).unwrap_or(0);
        let mut frames = Vec::with_capacity(st.live_leaves.len());
        let mut wframes = Vec::with_capacity(st.live_leaves.len());
        for (k, lv) in st.live_leaves.iter().enumerate() {
            // The boundary this frame shows: the recorded boundary while the playhead moved,
            // the last one during the post-horizon rounds.
            let j = k.min(n_b.saturating_sub(1));
            let projected: Vec<Vec<PixelOut>> = (0..st.pixels.len())
                .map(|i| st.pixels[i].iter().map(|p| scheduler::project_at(p, j)).collect())
                .collect();
            let (f, _) = adaptive::render_leaves(&tree, &projected, &cam, res, adaptive::TexelMode::Adaptive, rgb, lv);
            let mut wf = f.clone();
            wire::draw(&mut wf, res, res, &wire::boxes_from_leaves(&tree, &cam, res, lv), depth.max(1));
            frames.push(f);
            wframes.push(wf);
        }
        // Hold the finished frame, so the loop reads as an ending rather than a snap back.
        for _ in 0..6 {
            frames.push(frames.last().unwrap().clone());
            wframes.push(wframes.last().unwrap().clone());
        }
        let dup = apng::adjacent_duplicates(&frames);
        let stem = format!("{dir}/{}", t.name.replace(' ', "_"));
        let _ = apng::write(&format!("{stem}_live.png"), res, res, &frames, 1, 3);
        let _ = apng::write(&format!("{stem}_live_wire.png"), res, res, &wframes, 1, 3);
        let curve: Vec<String> = st.live.iter().map(|p| format!("{}:{:.2}:{}", p.j, p.t, p.computed)).collect();
        for suffix in ["_live", "_live_wire"] {
            let _ = prin_rs::output::provenance_sidecar(
                &format!("{stem}{suffix}.png"),
                &ens,
                &format!(
                    "chart={} animation=time_axis frames={} deliberate_hold=6 adjacent_duplicates={dup} \
                     veto=quiet vetoed_footprints={vetoed} of={} \
                     boundaries={n_b} live_stride={live_stride} scalar=ShapeSpread window=({lo:.4e},{hi:.4e}) \
                     res={res} viewport={res} budget={budget} tau_display={eps:e} k_frac={k_frac} \
                     policy={} quads={} leaves={} depth={depth} catchup_substeps={} total_substeps={steps} \
                     stop={} growth={}\n",
                    t.chart.name(), frames.len(), all_px.len(), cfg.policy.name(), st.quads_computed, leaves.len(),
                    st.catchup_substeps, tree.stop_breakdown(), curve.join(" ")
                ),
            );
        }
        println!("{:>18} {:>7} {:>7} {:>6} {:>7} {:>7.1}% {:>7.1}  {}   dup {dup} vetoed {vetoed}/{}",
                 t.name, st.quads_computed, leaves.len(), depth, frames.len(),
                 100.0 * st.catchup_substeps as f64 / steps.max(1) as f64,
                 t0.elapsed().as_secs_f64(), tree.stop_breakdown(), all_px.len());
    }
    println!("\n{dir}/<chart>_live.png and _live_wire.png: one frame per boundary while the playhead moves,");
    println!("then one per post-horizon round. Leaves split and merge; the footprints are as they stood.");
}
