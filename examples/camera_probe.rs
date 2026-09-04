//! **Does `camera_bias` change anything, and what does the margin buy?**
//!
//! §4.3 wants priority to be a **product** of camera relevance and structure, with the destination
//! model stated. §3 offers three: velocity extrapolation (rejected — it prefetches past a
//! flick-and-stop), the swept path, and *"simply widen the viewport margin and drop prediction
//! entirely"*, which is the honest baseline the others must beat. The margin is that baseline's
//! only parameter, and this measures it.
//!
//! **First, whether the knob bites at all** — the question A2 forced. `camera_bias` reaches only
//! `priority`, which reaches only `order_queue`, which changes the tree **only where the ranking is
//! truncated**. Under `k_frac < 1` a deferred quad is re-decided next round, so at a *non-binding*
//! budget everything the criterion wants eventually happens and the ordering washes out: the
//! balance census measured exactly that, with `k_frac` 0.25 and 1.0 giving identical trees at
//! budget 20000. So this runs **both regimes** and prints `budget_exhausted`, because a knob
//! measured only where it cannot act is the failure this project keeps meeting.
//!
//! The camera is zoomed **into a corner**, not framing the root: `Camera::framing` sets
//! `half_world` to the root half-width, every quad is then fully visible, `relevance` is 1.0
//! everywhere and the bias is an identity. `rel span` is printed so a dead arm is visible rather
//! than inferred.
//!
//! Run: `cargo run --release --example camera_probe -- [root=results] [charts] [res=512]`

use std::collections::HashMap;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::quad::{Decision, QuadTree};
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

/// Decisions keyed by box, so two trees that diverged are still comparable.
fn by_box(t: &QuadTree) -> HashMap<(u32, i64, i64), Decision> {
    t.nodes
        .iter()
        .map(|q| ((q.level, (q.cx / 1e-12).round() as i64, (q.cy / 1e-12).round() as i64), q.decision))
        .collect()
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let charts: Vec<String> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "near-field,deep interior,preset_shape_h1".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let res: usize = arg(3, 512);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    println!("camera_probe: {res}^2, policy {:?}, k_frac {}, t={} n_sync={}",
             SchedCfg::default().policy, SchedCfg::default().k_frac, ens.t_max, ens.n_sync);
    println!("  config: {}", ens.provenance());
    println!("  camera zoomed to a quarter width toward a corner, so `relevance` varies;");
    println!("  `rel span` is the liveness arm -- 0.000 means the bias is an identity.\n");
    println!("{:>18} {:>8} {:>9} {:>7} {:>7} {:>7} {:>8} {:>8} {:>7}",
             "chart", "budget", "margin", "quads", "leaves", "moved", "of", "exhausted", "rel span");

    for name in &charts {
        let Some(t) = target(name) else {
            println!("{name:>18}  unknown target, skipped");
            continue;
        };
        let cam = Camera::framing(t.cx + t.half * 0.4, t.cy + t.half * 0.4, t.half * 0.25, res);

        // **The binding budget is DERIVED, not chosen.** A constant does not bind: 400 was picked
        // against a remembered 125-quad near-field tree and the zoomed camera makes it 389, so the
        // "binding" arm reproduced the non-binding one exactly — the same regime measured twice,
        // which is the failure this whole probe exists to avoid. Measure the unconstrained tree
        // first and take a fraction of it.
        let free = {
            let cfg = SchedCfg { camera: Some(cam), chart: t.chart, budget: 20000, ..Default::default() };
            scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64).1.quads_computed
        };
        let tight = (free as f64 * 0.4).max(21.0) as usize;
        println!("{:>18}  unconstrained tree is {free} quads; the binding arm uses {tight}", t.name);

        for budget in [20000usize, tight] {
            let run = |margin: Option<f64>| {
                let cfg = SchedCfg {
                    camera: Some(cam),
                    camera_bias: margin,
                    chart: t.chart,
                    budget,
                    ..Default::default()
                };
                scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64)
            };
            let (base, sbase) = run(None);
            let a = by_box(&base);

            for margin in [0.0f64, 0.5, 2.0] {
                let (tr, st) = run(Some(margin));
                let b = by_box(&tr);
                let shared: Vec<_> = a.keys().filter(|k| b.contains_key(*k)).collect();
                let moved = shared.iter().filter(|k| a[**k] != b[**k]).count();

                // The liveness arm: relevance must vary over the quads actually ranked.
                let rels: Vec<f64> =
                    tr.leaves().map(|i| {
                        let q = &tr.nodes[i];
                        cam.relevance(q.cx, q.cy, q.half, margin)
                    }).collect();
                let span = rels.iter().cloned().fold(f64::MIN, f64::max)
                    - rels.iter().cloned().fold(f64::MAX, f64::min);

                println!("{:>18} {budget:>8} {margin:>9.1} {:>7} {:>7} {moved:>7} {:>8} {:>8} {span:>7.3}{}",
                         if margin == 0.0 { t.name.as_str() } else { "" },
                         st.quads_computed, tr.leaves().count(), shared.len(),
                         st.budget_exhausted,
                         if span <= 0.0 { "  <- DEAD: relevance constant" } else { "" });
            }
            if !sbase.budget_exhausted && budget == tight {
                println!("{:>18}  budget {budget} did NOT bind -- the truncation regime was not reached, \
                          so this arm is the non-binding one again", "");
            }
        }
        println!();
    }
    println!("`moved` counts boxes present in both trees whose decision differs. Zero at a");
    println!("non-binding budget is expected: the bias reaches only `order_queue`, and a deferred");
    println!("quad is re-decided next round, so the ordering washes out when nothing truncates.");
}
