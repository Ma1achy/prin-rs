//! **§3.4 in situ — the decode paths inside a real descent, and the silent failure mode.**
//!
//! The arithmetic ladder (`decode_ladder`) says where each path stops resolving its samples.
//! This asks what the *scheduler* does when it happens, and the answer is the seam worth having:
//!
//! **A collapsed decode reports a spread of exactly zero, which the criterion reads as
//! "perfectly resolved".** Every footprint in the quad is the same initial condition, every copy
//! is the same trajectory, the ensemble agrees completely — and the tree stops, confident, with
//! no data at all. That is the project's own standing pattern: *a statistic can report maximum
//! confidence precisely when it is least informed.*
//!
//! So the number to read here is `zero-spread quads`, not tree size. A small tree under a
//! collapsed path is not the criterion working.
//!
//! Run: `cargo run --release --example deep_zoom [budget]`

use prin_rs::camera::Camera;
use prin_rs::decode::{self, Path};
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::physics::Cart;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

const PATHS: [Path; 4] = [Path::DirectF64, Path::DirectF32, Path::LinNaiveF32, Path::LinSplitF32];

fn main() {
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(400);
    let root = grid::region("near-field", 2, 2, 0.05).unwrap();

    println!("near-field, N=8, E+1=8, t=13, tau=1e-4, alpha_hi=0.2, budget {budget} quads.");
    println!("The camera is zoomed to the stated depth about the region centre, so the root quad");
    println!("is a level-`depth` box: half = 0.05 / 2^depth.\n");
    println!("'distinct' is how many of the 64 footprint ICs in the root quad are actually");
    println!("different. Read it BEFORE the tree columns: a path that has collapsed builds a");
    println!("small confident tree out of nothing.\n");
    println!("{:>6} {:>15} {:>9} {:>13} {:>7} {:>7} {:>6} {:>6} {:>9}",
             "depth", "path", "distinct", "zero-spread", "quads", "leaves", "depth", "keep", "wall_s");

    for zoom_depth in [0u32, 20, 40, 46] {
        let half = 0.05 / (2f64).powi(zoom_depth as i32);
        for path in PATHS {
            // Distinctness of the root quad's own footprints, measured directly.
            let slice = grid::Slice::body_plane(8, 8, root.cx, root.cy, half, root.body);
            let lin = decode::linearise(&slice.chart, slice.body, slice.cx, slice.cy, slice.half);
            let ics: Vec<Cart<f64>> = (0..64)
                .map(|k| {
                    let (u, v) = slice.decode_pos(k);
                    decode::sample(path, &slice.chart, slice.body, slice.cx, slice.cy, slice.half,
                                   (u - slice.cx) / slice.half, (v - slice.cy) / slice.half, &lin)
                })
                .collect();
            let n_distinct = decode::distinct(&ics);

            let ens = EnsembleCfg { refine_flagged: false, decode_path: path, ..Default::default() };
            let cfg = SchedCfg {
                budget, tau_display: 1e-4, alpha_hi: 0.2, alpha_lo: 0.08,
                camera: Some(Camera::at_depth(root.cx, root.cy, 0.05, 512, zoom_depth)),
                ..Default::default()
            };
            let (t, st) = scheduler::descend(
                root.cx, root.cy, half, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let zero = t.nodes.iter()
                .filter(|q| q.red.n_footprints > 0 && q.red.spread_median == 0.0)
                .count();
            println!("{:>6} {:>15} {:>7}/64 {:>10}/{:<3} {:>7} {:>7} {:>6} {:>6} {:>9.1}",
                     zoom_depth, path.name(), n_distinct, zero, st.quads_computed,
                     st.quads_computed, leaves.len(),
                     t.depth_histogram().len().saturating_sub(1),
                     leaves.iter().filter(|&&i| t.nodes[i].decision == D::Keep).count(),
                     st.wall_seconds);
        }
        println!();
    }

    println!("A row with 1/64 distinct and a high zero-spread count is the failure: the tree is");
    println!("not small because the region is tame, it is small because there is nothing in it.");
    println!("Compare against decode_ladder.txt, which gives the same collapse depths without");
    println!("integrating anything.");
}
