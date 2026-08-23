//! **§3.4 in situ — the decode paths inside a real descent, and the silent failure mode.**
//!
//! The arithmetic ladder (`decode_ladder`) says where each path stops resolving its samples.
//! This asks what the *scheduler* does when it happens, and the answer is the seam worth having:
//!
//! **A collapsed decode reports a spread of ~5.6e-17, which the criterion reads as "perfectly
//! resolved".** Every footprint in the quad is the same initial condition, every copy is the same
//! trajectory, the ensemble agrees completely — and the tree stops, confident, having integrated
//! nothing distinguishable at all. That is the project's own standing pattern arriving from a new
//! direction: *a statistic can report maximum confidence precisely when it is least informed.*
//!
//! **And "zero spread" is not zero.** The first version of this measurement tested
//! `spread_median == 0.0` and reported no collapse anywhere, including where 1 of 64 initial
//! conditions was distinct. Eight identical `shape_vec`s summed and divided by eight do not
//! return the value bitwise, so the residual is ~5.55e-17 — twelve orders below `tau = 1e-4`, and
//! therefore indistinguishable from a genuinely resolved quad by any threshold anyone would set.
//! A collapse detector written as a spread comparison **cannot fire**.
//!
//! So collapse is detected **exactly**, by counting distinct initial conditions per quad
//! ([`decode::distinct`], a bitwise comparison of all 12 state components), and the spread is
//! reported beside it to show what the criterion would have believed.
//!
//! Run: `cargo run --release --example deep_zoom [budget]`

use prin_rs::camera::Camera;
use prin_rs::decode::{self, Path};
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::physics::Cart;
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
    println!("{:>6} {:>15} {:>9} {:>10} {:>13} {:>7} {:>7} {:>6} {:>17}",
             "depth", "path", "distinct", "collapsed", "spread(collap)", "quads", "leaves",
             "depth", "root decision");

    for zoom_depth in [0u32, 14, 18, 20, 22, 30, 40, 46] {
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

            // Collapse, counted EXACTLY: how many computed quads have fewer than N^2 distinct
            // initial conditions. Never a spread comparison — see the module note.
            let mut collapsed = 0usize;
            let mut computed = 0usize;
            let mut csp: Vec<f64> = Vec::new();
            for node in t.nodes.iter().filter(|q| q.red.n_footprints > 0) {
                computed += 1;
                let sl = node.slice(cfg.n, t.body, t.chart);
                let l = decode::linearise(&sl.chart, sl.body, sl.cx, sl.cy, sl.half);
                let ics: Vec<Cart<f64>> = (0..sl.npix())
                    .map(|k| {
                        let (u, v) = sl.decode_pos(k);
                        decode::sample(path, &sl.chart, sl.body, sl.cx, sl.cy, sl.half,
                                       (u - sl.cx) / sl.half, (v - sl.cy) / sl.half, &l)
                    })
                    .collect();
                if decode::distinct(&ics) < sl.npix() {
                    collapsed += 1;
                    if node.red.spread_median.is_finite() {
                        csp.push(node.red.spread_median);
                    }
                }
            }
            csp.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let csp_med = csp.get(csp.len() / 2).cloned();
            println!("{:>6} {:>15} {:>7}/64 {:>6}/{:<3} {:>13} {:>7} {:>7} {:>6} {:>17}",
                     zoom_depth, path.name(), n_distinct, collapsed, computed,
                     csp_med.map(|x| format!("{x:.3e}")).unwrap_or_else(|| "-".into()),
                     st.quads_computed, leaves.len(),
                     t.depth_histogram().len().saturating_sub(1),
                     t.nodes[0].decision.name());
        }
        println!();
    }

    println!("A row with 1/64 distinct and every quad collapsed is the failure: the tree is not");
    println!("small because the region is tame, it is small because there is nothing in it. Read");
    println!("`spread(collap)` beside it — that is the number the criterion saw, and at ~1e-17 it");
    println!("is twelve orders below any tau anyone would set. Nothing downstream can tell that");
    println!("apart from a perfectly resolved quad.");
    println!();
    println!("`root decision` explains the one-quad rows: at depth 40 the root quad's own cell");
    println!("width is already below the PRECISION floor (level ~36 at half0 = 0.05), so the");
    println!("descent stops for a numerical reason before any decode path is tested. That is a");
    println!("different limit from the collapse, and the column keeps them apart.");
}
