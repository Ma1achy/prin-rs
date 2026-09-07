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
//! **The stop-reason breakdown, added 2026-09-06 — and it is the column this harness lacked.**
//! It printed `root decision` alone, which names what stopped the *root* and says nothing about
//! the other leaves. The standing rule is *never quote a leaf count without the stop-reason
//! breakdown*, and this file's own concluding paragraph reads a mechanism off a leaf count. The
//! headroom is real: at depth 0 the descent computes **397 quads against a budget of 400**, three
//! spare, so a change anywhere upstream can make the tree budget-bound with nothing saying so.
//! `budget?` is printed per row for exactly that.
//!
//! **And the live arm is thin — 28 of the 32 rows are inert.** Only the depth-0 block has a real
//! tree; every block from depth 14 down reads 21 quads / 16 leaves / depth 2 for all four paths,
//! camera-vetoed. Quoting a 32-row diff as "reproduces" would be the standing `nf w/ hot nbr`
//! failure, 1.0000 in every cell. The `distinct` column is the arm with teeth and it is an
//! **identity**: it is computed from `Slice::body_plane -> decode::linearise -> decode::sample ->
//! decode::distinct`, a bitwise comparison, before `EnsembleCfg` is constructed at all. **No
//! integration enters it**, so no integrator, budget or step-control change can move it. If it
//! moves, the cause is the decoder or `Chart::default_half()`, and the response is to stop and
//! look there rather than to re-baseline.
//!
//! Run: `cargo run --release --example deep_zoom [budget] [root]`
//!
//! It writes `<root>/output/deep_zoom.txt` itself. The committed file was redirected stdout, so
//! the reproduction command lived only in shell history — *a documented reproduction command can
//! be wrong, and only running it finds out*, and one that is not written down cannot even be run.

use prin_rs::camera::Camera;
use prin_rs::decode::{self, Path};
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::physics::Cart;
use prin_rs::output::Log;
use prin_rs::quad::Decision;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::logln;

const PATHS: [Path; 4] = [Path::DirectF64, Path::DirectF32, Path::LinNaiveF32, Path::LinSplitF32];

fn main() {
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(400);
    let out: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let root = grid::region("near-field", 2, 2, 0.05).unwrap();
    let _ = std::fs::create_dir_all(format!("{out}/output"));
    let log = Log::tee(&format!("{out}/output/deep_zoom.txt"));
    let log = &log;

    // The kernel this ran on, absolutely and not as a diff. The decode ladder is an identity and
    // cannot move under it; the tree columns can, and this is what dates them.
    logln!(log, "config: {}",
           EnsembleCfg { refine_flagged: false, ..Default::default() }.provenance());
    logln!(log, "reproduce: cargo run --release --example deep_zoom {budget} {out}\n");
    logln!(log, "near-field, N=8, E+1=8, t=13, tau=1e-4, alpha_hi=0.2, budget {budget} quads.");
    logln!(log, "The camera is zoomed to the stated depth about the region centre, so the root quad");
    logln!(log, "is a level-`depth` box: half = 0.05 / 2^depth.\n");
    logln!(log, "'distinct' is how many of the 64 footprint ICs in the root quad are actually");
    logln!(log, "different. Read it BEFORE the tree columns: a path that has collapsed builds a");
    logln!(log, "small confident tree out of nothing.\n");
    logln!(log, "{:>6} {:>15} {:>9} {:>10} {:>13} {:>7} {:>7} {:>6} {:>7} {:>17} {}",
           "depth", "path", "distinct", "collapsed", "spread(collap)", "quads", "leaves",
           "depth", "budget?", "root decision", "leaf stop reasons");

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
            // **Every leaf's stop reason, from `Decision::ALL` rather than a hand-written
            // table.** The two `.prnq` decoders went stale by keeping their own list and one of
            // them dropped five variants silently; `ALL` is the table, so a new variant appears
            // here without an edit.
            let mut stops: Vec<String> = Decision::ALL
                .iter()
                .filter_map(|d| {
                    let n = leaves.iter().filter(|&&i| t.nodes[i].decision == *d).count();
                    (n > 0).then(|| format!("{}:{n}", d.name()))
                })
                .collect();
            if stops.is_empty() {
                stops.push("-".into());
            }
            logln!(log, "{:>6} {:>15} {:>7}/64 {:>6}/{:<3} {:>13} {:>7} {:>7} {:>6} {:>7} {:>17} {}",
                   zoom_depth, path.name(), n_distinct, collapsed, computed,
                   csp_med.map(|x| format!("{x:.3e}")).unwrap_or_else(|| "-".into()),
                   st.quads_computed, leaves.len(),
                   t.depth_histogram().len().saturating_sub(1),
                   if st.budget_exhausted { "BOUND" } else { "no" },
                   t.nodes[0].decision.name(), stops.join(" "));
        }
        logln!(log, "");
    }

    logln!(log, "A row with 1/64 distinct and every quad collapsed is the failure: the tree is not");
    logln!(log, "small because the region is tame, it is small because there is nothing in it. Read");
    logln!(log, "`spread(collap)` beside it — that is the number the criterion saw, and at ~1e-17 it");
    logln!(log, "is twelve orders below any tau anyone would set. Nothing downstream can tell that");
    logln!(log, "apart from a perfectly resolved quad.");
    logln!(log, "");
    logln!(log, "READ `budget?` AND THE STOP REASONS BEFORE THE COUNTS. At a bound budget `quads`");
    logln!(log, "and `leaves` are arithmetic in the split count and not evidence: 99 splits give");
    logln!(log, "1 + 4*99 = 397 computed and 397 - 99 = 298 leaves, whatever the tree did. Depth 0");
    logln!(log, "is bound here, so only `depth` and the stop reasons carry information in it. And");
    logln!(log, "28 of the 32 rows are inert as a live arm -- every block from depth 14 down reads");
    logln!(log, "21 quads / 16 leaves / depth 2 for all four paths, camera-vetoed. The `distinct`");
    logln!(log, "column is the arm with teeth, and it is an IDENTITY: no integration enters it.");
    logln!(log, "");
    logln!(log, "`root decision` explains the one-quad rows: at depth 40 the root quad's own cell");
    logln!(log, "width is already below the PRECISION floor (level ~36 at half0 = 0.05), so the");
    logln!(log, "descent stops for a numerical reason before any decode path is tested. That is a");
    logln!(log, "different limit from the collapse, and the column keeps them apart.");
}
