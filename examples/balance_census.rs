//! **What the 2:1 balance constraint costs under `Policy::Tolerance`, on the six charts.**
//!
//! `SchedCfg::balance` defaults to `false` and `descend_live_with` never calls the pass at all, so
//! every tree in `results/` is unbalanced. Whether to turn it on looks like a corpus-invalidation
//! decision and may be a non-decision: `tests/slippy.rs` had to drop to `n = 4` on `deep interior`
//! to produce a 2:1 violation at all, because at the production `N = 8` the trees it swept were
//! already gap 1 in **all twenty-four** cells of `alpha_hi x tau x n`. But that sweep is
//! `Policy::Alpha`-era; under the tolerance policy the trees are deeper and more varied, and the
//! balance situation is unmeasured. This measures it.
//!
//! **The number that decides it is `forced/split`** -- the share of splits the criterion did not
//! ask for. §4.4: if it is large, the budget is going on geometry rather than physics, and that
//! has to be visible rather than inferred. If it is zero the default question is moot.
//!
//! **`gap` is the control.** A balanced tree trivially satisfies 2:1 if the unbalanced one already
//! did, so the unbalanced column is printed beside it: where `gap off` is already 1, this chart
//! says nothing about the constraint and the row is marked inert. That is the arm that turned a
//! silently-vacuous test into a failing one twice in `tests/slippy.rs`.
//!
//! **And `moved` is not `forced`.** A forced split changes the tree below it, so the decisions of
//! quads that were never forced can move too. The diff is over the decision of every quad present
//! in both trees, keyed by box rather than by index, because a split renumbers nothing but adds.
//!
//! What an unbalanced tree costs *here* is a 16x jump in apparent resolution across one edge, not
//! a hole: texels are nearest-neighbour clipped to the quad box and quadtree leaves tile the root
//! exactly at any depth difference (`adaptive::coverage` reads zero gaps, zero overlaps). It
//! becomes a literal crack under interpolation across leaves or when quads are drawn as GPU
//! geometry. Cosmetic today, load-bearing later -- so the cost is what is being measured, not the
//! rule.
//!
//! Run: `cargo run --release --example balance_census -- [root=results] [charts] [res=512]`

use std::collections::HashMap;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::quad::{Decision, Dir, QuadTree};
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
    // Every named slice with its own window, from the one table in `grid`.
    if let Some((chart, cx, cy, half)) = grid::named_slice(name) {
        return Some(Target { name: name.into(), chart, cx, cy, half, body: 0 });
    }
    grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == name)
        .map(|(n, chart, cx, cy, half)| Target { name: n.into(), chart, cx, cy, half, body: 0 })
}

/// The worst adjacent-leaf level difference over the whole tree. `1` satisfies 2:1.
fn worst_gap(t: &QuadTree) -> u32 {
    let mut w = 0u32;
    for i in t.leaves() {
        for d in Dir::ALL {
            if let Some(j) = t.neighbour(i, d) {
                if t.nodes[j].children.is_none() {
                    w = w.max(t.nodes[i].level.abs_diff(t.nodes[j].level));
                }
            }
        }
    }
    w
}

/// **Adjacent leaf pairs differing by two or more levels** — what the pass actually fires on.
///
/// Depth *variance* was the first candidate and it points the wrong way on the first two charts:
/// `near-field` has the higher variance (2.098) and the lower tax (1.51x), `deep interior` the
/// lower variance (1.516) and the higher tax (1.76x). Variance is the spread of depths; the pass
/// fires on the **length of the depth discontinuity**, which is a perimeter and not a spread. A
/// tree can hold high variance across few boundaries (one deep blob in a coarse field) or low
/// variance across many (interleaved). This counts the boundaries.
///
/// It should predict `forced` closely but not exactly, because one split can repair several
/// adjacencies at once and splitting can create new ones — so it is an upper bound on the repairs
/// needed at the moment it is measured, not a forecast of the fixed point.
fn violating_adjacencies(t: &QuadTree) -> usize {
    let mut n = 0usize;
    for i in t.leaves() {
        for d in Dir::ALL {
            if let Some(j) = t.neighbour(i, d) {
                if t.nodes[j].children.is_none() && t.nodes[i].level.abs_diff(t.nodes[j].level) >= 2 {
                    n += 1;
                }
            }
        }
    }
    n / 2 // each pair is seen from both sides
}

/// Variance of leaf depth. **The direct measure of the thing 2:1 removes**, so it is what the
/// geometry tax should scale with -- and reasoning from tree size instead would be inferring the
/// mechanism from a proxy. A selective tree (coarse everywhere, deep on a filament) has high
/// variance; a tree sitting mostly at the screen floor has low variance whatever its size.
fn depth_variance(t: &QuadTree) -> f64 {
    let d: Vec<f64> = t.leaves().map(|i| t.nodes[i].level as f64).collect();
    if d.is_empty() {
        return f64::NAN;
    }
    let m = d.iter().sum::<f64>() / d.len() as f64;
    d.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / d.len() as f64
}

/// Decisions keyed by box, not by index: a split adds nodes, so indices are not comparable
/// between two trees that diverged.
fn by_box(t: &QuadTree) -> HashMap<(u32, i64, i64), Decision> {
    t.nodes
        .iter()
        .map(|q| {
            (
                (q.level, (q.cx / 1e-12).round() as i64, (q.cy / 1e-12).round() as i64),
                q.decision,
            )
        })
        .collect()
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let charts: Vec<String> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| {
            "near-field,deep interior,preset_prho,preset_shape,config_stability,preset_shape_h1"
                .into()
        })
        .split(',')
        .map(|s| s.trim().to_string())
        // **`-` runs the frame arm alone.** The static block is six charts at 512^2 and is the
        // expensive half; the frame arm sweeps its own throttle, so a third pass that re-ran the
        // static table would measure the same thing a third time and cost an hour doing it.
        .filter(|s| !s.is_empty() && s != "-")
        .collect();
    let res: usize = arg(3, 512);
    // **The throttle is a control arm, not a constant.** Under `k_frac < 1` the criterion splits
    // only the top slice of its want-list each round while `balance_pass` repairs *every*
    // violation, so a share measured at the production 0.25 may be a property of the throttle
    // rather than of balance. Argument four, so both arms run without an edit -- the standing
    // *"an argument hardcoded past is worse than an argument missing"*.
    let k_frac: f64 = arg(4, SchedCfg::default().k_frac);
    // **Argument five: frames under a BINDING FRAME QUOTA.** The throttle-invariance result below
    // was measured at a non-binding total budget (20000 against a largest tree of 4869), and its
    // own write-up says so: `k_frac` truncates per round, deferred quads are re-decided next round,
    // so everything the criterion wants eventually happens and only the ORDER changes. Under a
    // frame quota -- which is the whole point of the slippy map -- it binds by construction, and
    // the record marks that cell **unmeasured**. `0` skips the arm.
    let frames: usize = arg(5, 0);
    let frame_quads: usize = arg(6, 16);
    let frame_charts: Vec<String> = std::env::args()
        .nth(7)
        .unwrap_or_else(|| "near-field,deep interior".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    println!("balance_census: {res}^2 camera, policy {:?}, tau {:e}, k_frac {}, t={} n_sync={}",
             SchedCfg::default().policy, SchedCfg::default().tau_display,
             k_frac, ens.t_max, ens.n_sync);
    println!("  config: {}", ens.provenance());
    println!("  `forced/split` is the number that decides the default; `gap off` is the control --");
    println!("  where it is already 1 the chart says nothing about the constraint.\n");
    if charts.is_empty() {
        println!("  (static block skipped: chart list empty)");
    } else {
        println!("{:>18} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>9} {:>10} {:>8} {:>7} {:>7} {:>6} {:>6}",
                 "chart", "gap off", "gap on", "quads-", "quads+", "quad x", "forced", "splits",
                 "fr/split", "substeps+", "steps x", "crit-", "dvar-", "viol-", "moved");
    }

    for name in &charts {
        let Some(t) = target(name) else {
            println!("{name:>18}  unknown target, skipped");
            continue;
        };
        let run = |balance: bool| {
            let cfg = SchedCfg {
                camera: Some(Camera::framing(t.cx, t.cy, t.half, res)),
                chart: t.chart,
                balance,
                k_frac,
                budget: 20000,
                ..Default::default()
            };
            scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens, Precision::F64)
        };
        let (t_off, s_off) = run(false);
        let (t_on, s_on) = run(true);

        // Splits the criterion asked for, in the balanced tree: every internal node that is not
        // itself labelled BalanceForced.
        let splits = t_on
            .nodes
            .iter()
            .filter(|q| q.children.is_some() && q.decision != Decision::BalanceForced)
            .count();
        // **`SchedStats::balance_forced` counts NODES CREATED, not quads split** -- `balance_pass`
        // returns `made`, which is four children per split. Comparing it against a count of split
        // parents is a factor-of-four error in the direction that flatters the finding, and the
        // arithmetic that catches it is `1 + 4*(forced + splits) == nodes`, asserted below.
        let forced = t_on
            .nodes
            .iter()
            .filter(|q| q.decision == Decision::BalanceForced)
            .count();
        let frac = forced as f64 / (splits + forced).max(1) as f64;
        assert_eq!(
            1 + 4 * (splits + forced),
            t_on.nodes.len(),
            "every node is the root or a child of a split; if this fails the two counts are in \
             different units"
        );

        // **`forced` is not additive, and the sign goes both ways.** A forced split changes what
        // its descendants and its parent subsequently decide, so the criterion's own split count
        // moves too: measured, near-field 31 -> 33 and `deep interior` 56 -> 70 (the criterion
        // splits MORE), `preset_prho` 244 -> 200 (it splits FEWER, 44 of them). So `forced` alone
        // does not predict the tree growth and `quad x` is the honest cost column.
        let splits_off = t_off.nodes.iter().filter(|q| q.children.is_some()).count();

        // Decisions that moved, over the boxes present in BOTH trees.
        let (a, b) = (by_box(&t_off), by_box(&t_on));
        let shared: Vec<_> = a.keys().filter(|k| b.contains_key(*k)).collect();
        let moved = shared.iter().filter(|k| a[**k] != b[**k]).count();

        // **Quads are not the cost; substeps are.** A forced split integrates `N^2 * (E+1)` new
        // trajectories, and quads cost differently -- `results/payload/README.md`'s cost ledger
        // measures steps/trajectory spreading 9.2x within one level on the sea charts. Reporting
        // the geometry tax in quads alone would price a forced split in the smooth surroundings
        // the same as one at a close encounter. *Read `steps`, not `secs`* -- and not quads either.
        let steps = |t: &QuadTree| t.nodes.iter().map(|q| q.red.total_substeps as f64).sum::<f64>();
        let (sub_off, sub_on) = (steps(&t_off), steps(&t_on));

        let (g_off, g_on) = (worst_gap(&t_off), worst_gap(&t_on));
        println!("{:>18} {g_off:>7} {g_on:>7} {:>7} {:>7} {:>6.2}x {forced:>7} {splits:>7} {frac:>9.4} {:>10.3e} {:>7.2}x {:>7} {:>7.3} {:>6} {moved:>6}{}",
                 t.name, s_off.quads_computed, s_on.quads_computed,
                 s_on.quads_computed as f64 / s_off.quads_computed.max(1) as f64,
                 sub_on, sub_on / sub_off.max(1.0), splits_off, depth_variance(&t_off), violating_adjacencies(&t_off),
                 if g_off <= 1 { "   [inert: already 2:1 without the pass]" } else { "" });

        // **A budget-bound row is not a measurement.** `balance_pass` takes `room` and returns
        // early when it runs out, so a forced count read under an exhausted budget is low because
        // there was nowhere to put the splits, not because balance is cheap. Same shape as the
        // standing *"a control that is budget-bound is not a control"*.
        if s_off.budget_exhausted || s_on.budget_exhausted {
            println!("{:>18}  BUDGET EXHAUSTED (off {}, on {}) -- `forced` is a floor, not a count",
                     "", s_off.budget_exhausted, s_on.budget_exhausted);
        }

        if g_on > 1 {
            println!("{:>18}  BALANCED TREE STILL VIOLATES 2:1 (gap {g_on}) -- the pass did not \
                      reach a fixed point", "");
        }
    }

    if !charts.is_empty() {
        println!("\nThe pass runs after the criterion and can only add, so `quads+ >= quads-` always.");
        println!("`moved` counts boxes in both trees whose decision differs -- a forced split changes");
        println!("what its descendants decide, so it is not bounded by `forced`.");
    }

    // -------------------------------------------------------------------------------------
    // The cell the record marks unmeasured: `quad x` under a BINDING frame quota.
    // -------------------------------------------------------------------------------------
    if frames > 0 {
        use prin_rs::session::{FrameQuota, Session, SessionCfg};
        println!("\n== the balance tax under a BINDING FRAME QUOTA");
        println!("  {frames} frames x {frame_quads} quads, rounds 2 -- a hard ceiling of {} quads,",
                 frames * frame_quads);
        println!("  against unconstrained trees in the hundreds to thousands above. `bound` is the");
        println!("  count of frames the quota actually stopped: **a non-binding row proves nothing**");
        println!("  and is the same regime measured twice, which is how the constant 400 failed.\n");
        println!("{:>18} {:>6} {:>8} {:>8} {:>8} {:>8} {:>7} {:>7}",
                 "chart", "k_frac", "quads-", "quads+", "quad x", "forced", "bound-", "bound+");
        let mut pairs: Vec<(String, f64, f64)> = Vec::new();
        for name in &frame_charts {
            let Some(t) = target(name) else { continue };
            for &k in &[0.25f64, 1.0] {
                let run = |balance: bool| {
                    let cam = Camera::framing(t.cx, t.cy, t.half, res);
                    let cfg = SessionCfg {
                        sched: SchedCfg {
                            camera: Some(cam),
                            chart: t.chart,
                            balance,
                            k_frac: k,
                            budget: 20000,
                            ..Default::default()
                        },
                        quota: FrameQuota { quads: frame_quads, substeps: None, rounds: 2 },
                        ..Default::default()
                    };
                    let mut s = Session::new(t.cx, t.cy, t.half, t.body, cam, cfg, ens.t_max);
                    let mut bound = 0usize;
                    for _ in 0..frames {
                        let (hit, _) =
                            s.step(&|sl, i| prin_rs::ensemble::pixel::evaluate::<f64>(sl, i, &ens));
                        bound += usize::from(hit != prin_rs::session::QuotaHit::Drained);
                    }
                    let forced = s
                        .tree()
                        .nodes
                        .iter()
                        .filter(|q| q.decision == Decision::BalanceForced)
                        .count();
                    (s.stats().quads_computed, forced, bound)
                };
                let (q_off, _, b_off) = run(false);
                let (q_on, forced, b_on) = run(true);
                let x = q_on as f64 / q_off.max(1) as f64;
                println!("{name:>18} {k:>6.2} {q_off:>8} {q_on:>8} {x:>8.4} {forced:>8} {b_off:>7} {b_on:>7}");
                pairs.push((name.clone(), k, x));
            }
        }
        // **The claim under test, stated as a comparison rather than left to the eye.**
        for name in &frame_charts {
            let a = pairs.iter().find(|p| &p.0 == name && p.1 == 0.25).map(|p| p.2);
            let b = pairs.iter().find(|p| &p.0 == name && p.1 == 1.0).map(|p| p.2);
            if let (Some(a), Some(b)) = (a, b) {
                println!("  {name}: quad x is {a:.4} at k=0.25 and {b:.4} at k=1.0 -- {}",
                         if (a - b).abs() < 0.01 { "INVARIANT" } else { "THROTTLE-DEPENDENT" });
            }
        }
        println!("\n  The non-binding result says `quad x` is throttle-invariant because everything");
        println!("  the criterion wants eventually happens and only the order changes. A frame quota");
        println!("  removes that: work the throttle defers may never be reached at all.");
    }
}
