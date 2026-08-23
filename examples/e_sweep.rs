//! **§4.1 — the ensemble-copy sweep, reported in total trajectories.**
//!
//! The prediction, on record before the run: *low `E` biases toward refine, exactly as low `N`
//! did.* Same mechanism — a noisy spread estimate inflates apparent disagreement, and the
//! conservative failure direction turns that into extra splits. PR #11 measured it for `N`:
//! **N=4 spent 4x the quads of N=16.** If it holds for `E`, the cheap tier spends more quads
//! than the expensive one and the tier design partly cancels itself.
//!
//! **The seam that would hide it, and why both columns exist.** Under the screen floor a tree
//! saturates at `4^6 = 4096` leaves. If low `E` over-refines *into a cap*, the leaf counts come
//! out equal and the sweep reports a null that the veto manufactured. So every row runs twice.
//! A null with the veto on and an effect with it off is the veto, not `E`.
//!
//! The quantity that decides it is **total trajectories = leaves x N^2 x (E+1)**, not leaf count:
//! a cheap tier that needs 4x the leaves is not cheap.
//!
//! Run: `cargo run --release --example e_sweep [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn median(v: &mut Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let budget: usize = arg(1, 6000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    const N: usize = 8;
    let copies = [2usize, 4, 8, 16, 32];

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N={N}, t=13, f64");
    println!("viewport 512x512 when the veto is on; camera framing the root box.\n");
    println!("PREDICTION ON RECORD: low E biases toward refine, as low N did (N=4 spent 4x the");
    println!("quads of N=16). 'trajectories' is the quantity that decides the tier question.\n");
    println!("**A budget-limited row is not a measurement.** If `capped` is yes, the leaf count is");
    println!("a floor set by the budget, not a number the criterion chose, and it cannot be");
    println!("compared against an uncapped row in the same column.\n");
    println!("{:>14} {:>6} {:>5} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>14} {:>11} {:>9}",
             "region", "veto", "E+1", "quads", "leaves", "depth", "floor", "screen", "capped",
             "trajectories", "sib range", "wall_s");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        for veto in [true, false] {
            for &e1 in &copies {
                let ens = EnsembleCfg { n_extra: e1 - 1, refine_flagged: false, ..Default::default() };
                let cam = veto.then(|| Camera::framing(root.cx, root.cy, 0.05, 512));
                let cfg = SchedCfg {
                    n: N, budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4,
                    max_level: None, camera: cam, ..Default::default()
                };
                let (t, st) = scheduler::descend(
                    root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
                let leaves: Vec<usize> = t.leaves().collect();
                let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
                let depth = t.depth_histogram().len().saturating_sub(1);
                // The estimator's own noise, so the mechanism is visible rather than inferred:
                // the range of four siblings' alphas is a direct read on how well a quad knows
                // its own exponent, and it is already computed.
                let mut sib: Vec<f64> = t.nodes.iter()
                    .filter_map(|q| q.alpha_sibling_spread)
                    .filter(|x| x.is_finite())
                    .collect();
                let traj = leaves.len() * N * N * e1;
                println!("{:>14} {:>6} {:>5} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>14} {:>11.4} {:>9.1}",
                         region, if veto { "on" } else { "OFF" }, e1,
                         st.quads_computed, leaves.len(), depth, c(D::Floor), c(D::ScreenFloor),
                         if st.budget_exhausted { "YES" } else { "no" },
                         traj, median(&mut sib), st.wall_seconds);
            }
            println!();
        }
    }

    println!("Read `capped` first. Then read 'sib range' as the mechanism: if it falls with E while the leaf");
    println!("count also falls, the extra splits at low E were the estimator's noise, not");
    println!("structure. If it falls and the leaf count does not, the prediction is wrong.");
}
