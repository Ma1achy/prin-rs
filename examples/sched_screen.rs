//! **§4.2 — q1, q2 and q3 re-run with the screen floor.**
//!
//! PR #11 measured these without a camera, so without the veto the contract calls *"the everyday
//! refinement stop"*. At `N = 8` on a 512² viewport samples stop being displayable at **level 6**,
//! and that descent reached **level 12** — 4096x past it. So its answers describe the criterion
//! minus its principal stop condition.
//!
//! Both configurations run here, in one process, so the comparison is not across two write-ups.
//!
//! **A structural cap, stated before the numbers are read:** under the veto a tree cannot exceed
//! `4^6 = 4096` leaves. near-field's 4617 cannot recur, and its absence is arithmetic, not
//! improvement.
//!
//! Run: `cargo run --release --example sched_screen [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let budget: usize = arg(1, 50_000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, E+1=8, t=13, f64");
    println!("viewport 512x512, camera framing the whole root box -> camera_depth 0.");
    println!();
    println!("The screen floor caps the tree at 4^6 = 4096 leaves STRUCTURALLY. A smaller leaf");
    println!("count than PR #11's is arithmetic, not improvement, and is not read as one.");
    println!();
    println!("{:>14} {:>9} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>7} {:>7} {:>8}",
             "region", "veto", "quads", "leaves", "depth", "floor", "keep", "screen", "relcap",
             "budget", "wall_s");

    for region in ["far", "near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        for veto in [false, true] {
            let cam = veto.then(|| Camera::framing(root.cx, root.cy, 0.05, 512));
            let cfg = SchedCfg {
                budget,
                tau_display: tau,
                alpha_hi,
                alpha_lo: alpha_hi * 0.4,
                max_level: None,
                camera: cam,
                ..Default::default()
            };
            let (t, st) =
                scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
            let depth = t.depth_histogram().len().saturating_sub(1);
            println!("{:>14} {:>9} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>7} {:>7} {:>8.1}",
                     region, if veto { "on" } else { "OFF" },
                     st.quads_computed, leaves.len(), depth,
                     c(D::Floor), c(D::Keep), c(D::ScreenFloor), c(D::MaxRelDepth),
                     c(D::BudgetExhausted), st.wall_seconds);

            if veto {
                // q2: does the floor engage, and which floor is doing the work?
                let n = leaves.len().max(1) as f64;
                println!("{:>14}          q2: floor {:.1}%  keep {:.1}%  screen {:.1}%  relcap {:.1}%  precision {:.1}%",
                         "", 100.0 * c(D::Floor) as f64 / n, 100.0 * c(D::Keep) as f64 / n,
                         100.0 * c(D::ScreenFloor) as f64 / n,
                         100.0 * c(D::MaxRelDepth) as f64 / n,
                         100.0 * c(D::PrecisionFloor) as f64 / n);
                let hist = t.depth_histogram();
                print!("{:>14}          depth histogram:", "");
                for (l, k) in hist.iter().enumerate() {
                    if *k > 0 {
                        print!(" {l}:{k}");
                    }
                }
                println!();
                // Which stop actually terminated the descent: the criterion, or the view?
                let by_view = c(D::ScreenFloor) + c(D::MaxRelDepth);
                let by_crit = c(D::Floor) + c(D::Keep);
                println!("{:>14}          stopped by the VIEW: {by_view}, by the CRITERION: {by_crit}", "");
            }
        }
        println!();
    }

    println!("q1 under the veto is a different question from q1 without it: the descent cannot");
    println!("run away, so 'does it terminate' becomes 'what terminated it'. Read the last line");
    println!("of each region: if the view stopped most leaves, the criterion was never tested.");
}
