//! `r_coll` has no reference. Its default is therefore reported as a measurement.
//!
//! Two questions, one sweep:
//!
//! 1. **How much does the choice move the answer?** Outcome fractions across
//!    `r_coll/R ∈ {1e-4, 1e-3, 1e-2}` on a 64x64 near-field grid.
//!
//! 2. **Does the LC branch move them?** An outcome is a *classification boundary*, so a
//!    registration error near the branch cut could shift which side of `r_coll` a close
//!    approach lands on — flipping a discrete label rather than perturbing a continuous value.
//!    That is a different failure mode from the one measured for `spread_shape`, and the same
//!    reason for measuring it rather than assuming it is inert.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::outcome::State;

const SIZE: usize = 64;

fn render(r_coll_frac: f64, lc_stable: bool) -> Vec<PixelOut> {
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let cfg = EnsembleCfg { r_coll_frac, lc_stable, refine_flagged: false, ..Default::default() };
    (0..s.npix())
        .into_par_iter()
        .map(|i| evaluate::<f64>(&s, i, &cfg))
        .collect()
}

fn fracs(px: &[PixelOut]) -> [f64; 6] {
    let mut f = [0.0f64; 6];
    for p in px {
        f[(p.state as usize).min(5)] += 1.0;
    }
    f.map(|x| x / px.len() as f64)
}

fn triple_frac(px: &[PixelOut]) -> f64 {
    px.iter()
        .filter(|p| p.state == State::Collision as u8 && p.detail == 3)
        .count() as f64
        / px.len() as f64
}

fn main() {
    println!("near-field {SIZE}x{SIZE}, t=13, E+1=8, eta=0.01, f64");
    println!();
    println!("=== 1. r_coll sensitivity ===");
    println!("{:>10}{:>11}{:>11}{:>11}{:>11}{:>11}{:>12}{:>12}",
             "r_coll/R", "escape", "bounded", "collision", "running", "failed", "triples", "med t_end");

    let mut baseline: Option<Vec<PixelOut>> = None;
    for frac in [1e-4f64, 1e-3, 1e-2] {
        let px = render(frac, true);
        let f = fracs(&px);
        let mut te: Vec<f64> = px.iter().map(|p| p.t_end).filter(|x| x.is_finite()).collect();
        te.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("{frac:>10.0e}{:>11.4}{:>11.4}{:>11.4}{:>11.4}{:>11.4}{:>12.4}{:>12.4}",
                 f[State::Escape as usize], f[State::Bounded as usize],
                 f[State::Collision as usize], f[State::Running as usize],
                 f[State::SimFailed as usize], triple_frac(&px), te[te.len() / 2]);
        if frac == 1e-3 {
            baseline = Some(px);
        }
    }
    println!();
    println!("r_coll is a threshold on a monotone quantity, so raising it can only turn");
    println!("non-collisions into collisions. The table shows how fast - and the answer is");
    println!("that on this slice there is NO PLATEAU to sit on. The collision fraction goes");
    println!("0.0000 -> 0.0242 -> 1.0000 across three decades, because every pixel's d_min");
    println!("lands inside a single decade (see 1b). No choice of r_coll in that range is a");
    println!("physical event threshold; it is a readout of the d_min distribution. d_min_true");
    println!("is the primary quantity here and the collision label is derived from it.");
    println!();
    println!("Note also that `triples` is 0.0000 at every threshold: no pixel on this slice");
    println!("ever has two pairs below r_coll at once, so the >=2-pair rule is exercised only");
    println!("by construction in tests/outcome_encoding.rs, never by this data.");

    // The collision fractions are a threshold applied to d_min, so they must be reproducible
    // from the d_min distribution of an *unterminated* run. If they are not, the in-loop
    // sampling and the d_min tracking disagree and one of them is wrong.
    println!();
    println!("=== 1b. cross-check: the same fractions from d_min directly ===");
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let cfg = EnsembleCfg { stop_on_event: false, r_coll_frac: 0.0, refine_flagged: false, ..Default::default() };
    let free: Vec<PixelOut> = (0..s.npix())
        .into_par_iter()
        .map(|i| evaluate::<f64>(&s, i, &cfg))
        .collect();
    let mut dm: Vec<f64> = free.iter().map(|p| p.d_min_true).filter(|x| x.is_finite()).collect();
    dm.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let r0 = 2.2360679774997896f64; // hyperradius of the Burrau configuration at t=0
    println!("d_min_true/R over the grid: min {:.3e}  median {:.3e}  max {:.3e}",
             dm[0] / r0, dm[dm.len() / 2] / r0, dm[dm.len() - 1] / r0);
    for frac in [1e-4f64, 1e-3, 1e-2] {
        let below = dm.iter().filter(|d| **d / r0 < frac).count();
        println!("  fraction with d_min_true < {frac:.0e} R: {:.4}", below as f64 / dm.len() as f64);
    }
    println!("These are per-pixel minima over the ensemble, so they bound the nominal-copy");
    println!("fractions above from one side rather than matching them exactly.");

    println!();
    println!("=== 2. does the LC branch move the outcome fractions? ===");
    let base = baseline.unwrap();
    println!("{:>10}{:>11}{:>11}{:>11}{:>12}{:>14}",
             "r_coll/R", "escape", "bounded", "collision", "triples", "label flips");
    for frac in [1e-4f64, 1e-3, 1e-2] {
        let stable = if frac == 1e-3 { base.clone() } else { render(frac, true) };
        let unstable = render(frac, false);
        let f = fracs(&unstable);
        let flips = stable
            .iter()
            .zip(unstable.iter())
            .filter(|(a, b)| a.outcome != b.outcome)
            .count();
        println!("{frac:>10.0e}{:>11.4}{:>11.4}{:>11.4}{:>12.4}{:>10} / {}",
                 f[State::Escape as usize], f[State::Bounded as usize],
                 f[State::Collision as usize], triple_frac(&unstable), flips, stable.len());
    }
    println!();
    println!("Rows are the UNSTABLE branch; 'label flips' counts pixels whose packed outcome");
    println!("differs from the stable branch's at the same r_coll. A nonzero count means the");
    println!("branch cut reaches the outcome encoding, not only the continuous fields.");
}
