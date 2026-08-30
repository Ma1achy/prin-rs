//! **What are the pixels the predictive limit does NOT fix?**
//!
//! Under `StepLimit::Predictive` at `f = 0.02`, `config_stability` at 384^2 leaves
//! `err>10 = 0.0001` — of order fourteen pixels of 147456, in a thin diagonal streak rather than
//! scattered. A streak is structure; scattered pixels would be noise. This asks what they are.
//!
//! # The candidates, and the discriminator for each
//!
//! **MEASURED, AND IT IS NONE OF THESE.** All 18 residual pixels at 384^2 have **exactly one**
//! pair below `r_coll`, `d_min` between `4.43e-3` and `4.98e-3` against `r_coll = 5e-3` — a
//! *grazing* binary encounter, with the third body 350x further away — energy drift of only
//! `2e-4` to `9e-4`, no `ab_floored`, no non-finite copy, and terminal states split between
//! `Collision` and `Bounded`. They form a contiguous diagonal filament in the chart plane, not a
//! scatter. **`error_ratio` is a spread ACROSS COPIES, and here the copies straddle the collision
//! threshold** — some collide, some do not — so they stop in genuinely different states and their
//! energies differ for a physical reason rather than a numerical one. `n_outcome_disagree` is the
//! column that shows it, which is why it replaced the pair count.
//!
//! ```text
//!   triple collision      >=2 pairs below r_coll -- AZ regularises the two pairs adjacent to
//!                         the reference body and the triple singularity is NOT removable, so
//!                         no step size fixes it and it is a MEASUREMENT OUTCOME, not a defect
//!   binary collision      exactly one pair below r_coll and d_min very small
//!   TINY floor            ab_floored -- a fabricated denominator
//!   step budget           budget_exhausted / retry_exhausted
//!   residual overshoot    n_overshoot > 0 -- the limit failed to bound the step
//!   just chaos            none of the above, d_min unremarkable
//! ```
//!
//! **The `eta` arm is what separates the first from the last.** A triple collision does not
//! improve with a finer step — that is the standing diagnostic signature of this project, one
//! level down: *a quantity that does not converge under refinement is measuring the sampling
//! rather than the system*. So each residual pixel is re-run at `eta/4`, `eta/16` and `eta/64`,
//! and whether `error_ratio` falls is the answer.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, State, CLOSURE_TAU};
use prin_rs::physics::newton;

const WINDOW: f64 = 0.4;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 384);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    let (chart, cx, cy, half) = Chart::config_stability();
    let r_coll = 0.005;
    let cfg = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: r_coll,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", cfg.provenance());

    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    let bad: Vec<usize> = (0..px.len()).filter(|&i| px[i].error_ratio > 10.0).collect();
    println!(
        "residual err>10: {} of {} ({:.6})\n",
        bad.len(),
        px.len(),
        bad.len() as f64 / px.len() as f64
    );
    if bad.is_empty() {
        println!("nothing residual at this resolution -- the probe has no subject.");
        return;
    }

    println!(
        "  {:>5} {:>5} {:>11} {:>11} {:>10} {:>10} {:>10} {:>9} {:>7} {:>7} {:>6} {:>7}",
        "x", "y", "err_ratio", "drift", "d_min", "d_min/rcoll", "3rd sep", "state", "detail",
        "disagree", "nonfin", "floored"
    );
    for &i in &bad {
        let p = &px[i];
        let (x, y) = sl.decode_pos(i);
        let st = grid::decode_state(&chart, 0, x, y);
        let mut d = newton::pair_dists(&st.s.r);
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // How many pairs came within `r_coll` over the whole march -- the >=2-pair rule is the
        // project's triple test, and `d_min_true` is the tightest any pair reached.
        println!(
            "  {:>5} {:>5} {:>11.3e} {:>11.3e} {:>10.3e} {:>10.4} {:>10.3e} {:>9} {:>7} \
             {:>7} {:>6} {:>7}",
            i % res,
            i / res,
            p.error_ratio,
            p.energy_drift_max,
            p.d_min_true,
            p.d_min_true / r_coll,
            d[2],
            format!("{:?}", State::from_bits(p.state)),
            p.detail,
            p.n_outcome_disagree,
            p.n_nonfinite,
            p.ab_floored,
        );
    }

    // --- the discriminator ---------------------------------------------------------------
    println!(
        "\n== DOES IT CONVERGE? ==\n\
         A triple collision is a genuine singularity AZ does not remove, so `error_ratio` will\n\
         NOT fall with `eta`. Ordinary under-resolution will. **This is the whole question**, and\n\
         it is the project's own diagnostic signature: a quantity that does not converge under\n\
         refinement is measuring the sampling rather than the system.\n"
    );
    println!(
        "  {:>5} {:>5} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "x", "y", "eta", "eta/4", "eta/16", "eta/64", "verdict"
    );
    let mut n_conv = 0usize;
    let mut n_stuck = 0usize;
    for &i in &bad {
        let e: Vec<f64> = [1.0, 0.25, 0.0625, 0.015_625]
            .into_iter()
            .map(|s| pixel::evaluate_at::<f64>(&sl, i, &cfg, cfg.eta * s).error_ratio)
            .collect();
        // Converged means it reached the healthy value, not merely fell: `error_ratio` is
        // normalised to exactly 1.0 under exact dynamics, so it has an absolute target.
        let conv = e[3] <= 2.0;
        if conv {
            n_conv += 1;
        } else {
            n_stuck += 1;
        }
        println!(
            "  {:>5} {:>5} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>10}",
            i % res,
            i / res,
            e[0],
            e[1],
            e[2],
            e[3],
            if conv { "resolves" } else { "STUCK" }
        );
    }
    println!(
        "\n  resolves at eta/64: {n_conv}   stuck: {n_stuck}\n\n\
         A STUCK pixel is not a defect. `error_ratio` exists to say *this pixel is not data*, and\n\
         a trajectory that passes arbitrarily close to a triple collision is a measurement\n\
         outcome, not missing data -- the no-discard rule. A pixel that RESOLVES is ordinary\n\
         under-resolution and `f` could be tightened, at a cost the table in\n\
         `results/step_control/README.md` prices."
    );
}
