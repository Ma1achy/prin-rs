//! **Is `EnsembleCfg::production`'s `max_steps` still right, now that the default integrator is
//! Heggie?**
//!
//! `max_steps: 30_000` is on this record as *"sized for AZ"*, quoted from `integrator_gallery`,
//! which raises the budget to **400_000 for both arms** and says why in its own header: *"at the
//! production max_steps = 30_000 Heggie exhausts the budget on 8.6% of `config_stability` and its
//! drift panel comes back dominated by the magenta veto set"*. 27 diagnostic harnesses raise it;
//! **every render and scheduler harness silently took 30_000**. So a footprint that is
//! undetermined in a gallery panel is determined in an integrator panel, and the difference was
//! invisible because nothing took it as an argument. The default moved to Heggie on 2026-09-02
//! and its step budget did not follow -- the same shape as `refine_flagged` spreading by copy and
//! `k_frac = 1.0` shipping as the default.
//!
//! # The design, and each choice is a defect this project has already paid for
//!
//! **A FIXED UNIFORM GRID, NOT A TREE.** The first pass at this read the flagged count off
//! `chart_gallery`'s live march, where the tree changes with the budget -- so the denominators ran
//! 29632, 22912, 23296, 22144 and the "fraction" was over a moving population. Here every rung
//! evaluates the identical `res^2` grid, so the denominator is a constant and a difference in the
//! numerator is a difference in the physics.
//!
//! **`substeps`, NEVER `secs`.** *Read `steps`, not `secs`* is on this record because a run under
//! load timed faster while doing more work. `total_substeps` is the machine-independent cost and
//! it is the column that decides whether raising the budget is affordable.
//!
//! **THE CONTROL IS A CHART THAT FLAGS NOTHING.** If a chart reading 0 at 30_000 moves at all
//! across the ladder, the sweep is measuring something other than the budget and every other row
//! is suspect. `latent_shape` is that arm and it is asserted, not merely printed: its whole row
//! must be bitwise constant.
//!
//! **AND THE TREE IS MEASURED BESIDE THE FOOTPRINTS**, because the damage is not only cosmetic. A
//! truncated footprint reads as *undetermined*, so its quad can never resolve, so it splits -- and
//! a too-small budget therefore MANUFACTURES work rather than saving it. Measured on
//! `latent_mixed_h3`: the static tree runs 214 leaves at 30_000 against 142 at 480_000, and the
//! whole chart ran 2.5x faster at the larger budget. That inflation lands in leaf counts and
//! stop-reason breakdowns, which are the numbers the gallery table exists to report.
//!
//! Run: `cargo run --release --example step_budget -- [charts=all] [res=64] [ladder=30000,...]`
//! Writes nothing. `all` is the 26 gallery charts plus the 8 Burrau regions plus
//! `config_stability`.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

struct Target {
    name: String,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
}

fn targets(only: &Option<Vec<String>>) -> Vec<Target> {
    let mut v: Vec<Target> = Vec::new();
    for &(n, cx, cy, body) in grid::REGIONS.iter() {
        v.push(Target { name: n.into(), chart: Chart::BodyPlane, cx, cy, half: 0.05, body });
    }
    let (c, cx, cy, half) = Chart::config_stability();
    v.push(Target { name: "config_stability".into(), chart: c, cx, cy, half, body: 0 });
    for (n, c, cx, cy, half) in grid::gallery_cases() {
        v.push(Target { name: n.into(), chart: c, cx, cy, half, body: 0 });
    }
    match only {
        None => v,
        Some(want) => v.into_iter().filter(|t| want.iter().any(|w| w == &t.name)).collect(),
    }
}

/// One cell: the identical grid at one budget.
struct Row {
    flagged: usize,
    budget: usize,
    drift50: f64,
    drift99: f64,
    err99: f64,
    steps50: f64,
    steps_tot: u64,
}

fn cell(t: &Target, res: usize, ms: usize) -> Row {
    let ens = EnsembleCfg { max_steps: ms, ..EnsembleCfg::production() };
    let sl = grid::Slice::body_plane(res, res, t.cx, t.cy, t.half, t.body).with_chart(t.chart);
    let px: Vec<PixelOut> =
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();
    let flagged = px.iter().filter(|p| p.n_nonfinite > 0).count();
    let budget = px.iter().filter(|p| p.budget_exhausted).count();
    // Non-finite drift is an outcome, not missing data -- it is kept out of the quantile and
    // counted by `flagged` instead, which is the column that would otherwise silently shrink.
    let mut d: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
    let mut e: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
    let mut s: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
    Row {
        flagged,
        budget,
        drift50: q(&mut d, 0.50),
        drift99: q(&mut d, 0.99),
        err99: q(&mut e, 0.99),
        steps50: q(&mut s, 0.50),
        steps_tot: px.iter().map(|p| p.total_substeps).sum(),
    }
}

fn main() {
    let only: Option<Vec<String>> = std::env::args()
        .nth(1)
        .filter(|s| s != "all")
        .map(|s| s.split(',').map(|x| x.trim().to_string()).collect());
    let res: usize = arg(2, 64);
    let ladder: Vec<usize> = std::env::args()
        .nth(3)
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![30_000, 60_000, 120_000, 240_000, 480_000]);

    let ts = targets(&only);
    let n = res * res;
    println!(
        "step budget sweep: {} targets x {} rungs, a FIXED {res}^2 grid per rung ({n} footprints), \
         f64, production kernel otherwise.",
        ts.len(),
        ladder.len()
    );
    println!("  integrator: {:?}   ladder: {:?}", EnsembleCfg::production().integrator, ladder);
    println!(
        "  `flagged` is n_nonfinite>0 -- what the render used to paint magenta. `budget` is \
         budget_exhausted.\n  `steps` is total_substeps, the machine-independent cost: read steps, \
         not secs."
    );
    println!(
        "\n{:>18} {:>9} {:>8} {:>8} {:>10} {:>10} {:>10} {:>10} {:>11} {:>12} {:>8}",
        "case", "max_steps", "flagged", "budget", "flag%", "drift p50", "drift p99", "err p99",
        "steps p50", "steps TOTAL", "x base"
    );

    let mut control: Option<(String, Vec<(usize, usize, u64)>)> = None;
    for t in &ts {
        let mut sig: Vec<(usize, usize, u64)> = Vec::new();
        // **`steps TOTAL` against the coarsest rung is THE cost column, and `steps p50` is not.**
        // The budget caps only the tail, so the median trajectory is untouched by construction --
        // it reads identical at every rung on every target here, which is a fact about which
        // footprints the cap reaches and says nothing at all about what raising it costs. Printing
        // it alone would be *a statistic that cannot move being read as a statistic that did not
        // move*.
        let mut base: Option<u64> = None;
        for &ms in &ladder {
            let r = cell(t, res, ms);
            sig.push((r.flagged, r.budget, r.steps_tot));
            let b = *base.get_or_insert(r.steps_tot);
            println!(
                "{:>18} {ms:>9} {:>8} {:>8} {:>9.3}% {:>10.3e} {:>10.3e} {:>10.3e} {:>11.0} \
                 {:>12} {:>7.3}x",
                t.name,
                r.flagged,
                r.budget,
                100.0 * r.flagged as f64 / n as f64,
                r.drift50,
                r.drift99,
                r.err99,
                r.steps50,
                r.steps_tot,
                if b > 0 { r.steps_tot as f64 / b as f64 } else { f64::NAN }
            );
        }
        // The first target that flags nothing at the coarsest rung is the control.
        if control.is_none() && sig[0].0 == 0 {
            control = Some((t.name.clone(), sig));
        }
    }

    // **The arm that says the sweep measures the budget and not something else.** A chart the
    // budget never binds on must be bitwise constant across the whole ladder -- same flagged
    // count, same budget count, same total substeps. If it moves, every other row is suspect.
    match control {
        None => println!(
            "\n  NO CONTROL: every target flags at the coarsest rung, so nothing here says the \
             sweep is measuring the budget rather than the ladder. Add a tame chart."
        ),
        Some((name, sig)) => {
            let flat = sig.iter().all(|x| *x == sig[0]);
            println!(
                "\n  control {name}: {} across the ladder -- flagged/budget/steps {:?}",
                if flat { "CONSTANT" } else { "*** MOVED, the sweep is not clean" },
                sig[0]
            );
            assert!(flat, "the control moved across the ladder: {sig:?}");
        }
    }
}
