//! **Is A's reject-and-retry viable on a GPU, or is it a divergent branch?**
//!
//! A rejection is a per-lane branch. A warp executes in lockstep, so its cost is the **max over
//! its lanes**, not the mean: one lane retrying eight times stalls thirty-one others. A mode that
//! is cheap on CPU and divergent on GPU is not viable, and CPU wall clock hides that completely.
//!
//! # The control that makes the number mean anything
//!
//! Step counts already vary from lane to lane without any retry at all — a close encounter costs
//! more than a quiet one. So the absolute divergence factor is large for **every** mode, including
//! the baseline, and quoting it alone would condemn a mode for the field's own structure.
//! **What is reported is the factor under each candidate against the factor under `None`.** The
//! increase is A's cost; the level is the field's.
//!
//! # Two groupings, because dispatch order changes the answer
//!
//! Warps are formed **raster-linear** (32 consecutive pixels) and **8x4 tiled**, the usual
//! dispatch shape. Assuming one of them is how this measurement would go wrong: retries cluster
//! spatially, so a grouping that follows the structure and one that cuts across it give different
//! answers, and both are printed.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::StepLimit;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const WARP: usize = 32;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// Warp membership under the two dispatch shapes.
fn warps(res: usize, tiled: bool) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if !tiled {
        let mut cur = Vec::with_capacity(WARP);
        for i in 0..res * res {
            cur.push(i);
            if cur.len() == WARP {
                out.push(std::mem::take(&mut cur));
            }
        }
        if !cur.is_empty() {
            out.push(cur);
        }
        return out;
    }
    // 8 wide by 4 tall.
    for ty in (0..res).step_by(4) {
        for tx in (0..res).step_by(8) {
            let mut cur = Vec::with_capacity(WARP);
            for y in ty..(ty + 4).min(res) {
                for x in tx..(tx + 8).min(res) {
                    cur.push(y * res + x);
                }
            }
            if !cur.is_empty() {
                out.push(cur);
            }
        }
    }
    out
}

/// `mean(max per warp) / mean(per lane)`. 1.0 is perfect coherence; 32.0 is one lane doing all
/// the work.
fn divergence(cost: &[f64], groups: &[Vec<usize>]) -> f64 {
    let lane_mean = cost.iter().sum::<f64>() / cost.len() as f64;
    let warp_mean = groups
        .iter()
        .map(|g| g.iter().map(|&i| cost[i]).fold(0.0, f64::max))
        .sum::<f64>()
        / groups.len() as f64;
    warp_mean / lane_mean
}

fn main() {
    let res: usize = arg(1, 128);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    let (chart, cx, cy, half) = Chart::config_stability();
    let base = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);

    println!(
        "config_stability {res}^2, one lane per PIXEL (the copies are serial within a pixel here,\n\
         so this models a pixel-per-thread dispatch and says so rather than implying a lane is a\n\
         copy). Warp = {WARP}.\n\
         config: {}\n",
        base.provenance()
    );

    let lin = warps(res, false);
    let tiled = warps(res, true);

    println!(
        "  {:>22} {:>10} {:>10} {:>10} {:>11} {:>11} {:>10} {:>10}",
        "mode", "secs", "retry p90", "retry max", "div linear", "div tiled", "warps hit", "steps p50"
    );

    let mut base_div = (f64::NAN, f64::NAN);
    for (label, mode, f) in [
        ("None (control)", StepLimit::None, 0.0),
        ("Reject f=0.5", StepLimit::Reject, 0.5),
        ("Reject f=0.1", StepLimit::Reject, 0.1),
        ("Reject f=0.02", StepLimit::Reject, 0.02),
    ] {
        let cfg = EnsembleCfg { step_limit: mode, step_limit_f: f, ..base };
        let t = std::time::Instant::now();
        let px: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
            .collect();
        let secs = t.elapsed().as_secs_f64();

        // Cost is TOTAL STEPS, not retries: a retry's cost is the step it repeats, and a warp
        // stalls on total work rather than on how that work was labelled.
        let cost: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
        let mut rt: Vec<f64> = px.iter().map(|p| p.n_retry as f64).collect();
        let (dl, dt) = (divergence(&cost, &lin), divergence(&cost, &tiled));
        if mode == StepLimit::None {
            base_div = (dl, dt);
        }
        let hit = lin.iter().filter(|g| g.iter().any(|&i| px[i].n_retry > 0)).count() as f64
            / lin.len() as f64;
        let mut st: Vec<f64> = cost.clone();
        println!(
            "  {label:>22} {secs:>10.1} {:>10.0} {:>10.0} {dl:>11.3} {dt:>11.3} {hit:>10.4} {:>10.3e}",
            q(&mut rt.clone(), 0.9),
            q(&mut rt, 1.0),
            q(&mut st, 0.5)
        );
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **The absolute divergence factor is the FIELD's, not the mode's.** Step counts vary\n\
         lane to lane with or without retries, so the control row is large too. Read each\n\
         candidate's factor against the control's ({:.3} linear, {:.3} tiled): the INCREASE is\n\
         what reject-and-retry costs a warp.\n\n\
         `warps hit` is the fraction of warps containing at least one retrying lane. It is the\n\
         other half of the picture: rare retries that are perfectly scattered stall every warp,\n\
         and frequent retries that cluster stall few. A low retry rate with a high `warps hit` is\n\
         the bad case, and it is invisible in any per-pixel average.",
        base_div.0, base_div.1
    );
}
