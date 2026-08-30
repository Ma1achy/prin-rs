//! The boundary-overshoot measurement: the clamp's convergence order, its interaction with
//! `dtau_mode`, and whether this slice is resolved at all.
//!
//! # The mechanism
//!
//! The march exits an interval by **overshooting** it -- the loop condition is
//! `s.t >= dt_left` -- and only the *clock* was corrected (`t += s.t.min(dt_left)`). The
//! Cartesian state written back at the boundary is the overshot one. That is a **first-order**
//! error injected at every one of `n_sync` boundaries, inside an RK4 march.
//!
//! # Why it is the partner of `DtauMode::PerStepInterval` and not an independent knob
//!
//! Under `FixedPerInterval` `dtau` is constant across the interval, so the overshoot is a fixed
//! slice of fictitious time and neighbouring trajectories overshoot alike: the error is large
//! but spatially **smooth**, and it displaces the picture without breaking it. Under
//! `PerStepInterval` the final step's size is a function of the local `A*B`, so the overshoot
//! becomes a function of local state and neighbouring pixels overshoot by *different* amounts.
//! A spatially-varying error injected at every boundary is measurably worse than a smooth one,
//! and §1 and §3 below measure exactly that.
//!
//! **The nested-arc banding this was first proposed to explain is NOT caused by it.** All four
//! arms carry it, including the one predating both changes, and under outcome-class colouring it
//! vanishes -- a colouring artefact, per `RESULTS §21`. See `RESULTS §24.8`. The defect measured
//! here is real; it is not the cause of that appearance.
//!
//! # The four arms
//!
//! ```text
//!   A  dtau fixed     + overshoot present    the original committed behaviour
//!   B  dtau per-step  + overshoot present    the regression
//!   C  dtau fixed     + overshoot clamped
//!   D  dtau per-step  + overshoot clamped    the proposed default
//! ```
//!
//! # The order matters
//!
//! §1 is the **cheapest confirmation the clamp does what it claims** and it is a convergence
//! order, not an error: an error can fall for any number of reasons, and only the *order* says
//! the leading term changed. It runs in under a second on one closed orbit and needs no grid.
//!
//! §4 is the question underneath all of it. If the image keeps moving as `eta` falls, this slice
//! is not resolved at the current settings and no single bug is the answer. Horizon 50 is close
//! to the f64 predictability horizon for this system, so that is a live possibility.

use rayon::prelude::*;

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{integrate_az_opts, AzOpts, AzOut, DtauMode};
use prin_rs::physics::{shape::shape_vec, Cart};
use prin_rs::Vec2;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const MAX_STEPS: usize = 200_000;
const HOT: f64 = 1e-6;

/// The four arms, in the order the write-up quotes them.
const ARMS: [(&str, DtauMode, bool); 4] = [
    ("A fixed   +overshoot", DtauMode::FixedPerInterval, false),
    ("B perstep +overshoot", DtauMode::PerStepInterval, false),
    ("C fixed   +clamp    ", DtauMode::FixedPerInterval, true),
    ("D perstep +clamp    ", DtauMode::PerStepInterval, true),
];

fn opts<'a>(mode: DtauMode, clamp: bool, r_coll: f64) -> AzOpts<'a, f64> {
    AzOpts {
        step_limit: prin_rs::integrate::az::StepLimit::None,
        step_blend: prin_rs::integrate::az::StepBlend::Min,
        blend_p: 4.0,
        step_limit_f: 0.0,
        r_coll_frac: r_coll,
        // Nothing terminal. A run stopped early has a shorter step count and a smaller drift for
        // reasons that have nothing to do with the step control.
        stop_on_event: false,
        stop_on_escape: false,
        dtau_mode: mode,
        clamp_final_step: clamp,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------------------------
// §1  The figure-eight
// ---------------------------------------------------------------------------------------------

/// Chenciner-Montgomery, equal masses, `G = 1`. The orbit is exactly periodic, so the distance
/// between the state at `T` and the state at `0` is a pure error measure with no reference
/// trajectory to compute and no chaos to contaminate it.
fn figure_eight() -> (Cart<f64>, [f64; 3], f64) {
    let x = 0.97000436;
    let y = -0.24308753;
    let vx = -0.93240737;
    let vy = -0.86473146;
    let r = [
        Vec2::new(x, y),
        Vec2::new(0.0, 0.0),
        Vec2::new(-x, -y),
    ];
    let v = [
        Vec2::new(-vx / 2.0, -vy / 2.0),
        Vec2::new(vx, vy),
        Vec2::new(-vx / 2.0, -vy / 2.0),
    ];
    (Cart::new(r, v), [1.0, 1.0, 1.0], 6.325_913_98)
}

fn closure_err(a: &Cart<f64>, b: &Cart<f64>) -> f64 {
    let mut m = 0.0f64;
    for i in 0..3 {
        m = m.max((a.r[i].x - b.r[i].x).abs()).max((a.r[i].y - b.r[i].y).abs());
        m = m.max((a.v[i].x - b.v[i].x).abs()).max((a.v[i].y - b.v[i].y).abs());
    }
    m
}

fn section_1() {
    let (s0, m, period) = figure_eight();
    let etas = [0.02, 0.01, 0.005, 0.002, 0.001];
    let n_sync = 32usize;

    println!("== 1. THE FIGURE-EIGHT -- CONVERGENCE ORDER ==");
    println!(
        "Chenciner-Montgomery at equal masses, integrated over exactly one period\n\
         (T = {period}), n_sync = {n_sync} fixed so the number of boundaries -- and so the\n\
         number of overshoots -- is the same at every `eta`. `closure` is the max component\n\
         difference between the state at T and the state at 0.\n\n\
         **READ THE ORDER, NOT THE ERROR.** An error falls for many reasons; only the order says\n\
         the leading term changed. Without the clamp the overshoot is O(h) and dominates an RK4\n\
         march; with it, the boundary contributes at the stepper's own order.\n"
    );
    println!("{:<22} {:>8} {:>12} {:>12} {:>8} {:>9}", "arm", "eta", "closure", "drift", "order", "steps");
    for (label, mode, clamp) in ARMS {
        let mut prev: Option<(f64, f64)> = None;
        let mut ends: Vec<(f64, f64)> = Vec::new();
        for eta in etas {
            let o = opts(mode, clamp, 0.0);
            let out: AzOut<f64> =
                integrate_az_opts(s0, &m, period, n_sync, eta, MAX_STEPS, &o);
            let e = closure_err(&out.state, &s0);
            let ord = match prev {
                Some((pe, pv)) if e > 0.0 && pe > 0.0 => {
                    format!("{:8.2}", (pe / e).ln() / (pv / eta).ln())
                }
                _ => format!("{:>8}", "--"),
            };
            println!(
                "{label:<22} {eta:>8.4} {e:>12.4e} {:>12.3e} {ord} {:>9}",
                out.drift, out.steps
            );
            prev = Some((e, eta));
            ends.push((eta, e));
        }
        // The pairwise orders are noisy -- each is a two-point estimate over a factor of two.
        // The endpoint-to-endpoint slope over the whole decade is the number to read.
        let (e0, v0) = (ends[0].1, ends[0].0);
        let (e1, v1) = (ends[ends.len() - 1].1, ends[ends.len() - 1].0);
        println!("{label:<22} {:>8} {:>12} {:>12} {:>8.2} {:>9}",
                 "overall", "", "", (e0 / e1).ln() / (v0 / v1).ln(), "");
        println!();
    }
}

// ---------------------------------------------------------------------------------------------
// §2, §3  Regions
// ---------------------------------------------------------------------------------------------

struct Case {
    name: String,
    t_max: f64,
    n_sync: usize,
    r_coll: f64,
    n: usize,
    ics: Vec<(Cart<f64>, [f64; 3])>,
}

fn sample(chart: &Chart, body: usize, cx: f64, cy: f64, half: f64, n: usize)
    -> Vec<(Cart<f64>, [f64; 3])>
{
    let mut out = Vec::with_capacity(n * n);
    for j in 0..n {
        for i in 0..n {
            let u = cx - half + 2.0 * half * (i as f64 + 0.5) / n as f64;
            let v = cy - half + 2.0 * half * (j as f64 + 0.5) / n as f64;
            let ic = grid::decode_state(chart, body, u, v);
            out.push((ic.s, ic.m));
        }
    }
    out
}

fn run(c: &Case, mode: DtauMode, clamp: bool, eta: f64) -> Vec<AzOut<f64>> {
    let o = opts(mode, clamp, c.r_coll);
    c.ics
        .par_iter()
        .map(|(s, m)| integrate_az_opts(*s, m, c.t_max, c.n_sync, eta, MAX_STEPS, &o))
        .collect()
}

fn q(v: &[f64], p: f64) -> f64 {
    let mut w: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if w.is_empty() { f64::NAN } else { prin_rs::stats::quantile(&mut w, p) }
}

/// The shape vector of a final state, or `None` where it is not defined.
fn shape_of(o: &AzOut<f64>, m: &[f64; 3]) -> Option<[f64; 3]> {
    let n = shape_vec(&o.state.r, m);
    n.iter().all(|x| x.is_finite()).then_some(n)
}

/// Pixels moved, and the median and max chord between the two arms' shape vectors.
///
/// **This is the convergence red flag.** A converged integrator does not move most of a field on
/// a step-control change. Where both arms are non-finite the pixel is not counted as moved --
/// it has no value under either -- and that count is printed separately rather than folded in.
fn compare(a: &[AzOut<f64>], b: &[AzOut<f64>], ics: &[(Cart<f64>, [f64; 3])]) -> (usize, usize, f64, f64) {
    let mut d: Vec<f64> = Vec::with_capacity(a.len());
    let mut undef = 0usize;
    for k in 0..a.len() {
        match (shape_of(&a[k], &ics[k].1), shape_of(&b[k], &ics[k].1)) {
            (Some(p), Some(qv)) => {
                let s: f64 = (0..3).map(|i| (p[i] - qv[i]).powi(2)).sum();
                d.push(s.sqrt());
            }
            _ => undef += 1,
        }
    }
    let moved = d.iter().filter(|x| **x > 0.0).count();
    (moved, undef, q(&d, 0.5), d.iter().cloned().fold(0.0f64, f64::max))
}

fn cases(n: usize) -> Vec<Case> {
    let mut cases: Vec<Case> = Vec::new();
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            cases.push(Case { name: region.into(), t_max: 13.0, n_sync: 33, r_coll: 1e-3, n,
                              ics: sample(&Chart::BodyPlane, body, cx, cy, 0.05, n) });
        }
    }
    let (cs, sx, sy, sh) = Chart::config_stability();
    cases.push(Case { name: "config_stability".into(), t_max: 50.0, n_sync: 125, r_coll: 0.005, n,
                      ics: sample(&cs, 0, sx, sy, sh, n) });
    cases
}

fn main() {
    let n: usize = arg(1, 48);
    let eta: f64 = arg(2, 1e-2);

    section_1();

    let cs = cases(n);

    println!("== 2. DRIFT AND NON-FINITE, PER ARM ==");
    println!(
        "{n}x{n} nominal trajectories per case at eta = {eta:e}, nothing terminal.\n\
         `hot` is the fraction above {HOT:e}. `budget` is the step-budget exhaustion count and\n\
         is SEPARATE from `nonfin` on purpose: a change that swapped one failure for the other\n\
         would have moved the magenta count and fixed nothing.\n"
    );
    println!("{:<18} {:<22} {:>11} {:>11} {:>8} {:>8} {:>8} {:>10}",
             "case", "arm", "drift p50", "drift p99", "hot", "nonfin", "budget", "steps p50");
    let mut all: Vec<(String, Vec<Vec<AzOut<f64>>>)> = Vec::new();
    for c in &cs {
        let mut per_arm = Vec::new();
        for (label, mode, clamp) in ARMS {
            let out = run(c, mode, clamp, eta);
            let dr: Vec<f64> = out.iter().map(|o| o.drift).collect();
            let hot = dr.iter().filter(|x| !(**x <= HOT)).count() as f64 / dr.len() as f64;
            let nf = out.iter().filter(|o| !o.finite).count();
            let bud = out.iter().filter(|o| o.budget_exhausted).count();
            let st: Vec<f64> = out.iter().map(|o| o.steps as f64).collect();
            println!("{:<18} {label:<22} {:>11.3e} {:>11.3e} {hot:>8.4} {nf:>8} {bud:>8} {:>10.0}",
                     c.name, q(&dr, 0.5), q(&dr, 0.99), q(&st, 0.5));
            per_arm.push(out);
        }
        println!();
        all.push((c.name.clone(), per_arm));
    }

    println!("== 3. HOW MUCH OF THE FIELD MOVES BETWEEN ARMS ==");
    println!(
        "Chord between final shape vectors, over the unit sphere (diameter 2). **A CONVERGED\n\
         INTEGRATION DOES NOT MOVE MOST OF A FIELD ON A STEP-CONTROL CHANGE.** The A->B figure\n\
         is the red flag; if the clamp is doing what it claims, C->D is much smaller. `undef`\n\
         counts pixels with no shape under one arm or the other -- not moved, absent.\n"
    );
    println!("{:<18} {:<12} {:>10} {:>8} {:>8} {:>11} {:>11}",
             "case", "pair", "moved", "frac", "undef", "chord p50", "chord max");
    let pairs = [(0usize, 1usize, "A->B"), (2, 3, "C->D"), (1, 3, "B->D"), (0, 2, "A->C"),
                 (0, 3, "A->D")];
    for (ci, c) in cs.iter().enumerate() {
        for (i, j, lab) in pairs {
            let (moved, undef, p50, mx) = compare(&all[ci].1[i], &all[ci].1[j], &c.ics);
            println!("{:<18} {lab:<12} {moved:>10} {:>8.4} {undef:>8} {p50:>11.3e} {mx:>11.3e}",
                     c.name, moved as f64 / c.ics.len() as f64);
        }
        println!();
    }

    println!("== 4. IS THIS SLICE RESOLVED AT ALL? ==");
    println!(
        "Arm D at eta, eta/2, eta/4. **If the image keeps moving as the step falls, no single\n\
         bug is the answer** -- the slice is unresolved at these settings and every render from\n\
         it is a picture of the discretisation as much as of the physics. Horizon 50 is close to\n\
         the f64 predictability horizon for this system, so this is a live possibility rather\n\
         than a formality. The figure to watch is whether `chord p50` FALLS between the rungs.\n"
    );
    println!("{:<18} {:<16} {:>10} {:>8} {:>11} {:>11} {:>11}",
             "case", "pair", "moved", "frac", "chord p50", "chord max", "drift p50");
    for c in &cs {
        let mut prev: Option<Vec<AzOut<f64>>> = None;
        let mut pe = eta;
        for k in 0..3 {
            let e = eta / (1 << k) as f64;
            let out = run(c, DtauMode::PerStepInterval, true, e);
            if let Some(p) = prev {
                let (moved, _, p50, mx) = compare(&p, &out, &c.ics);
                let dr: Vec<f64> = out.iter().map(|o| o.drift).collect();
                println!("{:<18} {:<16} {moved:>10} {:>8.4} {p50:>11.3e} {mx:>11.3e} {:>11.3e}",
                         c.name, format!("{pe:.2e}->{e:.2e}"),
                         moved as f64 / c.ics.len() as f64, q(&dr, 0.5));
            }
            pe = e;
            prev = Some(out);
        }
        println!();
    }

    println!("DONE");
}
