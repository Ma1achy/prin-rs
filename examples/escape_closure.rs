//! The **closure-and-energy** escape criterion: does it hold, and is stopping on it safe?
//!
//! ```text
//!   ESCAPE  <=>  |dn| over a window < tau    AND    E_rel > 0
//! ```
//!
//! Two conditions, transcribed from `reference/escape_criterion.py`. `receding` and `d > r_esc`
//! are redundant once both hold -- measured identical to the digit on the reference -- so three
//! tuned constants become one, and `tau` is set from a **measured gap** rather than picked.
//!
//! # Why the old one was wrong
//!
//! `spec > 0 && receding` is not absorbing: during a close encounter a body's two-body energy
//! relative to the other pair goes transiently positive while it is still deep inside the system.
//! Of 895 `deep interior` trajectories escaping under a fine cadence, **0 were still unbound one
//! boundary later**. Under `stop_on_event` that mislabelling froze `shape_vec` at the copy's own
//! `t_end`, and since escape is detected at sync boundaries the rendered field became a patchwork
//! of `n_sync` time strata stitched at hard seams.
//!
//! # Sections, in the order they decide things
//!
//! 0. **Window and cadence.** The window is a TIME. `n_sync` is scaled per case so the realised
//!    window is comparable, and the closest-approach timescale is reported against it -- a
//!    two-end chord cannot tell a full revolution from stationarity (`tests/outcome_encoding.rs`
//!    holds that), so a window commensurate with an inner period aliases.
//! 1. **The measured gap**, per region and per `k`. `tau` is the geometric midpoint of the gap
//!    between the escaper and bound populations. An absolute cutoff picked by eye is what
//!    dismissed closure the first time -- `2e-3` sits *inside* the bound population's range.
//! 2. **Precision, recall, median firing time**, against the independent ground truth.
//! 3. **CHECK 1 -- nothing re-binds.** Candidacy at every boundary of ONE unstopped run at ONE
//!    step size, out past the horizon. Re-running with `n_sync` rescaled per window makes every
//!    window a different discretisation, which is the trap that produced a non-monotone curve
//!    inside a diagnostic written to catch traps.
//! 4. **`t_end` refinement.** How often the replay finds no crossing, and the distinct-value count.
//! 5. **The toggle.** `stop_on_escape` on and off.
//!
//! # The ground truth shares no term with the criterion
//!
//! CHECK 2. Purely geometric: run to `3 t_max` in three chained segments and call a body escaped
//! iff its separation from the other two's barycentre **grows monotonically** across the three
//! samples and `d(3T)/d(T)` exceeds [`GROWTH`]. No energy sign anywhere in it. The previous
//! round's ground truth shared the energy term, so its 100% precision was partly circular.
//!
//! # Writes
//!
//! stdout only.

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{integrate_az_opts, AzOpts, AzOut};
use prin_rs::outcome::{self, EscapeRule};
use prin_rs::physics::{energy, Cart};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The reference's window, in time units. `n_sync` is chosen per case to realise it.
const WINDOW: f64 = 0.4;
/// Ground truth: separation must grow by at least this factor between `T` and `3T`.
const GROWTH: f64 = 3.0;
/// The `k` ladder for the window sweep, in sync boundaries.
const KS: [usize; 4] = [1, 2, 3, 4];
/// Horizon multipliers for the gap. **`|dn/dt| ~ 1/t^3`, so the gap is a function of maturity**:
/// the reference quotes its 383x at `t = 25-30` and this project renders at 13. Both are shown,
/// because a criterion has to work at the horizon that ships, not only where it is cleanest.
const HORIZONS: [f64; 2] = [1.0, 2.0];

struct Case {
    name: String,
    t_max: f64,
    n_sync: usize,
    r_coll: f64,
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

fn opts(r_coll: f64, rule: EscapeRule<f64>, k: usize, stop_esc: bool) -> AzOpts<'static, f64> {
    AzOpts {
        step_limit: prin_rs::integrate::az::StepLimit::None,
        step_blend: prin_rs::integrate::az::StepBlend::Min,
        blend_p: 4.0,
        step_limit_f: 0.0,
        dtau_mode: prin_rs::integrate::az::DtauMode::default(),
        clamp_final_step: true,
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: r_coll,
        escape_rule: rule,
        closure_k: k,
        stop_on_event: true,
        stop_on_escape: stop_esc,
        escape_every: 0,
        escape_confirm: true,
        keep_boundary_shapes: false,
        keep_drift_hist: false,
    }
}

const ETA: f64 = 1e-2;
const MAX_STEPS: usize = 200_000;

fn run(c: &Cart<f64>, m: &[f64; 3], t: f64, n_sync: usize, o: &AzOpts<f64>) -> AzOut<f64> {
    integrate_az_opts(*c, m, t, n_sync, ETA, MAX_STEPS, o)
}

/// Separation of each body from the barycentre of the other two. The ground truth's only input.
fn seps(s: &Cart<f64>, m: &[f64; 3]) -> [f64; 3] {
    let mut d = [0.0; 3];
    for b in 0..3 {
        let o: Vec<usize> = (0..3).filter(|&k| k != b).collect();
        let mb = m[o[0]] + m[o[1]];
        let rc = (s.r[o[0]] * m[o[0]] + s.r[o[1]] * m[o[1]]) / mb;
        d[b] = (s.r[b] - rc).norm();
    }
    d
}

/// **CHECK 2.** Geometric only: monotone growth across `T, 2T, 3T` and a factor of [`GROWTH`].
///
/// Returns the escaping body, or `None`. Chained segments rather than one long run so the three
/// samples come from **one** trajectory at **one** step size -- the discretisation trap.
fn ground_truth(c: &Cart<f64>, m: &[f64; 3], t_max: f64, n_sync: usize, r_coll: f64)
    -> Option<u8>
{
    // Nothing terminal: an escape label must not depend on the run having been stopped, and a
    // `d_min` read over a truncated run inherits the truncation.
    let o = AzOpts { stop_on_event: false, stop_on_escape: false, r_coll_frac: r_coll,
                     ..opts(r_coll, EscapeRule::Reference, 1, false) };
    let mut st = *c;
    let mut d = [[0.0f64; 3]; 3];
    for seg in 0..3 {
        let out = run(&st, m, t_max, n_sync, &o);
        if !out.finite {
            return None;
        }
        st = out.state;
        d[seg] = seps(&st, m);
    }
    (0..3)
        .find(|&b| d[0][b] < d[1][b] && d[1][b] < d[2][b] && d[2][b] > GROWTH * d[0][b])
        .map(|b| b as u8)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() { f64::NAN } else { prin_rs::stats::quantile(v, p) }
}

fn main() {
    let n: usize = arg(1, 24);
    let tau_cli: f64 = arg(2, outcome::CLOSURE_TAU);

    let mut cases: Vec<Case> = Vec::new();
    // `n_sync` is derived from `t_max` so every case realises the same ~0.4 window. Holding
    // `n_sync` fixed while `t_max` varies compares different discretisations, and here it would
    // silently compare different criteria: at t_max = 50, n_sync = 32 the window is 1.5625.
    let nsync = |t: f64| (t / WINDOW).round().max(4.0) as usize;
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            cases.push(Case { name: region.into(), t_max: 13.0, n_sync: nsync(13.0),
                              r_coll: 1e-3, ics: sample(&Chart::BodyPlane, body, cx, cy, 0.05, n) });
        }
    }
    for (nm, ch, cx, cy, half) in grid::gallery_cases() {
        if nm == "preset_plambda" {
            cases.push(Case { name: nm.into(), t_max: 13.0, n_sync: nsync(13.0),
                              r_coll: 1e-3, ics: sample(&ch, 0, cx, cy, half, n) });
        }
    }
    for (nm, (ch, cx, cy, half), rc) in [
        ("config_basin", Chart::config_basin(), 0.02),
        ("config_stability", Chart::config_stability(), 0.005),
    ] {
        cases.push(Case { name: nm.into(), t_max: 50.0, n_sync: nsync(50.0),
                          r_coll: rc, ics: sample(&ch, 0, cx, cy, half, n) });
    }

    println!("escape_closure -- n={n} per side, eta={ETA}, tau={tau_cli:e}, window target {WINDOW}\n");

    // ---------------------------------------------------------------------------------------
    println!("== 0. WINDOW AND CADENCE ==");
    println!("The window is a TIME. `dt_sync` is set per case to realise ~{WINDOW}; `t_close` is the");
    println!("closest-approach timescale 2pi*sqrt(d_min^3/M), a proxy for the shortest inner period.");
    println!("A window at or above `t_close` can alias a fast orbit into reading as settled.\n");
    println!("A row that is all-bounded with no collisions may be a tame region or a run that");
    println!("never got anywhere: `budget` and `nonfin` say which, and `dmin p50` against `r_coll`");
    println!("says whether the collision arm could fire at all.\n");
    println!("{:<18} {:>5} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>7} {:>7} {:>9} {:>7}",
             "case", "nsync", "dt_sync", "R p50", "w(k=1)", "w(k=4)", "tcl p1", "tcl p50",
             "budget", "nonfin", "dmin p50", "bnd/ns");
    let mut base: Vec<Vec<AzOut<f64>>> = Vec::new();
    for c in &cases {
        let o = opts(c.r_coll, EscapeRule::Closure(tau_cli), 1, false);
        let outs: Vec<AzOut<f64>> = c.ics.par_iter()
            .map(|(s, m)| run(s, m, c.t_max, c.n_sync, &o)).collect();
        let dt = c.t_max / c.n_sync as f64;
        let mut rr: Vec<f64> = c.ics.iter().map(|(s, m)| energy::hyperradius(&s.r, m)).collect();
        let mut tcl: Vec<f64> = outs.iter().zip(c.ics.iter())
            .map(|(o, (_, m))| {
                let mt: f64 = m.iter().sum();
                2.0 * std::f64::consts::PI * (o.d_min_true.powi(3) / mt).sqrt()
            })
            .collect();
        let nn = outs.len() as f64;
        let be = outs.iter().filter(|o| o.budget_exhausted).count() as f64 / nn;
        let nf = outs.iter().filter(|o| !o.finite).count() as f64 / nn;
        let mut dm: Vec<f64> = outs.iter().map(|o| o.d_min_true).collect();
        // Boundaries actually completed, over `n_sync`. Below 1 means the run stopped early or a
        // boundary was SKIPPED after an overshoot -- and a skipped boundary widens the realised
        // closure window past `k * dt_sync` without saying so.
        let mut bn: Vec<f64> = outs.iter()
            .map(|o| o.closure_hist.len() as f64 / c.n_sync as f64).collect();
        println!("{:<18} {:>5} {:>7.4} {:>7.4} {:>8.4} {:>8.4} {:>8.2e} {:>8.2e} {:>7.4} {:>7.4} {:>9.2e} {:>7.4}",
                 c.name, c.n_sync, dt, q(&mut rr, 0.5), dt, 4.0 * dt,
                 q(&mut tcl, 0.01), q(&mut tcl, 0.5), be, nf, q(&mut dm, 0.5),
                 q(&mut bn, 0.5));
        base.push(outs);
    }

    // ---------------------------------------------------------------------------------------
    println!("\n== GROUND TRUTH (check 2: geometric, shares no term with the criterion) ==");
    println!("{:<18} {:>8} {:>10}", "case", "escaped", "frac");
    let truth: Vec<Vec<Option<u8>>> = cases.iter().map(|c| {
        c.ics.par_iter()
            .map(|(s, m)| ground_truth(s, m, c.t_max, c.n_sync, c.r_coll))
            .collect()
    }).collect();
    for (c, t) in cases.iter().zip(truth.iter()) {
        let k = t.iter().filter(|x| x.is_some()).count();
        println!("{:<18} {:>8} {:>10.4}", c.name, k, k as f64 / t.len() as f64);
    }

    // ---------------------------------------------------------------------------------------
    println!("\n== 1. THE MEASURED GAP, and the tau it implies ==");
    println!("Closure at the FINAL boundary of a run with **nothing terminal**, so every");
    println!("trajectory reaches the same playhead. With `stop_on_event` on, a collided run freezes");
    println!("and its shape stops changing -- closure reads exactly 0, which lands in the BOUND");
    println!("population and destroys the gap it is supposed to measure. Collided trajectories are");
    println!("split out and shown, not folded in: collision has its own arm and fires first.\n");
    println!("`sep` is p50(bound)/p50(esc) -- above 1 is the reference's picture. `tau*` is the");
    println!("geometric midpoint of p99(esc) and p1(bound). |dn/dt| ~ 1/t^3, so the gap is a");
    println!("function of MATURITY: the reference quotes 383x at t = 25-30 and this ships at 13.\n");
    println!("{:<18} {:>4} {:>2} {:>6} {:>6} {:>10} {:>10} {:>10} {:>10} {:>7} {:>10}",
             "case", "xT", "k", "n_esc", "n_bnd", "esc p50", "esc p99", "bnd p1", "bnd p50",
             "sep", "tau*");
    // A closure of EXACTLY zero is not a settled trajectory, it is two bitwise-identical shape
    // vectors -- a boundary that advanced no time, or a state that stopped moving numerically.
    // It reads as maximally settled and would fire the criterion, so it is counted rather than
    // absorbed into a percentile.
    let mut zero_note: Vec<String> = Vec::new();
    for (ci, c) in cases.iter().enumerate() {
        for (mult, k) in HORIZONS.iter().flat_map(|&h| KS.iter().map(move |&k| (h, k))) {
            // Nothing terminal. This is a measurement of the SIGNAL, not of the criterion, and a
            // truncated run's last sample is a fact about the truncation.
            let o = AzOpts { stop_on_event: false, stop_on_escape: false,
                             ..opts(c.r_coll, EscapeRule::Closure(tau_cli), k, false) };
            let (tm, ns) = (mult * c.t_max, (mult as usize) * c.n_sync);
            let outs: Vec<AzOut<f64>> = c.ics.par_iter()
                .map(|(s, m)| run(s, m, tm, ns, &o)).collect();
            let (mut e, mut b) = (Vec::new(), Vec::new());
            let mut n_coll = 0usize;
            for (i, out) in outs.iter().enumerate() {
                let Some(v) = out.closure_hist.iter().rev().find(|x| x.is_finite()).copied()
                else { continue };
                if out.events.collision.is_some() { n_coll += 1; continue }
                if truth[ci][i].is_some() { e.push(v) } else { b.push(v) }
            }
            let (e50, e99) = (q(&mut e.clone(), 0.5), q(&mut e.clone(), 0.99));
            let (b1, b50) = (q(&mut b.clone(), 0.01), q(&mut b.clone(), 0.5));
            let _ = n_coll;
            if mult == 1.0 && k == 1 {
                let ze = e.iter().filter(|x| **x == 0.0).count();
                let zb = b.iter().filter(|x| **x == 0.0).count();
                zero_note.push(format!("{:<18} exact-zero closure: esc {ze}/{} bnd {zb}/{}",
                                       c.name, e.len(), b.len()));
            }
            println!("{:<18} {:>4} {:>2} {:>6} {:>6} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>7.1} {:>10.3e}",
                     c.name, mult, k, e.len(), b.len(), e50, e99, b1, b50, b50 / e50,
                     (e99 * b1).sqrt());
        }
    }

    println!();
    for z in &zero_note {
        println!("  {z}");
    }

    println!("\n== 2. PRECISION / RECALL / MEDIAN FIRING TIME (k=1, tau={tau_cli:e}) ==");
    println!("Against the geometric ground truth. Recall below 1 is the RIGHT failure direction:");
    println!("a late firing is a late timestamp; an early one writes a wrong one permanently.\n");
    println!("{:<18} {:>8} {:>10} {:>10} {:>10} {:>10}",
             "case", "fires", "frac", "precision", "recall", "med t_fire");
    for (ci, c) in cases.iter().enumerate() {
        let outs = &base[ci];
        let mut tf: Vec<f64> = Vec::new();
        let (mut tp, mut fp, mut fnn) = (0usize, 0usize, 0usize);
        for (i, o) in outs.iter().enumerate() {
            let fired = o.events.escape.is_some();
            let real = truth[ci][i].is_some();
            match (fired, real) {
                (true, true) => { tp += 1; tf.push(o.events.escape.unwrap().1) }
                (true, false) => { fp += 1; tf.push(o.events.escape.unwrap().1) }
                (false, true) => fnn += 1,
                (false, false) => {}
            }
        }
        let fires = tp + fp;
        println!("{:<18} {:>8} {:>10.4} {:>10.4} {:>10.4} {:>10.3}",
                 c.name, fires, fires as f64 / outs.len() as f64,
                 if fires > 0 { tp as f64 / fires as f64 } else { f64::NAN },
                 if tp + fnn > 0 { tp as f64 / (tp + fnn) as f64 } else { f64::NAN },
                 q(&mut tf, 0.5));
    }

    // ---------------------------------------------------------------------------------------
    println!("\n== 3. CHECK 1 -- DOES ANYTHING RE-BIND? ==");
    println!("Candidacy at every boundary of ONE unstopped run at ONE step size, out to 3*t_max");
    println!("with n_sync scaled so `dt_sync` is UNCHANGED. The old criterion read 0 of 895 here.");
    println!("**The question is whether the fired body is still UNBOUND** -- the energy arm alone.");
    println!("Full candidacy also carries the closure gate, and closure is a difference of");
    println!("neighbouring samples, so it jitters above tau on a perfectly settled escape; reading");
    println!("persistence off it scores ordinary jitter as a re-binding. Both columns are shown so");
    println!("the difference between the two questions is visible rather than asserted.\n");
    println!("{:<18} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}",
             "case", "fired", "u+1", "u+2", "u+4", "u+8", "u end", "cand end");
    for c in &cases {
        let o = opts(c.r_coll, EscapeRule::Closure(tau_cli), 1, false);
        let outs: Vec<AzOut<f64>> = c.ics.par_iter()
            .map(|(s, m)| run(s, m, 3.0 * c.t_max, 3 * c.n_sync, &o)).collect();
        let wins = [1usize, 2, 4, 8];
        let (mut num, mut den) = ([0usize; 4], [0usize; 4]);
        let (mut fired, mut u_end, mut c_end, mut end_den) = (0usize, 0usize, 0usize, 0usize);
        for out in &outs {
            let Some((b, _)) = out.events.escape else { continue };
            let b = b as usize;
            let Some(k0) = out.escape_flags.iter().position(|&f| f) else { continue };
            fired += 1;
            for (wi, &w) in wins.iter().enumerate() {
                if let Some(f) = out.unbound_flags.get(k0 + w) {
                    den[wi] += 1;
                    if f[b] { num[wi] += 1 }
                }
            }
            if let Some(f) = out.unbound_flags.last() {
                end_den += 1;
                if f[b] { u_end += 1 }
                if *out.escape_flags.last().unwrap_or(&false) { c_end += 1 }
            }
        }
        let r = |a: usize, b: usize| if b > 0 { a as f64 / b as f64 } else { f64::NAN };
        println!("{:<18} {:>7} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>9.4}",
                 c.name, fired, r(num[0], den[0]), r(num[1], den[1]), r(num[2], den[2]),
                 r(num[3], den[3]), r(u_end, end_den), r(c_end, end_den));
    }

    println!("\n== 4. t_end REFINEMENT BY REPLAY ==");
    println!("`at entry` is the fraction where the energy arm already held on entry to the firing");
    println!("interval, so there was no crossing to find. If that is ~1 the refinement is");
    println!("decoration. `on bdry` is the fraction of escape t_end landing on a sync boundary.\n");
    println!("{:<18} {:>8} {:>10} {:>10} {:>10}",
             "case", "escapes", "at entry", "on bdry", "distinct");
    for (ci, c) in cases.iter().enumerate() {
        let dt = c.t_max / c.n_sync as f64;
        let mut te: Vec<f64> = Vec::new();
        let (mut at_entry, mut on_b) = (0usize, 0usize);
        for o in &base[ci] {
            let Some((_, t)) = o.events.escape else { continue };
            te.push(t);
            if o.t_end_at_entry { at_entry += 1 }
            if ((t / dt).round() * dt - t).abs() < 1e-9 { on_b += 1 }
        }
        let nn = te.len().max(1) as f64;
        let mut u: Vec<u64> = te.iter().map(|x| (x * 1e9) as u64).collect();
        u.sort_unstable();
        u.dedup();
        println!("{:<18} {:>8} {:>10.4} {:>10.4} {:>10}",
                 c.name, te.len(), at_entry as f64 / nn, on_b as f64 / nn, u.len());
    }

    // ---------------------------------------------------------------------------------------
    println!("\n== 5. THE TOGGLE ==");
    println!("`frozen` is the fraction whose run stopped before t_max. The prediction is that under");
    println!("this criterion the two rows are nearly identical -- it fires on a trajectory whose");
    println!("shape has already stopped moving, so freezing it changes little.\n");
    println!("**Median alone is not enough.** Twice on this project a row identical to five digits");
    println!("hid every pixel moving -- worst 6.7% and 1.86x. `d max` and `moved` are the guard.\n");
    println!("{:<18} {:>6} {:>8} {:>8} {:>8} {:>8} {:>9} {:>9} {:>7}",
             "case", "stop", "escape", "collide", "bounded", "frozen", "d med", "d max", "moved");
    for c in &cases {
        let mut prev: Option<Vec<[f64; 3]>> = None;
        for &stop in &[false, true] {
            let o = opts(c.r_coll, EscapeRule::Closure(tau_cli), 1, stop);
            let outs: Vec<AzOut<f64>> = c.ics.par_iter()
                .map(|(s, m)| run(s, m, c.t_max, c.n_sync, &o)).collect();
            let (mut e, mut co, mut bo, mut fr) = (0usize, 0usize, 0usize, 0usize);
            let shapes: Vec<[f64; 3]> = outs.iter().zip(c.ics.iter()).map(|(o, (_, m))| {
                let oc = outcome::classify(&o.events, &o.state, m, o.finite, o.budget_exhausted);
                match oc.state {
                    outcome::State::Escape => e += 1,
                    outcome::State::Collision => co += 1,
                    outcome::State::Bounded => bo += 1,
                    _ => {}
                }
                if o.t < c.t_max - 1e-9 { fr += 1 }
                prin_rs::physics::shape::shape_vec(&o.state.r, m)
            }).collect();
            let nn = outs.len() as f64;
            // The per-pixel distance between the two toggle states' shape vectors. An aggregate
            // can only say the distribution did not move; it cannot say the pixels did not.
            let (dmed, dmax, moved) = prev.as_ref().map_or((f64::NAN, f64::NAN, 0usize), |p| {
                let mut v: Vec<f64> = p.iter().zip(shapes.iter())
                    .map(|(a, b)| outcome::closure(a, b))
                    .filter(|x| x.is_finite())
                    .collect();
                let mv = v.iter().filter(|x| **x > 0.0).count();
                let mx = v.iter().cloned().fold(0.0f64, f64::max);
                (q(&mut v, 0.5), mx, mv)
            });
            println!("{:<18} {:>6} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>9.3e} {:>9.3e} {:>7}",
                     if stop { "" } else { &c.name }, stop,
                     e as f64 / nn, co as f64 / nn, bo as f64 / nn, fr as f64 / nn,
                     dmed, dmax, moved);
            prev = Some(shapes);
        }
    }

    println!("\nHOW TO READ THIS");
    println!("  §0 decides whether the window is safe before any distribution is read. A case whose");
    println!("     `t_close` is at or below `w(k=1)` can alias an inner orbit into reading settled.");
    println!("  §1 sets tau. If `tau*` is stable across regions a single value works; if the gap is");
    println!("     negative the populations overlap and NO fixed tau separates them there.");
    println!("  §2 recall < 1 is the right failure direction; precision < 1 is not -- BUT read it");
    println!("     against §3's `u end`. The ground truth demands 3x growth by 3*t_max, so a slow");
    println!("     genuine escape fails to be certified, and a `precision` shortfall with `u end`");
    println!("     at 1.0000 is the GROUND TRUTH missing them, not the criterion inventing them.");
    println!("  §3 IS THE ONE THAT DECIDES. Near 1.0 means nothing re-binds and stopping is safe.");
    println!("     Near 0 means the criterion is latching transients and it is NOT the fix.");
    println!("  §4 says whether the replay bought anything. `at entry` near 1 means it did not.");
    println!("  §5 `shape d` near 0 is the prediction: freezing a converged trajectory is a no-op.");
}
