//! The escape gate's **missing distance condition**: does adding it kill the transients?
//!
//! # The difference, read from source
//!
//! ```text
//!   frag.glsl:104        dist > r_esc && outward > 0.0 && E_out > 0.0     THREE conditions
//!   tb_all_az.py:59-75   spec > 0 && dot(dr, dv) > 0                      TWO
//! ```
//!
//! This port transcribed the numpy form, so it had **no distance gate at all** and declared
//! escape the instant the relative energy went positive and the body was receding -- at any
//! separation, including mid-encounter. That is exactly the population §21.7 measured: of the
//! 895 `deep interior` trajectories that escape under an in-loop test and not at the reference
//! cadence, **0 of 895** were still unbound one boundary later. `r_esc` is the same persistence
//! guard done **geometrically** rather than temporally, and distance is monotone on a real
//! escape where a time window is a heuristic.
//!
//! # What is measured, in the order it decides things
//!
//! 0. **Units.** `r_esc` must be canonical -- a fraction of the initial hyperradius `R`, fixed
//!    at `t = 0` -- not the reference's absolute literal. Whether the literal `5` transfers is
//!    an arithmetic question about `R` on each chart, and it is measured here rather than
//!    assumed.
//! 1. **The transient population.** Of the escapes that fire, how many are still unbound at
//!    +1, +2, +4, +8 boundaries? It was 0 of 895 ungated. **If the gate does not take that near
//!    100%, the gate is not the fix and this example says so.**
//! 2. **Escape fraction per region, before and after, at every cadence.**
//! 3. (renders: `examples/escape_gate_render.rs`)
//! 4. **`r_esc` sensitivity.** It is a new constant, so it is a reported measurement and not a
//!    picked number -- the same treatment `r_coll` gets in `examples/r_coll_sweep.rs`.
//!
//! Persistence is read off `AzOut::escape_flags` -- instantaneous candidacy at every boundary
//! of **one** run at **one** step size. Re-running to `t_e + w` with `n_sync` rescaled makes
//! every window a different discretisation, which produced a non-monotone curve that was the
//! instrument rather than the physics.
//!
//! # Writes
//!
//! stdout only.

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{integrate_az_opts, AzOpts};
use prin_rs::physics::{energy, Cart};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

struct Case {
    name: String,
    t_max: f64,
    r_coll: f64,
    ics: Vec<(Cart<f64>, [f64; 3])>,
}

fn sample(chart: &Chart, body: usize, cx: f64, cy: f64, half: f64, n: usize)
    -> Vec<(Cart<f64>, [f64; 3])>
{
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let u = cx - half + 2.0 * half * (i as f64 + 0.5) / n as f64;
            let v = cy - half + 2.0 * half * (j as f64 + 0.5) / n as f64;
            let ic = grid::decode_state(chart, body, u, v);
            out.push((ic.s, ic.m));
        }
    }
    out
}

// **The `all_bodies` axis no longer exists as a knob.** The body set is now a property of the
// rule -- `Reference` labels only the body outside the tightest pair, `Distance` tests all three
// -- so each rule matches its own reference exactly. The separated-axis measurement this example
// made is on record: the body arm alone moved near-field 0.0000 -> 1.0000, so it is *not* free.
// `all` here selects the rule: `false` gives `Reference` (one body, no gate) and `true` gives
// `Distance(r_esc)` (all three), so `Distance(0.0)` is still the ungated all-bodies cell.
fn opts(r_coll: f64, r_esc: f64, all: bool, ev: usize, stop: bool) -> AzOpts<'static, f64> {
    AzOpts {
        land_iterate: true,
        land_max_iters: 4,
        // Zero, the shipped default. These runs are about the escape criterion; hysteresis
        // changes WHICH reference body is chosen at a boundary and would vary a second thing.
        ref_hysteresis: 0.0,
        step_limit: prin_rs::integrate::az::StepLimit::None,
        step_blend: prin_rs::integrate::az::StepBlend::Min,
        blend_p: 4.0,
        step_limit_f: 0.0,
        dtau_mode: prin_rs::integrate::az::DtauMode::default(),
        clamp_final_step: true,
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: r_coll,
        escape_rule: if r_esc > 0.0 || all {
            prin_rs::outcome::EscapeRule::Distance(r_esc)
        } else {
            prin_rs::outcome::EscapeRule::Reference
        },
        closure_k: 1,
        stop_on_escape: stop,
        stop_on_event: stop,
        escape_every: ev,
        // **Deliberately off.** The question is whether the geometric guard works; running the
        // temporal one alongside would measure the pair.
        escape_confirm: false,
        keep_boundary_shapes: false,
        keep_drift_hist: false,
    }
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { f64::NAN } else { prin_rs::stats::quantile(v, p) }
}

fn main() {
    let n: usize = arg(1, 32);
    let n_sync: usize = arg(2, 32);
    let eta = 1e-2f64;
    let esc_ladder: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0];
    let strides: Vec<usize> = vec![0, 4];

    let mut cases: Vec<Case> = Vec::new();
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            cases.push(Case {
                name: region.into(),
                t_max: 13.0,
                r_coll: 1e-3,
                ics: sample(&Chart::BodyPlane, body, cx, cy, 0.05, n),
            });
        }
    }
    for (nm, ch, cx, cy, half) in grid::gallery_cases() {
        if nm == "preset_plambda" || nm == "preset_shape_pl_h1" {
            cases.push(Case { name: nm.into(), t_max: 13.0, r_coll: 1e-3,
                              ics: sample(&ch, 0, cx, cy, half, n) });
        }
    }
    // The user's saved slices, at **their own** settings -- horizon 50, their own `r_coll`, and
    // the `r_esc` the config was saved with. `config_basin` is the control: in basin mode the
    // colour IS the outcome, so the freezing artefact has no exposure there.
    for (nm, (ch, cx, cy, half), rc) in [
        ("config_basin", Chart::config_basin(), 0.02),
        ("config_stability", Chart::config_stability(), 0.005),
    ] {
        cases.push(Case { name: nm.into(), t_max: 50.0, r_coll: rc,
                          ics: sample(&ch, 0, cx, cy, half, n) });
    }

    println!(
        "{n}x{n} per case, n_sync={n_sync}, eta={eta:e}, escape_confirm=OFF, \
         escape_all_bodies=ON\nr_esc ladder {esc_ladder:?} (fractions of R), strides {strides:?}\n"
    );

    // ---- 0. units --------------------------------------------------------------------------
    //
    // `r_esc` is canonical: a fraction of R, fixed at t = 0. The GLSL's 5 and 12 are absolute
    // lengths in the latent decode's own units, so whether the literal transfers is a question
    // about R on that chart -- and the latent decode normalises M = 1 with I an algebraic
    // identity equal to 1, which predicts R = sqrt(I/M) = 1 exactly. Measured, not assumed.
    println!("=== 0. canonical units: R = sqrt(I/M) at t = 0 ===");
    println!("{:>22}  {:>12} {:>12} {:>12}   r_esc=5 means", "case", "R p1", "R median", "R p99");
    for c in &cases {
        let mut r: Vec<f64> = c.ics.iter()
            .map(|(s, m)| energy::hyperradius(&s.r, m))
            .filter(|x| x.is_finite())
            .collect();
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (lo, md, hi) = (q(&mut r, 0.01), q(&mut r, 0.5), q(&mut r, 0.99));
        println!("{:>22}  {lo:>12.6} {md:>12.6} {hi:>12.6}   {:.4} in absolute units",
                 c.name, 5.0 * md);
    }
    println!(
        "\n  A latent decode has M = 1 and I = 1 identically, so R = 1 and the GLSL's absolute\n\
         `rEsc = 5` IS `5 R` there -- the literal transfers on that family with no conversion.\n\
         Burrau's R = 2.2361, where an absolute 5 would mean 2.24 R instead: the same defect a\n\
         shared `half` default had across two coordinate systems.\n"
    );

    // ---- 1 + 2. persistence and escape fraction ---------------------------------------------
    let windows: [usize; 4] = [1, 2, 4, 8];
    println!("=== 1+2. transient population and escape fraction ===");
    for c in &cases {
        let dt_sync = c.t_max / n_sync as f64;
        println!("\n--- {} --- {} trajectories, t_max={}, r_coll={:e}, sync interval {dt_sync:.5}",
                 c.name, c.ics.len(), c.t_max, c.r_coll);
        println!("{:>7} {:>7} {:>6} {:>9} {:>9}   {:>8} {:>8} {:>8} {:>8}",
                 "r_esc", "bodies", "stride", "escape", "collision", "+1", "+2", "+4", "+8");
        // **Two changes, two axes.** The distance gate and the all-three-bodies test are
        // separate divergences from the numpy form, and `(0, one)` is the reference this port
        // has shipped. Sweeping them together would score their sum -- and the sum is not
        // small: the body arm alone moves near-field's escape fraction from 0 to 1.
        for &(r_esc, all) in &[(0.0f64, false), (5.0, true)] {
            for &ev in &strides {
                let rows: Vec<(u8, Option<usize>, [Option<bool>; 4])> = c.ics.par_iter().map(|(s0, m)| {
                    // One unstopped run: both arms accumulate and `escape_flags` carries the
                    // whole boundary history at a single discretisation.
                    let o = integrate_az_opts(*s0, m, c.t_max, n_sync, eta,
                                              200_000, &opts(c.r_coll, r_esc, all, ev, false));
                    let esc = o.events.escape.map(|(_, te)| (te / dt_sync).floor() as usize);
                    let mut w = [None; 4];
                    if let Some(k0) = esc {
                        for (i, ww) in windows.iter().enumerate() {
                            w[i] = o.escape_flags.get(k0 + ww).copied();
                        }
                    }
                    let cls = if esc.is_some() { 1 } else { 0 }
                        | if o.events.collision.is_some() { 2 } else { 0 };
                    (cls, esc, w)
                }).collect();

                let nesc = rows.iter().filter(|r| r.0 & 1 != 0).count();
                let ncol = rows.iter().filter(|r| r.0 & 2 != 0).count();
                let mut cells = String::new();
                for i in 0..4 {
                    // A trajectory with no boundary that far ahead is EXCLUDED from the
                    // denominator, not scored either way: never given a chance to re-bind is
                    // not evidence that it did not.
                    let t = rows.iter().filter(|r| r.2[i].is_some()).count();
                    let s = rows.iter().filter(|r| r.2[i] == Some(true)).count();
                    if t == 0 {
                        cells.push_str(&format!("{:>9}", "n/a"));
                    } else {
                        cells.push_str(&format!("{:>9.3}", s as f64 / t as f64));
                    }
                }
                println!("{r_esc:>7.1} {:>7} {ev:>6} {:>9.4} {:>9.4}  {cells}",
                         if all { "all" } else { "one" },
                         nesc as f64 / rows.len() as f64,
                         ncol as f64 / rows.len() as f64);
            }
        }
    }

    // ---- 4. r_esc sensitivity ----------------------------------------------------------------
    //
    // Terminal-outcome fractions under `stop_on_event`, which is production. A new constant is
    // reported as a curve, never picked to make a picture look right.
    println!("\n=== 4. r_esc sensitivity: terminal outcome fractions (stop_on_event, stride 0) ===");
    for c in &cases {
        println!("\n--- {} --- t_max={}, r_coll={:e}", c.name, c.t_max, c.r_coll);
        println!("{:>7} {:>7}  {:>8} {:>9} {:>8} {:>8} {:>8}   {:>11}",
                 "r_esc", "bodies", "escape", "collision", "bounded", "running", "failed",
                 "t_end p50");
        for &all in &[true] {
        for &r_esc in &esc_ladder {
            let rows: Vec<(u8, f64)> = c.ics.par_iter().map(|(s0, m)| {
                let o = integrate_az_opts(*s0, m, c.t_max, n_sync, eta,
                                          200_000, &opts(c.r_coll, r_esc, all, 0, true));
                let oc = prin_rs::outcome::classify(&o.events, &o.state, m,
                                                    o.finite, o.budget_exhausted);
                (oc.state as u8, o.t_end)
            }).collect();
            let f = |k: u8| rows.iter().filter(|r| r.0 == k).count() as f64 / rows.len() as f64;
            let mut te: Vec<f64> = rows.iter().map(|r| r.1).filter(|x| x.is_finite()).collect();
            println!("{r_esc:>7.1} {:>7}  {:>8.4} {:>9.4} {:>8.4} {:>8.4} {:>8.4}   {:>11.5}",
                     if all { "all" } else { "one" },
                     f(0), f(2), f(1), f(3), f(4) + f(5), q(&mut te, 0.5));
        }
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         [1] is the test the fix stands or falls on. Ungated, the escapes that only a fine\n\
         stride sees do not persist -- 0 of 895 in `deep interior`. If the gated rows do not\n\
         read near 1.000 at +1 through +8, the distance gate is NOT the fix and no render\n\
         should be quoted from it.\n\n\
         [2] the escape fraction is expected to FALL where transients were being latched\n\
         (`deep interior`) and to be essentially UNCHANGED where escape genuinely terminates\n\
         (`preset_plambda`). A guard that cuts both is cutting too much, and that arm is what\n\
         separates a guard from a refusal.\n\n\
         [4] `r_esc` is a new constant. The curve is the argument for whatever value ships, and\n\
         a plateau in it is what says the answer is not sensitive to the choice.\n\n\
         The `bodies` column is the SECOND divergence from the numpy form and it is not free.\n\
         `one` is what this port has shipped; `all` is what the GLSL does. Read the two rows at\n\
         `r_esc = 0` first -- their difference is the body arm on its own, with no gate."
    );
}
