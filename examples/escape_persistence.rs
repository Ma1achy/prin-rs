//! Does a fired escape PERSIST, and how much of the precedence bug was real?
//!
//! Two questions, both of which decide what the escape-test stride *means*.
//!
//! # 1. Persistence, and why "it converged" is not established
//!
//! `escape_candidate` is relative energy `> 0` and receding. During a close encounter a pair's
//! two-body energy can be **transiently positive and then re-bind**, so a test run every RK4
//! step catches transients that a boundary-only test mostly misses. If those transients are
//! being latched as terminal escapes, a finer stride does not resolve a quantity -- it invents
//! events.
//!
//! The escape fraction in `deep interior` ran **0.0945 -> 0.2153 -> 0.4423 -> 0.5494** across
//! strides `0, 32, 4, 1`. That was called converging. It is not: the step from 4 to 1 is
//! **+24% relative**, the largest relative move in the sequence, and an absorbing state
//! sampled ever more finely should settle. **Treat convergence as unproven** until the
//! escapes are shown to stick.
//!
//! The test: take the trajectories that escape only under the fine stride, integrate **past**
//! detection, and ask whether they are still unbound one sync interval later. All persist and
//! no guard is needed. A fraction re-bind and a **persistence guard is required** -- the same
//! shape as `spread_event_latched`'s `LATCH_RUN`, which exists for exactly this reason on a
//! different field.
//!
//! # 2. How much the precedence bug was actually costing
//!
//! `classify` used to rank collision above escape unconditionally, discarding both times.
//! Under `stop_on_event` that is nearly unobservable, because the loop breaks on the first
//! *detected* event and only one is ever recorded -- measured, the repair moved **one footprint
//! of 5440**. On the reference path (`stop_on_event = false`) both arms accumulate over the
//! whole run, and that is where the mis-ordering is visible. Counted here.
//!
//! # Writes
//!
//! stdout only.

use prin_rs::integrate::az::{integrate_az_opts, AzOpts};
use prin_rs::grid::{self, Chart};
use prin_rs::physics::Cart;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

struct Case {
    name: &'static str,
    /// `(u, v) -> latent`, so a whole slice can be swept from one closure.
    ics: Vec<(Cart<f64>, [f64; 3])>,
}

/// Sample a chart over `[cx +/- half] x [cy +/- half]` at cell centres.
///
/// **Through `grid::decode_state`, not a hand-rolled slice.** The `shape_pl` basis literal
/// already appeared three times on this project and was wrong in all three; a second definition
/// of what a region *is* would be a fourth place for it to disagree.
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

fn main() {
    let n: usize = arg(1, 48);
    let t_max: f64 = arg(2, 13.0);
    let n_sync: usize = arg(3, 32);
    let eta = 1e-2f64;
    let dt_sync = t_max / n_sync as f64;

    let opts = |ev: usize, stop: bool| AzOpts::<f64> {
        step_limit: prin_rs::integrate::az::StepLimit::None,
        step_blend: prin_rs::integrate::az::StepBlend::Min,
        blend_p: 4.0,
        step_limit_f: 0.0,
        dtau_mode: prin_rs::integrate::az::DtauMode::default(),
        clamp_final_step: true,
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: 1e-3,
        stop_on_event: stop,
        escape_every: ev,
        // **Deliberately OFF here.** This example is the diagnostic that measures whether the
        // guard is needed; running it with the guard on would measure the guard.
        escape_confirm: false,
        // The numpy reference's ungated escape test: every result in this diagnostic
        // predates the distance gate and is quoted against that form.
        escape_rule: prin_rs::outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_escape: stop,
        keep_boundary_shapes: false,
        keep_drift_hist: false,
    };

    let mut cases: Vec<Case> = Vec::new();
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            cases.push(Case {
                name: region,
                ics: sample(&Chart::BodyPlane, body, cx, cy, 0.05, n),
            });
        }
    }
    // `plambda`: momentum coordinates only, configuration fixed. Where escape terminates, and
    // therefore the only family the cadence can reach.
    for (nm, ch, cx, cy, half) in grid::gallery_cases() {
        if nm == "preset_plambda" {
            cases.push(Case { name: nm, ics: sample(&ch, 0, cx, cy, half, n) });
        }
    }

    println!(
        "{n}x{n} per case, t={t_max}, n_sync={n_sync} (interval {dt_sync:.6}), eta={eta:e}, \
         r_coll_frac=1e-3\n"
    );

    for c in &cases {
        println!("=== {} === {} trajectories", c.name, c.ics.len());

        // ---- 2. the precedence population, on the reference path ----
        //
        // `stop_on_event = false` lets BOTH arms accumulate over the whole run, which is the
        // only configuration in which the old fixed precedence is observable at all.
        let (mut both, mut esc_first) = (0usize, 0usize);
        let mut lead: Vec<f64> = Vec::new();
        for (s0, m) in &c.ics {
            let o = integrate_az_opts(*s0, m, t_max, n_sync, eta, 200_000, &opts(0, false));
            if let (Some((_, tc)), Some((_, te))) = (o.events.collision, o.events.escape) {
                both += 1;
                if te < tc {
                    esc_first += 1;
                    lead.push(tc - te);
                }
            }
        }
        let med = if lead.is_empty() {
            f64::NAN
        } else {
            prin_rs::stats::quantile(&mut lead.clone(), 0.5)
        };
        println!(
            "  [precedence] on the reference path (stop_on_event=false): {both} trajectories \
             fired BOTH arms;\n\
             \x20             {esc_first} of them escaped FIRST and were labelled `collision` by \
             the old fixed order\n\
             \x20             ({:.2}% of all, {:.1}% of those with both). Median lead: {med:.4} \
             ({:.2} sync intervals).",
            esc_first as f64 / c.ics.len() as f64 * 100.0,
            if both == 0 { f64::NAN } else { esc_first as f64 / both as f64 * 100.0 },
            med / dt_sync
        );

        // ---- 1. persistence of the escapes only the fine stride sees ----
        //
        // The population is the trajectories that escape at `escape_every = 1` and did NOT at
        // the reference cadence. For each, integrate past detection and ask whether it is still
        // an escape candidate a sync interval later. **This is the whole question**: a transient
        // that re-binds is not a terminal event, and latching it would invent one.
        // **One integration per IC, one discretisation, the whole history.** The first cut
        // re-ran to `t_e + w` with `n_sync` rescaled per window, which makes every window a
        // different discretisation -- and it produced a non-monotone curve (0.162, 0.219,
        // 0.011, 0.083, 0.335) that was the instrument, not the physics. `escape_flags` records
        // instantaneous candidacy at every boundary of a single run instead.
        //
        // The population is the trajectories that escape at `escape_every = 1` and did NOT at
        // the reference cadence -- the ones the finer stride invents or reveals.
        let windows: Vec<usize> = vec![1, 2, 3, 4, 8];
        let mut still: Vec<usize> = vec![0; windows.len()];
        let mut testable: Vec<usize> = vec![0; windows.len()];
        let (mut newly, mut ever_rebound) = (0usize, 0usize);
        let mut t_e_all: Vec<f64> = Vec::new();
        for (s0, m) in &c.ics {
            let coarse = integrate_az_opts(*s0, m, t_max, n_sync, eta, 200_000, &opts(0, true));
            let fine = integrate_az_opts(*s0, m, t_max, n_sync, eta, 200_000, &opts(1, true));
            let (Some((_, te)), None) = (fine.events.escape, coarse.events.escape) else {
                continue;
            };
            newly += 1;
            t_e_all.push(te);
            // The full history at the reference cadence, unstopped: one trajectory, one step
            // size, candidacy at every boundary.
            let hist = integrate_az_opts(*s0, m, t_max, n_sync, eta, 200_000, &opts(0, false));
            let k0 = (te / dt_sync).floor() as usize;
            if hist.escape_flags.iter().skip(k0).any(|f| !f) {
                ever_rebound += 1;
            }
            for (wi, w) in windows.iter().enumerate() {
                match hist.escape_flags.get(k0 + w) {
                    // No boundary that far ahead: the horizon, not a re-binding. Excluded from
                    // the denominator rather than scored either way -- a trajectory that was
                    // never given a chance to re-bind is not evidence that it did not.
                    None => {}
                    Some(f) => {
                        testable[wi] += 1;
                        if *f {
                            still[wi] += 1;
                        }
                    }
                }
            }
        }
        let mut tv = t_e_all.clone();
        let tmed = if tv.is_empty() { f64::NAN } else { prin_rs::stats::quantile(&mut tv, 0.5) };
        println!(
            "  [persistence] {newly} trajectories escape at stride 1 and NOT at the reference \
             cadence.\n\
             \x20              {ever_rebound} of them are NOT unbound at some later boundary \
             ({:.1}% flicker at least once).\n\
             \x20              median detection time {tmed:.4} ({:.2} sync intervals in).",
            if newly == 0 { f64::NAN } else { ever_rebound as f64 / newly as f64 * 100.0 },
            tmed / dt_sync
        );
        if newly > 0 {
            print!("               still unbound N boundaries later: ");
            for (wi, w) in windows.iter().enumerate() {
                let t = testable[wi];
                if t == 0 {
                    print!("+{w}: n/a  ");
                } else {
                    print!("+{w}: {:.3} (n {t})  ", still[wi] as f64 / t as f64);
                }
            }
            println!();
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         [persistence] decides what the stride MEANS. `escape_candidate` is relative energy > 0\n\
         and receding, and during a close encounter that can be transiently true. If a material\n\
         fraction RE-BIND, a finer stride is not resolving a quantity -- it is latching\n\
         transients as terminal events, and a persistence guard is required before the escape\n\
         fraction at stride 1 can be believed. If essentially all persist, the stride is a pure\n\
         resolution knob and 0.5494 stands.\n\n\
         [precedence] bounds the ordering bug. Under `stop_on_event` the loop breaks on the\n\
         first DETECTED event so only one is ever recorded, which is why repairing `classify`\n\
         moved one footprint of 5440 in production settings. The reference path accumulates\n\
         both, and that count is the size of the defect the repair addresses."
    );
}
