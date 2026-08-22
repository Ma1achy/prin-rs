//! What changes when `spread_event` is built on the event class instead of the terminal
//! outcome.
//!
//! BRIEF §4 defined it over the modal `(state, detail)`. That is terminal-grain: early in the
//! march nothing has terminated, so every copy agrees and the field reports maximum confidence
//! at exactly the playhead where least is known. The quantity it was meant to be is the
//! identity of the **currently tightest pair**, evaluated at every sync boundary and joined
//! with the terminal class where one exists.
//!
//! This prints both, side by side, and what `ensemble_spread` does as a result.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;

const SIZE: usize = 32;

fn render(t_max: f64) -> Vec<PixelOut> {
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let cfg = EnsembleCfg { t_max, ..Default::default() };
    (0..s.npix()).into_par_iter().map(|i| evaluate::<f64>(&s, i, &cfg)).collect()
}

fn stat(v: impl Iterator<Item = f64>) -> (f64, f64, f64) {
    let mut x: Vec<f64> = v.filter(|q| q.is_finite()).collect();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (x[0], x[x.len() / 2], x[x.len() - 1])
}

fn main() {
    println!("near-field {SIZE}x{SIZE}, E+1=8, eta=0.01, f64");
    println!();
    println!("{:>7}{:>34}{:>34}{:>26}", "t_max", "spread_event (event class)", "spread_event (terminal)", "nonzero pixels");
    println!("{:>7}{:>12}{:>11}{:>11}{:>12}{:>11}{:>11}{:>13}{:>13}",
             "", "min", "median", "max", "min", "median", "max", "event", "terminal");

    for t_max in [4.0f64, 8.0, 13.0, 20.0] {
        let px = render(t_max);
        let (a0, a1, a2) = stat(px.iter().map(|p| p.spread_event));
        let (b0, b1, b2) = stat(px.iter().map(|p| p.spread_event_terminal));
        let na = px.iter().filter(|p| p.spread_event > 0.0).count();
        let nb = px.iter().filter(|p| p.spread_event_terminal > 0.0).count();
        println!("{t_max:>7.1}{a0:>12.4}{a1:>11.4}{a2:>11.4}{b0:>12.4}{b1:>11.4}{b2:>11.4}{:>13}{:>13}",
                 format!("{na}/{}", px.len()), format!("{nb}/{}", px.len()));
    }

    println!();
    println!("=== at the project horizon, t = 13 ===");
    let px = render(13.0);
    let (_, ss_med, ss_max) = stat(px.iter().map(|p| p.spread_shape));
    let (_, se_med, se_max) = stat(px.iter().map(|p| p.spread_event));
    let (_, st_med, st_max) = stat(px.iter().map(|p| p.spread_event_terminal));
    let (_, es_med, es_max) = stat(px.iter().map(|p| p.ensemble_spread));

    // What ensemble_spread would have been under the old definition.
    let old_es: Vec<f64> = px.iter().map(|p| p.spread_shape.max(p.spread_event_terminal)).collect();
    let (_, oe_med, oe_max) = stat(old_es.iter().cloned());

    println!("{:>28}{:>14}{:>14}", "", "median", "max");
    println!("{:>28}{ss_med:>14.6}{ss_max:>14.6}", "spread_shape");
    println!("{:>28}{se_med:>14.6}{se_max:>14.6}", "spread_event (corrected)");
    println!("{:>28}{st_med:>14.6}{st_max:>14.6}", "spread_event (terminal)");
    println!("{:>28}{es_med:>14.6}{es_max:>14.6}", "ensemble_spread (corrected)");
    println!("{:>28}{oe_med:>14.6}{oe_max:>14.6}", "ensemble_spread (old)");

    let dominated = px.iter().filter(|p| p.spread_event > p.spread_shape).count();
    let dominated_old = px
        .iter()
        .filter(|p| p.spread_event_terminal > p.spread_shape)
        .count();
    println!();
    println!("pixels where spread_event exceeds spread_shape, so it is the one that sets");
    println!("ensemble_spread:  corrected {dominated}/{}   terminal {dominated_old}/{}",
             px.len(), px.len());

    let fired: Vec<f64> = px.iter().map(|p| p.t_spread_event).filter(|x| x.is_finite()).collect();
    let (t0, t1, t2) = stat(fired.iter().cloned());
    println!();
    println!("t_spread_event over the {} pixels that fire: min {t0:.4}  median {t1:.4}  max {t2:.4}",
             fired.len());
    println!("pixels whose copies never disagree before the horizon: {}/{}",
             px.len() - fired.len(), px.len());
    let (_, sm_med, sm_max) = stat(px.iter().map(|p| p.spread_event_max));
    println!();
    println!("running max over boundaries (monotone in the horizon): median {sm_med:.6}  max {sm_max:.6}");
    println!("pixels where it is nonzero: {}/{}",
             px.iter().filter(|p| p.spread_event_max > 0.0).count(), px.len());
    println!("The playhead value is a snapshot and can un-fire - the tightest pair fluctuates,");
    println!("so copies that disagreed at one boundary can agree again at the next. That is");
    println!("why the t_max=8 row above has MORE nonzero pixels than the t_max=13 row. The");
    println!("running max is monotone; both are dumped and the spec one is the default.");

    println!();
    println!("=== how much earlier, measured directly ===");
    first_disagreement();

    println!();
    println!("The terminal statistic cannot fire before something terminates; the event class");
    println!("is defined at every playhead. On this slice the two are strictly nested - every");
    println!("pixel the terminal statistic flags, the event class flags too, and 143 more.");
    println!();
    println!("But note WHERE the gain is. On the 22 pixels both flag, the lead time is exactly");
    println!("zero: they fire at the same boundary, because a collision is the tightest pair");
    println!("reaching r_coll and that usually settles the tightest-pair identity at the same");
    println!("boundary. The '~4 time units earlier' framing does not reproduce as a lead time");
    println!("here. It reproduces as COVERAGE - 165 pixels flagged against 22, a factor of");
    println!("7.5 - and as horizon-independence: at t_max = 8 it is 110 against 0.");
}

/// The claim the correction rests on: the event class fires earlier than the terminal class.
///
/// Measured directly rather than inferred. For each pixel, both statistics are evaluated at
/// every sync boundary and the first boundary at which each becomes nonzero is recorded. The
/// terminal class is "not yet terminated" for every copy until it terminates, so it cannot
/// fire before the first termination, by construction.
fn first_disagreement() {
    use prin_rs::ensemble::{jitter, stats};
    use prin_rs::integrate::az::{self, AzOpts};
    use prin_rs::outcome;
    use prin_rs::physics::burrau;

    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let (t_max, n_sync) = (13.0f64, 32usize);

    let rows: Vec<(f64, f64)> = (0..s.npix())
        .into_par_iter()
        .map(|i| {
            let copies = jitter::copies::<f64>(&s, i, 7, 0.5, 0);
            let outs: Vec<_> = copies
                .iter()
                .map(|c| {
                    az::integrate_az_opts(
                        *c, &m, t_max, n_sync, 0.01, 30_000,
                        // stop_on_event ON, matching the dumped field: a terminated copy's
                        // tight record ends there and the join with its terminal class takes
                        // over. With it off the join is never exercised and the comparison
                        // would be between pure tightest-pair and pure terminal.
                        &AzOpts { stop_on_event: true, r_coll_frac: 1e-3, ..Default::default() },
                    )
                })
                .collect();
            let term: Vec<u8> = outs
                .iter()
                .map(|o| outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted).pack())
                .collect();
            // Boundary index at which each copy terminated, if it did.
            let term_k: Vec<Option<usize>> = outs
                .iter()
                .map(|o| {
                    let te = o.events.collision.map(|(_, t)| t)
                        .or(o.events.escape.map(|(_, t)| t))?;
                    Some(((te / t_max * n_sync as f64).ceil() as usize).min(n_sync - 1))
                })
                .collect();

            let mut t_event = f64::NAN;
            let mut t_term = f64::NAN;
            for k in 0..n_sync {
                let t_now = (k + 1) as f64 * t_max / n_sync as f64;
                if t_event.is_nan() {
                    let cls: Vec<u8> = outs
                        .iter()
                        .zip(term.iter())
                        .map(|(o, &c)| stats::event_class_at(&o.tight, c, k))
                        .collect();
                    if stats::spread_event::<f64>(&cls) > 0.0 {
                        t_event = t_now;
                    }
                }
                if t_term.is_nan() {
                    // 0 = not yet terminated; otherwise the terminal class, offset.
                    let cls: Vec<u8> = term_k
                        .iter()
                        .zip(term.iter())
                        .map(|(tk, &c)| match tk {
                            Some(j) if *j <= k => 1 + c,
                            _ => 0,
                        })
                        .collect();
                    if stats::spread_event::<f64>(&cls) > 0.0 {
                        t_term = t_now;
                    }
                }
            }
            (t_event, t_term)
        })
        .collect();

    let both: Vec<(f64, f64)> = rows.iter().cloned().filter(|(a, b)| a.is_finite() && b.is_finite()).collect();
    let only_event = rows.iter().filter(|(a, b)| a.is_finite() && !b.is_finite()).count();
    let only_term = rows.iter().filter(|(a, b)| !a.is_finite() && b.is_finite()).count();
    println!("pixels where both statistics eventually fire: {}/{}", both.len(), rows.len());
    println!("  event class fires but terminal never does: {only_event}");
    println!("  terminal fires but event class never does: {only_term}");
    if !both.is_empty() {
        let mut d: Vec<f64> = both.iter().map(|(a, b)| b - a).collect();
        d.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let earlier = d.iter().filter(|x| **x > 0.0).count();
        println!("  event class fires earlier on {earlier} of {}", d.len());
        println!("  lead time (terminal - event): min {:.4}  median {:.4}  max {:.4}",
                 d[0], d[d.len() / 2], d[d.len() - 1]);
    }
}
