//! Should `spread_event` latch?
//!
//! The numpy work observed ensemble spread **falling 6x between t=6 and t=8** —
//! diverge-then-reconverge — which is why the divergence accumulator is a latching field. The
//! same shape appears here: the playhead `spread_event` fires on more pixels at `t_max = 8`
//! than at `t_max = 13`.
//!
//! **But a discrete label has a failure mode a continuous divergence measure does not.** If two
//! pairs are near-equal in separation, copies can disagree about which is *tightest* without
//! their trajectories having diverged at all. A running max would latch that artefact
//! permanently and it would never clear.
//!
//! So: at the boundary where disagreement first occurs, how close were the two tightest pairs?

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;

const SIZE: usize = 32;

fn qs(v: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| x[(((x.len() - 1) as f64) * f).round() as usize];
    (q(0.0), q(0.1), q(0.5), q(0.9), q(1.0))
}

fn main() {
    let s = grid::region("near-field", SIZE, SIZE, 0.05).unwrap();
    let cfg = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let px: Vec<PixelOut> = (0..s.npix())
        .into_par_iter()
        .map(|i| evaluate::<f64>(&s, i, &cfg))
        .collect();

    // Pixels that fired at some boundary, split by whether they are still firing at the
    // playhead. The un-firing ones are the whole question: latching would keep them lit.
    let fired: Vec<&PixelOut> = px.iter().filter(|p| p.spread_event_max > 0.0).collect();
    let unfired: Vec<&PixelOut> = fired.iter().filter(|p| p.spread_event == 0.0).cloned().collect();
    let still: Vec<&PixelOut> = fired.iter().filter(|p| p.spread_event > 0.0).cloned().collect();

    println!("near-field {SIZE}x{SIZE}, t=13, E+1=8, eta=0.01, f64");
    println!();
    println!("pixels that ever disagree: {}/{}", fired.len(), px.len());
    println!("  still disagreeing at the playhead: {}", still.len());
    println!("  UN-FIRED - disagreed, then re-agreed: {}", unfired.len());
    println!();
    println!("tie ratio = (second-tightest / tightest) pair separation, minimised over the 8");
    println!("copies, at the boundary where the copies FIRST disagree. 1.0 is an exact tie.");
    println!();
    println!("{:>26}{:>10}{:>10}{:>10}{:>10}{:>10}{:>8}",
             "population", "min", "p10", "median", "p90", "max", "n");

    for (name, set) in [
        ("all that ever disagree", &fired),
        ("un-fired (re-agreed)", &unfired),
        ("still disagreeing", &still),
    ] {
        if set.is_empty() {
            println!("{name:>26}{:>58}", "(empty)");
            continue;
        }
        let v: Vec<f64> = set.iter().map(|p| p.tie_ratio_at_disagree).collect();
        let (a, b, c, d, e) = qs(&v);
        println!("{name:>26}{a:>10.4}{b:>10.4}{c:>10.4}{d:>10.4}{e:>10.4}{:>8}", set.len());
    }

    // A control: the same ratio over the whole grid, sampled at the final boundary, so
    // "near 1" has something to be near-1 relative to.
    println!();
    let near_tie = |t: f64| fired.iter().filter(|p| p.tie_ratio_at_disagree < t).count();
    println!("of the {} pixels that ever disagree, how many were at a near-tie:", fired.len());
    for t in [1.01f64, 1.05, 1.1, 1.25, 1.5, 2.0] {
        println!("  tie ratio < {t:<5} : {:>4}", near_tie(t));
    }
    if !unfired.is_empty() {
        println!();
        let nt = unfired.iter().filter(|p| p.tie_ratio_at_disagree < 1.1).count();
        println!("of the {} UN-FIRED pixels, {nt} were at a near-tie (< 1.1) when they first",
                 unfired.len());
        println!("disagreed. That is the number the latching decision turns on.");
    }

    println!();
    println!("=== does persistence separate them where the tie ratio does not? ===");
    println!("{:>26}{:>16}{:>16}{:>16}{:>16}",
             "population", "n_disagree med", "n_disagree max", "run med", "run max");
    for (name, set) in [
        ("un-fired (re-agreed)", &unfired),
        ("still disagreeing", &still),
    ] {
        if set.is_empty() {
            continue;
        }
        let nd: Vec<f64> = set.iter().map(|p| p.n_disagree as f64).collect();
        let rn: Vec<f64> = set.iter().map(|p| p.longest_disagree_run as f64).collect();
        let (_, _, ndm, _, ndx) = qs(&nd);
        let (_, _, rm, _, rx) = qs(&rn);
        println!("{name:>26}{ndm:>16.1}{ndx:>16.1}{rm:>16.1}{rx:>16.1}");
    }
    println!();
    for k in [1u16, 2, 3, 4] {
        let a = unfired.iter().filter(|p| p.longest_disagree_run >= k).count();
        let b = still.iter().filter(|p| p.longest_disagree_run >= k).count();
        println!("  a latch requiring a run of >= {k}: keeps {b}/{} genuine, {a}/{} artefact",
                 still.len(), unfired.len());
    }

    println!();
    println!("=== the guarded latch, as implemented ===");
    let la = unfired.iter().filter(|p| p.spread_event_latched > 0.0).count();
    let lb = still.iter().filter(|p| p.spread_event_latched > 0.0).count();
    let lm = px.iter().filter(|p| p.spread_event_latched > 0.0).count();
    println!("spread_event_latched (run >= {}, joined with the playhead value):",
             prin_rs::ensemble::pixel::LATCH_RUN);
    println!("  lit on {lb}/{} genuine, {la}/{} artefact, {lm}/{} of the whole grid",
             still.len(), unfired.len(), px.len());
    println!("compare: an unguarded running max would be lit on {}/{}",
             fired.len(), px.len());

    println!();
    println!("=== monotonicity, measured WITHIN one run ===");
    println!("Varying t_max at fixed n_sync would not answer this: the sync grid changes with");
    println!("the horizon, so dtau changes and the rows are different discretisations, not one");
    println!("run truncated at different playheads. (Measured, for the record: sweeping t_max");
    println!("gives an unguarded max of 109, 297, 110, 488, 165 pixels at t = 4, 6, 8, 10, 13 -");
    println!("a running max cannot fall, so those rows are not nested.)");
    println!();
    println!("So: one run to t = 13 with n_sync = 32, evaluated at every boundary as if it were");
    println!("the playhead.");
    println!();
    println!("{:>5}{:>9}{:>14}{:>14}{:>14}", "k", "t", "playhead", "latched", "unguarded max");
    let curves = per_boundary_curves(&s);
    for (k, (t, a, b, c)) in curves.iter().enumerate() {
        if k % 4 == 3 || k + 1 == curves.len() {
            println!("{:>5}{t:>9.4}{a:>14}{b:>14}{c:>14}", k);
        }
    }

    println!();
    println!("Reading: near 1 means the copies disagreed about which pair is tightest while");
    println!("their trajectories had barely separated - a near-degeneracy, not decoherence,");
    println!("and latching it would light a pixel permanently for a labelling artefact.");
    println!("Well above 1 means the disagreement is genuine and latching is correct.");
}

/// Nonzero-pixel counts at every sync boundary of a single `t = 13`, `n_sync = 32` run:
/// the playhead value, the guarded latch, and the unguarded running max.
fn per_boundary_curves(s: &prin_rs::grid::Slice) -> Vec<(f64, usize, usize, usize)> {
    use prin_rs::ensemble::{jitter, stats};
    use prin_rs::integrate::az::{self, AzOpts};
    use prin_rs::outcome;
    use prin_rs::physics::burrau;

    const N_SYNC: usize = 32;
    let (t_max, m) = (13.0f64, burrau::masses::<f64>());
    let latch = prin_rs::ensemble::pixel::LATCH_RUN as usize;

    let per_pixel: Vec<Vec<f64>> = (0..s.npix())
        .into_par_iter()
        .map(|i| {
            let copies = jitter::copies::<f64>(s, i, 7, 0.5, 0);
            let outs: Vec<_> = copies
                .iter()
                .map(|c| {
                    az::integrate_az_opts(
                        *c, &m, t_max, N_SYNC, 0.01, 30_000,
                        &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
                    )
                })
                .collect();
            let term: Vec<u8> = outs
                .iter()
                .map(|o| outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted).pack())
                .collect();
            (0..N_SYNC)
                .map(|k| {
                    let cls: Vec<u8> = outs
                        .iter()
                        .zip(term.iter())
                        .map(|(o, &c)| stats::event_class_at(&o.tight, c, k))
                        .collect();
                    stats::spread_event::<f64>(&cls)
                })
                .collect()
        })
        .collect();

    (0..N_SYNC)
        .map(|k| {
            let t = (k + 1) as f64 * t_max / N_SYNC as f64;
            let mut a = 0;
            let mut b = 0;
            let mut c = 0;
            for row in &per_pixel {
                if row[k] > 0.0 {
                    a += 1;
                }
                if row[..=k].iter().any(|&x| x > 0.0) {
                    c += 1;
                }
                // The latch as implemented: any run of >= LATCH_RUN up to k, or the value at k.
                let mut run = 0usize;
                let mut lit = row[k] > 0.0;
                for &x in &row[..=k] {
                    run = if x > 0.0 { run + 1 } else { 0 };
                    if run >= latch {
                        lit = true;
                    }
                }
                if lit {
                    b += 1;
                }
            }
            (t, a, b, c)
        })
        .collect()
}
