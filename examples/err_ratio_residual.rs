//! **What `burrau_nu_k`'s `error_ratio` residual actually is.** Diagnosis only: nothing here
//! changes a production setting, and it writes to whatever root it is given (default `results`).
//!
//! Run: `cargo run --release --example err_ratio_residual [case] [res] [root]`
//!
//! ## The premise this was built on was wrong, twice
//!
//! **(1) "The p99 population is filtered."** `examples/step_budget.rs` quantiles
//! `.filter(|x| x.is_finite())`, and `stats::max_dev` returns `+inf` for a non-finite copy, so the
//! truncated footprints looked excluded — a different statistic in every row. Measured: **2304 of
//! 2304 finite at every rung, 0.0% dropped.** Nothing was ever filtered, because a *starved*
//! footprint does not go non-finite: every copy stops at the same early point and agrees perfectly,
//! which is the standing *a starved footprint reads exactly 1.0000* finding.
//!
//! **(2) "It is a tail."** Fixing the population moved the tail barely at all and moved the
//! **median 30x**: `err p50` runs `1.082e0 -> 3.189e1` across the budget ladder while `err p99`
//! runs `1.127e2 -> 1.244e2`. At 30_000, 22.6% of the frame was truncated and reading ~1.0, which
//! dragged the median under it. So this is not a few pixels that are not data — **more than half
//! the chart sits at `error_ratio` around 32**, and the p99 is the top of a broad bulk rather than
//! an outlier. `refine_threshold` is 10.
//!
//! ## What this harness asks, in the order that can dissolve the question
//!
//! 1. **`error_ratio` against `error_ratio_coarse`.** The bulk is above the flag threshold, so
//!    production's repair pass already re-integrated it at `eta/4` up to three times, i.e. to
//!    `eta/64`. `error_ratio_coarse` preserves the pre-repair value. If the two agree to the
//!    digit, **step size moved nothing**, and *insensitivity to step size is the tell*: an `eta`
//!    ladder is then guaranteed flat and must not be built.
//! 2. **Which arm is large.** `error_ratio = sigma_e_t / sigma_e_0`, both dumped. A large ratio
//!    from a *small denominator* is the standing cell-width artefact (`sigma_E(0)` is proportional
//!    to the jitter and so to cell width); a large ratio from a *large numerator* is energies
//!    genuinely spreading.
//! 3. **`n_outcome_disagree`.** Copies straddling `r_coll` terminate differently, so their
//!    energies differ **for a physical reason** and a large ratio is the instrument reporting.
//!    Conditioning the ratio on it is what separates that from a numerical fault.
//!
//! The comparison target is passed alongside so every number has a control: the ratio's
//! **converged** value is exactly 1.0, and a chart reading 1.0 is what says the instrument works.

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::stats;
use prin_rs::grid::Chart;
use prin_rs::physics::energy;
use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid;
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut [f64], p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// `(name, chart, cx, cy, half)`, taken from `grid::gallery_cases` rather than written out here.
///
/// **The window is the gallery's own or it is a different chart.** `half = 0.45` on
/// `burrau_nu_k` is a `(nu, K)` chart coordinate; `0.05` on `near-field` is a body position. Two
/// numbers in two coordinate systems under one name is how `preset_shape` shipped at a 3x crop,
/// so neither is typed in.
fn cases() -> Vec<(String, Chart, f64, f64, f64)> {
    let mut v: Vec<(String, Chart, f64, f64, f64)> = grid::gallery_cases()
        .into_iter()
        .filter(|(n, ..)| *n == "burrau_nu_k")
        .map(|(n, c, a, b, h)| (n.to_string(), c, a, b, h))
        .collect();
    // The control: a Burrau region whose whole budget ladder reads `err p99 = 1.000`, so a
    // decomposition that cannot separate it from the subject is not a decomposition.
    let r = grid::region("near-field", 2, 2, 0.05).unwrap();
    v.push(("near-field".into(), Chart::BodyPlane, r.cx, r.cy, 0.05));
    v
}

fn main() {
    let want: String = arg(1, "burrau_nu_k".into());
    let res: usize = arg(2, 48);
    let _root: String = arg(3, "results".into());
    // **Argument four is `stop_on_event`, and it is an arm rather than a setting.** The standing
    // result is that a trajectory stopped by an event is parked **at** a close approach, where the
    // Cartesian energy is a cancellation of two enormous terms — measured once already at
    // 27000x on `far` under Heggie. `error_ratio` is built from exactly that readout
    // (`pixel.rs:829`), so if the residual is the readout and not the integration, turning
    // termination off must collapse it while leaving `energy_drift_max` where it is.
    let stop: bool = arg(4, 1) != 0;
    let ens = EnsembleCfg { stop_on_event: stop, ..EnsembleCfg::production() };

    println!("err_ratio_residual: {res}^2 per case, f64.");
    println!("  config: {}", ens.provenance());
    println!(
        "  refine_flagged={} threshold={} eta_factor={} max_passes={}",
        ens.refine_flagged, ens.refine_threshold, ens.refine_eta_factor, ens.refine_max_passes
    );
    println!();

    for (name, chart, cx, cy, half) in cases() {
        if want != "all" && want != name {
            continue;
        }
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &ens)).collect();

        println!("=== {name} ({}, half {half}) ===", chart.name());

        // **Why this chart and not the control.** The parked-readout mechanism can only bite
        // where trajectories actually terminate before the horizon. `near-field` at `t = 13` is
        // 97.8% *Bounded* and sits at `t_end = t_max`, which is the standing result that its
        // escape arm is silent at this horizon — there is no parked state to read.
        let term = px.iter().filter(|p| p.t_end < ens.t_max - 1e-9).count();
        println!(
            "  terminated before the horizon: {term} of {} ({:.1}%)",
            px.len(),
            100.0 * term as f64 / px.len() as f64
        );

        // --- 1. the repair pass ---
        let refined = px.iter().filter(|p| p.refined).count();
        let mut r_now: Vec<f64> = px.iter().map(|p| p.error_ratio).collect();
        let mut r_pre: Vec<f64> = px.iter().map(|p| p.error_ratio_coarse).collect();
        let n_fin = r_now.iter().filter(|x| x.is_finite()).count();
        println!(
            "  {} footprints, {n_fin} with a finite ratio, {refined} re-integrated ({:.1}%)",
            px.len(),
            100.0 * refined as f64 / px.len() as f64
        );
        for p in [0.50f64, 0.90, 0.99] {
            println!(
                "    ratio p{:<3.0}  after repair {:>10.4e}   before {:>10.4e}   x {:>7.4}",
                p * 100.0,
                q(&mut r_now.clone(), p),
                q(&mut r_pre.clone(), p),
                q(&mut r_now.clone(), p) / q(&mut r_pre.clone(), p)
            );
        }
        // Paired, not quantile-against-quantile: two ladders can agree in distribution while every
        // pixel moved. *Never conclude "no effect" from an aggregate without the per-pixel
        // distribution* — four sites on this project already.
        let moved = px
            .iter()
            .filter(|p| p.error_ratio.is_finite() && p.error_ratio_coarse.is_finite())
            .filter(|p| (p.error_ratio - p.error_ratio_coarse).abs() > 1e-12 * p.error_ratio.abs())
            .count();
        let mut ch: Vec<f64> = px
            .iter()
            .filter(|p| p.error_ratio.is_finite() && p.error_ratio_coarse.is_finite())
            .map(|p| (p.error_ratio / p.error_ratio_coarse).ln().abs())
            .collect();
        println!(
            "    paired: {moved} of {} pixels moved at all; |ln ratio| p50 {:.3e} p99 {:.3e}",
            px.len(),
            q(&mut ch.clone(), 0.50),
            q(&mut ch, 0.99)
        );
        r_now.clear();
        r_pre.clear();

        // --- 2. which arm ---
        for (label, f) in [
            ("sigma_e_0", (|p: &PixelOut| p.sigma_e_0) as fn(&PixelOut) -> f64),
            ("sigma_e_t", |p: &PixelOut| p.sigma_e_t),
        ] {
            let mut v: Vec<f64> = px.iter().map(f).filter(|x| x.is_finite()).collect();
            println!(
                "    {label:<10} n {:>5}  p10 {:>10.3e}  p50 {:>10.3e}  p90 {:>10.3e}",
                v.len(),
                q(&mut v.clone(), 0.10),
                q(&mut v.clone(), 0.50),
                q(&mut v, 0.90)
            );
        }
        // `error_ratio` is NaN at `sigma_e_0 == 0`, and those pixels are silently absent from
        // every quantile above. Reported, because the population is the thing that misled here.
        let zero0 = px.iter().filter(|p| p.sigma_e_0 == 0.0).count();
        println!("    sigma_e_0 exactly zero on {zero0} pixels (ratio is 0/0 there)");

        // --- 3. the physical explanation ---
        let mut by_dis: [(usize, Vec<f64>); 2] = [(0, Vec::new()), (0, Vec::new())];
        for p in &px {
            let i = usize::from(p.n_outcome_disagree > 0);
            by_dis[i].0 += 1;
            if p.error_ratio.is_finite() {
                by_dis[i].1.push(p.error_ratio);
            }
        }
        for (i, label) in ["copies AGREE on outcome", "copies DISAGREE on outcome"].iter().enumerate()
        {
            let (n, v) = &mut by_dis[i];
            println!(
                "    {label:<28} {n:>5} px   ratio p50 {:>10.3e}  p99 {:>10.3e}",
                q(&mut v.clone(), 0.50),
                q(v, 0.99)
            );
        }

        // --- 4. THE MASS SEAM, evaluated at `t = 0` where the answer is known ---
        //
        // `pixel.rs:824-831` builds the two energy vectors the ratio is taken over:
        //
        // ```
        // let e0 = copies.iter().map(|c| energy(&c.s.r, &c.s.v, &c.m, 0));   //  own masses
        // let et = outs.iter().map(|o| energy(&o.state.r, &o.state.v, &m, 0)); //  NOMINAL masses
        // ```
        //
        // `m` is `copies[0].m`. Eleven lines above, the comment that introduces it states the
        // invariant this violates in as many words: *"per-copy masses are used where a copy is
        // integrated **or its energy taken**, because on a mass chart two jittered copies are two
        // different systems."* On `BodyPlane` every copy shares the nominal's masses, so the two
        // spellings are the same expression and the seam is invisible — the equal-mass blind spot
        // this project already has on record for the Jacobi round-trip.
        //
        // The arm below needs **no integration and no chaos**: it takes the same `t = 0` states
        // and evaluates their energies both ways. Under exact dynamics `sigma_e(t) == sigma_e(0)`,
        // so a ratio built from a mismatched pair reads its mismatch at `t = 0` too. A chart that
        // does not vary its masses must read exactly 1.0, which is what makes the number mean
        // something on the chart that does not.
        let ens0 = &ens;
        let mut mass_spread: Vec<f64> = Vec::new();
        let mut ratio0: Vec<f64> = Vec::new();
        // **The scale the ratio is a deviation OF.** `energy_drift_max` is relative and
        // `sigma_e_*` are absolute, so the two columns are only comparable through `|E|`: a
        // relative drift of `4e-9` on `|E| ~ 1e8` is an absolute move of `0.4`, which is the
        // order `sigma_e_t` reads. Without this column the drift and the ratio look like they
        // contradict each other, and one of them looks wrong.
        let mut escale: Vec<f64> = Vec::new();
        for k in 0..sl.npix() {
            let cs = jitter::copies_in_space::<f64>(
                &sl, k, ens0.n_extra, ens0.jitter_frac, ens0.seed, ens0.jitter_scheme,
                ens0.decode_path, ens0.sample_space,
            );
            let m_nom = cs[0].m;
            mass_spread.push(
                cs.iter()
                    .map(|c| (0..3).map(|i| (c.m[i] - m_nom[i]).abs()).fold(0.0, f64::max))
                    .fold(0.0, f64::max),
            );
            let own: Vec<f64> =
                cs.iter().map(|c| energy::energy(&c.s.r, &c.s.v, &c.m, 0.0)).collect();
            let nom: Vec<f64> =
                cs.iter().map(|c| energy::energy(&c.s.r, &c.s.v, &m_nom, 0.0)).collect();
            let (r, _, _) = stats::error_ratio(&own, &nom);
            ratio0.push(r);
            escale.push(own[0].abs());
        }
        println!(
            "    mass spread across copies  p50 {:>10.3e}  max {:>10.3e}",
            q(&mut mass_spread.clone(), 0.50),
            mass_spread.iter().cloned().fold(0.0, f64::max)
        );
        println!(
            "    |E| at t=0                 p50 {:>10.3e}  p99 {:>10.3e}   (drift is RELATIVE)",
            q(&mut escale.clone(), 0.50),
            q(&mut escale, 0.99)
        );
        println!(
            "    t=0 ratio, nominal masses over own  p50 {:>10.4e}  p90 {:>10.4e}  p99 {:>10.4e}",
            q(&mut ratio0.clone(), 0.50),
            q(&mut ratio0.clone(), 0.90),
            q(&mut ratio0, 0.99)
        );
        println!();
    }
}
