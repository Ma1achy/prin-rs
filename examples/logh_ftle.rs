//! **Prediction 3: does each occupant's drift field track a sensitivity measure?**
//!
//! `examples/independent_check.rs` established the standing pair, `config_stability` at 192^2,
//! Spearman with a half-frame shifted control on every row:
//!
//! ```text
//!         FTLE vs AZ drift        -0.0820      shifted -0.1022
//!     FTLE vs LEAPFROG drift      +0.3048      shifted +0.0240
//! ```
//!
//! # Two corrections before that pair can be used as a target
//!
//! **AZ's `-0.0820` is a null, not a negative correlation.** Its own shifted control reads
//! `-0.1022` — *larger in magnitude*. A correlation that its own displaced control matches is
//! about the two marginals, not about where the fields are, so the honest statement is "AZ's
//! drift does not track FTLE". Reading it as a negative correlation would be reading a control
//! as a signal.
//!
//! **And the row that works has the same integrator on both sides.** `src/physics/ftle.rs:26`:
//! the FTLE sits on the plain leapfrog. So `FTLE vs LEAPFROG drift` correlates two quantities
//! produced by one stepper, and `FTLE vs AZ drift` does not. That asymmetry is load-bearing for
//! a logH comparison, because `logh_lf` and `plain_lf` share a *stepper* with the FTLE while
//! `logh_rk4`, `plain_rk4`, `az` and `heggie` do not.
//!
//! **This harness does not fix that.** Fixing it means an FTLE per occupant, which is
//! propagating a tangent vector through a regularised chart, and is a separate build. What it
//! does instead is make the confound legible: one shared FTLE field, every occupant against it,
//! the shifted control on every row, and the occupants **grouped by whether they share the
//! FTLE's stepper**. If the KDK arms come out high and the RK4 arms low regardless of
//! regularisation, the table is measuring the shared stepper and says so on its face.
//!
//! `frac_ok` is printed because the FTLE's own march is unregularised and expected to fail on
//! close encounters — *a correlation over a mostly-failed field is a correlation with the
//! failure pattern*.
//!
//! Args: `res root case`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::StepLimit;
use prin_rs::integrate::Integrator;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::ftle::{self, FtleOpts};

const DLO: f64 = 1e-12;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0;
        for &k in &idx[i..=j] {
            r[k] = avg;
        }
        i = j + 1;
    }
    r
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    if a.len() < 3 {
        return f64::NAN;
    }
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..a.len() {
        let (x, y) = (a[i] - ma, b[i] - mb);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    num / (da.sqrt() * db.sqrt()).max(f64::MIN_POSITIVE)
}

/// Spearman over a selected subset, with the half-frame shifted control alongside.
///
/// **The control is the load-bearing half.** A correlation that survives displacing one field by
/// half a frame is about the marginals rather than the alignment, and that is exactly how AZ's
/// `-0.0820` resolves into a null.
fn rho_pair(a: &[f64], b: &[f64], sel: &[usize], res: usize) -> (f64, f64) {
    let (x, y): (Vec<f64>, Vec<f64>) = sel.iter().map(|&i| (a[i], b[i])).unzip();
    let straight = pearson(&ranks(&x), &ranks(&y));
    let shifted: Vec<f64> = sel
        .iter()
        .map(|&i| {
            let (px, py) = (i % res, i / res);
            a[((py + res / 2) % res) * res + (px + res / 2) % res]
        })
        .collect();
    (straight, pearson(&ranks(&shifted), &ranks(&y)))
}

/// The occupants, with the flag that decides how the table must be read.
fn arms() -> Vec<(&'static str, Integrator, f64, bool, bool)> {
    // label, integrator, eta scale, predictive limit, shares the FTLE's stepper
    vec![
        ("az", Integrator::Az, 1.0, true, false),
        ("heggie", Integrator::Heggie, 1.0, true, false),
        ("logh_rk4", Integrator::LogHRk4, 1.0, true, false),
        ("plain_rk4", Integrator::PlainRk4, 1.0, false, false),
        ("logh_lf", Integrator::LogHLeapfrog, 0.25, true, true),
        ("plain_lf", Integrator::PlainLeapfrog, 0.25, false, true),
    ]
}

fn main() {
    let res: usize = arg(1, 192);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let case: String = std::env::args().nth(3).unwrap_or_else(|| "config_stability".into());
    let dir = format!("{root}/logh_arms");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half, t_max, _r_coll) = if case == "config_stability" {
        let (c, x, y, h) = Chart::config_stability();
        (c, x, y, h, 50.0f64, 0.005f64)
    } else {
        let sl = grid::region(&case.replace('_', " "), 4, 4, 0.05)
            .unwrap_or_else(|| panic!("unknown case {case}"));
        (sl.chart, sl.cx, sl.cy, sl.half, 13.0f64, 0.001f64)
    };
    let n_sync = (t_max / 0.4).round().max(4.0) as usize;
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let n = sl.npix();

    println!("{case} {res}^2, t_max = {t_max}, n_sync = {n_sync}\n");

    // --- the sensitivity field, computed ONCE and shared by every row ---
    let t0 = std::time::Instant::now();
    let fo = FtleOpts::default();
    let pert = ftle::unit_perturbation::<f64>(0);
    let ftl: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = sl.decode_pos(k);
            let st = grid::decode_state(&chart, 0, x, y);
            let o = ftle::integrate_full::<f64>(st.s, &st.m, t_max, 1e-4, &fo, &pert);
            if o.n_renorm > 0 {
                o.ftle
            } else {
                f64::NAN
            }
        })
        .collect();
    let ok: Vec<usize> = (0..n).filter(|&i| ftl[i].is_finite()).collect();
    println!(
        "FTLE pass {:.1}s -- frac_ok {:.4}. **The FTLE's own march is the unregularised\n\
         leapfrog and is expected to fail on close encounters**, so every row below is over the\n\
         pixels where it completed: a correlation over a mostly-failed field is a correlation\n\
         with the failure pattern.\n",
        t0.elapsed().as_secs_f64(),
        ok.len() as f64 / n as f64
    );

    println!("== SPEARMAN against ONE shared FTLE field, shifted control on every row ==\n");
    println!(
        "  **Read the `same stepper` column first.** The FTLE sits on the plain leapfrog\n  \
         (`src/physics/ftle.rs:26`), so the two `*_lf` rows share a stepper with it and the\n  \
         other four do not. If the correlation sorts by that column rather than by\n  \
         regularisation, this table is measuring the shared stepper.\n"
    );
    println!(
        "  {:>10} {:>12} {:>10} {:>10} {:>11} {:>11}",
        "arm", "same stepper", "spearman", "shifted", "drift p50", "secs"
    );

    let lg = |x: f64| if x.is_finite() && x > 0.0 { x.log10() } else { DLO.log10() };
    for (label, integ, eta_scale, limit, same) in arms() {
        let t = std::time::Instant::now();
        let cfg = EnsembleCfg::production().with_overrides(&[
            Override::TMax(t_max),
            Override::NSync(n_sync),
            // The DIAGNOSTIC pass: a trajectory stopped by an event is parked at a close
            // approach where the Cartesian energy is a cancellation of enormous terms, and that
            // flag alone produced five wrong conclusions in a row once.
            Override::RCollFrac(0.0),
            Override::StopOnEvent(false),
            Override::EscapeRule(EscapeRule::Closure(CLOSURE_TAU)),
            Override::ClosureK(1),
            Override::Integrator(integ),
            Override::Eta(EnsembleCfg::production().eta * eta_scale),
            Override::MaxSteps(400_000),
            Override::RefineFlagged(false),
            Override::StepLimit(if limit { StepLimit::Predictive } else { StepLimit::None }),
        ]);
        let px: Vec<PixelOut> =
            (0..n).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let d: Vec<f64> = px.iter().map(|p| lg(p.energy_drift_max)).collect();
        let (r, s) = rho_pair(&ftl, &d, &ok, res);
        let mut dr: Vec<f64> =
            px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
        dr.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p50 = if dr.is_empty() { f64::NAN } else { dr[dr.len() / 2] };
        println!(
            "  {label:>10} {:>12} {r:>10.4} {s:>10.4} {p50:>11.3e} {:>11.1}",
            if same { "YES" } else { "no" },
            t.elapsed().as_secs_f64()
        );
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **A row whose `shifted` control matches its `spearman` in magnitude is a null.** That\n\
         is how the standing `FTLE vs AZ drift = -0.0820` resolves: its control reads -0.1022,\n\
         larger, so AZ's drift does not track FTLE and the number is not a negative correlation.\n\
         The leapfrog's +0.3048 against a control of +0.0240 is a real signal.\n\n\
         **The confound is stated, not fixed.** One FTLE field is shared by six rows and it is\n\
         computed on the unregularised leapfrog, so two rows share its stepper. An FTLE per\n\
         occupant means propagating a tangent vector through a regularised chart and is a\n\
         separate build; until it exists, no number here separates 'this integrator's error is\n\
         shaped like the physics' from 'this integrator is the one the FTLE was computed with'."
    );
}
