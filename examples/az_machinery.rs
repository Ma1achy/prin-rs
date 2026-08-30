//! **Which piece of Aarseth-Zare draws the wedges?** The re-registration count, or the LC branch.
//!
//! The argmax is ruled out four ways (depleted in the hot regions, transverse to the wedge edges,
//! `spread_shape`-invariant under hysteresis, and drift-invariant under it too). What has never
//! been touched is the machinery that runs at **every** boundary whether or not the reference
//! changes: the Cartesian -> regularised -> Cartesian round-trip, and the LC branch choice inside
//! it.
//!
//! # The confound that would eat this experiment, and the control for it
//!
//! `dt ~ eta * t_max / n_sync`, so raising `n_sync` at fixed `eta` also makes every step smaller.
//! "More boundaries" would then be inseparable from "finer stepping", and the result would be a
//! step-size result wearing a re-registration label. **`eta` is scaled with `n_sync`** to hold the
//! step size fixed, and `total_substeps` is printed as the check that it worked.
//!
//! A second confound, from this project's own record: the closure escape window is a **time**,
//! `closure_k * t_max / n_sync`. Changing `n_sync` at fixed `closure_k` silently changes the
//! escape criterion. `closure_k` is scaled too, which is why the ladder starts at the production
//! `n_sync` and doubles — below it the window cannot be held with an integer `k`.
//!
//! And a **deliberately confounded arm** is included, `n_sync` doubled with `eta` and `closure_k`
//! left alone, to show what the uncontrolled comparison would have claimed. Without it the
//! controlled rows are just numbers; with it they are a demonstration that the control mattered.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

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

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.015], [0.34, 0.06, 0.43], [0.72, 0.21, 0.33],
        [0.98, 0.55, 0.04], [0.99, 1.0, 0.64],
    ];
    let t = x.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    let mut o = [0u8; 3];
    for k in 0..3 {
        o[k] = (255.0 * (S[i][k] * (1.0 - f) + S[i + 1][k] * f)).clamp(0.0, 255.0) as u8;
    }
    o
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/az_machinery");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let n0 = (50.0f64 / WINDOW).round() as usize; // 125, the production cadence
    let base = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: n0,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", base.provenance());

    // (label, n_sync, eta, closure_k, lc_stable)
    let arms: Vec<(String, usize, f64, usize, bool)> = vec![
        ("baseline n=125".into(), n0, base.eta, 1, true),
        ("n=250 CONTROLLED".into(), n0 * 2, base.eta * 2.0, 2, true),
        ("n=500 CONTROLLED".into(), n0 * 4, base.eta * 4.0, 4, true),
        // The uncontrolled comparison, included so the controlled rows mean something.
        ("n=250 confounded".into(), n0 * 2, base.eta, 1, true),
        ("LC branch unconditioned".into(), n0, base.eta, 1, false),
    ];

    println!(
        "  {:>26} {:>7} {:>9} {:>4} {:>11} {:>11} {:>12} {:>12}",
        "arm", "n_sync", "eta", "lc", "steps p50", "drift p50", "hot lift", "chord p50"
    );

    let mut ref0: Option<(Vec<bool>, Vec<f64>)> = None;
    for (label, ns, eta, ck, lc) in arms {
        let cfg = EnsembleCfg {
            n_sync: ns,
            eta,
            closure_k: ck,
            lc_stable: lc,
            ..base
        };
        let px: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
            .collect();
        let n = px.len();

        let lg = |x: f64| if x.is_finite() && x > 0.0 { x.log10() } else { DLO.log10() };
        let d: Vec<f64> = px.iter().map(|p| lg(p.energy_drift_max)).collect();
        let mut dr: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
        let cut = q(&mut dr.clone(), 0.75);
        let hot: Vec<bool> = px.iter().map(|p| p.energy_drift_max > cut).collect();
        let mut st: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();

        let (mut lift, mut chord) = (f64::NAN, f64::NAN);
        if let Some((h0, d0)) = &ref0 {
            // Base rate is 0.25 by construction, so 1.0 is chance and 4.0 is perfect agreement.
            let n_h0 = h0.iter().filter(|x| **x).count().max(1) as f64;
            let bse = hot.iter().filter(|x| **x).count() as f64 / n as f64;
            let p = (0..n).filter(|&i| h0[i] && hot[i]).count() as f64 / n_h0;
            lift = p / bse.max(f64::MIN_POSITIVE);
            let mut c: Vec<f64> =
                (0..n).map(|i| (d[i] - d0[i]).abs()).filter(|x| x.is_finite()).collect();
            chord = q(&mut c, 0.5);
        }
        println!(
            "  {label:>26} {ns:>7} {eta:>9.4} {:>4} {:>11.3e} {:>11.3e} {lift:>12.3} {chord:>12.3e}",
            if lc { "on" } else { "OFF" },
            q(&mut st, 0.5),
            q(&mut dr, 0.5),
        );

        let buf: Vec<u8> = px
            .iter()
            .flat_map(|p| {
                let x = p.energy_drift_max;
                if x.is_finite() && x > 0.0 {
                    ramp((lg(x) - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    [255, 0, 255]
                }
            })
            .collect();
        let slug = label.replace(' ', "_").replace('=', "");
        let p = format!("{dir}/drift_{slug}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p, &cfg,
            &format!("res={res}x{res}\ncase=config_stability\narm={label}\n\
                      drift ramp=({DLO:e},{DHI:e}) FIXED and shared across arms\n"),
        );
        if ref0.is_none() {
            ref0 = Some((hot, d));
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **`steps p50` is the control and must be roughly flat across the CONTROLLED rows.** If\n\
         it moves, `eta` did not hold the step size and the comparison is a step-size comparison\n\
         wearing a re-registration label.\n\n\
         `hot lift` is against the baseline's hot set: 1.0 is chance, 4.0 is perfect agreement.\n\
         **A high lift means the wedges survive the change; a lift near 1 means they moved.**\n\
         `chord p50` is the median |delta log10 drift|, the magnitude.\n\n\
         The `confounded` row is what the uncontrolled experiment would have claimed. Compare it\n\
         against `n=250 CONTROLLED`: if they differ, the control was load-bearing."
    );
}
