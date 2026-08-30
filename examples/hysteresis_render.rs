//! **The intervention test, as pictures.** Move the chart-switching surfaces and see whether the
//! wedges move with them.
//!
//! The statistics say the reference argmax is a *symptom*: crossing it costs 1.25x the local step
//! error against 1.07x at a non-switching boundary. But the wedges are a **visual** complaint, and
//! a ratio is a poor answer to one. This is the direct test.
//!
//! `ref_hysteresis` holds the current reference until a rival beats its opposite side by `1+eps`.
//! That **displaces every switching surface in initial-condition space** while leaving every
//! trajectory a legitimate integration of the same Hamiltonian.
//!
//! ```text
//!   wedges MOVE with eps    ->  they are chart-selection artefacts after all, and the
//!                               branch-jump amplitude was measuring the wrong thing
//!   wedges STAY put         ->  they are dynamical; the itinerary partition merely reports
//!                               structure the trajectories already have
//! ```
//!
//! **The control that makes it readable:** the cell mask itself. Hysteresis changes the itinerary
//! *by construction*, so the mask MUST move — if it does not, `eps` is not doing anything and the
//! whole panel is a null run rather than a result. The mask moving while the field does not is
//! the finding; both moving, or neither, means something else.
//!
//! Every ramp is a fixed constant shared across arms. Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const SLO: f64 = 1e-6;
const SHI: f64 = 1e0;
/// The DRIFT window, fixed and shared across arms. The bright wedges live in this field, not in
/// `spread_shape`, and the first cut of this harness rendered only the latter.
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
    let dir = format!("{root}/step_control/hysteresis");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let base = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        keep_ref_path: true,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", base.provenance());
    println!(
        "  {:>8} {:>8} {:>11} {:>12} {:>12} {:>12} {:>12}",
        "eps", "secs", "switches", "cell frac", "mask moved", "shape chord", "DRIFT chord p50"
    );

    let mut ref0: Option<(Vec<bool>, Vec<f64>, Vec<f64>)> = None;
    for eps in [0.0f64, 0.02, 0.05, 0.20] {
        let cfg = EnsembleCfg { ref_hysteresis: eps, ..base };
        let t = std::time::Instant::now();
        let px: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
            .collect();
        let secs = t.elapsed().as_secs_f64();

        let cell: Vec<bool> = (0..px.len())
            .map(|i| {
                let (x, y) = (i % res, i / res);
                let d = |j: usize| px[i].ref_path != px[j].ref_path;
                (x + 1 < res && d(i + 1)) || (y + 1 < res && d(i + res))
            })
            .collect();
        let field: Vec<f64> = px.iter().map(|p| p.spread_shape).collect();
        let drift: Vec<f64> = px.iter().map(|p| p.energy_drift_max).collect();

        let (mut mask_moved, mut chord, mut drift_chord) = (f64::NAN, f64::NAN, f64::NAN);
        let _ = &drift;
        if let Some((c0, f0, _)) = &ref0 {
            mask_moved = (0..cell.len()).filter(|&i| cell[i] != c0[i]).count() as f64
                / cell.len() as f64;
            // Magnitude, not just a moved count: **a count of pixels differing in the last bit
            // is a fact about a chaotic field, not about the intervention.**
            let mut d: Vec<f64> = (0..field.len())
                .map(|i| (field[i] - f0[i]).abs())
                .filter(|x| x.is_finite())
                .collect();
            chord = q(&mut d, 0.5);
            // The drift field is the one the wedges are in, and it is ramped in LOG, so the
            // honest magnitude is a log difference: an absolute difference on a quantity
            // spanning fourteen decades is dominated by wherever the field happens to be large.
            let d0 = &ref0.as_ref().unwrap().2;
            let lg = |x: f64| if x.is_finite() && x > 0.0 { x.log10() } else { DLO.log10() };
            let mut dd: Vec<f64> = (0..drift.len())
                .map(|i| (lg(drift[i]) - lg(d0[i])).abs())
                .filter(|x| x.is_finite())
                .collect();
            drift_chord = q(&mut dd, 0.5);
        }
        println!(
            "  {eps:>8.2} {secs:>8.1} {:>11.3e} {:>12.4} {:>12.4} {:>12.4} {:>12.3e}",
            px.iter().map(|p| p.switches as f64).sum::<f64>() / px.len() as f64,
            cell.iter().filter(|x| **x).count() as f64 / cell.len() as f64,
            mask_moved,
            chord,
            drift_chord
        );

        let fbuf: Vec<u8> = field
            .iter()
            .flat_map(|&x| {
                if x.is_finite() && x > 0.0 {
                    ramp((x.log10() - SLO.log10()) / (SHI.log10() - SLO.log10()))
                } else if x == 0.0 {
                    ramp(0.0)
                } else {
                    [255, 0, 255]
                }
            })
            .collect();
        let mbuf: Vec<u8> =
            cell.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect();
        let dbuf: Vec<u8> = drift
            .iter()
            .flat_map(|&x| {
                if x.is_finite() && x > 0.0 {
                    ramp((x.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    [255, 0, 255]
                }
            })
            .collect();
        for (n, b) in [("field", &fbuf), ("drift", &dbuf), ("cells", &mbuf)] {
            let p = format!("{dir}/{n}_eps{eps:.2}.png");
            let _ = prin_rs::output::adaptive::save_rect(&p, res, res, b);
            let _ = prin_rs::output::provenance_sidecar(
                &p,
                &cfg,
                &format!(
                    "res={res}x{res}\ncase=config_stability\nref_hysteresis={eps}\n\
                     field=spread_shape  ramp=({SLO:e},{SHI:e})  <- FIXED, shared across arms\n"
                ),
            );
        }
        if ref0.is_none() {
            ref0 = Some((cell, field, drift));
        }
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **`mask moved` is the control and must be large**, because hysteresis changes the\n\
         itinerary by construction. If it is ~0 then `eps` did nothing and every other column is\n\
         a null run rather than a result.\n\n\
         Then compare the PICTURES, not the moved counts. `field moved` counts pixels differing\n\
         in the last bit, which on a chaotic field at t = 50 is a fact about the field -- it ran\n\
         at 0.87-0.93 for step-control changes that were plainly correct. `field chord p50` is\n\
         the magnitude, and the panels are the answer: **if the wedge OUTLINES sit in the same\n\
         places at every eps while the cell mask visibly rearranges, the wedges are not drawn by\n\
         the selector.**"
    );
}
