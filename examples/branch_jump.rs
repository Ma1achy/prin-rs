//! **How large a perturbation does crossing the argmax surface actually inject?**
//!
//! `ref_tie -> 1` locates the selector's switching surface and says nothing about its amplitude.
//! This measures the amplitude directly. At the first sync boundary where two neighbouring pixels'
//! reference itineraries diverge, take the incoming Cartesian state and integrate that one
//! interval under **both** the chosen chart and the runner-up:
//!
//! ```text
//!   delta_chart = || Phi_win(x_n) - Phi_alt(x_n) ||
//!   delta_step  = || Phi_win(x_n) at eta  -  Phi_win(x_n) at eta/2 ||
//! ```
//!
//! **The ratio is the finding.** `delta_step` is the ordinary local truncation error of that same
//! interval in that same chart, so `delta_chart / delta_step` answers: *is the chart jump larger
//! than the numerical noise the integrator already carries?* Without the normaliser, an absolute
//! `delta_chart` is unreadable — a large number on a violent interval and a small one on a quiet
//! interval mean the same thing.
//!
//! A ratio near 1 would say crossing the surface costs no more than an ordinary step, and the
//! selector is a **symptom** — chaotic divergence changed the geometry, which changed the chart.
//! A ratio far above 1 says the selector is a **seed**: it injects a jump the dynamics then
//! amplifies over the remaining `t_max - t_n`, which is the mechanism
//!
//! ```text
//!   d_i - d_j -> 0   =>   Phi_i - Phi_j != 0   =>   delta(t) ~ e^(lambda t) delta(0)   =>   wedge
//! ```
//!
//! # The control
//!
//! The same measurement at boundaries where the itinerary does **not** diverge, on the same
//! pixels. Without it, "the jump is 40x the step error at a switch" could be true at every
//! boundary of a chaotic slice and would mean nothing about switching.
//!
//! Args: `res cap root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts, StepLimit};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;

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

fn main() {
    let res: usize = arg(1, 192);
    let cap: usize = arg(2, 400);
    let root: String = std::env::args().nth(3).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    let (chart, cx, cy, half) = Chart::config_stability();
    let cfg = EnsembleCfg {
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
    println!("config_stability {res}^2\nconfig: {}\n", cfg.provenance());

    let opts = AzOpts::<f64> {
        step_limit: StepLimit::Predictive,
        step_limit_f: cfg.step_limit_f,
        step_blend: cfg.step_blend,
        blend_p: cfg.blend_p,
        dtau_mode: cfg.dtau_mode,
        clamp_final_step: cfg.clamp_final_step,
        lc_stable: cfg.lc_stable,
        ..Default::default()
    };

    let t0 = std::time::Instant::now();
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    println!("itinerary pass {:.1}s\n", t0.elapsed().as_secs_f64());

    // Neighbour pairs whose itineraries diverge, and where first.
    let mut pairs: Vec<(usize, usize)> = Vec::new(); // (pixel, first differing boundary)
    for i in 0..px.len() {
        let (x, y) = (i % res, i / res);
        for j in [
            if x + 1 < res { Some(i + 1) } else { None },
            if y + 1 < res { Some(i + res) } else { None },
        ]
        .into_iter()
        .flatten()
        {
            let (a, b) = (&px[i].ref_path, &px[j].ref_path);
            if let Some(k) = (0..a.len().min(b.len())).find(|&k| a[k] != b[k]) {
                pairs.push((i, k));
            }
        }
    }
    // Evenly spaced, never random: the same pairs every run and no seed to report.
    let take = |v: &[(usize, usize)], n: usize| -> Vec<(usize, usize)> {
        if v.len() <= n { v.to_vec() } else { (0..n).map(|k| v[k * v.len() / n]).collect() }
    };
    let sample = take(&pairs, cap);
    println!(
        "  {} neighbour pairs diverge; measuring {} (cap {cap}).\n  \
         **The cap is printed because a silently truncated set reads as full coverage.**\n",
        pairs.len(),
        sample.len()
    );

    let measure = |set: &[(usize, usize)]| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let out: Vec<(f64, f64, f64)> = set
            .par_iter()
            .filter_map(|&(i, k)| {
                let (x, y) = sl.decode_pos(i);
                let st = grid::decode_state(&chart, 0, x, y);
                let b = az::branch_jump::<f64>(
                    st.s, &st.m, cfg.t_max, cfg.n_sync, cfg.eta, cfg.max_steps, &opts, k,
                );
                if b.ok && b.delta_step > 0.0 && b.delta_chart.is_finite() {
                    Some((b.delta_chart, b.delta_step, b.ref_tie))
                } else {
                    None
                }
            })
            .collect();
        (
            out.iter().map(|o| o.0).collect(),
            out.iter().map(|o| o.1).collect(),
            out.iter().map(|o| o.2).collect(),
        )
    };

    // THE CONTROL: the same pixels, at a boundary where the itinerary does NOT diverge.
    // Without it, a large ratio at a switch could be true at every boundary and mean nothing.
    let control: Vec<(usize, usize)> = sample
        .iter()
        .map(|&(i, k)| (i, if k > 4 { k / 2 } else { k + 1 }))
        .filter(|&(i, k)| k < px[i].ref_path.len())
        .collect();

    println!(
        "  {:>26} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "population", "dchart p50", "dstep p50", "RATIO p50", "RATIO p90", "ref_tie p50"
    );
    for (name, set) in [("at the first divergence", &sample), ("control: a non-switch", &control)] {
        let (dc, ds, rt) = measure(set);
        let mut ratio: Vec<f64> =
            dc.iter().zip(ds.iter()).map(|(a, b)| a / b).filter(|x| x.is_finite()).collect();
        println!(
            "  {name:>26} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>10.4}",
            q(&mut dc.clone(), 0.5),
            q(&mut ds.clone(), 0.5),
            q(&mut ratio.clone(), 0.5),
            q(&mut ratio, 0.9),
            q(&mut rt.clone(), 0.5)
        );
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **RATIO is the finding.** Near 1 means crossing the selector's surface costs no more\n\
         than an ordinary step, and the argmax is a SYMPTOM -- the trajectories had already\n\
         diverged and the geometry followed. Far above 1 means it is a SEED: a jump the dynamics\n\
         then amplifies over the remaining horizon.\n\n\
         Read it against the control row, never alone. And `dchart` alone is unreadable: a large\n\
         number on a violent interval and a small one on a quiet interval say the same thing,\n\
         which is why the local step error is the normaliser."
    );
}
