//! **How large a perturbation does crossing the argmax surface actually inject?**
//!
//! `ref_tie -> 1` locates the selector's switching surface and says nothing about its amplitude.
//! This measures the amplitude directly. At the first sync boundary where two neighbouring pixels'
//! reference itineraries diverge, take the incoming Cartesian state and integrate that one
//! interval under **both** the chosen chart and the runner-up:
//!
//! ```text
//!   dr_chart = || r_win(x_n) - r_alt(x_n) ||      at a COMMON PHYSICAL TIME
//!   dr_step  = || r_win at eta - r_win at eta/2 ||   same treatment, same convention
//! ```
//!
//! **Position and velocity are reported separately**, never as one norm: they are dimensionally
//! different and a combined Euclidean norm is arbitrary unless phase space has been explicitly
//! non-dimensionalised, which it has not been.
//!
//! **And the arms are compared at the same physical time, not after equal fictitious-time
//! increments.** The landing residual is `O(h^2)` and `A*B` differs between charts, so the two
//! arms stop at different `t`; without correcting for it the time transformations manufacture a
//! branch discrepancy on their own. `dt_mismatch` and the uncorrected `dr_raw` are both printed,
//! so the size of the confound is visible rather than trusted.
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

    let measure = |set: &[(usize, usize)]| -> Vec<prin_rs::integrate::az::BranchJump<f64>> {
        set
            .par_iter()
            .filter_map(|&(i, k)| {
                let (x, y) = sl.decode_pos(i);
                let st = grid::decode_state(&chart, 0, x, y);
                let b = az::branch_jump::<f64>(
                    st.s, &st.m, cfg.t_max, cfg.n_sync, cfg.eta, cfg.max_steps, &opts, k,
                );
                if b.ok && b.dr_step > 0.0 && b.dr_chart.is_finite() {
                    Some(b)
                } else {
                    None
                }
            })
            .collect()
    };

    // THE CONTROL: the same pixels, at a boundary where the itinerary does NOT diverge.
    // Without it, a large ratio at a switch could be true at every boundary and mean nothing.
    let control: Vec<(usize, usize)> = sample
        .iter()
        .map(|&(i, k)| (i, if k > 4 { k / 2 } else { k + 1 }))
        .filter(|&(i, k)| k < px[i].ref_path.len())
        .collect();

    println!(
        "  {:>26} {:>7} {:>11} {:>11} {:>11} {:>11} {:>10} {:>11} {:>10}",
        "population", "n", "dr_chart", "dr_step", "R_pos p50", "R_pos p90", "R_vel p50",
        "dt_mismatch", "ref_tie"
    );
    for (name, set) in [("at the first divergence", &sample), ("control: a non-switch", &control)] {
        let b = measure(set);
        let col = |f: &dyn Fn(&prin_rs::integrate::az::BranchJump<f64>) -> f64| -> Vec<f64> {
            b.iter().map(|x| f(x)).filter(|x| x.is_finite()).collect()
        };
        let mut rp = col(&|x| x.dr_chart / x.dr_step);
        let mut rv = col(&|x| x.dv_chart / x.dv_step);
        println!(
            "  {name:>26} {:>7} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>10.3e} {:>11.3e} {:>10.4}",
            b.len(),
            q(&mut col(&|x| x.dr_chart), 0.5),
            q(&mut col(&|x| x.dr_step), 0.5),
            q(&mut rp.clone(), 0.5),
            q(&mut rp, 0.9),
            q(&mut rv, 0.5),
            q(&mut col(&|x| x.dt_mismatch), 0.5),
            q(&mut col(&|x| x.ref_tie), 0.5),
        );
        // The confound, made checkable: how much of the raw jump was the time mismatch?
        let mut expl = col(&|x| (x.speed * x.dt_mismatch) / x.dr_chart.max(f64::MIN_POSITIVE));
        let mut raw = col(&|x| x.dr_raw / x.dr_chart.max(f64::MIN_POSITIVE));
        println!(
            "  {:>26} time-mismatch displacement / dr_chart p50 {:.3e} p90 {:.3e};  \
             raw/corrected p50 {:.4}",
            "",
            q(&mut expl.clone(), 0.5),
            q(&mut expl, 0.9),
            q(&mut raw, 0.5)
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
         which is why the local step error is the normaliser.\n\n\
         **And read the confound line before the ratio.** If the time-mismatch displacement is a\n\
         material fraction of `dr_chart`, the arms were compared at different physical times and\n\
         the difference is partly the time transformations rather than the charts. `raw/corrected`\n\
         near 1 means the correction changed nothing and the comparison was clean anyway."
    );
}
