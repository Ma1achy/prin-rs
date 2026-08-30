//! **Softmax instead of argmax on the reference body — does it smooth the field?**
//!
//! `choose_reference` is a bare `argmax` over the three pair separations, re-evaluated at every
//! sync boundary with no hysteresis. Measured on `config_stability`: the regions where
//! neighbouring pixels take the **same** reference sequence over all 125 boundaries are exactly
//! the wedges in the drift field, and their boundaries are where the field has its edges.
//!
//! # Why this is well posed, and where the blend has to happen
//!
//! Not inside a step — you cannot be in two regularised charts at once. But at every sync
//! boundary the state is **Cartesian and chart-free**, and the reference choice governs only how
//! the next interval is integrated. All three choices approximate the same true trajectory, so a
//! convex combination of their endpoints is another approximation to it and converges to the same
//! limit. See [`prin_rs::integrate::az::integrate_softref`].
//!
//! # What would make this a failure, stated first
//!
//! - **Roughness does not fall.** Then the edges are not the argmax after all, and the
//!   `same-reference-path` correlation was coincidental.
//! - **Roughness falls and drift rises.** Blending two approximations is not itself a solution of
//!   the ODE, so it can trade smoothness for accuracy. Both columns are printed; a smoother field
//!   that conserves energy worse is not an improvement.
//! - **`arms/boundary` is near 2 or 3 everywhere.** Then it is not a thin-shell cost and it is
//!   three integrators, not one. The whole argument for cheapness is that away from a tie exactly
//!   one arm survives the prune.
//!
//! `temp = 0` is exact `argmax` including its first-maximum tie-break, so the first row is the
//! shipped path and every other row is measured against it.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts, StepLimit};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::energy;

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

fn grad(f: &[f64], res: usize) -> Vec<f64> {
    (0..f.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let mut g: f64 = 0.0;
            if x + 1 < res {
                g = g.max((f[i + 1] - f[i]).abs());
            }
            if y + 1 < res {
                g = g.max((f[i + res] - f[i]).abs());
            }
            g
        })
        .collect()
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/softref");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let cfg = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!(
        "config_stability {res}^2, NOMINAL COPY ONLY (this integrator carries no ensemble, no\n\
         events and no outcome -- a blended state has no single terminal class, and inventing one\n\
         would put a discontinuity back in by another door).\nconfig: {}\n",
        cfg.provenance()
    );

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

    println!(
        "  {:>8} {:>9} {:>12} {:>12} {:>12} {:>12} {:>12} {:>9}",
        "temp", "secs", "rough p50", "rough p90", "drift p50", "drift p99", "arms/bound", "nonfin"
    );
    let mut base_rough = f64::NAN;
    for temp in [0.0f64, 0.01, 0.03, 0.10, 0.30] {
        let t0 = std::time::Instant::now();
        let out: Vec<(f64, u64, bool)> = (0..sl.npix())
            .into_par_iter()
            .map(|k| {
                let (x, y) = sl.decode_pos(k);
                let st = grid::decode_state(&chart, 0, x, y);
                let o = az::integrate_softref::<f64>(
                    st.s, &st.m, cfg.t_max, cfg.n_sync, cfg.eta, cfg.max_steps, &opts, temp,
                );
                let _ = energy::energy(&st.s.r, &st.s.v, &st.m, 0.0);
                (o.drift, o.arms, o.finite)
            })
            .collect();
        let secs = t0.elapsed().as_secs_f64();

        let lg: Vec<f64> = out
            .iter()
            .map(|o| if o.0.is_finite() && o.0 > 0.0 { o.0.log10() } else { DLO.log10() })
            .collect();
        let mut g = grad(&lg, res);
        let mut d: Vec<f64> = out.iter().map(|o| o.0).filter(|x| x.is_finite()).collect();
        let arms = out.iter().map(|o| o.1).sum::<u64>() as f64
            / (out.len() * cfg.n_sync) as f64;
        let r50 = q(&mut g.clone(), 0.5);
        if temp == 0.0 {
            base_rough = r50;
        }
        println!(
            "  {temp:>8.2} {secs:>9.1} {r50:>12.4} {:>12.4} {:>12.3e} {:>12.3e} {arms:>12.4} {:>9.4}",
            q(&mut g, 0.9),
            q(&mut d.clone(), 0.5),
            q(&mut d, 0.99),
            out.iter().filter(|o| !o.2).count() as f64 / out.len() as f64
        );

        let buf: Vec<u8> = lg
            .iter()
            .flat_map(|&x| ramp((x - DLO.log10()) / (DHI.log10() - DLO.log10())).into_iter())
            .collect();
        let name = format!("{dir}/softref_temp{:.2}.png", temp);
        let _ = prin_rs::output::adaptive::save_rect(&name, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &name,
            &cfg,
            &format!("res={res}x{res}\ncase=config_stability\nsoftref_temp={temp}\n\
                      NOTE: nominal copy only, no events, no outcome\n"),
        );
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         `temp = 0` is exact argmax and is the reference row. **Roughness is the median\n\
         |grad log10 drift| -- a per-pixel number, not an image impression** (base {base_rough:.4}).\n\n\
         A win is roughness DOWN with drift FLAT and `arms/bound` near 1. Roughness down with\n\
         drift up is a trade, not a fix: blending two approximations is not itself a solution of\n\
         the ODE. `arms/bound` near 2 means the tie shell is not thin and this is two integrators."
    );
}
