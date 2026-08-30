//! **What ARE the bright wedges in the drift field?**
//!
//! An earlier analysis aligned the reference-cell mask against the top decile of `|grad drift|` —
//! the *edges* of the drift field. That is not the feature in question. The wedges are large,
//! bright, sharply-bounded regions of **high drift magnitude**, and a gradient mask and a
//! magnitude mask are different objects. Stated because it is the mistake this harness exists to
//! correct.
//!
//! # What is asked
//!
//! Take `hot = drift > p90` — the bright regions themselves — and ask what they coincide with,
//! each against its own base rate because a candidate covering half the frame has a lift of ~1
//! by arithmetic:
//!
//! ```text
//!   coherent cell        neighbours share the whole reference itinerary
//!   many ref switches    the trajectory changed chart often
//!   terminated early     t_end < t_max: collision or escape stopped it
//!   deep close approach  small d_min
//!   expensive            high step count
//! ```
//!
//! **The point is that these are rival explanations, not one hypothesis with controls.** If the
//! hot regions are simply where trajectories terminate, or where they pass closest, then the
//! reference partition is irrelevant to them and the whole chart-selection thread was about the
//! wrong feature.
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
    let dir = format!("{root}/step_control/wedge_id");
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
        keep_ref_path: true,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", cfg.provenance());

    let t0 = std::time::Instant::now();
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    println!("{:.1}s\n", t0.elapsed().as_secs_f64());
    let n = px.len();

    let mut d: Vec<f64> = px.iter().map(|p| p.energy_drift_max).collect();
    let p90 = q(&mut d.clone(), 0.90);
    let p75 = q(&mut d, 0.75);
    let hot: Vec<bool> = px.iter().map(|p| p.energy_drift_max > p90).collect();

    // The rival explanations, each a boolean over the frame.
    let coherent: Vec<bool> = (0..n)
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let same = |j: usize| px[i].ref_path == px[j].ref_path;
            !((x + 1 < res && !same(i + 1)) || (y + 1 < res && !same(i + res)))
        })
        .collect();
    let mut sw: Vec<f64> = px.iter().map(|p| p.switches as f64).collect();
    let sw_hi = q(&mut sw, 0.75);
    let mut dm: Vec<f64> = px.iter().map(|p| p.d_min_true).filter(|x| x.is_finite()).collect();
    let dm_lo = q(&mut dm, 0.25);
    let mut st: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
    let st_hi = q(&mut st, 0.75);

    let cands: Vec<(&str, Vec<bool>)> = vec![
        ("inside a coherent cell", coherent.clone()),
        ("many ref switches (>p75)", px.iter().map(|p| p.switches as f64 > sw_hi).collect()),
        ("terminated early", px.iter().map(|p| p.t_end < cfg.t_max - 1e-9).collect()),
        ("deep close approach (<p25)", px.iter().map(|p| p.d_min_true < dm_lo).collect()),
        ("expensive (steps > p75)", px.iter().map(|p| p.total_substeps as f64 > st_hi).collect()),
        ("collided", px.iter().map(|p| p.state == 2).collect()),
        ("escaped", px.iter().map(|p| p.state == 0).collect()),
    ];

    println!(
        "== WHAT COINCIDES WITH THE BRIGHT WEDGES (drift > p90 = {p90:.3e})? ==\n\
         Base rate first: **a candidate covering half the frame has a lift of ~1 by arithmetic.**\n"
    );
    println!("  {:>28} {:>11} {:>14} {:>9}", "candidate", "base rate", "P(cand | hot)", "lift");
    let n_hot = hot.iter().filter(|x| **x).count().max(1) as f64;
    for (name, m) in &cands {
        let base = m.iter().filter(|x| **x).count() as f64 / n as f64;
        let ph = (0..n).filter(|&i| hot[i] && m[i]).count() as f64 / n_hot;
        println!(
            "  {name:>28} {base:>11.4} {ph:>14.4} {:>9.3}",
            ph / base.max(f64::MIN_POSITIVE)
        );
    }

    // And the reverse question, which is the one that decides: is drift SMOOTH inside a cell?
    let mut din: Vec<f64> = (0..n).filter(|&i| coherent[i]).map(|i| px[i].energy_drift_max).collect();
    let mut dout: Vec<f64> =
        (0..n).filter(|&i| !coherent[i]).map(|i| px[i].energy_drift_max).collect();
    println!(
        "\n  drift inside a coherent cell:  p50 {:.3e}  p90 {:.3e}   (n = {})\n  \
           drift elsewhere:               p50 {:.3e}  p90 {:.3e}   (n = {})",
        q(&mut din.clone(), 0.5),
        q(&mut din.clone(), 0.9),
        din.len(),
        q(&mut dout.clone(), 0.5),
        q(&mut dout.clone(), 0.9),
        dout.len()
    );

    // --- panels, all at the same size so they can be compared directly --------------------
    let save = |name: &str, buf: &[u8]| {
        let p = format!("{dir}/{name}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p,
            &cfg,
            &format!("res={res}x{res}\ncase=config_stability\npanel={name}\n\
                      drift ramp=({DLO:e},{DHI:e})  <- FIXED\nhot cut = drift p90 = {p90:e}\n"),
        );
    };
    save(
        "drift",
        &px.iter()
            .flat_map(|p| {
                if p.energy_drift_max.is_finite() && p.energy_drift_max > 0.0 {
                    ramp((p.energy_drift_max.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    [255, 0, 255]
                }
            })
            .collect::<Vec<u8>>(),
    );
    let mask = |m: &[bool]| -> Vec<u8> {
        m.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect()
    };
    save("hot_p90", &mask(&hot));
    save(
        "hot_p75",
        &mask(&px.iter().map(|p| p.energy_drift_max > p75).collect::<Vec<bool>>()),
    );
    save("coherent_cells", &mask(&coherent));
    for (name, m) in &cands[1..] {
        save(&name.split(' ').next().unwrap().to_lowercase(), &mask(m));
    }
    println!("\nWrote {dir}/ -- drift, hot_p90, hot_p75, coherent_cells and one panel per rival.");
}
