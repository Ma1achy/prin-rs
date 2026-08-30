//! **Are the early switching lines the EDGES of the drift wedges?**
//!
//! Every previous test asked the wrong pairing. `edge_anatomy` compared `|grad drift|` edges
//! against the ref-path mask at *all* `k` — a mask covering 72% of the frame, saturated and
//! unreadable. `wedge_id` compared hot *interiors* against coherent cells — lift 1.49, weak.
//! Neither asked the question actually on the table: **does the boundary of the hot set lie on
//! the early, straight switching surfaces?**
//!
//! That pairing is the one with power, because the early set is *small*: `first_k <= 4` is ~5.5%
//! of the frame, where every earlier candidate was 25-72% and had a lift near 1 by arithmetic.
//!
//! # Why an early surface could draw an edge when `branch_jump` says the jump is small
//!
//! `branch_jump` measured **amplitude**: crossing a chart boundary displaces the state by ~1.25x
//! the local step error. It did **not** measure **coherence**, and that is the distinction that
//! matters here. Both perturbations are amplified by the same `e^(lambda t)` over the remaining
//! horizon and both saturate — but the chart jump is *systematic across the surface* (every pixel
//! on one side takes the same branch) while step error is incoherent pixel to pixel. A coherent
//! perturbation across a line produces a systematic difference between the two sides, which is an
//! **edge**; an incoherent one of the same size produces noise. So a small-amplitude crossing can
//! still draw a sharp boundary, and `branch_jump` alone cannot rule it out.
//!
//! # Dilation, because coincidence needs a tolerance
//!
//! The wedge outline and the switching line may sit a pixel or two apart — the hot cut is a
//! quantile on a continuous field. The lift is reported at dilation 0, 1 and 2 so the answer does
//! not hang on exact pixel registration, and the base rate rises with dilation so the comparison
//! stays honest.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn dilate(m: &[bool], res: usize, r: i64) -> Vec<bool> {
    if r == 0 {
        return m.to_vec();
    }
    (0..m.len())
        .map(|i| {
            let (cx, cy) = ((i % res) as i64, (i / res) as i64);
            for dy in -r..=r {
                for dx in -r..=r {
                    let (x, y) = (cx + dx, cy + dy);
                    if x >= 0 && y >= 0 && x < res as i64 && y < res as i64 && m[y as usize * res + x as usize] {
                        return true;
                    }
                }
            }
            false
        })
        .collect()
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/wedge_edge");
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

    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    let n = px.len();

    let mut d: Vec<f64> = px.iter().map(|p| p.energy_drift_max).collect();
    let p75 = q(&mut d, 0.75);
    let hot: Vec<bool> = px.iter().map(|p| p.energy_drift_max > p75).collect();
    // The wedge OUTLINE: hot on one side, not hot on the other.
    let hot_edge: Vec<bool> = (0..n)
        .map(|i| {
            let (x, y) = (i % res, i / res);
            (x + 1 < res && hot[i] != hot[i + 1]) || (y + 1 < res && hot[i] != hot[i + res])
        })
        .collect();

    let first: Vec<Option<usize>> = (0..n)
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let mut best: Option<usize> = None;
            for j in [
                if x + 1 < res { Some(i + 1) } else { None },
                if y + 1 < res { Some(i + res) } else { None },
            ]
            .into_iter()
            .flatten()
            {
                let (a, b) = (&px[i].ref_path, &px[j].ref_path);
                if let Some(k) = (0..a.len().min(b.len())).find(|&k| a[k] != b[k]) {
                    best = Some(best.map_or(k, |p: usize| p.min(k)));
                }
            }
            best
        })
        .collect();

    println!(
        "== DOES THE WEDGE OUTLINE LIE ON THE EARLY SWITCHING LINES? ==\n\
         hot = drift > p75; the outline is where hot meets not-hot. **The early set is small\n\
         ({:.4} of the frame), which is what gives this pairing power where every earlier\n\
         candidate sat at 25-72% and had a lift near 1 by arithmetic.**\n",
        first.iter().filter(|f| f.map_or(false, |k| k <= 4)).count() as f64 / n as f64
    );
    println!(
        "  {:>26} {:>6} {:>11} {:>15} {:>9}",
        "switching set", "dilate", "base rate", "P(set | outline)", "lift"
    );
    let n_edge = hot_edge.iter().filter(|x| **x).count().max(1) as f64;
    for (name, lo, hi) in [
        ("first_k <= 4 (straight)", 0usize, 5usize),
        ("first_k 5-9", 5, 10),
        ("first_k 10-24", 10, 25),
        ("first_k >= 25", 25, cfg.n_sync),
        ("any k (the old mask)", 0, cfg.n_sync),
    ] {
        let m: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k >= lo && k < hi)).collect();
        for r in [0i64, 1, 2] {
            let dm = dilate(&m, res, r);
            let base = dm.iter().filter(|x| **x).count() as f64 / n as f64;
            let pe = (0..n).filter(|&i| hot_edge[i] && dm[i]).count() as f64 / n_edge;
            println!(
                "  {:>26} {r:>6} {base:>11.4} {pe:>15.4} {:>9.3}",
                if r == 0 { name } else { "" },
                pe / base.max(f64::MIN_POSITIVE)
            );
        }
    }

    let mask = |m: &[bool]| -> Vec<u8> {
        m.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect()
    };
    // The two sets in one image, so coincidence is visible and not only tabulated.
    let early: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k <= 4)).collect();
    let both: Vec<u8> = (0..n)
        .flat_map(|i| match (hot_edge[i], early[i]) {
            (true, true) => [255u8, 255, 255],  // both
            (true, false) => [220, 60, 60],     // wedge outline only
            (false, true) => [60, 200, 255],    // early switching line only
            _ => [12, 12, 16],
        })
        .collect();
    // **The picture that actually answers the question**: the drift field itself with the early
    // switching lines drawn on it. A mask-against-mask overlay can only show where two derived
    // sets agree; this shows whether the lines bound the bright regions a reader is looking at.
    const DLO: f64 = 1e-12;
    const DHI: f64 = 1e2;
    let rmp = |x: f64| -> [u8; 3] {
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
    };
    let mut dl: Vec<u8> = Vec::with_capacity(n * 3);
    for p in &px {
        let x = p.energy_drift_max;
        dl.extend_from_slice(&if x.is_finite() && x > 0.0 {
            rmp((x.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
        } else {
            [255, 0, 255]
        });
    }
    let mut dl_lines = dl.clone();
    for i in 0..n {
        if early[i] {
            dl_lines[3 * i] = 0;
            dl_lines[3 * i + 1] = 255;
            dl_lines[3 * i + 2] = 255;
        }
    }
    for (nm, b) in [
        ("drift", dl),
        ("drift_with_early_lines", dl_lines),
        ("hot_outline", mask(&hot_edge)),
        ("early_switch_lines", mask(&early)),
        ("overlay", both),
    ] {
        let p = format!("{dir}/{nm}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &b);
        let _ = prin_rs::output::provenance_sidecar(
            &p,
            &cfg,
            &format!(
                "res={res}x{res}\ncase=config_stability\npanel={nm}\n\
                 overlay: WHITE = both, RED = wedge outline only, CYAN = early switching line only\n"
            ),
        );
    }
    println!("\nWrote {dir}/ -- hot_outline, early_switch_lines, overlay (white = both).");
}
