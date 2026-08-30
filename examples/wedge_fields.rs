//! **The cell boundaries drawn ON the fields, not just correlated against them.**
//!
//! `wedge_origin` measured the ratio of median `|grad log10 field|` across a reference-cell
//! boundary to inside one — 3.235 on `energy_drift_max`, **2.445 on `spread_shape`**, 2.769 on
//! `ensemble_spread`, and 0.000 on `t_end`. Those are numbers about an image nobody had drawn.
//! This draws it: each field on a fixed ramp, and the same field with the cell boundaries
//! overlaid, so "the partition draws this field's edges" can be checked by eye against the
//! statistic that claims it.
//!
//! `spread_shape` is the one that matters — it is what the shipping bivariate colouring reads for
//! lightness, so it is the difference between *the wedges are in the error diagnostic* and *the
//! wedges are in the picture*.
//!
//! **Every ramp is a fixed constant, never auto-ranged.** A per-field window would stretch each
//! field's own range to full scale and make a field with no structure look like one with plenty
//! — the failure mode this project has on record from the auto-ranged lightness ramp.
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
    let dir = format!("{root}/step_control/wedge_origin");
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

    let cell: Vec<bool> = (0..px.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let d = |j: usize| px[i].ref_path != px[j].ref_path;
            (x + 1 < res && d(i + 1)) || (y + 1 < res && d(i + res))
        })
        .collect();

    // (name, values, lo, hi) -- every window a FIXED constant. See the module docs.
    let fields: [(&str, Vec<f64>, f64, f64); 3] = [
        ("drift", px.iter().map(|p| p.energy_drift_max).collect(), 1e-12, 1e2),
        ("spread_shape", px.iter().map(|p| p.spread_shape).collect(), 1e-6, 1e0),
        ("ensemble_spread", px.iter().map(|p| p.ensemble_spread).collect(), 1e-6, 1e0),
    ];

    for (name, v, lo, hi) in &fields {
        let base: Vec<u8> = v
            .iter()
            .flat_map(|&x| {
                if x.is_finite() && x > 0.0 {
                    ramp((x.log10() - lo.log10()) / (hi.log10() - lo.log10()))
                } else if x == 0.0 {
                    ramp(0.0)
                } else {
                    [255, 0, 255]
                }
            })
            .collect();
        let mut over = base.clone();
        for (i, &c) in cell.iter().enumerate() {
            if c {
                // Cyan, and the boundary is INVERTED against the inferno ramp so it reads on
                // both the dark and the bright end. A boundary drawn in a colour the ramp also
                // uses would be invisible on half the frame.
                over[3 * i] = 0;
                over[3 * i + 1] = 255;
                over[3 * i + 2] = 255;
            }
        }
        for (suf, buf) in [("", &base), ("_cells", &over)] {
            let path = format!("{dir}/field_{name}{suf}.png");
            let _ = prin_rs::output::adaptive::save_rect(&path, res, res, buf);
            let _ = prin_rs::output::provenance_sidecar(
                &path,
                &cfg,
                &format!(
                    "res={res}x{res}\ncase=config_stability\nfield={name}\n\
                     ramp=({lo:e},{hi:e})  <- FIXED, never auto-ranged\n\
                     overlay=cyan where the reference itinerary differs from a neighbour\n"
                ),
            );
        }
        println!("  wrote field_{name}.png and field_{name}_cells.png");
    }

    println!(
        "\n  `spread_shape` is the one that matters: it is what the shipping colouring reads for\n\
         lightness, so it is the difference between *the wedges are in the error diagnostic* and\n\
         *the wedges are in the picture*. Measured ratio across a cell boundary to inside one:\n\
         drift 3.235, spread_shape 2.445, ensemble_spread 2.769."
    );
}
