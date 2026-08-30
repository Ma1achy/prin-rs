//! **The straight black lines cutting through the coherent ribbons — what are they, and why
//! straight?**
//!
//! They are loci where a *single* symbol of the reference itinerary flips: a codimension-1 surface
//! in the chart plane across which one of the `n_sync` argmax decisions changes, inside a region
//! that otherwise agrees on all the others.
//!
//! # The claim, stated so it can fail
//!
//! **Straightness encodes WHEN.** The map from initial condition to the state at boundary `k` is
//! nearly affine for small `k` — the flow has not had time to fold — so the pulled-back tie
//! surface `d_i = d_j` is nearly a straight line at early boundaries and progressively more
//! folded at later ones. `ic_class.png`, the reference partition at `t = 0`, is already on record
//! as a smooth six-sector pinwheel with straight edges; those are the `k = 0` members of the same
//! family.
//!
//! If that is right, the mask of "first differs at exactly `k`" is a set of clean lines for small
//! `k` and a fractal tangle for large `k`. If the early masks are *also* tangled, the claim is
//! wrong and straightness means something else.
//!
//! # Measured, not eyeballed
//!
//! Per `k`, the **local straightness** of the differing set: over a 9x9 window centred on each
//! differing pixel, the structure tensor's anisotropy `(l1 - l2) / (l1 + l2)`. 1 is a perfect
//! line, 0 is isotropic scatter. Reported with the pixel count, because anisotropy over three
//! pixels is not a measurement.
//!
//! # Polarity, stated because two panels in this project disagree on it
//!
//! Here **white = differs at this k**. In `wedge_id/coherent_cells.png` white was *coherent*, and
//! in `wedge_origin/cells_eta.png` white was *differs*. Every panel below says so in its sidecar.
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
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// Structure-tensor anisotropy of a binary set in a window: 1 = a line, 0 = isotropic.
fn straightness(m: &[bool], res: usize, i: usize, half: i64) -> Option<f64> {
    let (cx, cy) = ((i % res) as i64, (i / res) as i64);
    let (mut n, mut sx, mut sy) = (0.0f64, 0.0f64, 0.0f64);
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for dy in -half..=half {
        for dx in -half..=half {
            let (x, y) = (cx + dx, cy + dy);
            if x < 0 || y < 0 || x >= res as i64 || y >= res as i64 {
                continue;
            }
            if m[y as usize * res + x as usize] {
                pts.push((x as f64, y as f64));
                sx += x as f64;
                sy += y as f64;
                n += 1.0;
            }
        }
    }
    // Anisotropy over three pixels is not a measurement.
    if n < 6.0 {
        return None;
    }
    let (mx, my) = (sx / n, sy / n);
    let (mut a, mut b, mut c) = (0.0f64, 0.0f64, 0.0f64);
    for (x, y) in pts {
        let (u, v) = (x - mx, y - my);
        a += u * u;
        b += u * v;
        c += v * v;
    }
    let tr = a + c;
    let det = a * c - b * b;
    let disc = (tr * tr / 4.0 - det).max(0.0).sqrt();
    let (l1, l2) = (tr / 2.0 + disc, tr / 2.0 - disc);
    if l1 + l2 <= 0.0 {
        return None;
    }
    Some((l1 - l2) / (l1 + l2))
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/switch_depth");
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

    // First boundary at which this pixel's itinerary differs from a neighbour's; `None` if it
    // agrees with both all the way.
    let first: Vec<Option<usize>> = (0..px.len())
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

    let n_diff = first.iter().filter(|x| x.is_some()).count();
    println!(
        "  {} of {} pixels differ from a neighbour somewhere ({:.4})\n",
        n_diff,
        px.len(),
        n_diff as f64 / px.len() as f64
    );

    println!("== STRAIGHTNESS BY FIRST-DIVERGENCE BOUNDARY ==");
    println!(
        "  Structure-tensor anisotropy over a 9x9 window: 1 is a line, 0 is isotropic scatter.\n\
         **If straightness encodes WHEN, this falls with k.**\n"
    );
    println!("  {:>12} {:>10} {:>14} {:>14}", "k", "pixels", "anisotropy p50", "anisotropy p90");
    let bands: [(&str, std::ops::Range<usize>); 7] = [
        ("0", 0..1),
        ("1", 1..2),
        ("2", 2..3),
        ("3-4", 3..5),
        ("5-9", 5..10),
        ("10-24", 10..25),
        ("25+", 25..cfg.n_sync),
    ];
    for (name, r) in &bands {
        let m: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| r.contains(&k))).collect();
        let cnt = m.iter().filter(|x| **x).count();
        let mut a: Vec<f64> = (0..m.len())
            .filter(|&i| m[i])
            .filter_map(|i| straightness(&m, res, i, 4))
            .collect();
        if cnt > 0 {
            println!(
                "  {name:>12} {cnt:>10} {:>14.4} {:>14.4}",
                q(&mut a.clone(), 0.5),
                q(&mut a, 0.9)
            );
        }
        let buf: Vec<u8> =
            m.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect();
        let p = format!("{dir}/first_k_{}.png", name.replace('+', "plus").replace('-', "_"));
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p,
            &cfg,
            &format!(
                "res={res}x{res}\ncase=config_stability\nfirst-divergence boundary k in {name}\n\
                 POLARITY: white = differs from a neighbour first at this k\n"
            ),
        );
    }

    // The whole thing as one field, k on a ramp, so the layering is visible at once.
    let buf: Vec<u8> = first
        .iter()
        .flat_map(|f| match f {
            None => [12u8, 12, 16],
            Some(k) => {
                // log in k: the early boundaries are the interesting ones and a linear ramp
                // would compress them all into one colour.
                let t = ((*k as f64 + 1.0).ln() / ((cfg.n_sync as f64).ln())).clamp(0.0, 1.0);
                [
                    (255.0 * (1.0 - t)) as u8,
                    (120.0 + 100.0 * (1.0 - t) * t * 4.0) as u8,
                    (255.0 * t) as u8,
                ]
            }
        })
        .collect();
    let p = format!("{dir}/first_k_field.png");
    let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
    let _ = prin_rs::output::provenance_sidecar(
        &p,
        &cfg,
        &format!("res={res}x{res}\nfirst-divergence boundary k, LOG ramp\n\
                  RED = early (k small), BLUE = late, dark = never differs\n"),
    );

    println!(
        "\nWrote {dir}/ -- one mask per k band, plus first_k_field.png (RED early, BLUE late).\n\n\
         **If the early bands are clean lines and the late ones are tangles, the straightness of\n\
         a black line is a clock**: it says how many sync boundaries the two neighbouring\n\
         trajectories stayed together before one of them changed which body was farthest."
    );
}
