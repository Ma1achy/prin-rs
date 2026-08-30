//! **The lattice of black lines cutting through the coherent ribbons — which surfaces are they?**
//!
//! # The measurement that was wrong, and why
//!
//! `switch_depth` scored straightness as structure-tensor anisotropy over a 9x9 window and
//! concluded the surfaces "fold by `k = 5`". **That metric is confounded by density.** A window
//! placed on a *lattice* of straight lines contains several lines in different directions and
//! reads isotropic — indistinguishable from a tangle. So a falling anisotropy may mean the set
//! got *denser*, not curvier, and the visible lattice was measured past entirely.
//!
//! # The metric that cannot be fooled that way
//!
//! **Connected components.** For each component of the differing set: its size, its extent, and
//! the RMS deviation from a total-least-squares line through it. `rms / extent` is a per-object
//! straightness that does not care how many other components share the neighbourhood. A lattice
//! of straight lines is many long components each with a small ratio; a tangle is few large
//! components with a large one.
//!
//! Components below a size floor are excluded and **counted**: a two-pixel component is perfectly
//! straight by arithmetic, and letting those into the statistic would report a lattice wherever
//! there is dust.
//!
//! # Polarity
//!
//! Here **white/coloured = the pixel's reference itinerary DIFFERS from a neighbour's** — the
//! black lines of `wedge_id/coherent_cells.png`, shown positively.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;
/// Components smaller than this are dust, not lines, and are excluded and counted.
const MIN_COMP: usize = 12;

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
    let dir = format!("{root}/step_control/line_id");
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

    let differs: Vec<bool> = (0..n)
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let d = |j: usize| px[i].ref_path != px[j].ref_path;
            (x + 1 < res && d(i + 1)) || (y + 1 < res && d(i + res))
        })
        .collect();
    let first: Vec<Option<usize>> = (0..n)
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let mut b: Option<usize> = None;
            for j in [
                if x + 1 < res { Some(i + 1) } else { None },
                if y + 1 < res { Some(i + res) } else { None },
            ]
            .into_iter()
            .flatten()
            {
                let (a, c) = (&px[i].ref_path, &px[j].ref_path);
                if let Some(k) = (0..a.len().min(c.len())).find(|&k| a[k] != c[k]) {
                    b = Some(b.map_or(k, |p: usize| p.min(k)));
                }
            }
            b
        })
        .collect();

    // **THE DIRECT CHECK.** The lattice a reader sees in `coherent_cells.png` is black facets
    // INSIDE white ribbons. If those are early switching surfaces, then drawing the early set on
    // the coherent mask must land exactly on them. Only the segments crossing a ribbon are
    // visible there -- the same lines continue through the dense regions, where they are
    // invisible against a background of the same colour, which is why they read as short facets
    // rather than as the long sweeps they are.
    for kcut in [4usize, 9, 19] {
        let early: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k <= kcut)).collect();
        let buf: Vec<u8> = (0..n)
            .flat_map(|i| {
                if early[i] {
                    [255u8, 40, 40]
                } else if differs[i] {
                    [12, 12, 16]
                } else {
                    [235, 235, 245]
                }
            })
            .collect();
        let p = format!("{dir}/coherent_with_k{kcut}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p, &cfg,
            &format!("res={res}x{res}\nPALE = coherent ribbon, DARK = differs late,\n\
                      RED = differs first at k <= {kcut}\n\
                      If the red lands on the dark facets inside the pale ribbons, the lattice\n\
                      is the early switching surfaces.\n"),
        );
        println!(
            "  wrote coherent_with_k{kcut}.png -- early set is {:.4} of the frame",
            early.iter().filter(|x| **x).count() as f64 / n as f64
        );
    }

    // **Isolate the thin lines from the dense tangle.** The whole differing set is one 8-connected
    // blob (measured: 1 component covering 67.8% of the frame), so components of it say nothing.
    // A lattice line inside a ribbon is a differing pixel whose neighbourhood is mostly COHERENT;
    // a pixel in the tangle is surrounded by other differing pixels. That distinction is what
    // separates the two populations, and without it the component statistic is a null.
    let thin: Vec<bool> = (0..n)
        .map(|i| {
            if !differs[i] {
                return false;
            }
            let (cxp, cyp) = ((i % res) as i64, (i / res) as i64);
            let (mut tot, mut coh) = (0.0f64, 0.0f64);
            for dy in -2i64..=2 {
                for dx in -2i64..=2 {
                    let (x, y) = (cxp + dx, cyp + dy);
                    if x < 0 || y < 0 || x >= res as i64 || y >= res as i64 {
                        continue;
                    }
                    tot += 1.0;
                    if !differs[y as usize * res + x as usize] {
                        coh += 1.0;
                    }
                }
            }
            coh / tot > 0.6
        })
        .collect();
    println!(
        "\n  thin lines (differing pixel with >60% coherent neighbours): {:.4} of the frame\n",
        thin.iter().filter(|x| **x).count() as f64 / n as f64
    );

    // 8-connected components of the THIN set, which is what the lattice actually is.
    let differs = thin;
    // 8-connected components of the differing set.
    let mut lab = vec![usize::MAX; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for s in 0..n {
        if !differs[s] || lab[s] != usize::MAX {
            continue;
        }
        let id = comps.len();
        let mut stack = vec![s];
        let mut cur = Vec::new();
        lab[s] = id;
        while let Some(i) = stack.pop() {
            cur.push(i);
            let (cxp, cyp) = ((i % res) as i64, (i / res) as i64);
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let (x, y) = (cxp + dx, cyp + dy);
                    if x < 0 || y < 0 || x >= res as i64 || y >= res as i64 {
                        continue;
                    }
                    let j = y as usize * res + x as usize;
                    if differs[j] && lab[j] == usize::MAX {
                        lab[j] = id;
                        stack.push(j);
                    }
                }
            }
        }
        comps.push(cur);
    }

    // Per component: total-least-squares line fit. `rms / extent` is straightness, and it does
    // not care how many other components share the neighbourhood.
    let stats: Vec<(usize, f64, f64, f64)> = comps
        .iter()
        .map(|c| {
            let nn = c.len() as f64;
            let (mut mx, mut my) = (0.0, 0.0);
            for &i in c {
                mx += (i % res) as f64;
                my += (i / res) as f64;
            }
            mx /= nn;
            my /= nn;
            let (mut a, mut b, mut cc) = (0.0f64, 0.0f64, 0.0f64);
            for &i in c {
                let (u, v) = ((i % res) as f64 - mx, (i / res) as f64 - my);
                a += u * u;
                b += u * v;
                cc += v * v;
            }
            let tr = a + cc;
            let disc = (tr * tr / 4.0 - (a * cc - b * b)).max(0.0).sqrt();
            let (l1, l2) = (tr / 2.0 + disc, tr / 2.0 - disc);
            // l1 is variance along the line, l2 across it.
            let extent = (l1 / nn).max(0.0).sqrt();
            let rms = (l2 / nn).max(0.0).sqrt();
            let med_k = {
                let mut v: Vec<f64> =
                    c.iter().filter_map(|&i| first[i]).map(|k| k as f64).collect();
                q(&mut v, 0.5)
            };
            (c.len(), extent, rms, med_k)
        })
        .collect();

    let big: Vec<&(usize, f64, f64, f64)> = stats.iter().filter(|s| s.0 >= MIN_COMP).collect();
    println!(
        "  {} components; {} with >= {MIN_COMP} pixels. **The {} below the floor are excluded and\n\
         counted: a two-pixel component is perfectly straight by arithmetic, and admitting dust\n\
         would report a lattice wherever there is noise.**\n",
        comps.len(),
        big.len(),
        comps.len() - big.len()
    );

    println!("== PER-COMPONENT STRAIGHTNESS (rms across / extent along; small = a line) ==");
    println!(
        "  {:>16} {:>9} {:>11} {:>11} {:>11} {:>11}",
        "component size", "count", "extent p50", "rms/ext p50", "rms/ext p90", "median k"
    );
    for (name, lo, hi) in [
        ("12-50", 12usize, 51usize),
        ("51-200", 51, 201),
        ("201-1000", 201, 1001),
        ("1000+", 1001, usize::MAX),
    ] {
        let g: Vec<&&(usize, f64, f64, f64)> =
            big.iter().filter(|s| s.0 >= lo && s.0 < hi).collect();
        if g.is_empty() {
            continue;
        }
        let mut r: Vec<f64> = g.iter().map(|s| s.2 / s.1.max(1e-12)).collect();
        let mut e: Vec<f64> = g.iter().map(|s| s.1).collect();
        let mut k: Vec<f64> = g.iter().map(|s| s.3).filter(|x| x.is_finite()).collect();
        println!(
            "  {name:>16} {:>9} {:>11.2} {:>11.4} {:>11.4} {:>11.1}",
            g.len(),
            q(&mut e, 0.5),
            q(&mut r.clone(), 0.5),
            q(&mut r, 0.9),
            q(&mut k, 0.5)
        );
    }

    // Panels. The differing set coloured by `first_k`, and the components coloured by
    // straightness, so which surfaces form the lattice is visible rather than inferred.
    let mut by_k: Vec<u8> = Vec::with_capacity(n * 3);
    let mut by_s: Vec<u8> = Vec::with_capacity(n * 3);
    let mut over: Vec<u8> = Vec::with_capacity(n * 3);
    for i in 0..n {
        let base = {
            let x = px[i].energy_drift_max;
            if x.is_finite() && x > 0.0 {
                ramp((x.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
            } else {
                [255, 0, 255]
            }
        };
        if !differs[i] {
            by_k.extend_from_slice(&[12, 12, 16]);
            by_s.extend_from_slice(&[12, 12, 16]);
            over.extend_from_slice(&base);
            continue;
        }
        let k = first[i].unwrap_or(0) as f64;
        let t = ((k + 1.0).ln() / (cfg.n_sync as f64).ln()).clamp(0.0, 1.0);
        by_k.extend_from_slice(&[(255.0 * (1.0 - t)) as u8, 60, (255.0 * t) as u8]);
        let s = stats[lab[i]];
        let straight = s.0 >= MIN_COMP && s.2 / s.1.max(1e-12) < 0.15;
        by_s.extend_from_slice(&if straight { [255, 255, 255] } else { [70, 70, 90] });
        over.extend_from_slice(&if straight { [0, 255, 255] } else { base });
    }
    for (nm, b) in [
        ("differs_by_k", by_k),
        ("differs_by_straightness", by_s),
        ("drift_with_straight_lines", over),
    ] {
        let p = format!("{dir}/{nm}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &b);
        let _ = prin_rs::output::provenance_sidecar(
            &p,
            &cfg,
            &format!(
                "res={res}x{res}\ncase=config_stability\npanel={nm}\n\
                 POLARITY: coloured = itinerary DIFFERS from a neighbour\n\
                 by_k: RED early, BLUE late (log ramp)\n\
                 straight = component >= {MIN_COMP} px with rms/extent < 0.15\n"
            ),
        );
    }
    println!("\nWrote {dir}/ -- differs_by_k, differs_by_straightness, drift_with_straight_lines.");
}
