//! **Is the wedge structure where the early switching surfaces are DENSE?**
//!
//! Every previous test asked *membership*: is this pixel on a switching line. The lifts were mild
//! — 1.825 at best — and the lines looked transverse to the wedges. But the eye does not read
//! membership, it reads **density**: "there is a lot of red here". Those are different statistics
//! and only the first was ever measured.
//!
//! This is a hypothesis that came from looking at the pictures, not from the numbers, and it is
//! the one statistic in this whole investigation that was never tried.
//!
//! # What is measured
//!
//! `early_density[i]` — the fraction of pixels within a window whose reference itinerary first
//! differs from a neighbour's at `k <= K`. Then:
//!
//! - Spearman correlation of `early_density` against `log10 drift`, on ranks so a heavy-tailed
//!   field cannot dominate it.
//! - Median drift by density decile: **a monotone rise is the claim**, and a flat table is its
//!   refutation.
//! - The same for late surfaces (`k > K`) as a control — if BOTH correlate, the statistic is
//!   reading overall switching activity rather than anything about early surfaces.
//! - And a **shifted control**: the density field displaced by half a frame, which destroys any
//!   real spatial relationship while preserving the field's own distribution. A correlation that
//!   survives the shift is an artefact of the two fields' marginals, not of their alignment.
//!
//! Window radius is swept, because "dense" has a scale and picking one would be choosing the
//! answer.
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

/// Rank transform, ties given their average rank.
fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0;
        for &k in &idx[i..=j] {
            r[k] = avg;
        }
        i = j + 1;
    }
    r
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..a.len() {
        let (x, y) = (a[i] - ma, b[i] - mb);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    num / (da.sqrt() * db.sqrt()).max(f64::MIN_POSITIVE)
}

/// Box-filter density of a boolean set.
fn density(m: &[bool], res: usize, r: i64) -> Vec<f64> {
    (0..m.len())
        .into_par_iter()
        .map(|i| {
            let (cx, cy) = ((i % res) as i64, (i / res) as i64);
            let (mut tot, mut hit) = (0.0f64, 0.0f64);
            for dy in -r..=r {
                for dx in -r..=r {
                    let (x, y) = (cx + dx, cy + dy);
                    if x < 0 || y < 0 || x >= res as i64 || y >= res as i64 {
                        continue;
                    }
                    tot += 1.0;
                    if m[y as usize * res + x as usize] {
                        hit += 1.0;
                    }
                }
            }
            hit / tot
        })
        .collect()
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/line_density");
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

    let lg: Vec<f64> = px
        .iter()
        .map(|p| {
            let x = p.energy_drift_max;
            if x.is_finite() && x > 0.0 { x.log10() } else { DLO.log10() }
        })
        .collect();
    let rd = ranks(&lg);

    const K: usize = 9;
    let early: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k <= K)).collect();
    let late: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k > K)).collect();

    println!(
        "== SPEARMAN: local density of switching surfaces vs log drift ==\n\
         Ranks, so a heavy-tailed field cannot dominate. `shifted` displaces the density field by\n\
         half a frame: it preserves the marginal and destroys the alignment, so **a correlation\n\
         that survives it is an artefact of the two distributions rather than of where they are.**\n"
    );
    println!(
        "  {:>8} {:>16} {:>12} {:>12} {:>12}",
        "radius", "set", "spearman", "shifted", "density p50"
    );
    let mut best: Option<(i64, Vec<f64>)> = None;
    for r in [3i64, 7, 15, 31] {
        for (name, m) in [("early k<=9", &early), ("late k>9 (ctrl)", &late)] {
            let d = density(m, res, r);
            let rho = pearson(&ranks(&d), &rd);
            // Half-frame shift, wrapped.
            let sh: Vec<f64> = (0..n)
                .map(|i| {
                    let (x, y) = (i % res, i / res);
                    d[((y + res / 2) % res) * res + (x + res / 2) % res]
                })
                .collect();
            let rho_s = pearson(&ranks(&sh), &rd);
            println!(
                "  {r:>8} {name:>16} {rho:>12.4} {rho_s:>12.4} {:>12.4}",
                q(&mut d.clone(), 0.5)
            );
            if name.starts_with("early") && r == 15 {
                best = Some((r, d));
            }
        }
    }

    // The decile table: a monotone rise is the claim, a flat table refutes it.
    if let Some((r, d)) = &best {
        println!("\n== DRIFT BY EARLY-DENSITY DECILE (radius {r}) ==");
        println!(
            "  **A monotone rise is the claim.** A flat table means density explains nothing,\n\
             whatever the correlation says.\n"
        );
        println!("  {:>8} {:>14} {:>13} {:>13}", "decile", "density range", "drift p50", "n");
        let mut s = d.clone();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for k in 0..10 {
            let lo = s[(s.len() - 1) * k / 10];
            let hi = s[(s.len() - 1) * (k + 1) / 10];
            let sel: Vec<usize> =
                (0..n).filter(|&i| d[i] >= lo && (k == 9 || d[i] < hi)).collect();
            let mut dr: Vec<f64> = sel
                .iter()
                .map(|&i| px[i].energy_drift_max)
                .filter(|x| x.is_finite())
                .collect();
            println!(
                "  {:>8} {:>6.3}-{:>6.3} {:>13.3e} {:>13}",
                k + 1,
                lo,
                hi,
                q(&mut dr, 0.5),
                sel.len()
            );
        }

        // Panels: the density field beside the drift field, same size, both fixed-ramped.
        let ramp = |x: f64| -> [u8; 3] {
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
        let dmax = q(&mut d.clone(), 0.99).max(1e-9);
        for (nm, buf, note) in [
            (
                "early_density",
                d.iter().flat_map(|&x| ramp(x / dmax)).collect::<Vec<u8>>(),
                format!("early (k<=9) switching-line density, radius {r}, ramp 0..{dmax:.4} (p99)"),
            ),
            (
                "drift",
                px.iter()
                    .flat_map(|p| {
                        let x = p.energy_drift_max;
                        if x.is_finite() && x > 0.0 {
                            ramp((x.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
                        } else {
                            [255, 0, 255]
                        }
                    })
                    .collect::<Vec<u8>>(),
                format!("energy_drift_max, FIXED ramp ({DLO:e},{DHI:e})"),
            ),
        ] {
            let p = format!("{dir}/{nm}.png");
            let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
            let _ = prin_rs::output::provenance_sidecar(
                &p, &cfg, &format!("res={res}x{res}\ncase=config_stability\n{note}\n"),
            );
        }
        println!("\nWrote {dir}/early_density.png and drift.png -- same size, both fixed-ramped.");
    }
}
