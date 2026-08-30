//! **Is the wedge structure a property of the TRAJECTORIES, or of Aarseth-Zare's machinery?**
//!
//! Every test so far ran inside AZ, so every one of them could only compare AZ's own quantities
//! against each other. The one control that can separate "sensitivity structure" from "chart
//! machinery" is an **independent integrator on the same initial conditions**.
//!
//! `integrate::leapfrog` is exactly that: unregularised KDK, no reference body, no argmax, no
//! chart switching, no fictitious time. If the wedges appear in ITS drift field they are a
//! property of the trajectories; if they do not, they are AZ's.
//!
//! # The argument this is answering
//!
//! Against the sensitivity reading: **the predictive step limit reduced the error in exactly
//! these regions** — `err>10` 0.1114 -> 0.0001, overshoot 634 -> 0. Structure that a step-control
//! change can attenuate is not obviously physics. The counter-reading is that sensitivity sets
//! *where* error concentrates while stepping sets *how much*, so both can be true. Neither can be
//! settled from inside one integrator.
//!
//! # What would make this control useless, stated first
//!
//! Leapfrog is **expected to fail** on close encounters — that failure is why AZ exists — so its
//! drift field will be dominated by them and its budget exhaustion will be large. Two guards:
//! `frac_ok` is printed, and every correlation is computed **on the pixels where leapfrog
//! actually completed**. A correlation over a field that mostly failed is a correlation with the
//! failure pattern.
//!
//! # And the FTLE arm
//!
//! If the wedges are early dynamical sensitivity, early switching-line density should track a
//! genuine sensitivity measure. FTLE is computed on the nominal copy and correlated against both
//! the density and the AZ drift. **The shifted control is carried through**: a correlation that
//! survives displacing one field by half a frame is about the marginals, not the alignment.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::leapfrog;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::ftle::{self, FtleOpts};

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
    if a.len() < 3 {
        return f64::NAN;
    }
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

/// Spearman over a selected subset, with the shifted control alongside.
fn rho_pair(a: &[f64], b: &[f64], sel: &[usize], res: usize) -> (f64, f64) {
    let (x, y): (Vec<f64>, Vec<f64>) = sel.iter().map(|&i| (a[i], b[i])).unzip();
    let straight = pearson(&ranks(&x), &ranks(&y));
    let shifted: Vec<f64> = sel
        .iter()
        .map(|&i| {
            let (px, py) = (i % res, i / res);
            a[((py + res / 2) % res) * res + (px + res / 2) % res]
        })
        .collect();
    (straight, pearson(&ranks(&shifted), &ranks(&y)))
}

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
    let res: usize = arg(1, 192);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/independent");
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

    // --- AZ ------------------------------------------------------------------------------
    let t0 = std::time::Instant::now();
    let az: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    println!("AZ pass {:.1}s", t0.elapsed().as_secs_f64());
    let n = az.len();

    // --- leapfrog, the independent arm -----------------------------------------------------
    let t1 = std::time::Instant::now();
    let lf: Vec<(f64, bool)> = (0..n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = sl.decode_pos(k);
            let st = grid::decode_state(&chart, 0, x, y);
            let o = leapfrog::integrate::<f64>(
                st.s, &st.m, cfg.t_max, cfg.eta, 0.0, 2_000_000,
            );
            (o.drift, o.finite && o.reached(cfg.t_max))
        })
        .collect();
    let ok: Vec<usize> = (0..n).filter(|&i| lf[i].1 && lf[i].0.is_finite()).collect();
    println!(
        "leapfrog pass {:.1}s -- frac_ok {:.4}. **Leapfrog is EXPECTED to fail on close\n\
         encounters; that failure is why AZ exists. Every correlation below is over the pixels\n\
         where it completed, because a correlation over a mostly-failed field is a correlation\n\
         with the failure pattern.**\n",
        t1.elapsed().as_secs_f64(),
        ok.len() as f64 / n as f64
    );

    // --- FTLE, the sensitivity arm ---------------------------------------------------------
    let t2 = std::time::Instant::now();
    let fo = FtleOpts::default();
    let pert = ftle::unit_perturbation::<f64>(0);
    let ftl: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = sl.decode_pos(k);
            let st = grid::decode_state(&chart, 0, x, y);
            let o = ftle::integrate_full::<f64>(st.s, &st.m, cfg.t_max, cfg.ftle_dt, &fo, &pert);
            if o.n_renorm > 0 { o.ftle } else { f64::NAN }
        })
        .collect();
    println!("FTLE pass {:.1}s\n", t2.elapsed().as_secs_f64());

    // --- the fields --------------------------------------------------------------------
    let lg = |x: f64| if x.is_finite() && x > 0.0 { x.log10() } else { DLO.log10() };
    let d_az: Vec<f64> = az.iter().map(|p| lg(p.energy_drift_max)).collect();
    let d_lf: Vec<f64> = lf.iter().map(|o| lg(o.0)).collect();
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
                let (a, c) = (&az[i].ref_path, &az[j].ref_path);
                if let Some(k) = (0..a.len().min(c.len())).find(|&k| a[k] != c[k]) {
                    b = Some(b.map_or(k, |p: usize| p.min(k)));
                }
            }
            b
        })
        .collect();
    let early: Vec<bool> = first.iter().map(|f| f.map_or(false, |k| k <= 9)).collect();
    let dens = density(&early, res, 15);
    let all: Vec<usize> = (0..n).collect();

    println!("== SPEARMAN, with the half-frame shifted control on every row ==");
    println!(
        "  **A correlation that survives the shift is about the two marginals, not about where\n\
         the fields are.** `n` is the population; the leapfrog rows use only completed pixels.\n"
    );
    println!("  {:>34} {:>9} {:>11} {:>11}", "pair", "n", "spearman", "shifted");
    for (name, a, b, sel) in [
        ("early-line density vs AZ drift", &dens, &d_az, &all),
        ("early-line density vs FTLE", &dens, &ftl, &all),
        ("FTLE vs AZ drift", &ftl, &d_az, &all),
        ("AZ drift vs LEAPFROG drift", &d_az, &d_lf, &ok),
        ("early-line density vs LEAPFROG drift", &dens, &d_lf, &ok),
        ("FTLE vs LEAPFROG drift", &ftl, &d_lf, &ok),
    ] {
        let (r, s) = rho_pair(a, b, sel, res);
        println!("  {name:>34} {:>9} {r:>11.4} {s:>11.4}", sel.len());
    }

    // --- panels ---------------------------------------------------------------------------
    let save = |nm: &str, buf: Vec<u8>, note: &str| {
        let p = format!("{dir}/{nm}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p, &cfg, &format!("res={res}x{res}\ncase=config_stability\n{note}\n"),
        );
    };
    save(
        "drift_az",
        az.iter()
            .flat_map(|p| {
                let x = p.energy_drift_max;
                if x.is_finite() && x > 0.0 {
                    ramp((lg(x) - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    [255, 0, 255]
                }
            })
            .collect(),
        &format!("Aarseth-Zare drift, FIXED ramp ({DLO:e},{DHI:e})"),
    );
    save(
        "drift_leapfrog",
        lf.iter()
            .flat_map(|o| {
                if o.1 && o.0.is_finite() {
                    ramp((lg(o.0) - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else {
                    // Magenta = leapfrog did not complete. **Expected**, and it must be visible:
                    // a failed pixel painted like a low-drift one would fake agreement.
                    [255, 0, 255]
                }
            })
            .collect(),
        &format!("UNREGULARISED leapfrog drift, same FIXED ramp. MAGENTA = did not complete."),
    );
    let mut fs: Vec<f64> = ftl.iter().cloned().filter(|x| x.is_finite()).collect();
    let (flo, fhi) = (q(&mut fs.clone(), 0.02), q(&mut fs, 0.98));
    save(
        "ftle",
        ftl.iter()
            .flat_map(|&x| {
                if x.is_finite() { ramp((x - flo) / (fhi - flo).max(1e-12)) } else { [255, 0, 255] }
            })
            .collect(),
        &format!("FTLE, ramp p2-p98 ({flo:.4},{fhi:.4})"),
    );
    println!("\nWrote {dir}/ -- drift_az, drift_leapfrog, ftle.");
}
