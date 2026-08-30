//! **Are the wedges physics or numerics — and are they even in the field that ships?**
//!
//! Two questions, one harness, because they share a pass and the answer to either can make the
//! other moot.
//!
//! # 4. Is the structure in the SCIENCE field, or only in the error field?
//!
//! Everything so far is measured on `energy_drift_max`, which is an **error** field: it is
//! supposed to show where the hard regions are, and wedges in it may be the instrument working
//! rather than failing. The field the shipping colouring actually reads is `spread_shape`. This
//! computes the same inside-cell-versus-across-boundary gradient ratio on **both**, plus `t_end`
//! and `ensemble_spread`, so "the wedges are in the drift field" and "the wedges are in the
//! image" stop being the same claim.
//!
//! Precedent, and it points both ways: on this same slice the nested arcs vanished entirely under
//! outcome-class colouring while the crisp edges survived and sharpened.
//!
//! # 2. Do the cells CONVERGE under refinement?
//!
//! The reference sequence is a **symbolic itinerary** — at each of `n_sync` times, which body is
//! the odd one out. Regions of constant itinerary are cylinder sets of a coarse symbolic dynamics
//! on the three-body problem. If the cell boundaries hold still as `eta` falls they are a property
//! of the *system*, and removing them would be removing physics. If they wander, they are a
//! property of the *discretisation*.
//!
//! This is the project's own refinement signature applied to a **symbolic** structure rather than
//! to a number: *a quantity that does not converge under refinement is measuring the sampling
//! rather than the system* — read here in the affirmative.
//!
//! **The control that makes it readable:** the itinerary is compared to the FINEST rung, not to
//! the previous one. A chain of pairwise comparisons can shrink while the sequence walks steadily
//! away, which would read as convergence and be the opposite.
//!
//! # Resumable
//!
//! Every rung is checkpointed **per row**, so a kill costs at most one row rather than the whole
//! pass — the `eta/16` arm is ~80 minutes at 256^2 and was lost three times as one indivisible
//! block. Re-running the same command resumes; the checkpoint is keyed on the full config
//! provenance plus resolution, and **refuses to resume under different settings** rather than
//! mixing two experiments.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::ckpt::Ckpt;

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

/// The cell mask: true where a pixel's reference itinerary differs from a neighbour's.
fn cells(px: &[PixelOut], res: usize) -> Vec<bool> {
    (0..px.len())
        .map(|i| {
            let (x, y) = (i % res, i / res);
            let d = |j: usize| px[i].ref_path != px[j].ref_path;
            (x + 1 < res && d(i + 1)) || (y + 1 < res && d(i + res))
        })
        .collect()
}

fn lg(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 {
        x.log10()
    } else {
        -12.0
    }
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

    // The checkpoint key. **Everything that changes a number goes in it**: a key that omits a
    // swept parameter is the `criterion_sweep` filename bug at a new site, where six rows landed
    // on one stem and the last writer won.
    let key = format!("wedge_origin res={res} {}", cfg.provenance());
    let ck_path = format!("{root}/output/wedge_origin_{res}.ckpt");
    let ck = std::sync::Mutex::new(match Ckpt::open(&ck_path, &key) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("checkpoint: {e}");
            std::process::exit(1);
        }
    });
    {
        let n = ck.lock().unwrap().1.len();
        println!("checkpoint {ck_path}: {n} rows already done\n");
    }

    // Per-pixel record: the few scalars the analysis reads, then the reference itinerary.
    // `f64` first at a fixed offset so the byte layout needs no framing.
    let n_sync = cfg.n_sync;
    let enc = |p: &PixelOut| -> Vec<u8> {
        let mut v = Vec::with_capacity(4 * 8 + n_sync);
        for x in [p.energy_drift_max, p.spread_shape, p.ensemble_spread, p.t_end] {
            v.extend_from_slice(&x.to_le_bytes());
        }
        let mut path = p.ref_path.clone();
        path.resize(n_sync, 255); // 255 marks "this trajectory ended early", not a body index
        v.extend_from_slice(&path);
        v
    };
    let dec = |b: &[u8]| -> PixelOut {
        let g = |k: usize| f64::from_le_bytes(b[k * 8..k * 8 + 8].try_into().unwrap());
        PixelOut {
            energy_drift_max: g(0),
            spread_shape: g(1),
            ensemble_spread: g(2),
            t_end: g(3),
            ref_path: b[32..32 + n_sync].to_vec(),
            ..Default::default()
        }
    };

    // Rung index `r` and row `y` share one id space, so one file covers the whole experiment.
    let run_rung = |r: u64, eta: f64| -> Vec<PixelOut> {
        let mut out: Vec<Option<PixelOut>> = vec![None; sl.npix()];
        {
            let done = &ck.lock().unwrap().1;
            for y in 0..res {
                if let Some(b) = done.get(&(r * res as u64 + y as u64)) {
                    for x in 0..res {
                        let w = b.len() / res;
                        out[y * res + x] = Some(dec(&b[x * w..(x + 1) * w]));
                    }
                }
            }
        }
        let todo: Vec<usize> = (0..res).filter(|&y| out[y * res].is_none()).collect();
        // Rows are computed in parallel and written as they land; the checkpoint is append-only
        // and keyed by id, so out-of-order completion is not a hazard.
        let rows: Vec<(usize, Vec<PixelOut>)> = todo
            .par_iter()
            .map(|&y| {
                let v: Vec<PixelOut> = (0..res)
                    .map(|x| pixel::evaluate_at::<f64>(&sl, y * res + x, &cfg, eta))
                    .collect();
                (y, v)
            })
            .collect();
        for (y, v) in rows {
            let mut bytes = Vec::new();
            for p in &v {
                bytes.extend_from_slice(&enc(p));
            }
            let _ = ck.lock().unwrap().0.put(r * res as u64 + y as u64, &bytes);
            for (x, p) in v.into_iter().enumerate() {
                out[y * res + x] = Some(p);
            }
        }
        out.into_iter().map(|o| o.unwrap()).collect()
    };
    let run = |eta: f64| -> Vec<PixelOut> {
        let r = match eta {
            e if (e - cfg.eta).abs() < 1e-18 => 0u64,
            e if (e - cfg.eta / 4.0).abs() < 1e-18 => 1,
            _ => 2,
        };
        run_rung(r, eta)
    };

    // --- 4. which fields carry the structure? ---------------------------------------------
    let t0 = std::time::Instant::now();
    let px = run(cfg.eta);
    println!("nominal pass {:.1}s\n", t0.elapsed().as_secs_f64());
    let cell = cells(&px, res);
    let cell_frac = cell.iter().filter(|x| **x).count() as f64 / cell.len() as f64;

    println!("== 4. IS THE STRUCTURE IN THE FIELD THAT SHIPS? ==");
    println!(
        "  Ratio of median |grad log10 field| ACROSS a cell boundary to INSIDE a cell. A large\n\
         ratio means the reference partition draws that field's edges; ~1 means it does not.\n\
         `drift` is the error field and is the comparison, not the subject.\n\n\
         Cell boundaries cover {cell_frac:.4} of the frame.\n"
    );
    println!(
        "  {:>18} {:>13} {:>13} {:>8}",
        "field", "inside cell", "across bound", "ratio"
    );
    let fields: [(&str, Vec<f64>); 4] = [
        ("energy_drift_max", px.iter().map(|p| lg(p.energy_drift_max)).collect()),
        ("spread_shape", px.iter().map(|p| lg(p.spread_shape)).collect()),
        ("ensemble_spread", px.iter().map(|p| lg(p.ensemble_spread)).collect()),
        ("t_end", px.iter().map(|p| p.t_end).collect()),
    ];
    for (name, f) in &fields {
        let g = grad(f, res);
        let mut gi: Vec<f64> = (0..g.len()).filter(|&i| !cell[i]).map(|i| g[i]).collect();
        let mut ga: Vec<f64> = (0..g.len()).filter(|&i| cell[i]).map(|i| g[i]).collect();
        let (a, b) = (q(&mut gi, 0.5), q(&mut ga, 0.5));
        println!("  {name:>18} {a:>13.5} {b:>13.5} {:>8.3}", b / a.max(f64::MIN_POSITIVE));
    }

    // --- 2. do the cells converge? ----------------------------------------------------------
    println!("\n== 2. DO THE CELLS CONVERGE UNDER REFINEMENT? ==");
    println!(
        "  The itinerary is compared to the FINEST rung, never to the previous one: a chain of\n\
         pairwise comparisons can shrink while the sequence walks steadily away, which reads as\n\
         convergence and is its opposite.\n"
    );
    let fine_eta = cfg.eta / 16.0;
    let t1 = std::time::Instant::now();
    let fine = run(fine_eta);
    println!("  reference pass at eta/16 in {:.1}s\n", t1.elapsed().as_secs_f64());
    let cell_fine = cells(&fine, res);

    println!(
        "  {:>10} {:>14} {:>14} {:>14} {:>14}",
        "eta", "itin identical", "hamming p50", "cell agree", "cell frac"
    );
    for div in [1.0f64, 4.0, 16.0] {
        let v = if div == 16.0 { fine.clone() } else { run(cfg.eta / div) };
        let mut ham: Vec<f64> = (0..v.len())
            .map(|i| {
                v[i].ref_path
                    .iter()
                    .zip(fine[i].ref_path.iter())
                    .filter(|(a, b)| a != b)
                    .count() as f64
            })
            .collect();
        let ident = ham.iter().filter(|x| **x == 0.0).count() as f64 / ham.len() as f64;
        let c = cells(&v, res);
        let agree = (0..c.len()).filter(|&i| c[i] == cell_fine[i]).count() as f64 / c.len() as f64;
        println!(
            "  {:>10.3e} {ident:>14.4} {:>14.1} {agree:>14.4} {:>14.4}",
            cfg.eta / div,
            q(&mut ham, 0.5),
            c.iter().filter(|x| **x).count() as f64 / c.len() as f64
        );
    }

    println!(
        "\n  **`cell agree` is the number that decides.** If the cell mask is essentially the same\n\
         at every rung, the wedges are a property of the SYSTEM -- cylinder sets of the symbolic\n\
         itinerary -- and there is nothing to remove. If it moves with `eta`, they are a property\n\
         of the DISCRETISATION and event-driven registration is worth building.\n\n\
         `itin identical` will fall with `t_max` on any chaotic slice however real the cells are:\n\
         125 boundaries is 125 chances to differ, and one flip late in a trajectory does not move\n\
         a cell BOUNDARY. Read the mask agreement, not the sequence equality."
    );

    // --- panels -----------------------------------------------------------------------------
    let mask = |m: &[bool]| -> Vec<u8> {
        m.iter().flat_map(|&x| if x { [255u8, 255, 255] } else { [12, 12, 16] }).collect()
    };
    for (n, m) in [("cells_eta", &cell), ("cells_eta_over_16", &cell_fine)] {
        let p = format!("{dir}/{n}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, &mask(m));
        let _ = prin_rs::output::provenance_sidecar(&p, &cfg, &format!("res={res}x{res}\n"));
    }
    println!("\nWrote {dir}/cells_eta.png and cells_eta_over_16.png -- the two masks, to compare by eye\nas well as by the number.");
}
