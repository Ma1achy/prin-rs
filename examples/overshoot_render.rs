//! Four-arm renders for the boundary-overshoot clamp, so the **interaction** with `dtau_mode`
//! is visible rather than inferred.
//!
//! ```text
//!   A  dtau fixed     + overshoot present    the original committed behaviour
//!   B  dtau per-step  + overshoot present    the regression
//!   C  dtau fixed     + overshoot clamped
//!   D  dtau per-step  + overshoot clamped    the proposed default
//! ```
//!
//! Two knobs, four cells. Rendering only the diagonal (A and D) would show a difference and say
//! nothing about which knob produced it -- and the claim being made here is specifically about
//! the *cross* terms: that the step-control fix alone (B) is visually worse than the behaviour
//! it replaced (A), because it turns a spatially smooth error into a spatially varying one.
//!
//! # Writes
//!
//! `<root>/<sub>/`, both arguments, defaulting to `results/overshoot_fix`. **Never** into
//! `results/dtau_fix`: those renders are the "before" for this comparison, exactly as the
//! committed ones were the "before" for that one.
//!
//! Args are `res root only sub arms`. `arms` is a subset of `ABCD`, so an interrupted run resumes
//! instead of repeating ~15 minutes per arm. The pair table below prints only when all four arms
//! are computed in **one** run; a resumed run writes its panels and says so by omission.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::DtauMode;
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// (name, chart, cx, cy, half, t_max, r_coll, all_four_arms)
type Target = (String, Chart, f64, f64, f64, f64, f64, bool);

const ARMS: [(&str, DtauMode, bool); 4] = [
    ("A", DtauMode::FixedPerInterval, false),
    ("B", DtauMode::PerStepInterval, false),
    ("C", DtauMode::FixedPerInterval, true),
    ("D", DtauMode::PerStepInterval, true),
];

/// Just enough of a pixel to compare two arms. Keeping four whole `PixelOut` fields at 1024^2
/// would be most of a gigabyte for two numbers.
#[derive(Clone, Copy)]
struct Slim {
    n: [f64; 3],
    outcome: u8,
}

fn main() {
    let res: usize = arg(1, 1024);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(3).unwrap_or_else(|| "all".into());
    let sub: String = std::env::args().nth(4).unwrap_or_else(|| "overshoot_fix".into());
    // **An argument, not a constant.** Each arm is ~15 minutes at 1024^2, and a run interrupted
    // after two of them should resume rather than repeat -- the pair comparisons are printed only
    // when all four are present in one run, which is a limitation of the run and is stated below.
    let arms_arg: String = std::env::args().nth(5).unwrap_or_else(|| "ABCD".into());
    let dir = format!("{root}/{sub}");
    let _ = std::fs::create_dir_all(&dir);

    let mut targets: Vec<Target> = Vec::new();
    let (cs, sx, sy, sh) = Chart::config_stability();
    targets.push(("config_stability".into(), cs, sx, sy, sh, 50.0, 0.005, true));
    for w in ["preset_shape", "preset_prho", "preset_plambda", "preset_shape_pl"] {
        if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == w) {
            targets.push((c.0.into(), c.1, c.2, c.3, c.4, 13.0, 1e-3, false));
        }
    }

    println!(
        "{res}^2, one sample per pixel. TWO knobs, FOUR cells -- `dtau_mode` and\n\
         `clamp_final_step`. Everything else is `closure_render`'s, the settings the committed\n\
         images were made under.\n\n\
         `nonfin` is the magenta count; `simfail` the budget-exhausted count, separate on\n\
         purpose. `hot` is the fraction above 1e-6 energy drift -- the speckle. `ramp` is each\n\
         panel's own auto-range window and belongs with the image.\n"
    );

    for (name, chart, cx, cy, half, t_max, r_coll, four) in targets {
        let name = name.as_str();
        if only != "all" && only != name {
            continue;
        }
        let arms: Vec<(&str, DtauMode, bool)> = if four {
            ARMS.iter().filter(|a| arms_arg.contains(a.0)).cloned().collect()
        } else {
            vec![ARMS[3]]
        };
        let mut slim: Vec<(String, Vec<Slim>)> = Vec::new();
        for (tag, mode, clamp) in arms {
            let ens = EnsembleCfg {
                refine_flagged: false,
                t_max,
                r_coll_frac: r_coll,
                dtau_mode: mode,
                clamp_final_step: clamp,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let nonfin = px.iter().filter(|p| p.n_nonfinite > 0).count();
            let simfail = px.iter().filter(|p| p.state == 6).count();
            let hot = px.iter().filter(|p| !(p.energy_drift_max <= 1e-6)).count();
            let esc = px.iter().filter(|p| p.state == 0).count();
            let col = px.iter().filter(|p| p.state == 2).count();

            let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
            let (dlo, dhi) = colour::range_q(&px, Scalar::Drift, 0.02, 0.98);
            let m_here = grid::decode_state(&chart, 0, cx, cy).m;
            let sites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(px.len() * 3);
            let mut obuf = Vec::with_capacity(px.len() * 3);
            let mut dbuf = Vec::with_capacity(px.len() * 3);
            for p in &px {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
                obuf.extend_from_slice(&png::outcome_rgb(p));
                dbuf.extend_from_slice(&colour::drift_rgb(p, dlo, dhi));
            }
            let stem = format!("{dir}/{name}_arm{tag}");
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &buf);
            let _ = adaptive::save_rect(&format!("{stem}_outcome.png"), res, res, &obuf);
            let _ = adaptive::save_rect(&format!("{stem}_drift.png"), res, res, &dbuf);

            println!(
                "{name:>20} arm {tag} {secs:>7.1}s  escape {:.4} collision {:.4}  \
                 nonfin {:>7}  simfail {:>7}  hot {:.4}  spread ramp ({lo:.3e}, {hi:.3e})  \
                 drift ramp ({dlo:.3e}, {dhi:.3e})",
                esc as f64 / px.len() as f64,
                col as f64 / px.len() as f64,
                nonfin,
                simfail,
                hot as f64 / px.len() as f64,
            );

            slim.push((
                tag.to_string(),
                px.iter().map(|p| Slim { n: p.shape_vec, outcome: p.outcome }).collect(),
            ));
        }

        // Per-pixel, never the aggregate alone. **The A->B figure is the convergence red flag**:
        // a converged integration does not move most of a field on a step-control change.
        if slim.len() == 4 {
            println!(
                "{:>20}  pair comparisons (chord on the shape sphere, diameter 2):",
                ""
            );
            for (i, j) in [(0usize, 1usize), (2, 3), (1, 3), (0, 2), (0, 3)] {
                let (a, b) = (&slim[i].1, &slim[j].1);
                let mut d: Vec<f64> = Vec::with_capacity(a.len());
                let (mut moved, mut lab) = (0usize, 0usize);
                for (x, y) in a.iter().zip(b.iter()) {
                    let s: f64 =
                        (0..3).map(|k| (x.n[k] - y.n[k]).powi(2)).sum::<f64>().sqrt();
                    if s > 0.0 {
                        moved += 1;
                    }
                    if x.outcome != y.outcome {
                        lab += 1;
                    }
                    d.push(s);
                }
                let mut f: Vec<f64> = d.iter().cloned().filter(|x| x.is_finite()).collect();
                let med = prin_rs::stats::quantile(&mut f.clone(), 0.5);
                println!(
                    "{:>20}    {}->{}  moved {moved:>8} ({:.4})  median {med:.3e}  max {:.3e}  \
                     labels flipped {lab}",
                    "",
                    slim[i].0,
                    slim[j].0,
                    moved as f64 / a.len() as f64,
                    f.iter().cloned().fold(0.0, f64::max)
                );
            }
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         THE PREDICTION, stated so it can fail: the stacked-crescent banding in the green and\n\
         blue mid-field of `_uniform` is PRESENT in arm B and GONE in arm D; arm C is smooth too\n\
         and probably smoother than A; the white regions do not grow further and ideally shrink\n\
         back toward A; magenta stays at arm B's level, near zero, because that was the `dtau`\n\
         fix and it stands; and drift improves in D over B, not merely over A.\n\n\
         If the crescents survive in arm D, that is a THIRD mechanism and it is worth knowing\n\
         immediately rather than explaining.\n\n\
         MEASURED, AND THE FIRST CLAUSE FAILED. The crescents are present in ALL FOUR arms,\n\
         including A, which predates both changes -- so neither knob draws them. Under\n\
         outcome-class colouring arm D's arcs vanish entirely while the region boundaries\n\
         sharpen, which is RESULTS §21's standing result at a new site: the banding is a\n\
         colouring artefact, the crisp edges are not. Every other clause held, and the magenta\n\
         ran 30109 -> 2071 (clamp alone) -> 178 (both). See RESULTS §24.8.\n"
    );
}
