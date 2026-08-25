//! **§3.2 — does balanced mode hold a steady state, or drift to uniform depth?**
//!
//! March the playhead and plot depth variance across leaves against `t`.
//!
//! ```text
//! degenerating : variance -> 0   (everything converges to max depth)
//! balanced     : variance roughly constant, while individual quads still move
//! ```
//!
//! # Both curves, because one of them cannot tell the two apart
//!
//! **A tree that is stable because nothing moves is frozen, not balanced, and the two are
//! identical in a variance plot alone.** So churn is reported beside it: the fraction of quads
//! present at both playheads whose `Decision` changed. Frozen reads variance-flat and
//! churn-zero; balanced reads variance-flat and churn-nonzero.
//!
//! # The control
//!
//! `Mode::Uniform` runs in the same figure. It has the criterion **off** and must degenerate --
//! that is what proves the test discriminates rather than merely producing a line. If the uniform
//! arm does *not* degenerate, it is budget-bound rather than veto-bound and is no control at all;
//! `budget_exhausted` is printed per row so that is visible rather than assumed.
//!
//! # `n_sync` scales with `t_max`
//!
//! `dtau = eta*dt_left/(A0*B0)`, so holding `n_sync` fixed while `t_max` moves changes the step
//! size and the rows become different discretisations rather than one trajectory at several
//! playheads. Scaled here, as `criterion_metric` already does.

use std::collections::HashMap;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::output::plot::{Figure, Series};
use prin_rs::quad::Decision;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Mode, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// A quad's identity across playheads: its box, quantised. Two runs at different `t` build
/// different trees, so a node index means nothing between them.
fn key(cx: f64, cy: f64, half: f64) -> (i64, i64, i64) {
    let q = |x: f64| (x / (half * 1e-3)).round() as i64;
    (q(cx), q(cy), (half.log2() * 1e6).round() as i64)
}

fn main() {
    let budget: usize = arg(1, 800);
    let n: usize = arg(2, 4);
    let tau: f64 = arg(3, 1e-4);
    let viewport: usize = arg(4, 64);
    let ts = [4.0f64, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0];
    let base = EnsembleCfg::default();

    println!("balanced march. budget {budget}, N={n}, E+1=8, tau={tau:.0e}, viewport {viewport}²");
    println!("n_sync scales with t_max, or the rows are different discretisations.");
    println!();

    for region in ["near-field", "deep interior"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        println!("=== {region} ===");
        println!("{:>10}{:>6}{:>8}{:>10}{:>9}{:>9}{:>7}{:>7}{:>8}{:>12}{:>9}",
                 "mode", "t", "leaves", "depth var", "churn", "shared", "floor", "keep",
                 "screen", "spread med", "wall s");

        let mut fig_var: Vec<Series> = Vec::new();
        let mut fig_churn: Vec<Series> = Vec::new();

        for mode in [Mode::Balanced, Mode::Uniform] {
            let mut prev: Option<HashMap<(i64, i64, i64), Decision>> = None;
            let (mut vs, mut cs) = (Vec::new(), Vec::new());

            for &t_max in &ts {
                let n_sync =
                    ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
                let ens = EnsembleCfg {
                    refine_flagged: false,
                    t_max,
                    n_sync,
                    ..Default::default()
                };
                let cam = Camera::framing(root.cx, root.cy, 0.05, viewport);
                let cfg = SchedCfg {
                    n,
                    budget,
                    tau_display: tau,
                    alpha_hi: 0.2,
                    alpha_lo: 0.2,
                    mode,
                    camera: Some(cam),
                    ..Default::default()
                };
                let (tree, st) =
                    scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);

                let lv: Vec<f64> = tree.leaves().map(|i| tree.nodes[i].level as f64).collect();
                let m = lv.iter().sum::<f64>() / lv.len().max(1) as f64;
                let var =
                    lv.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / lv.len().max(1) as f64;

                let now: HashMap<(i64, i64, i64), Decision> = tree
                    .leaves()
                    .map(|i| {
                        let q = &tree.nodes[i];
                        (key(q.cx, q.cy, q.half), q.decision)
                    })
                    .collect();

                // Churn over the SHARED quads only. A quad that exists at one playhead and not
                // the other has not "changed decision"; counting it would fold the tree's size
                // change into a statistic about its stability, and the two move for different
                // reasons. `shared` is printed so a churn read over three quads is visible.
                let (churn, shared) = match &prev {
                    None => (f64::NAN, 0usize),
                    Some(p) => {
                        let common: Vec<_> =
                            now.keys().filter(|k| p.contains_key(*k)).collect();
                        let changed =
                            common.iter().filter(|k| p[**k] != now[**k]).count();
                        (
                            if common.is_empty() {
                                f64::NAN
                            } else {
                                changed as f64 / common.len() as f64
                            },
                            common.len(),
                        )
                    }
                };

                let c = |d: Decision| now.values().filter(|&&x| x == d).count();
                // **The premise being tested, measured rather than assumed.** Section 3 argues
                // that spread rises with `t` everywhere, so any fixed threshold must eventually
                // fire on every quad and balanced mode must degenerate to uniform. Print the
                // median leaf spread beside `keep` and the direction is a number, not an
                // inference from a decision count.
                let mut sp: Vec<f64> = tree
                    .leaves()
                    .map(|i| tree.nodes[i].red.spread_median)
                    .filter(|x| x.is_finite())
                    .collect();
                let spmed = prin_rs::stats::quantile(&mut sp, 0.5);
                println!("{:>10}{:>6.0}{:>8}{:>10.4}{:>9.4}{:>9}{:>7}{:>7}{:>8}{:>12.3e}{:>9.1}",
                         mode.name(), t_max, lv.len(), var, churn, shared,
                         c(Decision::Floor), c(Decision::Keep), c(Decision::ScreenFloor),
                         spmed, st.wall_seconds);
                if st.budget_exhausted {
                    println!("{:>10}       ^ BUDGET EXHAUSTED -- this row is budget-bound, not \
                              criterion-bound, and the uniform arm is no control if it fires.", "");
                }
                if shared > 0 && shared < 20 {
                    println!("{:>10}       ^ churn read over only {shared} shared quads. The tree \
                              changed size; this is thin.", "");
                }

                vs.push((t_max, var));
                if churn.is_finite() {
                    cs.push((t_max, churn));
                }
                prev = Some(now);
            }

            let rgb = if mode == Mode::Uniform { (220, 120, 90) } else { (120, 170, 230) };
            let mut sv = Series::new(mode.name().to_string(), vs, rgb);
            let mut sc = Series::new(mode.name().to_string(), cs, rgb);
            if mode == Mode::Uniform {
                sv = sv.dashed();
                sc = sc.dashed();
            }
            fig_var.push(sv);
            fig_churn.push(sc);
        }
        println!();

        let stem = region.replace(' ', "_");
        Figure {
            title: format!("depth variance against t -- {region}"),
            x_label: "t_max".into(),
            y_label: "depth variance across leaves".into(),
            series: fig_var,
            notes: vec![
                "uniform (dashed) is the CONTROL and must degenerate: criterion off, split to \
                 the veto."
                    .into(),
                "A flat variance curve alone cannot distinguish balanced from FROZEN. Read the \
                 churn figure beside it."
                    .into(),
            ],
            // Log panel with a zero band below it: the uniform control's exact 0 lands in the
            // band rather than being dropped, which is the whole point of showing it.
            y_lo: 1e-3,
            y_hi: 1.0,
        }
        .save(&format!("results/criterion/march_var_{stem}"))
        .unwrap();

        Figure {
            title: format!("per-quad churn against t -- {region}"),
            x_label: "t_max".into(),
            y_label: "fraction of shared quads whose decision changed".into(),
            series: fig_churn,
            notes: vec![
                "Churn is over quads present at BOTH playheads. A quad that appeared or vanished \
                 has not changed its decision."
                    .into(),
                "Balanced is a steady state, not a frozen one: variance flat AND churn nonzero."
                    .into(),
            ],
            y_lo: 1e-3,
            y_hi: 1.0,
        }
        .save(&format!("results/criterion/march_churn_{stem}"))
        .unwrap();
    }

    println!("Read the two figures together. Variance -> 0 is degeneration; variance flat with");
    println!("churn ~ 0 is a FROZEN tree, which a variance plot alone reports as success.");
}
