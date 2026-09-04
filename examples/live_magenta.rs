//! **Where the magenta in a live frame comes from.**
//!
//! `colour::rgb` paints a footprint magenta ([`colour::DEBUG_NAN`]) for four reasons, and under a
//! live playhead only three of them are things a viewer at boundary `j` could know:
//!
//! 1. `n_nonfinite > 0` -- **a whole-run count**. `evaluate` sets it from the driver's usability
//!    flag over the march to `t_max`, so a copy that diverges at `t = 12` makes the footprint
//!    magenta in the frame at `t = 0.8`. That is the future, painted into the past.
//! 2. a non-finite nominal `shape_vec` at this boundary -- live, and correct.
//! 3. a non-finite scalar (here `spread_shape`) at this boundary -- live, and correct.
//! 4. `State::SimFailed` / `DecodeFailed` -- `project_at` overwrites the state of a run that has
//!    not terminated by `j`, so this one is already handled.
//!
//! This harness separates them: one `descend_live` per chart, and per boundary the count of
//! footprints each route would paint. **Route 1 constant from the first boundary is the
//! signature**; a genuine divergence count starts at zero and grows.
//!
//! Run: `cargo run --release --example live_magenta -- [charts=preset_shape_h1,...] [res=64]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::State;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn target(name: &str) -> Option<(String, Chart, f64, f64, f64, usize)> {
    if let Some(&(n, cx, cy, body)) = grid::REGIONS.iter().find(|r| r.0 == name) {
        return Some((n.into(), Chart::BodyPlane, cx, cy, 0.05, body));
    }
    if name == "config_stability" {
        let (chart, cx, cy, half) = Chart::config_stability();
        return Some((name.into(), chart, cx, cy, half, 0));
    }
    grid::gallery_cases().into_iter().find(|c| c.0 == name).map(|(n, c, cx, cy, h)| (n.into(), c, cx, cy, h, 0))
}

fn main() {
    let charts: Vec<String> = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "preset_shape_h1,config_stability".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let res: usize = arg(2, 64);
    let ens = EnsembleCfg { refine_flagged: false, keep_live_series: true, live_stride: 2, ..Default::default() };
    println!("live magenta attribution: {res}^2, t={} n_sync={} stride 2", ens.t_max, ens.n_sync);
    println!("  config: {}", ens.provenance());

    for name in &charts {
        let Some((n, chart, cx, cy, half, body)) = target(name) else {
            println!("{name}: unknown target");
            continue;
        };
        let cfg = SchedCfg {
            budget: 400,
            camera: Some(Camera::framing(cx, cy, half, res)),
            keep_pixels: true,
            chart,
            ..Default::default()
        };
        let (tree, st) = scheduler::descend_live(cx, cy, half, body, &cfg, &ens, Precision::F64);
        let px: Vec<&PixelOut> = tree.leaves().flat_map(|i| st.pixels.get(i).into_iter().flatten()).collect();
        let n_b = px.first().map(|p| p.live_t.len()).unwrap_or(0);
        println!("\n{n}: {} leaves, {} footprints, {n_b} boundaries", tree.leaves().count(), px.len());
        println!("  {:>3} {:>6}  {:>8} {:>8} {:>8} {:>8} {:>8}   {:>8}",
                 "j", "t", "run-wide", "live-nf", "shape[j]", "scalar[j]", "state", "magenta");
        let mut mags = Vec::new();
        for j in 0..n_b {
            let (mut run, mut lnf, mut sh, mut sc, mut stt, mut mag) = (0, 0, 0, 0, 0, 0);
            for p in &px {
                let q = scheduler::project_at(p, j);
                if p.n_nonfinite > 0 {
                    run += 1;
                }
                let bad_nf = q.n_nonfinite > 0;
                let bad_state = matches!(State::from_bits(q.state), Some(State::SimFailed) | Some(State::DecodeFailed) | None);
                let bad_shape = !q.shape_vec.iter().all(|v| v.is_finite());
                let bad_scalar = !q.spread_shape.is_finite();
                lnf += bad_nf as usize;
                sh += bad_shape as usize;
                sc += bad_scalar as usize;
                stt += bad_state as usize;
                // Exactly `colour::rgb`'s veto set, in its order.
                mag += (bad_nf || bad_state || bad_shape || bad_scalar) as usize;
            }
            mags.push(mag);
            println!("  {j:>3} {:>6.2}  {run:>8} {lnf:>8} {sh:>8} {sc:>8} {stt:>8}   {mag:>8}", px[0].live_t[j]);
        }
        let run_wide = px.iter().filter(|p| p.n_nonfinite > 0).count();
        println!("  run-wide {run_wide} at every boundary; magenta {} -> {}, monotone {}",
                 mags[0], mags[n_b - 1], mags.windows(2).all(|w| w[0] <= w[1]));

        // **What the flagged footprints are.** The count says how many; this says what. Every
        // column is read against the healthy population, because a value is only diagnostic
        // beside the one it is not.
        let bad: Vec<&&PixelOut> = px.iter().filter(|p| p.n_nonfinite > 0).collect();
        let good: Vec<&&PixelOut> = px.iter().filter(|p| p.n_nonfinite == 0).collect();
        if bad.is_empty() {
            continue;
        }
        let med = |v: &mut Vec<f64>| -> f64 {
            v.retain(|x| x.is_finite());
            if v.is_empty() {
                return f64::NAN;
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v[v.len() / 2]
        };
        let col = |set: &[&&PixelOut], f: &dyn Fn(&PixelOut) -> f64| -> f64 {
            let mut v: Vec<f64> = set.iter().map(|p| f(p)).collect();
            med(&mut v)
        };
        let frac = |set: &[&&PixelOut], f: &dyn Fn(&PixelOut) -> bool| -> f64 {
            set.iter().filter(|p| f(p)).count() as f64 / set.len().max(1) as f64
        };
        println!("  flagged {} of {} footprints ({:.3}%). copies flagged: {:?}",
                 bad.len(), px.len(), 100.0 * bad.len() as f64 / px.len() as f64, {
                     let mut h = std::collections::BTreeMap::new();
                     for p in &bad { *h.entry(p.n_nonfinite).or_insert(0usize) += 1; }
                     h
                 });
        println!("  {:>22} {:>12} {:>12}", "", "flagged", "healthy");
        for (name, f) in [
            ("d_min_true (median)", &(|p: &PixelOut| p.d_min_true) as &dyn Fn(&PixelOut) -> f64),
            ("energy_drift_max", &|p: &PixelOut| p.energy_drift_max),
            ("error_ratio", &|p: &PixelOut| p.error_ratio),
            ("dt_max", &|p: &PixelOut| p.dt_max),
            ("ab_min", &|p: &PixelOut| p.ab_min),
            ("substeps", &|p: &PixelOut| p.total_substeps as f64),
            ("n_cap_hits", &|p: &PixelOut| p.n_cap_hits as f64),
        ] {
            println!("  {name:>22} {:>12.4e} {:>12.4e}", col(&bad, f), col(&good, f));
        }
        for (name, f) in [
            ("budget_exhausted", &(|p: &PixelOut| p.budget_exhausted) as &dyn Fn(&PixelOut) -> bool),
            ("ab_floored", &|p: &PixelOut| p.ab_floored),
            ("SimFailed/DecodeFailed", &|p: &PixelOut| matches!(State::from_bits(p.state), Some(State::SimFailed) | Some(State::DecodeFailed) | None)),
            ("nominal shape non-finite", &|p: &PixelOut| !p.shape_vec.iter().all(|v| v.is_finite())),
        ] {
            println!("  {name:>22} {:>11.1}% {:>11.1}%", 100.0 * frac(&bad, f), 100.0 * frac(&good, f));
        }
        // **Why the fix may move no tree.** A footprint is unresolved if its spread exceeds the
        // tolerance, its event arm disagrees, **or** a copy is unusable. The leak added the third
        // reason too early; it can only change a decision on a footprint that was *not* already
        // unresolved for one of the first two. This counts the leak window -- the boundaries at
        // which the run-wide verdict fired and the live count had not -- and how much of it is
        // decision-relevant.
        let tau = SchedCfg::default().tau_display;
        let (mut window, mut decisive) = (0usize, 0usize);
        for p in &px {
            if p.n_nonfinite == 0 {
                continue;
            }
            for j in 0..n_b {
                let live_nf = p.live_nonfinite.get(j).copied().unwrap_or(p.n_nonfinite);
                if live_nf > 0 {
                    continue; // no leak here: the live count already fired
                }
                window += 1;
                let otherwise = p.live_spread_shape[j] > tau || p.live_spread_event[j] > 0.0;
                if !otherwise {
                    decisive += 1;
                }
            }
        }
        println!("  leak window {window} footprint-boundaries; of those {decisive} were not already \
                  unresolved by spread or event ({:.2}%) -- only those could move a decision",
                 100.0 * decisive as f64 / window.max(1) as f64);

        // **The decision is per QUAD, not per footprint.** Under `Policy::Tolerance` a quad splits
        // if *any* footprint is unresolved, so a falsely-unresolved one only tips the decision
        // when every other footprint in its quad is resolved. This counts quad-boundaries whose
        // `any unresolved` verdict differs between the two rules -- the number that decides
        // whether a tree can move at all.
        let mut quad_bnd = 0usize;
        let mut quad_diff = 0usize;
        for i in tree.leaves() {
            let Some(group) = st.pixels.get(i) else { continue };
            if group.is_empty() {
                continue;
            }
            for j in 0..n_b {
                quad_bnd += 1;
                let any = |old: bool| {
                    group.iter().any(|p: &PixelOut| {
                        let nf = if old { p.n_nonfinite } else { p.live_nonfinite.get(j).copied().unwrap_or(p.n_nonfinite) };
                        nf > 0 || p.live_spread_shape[j] > tau || p.live_spread_event[j] > 0.0
                    })
                };
                if any(true) != any(false) {
                    quad_diff += 1;
                }
            }
        }
        println!("  quad-boundaries {quad_bnd}; `any unresolved` differs between the rules on {quad_diff} \
                  ({:.2}%) -- a tree can only move where this is nonzero",
                 100.0 * quad_diff as f64 / quad_bnd.max(1) as f64);

        // Terminal class of the flagged set: a triple collision is the instrument reporting.
        let mut cls = std::collections::BTreeMap::new();
        for p in &bad {
            *cls.entry(format!("{:?}/{}", State::from_bits(p.state), p.detail)).or_insert(0usize) += 1;
        }
        println!("  flagged terminal (state/detail): {cls:?}");
    }
    println!("\nRoute 1 is `n_nonfinite`, a whole-run count `project_at` clones unchanged.");
    println!("`live-true` is what a viewer at that boundary could actually know.");
}
