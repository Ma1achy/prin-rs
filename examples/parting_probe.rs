//! When have two runs of the same initial condition "parted"? Measured, to set the honesty
//! viewer's threshold.
//!
//! For each IC, four Heggie runs to the horizon (termination off, so every run reaches it):
//!   f64, f64 nudged by one ulp per component, f32 (from the f32-rounded IC), f32 nudged.
//! Each run records the shape-sphere vector n at a fine fixed boundary grid. Three pairs are
//! compared by the chord |n_a - n_b| (0..2; random uniform pairs average 4/3):
//!   f32 vs f64      what the atlas's precision costs
//!   f64 vs f64'     f64's own horizon
//!   f32 vs f32'     f32's own horizon
//!
//! Reported:
//!   1. the chord's growth: quantiles of log10(chord) at a ladder of times, and where it saturates
//!   2. parting time t_part(theta) = first time the chord exceeds theta, over a threshold ladder,
//!      and how much it moves with theta (exponential growth predicts ~ln(ratio)/lambda)
//!   3. whether the threshold MEANS anything: f32 and f64 first-event labels compared, split by
//!      whether the first event happened before or after t_part(theta). A useful theta gives
//!      near-perfect agreement before and near-chance after.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts, HgOut};
use prin_rs::physics::Ic;
use prin_rs::Real;
use rayon::prelude::*;

const T_MAX: f64 = 50.0;
const DT_SYNC: f64 = 0.1;

fn next_up64(x: f64) -> f64 {
    if x == 0.0 { f64::from_bits(1) } else if x > 0.0 { f64::from_bits(x.to_bits() + 1) } else { f64::from_bits(x.to_bits() - 1) }
}
fn next_up32(x: f32) -> f32 {
    if x == 0.0 { f32::from_bits(1) } else if x > 0.0 { f32::from_bits(x.to_bits() + 1) } else { f32::from_bits(x.to_bits() - 1) }
}

fn opts<T: Real>(cfg: &EnsembleCfg) -> HgOpts<T> {
    let n = (T_MAX / DT_SYNC).round() as usize;
    let _ = n;
    HgOpts::<T> {
        keep_boundary_shapes: true,
        stop_on_event: false,
        r_coll_frac: T::lit(cfg.r_coll_frac),
        escape_rule: cfg.escape_rule.lift(),
        closure_k: ((cfg.closure_k as f64) * (cfg.t_max / cfg.n_sync as f64) / DT_SYNC).round() as usize,
        step_limit_f: T::lit(cfg.step_limit_f),
        ..HgOpts::<T>::default()
    }
}

struct Run {
    n: Vec<[f64; 3]>,
    ev: Option<(u8, u8, f64)>, // (kind 0 collision / 1 escape, detail, time)
    ok: bool,
}

fn first_event<T: Real>(o: &HgOut<T>) -> Option<(u8, u8, f64)> {
    let c = o.events.collision.map(|(m, t)| (0u8, m, t.to_f64().unwrap()));
    let e = o.events.escape.map(|(b, t)| (1u8, b, t.to_f64().unwrap()));
    match (c, e) {
        (Some(a), Some(b)) => Some(if a.2 <= b.2 { a } else { b }),
        (a, b) => a.or(b),
    }
}

fn run<T: Real>(ic: &Ic<T>, cfg: &EnsembleCfg) -> Run {
    let n_sync = (T_MAX / DT_SYNC).round() as usize;
    // physical step held at production's eta * dt_sync
    let eta = cfg.eta * (cfg.t_max / cfg.n_sync as f64) / DT_SYNC;
    let o = integrate_hg(ic.s, &ic.m, T::lit(T_MAX), n_sync, T::lit(eta), 4_000_000, &opts::<T>(cfg));
    Run {
        n: o.boundary_shapes.iter().map(|v| [v[0].to_f64().unwrap(), v[1].to_f64().unwrap(), v[2].to_f64().unwrap()]).collect(),
        ev: first_event(&o),
        ok: !o.budget_exhausted,
    }
}

fn chord(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
    if d.is_finite() { d } else { 2.0 }
}
fn series(a: &Run, b: &Run) -> Vec<f64> {
    let m = a.n.len().min(b.n.len());
    (0..m).map(|k| chord(&a.n[k], &b.n[k])).collect()
}
fn t_part(s: &[f64], theta: f64) -> Option<f64> {
    s.iter().position(|&c| c > theta).map(|k| (k + 1) as f64 * DT_SYNC)
}
fn q(mut v: Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let per: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(120);
    let cfg = EnsembleCfg::production();
    let mut ics: Vec<(String, Ic<f64>)> = Vec::new();
    let mut sets: Vec<(String, Slice)> = Vec::new();
    let (c, cx, cy, h) = Chart::config_stability();
    sets.push(("config_stability".into(), Slice::body_plane(64, 64, cx, cy, h, 0).with_chart(c)));
    for &(n, cx, cy, body) in grid::REGIONS.iter() {
        if n == "near-field" || n == "deep interior" {
            sets.push((n.into(), Slice::body_plane(64, 64, cx, cy, 0.05, body)));
        }
    }
    for (n, c, cx, cy, h) in grid::gallery_cases() {
        if n == "preset_shape" {
            sets.push((n.into(), Slice::body_plane(64, 64, cx, cy, h, 0).with_chart(c)));
        }
    }
    // deterministic spread of pixels over each slice
    for (name, sl) in &sets {
        let np = sl.npix();
        for i in 0..per {
            let k = (i * 2_654_435_761usize) % np;
            ics.push((name.clone(), sl.nominal_ic::<f64>(k)));
        }
    }

    let results: Vec<(String, Run, Run, Run, Run)> = ics
        .par_iter()
        .map(|(name, ic)| {
            let mut icn = *ic;
            for b in 0..3 {
                icn.s.r[b].x = next_up64(icn.s.r[b].x); icn.s.r[b].y = next_up64(icn.s.r[b].y);
                icn.s.v[b].x = next_up64(icn.s.v[b].x); icn.s.v[b].y = next_up64(icn.s.v[b].y);
            }
            let ic32: Ic<f32> = ic.cast();
            let mut ic32n = ic32;
            for b in 0..3 {
                ic32n.s.r[b].x = next_up32(ic32n.s.r[b].x); ic32n.s.r[b].y = next_up32(ic32n.s.r[b].y);
                ic32n.s.v[b].x = next_up32(ic32n.s.v[b].x); ic32n.s.v[b].y = next_up32(ic32n.s.v[b].y);
            }
            (name.clone(), run(ic, &cfg), run(&icn, &cfg), run(&ic32, &cfg), run(&ic32n, &cfg))
        })
        .collect();

    let thetas = [1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0];
    for group in ["config_stability", "near-field", "deep interior", "preset_shape", "ALL"] {
        let rs: Vec<&(String, Run, Run, Run, Run)> =
            results.iter().filter(|r| group == "ALL" || r.0 == group).filter(|r| r.1.ok && r.2.ok && r.3.ok && r.4.ok).collect();
        let dropped = results.iter().filter(|r| group == "ALL" || r.0 == group).count() - rs.len();
        println!("\n===== {group}: {} ICs ({} with a run out of step budget, listed not used)", rs.len(), dropped);
        let pairs: [(&str, fn(&(String, Run, Run, Run, Run)) -> Vec<f64>); 3] = [
            ("f32 vs f64", |r| series(&r.3, &r.1)),
            ("f64 vs f64'", |r| series(&r.1, &r.2)),
            ("f32 vs f32'", |r| series(&r.3, &r.4)),
        ];
        // 1. growth
        println!("  chord quantiles (p10 / p50 / p90) at t =");
        let ts = [0.1, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0];
        print!("  {:>12}", "pair");
        for t in ts { print!(" {:>24}", format!("t={t}")); }
        println!();
        for (pn, f) in &pairs {
            let ss: Vec<Vec<f64>> = rs.iter().map(|r| f(r)).collect();
            print!("  {:>12}", pn);
            for t in ts {
                let k = ((t / DT_SYNC).round() as usize).saturating_sub(1);
                let v: Vec<f64> = ss.iter().filter_map(|s| s.get(k).copied()).collect();
                print!(" {:>24}", format!("{:.1e}/{:.1e}/{:.1e}", q(v.clone(), 0.1), q(v.clone(), 0.5), q(v, 0.9)));
            }
            println!();
        }
        // 2. parting time vs threshold
        println!("  parting time t_part(theta): median [p10, p90], fraction parted by t = 50");
        for (pn, f) in &pairs {
            let ss: Vec<Vec<f64>> = rs.iter().map(|r| f(r)).collect();
            print!("  {:>12}", pn);
            for &th in &thetas {
                let tp: Vec<f64> = ss.iter().filter_map(|s| t_part(s, th)).collect();
                let frac = tp.len() as f64 / ss.len() as f64;
                print!("  θ={:<6} {:>5.1} [{:>4.1},{:>4.1}] {:>4.2}", th, q(tp.clone(), 0.5), q(tp.clone(), 0.1), q(tp, 0.9), frac);
            }
            println!();
        }
        // 3. does theta mean anything? f32 vs f64 first-event agreement before / after t_part
        println!("  f32 vs f64 first-event label agreement, by whether the event came before or after t_part(θ)");
        let ss: Vec<Vec<f64>> = rs.iter().map(|r| series(&r.3, &r.1)).collect();
        for &th in &thetas {
            let (mut nb, mut ab, mut na, mut aa) = (0, 0, 0, 0);
            for (r, s) in rs.iter().zip(&ss) {
                let tp = t_part(s, th).unwrap_or(f64::INFINITY);
                let (e32, e64) = (r.3.ev, r.1.ev);
                // the earlier of the two first events decides "before" or "after"
                let te = match (e32, e64) { (Some(a), Some(b)) => a.2.min(b.2), (Some(a), None) | (None, Some(a)) => a.2, (None, None) => continue };
                let agree = match (e32, e64) { (Some(a), Some(b)) => a.0 == b.0 && a.1 == b.1, _ => false };
                if te < tp { nb += 1; if agree { ab += 1; } } else { na += 1; if agree { aa += 1; } }
            }
            println!("    θ={:<6} before: {:>4} events, agree {:>5.3}   after: {:>4} events, agree {:>5.3}",
                th, nb, ab as f64 / nb.max(1) as f64, na, aa as f64 / na.max(1) as f64);
        }
    }
}
