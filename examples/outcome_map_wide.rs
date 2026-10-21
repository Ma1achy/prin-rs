//! The event outcome per pixel, with its time, for rendering the outcome palette at several playheads.
//!
//! One march per pixel (the nominal IC only), Heggie at the production collision radius and escape
//! rule, stopping at the first collision or confirmed escape. Dumps `(state, detail)` packed as
//! `outcome::Outcome::pack` and `t_end` as f32, in slice order; a playhead shows a pixel's outcome
//! once `t_end <= playhead` and black before.
//!
//! Run: cargo run --release --example outcome_map -- [case] [n] [t_max] [out_dir]

use prin_rs::grid;
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::outcome::{self, EscapeRule, CLOSURE_TAU};
use rayon::prelude::*;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a
        .get(1)
        .cloned()
        .unwrap_or_else(|| "preset_shape_h1".into());
    let n: usize = a.get(2).and_then(|x| x.parse().ok()).unwrap_or(48);
    let t_max: f64 = a.get(3).and_then(|x| x.parse().ok()).unwrap_or(800.0);
    let out = a
        .get(4)
        .cloned()
        .unwrap_or_else(|| "results/outcome".into());
    let (_, chart, cx, cy, half0) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == case)
        .unwrap_or_else(|| panic!("no gallery case {case}"));
    let half = half0 * a.get(5).and_then(|x| x.parse::<f64>().ok()).unwrap_or(1.0); println!("# half={half}");
    let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
    // The closure window is a time: 0.4 per sync interval, as the reference.
    let n_sync = (t_max / 0.4).round() as usize;
    let opts = HgOpts::<f64> {
        r_coll_frac: 1e-3,
        stop_on_event: true,
        stop_on_escape: true,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let px: Vec<(u8, f32)> = (0..sl.npix())
        .into_par_iter()
        .map(|idx| {
            let ic = sl.nominal_ic::<f64>(idx);
            let o = integrate_hg(ic.s, &ic.m, t_max, n_sync, 1e-2, 20_000_000, &opts);
            let oc = outcome::classify(&o.events, &o.state, &ic.m, o.finite, o.budget_exhausted);
            (oc.pack(), o.t_end as f32)
        })
        .collect();
    let mut hist = std::collections::BTreeMap::new();
    for p in &px {
        *hist.entry(p.0).or_insert(0usize) += 1;
    }
    println!("# outcome_map case={case} n={n} t_max={t_max} n_sync={n_sync} heggie eta=1e-2 r_coll=1e-3 closure");
    println!("# {:.1}s", t0.elapsed().as_secs_f64());
    for (k, v) in &hist {
        let o = outcome::Outcome::unpack(*k).unwrap();
        println!("{:>10} detail {}  {:>6}", o.state.name(), o.detail, v);
    }
    let mut buf = Vec::with_capacity(px.len() * 5);
    for p in &px {
        buf.push(p.0);
    }
    for p in &px {
        buf.extend_from_slice(&p.1.to_le_bytes());
    }
    std::fs::create_dir_all(&out).unwrap();
    let path = format!("{out}/outcome_map_{case}_n{n}.bin");
    std::fs::write(&path, &buf).unwrap();
    println!("# wrote {path}");
}
