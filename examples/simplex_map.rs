//! Scratch: the outcome map of chart_reference §4.5's ternary mass plot (Principia's form), Burrau geometry at rest,
//! over the F3-extended view (u, v) in [-0.3, 1.3]^2. Negative-mass pixels are 255 (forbidden, never integrated).
//! Run: cargo run --release --example simplex_map -- [n] [t_max] [out]
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::outcome::{self, EscapeRule, CLOSURE_TAU};
use prin_rs::physics::decoder;
use prin_rs::vec2::Vec2;
use rayon::prelude::*;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).and_then(|x| x.parse().ok()).unwrap_or(320);
    let t_max: f64 = a.get(2).and_then(|x| x.parse().ok()).unwrap_or(20.0);
    let out = a.get(3).cloned().unwrap_or_else(|| "simplex_map.bin".into());
    let n_sync = (t_max / 0.4).round() as usize;
    let opts = HgOpts::<f64> { r_coll_frac: 1e-3, stop_on_event: true, stop_on_escape: true,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU), closure_k: 1, ..Default::default() };
    let eps = 1e-4;
    let px: Vec<u8> = (0..n * n).into_par_iter().map(|idx| {
        let (i, j) = (idx % n, idx / n);
        let u = -0.3 + 1.6 * (i as f64 + 0.5) / n as f64;
        let v = 1.3 - 1.6 * (j as f64 + 0.5) / n as f64; // row 0 = top = v high
        let raw = [(1.0 - u) * (1.0 - v), u, (1.0 - u) * v];
        if raw.iter().any(|&m| !(m >= 0.0)) { return 255u8; }
        let m = [0, 1, 2].map(|k| (1.0 - 3.0 * eps) * raw[k] + eps);
        let mut r = [Vec2::new(1.0, -1.0), Vec2::new(-2.0, -1.0), Vec2::new(1.0, 3.0)];
        let c = (r[0] * m[0] + r[1] * m[1] + r[2] * m[2]) / (m[0] + m[1] + m[2]);
        for k in 0..3 { r[k] = r[k] - c; }
        let mut p = [Vec2::new(0.0, 0.0); 3];
        let _ = decoder::scale_gauge(&mut r, &mut p, &m);
        let ic = decoder::to_ic(&r, &p, &m);
        let o = integrate_hg(ic.s, &ic.m, t_max, n_sync, 1e-2, 20_000_000, &opts);
        outcome::classify(&o.events, &o.state, &ic.m, o.finite, o.budget_exhausted).pack()
    }).collect();
    let mut hist = std::collections::BTreeMap::new();
    for p in &px { *hist.entry(*p).or_insert(0usize) += 1; }
    println!("{hist:?}");
    std::fs::write(&out, &px).unwrap();
}
