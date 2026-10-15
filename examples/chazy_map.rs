//! The Chazy class map as of time t, on one gallery chart, at several playheads.
//!
//! One march per pixel, classified at every snapshot by `chazy::classify` with the shipped
//! `ChazyCfg::default()`; the panels are the class at each playhead, left to right. Same
//! trajectories in every panel, so a pixel changing colour is the label evolving, not noise
//! between runs.
//!
//! Run: cargo run --release --example chazy_map -- [case] [n] [out_dir]

use prin_rs::chazy::{self, ChazyCfg, Latch, Motion};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::output::adaptive;
use prin_rs::physics::{energy, shape};
use prin_rs::grid;
use rayon::prelude::*;

const SNAP: f64 = 0.4;
const PLAYHEADS: [f64; 4] = [10.0, 50.0, 200.0, 800.0];

fn colour(c: Motion) -> [u8; 3] {
    match c {
        Motion::H => [0x34, 0x62, 0xE0],
        Motion::HP => [0x7F, 0xA0, 0xF0],
        Motion::P => [0xC0, 0xD0, 0xFA],
        Motion::HE => [0xDE, 0x2D, 0x2D],
        Motion::PE => [0xF2, 0x96, 0x20],
        Motion::BStable => [0x2E, 0xBC, 0x4E],
        Motion::BSoFar => [0x14, 0x14, 0x18],
        Motion::Transient => [0xE0, 0x34, 0xC6],
        Motion::OsCandidate => [0xF0, 0xDE, 0x32],
        Motion::Excursion => [0x30, 0xC8, 0xDC],
        Motion::Unresolved => [0x4E, 0x4E, 0x4E],
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or_else(|| "latent_shape".into());
    let n: usize = a.get(2).and_then(|x| x.parse().ok()).unwrap_or(64);
    let out = a.get(3).cloned().unwrap_or_else(|| "results/chazy".into());
    let (_, chart, cx, cy, half) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == case)
        .unwrap_or_else(|| panic!("no gallery case {case}"));
    let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
    let cfg = ChazyCfg::default();
    let opts = HgOpts::<f64> { stop_on_event: false, r_coll_frac: 0.0, ..Default::default() };
    let t0 = std::time::Instant::now();
    let t_end = *PLAYHEADS.last().unwrap();
    let px: Vec<(Vec<Motion>, bool)> = (0..sl.npix())
        .into_par_iter()
        .map(|idx| {
            let ic = sl.nominal_ic::<f64>(idx);
            let m = ic.m;
            let mut s = ic.s;
            let r0 = energy::hyperradius(&s.r, &m);
            let e = energy::energy(&s.r, &s.v, &m, 0.0);
            let mut latch = Latch::new(energy::inertia(&s.r, &m));
            let mut at = Vec::new();
            let mut ph = 0;
            let steps = (t_end / SNAP).round() as usize;
            for k in 1..=steps {
                let o = integrate_hg(s, &m, SNAP, 2, 1e-2, 2_000_000, &opts);
                if !o.finite || o.budget_exhausted {
                    while at.len() < PLAYHEADS.len() { at.push(Motion::Unresolved); }
                    return (at, false);
                }
                s = o.state;
                let t = k as f64 * SNAP;
                let i_now = energy::inertia(&s.r, &m);
                let nv = shape::shape_vec(&s.r, &m);
                let c = chazy::classify(&s, &m, e, r0, i_now, &nv, &latch, &cfg);
                latch.i_max = latch.i_max.max(i_now);
                latch.n_prev = Some(nv);
                while ph < PLAYHEADS.len() && t >= PLAYHEADS[ph] - 1e-9 {
                    at.push(c);
                    ph += 1;
                }
            }
            (at, true)
        })
        .collect();
    let failed = px.iter().filter(|p| !p.1).count();
    println!("# chazy_map case={case} n={n} cfg={cfg:?} snap={SNAP} heggie eta=1e-2 r_coll=0");
    println!("# {:.1}s, {failed} of {} trajectories failed (non-finite or budget), drawn grey", t0.elapsed().as_secs_f64(), px.len());
    print!("{:>12}", "class");
    for t in PLAYHEADS { print!(" {:>8}", format!("t={t}")); }
    println!();
    for c in Motion::ALL {
        print!("{:>12}", c.name());
        for p in 0..PLAYHEADS.len() {
            let k = px.iter().filter(|x| x.0[p] == c).count();
            print!(" {:>7.2}%", 100.0 * k as f64 / px.len() as f64);
        }
        println!();
    }
    // A strip of panels with a 4 px gutter.
    let g = 4;
    let w = PLAYHEADS.len() * n + (PLAYHEADS.len() - 1) * g;
    let mut buf = vec![255u8; w * n * 3];
    for p in 0..PLAYHEADS.len() {
        for idx in 0..sl.npix() {
            let (ix, iy) = (idx % n, idx / n);
            let o = (iy * w + p * (n + g) + ix) * 3;
            buf[o..o + 3].copy_from_slice(&colour(px[idx].0[p]));
        }
    }
    let path = format!("{out}/chazy_map_{case}_n{n}.png");
    adaptive::save_rect(&path, w, n, &buf).unwrap();
    println!("# wrote {path}");
}
