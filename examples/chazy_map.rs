//! The Chazy class map as of time t, on one gallery chart, at several playheads.
//!
//! One march per pixel, classified at every snapshot by `chazy::classify` with the shipped
//! `ChazyCfg::default()`; the panels are the class at each playhead, left to right. Same
//! trajectories in every panel, so a pixel changing colour is the label evolving, not noise
//! between runs.
//!
//! Besides the PNG strip it writes the per-pixel data the stain previews are recoloured from, so a
//! palette change needs no re-integration (all in slice order, playhead-major, beside the PNG):
//! `.cls` the class (`Motion as u8`); `.st` the structure (0 dispersal, 1 hierarchy, 2 interacting,
//! 3 oscillatory) and `.q` its `K/|U|` as f32; `.esc` the escaper index and `.qk` each body's own
//! `K/|U|` against the other two (3 x f32); `.m` the masses once per pixel (3 x f32).
//!
//! Run: cargo run --release --example chazy_map -- [case] [n] [out_dir]

use prin_rs::chazy::{self, ChazyCfg, Latch, Motion};
use prin_rs::grid;
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::output::adaptive;
use prin_rs::physics::{energy, shape};
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
    let opts = HgOpts::<f64> {
        stop_on_event: false,
        r_coll_frac: 0.0,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let t_end = *PLAYHEADS.last().unwrap();
    let px: Vec<(
        Vec<Motion>,
        bool,
        Vec<(u8, f32)>,
        Vec<(u8, [f32; 3])>,
        [f32; 3],
    )> = (0..sl.npix())
        .into_par_iter()
        .map(|idx| {
            let ic = sl.nominal_ic::<f64>(idx);
            let m = ic.m;
            let mut s = ic.s;
            let r0 = energy::hyperradius(&s.r, &m);
            let e = energy::energy(&s.r, &s.v, &m, 0.0);
            let mut latch = Latch::new(energy::inertia(&s.r, &m));
            let mut at = Vec::new();
            let mut sqs = Vec::new();
            let mut bods: Vec<(u8, [f32; 3])> = Vec::new();
            let mut ph = 0;
            let steps = (t_end / SNAP).round() as usize;
            for k in 1..=steps {
                let o = integrate_hg(s, &m, SNAP, 2, 1e-2, 2_000_000, &opts);
                if !o.finite || o.budget_exhausted {
                    while at.len() < PLAYHEADS.len() {
                        at.push(Motion::Unresolved);
                        sqs.push((255u8, f32::NAN));
                        bods.push((255u8, [f32::NAN; 3]));
                    }
                    return (at, false, sqs, bods, m.map(|x| x as f32));
                }
                s = o.state;
                let t = k as f64 * SNAP;
                let i_now = energy::inertia(&s.r, &m);
                let nv = shape::shape_vec(&s.r, &m);
                let c = chazy::classify(&s, &m, e, r0, i_now, &nv, &latch, &cfg);
                let sq = structure_q(&s, &m, c, r0, &cfg);
                let bod = bodies(&s, &m);
                latch.i_max = latch.i_max.max(i_now);
                latch.n_prev = Some(nv);
                while ph < PLAYHEADS.len() && t >= PLAYHEADS[ph] - 1e-9 {
                    at.push(c);
                    sqs.push(sq);
                    bods.push(bod);
                    ph += 1;
                }
            }
            (at, true, sqs, bods, m.map(|x| x as f32))
        })
        .collect();
    let failed = px.iter().filter(|p| !p.1).count();
    println!("# chazy_map case={case} n={n} cfg={cfg:?} snap={SNAP} heggie eta=1e-2 r_coll=0");
    println!(
        "# {:.1}s, {failed} of {} trajectories failed (non-finite or budget), drawn grey",
        t0.elapsed().as_secs_f64(),
        px.len()
    );
    print!("{:>12}", "class");
    for t in PLAYHEADS {
        print!(" {:>8}", format!("t={t}"));
    }
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
    // Class indices (`Motion as u8`), playhead-major then pixel in slice order, so a palette can
    // be re-applied without re-integrating.
    let mut cls = Vec::with_capacity(PLAYHEADS.len() * sl.npix());
    for p in 0..PLAYHEADS.len() {
        cls.extend(px.iter().map(|x| x.0[p] as u8));
    }
    let cpath = format!("{out}/chazy_map_{case}_n{n}.cls");
    std::fs::write(&cpath, &cls).unwrap();
    println!("# wrote {cpath}");
    // Structure (0 dispersal, 1 hierarchy, 2 interacting, 3 oscillatory) and the energy ratio
    // q = K/|U| of the separating motion, same layout as the class dump; q as little-endian f32.
    let mut st = Vec::new();
    let mut qb = Vec::new();
    for p in 0..PLAYHEADS.len() {
        for x in &px {
            st.push(x.2[p].0);
            qb.extend_from_slice(&x.2[p].1.to_le_bytes());
        }
    }
    std::fs::write(format!("{out}/chazy_map_{case}_n{n}.st"), &st).unwrap();
    std::fs::write(format!("{out}/chazy_map_{case}_n{n}.q"), &qb).unwrap();
    println!("# wrote structure and q dumps");
    // Per body: the escaper index (the body most separated from its pair) and each body's K/U
    // against the other two's centre of mass; plus the masses once per pixel.
    let mut esc = Vec::new();
    let mut qk = Vec::new();
    for p in 0..PLAYHEADS.len() {
        for x in &px {
            esc.push(x.3[p].0);
            for k in 0..3 {
                qk.extend_from_slice(&x.3[p].1[k].to_le_bytes());
            }
        }
    }
    let mut ms = Vec::new();
    for x in &px {
        for k in 0..3 {
            ms.extend_from_slice(&x.4[k].to_le_bytes());
        }
    }
    std::fs::write(format!("{out}/chazy_map_{case}_n{n}.esc"), &esc).unwrap();
    std::fs::write(format!("{out}/chazy_map_{case}_n{n}.qk"), &qk).unwrap();
    std::fs::write(format!("{out}/chazy_map_{case}_n{n}.m"), &ms).unwrap();
    println!("# wrote per-body dumps");
}

/// The bivariate view's two inputs at one snapshot: which bodies are separating from which
/// (the hue), and `q = K/|U|` of that separation (the lightness; 1 is parabolic).
fn structure_q(
    s: &prin_rs::physics::Cart<f64>,
    m: &[f64; 3],
    c: Motion,
    _r0: f64,
    cfg: &ChazyCfg,
) -> (u8, f32) {
    let pairs = [(0usize, 1usize), (0, 2), (1, 2)].map(|(i, j)| chazy::pair_kup(s, m, i, j));
    let total = || {
        let k: f64 = (0..3).map(|i| 0.5 * m[i] * s.v[i].norm_sq()).sum();
        let u: f64 = pairs.iter().map(|x| x.1).sum();
        (k / u) as f32
    };
    match c {
        Motion::H | Motion::HP | Motion::P => (
            0,
            pairs
                .iter()
                .map(|x| x.0 / x.1)
                .fold(f64::INFINITY, f64::min) as f32,
        ),
        Motion::OsCandidate => (3, total()),
        _ => {
            // The body most separated from its pair, as the classifier's hierarchy test.
            let (b, ratio) = (0..3)
                .map(|b| {
                    let (i, j) = chazy::other(b);
                    let k = [(0, 1), (0, 2), (1, 2)]
                        .iter()
                        .position(|&p| p == (i, j))
                        .unwrap();
                    (b, chazy::out_kup(s, m, b).2 / pairs[k].2)
                })
                .fold((0, 0.0), |a, x| if x.1 > a.1 { x } else { a });
            if ratio > cfg.h_far {
                let (ko, uo, _, _) = chazy::out_kup(s, m, b);
                (1, (ko / uo) as f32)
            } else {
                (2, total())
            }
        }
    }
}

/// The escaper index and every body's own `K/U` against the other two's centre of mass.
fn bodies(s: &prin_rs::physics::Cart<f64>, m: &[f64; 3]) -> (u8, [f32; 3]) {
    let pairs = [(0usize, 1usize), (0, 2), (1, 2)].map(|(i, j)| chazy::pair_kup(s, m, i, j));
    let (b, _) = (0..3)
        .map(|b| {
            let (i, j) = chazy::other(b);
            let k = [(0, 1), (0, 2), (1, 2)]
                .iter()
                .position(|&p| p == (i, j))
                .unwrap();
            (b, chazy::out_kup(s, m, b).2 / pairs[k].2)
        })
        .fold((0, 0.0), |a, x| if x.1 > a.1 { x } else { a });
    let q = [0, 1, 2].map(|k| {
        let (ko, uo, _, _) = chazy::out_kup(s, m, k);
        (ko / uo) as f32
    });
    (b as u8, q)
}
