//! Per-box diagnostics for the marked regions: what the panel shows there, and what the initial
//! conditions are. The IC half needs no integration, so ten boxes cost milliseconds.
use rayon::prelude::*;
use prin_rs::grid::{self, Chart};
use prin_rs::physics::{energy, newton, THIRD};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
/// `(name, u, v, half)` as fractions of the panel, origin top-left, digitised from the mark-up.
const BOXES: [(&str, f64, f64, f64); 16] = [
    ("B1", 0.890, 0.245, 0.0517), ("B2", 0.862, 0.332, 0.0571),
    ("B3", 0.379, 0.446, 0.0326), ("B4", 0.476, 0.497, 0.0294),
    ("B5", 0.590, 0.437, 0.0294), ("B6", 0.608, 0.528, 0.0337),
    ("B7", 0.419, 0.664, 0.0381), ("B8", 0.383, 0.748, 0.0403),
    ("B9", 0.807, 0.838, 0.0566), ("B10", 0.942, 0.789, 0.0522),
    // the second set, in purple
    ("P1", 0.539, 0.199, 0.0354), ("P2", 0.447, 0.426, 0.0354),
    ("P3", 0.533, 0.457, 0.0305), ("P4", 0.335, 0.682, 0.0408),
    ("P5", 0.428, 0.718, 0.0381), ("P6", 0.510, 0.742, 0.0397),
];

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0; q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx0, cy0, half0) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);

    println!("{n}x{n} decode samples per box. INITIAL CONDITIONS ONLY -- no integration.\n");
    println!("{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}  {:>22}  {:>22}",
        "box", "d(0,1)", "d(0,2)", "d(1,2)", "aspect", "alpha", "beta", "|Lz|",
        "ref body 0/1/2", "tightest (01)/(02)/(12)");
    // The frame, as the baseline row.
    let mut rows: Vec<(String, f64, f64, f64, f64, f64, f64, f64, [f64; 3], [f64; 3])> = Vec::new();
    for (name, u, v, h) in BOXES.iter().map(|b| (b.0, b.1, b.2, b.3)).chain([("FRAME", 0.5, 0.5, 0.5)]) {
        let cx = cx0 + (2.0 * u - 1.0) * half0;
        let cy = cy0 + (2.0 * v - 1.0) * half0;
        let half = 2.0 * h * half0;
        let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
        let d: Vec<[f64; 8]> = (0..sl.npix()).into_par_iter().map(|k| {
            let (x, y) = sl.decode_pos(k);
            let s = grid::decode_state(&chart, 0, x, y);
            let dd = newton::pair_dists(&s.s.r);
            let (mut lo, mut hi) = (0usize, 0usize);
            for j in 1..3 { if dd[j] > dd[hi] { hi = j; } if dd[j] < dd[lo] { lo = j; } }
            let (al, be) = prin_rs::physics::decoder::angles(Z[1] + y, Z[0] + x);
            let lz: f64 = (0..3).map(|i| s.m[i] * (s.s.r[i].x * s.s.v[i].y - s.s.r[i].y * s.s.v[i].x)).sum();
            let _ = energy::energy(&s.s.r, &s.s.v, &s.m, 0.0);
            [dd[0], dd[1], dd[2], al, be, lz.abs(), THIRD[hi] as f64, lo as f64]
        }).collect();
        let col = |i: usize| { let mut v: Vec<f64> = d.iter().map(|r| r[i]).collect(); q(&mut v, 0.5) };
        let fr = |i: usize, t: f64| d.iter().filter(|r| r[i] == t).count() as f64 / d.len() as f64;
        let (a, b, c) = (col(0), col(1), col(2));
        let asp = { let mut v: Vec<f64> = d.iter().map(|r| {
            let m = r[0].min(r[1]).min(r[2]); let x = r[0].max(r[1]).max(r[2]); x / m }).collect(); q(&mut v, 0.5) };
        rows.push((name.into(), a, b, c, asp, col(3), col(4), col(5),
            [fr(6, 0.0), fr(6, 1.0), fr(6, 2.0)], [fr(7, 0.0), fr(7, 1.0), fr(7, 2.0)]));
    }
    for (nm, a, b, c, asp, al, be, lz, rb, tp) in rows {
        println!("{nm:>4} {a:>8.4} {b:>8.4} {c:>8.4} {asp:>8.3} {al:>8.4} {be:>8.4} {lz:>8.4}  \
                  {:>6.3}{:>8.3}{:>8.3}  {:>6.3}{:>8.3}{:>8.3}",
            rb[0], rb[1], rb[2], tp[0], tp[1], tp[2]);
    }
}
