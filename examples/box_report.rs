//! **Full per-box report** for the marked regions: outcome, drift, cost, timing, switching, and
//! the initial conditions, side by side.
//!
//! The boxes are the ones drawn on `results/closure/config_stability_stop0_uniform.png`, in
//! fractions of the panel with the origin top-left. Every box is integrated at its OWN window at
//! `closure_render`'s settings -- `t_max = 50`, `r_coll = 0.005`, `n_sync = round(t_max/0.4)` --
//! so this is the same experiment as the panel, sampled finer inside each box.
//!
//! `detail` is decoded rather than printed as a number, because it means different things on the
//! two arms: under `Collision` it is the pair index, or `3` for two-or-more (the >=2-pair triple
//! rule); under `Escape` it is the **escaping body**, which the reference derives as the lowest
//! firing index and not the physically escaping one -- transcribed, and stated wherever it is read.
//!
//! `FRAME` is the whole panel window at the same sample count, as the baseline every row is
//! against. A per-box number without it is a number, not a comparison.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, State, CLOSURE_TAU};
use prin_rs::physics::{energy, newton, THIRD};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
const WINDOW: f64 = 0.4;

const BOXES: [(&str, f64, f64, f64); 17] = [
    ("B1", 0.890, 0.245, 0.0517), ("B2", 0.862, 0.332, 0.0571),
    ("B3", 0.379, 0.446, 0.0326), ("B4", 0.476, 0.497, 0.0294),
    ("B5", 0.590, 0.437, 0.0294), ("B6", 0.608, 0.528, 0.0337),
    ("B7", 0.419, 0.664, 0.0381), ("B8", 0.383, 0.748, 0.0403),
    ("B9", 0.807, 0.838, 0.0566), ("B10", 0.942, 0.789, 0.0522),
    ("P1", 0.539, 0.199, 0.0354), ("P2", 0.447, 0.426, 0.0354),
    ("P3", 0.533, 0.457, 0.0305), ("P4", 0.335, 0.682, 0.0408),
    ("P5", 0.428, 0.718, 0.0381), ("P6", 0.510, 0.742, 0.0397),
    ("FRAME", 0.5, 0.5, 0.5),
];

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn class_name(state: u8, detail: u8) -> String {
    const PAIR: [&str; 3] = ["(0,1)", "(0,2)", "(1,2)"];
    match State::from_bits(state) {
        Some(State::Collision) => if detail == 3 { "collision triple".into() }
            else { format!("collision {}", PAIR[(detail as usize).min(2)]) },
        Some(State::Escape) => if detail == 3 { "escape triple".into() }
            else { format!("escape body {detail}") },
        Some(State::Bounded) => "bounded".into(),
        Some(State::Running) => "running".into(),
        Some(State::SimFailed) => "SIM FAILED".into(),
        Some(State::DecodeFailed) => "DECODE FAILED".into(),
        None => format!("invalid {state}"),
    }
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(128);
    let only: String = std::env::args().nth(2).unwrap_or_else(|| "all".into());
    // The repair mechanism this project built for exactly this failure -- BRIEF's cliff, where
    // `eta = 1e-2` stops being sufficient -- re-integrates flagged pixels at finer `eta`.
    // `EnsembleCfg::default()` has it ON; `closure_render`, which made the committed panel, turns
    // it OFF. Argument 3 switches it so the two can be compared on the same window.
    let refine: bool = std::env::args().nth(3).map(|s| s == "on").unwrap_or(false);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0; q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx0, cy0, half0) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (t_max, r_coll) = (50.0f64, 0.005f64);
    let n_sync = (t_max / WINDOW).round().max(4.0) as usize;

    println!(
        "{n}x{n} per box, t_max {t_max}, r_coll {r_coll}, n_sync {n_sync}, eta 1e-2, \
         escape rule Closure({CLOSURE_TAU:e}),\n\
         refine_flagged {}.\n\
         `detail` decoded: on the collision arm it is the PAIR; on the escape arm it is the\n\
         reference's escaping-body index -- the LOWEST FIRING INDEX, not necessarily the body\n\
         that physically leaves. Transcribed, and said here because this table reads it.\n",
        if refine { "ON" } else { "OFF (as `closure_render` sets it)" }
    );

    for (name, u, v, h) in BOXES {
        if only != "all" && only != name { continue; }
        let cx = cx0 + (2.0 * u - 1.0) * half0;
        let cy = cy0 + (2.0 * v - 1.0) * half0;
        let half = 2.0 * h * half0;
        let ens = EnsembleCfg {
            refine_flagged: refine, t_max, n_sync, r_coll_frac: r_coll,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU), closure_k: 1,
            stop_on_escape: false, keep_drift_hist: true, ..Default::default()
        };
        let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
        let t0 = std::time::Instant::now();
        let px: Vec<PixelOut> = (0..sl.npix()).into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens)).collect();
        let secs = t0.elapsed().as_secs_f64();
        let np = px.len() as f64;
        let ncopy = (ens.n_extra + 1) as f64;

        // ---- initial conditions, no integration ------------------------------------------
        let ic: Vec<[f64; 9]> = (0..sl.npix()).into_par_iter().map(|k| {
            let (x, y) = sl.decode_pos(k);
            let s = grid::decode_state(&chart, 0, x, y);
            let d = newton::pair_dists(&s.s.r);
            let (mut lo, mut hi) = (0usize, 0usize);
            for j in 1..3 { if d[j] > d[hi] { hi = j; } if d[j] < d[lo] { lo = j; } }
            let (al, be) = prin_rs::physics::decoder::angles(Z[1] + y, Z[0] + x);
            let lz: f64 = (0..3).map(|i| s.m[i] * (s.s.r[i].x * s.s.v[i].y - s.s.r[i].y * s.s.v[i].x)).sum();
            let e0 = energy::energy(&s.s.r, &s.s.v, &s.m, 0.0);
            let kin: f64 = (0..3).map(|i| 0.5 * s.m[i] * s.s.v[i].norm_sq()).sum();
            [d[0], d[1], d[2], al, be, lz, e0, kin / (e0 - kin).abs(), (THIRD[hi] * 3 + lo) as f64]
        }).collect();

        let f = |g: &dyn Fn(&PixelOut) -> f64, p: f64| {
            let mut v: Vec<f64> = px.iter().map(g).filter(|x| x.is_finite()).collect();
            q(&mut v, p)
        };
        let fi = |i: usize, p: f64| { let mut v: Vec<f64> = ic.iter().map(|r| r[i]).collect(); q(&mut v, p) };

        println!("================ {name}  ({n}x{n} = {} px, {secs:.1}s) ================", px.len());
        println!(
            "  OUTCOME   escape {:.4}  collision {:.4}  bounded {:.4}  running {:.4}  \
             nonfinite {:.4}  simfail {:.4}",
            px.iter().filter(|p| p.state == 0).count() as f64 / np,
            px.iter().filter(|p| p.state == 2).count() as f64 / np,
            px.iter().filter(|p| p.state == 1).count() as f64 / np,
            px.iter().filter(|p| p.state == 3).count() as f64 / np,
            px.iter().filter(|p| p.n_nonfinite > 0).count() as f64 / np,
            px.iter().filter(|p| p.state == 4).count() as f64 / np,
        );
        let mut hist: std::collections::BTreeMap<(u8, u8), usize> = Default::default();
        for p in &px { *hist.entry((p.state, p.detail)).or_default() += 1; }
        let mut hv: Vec<_> = hist.into_iter().collect();
        hv.sort_by_key(|&(_, c)| std::cmp::Reverse(c));
        print!("  CLASS    ");
        for ((s, d), c) in hv.iter().take(6) {
            print!(" {} {:.4} |", class_name(*s, *d), *c as f64 / np);
        }
        println!();
        println!(
            "  TIGHTEST PAIR at the last boundary (binary_id)  (0,1) {:.4}  (0,2) {:.4}  (1,2) {:.4}",
            px.iter().filter(|p| p.binary_id == 0).count() as f64 / np,
            px.iter().filter(|p| p.binary_id == 1).count() as f64 / np,
            px.iter().filter(|p| p.binary_id == 2).count() as f64 / np,
        );
        println!(
            "  DRIFT     p50 {:.3e}  p90 {:.3e}  p99 {:.3e}  max {:.3e}   frac>1e-6 {:.4}   \
             gamma p50 {:.3e}",
            f(&|p| p.energy_drift_max, 0.5), f(&|p| p.energy_drift_max, 0.9),
            f(&|p| p.energy_drift_max, 0.99), f(&|p| p.energy_drift_max, 1.0),
            px.iter().filter(|p| !(p.energy_drift_max <= 1e-6)).count() as f64 / np,
            f(&|p| p.gamma_max, 0.5),
        );
        println!(
            "  COST      steps/copy p50 {:.0}  p90 {:.0}  max {:.0}   (nominal ~{:.0})",
            f(&|p| p.total_substeps as f64 / ncopy, 0.5),
            f(&|p| p.total_substeps as f64 / ncopy, 0.9),
            f(&|p| p.total_substeps as f64 / ncopy, 1.0),
            n_sync as f64 / ens.eta,
        );
        println!(
            "  TIMING    t_end p10 {:.3}  p50 {:.3}  p90 {:.3}   frozen (t_end<t_max) {:.4}   \
             on a sync boundary {:.4}",
            f(&|p| p.t_end, 0.1), f(&|p| p.t_end, 0.5), f(&|p| p.t_end, 0.9),
            px.iter().filter(|p| p.t_end < t_max - 1e-9).count() as f64 / np,
            px.iter().filter(|p| {
                let dt = t_max / n_sync as f64;
                (p.t_end / dt - (p.t_end / dt).round()).abs() < 1e-9
            }).count() as f64 / np,
        );
        println!(
            "  APPROACH  d_min_true p10 {:.3e}  p50 {:.3e}   d_min_ref p50 {:.3e}   \
             frac < r_coll {:.4}",
            f(&|p| p.d_min_true, 0.1), f(&|p| p.d_min_true, 0.5), f(&|p| p.d_min_ref, 0.5),
            px.iter().filter(|p| p.d_min_true < r_coll * 1.0).count() as f64 / np,
        );
        println!(
            "  SWITCHES  count p50 {:.0}  p90 {:.0}  max {:.0}   never {:.4}   \
             t_first p50 {:.3}   jump switch/hold p50 {:.3e} / {:.3e}",
            f(&|p| p.switches as f64, 0.5), f(&|p| p.switches as f64, 0.9),
            f(&|p| p.switches as f64, 1.0),
            px.iter().filter(|p| p.switches == 0).count() as f64 / np,
            f(&|p| p.t_first_switch, 0.5),
            f(&|p| p.switch_jump_med, 0.5), f(&|p| p.hold_jump_med, 0.5),
        );
        println!(
            "  SPREAD    shape p50 {:.3e}  p90 {:.3e}   event p50 {:.3e}   \
             error_ratio p50 {:.3}  p99 {:.3}",
            f(&|p| p.spread_shape, 0.5), f(&|p| p.spread_shape, 0.9),
            f(&|p| p.spread_event, 0.5), f(&|p| p.error_ratio, 0.5), f(&|p| p.error_ratio, 0.99),
        );
        let rb = |k: usize| ic.iter().filter(|r| (r[8] as usize) / 3 == k).count() as f64 / np;
        let tp = |k: usize| ic.iter().filter(|r| (r[8] as usize) % 3 == k).count() as f64 / np;
        println!(
            "  ICs       d(0,1) {:.4}  d(0,2) {:.4}  d(1,2) {:.4}   alpha {:.4}  beta {:.4}   \
             Lz {:+.4}  E0 {:.4}  K/|U| {:.4}",
            fi(0, 0.5), fi(1, 0.5), fi(2, 0.5), fi(3, 0.5), fi(4, 0.5),
            fi(5, 0.5), fi(6, 0.5), fi(7, 0.5),
        );
        println!(
            "            AZ reference at t=0  body0 {:.3} body1 {:.3} body2 {:.3}   \
             tightest at t=0  (0,1) {:.3} (0,2) {:.3} (1,2) {:.3}",
            rb(0), rb(1), rb(2), tp(0), tp(1), tp(2),
        );
        println!();
    }
}
