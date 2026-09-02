//! **The guard `band_survey` was missing: is each probe pair alive on both sides?**
//!
//! `config_basin` returned lag **exactly 0.0000** and correlation **exactly 1.0000** on all 16
//! probes, and several `preset_prho`/`preset_plambda` rows did the same. Perfect agreement is
//! what a regular island looks like AND what two identical inputs look like AND what two frozen
//! outputs look like. Three states, one number.
//!
//! Two things must hold before a correlation means anything:
//!
//! 1. **The ICs differ.** `decode::distinct` is the standing form of this. Note the keying trap
//!    already on record: positions in the latent decode do not depend on the momentum
//!    coordinates, so on `preset_prho`/`preset_plambda` every pixel is the SAME triangle at a
//!    different velocity and a position-keyed check reads "collapsed" when nothing has.
//!    Positions and velocities are reported separately for that reason.
//! 2. **The signal varies.** A trajectory whose pair has escaped has a frozen shape vector, so
//!    two of them correlate at 1.0000 while carrying no information. `sd_late` is that arm.
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};

fn main() {
    let nprobe: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(4);
    let mut cases: Vec<(String, Chart, f64, f64, f64)> = Vec::new();
    {
        let (c, x, y, h) = Chart::config_stability();
        cases.push(("config_stability".into(), c, x, y, h));
        let (c, x, y, h) = Chart::config_basin();
        cases.push(("config_basin".into(), c, x, y, h));
    }
    for (nm, ch) in [
        ("preset_prho", Chart::preset_prho()),
        ("preset_plambda", Chart::preset_plambda()),
    ] {
        let h = ch.default_half();
        cases.push((nm.into(), ch, 0.0, 0.0, h));
    }

    let hopts = HgOpts::<f64> {
        r_coll_frac: 0.0, stop_on_event: false, keep_boundary_shapes: true, ..Default::default()
    };
    println!("{:>18} {:>5} {:>5} {:>11} {:>11} {:>11} {:>11} {:>9}",
             "case", "fx", "fy", "|dr|", "|dv|", "sd_late_a", "sd_late_b", "verdict");
    for (nm, chart, cx, cy, half) in &cases {
        let res = 96usize;
        for i in 0..nprobe {
            for j in 0..nprobe {
                let (fx, fy) = ((i as f64 + 0.5) / nprobe as f64, (j as f64 + 0.5) / nprobe as f64);
                let (x, y) = (cx - half + 2.0 * half * fx, cy - half + 2.0 * half * fy);
                let d = 2.0 * half / res as f64;
                let a = grid::decode_state(chart, 0, x, y);
                let b = grid::decode_state(chart, 0, x + d, y);
                let dr: f64 = (0..3).map(|k| {
                    (a.s.r[k].x - b.s.r[k].x).powi(2) + (a.s.r[k].y - b.s.r[k].y).powi(2)
                }).sum::<f64>().sqrt();
                let dv: f64 = (0..3).map(|k| {
                    (a.s.v[k].x - b.s.v[k].x).powi(2) + (a.s.v[k].y - b.s.v[k].y).powi(2)
                }).sum::<f64>().sqrt();
                let oa = integrate_hg(a.s, &a.m, 50.0, 20_000, 1e-2, 4_000_000, &hopts);
                let ob = integrate_hg(b.s, &b.m, 50.0, 20_000, 1e-2, 4_000_000, &hopts);
                let sd = |o: &prin_rs::integrate::heggie::HgOut<f64>| {
                    let v: Vec<f64> = o.boundary_shapes.iter().map(|s| s[0]).collect();
                    if v.len() < 100 { return f64::NAN; }
                    let l = &v[v.len() / 2..];
                    let m = l.iter().sum::<f64>() / l.len() as f64;
                    (l.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / l.len() as f64).sqrt()
                };
                let (sa, sb) = (sd(&oa), sd(&ob));
                // A correlation is only readable when the ICs differ AND the signal moves.
                let verdict = if !(dr > 0.0 || dv > 0.0) { "IC IDENTICAL" }
                    else if !(sa > 1e-9 && sb > 1e-9) { "FROZEN" }
                    else { "live" };
                println!("{nm:>18} {fx:>5.2} {fy:>5.2} {dr:>11.3e} {dv:>11.3e} {sa:>11.3e} {sb:>11.3e} {verdict:>9}");
            }
        }
    }
}
