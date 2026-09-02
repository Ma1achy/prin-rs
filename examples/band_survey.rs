//! **Is the phase-beat structure local to one window, or general?**
//!
//! `config_stability`'s banded region was measured to be a **regular island**: two ICs a whole
//! band apart are the same trajectory shifted in orbital phase, the lag grows **linearly**
//! (0.025 -> 0.190 over `t = 5..45`, ratios exactly 2/1, 3/2, 4/3, ...) and the correlation holds
//! at **0.9999** to the end. That is a frequency beat between two slightly different pair
//! periods, not chaotic divergence -- and one full cycle of accumulated phase came to 79.8 px
//! against the raw field's 80.95 px tier.
//!
//! This surveys whether that holds elsewhere. The discriminator needs no spectrum:
//!
//! | | lag growth | corr at `t ~ 45` |
//! |---|---|---|
//! | **regular island** — phase beat | linear in `t` | ~1.0 |
//! | **chaotic** | exponential | collapses |
//!
//! **The correlation is the honest arm and the lag is the fragile one.** Once two trajectories
//! have decorrelated the "lag" is the argmax of noise -- a number, always, and meaningless. So
//! `corr` is printed beside every lag and a lag whose correlation has collapsed is marked, never
//! quietly tabulated. Same shape as this project's standing rule that a difference can be small
//! because both sides are dead.
//!
//! Probes are a GRID over each frame, not points picked by eye: the question is which colours
//! and which regions behave which way, and choosing the samples by their appearance would decide
//! the answer before measuring it. The rendered colour is reported per probe so the mapping from
//! appearance to behaviour is read off rather than assumed.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::outcome::State;
use prin_rs::output::colour::{self, Scalar};

const NSY: usize = 20_000;
const T_MAX: f64 = 50.0;

/// Best lag of `b` against `a` over a window, with the correlation achieved there.
fn lag_at(a: &[f64], b: &[f64], lo: usize, hi: usize, rng: isize) -> (f64, f64) {
    let n = a.len();
    let seg = &a[lo..hi];
    let mu = seg.iter().sum::<f64>() / seg.len() as f64;
    let sa: f64 = (seg.iter().map(|x| (x - mu) * (x - mu)).sum::<f64>() / seg.len() as f64).sqrt();
    let (mut best_l, mut best_c) = (0isize, -2.0f64);
    for l in -rng..=rng {
        let (s, e) = (lo as isize + l, hi as isize + l);
        if s < 0 || e > n as isize {
            continue;
        }
        let w = &b[s as usize..e as usize];
        let mw = w.iter().sum::<f64>() / w.len() as f64;
        let sw: f64 = (w.iter().map(|x| (x - mw) * (x - mw)).sum::<f64>() / w.len() as f64).sqrt();
        if sa <= 0.0 || sw <= 0.0 {
            continue;
        }
        let c = seg.iter().zip(w).map(|(x, y)| (x - mu) * (y - mw)).sum::<f64>()
            / (seg.len() as f64 * sa * sw);
        if c > best_c {
            best_c = c;
            best_l = l;
        }
    }
    (best_l as f64 * (T_MAX / n as f64), best_c)
}

struct Probe {
    case: String,
    fx: f64,
    fy: f64,
    rgb: [u8; 3],
    state: &'static str,
    lag_e: f64,
    corr_e: f64,
    lag_l: f64,
    corr_l: f64,
    period: f64,
}

fn main() {
    let out: String = std::env::args().nth(1).unwrap_or_else(|| "results/osc".into());
    let nprobe: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let _ = std::fs::create_dir_all(&out);

    // (name, chart, cx, cy, half). The Burrau regions use their own windows via `grid`.
    let mut cases: Vec<(String, Chart, f64, f64, f64)> = Vec::new();
    {
        let (c, x, y, h) = Chart::config_stability();
        cases.push(("config_stability".into(), c, x, y, h));
        let (c, x, y, h) = Chart::config_basin();
        cases.push(("config_basin".into(), c, x, y, h));
    }
    for (nm, ch) in [
        ("preset_shape", Chart::preset_shape()),
        ("preset_prho", Chart::preset_prho()),
        ("preset_plambda", Chart::preset_plambda()),
        ("preset_shape_pl", Chart::preset_shape_pl()),
    ] {
        let h = ch.default_half();
        cases.push((nm.into(), ch, 0.0, 0.0, h));
    }

    let cfg = EnsembleCfg::production().with_overrides(&[
        Override::TMax(T_MAX), Override::NSync(125),
        Override::RefineFlagged(false), Override::MaxSteps(4_000_000),
    ]);
    let hopts = HgOpts::<f64> {
        r_coll_frac: 0.0, stop_on_event: false, keep_boundary_shapes: true, ..Default::default()
    };

    println!("BAND SURVEY -- is the phase beat local to one window?");
    println!("`corr` is the arm to read. A lag whose correlation has collapsed is the argmax of");
    println!("noise; those rows are marked `--` rather than tabulated as a number.");
    println!();

    let mut rows: Vec<Probe> = Vec::new();
    for (nm, chart, cx, cy, half) in &cases {
        // One coarse render of the whole frame, for the colour window and each probe's colour.
        let res = 96usize;
        let sl = grid::Slice::body_plane(res, res, *cx, *cy, *half, 0).with_chart(*chart);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
        let m = grid::decode_state(chart, 0, *cx, *cy).m;
        let sites = colour::landmarks(&m);

        // Probes on an interior grid: the frame edges are where a window is least representative.
        let probes: Vec<(f64, f64)> = (0..nprobe)
            .flat_map(|i| (0..nprobe).map(move |j| {
                ((i as f64 + 0.5) / nprobe as f64, (j as f64 + 0.5) / nprobe as f64)
            }))
            .collect();

        let got: Vec<Probe> = probes.par_iter().map(|&(fx, fy)| {
            let (x, y) = (cx - half + 2.0 * half * fx, cy - half + 2.0 * half * fy);
            let ix = ((fx * res as f64) as usize).min(res - 1);
            let iy = ((fy * res as f64) as usize).min(res - 1);
            let p = &px[iy * res + ix];
            let rgb = colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);
            let state = State::from_bits(p.state).map(|s| s.name()).unwrap_or("?");

            // The pair: one probe and a neighbour ONE COARSE PIXEL away. The growth SHAPE is
            // what is being read, and it does not depend on the separation.
            let d = 2.0 * half / res as f64;
            let a = grid::decode_state(chart, 0, x, y);
            let b = grid::decode_state(chart, 0, x + d, y);
            let oa = integrate_hg(a.s, &a.m, T_MAX, NSY, 1e-2, 4_000_000, &hopts);
            let ob = integrate_hg(b.s, &b.m, T_MAX, NSY, 1e-2, 4_000_000, &hopts);
            let (sa, sb): (Vec<f64>, Vec<f64>) = (
                oa.boundary_shapes.iter().map(|v| v[0]).collect(),
                ob.boundary_shapes.iter().map(|v| v[0]).collect(),
            );
            let n = sa.len().min(sb.len());
            if n < 4000 || !oa.finite || !ob.finite {
                return Probe { case: nm.clone(), fx, fy, rgb, state,
                               lag_e: f64::NAN, corr_e: f64::NAN, lag_l: f64::NAN,
                               corr_l: f64::NAN, period: f64::NAN };
            }
            let (le, ce) = lag_at(&sa, &sb, (0.10 * n as f64) as usize, (0.20 * n as f64) as usize, 900);
            let (ll, cl) = lag_at(&sa, &sb, (0.80 * n as f64) as usize, (0.90 * n as f64) as usize, 900);

            // The pair period, from the dominant late oscillation. Zero crossings of the
            // mean-removed late half: cheap, and immune to the FFT's bin spacing.
            let late = &sa[n / 2..];
            let mu = late.iter().sum::<f64>() / late.len() as f64;
            let mut cross = 0usize;
            for w in late.windows(2) {
                if (w[0] - mu) * (w[1] - mu) < 0.0 { cross += 1; }
            }
            let period = if cross > 1 {
                2.0 * (T_MAX / 2.0) / cross as f64
            } else { f64::NAN };

            Probe { case: nm.clone(), fx, fy, rgb, state,
                    lag_e: le, corr_e: ce, lag_l: ll, corr_l: cl, period }
        }).collect();
        rows.extend(got);
    }

    println!("{:>18} {:>5} {:>5} {:>14} {:>10} {:>9} {:>8} {:>9} {:>8} {:>9} {:>8}",
             "case", "fx", "fy", "rgb", "state", "T_pair", "lag_e", "corr_e", "lag_l", "corr_l", "growth");
    let mut reg = 0usize;
    let mut cha = 0usize;
    for r in &rows {
        let live = r.corr_l > 0.90;
        let growth = if !r.corr_l.is_finite() { "n/a".to_string() }
            else if !live { "DECORRELATED".to_string() }
            else if r.lag_e.abs() > 0.0 { format!("x{:.2}", r.lag_l.abs() / r.lag_e.abs()) }
            else { "flat".to_string() };
        if live && r.corr_l > 0.99 { reg += 1; } else if r.corr_l.is_finite() { cha += 1; }
        println!("{:>18} {:>5.2} {:>5.2} {:>14} {:>10} {:>9.4} {:>8} {:>9.4} {:>8} {:>9.4} {:>8}",
                 r.case, r.fx, r.fy,
                 format!("{:3},{:3},{:3}", r.rgb[0], r.rgb[1], r.rgb[2]),
                 r.state, r.period,
                 if live { format!("{:.4}", r.lag_e.abs()) } else { "--".into() }, r.corr_e,
                 if live { format!("{:.4}", r.lag_l.abs()) } else { "--".into() }, r.corr_l,
                 growth);
    }
    println!();
    println!("  regular (corr_l > 0.99): {reg}    decorrelated or chaotic: {cha}    of {}", rows.len());
    println!();
    println!("`growth` is lag(t~42)/lag(t~7). LINEAR growth reads about x6 (42/7); exponential");
    println!("reads far higher and its correlation is gone before it gets there, so the two are");
    println!("told apart by `corr_l`, not by the ratio.");
}
