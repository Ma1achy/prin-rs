//! **Experiment 1 and 2: the fixed Halton (2,3) prefix, and the `alpha_E` control variate.**
//!
//! The spec calls for copy offsets from a fixed low-discrepancy Halton (2,3) prefix indexed by
//! copy index. The port inherited the reference's per-pixel PCG stream instead — pseudo-random
//! and different in every footprint. Two properties were lost: even coverage at small `E`, and
//! the common-random-numbers structure that makes sampling noise cancel in the parent/child
//! ratio the refinement exponent is built from.
//!
//! Everything here is measured both ways, same grid, same integrator, same estimator.

use rayon::prelude::*;

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::grid::{self, Slice};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::physics::{burrau, energy, shape};

const FINE: usize = 64;

struct Copies {
    e0: Vec<f64>,
    et: Vec<f64>,
    shapes: Vec<[f64; 3]>,
}

/// Root-mean-square deviation from the median. Its expectation does not grow with sample size,
/// unlike an order statistic — which matters because a parent pools 4x as many copies as a child.
fn rms_dev(v: &[f64]) -> f64 {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    if x.len() < 2 {
        return f64::NAN;
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if x.len() % 2 == 1 {
        x[x.len() / 2]
    } else {
        0.5 * (x[x.len() / 2 - 1] + x[x.len() / 2])
    };
    (x.iter().map(|q| (q - med).powi(2)).sum::<f64>() / x.len() as f64).sqrt()
}

fn qs(v: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    if x.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| x[(((x.len() - 1) as f64) * f).round() as usize];
    (q(0.0), q(0.1), q(0.5), q(0.9), q(1.0))
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let v: Vec<(f64, f64)> = a
        .iter()
        .zip(b)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .map(|(x, y)| (*x, *y))
        .collect();
    if v.len() < 3 {
        return f64::NAN;
    }
    let n = v.len() as f64;
    let (ma, mb) = (
        v.iter().map(|p| p.0).sum::<f64>() / n,
        v.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for (x, y) in &v {
        num += (x - ma) * (y - mb);
        da += (x - ma).powi(2);
        db += (y - mb).powi(2);
    }
    if da <= 0.0 || db <= 0.0 {
        return f64::NAN;
    }
    num / (da.sqrt() * db.sqrt())
}

fn quads() -> Vec<[usize; 4]> {
    let mut v = Vec::new();
    for py in 0..FINE / 2 {
        for px in 0..FINE / 2 {
            v.push([
                2 * py * FINE + 2 * px,
                2 * py * FINE + 2 * px + 1,
                (2 * py + 1) * FINE + 2 * px,
                (2 * py + 1) * FINE + 2 * px + 1,
            ]);
        }
    }
    v
}

/// `sigma_E(0)` only — no integration, so this is cheap enough to sweep `E`.
fn e0_grid(s: &Slice, n_copies: usize, scheme: Scheme) -> Vec<Vec<f64>> {
    (0..s.npix())
        .into_par_iter()
        .map(|i| {
            jitter::copies_with::<f64>(s, i, n_copies - 1, 0.5, 0, scheme)
                .iter()
                .map(|x| energy::energy(&x.s.r, &x.s.v, &x.m, 0.0))
                .collect()
        })
        .collect()
}

fn child_median(kids: &[usize; 4], f: &dyn Fn(usize) -> f64) -> f64 {
    let mut v: Vec<f64> = kids.iter().map(|&k| f(k)).filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let s = grid::region("near-field", FINE, FINE, 0.05).unwrap();
    let qd = quads();

    // ---------------------------------------------------------------- experiment 1
    println!("=== EXPERIMENT 1: the noise floor, Halton against PCG ===");
    println!("near-field {FINE}x{FINE} -> {}x{} parents. alpha for sigma_E(0), whose true value",
             FINE / 2, FINE / 2);
    println!("is exactly 1.0, so every departure below is measurement error.");
    println!();
    println!("{:>8}{:>10}{:>12}{:>12}{:>14}{:>16}",
             "E+1", "scheme", "median", "p90-p10", "|median-1|", "parent/child r");

    for n_copies in [4usize, 8, 16, 32, 64, 128] {
        for scheme in [Scheme::Halton, Scheme::Pcg] {
            let e0 = e0_grid(&s, n_copies, scheme);
            let mut alpha = Vec::with_capacity(qd.len());
            let mut par = Vec::with_capacity(qd.len());
            let mut chi = Vec::with_capacity(qd.len());
            for kids in &qd {
                let pool: Vec<f64> = kids.iter().flat_map(|&k| e0[k].iter().cloned()).collect();
                let p = rms_dev(&pool);
                let c = child_median(kids, &|k| rms_dev(&e0[k]));
                alpha.push((p / c).log2());
                par.push(p);
                chi.push(c);
            }
            let (_, p10, med, p90, _) = qs(&alpha);
            println!("{n_copies:>8}{:>10}{med:>12.4}{:>12.4}{:>14.4}{:>16.4}",
                     format!("{scheme:?}"), p90 - p10, (med - 1.0).abs(), pearson(&par, &chi));
        }
    }

    println!();
    println!("The last column is the correlation between the parent spread estimate and the");
    println!("child one across quads. Common random numbers should make it near 1; independent");
    println!("draws per footprint should not.");
    println!();
    println!("Halton's |median-1| against E — if this is a 1/E law, the bias is geometric:");
    println!("{:>8}{:>14}{:>14}", "E+1", "|median-1|", "x (E+1)");
    for n_copies in [4usize, 8, 16, 32, 64, 128] {
        let e0 = e0_grid(&s, n_copies, Scheme::Halton);
        let alpha: Vec<f64> = qd
            .iter()
            .map(|kids| {
                let pool: Vec<f64> = kids.iter().flat_map(|&k| e0[k].iter().cloned()).collect();
                (rms_dev(&pool) / child_median(kids, &|k| rms_dev(&e0[k]))).log2()
            })
            .collect();
        let (_, _, med, _, _) = qs(&alpha);
        let d = (med - 1.0).abs();
        println!("{n_copies:>8}{d:>14.4}{:>14.3}", d * n_copies as f64);
    }
    println!("A constant last column means the excess falls as 1/E, which is what a geometric");
    println!("surrogate error diluting as the offset set fills the footprint would do.");

    // ---------------------------------------------------------------- experiment 2
    println!();
    println!("=== EXPERIMENT 2: the alpha_E control variate ===");
    println!("Integrating {}x{} at t=13 for both schemes; this is the expensive part.",
             FINE, FINE);
    println!();

    let m = burrau::masses::<f64>();
    for scheme in [Scheme::Halton, Scheme::Pcg] {
        let all: Vec<Copies> = (0..s.npix())
            .into_par_iter()
            .map(|i| {
                let c = jitter::copies_with::<f64>(&s, i, 7, 0.5, 0, scheme);
                let mut e0 = Vec::new();
                let mut et = Vec::new();
                let mut shapes = Vec::new();
                for x in &c {
                    e0.push(energy::energy(&x.s.r, &x.s.v, &x.m, 0.0));
                    let o = az::integrate_az_opts(
                        x.s, &x.m, 13.0, 32, 0.01, 30_000,
                        &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
                    );
                    et.push(energy::energy(&o.state.r, &o.state.v, &m, 0.0));
                    shapes.push(shape::shape_vec(&o.state.r, &m));
                }
                Copies { e0, et, shapes }
            })
            .collect();

        let mut a_e = Vec::with_capacity(qd.len());
        let mut a_sh = Vec::with_capacity(qd.len());
        let mut a_et = Vec::with_capacity(qd.len());
        for kids in &qd {
            let pe0: Vec<f64> = kids.iter().flat_map(|&k| all[k].e0.iter().cloned()).collect();
            let pet: Vec<f64> = kids.iter().flat_map(|&k| all[k].et.iter().cloned()).collect();
            let psh: Vec<[f64; 3]> =
                kids.iter().flat_map(|&k| all[k].shapes.iter().cloned()).collect();
            let lg = |p: f64, c: f64| if p > 0.0 && c > 0.0 { (p / c).log2() } else { f64::NAN };
            a_e.push(lg(rms_dev(&pe0), child_median(kids, &|k| rms_dev(&all[k].e0))));
            a_et.push(lg(rms_dev(&pet), child_median(kids, &|k| rms_dev(&all[k].et))));
            a_sh.push(lg(
                shape::spread_shape(&psh),
                child_median(kids, &|k| shape::spread_shape(&all[k].shapes)),
            ));
        }

        // The control variate. alpha_E's true value is exactly 1, so (alpha_E - 1) is OBSERVED
        // noise sharing a source with alpha_shape: same copies, same offsets, same integration.
        let dev: Vec<f64> = a_e.iter().map(|x| x - 1.0).collect();
        let rho = pearson(&dev, &a_sh);
        let pairs: Vec<(f64, f64)> = dev
            .iter()
            .zip(&a_sh)
            .filter(|(a, b)| a.is_finite() && b.is_finite())
            .map(|(a, b)| (*a, *b))
            .collect();
        let n = pairs.len() as f64;
        let md = pairs.iter().map(|p| p.0).sum::<f64>() / n;
        let ms = pairs.iter().map(|p| p.1).sum::<f64>() / n;
        let cov = pairs.iter().map(|p| (p.0 - md) * (p.1 - ms)).sum::<f64>() / n;
        let var_d = pairs.iter().map(|p| (p.0 - md).powi(2)).sum::<f64>() / n;
        let beta = if var_d > 0.0 { cov / var_d } else { f64::NAN };
        let corrected: Vec<f64> = dev.iter().zip(&a_sh).map(|(d, x)| x - beta * d).collect();
        // The additive form, beta fixed at 1. Under Halton the regression form is degenerate —
        // var(alpha_E) is ~0, so beta is a ratio to nothing — but the geometric bias the two
        // exponents SHARE is removed exactly by subtracting (alpha_E - 1).
        let additive: Vec<f64> = dev.iter().zip(&a_sh).map(|(d, x)| x - d).collect();

        let (_, e10, emed, e90, _) = qs(&a_e);
        let (_, s10, smed, s90, _) = qs(&a_sh);
        let (_, c10, cmed, c90, _) = qs(&corrected);
        let (_, ad10, admed, ad90, _) = qs(&additive);
        let (_, t10, tmed, t90, _) = qs(&a_et);

        println!("--- {scheme:?} ---");
        println!("{:>26}{:>12}{:>14}", "quantity", "median", "p90-p10");
        println!("{:>26}{emed:>12.4}{:>14.4}", "alpha sigma_E(0) [truth 1]", e90 - e10);
        println!("{:>26}{tmed:>12.4}{:>14.4}", "alpha sigma_E(t)", t90 - t10);
        println!("{:>26}{smed:>12.4}{:>14.4}", "alpha spread_shape", s90 - s10);
        println!("{:>26}{cmed:>12.4}{:>14.4}", "shape, regression beta", c90 - c10);
        println!("{:>26}{admed:>12.4}{:>14.4}", "shape, additive beta=1", ad90 - ad10);
        println!("  rho(alpha_E - 1, alpha_shape) = {rho:+.4}");
        println!("  fitted beta = {beta:+.4}");
        println!("  predicted variance reduction 1 - rho^2 = {:.4}", 1.0 - rho * rho);
        // Variance and interdecile imply very different widths for alpha_shape. They can only
        // do that if the distribution is heavy-tailed, and which one a scheduler should read
        // depends on the answer: it decides per TYPICAL quad, so the interdecile is the right
        // measure — and the interdecile is the one that did not move between schemes.
        {
            let f: Vec<f64> = a_sh.iter().cloned().filter(|x| x.is_finite()).collect();
            let nn = f.len() as f64;
            let mean = f.iter().sum::<f64>() / nn;
            let var = f.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / nn;
            let sd = var.sqrt();
            let kurt = f.iter().map(|x| ((x - mean) / sd).powi(4)).sum::<f64>() / nn - 3.0;
            println!("  alpha_shape sd = {sd:.4}, interdecile = {:.4}, ratio = {:.3} (normal 2.563)",
                     s90 - s10, (s90 - s10) / sd);
            println!("  excess kurtosis = {kurt:.2} (normal 0) — the variance lives in the tails");
        }
        println!("  measured floor ratio, regression      = {:.4}", (c90 - c10) / (s90 - s10));
        println!("  measured floor ratio, additive        = {:.4}", (ad90 - ad10) / (s90 - s10));
        println!("  var(alpha_E) = {var_d:.3e}   var(alpha_shape) = {:.3e}",
                 pairs.iter().map(|p| (p.1 - ms).powi(2)).sum::<f64>() / n);
        println!();
    }
    println!("rho near zero means the control variate buys nothing and should be dropped.");
    println!("This is a per-quad correlation between two EXPONENTS across quads — not the");
    println!("per-pixel correlation between two within-footprint estimators measured earlier.");
}
