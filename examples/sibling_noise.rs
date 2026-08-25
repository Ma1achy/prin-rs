//! §3.4 — how noisy is `alpha_sibling_spread` itself?
//!
//! It is the range (max - min) of a parent's four children's `alpha`, flagged twice as needing
//! its own noise characterised before it is trusted. **The range of four samples is a noisy
//! statistic**, and if its scatter is comparable to the signal it cannot separate a reliable
//! exponent from an unreliable one. `sib_tau` ships at 0.5 and is compared against it.
//!
//! # Part 1 is a control that could not fail, and is kept as the demonstration
//!
//! The obvious control is `sigma_E(0)` — the spread of *initial* energies, proportional to the
//! jitter and therefore to the cell width, so its true `alpha` is exactly 1.0 and its true
//! sibling range exactly 0. It needs no integration.
//!
//! It reads **~0.003 and does not move with `N` or `E+1` at all**, which is the tell. With the
//! fixed Halton prefix there is **no sampling randomness anywhere in it**: the offsets are
//! fixed, the footprint positions are fixed, so `sigma_E(0)` is a deterministic function of the
//! box and its `alpha` is deterministic too. The residual is the geometric non-linearity of the
//! energy over the box and nothing else.
//!
//! So part 1 bounds the *geometric* term and is **structurally incapable of measuring sampling
//! noise**. It is retained, labelled, and read as a floor — not deleted, because the floor is
//! genuinely useful: anything above it in part 2 is not geometry.
//!
//! # Part 2 measures what part 1 cannot
//!
//! Sampling noise on the **real** signal needs an actual random ensemble, so part 2 runs under
//! `Scheme::Pcg` — the per-pixel stream — and varies the seed. The *same quad* is measured
//! several times; how far its `alpha` moves between seeds is sampling noise, with everything
//! else held fixed.
//!
//! `Scheme::Pcg` is used **here and nowhere else**, and never as a default: it is the only way
//! to get a second independent draw of an ensemble whose whole design is to have no draw.
//!
//! # How to misread this
//!
//! **`sd` is the wrong statistic.** Excess kurtosis on these distributions is 110; the variance
//! lives in the tail while a scheduler decides per typical quad. Every column here is an
//! interdecile.
//!
//! **A large sibling range is not automatically noise.** Compare `sib` against `seed move`: if
//! the range is large and the seed movement small, the siblings genuinely disagree and the
//! statistic is doing its job. If they are comparable, it is measuring its own draw.

use rayon::prelude::*;

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid::{self, Slice};
use prin_rs::physics::energy;
use prin_rs::quad::Agg;
use prin_rs::scheduler::reduce;
use prin_rs::stats;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// RMS deviation from the median: its expectation does not grow with sample size, unlike an
/// order statistic. The same estimator the `alpha_E` control variate uses.
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

fn sigma0(cx: f64, cy: f64, half: f64, body: usize, n: usize, e: usize, sc: Scheme) -> f64 {
    let s = Slice::body_plane(n, n, cx, cy, half, body);
    let e: Vec<f64> = (0..s.npix())
        .into_par_iter()
        .flat_map(|i| {
            jitter::copies_with::<f64>(&s, i, e, 0.5, 0, sc)
                .into_iter()
                .map(|x| energy::energy(&x.s.r, &x.s.v, &x.m, 0.0))
                .collect::<Vec<f64>>()
        })
        .collect();
    rms_dev(&e)
}

/// The real criterion signal for one quad: `ensemble_spread` aggregated by median.
fn signal(cx: f64, cy: f64, half: f64, body: usize, n: usize, ens: &EnsembleCfg, tau: f64) -> f64 {
    let s = Slice::body_plane(n, n, cx, cy, half, body);
    let px: Vec<_> = (0..s.npix()).into_par_iter().map(|i| evaluate::<f64>(&s, i, ens)).collect();
    reduce(&px, n, tau, ens.t_max).spread(Agg::Median)
}

/// `(alphas, sibling ranges)` over every parent above `levels`, under an arbitrary per-quad
/// scalar.
fn sweep(
    cx0: f64,
    cy0: f64,
    half0: f64,
    levels: u32,
    f: &(dyn Fn(f64, f64, f64) -> f64 + Sync),
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut alphas, mut ranges, mut parents) = (Vec::new(), Vec::new(), Vec::new());
    for l in 0..levels {
        let w = 1u32 << l;
        let h = half0 / (1u64 << l) as f64;
        for iy in 0..w {
            for ix in 0..w {
                let pcx = cx0 - half0 + (2 * ix + 1) as f64 * h;
                let pcy = cy0 - half0 + (2 * iy + 1) as f64 * h;
                let sp = f(pcx, pcy, h);
                let q = h / 2.0;
                let kid: Vec<f64> = [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0)]
                    .iter()
                    .map(|(dx, dy)| {
                        let sc = f(pcx + dx * q, pcy + dy * q, q);
                        if sp > 0.0 && sc > 0.0 { (sp / sc).log2() } else { f64::NAN }
                    })
                    .collect();
                if kid.iter().all(|x| x.is_finite()) {
                    alphas.extend_from_slice(&kid);
                    let hi = kid.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let lo = kid.iter().cloned().fold(f64::INFINITY, f64::min);
                    ranges.push(hi - lo);
                    parents.push(kid.iter().sum::<f64>() / 4.0);
                }
            }
        }
    }
    (alphas, ranges, parents)
}

fn main() {
    let levels: u32 = arg(1, 5);
    let sig_levels: u32 = arg(2, 3);

    // ---------------- part 1: the deterministic floor ----------------
    println!(
        "PART 1 -- control: sigma_E(0), true alpha EXACTLY 1.0, true sibling range EXACTLY 0.\n\
         near-field, levels 0..{levels}, no integration.\n\
         This measures the GEOMETRIC residual only. Under the fixed Halton prefix there is no\n\
         sampling randomness in it at all, so it CANNOT measure sampling noise -- it is a floor,\n\
         and anything above it in part 2 is not geometry.\n"
    );
    println!(
        "{:>7} {:>5} {:>7} {:>9} {:>9} {:>9} {:>9}",
        "N", "E+1", "parents", "a p50", "a idec", "sib p50", "sib p90"
    );
    let s = grid::region("near-field", 2, 2, 0.05).unwrap();
    for n in [4usize, 8, 16] {
        for e1 in [2usize, 8, 32] {
            let f = move |cx: f64, cy: f64, h: f64| sigma0(cx, cy, h, s.body, n, e1 - 1, Scheme::Halton);
            let (a, r, _) = sweep(s.cx, s.cy, 0.05, levels, &f);
            let (_, ap50, _, aid) = stats::interdecile(&a);
            let (_, rp50, rp90, _) = stats::interdecile(&r);
            println!("{n:>7} {e1:>5} {:>7} {ap50:>9.4} {aid:>9.4} {rp50:>9.4} {rp90:>9.4}", r.len());
        }
    }

    // ---------------- part 2: the real signal, with a real draw ----------------
    println!(
        "\nPART 2 -- the REAL signal (ensemble_spread, median), under Scheme::Pcg so that\n\
         varying the seed is an independent draw. Levels 0..{sig_levels}, t=13.\n\
         Pcg is used HERE ONLY and is never a default: the shipped Halton prefix is fixed by\n\
         design, so there is no second draw of it to compare against.\n"
    );
    println!(
        "{:>14} {:>5} {:>7} {:>9} {:>9} {:>9} {:>9} {:>10}",
        "region", "seed", "parents", "a p50", "a idec", "sib p50", "sib p90", "seed move"
    );

    for region in ["near-field", "deep interior"] {
        let s = grid::region(region, 2, 2, 0.05).unwrap();
        let mut per_seed: Vec<Vec<f64>> = Vec::new();
        for seed in 0..3u64 {
            let ens = EnsembleCfg {
                jitter_scheme: Scheme::Pcg,
                seed,
                refine_flagged: false,
                ..Default::default()
            };
            let f = move |cx: f64, cy: f64, h: f64| signal(cx, cy, h, s.body, 8, &ens, 1e-4);
            let (a, r, _) = sweep(s.cx, s.cy, 0.05, sig_levels, &f);
            let (_, ap50, _, aid) = stats::interdecile(&a);
            let (_, rp50, rp90, _) = stats::interdecile(&r);
            per_seed.push(a.clone());
            // Seed-to-seed movement of the SAME quad's alpha: sampling noise with everything
            // else held fixed. Only defined from the second seed on.
            let mv = if per_seed.len() > 1 {
                let base = &per_seed[0];
                let d: Vec<f64> = base
                    .iter()
                    .zip(&a)
                    .filter(|(x, y)| x.is_finite() && y.is_finite())
                    .map(|(x, y)| (x - y).abs())
                    .collect();
                let (_, _, p90, _) = stats::interdecile(&d);
                p90
            } else {
                f64::NAN
            };
            println!(
                "{region:>14} {seed:>5} {:>7} {ap50:>9.4} {aid:>9.4} {rp50:>9.4} {rp90:>9.4} {mv:>10.4}",
                r.len()
            );
        }
    }

    println!(
        "\nPart 1's columns have TRUE values of 1.0000 and 0.0000. Part 2's do not -- there is no\n\
         known answer for the real signal, which is exactly why part 1 exists as a floor.\n\
         \n\
         `seed move` is the p90 of |alpha(seed k) - alpha(seed 0)| over the same quads: sampling\n\
         noise on the real signal, everything else fixed. Compare it against `sib p90`.\n\
         \n\
         If seed move ~ sib p90, the sibling range is measuring its own draw and the Sibling\n\
         policy is thresholding on sampling error. If seed move << sib p90, the siblings\n\
         genuinely disagree and the statistic carries signal -- though `alpha` is a chaotic\n\
         quantity, so part of any residual is divergence that no ensemble size reduces.\n\
         \n\
         Both are compared against the shipped `sib_tau = 0.5`."
    );
}
