//! Does the fixed offset prefix change the f32 answer, and if so why?
//!
//! Switching to the spec's fixed Halton (2,3) prefix took the f32 `spread_shape` from 0.5% off
//! the f64 truth to **30%** off. That is a regression, not a tolerance to widen, and this is the
//! measurement that says what it is.
//!
//! The hypothesis: with a fixed prefix, every footprint uses the *same* offsets, so if one of
//! them places a configuration near the Levi-Civita branch cut it does so in every pixel at
//! once. Under PCG the offsets varied per footprint, so the same effect averaged out.

use prin_rs::ensemble::jitter::Scheme;
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

fn mean_spread<T: prin_rs::Real>(region: &str, lc_stable: bool, scheme: Scheme) -> (f64, usize) {
    let s = grid::region(region, 5, 5, 0.05).unwrap();
    let cfg = EnsembleCfg {
        t_max: 13.0,
        lc_stable,
        jitter_scheme: scheme,
        refine_flagged: false,
        ..Default::default()
    };
    let mut acc = 0.0;
    let mut bad = 0usize;
    for i in 0..s.npix() {
        let p = evaluate::<T>(&s, i, &cfg);
        if p.spread_shape.is_finite() {
            acc += p.spread_shape;
        } else {
            bad += 1;
        }
    }
    (acc / s.npix() as f64, bad)
}

fn main() {
    println!("=== the f32 spread_shape answer, by offset scheme ===");
    println!("{:>14}{:>10}{:>14}{:>14}{:>14}{:>10}",
             "region", "scheme", "f64 truth", "f32 stable", "rel err", "NaN pix");
    for region in ["near-field", "body2 core", "body1 slice", "mid-field"] {
        for scheme in [Scheme::Halton, Scheme::Pcg] {
            let (truth, _) = mean_spread::<f64>(region, true, scheme);
            let (s32, bad) = mean_spread::<f32>(region, true, scheme);
            println!("{region:>14}{:>10}{truth:>14.4e}{s32:>14.4e}{:>14.2}{bad:>10}",
                     format!("{scheme:?}"), (s32 - truth).abs() / truth);
        }
    }

    println!();
    println!("The f64 truths differ between schemes by up to 1.8x, which is expected: different");
    println!("offsets measure a different ensemble. What matters is whether f32 tracks f64.");

    println!();
    println!("=== per pixel, not per mean — near-field, f32 against f64 ===");
    println!("A mean over 25 pixels can be moved by one of them, so the distribution is what");
    println!("says whether this is a systematic loss or a single bad pixel.");
    println!();
    let s = grid::region("near-field", 5, 5, 0.05).unwrap();
    println!("{:>10}{:>12}{:>12}{:>12}{:>12}{:>14}",
             "scheme", "rel p50", "rel p90", "rel max", "argmax px", "worst f64/f32");
    for scheme in [Scheme::Halton, Scheme::Pcg] {
        let cfg = |lc_stable: bool| EnsembleCfg {
            t_max: 13.0,
            lc_stable,
            jitter_scheme: scheme,
            refine_flagged: false,
            ..Default::default()
        };
        let mut rel: Vec<(f64, usize, f64, f64)> = (0..s.npix())
            .map(|i| {
                let a = evaluate::<f64>(&s, i, &cfg(true)).spread_shape;
                let b = evaluate::<f32>(&s, i, &cfg(true)).spread_shape;
                ((b - a).abs() / a.abs().max(1e-300), i, a, b)
            })
            .collect();
        rel.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        let q = |f: f64| rel[(((rel.len() - 1) as f64) * f).round() as usize];
        let w = rel[rel.len() - 1];
        println!("{:>10}{:>12.4}{:>12.4}{:>12.4}{:>12}{:>14}",
                 format!("{scheme:?}"), q(0.5).0, q(0.9).0, w.0, w.1,
                 format!("{:.3e}/{:.3e}", w.2, w.3));
    }
    println!();
    println!("If the median is small and only the max is large, this is one pixel where the two");
    println!("precisions took different branches at a close approach - a chaotic divergence, not");
    println!("a systematic loss of accuracy from the offset scheme.");
}
