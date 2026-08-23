//! **Why low `E` behaves the opposite way to low `N`, measured rather than inferred.**
//!
//! The prediction on record was that low `E` biases toward *refine*, as low `N` did. The sweep
//! says otherwise: without the veto, near-field's leaf count **rises** with `E` (742 -> 2713 ->
//! 3463 at `E+1 = 2, 4, 8`). Low `E` **under**-refines.
//!
//! The candidate mechanism is a small-sample bias, not noise. `spread_shape` is the mean distance
//! of the copies' `shape_vec` from their centroid: with two points the centroid sits exactly
//! between them, so the statistic measures half the separation of a single pair, and it is
//! systematically **smaller** than the same quantity over eight copies. A low spread falls below
//! `tau_display` and the quad is *kept*.
//!
//! That is the opposite failure direction from `N`. `N` controls how well a quad knows its own
//! **area**, and undersampling it inflates the between-footprint variation that drives `alpha`.
//! `E` controls how well a footprint knows its own **value**, and undersampling it deflates the
//! within-footprint spread that is compared against `tau`. **One over-refines, the other
//! under-refines, and they are not interchangeable knobs.**
//!
//! This measures the bias directly, on identical footprints, so the explanation is a number.
//!
//! Run: `cargo run --release --example spread_bias_e [n]`

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::render::{self, Precision};

fn q(v: &mut Vec<f64>, f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * f).round() as usize]
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(48);
    println!("uniform {n}x{n} renders, identical nominal footprints, only E+1 varying.");
    println!("Halton offsets are a fixed prefix, so the E+1 = 2 copies are a SUBSET of the");
    println!("E+1 = 32 copies. Any movement is the estimator, not a different sample.\n");
    println!("{:>14} {:>5} {:>12} {:>12} {:>12} {:>12} {:>10}",
             "region", "E+1", "spread p10", "spread med", "spread p90", "vs E+1=32", "wall_s");

    for region in ["near-field", "deep interior", "far"] {
        let slice = grid::region(region, n, n, 0.05).unwrap();
        let mut ref_med = f64::NAN;
        let mut rows = Vec::new();
        for e1 in [2usize, 4, 8, 16, 32] {
            let cfg = EnsembleCfg { n_extra: e1 - 1, refine_flagged: false, ..Default::default() };
            let t0 = std::time::Instant::now();
            let px: Vec<PixelOut> = render::render(&slice, &cfg, Precision::F64);
            let mut sp: Vec<f64> =
                px.iter().map(|p| p.ensemble_spread).filter(|x| x.is_finite()).collect();
            let (a, m, b) = (q(&mut sp.clone(), 0.1), q(&mut sp.clone(), 0.5), q(&mut sp, 0.9));
            rows.push((e1, a, m, b, t0.elapsed().as_secs_f64()));
            ref_med = m;
        }
        for (e1, a, m, b, w) in rows {
            println!("{:>14} {:>5} {:>12.4e} {:>12.4e} {:>12.4e} {:>12.4} {:>10.1}",
                     region, e1, a, m, b, m / ref_med, w);
        }
        println!();
    }

    println!("If the median spread rises monotonically with E on identical footprints, the low-E");
    println!("under-refinement is a SMALL-SAMPLE BIAS in the estimator, not sampling noise, and");
    println!("more copies is the only fix. If it is flat, the explanation is wrong and the leaf-");
    println!("count trend needs another one.");
}
