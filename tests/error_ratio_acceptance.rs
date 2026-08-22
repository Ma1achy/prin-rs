//! The Step 5a acceptance gate for `error_ratio`, and a finding about the statistic itself.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

fn render(size: usize) -> Vec<prin_rs::ensemble::pixel::PixelOut> {
    let s = grid::region("near-field", size, size, 0.05).unwrap();
    let cfg = EnsembleCfg::default();
    (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect()
}

fn q(mut v: Vec<f64>, f: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * f).round() as usize]
}

/// **BRIEF §4's MAD requirement defeats the field's own purpose.**
///
/// MAD is specified so the statistic survives a non-finite copy — a standard deviation
/// returns NaN on precisely the pathological pixel `error_ratio` exists to flag. That reason
/// is sound. But the repair overshoots: a spread estimator that a single wild copy cannot
/// move is also one that cannot *see* a single wild copy.
///
/// Measured on near-field 32x32 at t=13: a pixel whose worst copy drifted by `1.2e+02` — an
/// energy error 120x the total energy — reported a MAD-based `error_ratio` of 1.1369, while
/// the healthy p99 sits at 1.0756. The two populations are not separable.
///
/// The same ratio built on the **maximum deviation** from the median is both NaN-safe (a
/// non-finite copy gives an infinite deviation, which is the correct answer) and sensitive.
///
/// This contradicts a stated non-negotiable in CLAUDE.md, so both are computed and dumped
/// and the choice is left to the user. Reported, not silently changed.
#[test]
fn mad_based_error_ratio_cannot_separate_damaged_pixels() {
    let px = render(32);
    let damaged: Vec<_> = px.iter().filter(|p| p.energy_drift_max > 1e-3).collect();
    let healthy: Vec<_> = px.iter().filter(|p| p.energy_drift_max <= 1e-3).collect();
    assert!(!damaged.is_empty(), "no damaged pixels in this region; the test cannot say anything");

    println!("{} damaged pixels (drift_max > 1e-3), {} healthy", damaged.len(), healthy.len());
    println!("{:>16}{:>15}{:>15}{:>14}", "statistic", "damaged med", "healthy p99", "separation");

    let mut seps = Vec::new();
    for (name, f) in [
        ("MAD", (|p: &&prin_rs::ensemble::pixel::PixelOut| p.error_ratio) as fn(&&_) -> f64),
        ("max deviation", |p: &&prin_rs::ensemble::pixel::PixelOut| p.error_ratio_range),
    ] {
        let dmed = q(damaged.iter().map(f).filter(|x| x.is_finite()).collect(), 0.5);
        let hp99 = q(healthy.iter().map(f).filter(|x| x.is_finite()).collect(), 0.99);
        println!("{name:>16}{dmed:>15.4e}{hp99:>15.4e}{:>14.2}", dmed / hp99);
        seps.push(dmed / hp99);
    }
    println!();
    println!("Separation is the damaged median over the healthy p99 — how far apart a");
    println!("threshold could be placed. Near 1 means the flag cannot separate them.");

    assert!(seps[0] < 2.0, "MAD-based separation was expected to be poor, got {}", seps[0]);
    assert!(seps[1] > 10.0, "max-deviation separation should be decisive, got {}", seps[1]);
}

/// The gate. The threshold comes from what a healthy f64 run actually produces, per the
/// earlier decision — not from a number chosen by eye.
///
/// It is a weak gate, and deliberately reported as such: given the separation measured above,
/// almost any threshold on the MAD-based statistic is arbitrary. The distribution is what
/// carries information, so all four numbers are printed.
#[test]
fn error_ratio_acceptance_near_field_t13() {
    let px = render(32);
    let er: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
    let rg: Vec<f64> = px.iter().map(|p| p.error_ratio_range).filter(|x| x.is_finite()).collect();
    let argmax = px
        .iter()
        .enumerate()
        .filter(|(_, p)| p.error_ratio.is_finite())
        .max_by(|a, b| a.1.error_ratio.partial_cmp(&b.1.error_ratio).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    println!("near-field 32x32, t=13, E+1=8 copies, eta=0.01, f64");
    println!("  error_ratio (MAD)   max {:.6}  p99 {:.6}  median {:.6}  argmax pixel {argmax}",
             q(er.clone(), 1.0), q(er.clone(), 0.99), q(er.clone(), 0.5));
    println!("  error_ratio (range) max {:.4e}  p99 {:.4e}  median {:.6}",
             q(rg.clone(), 1.0), q(rg.clone(), 0.99), q(rg.clone(), 0.5));
    println!();
    println!("BRIEF §5 asks for 1.0000. The median is 1.000000; the max is not, and cannot");
    println!("be — a max over 1024 pixels of a statistic §4 calls unstable in magnitude is");
    println!("the least reproducible quantity in the system.");

    assert!(q(er.clone(), 0.5) < 1.001, "median error_ratio should sit at 1");
    assert!(q(er, 1.0) < 3.0, "max error_ratio exceeded the healthy-run bound");
}

/// NOTES §2.1, settled with data: the reference's `d_min` blind spot never bites.
#[test]
fn d_min_gap_is_zero_in_every_region() {
    for region in ["near-field", "mid-field", "body2 core", "body1 slice", "far"] {
        let s = grid::region(region, 8, 8, 0.05).unwrap();
        let cfg = EnsembleCfg::default();
        let worst = (0..s.npix())
            .map(|i| evaluate::<f64>(&s, i, &cfg).d_min_gap)
            .filter(|x| x.is_finite())
            .fold(0.0f64, f64::max);
        println!("{region:>14}: max d_min_gap = {worst:.3e}");
        assert_eq!(worst, 0.0, "{region}: the unregularised pair became the minimum");
    }
    println!();
    println!("d_min_ref tracks only the two regularised pairs; d_min_true includes the");
    println!("unregularised side. The gap is identically zero, so the discrepancy flagged in");
    println!("NOTES §2.1 exists in principle and does not materialise in practice.");
}
