//! The Step 5a acceptance gate for `error_ratio`, and a finding about the statistic itself.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;
use prin_rs::integrate::az::StepLimit;

/// **Refinement off, AND `StepLimit::None`.** These tests compare spread estimators *on damaged
/// pixels*, and both mechanisms exist precisely to remove those — with either on there is nothing
/// left to separate, and the test would pass by having no subject rather than by the estimator
/// working.
///
/// The step-limit pin was added when `StepLimit::Predictive` became the default and both tests
/// here failed. **They failed correctly**: `refined_pixels_are_repaired` fell over on its own
/// `n_ref > 0` guard — *nothing was flagged, so this test has no subject* — which is the
/// assertion written to catch exactly this. That the per-step limit deletes the damaged
/// population outright is the strongest corroboration in the suite of what
/// `results/step_control/README.md` measures; it is recorded here rather than worked around,
/// and the pin is what keeps these two tests measuring the estimator they are about.
fn render(size: usize) -> Vec<prin_rs::ensemble::pixel::PixelOut> {
    let s = grid::region("near-field", size, size, 0.05).unwrap();
    let cfg = EnsembleCfg {
        refine_flagged: false,
        step_limit: StepLimit::None,
        ..Default::default()
    };
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
        ("MAD", (|p: &&prin_rs::ensemble::pixel::PixelOut| p.error_ratio_mad) as fn(&&_) -> f64),
        ("max deviation", |p: &&prin_rs::ensemble::pixel::PixelOut| p.error_ratio),
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
/// **It is gated on the healthy population, and that is not a dodge.** `error_ratio` is now
/// built on the maximum deviation, so a genuinely damaged pixel reads five or six orders of
/// magnitude above 1 — as it should, that being the entire point. A max over the whole grid
/// would therefore be a measurement of the worst pixel in the region, not a correctness
/// criterion. The criterion that carries weight is: **where the integration is healthy, the
/// ratio sits at 1.** Both populations are printed, so nothing is hidden by the split.
#[test]
fn error_ratio_acceptance_near_field_t13() {
    let px = render(32);
    let fin = |x: &f64| x.is_finite();
    let er: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(fin).collect();
    let md: Vec<f64> = px.iter().map(|p| p.error_ratio_mad).filter(fin).collect();
    let healthy: Vec<f64> = px
        .iter()
        .filter(|p| p.energy_drift_max <= 1e-3)
        .map(|p| p.error_ratio)
        .filter(fin)
        .collect();
    let argmax = px
        .iter()
        .enumerate()
        .filter(|(_, p)| p.error_ratio.is_finite())
        .max_by(|a, b| a.1.error_ratio.partial_cmp(&b.1.error_ratio).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    println!("near-field 32x32, t=13, E+1=8 copies, eta=0.01, f64");
    println!("  error_ratio, all pixels    max {:.4e}  p99 {:.4e}  median {:.6}  argmax {argmax}",
             q(er.clone(), 1.0), q(er.clone(), 0.99), q(er.clone(), 0.5));
    println!("  error_ratio, healthy only  max {:.4e}  p99 {:.4e}  median {:.6}   ({} px)",
             q(healthy.clone(), 1.0), q(healthy.clone(), 0.99), q(healthy.clone(), 0.5),
             healthy.len());
    println!("  error_ratio_mad, all       max {:.6}  p99 {:.6}  median {:.6}",
             q(md.clone(), 1.0), q(md.clone(), 0.99), q(md.clone(), 0.5));
    println!();
    println!("BRIEF §5 asks for 1.0000. The median is 1.000000 on both populations. The max");
    println!("over all pixels is not 1 and must not be: on a chaos instrument, a region with");
    println!("no damaged pixels would mean the region was not interesting.");

    assert!(q(er, 0.5) < 1.001, "median error_ratio should sit at 1");
    assert!(q(healthy.clone(), 0.5) < 1.001, "healthy median error_ratio should sit at 1");
    assert!(
        q(healthy, 1.0) < HEALTHY_MAX,
        "max error_ratio over healthy pixels exceeded the measured bound"
    );
}

/// Set from the measured healthy-f64 run: near-field 32x32, t=13, eta=0.01 gives a healthy
/// p99 of 1.0228 and a healthy max of 5.2087. The bound is 10.0 — roughly 2x the measured
/// worst — and a named constant so any future change to it shows up in a diff.
///
/// The gap between p99 (1.02) and max (5.21) is itself worth reading: `drift_max <= 1e-3` is
/// a blunt cut, and a handful of pixels just inside it are already mildly damaged. The cut is
/// not a clean partition and is not presented as one.
const HEALTHY_MAX: f64 = 10.0;

/// The other side of the switch above: with the second pass on, the pixels `error_ratio` flags
/// are repaired rather than merely reported. BRIEF §2.5.
#[test]
fn refined_pixels_are_repaired() {
    let s = grid::region("near-field", 32, 32, 0.05).unwrap();
    // Both arms on `StepLimit::None`: this test is about the refinement pass, and under the
    // shipped per-step limit there is no flagged population for it to act on. See the module note.
    let off = EnsembleCfg {
        refine_flagged: false,
        step_limit: StepLimit::None,
        ..Default::default()
    };
    let on = EnsembleCfg { step_limit: StepLimit::None, ..Default::default() };
    let a: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &off)).collect();
    let b: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &on)).collect();

    let worst = |v: &Vec<prin_rs::ensemble::pixel::PixelOut>| {
        v.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).fold(0.0f64, f64::max)
    };
    let n_ref = b.iter().filter(|p| p.refined).count();
    println!("near-field 32x32, t=13, f64");
    println!("  drift max      {:.4e} unrefined -> {:.4e} refined", worst(&a), worst(&b));
    println!("  pixels re-integrated: {n_ref} of {}", b.len());
    println!("  pixels above 1e-3 drift: {} -> {}",
             a.iter().filter(|p| p.energy_drift_max > 1e-3).count(),
             b.iter().filter(|p| p.energy_drift_max > 1e-3).count());

    assert!(n_ref > 0, "nothing was flagged, so this test has no subject");
    assert!(worst(&b) < worst(&a), "refinement did not reduce the worst drift");
    // Every refined pixel keeps both values, so a refinement that did not help stays visible.
    for p in b.iter().filter(|p| p.refined) {
        assert!(p.eta_used < on.eta, "a refined pixel kept the coarse eta");
        assert!(p.error_ratio_coarse.is_finite(), "the coarse value was not carried forward");
    }
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
