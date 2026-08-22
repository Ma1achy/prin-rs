//! Arithmetic of the robust statistics, including the conventions that matter.

use prin_rs::ensemble::stats;

#[test]
fn median_uses_the_even_count_convention() {
    assert_eq!(stats::median(&[3.0f64, 1.0, 2.0]), 2.0);
    // Even count: the mean of the two central order statistics, matching numpy.
    assert_eq!(stats::median(&[1.0f64, 2.0, 3.0, 4.0]), 2.5);
    assert_eq!(stats::median(&[4.0f64, 1.0, 3.0, 2.0]), 2.5);
}

#[test]
fn mad_is_the_scaled_median_absolute_deviation() {
    // median = 3; deviations 2,1,0,1,2; median deviation = 1; MAD = 1.4826
    let v = [1.0f64, 2.0, 3.0, 4.0, 5.0];
    assert!((stats::mad(&v) - 1.4826).abs() < 1e-12, "{}", stats::mad(&v));
    // Constant data has zero spread.
    assert_eq!(stats::mad(&[7.0f64; 8]), 0.0);
}

/// The reason BRIEF §4 specifies MAD rather than a standard deviation: a std returns NaN the
/// moment one copy is non-finite — precisely the pathological pixel it exists to flag.
///
/// A non-finite value is an *infinitely bad* outcome, not missing data, so the copy keeps its
/// slot in the ordering and pushes the median. That is how "never discard a copy" survives
/// contact with a median. Fewer than half non-finite and the statistic is still meaningful.
#[test]
fn mad_survives_a_non_finite_copy_where_a_std_would_not() {
    let v = [1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, f64::NAN];

    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    let std = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
    assert!(std.is_nan(), "the std was expected to be poisoned");

    let m = stats::mad(&v);
    println!("with one NaN of 8 copies: std = {std}, MAD = {m}");
    assert!(m.is_finite() && m > 0.0, "MAD was poisoned too: {m}");

    // It does NOT shift the statistic, and that is the point rather than a defect:
    // insensitivity to a single extreme value is what robustness means.
    let clean = [1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    assert_eq!(m, stats::mad(&clean), "one bad copy of eight should not move a robust estimator");

    // But the copy is COUNTED, not discarded — it holds its slot in the ordering. Once a
    // majority is non-finite the statistic responds, which is the correct behaviour: at that
    // point the pixel genuinely is undetermined.
    let mostly_bad = [1.0f64, 2.0, 3.0, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN];
    let mb = stats::mad(&mostly_bad);
    println!("with five NaN of 8 copies: MAD = {mb}");
    assert!(!mb.is_finite(), "a majority-bad ensemble should not report a finite spread: {mb}");
}

#[test]
fn spread_event_is_normalised_to_one_when_maximally_split() {
    // All agree.
    assert_eq!(stats::spread_event::<f64>(&[1, 1, 1, 1]), 0.0);
    // Maximally split over 4 copies: modal count 1, raw fraction 3/4, normalised by
    // (1 - 1/4) = 3/4, giving exactly 1.
    let s = stats::spread_event::<f64>(&[0, 1, 2, 3]);
    assert!((s - 1.0).abs() < 1e-15, "{s}");
    // Half and half.
    let s = stats::spread_event::<f64>(&[0, 0, 1, 1]);
    assert!((s - (0.5 / 0.75)).abs() < 1e-15, "{s}");
}

/// `spread_event` is `refine_test.disagree` divided by `1 - 1/(E+1)`. That relation is the
/// only reference check available for this field.
#[test]
fn spread_event_matches_the_reference_statistic_up_to_its_normalisation() {
    for classes in [
        vec![0u8, 0, 1, 1, 2, 3, 3, 3],
        vec![3u8, 3, 3, 3, 3, 3, 3, 0],
        vec![0u8, 1, 2, 3, 0, 1, 2, 3],
    ] {
        let n = classes.len() as f64;
        let mut best = 0usize;
        for c in &classes {
            best = best.max(classes.iter().filter(|x| *x == c).count());
        }
        let disagree = 1.0 - best as f64 / n; // refine_test.disagree
        let want = disagree / (1.0 - 1.0 / n);
        let got = stats::spread_event::<f64>(&classes);
        assert!((got - want).abs() < 1e-15, "{got} vs {want}");
    }
}
