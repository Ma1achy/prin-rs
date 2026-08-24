//! Stage 1: the between-footprint arm, the matched-count controls, and the hot-set layout.
//!
//! Every assertion here is written to be able to fail. Where a property is *structural* it is
//! stated as structural in a comment and not dressed as a measurement.

use prin_rs::ensemble::pixel::PixelOut;
use prin_rs::quad::{Criterion, QuadReduction};
use prin_rs::scheduler::reduce;
use prin_rs::spatial::{self, Layout};
use prin_rs::stats;

// ---------------------------------------------------------------------------------------
// §3.2 — the hot-set layout
// ---------------------------------------------------------------------------------------

fn mask(n: usize, cells: &[(usize, usize)]) -> Vec<bool> {
    let mut m = vec![false; n * n];
    for &(jx, jy) in cells {
        m[jy * n + jx] = true;
    }
    m
}

#[test]
fn an_empty_hot_set_reports_nan_and_not_zero() {
    let l = spatial::layout(&vec![false; 64], 8);
    assert_eq!(l.n_hot, 0);
    assert_eq!(l.n_components, 0);
    // NaN, deliberately. A perimeter ratio of 0 is what a *featureless fully hot* quad reads,
    // and "nothing is hot" must not be arithmetically indistinguishable from it.
    assert!(l.perimeter_ratio.is_nan(), "empty must be NaN, got {}", l.perimeter_ratio);
}

#[test]
fn a_fully_hot_quad_has_no_perimeter_and_one_component() {
    let l = spatial::layout(&vec![true; 64], 8);
    assert_eq!(l.n_hot, 64);
    assert_eq!(l.n_components, 1);
    assert_eq!(l.largest_component, 64);
    // Internal edges only: a uniformly hot quad has no hot/cold edge anywhere inside it.
    assert_eq!(l.perimeter_ratio, 0.0);
}

#[test]
fn an_isolated_cell_reads_four_and_a_filament_reads_about_two() {
    let one = spatial::layout(&mask(8, &[(3, 3)]), 8);
    assert_eq!((one.n_hot, one.n_components, one.largest_component), (1, 1, 1));
    assert_eq!(one.perimeter_ratio, 4.0);

    // A straight one-cell-wide filament spanning the quad: 2 cold neighbours per cell, and the
    // two ends touch the border, which is not counted.
    let row: Vec<(usize, usize)> = (0..8).map(|jx| (jx, 4)).collect();
    let fil = spatial::layout(&mask(8, &row), 8);
    assert_eq!((fil.n_hot, fil.n_components, fil.largest_component), (8, 1, 8));
    assert_eq!(fil.perimeter_ratio, 2.0);

    // A compact blob of the same count is markedly less thin — the separation the signal rests
    // on, asserted rather than asserted-about.
    let blob: Vec<(usize, usize)> =
        (2..5).flat_map(|jx| (2..5).map(move |jy| (jx, jy))).take(8).collect();
    let b = spatial::layout(&mask(8, &blob), 8);
    assert!(
        b.perimeter_ratio < fil.perimeter_ratio,
        "blob {} should be less thin than filament {}",
        b.perimeter_ratio,
        fil.perimeter_ratio
    );
}

#[test]
fn two_blobs_are_two_components_and_a_checkerboard_is_scatter() {
    let two = spatial::layout(&mask(8, &[(1, 1), (2, 1), (6, 6), (6, 5)]), 8);
    assert_eq!(two.n_components, 2);
    assert_eq!(two.largest_component, 2);

    // 4-connectivity, chosen so that diagonal-only contact is *not* connection: a checkerboard
    // is scatter, which is the reading that separates chaos from a boundary.
    let checks: Vec<(usize, usize)> =
        (0..8).flat_map(|jx| (0..8).map(move |jy| (jx, jy))).filter(|(x, y)| (x + y) % 2 == 0).collect();
    let c = spatial::layout(&mask(8, &checks), 8);
    assert_eq!(c.n_hot, 32);
    assert_eq!(c.n_components, 32, "a checkerboard must read as 32 separate components");
    assert_eq!(c.largest_component, 1);
}

#[test]
fn a_strictly_diagonal_mask_reads_as_scatter_and_that_is_a_known_limit() {
    // **A negative result, asserted so it cannot drift silently.** Under 4-connectivity an
    // exactly-diagonal one-cell filament is n_components == length, i.e. indistinguishable from
    // scatter by the connectivity test alone.
    //
    // It is left this way on purpose. The hot mask is "footprints whose cell straddles the
    // boundary", which for a generic line is the *supercover* of that line and is 4-connected;
    // the exact 45-degree single-cell diagonal is the measure-zero case that a real basin
    // boundary does not produce. Switching to 8-connectivity would fix this case and break the
    // checkerboard, which is the case that actually occurs — chaos.
    let diag: Vec<(usize, usize)> = (0..8).map(|k| (k, k)).collect();
    let d = spatial::layout(&mask(8, &diag), 8);
    assert_eq!(d.n_hot, 8);
    assert_eq!(d.n_components, 8);
    assert_eq!(d.largest_component, 1);
}

#[test]
fn looks_like_a_boundary_separates_the_filament_from_the_checkerboard() {
    let row: Vec<(usize, usize)> = (0..8).map(|jx| (jx, 4)).collect();
    let fil = spatial::layout(&mask(8, &row), 8);
    let checks: Vec<(usize, usize)> =
        (0..8).flat_map(|jx| (0..8).map(move |jy| (jx, jy))).filter(|(x, y)| (x + y) % 2 == 0).collect();
    let c = spatial::layout(&mask(8, &checks), 8);

    assert!(fil.looks_like_boundary(8, 1.5), "a filament is a boundary");
    assert!(!c.looks_like_boundary(8, 1.5), "a checkerboard is not");
    assert!(!Layout::default().looks_like_boundary(8, 1.5), "nothing hot is not a boundary");
}

// ---------------------------------------------------------------------------------------
// §1 — the between arm
// ---------------------------------------------------------------------------------------

/// A footprint with a chosen nominal shape vector, event class and terminal outcome.
fn fp(shape: [f64; 3], event_class: u8, outcome: u8) -> PixelOut {
    PixelOut {
        shape_vec: shape,
        event_class,
        outcome,
        state: outcome >> 2,
        detail: outcome & 3,
        ensemble_spread: 0.0,
        spread_shape: 0.0,
        spread_event: 0.0,
        error_ratio: 1.0,
        t_end: 13.0,
        censored: true,
        ..Default::default()
    }
}

#[test]
fn between_event_reads_the_event_class_and_never_the_terminal_outcome() {
    // The regression guard for a non-negotiable. `spread_event` is over the EVENT CLASS; the
    // terminal outcome is terminal-grain and inverts under lockstep, so a between-footprint arm
    // built on `outcome` would reinstate exactly what was removed, one level up.
    //
    // Sixteen footprints whose terminal outcomes are maximally split but whose event classes
    // all agree. If the arm read `outcome`, this would be ~1.0.
    let px: Vec<PixelOut> = (0..16).map(|k| fp([1.0, 0.0, 0.0], 2, k as u8)).collect();
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert_eq!(r.between_event, 0.0, "agreeing event classes must give zero disagreement");

    // And the converse: identical terminal outcomes, maximally split event classes.
    let px: Vec<PixelOut> = (0..16).map(|k| fp([1.0, 0.0, 0.0], k as u8, 7)).collect();
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert!(r.between_event > 0.99, "split event classes must fire, got {}", r.between_event);
}

#[test]
fn between_shape_is_the_same_estimator_over_the_nominals() {
    // Two clusters on the sphere, eight footprints each: the clean-boundary case §1.3 predicts
    // the within arm is blind to. The between arm must see it.
    let a = [1.0, 0.0, 0.0];
    let b = [-1.0, 0.0, 0.0];
    let px: Vec<PixelOut> = (0..16)
        .map(|k| fp(if k < 8 { a } else { b }, 1, 4))
        .collect();
    let r = reduce(&px, 4, 1e-3, 13.0);

    // Every footprint is exactly 1 from the centroid at the origin, halved by the chord
    // convention.
    assert!((r.between_shape - 0.5).abs() < 1e-12, "got {}", r.between_shape);
    // The within arm is identically zero here by construction — each footprint's own copies
    // were never populated — which is the whole shape of §1.3's prediction.
    assert_eq!(r.spread_median, 0.0);
    assert!(r.between_spread > r.spread_median);
}

#[test]
fn between_matched_equals_between_full_when_the_counts_already_agree() {
    // The matched-count control must be the *same estimator*, differing only in how many
    // samples it is handed. Give it every footprint and it must return the full value exactly.
    let n = 4;
    let px: Vec<PixelOut> = (0..n * n)
        .map(|k| {
            let t = k as f64 * 0.1;
            let mut p = fp([t.cos(), t.sin(), 0.0], 1, 4);
            // `E+1` is read from the copies the footprints carry; hand it N^2 of them.
            p.copy_shapes = vec![[t.cos(), t.sin(), 0.0]; n * n];
            p
        })
        .collect();
    let r = reduce(&px, n, 1e-3, 13.0);
    assert_eq!(
        r.between_matched, r.between_shape,
        "matched at full count must be bitwise the full value"
    );
}

#[test]
fn within_pooled_is_nan_when_the_copies_were_not_kept() {
    // Reported as "not measured", never as zero. A zero here would read as perfect agreement
    // among copies that were never looked at.
    let px: Vec<PixelOut> = (0..16).map(|_| fp([1.0, 0.0, 0.0], 1, 4)).collect();
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert!(r.within_pooled.is_nan(), "got {}", r.within_pooled);
}

#[test]
fn a_collapsed_decode_is_undetermined_and_a_uniform_region_is_not() {
    // Tested on initial conditions, never on the spread being zero. A genuinely uniform region
    // has zero spread over perfectly distinct ICs, and flagging that as a numerical failure
    // would be reporting the physics as a bug.
    let mut r = QuadReduction { n_footprints: 64, n_distinct_ic: 64, between_shape: 0.0, ..Default::default() };
    assert!(!r.between_collapsed(), "zero spread over distinct ICs is a uniform region");

    r.n_distinct_ic = 63;
    assert!(r.between_collapsed(), "one repeated IC is a collapsed decode");

    // And the trap it exists to catch: a collapsed quad's spread is exactly zero, which the
    // criterion would otherwise read as maximum confidence.
    assert_eq!(r.signal(Criterion::Between, Default::default()), 0.0);
}

#[test]
fn the_escape_gradient_refuses_to_answer_when_nothing_escaped() {
    // At t_max = 13 zero of 1024 near-field pixels escape; 109 do at t_max = 20. A gradient
    // returned over an empty set would be a null that could not have failed.
    let px: Vec<PixelOut> = (0..16).map(|_| fp([1.0, 0.0, 0.0], 1, 4)).collect();
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert_eq!(r.escaped_fraction, 0.0);
    assert!(r.t_end_gradient.is_nan(), "must decline, got {}", r.t_end_gradient);

    // With escapes present it answers, and the value is the mean adjacent difference.
    let px: Vec<PixelOut> = (0..16)
        .map(|k| {
            let mut p = fp([1.0, 0.0, 0.0], 1, 4);
            p.censored = false;
            p.t_end = 1.0 + (k % 4) as f64; // varies along x only
            p
        })
        .collect();
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert_eq!(r.escaped_fraction, 1.0);
    // 12 x-pairs differ by 1, 12 y-pairs differ by 0.
    assert!((r.t_end_gradient - 0.5).abs() < 1e-12, "got {}", r.t_end_gradient);
}

#[test]
fn a_non_finite_footprint_is_hot_and_not_calm() {
    // "Never discard a copy" at quad level: a footprint that could not be determined is not
    // evidence of calm, and treating it as cold would hide it from the one statistic built to
    // find structure.
    let mut px: Vec<PixelOut> = (0..16).map(|_| fp([1.0, 0.0, 0.0], 1, 4)).collect();
    px[5].ensemble_spread = f64::NAN;
    let r = reduce(&px, 4, 1e-3, 13.0);
    assert_eq!(r.layout_within.n_hot, 1, "the non-finite footprint must count as hot");
}

// ---------------------------------------------------------------------------------------
// Rank statistics
// ---------------------------------------------------------------------------------------

#[test]
fn spearman_declines_rather_than_returning_zero_on_a_dead_input() {
    // A difference can be small because both sides are right or because both are dead. A
    // constant input has no ordering, and 0.0 would read as "no relationship".
    assert!(stats::spearman(&[1.0, 1.0, 1.0, 1.0], &[1.0, 2.0, 3.0, 4.0]).is_nan());
    assert!(stats::spearman(&[1.0, 2.0], &[1.0, 2.0]).is_nan(), "under three points");

    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [5.0, 4.0, 3.0, 2.0, 1.0];
    assert!((stats::spearman(&x, &x) - 1.0).abs() < 1e-12);
    assert!((stats::spearman(&x, &y) + 1.0).abs() < 1e-12);
}

#[test]
fn rank_displacement_shows_movement_a_high_rho_conceals() {
    // The guard against "never conclude no effect from an aggregate". Two orderings agreeing
    // strongly overall, with one item crossing the whole list.
    let x: Vec<f64> = (0..50).map(|k| k as f64).collect();
    let mut y = x.clone();
    y[0] = 100.0; // the first item moves to last
    let rho = stats::spearman(&x, &y);
    let d = stats::rank_displacement(&x, &y);
    assert!(rho > 0.85, "rho should still read high, got {rho}");
    let worst = d.iter().cloned().fold(0.0f64, f64::max);
    assert!(worst > 0.9, "and one item moved almost the whole list: {worst}");
}
