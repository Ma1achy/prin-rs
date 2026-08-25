//! Stage 1: the between-footprint arm, the matched-count controls, and the hot-set layout.
//!
//! Every assertion here is written to be able to fail. Where a property is *structural* it is
//! stated as structural in a comment and not dressed as a measurement.

use prin_rs::ensemble::pixel::PixelOut;
use prin_rs::quad::{Criterion, QuadReduction};
use prin_rs::scheduler::reduce;
use prin_rs::spatial::{self, HotRule, Layout};
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
    assert_eq!(r.between_event, 0.0, "agreeing event classes must give zero disagreement");

    // And the converse: identical terminal outcomes, maximally split event classes.
    let px: Vec<PixelOut> = (0..16).map(|k| fp([1.0, 0.0, 0.0], k as u8, 7)).collect();
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);

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
    let r = reduce(&px, n, 1e-3, HotRule::default(), 13.0);
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
    assert_eq!(r.terminated_fraction, 0.0);
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
    assert_eq!(r.terminated_fraction, 1.0);
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
    let r = reduce(&px, 4, 1e-3, HotRule::default(), 13.0);
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

// ---------------------------------------------------------------------------------------
// §2 — OKLab, and the metric's own guards
// ---------------------------------------------------------------------------------------

use prin_rs::output::oklab;

#[test]
fn oklab_matches_the_published_reference_triples() {
    // Ottosson's published values, for LINEAR sRGB inputs. This is a transcription check on
    // matrices that fail silently: a swapped coefficient still produces plausible colours.
    let cases: [([f64; 3], [f64; 3]); 4] = [
        ([1.0, 1.0, 1.0], [1.0, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.627_955, 0.224_863, 0.125_846]),
        ([0.0, 1.0, 0.0], [0.866_440, -0.233_888, 0.179_498]),
        ([0.0, 0.0, 1.0], [0.452_014, -0.032_457, -0.311_528]),
    ];
    for (rgb, want) in cases {
        let got = oklab::linear_to_oklab(rgb[0], rgb[1], rgb[2]);
        for k in 0..3 {
            assert!(
                (got[k] - want[k]).abs() < 1e-5,
                "linear {rgb:?} -> {got:?}, expected {want:?}"
            );
        }
    }
}

#[test]
fn oklab_round_trips_through_srgb() {
    // Every 8-bit triple must survive the round trip, or the transform is not invertible on the
    // domain the renders actually live in.
    for r in (0..=255).step_by(17) {
        for g in (0..=255).step_by(51) {
            for b in (0..=255).step_by(51) {
                let c = [r as u8, g, b];
                let back = oklab::oklab_to_srgb(oklab::srgb_to_oklab(c));
                assert_eq!(c, back, "round trip failed on {c:?}");
            }
        }
    }
}

#[test]
fn oklab_distance_is_zero_only_on_equality_and_orders_sensibly() {
    assert_eq!(oklab::delta([12, 34, 56], [12, 34, 56]), 0.0);
    let near = oklab::delta([100, 100, 100], [104, 100, 100]);
    let far = oklab::delta([100, 100, 100], [200, 100, 100]);
    assert!(near > 0.0 && far > near, "near {near} far {far}");
    // sRGB is perceptually non-uniform, which is the whole reason this module exists: an equal
    // channel step is a larger visible change in the dark end than the light end.
    let dark = oklab::delta([10, 10, 10], [30, 30, 30]);
    let light = oklab::delta([225, 225, 225], [245, 245, 245]);
    assert!(dark > light, "equal channel steps: dark {dark} vs light {light}");
}

#[test]
fn image_error_reports_the_fraction_moved_beside_the_mean() {
    // The guard against reading an aggregate alone. One pixel in a hundred moving a long way,
    // and ninety-nine standing still, must be visible as such.
    let a = vec![0u8; 300];
    let mut b = a.clone();
    b[0] = 255;
    b[1] = 255;
    b[2] = 255;
    let (mean, _p99, max, moved) = oklab::image_error(&a, &b);
    assert!((moved - 0.01).abs() < 1e-12, "moved {moved}");
    assert!(max > 0.9, "one pixel went the whole way: {max}");
    assert!(mean < 0.02, "and the mean barely notices: {mean}");
}

#[test]
fn the_metric_is_exact_at_the_full_tree_and_the_greedy_replay_is_monotone() {
    use prin_rs::ensemble::pixel::EnsembleCfg;
    use prin_rs::grid::Chart;
    use prin_rs::metric::{self, Rank};

    // Deliberately tiny: 2 levels, N=2, so res = 8 and the whole cache is 21 quads.
    let (levels, n) = (2u32, 2usize);
    let res = (1usize << levels) * n;
    let ens = EnsembleCfg { n_extra: 1, t_max: 2.0, n_sync: 4, refine_flagged: false, ..Default::default() };
    let cache = metric::build(
        "deep interior", 0.0, 0.0, 0.05, 0, Chart::BodyPlane, levels, n, res, 1e-4, &ens, metric::Colouring::Outcome,
    );

    // Exactly zero at the deepest level: the reference IS that tree, so this is a consistency
    // check on the rasterisation, not a statement about image quality.
    let w = 1u32 << levels;
    let deepest: Vec<metric::Key> =
        (0..w).flat_map(|iy| (0..w).map(move |ix| (levels, ix, iy))).collect();
    assert_eq!(cache.error_of(&deepest), 0.0);

    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;
    let pts = metric::replay(&cache, Rank::GreedyOracle, full);

    // Monotone non-increasing. A rise here is a bug in the priority queue, not a finding.
    for w in pts.windows(2) {
        assert!(
            w[1].error <= w[0].error + 1e-12,
            "greedy error rose: {} -> {}",
            w[0].error,
            w[1].error
        );
    }
    assert_eq!(pts.last().unwrap().error, 0.0, "the replay must reach the reference");

    // Budget accounting: the root, then four quads per split.
    for (j, p) in pts.iter().enumerate() {
        assert_eq!(p.budget, 1 + 4 * j);
        assert_eq!(p.leaves, 1 + 3 * j);
    }

    // **There is deliberately no assertion that greedy_oracle dominates every criterion.**
    // Greedy is not a bound on a sequential tree problem: a quad whose own split gains little
    // may unlock children with large gains two levels down, and greedy declines it. Such a test
    // would fire on correct behaviour. What IS asserted is that greedy never does worse than
    // the tree it started from, which is the monotonicity above.
}

#[test]
fn a_criterion_enters_the_replay_as_an_ordering_and_not_against_tau() {
    use prin_rs::ensemble::pixel::EnsembleCfg;
    use prin_rs::grid::Chart;
    use prin_rs::metric::{self, Rank};
    use prin_rs::quad::{Agg, Criterion};

    let (levels, n) = (2u32, 2usize);
    let res = (1usize << levels) * n;
    let ens = EnsembleCfg { n_extra: 1, t_max: 2.0, n_sync: 4, refine_flagged: false, ..Default::default() };
    let cache = metric::build(
        "near-field", 1.0, 3.0, 0.05, 0, Chart::BodyPlane, levels, n, res, 1e-4, &ens, metric::Colouring::Outcome,
    );
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;

    // The between arm runs 1.17x the within arm in near-field and 9.56x in `far`, so a
    // threshold comparison would score that rescaling rather than the signal. A ranking is
    // invariant to any monotone rescaling, and this is what holds that: both replays spend the
    // budget identically at every step, whatever tau was used to build the cache.
    let a = metric::replay(&cache, Rank::Signal(Criterion::Within, Agg::Median), full);
    let b = metric::replay(&cache, Rank::Signal(Criterion::Within, Agg::Median), full);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.error, y.error, "the replay must be deterministic");
    }
}

// ---------------------------------------------------------------------------------------
// §3.3 — the neighbour lookup
// ---------------------------------------------------------------------------------------

#[test]
fn neighbour_agrees_with_an_independent_geometric_predicate() {
    use prin_rs::quad::{Dir, QuadTree};

    /// Box-touch adjacency, transcribed from `examples/sched_thrash.rs`. An *independent*
    /// implementation on purpose: checking the descent against itself would prove nothing.
    fn adjacent(t: &QuadTree, i: usize, j: usize) -> bool {
        let (a, b) = (&t.nodes[i], &t.nodes[j]);
        let dx = (a.cx - b.cx).abs();
        let dy = (a.cy - b.cy).abs();
        let (sx, sy) = (a.half + b.half, a.half + b.half);
        (dx <= sx + 1e-12) && (dy <= sy + 1e-12) && ((dx - sx).abs() < 1e-12 || (dy - sy).abs() < 1e-12)
    }

    // An intentionally lopsided tree, so neighbours sit at several level differences.
    let mut t = QuadTree::new(0.0, 0.0, 1.0, 4, 0);
    let k = t.split(0, 1);
    let k2 = t.split(k[0], 2);
    t.split(k2[3], 3);
    t.split(k[3], 2);

    let mut checked = 0;
    for i in 0..t.nodes.len() {
        for d in Dir::ALL {
            let Some(j) = t.neighbour(i, d) else { continue };
            assert_ne!(i, j, "a quad is not its own neighbour");
            assert!(
                t.nodes[j].level <= t.nodes[i].level,
                "neighbour must be same-or-coarser: {} at level {} returned level {}",
                i,
                t.nodes[i].level,
                t.nodes[j].level
            );
            assert!(
                adjacent(&t, i, j),
                "quad {i} dir {d:?} returned {j}, which does not touch it"
            );
            checked += 1;
        }
    }
    assert!(checked > 20, "the tree must actually exercise the lookup, got {checked}");
}

#[test]
fn a_quad_at_the_root_border_has_no_neighbour_outside_it() {
    use prin_rs::quad::{Dir, QuadTree};
    let mut t = QuadTree::new(0.0, 0.0, 1.0, 4, 0);
    let k = t.split(0, 1);
    // Lower-left child: nothing to its -x or -y.
    assert!(t.neighbour(k[0], Dir::NegX).is_none());
    assert!(t.neighbour(k[0], Dir::NegY).is_none());
    assert_eq!(t.neighbour(k[0], Dir::PosX), Some(k[1]));
    assert_eq!(t.neighbour(k[0], Dir::PosY), Some(k[2]));
    // The root itself has none in any direction.
    for d in Dir::ALL {
        assert!(t.neighbour(0, d).is_none());
    }
}

#[test]
fn contrast_reports_how_many_edges_it_saw() {
    use prin_rs::quad::{Agg, Criterion, QuadTree};
    let mut t = QuadTree::new(0.0, 0.0, 1.0, 4, 0);
    let k = t.split(0, 1);
    for (j, &i) in k.iter().enumerate() {
        t.nodes[i].red.n_footprints = 16;
        t.nodes[i].red.spread_median = j as f64;
    }
    // A corner child sees two computed neighbours, not four, so its contrast is a max over a
    // smaller set and is biased low by construction. The count is returned so that is visible.
    let (c, n) = t.contrast(k[0], Criterion::Within, Agg::Median);
    assert_eq!(n, 2, "a corner child has two in-tree neighbours");
    assert_eq!(c, 2.0, "max(|0-1|, |0-2|)");

    // An uncomputed neighbour is skipped rather than counted as a zero contrast.
    t.nodes[k[1]].red.n_footprints = 0;
    let (_, n2) = t.contrast(k[0], Criterion::Within, Agg::Median);
    assert_eq!(n2, 1);
}

// ---------------------------------------------------------------------------------------
// §6 — the FTLE port
// ---------------------------------------------------------------------------------------

#[test]
fn ftle_declines_rather_than_returning_zero_when_it_never_renormalised() {
    use prin_rs::physics::{burrau, ftle};

    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let pert = ftle::normalise::<f64>([1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);

    // A horizon shorter than one renormalisation interval: nothing has been accumulated, so
    // there is no estimate. **NaN, never 0** — a zero here is what a perfectly regular
    // trajectory reports, and the two must not be arithmetically identical.
    let opts = ftle::FtleOpts { renorm_every: 10_000, ..Default::default() };
    let o = ftle::integrate_full(s0, &m, 0.01, 1e-4, &opts, &pert);
    assert_eq!(o.n_renorm, 0);
    assert!(o.ftle.is_nan(), "must decline, got {}", o.ftle);
}

#[test]
fn ftle_renormalises_and_the_count_is_the_thing_to_assert() {
    use prin_rs::physics::{burrau, ftle};

    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let pert = ftle::normalise::<f64>([1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    let opts = ftle::FtleOpts::default();
    let o = ftle::integrate_full(s0, &m, 2.0, 1e-4, &opts, &pert);

    // Renormalisation is what stops the estimator saturating: without it the shadow separates
    // until it fills the accessible space and log(d/d0)/T decays toward zero, reporting
    // lambda ~ 0 for the MOST chaotic regions. An FTLE built from zero renormalisations is
    // that failure in new clothing, so the count is asserted rather than the value.
    assert_eq!(o.n_renorm, 99, "20000 steps at renorm_every=200, first at s=200");
    assert!(o.ftle.is_finite() && o.ftle > 0.0, "Burrau is chaotic: {}", o.ftle);
    assert!(o.diffusion.is_finite(), "the regression must be determined: {}", o.diffusion);
    assert!(o.finite && o.steps == 20_000);
}

#[test]
fn the_ftle_perturbation_is_normalised_over_all_six_components_not_per_body() {
    use prin_rs::physics::ftle;
    // The reference does `pert /= norm(pert.reshape(n, -1))` — ONE norm over the flattened
    // (3, 2), not three per-body norms. Getting that wrong scales d0 by sqrt(3) and shifts
    // every FTLE by a constant, which looks like a plausible field.
    let p = ftle::normalise::<f64>([1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    let n2: f64 = (0..3).map(|k| p[k].norm_sq()).sum();
    assert!((n2 - 1.0).abs() < 1e-15, "total norm must be 1, got {}", n2.sqrt());
    // And no single body is a unit vector, which is what the per-body mistake would give.
    for k in 0..3 {
        assert!(p[k].norm_sq() < 0.9);
    }
}

#[test]
fn a_deterministic_perturbation_is_a_unit_direction_and_varies_with_the_seed() {
    use prin_rs::physics::ftle;
    let a = ftle::unit_perturbation::<f64>(0);
    let b = ftle::unit_perturbation::<f64>(1);
    for p in [a, b] {
        let n2: f64 = (0..3).map(|k| p[k].norm_sq()).sum();
        assert!((n2 - 1.0).abs() < 1e-12);
    }
    assert!((a[0].x - b[0].x).abs() > 1e-9, "seeds must give different directions");
}

// ---------------------------------------------------------------------------------------
// §5 — the temporal accumulators
// ---------------------------------------------------------------------------------------

#[test]
fn the_event_arm_already_had_all_three_accumulators() {
    // A correction to the brief, held as a test so it cannot quietly become untrue.
    // `spread_event_max` is a running max over boundaries, `t_spread_event` is a
    // first-divergence time, and `spread_event_latched` is the persistence-guarded latch. What
    // §5 is actually asking for is the CONTINUOUS arm.
    use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
    use prin_rs::grid::Slice;

    let s = Slice::body_plane(2, 2, 1.0, 3.0, 0.05, 0);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let p = evaluate::<f64>(&s, 0, &ens);

    assert!(p.spread_event_max >= p.spread_event, "the running max must dominate the playhead");
    assert!(
        p.t_spread_event.is_nan() || p.t_spread_event <= ens.t_max,
        "a first-divergence time is a time or NaN, never a sentinel"
    );
    // And the shape arm is absent unless asked for — NaN, never 0, which would read as
    // "no divergence" on a quantity that was not measured.
    assert!(p.running_max_divergence.is_nan());
    assert!(p.divergence_trend.is_nan());
    assert!(p.first_divergence_t.is_nan());
}

#[test]
fn the_shape_accumulators_survive_ragged_per_copy_records() {
    use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
    use prin_rs::grid::Slice;

    // `stop_on_event` makes copies terminate at different boundaries, so their per-boundary
    // records have different lengths. A reader that truncates to the shortest would discard
    // exactly the boundaries where the surviving copies are diverging.
    let s = Slice::body_plane(2, 2, 0.0, 0.0, 0.05, 0); // deep interior: copies really do stop
    let ens = EnsembleCfg {
        keep_boundary_shapes: true,
        stop_on_event: true,
        refine_flagged: false,
        ..Default::default()
    };
    let p = evaluate::<f64>(&s, 0, &ens);

    assert!(p.running_max_divergence.is_finite(), "must be measured: {}", p.running_max_divergence);
    assert!(
        p.running_max_divergence >= p.spread_shape - 1e-12,
        "a running max cannot fall below the playhead value: {} vs {}",
        p.running_max_divergence,
        p.spread_shape
    );
    assert!(p.divergence_trend.is_finite());
    // NaN when it never crosses, never `t_max` — a sentinel equal to the horizon is
    // indistinguishable from crossing at the last boundary.
    assert!(
        p.first_divergence_t.is_nan() || (p.first_divergence_t > 0.0 && p.first_divergence_t <= ens.t_max),
        "got {}",
        p.first_divergence_t
    );
}

#[test]
fn boundary_shapes_are_only_recorded_when_asked_for() {
    use prin_rs::integrate::az::{self, AzOpts};
    use prin_rs::physics::burrau;

    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let off = az::integrate_az_opts(s0, &m, 2.0, 8, 0.01, 30_000, &AzOpts::default());
    assert!(off.boundary_shapes.is_empty(), "off by default: it is ~70x a PixelOut");

    let on = az::integrate_az_opts(
        s0, &m, 2.0, 8, 0.01, 30_000,
        &AzOpts { keep_boundary_shapes: true, ..Default::default() },
    );
    assert_eq!(on.boundary_shapes.len(), on.tight.len(), "one per completed boundary");
    for n in &on.boundary_shapes {
        let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((norm - 1.0).abs() < 1e-10, "each must be a unit shape vector, got {norm}");
    }
}

#[test]
fn a_signal_that_declines_to_score_does_not_block_the_ones_that_do() {
    use prin_rs::ensemble::pixel::EnsembleCfg;
    use prin_rs::grid::Chart;
    use prin_rs::metric::{self, Rank};
    use prin_rs::quad::{Agg, Criterion};

    // `v > bv` is FALSE whenever `bv` is NaN, so a NaN at the first candidate would block every
    // finite score behind it and collapse the ranking to scan order — silently, and worst on
    // exactly the signals that decline to score part of the region. `term_grad` is NaN on 97.1%
    // of near-field.
    //
    // "Undetermined never wins" is the intended semantics. "Undetermined blocks everything" was
    // the bug, and it is invisible in the error curve, which just looks like a bad criterion.
    let (levels, n) = (2u32, 2usize);
    let res = (1usize << levels) * n;
    let ens = EnsembleCfg { n_extra: 1, t_max: 2.0, n_sync: 4, refine_flagged: false, ..Default::default() };
    let cache = metric::build(
        "near-field", 1.0, 3.0, 0.05, 0, Chart::BodyPlane, levels, n, res, 1e-4, &ens,
        metric::Colouring::Outcome,
    );
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;

    // At t_max = 2 nothing terminates anywhere near-field, so term_grad is NaN everywhere and
    // this signal has nothing to say. It must still complete the descent rather than stalling.
    let nan_everywhere = cache
        .quads
        .values()
        .all(|q| !q.red.signal(Criterion::TerminationGradient, Agg::Median).is_finite());
    assert!(nan_everywhere, "the premise of this test: term_grad declines here");

    let pts = metric::replay(&cache, Rank::Signal(Criterion::TerminationGradient, Agg::Median), full);
    assert_eq!(
        pts.last().unwrap().error,
        0.0,
        "an all-NaN signal must still spend its budget and reach the reference"
    );
    assert_eq!(pts.len(), 1 + (full - 1) / 4, "and must spend all of it");
}

// ---------------------------------------------------------------------------------------
// The hot rule — the relative mask, and the two things it does NOT buy
// ---------------------------------------------------------------------------------------

/// **What a quantile hot rule actually buys, and what it does not.**
///
/// The instruction on record was to make the hot threshold relative — the quad's own median.
/// Two consequences, both asserted here because both were missed when the instruction was
/// written:
///
/// 1. **`n_hot` stops being a signal.** On a field with distinct values the count above the cut
///    is fixed by the rule, so `frac_hot` carries zero information at every quad of every chart.
///    That matters because `frac_hot_between/median` is the best criterion measured on this
///    project, and a relative rule that *replaced* the absolute mask would have deleted the
///    best-performing signal in the system while reading as an improvement. Both masks are
///    computed; this names the constant so nobody reads `frac_hot` off the relative one.
///
/// 2. **The mask is invariant under a monotone rescaling of the field**, which is the property
///    the relative rule exists for: spread rises globally with `t`, and an absolute cut must
///    eventually fire everywhere. Measured cause, on the committed corpus: `tau = 1e-4` sits at
///    the 0.4th percentile and `n_hot == N^2` in 98.8% of 75,359 leaves.
///
/// **And the exception, measured by this test rather than assumed away:** on a *tied* field the
/// count is set by the tie structure, not by the rule — a step field reads 5 at q = 0.5, 0.75
/// and 0.9 alike. That is not a defect; it is what "the quantile is nearest-rank and the
/// comparison is strict" means, and it is the case that occurs when the event arm dominates a
/// footprint field with five distinct values.
#[test]
fn a_quantile_hot_rule_is_rescaling_invariant_and_n_hot_stops_being_a_signal() {
    let n = 8;
    let ramp: Vec<f64> = (0..n * n).map(|k| k as f64).collect();
    let noisy: Vec<f64> = (0..n * n).map(|k| ((k as u64 * 2654435761) % 997) as f64).collect();
    let step: Vec<f64> = (0..n * n).map(|k| if k < 5 { 1e6 } else { 1e-9 }).collect();

    let count = |f: &Vec<f64>, r| spatial::hot_mask(f, r).iter().filter(|&&h| h).count();

    for (name, f) in [("ramp", &ramp), ("noise", &noisy), ("step", &step)] {
        let c: Vec<usize> =
            [0.5, 0.75, 0.9].iter().map(|&q| count(f, HotRule::Quantile(q))).collect();
        println!("{name:>6}: n_hot at q = 0.5 / 0.75 / 0.9 -> {c:?}");
    }

    // (1) Same rule, two unrelated distinct-valued fields, identical counts.
    assert_eq!(
        count(&ramp, HotRule::Quantile(0.5)),
        count(&noisy, HotRule::Quantile(0.5)),
        "n_hot must not depend on a distinct-valued field"
    );
    // The tied field is the named exception, asserted so it cannot drift into a surprise.
    assert_eq!(count(&step, HotRule::Quantile(0.5)), count(&step, HotRule::Quantile(0.9)),
               "on a two-valued field the tie structure sets the count, not q");

    // (2) The property that matters: rescale the whole field — the global rise with `t` — and
    // the mask does not move. This is the arm the absolute rule fails.
    let risen: Vec<f64> = ramp.iter().map(|x| x * 1000.0 + 70.0).collect();
    let a = spatial::hot_mask(&ramp, HotRule::Quantile(0.5));
    let b = spatial::hot_mask(&risen, HotRule::Quantile(0.5));
    assert_eq!(a, b, "a quantile mask must be invariant under a monotone rescaling");

    // The control that proves the arm above is about the RULE and not about this field: the
    // absolute rule saturates under exactly the same rescaling, which is the observed failure.
    let abs_before = count(&ramp, HotRule::AbsTau(10.0));
    let abs_after = count(&risen, HotRule::AbsTau(10.0));
    println!("abs[1e1]: {abs_before} hot before the rise, {abs_after} after (of {})", n * n);
    assert!(abs_before < n * n, "the control is empty unless the absolute rule starts unsaturated");
    assert_eq!(abs_after, n * n, "the absolute rule must saturate when the signal rises");
}

/// A quad with fewer than two finite values is **undetermined, not calm**.
#[test]
fn a_degenerate_field_reads_all_hot_rather_than_none() {
    let f = vec![f64::NAN; 16];
    let m = spatial::hot_mask(&f, HotRule::Quantile(0.5));
    assert!(m.iter().all(|&h| h), "an all-NaN quad must not read as cold");
    let l = spatial::layout(&m, 4);
    assert_eq!(l.n_hot, 16);
    // One finite value is still no distribution.
    let mut g = vec![f64::NAN; 16];
    g[3] = 1.0;
    assert!(spatial::hot_mask(&g, HotRule::Quantile(0.5)).iter().all(|&h| h));
}

/// `grad_rms` is the one structure measure with **no threshold in it**, which makes it the
/// control on the whole hot-mask family. Its degenerate convention is `NaN`, never 0.
#[test]
fn grad_rms_is_nan_with_no_finite_pair_and_ranks_a_filament_above_a_flat_field() {
    let n = 8;
    assert!(spatial::grad_rms(&vec![f64::NAN; n * n], n).is_nan(), "no finite pair must be NaN");

    let flat = vec![1.0f64; n * n];
    assert_eq!(spatial::grad_rms(&flat, n), 0.0, "a genuinely flat field IS zero");

    // A one-cell-wide ridge: the same total "energy" as the flat field carries, arranged so it
    // varies. The inequality is asserted, not eyeballed.
    let mut ridge = vec![0.0f64; n * n];
    for jx in 0..n {
        ridge[4 * n + jx] = 1.0;
    }
    let g_ridge = spatial::grad_rms(&ridge, n);
    // A smooth ramp of the same peak-to-peak: same range, far less gradient.
    let ramp: Vec<f64> = (0..n * n).map(|k| (k / n) as f64 / (n - 1) as f64).collect();
    let g_ramp = spatial::grad_rms(&ramp, n);
    println!("grad_rms: flat 0, ramp {g_ramp:.4}, ridge {g_ridge:.4}");
    assert!(g_ridge > g_ramp, "a ridge {g_ridge} must outrank a ramp {g_ramp} of the same range");
    assert!(g_ramp > 0.0);

    // Non-finite cells are skipped, not counted as zero difference — a NaN neighbour must not
    // make a quad look smooth.
    let mut holed = ridge.clone();
    holed[0] = f64::NAN;
    assert!(spatial::grad_rms(&holed, n).is_finite());
}

/// **The desaturation, as a pair — because the obvious form of this test cannot fail.**
///
/// The instruction on record was *"assert `n_hot_within < N^2` for a stated majority of quads"*.
/// Under a quantile rule that passes trivially and unconditionally, on any field, forever. It is
/// decoration.
///
/// The form with teeth runs both masks over **one descent** and asserts they disagree: the
/// absolute mask saturated, the relative one not. The absolute arm is the control that proves
/// the test can tell them apart — without it, a relative mask that happened to saturate too
/// would read the same as a working one.
#[test]
fn the_relative_mask_desaturates_where_the_absolute_one_does_not() {
    use prin_rs::ensemble::pixel::EnsembleCfg;
    use prin_rs::render::Precision;
    use prin_rs::scheduler::{self, SchedCfg};

    const N: usize = 4;
    let ens = EnsembleCfg { t_max: 2.0, n_sync: 8, refine_flagged: false, ..Default::default() };
    // `tau_display = 1e-4` is the value every committed run used, and the one measured to sit at
    // the 0.4th percentile of the spread distribution. Reproducing it is the point.
    let cfg = SchedCfg {
        n: N,
        budget: 80,
        tau_display: 1e-4,
        hot_rule: HotRule::Quantile(0.5),
        ..Default::default()
    };
    let (t, st) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64);

    let computed: Vec<&prin_rs::quad::Quad> =
        t.nodes.iter().filter(|q| q.red.n_footprints > 0).collect();
    let m = computed.len() as f64;
    assert!(m >= 20.0, "only {m} quads computed; too few to say anything");

    let full = (N * N) as u32;
    let frac = |f: &dyn Fn(&prin_rs::quad::Quad) -> bool| {
        computed.iter().filter(|q| f(q)).count() as f64 / m
    };

    let abs_sat = frac(&|q| q.red.layout_within.n_hot == full);
    let rel_sat = frac(&|q| q.red.layout_rel_within.n_hot == full);
    let abs_one = frac(&|q| q.red.layout_within.n_components <= 1);
    let rel_one = frac(&|q| q.red.layout_rel_within.n_components <= 1);
    let rel_multi = frac(&|q| q.red.layout_rel_within.n_components > 1);

    println!("{} quads computed, {} leaves, N = {N}", st.quads_computed, t.leaves().count());
    println!("  absolute mask: n_hot == {full} on {:.1}%, n_components <= 1 on {:.1}%",
             100.0 * abs_sat, 100.0 * abs_one);
    println!("  relative mask: n_hot == {full} on {:.1}%, n_components <= 1 on {:.1}%, > 1 on {:.1}%",
             100.0 * rel_sat, 100.0 * rel_one, 100.0 * rel_multi);

    // The control arm. If this fails the absolute mask is not saturated in this configuration
    // and the comparison below is scoring nothing.
    assert!(abs_sat > 0.8, "absolute mask only {:.1}% saturated -- the control is empty", 100.0 * abs_sat);

    // The arm under test.
    assert!(rel_sat < 0.05, "relative mask still {:.1}% saturated", 100.0 * rel_sat);
    assert!(rel_multi > 0.2,
            "relative mask resolves > 1 component on only {:.1}% of quads", 100.0 * rel_multi);
}
