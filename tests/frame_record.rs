//! Tests for the per-frame record and the fill fractions.
//!
//! The claims: a stage this build does not have reads `NaN` rather than `0.0`; the headline
//! statistic is measured during motion and says so when nothing moved; percentiles are not means;
//! and the fill fractions sum to one and separate the leaf's own paint from an ancestor's.

use prin_rs::camera::Camera;
use prin_rs::output::frame::{
    frac_over_floor_moving, percentiles, record, FrameRecord, StageMs, BUDGET_24, BUDGET_60,
    FIELDS, N_FIELDS,
};

/// **A stage this build does not have is `NaN`, never `0.0`.**
///
/// There is no GPU and no window here, so `upload` and `present` are not measured — and a zero
/// there reads as *instant* where the truth is *absent*. The same conflation as an empty mask
/// reading as "no structure found", or `Absent` pooled with `Evicted`.
#[test]
fn an_unmeasured_stage_is_nan_and_not_zero() {
    let s = StageMs::default();
    assert!(s.upload.is_nan(), "upload must be NaN: there is no GPU in this build");
    assert!(s.present.is_nan(), "present must be NaN: there is no window in this build");
    assert_eq!(s.integrate, 0.0, "a stage that IS measured starts at zero");

    // And it survives the record, rather than being flattened on the way out.
    let r = FrameRecord { stage_ms: s, ..Default::default() };
    let row = record(&r);
    let up = FIELDS.iter().position(|&f| f == "upload_ms").unwrap();
    let pr = FIELDS.iter().position(|&f| f == "present_ms").unwrap();
    assert!(row[up].is_nan() && row[pr].is_nan(), "the record flattened an absent stage to a number");

    // A frontier audit that did not run must not report a pass.
    assert!(FrameRecord::default().frontier_agrees.is_nan() || FrameRecord::default().frontier_agrees == 0.0);
}

/// **The field list and the record are tied at compile time**, so adding one alone breaks the
/// build rather than shifting every column at read time.
#[test]
fn the_field_names_and_the_record_agree() {
    assert_eq!(FIELDS.len(), N_FIELDS);
    assert_eq!(record(&FrameRecord::default()).len(), N_FIELDS);
    let mut names: Vec<&str> = FIELDS.to_vec();
    names.sort_unstable();
    let n = names.len();
    names.dedup();
    assert_eq!(names.len(), n, "two fields share a name");
}

/// **The headline is measured during motion, and says nothing when nothing moved.**
///
/// Static frames may take longer without anyone minding; a dropped frame mid-pan is immediately
/// visible. A harness that reported `0.0` for a log with no motion would be claiming a pass it
/// never earned, so the empty case is `NaN`.
#[test]
fn the_headline_is_motion_only_and_nan_when_nothing_moved() {
    let still = |ms: f64| FrameRecord { frame_ms: ms, ..Default::default() };
    let moving = |ms: f64| FrameRecord { frame_ms: ms, camera_delta: 0.01, ..Default::default() };

    // All static, several of them slow: the question was not asked.
    let log = vec![still(100.0), still(200.0), still(1.0)];
    assert!(frac_over_floor_moving(&log).is_nan(), "a log with no motion must not score");

    // Static frames must not dilute the moving ones either.
    let log = vec![still(999.0), moving(10.0), moving(50.0), still(999.0)];
    let f = frac_over_floor_moving(&log);
    assert!((f - 0.5).abs() < 1e-12, "expected 1 of 2 moving frames over the floor, got {f}");

    assert!(moving(50.0).over_floor() && !moving(10.0).over_floor());
    assert!(moving(20.0).over_goal() && !moving(10.0).over_goal());
    assert!(BUDGET_60 < BUDGET_24, "the goal must be tighter than the floor");
}

/// **Percentiles, never means** — a mean hides the stutter that makes a thing feel broken.
#[test]
fn the_percentiles_are_not_a_mean() {
    // Ninety-nine fast frames and one catastrophic one: the mean barely moves, p99 and max do.
    let mut v: Vec<f64> = vec![10.0; 99];
    v.push(1000.0);
    let p = percentiles(&v);
    let mean = v.iter().sum::<f64>() / v.len() as f64;

    assert_eq!(p.p50, 10.0);
    assert_eq!(p.max, 1000.0);
    assert!(mean < BUDGET_24, "the mean reads {mean} -- inside budget, which is the whole problem");
    assert!(p.max > BUDGET_24, "the max is what says a frame was dropped");
    assert_eq!(p.n, 100);

    assert_eq!(percentiles(&[]).n, 0);
    assert!(percentiles(&[]).p50.is_nan(), "an empty log has no p50, and must not report 0");
}

/// **The fill fractions sum to one and separate a leaf's own paint from an ancestor's.**
///
/// §4.5 takes the coarse-ancestor fill as the option that never lies, and big texels during motion
/// as deliberate — but the share has to be visible, or a scheduler that never converges looks the
/// same as one that has. The control is the complete tree, where nothing is served by an ancestor.
#[test]
fn the_fill_fractions_separate_own_paint_from_ancestor_paint() {
    use prin_rs::output::adaptive;
    use prin_rs::quad::QuadTree;

    let mut tree = QuadTree::with_chart(0.0, 0.0, 1.0, 4, 0, prin_rs::grid::Chart::BodyPlane);
    tree.split(0, 0);
    let kids = tree.nodes[0].children.unwrap();
    for k in kids {
        tree.split(k, 1);
    }
    let cam = Camera::framing(0.0, 0.0, 1.0, 64);
    let leaves: Vec<usize> = tree.leaves().collect();

    // Everything resident: no ancestor is ever needed.
    let all = adaptive::fill_fractions(&tree, &cam, 64, &leaves, &|_| true);
    assert!((all.own + all.ancestor + all.background - 1.0).abs() < 1e-9, "shares must sum to 1");
    assert!(all.own > 0.99, "a complete resident tree is painted by its own leaves: {all:?}");
    assert_eq!(all.ancestor, 0.0, "nothing should fall back when everything is resident");

    // One quadrant's leaves released: their region is served by the ancestor, and the SHARE says so.
    let gone: std::collections::HashSet<usize> =
        tree.nodes[kids[0]].children.unwrap().iter().cloned().collect();
    let some = adaptive::fill_fractions(&tree, &cam, 64, &leaves, &|i| !gone.contains(&i));
    assert!((some.own + some.ancestor + some.background - 1.0).abs() < 1e-9);
    assert!(some.ancestor > 0.2, "a released quadrant must show as ancestor fill: {some:?}");
    assert!(some.own < all.own, "own paint must fall when leaves are released");
    assert_eq!(some.background, 0.0, "the root is resident, so nothing is bare background");

    // **The whole chain gone: now there IS bare background** — a third statement again, and the
    // one the renderer must never show without saying so. Releasing the leaves and the root is not
    // enough: the level-1 quadrant between them is still resident and paints the region as an
    // ancestor, which is exactly the fallback working.
    let chain: std::collections::HashSet<usize> =
        gone.iter().cloned().chain([kids[0], 0]).collect();
    let bare = adaptive::fill_fractions(&tree, &cam, 64, &leaves, &|i| !chain.contains(&i));
    assert!(bare.background > 0.2, "with no ancestor left in the chain it is background: {bare:?}");
    assert!(bare.ancestor < some.ancestor, "removing the fallback must reduce ancestor paint");
}
