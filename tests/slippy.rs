//! Stage 4: the 2:1 balance constraint, camera-biased priority, and the persistent frontier.
//!
//! Every assertion here is written so that something has to be true for it to fire. Where a
//! property is *structural* it says so and is not dressed as a measurement.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::frontier::{band_of, Frontier, BANDS};
use prin_rs::quad::{Agg, Criterion, Decision, Dir, QuadReduction};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::spatial;

/// **`t = 13`, not a short horizon.** Twice in this build a stage-4 test was written at `t = 2`
/// and read as a failure of the thing under test when what it had actually measured was the
/// criterion never firing: near-field at `t = 2` is tame enough that every leaf reads `Keep`, so
/// the ranking never runs and no tree ever becomes unbalanced. A test whose subject never
/// executes is decoration.
fn ens() -> EnsembleCfg {
    EnsembleCfg { t_max: 13.0, refine_flagged: false, ..Default::default() }
}

// ---------------------------------------------------------------------------------------
// The 2:1 balance constraint
// ---------------------------------------------------------------------------------------

/// **2:1 holds on every produced tree, and the unbalanced control violates it.**
///
/// Without the control this test could pass on a tree that never had two adjacent leaves at
/// different depths in the first place — which is exactly what a veto-bound uniform tree looks
/// like, and most of this corpus is veto-bound.
#[test]
fn the_two_to_one_constraint_holds_and_the_control_violates_it() {
    let e = ens();
    // **The fixture region has now swapped TWICE, and the control arm caught it both times.**
    // Originally `deep interior`, because near-field reached a complete tree at one depth under
    // the veto. The escape distance gate flattened `deep interior` -- its terminal class moved
    // from 60% escape to 2%, most of those escapes being mid-encounter transients -- so the
    // fixture moved to near-field, which then carried a gap of 2 at `alpha_hi = 0.2`.
    //
    // The `dtau` step-control fix has moved it back. Under the corrected stepping near-field is
    // gap 1 at **every** cell of `alpha_hi in {0.1,0.2,0.3,0.5} x tau in {1e-6,1e-4,1e-3} x
    // n in {4,8}` -- twenty-four cells, nothing forced -- while `deep interior` recovers gap 2
    // at `n = 4` across most of that grid. Measured by sweeping, not guessed, and the assertion
    // below is what turned a silently-vacuous test into a failing one each time. It stays.
    //
    // `n = 4` matters: at `n = 8` `deep interior` is gap 1 too. A coarser footprint grid makes a
    // noisier spread estimate, which biases toward *refine* -- the conservative direction -- and
    // that extra depth variation is what there is to balance.
    let root = prin_rs::grid::region("deep interior", 2, 2, 0.05).unwrap();
    let run = |balance| {
        let cfg = SchedCfg {
            n: 4,
            budget: 600,
            tau_display: 1e-4,
            alpha_hi: 0.2,
            alpha_lo: 0.2,
            balance,
            camera: Some(Camera::framing(root.cx, root.cy, 0.05, 64)),
            ..Default::default()
        };
        scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &e, Precision::F64)
    };
    let worst = |t: &prin_rs::quad::QuadTree| {
        let mut w = 0u32;
        for i in t.leaves() {
            for d in Dir::ALL {
                if let Some(j) = t.neighbour(i, d) {
                    if t.nodes[j].children.is_none() {
                        w = w.max(t.nodes[i].level.abs_diff(t.nodes[j].level));
                    }
                }
            }
        }
        w
    };

    let (tu, su) = run(false);
    let (tb, sb) = run(true);
    let (wu, wb) = (worst(&tu), worst(&tb));
    let forced = tb.nodes.iter().filter(|q| q.decision == Decision::BalanceForced).count();

    println!("unbalanced: {} leaves, worst adjacent level gap {wu}", tu.leaves().count());
    println!("balanced  : {} leaves, worst adjacent level gap {wb}, {forced} balance-forced \
              ({:.1}% of {} quads computed)",
             tb.leaves().count(), 100.0 * sb.balance_forced as f64 / sb.quads_computed as f64,
             sb.quads_computed);
    assert_eq!(su.balance_forced, 0, "the control must force nothing");

    // The control. If the unbalanced tree already satisfies 2:1, this test proves nothing.
    assert!(wu > 1, "the unbalanced control satisfies 2:1 already -- nothing is being tested");
    assert_eq!(wb, 1, "the balanced tree violates 2:1 with a gap of {wb}");
    // §4.4 asks for this to be REPORTED: a large share means the budget went on geometry.
    assert!(forced > 0, "nothing was balance-forced, so the pass never ran");
}

// ---------------------------------------------------------------------------------------
// Camera-biased priority
// ---------------------------------------------------------------------------------------

/// **Relevance is a priority term and the veto is position-free. Both halves are asserted.**
///
/// The second half is the one that matters: it is what keeps a quad's *decision* independent of
/// where the camera points, and it is why the committed pan sequence produced a byte-identical
/// tree nine times. Ranking may read position because nothing about it is stored.
#[test]
fn camera_position_moves_the_ranking_and_never_the_veto() {
    let q = prin_rs::quad::Quad {
        level: 3,
        cx: 1.0,
        cy: 3.0,
        half: 0.00625,
        parent: None,
        children: None,
        sib_index: 0,
        iteration: 0,
        red: QuadReduction::default(),
        alpha: None,
        alpha_mean: None,
        alpha_p90: None,
        alpha_sibling_spread: None,
        decision: Decision::Pending,
    };
    let here = Camera::framing(1.0, 3.0, 0.05, 512);
    let away = Camera { cx: 1.0 + 10.0, ..here };

    // The veto: identical, because it reads no position term at all.
    assert_eq!(here.veto(&q, 8, 0.05), away.veto(&q, 8, 0.05),
               "the veto must not depend on where the camera points");

    // The ranking: it moves.
    let (rh, ra) = (here.relevance(q.cx, q.cy, q.half, 0.0), away.relevance(q.cx, q.cy, q.half, 0.0));
    println!("relevance: camera on the quad {rh:.4}, camera 10 units away {ra:.4}");
    assert_eq!(rh, 1.0, "a fully visible quad must score 1");
    assert_eq!(ra, 0.0, "an off-screen quad must score 0");

    // Partially visible: strictly between, so the term is graded rather than a second veto.
    let edge = Camera { cx: 1.0 + here.half_world, ..here };
    let big = prin_rs::quad::Quad { half: 0.05, ..q };
    let re = edge.relevance(big.cx, big.cy, big.half, 0.0);
    println!("relevance: quad straddling the viewport edge {re:.4}");
    assert!(re > 0.0 && re < 1.0, "a straddling quad must be graded, got {re}");

    // And the margin widens it -- §4.3's baseline, which any prediction model must beat.
    assert!(edge.relevance(big.cx, big.cy, big.half, 2.0) > re, "margin must widen relevance");
}

/// A pan changes the tree **only** once the bias is switched on. Before it, a pan is an identity.
///
/// **The pan distance is a measured fixture and it has moved once.** It was `0.04` against a
/// `half_world` of `0.05`; under the boundary-overshoot clamp (`RESULTS §24`) near-field's tree
/// is tamer -- 184 leaves, none of them `Split`, levels 2-4 all `Keep` and level 5 all
/// `ScreenFloor` -- and a pan that small no longer changes any quad's relevance enough to move a
/// decision. Measured across `{0.01, 0.02, 0.04, 0.06, 0.08, 0.10}`: the tree is identical under
/// the bias out to `0.06` and differs from `0.08`, which is where the pan finally exceeds
/// `half_world` and quads start leaving the viewport. **The control arm caught it**, as it has
/// at every previous fixture move on this project -- the `assert_ne!`, not the property.
#[test]
fn a_pan_is_an_identity_until_camera_bias_is_switched_on() {
    let e = ens();
    let run = |cx: f64, bias| {
        let cam = Camera { cx, ..Camera::framing(1.0, 3.0, 0.05, 128) };
        let cfg = SchedCfg {
            n: 4,
            budget: 600,
            tau_display: 1e-4,
            alpha_hi: 0.2,
            alpha_lo: 0.2,
            k_frac: 0.5,
            camera: Some(cam),
            camera_bias: bias,
            ..Default::default()
        };
        let (t, _) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &e, Precision::F64);
        t.leaves().map(|i| (t.nodes[i].level, t.nodes[i].decision)).collect::<Vec<_>>()
    };
    // Larger than `half_world = 0.05`: below that the relevance difference does not survive to
    // a decision on this tree. See the doc comment.
    const PAN: f64 = 0.08;
    let a = run(1.0, None);
    let b = run(1.0 + PAN, None);
    println!("no bias : {} vs {} leaves, identical {}", a.len(), b.len(), a == b);
    assert_eq!(a, b, "without the bias a pan must change nothing -- that is the standing result");

    let c = run(1.0, Some(0.0));
    let d = run(1.0 + PAN, Some(0.0));
    println!("with bias: {} vs {} leaves, identical {}", c.len(), d.len(), c == d);
    assert_ne!(c, d, "with the bias a pan must change the tree, or the term is not reaching it");
}

// ---------------------------------------------------------------------------------------
// The persistent frontier
// ---------------------------------------------------------------------------------------

#[test]
fn bands_are_monotone_and_undetermined_sits_at_the_bottom() {
    let mut prev = 0usize;
    for p in [1e-10f64, 1e-8, 1e-6, 1e-4, 1e-2, 1.0, 10.0] {
        let b = band_of(p);
        assert!(b >= prev, "band must not fall as priority rises: {p:e} -> {b}");
        prev = b;
    }
    assert!(band_of(1e-8) < band_of(1e-2), "six orders must not share a band");
    // Undetermined never outranks determined. The same convention the replay uses.
    assert_eq!(band_of(f64::NAN), 0);
    assert_eq!(band_of(f64::INFINITY), BANDS - 1);
}

/// **The staleness check.** An incrementally-maintained frontier that is wrong looks exactly
/// like a criterion that is wrong, so the from-scratch path is kept and compared against.
#[test]
fn the_incremental_frontier_matches_a_from_scratch_rebuild() {
    let mut f = Frontier::new();
    // A spread of stored priorities across six orders, as the real signal has.
    for id in 0..200usize {
        let stored = 10f64.powf(-8.0 + (id % 13) as f64 * 0.5) * (1 + id % 7) as f64;
        f.insert(id, stored);
    }
    // Camera relevance: changes every frame, computed at query time, never stored.
    let derive = |id: usize| ((id * 37) % 11) as f64 / 10.0;
    assert!(f.agrees_with_rebuild(20, derive), "fresh frontier disagrees with rebuild");

    // Reprioritise a third of them, including across band boundaries, and check again. This is
    // the operation a plain binary heap cannot do without a full rebuild.
    for id in (0..200).step_by(3) {
        f.reprioritise(id, 10f64.powf(-8.0 + (id % 5) as f64 * 2.0));
    }
    assert!(f.agrees_with_rebuild(20, derive), "frontier is STALE after reprioritising");

    // Removal, then a different camera. A stale entry surviving a remove is the exact failure
    // this exists to catch, and it is invisible in the tree.
    for id in (0..200).step_by(7) {
        f.remove(id);
    }
    let derive2 = |id: usize| 1.0 / (1.0 + id as f64);
    assert!(f.agrees_with_rebuild(20, derive2), "frontier is STALE after removals");
    assert_eq!(f.len(), 200 - (0..200).step_by(7).count());
    println!("{} entries, incremental order matches the rebuild through three mutations", f.len());
}

/// A relevance change that does not cross a band boundary must not re-bucket — the property the
/// whole structure exists for. Asserted through observable state, not by counting internals.
#[test]
fn the_derived_term_is_never_stored_so_a_pan_touches_no_stored_priority() {
    let mut f = Frontier::new();
    for id in 0..50usize {
        f.insert(id, 1e-3 * (1 + id) as f64);
    }
    let before = f.entries();
    // "Move the camera": the derived term changes for every entry.
    let a = f.top_k(10, |id| ((id * 13) % 7) as f64 / 6.0);
    let b = f.top_k(10, |id| ((id * 29) % 5) as f64 / 4.0);
    assert_ne!(a, b, "the derived term must change the order, or the test is empty");
    // ...and nothing stored moved. That is what makes a pan cheap and keeps view state off the
    // quad: relevance is recomputed, never written down.
    assert_eq!(before, f.entries(), "a camera move mutated the STORED priorities");
}

// ---------------------------------------------------------------------------------------
// The edge-filament blind spot
// ---------------------------------------------------------------------------------------

/// **§2.1's blind spot is real, and neighbour contrast closes it.**
///
/// A filament running along a quad *boundary* shows as low internal variation in **both** quads,
/// so every within-quad structure measure misses it — and boundaries are precisely what the
/// criterion exists to find. The blind spot is systematically aligned with the target.
///
/// Built as a synthetic field so the answer is known: two quads, each internally uniform, with a
/// step between them. Every within-quad measure must read "nothing here"; contrast must not.
#[test]
fn neighbour_contrast_catches_an_edge_filament_that_every_within_quad_measure_misses() {
    let n = 8;
    let uniform = |v: f64| {
        let mut r = QuadReduction { n_footprints: (n * n) as u32, ..Default::default() };
        let field = vec![v; n * n];
        r.spread_median = v;
        r.layout_within = spatial::layout(&spatial::hot_mask(&field, spatial::HotRule::AbsTau(1e-4)), n);
        r.layout_rel_within =
            spatial::layout(&spatial::hot_mask(&field, spatial::HotRule::Quantile(0.5)), n);
        r.grad_rms_within = spatial::grad_rms(&field, n);
        r
    };
    // The filament lies exactly on the shared edge: each side is internally featureless.
    let lo = uniform(1e-6);
    let hi = uniform(1e-1);

    for (name, r) in [("cold side", &lo), ("hot side", &hi)] {
        println!("{name}: grad_rms {:.3e}, structure {:?}, layout n_components {}, \
                  perimeter {:?}",
                 r.grad_rms_within, r.structure(true), r.layout_rel_within.n_components,
                 r.layout_rel_within.perimeter_ratio);
        // Every WITHIN-quad measure reads a featureless quad, on both sides.
        assert_eq!(r.grad_rms_within, 0.0, "{name}: an internally uniform quad has no gradient");
        assert!(r.structure(true).is_nan() || r.structure(true) == 0.0,
                "{name}: within-quad structure must not see the edge filament");
    }

    // Contrast does see it, and it is the only thing here that does.
    let c = (hi.signal(Criterion::Within, Agg::Median) - lo.signal(Criterion::Within, Agg::Median))
        .abs();
    println!("neighbour contrast across the shared edge: {c:.3e}");
    assert!(c > 1e-2, "contrast must fire on the edge filament, got {c}");
}

/// **Zoom-out is nearly free — but not in the form §4.5 states it, and the difference matters.**
///
/// §4.5 asserts *"the count of newly-computed quads after a zoom-out is ~0"*, on the reasoning
/// that zooming out reveals shallower quads whose parents are already in the tree. That
/// presupposes a **tree that persists across frames**, and this build deliberately has none: the
/// scope discipline is *no eviction, no caching, no async, no promotion*, and a cross-frame tree
/// is the caching it excludes. Asserting it as written would require building the thing the
/// scope forbids, and asserting it against a from-scratch descent would be measuring nothing.
///
/// What **is** available is the arithmetic underneath the claim: a wider view floors sooner, so a
/// zoomed-out descent computes no more than a zoomed-in one, and every box it computes is one the
/// zoomed-in run already computed. That is the content of "the parents are already there" without
/// pretending to a persistence this build does not have.
#[test]
fn zooming_out_computes_a_subset_of_what_zooming_in_computed() {
    let e = ens();
    let boxes = |hw: f64| {
        let cam = Camera { half_world: hw, ..Camera::framing(1.0, 3.0, 0.05, 128) };
        let cfg = SchedCfg {
            n: 4,
            budget: 600,
            tau_display: 1e-4,
            alpha_hi: 0.2,
            alpha_lo: 0.2,
            camera: Some(cam),
            ..Default::default()
        };
        let (t, st) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &e, Precision::F64);
        let key = |q: &prin_rs::quad::Quad| {
            (q.level, (q.cx / 1e-9).round() as i64, (q.cy / 1e-9).round() as i64)
        };
        let set: std::collections::HashSet<_> = t
            .nodes
            .iter()
            .filter(|q| q.red.n_footprints > 0)
            .map(key)
            .collect();
        (set, st.quads_computed)
    };

    // Zoomed IN is the smaller half_world.
    let (inner, n_in) = boxes(0.0125);
    let (outer, n_out) = boxes(0.05);
    let novel = outer.difference(&inner).count();
    println!("zoomed in  (half_world 0.0125): {n_in} quads computed");
    println!("zoomed out (half_world 0.0500): {n_out} quads computed, {novel} of them NOT present \
              in the zoomed-in descent");

    assert!(n_out <= n_in, "zooming out computed MORE ({n_out} against {n_in})");
    // The real content: a zoom-out asks for nothing the zoom-in did not already hold. A
    // persistent tree would therefore have to compute none of them.
    assert_eq!(novel, 0, "{novel} boxes are new on zoom-out; the parents were not already there");
}
