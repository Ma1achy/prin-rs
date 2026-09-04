//! Tests for the frame loop.
//!
//! The claims: a frame spends at most its quota; the tree **persists** across frames rather than
//! being rebuilt; the quota is a count and not a clock; a pan touches nothing stored and a zoom
//! touches every stored term; and the degenerate cell is refused rather than trusted.

use prin_rs::camera::Camera;
use prin_rs::grid::Slice;
use prin_rs::scheduler::{Mode, Policy, SchedCfg};
use prin_rs::session::{assert_session_engages, FrameQuota, Regrow, Session, SessionCfg, QuotaHit};

fn cfg(quads: usize, rounds: u32) -> SessionCfg {
    SessionCfg {
        sched: SchedCfg {
            n: 8,
            budget: 20000,
            tau_display: 1e-2,
            policy: Policy::Tolerance,
            mode: Mode::Balanced,
            max_level: Some(5),
            camera: None,
            ..Default::default()
        },
        quota: FrameQuota { quads, substeps: None, rounds },
        ..Default::default()
    }
}

fn field() -> impl Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync {
    // The step, offset off the midline: at `x0 = 0` on a root spanning [-1,1] the discontinuity
    // falls exactly on the quad boundary at every level, no quad straddles it, and the field is
    // featureless to the criterion. That cost the golden matrix four rounds to notice.
    prin_rs::testing::step(0.137, 13.0)
}

/// **A frame spends at most its quota, and says which one bound.**
#[test]
fn a_frame_is_bounded_by_its_quota_and_reports_which() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(8, 32), 13.0);

    // 40 frames at 8 quads each, against a tree of ~133: enough to drain with room to spare.
    // The first cut used 12 and failed on its own arithmetic, not on the code.
    let mut hits = Vec::new();
    for _ in 0..40 {
        let (hit, spend) = s.step(&f);
        assert!(spend.quads <= 8, "a frame computed {} quads against a quota of 8", spend.quads);
        hits.push(hit);
        if hit == QuotaHit::Drained {
            break;
        }
    }
    assert!(hits.iter().any(|&h| h == QuotaHit::Quads), "the quad quota never bound: {hits:?}");
    assert!(
        hits.last() == Some(&QuotaHit::Drained),
        "the descent never drained in 40 frames: {hits:?}"
    );
}

/// **The tree persists across frames — it is not rebuilt.**
///
/// This is the property the whole phase turns on, and the arm that makes it non-trivial is that
/// the quad count must be **monotone**: a rebuilt tree would restart from the root every frame and
/// the count would sawtooth rather than climb.
#[test]
fn the_tree_persists_and_grows_across_frames() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(8, 2), 13.0);

    let mut counts = Vec::new();
    for _ in 0..10 {
        s.step(&f);
        counts.push(s.stats().quads_computed);
    }
    assert!(counts.windows(2).all(|w| w[1] >= w[0]), "quad count is not monotone: {counts:?}");
    assert!(counts.last().unwrap() > &counts[0], "the tree never grew past frame 1: {counts:?}");
    assert!(counts[0] > 0, "frame 1 computed nothing");

    // And the whole run reaches the same place a one-shot descent does — same code, same rounds.
    let one_shot = prin_rs::scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(8, 2).sched, 13.0, &f);
    let mut s2 = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(8, 2), 13.0);
    for _ in 0..200 {
        if s2.step(&f).0 == QuotaHit::Drained {
            break;
        }
    }
    assert_eq!(
        s2.tree().leaves().count(),
        one_shot.0.leaves().count(),
        "run to exhaustion, the frame loop must reach the one-shot tree"
    );
    assert_eq!(s2.stats().quads_computed, one_shot.1.quads_computed);
}

/// **The quota is a count, not a clock.**
///
/// A sampler that sleeps must not change a single decision. A wall-clock deadline would make the
/// tree a function of machine load, and every result in `results/` is reproducible precisely
/// because nothing in the descent reads a clock.
#[test]
fn the_quota_is_not_a_deadline() {
    let plain = field();
    let slow = move |sl: &Slice, k: usize| {
        std::thread::sleep(std::time::Duration::from_micros(50));
        plain(sl, k)
    };
    let quick = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);

    let run = |samp: &(dyn Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync)| {
        let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(8, 4), 13.0);
        let mut trace = Vec::new();
        for _ in 0..6 {
            let (hit, spend) = s.step(samp);
            trace.push((hit, spend.quads, spend.rounds));
        }
        (trace, s.tree().leaves().count(), s.stats().quads_computed)
    };
    let a = run(&quick);
    let b = run(&slow);
    assert_eq!(a, b, "a slow sampler changed the frame trace -- the quota is reading a clock");
}

/// **A pan touches nothing stored; a zoom touches every stored term.**
///
/// The asymmetry §4.5 asks to be measured separately, and the reason the frontier splits stored
/// from derived at all. `restored == 0` on a pan is the measurement, not an omission — the stored
/// term is position-free by construction.
#[test]
fn a_pan_restores_nothing_and_a_zoom_restores_everything() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(64, 8), 13.0);
    for _ in 0..4 {
        s.step(&f);
    }
    let leaves = s.tree().leaves().count();
    assert!(leaves > 1, "the tree is a single leaf, so there is nothing to restore either way");

    let panned = s.set_camera(Camera::framing(0.3, 0.0, 1.0, 512));
    assert_eq!(panned.restored, 0, "a pan invalidated a stored term");
    assert!(panned.d_centre > 0.0 && panned.d_zoom_octaves == 0.0, "the pan arm is not live");

    let zoomed = s.set_camera(Camera::framing(0.3, 0.0, 0.5, 512));
    assert_eq!(zoomed.restored, leaves, "a zoom must invalidate every stored ranking term");
    assert!(zoomed.d_zoom_octaves > 0.0, "zooming in is a positive octave count");

    // And neither touched the physics: the tree is the same size after both moves.
    assert_eq!(s.tree().leaves().count(), leaves, "a camera move recomputed something");
}

/// **The degenerate cell is refused, and all four cells are tested** — a guard that always fires
/// passes as easily as one that never does.
#[test]
fn the_session_guard_fires_only_on_the_inert_cell() {
    let base = cfg(8, 4);
    let cell = |cap: Option<usize>, regrow: Regrow| SessionCfg {
        payload_cap: cap,
        regrow,
        ..SessionCfg { sched: base.sched.clone(), ..Default::default() }
    };

    // The inert cell, writing to results: refused.
    let inert = cell(None, Regrow::Off);
    assert!(std::panic::catch_unwind(|| {
        assert_session_engages(&inert, "results/session/x.prnf", false);
    })
    .is_err(), "the inert cell was allowed");

    // ...and permitted when named as the control.
    assert_session_engages(&inert, "results/session/x.prnf", true);
    // ...and outside `results/`, where nothing is being claimed.
    assert_session_engages(&inert, "/tmp/scratch/x.prnf", false);

    // The three engaged cells: all allowed.
    for c in [cell(Some(64), Regrow::Off), cell(None, Regrow::Upto(2)), cell(Some(64), Regrow::Upto(2))] {
        assert_session_engages(&c, "results/session/x.prnf", false);
    }
}

/// **Eviction changes no decision — with the arm that says it fired.**
///
/// Under a frozen playhead the argument is structural rather than empirical: `decide` does not
/// take pixels and nothing it calls does, so once `reduce` has run there is no path from a payload
/// to a decision at all. This asserts the consequence *and* prints the evicted count, because a
/// null read off an arm that never engaged is this project's recorded failure mode — a cap that
/// evicted nothing would pass the equality trivially.
#[test]
fn eviction_changes_no_decision_and_the_arm_fires() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);

    let run = |cap: Option<usize>| {
        let mut c = cfg(32, 8);
        c.sched.keep_pixels = true;
        c.payload_cap = cap;
        let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
        let mut evicted = 0usize;
        for _ in 0..40 {
            let (hit, _) = s.step(&f);
            evicted += s.evict_after_frame(&|_| false).len();
            if hit == QuotaHit::Drained {
                break;
            }
        }
        let decisions: Vec<u8> = s.tree().nodes.iter().map(|q| q.decision.code()).collect();
        let boxes: Vec<(u32, u64, u64)> = s
            .tree()
            .nodes
            .iter()
            .map(|q| (q.level, q.cx.to_bits(), q.cy.to_bits()))
            .collect();
        (decisions, boxes, evicted, s.store().resident(), s.stats().quads_computed)
    };

    let (d_keep, b_keep, e_keep, r_keep, q_keep) = run(None);
    let (d_eyes, b_eyes, e_eyes, r_eyes, q_eyes) = run(Some(0));

    // The arm: eviction actually happened, and the two caps really differ.
    assert!(e_eyes > 0, "cap 0 evicted nothing, so the equality below is vacuous");
    assert_eq!(e_keep, 0, "the uncapped arm evicted something");
    assert!(r_keep > r_eyes, "residency did not differ: {r_keep} against {r_eyes}");
    println!("evicted {e_eyes} payloads, resident {r_eyes} against {r_keep}");

    // The claim.
    assert_eq!(d_keep, d_eyes, "eviction changed a decision");
    assert_eq!(b_keep, b_eyes, "eviction changed the tree's shape");
    assert_eq!(q_keep, q_eyes, "eviction changed how many quads were computed");
}

// -------------------------------------------------------------------------------------------
// The frontier, wired into the frame loop.
// -------------------------------------------------------------------------------------------

/// **`priority` IS `stored_term × derived_term`, exactly** — so the frontier and `order_queue`
/// cannot be ranking on two different functions.
///
/// The plan named this as a weak joint with no check: *"the frontier selects on `stored × derived`;
/// `order_queue` then sorts the selected set on `scheduler::priority`"*, and where those differ a
/// frame refines the top-`k` of one ordering in the order of another — an invisible staleness that
/// `agrees_with_rebuild` sits one level below and cannot see. Making `priority` the product by
/// construction removes the class; this asserts it bitwise so a future edit to either half cannot
/// reintroduce it.
#[test]
fn the_two_halves_of_priority_multiply_back_to_it() {
    use prin_rs::scheduler::{derived_term, priority, stored_term};

    let f = field();
    // **Zoomed to a corner, not framing the root.** `Camera::framing` sets `half_world` to the
    // root half-width, which makes every quad fully visible and `relevance` identically 1.0 --
    // so the product identity would hold trivially and the bound would be untested. That fixture
    // defect has now appeared four times on this project; the `nonunit` control below is what
    // catches it rather than a reading of the constructor.
    let base = Camera::framing(0.0, 0.0, 1.0, 512);
    let cam = Camera { cx: 0.5, cy: 0.5, half_world: 0.25, ..base };
    let mut c = cfg(400, 8);
    c.sched.camera = Some(cam);
    c.sched.camera_bias = Some(0.5);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
    for _ in 0..6 {
        s.step(&f);
    }

    let (mut checked, mut nonunit) = (0, 0);
    for i in 0..s.tree().nodes.len() {
        let (st, de) = (
            stored_term(s.tree(), i, s.sched_cfg()),
            derived_term(s.tree(), i, s.sched_cfg()),
        );
        let p = priority(s.tree(), i, s.sched_cfg());
        assert_eq!(st * de, p, "node {i}: {st} * {de} != {p}");
        // **The derived term's bound is the soundness argument for `top_k_bounded`.**
        assert!((0.0..=1.0).contains(&de) || de.is_nan(), "node {i}: derived {de} is outside [0,1]");
        checked += 1;
        nonunit += usize::from(de < 1.0);
    }
    assert!(checked > 20, "only {checked} nodes");
    // The control: if every derived term were 1.0 the product identity would be vacuous.
    assert!(nonunit > 0, "no quad had a derived term below 1, so the camera arm is inert here");
}

/// **A pending quad inherits its PARENT's stored term, and its own would be zero.**
///
/// A freshly-split child carries `QuadReduction::default()`. Ranking on that ranks every child of
/// every parent at exactly zero — *no* ordering, not a weak one, which this project has twice been
/// caught reading a flat curve off. The arm that makes this fire is the comparison against the
/// parent, not the mere fact that the value is non-zero.
#[test]
fn an_uncomputed_quad_ranks_on_its_parent() {
    use prin_rs::scheduler::stored_term;

    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(4, 1), 13.0);
    // Enough frames to have split something and left children pending under the tight quota.
    for _ in 0..4 {
        s.step(&f);
    }
    let pending: Vec<usize> = s.pending().to_vec();
    assert!(!pending.is_empty(), "nothing is pending, so this test has no subject");

    let mut checked = 0;
    for &i in &pending {
        assert_eq!(s.tree().nodes[i].red.n_footprints, 0, "quad {i} is not actually uncomputed");
        let Some(p) = s.tree().nodes[i].parent else { continue };
        assert_eq!(
            stored_term(s.tree(), i, s.sched_cfg()),
            stored_term(s.tree(), p, s.sched_cfg()),
            "pending quad {i} did not inherit parent {p}"
        );
        // And the parent is computed, so the inherited value is real rather than zero twice over.
        assert!(s.tree().nodes[p].red.n_footprints > 0, "parent {p} is uncomputed too");
        checked += 1;
    }
    assert!(checked > 0, "no pending quad had a parent");
}

/// **The frontier walks a fraction of its entries, and reports which fraction.**
///
/// Whether the bucketing earns its place is empirical: if the signal piles into two or three bands
/// the walk degenerates to a full scan and the frontier is a `HashMap` with extra steps. The
/// measurement is committed in `results/frontier/`; this pins that the wiring actually reports it
/// rather than leaving the column at zero.
#[test]
fn the_frame_frontier_reports_its_scan() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    // A quota well below the frontier, or nothing is truncated and the ranking is inert by
    // construction -- the `k_frac = 1.0` cell, at a new site.
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(8, 1), 13.0);
    let mut ever_ranked = false;
    for _ in 0..12 {
        s.step(&f);
        let (scan, len, agrees) = s.frontier_telemetry();
        assert!(scan <= len, "scanned {scan} of {len}");
        if scan > 0 {
            ever_ranked = true;
            assert!(len > 0);
        }
        assert!(agrees.is_nan() || agrees == 1.0, "the frontier disagreed with its rebuild");
    }
    assert!(ever_ranked, "the quota never bound, so the frontier never ranked anything");
}

/// **The audit runs on schedule and is `NaN` otherwise — never `1.0` by default.**
#[test]
fn the_frontier_audit_is_nan_when_it_did_not_run() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut c = cfg(8, 1);
    c.audit_every = 3;
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
    let (mut ran, mut skipped) = (0, 0);
    for _ in 0..9 {
        s.step(&f);
        if s.frontier_telemetry().2.is_nan() { skipped += 1 } else { ran += 1 }
    }
    assert_eq!((ran, skipped), (3, 6), "the audit did not run on exactly every third frame");

    // And `audit_every = 0` means never, said in the type rather than by an accident of modulo.
    let mut c0 = cfg(8, 1);
    c0.audit_every = 0;
    let mut s0 = Session::new(0.0, 0.0, 1.0, 0, cam, c0, 13.0);
    for _ in 0..4 {
        s0.step(&f);
        assert!(s0.frontier_telemetry().2.is_nan(), "audit_every = 0 ran an audit");
    }
}

// -------------------------------------------------------------------------------------------
// Re-rooting, driven by the camera.
// -------------------------------------------------------------------------------------------

/// **A zoom-out past the root box grows the root; a zoom-in never does.**
///
/// The asymmetry is the test. A regrow triggered by a zoom-in would mean the containment test is
/// reading something other than containment, and the count would still look plausible.
#[test]
fn a_zoom_out_grows_the_root_and_a_zoom_in_does_not() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut c = cfg(64, 4);
    c.regrow = Regrow::Upto(3);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
    s.step(&f);
    let root_before = s.tree().root_node().half;

    // Zoom IN: the camera is well inside the root box.
    let d_in = s.set_camera(Camera { half_world: 0.25, ..cam });
    assert_eq!(d_in.regrown, 0, "a zoom-in grew the root");
    assert_eq!(s.tree().root_node().half, root_before);

    // Zoom OUT past the root: 2.5 against a root half of 1.0 needs two doublings.
    let d_out = s.set_camera(Camera { half_world: 2.5, ..cam });
    assert_eq!(d_out.regrown, 2, "expected exactly two doublings, got {}", d_out.regrown);
    assert_eq!(s.tree().root_node().half, root_before * 4.0);
    assert_eq!(s.grown(), 2);

    // And the bound holds: a further zoom-out can add only the one level left of `Upto(3)`.
    let d_far = s.set_camera(Camera { half_world: 40.0, ..cam });
    assert_eq!(d_far.regrown, 1, "the bound did not clamp the growth");
    assert_eq!(s.grown(), 3);
    let d_more = s.set_camera(Camera { half_world: 80.0, ..cam });
    assert_eq!(d_more.regrown, 0, "growth continued past Upto(3)");
}

/// **`Regrow::Off` is the named control, and it must be genuinely off.**
#[test]
fn regrow_off_never_grows() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, cfg(64, 4), 13.0);
    s.step(&f);
    let before = s.tree().root_node().half;
    let d = s.set_camera(Camera { half_world: 32.0, ..cam });
    assert_eq!(d.regrown, 0);
    assert_eq!(s.tree().root_node().half, before, "the root grew with Regrow::Off");
    // The control: the same camera under `Upto` must actually grow, or "off" proves nothing.
    let mut c = cfg(64, 4);
    c.regrow = Regrow::Upto(6);
    let mut s2 = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
    s2.step(&f);
    assert!(s2.set_camera(Camera { half_world: 32.0, ..cam }).regrown > 0,
            "the control did not grow either, so `Off` is untested");
}

/// **Growing the root preserves every existing quad's box, decision and payload.**
///
/// `grow_root` pushes the new root at the end, so no index moves — which is what lets the store
/// and the frontier stay untouched. This asserts the property through the *session*, because the
/// unit test in `tests/regrow.rs` cannot see the store.
#[test]
fn a_regrow_disturbs_neither_the_store_nor_any_decision() {
    let f = field();
    let cam = Camera::framing(0.0, 0.0, 1.0, 512);
    let mut c = cfg(64, 4);
    c.regrow = Regrow::Upto(2);
    c.payload_cap = Some(4096);
    c.sched.keep_pixels = true;
    let mut s = Session::new(0.0, 0.0, 1.0, 0, cam, c, 13.0);
    for _ in 0..4 {
        s.step(&f);
    }
    let n = s.tree().nodes.len();
    let boxes: Vec<(f64, f64, f64)> =
        s.tree().nodes.iter().map(|q| (q.cx, q.cy, q.half)).collect();
    let decisions: Vec<_> = s.tree().nodes.iter().map(|q| q.decision).collect();
    let resident = s.store().resident();

    let d = s.set_camera(Camera { half_world: 3.0, ..cam });
    assert!(d.regrown > 0, "nothing grew, so this test has no subject");
    for i in 0..n {
        assert_eq!((s.tree().nodes[i].cx, s.tree().nodes[i].cy, s.tree().nodes[i].half), boxes[i],
                   "node {i}'s box moved under a regrow");
        assert_eq!(s.tree().nodes[i].decision, decisions[i], "node {i}'s decision moved");
    }
    assert_eq!(s.store().resident(), resident, "a regrow evicted or added a payload");
}
