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
