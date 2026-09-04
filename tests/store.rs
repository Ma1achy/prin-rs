//! Tests for payload residency and eviction.
//!
//! The claims: an evicted quad is distinguishable from one never computed; eviction is a total,
//! deterministic order; the exempt class is honoured; and `Retain::Never` — the one-shot descent's
//! setting — allocates nothing and returns the same sparse shape it always has.

use prin_rs::ensemble::pixel::PixelOut;
use prin_rs::store::{PixelStore, Residency, Retain};

fn px(n: usize) -> Vec<PixelOut> {
    vec![PixelOut::default(); n]
}

/// **`Absent` and `Evicted` must not be the same state, and both read empty to the renderer.**
///
/// A hole from an eviction is a budget fact; a hole from an unreached quad is a coverage fact.
/// Pooling them is the stop-reason conflation one level down. The renderer is *supposed* to treat
/// them alike — both fall through to the coarse ancestor — so the distinction lives in
/// `residency()` and not in `get()`, and this pins both halves.
#[test]
fn an_evicted_quad_is_not_an_absent_one() {
    let mut s = PixelStore::new(Retain::All);
    s.put(3, px(64), 1);

    assert_eq!(s.residency(3), Residency::Resident);
    assert_eq!(s.residency(9), Residency::Absent, "never written");
    assert_eq!(s.get(3).len(), 64);

    let freed = s.evict_to(0, 7, &|_| false, &|_| 0.0, &|_| 4242);
    assert_eq!(freed, vec![3]);

    match s.residency(3) {
        Residency::Evicted { frame, substeps } => {
            assert_eq!(frame, 7, "the frame it was released on");
            assert_eq!(substeps, 4242, "what recomputing it would cost, without recomputing it");
        }
        other => panic!("expected Evicted, got {other:?}"),
    }
    assert_ne!(s.residency(3), s.residency(9), "evicted and absent must be distinguishable");
    // ...and the renderer sees the same thing for both, which is why the fill works unchanged.
    assert!(s.get(3).is_empty() && s.get(9).is_empty());
    assert_eq!(s.resident(), 0);
}

/// **The exempt class is honoured, and the cap can be unsatisfiable.**
///
/// The caching contract pins the baseline cover and the visible ancestor chain: *breaking the
/// visible fallback chain is the one eviction that visibly hurts.* A cap below the size of that
/// class cannot be met, and the store must leave it unmet rather than evict into it — a cap that
/// silently breaks its own exemption reads as a cap that is being honoured.
#[test]
fn the_exempt_class_survives_a_cap_it_cannot_satisfy() {
    let mut s = PixelStore::new(Retain::All);
    for i in 0..10 {
        s.put(i, px(4), 1);
    }
    let pinned = |i: usize| i < 6;

    let freed = s.evict_to(0, 2, &pinned, &|_| 0.0, &|_| 0);
    assert_eq!(freed.len(), 4, "only the four unpinned may go");
    assert_eq!(s.resident(), 6, "the cap of 0 is unsatisfiable and is left unmet");
    for i in 0..6 {
        assert_eq!(s.residency(i), Residency::Resident, "pinned {i} was evicted");
    }
}

/// **Eviction is a total order, so a frame's releases are a function of (tree, camera, frame).**
///
/// Score first, LRU second, index last. Without the final key the order would depend on the
/// enumeration, and the session's memory behaviour would vary run to run — the scheduler-firewall
/// violation one level down. Tested by giving every candidate the *same* score and the same touch
/// time, where only the index can break the tie.
#[test]
fn eviction_is_deterministic_under_ties() {
    let build = || {
        let mut s = PixelStore::new(Retain::All);
        for i in 0..8 {
            s.put(i, px(4), 5);
        }
        s
    };
    let a = build().evict_to(3, 9, &|_| false, &|_| 1.0, &|_| 0);
    let b = build().evict_to(3, 9, &|_| false, &|_| 1.0, &|_| 0);
    assert_eq!(a, b, "two identical stores evicted differently");
    assert_eq!(a, vec![0, 1, 2, 3, 4], "an exact tie must break by index ascending");

    // And the score dominates the index when it differs: high score survives.
    let freed = build().evict_to(3, 9, &|_| false, &|i| i as f64, &|_| 0);
    assert_eq!(freed, vec![0, 1, 2, 3, 4], "the lowest-scoring five must go");

    // LRU is the tie-break under equal score, ahead of the index.
    let mut s = PixelStore::new(Retain::All);
    for i in 0..4 {
        s.put(i, px(4), 100 - i as u64); // slot 3 is the oldest touch
    }
    let freed = s.evict_to(3, 9, &|_| false, &|_| 1.0, &|_| 0);
    assert_eq!(freed, vec![3], "under equal score the least recently touched goes first");
}

/// **`Retain::Never` allocates nothing and returns the shape the one-shot descent always did.**
///
/// This is the batch descent's setting, and it is why `descend_with`'s memory does not move: the
/// store records nothing at all. And `into_dense` must reproduce `resize(i + 1, Vec::new())`
/// exactly — same length, same empty slots — or `SchedStats::pixels` changes shape under every
/// consumer that indexes it.
#[test]
fn retain_never_records_nothing_and_dense_keeps_the_sparse_shape() {
    let mut never = PixelStore::new(Retain::Never);
    never.put(5, px(64), 1);
    assert_eq!(never.residency(5), Residency::Absent, "Never must not record");
    assert_eq!(never.resident(), 0);
    assert!(never.into_dense().is_empty(), "Never must not allocate");

    let mut all = PixelStore::new(Retain::All);
    all.put(2, px(3), 1);
    all.put(5, px(3), 1);
    let dense = all.into_dense();
    assert_eq!(dense.len(), 6, "length is the highest index + 1, as `resize` gave");
    assert!(dense[0].is_empty() && dense[1].is_empty() && dense[3].is_empty() && dense[4].is_empty());
    assert_eq!(dense[2].len(), 3);
    assert_eq!(dense[5].len(), 3);
}

/// A re-put after eviction restores residency and the resident count — the resume path.
#[test]
fn a_recomputed_quad_becomes_resident_again() {
    let mut s = PixelStore::new(Retain::All);
    s.put(1, px(4), 1);
    s.evict_to(0, 2, &|_| false, &|_| 0.0, &|_| 99);
    assert!(matches!(s.residency(1), Residency::Evicted { .. }));
    assert_eq!(s.resident(), 0);

    s.put(1, px(4), 3);
    assert_eq!(s.residency(1), Residency::Resident);
    assert_eq!(s.resident(), 1, "the resident count must not double-count a re-put");
    assert_eq!(s.last_touched(1), 3);
}
