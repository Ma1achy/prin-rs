//! **Growing the root upward: the old subtree survives bit for bit, and the veto does not move.**
//!
//! A zoom-out past the root box needs area the tree does not contain. Re-rooting and discarding
//! would throw away exactly the quads the zoom-out is about to display — turning the one gesture
//! the caching contract calls *nearly free* into the most expensive in the system. So the root
//! grows, and the two things that could go silently wrong are the geometry and the depth cap.

use prin_rs::camera::Camera;
use prin_rs::grid::Chart;
use prin_rs::quad::QuadTree;

fn built() -> QuadTree {
    let mut t = QuadTree::with_chart(1.0, 3.0, 0.05, 8, 0, Chart::BodyPlane);
    t.split(0, 0);
    let kids = t.nodes[0].children.unwrap();
    t.split(kids[0], 1);
    t.split(kids[3], 1);
    t
}

/// **Every existing box is preserved exactly, and no index moves.**
///
/// Re-deriving the old subtree's boxes from the new root would shift the whole tree by an ulp —
/// `old_cx + old_half - old_half` is not `old_cx` in f64 — which is the half-cell class of defect
/// `Cache::key_of` carries a paragraph about. `grow_root` does not touch them at all, and the
/// constructor asserts `child_boxes()[quadrant]` reproduces the old root's box before committing.
#[test]
fn growing_the_root_preserves_every_box_and_index() {
    for quadrant in 0..4 {
        let before = built();
        let boxes: Vec<(f64, f64, f64)> =
            before.nodes.iter().map(|q| (q.cx, q.cy, q.half)).collect();
        let kids: Vec<Option<[usize; 4]>> = before.nodes.iter().map(|q| q.children).collect();
        let n_before = before.nodes.len();

        let mut t = built();
        let made = t.grow_root(quadrant, 7);

        assert_eq!(made.len(), 3, "three new siblings, not four");
        assert_eq!(t.root, n_before, "the new root is pushed at the END, so no index moves");
        for (i, &(cx, cy, half)) in boxes.iter().enumerate() {
            assert_eq!((t.nodes[i].cx, t.nodes[i].cy, t.nodes[i].half), (cx, cy, half),
                       "quadrant {quadrant}: node {i}'s box moved -- bitwise equality is the bar");
            assert_eq!(t.nodes[i].children, kids[i], "node {i}'s children were remapped");
        }

        // The old root is the named child, and the new root's own box halves back onto it.
        let old_root_box = boxes[0];
        let rb = t.root_node().child_boxes()[quadrant];
        assert_eq!(rb, old_root_box, "the new root does not reproduce the old root's box");
        assert_eq!(t.nodes[t.root].children.unwrap()[quadrant], 0);
        assert_eq!(t.nodes[0].parent, Some(t.root));
        assert_eq!(t.nodes[0].sib_index, quadrant as u8);
    }
}

/// **Every level rises by one, which is a semantic change and not bookkeeping.**
///
/// A leaf that was inside `bootstrap_levels` may no longer be. Recorded rather than absorbed —
/// and the tree's *shape* is unchanged, which is what says the shift is uniform.
#[test]
fn every_level_rises_by_exactly_one_and_the_shape_holds() {
    let before = built();
    let levels: Vec<u32> = before.nodes.iter().map(|q| q.level).collect();
    let leaves_before = before.leaves().count();

    let mut t = built();
    t.grow_root(2, 0);

    for (i, &l) in levels.iter().enumerate() {
        assert_eq!(t.nodes[i].level, l + 1, "node {i}'s level did not rise by exactly one");
    }
    assert_eq!(t.nodes[t.root].level, 0, "the new root is level 0");
    // Three new leaves join; the old leaves are all still leaves.
    assert_eq!(t.leaves().count(), leaves_before + 3);
}

/// **`Camera::veto` is invariant under a re-root**, and that is asserted rather than assumed.
///
/// Both `q.level` and `floor(camera_depth)` rise together, because the root box doubles — so
/// `q.level >= camera_depth + max_rel_depth` is preserved. If it were not, a zoom-out would change
/// the depth cap of quads whose physics did not move, which is precisely what the position-free
/// veto exists to prevent.
#[test]
fn the_veto_is_invariant_under_a_regrow() {
    let before = built();
    // **A 32-pixel viewport, so the screen floor actually bites on this three-level tree.** At
    // 512 nothing is vetoed at all and the invariance below passes on an empty set — which the
    // control at the end of this test caught on the first run.
    let cam = Camera::framing(1.0, 3.0, 0.05, 32);
    let root_half_before = before.root_node().half;
    let vetoes_before: Vec<Option<u8>> = before
        .nodes
        .iter()
        .map(|q| cam.veto(q, before.n, root_half_before).map(|d| d.code()))
        .collect();

    let mut t = built();
    t.grow_root(1, 0);
    let root_half_after = t.root_node().half;
    assert_eq!(root_half_after, root_half_before * 2.0, "the root box must double");

    let mut checked = 0;
    for (i, want) in vetoes_before.iter().enumerate() {
        let got = cam.veto(&t.nodes[i], t.n, root_half_after).map(|d| d.code());
        assert_eq!(&got, want, "node {i}'s veto moved under a re-root");
        checked += 1;
    }
    assert!(checked > 4, "only {checked} nodes were checked; the fixture is too small");

    // The control: the veto is NOT trivially constant — some node must actually be vetoed, or
    // this test would pass on a camera that vetoes nothing.
    assert!(
        vetoes_before.iter().any(|v| v.is_some()),
        "no node was vetoed at all, so the invariance above is vacuous"
    );
}

/// **`neighbour` still descends from the right place**, and finds across the new siblings.
#[test]
fn neighbour_follows_the_new_root() {
    use prin_rs::quad::Dir;

    let mut t = built();
    t.grow_root(0, 0);

    // The old root sits in quadrant 0 of the new root; its siblings are new leaves. A leaf on the
    // old root's outer edge must now find one of them rather than falling off the tree.
    let found = t
        .leaves()
        .filter(|&i| Dir::ALL.iter().any(|&d| t.neighbour(i, d).is_some()))
        .count();
    assert!(found > 0, "no leaf found any neighbour after a re-root");

    // And every neighbour is same-or-coarser and genuinely adjacent, as before.
    for i in t.leaves() {
        for d in Dir::ALL {
            if let Some(j) = t.neighbour(i, d) {
                assert_ne!(i, j);
                assert!(t.nodes[j].level <= t.nodes[i].level, "neighbour is finer than the query");
                let (a, b) = (&t.nodes[i], &t.nodes[j]);
                let gap = (a.cx - b.cx).abs().max((a.cy - b.cy).abs());
                assert!(gap <= a.half + b.half + 1e-9, "neighbour {j} does not touch {i}");
            }
        }
    }
}
