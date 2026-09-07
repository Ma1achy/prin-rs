//! **UV addressing: the exactness claim that is true here and false in chart space.**
//!
//! §12's *"`h` is an exact power of two at every depth"* was read as an arithmetic claim about
//! this repo and objected to; the objection was in the wrong coordinate system. The claim is about
//! UV, where the root half is exactly `1/2`. What is true of chart space is the weaker
//! `w == prev / 2.0`, which `tests/quadtree.rs` already pins. Both are asserted here, **side by
//! side**, because the pair is the finding.

use prin_rs::uv::{axis_local, uv_half, QuadId, UvFrame};

/// **In UV a cell half-width is a strict power of two; in chart space it is not.**
#[test]
fn uv_widths_are_powers_of_two_and_chart_widths_are_not() {
    for d in 0..40i32 {
        let h = uv_half(d);
        assert_eq!(h, 2f64.powi(-(d + 1)), "uv_half({d}) is not 2^-(d+1)");
        // Strict: the exponent is an integer AND the mantissa is 1.
        assert_eq!(h.log2(), f64::from(-(d + 1)), "uv_half({d}).log2() is not an integer");
        assert_eq!(h.to_bits() & 0x000F_FFFF_FFFF_FFFF, 0, "uv_half({d}) has a mantissa");
    }

    // The chart-space control, on the tree's real root half. `0.05` is not a power of two, so
    // halving it stays exact and never becomes one -- which is exactly why an integer address
    // needs the other space.
    let mut w = 0.05f64;
    let mut ever_pow2 = false;
    for _ in 0..40 {
        let prev = w;
        w /= 2.0;
        assert_eq!(w, prev / 2.0, "chart halving is not exact");
        ever_pow2 |= w.to_bits() & 0x000F_FFFF_FFFF_FFFF == 0;
    }
    assert!(!ever_pow2, "a chart width became a power of two, so this control proves nothing");
}

/// **The address round-trips, at every depth and off-centre.**
#[test]
fn a_box_round_trips_through_its_address() {
    let f = UvFrame::new(1.0, 3.0, 0.05);
    let mut checked = 0;
    for depth in 0..24i32 {
        for &(tx, ty) in &[(0i64, 0i64), (1, 0), (0, 1), ((1i64 << depth) - 1, (1i64 << depth) / 2)]
        {
            let id = QuadId { depth, tx, ty };
            let (cx, cy, half) = f.box_of(id);
            assert_eq!(f.id_of(cx, cy, half), Some(id), "round trip failed at {id:?}");
            checked += 1;
        }
    }
    assert!(checked > 80, "only {checked} addresses checked");
}

/// **The half-cell trap is REFUSED, not rounded.**
///
/// `metric::Cache::key_of` scored a perfectly coherent leaf set belonging to a *shifted* tree
/// because `(2i+1)h / 2h = i + 0.5` rounds to `i + 1`. A box offset by half a cell is not on the
/// lattice and must come back `None`, and the arm that makes that mean something is the aligned
/// box beside it returning `Some`.
#[test]
fn a_half_cell_offset_is_refused_rather_than_rounded() {
    let f = UvFrame::new(1.0, 3.0, 0.05);
    let id = QuadId { depth: 5, tx: 9, ty: 17 };
    let (cx, cy, half) = f.box_of(id);
    assert_eq!(f.id_of(cx, cy, half), Some(id), "the aligned control must resolve");

    // Half a cell right: the classic wrong answer is `tx + 1` with everything else intact.
    assert_eq!(f.id_of(cx + half, cy, half), None, "a half-cell offset resolved to an address");
    assert_eq!(f.id_of(cx, cy + half, half), None);
    // A width off the ladder.
    assert_eq!(f.id_of(cx, cy, half * 0.75), None, "an off-ladder width resolved");
    // And the degenerate inputs.
    assert_eq!(f.id_of(cx, cy, 0.0), None);
    assert_eq!(f.id_of(f64::NAN, cy, half), None);
}

/// **Parent and children compose, including across zero where a truncating division is wrong.**
#[test]
fn parent_and_children_compose_and_negative_indices_work() {
    for &id in &[
        QuadId { depth: 3, tx: 5, ty: 2 },
        QuadId { depth: 0, tx: 0, ty: 0 },
        // A grown root: depth -1 is twice the frame, and -1 is the cell left of it.
        QuadId { depth: -1, tx: -1, ty: 0 },
        QuadId { depth: 2, tx: -3, ty: -1 },
    ] {
        for (j, c) in id.children().iter().enumerate() {
            assert_eq!(c.parent(), id, "child {j} of {id:?} does not name it back");
        }
    }
    // The arm that catches a truncating `/2`: `-1 / 2 == 0` in Rust, so the cell left of the
    // frame would claim the frame's own parent and two disjoint subtrees would merge.
    let left = QuadId { depth: 1, tx: -1, ty: 0 };
    assert_eq!(left.parent().tx, -1, "a truncating division would give 0 here");
    assert_ne!(left.parent(), QuadId { depth: 0, tx: 0, ty: 0 });
}

/// **A GROWN ROOT LEAVES THE ABSOLUTE LATTICE IN THREE DIRECTIONS OF FOUR, AND THAT IS THE SEAM
/// BETWEEN THE ROOTED TREE AND THE FLAT STORE.**
///
/// The caching contract wants a flat store keyed by an *absolute* `QuadID`; this repo has a
/// *rooted* tree that `grow_root` doubles toward the camera. A cell has exactly **one** parent in
/// an absolute lattice, so only one of the four growth directions produces a box that is on it.
/// The tree's initial root is the frame, address `(0, 0, 0)`, whose lattice parent is the cell to
/// its **upper-right** — reached by putting the old root at `child_boxes[0]`, the lower-left.
///
/// So a zoom-out toward the lower-left grows a perfectly correct tree whose new root **has no
/// address**. That is not a defect in either construction: it is the incompatibility the plan
/// named as *"keep the rooted tree and grow it, or build UV addressing first and let the store be
/// flat"*, as a number. `id_of` returns `None` rather than a plausible wrong answer, which is the
/// same refusal the half-cell guard makes.
#[test]
fn a_grown_root_addresses_only_in_the_lattice_direction() {
    use prin_rs::grid::Chart;
    use prin_rs::quad::QuadTree;

    let f = UvFrame::new(1.0, 3.0, 0.05);
    let t = QuadTree::with_chart(1.0, 3.0, 0.05, 8, 0, Chart::BodyPlane);
    assert_eq!(
        f.id_of(t.nodes[0].cx, t.nodes[0].cy, t.nodes[0].half),
        Some(QuadId { depth: 0, tx: 0, ty: 0 }),
        "the frame's own box must be depth 0 at the origin"
    );

    let mut on = 0;
    let mut off = 0;
    for quadrant in 0..4 {
        let mut t2 = t.clone();
        t2.grow_root(quadrant, 0);
        let r = t2.root_node();
        match f.id_of(r.cx, r.cy, r.half) {
            Some(id) => {
                assert_eq!(quadrant, 0, "quadrant {quadrant} resolved; only the lattice one may");
                assert_eq!(id.depth, -1, "a doubled root is depth -1");
                assert_eq!(id, QuadId { depth: -1, tx: 0, ty: 0 });
                assert_eq!(
                    id.children()[quadrant],
                    QuadId { depth: 0, tx: 0, ty: 0 },
                    "the old root is not the child grow_root named"
                );
                on += 1;
            }
            None => {
                assert_ne!(quadrant, 0, "the lattice direction must resolve");
                off += 1;
            }
        }
    }
    // Both arms must have fired, or this is a test that cannot fail.
    assert_eq!((on, off), (1, 3), "expected exactly one lattice direction of four");
}

/// **The seed depends on the address and nothing else**, so it survives eviction and a re-root.
#[test]
fn the_seed_is_a_function_of_the_address_alone() {
    let a = QuadId { depth: 7, tx: 33, ty: 91 };
    assert_eq!(a.seed(0), a.seed(0), "not deterministic");
    assert_ne!(a.seed(0), a.seed(1), "the salt does nothing");

    // Neighbouring addresses must not collide, and a swap of tx/ty must not either — a seed that
    // hashed `tx + ty` would pass a determinism check and alias the whole diagonal.
    let mut seen = std::collections::HashSet::new();
    for d in 0..8i32 {
        for tx in -4..12i64 {
            for ty in -4..12i64 {
                assert!(seen.insert(QuadId { depth: d, tx, ty }.seed(0xA5)),
                        "seed collision at ({d}, {tx}, {ty})");
            }
        }
    }
    assert_eq!(seen.len(), 8 * 16 * 16);
}

/// **`axis_local` spans `[-1, 1]` exactly at both ends and is monotone.**
///
/// And it is **not** bitwise `(axis(c, h, n, i) - c) / h`, which is the point rather than a defect:
/// the two orders of operation round differently, and the local one is the one that has not yet
/// thrown the information away. `examples/sample_space.rs` measures the gap.
#[test]
fn the_local_axis_is_exact_at_both_ends() {
    for n in [2usize, 3, 8, 64, 1024] {
        assert_eq!(axis_local(n, 0), -1.0, "n={n} does not start at -1");
        assert_eq!(axis_local(n, n - 1), 1.0, "n={n} does not end at +1");
        let mut prev = f64::NEG_INFINITY;
        for i in 0..n {
            let x = axis_local(n, i);
            assert!(x > prev, "n={n} is not strictly monotone at {i}");
            assert!((-1.0..=1.0).contains(&x));
            prev = x;
        }
    }
    assert_eq!(axis_local(1, 0), -1.0);
}
