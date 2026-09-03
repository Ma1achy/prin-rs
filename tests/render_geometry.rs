//! The adaptive render's geometry, pinned: orientation, tiling, the coarse fill, and truncation
//! without a mask.
//!
//! Four faults sat in the render path unmeasured. The adaptive render flipped `y` and the
//! uniform panels did not, so every `results/charts/<case>.png` was a vertical mirror of the
//! `<case>_uniform.png` beside it. Sample-centred tiles overhung the quad box by half a cell,
//! so same-level seams were decided by node index and the wire sat half a texel off the colour.
//! The only way to truncate a frame was to empty a node's samples, and that same key was what
//! the coarse-ancestor fill painted from, so a truncated frame could not show a fill. And no
//! test looked at any of it: `tests/colour.rs` asserted frames differ, which they did, in the
//! wrong orientation with cracked seams. Every test here has the arm that shows it can fire.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::PixelOut;
use prin_rs::grid::Slice;
use prin_rs::output::adaptive::{self, TexelMode};
use prin_rs::output::colour::BACKGROUND;
use prin_rs::quad::QuadTree;

const CX: f64 = 1.0;
const CY: f64 = 3.0;
const HALF: f64 = 0.05;

/// A complete tree to `levels`, with `n` samples per axis.
fn complete(levels: u32, n: usize) -> QuadTree {
    let mut t = QuadTree::new(CX, CY, HALF, n, 0);
    let mut frontier = vec![0usize];
    for l in 0..levels {
        let mut next = Vec::new();
        for i in frontier {
            next.extend(t.split(i, l + 1));
        }
        frontier = next;
    }
    t
}

/// Every node carries `N²` footprints whose `ensemble_spread` is its own index, so a colouring
/// keyed on that field paints each quad flat in a colour that names the node.
fn node_pixels(t: &QuadTree) -> Vec<Vec<PixelOut>> {
    (0..t.nodes.len())
        .map(|i| {
            let mut p = PixelOut::default();
            p.ensemble_spread = i as f64;
            vec![p; t.n * t.n]
        })
        .collect()
}

fn node_rgb(p: &PixelOut) -> [u8; 3] {
    let i = p.ensemble_spread as usize;
    [(i * 37 % 251 + 1) as u8, (i * 59 % 251 + 1) as u8, (i * 83 % 251 + 1) as u8]
}

/// The leaf whose box contains the world point: at each level, the child whose centre is
/// nearest. Nearest rather than a containment test, because `Slice::axis` is endpoint-inclusive
/// and its last sample sits exactly on the root edge, where `<= half` against a child edge
/// computed as `cx + half/2 + half/2` can lose by an ulp. Edge pixels are excluded from the
/// comparison anyway; this only has to be total.
fn leaf_at(t: &QuadTree, x: f64, y: f64) -> usize {
    let mut i = 0usize;
    while let Some(kids) = t.nodes[i].children {
        i = *kids
            .iter()
            .min_by(|&&a, &&b| {
                let d = |k: usize| {
                    let q = &t.nodes[k];
                    (x - q.cx).abs() + (y - q.cy).abs()
                };
                d(a).partial_cmp(&d(b)).unwrap()
            })
            .unwrap();
    }
    i
}

/// **Row 0 is the minimum `y`, through the render exactly as through a `Slice` buffer.**
///
/// A per-leaf colouring makes the expected image independent of how the tiles inside a leaf are
/// cut: each pixel's colour is the colour of the leaf holding its centre, looked up in the same
/// `Slice` index order every uniform panel is written in. The render has to match it bitwise.
/// Before the fix it matched the row-reversed image instead — kept as the negative arm.
#[test]
fn the_render_is_the_slice_buffer_row_for_row_and_not_its_mirror() {
    for (levels, n, res) in [(1u32, 4usize, 64usize), (2, 4, 64), (3, 8, 256)] {
        let t = complete(levels, n);
        let px = node_pixels(&t);
        let cam = Camera::framing(CX, CY, HALF, res);
        let (img, _) = adaptive::render(&t, &px, &cam, res, TexelMode::Adaptive, node_rgb);

        // The uniform panel's own path: pixel centres in `Slice` order, written row by row.
        let sl = Slice::body_plane(res, res, CX, CY, HALF, 0);
        let mut expect = Vec::with_capacity(res * res * 3);
        for k in 0..sl.npix() {
            let (x, y) = sl.decode_pos(k);
            expect.extend_from_slice(&node_rgb(&px[leaf_at(&t, x, y)][0]));
        }
        // `Slice::axis` puts its samples at the corners, so a pixel centre on a leaf edge is
        // ambiguous; compare away from edges by testing the interior of every leaf instead of
        // every pixel.
        let mut mism = 0usize;
        let mut mirror = 0usize;
        let mut checked = 0usize;
        for py in 0..res {
            for pxx in 0..res {
                let (x, y) = sl.decode_pos(py * res + pxx);
                let q = &t.nodes[leaf_at(&t, x, y)];
                let margin = 2.0 * cam.pixel_size();
                if (x - q.cx).abs() > q.half - margin || (y - q.cy).abs() > q.half - margin {
                    continue;
                }
                checked += 1;
                let o = (py * res + pxx) * 3;
                let om = ((res - 1 - py) * res + pxx) * 3;
                if img[o..o + 3] != expect[o..o + 3] {
                    mism += 1;
                }
                if img[om..om + 3] == expect[o..o + 3] {
                    mirror += 1;
                }
            }
        }
        assert!(checked > res * res / 2, "the interior check covered too little: {checked}");
        assert_eq!(mism, 0, "levels {levels}, N {n}: {mism} of {checked} interior pixels differ from the slice buffer");
        // The negative arm: the mirrored image must NOT match, or the test could not tell the
        // two conventions apart. (A tree with a symmetric colouring would pass both.)
        assert!(
            mirror < checked / 2,
            "levels {levels}: the row-reversed image also matches ({mirror} of {checked}); the test cannot see orientation"
        );
    }
}

/// **The leaves of a complete tree tile the root: every pixel covered exactly once.**
///
/// At `N = 8` the interior tile edges fall at fractional pixels, which is where the old
/// sample-centred tiles overlapped by up to a pixel and the frame edge painted half a cell of
/// nothing. A mixed-depth tree is the case that matters: a coarse leaf beside four fine ones.
#[test]
fn leaves_cover_every_pixel_exactly_once() {
    for (levels, n, res) in [(2u32, 4usize, 64usize), (3, 8, 128), (3, 8, 200)] {
        let t = complete(levels, n);
        let leaves: Vec<usize> = t.leaves().collect();
        let cam = Camera::framing(CX, CY, HALF, res);
        let cov = adaptive::coverage(&t, &cam, res, &leaves);
        let gaps = cov.iter().filter(|&&c| c == 0).count();
        let overlaps = cov.iter().filter(|&&c| c > 1).count();
        assert_eq!(gaps, 0, "levels {levels} N {n} res {res}: {gaps} pixels uncovered");
        assert_eq!(overlaps, 0, "levels {levels} N {n} res {res}: {overlaps} pixels covered twice");
    }

    // Mixed depth: split only the first child of the root two levels further.
    let mut t = QuadTree::new(CX, CY, HALF, 8, 0);
    let kids = t.split(0, 1);
    let deeper = t.split(kids[0], 2);
    let _ = t.split(deeper[3], 3);
    let leaves: Vec<usize> = t.leaves().collect();
    let cam = Camera::framing(CX, CY, HALF, 160);
    let cov = adaptive::coverage(&t, &cam, 160, &leaves);
    assert_eq!(cov.iter().filter(|&&c| c != 1).count(), 0, "mixed-depth tree does not tile");

    // The negative arm: drop one leaf and the gap appears; duplicate one and the overlap does.
    let cov = adaptive::coverage(&t, &cam, 160, &leaves[1..]);
    assert!(cov.iter().any(|&c| c == 0), "removing a leaf must leave a gap");
    let mut dup = leaves.clone();
    dup.push(leaves[0]);
    let cov = adaptive::coverage(&t, &cam, 160, &dup);
    assert!(cov.iter().any(|&c| c == 2), "repeating a leaf must show as an overlap");
}

/// **A leaf with nothing computed shows its ancestor, never the background — and a leaf with no
/// computed ancestor shows the background, never a fabrication.**
#[test]
fn an_uncomputed_leaf_is_filled_by_its_ancestor_and_only_by_its_ancestor() {
    let t = complete(1, 4);
    let res = 64;
    let cam = Camera::framing(CX, CY, HALF, res);
    let leaves: Vec<usize> = t.leaves().collect();
    let mut px = node_pixels(&t);
    px[leaves[0]].clear(); // one child never computed

    let (img, _) = adaptive::render_leaves(&t, &px, &cam, res, TexelMode::Adaptive, node_rgb, &leaves);
    let root_colour = node_rgb(&px[0][0]);
    let n_root = img.chunks(3).filter(|c| *c == root_colour).count();
    let n_bg = img.chunks(3).filter(|c| *c == BACKGROUND).count();
    assert_eq!(n_bg, 0, "an uncomputed leaf under a computed root left {n_bg} background pixels");
    assert_eq!(n_root, res * res / 4, "the fill covers exactly the uncomputed quarter");

    // Control: with the root uncomputed too, the quarter is background and nothing invents it.
    px[0].clear();
    let (img, _) = adaptive::render_leaves(&t, &px, &cam, res, TexelMode::Adaptive, node_rgb, &leaves);
    let n_bg = img.chunks(3).filter(|c| *c == BACKGROUND).count();
    assert_eq!(n_bg, res * res / 4, "with no ancestor computed the quarter must read as background");
}

/// **A truncated frame needs no mask: the leaf set is the argument.**
///
/// Three frames of a two-level tree — root, the four level-1 quads, the sixteen leaves — must
/// all differ, without emptying anyone's samples; the finished render is the last frame; and a
/// capped frame carries **no texel from below its cap**, which is the fault every animation
/// used to have. (A "shadow tree" with `children = None` at the cap is not a route: the deeper
/// nodes stay in the arena as childless nodes and `leaves()` lists them.)
#[test]
fn truncation_is_the_leaf_set_and_not_a_mask() {
    let t = complete(2, 4);
    let px = node_pixels(&t);
    let res = 64;
    let cam = Camera::framing(CX, CY, HALF, res);
    let at_level = |l: u32| -> Vec<usize> {
        (0..t.nodes.len()).filter(|&i| t.nodes[i].level == l).collect()
    };
    let frames: Vec<Vec<u8>> = (0..=2u32)
        .map(|l| {
            adaptive::render_leaves(&t, &px, &cam, res, TexelMode::Adaptive, node_rgb, &at_level(l)).0
        })
        .collect();
    assert_ne!(frames[0], frames[1], "cap 0 and cap 1 render alike: truncation is not restricting");
    assert_ne!(frames[1], frames[2], "cap 1 and cap 2 render alike: truncation is not restricting");
    let (full, _) = adaptive::render(&t, &px, &cam, res, TexelMode::Adaptive, node_rgb);
    assert_eq!(full, frames[2], "the finished render is the deepest frame");

    // No texel from below the cap, and every texel from the cap: the frame at cap `l` is made
    // of exactly the colours of level `l`.
    for l in 0..=2u32 {
        let colours_at = |lv: u32| -> Vec<[u8; 3]> {
            at_level(lv).iter().map(|&i| node_rgb(&px[i][0])).collect()
        };
        let own = colours_at(l);
        let below: Vec<[u8; 3]> = (l + 1..=2).flat_map(colours_at).collect();
        let (mut n_own, mut n_below) = (0usize, 0usize);
        for c in frames[l as usize].chunks(3) {
            let c = [c[0], c[1], c[2]];
            if own.contains(&c) {
                n_own += 1;
            }
            if below.contains(&c) {
                n_below += 1;
            }
        }
        assert_eq!(n_below, 0, "cap {l} frame carries {n_below} pixels from below the cap");
        assert_eq!(n_own, res * res, "cap {l} frame is not made only of level-{l} texels");
    }
}
