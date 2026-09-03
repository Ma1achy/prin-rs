//! Quadtree wireframe: leaf boundaries drawn over a render.
//!
//! The adaptive render shows **what the system displays**. The wireframe shows **where the tree
//! put its boundaries**. They are different questions and for a criterion study both are wanted:
//! a coarse texel tells you a leaf is coarse, and only the wire tells you whether the
//! surrounding structure was subdivided *around* it or straight *through* it.
//!
//! PR #11 drew boundaries over a **uniform** base, which conflated the two and hid the failure it
//! was supposed to expose. These draw over the adaptive render, so the texel size and the box it
//! belongs to are visible together, and every image has a `_wire` twin rather than the wire
//! replacing it — one is not a substitute for the other.
//!
//! # Depth grading, and why it is not decoration
//!
//! Opacity rises with level. Without it a deep tree is a uniform mesh in which no boundary can
//! be attributed to a level, and the whole point of looking is to see *which* levels the budget
//! went to. `max_level` must be the deepest level of the **finished** tree, held fixed across
//! the frames of an animation: two harnesses passed a literal `1`, which clamps every level to
//! full brightness and deletes the grading, and one passed the frame's own cap, which regrades
//! every frame so a level-3 box is bright in frame 3 and dim in frame 6 — the ramp moving
//! rather than the tree.
//!
//! # Contrast is taken from the pixel underneath
//!
//! A fixed grey-blue line has no contrast on a pale field: `body_plane` at `t = 13` is flat cyan
//! at lightness 0.9 and its wire was invisible. The line is dark on a light pixel and light on a
//! dark one, decided per pixel from the render's own luma, with the level carried by opacity and
//! a small tint so deep and shallow still read apart.
//!
//! Boxes use the crate's one projection, [`crate::camera::Camera::to_px`], so the wire sits on
//! the tile edges it describes — the old copies here flipped `y` and the render did not.

/// One leaf's box in **pixel** coordinates, `x0 < x1`, `y0 < y1`, with the level it sits at.
#[derive(Clone, Copy, Debug)]
pub struct Box2 {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub level: u32,
}

/// Blend a level-`t` line pixel into `(x, y)`, choosing dark or light against what is there.
///
/// Blended rather than overwritten so the render underneath stays readable: an opaque wire on a
/// deep tree erases most of the image it is annotating.
fn mark(img: &mut [u8], w: usize, h: usize, x: isize, y: isize, t: f64) {
    if x < 0 || y < 0 || x as usize >= w || y as usize >= h {
        return;
    }
    let o = (y as usize * w + x as usize) * 3;
    let luma = 0.299 * img[o] as f64 + 0.587 * img[o + 1] as f64 + 0.114 * img[o + 2] as f64;
    // Deep boxes more opaque, so the fine mesh reads as fine rather than as noise. The tint is
    // the same hue either way — bluer when deep — so a level can be read off either polarity.
    let a = 0.40 + 0.45 * t;
    let rgb: [f64; 3] = if luma > 128.0 {
        [30.0 + 30.0 * (1.0 - t), 30.0 + 25.0 * (1.0 - t), 70.0 + 40.0 * t]
    } else {
        [215.0 + 40.0 * t, 220.0 + 35.0 * t, 235.0 + 20.0 * t]
    };
    for k in 0..3 {
        img[o + k] = (img[o + k] as f64 * (1.0 - a) + rgb[k] * a).round().clamp(0.0, 255.0) as u8;
    }
}

/// Draw leaf boundaries over an RGB8 image, opacity graded by level.
///
/// `max_level` sets the top of the grade; pass the deepest leaf of the finished tree, and the
/// same value for every frame of an animation.
pub fn draw(img: &mut [u8], w: usize, h: usize, boxes: &[Box2], max_level: u32) {
    let hi = max_level.max(1) as f64;
    for b in boxes {
        let t = (b.level as f64 / hi).clamp(0.0, 1.0);
        let (x0, y0) = (b.x0.round() as isize, b.y0.round() as isize);
        let (x1, y1) = ((b.x1.round() - 1.0) as isize, (b.y1.round() - 1.0) as isize);
        if x1 < x0 || y1 < y0 {
            continue;
        }
        for x in x0..=x1 {
            mark(img, w, h, x, y0, t);
            mark(img, w, h, x, y1, t);
        }
        for y in y0..=y1 {
            mark(img, w, h, x0, y, t);
            mark(img, w, h, x1, y, t);
        }
    }
}

/// Boxes for an **explicit** leaf set, projected through a camera.
///
/// [`boxes_from_tree`] walks `tree.leaves()`, which is every node whose `children` is `None`.
/// For a finished tree that is the leaf set. For a **truncated** one — an animation frame built
/// by capping a completed descent — it is not: the deep quads outside the cap were already
/// leaves, so their boxes come along and the wireframe shows the finished tree in every frame.
/// A truncated view of a finished tree has to name its own leaf set rather than infer one from
/// `children`, and so does the colour frame beside it (`adaptive::render_leaves`).
pub fn boxes_from_leaves(
    tree: &crate::quad::QuadTree,
    cam: &crate::camera::Camera,
    res: usize,
    leaves: &[usize],
) -> Vec<Box2> {
    leaves
        .iter()
        .map(|&i| {
            let q = &tree.nodes[i];
            let (ax, ay) = cam.to_px(res, q.cx - q.half, q.cy - q.half);
            let (bx, by) = cam.to_px(res, q.cx + q.half, q.cy + q.half);
            Box2 { x0: ax, y0: ay, x1: bx, y1: by, level: q.level }
        })
        .collect()
}

/// Leaf boxes of a finished [`crate::quad::QuadTree`], projected through a camera.
pub fn boxes_from_tree(
    tree: &crate::quad::QuadTree,
    cam: &crate::camera::Camera,
    res: usize,
) -> Vec<Box2> {
    let leaves: Vec<usize> = tree.leaves().collect();
    boxes_from_leaves(tree, cam, res, &leaves)
}
