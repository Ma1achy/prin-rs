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
//! Brightness rises with level. Without it a deep tree is a uniform mesh in which no boundary
//! can be attributed to a level, and the whole point of looking is to see *which* levels the
//! budget went to.

/// One leaf's box in **pixel** coordinates, with the level it sits at.
#[derive(Clone, Copy, Debug)]
pub struct Box2 {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub level: u32,
}

/// Blend `rgb` into the pixel at `(x, y)` with weight `a`.
///
/// Blended rather than overwritten so the render underneath stays readable: an opaque wire on a
/// deep tree erases most of the image it is annotating.
fn blend(img: &mut [u8], w: usize, h: usize, x: isize, y: isize, rgb: [u8; 3], a: f64) {
    if x < 0 || y < 0 || x as usize >= w || y as usize >= h {
        return;
    }
    let o = (y as usize * w + x as usize) * 3;
    for k in 0..3 {
        img[o + k] = (img[o + k] as f64 * (1.0 - a) + rgb[k] as f64 * a).round().clamp(0.0, 255.0) as u8;
    }
}

/// Draw leaf boundaries over an RGB8 image, brightness graded by level.
///
/// `max_level` sets the top of the grade; pass the deepest leaf present so the grading uses the
/// full range rather than a fixed one that washes out on a shallow tree.
pub fn draw(img: &mut [u8], w: usize, h: usize, boxes: &[Box2], max_level: u32) {
    let hi = max_level.max(1) as f64;
    for b in boxes {
        let t = (b.level as f64 / hi).clamp(0.0, 1.0);
        // Deep boxes brighter, and slightly more opaque, so the fine mesh reads as fine rather
        // than as noise.
        let rgb = [
            (70.0 + 185.0 * t) as u8,
            (80.0 + 175.0 * t) as u8,
            (110.0 + 145.0 * t) as u8,
        ];
        let a = 0.35 + 0.45 * t;

        let (x0, y0) = (b.x0.round() as isize, b.y0.round() as isize);
        let (x1, y1) = ((b.x1 - 1.0).round() as isize, (b.y1 - 1.0).round() as isize);
        if x1 < x0 || y1 < y0 {
            continue;
        }
        for x in x0..=x1 {
            blend(img, w, h, x, y0, rgb, a);
            blend(img, w, h, x, y1, rgb, a);
        }
        for y in y0..=y1 {
            blend(img, w, h, x0, y, rgb, a);
            blend(img, w, h, x1, y, rgb, a);
        }
    }
}

/// Leaf boxes of a live [`crate::quad::QuadTree`], projected through a camera.
/// Boxes for an **explicit** leaf set.
///
/// [`boxes_from_tree`] walks `tree.leaves()`, which is every node whose `children` is `None`.
/// For a finished tree that is the leaf set. For a **truncated** one — an animation frame built
/// by capping a completed descent — it is not: the deep quads outside the cap were already
/// leaves, so their boxes come along and the wireframe shows the finished tree in every frame.
///
/// That is exactly what was shipping. The colour frames had the same fault one layer down, in
/// `adaptive::render`; both are the same mistake, that a truncated view of a finished tree has to
/// name its own leaf set rather than infer one from `children`.
pub fn boxes_from_leaves(
    tree: &crate::quad::QuadTree,
    cam: &crate::camera::Camera,
    res: usize,
    leaves: &[usize],
) -> Vec<Box2> {
    let px = cam.pixel_size();
    let to = |x: f64, y: f64| -> (f64, f64) {
        (
            (x - cam.cx) / px + res as f64 / 2.0,
            res as f64 / 2.0 - (y - cam.cy) / px,
        )
    };
    leaves
        .iter()
        .map(|&i| {
            let q = &tree.nodes[i];
            let (ax, ay) = to(q.cx - q.half, q.cy + q.half);
            let (bx, by) = to(q.cx + q.half, q.cy - q.half);
            Box2 { x0: ax, y0: ay, x1: bx, y1: by, level: q.level }
        })
        .collect()
}

pub fn boxes_from_tree(
    tree: &crate::quad::QuadTree,
    cam: &crate::camera::Camera,
    res: usize,
) -> Vec<Box2> {
    let px = cam.pixel_size();
    let to = |x: f64, y: f64| -> (f64, f64) {
        (
            (x - cam.cx) / px + res as f64 / 2.0,
            res as f64 / 2.0 - (y - cam.cy) / px,
        )
    };
    tree.leaves()
        .map(|i| {
            let q = &tree.nodes[i];
            let (ax, ay) = to(q.cx - q.half, q.cy + q.half);
            let (bx, by) = to(q.cx + q.half, q.cy - q.half);
            Box2 { x0: ax, y0: ay, x1: bx, y1: by, level: q.level }
        })
        .collect()
}
