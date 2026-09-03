//! The adaptive render: **texels at their true per-quad sizes**.
//!
//! PR #11's overlay drew leaf boundaries over a *uniform* render, so every texel was the same
//! size. That is the wrong instrument — it shows where boundaries fell, not what the system
//! displays. A leaf at level 3 must be drawn with 4x the linear texel size of a leaf at
//! level 5, or the tree's quality cannot be judged by eye at all.
//!
//! **One sample, one tile, no interpolation.** A coarse quad is never upsampled to fill pixels
//! smoothly; that fabricates structure it does not have, which on a chaos instrument is the
//! one thing a picture must not do.
//!
//! # Orientation: the image is the array
//!
//! Row 0 of every image this crate writes is the **minimum** `y` of the slice, because
//! `Slice::decode_pos` runs `idx = jy*nx + jx` with `y` increasing with `jy`, and thirty-odd
//! harnesses hand that buffer to [`save_rect`] in index order. This render, the wireframe and
//! `tree::overlay` used to flip — world-up as image-top — so every adaptive panel in
//! `results/charts/` was a vertical mirror of the `_uniform` panel beside it, and the two
//! conventions sat in one directory as a comparison pair. There is now one projection,
//! [`Camera::to_px`], and `tests/render_geometry.rs` pins that a leaf at the minimum `y` lands
//! on row 0 through every path.
//!
//! # Tiles tile the quad, and the quad tiles its parent
//!
//! `Slice::axis` is endpoint-inclusive, so a quad's `N` samples run corner to corner at spacing
//! `2h/(N-1)`. A sample-centred tile of that width would overhang the box by half a cell on each
//! side, and the old render painted exactly that, coarsest first, relying on a finer neighbour
//! to overwrite the overhang: at a same-level seam the winner was whichever node had the higher
//! index, at the frame edge nothing overwrote it, and the wireframe (which draws the exact box)
//! sat half a texel off the colour everywhere. Each sample now paints its cell **clipped to the
//! quad box** — a full cell in the interior, a half cell at the edge — with every edge computed
//! by one expression and rounded once, so neighbouring tiles share their boundary exactly and
//! the leaves of a complete tree cover the root with no gap, no overlap, and no ancestor showing
//! through. [`coverage`] is the instrument that measures it.
//!
//! **The screen floor and this raster disagree by `N/(N-1)`** — 14.3% at `N = 8`.
//! `Camera::tile_size_px` is `2h/N` per pixel, the nominal tile; the painted cell is `2h/(N-1)`.
//! Left as it is, stated: moving the floor to the painted width would push the everyday stop one
//! level deeper at the standard viewport (level 6 to 7 at `N = 8` on 512²), which is a regime
//! change to every committed tree and not a rendering fix. The floor firing slightly early is
//! the conservative direction.
//!
//! # The visible set is named, never inferred
//!
//! [`render`] draws the finished tree: every leaf, with every ancestor painted underneath it,
//! coarsest first — the coarse-ancestor fill of §4.5, which is bitwise inert wherever the tree
//! is complete and is what stops an uncomputed leaf from reading as a hole. [`render_leaves`]
//! draws an **explicit** leaf set the same way, and is what an animation frame must use. The
//! harnesses used to truncate a frame by emptying the samples of every node outside the set,
//! because the old render keyed painting on "has samples" — and that same key was the fill, so
//! a truncated frame could not show a fill. The two were never separable until the set was an
//! argument.

use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use crate::camera::Camera;
use crate::ensemble::pixel::PixelOut;
use crate::quad::QuadTree;

/// What one leaf contributed to the image, so the render can be *measured* and not just looked
/// at. `texel_px` is the linear size actually used, not one derived from the level — which is
/// what lets [`texel_scaling`] tell an adaptive render from a uniform one.
#[derive(Clone, Copy, Debug)]
pub struct LeafTexel {
    pub node: usize,
    pub level: u32,
    pub texel_px: f64,
    pub tiles_drawn: usize,
}

/// Fitted exponent of texel size against level: `log2(texel_px) = a + slope * level`.
///
/// **This is the acceptance test, and it can fail.** An adaptive render gives exactly `-1`: a
/// level-3 leaf's texels are 4x a level-5 leaf's. A uniform render gives `0`, because every
/// texel is the same size whatever the level — the PR #11 failure, which now has an assertion
/// that fires on it. Returns `None` when every leaf is at one level, where the fit is
/// undefined and a number would be an invention.
pub fn texel_scaling(t: &[LeafTexel]) -> Option<f64> {
    let pts: Vec<(f64, f64)> = t
        .iter()
        .filter(|x| x.texel_px > 0.0 && x.texel_px.is_finite())
        .map(|x| (x.level as f64, x.texel_px.log2()))
        .collect();
    if pts.len() < 2 {
        return None;
    }
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    if sxx <= 0.0 {
        return None; // one level only: no variation to fit against
    }
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    Some(sxy / sxx)
}

/// Whether each leaf rasterises at its own texel size (`true`) or all at one size (`false`).
///
/// `Uniform` is not a rendering mode anyone wants; it exists so the acceptance test has
/// something to reject, and so PR #11's instrument stays reproducible.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TexelMode {
    Adaptive,
    /// Every leaf drawn at the finest leaf's texel size — PR #11's uniform render.
    Uniform,
}

/// The nodes to paint for a leaf set: the set and every ancestor of it, coarsest first.
///
/// Ancestors are what the coarse fill draws under a leaf that has no samples yet. Descendants
/// of the set are *not* included — that is the whole difference from the old "every node with
/// samples" rule, and it is what makes a truncated frame a truncated frame.
fn paint_order(tree: &QuadTree, leaves: &[usize]) -> Vec<usize> {
    let mut set: HashSet<usize> = HashSet::with_capacity(leaves.len() * 2);
    for &i in leaves {
        let mut j = Some(i);
        while let Some(k) = j {
            if !set.insert(k) {
                break;
            }
            j = tree.nodes[k].parent;
        }
    }
    let mut v: Vec<usize> = set.into_iter().collect();
    v.sort_by_key(|&i| (tree.nodes[i].level, i));
    v
}

/// One tile's pixel box, `[x0, x1) x [y0, y1)`, clamped to the raster; `None` when it is empty.
///
/// The tile is sample `(jx, jy)`'s cell **clipped to the quad box**. Every edge between two
/// tiles is one expression — `c - h + (k - 0.5) * cell` for the `k`th interior edge — evaluated
/// once per edge and rounded once, so the tiles of one quad partition its pixel box exactly.
/// Two sibling quads share their common edge as `c ± h` of each, which agree to an ulp; a
/// half-pixel landing within that ulp is the one case a seam can move by a pixel.
#[allow(clippy::too_many_arguments)]
fn tile_px(
    q: &crate::quad::Quad,
    n: usize,
    cam: &Camera,
    res: usize,
    mode: TexelMode,
    deepest: u32,
    jx: usize,
    jy: usize,
) -> Option<(usize, usize, usize, usize)> {
    let cell = 2.0 * q.half / (n - 1) as f64;
    let bounds = |c: f64, j: usize| -> (f64, f64) {
        match mode {
            TexelMode::Adaptive => {
                let lo = if j == 0 { c - q.half } else { c - q.half + (j as f64 - 0.5) * cell };
                let hi =
                    if j + 1 == n { c + q.half } else { c - q.half + (j as f64 + 0.5) * cell };
                (lo, hi)
            }
            // Every leaf at the finest leaf's texel size, centred on the sample — the wrong
            // instrument, kept reproducible so the acceptance test has a negative case.
            TexelMode::Uniform => {
                let dc = 2.0 * (q.half * (2f64).powi(q.level as i32 - deepest as i32))
                    / (n - 1) as f64;
                let s = c - q.half + j as f64 * cell;
                ((s - dc / 2.0).max(c - q.half), (s + dc / 2.0).min(c + q.half))
            }
        }
    };
    let (wx0, wx1) = bounds(q.cx, jx);
    let (wy0, wy1) = bounds(q.cy, jy);
    let (ax, ay) = cam.to_px(res, wx0, wy0);
    let (bx, by) = cam.to_px(res, wx1, wy1);
    let clamp = |v: f64| v.round().clamp(0.0, res as f64) as usize;
    let (x0, x1, y0, y1) = (clamp(ax), clamp(bx), clamp(ay), clamp(by));
    if x1 <= x0 || y1 <= y0 {
        None
    } else {
        Some((x0, x1, y0, y1))
    }
}

/// Rasterise the finished tree into an RGB buffer of `res × res`: every leaf, with its
/// ancestors filled in underneath.
///
/// `pixels[node]` holds that quad's `N²` footprints in `idx = jy*n + jx` order.
pub fn render(
    tree: &QuadTree,
    pixels: &[Vec<PixelOut>],
    cam: &Camera,
    res: usize,
    mode: TexelMode,
    rgb: impl Fn(&PixelOut) -> [u8; 3],
) -> (Vec<u8>, Vec<LeafTexel>) {
    let leaves: Vec<usize> = tree.leaves().collect();
    render_leaves(tree, pixels, cam, res, mode, rgb, &leaves)
}

/// Rasterise an **explicit** leaf set — an animation frame, a depth cap, a budget round — with
/// the set's ancestors filled in underneath and nothing outside the set painted at all.
///
/// A leaf with no samples in `pixels` paints nothing and shows its nearest painted ancestor;
/// a leaf whose ancestors have no samples either shows [`crate::output::colour::BACKGROUND`],
/// which is the honest picture of a region nothing has computed.
#[allow(clippy::too_many_arguments)]
pub fn render_leaves(
    tree: &QuadTree,
    pixels: &[Vec<PixelOut>],
    cam: &Camera,
    res: usize,
    mode: TexelMode,
    rgb: impl Fn(&PixelOut) -> [u8; 3],
    leaves: &[usize],
) -> (Vec<u8>, Vec<LeafTexel>) {
    let n = tree.n;
    let mut img: Vec<u8> = crate::output::colour::BACKGROUND
        .iter()
        .cloned()
        .cycle()
        .take(res * res * 3)
        .collect();
    let leaf_set: HashSet<usize> = leaves.iter().cloned().collect();
    let deepest = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
    let px_size = cam.pixel_size();

    let mut info = Vec::with_capacity(leaves.len());
    for i in paint_order(tree, leaves) {
        let Some(px) = pixels.get(i).filter(|p| !p.is_empty()) else {
            continue; // nothing computed here; whatever was painted underneath stays
        };
        let q = &tree.nodes[i];
        let mut drawn = 0usize;
        for (k, p) in px.iter().enumerate() {
            let (jx, jy) = (k % n, k / n);
            let Some((x0, x1, y0, y1)) = tile_px(q, n, cam, res, mode, deepest, jx, jy) else {
                continue;
            };
            let c = rgb(p);
            for y in y0..y1 {
                for x in x0..x1 {
                    let o = (y * res + x) * 3;
                    img[o..o + 3].copy_from_slice(&c);
                }
            }
            drawn += 1;
        }
        // Leaves only: `info` is the texel-scaling instrument, and an ancestor drawn as fill is
        // not a texel anyone is measuring. Including them would silently double the rows and
        // halve the apparent texel size at every level.
        if !leaf_set.contains(&i) {
            continue;
        }
        let dc = match mode {
            TexelMode::Adaptive => 2.0 * q.half / (n - 1) as f64,
            TexelMode::Uniform => {
                2.0 * (q.half * (2f64).powi(q.level as i32 - deepest as i32)) / (n - 1) as f64
            }
        };
        info.push(LeafTexel { node: i, level: q.level, texel_px: dc / px_size, tiles_drawn: drawn });
    }
    (img, info)
}

/// How many leaf tiles cover each pixel — the tiling instrument, geometry only.
///
/// For the leaves of a complete tree over a camera framing the root, every pixel reads exactly
/// `1`: no gap (a `0`), no overlap (a `2`). Samples are not consulted, so a leaf with nothing
/// computed still counts; this measures where the tiles *would* land.
pub fn coverage(tree: &QuadTree, cam: &Camera, res: usize, leaves: &[usize]) -> Vec<u16> {
    let n = tree.n;
    let mut cov = vec![0u16; res * res];
    let deepest = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
    for &i in leaves {
        let q = &tree.nodes[i];
        for jy in 0..n {
            for jx in 0..n {
                let Some((x0, x1, y0, y1)) =
                    tile_px(q, n, cam, res, TexelMode::Adaptive, deepest, jx, jy)
                else {
                    continue;
                };
                for y in y0..y1 {
                    for x in x0..x1 {
                        cov[y * res + x] = cov[y * res + x].saturating_add(1);
                    }
                }
            }
        }
    }
    cov
}

pub fn save(path: &str, res: usize, data: &[u8]) -> std::io::Result<()> {
    save_rect(path, res, res, data)
}

/// As [`save`], for a non-square image. Plots and side-by-side frames are not square.
///
/// Rows are written in buffer order: row 0 is `data[0..w*3]`. For a `Slice`-ordered buffer that
/// is the minimum `y` — the convention every image in this crate now shares.
pub fn save_rect(path: &str, w: usize, h: usize, data: &[u8]) -> std::io::Result<()> {
    assert_eq!(data.len(), w * h * 3, "buffer is not {w}x{h} RGB8");
    if let Some(dir) = Path::new(path).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let file = File::create(Path::new(path))?;
    let mut enc = png::Encoder::new(BufWriter::new(file), w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header()?.write_image_data(data)?;
    Ok(())
}
