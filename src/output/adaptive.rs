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
//! **A known geometric cost, stated rather than hidden.** `Slice::axis` is endpoint-inclusive,
//! so a quad's `N` samples run corner to corner and a sample-centred tile overhangs the quad
//! box by half a cell on each side. Leaves are therefore painted **coarsest first**, so a finer
//! neighbour overwrites the overhang. This is the same endpoint-inclusive duplication already
//! recorded at sibling edges (1/N of a quad, 12.5% at N=8), seen from the render side.

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

/// Rasterise the leaves into an RGB buffer of `res × res`.
///
/// `pixels[node]` holds that quad's `N²` footprints in `idx = jy*n + jx` order.
#[allow(clippy::too_many_arguments)]
pub fn render(
    tree: &QuadTree,
    pixels: &[Vec<PixelOut>],
    cam: &Camera,
    res: usize,
    mode: TexelMode,
    rgb: impl Fn(&PixelOut) -> [u8; 3],
) -> (Vec<u8>, Vec<LeafTexel>) {
    let n = tree.n;
    let mut img = vec![18u8; res * res * 3];

    // Coarsest first, so a finer neighbour overwrites the half-cell overhang.
    let mut leaves: Vec<usize> = tree.leaves().filter(|&i| pixels.get(i).is_some_and(|p| !p.is_empty())).collect();
    leaves.sort_by_key(|&i| tree.nodes[i].level);

    let deepest = leaves.iter().map(|&i| tree.nodes[i].level).max().unwrap_or(0);
    let px_size = cam.pixel_size();
    // World -> screen, with the camera centre at the image centre.
    let to_px = |x: f64, y: f64| -> (f64, f64) {
        (
            (x - cam.cx) / px_size + res as f64 / 2.0,
            res as f64 / 2.0 - (y - cam.cy) / px_size,
        )
    };

    let mut info = Vec::with_capacity(leaves.len());
    for &i in &leaves {
        let q = &tree.nodes[i];
        let cell = 2.0 * q.half / (n - 1) as f64;
        let draw_cell = match mode {
            TexelMode::Adaptive => cell,
            // Every leaf at the finest leaf's texel size — the wrong instrument, kept
            // reproducible so the acceptance test has a negative case.
            TexelMode::Uniform => 2.0 * (q.half * (2f64).powi(q.level as i32 - deepest as i32))
                / (n - 1) as f64,
        };
        let mut drawn = 0usize;
        for (k, p) in pixels[i].iter().enumerate() {
            let (jx, jy) = (k % n, k / n);
            let sx = q.cx - q.half + jx as f64 * cell;
            let sy = q.cy - q.half + jy as f64 * cell;
            let (ax, ay) = to_px(sx - draw_cell / 2.0, sy + draw_cell / 2.0);
            let (bx, by) = to_px(sx + draw_cell / 2.0, sy - draw_cell / 2.0);
            let c = rgb(p);
            let (x0, x1) = (ax.floor().max(0.0) as usize, bx.ceil().min(res as f64) as usize);
            let (y0, y1) = (ay.floor().max(0.0) as usize, by.ceil().min(res as f64) as usize);
            if x1 <= x0 || y1 <= y0 {
                continue;
            }
            for y in y0..y1 {
                for x in x0..x1 {
                    let o = (y * res + x) * 3;
                    img[o..o + 3].copy_from_slice(&c);
                }
            }
            drawn += 1;
        }
        info.push(LeafTexel {
            node: i,
            level: q.level,
            texel_px: draw_cell / px_size,
            tiles_drawn: drawn,
        });
    }
    (img, info)
}

pub fn save(path: &str, res: usize, data: &[u8]) -> std::io::Result<()> {
    let file = File::create(Path::new(path))?;
    let mut enc = png::Encoder::new(BufWriter::new(file), res as u32, res as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header()?.write_image_data(data)?;
    Ok(())
}
