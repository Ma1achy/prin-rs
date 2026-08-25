//! The viridis colourmap, for **categorical** class images matched to the reference.
//!
//! Kept separate from [`crate::output::plot::palette`], which spreads hue evenly for line-chart
//! series and is not a colourmap: it has no ordering and no perceptual monotonicity, and using it
//! for a class field would imply neither.
//!
//! It exists because the `Ma1achy/principia-ii` WebGPU panel reads
//! `Colour mode: Event class, Palette: viridis`, and this project had no categorical colourmap at
//! all. **Comparing a continuous field against a categorical map is how a rendering choice gets
//! mistaken for a physics bug**, which is most of what went wrong with the preset port: the
//! smooth-rainbow GLSL image being compared against was a continuous field, and the panel's own
//! event-class render is discrete and looks far closer.

/// Matplotlib's viridis at nine evenly spaced control points, interpolated linearly in sRGB.
///
/// Nine rather than 256: the map is smooth enough that the interpolation error is under a code
/// value, and a table small enough to read is a table that can be checked against the source.
const ANCHORS: [[f64; 3]; 9] = [
    [68.0, 1.0, 84.0],
    [72.0, 40.0, 120.0],
    [62.0, 74.0, 137.0],
    [49.0, 104.0, 142.0],
    [38.0, 130.0, 142.0],
    [31.0, 158.0, 137.0],
    [53.0, 183.0, 121.0],
    [109.0, 205.0, 89.0],
    [253.0, 231.0, 37.0],
];

/// Viridis at `t`, clamped to `[0, 1]`. A non-finite `t` returns the low end rather than casting
/// NaN into a `u8` and coming out as black -- but callers should not be handing this a NaN: an
/// undetermined pixel takes [`crate::output::colour::DEBUG_NAN`] and never a colourmap entry.
pub fn viridis(t: f64) -> [u8; 3] {
    let t = if t.is_finite() { t.clamp(0.0, 1.0) } else { 0.0 };
    let n = ANCHORS.len() - 1;
    let x = t * n as f64;
    let i = (x.floor() as usize).min(n - 1);
    let f = x - i as f64;
    let (a, b) = (ANCHORS[i], ANCHORS[i + 1]);
    [
        (a[0] + f * (b[0] - a[0])).round() as u8,
        (a[1] + f * (b[1] - a[1])).round() as u8,
        (a[2] + f * (b[2] - a[2])).round() as u8,
    ]
}
