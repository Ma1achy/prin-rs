//! Two images, per BRIEF §7: outcome, and ensemble spread.
//!
//! Both are diagnostics. Anything read off them should be confirmed against the raw dump.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use crate::ensemble::pixel::PixelOut;
use crate::grid::Slice;

/// Outcome colouring. The legacy classifier emits body 0/1/2 escaping, or 3 for bound;
/// a pixel with any non-finite copy is flagged separately, because "undetermined" is a
/// distinct answer from "bound" and must not be painted as one.
fn outcome_rgb(p: &PixelOut) -> [u8; 3] {
    if p.n_nonfinite > 0 {
        return [255, 0, 255]; // magenta: undetermined, deliberately loud
    }
    match p.legacy_class {
        0 => [220, 80, 60],
        1 => [70, 150, 220],
        2 => [110, 190, 110],
        _ => [40, 40, 48],
    }
}

/// Perceptually monotone ramp for a value in `[0, 1]`. Not a scientific colourmap; it is
/// only here to make structure visible at a glance.
fn ramp(x: f64) -> [u8; 3] {
    let t = x.clamp(0.0, 1.0);
    let r = (255.0 * t.powf(0.6)) as u8;
    let g = (255.0 * (t * (1.0 - t) * 4.0).powf(0.8)) as u8;
    let b = (255.0 * (1.0 - t).powf(0.6)) as u8;
    [r, g, b]
}

fn save(path: &Path, w: u32, h: u32, data: &[u8]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), w, h);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    writer.write_image_data(data)?;
    Ok(())
}

pub fn write_pair(stem: &str, slice: &Slice, pixels: &[PixelOut]) -> std::io::Result<()> {
    let (w, h) = (slice.nx as u32, slice.ny as u32);

    let mut a = Vec::with_capacity(pixels.len() * 3);
    for p in pixels {
        a.extend_from_slice(&outcome_rgb(p));
    }
    save(Path::new(&format!("{stem}_outcome.png")), w, h, &a)?;

    let mut b = Vec::with_capacity(pixels.len() * 3);
    for p in pixels {
        let v = if p.ensemble_spread.is_finite() { p.ensemble_spread } else { 1.0 };
        b.extend_from_slice(&ramp(v));
    }
    save(Path::new(&format!("{stem}_spread.png")), w, h, &b)?;
    Ok(())
}
