//! Animated PNG, written natively by the `png` crate.
//!
//! Promoted out of `examples/zoom_sequence.rs`, where it was inline. There is no GIF dependency
//! in this project and adding one to animate a few frames was not worth it — an APNG **is** a
//! PNG, so viewers that do not animate show the first frame rather than nothing, and every
//! frame is written as an ordinary PNG beside it so nothing here depends on APNG support to be
//! readable.
//!
//! The file is named `.png` deliberately, for the same reason.

use std::fs::File;
use std::io::BufWriter;

/// Write `frames` (each `w*h*3` bytes, RGB8) as an animated PNG.
///
/// `delay_num/delay_den` is the per-frame delay in seconds. All frames must be the same size —
/// asserted, because a mismatched frame is a silent corruption in the APNG chunk stream rather
/// than an error at write time.
pub fn write(
    path: &str,
    w: usize,
    h: usize,
    frames: &[Vec<u8>],
    delay_num: u16,
    delay_den: u16,
) -> std::io::Result<()> {
    assert!(!frames.is_empty(), "an animation needs at least one frame");
    for (i, f) in frames.iter().enumerate() {
        assert_eq!(f.len(), w * h * 3, "frame {i} is the wrong size for {w}x{h}");
    }
    if let Some(dir) = std::path::Path::new(path).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.set_animated(frames.len() as u32, 0)
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    enc.set_frame_delay(delay_num, delay_den)
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    let mut w2 = enc.write_header()?;
    for f in frames {
        w2.write_image_data(f)?;
    }
    w2.finish()?;
    Ok(())
}

/// Lay two same-sized RGB8 images side by side into one frame.
///
/// Used for before/after comparisons where the interesting thing is *which* quads each side
/// spent its budget on. Two separate animations would make the reader hold one in memory while
/// watching the other, which is precisely the comparison the picture exists to remove.
pub fn side_by_side(a: &[u8], b: &[u8], w: usize, h: usize) -> Vec<u8> {
    let mut out = vec![0u8; w * 2 * h * 3];
    for y in 0..h {
        let src = y * w * 3;
        let dst = y * w * 2 * 3;
        out[dst..dst + w * 3].copy_from_slice(&a[src..src + w * 3]);
        out[dst + w * 3..dst + w * 6].copy_from_slice(&b[src..src + w * 3]);
    }
    out
}

/// Draw a 1px vertical divider down the middle of a side-by-side frame, so the seam is not
/// mistaken for structure.
pub fn divide(frame: &mut [u8], w: usize, h: usize, rgb: [u8; 3]) {
    for y in 0..h {
        let o = (y * w * 2 + w) * 3;
        frame[o] = rgb[0];
        frame[o + 1] = rgb[1];
        frame[o + 2] = rgb[2];
    }
}
