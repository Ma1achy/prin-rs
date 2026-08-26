//! **GIF output, because APNG does not animate where these get looked at.**
//!
//! `apng.rs` writes a structurally valid animation — `acTL`, one `fcTL` per frame, `fdAT` for
//! every frame after the first. It is the better format on the merits: 24-bit colour, no palette,
//! and an APNG *is* a PNG so a viewer that cannot animate shows the first frame rather than
//! refusing the file.
//!
//! **And GitHub's blob viewer does not animate it.** Neither do a lot of image viewers. A
//! diagnostic nobody can see move is not a diagnostic, so these ship as GIF as well.
//!
//! # What the palette costs
//!
//! GIF is 256 colours per frame. The shipping colouring is a continuous OKLCh field, so it is
//! quantised — with **one palette for the whole animation**, computed from a sample across all
//! frames rather than per frame. Per-frame palettes make flat regions shimmer between frames as
//! the quantiser makes different choices about a colour that did not change, which reads as noise
//! in exactly the still areas the eye uses to judge that something else moved.
//!
//! The APNG is kept beside it as the lossless record.

use std::fs::File;
use std::io::BufWriter;

/// Write `frames` as an animated GIF at `delay_cs` hundredths of a second per frame.
///
/// Every frame must be `w * h * 3` bytes, RGB.
pub fn write(
    path: &str,
    w: usize,
    h: usize,
    frames: &[Vec<u8>],
    delay_cs: u16,
) -> std::io::Result<()> {
    assert!(!frames.is_empty(), "an animation needs at least one frame");
    for (i, f) in frames.iter().enumerate() {
        assert_eq!(f.len(), w * h * 3, "frame {i} is the wrong size for {w}x{h}");
    }
    if let Some(dir) = std::path::Path::new(path).parent() {
        std::fs::create_dir_all(dir)?;
    }

    // One palette for the whole animation. Sampled across every frame so a colour that appears
    // only late still gets represented -- taking the palette from frame 0 alone would quantise
    // the finished picture against the coarsest one, which is the frame with the least colour in
    // it and therefore the worst possible reference.
    let stride = (frames.len() * w * h / 40_000).max(1);
    let sample: Vec<u8> = frames
        .iter()
        .flat_map(|f| f.chunks_exact(3).step_by(stride).flat_map(|p| [p[0], p[1], p[2], 255]))
        .collect();
    let nq = color_quant::NeuQuant::new(10, 256, &sample);
    let palette: Vec<u8> = nq.color_map_rgb();

    let file = File::create(path)?;
    let mut enc = gif::Encoder::new(BufWriter::new(file), w as u16, h as u16, &palette)
        .map_err(std::io::Error::other)?;
    enc.set_repeat(gif::Repeat::Infinite).map_err(std::io::Error::other)?;

    for f in frames {
        let idx: Vec<u8> =
            f.chunks_exact(3).map(|p| nq.index_of(&[p[0], p[1], p[2], 255]) as u8).collect();
        let mut frame = gif::Frame::from_indexed_pixels(w as u16, h as u16, idx, None);
        frame.delay = delay_cs;
        enc.write_frame(&frame).map_err(std::io::Error::other)?;
    }
    Ok(())
}
