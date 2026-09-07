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

/// **The same file, written a frame at a time, against a palette chosen in advance.**
///
/// [`write`] holds every frame so it can sample them for the palette and then encode them. At
/// 1024² that is 3 MB a frame, and a thousand-frame animation is 3 GB held for no other reason —
/// the same failure [`crate::output::apng::Stream`] exists for. The split is where it has to be:
/// the palette must be decided **before** the first frame is encoded, so a caller that cannot
/// hold the frames has to sample in one pass and encode in another. [`build_palette`] is that
/// first pass, and it is the same sampling [`write`] does inline.
///
/// `write` is implemented on top of this, so there is one encoder and the two cannot drift.
pub struct Stream {
    enc: gif::Encoder<BufWriter<File>>,
    nq: color_quant::NeuQuant,
    dims: (usize, usize),
    delay_cs: u16,
    written: usize,
}

impl Stream {
    pub fn begin(
        path: &str,
        w: usize,
        h: usize,
        nq: color_quant::NeuQuant,
        delay_cs: u16,
    ) -> std::io::Result<Self> {
        if let Some(dir) = std::path::Path::new(path).parent() {
            std::fs::create_dir_all(dir)?;
        }
        let palette: Vec<u8> = nq.color_map_rgb();
        let file = File::create(path)?;
        let mut enc = gif::Encoder::new(BufWriter::new(file), w as u16, h as u16, &palette)
            .map_err(std::io::Error::other)?;
        enc.set_repeat(gif::Repeat::Infinite).map_err(std::io::Error::other)?;
        Ok(Self { enc, nq, dims: (w, h), delay_cs, written: 0 })
    }

    pub fn push_rgb(&mut self, f: &[u8]) -> std::io::Result<()> {
        let (w, h) = self.dims;
        assert_eq!(f.len(), w * h * 3, "frame {} is the wrong size for {w}x{h}", self.written);
        let idx: Vec<u8> =
            f.chunks_exact(3).map(|p| self.nq.index_of(&[p[0], p[1], p[2], 255]) as u8).collect();
        let mut frame = gif::Frame::from_indexed_pixels(w as u16, h as u16, idx, None);
        frame.delay = self.delay_cs;
        self.enc.write_frame(&frame).map_err(std::io::Error::other)?;
        self.written += 1;
        Ok(())
    }

    /// Write an already-quantised frame. For a caller that quantises in parallel — `index_of` is
    /// a nearest-colour search per pixel, which at 1024² over a thousand frames is a billion of
    /// them, and it is the only part of this that is worth threading.
    pub fn push_indexed(&mut self, idx: Vec<u8>) -> std::io::Result<()> {
        let (w, h) = self.dims;
        assert_eq!(idx.len(), w * h, "indexed frame {} is the wrong size", self.written);
        let mut frame = gif::Frame::from_indexed_pixels(w as u16, h as u16, idx, None);
        frame.delay = self.delay_cs;
        self.enc.write_frame(&frame).map_err(std::io::Error::other)?;
        self.written += 1;
        Ok(())
    }

    pub fn quantiser(&self) -> &color_quant::NeuQuant {
        &self.nq
    }

    pub fn finish(self) -> std::io::Result<usize> {
        assert!(self.written > 0, "an animation needs at least one frame");
        Ok(self.written)
    }
}

/// The palette pass: sample `n_frames` frames' worth of pixels and train the quantiser.
///
/// The sample must span **every** frame, not frame 0: taking the palette from the first frame
/// quantises the finished picture against the emptiest one, which is the worst possible
/// reference. `total_px` is the whole animation's pixel count, and the stride is derived from it
/// so the sample size is the same however the caller feeds the pixels in.
pub fn sample_stride(total_px: usize) -> usize {
    (total_px / 40_000).max(1)
}

pub fn build_palette(sample_rgba: &[u8]) -> color_quant::NeuQuant {
    color_quant::NeuQuant::new(10, 256, sample_rgba)
}

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
    let stride = sample_stride(frames.len() * w * h);
    let sample: Vec<u8> = frames
        .iter()
        .flat_map(|f| f.chunks_exact(3).step_by(stride).flat_map(|p| [p[0], p[1], p[2], 255]))
        .collect();
    let mut s = Stream::begin(path, w, h, build_palette(&sample), delay_cs)?;
    for f in frames {
        s.push_rgb(f)?;
    }
    s.finish()?;
    Ok(())
}
