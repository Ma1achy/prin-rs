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

/// **The same file, written a frame at a time.**
///
/// [`write`] takes every frame at once, which at 1024² is 3 MB each: a thousand-frame animation
/// is 3 GB of `Vec<u8>` held only so the encoder can be handed a slice. This owns the encoder
/// instead and takes frames as they are painted, so the caller holds one frame.
///
/// **The saving is in the caller, not in the encoder.** The first cut reached for
/// `into_stream_writer`, which is for a *single* image: it concatenated every frame into one and
/// the file decoded to **zero** frames while still carrying plausible `fcTL`/`fdAT` chunks and
/// being three times smaller — a corruption that reads as a compression win. `write_image_data`
/// already takes one frame at a time; what had to change was the caller holding a `Vec` of them.
///
/// **It keeps the previous frame, and that is not an optimisation to remove.** The
/// adjacent-duplicate count is the arm that says the playhead moved — *every animation this
/// project produced before the truncation fix was one image repeated N times* — and a streaming
/// writer that dropped the previous frame could not compute it. So the check survives the
/// restructure at a cost of exactly one frame.
///
/// The frame count is declared to the encoder up front and asserted at [`Stream::finish`]: an
/// APNG whose `acTL` count disagrees with the frames written is a corrupt file that most viewers
/// show as a still, which is the same failure the duplicate count exists to catch wearing a
/// different cause.
pub struct Stream {
    w: png::Writer<BufWriter<File>>,
    frame_len: usize,
    prev: Option<Vec<u8>>,
    dup: usize,
    written: usize,
    declared: usize,
}

impl Stream {
    pub fn begin(
        path: &str,
        w: usize,
        h: usize,
        n_frames: usize,
        delay_num: u16,
        delay_den: u16,
    ) -> std::io::Result<Self> {
        assert!(n_frames > 0, "an animation needs at least one frame");
        if let Some(dir) = std::path::Path::new(path).parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let file = File::create(path)?;
        let mut enc = png::Encoder::new(BufWriter::new(file), w as u32, h as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        enc.set_animated(n_frames as u32, 0)
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        enc.set_frame_delay(delay_num, delay_den)
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        Ok(Self {
            w: enc.write_header()?,
            frame_len: w * h * 3,
            prev: None,
            dup: 0,
            written: 0,
            declared: n_frames,
        })
    }

    pub fn push(&mut self, f: &[u8]) -> std::io::Result<()> {
        assert_eq!(f.len(), self.frame_len, "frame {} is the wrong size", self.written);
        assert!(self.written < self.declared, "more frames than were declared");
        if self.prev.as_deref() == Some(f) {
            self.dup += 1;
        }
        self.w.write_image_data(f)?;
        match &mut self.prev {
            Some(p) => p.copy_from_slice(f),
            None => self.prev = Some(f.to_vec()),
        }
        self.written += 1;
        Ok(())
    }

    /// How many adjacent frame pairs written so far were byte-identical.
    pub fn adjacent_duplicates(&self) -> usize {
        self.dup
    }

    pub fn finish(self) -> std::io::Result<usize> {
        assert_eq!(
            self.written, self.declared,
            "declared {} frames and wrote {}: the acTL count would be wrong",
            self.declared, self.written
        );
        self.w.finish()?;
        Ok(self.dup)
    }
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

/// How many adjacent frame pairs are byte-identical.
///
/// **Print it before writing any animation.** Every animation this project produced before
/// the truncation fix was one image repeated N times, and the check that would have caught it
/// is this one line. A deliberate hold on the final frame is the only legitimate source of
/// duplicates; anything past the hold is a still wearing an animation's name.
pub fn adjacent_duplicates(frames: &[Vec<u8>]) -> usize {
    frames.windows(2).filter(|w| w[0] == w[1]).count()
}
