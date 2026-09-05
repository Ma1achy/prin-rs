//! **Rebuild the gallery flip-books from the committed stills.**
//!
//! `chart_gallery` writes `animated/gallery.png` and `gallery_wire.png` from the charts of *that
//! run*, so a gallery rendered in batches -- which is what a memory-constrained machine forces --
//! leaves the two APNGs holding only the last batch. That is not a small cosmetic loss: a
//! flip-book claiming to be the gallery while showing four of twenty-six charts is an artefact
//! that lies about its own scope, and nothing in the file says so.
//!
//! This reads the per-chart stills back and re-encodes both APNGs over the full case list, in
//! `grid::gallery_cases()` order, so the result is independent of how the run was batched.
//!
//! # The guards, and what each one would catch
//!
//! - **Frame count equals the case count.** A missing chart is a silent short flip-book
//!   otherwise; this refuses rather than writing one.
//! - **Every frame the same dimensions.** Mixing rasters would encode, and the record already
//!   carries *softness in an image is a raster size* four times over.
//! - **Not every adjacent pair identical.** If the decode silently produced blank or repeated
//!   buffers the APNG would still write and still look like a gallery. A flip-book whose frames
//!   all match is a still, and that is the arm that says the reads worked at all.
//!
//! Run: `cargo run --release --example gallery_apng -- [root=results]`

use prin_rs::grid;
use prin_rs::output::apng;

fn read_rgb(path: &str) -> Option<(usize, usize, Vec<u8>)> {
    let f = std::fs::File::open(path).ok()?;
    let dec = png::Decoder::new(std::io::BufReader::new(f));
    let mut rd = dec.read_info().ok()?;
    let mut buf = vec![0u8; rd.output_buffer_size()];
    let info = rd.next_frame(&mut buf).ok()?;
    let (w, h) = (info.width as usize, info.height as usize);
    // `adaptive::save` writes RGB8; RGBA is accepted and its alpha dropped rather than refused,
    // because a frame that decodes fine should not fail on a channel the flip-book never uses.
    let rgb = match info.color_type {
        png::ColorType::Rgb => buf[..w * h * 3].to_vec(),
        png::ColorType::Rgba => {
            buf[..w * h * 4].chunks_exact(4).flat_map(|p| [p[0], p[1], p[2]]).collect()
        }
        other => {
            eprintln!("  {path}: unsupported colour type {other:?}");
            return None;
        }
    };
    Some((w, h, rgb))
}

/// Collect one flip-book's frames. **Collect only -- nothing is written here**, because the two
/// APNGs are a matched pair: writing the wire book while the colour book was refused leaves a
/// complete flip-book beside a stale one and no file saying which is which. The first cut did
/// exactly that and its own missing-chart test caught it.
fn collect(root: &str, suffix: &str) -> Option<(usize, usize, Vec<Vec<u8>>)> {
    let cases = grid::gallery_cases();
    let mut frames: Vec<Vec<u8>> = Vec::with_capacity(cases.len());
    let mut dims: Option<(usize, usize)> = None;
    for (name, ..) in &cases {
        let p = format!("{root}/charts/{name}{suffix}.png");
        let Some((w, h, rgb)) = read_rgb(&p) else {
            eprintln!("  MISSING or unreadable: {p}");
            return None;
        };
        match dims {
            None => dims = Some((w, h)),
            Some((w0, h0)) if (w0, h0) != (w, h) => {
                eprintln!("  {p}: {w}x{h} against {w0}x{h0} -- mixed rasters, refusing");
                return None;
            }
            _ => {}
        }
        frames.push(rgb);
    }
    assert_eq!(frames.len(), cases.len(), "one frame per gallery case");
    let (w, h) = dims?;
    Some((w, h, frames))
}

fn write_book(out: &str, w: usize, h: usize, frames: &[Vec<u8>]) {
    let dup = apng::adjacent_duplicates(frames);
    // Every adjacent pair identical means the decode produced the same buffer every time --
    // which writes a perfectly valid APNG that is a still wearing a flip-book's name.
    assert!(
        dup < frames.len().saturating_sub(1),
        "all {} frames identical: the reads did not work",
        frames.len()
    );
    let _ = apng::write(out, w, h, frames, 1, 2);
    // One pair is expected and explained: `body_plane` and `plane_00deg` are the same chart under
    // two names, which the gallery's own control asserts at `max |dIC| = 0e0`.
    println!("  {out}: {} frames at {w}x{h}, {dup} identical adjacent pairs", frames.len());
}

fn main() {
    let root = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    println!("rebuilding the gallery flip-books from {root}/charts over all {} cases",
             grid::gallery_cases().len());
    let (a, b) = (collect(&root, ""), collect(&root, "_wire"));
    let (Some((aw, ah, af)), Some((bw, bh, bf))) = (a, b) else {
        eprintln!(
            "\nNEITHER written: a chart is missing or the rasters disagree. Run the gallery for \
             the missing charts first -- a short flip-book is worse than none, and a complete \
             one beside a stale one is worse than both."
        );
        std::process::exit(1);
    };
    write_book(&format!("{root}/animated/gallery.png"), aw, ah, &af);
    write_book(&format!("{root}/animated/gallery_wire.png"), bw, bh, &bf);
}
