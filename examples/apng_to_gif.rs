//! **An APNG converted to GIF, a frame at a time, with one palette for the whole animation.**
//!
//! Run: `cargo run --release --example apng_to_gif <in.png> <out.gif> [delay_cs] [div] [stride] [cap]`
//!
//! The APNG is the lossless record and is never touched. GIF exists because APNG does not animate
//! in a lot of viewers, and a diagnostic nobody can see move is not a diagnostic.
//!
//! # Two passes, because the palette has to be decided before the first frame is encoded
//!
//! Pass one decodes every frame and samples pixels across **all** of them; pass two decodes again
//! and encodes. Taking the palette from frame 0 would quantise the finished picture against the
//! emptiest frame in the animation. Neither pass holds more than one frame — the input here is
//! 2.57 GB over 1040 frames at 1024², and holding them is the failure this streams to avoid.
//!
//! # What the conversion costs, stated rather than discovered
//!
//! GIF is **256 colours** against the APNG's 24-bit, and this field is a continuous OKLCh ramp,
//! so banding in a smooth region is the format and not the physics. GIF delays are in hundredths
//! of a second, so 30 fps is not representable: `delay_cs = 3` is 33.3 fps and `4` is 25.
//! `div` box-averages the raster down (never nearest — a downscale that drops pixels aliases the
//! fine structure this field is full of), and `stride` keeps every `stride`th frame.
//!
//! `cap` stops after that many kept frames, which truncates the **horizon**: `stride` lowers the
//! frame rate over the whole animation, `cap` shortens it. They are separate because a GIF small
//! enough to sit in a README can be bought either way and they cost different things -- a coarser
//! stride aliases the motion, a shorter cap simply shows less of it.

use prin_rs::output::gifout;
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// Box-average by an integer divisor. `div == 1` is a move, not a copy of the arithmetic.
fn shrink(src: &[u8], w: usize, h: usize, div: usize) -> (Vec<u8>, usize, usize) {
    if div <= 1 {
        return (src.to_vec(), w, h);
    }
    let (nw, nh) = (w / div, h / div);
    let mut out = vec![0u8; nw * nh * 3];
    for y in 0..nh {
        for x in 0..nw {
            let mut acc = [0u32; 3];
            for dy in 0..div {
                for dx in 0..div {
                    let o = ((y * div + dy) * w + x * div + dx) * 3;
                    for c in 0..3 {
                        acc[c] += src[o + c] as u32;
                    }
                }
            }
            let n = (div * div) as u32;
            for c in 0..3 {
                out[(y * nw + x) * 3 + c] = (acc[c] / n) as u8;
            }
        }
    }
    (out, nw, nh)
}

fn each_frame(path: &str, mut f: impl FnMut(usize, &[u8], usize, usize)) -> usize {
    let dec = png::Decoder::new(std::fs::File::open(path).expect("cannot open the APNG"));
    let mut rdr = dec.read_info().expect("not a readable PNG");
    let mut buf = vec![0u8; rdr.output_buffer_size()];
    let mut n = 0usize;
    while let Ok(info) = rdr.next_frame(&mut buf) {
        let (w, h) = (info.width as usize, info.height as usize);
        assert_eq!(info.color_type, png::ColorType::Rgb, "expected an RGB8 APNG");
        f(n, &buf[..info.buffer_size()], w, h);
        n += 1;
    }
    n
}

fn main() {
    let src: String = std::env::args().nth(1).expect("argument one is the input .png");
    let dst: String = std::env::args().nth(2).expect("argument two is the output .gif");
    let delay_cs: u16 = arg(3, 3);
    let div: usize = arg(4, 1);
    let stride: usize = arg(5, 1);
    let cap: usize = arg(6, usize::MAX);

    println!("apng -> gif: {src}");
    println!("  delay {delay_cs} cs ({:.1} fps), raster /{div}, every {stride} frame(s), cap {}",
             100.0 / delay_cs as f64,
             if cap == usize::MAX { "none".to_string() } else { cap.to_string() });

    // Pass one: the palette, sampled across every frame that will be written.
    let t0 = std::time::Instant::now();
    let mut sample: Vec<u8> = Vec::new();
    let mut dims = (0usize, 0usize);
    let mut kept = 0usize;
    // The stride is derived from the whole animation's pixel count, so the sample size does not
    // depend on how the pixels are fed in. It needs the frame count, which needs a first decode;
    // rather than decode three times, sample generously here and let `build_palette` see more.
    let total = each_frame(&src, |i, f, w, h| {
        if i % stride != 0 || kept >= cap {
            return;
        }
        let (small, sw, sh) = shrink(f, w, h, div);
        dims = (sw, sh);
        kept += 1;
        let st = gifout::sample_stride(sw * sh * 64);
        sample.extend(small.chunks_exact(3).step_by(st).flat_map(|p| [p[0], p[1], p[2], 255]));
    });
    let (w, h) = dims;
    println!("  {total} frames in, {kept} kept, {w}x{h}, palette sample {} px, {:.1}s",
             sample.len() / 4, t0.elapsed().as_secs_f64());
    assert!(kept > 1, "only {kept} frames survive the stride: nothing to animate");

    let nq = gifout::build_palette(&sample);
    drop(sample);
    println!("  palette trained, {:.1}s", t0.elapsed().as_secs_f64());

    // Pass two: encode. `index_of` is a nearest-colour search per pixel — a billion of them at
    // these settings — so it is threaded per row and the encoder itself stays sequential.
    let mut out = gifout::Stream::begin(&dst, w, h, nq, delay_cs).expect("cannot open the GIF");
    let mut written = 0usize;
    each_frame(&src, |i, f, fw, fh| {
        if i % stride != 0 || written >= kept {
            return;
        }
        let (small, sw, _) = shrink(f, fw, fh, div);
        let q = out.quantiser();
        let idx: Vec<u8> = small
            .par_chunks_exact(sw * 3)
            .flat_map_iter(|row| {
                row.chunks_exact(3).map(|p| q.index_of(&[p[0], p[1], p[2], 255]) as u8)
            })
            .collect();
        out.push_indexed(idx).expect("gif frame write failed");
        written += 1;
        if written % 50 == 0 || written == kept {
            println!("    {written}/{kept} frames, {:.1}s", t0.elapsed().as_secs_f64());
        }
    });
    let n = out.finish().expect("gif finish failed");
    assert_eq!(n, kept, "wrote {n} frames against {kept} sampled: the two passes disagree");

    let bytes = std::fs::metadata(&dst).map(|m| m.len()).unwrap_or(0);
    println!("wrote {dst}, {:.2} GB, {n} frames at {delay_cs} cs", bytes as f64 / 1e9);
}
