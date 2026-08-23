//! **The screen floor is view-relative, and a still image cannot show it.**
//!
//! The contract: *"Do not cache the screen floor as a quad fact. It is view-relative and
//! evaluated live. A quad floored at one zoom must refine again when zoomed into."* In a single
//! frame that is indistinguishable from a quad that was simply never interesting. Across a zoom
//! ladder it is the whole behaviour: each frame descends six levels below its own camera depth,
//! the patch that was one flat texel becomes 4096 of them, and the samples are **new**, not
//! upsampled.
//!
//! Writes one PNG per frame, plus an animated APNG of the whole ladder (named `.png`, since an
//! APNG is a PNG and a viewer that does not animate shows the first frame), plus the raw tree dump
//! for every frame so the pictures can be checked against numbers.
//!
//! Run: `cargo run --release --example zoom_sequence [region] [frames] [budget]`

use std::fs::File;
use std::io::BufWriter;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::adaptive::{self, TexelMode};
use prin_rs::output::tree as treedump;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn spread_rgb(p: &PixelOut) -> [u8; 3] {
    let (lo, hi) = (1e-8f64, 1e-1f64);
    let v = p.ensemble_spread;
    let t = if !v.is_finite() {
        1.0
    } else if v <= lo {
        0.0
    } else {
        ((v.ln() - lo.ln()) / (hi.ln() - lo.ln())).clamp(0.0, 1.0)
    };
    [
        (255.0 * t.powf(0.6)) as u8,
        (255.0 * (t * (1.0 - t) * 4.0).powf(0.8)) as u8,
        (255.0 * (1.0 - t).powf(0.6)) as u8,
    ]
}

fn main() -> std::io::Result<()> {
    let region: String = arg(1, "near-field".to_string());
    let frames: u32 = arg(2, 9);
    let budget: usize = arg(3, 2000);
    const RES: usize = 384;

    let root = grid::region(&region, 2, 2, 0.05).expect("unknown region");
    let ens = EnsembleCfg { refine_flagged: false, keep_copy_outcomes: true, ..Default::default() };
    let stem = format!("results/vertical/zoom_{}", region.replace(' ', "_"));

    println!("zoom ladder over {region}: {frames} frames, {RES}x{RES}, budget {budget} quads each.");
    println!("Each frame re-descends from a root box of half = 0.05 / 2^depth with the camera");
    println!("framing it, so `camera_depth` is 0 in every frame and the screen floor always sits");
    println!("six levels below the view. The samples in frame k+1 are NEW, not upsampled.\n");
    println!("{:>6} {:>12} {:>7} {:>7} {:>6} {:>6} {:>6} {:>11} {:>9}",
             "frame", "half", "quads", "leaves", "depth", "screen", "keep", "spread med", "wall_s");

    let mut buffers: Vec<Vec<u8>> = Vec::new();
    for k in 0..frames {
        let half = 0.05 / (2f64).powi(k as i32);
        let cam = Camera::framing(root.cx, root.cy, half, RES);
        let cfg = SchedCfg {
            budget, tau_display: 1e-4, alpha_hi: 0.2, alpha_lo: 0.08,
            camera: Some(cam), keep_pixels: true, ..Default::default()
        };
        let (t, st) = scheduler::descend(
            root.cx, root.cy, half, root.body, &cfg, &ens, Precision::F64);
        let leaves: Vec<usize> = t.leaves().collect();
        let mut sp: Vec<f64> = leaves.iter().map(|&i| t.nodes[i].red.spread_median)
            .filter(|x| x.is_finite()).collect();
        sp.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("{:>6} {half:>12.4e} {:>7} {:>7} {:>6} {:>6} {:>6} {:>11.3e} {:>9.1}",
                 k, st.quads_computed, leaves.len(),
                 t.depth_histogram().len().saturating_sub(1),
                 leaves.iter().filter(|&&i| t.nodes[i].decision == D::ScreenFloor).count(),
                 leaves.iter().filter(|&&i| t.nodes[i].decision == D::Keep).count(),
                 sp.get(sp.len() / 2).cloned().unwrap_or(f64::NAN), st.wall_seconds);

        let (img, _) = adaptive::render(&t, &st.pixels, &cam, RES, TexelMode::Adaptive, spread_rgb);
        adaptive::save(&format!("{stem}_{k:02}.png"), RES, &img)?;
        buffers.push(img);

        // The raw dump beside every picture, so nothing here has to be taken on trust.
        let mut w = BufWriter::new(File::create(format!("{stem}_{k:02}.prnq"))?);
        treedump::write(&mut w, &t, &cfg, &ens, &st, &region, "f64")?;
    }

    // APNG of the ladder. No new dependency: png 0.17 writes animated PNG directly, and every
    // frame is also on disk as a still, so nothing depends on APNG support to be readable.
    // Named .png, not .apng: an APNG *is* a PNG, and viewers that do not animate show
    // the first frame rather than refusing the file. The .apng extension mostly stops them trying.
    let file = File::create(format!("{stem}_animated.png"))?;
    let mut enc = png::Encoder::new(BufWriter::new(file), RES as u32, RES as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.set_animated(buffers.len() as u32, 0)?;
    enc.set_frame_delay(1, 2)?;
    let mut writer = enc.write_header()?;
    for b in &buffers {
        writer.write_image_data(b)?;
    }
    writer.finish()?;

    println!("\nwrote {frames} frames, {stem}_animated.png (APNG), and one .prnq per frame.");
    println!("The `screen` column is the point: at every zoom level a fresh population of leaves");
    println!("is stopped by the view, and the previous frame's floored quads have refined again.");
    Ok(())
}
