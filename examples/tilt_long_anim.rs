//! **`tilt_plambda` at the chosen window, 1024², a long horizon and a frame per sync boundary.**
//!
//! Run: `cargo run --release --example tilt_long_anim <out_root> [res] [t_max] [frames] [fps]
//!       [zoom_mult] [u_shift]`
//!
//! The output root is argument one and has **no default** — this writes several gigabytes.
//!
//! # The store is on disk, because at this frame count it cannot be in memory
//!
//! `tilt_slice_render`'s animation holds the projected field as `PixelSlim`, 48 bytes a
//! footprint-frame. That is 1.6 GB at 1024² over 32 frames and **15 GB over 300** — the frame
//! count, not the raster, is what breaks it. So each strip's projections are appended to a single
//! file laid out `[strip][frame][pixel]`, and painting frame `j` reads one contiguous block per
//! strip. The record is **`f64`**, not the `f32` that would halve the file: the colouring is a
//! pure function of `shape_vec` and `spread_shape`, and narrowing them would make this animation
//! differ from every other panel in the corpus for a reason that has nothing to do with physics.
//!
//! Peak resident is one strip of `PixelOut` with the live series attached — about `40 * n_sync`
//! bytes a footprint on top of the 656 — plus **two** painted frames, because the animation is
//! written through `apng::Stream` rather than `apng::write`. A thousand 1024² frames held for the
//! encoder is 3 GB, which is the same failure as the store one level up.
//!
//! # The window is the one picked off the sweep, and it is not the shipped one
//!
//! `zoom_mult` and `u_shift` are applied to `Chart::tilt_plambda`'s window: `half = zoom * mult`
//! and `cx` shifted right by `u_shift` half-windows, `cy` held. Both are printed and land in the
//! sidecar, because a panel at a window that is not the chart's own default is otherwise
//! indistinguishable from one at a window that is.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelSlim};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng};
use prin_rs::scheduler;
use rayon::prelude::*;
use std::io::{Read, Seek, SeekFrom, Write};

const Z10: [f64; 10] = [2.23, -0.56, -0.05, 0.0, 0.05, -0.04, -0.02, 0.12, -0.1, 0.02];
const ZOOM: f64 = 0.167_798_447_231_782_42;
const PAN: (f64, f64) = (0.516_965_993_956_604_7, 0.538_322_670_350_332_3);
const GAMMA: f64 = 2.0;
const REC: usize = 40; // 4 x f64 + 4 x u8, padded

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn put(buf: &mut [u8], p: &PixelSlim) {
    buf[0..8].copy_from_slice(&p.spread_shape.to_le_bytes());
    for k in 0..3 {
        buf[8 + k * 8..16 + k * 8].copy_from_slice(&p.shape_vec[k].to_le_bytes());
    }
    buf[32] = p.n_nonfinite;
    buf[33] = p.state;
    buf[34] = p.detail;
    buf[35] = p.event_class;
}

fn get(buf: &[u8]) -> PixelSlim {
    let f = |o: usize| f64::from_le_bytes(buf[o..o + 8].try_into().unwrap());
    PixelSlim {
        spread_shape: f(0),
        shape_vec: [f(8), f(16), f(24)],
        n_nonfinite: buf[32],
        state: buf[33],
        detail: buf[34],
        event_class: buf[35],
    }
}

fn main() {
    let root: String = std::env::args().nth(1).expect("argument one is the output root, no default");
    let res: usize = arg(2, 1024);
    let t_max: f64 = arg(3, 50.0);
    let frames: usize = arg(4, 300);
    let fps: u16 = arg(5, 30);
    let zmult: f64 = arg(6, 0.55);
    let ushift: f64 = arg(7, 0.8);

    let dir = format!("{root}/tilt_long");
    let _ = std::fs::create_dir_all(&dir);

    // One recorded boundary per frame: `n_sync = frames` at `live_stride = 1`. `n_sync` scales
    // with `t_max` by construction here rather than being inherited — a fixed `n_sync` while
    // `t_max` varies compares different discretisations, which is on this project's record.
    let base = EnsembleCfg::production();
    let ens = EnsembleCfg {
        t_max,
        n_sync: frames,
        keep_live_series: true,
        live_stride: 1,
        ..base
    };
    scheduler::assert_production_kernel(&ens, &dir);

    let (chart, _, _, _) = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], GAMMA, ZOOM, PAN);
    let half = ZOOM * zmult;
    let (cx, cy) = (2.0 * PAN.0 - 1.0 + ushift * half, 2.0 * PAN.1 - 1.0);
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let npix = sl.npix();

    println!("tilt_plambda long animation: {res}^2, t_max {t_max}, {frames} frames at {fps} fps");
    println!("config: {}", ens.provenance());
    println!("window: cx {cx:.9} cy {cy:.9} half {half:.9}  (zoom x{zmult}, u shift +{ushift} \
              half-windows; the SHIPPED window is cx {:.6} cy {:.6} half {ZOOM:.6})",
             2.0 * PAN.0 - 1.0, 2.0 * PAN.1 - 1.0);
    println!("reproduce: cargo run --release --example tilt_long_anim {root} {res} {t_max} \
              {frames} {fps} {zmult} {ushift}");

    let strip = 8usize;
    let nstrip = sl.ny.div_ceil(strip);
    let store_path = format!("{dir}/store.bin");
    let bytes = (frames * npix * REC) as f64 / 1e9;
    println!("store: {store_path}, {bytes:.2} GB, laid out [strip][frame][pixel]");
    let mut store = std::fs::File::create(&store_path).expect("cannot create the store");

    let t0 = std::time::Instant::now();
    let mut nb = frames;
    for s in 0..nstrip {
        let (r0, r1) = (s * strip, ((s + 1) * strip).min(sl.ny));
        let band: Vec<pixel::PixelOut> = (r0 * sl.nx..r1 * sl.nx)
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &ens))
            .collect();
        let here = band.iter().map(|p| p.live_t.len()).min().unwrap_or(0);
        assert!(here > 1, "strip {s} recorded {here} boundaries: nothing to animate");
        nb = nb.min(here);
        // Frame-major within the strip, so painting reads one contiguous run per strip.
        let rows = r1 - r0;
        let mut blk = vec![0u8; frames * rows * sl.nx * REC];
        for j in 0..frames.min(here) {
            for (i, p) in band.iter().enumerate() {
                let o = (j * rows * sl.nx + i) * REC;
                put(&mut blk[o..o + REC], &PixelSlim::thin(&scheduler::project_at(p, j)));
            }
        }
        store.write_all(&blk).expect("store write failed");
        if s % 16 == 0 || s + 1 == nstrip {
            println!("  strip {}/{nstrip}  rows {r0}..{r1}  {:.1}s elapsed",
                     s + 1, t0.elapsed().as_secs_f64());
        }
    }
    store.flush().unwrap();
    drop(store);
    println!("{npix} footprints in {:.1}s, {nb} boundaries usable", t0.elapsed().as_secs_f64());

    // **One window over every frame**, read back from the store. A per-frame range would rescale
    // the ramp as the field develops and turn a growing signal into a constant one.
    let mut f = std::fs::File::open(&store_path).unwrap();
    let mut rec = vec![0u8; REC];
    let mut spreads: Vec<f64> = Vec::with_capacity(nb * npix / 16);
    // Every 16th footprint of every frame: 6.5 M values at these settings, which pins a p1/p99
    // far past the precision the ramp has, and reads a sixteenth of the file.
    for s in 0..nstrip {
        let rows = ((s + 1) * strip).min(sl.ny) - s * strip;
        let base_off = (0..s).map(|q| {
            let r = ((q + 1) * strip).min(sl.ny) - q * strip;
            frames * r * sl.nx * REC
        }).sum::<usize>();
        for j in 0..nb {
            for i in (0..rows * sl.nx).step_by(16) {
                let o = base_off + (j * rows * sl.nx + i) * REC;
                f.seek(SeekFrom::Start(o as u64)).unwrap();
                f.read_exact(&mut rec).unwrap();
                spreads.push(get(&rec).spread_shape);
            }
        }
    }
    let (alo, ahi) = colour::range_q_of(spreads.iter().copied(), 0.01, 0.99);
    println!("fixed ramp window over all {nb} frames ({} sampled values): ({alo:.4e}, {ahi:.4e})",
             spreads.len());
    drop(spreads);

    let sites = colour::landmarks(&grid::decode_state(&chart, 0, cx, cy).m);
    let path = format!("{dir}/tilt_plambda_long_{res}_t{t_max}_{nb}f.png");
    let mut anim = apng::Stream::begin(&path, res, res, nb, 1, fps).expect("cannot open the APNG");
    let mut last_frame: Vec<u8> = Vec::new();
    let mut vetoed_any = 0usize;
    for j in 0..nb {
        let mut buf = vec![0u8; npix * 3];
        for s in 0..nstrip {
            let (r0, r1) = (s * strip, ((s + 1) * strip).min(sl.ny));
            let rows = r1 - r0;
            let base_off = (0..s).map(|q| {
                let r = ((q + 1) * strip).min(sl.ny) - q * strip;
                frames * r * sl.nx * REC
            }).sum::<usize>();
            let mut blk = vec![0u8; rows * sl.nx * REC];
            f.seek(SeekFrom::Start((base_off + j * rows * sl.nx * REC) as u64)).unwrap();
            f.read_exact(&mut blk).unwrap();
            for i in 0..rows * sl.nx {
                let p = get(&blk[i * REC..(i + 1) * REC]).fat();
                if colour::vetoed(&p, Scalar::ShapeSpread, &sites, alo, ahi) {
                    vetoed_any += 1;
                }
                let px = colour::rgb_veto(&p, Scalar::ShapeSpread, &sites, alo, ahi, colour::Veto::None);
                let o = (r0 * sl.nx + i) * 3;
                buf[o..o + 3].copy_from_slice(&px);
            }
        }
        anim.push(&buf).expect("frame write failed");
        if j + 1 == nb {
            last_frame = buf;
        }
        if j % 25 == 0 || j + 1 == nb {
            println!("  painted frame {}/{nb}, {} duplicates so far, {:.1}s elapsed",
                     j + 1, anim.adjacent_duplicates(), t0.elapsed().as_secs_f64());
        }
    }

    let dup = anim.finish().expect("APNG finish failed");
    println!("{nb} frames, {dup} adjacent duplicates, {vetoed_any} flagged footprint-frames");
    // The duplicate check runs AFTER the file is closed rather than instead of writing it: a
    // thousand-frame run that trips this is still worth looking at, and the assertion is what
    // says so out loud.
    assert!(dup + 1 < nb, "every frame duplicates its neighbour: the playhead never moved");
    let _ = prin_rs::output::provenance_sidecar(
        &path,
        &ens,
        &format!(
            "chart=latent_tilt_plambda panel=uniform_over_time scalar=ShapeSpread \
             window=({alo:.4e},{ahi:.4e}) window_from=all_frames res={res} frames={nb} \
             t_max={t_max} fps={fps} cx={cx} cy={cy} half={half} zoom_mult={zmult} \
             u_shift={ushift} shipped_window=false adjacent_duplicates={dup} tree=none \
             veto=none\n"
        ),
    );
    println!("wrote {path}");

    // A still at the same window and the same ramp, so the animation has a frame that can be
    // looked at without a player.
    let last = format!("{dir}/tilt_plambda_long_{res}_final.png");
    let _ = adaptive::save_rect(&last, res, res, &last_frame);
    println!("wrote {last}");

    let _ = std::fs::remove_file(&store_path);
    println!("store removed ({bytes:.2} GB reclaimed)");
}
