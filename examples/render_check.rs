//! **One region through the adaptive render and the wire, for looking at.**
//!
//! A small descent at a small viewport, rendered three ways into a root you name: the adaptive
//! panel, its wire twin, and the leaf-coverage map (`adaptive::coverage`, painted so that a gap
//! is black and an overlap is red). Beside them, the same window on a uniform grid at the same
//! raster, so the pair can be read as a pair — same orientation, same colour window.
//!
//! **The root is required and must not be `results/`.** This is a validation run; a small raster
//! under `results/` reads as a rendering fault rather than a stale file, which has cost a round
//! trip twice on this project.
//!
//! Run: `cargo run --release --example render_check -- <root> [region] [res] [budget]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, wire};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let Some(root) = std::env::args().nth(1) else {
        eprintln!("usage: render_check <root> [region] [res] [budget]  -- root must not be results/");
        std::process::exit(2);
    };
    assert!(
        !root.replace('\\', "/").split('/').any(|c| c == "results"),
        "render_check is a validation run and must not write under results/"
    );
    let region: String = arg(2, "near-field".to_string());
    let res: usize = arg(3, 256);
    let budget: usize = arg(4, 2000);
    let _ = std::fs::create_dir_all(&root);

    let r = grid::region(&region, 2, 2, 0.05).expect("unknown region");
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let cam = Camera::framing(r.cx, r.cy, 0.05, res);
    let cfg = SchedCfg { budget, camera: Some(cam), keep_pixels: true, ..Default::default() };
    println!("config: {}   sched: tau_display={:e} alpha_hi={} alpha_lo={} policy={:?} k_frac={}",
             ens.provenance(), cfg.tau_display, cfg.alpha_hi, cfg.alpha_lo, cfg.policy, cfg.k_frac);
    let (t, st) = scheduler::descend(r.cx, r.cy, 0.05, r.body, &cfg, &ens, Precision::F64);
    let leaves: Vec<usize> = t.leaves().collect();
    let depth = leaves.iter().map(|&i| t.nodes[i].level).max().unwrap_or(0);
    println!("{region}: {} quads, {} leaves, depth {depth}, stop {}, {:.1} s",
             st.quads_computed, leaves.len(), t.stop_breakdown(), st.wall_seconds);

    // The uniform grid at the same raster, evaluated first so its window colours both panels.
    let sl = grid::Slice::body_plane(res, res, r.cx, r.cy, 0.05, r.body);
    let upx: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();
    let (lo, hi) = colour::range(&upx, Scalar::ShapeSpread);
    let sites = colour::landmarks(&grid::decode_state(&t.chart, r.body, r.cx, r.cy).m);
    let rgb = |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);
    println!("window ({lo:.3e}, {hi:.3e}) from the uniform grid; {} of {} uniform pixels non-finite",
             upx.iter().filter(|p| p.n_nonfinite > 0).count(), upx.len());

    let stem = format!("{root}/{}", region.replace(' ', "_"));
    let (img, tex) = adaptive::render(&t, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, rgb);
    let mut wimg = img.clone();
    wire::draw(&mut wimg, res, res, &wire::boxes_from_tree(&t, &cam, res), depth.max(1));
    let _ = adaptive::save(&format!("{stem}_adaptive.png"), res, &img);
    let _ = adaptive::save(&format!("{stem}_adaptive_wire.png"), res, &wimg);
    println!("texel scaling slope: {:?} (adaptive wants -1)", adaptive::texel_scaling(&tex));

    let cov = adaptive::coverage(&t, &cam, res, &leaves);
    let mut cimg = vec![0u8; res * res * 3];
    for (k, &c) in cov.iter().enumerate() {
        cimg[k * 3..k * 3 + 3].copy_from_slice(&match c {
            0 => [0, 0, 0],
            1 => [90, 200, 120],
            _ => [230, 40, 40],
        });
    }
    let _ = adaptive::save(&format!("{stem}_coverage.png"), res, &cimg);
    println!("coverage: {} gaps, {} overlaps of {} pixels",
             cov.iter().filter(|&&c| c == 0).count(),
             cov.iter().filter(|&&c| c > 1).count(), cov.len());

    let mut ubuf = Vec::with_capacity(upx.len() * 3);
    for p in &upx {
        ubuf.extend_from_slice(&rgb(p));
    }
    let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), res, res, &ubuf);
    println!("wrote {stem}_{{adaptive,adaptive_wire,coverage,uniform}}.png");
}
