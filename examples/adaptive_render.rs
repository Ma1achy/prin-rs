//! **§3.2 and §4.2 q3 — the adaptive render, at true per-quad texel sizes.**
//!
//! PR #11's overlay drew leaf boundaries over a *uniform* render, so every texel was the same
//! size. It answered "where did the boundaries fall", not "what does the system display", and
//! the tree's quality could not be judged by eye at all. Both are produced here, and the texel
//! scaling of each is measured: **an adaptive render fits `2^-level` exactly; a uniform one does
//! not, and the same assertion has to reject it.**
//!
//! Also produced: the resolved (SSAA) adaptive render, so §3.3's effect is visible in the same
//! frame rather than described.
//!
//! Run: `cargo run --release --example adaptive_render [region] [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::output::adaptive::{self, TexelMode};
use prin_rs::output::png::outcome_rgb;
use prin_rs::output::ssaa;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn spread_rgb(p: &prin_rs::ensemble::pixel::PixelOut) -> [u8; 3] {
    // Log ramp over a fixed window, so panels from different runs are comparable. The window is
    // printed; a false-colour image without its scale is decoration.
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
    let budget: usize = arg(2, 6000);
    let tau: f64 = arg(3, 1e-4);
    let alpha_hi: f64 = arg(4, 0.2);
    const RES: usize = 512;

    let root = grid::region(&region, 2, 2, 0.05).expect("unknown region");
    let ens = EnsembleCfg { refine_flagged: false, keep_copy_outcomes: true, ..Default::default() };
    let cam = Camera::framing(root.cx, root.cy, 0.05, RES);
    let cfg = SchedCfg {
        budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4,
        camera: Some(cam), keep_pixels: true, ..Default::default()
    };

    let (t, st) = scheduler::descend(root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
    let leaves: Vec<usize> = t.leaves().collect();
    let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
    println!("region {region}, budget {budget}, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, E+1=8");
    println!("viewport {RES}x{RES}, camera framing the root box, screen floor ON");
    println!("{} quads, {} leaves, depth {}, floor {} keep {} screen {}, {:.1} s",
             st.quads_computed, leaves.len(), t.depth_histogram().len() - 1,
             c(D::Floor), c(D::Keep), c(D::ScreenFloor), st.wall_seconds);

    let stem = format!("results/vertical/{}", region.replace(' ', "_"));
    let report = |suffix: &str, mode: TexelMode, f: &dyn Fn(&prin_rs::ensemble::pixel::PixelOut) -> [u8; 3]| -> std::io::Result<Option<f64>> {
        let (img, info) = adaptive::render(&t, &st.pixels, &cam, RES, mode, f);
        adaptive::save(&format!("{stem}_{suffix}.png"), RES, &img)?;
        let slope = adaptive::texel_scaling(&info);
        let tiles: usize = info.iter().map(|x| x.tiles_drawn).sum();
        let (lo, hi) = info.iter().fold((f64::MAX, f64::MIN), |(a, b), x| {
            (a.min(x.texel_px), b.max(x.texel_px))
        });
        println!("  {suffix:<22} texel px {lo:.3}..{hi:.3}  tiles {tiles}  slope {}",
                 slope.map(|s| format!("{s:+.6}")).unwrap_or_else(|| "n/a (one level)".into()));
        Ok(slope)
    };

    println!("\nspread window for the false colour: 1e-8 .. 1e-1, log.");
    let ad = report("adaptive_spread", TexelMode::Adaptive, &spread_rgb)?;
    let un = report("uniform_spread", TexelMode::Uniform, &spread_rgb)?;
    report("adaptive_outcome", TexelMode::Adaptive, &outcome_rgb)?;
    report("adaptive_resolved", TexelMode::Adaptive, &ssaa::resolve_rgb)?;

    println!("\n**The acceptance test.** Texel size must vary as 2^-level, i.e. slope -1.");
    match (ad, un) {
        (Some(a), Some(u)) => {
            println!("  adaptive slope {a:+.6}  -> {}", if (a + 1.0).abs() < 1e-12 { "PASS" } else { "FAIL" });
            println!("  uniform  slope {u:+.6}  -> {} (it MUST fail; a render where every texel",
                     if (u + 1.0).abs() < 1e-12 { "FAIL: it passed" } else { "correctly rejected" });
            println!("     is the same size is the PR #11 instrument, not the system)");
        }
        _ => println!("  only one leaf level in this tree; the scaling fit is undefined and no"),
    }

    // SSAA, in the same frame.
    let all: Vec<prin_rs::ensemble::pixel::PixelOut> =
        leaves.iter().flat_map(|&i| st.pixels[i].iter().cloned()).collect();
    let (mean, worst, frac) = ssaa::resolve_effect(&all);
    println!("\nSSAA resolve over {} leaf footprints:", all.len());
    println!("  mean |dRGB| {mean:.3}  max {worst:.0}  fraction of footprints that moved {:.4}", frac);
    println!("  (per-footprint, not an aggregate: an aggregate can only say the distribution did");
    println!("   not move, never that the footprints did not.)");
    Ok(())
}
