//! **§3.3 — SSAA resolve, and what it does to the image at each `E`.**
//!
//! The ensemble has been used for exactly one thing: `ensemble_spread`, a **disagreement**
//! statistic that drives scheduling. Its other job is **resolve** — many sub-pixel samples to one
//! pixel colour — an *average*, and that path had never run.
//!
//! The two must not be confused. A footprint whose copies split 4/4 has a large spread *and* a
//! blended colour; the spread says "refine here", the blend says "this is what the pixel looks
//! like". Neither substitutes for the other.
//!
//! **Reported per footprint, not only as an aggregate.** An aggregate can only say the
//! distribution did not move; it cannot say the footprints did not. That mistake has been made
//! twice on this project in one PR.
//!
//! Run: `cargo run --release --example ssaa_resolve [res]`

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::png::outcome_rgb;
use prin_rs::output::ssaa;
use prin_rs::render::{self, Precision};

fn main() -> std::io::Result<()> {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);

    println!("uniform render {res}x{res}, t=13, eta=1e-2, f64, Halton offsets, jitter_frac 0.5.");
    println!("Resolve = mean of the copies' outcome colours. Spread is NOT resolve.\n");
    println!("{:>14} {:>5} {:>12} {:>12} {:>12} {:>12} {:>10}",
             "region", "E+1", "moved frac", "mean |dRGB|", "max |dRGB|", "spread med", "wall_s");

    for region in ["near-field", "deep interior"] {
        let slice = grid::region(region, res, res, 0.05).unwrap();
        for e1 in [1usize, 2, 4, 8, 16, 32] {
            let cfg = EnsembleCfg {
                n_extra: e1 - 1,
                keep_copy_outcomes: true,
                refine_flagged: false,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let px: Vec<PixelOut> = render::render(&slice, &cfg, Precision::F64);
            let (mean, worst, frac) = ssaa::resolve_effect(&px);
            let mut sp: Vec<f64> =
                px.iter().map(|p| p.ensemble_spread).filter(|x| x.is_finite()).collect();
            sp.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let med = sp.get(sp.len() / 2).cloned().unwrap_or(f64::NAN);
            println!("{:>14} {:>5} {:>12.4} {:>12.3} {:>12.0} {:>12.3e} {:>10.1}",
                     region, e1, frac, mean, worst, med, t0.elapsed().as_secs_f64());

            // The images, at the two ends and the default, so the effect is visible not described.
            if e1 == 1 || e1 == 8 || e1 == 32 {
                let stem = format!("results/vertical/ssaa_{}_{e1}", region.replace(' ', "_"));
                let mut a = Vec::with_capacity(px.len() * 3);
                let mut b = Vec::with_capacity(px.len() * 3);
                for p in &px {
                    a.extend_from_slice(&outcome_rgb(p));
                    b.extend_from_slice(&ssaa::resolve_rgb(p));
                }
                prin_rs::output::adaptive::save(&format!("{stem}_nominal.png"), res, &a)?;
                prin_rs::output::adaptive::save(&format!("{stem}_resolved.png"), res, &b)?;
            }
        }
        println!();
    }

    println!("At E+1 = 1 the resolve is the nominal copy by definition, so 'moved frac' must be");
    println!("0.0000 there. It is the control: a nonzero value would mean the two paths disagree");
    println!("about the same single copy, which would be a bug in the resolve and not a finding.");
    Ok(())
}
