//! **Is the deep-zoom cross-hatch an ALIAS of real structure, or structure at the displayed scale?**
//!
//! The `z4` panel's lightness carries a sharp coherent beat at **2.44 px** -- the Nyquist floor.
//! A beat there is what aliasing looks like, and aliasing has one honest test: sample the field
//! finer and average down. If the cross-hatch is the pixel grid folding real sub-pixel structure,
//! it weakens with the supersampling factor. If it survives, it is in the field at the scale
//! being displayed and no amount of sampling removes it.
//!
//! **This is a different supersample from `colour::rgb_resolved`, and the difference is the
//! point.** `rgb_resolved` averages the `E+1` copies' **hue** and holds `l` fixed -- `spread_shape`
//! is a footprint statistic with no per-copy analogue, so the shipped supersampler structurally
//! cannot touch a beat in the lightness channel. This one resamples the *field*: it renders at
//! `f x` the resolution and box-averages, which is ground truth for what the panel would look
//! like if the grid could carry it.
//!
//! Both are reported per OKLab channel, because the two channels come from different arms --
//! `l` from the ensemble spread, `a`/`b` from the nominal copy's shape vector.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, colour as c2};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn main() {
    let out: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(256);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/osc".into());
    let zf: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.015);
    let facs: Vec<usize> = std::env::args().nth(4)
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_else(|| vec![1, 2, 3]);
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0], z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (fcx, fcy, fh) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (cx, cy) = (fcx + fh * (2.0 * 0.10 - 1.0), fcy + fh * (2.0 * 0.45 - 1.0));
    let half = fh * zf;
    let m = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = c2::landmarks(&m);
    let cfg = EnsembleCfg::production().with_overrides(&[
        Override::TMax(50.0), Override::NSync(125), Override::Eta(1e-2),
        Override::RefineFlagged(false), Override::MaxSteps(4_000_000),
        // WITHOUT THIS THE RESOLVED ARM IS INERT AND LOOKS LIKE A NULL. `copy_shapes` is only
        // filled when `keep_copy_shapes` is set, and `rgb_resolved` returns `rgb` unchanged on
        // `copy_shapes.len() < 2`. The first run of this harness produced BITWISE IDENTICAL
        // files for both arms and read as "the supersampler changes nothing".
        Override::KeepCopyShapes(true),
    ]);

    println!("SUPERSAMPLE at zf = {zf}, output {out}^2. `f` is samples per output pixel edge.");
    println!("The colour window is taken ONCE at f = 1 and shared -- a per-factor window would");
    println!("rescale each panel and the amplitudes would stop being comparable.");
    println!();
    let mut window: Option<(f64, f64)> = None;

    for f in facs {
        let res = out * f;
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let (lo, hi) = *window.get_or_insert_with(|| colour::range(&px, Scalar::ShapeSpread));

        // THE ARM IS NOT INERT. Asserted, not assumed: a harness whose two arms are the same
        // code path reports a clean null and reads as a finding.
        let ncs = px.iter().filter(|p| p.copy_shapes.len() >= 2).count();
        let differ = px.iter().filter(|p| {
            colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi)
                != colour::rgb_resolved(p, Scalar::ShapeSpread, &sites, lo, hi)
        }).count();
        println!("  f = {f}  copies kept on {ncs} of {} px;  plain != resolved on {differ} px  -> {}",
                 px.len(), if differ > 0 { "NOT INERT" } else { "INERT -- the comparison is dead" });

        // Two renders of the SAME march. `plain` is the production map; `res` averages the E+1
        // copies' hue. Both are then box-averaged to the output grid, in LINEAR sRGB -- averaging
        // gamma-encoded bytes would darken every edge and put a bias in the comparison.
        for (nm, resolved) in [("plain", false), ("ssres", true)] {
            let mut lin = vec![0.0f64; res * res * 3];
            for (i, p) in px.iter().enumerate() {
                let c = if resolved {
                    colour::rgb_resolved(p, Scalar::ShapeSpread, &sites, lo, hi)
                } else {
                    colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi)
                };
                for k in 0..3 {
                    let v = c[k] as f64 / 255.0;
                    lin[i * 3 + k] =
                        if v <= 0.04045 { v / 12.92 } else { ((v + 0.055) / 1.055).powf(2.4) };
                }
            }
            let mut rgb = Vec::with_capacity(out * out * 3);
            for y in 0..out {
                for x in 0..out {
                    for k in 0..3 {
                        let mut a = 0.0;
                        for dy in 0..f {
                            for dx in 0..f {
                                a += lin[((y * f + dy) * res + x * f + dx) * 3 + k];
                            }
                        }
                        a /= (f * f) as f64;
                        let s = if a <= 0.003_130_8 { 12.92 * a } else { 1.055 * a.powf(1.0 / 2.4) - 0.055 };
                        rgb.push((s.clamp(0.0, 1.0) * 255.0).round() as u8);
                    }
                }
            }
            let _ = adaptive::save_rect(&format!("{dir}/ss_{nm}_f{f}.png"), out, out, &rgb);
        }
        println!("  f = {f}   marched {res}^2 = {} pixels", res * res);
    }
    println!();
    println!("`plain` isolates the FIELD sampling; `ssres` adds the shipped hue supersampler on");
    println!("top. If the cross-hatch falls with `f` in `plain` it is an alias of real structure.");
    println!("If `ssres` moves it and `plain` does not, it was the hue channel; the lightness");
    println!("channel has no per-copy analogue and `rgb_resolved` cannot reach it.");
}
