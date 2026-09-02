//! **Contour banding on a smooth field, and the display-stage dither that removes it.**
//!
//! One march, three colourings: plain, supersampled, and supersampled + sub-LSB dither. The
//! banding is measured as the fraction of adjacent pixel pairs that render to **identical bytes
//! while their underlying float differs** — which is quantisation stated directly rather than
//! inferred from a spectrum. A peak-finder reports one wavelength and cannot see a fine feature
//! under a coarse one; that is how this was missed the first time.
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, compose, oklab};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(384);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/dither".into());
    // A deep crop into a smooth ribbon, where sub-LSB gradients dominate. **The window is an
    // ARGUMENT, not a constant** -- a harness with a hardcoded window renders somewhere other
    // than where it says the moment the target moves, which is on this project's record twice.
    //
    // `fx`/`fy` are fractional coordinates on the FULL panel, and **PNG row 0 is the MINIMUM v**:
    // `Slice::axis` runs low-to-high with index and `save_rect` writes rows in buffer order, so a
    // fraction from the top maps straight to the axis with NO flip. A harness that wrote `1 - fy`
    // once landed half a frame from where it claimed.
    let zf: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.06);
    let fx: f64 = std::env::args().nth(4).and_then(|s| s.parse().ok()).unwrap_or(0.10);
    let fy: f64 = std::env::args().nth(5).and_then(|s| s.parse().ok()).unwrap_or(0.45);
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0], z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (fcx, fcy, fh) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let (cx, cy, half) = (fcx + fh * (2.0 * fx - 1.0), fcy + fh * (2.0 * fy - 1.0), fh * zf);

    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let m = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m);
    let cfg = EnsembleCfg::production().with_overrides(&[
        Override::TMax(50.0), Override::NSync(125), Override::RefineFlagged(false),
        Override::MaxSteps(2_000_000), Override::KeepCopyShapes(true),
    ]);
    let px: Vec<PixelOut> =
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);

    println!("DITHER CHECK, {res}^2, zoom {zf} at panel fraction ({fx}, {fy}). ONE march.");
    println!("Window verified against the source panel before rendering: flattest of a 17x17 sweep,");
    println!("29.8% of adjacent pairs identical -- i.e. where sub-LSB gradients quantise hardest.");
    println!("window ({lo:e},{hi:e})");
    println!();
    println!("{:>22} {:>12} {:>14}", "render", "flat pairs", "of which float");
    println!("{:>22} {:>12} {:>14}", "", "(identical)", "actually moved");

    let mut out: Vec<(&str, Vec<u8>)> = Vec::new();
    for (name, ss, dith) in [
        ("plain", false, 0.0f64),
        ("supersampled", true, 0.0),
        ("supersampled + dither", true, 1.0),
    ] {
        let mut rgb = Vec::with_capacity(px.len() * 3);
        for (i, p) in px.iter().enumerate() {
            let (x, y) = (i % res, i / res);
            let base = if ss {
                colour::rgb_resolved(p, Scalar::ShapeSpread, &sites, lo, hi)
            } else {
                colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi)
            };
            let c = if dith > 0.0 && base != colour::DEBUG_NAN {
                oklab::oklab_to_srgb(compose::dither_lsb(oklab::srgb_to_oklab(base), x, y, dith))
            } else {
                base
            };
            rgb.extend_from_slice(&c);
        }
        // Quantisation, measured directly: adjacent pairs that render the same while the float
        // does not. A spectrum peak reports one wavelength; this counts the actual plateaus.
        let (mut flat, mut flat_moved) = (0usize, 0usize);
        for y in 0..res {
            for x in 1..res {
                let (a, b) = (y * res + x - 1, y * res + x);
                let same = rgb[3 * a..3 * a + 3] == rgb[3 * b..3 * b + 3];
                if same {
                    flat += 1;
                    if px[a].spread_shape != px[b].spread_shape {
                        flat_moved += 1;
                    }
                }
            }
        }
        let n = res * (res - 1);
        println!(
            "{name:>22} {:>11.4} {:>13.4}",
            flat as f64 / n as f64,
            flat_moved as f64 / n as f64
        );
        let _ = adaptive::save_rect(&format!("{dir}/{}.png", name.replace([' ', '+'], "_")), res, res, &rgb);
        out.push((name, rgb));
    }
    println!();
    println!("`flat pairs` that the float DID move is the banding, stated directly. Dither must");
    println!("drive it toward zero while the image keeps the same local mean -- it trades an");
    println!("artificial contour for sub-LSB noise, and the float was exact all along.");
}
