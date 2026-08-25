//! **The negative control on the preset fix**, and it runs before the gallery is regenerated.
//!
//! Two comparisons, both at one sample per pixel with no tree in the way:
//!
//! 1. **`shape_pl`, correct pairing against the crossed one, at the same window.** If these two
//!    images do not differ, the cross-coupling diagnosis is wrong and nothing downstream is worth
//!    paying for. The pairing is by GLSL slot -- `alpha` with `pLambda.y` -- and pairing it with
//!    `pLambda.x` instead is a genuinely different 2-plane through the 8D space, not a
//!    reorientation of the same one. Transposing `q1`/`q2` does not recover it, so that is the
//!    third image here rather than an argument in a comment.
//!
//! 2. **The crop: `half = 3.0` against `half = 1.0`, same chart, same basis.** The reference UI
//!    reads `Slice +/- 3.0e+0`; the port shipped 1.0, a 3x zoom on the middle of the picture.
//!
//! **Nothing here writes into `results/`.** A validation run once overwrote committed 1024^2
//! artefacts with small ones, and a small raster reads as a rendering fault rather than a stale
//! file. The output directory is the first argument and defaults to a scratch path.
//!
//! Rendered under **two** colour modes: the shipping bivariate scheme and `event_class/viridis`,
//! which is the mode the reference's WebGPU panel uses. Comparing across colour modes is how a
//! rendering choice gets mistaken for a physics bug, and it is most of what went wrong here.

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, png};
use prin_rs::physics::decoder::Latent;
use rayon::prelude::*;

fn main() {
    let dir: String = std::env::args().nth(1).unwrap_or_else(|| "/tmp/preset_control".into());
    let res: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(512);
    assert!(
        !dir.trim_end_matches('/').ends_with("results"),
        "this is a validation run and must not write into results/"
    );
    let _ = std::fs::create_dir_all(&dir);

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let wide = Chart::preset_shape_pl().default_half();

    // The shipped-wrong basis, and the "fix" that is not one.
    let crossed = Chart::Latent {
        z0: Latent::default(),
        q1: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        q2: [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    };
    let transposed = Chart::Latent {
        z0: Latent::default(),
        q1: [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        q2: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    };

    let cases: Vec<(&str, Chart, f64)> = vec![
        ("shape_pl_correct_h3", Chart::preset_shape_pl(), wide),
        ("shape_pl_crossed_h3", crossed, wide),
        ("shape_pl_crossed_transposed_h3", transposed, wide),
        ("shape_pl_correct_h1", Chart::preset_shape_pl(), 1.0),
        ("shape_correct_h3", Chart::preset_shape(), wide),
        ("shape_correct_h1", Chart::preset_shape(), 1.0),
    ];

    println!(
        "preset control, {res}^2, one sample per pixel, E+1 = {}, t = {}, f64.\n\
         Colour modes: bivariate/spread_shape (shipping) and event_class/viridis (the mode the\n\
         reference panel renders). Comparing ACROSS modes is not a comparison.\n",
        ens.n_extra + 1,
        ens.t_max
    );
    println!(
        "{:>32} {:>6} {:>9} {:>9} {:>7}   event classes",
        "case", "half", "distinct", "nan px", "classes"
    );

    let mut rendered: Vec<(String, Vec<u8>)> = Vec::new();
    for (name, chart, half) in &cases {
        let sl = grid::Slice::body_plane(res, res, 0.0, 0.0, *half, 0).with_chart(*chart);
        let px: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
            .collect();

        let m_here = grid::decode_state(chart, 0, 0.0, 0.0).m;
        let sites = colour::landmarks(&m_here);
        let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
        let (distinct, _, _) = colour::quantisation(&px, Scalar::ShapeSpread);

        let mut biv = Vec::with_capacity(px.len() * 3);
        let mut ev = Vec::with_capacity(px.len() * 3);
        for p in &px {
            biv.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
            ev.extend_from_slice(&png::event_class_rgb(p));
        }
        let _ = adaptive::save_rect(&format!("{dir}/{name}_spread.png"), res, res, &biv);
        let _ = adaptive::save_rect(&format!("{dir}/{name}_event.png"), res, res, &ev);

        let (rows, undet) = png::event_class_histogram(&px);
        let live: Vec<String> = rows
            .iter()
            .filter(|&&(_, n)| n > 0)
            .map(|&(c, n)| format!("{}={n}", png::event_class_name(c)))
            .collect();
        println!(
            "{name:>32} {half:>6.1} {distinct:>9} {undet:>9} {:>7}   {}",
            live.len(),
            live.join(", ")
        );
        rendered.push((format!("{name}_spread"), biv));
        rendered.push((format!("{name}_event"), ev));
    }

    // The gate. A pixel-level difference, not a visual impression -- "they look different" is
    // not a measurement, and the whole point of this run is that it can fail.
    let find = |k: &str| -> &Vec<u8> { &rendered.iter().find(|(n, _)| n == k).unwrap().1 };

    // **The mirror symmetry about the horizontal midline -- a NEGATIVE result, kept.**
    //
    // The eye reads a clear mirror symmetry in `preset_shape` and in the correct `shape_pl`, and
    // not in the crossed one. That looked like a far stronger discriminator than a pixel diff, so
    // it was measured. **It is not one, in either of the two forms tried.**
    //
    // Not at the IC level: `preset_shape` itself does not satisfy an IC-level mirror
    // (`worst |IC(u,v) - mirror_x IC(u,-v)| = 2.69`, the same as the crossed form), so the
    // symmetry is a property of the rendered field rather than of the chart.
    //
    // And not at the pixel level, below. Exact pixel equality is the wrong instrument for a
    // chaotic field: a filament landing one pixel to the other side reads as fully broken, so the
    // statistic is dominated by the fine structure it should be looking past. It **inverts** --
    // the crossed plane scores 0.8007 against the correct one's 0.7112 under event class, and
    // `shape_correct_h3`, the most obviously symmetric image of the set by eye, scores 0.4118.
    //
    // Kept and labelled rather than deleted, because the inference is one a reader will make from
    // the pictures and the number is what says not to. **Nothing is asserted on it**, and the
    // gate below is the pixel diff.
    let mirror_frac = |img: &Vec<u8>| -> f64 {
        let mut same = 0usize;
        for y in 0..res {
            let ym = res - 1 - y;
            for x in 0..res {
                let (a, b) = (3 * (y * res + x), 3 * (ym * res + x));
                if img[a..a + 3] == img[b..b + 3] {
                    same += 1;
                }
            }
        }
        same as f64 / (res * res) as f64
    };
    println!("\nMirror symmetry about the horizontal midline (fraction of pixels equal to their\n\
              vertical reflection). **This does NOT discriminate the three planes** -- it inverts,\n\
              scoring the crossed plane above the correct one. Exact pixel equality is the wrong\n\
              instrument for a chaotic field. Printed so the eyeball inference is refuted in the\n\
              output rather than made silently; nothing is asserted on it.");
    for (name, _, _) in &cases {
        for mode in ["spread", "event"] {
            let k = format!("{name}_{mode}");
            println!("  {k:>38}  {:.4}", mirror_frac(find(&k)));
        }
    }

    let diff = |a: &Vec<u8>, b: &Vec<u8>| -> (usize, f64) {
        let n = a.len() / 3;
        let d = (0..n).filter(|&i| a[3 * i..3 * i + 3] != b[3 * i..3 * i + 3]).count();
        (d, 100.0 * d as f64 / n as f64)
    };
    println!();
    for mode in ["spread", "event"] {
        for (label, a, b) in [
            (
                "shape_pl: correct vs crossed",
                format!("shape_pl_correct_h3_{mode}"),
                format!("shape_pl_crossed_h3_{mode}"),
            ),
            (
                "shape_pl: correct vs crossed-transposed",
                format!("shape_pl_correct_h3_{mode}"),
                format!("shape_pl_crossed_transposed_h3_{mode}"),
            ),
            (
                "shape_pl: half 3.0 vs 1.0",
                format!("shape_pl_correct_h3_{mode}"),
                format!("shape_pl_correct_h1_{mode}"),
            ),
            (
                "shape:    half 3.0 vs 1.0",
                format!("shape_correct_h3_{mode}"),
                format!("shape_correct_h1_{mode}"),
            ),
        ] {
            let (n, pct) = diff(find(&a), find(&b));
            println!("  [{mode:>6}] {label:<42} {n:>9} px differ  ({pct:6.2}%)");
            assert!(
                pct > 1.0,
                "{label} under {mode}: only {pct:.2}% of pixels differ -- the control cannot \
                 fire, and the diagnosis behind this fix is not supported"
            );
        }
    }
    println!("\nWrote {} images to {dir}/", rendered.len());
}
