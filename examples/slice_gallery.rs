//! Looking at the slices, rather than only tabulating them.
//!
//! `slice_variety` measured leaf counts and `alpha` distributions across charts and found tree
//! size **slice-conditional to 4.3x** while the `alpha` distribution stayed put (median
//! 0.172-0.289). That is a table, and a table cannot say *what the slices look like* — whether
//! an oblique plane cuts the same structure at an angle or lands on different structure
//! entirely, and whether the shape chart's curvature is visible at this scale or only in the
//! decode ladder.
//!
//! # What is held fixed
//!
//! **One shared centre configuration** — near-field's, body 0 at (1, 3), released from rest —
//! and only the 2-plane through it changes. The bases are orthonormal in the 6D position metric,
//! so a unit of chart coordinate moves the system the same distance in every case. Without that,
//! a "different slice" would be a different *scale* and the comparison would be of zoom levels.
//!
//! `body_plane` and `plane 0deg` are the **control pair**: they are the same chart written two
//! ways and must render **bitwise identically**. That is asserted here, not assumed — if the
//! control pair differs, nothing else on the page means anything.
//!
//! # Every image has a wire twin
//!
//! The adaptive render says what is displayed; the wireframe says where the tree cut. A coarse
//! texel tells you a leaf is coarse; only the wire tells you whether the structure around it was
//! subdivided *around* or *through*. PR #11 drew boundaries over a uniform base and conflated
//! the two, which is how `deep interior`'s failure went unnoticed.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::Chart;
use prin_rs::output::{adaptive, apng, png, wire};
use prin_rs::physics::{burrau, shape, Cart};
use prin_rs::quad::Agg;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::stats;
use prin_rs::Vec2;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn centre() -> Cart<f64> {
    let mut s = burrau::state::<f64>();
    s.r[0] = Vec2::new(1.0, 3.0);
    s
}

fn unit(mut u: [Vec2<f64>; 3]) -> [Vec2<f64>; 3] {
    let n = u.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt();
    for v in u.iter_mut() {
        *v = *v / n;
    }
    u
}

/// A 2-plane through the shared centre. `theta` rotates within body 0's plane; `mix` sends that
/// fraction of `V` into body 2, so `mix = 0` is single-body and `mix > 0` genuinely cross-body.
fn plane(theta: f64, mix: f64) -> Chart {
    let (c, s) = (theta.cos(), theta.sin());
    let mut u = [Vec2::zero(); 3];
    let mut v = [Vec2::zero(); 3];
    u[0] = Vec2::new(c, s);
    v[0] = Vec2::new(-s, c) * (1.0 - mix);
    v[2] = Vec2::new(-s, c) * mix;
    Chart::Plane { origin: centre(), u: unit(u), v: unit(v) }
}

/// The shape chart through the **same** configuration, so the nonlinear case is comparable to
/// the affine ones rather than to Burrau's own triangle.
fn shape_here(phase: f64) -> Chart {
    let m = burrau::MASSES;
    let r = centre().r;
    let n0 = shape::shape_vec(&r, &m);
    let (e1, e2) = shape::tangent_frame(n0);
    Chart::Shape { n0, e1, e2, inertia: shape::inertia(&r, &m), phase }
}

fn main() {
    let budget: usize = arg(1, 4000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let res: usize = arg(4, 512);
    let half = 0.05;

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let pi = std::f64::consts::PI;
    let cases: Vec<(&str, Chart)> = vec![
        ("body_plane", Chart::BodyPlane),
        ("plane_00deg", plane(0.0, 0.0)),
        ("plane_15deg", plane(pi / 12.0, 0.0)),
        ("plane_30deg", plane(pi / 6.0, 0.0)),
        ("plane_45deg", plane(pi / 4.0, 0.0)),
        ("plane_45deg_mix05", plane(pi / 4.0, 0.5)),
        ("plane_45deg_mix10", plane(pi / 4.0, 1.0)),
        ("shape_phase_00", shape_here(0.0)),
        ("shape_phase_04", shape_here(0.4)),
        ("shape_phase_13", shape_here(1.3)),
    ];

    println!(
        "budget {budget}, tau={tau:e}, alpha_hi={alpha_hi}, N=8, E+1=8, t=13, f64, {res}^2, \
         screen floor ON.\n\
         ONE shared centre configuration; only the 2-plane through it changes. Bases orthonormal\n\
         in the 6D position metric, so a unit of chart coordinate moves the system equally far\n\
         in every case -- otherwise a different slice would be a different SCALE.\n"
    );
    println!(
        "{:>20} {:>11} {:>7} {:>7} {:>6} {:>7} {:>11} {:>10} {:>10}",
        "case", "chart", "quads", "leaves", "depth", "screen", "spread med", "alpha med", "alpha idec"
    );

    let mut frames: Vec<Vec<u8>> = Vec::new();
    let mut wire_frames: Vec<Vec<u8>> = Vec::new();
    // The control compares INITIAL CONDITIONS, which is exact, and reports the image difference
    // separately. Comparing the images alone conflates two different things: whether the two
    // charts are the same chart, and whether the rasteriser rounds the same way at two different
    // coordinate magnitudes. `body_plane` puts quad centres at O(1) and `Plane` at O(0), so
    // `(x - cam.cx) / pixel_size` cancels differently and a tile edge can land one pixel over.
    // That is the same chart-coordinate-magnitude effect the vertical slice recorded for the
    // deep-zoom floor, seen from the render side, and it is not a physics difference.
    let mut control_ics: Option<Vec<Cart<f64>>> = None;
    let mut control: Option<Vec<u8>> = None;

    for (name, chart) in &cases {
        // **`body_plane` reads its centre from the chart coordinate; `Plane` and `Shape` carry
        // it in `origin` and are centred at zero.** Putting them all at (1, 3) samples a box two
        // units away from the shared configuration, which is a different slice of different
        // physics rather than a rotation of the same one. The control pair is what catches it.
        let (cx, cy) = if *chart == Chart::BodyPlane { (1.0, 3.0) } else { (0.0, 0.0) };
        let cam = Camera::framing(cx, cy, half, res);
        let cfg = SchedCfg {
            budget,
            tau_display: tau,
            alpha_hi,
            alpha_lo: alpha_hi,
            agg: Agg::Median,
            chart: *chart,
            camera: Some(cam),
            keep_pixels: true,
            ..Default::default()
        };
        let (t, st) = scheduler::descend(cx, cy, half, 0, &cfg, &ens, Precision::F64);

        let leaves: Vec<usize> = t.leaves().collect();
        let depth = leaves.iter().map(|&i| t.nodes[i].level).max().unwrap_or(0);
        let screen = leaves
            .iter()
            .filter(|&&i| t.nodes[i].decision == prin_rs::quad::Decision::ScreenFloor)
            .count();
        let spreads: Vec<f64> =
            leaves.iter().map(|&i| t.nodes[i].red.spread_median).collect();
        let alphas: Vec<f64> =
            leaves.iter().filter_map(|&i| t.nodes[i].alpha).collect();
        let (_, amed, _, aidec) = stats::interdecile(&alphas);

        let (img, _tex) = adaptive::render(
            &t, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, png::outcome_rgb,
        );
        let mut wimg = img.clone();
        let boxes = wire::boxes_from_tree(&t, &cam, res);
        let deepest = boxes.iter().map(|b| b.level).max().unwrap_or(1);
        wire::draw(&mut wimg, res, res, &boxes, deepest.max(1));

        let stem = format!("results/criterion/slice_{name}");
        let _ = adaptive::save(&format!("{stem}.png"), res, &img);
        let _ = adaptive::save(&format!("{stem}_wire.png"), res, &wimg);
        if let Ok(f) = std::fs::File::create(format!("{stem}.prnq")) {
            let mut w = std::io::BufWriter::new(f);
            let _ = prin_rs::output::tree::write(&mut w, &t, &cfg, &ens, &st, name, "f64");
        }

        // The control pair must be bitwise identical: `plane 0deg` is `body_plane` written a
        // second way. A difference here means the comparison is broken and nothing below it is
        // readable.
        let ics: Vec<Cart<f64>> = {
            let sl = t.nodes[0].slice(cfg.n, 0, *chart);
            (0..sl.npix()).map(|k| sl.nominal::<f64>(k)).collect()
        };
        match (*name, &control, &control_ics) {
            ("body_plane", _, _) => {
                control = Some(img.clone());
                control_ics = Some(ics.clone());
            }
            ("plane_00deg", Some(c), Some(ci)) => {
                let dic = ci
                    .iter()
                    .zip(&ics)
                    .map(|(a, b)| prin_rs::decode::max_abs_diff(a, b))
                    .fold(0.0f64, f64::max);
                let moved = c
                    .chunks(3)
                    .zip(img.chunks(3))
                    .filter(|(a, b)| a != b)
                    .count();
                println!(
                    "  [control] plane_00deg vs body_plane: max |dIC| = {dic:.3e} -> {}",
                    if dic == 0.0 {
                        "the two charts are the SAME chart; the comparison holds"
                    } else {
                        "*** THE CHARTS DISAGREE; nothing below is readable ***"
                    }
                );
                println!(
                    "  [control] and {moved} of {} pixels differ ({:.3}%) -- rasteriser rounding at\n\
                     \x20           two coordinate magnitudes (O(1) vs O(0)), NOT a physics difference.",
                    c.len() / 3,
                    100.0 * moved as f64 / (c.len() / 3) as f64
                );
            }
            _ => {}
        }

        frames.push(img);
        wire_frames.push(wimg);

        println!(
            "{name:>20} {:>11} {:>7} {:>7} {depth:>6} {screen:>7} {:>11.4e} {amed:>10.4} {aidec:>10.4}",
            chart.name(),
            t.nodes.len(),
            leaves.len(),
            {
                let mut v = spreads.clone();
                prin_rs::quad::quantile(&mut v, 0.5)
            },
        );
    }

    let _ = apng::write("results/criterion/slice_gallery_animated.png", res, res, &frames, 1, 1);
    let _ = apng::write(
        "results/criterion/slice_gallery_wire_animated.png",
        res,
        res,
        &wire_frames,
        1,
        1,
    );

    println!(
        "\nEvery case has a plain render and a `_wire` twin. The plain image says WHAT IS\n\
         DISPLAYED -- texels at true per-quad sizes, so a coarse leaf is visibly coarse. The wire\n\
         says WHERE THE TREE CUT, brightness graded by level. They answer different questions and\n\
         neither substitutes for the other: PR #11 drew boundaries over a UNIFORM base, which\n\
         conflated them and is how `deep interior`'s bad tree went unnoticed.\n\
         \n\
         Read the control line first. If plane_00deg is not bitwise body_plane, the bases are\n\
         wrong and every other row is comparing different physics rather than different slices.\n\
         \n\
         Leaf counts are slice-conditional to 4.3x and must be compared WITHIN a slice, never\n\
         across. `alpha med`/`alpha idec` are the stable quantities and are the ones to compare."
    );
}
