//! **§3.5 — slice variety: does any of this depend on the slice we happened to pick?**
//!
//! Every experiment before this one is three regions of *one* slice family: one body's position
//! over a 2D box, axis-aligned. If oblique or nonlinear slices behave differently, every prior
//! conclusion is slice-conditional and has to be labelled as such.
//!
//! Three chart types, matched budget and thresholds:
//!
//! - `body_plane` — the historical slice. Affine, axis-aligned.
//! - `plane` — an **oblique** 2-plane in the 6D position space, mixing bodies. Still affine.
//! - `shape` — the shape sphere through Burrau's own configuration. **Nonlinear.**
//!
//! The oblique planes are normalised to the same step length per unit chart coordinate as the
//! axis-aligned one, so a difference in tree shape is not just a difference in how far a unit of
//! `u` moves the system. Without that they would not be comparable at all.
//!
//! Run: `cargo run --release --example slice_variety [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::Chart;
use prin_rs::physics::burrau;
use prin_rs::quad::Decision as D;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::Vec2;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * f).round() as usize]
}

/// An oblique 2-plane: `u` rotates body 0 in its own plane by `theta`, `v` moves body 0 and
/// body 2 together. Unit-normalised so one unit of chart coordinate moves the configuration the
/// same total distance as the axis-aligned chart does.
fn oblique(theta: f64) -> Chart {
    let mut origin = burrau::state::<f64>();
    origin.r[0] = Vec2::zero();
    let (c, s) = (theta.cos(), theta.sin());
    let mut u = [Vec2::zero(); 3];
    let mut v = [Vec2::zero(); 3];
    u[0] = Vec2::new(c, s);
    // Split across two bodies, then renormalise so |V| = 1 in the 6D position metric.
    v[0] = Vec2::new(-s, c) / 2f64.sqrt();
    v[2] = Vec2::new(0.0, 1.0) / 2f64.sqrt();
    Chart::Plane { origin, u, v }
}

fn main() {
    let budget: usize = arg(1, 4000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, E+1=8, t=13, f64");
    println!("screen floor ON, viewport 512x512. Oblique planes are unit-normalised in the 6D");
    println!("position metric, so a tree-shape difference is not a step-length difference.\n");

    // Centres chosen so every chart is looking at a comparable neighbourhood: the body_plane
    // regions, and the same coordinates read as chart coordinates for the other two.
    let cases: Vec<(&str, Chart, f64, f64, usize)> = vec![
        ("body_plane near-field", Chart::BodyPlane, 1.0, 3.0, 0),
        ("body_plane deep int.", Chart::BodyPlane, 0.0, 0.0, 0),
        ("plane obl 30deg", oblique(std::f64::consts::FRAC_PI_6), 1.0, 3.0, 0),
        ("plane obl 45deg", oblique(std::f64::consts::FRAC_PI_4), 1.0, 3.0, 0),
        ("plane obl 45 @origin", oblique(std::f64::consts::FRAC_PI_4), 0.0, 0.0, 0),
        ("shape @burrau", Chart::shape_at_burrau(0.4), 0.0, 0.0, 0),
        ("shape @burrau off", Chart::shape_at_burrau(0.4), 0.1, -0.1, 0),
        ("shape phase 1.3", Chart::shape_at_burrau(1.3), 0.0, 0.0, 0),
    ];

    println!("{:>22} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>11} {:>9} {:>9} {:>9}",
             "case", "chart", "quads", "leaves", "depth", "floor", "keep", "screen",
             "spread med", "alpha p10", "alpha med", "alpha p90");

    for (label, chart, cx, cy, body) in cases {
        let cfg = SchedCfg {
            budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4,
            camera: Some(Camera::framing(cx, cy, 0.05, 512)),
            chart, ..Default::default()
        };
        let (t, st) = scheduler::descend(cx, cy, 0.05, body, &cfg, &ens, Precision::F64);
        let leaves: Vec<usize> = t.leaves().collect();
        let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
        let mut sp: Vec<f64> = leaves.iter().map(|&i| t.nodes[i].red.spread_median)
            .filter(|x| x.is_finite()).collect();
        let mut al: Vec<f64> = t.nodes.iter().filter_map(|n| n.alpha)
            .filter(|x| x.is_finite()).collect();
        println!("{:>22} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>11.3e} {:>9.3} {:>9.3} {:>9.3}",
                 label, chart.name(), st.quads_computed, leaves.len(),
                 t.depth_histogram().len().saturating_sub(1),
                 c(D::Floor), c(D::Keep), c(D::ScreenFloor),
                 q(&mut sp, 0.5), q(&mut al.clone(), 0.1), q(&mut al.clone(), 0.5),
                 q(&mut al, 0.9));
    }

    println!();
    println!("If tree shape, leaf count and the alpha distribution differ by chart TYPE rather");
    println!("than by where the chart is centred, every prior conclusion is slice-conditional.");
    println!("The two 'body_plane' rows and the '@origin' rows are the controls for that: they");
    println!("separate 'which chart' from 'which neighbourhood'.");
}
