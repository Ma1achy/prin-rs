//! **§3.5 — slice variety: does any of this depend on the slice we happened to pick?**
//!
//! Every experiment before this one is three regions of *one* slice family: one body's position
//! over a 2D box, axis-aligned. If oblique or nonlinear slices behave differently, every prior
//! conclusion is slice-conditional and has to be labelled as such.
//!
//! # The comparison has to hold the configuration fixed, not the coordinates
//!
//! A first version of this ran oblique planes at the same *chart coordinates* `(1.0, 3.0)` as
//! near-field. That is not an orientation test: an oblique plane evaluated at those coordinates
//! lands on a completely different configuration, so the rows differed because they were looking
//! at different physics, not because the slice was rotated. The tamer spreads it reported were
//! about the neighbourhood, not the orientation.
//!
//! Here **every case shares one centre configuration** — near-field's, body 0 at `(1, 3)`, which
//! is Burrau's own — and the chart centre is `(0, 0)` in all of them. Only the 2-plane through
//! that configuration changes. `theta = 0` with `U = e(0,x)`, `V = e(0,y)` is the control.
//!
//! **The control check is on the initial conditions, not on the tree.** The initial conditions
//! are what the two charts must agree about; the tree is downstream of a chaotic integration, so
//! a check on the tree would be testing chaos rather than the charts. Measured, the agreement is
//! **exact** — `max |dIC| = 0`, and the two trees are identical quad for quad — so the
//! amplification line below reads zero. It is printed anyway, because if a future change makes
//! the charts differ in the last bit, that line is what says how far the difference travels.
//!
//! **The `shape phase` rows are a gauge check that could have failed.** The fibre phase is a
//! global rotation of the configuration, and the three-body problem is rotation-invariant, so
//! every phase must give a bitwise identical tree. If the Hopf inverse or the AZ port broke
//! rotational invariance, these rows would separate.
//!
//! Basis vectors are orthonormal in the 6D position metric, so one unit of chart coordinate moves
//! the configuration the same distance in every case.
//!
//! Run: `cargo run --release --example slice_variety [budget] [tau] [alpha_hi]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::Chart;
use prin_rs::physics::{burrau, shape, Cart};
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

/// The shared centre: near-field's, body 0 at (1, 3), released from rest.
fn centre() -> Cart<f64> {
    let mut s = burrau::state::<f64>();
    s.r[0] = Vec2::new(1.0, 3.0);
    s
}

fn norm6(u: &[Vec2<f64>; 3]) -> f64 {
    u.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt()
}

fn unit(mut u: [Vec2<f64>; 3]) -> [Vec2<f64>; 3] {
    let n = norm6(&u);
    for v in u.iter_mut() {
        *v = *v / n;
    }
    u
}

/// A 2-plane through the shared centre. `theta` rotates within body 0's plane; `mix` sends that
/// fraction of `V` into body 2 instead, so `mix = 0` is a single-body plane and `mix > 0` is a
/// genuinely cross-body one.
fn plane(theta: f64, mix: f64) -> Chart {
    let (c, s) = (theta.cos(), theta.sin());
    let mut u = [Vec2::zero(); 3];
    let mut v = [Vec2::zero(); 3];
    u[0] = Vec2::new(c, s);
    v[0] = Vec2::new(-s, c) * (1.0 - mix);
    v[2] = Vec2::new(-s, c) * mix;
    Chart::Plane { origin: centre(), u: unit(u), v: unit(v) }
}

/// The shape chart through the **same** configuration, so the nonlinear row is comparable to the
/// affine ones rather than to Burrau's own triangle.
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
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let half = 0.05;

    println!("budget {budget} quads, tau={tau:.0e}, alpha_hi={alpha_hi}, N=8, E+1=8, t=13, f64");
    println!("screen floor ON, viewport 512x512. half={half}.\n");
    println!("ONE shared centre configuration (near-field's: body 0 at (1,3), released from rest).");
    println!("Only the 2-plane through it changes. Bases are orthonormal in the 6D position");
    println!("metric, so a unit of chart coordinate moves the system the same distance in each.\n");
    println!("The 'plane 0deg' row is the CONTROL: it must reproduce 'body_plane' exactly. If it");
    println!("does not, the comparison is broken and nothing below it means anything.\n");

    let pi = std::f64::consts::PI;
    let cases: Vec<(&str, Chart)> = vec![
        ("body_plane (control)", Chart::BodyPlane),
        ("plane 0deg (control)", plane(0.0, 0.0)),
        ("plane 15deg", plane(pi / 12.0, 0.0)),
        ("plane 30deg", plane(pi / 6.0, 0.0)),
        ("plane 45deg", plane(pi / 4.0, 0.0)),
        ("plane 45deg mix 0.5", plane(pi / 4.0, 0.5)),
        ("plane 45deg mix 1.0", plane(pi / 4.0, 1.0)),
        ("shape phase 0.0", shape_here(0.0)),
        ("shape phase 0.4", shape_here(0.4)),
        ("shape phase 1.3", shape_here(1.3)),
    ];

    println!("{:>22} {:>11} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>11} {:>9} {:>9} {:>9}",
             "case", "chart", "quads", "leaves", "depth", "floor", "keep", "screen",
             "spread med", "alpha p10", "alpha med", "alpha p90");

    let mut control: Option<(usize, usize, f64)> = None;
    let mut control_ics: Option<Vec<Cart<f64>>> = None;
    for (label, chart) in cases {
        // body_plane reads its centre from the chart coordinate; the others carry it in `origin`
        // and are centred at zero.
        let (cx, cy) = if chart == Chart::BodyPlane { (1.0, 3.0) } else { (0.0, 0.0) };
        let cfg = SchedCfg {
            budget, tau_display: tau, alpha_hi, alpha_lo: alpha_hi * 0.4,
            camera: Some(Camera::framing(cx, cy, half, 512)),
            chart, ..Default::default()
        };
        let (t, st) = scheduler::descend(cx, cy, half, 0, &cfg, &ens, Precision::F64);
        let leaves: Vec<usize> = t.leaves().collect();
        let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
        let mut sp: Vec<f64> = leaves.iter().map(|&i| t.nodes[i].red.spread_median)
            .filter(|x| x.is_finite()).collect();
        let mut al: Vec<f64> = t.nodes.iter().filter_map(|n| n.alpha)
            .filter(|x| x.is_finite()).collect();
        let med = q(&mut sp, 0.5);
        println!("{:>22} {:>11} {:>7} {:>7} {:>6} {:>6} {:>6} {:>6} {:>11.3e} {:>9.3} {:>9.3} {:>9.3}",
                 label, chart.name(), st.quads_computed, leaves.len(),
                 t.depth_histogram().len().saturating_sub(1),
                 c(D::Floor), c(D::Keep), c(D::ScreenFloor),
                 med, q(&mut al.clone(), 0.1), q(&mut al.clone(), 0.5), q(&mut al, 0.9));
        // The root quad's own initial conditions, which is what the control has to agree on.
        let ics: Vec<Cart<f64>> = {
            let sl = t.nodes[0].slice(cfg.n, 0, chart);
            (0..sl.npix()).map(|k| sl.nominal::<f64>(k)).collect()
        };
        if label == "body_plane (control)" {
            control = Some((st.quads_computed, leaves.len(), med));
            control_ics = Some(ics.clone());
        }
        if label == "plane 0deg (control)" {
            if let (Some((cq, cl, cm)), Some(ci)) = (control, control_ics.as_ref()) {
                let dic = ci.iter().zip(&ics)
                    .map(|(a, b)| prin_rs::decode::max_abs_diff(a, b))
                    .fold(0.0f64, f64::max);
                let ok = dic < 1e-15;
                println!("{:>22} CONTROL: max |dIC| vs body_plane = {dic:.3e} -> {}",
                         "", if ok { "the two charts agree; the comparison holds" }
                             else { "THE CHARTS DISAGREE; nothing below is readable" });
                println!("{:>22}          that IC difference moves the tree by {} quads, {} leaves,",
                         "", (st.quads_computed as i64 - cq as i64).abs(),
                         (leaves.len() as i64 - cl as i64).abs());
                println!("{:>22}          and the median spread by {:.3}% — the chaotic amplification",
                         "", 100.0 * (med - cm).abs() / cm.abs());
                println!("{:>22}          of a last-bit change, and the floor under every row below.",
                         "");
            }
        }
    }

    println!();
    println!("Read the two control rows first. Then: if tree shape, leaf count and the alpha");
    println!("distribution move with ORIENTATION at a fixed centre, every prior conclusion on");
    println!("this project is slice-conditional and has to be labelled as such.");
}
