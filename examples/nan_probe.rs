//! **Where the magenta comes from, and how many trajectories are actually behind it.**
//!
//! `results/glsl/shape.png` carries visible magenta in the structured band. `DEBUG_NAN` is the
//! deliberately loud colour for "this pixel has no value", so the marks are the render working;
//! the question is what is undetermined and whether it is physics or a fault.
//!
//! **Read the block structure before the pixel count.** The adaptive render is
//! nearest-neighbour: one footprint of a level-2 quad paints a texel roughly `res / (4 * N)`
//! pixels on a side, so at `res = 512, N = 8` a single undetermined trajectory paints ~16x16 =
//! 256 pixels. Measured on the committed frames, the magenta is **3 axis-aligned blocks** of
//! 18x19, 18x18 and 20x19 in `shape` and **1** of 20x19 in `plambda` -- four footprints, not
//! 1426 scattered pixels. A pixel count there is a fact about the texel size.
//!
//! [`crate::output::colour::rgb`] has four exits to `DEBUG_NAN` and they mean different things:
//! a non-finite copy anywhere in the ensemble, the two failure states, a non-finite `shape_vec`,
//! and a non-finite scalar. This attributes each one. **A non-finite copy is a measurement
//! outcome -- "this could not be determined" -- and is never discarded**; a `DecodeFailed` is a
//! chart fault and would be.
//!
//! Run: `cargo run --release --example nan_probe [n] [chart]`

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid::Chart;
use prin_rs::outcome::State;
use prin_rs::output::colour::Scalar;
use prin_rs::quad::QuadTree;
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let n: usize = arg(1, 64);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("undetermined-footprint census. N={n} per axis, E+1={}, t={}, f64.",
             ens.n_extra + 1, ens.t_max);
    println!("The four exits of `colour::rgb` to DEBUG_NAN, attributed separately. They are not");
    println!("the same finding: a non-finite COPY is a measurement outcome and is never");
    println!("discarded; a DECODE failure is a chart fault.\n");
    println!("{:>16} {:>8} {:>10} {:>10} {:>10} {:>10} {:>9}",
             "chart", "samples", "nonfin cp", "sim fail", "decode f", "nonfin sh", "total %");

    for (name, chart) in [
        ("preset_shape", Chart::preset_shape()),
        ("preset_plambda", Chart::preset_plambda()),
        ("preset_prho", Chart::preset_prho()),
        ("preset_shape_pl", Chart::preset_shape_pl()),
    ] {
        let half = chart.default_half();
        let t = QuadTree::new(0.0, 0.0, half, n, 0);
        let slice = t.nodes[0].slice(n, 0, chart);
        let px: Vec<_> = (0..n * n)
            .into_par_iter()
            .map(|i| evaluate::<f64>(&slice, i, &ens))
            .collect();

        let (mut cp, mut sf, mut df, mut ns) = (0usize, 0, 0, 0);
        for p in &px {
            // Attribution follows `colour::rgb`'s own order, so the columns partition rather
            // than overlap: a footprint with two faults is counted at the first exit it takes,
            // which is the one that decided the pixel.
            if p.n_nonfinite > 0 {
                cp += 1;
            } else {
                match State::from_bits(p.state) {
                    Some(State::SimFailed) => sf += 1,
                    Some(State::DecodeFailed) | None => df += 1,
                    _ if !p.shape_vec.iter().all(|x| x.is_finite())
                        || !Scalar::ShapeSpread.value(p).is_finite() => ns += 1,
                    _ => {}
                }
            }
        }
        let tot = cp + sf + df + ns;
        println!("{name:>16} {:>8} {cp:>10} {sf:>10} {df:>10} {ns:>10} {:>8.3}%",
                 px.len(), 100.0 * tot as f64 / px.len() as f64);
    }

    println!();
    println!("`preset_shape` is the one chart in this set whose terminations are COLLISIONS --");
    println!("escape_fraction 0.0547 against 0.9894-1.0000 for the momentum slices -- so it is");
    println!("also the one that passes through collision-adjacent shapes, and an undetermined");
    println!("footprint there is the instrument reporting rather than failing. What would be a");
    println!("fault is a DECODE failure, which is the chart handing back a configuration that is");
    println!("not a three-body state at all.");
}
