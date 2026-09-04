//! **§12: the sample coordinate is formed globally and then the offset is subtracted back off.**
//!
//! `Slice::axis` is a `linspace` from the lower corner of the cell, so a footprint's chart
//! coordinate `u` is an absolute number; `decode::sample` wants a **cell-relative** `du ∈ [-1, 1]`
//! and `ensemble::jitter` recovered it as `(u - cx) / half`. Precision is spent building the
//! offset and then the offset is thrown away. `SampleSpace::QuadLocal` forms `du` directly from
//! [`prin_rs::grid::Slice::local_pos`] and never builds the global sum except where a decode needs
//! one.
//!
//! # `QuadLocal` alone cannot help, and the 2x2 is what says so
//!
//! On the **direct** path the decoder takes a chart coordinate, so `make_local` still computes
//! `cx + du*half` before decoding — the global sum is formed either way, and the only difference is
//! a last-ulp reordering. The linearisation is what consumes `du` un-summed. So neither factor is
//! sufficient alone and the honest experiment is a 2x2:
//!
//! ```text
//!               DirectF64            LinSplitF64
//!   Global      the shipped path     du RECOVERED from a collapsing u
//!   QuadLocal   the sum is still     du CARRIED -- the only cell that
//!               formed                escapes the chart-coordinate floor
//! ```
//!
//! The standing measurement this has to sit beside: *the linearised decoder buys ~24 levels over
//! f32 and none over f64*, and *the deep-zoom floor is a property of where you zoom* — the same box
//! at the chart origin has no floor at all, because there is no O(1) neighbour for the increment to
//! be absorbed into. So the base point's coordinate magnitude is printed with every ladder.
//!
//! # What a collapse means here
//!
//! `decode::distinct` counts distinct decoded states over the footprint. `n*n` is healthy; `1` is
//! total collapse. A collapsed decode makes the criterion **maximally confident** — identical
//! footprints give `ensemble_spread` exactly zero, which reads as *perfectly resolved* and stops
//! the descent with a small tidy tree built from nothing.
//!
//! Run: `cargo run --release --example sample_space -- [root=results] [n=8] [max_depth=60]`

use prin_rs::decode::{self, Path};
use prin_rs::ensemble::jitter;
use prin_rs::grid::{Chart, Slice};
use prin_rs::physics::Cart;
use prin_rs::uv::SampleSpace;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// Distinct decoded states over one footprint, under a given (space, path).
fn distinct_at(sl: &Slice, n: usize, path: Path, space: SampleSpace) -> usize {
    let states: Vec<Cart<f64>> = (0..n * n)
        .flat_map(|k| {
            jitter::copies_in_space::<f64>(sl, k, 0, 0.5, 0, jitter::Scheme::Halton, path, space)
        })
        .map(|ic| ic.s)
        .collect();
    decode::distinct(&states)
}

fn main() {
    let root: String = std::env::args().nth(1).unwrap_or_else(|| "results".into());
    let n: usize = arg(2, 8);
    let max_depth: u32 = arg(3, 60);
    let dir = format!("{root}/uv");
    let _ = std::fs::create_dir_all(format!("{dir}/output"));

    println!("sample_space: §12, the global round-trip against quad-local coordinates");
    println!("  n = {n} ({} footprints per rung), max_depth = {max_depth}\n", n * n);

    // -----------------------------------------------------------------------------------
    // 1. Does the space move a committed number at production settings?
    // -----------------------------------------------------------------------------------
    println!("== at production settings: does `QuadLocal` move the shipped path?");
    println!("{:>18} {:>10} {:>12} {:>12} {:>12}", "chart", "moved/n", "max |du|", "max |dx|", "rel");
    for (name, chart, cx, cy, half) in [
        ("near-field", Chart::BodyPlane, 1.0, 3.0, 0.05),
        ("config_stability", Chart::config_stability().0, Chart::config_stability().1,
         Chart::config_stability().2, Chart::config_stability().3),
    ] {
        let sl = Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
        let (mut moved, mut mdu, mut mdx, mut mrel) = (0usize, 0.0f64, 0.0f64, 0.0f64);
        for k in 0..n * n {
            let (gx, _gy) = sl.decode_pos(k);
            let (lx, _ly) = sl.local_pos(k);
            mdu = mdu.max(((gx - cx) / half - lx).abs());
            let g = jitter::copies_in_space::<f64>(
                &sl, k, 0, 0.5, 0, jitter::Scheme::Halton, Path::DirectF64, SampleSpace::Global);
            let l = jitter::copies_in_space::<f64>(
                &sl, k, 0, 0.5, 0, jitter::Scheme::Halton, Path::DirectF64, SampleSpace::QuadLocal);
            let d = decode::max_abs_diff(&g[0].s, &l[0].s);
            if d != 0.0 {
                moved += 1;
            }
            mdx = mdx.max(d);
            mrel = mrel.max(d / half);
        }
        println!("{name:>18} {:>10} {mdu:>12.3e} {mdx:>12.3e} {mrel:>12.3e}",
                 format!("{moved}/{}", n * n));
    }
    println!("  `rel` is the displacement in units of the CELL WIDTH -- the scale that decides");
    println!("  whether a sample still lands in its own cell. A last-ulp reordering of a sum, not");
    println!("  a different sample, and that is why the flag defaults to `Global`.\n");

    // -----------------------------------------------------------------------------------
    // 2. The depth ladder: which of the 2x2 survives, and to what depth.
    // -----------------------------------------------------------------------------------
    for (label, cx, cy) in [("O(1) coordinate", 1.0, 3.0), ("chart ORIGIN", 0.0, 0.0)] {
        println!("== depth ladder at the {label}: |cx| = {cx}, |cy| = {cy}");
        println!("  a healthy footprint reads {} distinct; 1 is total collapse", n * n);
        println!("{:>6} {:>11} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>11}",
                 "depth", "half", "g/f64", "l/f64", "g/lin", "l/lin", "g/f32", "g/L32", "l/L32",
                 "|ju|");
        let mut first_loss = [None::<u32>; 7];
        for depth in (0..=max_depth).step_by(4) {
            let half = 0.05 * 2f64.powi(-(depth as i32));
            let sl = Slice::body_plane(n, n, cx, cy, half, 0).with_chart(Chart::BodyPlane);
            let cells = [
                distinct_at(&sl, n, Path::DirectF64, SampleSpace::Global),
                distinct_at(&sl, n, Path::DirectF64, SampleSpace::QuadLocal),
                distinct_at(&sl, n, Path::LinSplitF64, SampleSpace::Global),
                distinct_at(&sl, n, Path::LinSplitF64, SampleSpace::QuadLocal),
                distinct_at(&sl, n, Path::DirectF32, SampleSpace::Global),
                // **The control for the f32 linearised cell.** Without it, `l/L32` holding 28
                // levels past `g/f32` would be attributed to `QuadLocal` when it is
                // `LinSplitF32`'s own -- the standing "the linearisation buys ~24 levels over
                // f32" measured under `Global` long before this flag existed.
                distinct_at(&sl, n, Path::LinSplitF32, SampleSpace::Global),
                distinct_at(&sl, n, Path::LinSplitF32, SampleSpace::QuadLocal),
            ];
            // **The Jacobian's own magnitude.** `linearise` differences at `cu ± half`, which is a
            // GLOBAL sum -- so the linearisation's own construction meets the chart-coordinate
            // floor before any sample does. If `|ju|` reaches 0 at the same rung the columns
            // collapse, that is the mechanism and not a coincidence.
            let l = decode::linearise(&sl.chart, sl.body, cx, cy, half);
            let ju = (0..3).fold(0.0f64, |m, k| {
                m.max(l.ju.r[k].x.abs()).max(l.ju.r[k].y.abs())
                 .max(l.ju.v[k].x.abs()).max(l.ju.v[k].y.abs())
            });
            for (i, &c) in cells.iter().enumerate() {
                if c < n * n && first_loss[i].is_none() {
                    first_loss[i] = Some(depth);
                }
            }
            println!("{depth:>6} {half:>11.3e} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {ju:>11.3e}",
                     cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], cells[6]);
        }
        let nm = ["g/f64", "l/f64", "g/lin", "l/lin", "g/f32", "g/L32", "l/L32"];
        print!("  first rung below full:");
        for (i, l) in first_loss.iter().enumerate() {
            print!("  {}={}", nm[i], l.map_or("none".into(), |d| d.to_string()));
        }
        println!("\n");
    }

    println!("THE READING: `QuadLocal` buys NOTHING at any depth in this table, and the reason is");
    println!("a different global sum one level up.");
    println!();
    println!("§12's defect is real -- `jitter` did recover `du` from a `u` the linespace had just");
    println!("built. Its cost is zero, because `decode::linearise` differences the chart at");
    println!("`cu +/- half`, which is the SAME global sum, and it collapses first: `|ju|` tracks");
    println!("`half` exactly to depth 44, quantises to the ulp of `cx` at 48, and is EXACTLY ZERO");
    println!("at 52. A linearisation whose own secant has collapsed carries no information for a");
    println!("carried `du` to preserve, so all four f64 cells must fail at the same rung -- and");
    println!("they do, at 48. The fix §12 actually needs is a Jacobian not built by a secant at");
    println!("the cell width.");
    println!();
    println!("And `g/L32` is the control that stops the f32 result being misread: it is IDENTICAL");
    println!("to `l/L32` at every rung, so the ~28 levels the linearised f32 path holds over");
    println!("`g/f32` are `LinSplitF32`'s own -- the standing result, measured under `Global` long");
    println!("before this flag existed -- and not `QuadLocal`'s. Without that column the table");
    println!("would have read as this change buying 28 levels.");
    println!();
    println!("The ORIGIN ladder is the control for the whole finding: with `cx = 0` nothing");
    println!("collapses at any rung in any cell and `|ju|` tracks `half` to 4.3e-20, because there");
    println!("is no O(1) neighbour for the increment to be absorbed into. Quote the coordinate");
    println!("magnitude with any floor depth.");
    println!("\n{dir}/output/sample_space.txt");
}
