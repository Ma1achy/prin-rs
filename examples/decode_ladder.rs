//! §3.4 — the two decode paths against depth, with **distinctness read before divergence**.
//!
//! Run: `cargo run --release --example decode_ladder`

use prin_rs::decode::{self, Path, DECODES_PER_LINEARISATION};
use prin_rs::grid::Chart;
use prin_rs::physics::Cart;

const N: usize = 8;
const PATHS: [Path; 5] = [
    Path::DirectF64,
    Path::DirectF32,
    Path::LinNaiveF32,
    Path::LinSplitF32,
    Path::LinSplitF64,
];

fn samples(p: Path, chart: &Chart, cu: f64, cv: f64, half: f64) -> Vec<Cart<f64>> {
    let lin = decode::linearise(chart, 0, cu, cv, half);
    let d = decode::deltas(N);
    let mut v = Vec::with_capacity(N * N);
    for &dv in &d {
        for &du in &d {
            v.push(decode::sample(p, chart, 0, cu, cv, half, du, dv, &lin));
        }
    }
    v
}

fn worst(a: &[Cart<f64>], b: &[Cart<f64>]) -> f64 {
    a.iter().zip(b).map(|(x, y)| decode::max_abs_diff(x, y)).fold(0.0, f64::max)
}

fn main() {
    let depths: Vec<i32> = vec![0, 8, 14, 16, 18, 20, 22, 24, 26, 30, 35, 40, 44, 45, 46, 47, 48, 50, 52, 55];
    // **Matched coordinate magnitudes.** The floor is set by the magnitude of the coordinate
    // the increment is added to, not by the chart type — a quad centred at the chart origin has
    // no O(1) neighbour to be absorbed into and therefore no such floor at all. Comparing a
    // body_plane quad at (1, 3) against a shape quad at (0, 0) would attribute that difference
    // to the chart. Both are run at both centres instead.
    let charts: [(&str, Chart, f64, f64); 4] = [
        ("body_plane (affine)   @ centre |c| ~ 3", Chart::BodyPlane, 1.0, 3.0),
        ("body_plane (affine)   @ centre 0", Chart::plane_for_body(0), 0.0, 0.0),
        ("shape (nonlinear)     @ centre |c| ~ 3", Chart::shape_at_burrau(0.4), 1.0, 3.0),
        ("shape (nonlinear)     @ centre 0", Chart::shape_at_burrau(0.4), 0.0, 0.0),
    ];

    println!("N = {N}, so {} samples per quad. half0 = 0.05.", N * N);
    println!(
        "Jacobian cost: {DECODES_PER_LINEARISATION} f64 decodes per quad (2 per axis + centre), \
         against {} trajectories per quad at E+1 = 8.\n",
        N * N * 8
    );

    for (label, chart, cu, cv) in charts {
        println!("=== {label} ===");
        println!(
            "**Distinctness first.** A collapsed path agrees perfectly with another collapsed\n\
             path; the agreement means nothing. Divergence columns are printed as '-' wherever\n\
             either side has lost samples.\n"
        );
        print!("{:>5} {:>11}", "depth", "half");
        for p in PATHS {
            print!(" {:>13}", p.name());
        }
        // NOT "curvature": on an affine chart the curvature term is structurally zero, so this
        // column is pure secant/accumulation rounding there. It is (curvature + rounding).
        println!("   {:>11} {:>11} {:>11}", "lin/space", "|d64-Lf32|", "|d64-naive|");

        for &d in &depths {
            let half = 0.05 / (2f64).powi(d);
            let s: Vec<Vec<Cart<f64>>> = PATHS.iter().map(|&p| samples(p, &chart, cu, cv, half)).collect();
            let n: Vec<usize> = s.iter().map(|x| decode::distinct(x)).collect();
            print!("{d:>5} {half:>11.3e}");
            for k in 0..PATHS.len() {
                print!(" {:>10}/{:<2}", n[k], N * N);
            }
            let full = |k: usize| n[k] == N * N;
            // curvature, relative to the sample spacing it would have to exceed to matter
            let spacing = 2.0 * half / (N - 1) as f64;
            let curv = if full(0) && full(4) {
                format!("{:>11.3e}", worst(&s[0], &s[4]) / spacing)
            } else {
                format!("{:>11}", "-")
            };
            let lf = if full(0) && full(3) {
                format!("{:>11.3e}", worst(&s[0], &s[3]))
            } else {
                format!("{:>11}", "-")
            };
            let nv = if full(0) && full(2) {
                format!("{:>11.3e}", worst(&s[0], &s[2]))
            } else {
                format!("{:>11}", "-")
            };
            println!("   {curv} {lf} {nv}");
        }
        println!();
    }

    // Where does the curvature actually matter? Only meaningful on the nonlinear chart.
    println!("=== where the linearisation starts to matter ===");
    println!("|direct - linearised| against the sample spacing. Above 1 the approximation moves a");
    println!("sample by more than the distance to its neighbour. On an affine chart this is");
    println!("STRUCTURALLY zero at every depth — reported as structural, never as a measurement.\n");
    for (label, chart, cu, cv) in [
        ("body_plane (affine)", Chart::BodyPlane, 1.0, 3.0),
        ("shape (nonlinear)", Chart::shape_at_burrau(0.4), 0.0, 0.0),
        ("shape, offset centre", Chart::shape_at_burrau(0.4), 1.0, 3.0),
    ] {
        print!("{label:>22}: ");
        let mut crossed = None;
        for d in 0..40 {
            let half = 0.5 / (2f64).powi(d);
            let a = samples(Path::DirectF64, &chart, cu, cv, half);
            let b = samples(Path::LinSplitF64, &chart, cu, cv, half);
            let r = worst(&a, &b) / (2.0 * half / (N - 1) as f64);
            if d == 0 {
                print!("ratio at half=0.5 is {r:.4e}; ");
            }
            if r < 1.0 && crossed.is_none() {
                crossed = Some(d);
            }
        }
        match crossed {
            Some(0) => println!("below the spacing at every depth tested"),
            Some(d) => println!("falls below the spacing at depth {d}"),
            None => println!("never falls below the spacing (would be a design problem)"),
        }
    }
}
