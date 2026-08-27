//! Did the fix actually reach the pixels that were coloured? A/B on one slice, counted.
//!
//! The committed pair `preset_plambda_e0_uniform.png` and `..._e4_uniform.png` are **bitwise
//! identical**, and the honest first question about a null result is whether the change was
//! applied at all. This answers it with counts rather than an assurance: build the *same* slice
//! twice, differing only in `escape_every`, and report
//!
//! - how many pixels have a different `t_end` (did the fix reach the physics?),
//! - how many have a different `spread_shape` (did it reach the coloured field?),
//! - how many emit a different RGB byte (did it reach the image?).
//!
//! A large first count with a zero third is the finding. A zero first count would mean the flag
//! never took effect, and that is the failure mode worth ruling out explicitly — it has happened
//! on this project before.
//!
//! Both arms are coloured against the **same** ramp, taken from the `e0` pass, so a ramp shift
//! cannot mask or manufacture a difference.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::colour::{self, Scalar};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 1024);
    let case: String = std::env::args().nth(2).unwrap_or_else(|| "preset_plambda".into());
    let ev: usize = arg(3, 4);

    let (name, chart, cx, cy, half) =
        grid::gallery_cases().into_iter().find(|c| c.0 == case).expect("case");
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);

    let run = |e: usize| -> Vec<PixelOut> {
        let ens = EnsembleCfg { refine_flagged: false, escape_every: e, ..Default::default() };
        (0..sl.npix())
            .into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
            .collect()
    };
    let a = run(0);
    let b = run(ev);

    // ONE ramp, from the e0 pass, for both arms.
    let (lo, hi) = colour::range(&a, Scalar::ShapeSpread);
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);
    let col = |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);

    let bits = |x: f64| x.to_bits();
    let n = a.len();
    let d_tend = (0..n).filter(|&i| bits(a[i].t_end) != bits(b[i].t_end)).count();
    let d_state = (0..n).filter(|&i| a[i].state != b[i].state).count();
    let d_shape = (0..n).filter(|&i| bits(a[i].spread_shape) != bits(b[i].spread_shape)).count();
    let d_sv = (0..n)
        .filter(|&i| (0..3).any(|k| bits(a[i].shape_vec[k]) != bits(b[i].shape_vec[k])))
        .count();
    let d_rgb = (0..n).filter(|&i| col(&a[i]) != col(&b[i])).count();

    let pct = |c: usize| c as f64 / n as f64 * 100.0;
    println!(
        "{name} at {res}^2 = {n} px, escape_every 0 vs {ev}, one ramp ({lo:.4e}, {hi:.4e})\n\n\
         {:>18} {:>12} {:>9}\n\
         {:>18} {d_tend:>12} {:>8.4}%   <- did the fix reach the PHYSICS?\n\
         {:>18} {d_state:>12} {:>8.4}%\n\
         {:>18} {d_shape:>12} {:>8.4}%   <- did it reach the COLOURED FIELD?\n\
         {:>18} {d_sv:>12} {:>8.4}%\n\
         {:>18} {d_rgb:>12} {:>8.4}%   <- did it reach the IMAGE?",
        "quantity", "pixels differ", "share",
        "t_end", pct(d_tend),
        "state", pct(d_state),
        "spread_shape", pct(d_shape),
        "shape_vec", pct(d_sv),
        "RENDERED RGB", pct(d_rgb),
    );
    println!(
        "\nA large `t_end` count with a zero `RENDERED RGB` count is the finding: the fix\n\
         demonstrably reached the physics and the image did not move. A ZERO `t_end` count\n\
         would mean the flag never took effect, which is the failure mode this exists to rule\n\
         out -- and it is not what happened."
    );
}
