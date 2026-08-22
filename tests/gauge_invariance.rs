//! BRIEF §5: rescale the initial conditions by `alpha` and time by `alpha^{3/2}`, and every
//! shape-space measurement must be identical.
//!
//! This is the check that catches an absolute length leaking into the kernel — the failure
//! that measured a 1.66x discrepancy in prior work, purely from an arbitrary choice of
//! overall size. Newtonian gravity has the exact symmetry `r -> alpha r`,
//! `t -> alpha^{3/2} t`, and the whole project is built on quotienting it out.

use prin_rs::ensemble::{jitter, stats};
use prin_rs::grid;
use prin_rs::integrate::az;
use prin_rs::physics::{burrau, shape, Cart};

const N_EXTRA: usize = 7;
const ETA: f64 = 0.01;

/// Spread statistics for one pixel's ensemble, with every length scaled by `alpha` and the
/// horizon by `alpha^{3/2}`.
fn measure(idx: usize, alpha: f64, t_max: f64) -> (f64, f64, f64) {
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let copies = jitter::copies::<f64>(&s, idx, N_EXTRA, 0.5, 0);
    let n_sync = 10;

    let shapes: Vec<[f64; 3]> = copies
        .iter()
        .map(|c| {
            let scaled = Cart::new(
                [c.r[0] * alpha, c.r[1] * alpha, c.r[2] * alpha],
                // Released from rest, so velocities are zero; the alpha^{-1/2} factor that
                // would apply to a moving configuration is written out anyway so the test
                // stays correct if the initial condition ever changes.
                [
                    c.v[0] / alpha.sqrt(),
                    c.v[1] / alpha.sqrt(),
                    c.v[2] / alpha.sqrt(),
                ],
            );
            let o = az::integrate_az(
                scaled,
                &m,
                t_max * alpha.powf(1.5),
                n_sync,
                ETA,
                30_000,
                None,
            );
            shape::shape_vec(&o.state.r, &m)
        })
        .collect();

    let classes: Vec<u8> = vec![0; shapes.len()];
    (
        shape::spread_shape(&shapes),
        shape::svar(&shapes),
        stats::spread_event::<f64>(&classes),
    )
}

#[test]
fn shape_spread_is_invariant_under_the_scale_symmetry() {
    let t_max = 4.0;
    println!("{:>8}{:>22}{:>22}{:>14}", "alpha", "spread_shape", "svar", "max |d| vs a=1");
    let (base_ss, base_sv, _) = measure(4, 1.0, t_max);
    println!("{:>8}{base_ss:>22.16e}{base_sv:>22.16e}{:>14}", 1.0, "-");

    // 0.25 and 4 are exact in binary, so they test scale-covariance with no rounding at all.
    // 3.7 and 1/3 are not, so they also exercise the roundoff path — otherwise a bitwise
    // result would partly be measuring binary exactness rather than the physics.
    let mut worst = 0.0f64;
    let mut worst_inexact = 0.0f64;
    for alpha in [0.25f64, 4.0, 3.7, 1.0 / 3.0] {
        let (ss, sv, _) = measure(4, alpha, t_max);
        let d = (ss - base_ss).abs().max((sv - base_sv).abs());
        println!("{alpha:>8.4}{ss:>22.16e}{sv:>22.16e}{d:>14.3e}");
        worst = worst.max(d);
        if alpha == 0.25 || alpha == 4.0 {
            assert_eq!(d, 0.0, "alpha {alpha} is exact in binary; the result must be bitwise identical");
        } else {
            worst_inexact = worst_inexact.max(d);
        }
    }
    println!();
    println!("BRIEF §5 asks for agreement to ~10 decimals. Powers of two come out bitwise");
    println!("identical; inexact factors land at {worst_inexact:.3e}, which is roundoff.");
    assert!(worst < 1e-10, "gauge invariance broken at {worst:e} — an absolute length has leaked in");
}

/// The same symmetry, applied across the whole grid rather than one pixel, so a leak that
/// only shows in particular geometries is not missed.
#[test]
fn shape_spread_is_invariant_across_the_grid() {
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();
    let mut worst = 0.0f64;
    let mut worst_pix = 0usize;
    for i in 0..s.npix() {
        let (b, _, _) = measure(i, 1.0, 2.0);
        for alpha in [0.25f64, 4.0] {
            let (v, _, _) = measure(i, alpha, 2.0);
            let d = (v - b).abs();
            if d > worst {
                worst = d;
                worst_pix = i;
            }
        }
    }
    println!("worst spread_shape deviation over the grid: {worst:.3e} at pixel {worst_pix}");
    assert!(worst < 1e-10, "gauge invariance broken at {worst:e}");
}
