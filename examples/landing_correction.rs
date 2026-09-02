//! **The landing predictor, not the stepper, is what caps marched accuracy in this codebase.**
//!
//! `clamp_final_step` sizes the last step of each sync interval as `(dt_left - s.t)/(dt/ds)` with
//! `dt/ds` read *before* the step. That is a first-order predictor of the time increment, so the
//! step overshoots, and the clock is then clamped over the top of it: the state is at one time
//! and the clock says another. Measured on the figure-eight, `|s.t - dt_left|` runs `~5e-6` while
//! the GBS macro-step it follows is accurate to `1e-13`.
//!
//! `LhOpts::land_iterate` re-takes the landing step with a secant on `t(ds)`, using the step just
//! taken as the second point and `t(0) = 0` as the first. One extra step per correction, no new
//! machinery.
//!
//! # What it buys, and the control that says what it is
//!
//! ```text
//!   stepper  land    eta      closure    land resid    evals   corrections
//!   Kdk      false   8.0e-2   2.9396e-4    5.278e-6      446         0
//!   Kdk      true    8.0e-2   3.6399e-4    1.721e-15     516       100
//!   Rk4      false   5.0e-3   8.1543e-7    8.444e-8    25856         0
//!   Rk4      true    4.0e-2   2.7692e-8    1.832e-15    3608        84
//!   Gbs      false   8.0e-2   7.1387e-5    5.281e-6     8388         0
//!   Gbs      true    8.0e-2   4.1022e-8    1.971e-15   10108        99
//! ```
//!
//! The residual falls **nine orders**, to round-off. RK4 reaches `2.77e-8` for **3608**
//! evaluations where it previously needed 25856 to reach `8.15e-7` — twenty times better for
//! seven times less work. GBS reaches the same place at the **coarsest step tested**.
//!
//! **And KDK is unaffected — that is the control.** It is second order, so an `O(h^2)` landing
//! was never its binding constraint, and the correction changes nothing it could not already do
//! (it is marginally *worse*, paying for corrections it does not need). If `land_iterate` had
//! improved KDK too, it would be doing something other than removing an order-two cap, and the
//! account above would be wrong.
//!
//! # The floor this now runs into
//!
//! Both corrected arms stop at **`4.1022e-8` at every step size**. That is the figure-eight
//! fixture's own floor — its initial conditions and period carry nine significant digits — and
//! not the method. So this measures that the cap is *gone*; it does not measure how good GBS is,
//! and no order may be quoted from the corrected rows. Better initial conditions are what that
//! would need.
//!
//! # Not switched on anywhere
//!
//! `land_iterate` defaults **off** and is off in `pixel.rs`. Every committed number in `results/`
//! was taken without it, and AZ and Heggie have no landing correction at all — so enabling it for
//! logH alone would give one arm a landing the others lack, inside a comparison whose whole point
//! is that the arms differ in one named way. **Whether to port it to AZ and Heggie is a
//! corpus-invalidating decision and is not taken here.** What is on record is what it is worth:
//! their measured orders, 2.08 and 2.40, are what an `O(h^2)` landing allows, on a fixture where
//! an exact predictor reaches 4.52.
//!
//! Args: `n_sync`.

use prin_rs::integrate::logh::{integrate_lh, LhOpts, Stepper};
use prin_rs::physics::Cart;
use prin_rs::Vec2;

fn fig8() -> (Cart<f64>, [f64; 3], f64) {
    let (x, y) = (0.970_004_36, -0.243_087_53);
    let (vx, vy) = (-0.932_407_37, -0.864_731_46);
    (
        Cart::new(
            [Vec2::new(x, y), Vec2::zero(), Vec2::new(-x, -y)],
            [Vec2::new(-vx / 2.0, -vy / 2.0), Vec2::new(vx, vy), Vec2::new(-vx / 2.0, -vy / 2.0)],
        ),
        [1.0; 3],
        6.325_913_98,
    )
}

fn closure(a: &Cart<f64>, b: &Cart<f64>) -> f64 {
    (0..3).fold(0.0f64, |w, i| w.max((a.r[i] - b.r[i]).norm()).max((a.v[i] - b.v[i]).norm()))
}

fn main() {
    let n_sync: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(32);
    let (s0, m, t) = fig8();
    println!("figure-eight closure over {n_sync} sync intervals, predictive limit off\n");
    println!(
        "  {:>6} {:>7} {:>9} {:>11} {:>11} {:>9} {:>6} {:>10}",
        "step", "land", "eta", "closure", "land resid", "evals", "corr", "order"
    );
    for stepper in [Stepper::Kdk, Stepper::Rk4, Stepper::Gbs] {
        for land in [false, true] {
            let mut prev: Option<(f64, f64)> = None;
            for eta in [8e-2, 4e-2, 2e-2, 1e-2, 5e-3] {
                let o = integrate_lh(
                    s0, &m, t, n_sync, eta, 40_000_000,
                    &LhOpts {
                        stepper,
                        land_iterate: land,
                        r_coll_frac: 0.0,
                        stop_on_event: false,
                        step_limit_f: 0.0,
                        gbs_tol: 1e-13,
                        gbs_k_max: 6,
                        ..Default::default()
                    },
                );
                let e = closure(&o.state, &s0);
                let ord = prev.map(|(pe, px): (f64, f64)| (px / e).ln() / (pe / eta).ln());
                println!(
                    "  {:>6} {:>7} {eta:>9.1e} {e:>11.4e} {:>11.3e} {:>9} {:>6} {:>10}",
                    format!("{stepper:?}"), land, o.land_residual_max, o.force_evals, o.land_iters,
                    ord.map(|x| format!("{x:.2}")).unwrap_or_else(|| "-".into())
                );
                prev = Some((eta, e));
            }
            println!();
        }
    }
    println!(
        "HOW TO READ THIS\n\n\
         **`land resid` is the quantity, `closure` is the consequence.** The residual falls nine\n\
         orders to round-off; the closure then falls until it meets the FIXTURE's own ~4.10e-8\n\
         floor, which the corrected rows sit on at every step size. **No order may be read off a\n\
         corrected row** -- that flat 4.1022e-8 is nine-significant-digit initial conditions, not\n\
         a method.\n\n\
         **KDK is the control.** Second order, so an O(h^2) landing was never its constraint, and\n\
         the correction buys it nothing -- it is marginally worse, paying for corrections it does\n\
         not need. A fix that improved every arm would not be removing an order-two cap."
    );
}
