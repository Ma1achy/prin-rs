//! The boundary-overshoot clamp: does it change the order, and is it separable from `dtau_mode`?
//!
//! The march exits each sync interval by **overshooting** (`s.t >= dt_left`) and the clock alone
//! was corrected -- the Cartesian state written back was the overshot one. A first-order error at
//! every boundary, inside an RK4 march.
//!
//! **The order is the test, not the error.** An error falls for many reasons; only the order says
//! the leading term changed. And every assertion here carries the arm that says it could have
//! failed: an order test that only checks the good arm passes equally well if the measurement is
//! broken.

use prin_rs::integrate::az::{integrate_az_opts, AzOpts, DtauMode};
use prin_rs::physics::Cart;
use prin_rs::Vec2;

/// Chenciner-Montgomery, equal masses, `G = 1`. Exactly periodic, so the distance between the
/// state at `T` and the state at `0` is an error measure with no reference trajectory to compute
/// and no chaotic amplification to contaminate it.
fn figure_eight() -> (Cart<f64>, [f64; 3], f64) {
    let (x, y) = (0.97000436, -0.24308753);
    let (vx, vy) = (-0.93240737, -0.86473146);
    (
        Cart::new(
            [Vec2::new(x, y), Vec2::new(0.0, 0.0), Vec2::new(-x, -y)],
            [
                Vec2::new(-vx / 2.0, -vy / 2.0),
                Vec2::new(vx, vy),
                Vec2::new(-vx / 2.0, -vy / 2.0),
            ],
        ),
        [1.0, 1.0, 1.0],
        6.325_913_98,
    )
}

fn closure(mode: DtauMode, clamp: bool, eta: f64) -> f64 {
    let (s0, m, period) = figure_eight();
    let o = AzOpts::<f64> {
        r_coll_frac: 0.0,
        stop_on_event: false,
        stop_on_escape: false,
        dtau_mode: mode,
        clamp_final_step: clamp,
        ..Default::default()
    };
    let out = integrate_az_opts(s0, &m, period, 32, eta, 200_000, &o);
    let mut d = 0.0f64;
    for i in 0..3 {
        d = d
            .max((out.state.r[i].x - s0.r[i].x).abs())
            .max((out.state.r[i].y - s0.r[i].y).abs())
            .max((out.state.v[i].x - s0.v[i].x).abs())
            .max((out.state.v[i].y - s0.v[i].y).abs());
    }
    d
}

/// The slope of `log(closure)` against `log(eta)` across a decade. Two-point estimates over a
/// factor of two are noisy enough to straddle any threshold; the endpoints are the measurement.
fn order(mode: DtauMode, clamp: bool) -> f64 {
    let (hi, lo) = (0.02, 0.001);
    (closure(mode, clamp, hi) / closure(mode, clamp, lo)).ln() / (hi / lo).ln()
}

#[test]
fn the_overshoot_is_first_order_and_the_clamp_removes_it() {
    // **The control arm.** Without it this is a test that the good number is good, which passes
    // just as well when the harness is broken. The overshoot must actually read first-order.
    let without = order(DtauMode::FixedPerInterval, false);
    assert!(
        without < 1.5,
        "the unclamped march must be first-order for this test to mean anything; got {without:.2}"
    );
    let with = order(DtauMode::FixedPerInterval, true);
    assert!(
        with > 2.5,
        "the clamp must raise the order well clear of first; got {with:.2} against {without:.2}"
    );
}

/// `PerStepInterval` is the *harder* arm and it is asserted separately, because the clamp's
/// interaction with a varying step is the whole reason the two changes ship together.
///
/// It lands at **~2**, not at `FixedPerInterval`'s ~3: the clamp sizes the last step from the
/// instantaneous `A*B`, which predicts the time increment to first order, so the landing residual
/// is `O(h^2)` per boundary -- and where that residual overshoots, the tolerance accepts it
/// rather than paying another step. Recorded as measured rather than asserted at the value hoped
/// for.
#[test]
fn the_clamp_raises_the_order_under_per_step_sizing_too() {
    let without = order(DtauMode::PerStepInterval, false);
    let with = order(DtauMode::PerStepInterval, true);
    assert!(without < 1.5, "control: unclamped must be first-order, got {without:.2}");
    assert!(with > 1.8, "clamped per-step order {with:.2} against unclamped {without:.2}");
    assert!(
        closure(DtauMode::PerStepInterval, true, 0.001)
            < closure(DtauMode::PerStepInterval, false, 0.001) / 1e3,
        "the clamp must buy orders of magnitude, not a factor"
    );
}

/// **The interaction claim, which is the reason neither change ships alone.**
///
/// Under the overshoot, switching `dtau_mode` moves the field a great deal -- the last step's
/// size becomes a function of local state, so neighbouring trajectories land at different times.
/// Under the clamp, both arms land on the boundary and the same switch moves it far less.
///
/// **The 10x here is a discrimination threshold on an 8x8 near-field grid, NOT the shipping
/// number.** At 1024^2 on `config_stability` the same ratio is **2.5x** (`RESULTS §24.6`): a
/// coarse grid over that window is dominated by the tame majority, while a million samples land
/// in the chaotic population where any step-control change diverges regardless. This test asks
/// whether the two knobs interact *at all*; it does not size the interaction, and no figure from
/// it should be quoted as if it did.
#[test]
fn the_two_knobs_are_not_independent() {
    use prin_rs::grid::{self, Chart};
    let (cx, cy, body) = grid::REGIONS
        .iter()
        .find(|r| r.0 == "near-field")
        .map(|r| (r.1, r.2, r.3))
        .expect("near-field");
    let n = 8usize;
    let half = 0.05;
    let ics: Vec<_> = (0..n * n)
        .map(|k| {
            let (i, j) = (k % n, k / n);
            let u = cx - half + 2.0 * half * (i as f64 + 0.5) / n as f64;
            let v = cy - half + 2.0 * half * (j as f64 + 0.5) / n as f64;
            let ic = grid::decode_state(&Chart::BodyPlane, body, u, v);
            (ic.s, ic.m)
        })
        .collect();
    let run = |mode, clamp| -> Vec<[f64; 3]> {
        let o = AzOpts::<f64> {
            r_coll_frac: 1e-3,
            stop_on_event: false,
            stop_on_escape: false,
            dtau_mode: mode,
            clamp_final_step: clamp,
            ..Default::default()
        };
        ics.iter()
            .map(|(s, m)| {
                let out = integrate_az_opts(*s, m, 13.0, 33, 1e-2, 30_000, &o);
                prin_rs::physics::shape::shape_vec(&out.state.r, m)
            })
            .collect()
    };
    let med = |a: &[[f64; 3]], b: &[[f64; 3]]| {
        let mut d: Vec<f64> = a
            .iter()
            .zip(b.iter())
            .map(|(p, q)| (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt())
            .filter(|x| x.is_finite())
            .collect();
        prin_rs::stats::quantile(&mut d, 0.5)
    };
    let ab = med(
        &run(DtauMode::FixedPerInterval, false),
        &run(DtauMode::PerStepInterval, false),
    );
    let cd = med(
        &run(DtauMode::FixedPerInterval, true),
        &run(DtauMode::PerStepInterval, true),
    );
    // The control: if switching `dtau_mode` moved nothing under the overshoot either, the
    // comparison below would be two zeros and would pass without measuring anything.
    assert!(ab > 1e-4, "control: the mode switch must move the field under the overshoot, got {ab:.3e}");
    assert!(
        cd < ab / 10.0,
        "under the clamp the same switch must move the field far less: {cd:.3e} against {ab:.3e}"
    );
}
