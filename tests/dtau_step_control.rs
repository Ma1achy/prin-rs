//! The `dtau` step-sizing modes: that they differ where the mechanism says they must, and agree
//! where it says they must not.
//!
//! **A test that cannot fail is indistinguishable from a test that passes.** "The modes produce
//! some number" would pass under a `dtau_mode` that was ignored entirely, which is exactly the
//! failure worth catching -- the field is threaded through four structs and a closure. So each
//! test below names the configuration in which it fires and asserts a *direction*, not a value.

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{integrate_az_opts, AzOpts, DtauMode};
use prin_rs::physics::{burrau, Cart};

const ETA: f64 = 1e-2;
const MAX_STEPS: usize = 30_000;

fn run(s: Cart<f64>, m: &[f64; 3], t_max: f64, n_sync: usize, mode: DtauMode)
    -> prin_rs::integrate::az::AzOut<f64>
{
    let o = AzOpts::<f64> {
        stop_on_event: false,
        stop_on_escape: false,
        dtau_mode: mode,
        ..Default::default()
    };
    integrate_az_opts(s, m, t_max, n_sync, ETA, MAX_STEPS, &o)
}

/// Burrau's `deep interior`: the region whose intervals routinely open on a close encounter, and
/// where `A*B` therefore grows by orders across an interval.
fn deep_interior(i: usize, j: usize, n: usize) -> (Cart<f64>, [f64; 3]) {
    let half = 0.05;
    let u = -half + 2.0 * half * (i as f64 + 0.5) / n as f64;
    let v = -half + 2.0 * half * (j as f64 + 0.5) / n as f64;
    let ic = grid::decode_state(&Chart::BodyPlane, 0, u, v);
    (ic.s, ic.m)
}

/// The modes are not the same algorithm, and this is where the difference lives.
///
/// `dt = A*B*dtau`. Sizing `dtau` once at entry pins the physical step to `eta*dt_left` only
/// while `A*B` stays near its entry value; where an interval opens at a close encounter and the
/// bodies then separate, `A*B` grows by orders and `dt` grows with it. So on a grid that
/// contains such trajectories the two modes must give materially different answers, and if they
/// do not, `dtau_mode` is not reaching the stepper.
#[test]
fn the_modes_disagree_where_ab_grows_across_an_interval() {
    let n = 8;
    let mut differing = 0usize;
    let mut worst = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let (s, m) = deep_interior(i, j, n);
            let a = run(s, &m, 13.0, 33, DtauMode::FixedPerInterval);
            let b = run(s, &m, 13.0, 33, DtauMode::PerStepInterval);
            let d: f64 = (0..3)
                .map(|k| (a.state.r[k] - b.state.r[k]).norm())
                .fold(0.0, f64::max);
            if d > 0.0 {
                differing += 1;
            }
            if d.is_finite() {
                worst = worst.max(d);
            }
        }
    }
    assert!(
        differing > n * n / 2,
        "only {differing} of {} deep-interior trajectories differ between the two modes; \
         `dtau_mode` is probably not reaching the stepper",
        n * n
    );
    assert!(worst > 1e-6, "the modes differ, but only at round-off (worst {worst:e})");
}

/// And they agree where `A*B` is flat -- which is what says the difference above is the
/// mechanism and not merely two arbitrary steppers.
///
/// A two-body-dominated Burrau trajectory over a **short** horizon never opens an interval at a
/// close encounter, so the recomputed `A*B` stays within a whisker of its entry value, the cap
/// binds throughout, and `PerStepInterval` degenerates to `FixedPerInterval`.
#[test]
fn the_modes_agree_to_roundoff_where_ab_is_flat() {
    let m = burrau::masses::<f64>();
    let s = grid::Slice::body_plane(1, 1, 1.0, 3.0, 0.0, 0).nominal::<f64>(0);
    let a = run(s, &m, 0.5, 2, DtauMode::FixedPerInterval);
    let b = run(s, &m, 0.5, 2, DtauMode::PerStepInterval);
    let d: f64 = (0..3).map(|k| (a.state.r[k] - b.state.r[k]).norm()).fold(0.0, f64::max);
    assert!(
        d < 1e-8,
        "over a short horizon with no encounter the two modes should be within round-off, got {d:e}"
    );
}

/// **Zeno by arithmetic**, asserted rather than argued.
///
/// `PerStepRemaining` sets `dt ~ eta*rem`, so `rem_{n+1} = rem_n (1 - eta)`: the interval is
/// approached geometrically and never completed. The tell is not the drift -- a stalled
/// trajectory has a *beautiful* drift, five orders better than either real mode -- it is that
/// `t` never advances. This is the test that stops that number being read as accuracy.
#[test]
fn per_step_remaining_stalls_rather_than_integrating() {
    let m = burrau::masses::<f64>();
    let s = grid::Slice::body_plane(1, 1, 1.0, 3.0, 0.0, 0).nominal::<f64>(0);
    let n_sync = 33;
    let o = run(s, &m, 13.0, n_sync, DtauMode::PerStepRemaining);
    assert!(o.budget_exhausted, "PerStepRemaining should exhaust its step budget");
    // It completes at most the first interval and then never satisfies `s.t >= dt_left`.
    assert!(
        o.t < 13.0 * 2.0 / n_sync as f64,
        "PerStepRemaining reached t = {} of 13; the geometric decay should stall it inside the \
         first interval",
        o.t
    );
    let ok = run(s, &m, 13.0, n_sync, DtauMode::PerStepInterval);
    assert!(!ok.budget_exhausted && (ok.t - 13.0).abs() < 1e-9, "the control must complete");
    // And the trap it sets: its drift looks *better* than the mode that actually ran.
    assert!(
        o.drift < ok.drift,
        "the stalled mode is expected to show a smaller drift than the completed one -- that is \
         the whole reason `t/t_max` is printed beside it"
    );
}

/// The floor is not the guard it is named for, and this records which one is doing the work.
///
/// `ab_min` is the raw `A*B` before the `T::TINY` clamp. At f64 the clamp is `1e-300` and the
/// measured minimum on `deep interior` sits around `1e-215` -- so the floor never binds and the
/// explicit `is_finite` test is what catches a degenerate state. At f32 `TINY*TINY` underflows
/// and `dtau` comes out `inf`, which the floor also does not catch.
#[test]
fn ab_min_is_recorded_and_the_f64_floor_never_binds() {
    let n = 6;
    let mut smallest = f64::INFINITY;
    let mut floored = 0usize;
    for j in 0..n {
        for i in 0..n {
            let (s, m) = deep_interior(i, j, n);
            let o = run(s, &m, 13.0, 33, DtauMode::PerStepInterval);
            if o.ab_min.is_finite() {
                smallest = smallest.min(o.ab_min);
            }
            floored += o.ab_floored as usize;
        }
    }
    assert!(smallest.is_finite() && smallest > 0.0, "ab_min was never recorded");
    assert!(
        smallest < 1e-3,
        "deep interior should drive A*B far below O(1); got {smallest:e}, which suggests ab_min \
         is being written from the entry state rather than per step"
    );
    assert_eq!(floored, 0, "the f64 TINY floor bound on {floored} trajectories, which it has not \
                            done before -- the report in RESULTS §23 needs updating");
}
