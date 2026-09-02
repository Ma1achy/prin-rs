//! The saturation telemetry reaches the payload, and carries values that could not be defaults.
//!
//! `AzOut::ab_floored` and `ab_min` were written on every march since they were added and read by
//! **nothing** — `pixel.rs` never touched either, so they stopped one layer below `PixelOut` and
//! no render, dump, criterion or test could see the floor fire. A sticky bit that nothing reads
//! is indistinguishable from one that never fires.
//!
//! **What each assertion could fail on**, because a plumbing test that only checks a field exists
//! is decoration: an unplumbed `dt_max` is exactly `0.0`, an unplumbed `ab_min` is `INFINITY`, and
//! an unplumbed `n_cap_hits` is `0`. Every assertion below is chosen to be false in that state.

// **Every `EnsembleCfg` in this file pins `Integrator::Az`, and the pin is not a convenience.**
//
// The subject of this file is AZ machinery -- `A*B` and its `TINY` floor, `DtauMode`, the LC
// branch, the reference body, `StepLimit::Reject`. None of it exists under Heggie, which has no
// reference body, no `A*B` and three KS charts. These configs took the integrator from the
// default; when the default moved to `Heggie` on 2026-09-02 they ran a kernel with no such
// machinery, and **every one of them failed on its own guard or control arm rather than on a
// wrong number** -- *"an unplumbed ab_min is INFINITY"*, *"the cap never fired -- either
// unplumbed or the premise is wrong"*, *"Reject at a strict parameter changed nothing"*,
// *"the hold arm is not empty and must be finite"*. The arms written to prove each test had a
// subject are what announced that it no longer did.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::DtauMode;
use prin_rs::integrate::Integrator;

fn near_field(eta: f64, mode: DtauMode) -> Vec<PixelOut> {
    let cfg = EnsembleCfg { eta, dtau_mode: mode, integrator: Integrator::Az, ..Default::default() };
    let (chart, cx, cy, half) = (Chart::BodyPlane, 0.0, 0.0, 0.05);
    let sl = grid::Slice::body_plane(6, 6, cx, cy, half, 2).with_chart(chart);
    (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect()
}

#[test]
fn dt_max_is_populated_and_sits_near_the_nominal_step() {
    let cfg = EnsembleCfg { integrator: Integrator::Az, ..Default::default() };
    let px = near_field(cfg.eta, cfg.dtau_mode);
    let nominal = cfg.eta * cfg.t_max / cfg.n_sync as f64;
    for p in &px {
        assert!(p.dt_max.is_finite(), "unplumbed or non-finite");
        assert!(p.dt_max > 0.0, "an unplumbed dt_max is exactly 0.0");
        // Under the landing clamp no step may exceed the interval it is inside.
        assert!(
            p.dt_max <= cfg.t_max / cfg.n_sync as f64 * 1.000_001,
            "a step larger than its own sync interval: {} > {nominal}",
            p.dt_max
        );
    }
    // And it is not all one value -- a constant would mean it is reading the setting, not the run.
    let mut v: Vec<u64> = px.iter().map(|p| p.dt_max.to_bits()).collect();
    v.sort_unstable();
    v.dedup();
    assert!(v.len() > 1, "dt_max took a single value over the footprint set");
}

#[test]
fn ab_min_reaches_the_payload_finite_and_positive() {
    let px = near_field(0.01, DtauMode::PerStepInterval);
    for p in &px {
        assert!(p.ab_min.is_finite(), "an unplumbed ab_min is INFINITY");
        assert!(p.ab_min > 0.0);
    }
}

/// **The cap is routine, not exceptional**, and that is worth a test rather than a comment.
///
/// `capped` fires whenever `A*B` falls below its value at the interval's entry — which is what
/// happens every time bodies approach mid-interval. So a near-field footprint set must show it
/// firing under `PerStepInterval` and **never** under `FixedPerInterval`, where there is no cap.
/// Read a nonzero `n_cap_hits` as "the mode was refused the step it asked for", not as a fault.
#[test]
fn the_cap_counter_fires_under_per_step_and_never_under_fixed() {
    let per: u64 = near_field(0.01, DtauMode::PerStepInterval).iter().map(|p| p.n_cap_hits).sum();
    let fixed: u64 = near_field(0.01, DtauMode::FixedPerInterval).iter().map(|p| p.n_cap_hits).sum();
    assert!(per > 0, "the cap never fired -- either unplumbed or the premise is wrong");
    assert_eq!(fixed, 0, "FixedPerInterval has no cap to hit");
}

/// The healthy negative control: a tame region reports no saturation of any kind. Without it a
/// mask that fired everywhere would pass every test above.
#[test]
fn a_tame_region_reports_no_saturation() {
    let px = near_field(0.01, DtauMode::PerStepInterval);
    assert!(!px.iter().any(|p| p.ab_floored), "the TINY floor should not fire in near-field");
    assert!(!px.iter().any(|p| p.budget_exhausted));
}
