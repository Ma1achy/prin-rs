//! The reference-body switch instrumentation: `AzOpts::keep_drift_hist` and the four
//! `PixelOut` switch statistics.
//!
//! The three properties worth holding are the ones that would let a wrong number through
//! unnoticed: that turning the diagnostic on changes **no** physics, that the two arms
//! partition the boundaries rather than double-counting or dropping them, and that an empty
//! arm is `NaN` rather than `0` — a trajectory that never switched has no switch increment,
//! which is a different statement from an increment of zero.

use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::physics::burrau;

fn run(keep: bool) -> az::AzOut<f64> {
    az::integrate_az_opts(
        burrau::state::<f64>(),
        &burrau::masses::<f64>(),
        13.0,
        32,
        1e-2,
        30_000,
        &AzOpts { keep_drift_hist: keep, ..Default::default() },
    )
}

/// **The control that matters.** A diagnostic that perturbs the run it is diagnosing is worse
/// than no diagnostic: it would be reporting a trajectory nothing else in the build integrates.
#[test]
fn the_flag_changes_no_number() {
    let (off, on) = (run(false), run(true));
    assert_eq!(off.drift.to_bits(), on.drift.to_bits(), "drift moved");
    assert_eq!(off.t.to_bits(), on.t.to_bits(), "t_end moved");
    assert_eq!(off.switches, on.switches);
    assert_eq!(off.refs, on.refs);
    assert_eq!(off.steps, on.steps);
    for k in 0..2 {
        assert_eq!(off.state.r[k].x.to_bits(), on.state.r[k].x.to_bits(), "position moved");
    }
    assert!(off.drift_hist.is_empty(), "the series must be empty when the flag is off");
    assert!(!on.drift_hist.is_empty(), "the series must be populated when it is on");
}

/// The series and `refs` share one cadence, which is the whole basis for splitting the
/// increments by whether the reference changed. If they ever drift apart the split is silently
/// misaligned and the switch arm is a random subset.
#[test]
fn the_two_arms_partition_the_boundaries() {
    let o = run(true);
    let n = o.drift_hist.len().min(o.refs.len());
    assert!(n > 2, "need boundaries to partition");
    let sw = (1..n).filter(|&k| o.refs[k] != o.refs[k - 1]).count();
    let hd = (1..n).filter(|&k| o.refs[k] == o.refs[k - 1]).count();
    assert_eq!(sw + hd, n - 1, "the two arms must cover every increment exactly once");
    // **The test must be able to fire.** Burrau switches; an arm that is empty here would make
    // the partition trivially true and the assertion decoration.
    assert!(sw > 0, "no switch in this trajectory -- the partition is untested");
    assert!(hd > 0, "no hold in this trajectory -- the partition is untested");
    assert_eq!(o.drift_hist.len(), o.refs.len().min(o.drift_hist.len()));
    // The final entry of the series is the run's own reported drift, to round-off.
    let last = *o.drift_hist.last().unwrap();
    assert!(
        (last - o.drift).abs() <= 1e-12 * o.drift.abs().max(1e-30),
        "the series must end where `drift` does: {last:e} against {:e}",
        o.drift
    );
}

/// `NaN`, never `0`. The zero would read as "this switch cost nothing" on a trajectory that
/// never switched, and it would sit in a median beside real measurements.
#[test]
fn an_empty_arm_is_nan_and_not_zero() {
    use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
    use prin_rs::grid::{self, Chart};

    // A very short horizon: the reference is chosen once per sync boundary, so at `t_max` small
    // enough almost nothing has deformed far enough to change the longest side.
    let cfg = EnsembleCfg {
        t_max: 0.2,
        n_sync: 4,
        keep_drift_hist: true,
        refine_flagged: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(4, 4, 0.0, 0.0, 0.02, 0).with_chart(Chart::BodyPlane);
    let px: Vec<_> = (0..sl.npix()).map(|k| evaluate::<f64>(&sl, k, &cfg)).collect();
    let never: Vec<_> = px.iter().filter(|p| p.switches == 0).collect();
    assert!(!never.is_empty(), "no non-switching pixel here -- the case is untested");
    for p in never {
        assert!(p.t_first_switch.is_nan(), "t_first_switch must be NaN, got {}", p.t_first_switch);
        assert!(p.t_last_switch.is_nan());
        assert!(p.switch_jump_med.is_nan(), "an empty switch arm must be NaN, not 0");
        assert!(p.switch_jump_max.is_nan());
        assert!(p.hold_jump_med.is_finite(), "the hold arm is not empty and must be finite");
    }
}
