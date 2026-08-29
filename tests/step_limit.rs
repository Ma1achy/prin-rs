//! The four step-control candidates: inert where they should be, active where they should be.
//!
//! `StepLimit::None` is the default and **every committed number in the corpus was taken under
//! it**, so the first property is that nothing moves. The second is that each mode *can* move
//! something — a limit that never engages passes an inertness test perfectly.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{DtauMode, StepLimit};

const MODES: [StepLimit; 4] =
    [StepLimit::Reject, StepLimit::Predictive, StepLimit::AbGrowth, StepLimit::Global];

/// A permissive parameter per mode: one that must not bind on anything.
fn permissive(m: StepLimit) -> f64 {
    match m {
        // A fraction of `d_min` a step may move: 1e30 accepts everything.
        StepLimit::Reject => 1e30,
        // A crossing-time fraction: 1e30 crossing times is no bound.
        StepLimit::Predictive => 1e30,
        // An `A*B` growth factor: 1e30x growth is never reached.
        StepLimit::AbGrowth => 1e30,
        // An `eta` multiplier: 1.0 is the identity.
        StepLimit::Global => 1.0,
        StepLimit::None => 0.0,
    }
}

/// And a strict one, chosen to bind hard rather than plausibly.
fn strict(m: StepLimit) -> f64 {
    match m {
        StepLimit::Reject => 1e-3,
        StepLimit::Predictive => 1e-3,
        StepLimit::AbGrowth => 1.000_001,
        StepLimit::Global => 0.25,
        StepLimit::None => 0.0,
    }
}

fn run(region: &str, limit: StepLimit, f: f64, n: usize) -> Vec<PixelOut> {
    run_mode(region, limit, f, n, EnsembleCfg::production().dtau_mode)
}

fn run_mode(
    region: &str,
    limit: StepLimit,
    f: f64,
    n: usize,
    dtau_mode: DtauMode,
) -> Vec<PixelOut> {
    let cfg =
        EnsembleCfg { step_limit: limit, step_limit_f: f, dtau_mode, ..Default::default() };
    let sl = grid::region(region, n, n, 0.05).unwrap().with_chart(Chart::BodyPlane);
    (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect()
}

fn bits(px: &[PixelOut]) -> Vec<u64> {
    px.iter().flat_map(|p| p.shape_vec.iter().map(|x| x.to_bits()).collect::<Vec<_>>()).collect()
}

/// **Nothing moves when nothing binds — bitwise.**
///
/// Not "close to": bitwise, because a limit that is supposed to be absent and instead perturbs
/// the last bit would change the whole corpus for no reason and read as noise.
#[test]
fn a_permissive_parameter_is_bitwise_inert() {
    let base = run("near-field", StepLimit::None, 0.0, 4);
    for m in MODES {
        let got = run("near-field", m, permissive(m), 4);
        assert_eq!(bits(&base), bits(&got), "{m:?} at a permissive parameter moved the field");
    }
}

/// **And each mode can bind.** A limit that never engages would pass the test above perfectly.
#[test]
fn a_strict_parameter_engages_every_mode() {
    let base = run("near-field", StepLimit::None, 0.0, 4);
    for m in MODES {
        // C is the exception, and the reason is the finding below rather than a weaker test.
        if m == StepLimit::AbGrowth {
            continue;
        }
        let got = run("near-field", m, strict(m), 4);
        assert_ne!(bits(&base), bits(&got), "{m:?} at a strict parameter changed nothing");
    }
}

/// **C IS ALREADY SHIPPED, UNDER ANOTHER NAME — `DtauMode::PerStepInterval` IS AN `A*B` GROWTH
/// CLAMP AT `C = 1`.**
///
/// The brief's formula, `dt = min(A*B, ab_entry*C) * dtau`, assumes `dtau` is fixed across the
/// interval — which is `FixedPerInterval`. Under the shipped `PerStepInterval`,
/// `dtau = eta*dt_left/(A*B)` is recomputed every step, so `dt ~ eta*dt_left` *however much `A*B`
/// grows*, and C's bound works out to `C` times the step the mode already chose. At `C = 1.000001`
/// it is strictly larger, so the `min` never binds and the mode is **bitwise inert**.
///
/// This is the subsumption question answered before the measurement rather than after: two
/// mechanisms doing one job, and the older one already does it. Held as a test so that if
/// `PerStepInterval` is ever changed, the thing that silently starts mattering announces itself.
#[test]
fn the_ab_growth_clamp_is_subsumed_by_per_step_interval_and_bites_under_fixed() {
    for mode in [DtauMode::PerStepInterval, DtauMode::FixedPerInterval] {
        let base = run_mode("near-field", StepLimit::None, 0.0, 4, mode);
        let got = run_mode("near-field", StepLimit::AbGrowth, 1.000_001, 4, mode);
        if mode == DtauMode::PerStepInterval {
            assert_eq!(bits(&base), bits(&got), "C must be inert where the mode already clamps");
        } else {
            assert_ne!(bits(&base), bits(&got), "C must bite where `dtau` is held fixed");
        }
    }
}

/// The retry counter is `Reject`'s and nobody else's, and `Reject` at a strict parameter must
/// actually retry — otherwise the mode is untested wherever it is measured.
#[test]
fn only_reject_retries_and_it_does_retry() {
    for m in MODES {
        let px = run("near-field", m, strict(m), 4);
        let retries: u64 = px.iter().map(|p| p.n_retry).sum();
        if m == StepLimit::Reject {
            assert!(retries > 0, "Reject at f = {} never retried", strict(m));
        } else {
            assert_eq!(retries, 0, "{m:?} reported retries");
        }
    }
}

/// **The tripwire reads zero on a healthy region**, which is what makes a nonzero reading mean
/// something. It counts a step that carried the interval clock past twice its own interval —
/// `dt > dt_left` is a bug, not a condition to handle.
#[test]
fn the_tripwire_is_silent_on_a_tame_region() {
    let px = run("near-field", StepLimit::None, 0.0, 4);
    assert_eq!(px.iter().map(|p| p.n_overshoot).sum::<u64>(), 0);
    assert!(!px.iter().any(|p| p.retry_exhausted));
}

/// And it is **conditioned on `clamp_final_step`**: with the clamp off, overshooting the boundary
/// is the expected behaviour of a named measurement axis, so the tripwire must not fire there.
/// An assert that fires on a deliberate mode is a broken assert.
#[test]
fn the_tripwire_does_not_fire_on_the_deliberate_overshoot_mode() {
    let cfg = EnsembleCfg { clamp_final_step: false, ..Default::default() };
    let sl = grid::region("near-field", 4, 4, 0.05).unwrap().with_chart(Chart::BodyPlane);
    let px: Vec<PixelOut> =
        (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
    assert_eq!(px.iter().map(|p| p.n_overshoot).sum::<u64>(), 0);
}
