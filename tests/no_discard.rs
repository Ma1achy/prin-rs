//! **A failed copy is a measurement outcome, not missing data.**
//!
//! `energy_drift_max` and `gamma_max` reduce over the `E+1` copies of a pixel. Both used to read
//! `if x.is_finite() { max }`, which **discarded** a copy that had diverged or been truncated and
//! left the pixel reporting a finite, healthy-looking maximum over the survivors.
//!
//! That is the no-discard rule broken, and it breaks **chaos-selectively**: a copy goes non-finite
//! because its integration was hard, integration is hard at a close encounter, and close
//! encounters are what this instrument exists to measure. The statistic was biased against exactly
//! the regions it is pointed at, in the direction that makes them look tame.
//!
//! # What makes this test able to fail
//!
//! A starved step budget is the lever: it drives `budget_exhausted`, which sets `finite = false`
//! while leaving a **perfectly finite drift at the point the march stopped**. So it exercises the
//! part of the fix that a plain `is_finite` check would still get wrong — measured on
//! `deep_interior`, 199 pixel-outs carried `budget_exhausted` while `nonfin` read 0.
//!
//! And the **control arm is the half that keeps it honest**: the same pixels at a generous budget
//! must come back finite. Without it, a reduction hard-wired to return `inf` would pass the
//! property arm exactly as well as the correct one.

use prin_rs::ensemble::pixel::{self, EnsembleCfg};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;

fn slice() -> grid::Slice {
    // `deep interior` at `t = 13` -- the region that actually exhausts budgets in production.
    grid::region("deep interior", 8, 8, 0.05).unwrap()
}

#[test]
fn a_truncated_copy_makes_the_pixel_undetermined_and_a_healthy_one_does_not() {
    let sl = slice();
    // Starved: enough steps to start, nowhere near enough to finish.
    let starved = EnsembleCfg::production().with_overrides(&[Override::MaxSteps(64)]);
    let healthy = EnsembleCfg::production();

    let (mut starved_exhausted, mut starved_undetermined) = (0usize, 0usize);
    let (mut healthy_finite, mut healthy_n) = (0usize, 0usize);

    for k in 0..sl.npix() {
        let s = pixel::evaluate::<f64>(&sl, k, &starved);
        if s.budget_exhausted {
            starved_exhausted += 1;
            if !s.energy_drift_max.is_finite() {
                starved_undetermined += 1;
            }
        }
        let h = pixel::evaluate::<f64>(&sl, k, &healthy);
        healthy_n += 1;
        if h.energy_drift_max.is_finite() {
            healthy_finite += 1;
        }
    }

    // The test has a subject: the starved arm really did exhaust budgets.
    assert!(
        starved_exhausted > 0,
        "NO SUBJECT: nothing was budget-exhausted, so this test asserts nothing. \
         Lower `MaxSteps` or pick a harder region."
    );
    // The property: every truncated pixel reports undetermined rather than a survivor's value.
    assert_eq!(
        starved_undetermined, starved_exhausted,
        "{} of {} budget-exhausted pixels reported a FINITE drift -- a discarded copy is being \
         silently dropped from the reduction",
        starved_exhausted - starved_undetermined,
        starved_exhausted
    );
    // The control: a reduction that always returned `inf` would pass the arm above.
    assert_eq!(
        healthy_finite, healthy_n,
        "CONTROL FAILED: {} of {} pixels at a generous budget are non-finite, so the property arm \
         above proves nothing",
        healthy_n - healthy_finite,
        healthy_n
    );
}

// ---------------------------------------------------------------------------------------------
// The same rule one layer up: `N x N` footprints reduced to one quad.
// ---------------------------------------------------------------------------------------------

/// `scheduler::reduce` had the identical defect, and **the fix above made it strictly worse.**
///
/// Both `error_ratio_max` and `worst_energy_drift` reduced with
/// `.filter(|x| x.is_finite()).fold(0.0, f64::max)`. Before the pixel-level fix a
/// budget-truncated pixel contributed a finite, healthy-looking drift and at least reached the
/// `max`; after it that pixel is `+inf`, the filter dropped it, and a quad whose footprints were
/// *all* undetermined folded to `0.0` — the best value the scale admits, from nothing.
///
/// # What makes this test able to fail
///
/// Three arms, and the first two are the ones a wrong implementation passes:
///
/// - **subject** — the starved quad really does contain undetermined footprints. Without it the
///   property arm is asserting over an empty set.
/// - **control** — the same quad at a generous budget reports a finite drift. A reduction
///   hard-wired to `inf` passes the property arm exactly as well as a correct one.
/// - **the old form, computed side by side** — the filtered reduction must give a *different*
///   answer on the same data. If it does not, this test cannot distinguish the fix from the bug
///   and is decoration.
#[test]
fn an_undetermined_footprint_makes_the_whole_quad_undetermined() {
    use prin_rs::scheduler;
    use prin_rs::spatial::HotRule;

    let sl = slice();
    let n = 8usize;
    let starved = EnsembleCfg::production().with_overrides(&[Override::MaxSteps(64)]);
    let healthy = EnsembleCfg::production();

    let px_s: Vec<_> = (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &starved)).collect();
    let px_h: Vec<_> = (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &healthy)).collect();

    let undetermined = px_s.iter().filter(|p| !p.energy_drift_max.is_finite()).count();
    assert!(
        undetermined > 0,
        "NO SUBJECT: the starved quad has no undetermined footprint, so nothing is asserted here."
    );

    let red_s = scheduler::reduce(&px_s, n, 1e-4, HotRule::AbsTau(1e-4), starved.t_max);
    let red_h = scheduler::reduce(&px_h, n, 1e-4, HotRule::AbsTau(1e-4), healthy.t_max);

    // The property: the quad says undetermined, it does not report the survivors' maximum.
    assert!(
        !red_s.worst_energy_drift.is_finite(),
        "{} of {} footprints are undetermined and the quad still reports a finite \
         worst_energy_drift = {:e} -- the reduction is discarding them",
        undetermined,
        px_s.len(),
        red_s.worst_energy_drift
    );

    // The control: a reduction hard-wired to `inf` would pass the arm above.
    assert!(
        red_h.worst_energy_drift.is_finite(),
        "CONTROL FAILED: the healthy quad reports {:e}, so the property arm proves nothing",
        red_h.worst_energy_drift
    );

    // The discriminator: the old, filtered form must disagree on this very data. Without this
    // the two implementations could be indistinguishable and the test would be decoration.
    let old_form = px_s
        .iter()
        .map(|p| p.energy_drift_max)
        .filter(|x| x.is_finite())
        .fold(0.0f64, f64::max);
    assert!(
        old_form.is_finite() && old_form != red_s.worst_energy_drift,
        "NOT DISCRIMINATING: the filtered form gives {:e} and the fixed form {:e} -- if these \
         agree, this test cannot tell the bug from the fix",
        old_form,
        red_s.worst_energy_drift
    );
}

/// `NaN` and `+inf` are different answers and the reduction must not collapse them.
///
/// `error_ratio` is `0/0` when `sigma_E(0) == 0` — a collapsed decode, or a family where the
/// statistic is structurally undefined. That is **not** evidence of damage, and a quad of nothing
/// but such footprints must read `NaN`, never `0.0`.
///
/// Note why the obvious minimal fix is not enough, which is the whole reason this arm exists:
/// `f64::max` already ignores `NaN`, so merely deleting the `.filter` would look correct while
/// still folding an all-`NaN` quad to `0.0` — maximum confidence from no information.
#[test]
fn a_quad_with_nothing_determinable_reads_nan_and_never_zero() {
    // Exercised through the same code path the reduction uses, on values it genuinely produces.
    let all_nan = [f64::NAN, f64::NAN, f64::NAN];
    let mixed = [1.0f64, f64::NAN, 3.0];
    let with_inf = [1.0f64, f64::INFINITY, 3.0];

    let naive = |v: &[f64]| v.iter().copied().fold(0.0f64, f64::max);
    assert_eq!(naive(&all_nan), 0.0, "premise: the naive fold really does return 0.0 here");

    // The REAL reduction, not a copy of it -- a test that re-implements its subject passes
    // whatever the subject does.
    let fold = |v: &[f64]| prin_rs::scheduler::max_no_discard(v.iter().copied());
    assert!(fold(&all_nan).is_nan(), "all-undefined must be NaN, not the best value on the scale");
    assert_eq!(fold(&mixed), 3.0, "NaN must not poison a quad that has determinable footprints");
    assert_eq!(fold(&with_inf), f64::INFINITY, "an undetermined footprint must propagate");
}
