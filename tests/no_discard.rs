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
