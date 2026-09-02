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

//! # The full audit: every `filter(is_finite)` in the tree, classified
//!
//! Fixing the instance that was pointed at is how this defect survived in the first place, so
//! all 26 sites were read. **Three classes, and only one is the bug** — a blanket repair would
//! have broken the other two.
//!
//! **A — silent discard in a reduction that feeds a number. THE BUG.**
//!
//! | site | old | all-unusable read as |
//! |---|---|---|
//! | `pixel.rs` `energy_drift_max`, `gamma_max` | `if is_finite { max }` | a survivor's finite max |
//! | `scheduler.rs` `error_ratio_max`, `worst_energy_drift` | `.filter(..).fold(0.0, max)` | **`0.0`** — perfectly clean |
//! | `pixel.rs` `dt_max` | `.filter(..).fold(0.0, max)` | **`0.0`** — "the largest step was zero" |
//! | `pixel.rs` `ab_min` | `.filter(..).fold(INFINITY, min)` | **`+inf`** — "never came close" |
//!
//! All four now report undetermined. `dt_max` is the sharpest: it is the diagnostic built to
//! catch a step of `2.209e128`, and it folded from `0.0`.
//!
//! **B — a ramp or axis window. CORRECT, and a "fix" would break it.** `render.rs:87-91,126`,
//! `colour.rs:530`, `prinq.rs:198`, `png.rs:184`, `plot.rs:87,279,423`. These set a colour or
//! plot range, and `colour::drift_rgb` paints the non-finite veto set **magenta separately**, so
//! an undetermined pixel renders as one rather than being mixed into the ramp. `range_q`'s own
//! doc gives the reason percentiles are used at all: one footprint at `1e12` would compress
//! every other pixel into the bottom of the range.
//!
//! **C — a guard that declines to compute, or a value whose non-finite case is meaningful.
//! CORRECT.** `scheduler.rs:726` requires all four children finite before reading an exponent;
//! `spatial.rs:196` returns an **all-hot** mask (undetermined) rather than an empty one;
//! `az/driver.rs:640`'s `+inf` *is* how "this constraint is not in force" is spelled;
//! `adaptive.rs:47` returns `None` below two points; `scheduler.rs:449`'s `first_divergence_t`
//! is `NaN` for "never diverged" and `frac_diverged` beside it carries that half explicitly,
//! counting non-crossers **in the denominator**.
//!
//! # Still open, and deliberately not fixed here because it moves every tree
//!
//! `scheduler.rs:378-380` filters `ensemble_spread` before the quantiles that feed `signal()`.
//! `quantile` returns `NaN` on empty, and `decide`'s `!(spread > tau)` sends `NaN` to
//! **`Decision::Keep`** — a quad where nothing could be integrated reports *refinement does not
//! pay*. `Decision::Collapsed` exists for undetermined quads but `between_collapsed()` tests
//! `n_distinct_ic < n_footprints`, a **decode** collapse; a quad with distinct ICs whose every
//! footprint diverged does not reach it. Two ways to be undetermined, one `Decision`. Changing
//! it is corpus-invalidating and wants its own measurement first.
//!

use prin_rs::ensemble::pixel::{self, EnsembleCfg};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::integrate::Integrator;

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

/// **The block-wide form of the property, and the reason it exists.**
///
/// The first version of this file checked `energy_drift_max` alone. `gamma_max` was fixed beside
/// it because it was in the same expression — and `dt_max` and `ab_min`, eight lines further down
/// the *same reduction*, were left filtering. They survived a fix, a full corpus regeneration and
/// a passing test suite. **A field-specific test is what let a block-wide defect through**, so
/// this arm asserts over every no-discard field at once and a new one has to be added here to
/// pass.
///
/// The two directions are not the same repair and the test knows it:
///
/// - `energy_drift_max`, `gamma_max`, `dt_max` are **max**-reductions — `+inf` is the absorbing
///   element and the pixel saturates to undetermined.
/// - `ab_min` is a **min** and takes `NaN`, because a min has no safe saturating value. That is
///   the same reasoning `d_min` is documented under, and `d_min` is deliberately left alone: it
///   is the quantity the collision labels are derived from, and `-inf` there would rewrite them.
#[test]
fn every_no_discard_field_in_the_reduction_reports_undetermined_together() {
    // **Two arms, because `ab_min` is AZ machinery and the other three are not.** Heggie has no
    // `A*B` at all, so its `ab_min` is non-finite on every pixel including healthy ones -- which
    // is correct and is asserted separately by `integrator_seam`. Pinning the whole test to `Az`
    // would have left the PRODUCTION integrator untested for a rule that applies to it, so the
    // Heggie arm runs over the three fields it actually has.
    for (integrator, fields) in [
        (Integrator::Az, &["energy_drift_max", "gamma_max", "dt_max", "ab_min"][..]),
        (Integrator::Heggie, &["energy_drift_max", "gamma_max", "dt_max"][..]),
    ] {
        let sl = slice();
        let starved = EnsembleCfg::production()
            .with_overrides(&[Override::MaxSteps(64), Override::Integrator(integrator)]);
        let healthy = EnsembleCfg::production()
            .with_overrides(&[Override::Integrator(integrator)]);

        let pick = |p: &prin_rs::ensemble::pixel::PixelOut, name: &str| match name {
            "energy_drift_max" => p.energy_drift_max,
            "gamma_max" => p.gamma_max,
            "dt_max" => p.dt_max,
            "ab_min" => p.ab_min,
            _ => unreachable!(),
        };

        let mut subject = 0usize;
        let mut bad: Vec<String> = Vec::new();
        let mut control_bad: Vec<String> = Vec::new();

        for k in 0..sl.npix() {
            let s = pixel::evaluate::<f64>(&sl, k, &starved);
            if s.budget_exhausted {
                subject += 1;
                for name in fields {
                    let v = pick(&s, name);
                    if v.is_finite() {
                        bad.push(format!("{name} = {v:e} on a truncated pixel"));
                    }
                }
                // The DIRECTION, not merely non-finiteness. `dt_max = 0.0` and `ab_min = +inf`
                // were the old readings and both are finite, so the arm above would have caught
                // them -- but `dt_max = NaN` would not be wrong and would not be the design.
                if s.dt_max != f64::INFINITY {
                    bad.push(format!("dt_max should saturate to +inf, got {:e}", s.dt_max));
                }
                if fields.contains(&"ab_min") && !s.ab_min.is_nan() {
                    bad.push(format!("ab_min should be NaN, got {:e}", s.ab_min));
                }
            }
            let h = pixel::evaluate::<f64>(&sl, k, &healthy);
            for name in fields {
                let v = pick(&h, name);
                if !v.is_finite() {
                    control_bad.push(format!("{name} = {v:e} at a generous budget"));
                }
            }
        }

        let who = integrator.name();
        assert!(subject > 0, "{who}: NO SUBJECT -- nothing was budget-exhausted.");
        assert!(
            bad.is_empty(),
            "{who}: {} findings across {subject} truncated pixels -- a pixel can contribute more than one, so this counts FINDINGS and not pixels:\n  {}",
            bad.len(),
            bad.join("\n  ")
        );
        assert!(
            control_bad.is_empty(),
            "{who}: CONTROL FAILED -- healthy pixels are non-finite, so the property arm proves nothing:\n  {}",
            control_bad.join("\n  ")
        );
    }
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
