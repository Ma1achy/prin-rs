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
//! # Closed: the quad that integrated nothing — and the stated mechanism was the smaller half
//!
//! The open item was that `scheduler.rs` filters `ensemble_spread` before the within-arm
//! quantiles, so an empty vector gives `quantile -> NaN` and `decide`'s `!(spread > tau)` sends
//! `NaN` to **`Decision::Keep`**: a quad where nothing could be integrated reporting *refinement
//! does not pay*.
//!
//! **That path is real and, at the shipped step control, essentially never taken.** Measured at
//! `N = 8` on `deep interior`, `near-field` and `far`, a quad with every footprint
//! budget-exhausted and **all 512 copies flagged unusable** reads `ensemble_spread` finite on
//! every one of them: a truncated state is a perfectly good number, it is simply not the number
//! the statistic claims. So the live failure is an **ordinary finite spread over a contaminated
//! sample**, and every field that existed before it reads clean:
//!
//! | quad | budget | nonfin copies | `spread_median` | `error_ratio_max` | `worst_drift` |
//! |---|---|---|---|---|---|
//! | near-field, production | 0/64 | 0/512 | 2.584e-3 | 1.0020 | 2.500e-6 |
//! | near-field, `MaxSteps=64` | 64/64 | 512/512 | **4.580e-4** | **1.0000** | inf |
//! | deep interior, `MaxSteps=64` | 64/64 | 512/512 | **1.221e-3** | **1.0000** | inf |
//!
//! `near-field`'s starved spread is **5.6x smaller** than its healthy one — it reads as *better*
//! resolved — while `deep interior`'s moves the other way, so it is not a bias that could be
//! corrected for. And `error_ratio_max` is **1.0000, exactly its converged value**, because every
//! copy stopped at the same early point and so agrees perfectly with the others: the statistic
//! whose job is to say *this pixel is not data* reports the ideal. Only `worst_energy_drift`
//! (`+inf`, from the fix above) and `n_nonfinite` (512/512) told the truth, and `decide` read
//! neither.
//!
//! Closed by `QuadReduction::n_undetermined` + `Decision::Undetermined`, keyed on
//! `scheduler::footprint_undetermined` — *not* on the spread being `NaN`, which would have been
//! a guard that could not fire. **`error_ratio`'s blindness is recorded and not repaired here**:
//! it is computed from the energy arrays and never consults the driver's usability flag, and
//! moving it moves every `error_ratio` in the corpus, which wants its own attribution.
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

// ---------------------------------------------------------------------------------------
// The quad that integrated nothing
// ---------------------------------------------------------------------------------------

/// Both arms of the guard, and the discriminating arm that says the defect was real.
///
/// The **property**: a quad none of whose footprints could be read stops as
/// `Decision::Undetermined`. The **control**: the same quad, integrated properly, does not — a
/// guard that always fires passes as easily as one that never does. The **discriminating arm**:
/// on the starved quad every pre-existing field reads clean, so nothing before `n_undetermined`
/// could have told the two apart.
#[test]
fn a_quad_that_integrated_nothing_stops_as_undetermined_and_a_healthy_one_does_not() {
    use prin_rs::quad::Decision;
    use prin_rs::render::Precision;
    use prin_rs::scheduler::{self, SchedCfg};

    let sched = SchedCfg { n: 8, budget: 24, tau_display: 1e-4, ..Default::default() };
    // `near-field` at t = 13, the region the scheduler tests are pinned to.
    let (cx, cy, half) = (1.0, 3.0, 0.05);

    let starved = EnsembleCfg::production().with_overrides(&[Override::MaxSteps(64)]);
    let healthy = EnsembleCfg::production();

    let (t_starved, _) = scheduler::descend(cx, cy, half, 0, &sched, &starved, Precision::F64);
    let (t_healthy, _) = scheduler::descend(cx, cy, half, 0, &sched, &healthy, Precision::F64);

    let computed = |t: &prin_rs::quad::QuadTree| -> Vec<prin_rs::quad::Quad> {
        t.nodes.iter().filter(|q| q.red.n_footprints > 0).cloned().collect()
    };
    let (cs, ch) = (computed(&t_starved), computed(&t_healthy));

    // The test has a subject: the starved arm really did fail to integrate.
    assert!(cs.len() >= 4 && ch.len() >= 4, "too few quads computed: {} / {}", cs.len(), ch.len());
    assert!(
        cs.iter().all(|q| q.red.n_undetermined == q.red.n_footprints),
        "NO SUBJECT: the starved arm produced readable footprints, so this test asserts nothing"
    );

    // The property.
    let n_undet = cs.iter().filter(|q| q.decision == Decision::Undetermined).count();
    assert!(
        n_undet > 0,
        "starved arm produced no Undetermined decision; stop reasons were {:?}",
        cs.iter().map(|q| q.decision.name()).collect::<Vec<_>>()
    );

    // The control: the healthy arm must reach none of it, and must reach it by having readable
    // footprints rather than by never being asked.
    assert!(
        ch.iter().all(|q| q.red.n_undetermined == 0),
        "the healthy arm has undetermined footprints, so the control is not clean"
    );
    assert!(
        ch.iter().all(|q| q.decision != Decision::Undetermined),
        "the healthy arm stopped as Undetermined; the guard fires on good data"
    );

    // The discriminating arm. On a starved quad every field that existed before this change
    // reads like ordinary data -- which is why the defect survived.
    let q = cs.iter().find(|q| q.decision == Decision::Undetermined).unwrap();
    assert!(
        q.red.spread_median.is_finite(),
        "the NaN route was taken, so this quad is not the case the fix is about"
    );
    assert!(
        (q.red.error_ratio_max - 1.0).abs() < 1e-6,
        "error_ratio_max reads {:.6}, not its converged 1.0 -- the table in the header is stale \
         and the discriminating arm needs re-measuring",
        q.red.error_ratio_max
    );

    // And it is NOT the other way of being undetermined: the decode is fine, the integration is
    // not. Without this the new variant could be a second name for `Collapsed`.
    assert!(
        !q.red.between_collapsed(),
        "the starved quad's decode collapsed too, so this does not separate the two causes"
    );
}

/// `ensemble_spread` swallows a `NaN` shape spread, and the guard must not depend on luck.
///
/// `ensemble_spread = sp_shape.max(sp_event)` and Rust's `f64::max` **ignores `NaN`**, so a
/// footprint whose shape spread is undetermined reports its *event* spread as an ordinary number.
/// Measured on `deep interior` under the pre-fix kernel: 11 footprints carry a `NaN`
/// `spread_shape` and **all 11** report a finite `ensemble_spread`.
///
/// On that corpus every one of them also carried an unusable copy, so
/// `footprint_undetermined`'s `n_nonfinite` arm caught them anyway. **Coincidence is not
/// coverage** — a triple collision reaches this with every copy still flagged usable, because
/// `shape_vec` is `NaN` at `I = 0` while the state stays finite. So the arm is asserted directly,
/// on a footprint constructed to have exactly that shape.
#[test]
fn a_nan_shape_spread_is_undetermined_even_though_ensemble_spread_swallows_it() {
    use prin_rs::ensemble::pixel::PixelOut;
    use prin_rs::scheduler::footprint_undetermined;

    // The premise, asserted rather than assumed: this is a property of `f64::max`, not of a
    // particular value, and if it ever stops holding the guard's second arm loses its subject.
    assert_eq!(f64::NAN.max(4.0), 4.0, "f64::max no longer ignores NaN; re-read this whole test");

    // A triple collision: shape spread undetermined, event spread ordinary, every copy usable.
    let mut p = PixelOut { n_nonfinite: 0, ..Default::default() };
    p.spread_shape = f64::NAN;
    p.spread_event = 4.0 / 7.0;
    p.ensemble_spread = p.spread_shape.max(p.spread_event);

    assert!(
        p.ensemble_spread.is_finite(),
        "NO SUBJECT: ensemble_spread is already non-finite, so the first arm covers this and the          test asserts nothing"
    );
    assert!(
        footprint_undetermined(&p),
        "a NaN shape spread reads as determined; the guard depends on n_nonfinite happening to fire"
    );

    // The control: the same footprint with a real shape spread is NOT undetermined. Without it a
    // predicate hard-wired to `true` would pass the arm above exactly as well.
    p.spread_shape = 1e-4;
    p.ensemble_spread = p.spread_shape.max(p.spread_event);
    assert!(!footprint_undetermined(&p), "the guard fires on a perfectly ordinary footprint");
}

/// `Collapsed` and `Undetermined` are independent, and `decide` reports the cause when both hold.
///
/// A decode collapse is the *cause* -- identical ICs -- and divergence is downstream of it, so
/// the more fundamental label wins. Asserted rather than left to the reading order of the source.
#[test]
fn a_collapsed_decode_outranks_a_failed_integration_when_both_hold() {
    use prin_rs::quad::{Decision, QuadReduction, QuadTree};
    use prin_rs::scheduler::{decide, SchedCfg};

    // `bootstrap_levels = 0` so the root is actually decided rather than split blind.
    let cfg = SchedCfg { n: 8, bootstrap_levels: 0, tau_display: 1e-4, ..Default::default() };
    let mut tree = QuadTree::new(1.0, 3.0, 0.05, cfg.n, 0);

    let base = QuadReduction {
        n_footprints: 64,
        n_distinct_ic: 64,
        n_undetermined: 0,
        spread_median: 1.0,
        ..Default::default()
    };

    // Neither: an ordinary quad reaches the policy branches.
    tree.nodes[0].red = base;
    assert!(!matches!(decide(&tree, 0, &cfg), Decision::Collapsed | Decision::Undetermined));

    // Integration failed, decode fine.
    tree.nodes[0].red = QuadReduction { n_undetermined: 64, ..base };
    assert_eq!(decide(&tree, 0, &cfg), Decision::Undetermined);

    // Decode collapsed, integration fine.
    tree.nodes[0].red = QuadReduction { n_distinct_ic: 63, ..base };
    assert_eq!(decide(&tree, 0, &cfg), Decision::Collapsed);

    // Both: the cause wins.
    tree.nodes[0].red = QuadReduction { n_distinct_ic: 63, n_undetermined: 64, ..base };
    assert_eq!(decide(&tree, 0, &cfg), Decision::Collapsed);

    // Partial is deliberately NOT a decision -- 63 of 64 footprints unreadable still reaches the
    // policy branches. If that ever becomes a threshold this assertion is what has to be argued
    // with, rather than the change landing silently.
    tree.nodes[0].red = QuadReduction { n_undetermined: 63, ..base };
    assert!(!matches!(decide(&tree, 0, &cfg), Decision::Collapsed | Decision::Undetermined));
}
