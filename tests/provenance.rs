//! One source of truth, and a record of every departure from it.
//!
//! `refine_flagged: false` propagated from experiment harnesses into render harnesses by copy
//! through five commits and six days. Nothing fired, because nothing recorded the choice. These
//! tests hold the mechanism that makes that impossible.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::ensemble::provenance::Override;

#[test]
fn default_is_production_and_says_so() {
    let c = EnsembleCfg::default();
    assert!(c.is_production());
    assert!(c.overrides_vs_production().is_empty());
    // Not an empty string. A blank field and an absent field look the same in a log.
    assert_eq!(c.provenance(), "production");
}

#[test]
fn a_named_override_is_recorded_with_the_value_it_replaced() {
    let c = EnsembleCfg::production().with_overrides(&[Override::RefineFlagged(false)]);
    assert!(!c.refine_flagged);
    let ov = c.overrides_vs_production();
    assert_eq!(ov.len(), 1, "exactly one departure, not all of them and not none");
    assert_eq!(ov[0].0, "refine_flagged");
    assert_eq!(ov[0].1, "false");
    assert_eq!(ov[0].2, "true", "production's value is carried, not merely the new one");
    assert!(c.provenance().contains("refine_flagged=false"));
}

/// **The case that actually matters.** The 111 existing struct literals were never going to be
/// rewritten in one go, and a mechanism that only declares configs someone remembered to annotate
/// would have missed every one of them — which is the original failure, one level up. The diff is
/// derived, so a legacy literal declares itself for free.
#[test]
fn a_plain_struct_literal_still_declares_itself() {
    let c = EnsembleCfg { refine_flagged: false, t_max: 50.0, ..Default::default() };
    let ov = c.overrides_vs_production();
    assert_eq!(ov.len(), 2);
    let names: Vec<&str> = ov.iter().map(|o| o.0).collect();
    assert!(names.contains(&"refine_flagged"));
    assert!(names.contains(&"t_max"));
}

/// A diff that reported *everything* would pass the tests above just as well as a correct one.
/// Walk one field at a time and assert each is detected **alone**.
#[test]
fn each_field_is_detected_and_only_that_field() {
    let cases: Vec<(&str, EnsembleCfg)> = vec![
        ("n_extra", EnsembleCfg { n_extra: 3, ..Default::default() }),
        ("eta", EnsembleCfg { eta: 0.005, ..Default::default() }),
        ("n_sync", EnsembleCfg { n_sync: 125, ..Default::default() }),
        ("max_steps", EnsembleCfg { max_steps: 20_000, ..Default::default() }),
        ("stop_on_escape", EnsembleCfg { stop_on_escape: true, ..Default::default() }),
        ("r_coll_frac", EnsembleCfg { r_coll_frac: 5e-3, ..Default::default() }),
        ("refine_max_passes", EnsembleCfg { refine_max_passes: 10, ..Default::default() }),
        ("keep_drift_hist", EnsembleCfg { keep_drift_hist: true, ..Default::default() }),
    ];
    for (name, c) in cases {
        let ov = c.overrides_vs_production();
        assert_eq!(ov.len(), 1, "{name}: expected exactly one departure, got {ov:?}");
        assert_eq!(ov[0].0, name);
    }
}

/// The named path and the literal path must agree, or the enum is a second source of truth.
#[test]
fn named_overrides_and_literals_agree() {
    let a = EnsembleCfg::production()
        .with_overrides(&[Override::RefineFlagged(false), Override::TMax(50.0)]);
    let b = EnsembleCfg { refine_flagged: false, t_max: 50.0, ..Default::default() };
    assert_eq!(a.provenance(), b.provenance());
}

/// The new step-control fields declare themselves like every other. Added when `StepLimit`
/// landed — and the exhaustive destructure in `overrides_vs_production` **failed the build**
/// until they were, which is the guard doing exactly what it exists for rather than a
/// convention being remembered.
#[test]
fn the_step_control_fields_are_declared() {
    use prin_rs::integrate::az::StepLimit;
    // `Predictive` at `f = 0.02` IS production, so overriding to it must report **nothing** --
    // and the first cut of this test asserted two overrides and failed the moment the default
    // changed. That is the mechanism working: a declaration is a difference from production, not
    // a restatement of it.
    let same = EnsembleCfg::production()
        .with_overrides(&[Override::StepLimit(StepLimit::Predictive), Override::StepLimitF(0.02)]);
    assert!(same.is_production(), "overriding to the production value is not an override");

    let c = EnsembleCfg::production()
        .with_overrides(&[Override::StepLimit(StepLimit::None), Override::StepLimitF(0.1)]);
    let ov = c.overrides_vs_production();
    assert_eq!(ov.len(), 2, "got {ov:?}");
    assert!(c.provenance().contains("step_limit=None"));
    assert!(c.provenance().contains("step_limit_f=0.1"));
    assert!(c.provenance().contains("production Predictive"), "the replaced value is carried");
}
