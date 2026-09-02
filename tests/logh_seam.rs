//! The four chartless occupants at the `Integrator` seam.
//!
//! `tests/integrator_seam.rs` holds the same guards for Heggie. They are repeated rather than
//! generalised because the *reasons* differ per occupant — Heggie's `ab_min` is absent because it
//! has no `A*B`, logH's because it has no chart at all — and a guard whose message no longer
//! names why it is asserting something is a guard nobody will trust when it fires.
//!
//! Two hazards this project has actually shipped: `refine_flagged: false` copied into every
//! render harness, and `k_frac = 1.0` making the ranked frontier take the top 100%. Both were
//! invisible because nothing printed the setting and no test asserted the mechanism fired.

use prin_rs::ensemble::pixel::{self, EnsembleCfg};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::integrate::az::RefPolicy;
use prin_rs::integrate::Integrator;

/// The four occupants added by the logH work, in the order the tables print them.
const NEW: [Integrator; 4] = [
    Integrator::LogHLeapfrog,
    Integrator::LogHRk4,
    Integrator::PlainLeapfrog,
    Integrator::PlainRk4,
];

fn slice() -> grid::Slice {
    // `deep interior` at `t = 13`, the same fixture the Heggie seam test uses: chaotic enough
    // that the integrators must disagree. A tame region would let an inert flag pass.
    grid::region("deep interior", 8, 8, 0.05).unwrap()
}

fn cfg_at(i: Integrator) -> EnsembleCfg {
    EnsembleCfg::production().with_overrides(&[Override::Integrator(i)])
}

fn eval(i: Integrator) -> Vec<pixel::PixelOut> {
    let sl = slice();
    let c = cfg_at(i);
    (0..sl.npix()).map(|k| pixel::evaluate::<f64>(&sl, k, &c)).collect()
}

/// The registry `principia_integrator_contract.md` describes, asserted in this codebase's terms.
///
/// That document is the **GLSL app's** contract and is not in this repo — its `substep_bucket`,
/// `N_sub`, `N_max` and descriptor bit 5 appear nowhere in `src/`. So the profile is built here
/// and pinned here, rather than borrowed from a spec this port does not implement.
///
/// The load-bearing row is `re_registers`: it is `true` for AZ alone, and **that single boolean
/// is the whole logH experiment**. If it ever reads `true` for a chartless occupant, the
/// falsification test has no subject.
#[test]
fn the_profile_registry_says_what_each_occupant_is() {
    println!("  {:>10}  {:>6} {:>8} {:>10} {:>6}", "occupant", "charts", "re-reg", "owns dt/ds", "evals");
    for i in [Integrator::Az, Integrator::Heggie].into_iter().chain(NEW) {
        let p = i.profile();
        println!(
            "  {:>10}  {:>6} {:>8} {:>10} {:>6}",
            i.name(), p.charts, p.re_registers, p.owns_time_mapping, p.evals_per_step
        );
        assert_eq!(Integrator::parse(i.name()), Some(i), "{} does not round-trip", i.name());
    }
    assert!(Integrator::Az.profile().re_registers, "AZ stopped re-registering");
    for i in NEW {
        let p = i.profile();
        assert_eq!(p.charts, 0, "{} claims a chart; it is algorithmic regularisation", i.name());
        assert!(!p.re_registers, "{} claims to re-register, which is the property under test", i.name());
    }
    // The controls are the only occupants that march in physical time.
    assert!(!Integrator::PlainLeapfrog.profile().owns_time_mapping);
    assert!(!Integrator::PlainRk4.profile().owns_time_mapping);
    assert!(Integrator::LogHLeapfrog.profile().owns_time_mapping);
}

/// **Force evaluations, not steps.** The relation is exact — nothing retries on these paths — so
/// it can be asserted rather than trusted, which is what makes the counter usable in a table
/// where one arm spends four evaluations per step and another spends one.
#[test]
fn the_evaluation_count_matches_each_occupants_profile() {
    for i in [Integrator::Az, Integrator::Heggie].into_iter().chain(NEW) {
        let px = eval(i);
        let per = i.profile().evals_per_step as u64;
        let (mut steps, mut evals) = (0u64, 0u64);
        for p in &px {
            assert_eq!(
                p.total_force_evals, p.total_substeps * per,
                "{}: evals {} != steps {} * {per}", i.name(), p.total_force_evals, p.total_substeps
            );
            steps += p.total_substeps;
            evals += p.total_force_evals;
        }
        println!("  {:>10}: steps {steps:>10}  evals {evals:>11}", i.name());
    }
}

/// Every new variant must move the pixels. A wired-but-inert flag is this project's recorded
/// failure mode, and here it would make the whole falsification test a null.
#[test]
fn every_new_occupant_changes_the_pixels() {
    let az = eval(Integrator::Az);
    for i in NEW {
        let o = eval(i);
        let moved = (0..az.len())
            .filter(|&k| {
                (0..3).any(|j| (az[k].shape_vec[j] - o[k].shape_vec[j]).abs() > 0.0)
            })
            .count();
        println!("  {:>10}: {moved} of {} pixels differ from AZ", i.name(), az.len());
        assert!(moved > 0, "{} changed nothing — it is wired but inert", i.name());
    }
}

/// **The two knobs behind the four occupants both reach the field.**
///
/// `LogH*` against `Plain*` is the time transformation; leapfrog against RK4 is the stepper.
/// If either pair agreed, two rows of every comparison table would be one row printed twice —
/// the failure that made three unrelated charts "agree" once before.
#[test]
fn the_transformation_and_the_stepper_are_both_live_through_the_seam() {
    let d = |a: &[pixel::PixelOut], b: &[pixel::PixelOut]| {
        (0..a.len())
            .filter(|&k| (0..3).any(|j| (a[k].shape_vec[j] - b[k].shape_vec[j]).abs() > 0.0))
            .count()
    };
    let lf = eval(Integrator::LogHLeapfrog);
    let rk = eval(Integrator::LogHRk4);
    let plf = eval(Integrator::PlainLeapfrog);
    let prk = eval(Integrator::PlainRk4);
    println!("  logH_lf vs plain_lf (transformation, KDK) : {}", d(&lf, &plf));
    println!("  logH_rk4 vs plain_rk4 (transformation, RK4): {}", d(&rk, &prk));
    println!("  logH_lf vs logH_rk4 (stepper)             : {}", d(&lf, &rk));
    println!("  plain_lf vs plain_rk4 (stepper, control)  : {}", d(&plf, &prk));
    assert!(d(&lf, &plf) > 0, "the time transformation is inert under KDK");
    assert!(d(&rk, &prk) > 0, "the time transformation is inert under RK4");
    assert!(d(&lf, &rk) > 0, "the stepper is inert under logH");
    assert!(d(&plf, &prk) > 0, "the stepper is inert under the control");
}

/// The AZ-only fields must read **absent**, never a plausible zero — and the AZ control arm is
/// what keeps each assertion from being vacuous.
#[test]
fn the_fields_a_chartless_integrator_does_not_have_read_as_absent() {
    let sl = slice();
    let q = pixel::evaluate::<f64>(
        &sl, 0,
        &EnsembleCfg::production()
            .with_overrides(&[Override::Integrator(Integrator::Az), Override::KeepRefPath(true)]),
    );
    assert!(q.ab_min.is_finite(), "AZ's ab_min is non-finite, so every assertion below is vacuous");
    assert!(!q.ref_path.is_empty(), "AZ's ref_path is empty, so every assertion below is vacuous");

    for i in NEW {
        let c = EnsembleCfg::production()
            .with_overrides(&[Override::Integrator(i), Override::KeepRefPath(true)]);
        let p = pixel::evaluate::<f64>(&sl, 0, &c);
        assert!(
            !p.ab_min.is_finite(),
            "{}: ab_min read {}. There is no chart here at all, so a finite value is fabricated — \
             and a 0.0 would read as 'the product hit its floor on every step'",
            i.name(), p.ab_min
        );
        assert!(p.ref_path.is_empty(), "{} reported a reference-body path it does not have", i.name());
        assert_eq!(p.n_retry, 0, "{} reported a retry; it has no reject-and-retry arm", i.name());
    }
}

/// `RefPolicy::Shared` shares the nominal copy's reference-body choices. There are none here, so
/// it must be **bitwise** inert — and the AZ arm is what says the fixture can tell the difference.
#[test]
fn the_shared_reference_policy_is_a_no_op_for_every_chartless_occupant() {
    let sl = slice();
    let moved = |i: Integrator| {
        let a = EnsembleCfg::production()
            .with_overrides(&[Override::Integrator(i), Override::RefPolicy(RefPolicy::PerCopy)]);
        let b = EnsembleCfg::production()
            .with_overrides(&[Override::Integrator(i), Override::RefPolicy(RefPolicy::Shared)]);
        (0..sl.npix())
            .filter(|&k| {
                let (x, y) =
                    (pixel::evaluate::<f64>(&sl, k, &a), pixel::evaluate::<f64>(&sl, k, &b));
                (0..3).any(|j| x.shape_vec[j] != y.shape_vec[j]) || x.spread_shape != y.spread_shape
            })
            .count()
    };
    let az = moved(Integrator::Az);
    println!("  az (control): {az} pixels move under RefPolicy::Shared");
    assert!(az > 0, "the shared policy is inert for AZ on this fixture, so the assertions below \
                     say nothing about anything");
    for i in NEW {
        let n = moved(i);
        println!("  {:>10}: {n}", i.name());
        assert_eq!(n, 0, "{} moved under the shared reference policy; it has no reference", i.name());
    }
}

/// Every occupant must name itself in the provenance line, and production must not name one at
/// all. `output::provenance_sidecar` puts this beside every panel — the blind spot that hid
/// `refine_flagged: false` for six days was PNGs carrying no settings.
#[test]
fn every_occupant_appears_in_the_provenance() {
    assert!(!EnsembleCfg::production().provenance().contains("integrator"));
    for i in NEW {
        let c = cfg_at(i);
        let p = c.provenance();
        println!("  {:>10}: {p}", i.name());
        assert!(p.contains("integrator"), "{}: integrator missing from provenance: {p}", i.name());
        assert!(!c.is_production(), "{} reads as production", i.name());
    }
}
