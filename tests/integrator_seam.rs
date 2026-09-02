//! The `Integrator` seam: does the flag do anything, does it declare itself, and is `Az` still
//! what it was?
//!
//! This project has shipped a mechanism **disabled** twice — `refine_flagged: false` copied into
//! every render harness, and `k_frac = 1.0` making the ranked frontier take the top 100%. Both
//! were invisible because nothing printed the setting and no test asserted the mechanism fired.
//! A new selectable integrator is the same hazard in a new place, so the guards come with it
//! rather than after it.

use prin_rs::ensemble::pixel::{self, EnsembleCfg};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::integrate::az::RefPolicy;
use prin_rs::integrate::Integrator;

fn slice() -> grid::Slice {
    // `deep interior` at `t = 13`: chaotic enough that the two integrators must disagree, and a
    // fixture the rest of the suite already uses. A tame region would let an inert flag pass.
    grid::region("deep interior", 8, 8, 0.05).unwrap()
}

fn cfg_at(i: Integrator) -> EnsembleCfg {
    EnsembleCfg::production().with_overrides(&[Override::Integrator(i)])
}

/// **`Heggie` is the default from 2026-09-02, and an `Az` run must declare itself.**
///
/// The guard this replaces said *"`Az` is the default and it must stay so: every committed number
/// in `results/` was taken under it, and a default that quietly changed the corpus is this
/// project's recorded failure."* That concern is real and is **not** answered by changing the
/// assertion — it is answered by the second arm below. The corpus does not become wrong when the
/// default moves; it becomes wrong when the corpus stops *saying which integrator it was*. So
/// this test asserts the property that keeps it honest: `overrides_vs_production` **derives** the
/// declaration by diffing, so every AZ artefact from here on carries `integrator=Az` in its
/// provenance sidecar, and the pre-2026-09-02 corpus is labelled by its own committed headers.
///
/// The reason for the default is on `Integrator::Heggie` and is the **cadence measurement**, not
/// the win count: a default with a scoreboard does not survive the next investigation.
#[test]
fn production_is_heggie_and_an_az_run_declares_itself() {
    assert_eq!(EnsembleCfg::production().integrator, Integrator::Heggie);
    assert!(EnsembleCfg::production().is_production());

    // The arm that makes the change safe. Without it this test is a restatement of the default.
    let az = cfg_at(Integrator::Az);
    assert!(!az.is_production(), "an AZ run reads as production, so it would declare nothing");
    let p = az.provenance();
    assert!(p.contains("integrator") && p.contains("Az"), "an AZ run does not name itself: {p}");
}

/// **The flag is not inert.** A selectable integrator that produces the same pixels as the old
/// one would pass every other test in this file.
#[test]
fn switching_the_integrator_changes_the_pixels() {
    let sl = slice();
    let (az, hg) = (cfg_at(Integrator::Az), cfg_at(Integrator::Heggie));
    let mut moved = 0usize;
    let mut worst = 0.0f64;
    for k in 0..sl.npix() {
        let a = pixel::evaluate::<f64>(&sl, k, &az);
        let h = pixel::evaluate::<f64>(&sl, k, &hg);
        let d = (0..3).fold(0.0f64, |w, i| w.max((a.shape_vec[i] - h.shape_vec[i]).abs()));
        if d > 0.0 {
            moved += 1;
        }
        worst = worst.max(d);
    }
    println!("deep interior 8x8: {moved} of {} pixels move, worst |d shape| = {worst:.3e}", sl.npix());
    assert!(moved > 0, "the integrator flag changed nothing — it is wired but inert");
}

/// **Absent is not zero.** Six of `MarchOut`'s fields have no Heggie analogue, and a plausible
/// zero would be read as a value: `ab_min = 0.0` says "the product hit its floor on every step",
/// which is the opposite of "there is no such product here".
#[test]
fn the_fields_heggie_does_not_have_read_as_absent() {
    let sl = slice();
    let cfg = EnsembleCfg::production().with_overrides(&[
        Override::Integrator(Integrator::Heggie),
        Override::KeepRefPath(true),
    ]);
    let p = pixel::evaluate::<f64>(&sl, 0, &cfg);
    println!("heggie: ab_min = {}, ref_path len = {}, n_retry = {}", p.ab_min, p.ref_path.len(), p.n_retry);
    assert!(
        !p.ab_min.is_finite(),
        "ab_min read {} under Heggie, which has no A*B — a finite value here is a fabricated one",
        p.ab_min
    );
    assert!(p.ref_path.is_empty(), "Heggie reported a reference-body path it does not have");
    assert_eq!(p.n_retry, 0);

    // ...and AZ still fills them, so the assertions above are about Heggie and not about the
    // fields being dead everywhere.
    //
    // **`Integrator::Az` is named explicitly and must stay named.** This arm relied on AZ being
    // the default, and when the default moved to Heggie it silently became Heggie-against-Heggie
    // — the control comparing the thing to itself. It failed loudly (`AZ's ab_min is non-finite`)
    // and that failure was correct: *the control arm caught that the control was no longer the
    // control*, which is now the fifth time on this project that a fixture moved under a test and
    // the control, not the property, is what noticed. Measured while diagnosing it: AZ fills
    // `ab_min` on **64 of 64** pixels of this fixture, so the arm has a subject.
    let a = EnsembleCfg::production()
        .with_overrides(&[Override::Integrator(Integrator::Az), Override::KeepRefPath(true)]);
    let q = pixel::evaluate::<f64>(&sl, 0, &a);
    assert!(q.ab_min.is_finite(), "AZ's ab_min is non-finite, so the Heggie assertion is vacuous");
    assert!(!q.ref_path.is_empty(), "AZ's ref_path is empty, so the Heggie assertion is vacuous");
}

/// `RefPolicy::Shared` is a **no-op** under Heggie, and that is a property of the method rather
/// than an unimplemented feature: there is no reference body to share.
///
/// The control arm is AZ, where the same switch must *change* something — otherwise this test
/// would pass on a fixture where sharing happens to be inert for both.
#[test]
fn the_shared_reference_policy_is_a_no_op_under_heggie() {
    let sl = slice();
    let mut hg_moved = 0usize;
    let mut az_moved = 0usize;
    for (i, moved) in [(Integrator::Heggie, &mut hg_moved), (Integrator::Az, &mut az_moved)] {
        let per = EnsembleCfg::production().with_overrides(&[
            Override::Integrator(i),
            Override::RefPolicy(RefPolicy::PerCopy),
        ]);
        let shared = EnsembleCfg::production().with_overrides(&[
            Override::Integrator(i),
            Override::RefPolicy(RefPolicy::Shared),
        ]);
        for k in 0..sl.npix() {
            let a = pixel::evaluate::<f64>(&sl, k, &per);
            let b = pixel::evaluate::<f64>(&sl, k, &shared);
            if a.shape_vec != b.shape_vec || a.spread_shape != b.spread_shape {
                *moved += 1;
            }
        }
    }
    println!("RefPolicy Shared vs PerCopy: heggie moves {hg_moved}, az moves {az_moved}");
    assert_eq!(hg_moved, 0, "the shared reference policy changed a Heggie result; it cannot");
    assert!(
        az_moved > 0,
        "the shared policy is inert for AZ too on this fixture, so the Heggie assertion says \
         nothing about Heggie"
    );
}

/// The setting declares itself. `overrides_vs_production` **derives** the declaration by diffing,
/// so a config says what it is however it was built — a hand-maintained list goes stale, which is
/// the same failure one level up.
///
/// **Both directions, and the roles swapped when the default moved.** The non-production value
/// must name itself and the production value must stay silent — a header line that is always true
/// carries no information. Asserting only the first half would pass on a `provenance` that named
/// the integrator unconditionally.
#[test]
fn the_integrator_appears_in_the_provenance() {
    // `Az` is now the override.
    let c = cfg_at(Integrator::Az);
    let p = c.provenance();
    println!("{p}");
    assert!(p.contains("integrator"), "the integrator is missing from the provenance: {p}");
    assert!(p.contains("Az"), "the provenance does not name the integrator: {p}");
    assert!(!c.is_production());

    // `Heggie` is production and declares nothing — reached by the default AND by an explicit
    // override to the same value, because *overriding to the production value is not an override*.
    assert!(!EnsembleCfg::production().provenance().contains("integrator"));
    assert!(!cfg_at(Integrator::Heggie).provenance().contains("integrator"));
    assert!(cfg_at(Integrator::Heggie).is_production());
}
