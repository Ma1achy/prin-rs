//! **The production step budget, with the arm that says the fixture is exercised.**
//!
//! `EnsembleCfg::production().max_steps` was `30_000` from before the integrator default moved to
//! Heggie (`84830a1`, 2026-09-02). It does not bind on any of the eight Burrau regions -- bitwise
//! zero flagged at every rung of `results/step_budget/` -- which is why it survived a whole corpus
//! and every scheduler test. It binds on the **latent charts**, where the gallery lives.
//!
//! The measurement is in `results/step_budget/README.md`. What is asserted here is the operative
//! claim rather than a fitted constant: on a chart where the old value truncated footprints, the
//! shipped value **determines** them.
//!
//! **The control is the half with teeth.** A test that only asserted "the production budget
//! truncates nothing" would pass just as well on a chart where nothing was ever near the budget --
//! *a test that cannot fail is indistinguishable from a test that passes*. So the old-budget arm
//! must fire, and the assertion says so in its own message.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid::{self};

/// 16^2 on `latent_mixed_h3`: 30 of 256 footprints exhaust `30_000` and every one of them clears
/// at the shipped budget. Runs in about a second.
#[test]
fn the_shipped_budget_determines_footprints_the_old_one_truncated() {
    let (_, chart, cx, cy, half) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == "latent_mixed_h3")
        .expect("latent_mixed_h3 is in the gallery");
    let n = 16;
    let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
    let prod = EnsembleCfg::production();
    let old = EnsembleCfg { max_steps: 30_000, ..prod.clone() };

    let mut truncated_old = 0usize;
    let mut truncated_now = 0usize;
    let mut undetermined_old = 0usize;
    let mut undetermined_now = 0usize;
    for k in 0..sl.npix() {
        let a = evaluate::<f64>(&sl, k, &old);
        let b = evaluate::<f64>(&sl, k, &prod);
        if a.budget_exhausted {
            truncated_old += 1;
        }
        if b.budget_exhausted {
            truncated_now += 1;
        }
        if a.n_nonfinite > 0 {
            undetermined_old += 1;
        }
        if b.n_nonfinite > 0 {
            undetermined_now += 1;
        }
    }
    println!(
        "latent_mixed_h3 {n}^2: budget_exhausted {truncated_old} -> {truncated_now}, \
         undetermined {undetermined_old} -> {undetermined_now}  (30_000 -> {})",
        prod.max_steps
    );

    // The control: without this the assertion below is about a chart that never approached the
    // budget, and it would pass on any value at all.
    assert!(
        truncated_old > 0,
        "the old budget truncates nothing on this fixture, so this test has no subject -- \
         re-derive it against `results/step_budget/` rather than loosening the bar"
    );
    assert_eq!(truncated_now, 0, "the shipped budget still truncates {truncated_now} footprints");
    assert!(
        undetermined_now < undetermined_old,
        "the budget cleared {truncated_old} truncations without determining any footprint, \
         which means the flag is not what `n_nonfinite` was counting here"
    );
}
