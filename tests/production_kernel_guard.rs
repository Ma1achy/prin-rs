//! **The guard on the guard.** `scheduler::assert_production_kernel` refuses a `results/` write
//! from a config whose integration kernel is not production's. Its `max_steps` arm was missing
//! until 2026-09-06, so a tree built at the superseded 30_000 budget would have been written into
//! the corpus with nothing saying so — and 15 of the 26 gallery trees moved on that field alone.
//!
//! *A guard that always fires passes just as easily as one that never does*, so both arms are
//! asserted: production must be **accepted** on the same path the departure is refused on, and
//! each departure is refused **one field at a time** from an otherwise-production config. A single
//! all-fields-wrong arm would pass on a guard that only checks `integrator`.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::integrate::Integrator;
use prin_rs::integrate::az::driver::{DtauMode, StepLimit};
use prin_rs::scheduler::assert_production_kernel;

fn refuses(mutate: impl FnOnce(&mut EnsembleCfg)) -> bool {
    let mut e = EnsembleCfg::production();
    mutate(&mut e);
    std::panic::catch_unwind(move || assert_production_kernel(&e, "results/sweep")).is_err()
}

#[test]
fn production_is_accepted_and_each_kernel_field_is_refused_on_its_own() {
    let p = EnsembleCfg::production();

    // The control. Without it every arm below could be satisfied by a guard that refuses
    // everything, which is the same failure as one that refuses nothing.
    assert_production_kernel(&p, "results/sweep");

    // And the path arm: a scratch root must be waved through whatever the kernel is, or the
    // "write it to a scratch root" the panic message offers would be a dead end.
    let mut old = p.clone();
    old.max_steps = 30_000;
    assert_production_kernel(&old, "/tmp/scratch/sweep");

    assert!(refuses(|e| e.max_steps = 30_000), "the superseded 30_000 step budget reached results/");
    assert!(refuses(|e| e.integrator = Integrator::Az));
    assert!(refuses(|e| e.step_limit = StepLimit::None));
    assert!(refuses(|e| e.step_limit_f *= 2.0));
    assert!(refuses(|e| e.dtau_mode = DtauMode::FixedPerInterval));
    assert!(refuses(|e| e.clamp_final_step = !e.clamp_final_step));

    // The named non-kernel fields stay waved through: they are experiment axes and arguments,
    // not the kernel, and a guard that refused them would make every legitimate harness illegal.
    assert!(!refuses(|e| e.refine_flagged = !e.refine_flagged));
    assert!(!refuses(|e| e.t_max *= 2.0));
}

/// The absolute kernel stamp is present unconditionally — including on a config with **no**
/// overrides, which is the case that printed nothing distinguishing before.
#[test]
fn the_provenance_line_carries_the_absolute_kernel_even_with_no_overrides() {
    let p = EnsembleCfg::production();
    let s = p.provenance();
    assert!(s.starts_with("production ["), "no override list, so the stamp must be the whole line: {s}");
    assert!(s.contains(&format!("max_steps={}", p.max_steps)), "{s}");
    assert!(s.contains("integrator="), "{s}");

    // And it must actually SEPARATE the two kernels — the defect was that it did not.
    let mut old = p.clone();
    old.max_steps = 30_000;
    assert_ne!(p.kernel_stamp(), old.kernel_stamp());
}
