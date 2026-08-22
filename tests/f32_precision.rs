//! Step 6: where f32 and f64 are not the same algorithm, and what that costs.
//!
//! Acceptance thresholds here are **parameterised by precision**. `|dE/E| < 1e-12` is
//! structurally impossible at f32 — eps is ~1.19e-7 — so asserting it would fail for a
//! property of the type rather than of the port. f32 numbers are reported and gated at f32
//! tolerances; they are never asserted at f64 tolerance.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::physics::{burrau, Cart};
use prin_rs::{grid, Real, Vec2};

/// **The floor divergence, stated as measurement rather than as a comment.**
///
/// The reference's guard literals are chosen for f64 and two of them do not survive the cast.
/// This is the one place f32 and f64 run different algorithms, and it lives exactly where
/// degenerate geometry lives — which is why it is asserted rather than trusted.
#[test]
fn the_reference_floors_do_not_survive_the_cast_to_f32() {
    println!("{:>16}{:>16}{:>16}{:>16}", "constant", "f64", "as f32", "f32 in use");
    println!("{:>16}{:>16.3e}{:>16.3e}{:>16.3e}", "TINY", f64::TINY, 1e-300f64 as f32, f32::TINY);
    println!("{:>16}{:>16.3e}{:>16.3e}{:>16.3e}", "SYNC_EPS", f64::SYNC_EPS, 1e-15f64 as f32, f32::SYNC_EPS);
    println!("{:>16}{:>16.3e}{:>16.3e}{:>16.3e}", "DRIFT_FLOOR", f64::DRIFT_FLOOR, 1e-30f64 as f32, f32::DRIFT_FLOOR);
    println!("{:>16}{:>16.3e}{:>16.3e}{:>16.3e}", "DIST_FLOOR", f64::DIST_FLOOR, 1e-12f64 as f32, f32::DIST_FLOOR);
    println!();

    // 1. TINY: the reference's 1e-300 is exactly zero at f32, so the guard stops guarding.
    assert_eq!(1e-300f64 as f32, 0.0, "if this ever stops holding the guard is unnecessary");
    assert!(f32::TINY > 0.0 && f32::TINY.is_normal(), "the f32 floor must be a normal number");
    println!("1e-300 casts to exactly {} at f32 — a floor of zero floors nothing.", 1e-300f64 as f32);
    println!("f32 uses {:.3e}, which is normal and leaves headroom above the {:.3e} minimum.",
             f32::TINY, f32::MIN_POSITIVE);

    // But TINY squared does underflow, and A*B is a product of two floored quantities.
    let sq = f32::TINY * f32::TINY;
    println!("TINY*TINY at f32 = {sq:.3e}. A*B is a product of two floored quantities, so a");
    println!("doubly-degenerate state gives dtau = eta*dt_left/(A*B) = inf rather than a large");
    println!("finite step. That is caught by the explicit is_finite test in the RK4 loop, not");
    println!("by the floor — worth knowing which guard is actually doing the work.");
    assert_eq!(sq, 0.0, "documented behaviour: the squared floor underflows at f32");

    // 2. SYNC_EPS: 1e-15 is below the ulp of t ~ 13 at f32, so the slack is not slack.
    let t = 13.0f32;
    assert_eq!(t - (1e-15f64 as f32), t, "1e-15 is a no-op at f32 near t = 13");
    assert_ne!(t - f32::SYNC_EPS, t, "the f32 sync epsilon must actually subtract");
    let ulp = f32::from_bits(t.to_bits() + 1) - t;
    println!();
    println!("ulp(13) at f32 = {ulp:.3e}; the reference's 1e-15 slack is {:.0}x below it, so",
             ulp as f64 / 1e-15);
    println!("`t < t_target - 1e-15` degenerates to `t < t_target`. f32 uses {:.3e}.", f32::SYNC_EPS);
}

/// Gate (b) at both precisions, with the threshold set by the type.
///
/// BRIEF §5 asks for `d_min < 1e-10` and `|dE/E| < 1e-12`. The energy bound is five orders
/// below f32 epsilon, so it is not a statement about the port at all at that precision.
#[test]
fn the_two_body_collision_gate_is_parameterised_by_precision() {
    // Same parameters as the f64 gate in tests/two_body_collision.rs: eta = 1e-4, n_sync = 1.
    // d_min at a genuine collision is sampling-limited, so the threshold is a statement about
    // resolution as much as about correctness, and comparing precisions requires the same
    // resolution on both sides.
    fn run<T: Real>(name: &str, d_tol: f64, e_tol: f64) {
        let m = [T::one(); 3];
        let s = Cart::new(
            [
                Vec2::new(T::lit(-0.5), T::zero()),
                Vec2::new(T::lit(0.5), T::zero()),
                Vec2::new(T::lit(0.1), T::lit(1000.0)),
            ],
            [Vec2::zero(); 3],
        );
        let o = az::integrate_az(s, &m, T::lit(1.0), 1, T::lit(1e-4), 20_000_000, None);
        let (d, e) = (o.d_min_ref.to_f64().unwrap(), o.drift.to_f64().unwrap());
        println!("{name:>5}: d_min {d:.4e} (tol {d_tol:.0e})   |dE/E| {e:.4e} (tol {e_tol:.0e})   \
                  steps {}  finite {}", o.steps, o.finite);
        assert!(d < d_tol, "{name}: d_min {d:e} exceeds {d_tol:e}");
        assert!(e < e_tol, "{name}: drift {e:e} exceeds {e_tol:e}");
    }
    println!("BRIEF §5's gate (b), with the tolerance set by the type rather than by the spec:");
    run::<f64>("f64", 1e-10, 1e-12);
    // f32 eps is ~1.19e-7. The d_min bound is a *distance* and regularisation still carries
    // the trajectory through the collision, so it survives the cast; the energy bound cannot.
    run::<f32>("f32", 1e-9, 1e-4);
    println!();
    println!("The d_min half survives the cast outright: f32 measures 1.2218e-11 against f64's");
    println!("1.2881e-11, so it meets BRIEF §5's 1e-10 bound as written. It is gated here at");
    println!("1e-9 only to keep a decade of headroom on a sampling-limited quantity. The");
    println!("energy half cannot: 1e-12 is five orders below f32 eps, and f32 delivers");
    println!("2.8553e-6 - a statement about the type, not about the port.");
}

/// The whole pipeline at f32, with f64 alongside. Reported, and gated only on what f32 can
/// structurally deliver.
#[test]
fn f32_renders_the_grid_and_agrees_with_f64_on_the_labels() {
    let s = grid::region("near-field", 16, 16, 0.05).unwrap();
    let cfg = EnsembleCfg::default();
    let a: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect();
    let b: Vec<_> = (0..s.npix()).map(|i| evaluate::<f32>(&s, i, &cfg)).collect();

    let med = |v: &Vec<prin_rs::ensemble::pixel::PixelOut>, f: fn(&prin_rs::ensemble::pixel::PixelOut) -> f64| {
        let mut x: Vec<f64> = v.iter().map(f).filter(|q| q.is_finite()).collect();
        x.sort_by(|p, q| p.partial_cmp(q).unwrap());
        x[x.len() / 2]
    };
    let flips = a.iter().zip(b.iter()).filter(|(p, q)| p.outcome != q.outcome).count();
    println!("16x16 near-field, t=13, conditioned LC branch:");
    println!("  median |dE/E|      f64 {:.3e}   f32 {:.3e}",
             med(&a, |p| p.energy_drift_max), med(&b, |p| p.energy_drift_max));
    println!("  median spread_shape f64 {:.4e}  f32 {:.4e}",
             med(&a, |p| p.spread_shape), med(&b, |p| p.spread_shape));
    println!("  outcome label flips f32 vs f64: {flips} of {}", a.len());
    println!("  non-finite pixels   f64 {}  f32 {}",
             a.iter().filter(|p| p.n_nonfinite > 0).count(),
             b.iter().filter(|p| p.n_nonfinite > 0).count());

    // Gated on the label, not on the arithmetic: f32 must not change *what happened*, even
    // though it certainly changes the digits.
    let rel = (med(&b, |p| p.spread_shape) / med(&a, |p| p.spread_shape) - 1.0).abs();
    println!("  spread_shape median differs from f64 by {:.2}%", 100.0 * rel);
    assert!(rel < 0.1, "f32 spread_shape median is {rel:.3} off f64 — the branch cut is back");
    assert!(
        flips * 100 <= a.len(),
        "{flips} of {} outcome labels differ between precisions", a.len()
    );
}

/// The shared-reference flag governs **cross-copy sharing only**, never freezing across time.
/// Since Step 5b the nominal copy can terminate early, so its `refs` record is shorter than
/// `n_sync`; the shared policy must fall back rather than index past the end.
#[test]
fn the_shared_reference_policy_survives_an_early_terminating_nominal() {
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let mut short = 0usize;
    for i in 0..s.npix() {
        let nominal = az::integrate_az_opts(
            s.nominal::<f64>(i), &m, 13.0, 32, 0.01, 30_000,
            &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
        );
        if nominal.refs.len() < 32 {
            short += 1;
        }
        // The copy runs the full horizon even though the nominal stopped: the flag shares a
        // choice, it does not share a stopping time.
        let copy = az::integrate_az_opts(
            s.nominal::<f64>(i), &m, 13.0, 32, 0.01, 30_000,
            &AzOpts {
                forced_refs: Some(&nominal.refs),
                r_coll_frac: 0.0,
                stop_on_event: false,
                ..Default::default()
            },
        );
        assert_eq!(copy.refs.len(), 32, "pixel {i}: the copy did not run every boundary");
    }
    println!("{short} of {} pixels have a nominal copy that terminates before the horizon,",
             s.npix());
    println!("so its refs record is shorter than n_sync. The shared policy falls back to the");
    println!("per-copy choice past that point rather than indexing off the end.");
    assert!(short > 0, "no nominal terminated early, so this test proves nothing");
}
