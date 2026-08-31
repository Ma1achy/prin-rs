//! The logH march, against the instruments this project already has.
//!
//! Every fixture here is one AZ or Heggie has already been graded on, so the numbers are
//! comparable rather than merely present: the figure-eight convergence order (AZ 2.08, Heggie
//! 2.40), the three-pair radial collision (Heggie 5.422e-27, identical step count for all three),
//! the scale gauge, and Burrau against both other integrators.
//!
//! Two of these tests exist only to say the flags are not inert. *A difference can be small
//! because both sides are right or because one side is dead*, and a control arm that silently
//! stopped varying has cost this project a fixture twice.

use prin_rs::integrate::logh::driver::LhDsMode;
use prin_rs::integrate::logh::{integrate_lh, LhOpts, LhTime, Stepper};
use prin_rs::integrate::{az, heggie};
use prin_rs::physics::{burrau, energy, Cart};
use prin_rs::Vec2;

fn opts() -> LhOpts<f64> {
    LhOpts { r_coll_frac: 0.0, stop_on_event: false, ..Default::default() }
}

fn kdk() -> LhOpts<f64> {
    LhOpts { stepper: Stepper::Kdk, ..opts() }
}

// ---------------------------------------------------------------------------------------------
// The globality claim, and the flags not being inert
// ---------------------------------------------------------------------------------------------

fn collision_setup(pair: (usize, usize)) -> (Cart<f64>, [f64; 3]) {
    let third = 3 - pair.0 - pair.1;
    let mut r = [Vec2::zero(); 3];
    r[pair.0] = Vec2::new(-0.5, 0.0);
    r[pair.1] = Vec2::new(0.5, 0.0);
    r[third] = Vec2::new(0.1, 1000.0);
    let m = [1.0f64; 3];
    let c: Vec2<f64> = (r[0] + r[1] + r[2]) / 3.0;
    for x in r.iter_mut() {
        *x -= c;
    }
    (Cart::new(r, [Vec2::zero(); 3]), m)
}

/// **The globality claim, and the first hard evidence that the regularisation is the leapfrog.**
///
/// AZ regularises the two pairs sharing its reference body and leaves the third alone, so its
/// version of this test must first assert the colliding pair is one of the regularised ones.
/// Heggie has no third side. logH has no *side*: it never distinguishes a pair at all, so at
/// equal masses the three runs are one computation with the labels permuted and the step counts
/// must come out **identical**.
///
/// # KDK traverses it; RK4 does not, at any step size tried
///
/// ```text
///   KDK, limit off:  eta 1e-3   d_min 3.600e-8   drift 1.123e-9   steps    34034   finished
///                    eta 1e-4   d_min 1.724e-10  drift 3.746e-8   steps   339986   finished
///                    eta 3e-5   d_min 5.248e-11  drift 1.996e-7   steps  1133204   finished
///                    eta 1e-5   d_min 1.852e-12  drift 2.332e-6   steps  3399527   finished
///   RK4, limit off:  eta 1e-3 .. 1e-4            budget exhausted at 8e6-40e6 steps
/// ```
///
/// That is Mikkola & Merritt's *"the regularization is achieved by using the leapfrog"* as a
/// number, on the sharpest fixture this project has. **Under RK4 the time transformation alone
/// does not survive an exact collision**, which is the prediction the two-stepper design exists
/// to test, landing before any field is rendered.
///
/// # And logH does not meet BRIEF §5's collision gate, which is a difference of kind
///
/// The gate is `d_min < 1e-10` **with** `|dE/E| < 1e-12`. KDK reaches the first at `eta <= 3e-5`
/// and never the second: its drift *rises* with penetration depth, 1.1e-9 to 2.3e-6. Heggie
/// reaches `d_min = 5.422e-27` with `drift_reg = 4.4e-15`, flat under refinement.
///
/// The reason is structural rather than a tuning gap. Heggie's KS map removes the `1/r`
/// singularity from the Hamiltonian, so there is a regularised energy that stays at round-off
/// through the collision and only the Cartesian *readout* degrades. logH leaves the coordinates
/// alone and only slows the clock, so the encounter still has to be **resolved** rather than
/// removed, and there is no second energy that is better conditioned — `rho` tracks the drift
/// (2.4e-8 against 3.7e-8 at `eta = 1e-4`) instead of staying flat.
///
/// **A chartless method is not a coordinate regularisation with the coordinates left out.**
#[test]
fn radial_collision_is_traversed_by_the_leapfrog_and_not_by_rk4() {
    let mut steps = Vec::new();
    for pair in [(0usize, 1usize), (0, 2), (1, 2)] {
        let (s, m) = collision_setup(pair);
        let o = integrate_lh(s, &m, 1.0, 32, 3e-5, 8_000_000, &kdk());
        println!(
            "KDK pair {pair:?}: d_min = {:.3e}  |dE/E| = {:.3e}  rho = {:.2e}  steps = {}  \
             evals = {}",
            o.d_min, o.drift, o.gamma_max, o.steps, o.force_evals
        );
        assert!(o.finite, "KDK pair {pair:?} went non-finite");
        assert!(o.d_min < 1e-10, "KDK pair {pair:?}: d_min = {:e}", o.d_min);
        steps.push(o.steps);
    }
    assert!(
        steps[0] == steps[1] && steps[1] == steps[2],
        "the three pairs cost different step counts {steps:?}, so something is distinguishing a \
         pair in a method that has no pairs"
    );

    // The RK4 arm, pinned rather than skipped. If this ever starts completing, the two-stepper
    // framing needs revisiting -- so it fails loudly instead of quietly becoming true.
    let (s, m) = collision_setup((0, 1));
    let r = integrate_lh(s, &m, 1.0, 32, 3e-5, 2_000_000, &LhOpts { step_limit_f: 0.0, ..opts() });
    println!(
        "RK4 pair (0, 1): finite = {}  budget_exhausted = {}  d_min = {:.3e}  steps = {}",
        r.finite, r.budget_exhausted, r.d_min, r.steps
    );
    assert!(
        r.budget_exhausted,
        "RK4 completed the collision. That contradicts the measured behaviour this test records \
         and weakens the leapfrog-is-the-regularisation reading; re-measure before relaxing it."
    );
}

/// **Force evaluations are counted, not derived.**
///
/// RK4 spends four per step and KDK one. Nothing retries here, so the relation is exact and can
/// be asserted — which is what makes the counter trustworthy in the comparison harness, where
/// `steps` stops being commensurable the moment two steppers share a table.
#[test]
fn force_evaluations_are_counted_and_the_two_steppers_differ() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    for (stepper, per) in [(Stepper::Rk4, 4), (Stepper::Kdk, 1)] {
        let o = integrate_lh(s0, &m, 2.0, 16, 1e-3, 4_000_000, &LhOpts { stepper, ..opts() });
        println!("{stepper:?}: steps = {}  evals = {}  drift = {:.3e}", o.steps, o.force_evals, o.drift);
        assert_eq!(o.force_evals, o.steps * per, "{stepper:?}: evals != steps * {per}");
    }
}

/// The two things that would make this whole module a null: an inert `LhTime` and an inert
/// `Stepper`.
///
/// Neither assertion is about quality. They say only that the knobs reach the arithmetic, which
/// is the precondition for every comparison downstream.
#[test]
fn the_time_transformation_and_the_stepper_both_reach_the_arithmetic() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let run = |o: LhOpts<f64>| integrate_lh(s0, &m, 2.0, 16, 1e-3, 4_000_000, &o);

    let logh = run(opts());
    let plain = run(LhOpts { time: LhTime::None, ..opts() });
    let leap = run(kdk());

    let d = |a: &Cart<f64>, b: &Cart<f64>| {
        (0..3).fold(0.0f64, |w, i| w.max((a.r[i] - b.r[i]).norm()))
    };
    println!("Burrau to t = 2, |dr| against logH+RK4:");
    println!("  LhTime::None (unregularised) : {:.3e}   drift {:.3e}", d(&logh.state, &plain.state), plain.drift);
    println!("  Stepper::Kdk                 : {:.3e}   drift {:.3e}", d(&logh.state, &leap.state), leap.drift);
    println!("  logH + RK4                                        drift {:.3e}", logh.drift);
    assert!(d(&logh.state, &plain.state) > 0.0, "LhTime is inert: the transformation changed nothing");
    assert!(d(&logh.state, &leap.state) > 0.0, "Stepper is inert: KDK and RK4 agree bitwise");
}

// ---------------------------------------------------------------------------------------------
// Convergence, the gauge, and the clamp
// ---------------------------------------------------------------------------------------------

fn figure_eight() -> (Cart<f64>, [f64; 3], f64) {
    let (x, y) = (0.970_004_36, -0.243_087_53);
    let (vx, vy) = (-0.932_407_37, -0.864_731_46);
    (
        Cart::new(
            [Vec2::new(x, y), Vec2::zero(), Vec2::new(-x, -y)],
            [Vec2::new(-vx / 2.0, -vy / 2.0), Vec2::new(vx, vy), Vec2::new(-vx / 2.0, -vy / 2.0)],
        ),
        [1.0; 3],
        6.325_913_98,
    )
}

fn closure_err(a: &Cart<f64>, b: &Cart<f64>) -> f64 {
    (0..3).fold(0.0f64, |w, i| w.max((a.r[i] - b.r[i]).norm()).max((a.v[i] - b.v[i]).norm()))
}

/// **Read the order, not the error**, and read it against the unclamped control.
///
/// The figure-eight is exactly periodic, so `|state(T) - state(0)|` is a pure error with no
/// reference trajectory and no chaos. Measured orders on this same fixture: AZ 2.08 clamped
/// against 1.06 unclamped, Heggie 2.40 against 1.03. The unclamped arm is what makes the clamped
/// number mean anything — without it, second order could be the stepper rather than the landing.
///
/// # The ladder stops at `eta = 2.5e-3` because **the fixture has a floor at about 5e-8**
///
/// Measured, KDK clamped, per-rung orders across `eta` from `8e-2` down:
///
/// ```text
///   2.87  1.99  2.06  1.87  1.92   then   1.08  0.51
///   closure 2.94e-4 ... 1.75e-7    then   6.50e-8  4.55e-8
/// ```
///
/// The turn is not the integrator: KDK's energy drift keeps falling cleanly through it
/// (3.3e-10 to 1.3e-11), so what stops improving is the **position**, not the conservation. The
/// cause is the fixture — the Chenciner-Montgomery initial conditions and period here carry nine
/// significant digits, so a trajectory error of order `1e-8` is built into the closure before any
/// integrator runs.
///
/// **This bounds the AZ and Heggie numbers too.** Their ladders end at `eta = 1e-3`, where AZ
/// reads `9.2e-9` and Heggie `9.22e-9` — at the floor. So the 2.40-against-2.08 gap between them
/// is partly a reading of this constant, and the honest range for any endpoint slope on this
/// fixture stops above it. Extending the ICs is the fix and is not done here; stating the bound
/// is.
#[test]
fn convergence_order_on_the_figure_eight() {
    let (s0, m, t) = figure_eight();
    let etas = [8e-2, 4e-2, 2e-2, 1e-2, 5e-3, 2.5e-3];
    for stepper in [Stepper::Rk4, Stepper::Kdk] {
        for clamp in [true, false] {
            let mut errs = Vec::new();
            for &eta in &etas {
                let o = integrate_lh(
                    s0, &m, t, 32, eta, 40_000_000,
                    // **The predictive limit is OFF here**, as it is in the Heggie and AZ
                    // versions of this test. It resizes steps from local geometry, so leaving it
                    // on would put a third thing in a measurement of the stepper and the landing.
                    &LhOpts { stepper, clamp_final_step: clamp, step_limit_f: 0.0, ..opts() },
                );
                assert!(o.finite, "{stepper:?} went non-finite at eta = {eta}");
                errs.push((eta, closure_err(&o.state, &s0), o.steps, o.force_evals));
            }
            // The slope is quoted over the ASYMPTOTIC WINDOW, rungs 2..=5, not end to end.
            // Above it the coarsest rung is pre-asymptotic — RK4 moves only 1.32x from `8e-2` to
            // `4e-2` — and below it the fixture's own `~5e-8` floor takes over. An endpoint slope
            // across either reads the boundary rather than the method: end-to-end here gives RK4
            // 1.71 and KDK 2.14, and over the window they are 2.04 and 1.95.
            const LO: usize = 2;
            let hi = errs.len() - 1;
            let order = (errs[LO].1.ln() - errs[hi].1.ln()) / (etas[LO].ln() - etas[hi].ln());
            println!("{stepper:?}, clamp_final_step = {clamp}:");
            for (i, (eta, e, st, ev)) in errs.iter().enumerate() {
                let mark = if (LO..=hi).contains(&i) { "*" } else { " " };
                println!("  {mark}eta {eta:.1e}   closure {e:.4e}   steps {st:>8}   evals {ev:>9}");
            }
            println!("   order over the * window = {order:.2}");
            if clamp {
                assert!(order > 1.7, "{stepper:?} clamped came out at order {order:.2}");
            } else {
                assert!(
                    order < 1.5,
                    "{stepper:?} UNCLAMPED came out at order {order:.2}, so the clamp is not what \
                     the clamped arm is measuring and that number means nothing"
                );
            }
        }
    }
}

/// The gauge, asserted the same way it is for AZ and Heggie.
///
/// Under `r -> alpha r`, `t -> alpha^{3/2} t`: `K + B ~ alpha^{-1}` and `dt_left ~ alpha^{3/2}`,
/// so `ds ~ alpha^{1/2}` and both halves of the step come out covariant. An **absolute**
/// tolerance anywhere in the landing test would break this — that bug is on record, caught at
/// `4.24e-15` by the bitwise version of the AZ gauge test, and `Real::LAND_EPS_REL` is the fix.
#[test]
fn the_march_respects_the_scale_gauge() {
    let (s0, m, t) = figure_eight();
    for stepper in [Stepper::Rk4, Stepper::Kdk] {
        let base = integrate_lh(s0, &m, t, 32, 5e-3, 4_000_000, &LhOpts { stepper, ..opts() });
        for alpha in [0.25f64, 4.0] {
            let scaled = Cart::new(
                std::array::from_fn(|i| s0.r[i] * alpha),
                std::array::from_fn(|i| s0.v[i] / alpha.sqrt()),
            );
            let o = integrate_lh(
                scaled, &m, t * alpha.powf(1.5), 32, 5e-3, 4_000_000,
                &LhOpts { stepper, ..opts() },
            );
            let mut worst = 0.0f64;
            for i in 0..3 {
                worst = worst
                    .max((o.state.r[i] / alpha - base.state.r[i]).norm())
                    .max((o.state.v[i] * alpha.sqrt() - base.state.v[i]).norm());
            }
            println!("{stepper:?}  alpha = {alpha}: max rescaled difference = {worst:.3e}  \
                      (steps {} against {})", o.steps, base.steps);
            assert_eq!(o.steps, base.steps, "{stepper:?} at alpha = {alpha}: step counts differ, \
                       so something in the sizing carries a scale");
            assert!(worst < 1e-11, "{stepper:?} at alpha = {alpha}: {worst:e}");
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Against the other two integrators
// ---------------------------------------------------------------------------------------------

/// **Three independently derived methods on one trajectory.** The strongest correctness gate
/// available here: AZ regularises around a reference body, Heggie regularises three vectors
/// symmetrically, logH does not transform coordinates at all, and at tight `eta` all three must
/// describe the same Burrau trajectory.
#[test]
fn burrau_agrees_with_az_and_heggie() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let (t, n, eta) = (6.0, 64, 1e-4);

    let a = az::integrate_az(s0, &m, t, n, eta, 40_000_000, None);
    let h = heggie::integrate_hg(
        s0, &m, t, n, eta, 40_000_000,
        &heggie::HgOpts { step_limit_f: 0.0, r_coll_frac: 0.0, stop_on_event: false, ..Default::default() },
    );
    let l = integrate_lh(s0, &m, t, n, eta, 40_000_000, &LhOpts { step_limit_f: 0.0, ..opts() });
    let lk = integrate_lh(s0, &m, t, n, eta * 0.25, 40_000_000, &LhOpts { step_limit_f: 0.0, ..kdk() });

    let d = |x: &Cart<f64>, y: &Cart<f64>| {
        (0..3).fold(0.0f64, |w, i| w.max((x.r[i] - y.r[i]).norm()))
    };
    println!("Burrau to t = {t} at eta = {eta:e}, max |dr| between methods:");
    println!("  AZ                 drift {:.3e}  steps {:>9}  evals {:>10}", a.drift, a.steps, a.steps * 4);
    println!("  Heggie             drift {:.3e}  steps {:>9}  evals {:>10}   vs AZ {:.3e}", h.drift, h.steps, h.steps * 4, d(&a.state, &h.state));
    println!("  logH + RK4         drift {:.3e}  steps {:>9}  evals {:>10}   vs AZ {:.3e}", l.drift, l.steps, l.force_evals, d(&a.state, &l.state));
    println!("  logH + KDK, eta/4  drift {:.3e}  steps {:>9}  evals {:>10}   vs AZ {:.3e}", lk.drift, lk.steps, lk.force_evals, d(&a.state, &lk.state));
    println!("  logH rho: RK4 {:.3e}   KDK {:.3e}", l.gamma_max, lk.gamma_max);
    assert!(l.finite && lk.finite);
    assert!(d(&a.state, &l.state) < 1e-8, "logH+RK4 disagrees with AZ: {:e}", d(&a.state, &l.state));
    assert!(d(&a.state, &lk.state) < 1e-6, "logH+KDK disagrees with AZ: {:e}", d(&a.state, &lk.state));
}

#[test]
fn burrau_constants_are_unchanged() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let mt: f64 = m.iter().sum();
    let r = energy::hyperradius(&s0.r, &m);
    let e = energy::energy(&s0.r, &s0.v, &m, 0.0);
    println!("M = {mt}   R = {r:.4}   E = {e:.4}");
    assert!((mt - 12.0).abs() < 1e-12);
    assert!((r - 2.2361).abs() < 1e-4);
    assert!((e + 12.8167).abs() < 1e-4);
}

/// Trivially exact, and **that is why it discriminates nothing**.
///
/// AZ is label-covariant only because it re-chooses its reference body: frozen at index 0 it
/// reads `3.41e-6` against `3.23e-15` free, a factor of 1.06e9. logH has no labels to be
/// sensitive to, so it cannot separate those two readings. Kept as a cheap positive control that
/// the march is symmetric, and labelled as such rather than quoted as a win.
#[test]
fn the_march_is_independent_of_the_body_labelling() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let base = integrate_lh(s0, &m, 4.0, 32, 1e-4, 40_000_000, &opts());
    let mut worst = 0.0f64;
    for p in [[0usize, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
        let sp = Cart::new(
            std::array::from_fn(|i| s0.r[p[i]]),
            std::array::from_fn(|i| s0.v[p[i]]),
        );
        let mp: [f64; 3] = std::array::from_fn(|i| m[p[i]]);
        let o = integrate_lh(sp, &mp, 4.0, 32, 1e-4, 40_000_000, &opts());
        for i in 0..3 {
            worst = worst.max((o.state.r[i] - base.state.r[p[i]]).norm());
        }
    }
    println!("Burrau to t = 4 over all five non-identity label permutations: {worst:.3e}");
    println!(
        "  Not exactly zero, and it cannot be: `PAIRS` is a fixed index order, so permuting the\n           labels permutes the summation order inside `accel`, `kinetic` and `potential_pos`. What\n           is asserted is round-off over four time units of Burrau, not identity."
    );
    assert!(worst < 1e-11, "logH moved under relabelling: {worst:e}");
}

// ---------------------------------------------------------------------------------------------
// The step limit, and the mode that is an axis rather than a candidate
// ---------------------------------------------------------------------------------------------

/// **The predictive limit was predicted inert. It is not — it inverts.**
///
/// Prediction: `dt = ds/(K+B)` and `K+B` is `U` on shell, which diverges at a close approach, so
/// the physical step already shrinks with the separation and the limit should buy nothing.
///
/// Measured, it does two opposite things:
///
/// - On **Burrau at `t = 13`** it is a large *improvement* — drift `2.3e-6 -> 1.2e-9` under RK4
///   for 36% more steps.
/// - On the **two-body radial collision** it is fatal: with `f = 0.02` neither stepper finishes,
///   both burning the whole budget, where KDK with the limit off finishes in 34034 steps. The
///   bound is `ds <= f d_min U/|v_rel|`, which in free fall tends to `f sqrt(d)` while the
///   unbounded step wants to grow as `1/d` — *do not shrink the fictitious step at close
///   approach*, at a third site.
///
/// So `step_limit_f` defaults to `0.0` here and to `0.02` for AZ and Heggie. Both arms are
/// carried and the comparison harness runs both, because a knob held fixed for fairness is a
/// knob whose effect is unattributed.
#[test]
fn the_predictive_limit_helps_on_burrau_and_is_fatal_at_a_collision() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    for stepper in [Stepper::Rk4, Stepper::Kdk] {
        let off = integrate_lh(s0, &m, 13.0, 32, 1e-2, 4_000_000, &LhOpts { stepper, ..opts() });
        let on = integrate_lh(
            s0, &m, 13.0, 32, 1e-2, 4_000_000,
            &LhOpts { stepper, step_limit_f: 0.02, ..opts() },
        );
        println!(
            "{stepper:?} Burrau t = 13: limit off steps {:>8} drift {:.3e} overshoot {}  |  \
             f = 0.02 steps {:>8} drift {:.3e} overshoot {}",
            off.steps, off.drift, off.n_overshoot, on.steps, on.drift, on.n_overshoot
        );
        assert_eq!(
            off.n_overshoot, 0,
            "{stepper:?} overshot an interval with the limit OFF, so the time transformation is \
             not bounding the physical step on its own"
        );
        assert_eq!(on.n_overshoot, 0);
        assert!(!off.den_degenerate, "{stepper:?}: a denominator went non-positive");
        assert!(
            on.drift < off.drift,
            "{stepper:?}: the limit did not improve Burrau's drift ({:e} against {:e}), so the \
             trade this test records has changed sign and the default needs re-deciding",
            on.drift, off.drift
        );
    }

    // The other half of the inversion, on the fixture where it is fatal.
    let (s, mm) = collision_setup((0, 1));
    let capped = integrate_lh(
        s, &mm, 1.0, 32, 1e-3, 2_000_000,
        &LhOpts { step_limit_f: 0.02, ..kdk() },
    );
    let (s, mm) = collision_setup((0, 1));
    let free = integrate_lh(s, &mm, 1.0, 32, 1e-3, 2_000_000, &kdk());
    println!(
        "KDK collision: f = 0.02 -> budget_exhausted {} at {} steps  |  limit off -> finished {} \
         at {} steps, drift {:.3e}",
        capped.budget_exhausted, capped.steps, free.finite, free.steps, free.drift
    );
    assert!(capped.budget_exhausted, "the limit no longer defeats the collision");
    assert!(free.finite && !free.budget_exhausted, "the unlimited arm failed to finish");
}

/// `PerStepRemaining` is Zeno by arithmetic, here as in the other two integrators.
///
/// `dt ~ eta * rem` gives `rem_{n+1} = rem_n (1 - eta)`, so the interval is approached
/// geometrically. Under `clamp_final_step` the relative landing tolerance gives it a floor and it
/// does complete, after `ln(1/eps)/eta` steps against a nominal `1/eta` — so the test is pinned
/// to the unclamped arm, which is the one the property is about.
#[test]
fn the_remaining_time_mode_is_zeno() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let o = integrate_lh(
        s0, &m, 13.0, 32, 1e-2, 200_000,
        &LhOpts { ds_mode: LhDsMode::PerStepRemaining, clamp_final_step: false, ..opts() },
    );
    println!(
        "PerStepRemaining, unclamped: t/t_max = {:.4}  steps = {}  budget_exhausted = {}  \
         drift = {:.3e}",
        o.t / 13.0, o.steps, o.budget_exhausted, o.drift
    );
    println!(
        "  **Print t/t_max before any drift column.** This mode's drift is excellent because it\n  \
         went nowhere."
    );
    assert!(o.budget_exhausted, "the remaining-time mode completed, so it is not Zeno here");
    assert!(o.t / 13.0 < 0.5, "it reached t/t_max = {:.4}", o.t / 13.0);
}
