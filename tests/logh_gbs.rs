//! Gragg-Bulirsch-Stoer over the logH leapfrog — the configuration Mikkola & Merritt recommend.
//!
//! Everything in `results/output/logh_arms.txt` runs the leapfrog **bare**, which is not how the
//! method is meant to be used, and the plan said the rematch had to be a separate experiment
//! rather than something folded into the RK4 comparison. This is that experiment's correctness
//! half.
//!
//! # The one test that decides whether any of it is real
//!
//! Extrapolation in `h^2` gains **two orders per level**, and only because the base method is
//! time-symmetric. If [`gbs::macro_step`]'s observed order does not rise as `2k`, the
//! extrapolation is not working and every number produced under it is a slower leapfrog wearing
//! a better name. That is `the_macro_step_order_rises_by_two_per_level` below, and it is measured
//! on the **macro-step alone**, deliberately.
//!
//! # Why not on the figure-eight, like every other convergence test here
//!
//! Because the sync-boundary landing would cap it at order two and hide the whole effect.
//! `clamp_final_step` sizes the last step of each interval from `dt/ds` evaluated *before* the
//! step — a first-order predictor whose residual is `O(h^2)` — which is exactly why the LogH arm
//! measures 2.04 on that fixture while `LhTime::None`, whose predictor is exact, measures 4.52.
//! A high-order stepper under an `O(h^2)` landing is an `O(h^2)` march.
//!
//! **That is a real limitation of GBS as wired here and not a testing inconvenience**, so it is
//! measured rather than stepped around: `the_landing_residual_dwarfs_the_macro_step_accuracy`
//! records `|s.t - dt_left|` directly and finds it **nine orders above** the macro-step's own
//! error, with the `LhTime::None` arm — whose predictor is exact — at zero.

use prin_rs::integrate::logh::driver::LhDsMode;
use prin_rs::integrate::logh::{gbs, integrate_lh, LhOpts, LhState, LhSystem, LhTime, Stepper};
use prin_rs::physics::{burrau, Cart};
use prin_rs::Vec2;

fn opts(stepper: Stepper) -> LhOpts<f64> {
    LhOpts { stepper, r_coll_frac: 0.0, stop_on_event: false, step_limit_f: 0.0, ..Default::default() }
}

fn fig8() -> (Cart<f64>, [f64; 3], f64) {
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

fn diff13(a: &LhState<f64>, b: &LhState<f64>) -> f64 {
    let (x, y) = (a.to_array14(), b.to_array14());
    (0..13).fold(0.0f64, |w, i| w.max((x[i] - y[i]).abs()))
}

/// A reference for one macro-step: the same map at a resolution far beyond anything under test.
///
/// Not a different integrator and not a different composition — literally `n` plain `kdk` steps,
/// which is what `macro_step`'s levels are made of. So the comparison isolates the
/// **extrapolation** and nothing else.
fn reference(sys: &LhSystem<f64>, s: &LhState<f64>, b: f64, h: f64, n: usize) -> LhState<f64> {
    let hs = h / n as f64;
    let mut cur = *s;
    for _ in 0..n {
        cur = prin_rs::integrate::logh::step::kdk(sys, &cur, b, LhTime::LogH, hs).0;
    }
    cur
}

/// **The decisive test**, and the answer is stronger than an order table.
///
/// A method of global order `p` has local order `p + 1`, so `k` levels over a second-order
/// symmetric base should read `2k + 1`: 3, 5, 7, 9. Measured against a level-10 reference at the
/// same `h`, so no independent integrator's accuracy enters:
///
/// ```text
///        h        k=1          k=2          k=3          k=4
///   8.0e-1     2.548e-6    8.774e-11    1.088e-13    1.086e-13
///   4.0e-1     3.182e-7    2.738e-12    2.838e-13    2.767e-13
///   2.0e-1     3.977e-8    1.210e-13    1.206e-13    1.199e-13
///   1.0e-1     4.972e-9    2.127e-13    2.141e-13    2.127e-13
///
///   per-rung order:  k=1 -> 3.00 at every rung.  k=2 -> 5.00 at the top rung.
/// ```
///
/// `k = 1` reads **3.00 exactly** everywhere and `k = 2` reads **5.00** on the one rung above the
/// floor. `k >= 3` cannot be given an order at all, because **it is already at double-precision
/// round-off at the largest step tested** — `1.09e-13` at `h = 0.8`, for twelve force
/// evaluations. That is not a weaker result than "order 7"; it is a stronger one, and quoting a
/// slope fitted through it would be quoting round-off, which is the same defect as reading an
/// endpoint slope across the figure-eight's floor.
///
/// So the assertions are: the two orders that are measurable, and a **floor** assertion for the
/// levels that are not.
#[test]
fn the_macro_step_order_rises_by_two_per_level() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let sys = LhSystem::new(m);
    let b = sys.b_of(&s0);
    let st = LhState::from_cart(&s0);

    // Reference at level 10, same `h`: far past the floor, and the same map, so this isolates
    // the extrapolation rather than importing another integrator's error.
    let err = |h: f64, k: usize| {
        let r = gbs::macro_step(&sys, &st, b, LhTime::LogH, h, 0.0, 10).state;
        diff13(&gbs::macro_step(&sys, &st, b, LhTime::LogH, h, 0.0, k).state, &r)
    };

    let hs = [8e-1, 4e-1, 2e-1, 1e-1];
    println!("local error of one GBS macro-step against a level-10 reference at the same h:");
    println!("  {:>8} {:>12} {:>12} {:>12} {:>12}", "h", "k=1", "k=2", "k=3", "k=4");
    for &h in &hs {
        println!(
            "  {h:>8.1e} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e}",
            err(h, 1), err(h, 2), err(h, 3), err(h, 4)
        );
    }

    // k = 1: order 3 over the whole range.
    let o1 = (err(hs[0], 1).ln() - err(hs[3], 1).ln()) / (hs[0].ln() - hs[3].ln());
    // k = 2: order 5, measurable on the top rung only -- below it the error is round-off.
    let o2 = (err(hs[0], 2).ln() - err(hs[1], 2).ln()) / (hs[0].ln() - hs[1].ln());
    println!("\n  k=1 order over the full range: {o1:.2}   (want 3)");
    println!("  k=2 order over the top rung  : {o2:.2}   (want 5)");
    assert!((o1 - 3.0).abs() < 0.3, "k=1 local order {o1:.2}, want 3");
    assert!(
        (o2 - 5.0).abs() < 0.6,
        "k=2 local order {o2:.2}, want 5. Each extrapolation level must buy TWO orders; if it \
         buys one, the base method is not time-symmetric and extrapolating in h^2 is invalid."
    );

    // k >= 3 has no order to measure -- it is at the floor. Assert the floor instead, and assert
    // that k = 1 is NOT at it, or this says nothing.
    let floor = err(hs[0], 3);
    println!("  k=3 at h = {:.1e}: {floor:.3e} -- round-off, for {} evaluations", hs[0], 3 * 4);
    assert!(floor < 1e-12, "k=3 at h={:.1e} is {floor:e}, not at round-off", hs[0]);
    assert!(
        err(hs[0], 1) > 1e-8,
        "k=1 is also at the floor, so the k=3 assertion is about the reference and not about \
         the extrapolation"
    );
}

/// Each level `j` costs `SEQ[j]` evaluations, so level `k` has spent `sum = k(k+1)`.
///
/// Asserted rather than trusted because `Stepper::Gbs.evals_per_step()` is deliberately `0`:
/// there is no fixed multiplier, and this is the only thing that makes the cost column honest.
#[test]
fn the_evaluation_count_is_the_sum_of_the_sequence() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let sys = LhSystem::new(m);
    let b = sys.b_of(&s0);
    let st = LhState::from_cart(&s0);
    assert_eq!(Stepper::Gbs.evals_per_step(), 0, "a fixed cost was claimed for GBS");
    assert!(!Stepper::Gbs.has_fixed_cost());
    for k in 1..=6usize {
        let g = gbs::macro_step(&sys, &st, b, LhTime::LogH, 1e-2, 0.0, k);
        let want: usize = gbs::SEQ[..k].iter().sum();
        println!("  k = {k}: evals {} against sum(SEQ[..{k}]) = {want} = k(k+1) = {}", g.evals, k * (k + 1));
        assert_eq!(g.evals, want);
        assert_eq!(g.evals, k * (k + 1));
        assert_eq!(g.k_used, k);
        assert!(!g.converged, "tol = 0 must never be met");
    }
}

/// The symmetry the `h^2` expansion rests on: `n` equal KDK steps forward then reversed must
/// return the start.
///
/// **This is a property of the base method, not of the extrapolation**, and if it fails then
/// extrapolating in `h^2` is invalid however good the order table looks — the expansion would
/// carry odd powers and each level would buy one order, not two, with the fit picking up
/// whatever the mixture happens to be.
#[test]
fn the_leapfrog_base_is_time_symmetric() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let sys = LhSystem::new(m);
    let b = sys.b_of(&s0);
    let st = LhState::from_cart(&s0);
    for n in [2usize, 8, 32] {
        let fwd = reference(&sys, &st, b, 1e-2, n);
        // Reversing `s` reverses the velocities and marches back; `t` runs backwards with it.
        let flipped =
            LhState { r: fwd.r, v: std::array::from_fn(|i| -fwd.v[i]), t: 0.0, w: fwd.w };
        let back = reference(&sys, &flipped, b, 1e-2, n);
        let err = (0..3).fold(0.0f64, |w, i| w.max((back.r[i] - st.r[i]).norm()));
        println!("  n = {n:>2}: reversal recovery {err:.3e}");
        assert!(err < 1e-13, "the KDK base is not time-symmetric at n = {n}: {err:e}");
    }
}

/// `gbs_unconverged` must fire when the tolerance cannot be met, and **must not** when it can.
///
/// A counter that always fires passes as easily as one that never does, so both arms are here.
#[test]
fn reaching_the_level_cap_without_converging_is_recorded() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let hard = integrate_lh(
        s0, &m, 2.0, 16, 1e-2, 4_000_000,
        &LhOpts { gbs_tol: 0.0, gbs_k_max: 3, ..opts(Stepper::Gbs) },
    );
    let easy = integrate_lh(
        s0, &m, 2.0, 16, 1e-2, 4_000_000,
        &LhOpts { gbs_tol: 1e-8, gbs_k_max: 8, ..opts(Stepper::Gbs) },
    );
    println!(
        "  tol = 0   : unconverged {} of {} steps, mean level {:.2}, evals {}",
        hard.gbs_unconverged, hard.steps, hard.gbs_levels as f64 / hard.steps as f64, hard.force_evals
    );
    println!(
        "  tol = 1e-8: unconverged {} of {} steps, mean level {:.2}, evals {}",
        easy.gbs_unconverged, easy.steps, easy.gbs_levels as f64 / easy.steps as f64, easy.force_evals
    );
    // **Not `== steps`.** `tol = 0` still admits `err <= 0` when two consecutive levels agree
    // **bitwise**, which is a real convergence and happens on 162 of 1764 macro-steps here --
    // the extrapolation reaching round-off, which is exactly what the order test shows it doing
    // in three levels. Asserting equality would have been asserting that never happens.
    assert!(
        hard.gbs_unconverged as f64 > 0.8 * hard.steps as f64,
        "an impossible tolerance recorded only {} misses in {} steps",
        hard.gbs_unconverged, hard.steps
    );
    assert_eq!(
        easy.gbs_unconverged, 0,
        "a reachable tolerance recorded misses, so the counter fires regardless and says nothing"
    );
    assert!(easy.force_evals < hard.force_evals, "the adaptive arm cost no less than the capped one");
}

/// **The landing residual dwarfs the macro-step's accuracy, and that is the whole story for GBS
/// here.**
///
/// The macro-step reaches round-off in three levels — `1.09e-13` for twelve force evaluations —
/// and the march is nonetheless *worse* than a bare leapfrog at matched cost. Something between
/// the two throws the accuracy away, and naming it needs a measurement.
///
/// `LhOut::land_residual_max` records `|s.t - dt_left|` at each boundary **before the clock is
/// clamped**, since clamping is what hides it. Measured on the figure-eight, 32 intervals:
///
/// ```text
///   LogH          eta 8e-1  1.2960e-4     eta 4e-1  1.3465e-4
///                 eta 2e-1  1.3748e-4     eta 1e-1  3.4175e-5
///   LhTime::None  zero to round-off at every step size
/// ```
///
/// **Nine orders of magnitude above the macro-step's own error.** A step accurate to `1e-13`
/// followed by a landing that misses the boundary by `1.3e-4` in time is a march with a `1.3e-4`
/// error, and no amount of extrapolation touches it.
///
/// # Two things this test deliberately does not claim
///
/// The residual was expected to fall as `h^2` — a first-order predictor for the time increment
/// missing by `O(h^2)` — and **it does not**: 1.296e-4, 1.3465e-4, 1.3748e-4, 3.4175e-5, a fitted
/// slope of 0.64. *A quantity that does not fall when the step shrinks is not a step-size
/// problem*, which is this project's own oldest diagnostic and applies to its own landing here.
/// The likely reason is that the **final** step's size is set by whatever is left of the interval
/// rather than by `eta`, so it moves with how `1/eta` happens to divide — the same
/// round-off-in-the-division effect that makes the unclamped fixed-step control non-monotone.
/// That is a hypothesis and is **not** measured, so the assertion is on the magnitude, which is
/// unambiguous, and not on a scaling that is not there.
///
/// Nor is an order fitted on the figure-eight closure: its nine-significant-digit initial
/// conditions floor it at `~5e-8`, which pins the `LhTime::None` arm at `4.102e-8` at *every*
/// step size and makes it useless as a control for that.
///
/// The `LhTime::None` arm **is** the control here, and a good one: `dt/ds` is exactly 1, the
/// prediction is exact, and the residual is zero. If it were not, the cap would be something
/// other than the predictor.
#[test]
fn the_landing_residual_dwarfs_the_macro_step_accuracy() {
    let (s0, m, t) = fig8();
    let etas = [8e-1, 4e-1, 2e-1, 1e-1];
    for time in [LhTime::LogH, LhTime::None] {
        let r: Vec<f64> = etas
            .iter()
            .map(|&eta| {
                // **`land_iterate: false`, and the paired arm below is why that is honest.**
                // This test characterises the residual left by `clamp_final_step` alone -- the
                // first-order predictor. The secant correction removes it, which deletes the
                // subject rather than contradicting the finding.
                let o = integrate_lh(
                    s0, &m, t, 32, eta, 40_000_000,
                    &LhOpts {
                        time,
                        land_iterate: false,
                        gbs_tol: 1e-13,
                        gbs_k_max: 6,
                        ..opts(Stepper::Gbs)
                    },
                );
                assert!(o.finite, "{time:?} went non-finite at eta = {eta}");
                o.land_residual_max
            })
            .collect();
        println!("{time:?}: landing residual |s.t - dt_left| at the boundary");
        for (eta, x) in etas.iter().zip(&r) {
            println!("   eta {eta:.1e}   {x:.4e}");
        }
        let worst = r.iter().cloned().fold(0.0f64, f64::max);
        match time {
            // TTL carries `W` as a fourteenth component and its landing behaviour is a separate
            // question from the one this test asks; it is excluded explicitly rather than folded
            // into the `LogH` arm, where it would be measured under an assertion written for a
            // different transformation.
            LhTime::Ttl => {}
            LhTime::LogH => {
                let order = (r[0].ln() - r[3].ln()) / (etas[0].ln() - etas[3].ln());
                println!(
                    "   worst {worst:.3e}, fitted slope {order:.2} -- NOT the h^2 expected, and \
                     the assertion is on the magnitude for that reason"
                );
                assert!(
                    worst > 1e-8,
                    "the landing residual is {worst:e}, which is near the macro-step's own \
                     accuracy -- the account of why GBS buys nothing here would then be wrong"
                );
                // **THE PAIRED ARM.** The same ladder with the secant correction on must drive
                // the residual down by orders. Without this the `land_iterate: false` above would
                // be an unexplained pin, indistinguishable from one added to keep a stale
                // assertion green.
                let corrected = etas
                    .iter()
                    .map(|&eta| {
                        let o = integrate_lh(
                            s0, &m, t, 32, eta, 40_000_000,
                            &LhOpts {
                                time,
                                land_iterate: true,
                                gbs_tol: 1e-13,
                                gbs_k_max: 6,
                                ..opts(Stepper::Gbs)
                            },
                        );
                        o.land_residual_max
                    })
                    .fold(0.0f64, f64::max);
                println!("   with the secant landing: worst {corrected:.3e}");
                assert!(
                    corrected < worst * 1e-3,
                    "the secant landing left {corrected:e} against the clamp's {worst:e} -- it is \
                     not removing the residual, so pinning this test away from it is hiding that"
                );
            }
            LhTime::None => {
                println!("   worst {worst:.3e}  (control: dt/ds is exactly 1, prediction exact)");
                assert!(
                    worst < 1e-14,
                    "the control also carries a landing residual ({worst:e}), so the cap is NOT \
                     the predictor and the LogH assertion above describes a coincidence"
                );
            }
        }
    }
}

/// Does GBS actually beat the bare leapfrog **at matched force evaluations**?
///
/// The practical question, and the one the whole rematch is for. Matched on evaluations rather
/// than steps, because a GBS macro-step costs `k(k+1)` of them and a KDK step costs one — a
/// comparison at equal `steps` would be a factor of forty in disguise.
#[test]
fn gbs_against_the_bare_leapfrog_at_matched_evaluations() {
    let (s0, m, t) = fig8();
    let run = |stepper: Stepper, eta: f64, tol: f64| {
        let o = integrate_lh(
            s0, &m, t, 32, eta, 40_000_000,
            &LhOpts { gbs_tol: tol, gbs_k_max: 8, ds_mode: LhDsMode::PerStepInterval, ..opts(stepper) },
        );
        let err = (0..3).fold(0.0f64, |w, i| {
            w.max((o.state.r[i] - s0.r[i]).norm()).max((o.state.v[i] - s0.v[i]).norm())
        });
        (err, o.force_evals, o.steps)
    };
    println!("figure-eight closure at matched evaluations:");
    let (ek, vk, sk) = run(Stepper::Kdk, 2.5e-3, 0.0);
    println!("  KDK  eta 2.5e-3       closure {ek:.4e}  evals {vk:>8}  steps {sk:>8}");
    for eta in [8e-2, 4e-2, 2e-2] {
        let (eg, vg, sg) = run(Stepper::Gbs, eta, 1e-13);
        println!("  GBS  eta {eta:.1e}       closure {eg:.4e}  evals {vg:>8}  steps {sg:>8}");
    }
    println!(
        "\n  **Read the evals column, not eta.** A GBS macro-step costs k(k+1) evaluations and a\n  \
         KDK step costs one, so equal `eta` is a factor of tens in disguise."
    );
}

/// **The secant landing removes the cap, and KDK is the control that says that is what it does.**
///
/// `land_iterate` re-takes the landing step with a secant on `t(ds)`. Three things must hold
/// together, and the third is what makes the first two mean anything:
///
///   - the landing residual falls to **round-off** — this is the quantity being fixed;
///   - a stepper better than second order gets **much** more accurate — this is the consequence;
///   - **KDK does not improve** — this is the control. It is second order, so an `O(h^2)` landing
///     was never its binding constraint. A correction that improved every arm would be doing
///     something other than removing an order-two cap, and the whole account would be wrong.
#[test]
fn the_secant_landing_removes_the_order_cap_and_kdk_is_the_control() {
    let (s0, m, t) = fig8();
    let run = |stepper: Stepper, land: bool, eta: f64| {
        let o = integrate_lh(
            s0, &m, t, 32, eta, 40_000_000,
            &LhOpts { land_iterate: land, gbs_tol: 1e-13, gbs_k_max: 6, ..opts(stepper) },
        );
        assert!(o.finite);
        let e = (0..3).fold(0.0f64, |w, i| {
            w.max((o.state.r[i] - s0.r[i]).norm()).max((o.state.v[i] - s0.v[i]).norm())
        });
        (e, o.land_residual_max, o.force_evals, o.land_iters)
    };

    for stepper in [Stepper::Kdk, Stepper::Rk4, Stepper::Gbs] {
        let (e0, r0, v0, _) = run(stepper, false, 4e-2);
        let (e1, r1, v1, c1) = run(stepper, true, 4e-2);
        println!(
            "  {stepper:?}: closure {e0:.4e} -> {e1:.4e}   resid {r0:.3e} -> {r1:.3e}   \
             evals {v0} -> {v1} ({c1} corrections)"
        );
        assert!(r1 < 1e-13, "{stepper:?}: the corrected landing residual is {r1:e}, not round-off");
        assert!(r0 > 1e-8, "{stepper:?}: the UNcorrected residual is {r0:e}, so there was nothing \
                            to correct and this comparison says nothing");
        assert!(c1 > 0, "{stepper:?}: no corrections were taken");
        match stepper {
            // The control. Second order: the landing was never its constraint, so it must not
            // improve -- and it is in fact slightly worse, paying for corrections it cannot use.
            Stepper::Kdk => assert!(
                e1 > e0 * 0.9,
                "KDK improved from {e0:e} to {e1:e}. A second-order stepper cannot be limited by \
                 an O(h^2) landing, so the correction is doing something other than what this \
                 test claims."
            ),
            // Better than second order: the cap was binding and removing it must show.
            _ => assert!(
                e1 < e0 * 0.1,
                "{stepper:?} improved only {e0:e} -> {e1:e}; the landing was expected to be its \
                 binding constraint"
            ),
        }
    }
}
