//! Phase 2: the Heggie march, against the instruments this project already trusts.
//!
//! Everything here has an AZ counterpart, run on the same fixture, so a Heggie number can be read
//! against a known one rather than against an expectation.

use prin_rs::integrate::az;
use prin_rs::integrate::heggie::driver::{q_dot, HgDtauMode};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts, HgSystem, HgTime};
use prin_rs::physics::{burrau, energy, Cart};
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

fn opts() -> HgOpts<f64> {
    HgOpts::default()
}

// ---------------------------------------------------------------------------------------------
// Eq. (7): the enlarged equations of motion.

/// `q_dot` is `dH/dp` for Heggie's Eq. (6), finite-differenced.
///
/// The index pattern is the trap: the momentum with the **further** cyclic index carries the
/// reciprocal of the **nearer** mass. The mutation arm swaps them, which is a change no
/// dimensional or symmetry argument would catch.
///
/// **The step is large — `h = 1` — and that is deliberate.** Eq. (6) is exactly quadratic in the
/// momenta, so central differencing has no truncation error at any step size; the only error is
/// roundoff. At `h = 1e-6` the difference `H(p+h) - H(p-h)` cancels a potential term
/// `-m_j m_k/|q_i|` that is *independent of p* and can reach `1e3` when a separation is small,
/// and the test read **1.6e-7** — entirely that cancellation, amplified by `1/h`. Same shape as
/// `Gamma*` being exactly differenced: ask what order the function is in the variable before
/// choosing the step.
#[test]
fn eq7_is_the_gradient_of_the_enlarged_hamiltonian() {
    let m = [3.0, 4.0, 5.0];
    let sys = HgSystem::new(m);
    let mut rng = SplitMix64::new(0x077);
    let mut worst = 0.0f64;
    let mut worst_swapped = 0.0f64;

    for _ in 0..256 {
        let q: [Vec2<f64>; 3] =
            std::array::from_fn(|_| Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)));
        let p: [Vec2<f64>; 3] =
            std::array::from_fn(|_| Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)));
        let an = q_dot(&sys, &p);

        for i in 0..3 {
            for c in 0..2 {
                let h = 1.0;
                let (mut hi, mut lo) = (p, p);
                if c == 0 {
                    hi[i].x += h;
                    lo[i].x -= h;
                } else {
                    hi[i].y += h;
                    lo[i].y -= h;
                }
                let fd = (sys.energy_enlarged(&q, &hi) - sys.energy_enlarged(&q, &lo)) / (2.0 * h);
                let a = if c == 0 { an[i].x } else { an[i].y };
                let scale = a.abs().max(fd.abs()).max(1.0);
                worst = worst.max((a - fd).abs() / scale);

                // The mutation arm: reciprocals swapped between the two coupling partners.
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                let bad = p[i] / sys.mu[i] - p[k] * sys.inv_m[k] - p[j] * sys.inv_m[j];
                let b = if c == 0 { bad.x } else { bad.y };
                worst_swapped = worst_swapped.max((b - fd).abs() / scale);
            }
        }
    }
    println!("Eq. (7) against the finite-differenced Eq. (6): {worst:.3e}");
    println!("  reciprocals swapped between partners:         {worst_swapped:.3e}");
    assert!(worst < 1e-12, "Eq. (7) is not dH/dp: {worst:e}");
    assert!(worst_swapped > 1e-2, "swapping the reciprocals still agreed at {worst_swapped:e}");
}

// ---------------------------------------------------------------------------------------------
// The two-body radial collision — for ALL THREE pairs, which is the globality claim.

/// Equal masses, unit separation, at rest, with the third body 1000 away and off axis.
fn collision_setup(pair: (usize, usize)) -> (Cart<f64>, [f64; 3]) {
    let third = 3 - pair.0 - pair.1;
    let mut r = [Vec2::zero(); 3];
    r[pair.0] = Vec2::new(-0.5, 0.0);
    r[pair.1] = Vec2::new(0.5, 0.0);
    r[third] = Vec2::new(0.1, 1000.0);
    // Centre it: `to_reg`'s Eq. (8b) reduction and `cart_from_enlarged` both assume the COM
    // frame, and an uncentred input would fail for a reason that is not the physics.
    let m = [1.0f64; 3];
    let c: Vec2<f64> = (r[0] + r[1] + r[2]) / 3.0;
    for x in r.iter_mut() {
        *x -= c;
    }
    (Cart::new(r, [Vec2::zero(); 3]), m)
}

/// **The globality claim, stated as a test.** AZ regularises the two pairs sharing its reference
/// body and leaves the third unregularised, so `tests/two_body_collision.rs` has to first assert
/// that the colliding pair is one of the regularised ones — if the geometry put it on the third
/// side the test would measure nothing.
///
/// Heggie has no third side. The same collision is run **three times, once per pair**, and all
/// three must pass with no configuration-dependent preamble. A version that passed for one pair
/// and not the others would be AZ wearing a different name.
#[test]
fn radial_collision_passes_through_for_every_pair() {
    for pair in [(0usize, 1usize), (0, 2), (1, 2)] {
        let (s, m) = collision_setup(pair);
        let o = integrate_hg(s, &m, 1.0, 32, 1e-3, 4_000_000, &opts());
        println!(
            "pair {pair:?}: d_min = {:.3e}  |dE/E| = {:.3e}  steps = {}  \
             gamma = {:.2e}  sum_q = {:.2e}",
            o.d_min, o.drift, o.steps, o.gamma_max, o.sum_q_max
        );
        assert!(o.finite, "pair {pair:?} went non-finite");
        assert!(o.d_min < 1e-10, "pair {pair:?}: d_min = {:e}", o.d_min);
        assert!(o.drift < 1e-12, "pair {pair:?}: drift = {:e}", o.drift);
    }
}

// ---------------------------------------------------------------------------------------------
// Burrau, against AZ.

/// Heggie and AZ integrate the same system, so at a tight step they must agree on the trajectory.
///
/// This is the cross-check that does not exist in the paper: Heggie compared his method against
/// Aarseth-Zare on *his* problems, and this compares the two ports on *this* project's.
#[test]
fn burrau_agrees_with_az() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let t_max = 6.0;

    let hg = integrate_hg(
        s0,
        &m,
        t_max,
        64,
        1e-4,
        20_000_000,
        &HgOpts { step_limit_f: 0.0, ..opts() },
    );
    let azo = az::integrate_az(s0, &m, t_max, 64, 1e-4, 20_000_000, None);

    let mut worst = 0.0f64;
    for i in 0..3 {
        worst = worst
            .max((hg.state.r[i] - azo.state.r[i]).norm())
            .max((hg.state.v[i] - azo.state.v[i]).norm());
    }
    println!("Burrau to t = {t_max}, eta = 1e-4:");
    println!("  Heggie  drift {:.3e}  steps {:>9}  gamma {:.2e}  sum_q {:.2e}",
        hg.drift, hg.steps, hg.gamma_max, hg.sum_q_max);
    println!("  AZ      drift {:.3e}  steps {:>9}", azo.drift, azo.steps);
    println!("  max |state difference| = {worst:.3e}");
    println!("  cost ratio (Heggie/AZ steps) = {:.2}", hg.steps as f64 / azo.steps as f64);
    assert!(hg.finite && azo.finite);
    assert!(worst < 1e-6, "the two integrators disagree by {worst:e}");
}

/// The Burrau constants, through Heggie's own reconstruction rather than through `Cart`.
#[test]
fn burrau_constants_survive_the_enlarged_construction() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let sys = HgSystem::new(m);
    let (q, p) = sys.enlarged_from_cart(&s0);
    let back = sys.cart_from_enlarged(&q, &p);

    let mtot: f64 = m.iter().sum();
    let r = energy::hyperradius(&back.r, &m);
    let e = sys.energy_enlarged(&q, &p);
    println!("M = {mtot}  R = {r:.4}  E = {e:.4}");
    assert!((mtot - 12.0).abs() < 1e-12, "M = {mtot}");
    assert!((r - 2.2361).abs() < 1e-4, "R = {r}");
    assert!((e - (-12.8167)).abs() < 1e-4, "E = {e}");
}

// ---------------------------------------------------------------------------------------------
// Convergence order, on the figure eight.

/// Chenciner-Montgomery, equal masses, `G = 1`. Exactly periodic, so `|state(T) - state(0)|` is a
/// pure error with no reference trajectory and no chaos.
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

/// **Read the order, not the error.**
///
/// The per-rung two-point estimates are noise — AZ's `fixed+clamp` arm runs 2.34, 2.58, 1.36,
/// 6.49 — so the endpoint slope is what is quoted, exactly as it is for AZ.
///
/// AZ's measured orders on this same fixture: `fixed+overshoot` **1.13**, `perstep+overshoot`
/// **1.06**, `fixed+clamp` **3.06**, `perstep+clamp` **2.08**. The clamp is a correctness
/// property and this asserts it here too, through its own control arm — an unclamped Heggie march
/// must come out near first order, or the clamp is not the thing being measured.
#[test]
fn convergence_order_on_the_figure_eight() {
    let (s0, m, t) = figure_eight();
    let etas = [2e-2, 1e-2, 5e-3, 2.5e-3, 1e-3];

    for clamp in [true, false] {
        let mut errs = Vec::new();
        for &eta in &etas {
            let o = integrate_hg(
                s0,
                &m,
                t,
                32,
                eta,
                40_000_000,
                &HgOpts { clamp_final_step: clamp, step_limit_f: 0.0, ..opts() },
            );
            assert!(o.finite, "figure eight went non-finite at eta = {eta}");
            errs.push((eta, closure_err(&o.state, &s0), o.steps, o.drift));
        }
        let slope = (errs[0].1.ln() - errs[4].1.ln()) / (etas[0].ln() - etas[4].ln());
        println!("clamp_final_step = {clamp}:");
        for (eta, e, st, d) in &errs {
            println!("   eta {eta:.1e}   closure {e:.4e}   steps {st:>8}   drift {d:.2e}");
        }
        println!("   endpoint slope = {slope:.2}");
        if clamp {
            assert!(slope > 1.7, "the clamped march is only order {slope:.2}");
        } else {
            assert!(
                slope < 1.5,
                "the UNCLAMPED march is order {slope:.2}; if it is already second order the \
                 clamp is not what the clamped result is measuring"
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// The scale gauge, and the residual with no AZ analogue.

/// Under `r -> alpha r`, `t -> alpha^{3/2} t`, the physics is identical. Eq. (22) makes `tau`
/// itself scale-invariant, which Eq. (20) does not — so both are checked, and the point is that
/// **neither may drift**, not that one is prettier.
///
/// `Real::LAND_EPS_REL` is relative to `dt_left` precisely so this holds; an absolute tolerance
/// broke the AZ version of this test at 4.24e-15.
#[test]
fn the_march_respects_the_scale_gauge() {
    let (s0, m, t) = figure_eight();
    for time in [HgTime::Product, HgTime::default()] {
        let base = integrate_hg(s0, &m, t, 32, 5e-3, 4_000_000, &HgOpts { time, ..opts() });
        for alpha in [0.25f64, 4.0] {
            let scaled = Cart::new(
                std::array::from_fn(|i| s0.r[i] * alpha),
                std::array::from_fn(|i| s0.v[i] / alpha.sqrt()),
            );
            let o = integrate_hg(
                scaled,
                &m,
                t * alpha.powf(1.5),
                32,
                5e-3,
                4_000_000,
                &HgOpts { time, ..opts() },
            );
            let mut worst = 0.0f64;
            for i in 0..3 {
                worst = worst
                    .max((o.state.r[i] / alpha - base.state.r[i]).norm())
                    .max((o.state.v[i] * alpha.sqrt() - base.state.v[i]).norm());
            }
            println!("{time:?}  alpha = {alpha}: max rescaled difference = {worst:.3e}");
            assert!(worst < 1e-11, "{time:?} at alpha = {alpha}: {worst:e}");
        }
    }
}

/// `sum q_i = 0` is Heggie's Eq. (9), an integral of the enlarged motion with **no AZ analogue**.
///
/// It is not conserved by construction — nothing projects onto the constraint surface — so it is
/// a genuine free residual and the march has to be shown to hold it. The mutation arm is a state
/// deliberately pushed off the constraint surface, which must read large: a residual that is
/// small for every input is not measuring the constraint.
#[test]
fn the_sum_q_constraint_survives_the_march() {
    let (s0, m, t) = figure_eight();
    let o = integrate_hg(s0, &m, t, 32, 5e-3, 4_000_000, &opts());
    println!("figure eight: sum_q residual max = {:.3e}, gamma max = {:.3e}", o.sum_q_max, o.gamma_max);
    assert!(o.finite);
    assert!(o.sum_q_max < 1e-10, "sum q_i drifted to {:e}", o.sum_q_max);

    // Off the constraint surface on purpose.
    let sys = HgSystem::new(m);
    let (mut s, _) = sys.to_reg(&s0);
    s.u[0] = s.u[0] * 1.5;
    let bad = sys.sum_q_residual(&s);
    println!("  a state pushed off the constraint surface reads {bad:.3e}");
    assert!(bad > 1e-2, "the residual is blind to a violated constraint: {bad:e}");
}

/// Heggie's two time transformations integrate the **same trajectory**; Eqs. (22)-(24) are a
/// reparameterisation, not different physics. So they must land in the same place.
#[test]
fn the_two_time_transformations_reach_the_same_state() {
    let (s0, m, t) = figure_eight();
    let a = integrate_hg(
        s0, &m, t, 32, 2e-3, 20_000_000,
        &HgOpts { time: HgTime::Product, step_limit_f: 0.0, ..opts() },
    );
    let b = integrate_hg(
        s0, &m, t, 32, 2e-3, 20_000_000,
        &HgOpts { time: HgTime::default(), step_limit_f: 0.0, ..opts() },
    );
    let worst = closure_err(&a.state, &b.state);
    println!("Eq. (20) vs Eqs. (22)-(24) at eta = 2e-3:");
    println!("  Product     drift {:.3e}  steps {:>8}  closure {:.3e}",
        a.drift, a.steps, closure_err(&a.state, &s0));
    println!("  SumPow32    drift {:.3e}  steps {:>8}  closure {:.3e}",
        b.drift, b.steps, closure_err(&b.state, &s0));
    println!("  difference between them = {worst:.3e}");
    assert!(a.finite && b.finite);
    assert!(worst < 1e-5, "the two time transformations disagree by {worst:e}");
}

/// `HgDtauMode::PerStepRemaining` is Zeno by arithmetic, and it must be **shown** to be, or it is
/// a mode nobody has checked. The clamp's relative landing tolerance is what lets it complete at
/// all; without it AZ's analogue reaches `t/t_max = 0.008` on its whole budget.
#[test]
fn the_remaining_time_mode_is_zeno() {
    let (s0, m, t) = figure_eight();
    let zeno = integrate_hg(
        s0, &m, t, 8, 1e-2, 2_000_000,
        &HgOpts {
            dtau_mode: HgDtauMode::PerStepRemaining,
            clamp_final_step: false,
            step_limit_f: 0.0,
            ..opts()
        },
    );
    let normal = integrate_hg(
        s0, &m, t, 8, 1e-2, 2_000_000,
        &HgOpts { clamp_final_step: false, step_limit_f: 0.0, ..opts() },
    );
    println!("PerStepRemaining, no clamp:  t/t_max = {:.4}  steps {}  budget_exhausted {}",
        zeno.t / t, zeno.steps, zeno.budget_exhausted);
    println!("PerStepInterval, no clamp:   t/t_max = {:.4}  steps {}",
        normal.t / t, normal.steps);
    assert!(normal.t / t > 0.99, "the control did not complete: {:.4}", normal.t / t);
    assert!(
        zeno.t / t < 0.5,
        "PerStepRemaining completed ({:.4}); the Zeno mode is not Zeno here and the AZ result \
         does not transfer",
        zeno.t / t
    );
}

// ---------------------------------------------------------------------------------------------
// Heggie's own §3 claim, on this project's code.

/// Permute the body labels: physically the same system, different index order.
fn permute(s: &Cart<f64>, m: &[f64; 3], p: [usize; 3]) -> (Cart<f64>, [f64; 3]) {
    (
        Cart::new(std::array::from_fn(|i| s.r[p[i]]), std::array::from_fn(|i| s.v[p[i]])),
        std::array::from_fn(|i| m[p[i]]),
    )
}

/// **Is the march independent of the body labelling?**
///
/// Heggie's §3, on his close-triple-encounter time reversals: *"the success of the time-reversals
/// recorded in Table III does not depend on any judicious choice for the initial labelling of the
/// bodies, which is the case with the method of Aarseth and Zare."* That is the one remark in the
/// paper that bears on why this port exists, so it is tested here rather than quoted.
///
/// Relabelling is a physically empty operation, so a *correct* integrator of either kind should
/// be covariant under it. Heggie is covariant **by construction** — nothing in `HgSystem` reads
/// an index for anything but cyclic bookkeeping. AZ's reference body is chosen from geometry, so
/// it is covariant too *except* where `longest_index`'s strict-`>` first-max tie-break decides,
/// and where relabelling changes the LC branch cut's alignment.
///
/// **Both numbers are printed and only Heggie's is asserted.** Asserting AZ fails would be
/// asserting a tie occurs, which is a property of the fixture and not of the method.
///
/// **Measured, AZ is label-covariant too — 3.2e-15, better than Heggie's 1.8e-14.** So the
/// paper's remark does NOT reproduce as a label-permutation property on Burrau at this horizon,
/// and quoting it as one would be wrong. What Heggie is describing is his own Table III
/// comparison, where AZ needs the bodies labelled so the regularised pairs are the ones that
/// close; that is a statement about setting a run up, not about covariance under relabelling.
/// Recorded here because it is a motivating quote for this whole port and it means less than it
/// appears to.
#[test]
fn the_march_is_independent_of_the_body_labelling() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let (t_max, n_sync, eta) = (6.0, 64, 1e-4);
    let perms = [[0usize, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

    let base_h = integrate_hg(
        s0, &m, t_max, n_sync, eta, 20_000_000,
        &HgOpts { step_limit_f: 0.0, ..opts() },
    );
    let base_a = az::integrate_az(s0, &m, t_max, n_sync, eta, 20_000_000, None);

    let (mut worst_h, mut worst_a) = (0.0f64, 0.0f64);
    for p in perms {
        let (sp, mp) = permute(&s0, &m, p);
        let h = integrate_hg(
            sp, &mp, t_max, n_sync, eta, 20_000_000,
            &HgOpts { step_limit_f: 0.0, ..opts() },
        );
        let a = az::integrate_az(sp, &mp, t_max, n_sync, eta, 20_000_000, None);
        for i in 0..3 {
            worst_h = worst_h.max((h.state.r[i] - base_h.state.r[p[i]]).norm());
            worst_a = worst_a.max((a.state.r[i] - base_a.state.r[p[i]]).norm());
        }
    }
    println!("Burrau to t = {t_max}, over all five non-identity label permutations:");
    println!("  Heggie  max |difference from the relabelled reference| = {worst_h:.3e}");
    println!("  AZ      max |difference from the relabelled reference| = {worst_a:.3e}");
    println!("  Heggie is covariant by construction -- no index is read for anything but cyclic");
    println!("  bookkeeping. AZ's is a geometric choice, so it is covariant except where the");
    println!("  argmax tie-break or the LC branch alignment decides. Only Heggie's is asserted.");
    assert!(
        worst_h < 1e-12,
        "Heggie is not label-covariant ({worst_h:e}); an index is being read that should not be"
    );
}

/// **Time reversal**, Heggie's Table III in miniature: integrate, negate the momenta, integrate
/// the same physical time, and see how well the initial conditions come back.
///
/// A pure error measure with no reference trajectory — the same virtue as the figure eight, but
/// on a chaotic system, so it also bounds how far the trajectory has been contaminated. His
/// Table III reports `3e-11` after 250 RKF7(8) steps at tolerance `1e-12`; this is a different
/// stepper with a different error control, so the step counts are **not** comparable and no
/// assertion is made against his figure. What transfers is the shape of the test.
#[test]
fn the_march_reverses() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let (t_max, eta) = (4.0, 1e-5);

    let fwd = integrate_hg(
        s0, &m, t_max, 32, eta, 40_000_000,
        &HgOpts { step_limit_f: 0.0, ..opts() },
    );
    assert!(fwd.finite);
    let flipped = Cart::new(fwd.state.r, std::array::from_fn(|i| -fwd.state.v[i]));
    let back = integrate_hg(
        flipped, &m, t_max, 32, eta, 40_000_000,
        &HgOpts { step_limit_f: 0.0, ..opts() },
    );
    assert!(back.finite);

    let err = (0..3).fold(0.0f64, |w, i| w.max((back.state.r[i] - s0.r[i]).norm()));
    println!("Burrau, forward to t = {t_max} at eta = {eta:.0e} then reversed:");
    println!("  d_min along the way = {:.3e}  steps = {}", fwd.d_min, fwd.steps);
    println!("  recovery error in position = {err:.3e}");
    println!("  drift out {:.2e}, back {:.2e}", fwd.drift, back.drift);
    assert!(err < 1e-6, "the reversal did not recover the initial conditions: {err:e}");
}

/// **Does a FROZEN reference body make AZ label-dependent?** The discriminator for the negative
/// above.
///
/// The covariance test found AZ label-covariant, so Heggie's §3 remark does not reproduce in that
/// form. The candidate explanation is that this AZ **re-chooses its reference body at every sync
/// boundary** — so a poor initial labelling is discarded at the first boundary and never costs
/// anything again, and the contrast he drew would only be visible against an AZ whose reference
/// is fixed at the start.
///
/// `forced_refs` freezes it. Under a label permutation, forcing the same *index* selects a
/// different *physical* body, so this asks precisely what a bad initial choice costs when the
/// march cannot correct it. Heggie is run through the identical permutation as the control: it
/// has no reference to freeze, so it must not move whatever is done here.
///
/// **Measured: free 3.23e-15, frozen 3.41e-6 — a factor of 1.06e9.** Confirmed. The irony is now
/// measured rather than argued: **the re-registration that causes the wedges is the same
/// mechanism that makes this AZ insensitive to its initial labelling.** Heggie's contrast is
/// real and is against a fixed-reference AZ; this port already spends the cost that buys it off,
/// and the wedges are what it spends.
#[test]
fn a_frozen_reference_body_is_what_makes_az_label_dependent() {
    let m = burrau::masses::<f64>();
    let s0 = burrau::state::<f64>();
    let (t_max, n_sync, eta) = (6.0, 64, 1e-4);
    let perms = [[0usize, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let frozen = vec![0u8; n_sync];

    let base_free = az::integrate_az(s0, &m, t_max, n_sync, eta, 20_000_000, None);
    let base_frozen =
        az::integrate_az(s0, &m, t_max, n_sync, eta, 20_000_000, Some(&frozen));
    let base_hg = integrate_hg(
        s0, &m, t_max, n_sync, eta, 20_000_000,
        &HgOpts { step_limit_f: 0.0, ..opts() },
    );

    let (mut free, mut froz, mut hgw) = (0.0f64, 0.0f64, 0.0f64);
    for p in perms {
        let (sp, mp) = permute(&s0, &m, p);
        let a_free = az::integrate_az(sp, &mp, t_max, n_sync, eta, 20_000_000, None);
        let a_froz =
            az::integrate_az(sp, &mp, t_max, n_sync, eta, 20_000_000, Some(&frozen));
        let h = integrate_hg(
            sp, &mp, t_max, n_sync, eta, 20_000_000,
            &HgOpts { step_limit_f: 0.0, ..opts() },
        );
        for i in 0..3 {
            free = free.max((a_free.state.r[i] - base_free.state.r[p[i]]).norm());
            froz = froz.max((a_froz.state.r[i] - base_frozen.state.r[p[i]]).norm());
            hgw = hgw.max((h.state.r[i] - base_hg.state.r[p[i]]).norm());
        }
    }
    println!("Burrau to t = {t_max}, over all five non-identity label permutations:");
    println!("  AZ, reference re-chosen every boundary : {free:.3e}   switches {}", base_free.switches);
    println!("  AZ, reference FROZEN at index 0        : {froz:.3e}");
    println!("  Heggie (no reference to freeze)        : {hgw:.3e}");
    println!("  frozen/free = {:.3e}", froz / free.max(1e-300));
    assert!(hgw < 1e-12, "Heggie moved under relabelling: {hgw:e}");
    assert!(
        froz > free * 1e3,
        "freezing the reference did NOT make AZ label-dependent (frozen {froz:e} against free \
         {free:e}), so re-registration is not what erases Heggie's contrast and the explanation \
         on record for that negative is wrong"
    );
}
