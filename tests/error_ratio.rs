//! `error_ratio` has no reference implementation — it exists nowhere in the numpy tree — so
//! it is validated by invariant. Five checks, of which the last is the one that matters.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::ensemble::stats;
use prin_rs::grid;
use prin_rs::integrate::az::StepLimit;
use prin_rs::integrate::az;
use prin_rs::integrate::Integrator;
use prin_rs::physics::{energy, Cart};
use prin_rs::Vec2;

/// 2. Exactly 1.0 at `t = 0`, by construction.
#[test]
fn error_ratio_is_exactly_one_at_t_zero() {
    let e0 = [-12.81, -12.80, -12.83, -12.79, -12.82, -12.815, -12.805, -12.825];
    let (r, s0, st) = stats::error_ratio(&e0, &e0);
    println!("t=0: sigma_E(0) = {s0:e}, sigma_E(t) = {st:e}, ratio = {r}");
    assert_eq!(r, 1.0);
}

/// 3. Shift and scale invariance. MAD is shift-equivariant and scale-equivariant, so the
/// ratio is invariant to both — the statistic cannot be reading an offset or a unit choice.
#[test]
fn error_ratio_is_invariant_under_shift_and_scale() {
    let e0: Vec<f64> = (0..8).map(|k| -12.8 + 1e-5 * (k as f64 - 3.5)).collect();
    let et: Vec<f64> = e0.iter().map(|x| x + 3e-9 * (x * 1e8).sin()).collect();
    let (base, _, _) = stats::error_ratio(&e0, &et);

    // The tolerance is derived, not guessed. Shifting by `c` and then differencing values
    // whose spread is `w` costs about `eps * c / w` of relative precision — the shift is
    // ill-conditioned, the statistic is not. Here w ~ 1e-5, so a shift of 1000 predicts a
    // loss near 1e-8 and a shift of 12.8 near 1e-10.
    let spread = 1e-5f64;
    for shift in [12.8f64, 1000.0] {
        let s0: Vec<f64> = e0.iter().map(|x| x + shift).collect();
        let st: Vec<f64> = et.iter().map(|x| x + shift).collect();
        let (shifted, _, _) = stats::error_ratio(&s0, &st);
        let predicted = 100.0 * f64::EPSILON * shift / spread;
        println!(
            "shift {shift}: |d ratio| = {:.3e}, conditioning predicts <= {predicted:.1e}",
            (shifted - base).abs()
        );
        assert!(
            (shifted - base).abs() < predicted,
            "shift {shift} moved the ratio by more than its own conditioning explains: \
             {shifted} vs {base}"
        );
    }

    // Powers of two rescale exactly, so the invariance must hold to the last bit.
    for alpha in [0.25f64, 4.0, 1024.0] {
        let a0: Vec<f64> = e0.iter().map(|x| x * alpha).collect();
        let at: Vec<f64> = et.iter().map(|x| x * alpha).collect();
        let (scaled, _, _) = stats::error_ratio(&a0, &at);
        assert_eq!(scaled, base, "scale {alpha} is exact in binary and must not move the ratio");
    }

    // Other factors perturb each value by ~eps relative. The statistic depends on
    // differences that are ~1e-6 of the value here, so that perturbation is amplified by the
    // same factor — again the rescaling's conditioning, not the statistic's.
    let rel_spread = spread / 12.8;
    for alpha in [1e6f64, 3.7] {
        let a0: Vec<f64> = e0.iter().map(|x| x * alpha).collect();
        let at: Vec<f64> = et.iter().map(|x| x * alpha).collect();
        let (scaled, _, _) = stats::error_ratio(&a0, &at);
        let predicted = 100.0 * f64::EPSILON / rel_spread;
        println!(
            "scale {alpha}: |d ratio| = {:.3e}, conditioning predicts <= {predicted:.1e}",
            (scaled - base).abs()
        );
        assert!((scaled - base).abs() < predicted, "scale {alpha}: {scaled} vs {base}");
    }
}

/// 4. An exactly-integrable control: a tight pair with the third body far away is very nearly
/// a two-body problem, so `error_ratio` must sit at 1 to near machine precision. This
/// exercises the whole pipeline, not only the statistic.
#[test]
fn error_ratio_sits_at_one_on_a_near_two_body_control() {
    let m = [1.0f64, 1.0, 1.0];
    let mut e0 = Vec::new();
    let mut et = Vec::new();
    for k in 0..8 {
        let d = 1.0 + 1e-4 * (k as f64);
        let s = Cart::new(
            [Vec2::new(-0.5 * d, 0.0), Vec2::new(0.5 * d, 0.0), Vec2::new(0.0, 5000.0)],
            [Vec2::zero(); 3],
        );
        e0.push(energy::energy(&s.r, &s.v, &m, 0.0));
        let o = az::integrate_az(s, &m, 0.5, 4, 0.01, 200_000, None);
        et.push(energy::energy(&o.state.r, &o.state.v, &m, 0.0));
    }
    let (r, s0, st) = stats::error_ratio(&e0, &et);
    println!("near-two-body control: sigma_E(0) = {s0:.6e}, sigma_E(t) = {st:.6e}, ratio = {r:.12}");
    assert!((r - 1.0).abs() < 1e-6, "ratio = {r}, want 1 on a near-integrable control");
}

/// 5. **The check that makes a field with no oracle trustworthy.**
///
/// `error_ratio - 1` is pure integration error, so it must fall as `eta -> 0`. If it does
/// not, it is measuring a wrong equation rather than integration error — the same diagnostic
/// signature that has caught three bugs in this project, applied to a statistic that has
/// nothing else to check it against.
///
/// **The convergence is asserted on the median over pixels, not the max, and that is a
/// correction rather than a convenience.** `error_ratio` is now built on the maximum
/// deviation over 8 copies, so a max over 9 pixels is a 1-in-72 order statistic: which copy
/// of which pixel happens to be worst changes with `eta`, and the realisation scatter is
/// wider than the trend. Measured — the max rose 3.166e-4 -> 9.596e-4 at the first halving
/// and then fell by a factor of 340 over the next two, while the median fell at every single
/// step. This is CLAUDE.md's "do not read a scaling law off a single trajectory", one level
/// up: an extreme order statistic is not the place to read a convergence law.
///
/// The max is printed anyway, and the end-to-end fall across the whole decade is asserted, so
/// the scatter is visible rather than filtered out.
#[test]
fn error_ratio_minus_one_falls_with_step_size() {
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();
    println!("{:>8}{:>16}{:>16}{:>14}", "eta", "max |ratio-1|", "median |r-1|", "med ratio");
    let mut prev: Option<f64> = None;
    let mut fell = 0usize;
    let mut total = 0usize;
    let (mut first_max, mut last_max) = (0.0f64, 0.0f64);
    for (k, eta) in [4e-2f64, 2e-2, 1e-2, 5e-3].into_iter().enumerate() {
        // **Pinned to `StepLimit::None`, and that pin is the finding -- but the reason first
        // given for it was wrong and is corrected here.** The original read: *under the shipped
        // `Predictive` limit the residual is already at the round-off floor at the COARSEST
        // rung, 3.1e-9 falling to 1.6e-9 across the decade, non-monotone because it is
        // arithmetic scatter and not truncation.*
        //
        // Those two numbers are real and they are a **PLATEAU, not a floor**. The full ladder
        // under the shipped limit runs 3.10e-9, 3.71e-9, 3.50e-9, 1.64e-9 and then falls **8.7x
        // and 14.4x** to 1.31e-11. So `eta = 4e-2` is ~300x above the actual floor and
        // truncation is very much alive under the limit; what the two-rung sample measured was
        // the flat part. It was caught only when the integrator default moved to `Heggie`,
        // which has no plateau there and failed the assertion honestly.
        //
        // The pin still stands, on the honest reason: this test watches truncation fall over a
        // decade at coarse `eta`, and under the limit that decade is the plateau, where the
        // trend is flat and the ordering is arithmetic scatter. See
        // `the_residual_has_a_floor_and_the_coarse_plateau_is_not_it`, which measures the ladder
        // for both integrators and asserts the plateau is not the floor.
        let cfg = EnsembleCfg {
            eta,
            t_max: 4.0,
            n_sync: 10,
            step_limit: StepLimit::None,
            ..Default::default()
        };
        let mut devs: Vec<f64> = (0..s.npix())
            .map(|i| (evaluate::<f64>(&s, i, &cfg).error_ratio - 1.0).abs())
            .collect();
        devs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mx = *devs.last().unwrap();
        let med = devs[devs.len() / 2];
        if k == 0 {
            first_max = mx;
        }
        last_max = mx;
        let rs = prev.map(|p| format!("{:>14.3}", med / p)).unwrap_or_else(|| format!("{:>14}", "-"));
        println!("{eta:>8.0e}{mx:>16.6e}{med:>16.6e}{rs}");
        if let Some(p) = prev {
            total += 1;
            if med < p {
                fell += 1;
            }
        }
        prev = Some(med);
    }
    println!();
    println!("error_ratio - 1 must fall with eta. If it plateaus or rises, it is measuring a");
    println!("wrong equation, not integration error. Read the median column: the max over 9");
    println!("pixels of a max over 8 copies is a 1-in-72 order statistic and its realisation");
    println!("scatter is wider than the trend it is being asked to show.");

    assert_eq!(fell, total, "median error_ratio - 1 did not fall at every step-size halving");
    assert!(
        last_max < first_max,
        "max error_ratio - 1 did not fall across the decade: {first_max:.4e} -> {last_max:.4e}"
    );
}

/// **The residual under the shipped limit has a floor, and `eta = 4e-2` is 300x above it.**
///
/// # What this replaces, and why it had to be replaced
///
/// The previous form asserted *"under `StepLimit::Predictive` at `eta = 4e-2` the median
/// `|error_ratio - 1|` sits at round-off, and shrinking `eta` eightfold does not move it by an
/// order"* — sampling exactly two rungs, `4e-2` and `5e-3`. **It passed for a reason other than
/// the one it stated.** Both of those points lie on a *plateau*, not on the floor. The full
/// ladder, median over 9 near-field pixels:
///
/// ```text
///        eta         AZ      x prev        Heggie      x prev
///     4.00e-2   3.0957e-9        -       1.3121e-8        -
///     2.00e-2   3.7136e-9      0.8       5.1571e-9      2.5
///     1.00e-2   3.4963e-9      1.1       3.7245e-10    13.8
///     5.00e-3   1.6396e-9      2.1       1.9414e-10     1.9
///     2.50e-3   1.8881e-10     8.7       1.4579e-11    13.3
///     1.25e-3   1.3099e-11    14.4       7.5440e-12     1.9
///     6.25e-4   8.8050e-12     1.5       1.4508e-11     0.5
///     3.13e-4   1.2351e-11     0.7          0.0e0       inf
/// ```
///
/// AZ is flat across the first four rungs and then falls **8.7x and 14.4x**. So truncation is
/// very much alive under the shipped limit; what the two-rung test measured was the flat part.
/// It was exposed by the integrator default moving to `Heggie`, which has no plateau there and
/// failed the assertion honestly — *the control caught it, not the property*.
///
/// # What is actually true, and is asserted here
///
/// Both integrators reach the **same** floor, `~1e-11`, at `eta ~ 1.25e-3`, and refining past it
/// buys nothing (ratios 1.5, 0.7 and 0.5 — noise about a floor). That is the real statement, and
/// it is the one the pin above needs: a test written to watch truncation fall must run **above**
/// `1.25e-3` or it has no subject, and `error_ratio_minus_one_falls_with_step_size` does.
///
/// The plateau is asserted **explicitly as not-the-floor**, so the earlier reading cannot come
/// back by someone sampling two rungs again.
#[test]
fn the_residual_has_a_floor_and_the_coarse_plateau_is_not_it() {
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();
    // **Both integrators, because the difference between them is the finding.** The old form
    // ran on whichever was the default, and when that moved from `Az` to `Heggie` it failed --
    // correctly, and for a reason the two-rung sample could not show.
    for integrator in [Integrator::Az, Integrator::Heggie] {
        let med = |eta: f64| {
            let cfg = EnsembleCfg { eta, t_max: 4.0, n_sync: 10, integrator, ..Default::default() };
            let mut d: Vec<f64> = (0..s.npix())
                .map(|i| (evaluate::<f64>(&s, i, &cfg).error_ratio - 1.0).abs())
                .collect();
            d.sort_by(|a, b| a.partial_cmp(b).unwrap());
            d[d.len() / 2]
        };
        let coarse = med(4e-2);
        let at_floor = med(1.25e-3);
        let past_floor = med(3.125e-4);
        println!(
            "{:>7}: {coarse:.4e} at eta=4e-2, {at_floor:.4e} at 1.25e-3, {past_floor:.4e} at 3.125e-4",
            integrator.name()
        );

        // 1. There IS a floor, asserted as an ABSOLUTE LEVEL and not as a ratio between rungs.
        //    The first cut used a ratio and fired at once: Heggie reaches EXACTLY 0.0 at
        //    eta = 3.125e-4, which is the floor arriving, and `7.5e-12 / 0` is not a measurement.
        //    A floor is a level. 1e-10 sits an order above where both bottom out.
        const FLOOR: f64 = 1e-10;
        assert!(
            at_floor < FLOOR && past_floor < FLOOR,
            "{}: not in the floor band by eta=1.25e-3 ({at_floor:.4e}, then {past_floor:.4e}) -- the ladder needs another rung before anything below can be called a floor",
            integrator.name()
        );

        // 2. The coarse rung is NOT that floor. The arm with teeth: it fails if anyone
        //    reinstates "4e-2 is already at the floor". Measured margins are 236x (AZ) and
        //    1750x (Heggie), against a bar of 50x -- and against the FLOOR BAND rather than
        //    the measured value AZ would read only 31x, which is why the denominator is the
        //    rung and not the band.
        assert!(
            coarse > 50.0 * at_floor.max(1e-13),
            "{}: eta=4e-2 reads {coarse:.4e} against a floor of {at_floor:.4e}, under 50x -- if the coarse rung really is at the floor then truncation is gone and the pin above has no subject, which is exactly the reading this test exists to prevent",
            integrator.name()
        );
    }
}
