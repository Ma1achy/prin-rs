//! `error_ratio` has no reference implementation — it exists nowhere in the numpy tree — so
//! it is validated by invariant. Five checks, of which the last is the one that matters.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::ensemble::stats;
use prin_rs::grid;
use prin_rs::integrate::az;
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
#[test]
fn error_ratio_minus_one_falls_with_step_size() {
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();
    println!("{:>8}{:>16}{:>16}{:>14}", "eta", "max |ratio-1|", "median |r-1|", "ratio");
    let mut prev: Option<f64> = None;
    let mut fell = 0usize;
    let mut total = 0usize;
    for eta in [4e-2f64, 2e-2, 1e-2, 5e-3] {
        let cfg = EnsembleCfg { eta, t_max: 4.0, n_sync: 10, ..Default::default() };
        let mut devs: Vec<f64> = (0..s.npix())
            .map(|i| (evaluate::<f64>(&s, i, &cfg).error_ratio - 1.0).abs())
            .collect();
        devs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mx = *devs.last().unwrap();
        let med = devs[devs.len() / 2];
        let rs = prev.map(|p| format!("{:>14.3}", mx / p)).unwrap_or_else(|| format!("{:>14}", "-"));
        println!("{eta:>8.0e}{mx:>16.6e}{med:>16.6e}{rs}");
        if let Some(p) = prev {
            total += 1;
            if mx < p {
                fell += 1;
            }
        }
        prev = Some(mx);
    }
    println!();
    println!("error_ratio - 1 must fall with eta. If it plateaus or rises, it is measuring a");
    println!("wrong equation, not integration error.");
    assert_eq!(fell, total, "error_ratio - 1 did not fall at every step-size halving");
}
