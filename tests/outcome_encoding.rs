//! Step 5b's property tests.
//!
//! `r_coll`, the >=2-pair rule and the 3+2-bit encoding have **no reference implementation** —
//! none of them appear anywhere in the numpy tree — so there is nothing to compare against and
//! the only available validation is by property. Each test below states the property it is
//! protecting rather than pinning a number.
//!
//! The escape arm is the exception: it is ported from `reference/tb_all_az.py:59-75`, and
//! `escape_matches_the_legacy_classifier` checks the transcription against the arm of
//! `classify_legacy` that has the same origin.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::outcome::{self, Outcome, State};
use prin_rs::physics::{burrau, energy, newton, Cart};
use prin_rs::Vec2;

/// Every reachable `(state, detail)` survives a round trip through one byte.
#[test]
fn the_encoding_packs_into_one_byte_without_loss() {
    for bits in 0..6u8 {
        let st = State::from_bits(bits).unwrap();
        for detail in 0..4u8 {
            let o = Outcome::new(st, detail);
            let back = Outcome::unpack(o.pack()).unwrap();
            assert_eq!(back, o, "{}/{detail} did not survive packing", st.name());
        }
    }
    assert_eq!(State::from_bits(6), None, "6 and 7 are unassigned and must decode as invalid");
    assert_eq!(State::from_bits(7), None);
}

/// **The rule the >=2-pair count exists to enforce.**
///
/// By the triangle inequality, `|AB| < r_coll` and `|AC| < r_coll` force `|BC| < 2 r_coll`, so
/// "exactly two pairs below threshold" is a reachable state and is already a near-triple.
/// Requiring all three would label it an ordinary binary collision. The test builds exactly
/// that configuration and asserts it is never given a single-pair detail.
#[test]
fn two_pairs_below_r_coll_is_never_an_ordinary_binary_collision() {
    let r_coll = 1e-2f64;
    // A and B close, A and C close, so |BC| is between the two and below 2 r_coll but above
    // r_coll — the case that is reachable and that "all three" would miss.
    let r = [
        Vec2::new(0.0, 0.0),
        Vec2::new(0.6 * r_coll, 0.0),
        Vec2::new(-0.6 * r_coll, 0.0),
    ];
    let d = newton::pair_dists(&r);
    println!("|AB| = {:.4e}  |AC| = {:.4e}  |BC| = {:.4e}   r_coll = {r_coll:.4e}", d[0], d[1], d[2]);
    assert!(d[0] < r_coll && d[1] < r_coll, "setup: two pairs must be below r_coll");
    assert!(d[2] > r_coll, "setup: the third pair must be above, or this is a plain triple");
    assert!(d[2] < 2.0 * r_coll, "the triangle inequality bound must hold");

    let mask = outcome::collision_pairs(&r, r_coll);
    assert_eq!(mask.count_ones(), 2, "exactly two pairs should register");
    let detail = outcome::collision_detail(mask);
    assert_eq!(detail, 3, "two pairs below r_coll must encode as detail = 3, all three");

    // And the opposite direction: exactly one pair below must keep its own index.
    for k in 0..3u8 {
        assert_eq!(outcome::collision_detail(1 << k), k, "a single pair keeps its own index");
    }
    println!();
    println!("Two pairs below r_coll encodes as detail = 3 (triple), not as pair {} — which",
             mask.trailing_zeros());
    println!("is what an 'all three' rule would have produced.");
}

/// The mapping from AZ's `(R1, R2, R3)` back to pair indices, over every reference body.
///
/// `R1 = r_b - r_a`, `R2 = r_c - r_a`, `R3 = r_c - r_b`, and the reference body changes
/// between sync boundaries. Getting this wrong attributes a collision to the wrong pair while
/// every magnitude stays correct — silent, and exactly the failure class this project keeps
/// running into.
#[test]
fn the_collision_mask_labels_the_right_pair_for_every_reference_body() {
    let r = [Vec2::new(0.0, 0.0), Vec2::new(1.0, 0.0), Vec2::new(0.0, 1.0)];
    let d_true = newton::pair_dists(&r);
    for a in 0..3usize {
        let (b, c) = {
            let o: Vec<usize> = (0..3).filter(|&k| k != a).collect();
            (o[0], o[1])
        };
        let d1 = (r[b] - r[a]).norm();
        let d2 = (r[c] - r[a]).norm();
        let d3 = (r[c] - r[b]).norm();
        // A threshold just above the shortest pair picks out exactly that pair.
        let shortest = d_true.iter().cloned().fold(f64::INFINITY, f64::min);
        let mask = outcome::collision_pairs_from(a, b, c, d1, d2, d3, shortest * 1.001);
        let expect = outcome::collision_pairs(&r, shortest * 1.001);
        assert_eq!(mask, expect, "reference body {a}: AZ-labelled mask disagrees with the direct one");
        println!("reference body {a}: mask {mask:03b} matches the direct pair_dists mask");
    }
}

/// The ported escape arm and `classify_legacy`'s escape arm come from the same reference
/// lines, so they must agree on every state they are both shown.
///
/// **The arm does not fire at `t = 13`, and under the corrected step control it does not fire at
/// `t = 20` either.** An arm that is never exercised is not tested by agreeing with another arm
/// that is also never exercised, so the horizon here is chosen to make it fire and the assertion
/// at the bottom is what enforces that.
///
/// **The `t = 20` figure on record was measuring the `dtau` bug.** With `dtau` sized once per
/// interval, 109 of 1024 near-field pixels fired at `t_max = 20`; with it recomputed per step,
/// **zero** do. Those 109 are not a chaotic reshuffle -- their median energy drift is **1.147**,
/// which is 115% of the total energy, against **6.2e-5** for the 915 that stay silent, and under
/// the fix the same pixels sit at 1.6e-3 and do not escape. A giant post-encounter step throws a
/// body outward, it reads as unbound and receding at the next boundary, and the arm latches.
/// Genuine escape is not suppressed: at `t = 40` the two modes give 280 and 308 of 1024.
///
/// So this runs at `t = 40`, where 7 of 25 fire with a worst drift of 7.5e-6.
///
/// **`n_sync` is scaled with `t_max`**, not held at 32. `dtau = eta*dt_left/(A0*B0)`, so a fixed
/// `n_sync` across two horizons compares two discretisations; the old form ran `t = 13` at a
/// 0.406 interval and `t = 20` at 0.625.
///
/// **Two properties of the ported arm, transcribed rather than fixed.** It is sampled at sync
/// boundaries, so an escape's `t_end` has the resolution of the sync grid. And it latches on
/// first firing, so a body can satisfy "unbound and receding" at one boundary and be recaptured
/// by the next.
#[test]
fn escape_matches_the_legacy_classifier() {
    let s = grid::region("near-field", 5, 5, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    // The sync interval every horizon is held to, so the rows are one discretisation.
    const DT_SYNC: f64 = 13.0 / 32.0;
    let mut fired = 0usize;
    let mut checked = 0usize;
    for t_max in [13.0f64, 40.0] {
        let n_sync = (t_max / DT_SYNC).round() as usize;
        for i in 0..s.npix() {
            let o = az::integrate_az_opts(
                s.nominal::<f64>(i), &m, t_max, n_sync, 0.01, 30_000,
                &AzOpts { stop_on_event: false, ..Default::default() },
            );
            let legacy = outcome::classify_legacy(&o.state, &m);
            let ported = outcome::escape_candidate(&o.state, &m).unwrap_or(3);
            assert_eq!(ported, legacy, "pixel {i} at t={t_max}: ported escape arm disagrees with tb.classify");
            if o.events.escape.is_some() {
                fired += 1;
            }
            checked += 1;
        }
    }
    println!("{checked} states: the ported escape test and tb.classify agree everywhere");
    println!("the arm fired at a sync boundary on {fired} of them");
    assert!(fired > 0, "the escape arm never fired, so the agreement above proves nothing");
}

/// The two `dtau` modes must **disagree** where `A*B` moves across an interval and **agree to
/// round-off** where it does not.
///
/// "The modes differ somewhere" is the test that cannot fail: they are different expressions, so
/// of course they do. The form with teeth is the pair -- a case that separates them and a case
/// that does not -- because a fix that changed every trajectory would pass the first arm exactly
/// as well as a correct one, and so would a no-op that changed none.
///
/// The separating case is `deep interior`, which is collision-rich, so encounters land inside
/// intervals and `A*B` swings by orders. The control is a **long sync interval on a tame
/// trajectory**: near-field at `t = 1` with `n_sync = 1`, where the geometry barely moves, the
/// recomputed `A*B` stays near its entry value, and the cap in `PerStepInterval` is active
/// throughout -- so the two modes are the same arithmetic and the states agree to round-off.
#[test]
fn the_dtau_modes_separate_where_ab_moves_and_agree_where_it_does_not() {
    use prin_rs::integrate::az::DtauMode;
    let m = burrau::masses::<f64>();
    let run = |s, t: f64, n_sync, mode| {
        az::integrate_az_opts(
            s, &m, t, n_sync, 0.01, 30_000,
            &AzOpts { stop_on_event: false, dtau_mode: mode, ..Default::default() },
        )
    };
    let dev = |a: &az::AzOut<f64>, b: &az::AzOut<f64>| {
        (0..3).map(|k| (a.state.r[k] - b.state.r[k]).norm()).fold(0.0f64, f64::max)
    };

    // Separating: encounters inside intervals, so `A*B` swings and the step control matters.
    let di = grid::region("deep interior", 4, 4, 0.05).unwrap();
    let mut worst = 0.0f64;
    for i in 0..di.npix() {
        let s0 = di.nominal::<f64>(i);
        let a = run(s0, 13.0, 32, DtauMode::FixedPerInterval);
        let b = run(s0, 13.0, 32, DtauMode::PerStepInterval);
        if a.finite && b.finite {
            worst = worst.max(dev(&a, &b));
        }
    }
    assert!(worst > 1e-3, "the modes agree on deep interior (worst dev {worst:e}); the step \
                           control cannot be doing anything there and the fix is a no-op");

    // Control: one long interval on a tame trajectory. `A*B` barely moves, the cap binds, and
    // the two modes are the same arithmetic.
    let nf = grid::region("near-field", 3, 3, 0.05).unwrap();
    for i in 0..nf.npix() {
        let s0 = nf.nominal::<f64>(i);
        let a = run(s0, 1.0, 1, DtauMode::FixedPerInterval);
        let b = run(s0, 1.0, 1, DtauMode::PerStepInterval);
        let d = dev(&a, &b);
        assert!(d < 1e-9, "pixel {i}: the modes disagree by {d:e} on a tame interval, where \
                           `A*B` does not move and they should be the same arithmetic");
    }
}

/// Rescale by `alpha`, rescale time by `alpha^{3/2}`, and both the **outcome** and `t_end/
/// alpha^{3/2}` must come back the same.
///
/// This is the check that catches an absolute length leaking into `r_coll` — the failure mode
/// that measured 1.66x in prior work. It is stronger than the `shape_vec` version because an
/// outcome is a *discrete* label: a leak does not perturb it slightly, it flips it.
#[test]
fn outcomes_and_t_end_are_invariant_under_the_scale_symmetry() {
    let s = grid::region("near-field", 5, 5, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    // t_max = 13 with n_sync = 32, not a shorter horizon: at t_max = 6 nothing on this grid
    // fires, every t_end is exactly the horizon, and the test passes while measuring nothing.
    // The assertion below on `terminated` is there so that can never silently recur.
    let t_max = 13.0f64;
    println!("{:>8}{:>12}{:>24}{:>16}", "alpha", "outcomes", "t_end/alpha^1.5", "max |dt| rel");
    let run = |i: usize, alpha: f64| -> (u8, f64) {
        let c = s.nominal::<f64>(i);
        let scaled = Cart::new(
            [c.r[0] * alpha, c.r[1] * alpha, c.r[2] * alpha],
            [c.v[0] / alpha.sqrt(), c.v[1] / alpha.sqrt(), c.v[2] / alpha.sqrt()],
        );
        let o = az::integrate_az_opts(
            scaled, &m, t_max * alpha.powf(1.5), 16, 0.01, 30_000,
            &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
        );
        let out = outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted);
        (out.pack(), o.t_end / alpha.powf(1.5))
    };

    let base: Vec<(u8, f64)> = (0..s.npix()).map(|i| run(i, 1.0)).collect();
    let terminated = base.iter().filter(|(_, te)| *te < t_max * (1.0 - 1e-12)).count();
    // Report the t_end column for a pixel that actually terminates; a censored pixel's t_end
    // is just the horizon and would show invariance that the rescaling gives for free.
    let probe = base
        .iter()
        .position(|(_, te)| *te < t_max * (1.0 - 1e-12))
        .unwrap_or(0);
    println!("{terminated} of {} pixels terminate before the horizon; showing pixel {probe}",
             s.npix());
    assert!(
        terminated > 0,
        "no pixel terminated early, so this test would pass without exercising t_end at all"
    );
    let mut worst_rel = 0.0f64;
    for alpha in [0.25f64, 4.0, 3.7, 1.0 / 3.0] {
        let mut flips = 0usize;
        let mut wr = 0.0f64;
        for i in 0..s.npix() {
            let (o, te) = run(i, alpha);
            if o != base[i].0 {
                flips += 1;
            }
            let rel = (te - base[i].1).abs() / base[i].1.abs().max(1e-30);
            wr = wr.max(rel);
        }
        println!("{alpha:>8.4}{:>12}{:>24.14e}{wr:>16.3e}",
                 format!("{}/{} same", s.npix() - flips, s.npix()), run(probe, alpha).1);
        assert_eq!(flips, 0, "alpha {alpha}: {flips} outcome labels flipped — a length has leaked in");
        if alpha == 0.25 || alpha == 4.0 {
            assert_eq!(wr, 0.0, "alpha {alpha} is exact in binary; t_end must be bitwise identical");
        }
        worst_rel = worst_rel.max(wr);
    }
    println!();
    println!("No outcome label moved. t_end agrees bitwise for the binary-exact factors and to");
    println!("{worst_rel:.3e} relative for the inexact ones, which is roundoff in the rescaling.");
    assert!(worst_rel < 1e-10, "t_end scale covariance broken at {worst_rel:e}");
}

/// BRIEF §2.6, and **the brief is wrong about this pixel.**
///
/// §2.6 says `deep interior` "drives all three bodies together", is not regularisable, will
/// fail however well the integrator is built, and should hit the triple-collision outcome —
/// citing 190 s per probe and still failing.
///
/// Measured, both here and in the numpy reference: it is an **ordinary binary encounter
/// between bodies 0 and 2**. The initial separations are 2.236, 1.414 and 3.0; pair (0,2)
/// reaches `1.67e-5 R` while pairs (0,1) and (1,2) never come within `R` of each other at
/// all — checked by sweeping `r_coll` up to `R` itself, where they still do not register. It
/// integrates to `t = 13` in about a second with `|dE/E| ~ 1.4e-7` and two reference
/// switches; the reference gives `dmin = 2.2976e-5`, `drift = 1.3936e-7`, `switches = 2` on
/// the same initial condition.
///
/// The 190 s failure §2.6 records is almost certainly the **unregularised** integrator. A
/// close binary approach with a distant third body is the exact case AZ exists to handle.
///
/// So this test asserts the hazard that is real — bounded wall-clock and a finite result,
/// BRIEF §5's third signature — and *records* the label rather than asserting a triple that
/// the geometry does not support. See `examples/deep_interior.rs` for the full probe.
#[test]
fn deep_interior_terminates_in_bounded_wall_clock() {
    let s = grid::region("deep interior", 1, 1, 1e-6).unwrap();
    let m = burrau::masses::<f64>();
    let t0 = std::time::Instant::now();
    let o = az::integrate_az_opts(
        s.nominal::<f64>(0), &m, 13.0, 32, 0.01, 30_000,
        &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
    );
    let out = outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted);
    let dt = t0.elapsed().as_secs_f64();
    println!("deep interior: {} detail {} at t_end {:.6}, {:.3} s wall, {} steps",
             out.state.name(), out.detail, o.t_end, dt, o.steps);
    if let Some((mask, tc)) = o.events.collision {
        let pairs: Vec<String> = (0..3)
            .filter(|k| mask & (1 << k) != 0)
            .map(|k| format!("{:?}", prin_rs::physics::PAIRS[k]))
            .collect();
        println!("  fired on pair(s) {} at t = {tc:.6}", pairs.join(", "));
    }
    println!();
    println!("BRIEF §2.6 predicts a triple collision here. It is a binary encounter between");
    println!("bodies 0 and 2; the other two pairs never register even at r_coll = R.");

    assert!(dt < 30.0, "took {dt:.1} s — this is the hang BRIEF §5 signature 3 warns about");
    assert!(o.finite, "the trajectory went non-finite");
    assert_eq!(out.state, State::Collision, "expected a collision outcome");
    assert_eq!(out.detail, 1, "the colliding pair is (0,2), index 1");
}

/// `r_coll` is a fraction of the **initial** hyperradius, evaluated once. A co-moving length
/// would make the threshold move with the system and destroy the fixed-at-t=0 guarantee.
#[test]
fn r_coll_is_fixed_at_t_zero_and_not_co_moving() {
    let s = grid::region("near-field", 1, 1, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let c = s.nominal::<f64>(0);
    let r0 = energy::hyperradius(&c.r, &m);
    let frac = 1e-3f64;

    let o = az::integrate_az_opts(
        c, &m, 13.0, 32, 0.01, 30_000,
        &AzOpts { r_coll_frac: frac, stop_on_event: false, ..Default::default() },
    );
    let r_final = energy::hyperradius(&o.state.r, &m);
    println!("R(0) = {r0:.6}, R(t_max) = {r_final:.6}, ratio {:.4}", r_final / r0);
    println!("r_coll = {:.6e} throughout, never {:.6e}", frac * r0, frac * r_final);
    assert!(
        (r_final / r0 - 1.0).abs() > 1e-3,
        "the system barely moved; this test cannot distinguish fixed from co-moving"
    );
    // The threshold that fired is the t=0 one: d_min_true must be compared against frac*R(0).
    // If it were co-moving the comparison below would be against a different number entirely.
    println!("d_min_true = {:.6e}", o.d_min_true);
}

/// The whole pipeline, at the pixel level: outcomes are produced, packed and spread.
#[test]
fn outcome_fields_populate_across_a_grid() {
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    let cfg = EnsembleCfg::default();
    let px: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect();
    let mut counts = [0usize; 6];
    for p in &px {
        counts[p.state as usize] += 1;
        assert_eq!(p.outcome, Outcome::new(State::from_bits(p.state).unwrap(), p.detail).pack());
        assert!(p.t_end > 0.0 && p.t_end <= cfg.t_max, "t_end {} out of range", p.t_end);
    }
    for k in 0..6 {
        if counts[k] > 0 {
            println!("{:>14}: {:>4} of {}", State::from_bits(k as u8).unwrap().name(), counts[k], px.len());
        }
    }
    let disagree = px.iter().filter(|p| p.n_outcome_disagree > 0).count();
    println!("pixels whose copies disagree on the outcome: {disagree} of {}", px.len());
}

/// `escape_every` must change an escape-terminated `t_end` and leave everything else alone.
///
/// **A flag that does nothing passes as easily as one that works**, so both arms are asserted:
/// the reference cadence must be reproduced bitwise at `0`, and the fine cadence must actually
/// move `t_end` off the sync boundary on a trajectory that escapes.
///
/// The mechanism it exists for: collision is sampled inside the RK4 loop and already carries
/// step resolution; **escape is sampled only at sync boundaries**, so an escape-terminated
/// `t_end` takes at most `n_sync` values across a whole chart and renders as concentric contour
/// bands. Measured on `preset_plambda` at 64²: `t_end` distinct **16 -> 2623**, and the fraction
/// landing exactly on a boundary **99.52% -> 0.26%**.
#[test]
fn escape_every_moves_t_end_off_the_sync_boundary_and_is_inert_at_zero() {
    use prin_rs::integrate::az::{integrate_az_opts, AzOpts};
    use prin_rs::physics::decoder;

    let n_sync = 32usize;
    let t_max = 13.0f64;
    let dt_sync = t_max / n_sync as f64;
    let opts = |ev: usize| AzOpts::<f64> {
        dtau_mode: prin_rs::integrate::az::DtauMode::default(),
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: 1e-3,
        stop_on_event: true,
        escape_every: ev,
        // Off: this test asserts the CADENCE moves `t_end`, and the guard would suppress the
        // in-loop detections it is asserting on. The guard has its own test.
        escape_confirm: false,
        // The numpy reference's ungated escape test: every result in this diagnostic
        // predates the distance gate and is quoted against that form.
        escape_rule: outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_escape: true,
        keep_boundary_shapes: false,
    };

    // A latent-chart configuration, because that is where escape terminates: Burrau's near-field
    // has a silent escape arm at t = 13 and could not exercise this at all.
    let mut moved = 0usize;
    let mut escaped = 0usize;
    let mut checked = 0usize;
    for i in 0..9usize {
        for j in 0..9usize {
            let u = -1.0 + 2.0 * i as f64 / 8.0;
            let v = -1.0 + 2.0 * j as f64 / 8.0;
            // A momentum plane: the configuration coordinates are held at zero and only the
            // momentum ones move, which is what makes essentially every trajectory escape --
            // and escape is the only arm this flag touches.
            let z = decoder::Latent {
                z_alpha: 0.0,
                z_beta: 0.0,
                z_q: [u, v, 0.0, 0.0],
                z_mu: [0.0, 0.0],
            };
            let d = decoder::decode(&z);
            let s0 = d.ic.s;
            let ic = d.ic;
            checked += 1;

            let a = integrate_az_opts(s0, &ic.m, t_max, n_sync, 1e-2, 200_000, &opts(0));
            let b = integrate_az_opts(s0, &ic.m, t_max, n_sync, 1e-2, 200_000, &opts(1));

            if a.events.escape.is_some() {
                escaped += 1;
                // At the reference cadence an escape time IS a boundary time, exactly.
                let t = a.events.escape.unwrap().1;
                let k = (t / dt_sync).round();
                assert!(
                    (t - k * dt_sync).abs() <= 1e-9 * t_max,
                    "escape at the reference cadence must land on a sync boundary, got {t}"
                );
                if (b.t_end - a.t_end).abs() > 1e-9 {
                    moved += 1;
                }
            }
        }
    }
    assert!(checked > 0, "no configuration decoded, so nothing was tested");
    assert!(
        escaped > 0,
        "no trajectory escaped, so the flag's only subject never executed -- this test would \
         pass without the feature existing"
    );
    assert!(
        moved > 0,
        "{escaped} trajectories escaped and not one had its t_end moved by the fine cadence, \
         so the flag is inert"
    );
}


/// The persistence guard must reject transients and keep genuine escapes, and BOTH arms matter.
///
/// **Measured basis, not a precaution.** `escape_candidate` is relative energy `> 0` and
/// receding, which during a close encounter is transiently true. In `deep interior`, of the 895
/// trajectories that escape under `escape_every = 1` and not at the reference cadence, **0.000
/// are still unbound one boundary later** — and 0.000 at +2, +3, +4 and +8. Latching them took
/// the escape fraction from 0.0947 to 0.5494.
///
/// So the guard must **cut** the count somewhere collisions are common, and must **not** cut it
/// where escape genuinely terminates. A guard that rejects everything passes the first arm as
/// easily as a correct one; the second arm is what distinguishes them.
#[test]
fn escape_confirm_cuts_transients_and_keeps_genuine_escapes() {
    use prin_rs::grid::{self, Chart};
    use prin_rs::integrate::az::{integrate_az_opts, AzOpts};

    let (t_max, n_sync) = (13.0f64, 32usize);
    let opts = |confirm: bool| AzOpts::<f64> {
        dtau_mode: prin_rs::integrate::az::DtauMode::default(),
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: 1e-3,
        stop_on_event: true,
        escape_every: 1,
        escape_confirm: confirm,
        // The numpy reference's ungated escape test: every result in this diagnostic
        // predates the distance gate and is quoted against that form.
        escape_rule: outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_escape: true,
        keep_boundary_shapes: false,
    };
    let count = |chart: &Chart, body: usize, cx: f64, cy: f64, half: f64, confirm: bool| {
        let n = 16usize;
        let mut esc = 0usize;
        for i in 0..n {
            for j in 0..n {
                let u = cx - half + 2.0 * half * (i as f64 + 0.5) / n as f64;
                let v = cy - half + 2.0 * half * (j as f64 + 0.5) / n as f64;
                let ic = grid::decode_state(chart, body, u, v);
                let o = integrate_az_opts(
                    ic.s, &ic.m, t_max, n_sync, 1e-2, 200_000, &opts(confirm),
                );
                if o.events.escape.is_some() {
                    esc += 1;
                }
            }
        }
        esc
    };

    // Arm 1: collision-rich. The guard must cut the count.
    let (dirty, clean) = (
        count(&Chart::BodyPlane, 0, 0.0, 0.0, 0.05, false),
        count(&Chart::BodyPlane, 0, 0.0, 0.0, 0.05, true),
    );
    assert!(
        clean < dirty,
        "`deep interior` unguarded {dirty} escapes, guarded {clean} -- the guard cut nothing, \
         so either the transients are gone or it is not running"
    );

    // Arm 2: escape genuinely terminates. The guard must NOT cut the count, or it is simply
    // rejecting escapes rather than rejecting transients.
    let plambda = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == "preset_plambda")
        .expect("preset_plambda is in the gallery");
    let (pd, pc) = (
        count(&plambda.1, 0, plambda.2, plambda.3, plambda.4, false),
        count(&plambda.1, 0, plambda.2, plambda.3, plambda.4, true),
    );
    assert!(
        pd > 0,
        "no escapes on `preset_plambda` at all, so arm 2's subject never executed"
    );
    assert_eq!(
        pd, pc,
        "the guard removed {} genuine escapes on `preset_plambda`, where escape is the \
         terminating event and the finer stride adds none -- it is rejecting escapes, not \
         transients",
        pd - pc
    );
}

/// The escape **distance gate**, with the arm that says it does not reject real escapes.
///
/// Two hand-built states, one masses. A body sitting close to a tight pair but momentarily
/// unbound and receding — the mid-encounter transient the ungated test latches — and the same
/// body far away on the same kind of orbit. The ungated test cannot tell them apart; that is the
/// defect. The gated one must reject the first and **keep the second**.
///
/// *A guard needs the arm that says it did not cut too much.* A gate that refused everything
/// would pass the first assertion exactly as well as a correct one.
#[test]
fn escape_distance_gate_rejects_transients_and_keeps_real_escapes() {
    let m = [3.0f64, 4.0, 5.0];
    // Bodies 0 and 1 are the tight pair, so body 2 is the candidate under either rule.
    let build = |d: f64| {
        let mut s = Cart::<f64>::default();
        s.r = [
            Vec2::new(-0.1, 0.0),
            Vec2::new(0.1, 0.0),
            Vec2::new(d, 0.0),
        ];
        // Fast and outward: unbound relative to the pair's barycentre at either separation.
        s.v = [Vec2::new(0.0, 0.0), Vec2::new(0.0, 0.0), Vec2::new(8.0, 0.0)];
        s
    };
    let near = build(0.5);
    let far = build(40.0);

    // Ungated — the numpy reference's form. It fires on BOTH, which is the fault.
    assert_eq!(outcome::escape_candidate(&near, &m), Some(2),
               "the ungated test is supposed to fire mid-encounter; if it does not, this test \
                is not exercising the defect it exists for");
    assert_eq!(outcome::escape_candidate(&far, &m), Some(2));

    // Gated at 5R. **`R` is fixed at `t = 0`**, from the compact configuration both states
    // start in — the driver forms it once from `s0` and never recomputes it. Taking it from the
    // instantaneous state instead makes the gate co-moving, and a co-moving length grows with
    // the very separation it is meant to bound: at `d = 40` its own `R` is 12.5, so `5R = 62.5`
    // and the gate rejects the escape it exists to admit. That is the same defect
    // `r_coll`/`epsilon` are canonical to avoid, and it fired here first.
    let r_scale = energy::hyperradius(&near.r, &m);
    let gated = outcome::EscapeRule::Distance(5.0);
    assert_eq!(
        outcome::escape_candidate_rule(&near, &m, gated, r_scale, None),
        None,
        "the distance gate must reject a body still deep inside the system"
    );
    assert_eq!(
        outcome::escape_candidate_rule(&far, &m, gated, r_scale, None),
        Some(2),
        "the distance gate must NOT reject a genuine, far escape — without this arm a gate that \
         refuses everything passes"
    );

    // `r_esc = 0` is the reference path and must be bit-for-bit the ungated test.
    for s in [&near, &far] {
        assert_eq!(
            outcome::escape_candidate_rule(s, &m, outcome::EscapeRule::Reference, r_scale, None),
            outcome::escape_candidate(s, &m),
            "EscapeRule::Reference is the numpy form and must not diverge"
        );
    }
}


// ---------------------------------------------------------------------------------------
// The closure-and-energy escape criterion.
// ---------------------------------------------------------------------------------------

/// Neither arm alone is sufficient, and **each fixture is scored at maximum by the other arm**.
///
/// The two failure modes are complementary and both are on record: closure alone reads 82.8%
/// precision because *settling* also happens for bound hierarchies; energy alone reads 97.9%
/// because it **flickers** during encounters. So a test that only showed a real escape firing
/// would pass under either arm on its own and prove nothing about the conjunction.
///
/// - `escaper` — both hold, and only body 2 is unbound. Fires, on body 2.
/// - `bound` — closure holds, energy does not. The closure-alone failure mode.
/// - `transient` — energy holds, closure does not. The **old criterion's** failure mode, and the
///   test asserts the old criterion fires on it, so the fixture is known to be live.
/// - `dispersing` — see below. The reference's body ordering, transcribed and recorded.
#[test]
fn closure_needs_both_arms_and_each_fixture_isolates_one() {
    let rule = outcome::EscapeRule::Closure(outcome::CLOSURE_TAU);
    let settled = Some(1e-6);
    let moving = Some(1e-1);

    // Bodies 0 and 1 are a tight bound pair in mutual orbit; body 2 is light, far and receding.
    // Chosen so that body 2 is the ONLY unbound body -- with a heavy escaper the (1,2) barycentre
    // runs off with it and body 0 comes out unbound too, which is the `dispersing` case below.
    let hier_m = [3.0f64, 4.0, 0.05];
    let mut escaper = Cart::<f64>::default();
    escaper.r = [Vec2::new(-0.4, 0.0), Vec2::new(0.3, 0.0), Vec2::new(8.0, 0.0)];
    escaper.v = [
        Vec2::new(0.0, -4.0 / 7.0),
        Vec2::new(0.0, 3.0 / 7.0),
        Vec2::new(2.0, 0.0),
    ];
    assert!(outcome::unbound(&escaper, &hier_m, 2), "fixture is not exercising the energy arm");
    for b in [0, 1] {
        assert!(!outcome::unbound(&escaper, &hier_m, b),
                "body {b} must be BOUND here, or the returned label says nothing about ordering");
    }
    assert_eq!(outcome::escape_candidate_rule(&escaper, &hier_m, rule, 1.0, settled), Some(2));

    let m = [3.0f64, 4.0, 5.0];
    // A wide but BOUND hierarchy whose shape is settling just as hard as a real escape's.
    let mut bound = Cart::<f64>::default();
    bound.r = [Vec2::new(-0.1, 0.0), Vec2::new(0.1, 0.0), Vec2::new(40.0, 0.0)];
    bound.v = [Vec2::new(0.0, 0.0), Vec2::new(0.0, 0.0), Vec2::new(0.02, 0.0)];
    for b in 0..3 {
        assert!(!outcome::unbound(&bound, &m, b),
                "the bound fixture must be bound on every body, or it tests nothing");
    }
    assert_eq!(outcome::escape_candidate_rule(&bound, &m, rule, 1.0, settled), None,
               "a settled BOUND hierarchy must not fire -- this is closure-alone's 82.8%");

    // Deep inside the system and moving fast: unbound this instant, shape still swinging.
    let mut transient = bound;
    transient.r[2] = Vec2::new(0.5, 0.0);
    transient.v[2] = Vec2::new(8.0, 0.0);
    assert!(outcome::unbound(&transient, &m, 2), "fixture is not exercising the energy arm");
    assert_eq!(outcome::escape_candidate_rule(&transient, &m, rule, 1.0, moving), None,
               "an unbound instant with a moving shape must not fire -- energy-alone's 97.9%");

    // ...and the fixture the closure arm rejects is one the OLD criterion accepts. Without this
    // the transient case could be inert geometry rather than the defect.
    assert_eq!(outcome::escape_candidate(&transient, &m), Some(2),
               "the mid-encounter fixture must fire under the ungated test, or this test is not \
                exercising the defect the criterion exists to fix");

    // **THE LABEL IS THE LOWEST FIRING INDEX, NOT THE ESCAPING BODY, AND ON A DISPERSING SYSTEM
    // THOSE DIFFER.** With a heavy third body receding, the (1,2) barycentre travels with it and
    // body 0 comes out unbound relative to that barycentre as well -- all three fire, and the
    // reference's `b = np.argmax(fire, -1)` returns 0. Checked against the reference itself:
    // `E = [9.474, 12.178, 31.825]`, `argmax = 0`.
    //
    // Transcribed rather than corrected. It matters because `detail` is the escaping-body label
    // classification and rendering read, so how often it disagrees with the physically escaping
    // body is a **measurement**, not something to quietly patch here.
    let mut dispersing = bound;
    dispersing.v[2] = Vec2::new(8.0, 0.0);
    for b in 0..3 {
        assert!(outcome::unbound(&dispersing, &m, b), "every body should be unbound here");
    }
    assert_eq!(outcome::escape_candidate_rule(&dispersing, &m, rule, 1.0, settled), Some(0),
               "the reference returns the lowest firing index; a tightest-pair-first ordering \
                would say 2 here");

    // An unfilled window is `None`, and **nothing can fire** -- not even the genuine escaper.
    // NaN, never 0: a zero would read as perfectly settled and fire on everything at t ~ 0.
    assert_eq!(outcome::escape_candidate_rule(&escaper, &hier_m, rule, 1.0, None), None);
    assert_eq!(outcome::escape_candidate_rule(&escaper, &hier_m, rule, 1.0, Some(f64::NAN)), None,
               "NaN closure is undetermined, and undetermined has not settled");
}

/// **The two-end chord cannot tell a full revolution from stationarity.** Transcribed, and a
/// limitation to report rather than a bug to fix.
///
/// The reference buffers `nbuf` samples and reads only `buf[-1]` and `buf[0]`, so a shape vector
/// that goes all the way round and comes back reads as perfectly settled. A circular inner binary
/// rotates `(2p/I, 2q/I)` at its orbital frequency while `n_0` stays put, so a window
/// commensurate with that period aliases the orbit into a closure of zero.
///
/// This is why `examples/escape_closure.rs` §0 measures the tightest pair's period against the
/// realised window in every region **before** any distribution is read. The energy arm covers
/// much of it -- a bound hierarchy has no unbound body -- but not the transient case.
///
/// Would fire if `closure` ever became a max over the window instead of a chord.
#[test]
fn a_full_revolution_aliases_to_zero_closure() {
    // A shape vector on the (n1, n2) circle at fixed n0, sampled at phase 0, half a turn and a
    // full turn -- what a circular inner binary does to the Hopf map.
    let at = |phase: f64| {
        let (s, c) = phase.sin_cos();
        [0.30, 0.954_f64 * c, 0.954 * s]
    };
    let start = at(0.0);
    let half = outcome::closure(&at(std::f64::consts::PI), &start);
    let full = outcome::closure(&at(2.0 * std::f64::consts::PI), &start);

    assert!(half > 1.0, "half a revolution must read as a large excursion, got {half:e}");
    assert!(full < 1e-12,
            "a FULL revolution reads as settled -- this is the aliasing the window has, and the \
             period is measured against the window for that reason. got {full:e}");
}

// generated from `reference/escape_criterion.py` (`_rel` plus the criterion line) over 40 random
// states at tau = 1e-3; do not hand-edit. `v = p/m`, matching the reference's own conversion.
const GOLDEN: [([f64; 3], [f64; 6], [f64; 6], f64, i8); 40] = [
    ([8.21425506922998983e-01, 1.74819465610028746e+00, 2.00374589405839387e+00], [-3.77048793302444363e+00, -2.81659132338035256e+00, 3.42568818368295602e+00, -3.43663539076642532e+00, -2.96180840480561614e+00, 3.58662762633420051e+00], [1.18704463661421977e+00, -1.27589781584408191e+00, 5.21224419192332644e-02, 7.45193686291454105e-01, -8.97084545121758503e-01, -1.44542051247739378e+00], 2.00000000000000010e-04, 0),
    ([2.17590146025620967e+00, 1.78095578370790109e+00, 2.54184108992414526e+00], [3.92602150960210672e-01, 3.84730911437844370e+00, -2.36392430935964093e+00, 4.29842902921822478e-01, -1.31002424612898238e-01, -1.17380115951641528e+00], [3.36762691223188915e-01, -9.73201307751322520e-01, 1.35748539764106124e+00, 1.65005130482301499e+00, -1.16841396653392726e+00, -1.03631319494363786e-01], 2.00000000000000004e-03, -1),
    ([1.19286223133426095e+00, 7.07792494338105982e-01, 2.73986077062591882e+00], [-5.60410463617317056e-01, -2.81846960030324745e+00, 1.38689885820728609e+00, -2.38227177769569476e+00, 3.21144862950541565e+00, -2.26281394065177466e+00], [-3.13146179236673738e+00, -2.00681042942623522e+00, -1.74347286786203903e+00, -3.51423184916262399e-01, 1.18585394776181885e+00, 5.76266037609526549e-01], 2.00000000000000004e-03, -1),
    ([5.42193037743744610e-01, 8.99559235262452406e-01, 2.99108968884917648e+00], [-3.22272160854523371e-01, 1.52831933055007418e+00, -3.56265550878529602e+00, -3.72759776842369916e+00, 2.76712085156045973e+00, 7.03055525334884024e-01], [-2.82246718019786957e+00, -2.69458807480999640e+00, -3.65301342684705777e+00, -2.91102918905322561e+00, -1.27154700658335518e+00, 9.07026892940171692e-01], 2.00000000000000010e-04, 0),
    ([1.66575799300791294e+00, 8.18007290146326005e-01, 2.34811718508422995e+00], [-2.43477736043743231e+00, -3.50463811881237941e+00, 7.87136858592305089e-01, 3.16606201393025266e+00, -3.78445270892237762e+00, 2.44108791913335388e+00], [-1.48799518091911653e+00, -1.95513912619847896e+00, -4.71426570549234469e+00, -2.02467567692098394e+00, 7.73766302379057302e-01, -2.32392115187842023e-02], 2.00000000000000010e-04, 0),
    ([1.04302825880372585e+00, 1.28795686785857910e+00, 1.14535212963647082e+00], [3.82640909851373578e+00, 3.52804795721147890e+00, -1.27451059314792747e+00, -5.11987971261351760e-01, -1.48544297036917694e+00, 1.97207133844328464e+00], [-3.52808770856797116e+00, -3.31773229158019989e+00, -5.96059790052725136e-01, -1.58331978498475801e+00, 2.41131070711537099e+00, 1.68896247865388593e+00], 1.00000000000000006e-01, -1),
    ([1.86448495630739197e+00, 2.15368878578375167e+00, 2.23069776671230979e+00], [2.24843858718673673e+00, 3.42002077692546180e+00, -2.80211860969813653e+00, 1.00904126196305199e+00, -2.85102845509315905e+00, -4.54952469722424979e-01], [1.22837049515580232e+00, 1.69354133902090997e+00, 9.62917783983143760e-01, -1.72575276051070148e+00, -5.04193533213340039e-01, -1.20838662272055331e+00], 5.00000000000000010e-04, 0),
    ([8.60038394075093415e-01, 1.11078610711620751e+00, 1.39304928194881761e+00], [-3.51290637672649630e+00, 2.96307933634667009e+00, 1.09089137490499422e+00, -2.72201408033264070e+00, -1.38515382995985092e-02, -3.37064193092033904e+00], [1.03232402761968389e+00, -2.49604458360134807e+00, -3.32277160601523525e+00, -2.77110053027225822e+00, 3.17362053654905774e-01, 7.86651916673339890e-01], 1.00000000000000006e-01, -1),
    ([1.31200161783549274e+00, 2.10863107390271720e+00, 1.38016524977153177e+00], [-2.95472794567813768e+00, -1.47831774461306242e+00, -8.37854189969443830e-01, 3.30111067762475852e+00, -3.07471701024194921e+00, -3.31068887546197210e+00], [3.75011897763829372e-01, 2.82556229913386625e+00, 1.54509769412828746e+00, 7.59712062533375088e-01, -2.51128915715600565e+00, 1.77602746819837676e+00], 1.00000000000000006e-01, -1),
    ([8.59825598127321067e-01, 1.66235931506842838e+00, 6.23306638448835160e-01], [2.41500772227719640e+00, 1.74907373750626327e+00, 2.43702479934442273e+00, 2.08332512496141753e+00, -1.86193938429842110e+00, 2.31792411594808900e+00], [-2.33376958886447428e+00, -3.36909738767733069e+00, -5.27491092436127240e-01, -1.03931663360419778e-02, -2.74892438682871987e+00, 1.35769970064390266e+00], 2.00000000000000004e-03, -1),
    ([2.00625932302348620e+00, 1.09954414553579105e+00, 2.05641871958162215e+00], [-1.14197665956658856e+00, 1.87751678923000220e+00, -1.67665638704861042e+00, 2.39034743072296507e+00, -6.79116152916724225e-01, 4.25921875410980810e-01], [6.91221995875577866e-01, 7.32801599097410855e-02, -1.76350779945339853e+00, 3.48785695193217204e+00, -1.57144880764298156e+00, -6.81092122791280397e-01], 1.00000000000000008e-05, 1),
    ([5.71751180657425029e-01, 8.87282715230064190e-01, 2.49816104640149472e+00], [2.27051505596052472e+00, 1.09573117698485323e+00, 3.21439728075510267e+00, 2.04418722669644559e+00, -1.60972426440863003e+00, 1.15893220970114186e+00], [-2.23758363499190827e+00, 3.49584842325144551e+00, -1.03873355518404153e+00, -3.12645911198210102e+00, 1.20554792500021346e+00, 6.09745949174486213e-01], 5.00000000000000010e-04, 0),
    ([2.36169684788578094e+00, 1.89899478625873286e+00, 2.45735404449246486e+00], [-4.17012828970423577e-01, 5.26344325305945659e-01, -3.49948974085998810e+00, 4.40549713107427188e-01, 2.51682879934244763e+00, 1.64436418975175602e+00], [1.02689621862548575e+00, -1.32138224169707583e-02, 1.60647154535486925e+00, -1.65012632857514063e+00, 1.22310590327575786e+00, -4.19080044410062191e-01], 1.00000000000000006e-01, -1),
    ([2.04833254622550287e+00, 1.63723151269854394e+00, 1.57093150633183276e+00], [2.81623232971644377e+00, -2.89766906047984918e+00, 9.35478474503745439e-01, -6.89998812628288327e-01, 2.26070039495249731e-01, -4.72495494284963513e-03], [-1.42895127276584666e+00, 4.69274661410573335e-02, 1.76935509811815384e+00, -1.60644723580109416e+00, -2.48785440384327261e+00, -2.20218234337396712e+00], 1.00000000000000008e-05, 0),
    ([1.64948925783775224e+00, 2.93613013422899316e+00, 6.10592247468807026e-01], [3.92743610724915992e+00, 2.88424657947305541e-01, -3.03887782445258470e+00, -6.51684999858151315e-01, -2.34111943644578258e+00, 1.71321174475834503e+00], [2.01340822407103881e-01, -1.02839128813448610e+00, -6.66393796229425406e-01, 9.99892664861131264e-01, 3.48831188274834281e+00, -8.36105708052499397e-01], 2.00000000000000004e-03, -1),
    ([2.34368226591409412e+00, 2.92683657238947736e+00, 6.98981579066859426e-01], [-2.72733106854575968e+00, -1.10138640141021416e+00, 6.17183499347175868e-02, -3.42873002447842179e+00, -3.57208310114243677e+00, -2.24430198735214592e+00], [-3.89269407381562904e-01, 8.16772562077782815e-01, 3.00431709601991226e-01, -1.28683892771160813e+00, -5.20742932437557737e+00, -5.49030478492554930e-01], 5.00000000000000010e-04, 0),
    ([2.68712663334609747e+00, 2.78749279401897754e+00, 1.41276357077400583e+00], [3.09898127342397167e+00, 3.43299446691703292e+00, -7.47087889728371479e-01, 1.96858529360819468e+00, -1.07726122871137164e+00, -2.77138390204214069e+00], [2.23876178420829075e-01, -1.23128603644164092e+00, 4.68817449981431933e-01, 8.87392159453115559e-01, 2.35230268704168344e+00, -2.94316847095695999e-01], 2.00000000000000004e-03, -1),
    ([2.75622878480601852e+00, 2.67592480076779715e+00, 2.91902554506579781e+00], [7.51830794422104987e-01, 1.38722843312086930e+00, -1.01376417641774808e+00, -2.53491217954877079e+00, -1.66627989794747933e+00, 1.76525231467294130e+00], [-5.07987919242874675e-01, 5.56082951241276713e-01, 6.45521625562663764e-03, -1.52482336543070085e-01, 1.28160002775486137e+00, -9.15125895948652790e-01], 1.00000000000000008e-05, -1),
    ([1.71790198449066556e+00, 8.99286601609012193e-01, 2.84364242829696323e+00], [-9.10447040453838952e-02, -1.34427187488138777e+00, -2.65343765439635071e+00, 7.52805594468277661e-01, -8.08328595834308494e-01, -2.98255891640183091e+00], [-2.17766206475160295e+00, -5.06186737383407359e-01, 7.39861364138662325e-01, 1.67467026112784301e-01, 1.12539535470418417e+00, 1.15800354747867007e+00], 2.00000000000000010e-04, 0),
    ([2.49923162403520660e+00, 1.69403945040702753e+00, 1.80948779518857883e+00], [-8.21479520590175483e-01, -8.03649660980040537e-01, -3.66736260509332990e+00, 2.39057976764618285e+00, -2.03183754001181871e+00, -3.76528764885022049e+00], [-1.07299049261853074e-01, 6.97354172847024056e-01, -3.02086045890237687e-01, -6.41622426087977926e-01, 7.00870554202079110e-01, -1.49056653707009401e+00], 1.00000000000000006e-01, -1),
    ([5.28329696295233409e-01, 1.98105753701153109e+00, 1.82612634026149512e+00], [2.92656403844747359e+00, -6.73038389494545974e-01, 2.28892239692227495e+00, -3.85659318557231767e+00, -3.78177719939714230e+00, 7.88523578655532908e-01], [-4.14529553788396576e+00, -6.55249214561848170e+00, -1.50323925689881865e+00, -5.05818431085650944e-01, -7.13479957881156146e-01, 2.79465830910195112e-01], 5.00000000000000010e-04, 0),
    ([1.41419392335346772e+00, 9.60391653279815238e-01, 1.94507561944581897e+00], [-3.98273777456413214e+00, -9.21453069167109895e-01, -1.46210489921412812e+00, -9.00869865265808478e-01, -3.59359175227228622e+00, 4.33735088941864078e-01], [7.02408418175887173e-01, 1.41120938902962068e+00, 3.81055938584671994e+00, -1.20998450849569927e+00, 6.42813780354990372e-01, -2.38377494658656236e-01], 1.00000000000000006e-01, -1),
    ([1.03761928517418922e+00, 2.72323803820993326e+00, 2.36476595578915427e+00], [3.54231118090015507e+00, 2.57233621038073501e+00, 3.33974411667288162e+00, -2.97561669771463944e+00, -3.87428371419795958e+00, -2.41130216424539157e+00], [-6.77511738372564509e-01, 1.44671956973171678e+00, -1.80781307498824789e-01, 3.92070907927856843e-02, 8.82935960625005012e-01, -1.45264798893149583e+00], 1.00000000000000008e-05, 0),
    ([6.21821645951936675e-01, 1.12162132896884015e+00, 2.30568443399346412e+00], [1.46454607924071656e+00, -7.57708082744532163e-01, -1.64353860334828283e+00, 2.62486209603672904e+00, -1.80197288840862324e+00, -3.75030921446738752e+00], [-3.25151381945897189e+00, -6.16005082970738815e+00, 2.21601240957638312e+00, -2.69417246339282412e+00, 1.64452701031124127e+00, -9.95856360474395297e-01], 2.00000000000000004e-03, -1),
    ([1.26517853749414266e+00, 1.96298071563836229e+00, 1.27052026762706283e+00], [-9.65648917588604938e-01, -3.16631844894755421e+00, -8.87432100771267329e-01, -3.79813648089691291e+00, -2.20826462063207263e+00, 3.27576263645369359e+00], [-1.38619819707015934e+00, -1.44865976288366571e+00, 7.30572392683087712e-01, 1.56049733223908826e+00, -1.28549326478893278e+00, -7.78712090981488747e-01], 1.00000000000000006e-01, -1),
    ([2.48036493875918884e+00, 2.03448478064020399e+00, 9.47493516683794557e-01], [-2.55411743809314018e+00, 2.06362149323945321e+00, 1.65463791608571231e+00, 3.17023139486791194e+00, 2.23840066611023669e+00, -6.30266752250897788e-01], [-1.18281151924340433e+00, 1.25940199379641582e+00, 7.68914591620714960e-01, 3.00323261171252354e-01, 1.53558178545708368e+00, -9.69758693412796302e-01], 2.00000000000000004e-03, -1),
    ([1.35820032764976850e+00, 1.93859663803132976e+00, 1.70686809144746121e+00], [1.67175805329135763e+00, 3.48874516002833257e+00, 2.77835867161819916e+00, -9.00236504873722865e-01, -1.02018678861710388e+00, 3.41465815240791315e+00], [-6.20692477071483006e-01, 1.76845185151478757e+00, -1.68063054399699618e+00, -2.01762789778448082e-01, 7.61097791145553260e-02, 1.44750820984447937e+00], 5.00000000000000010e-04, 1),
    ([8.81648067529461787e-01, 2.04288954842241921e+00, 2.07284894520525143e+00], [-2.78767819289293328e+00, 2.66857566373372102e+00, 4.27709025922266228e-01, 3.69373131779822206e+00, -2.50950972927053417e+00, 1.99760820133633388e-01], [-1.60198407965998779e+00, 1.29280671081680421e+00, 4.23222635270045611e-01, 4.86637925976622782e-02, -6.75953594757319087e-01, 1.52188724175742612e-02], 1.00000000000000006e-01, -1),
    ([1.43063970821897612e+00, 2.72695024567507049e+00, 1.65380837704583961e+00], [-1.98131522833173257e-01, -1.93152790149288833e+00, -1.96931834904708492e+00, 3.63110073820083379e+00, -1.48796626318837788e+00, 2.94673316199157753e+00], [5.29982555431854796e-01, 2.68701997359779110e+00, -9.85497750596677546e-01, 1.29543667708889232e+00, 1.31802318710232358e+00, 3.28996760246946740e-01], 1.00000000000000006e-01, -1),
    ([2.30609044745737979e+00, 1.62096072556320903e+00, 2.08619539905299467e+00], [-2.17151953815658949e+00, -3.82911875886016784e+00, 1.41204588136675291e+00, -2.04343819627061052e-01, 3.24452471881129423e+00, -1.88769912007179830e+00], [1.10303094634179999e+00, 1.68710603203480436e+00, 1.96467614566644766e+00, 2.29200149003061782e+00, -1.30739057049210983e+00, -1.28432892073955007e+00], 2.00000000000000004e-03, -1),
    ([1.07076591476969996e+00, 1.40370731217299705e+00, 7.88502427293758079e-01], [-2.79447255261550076e+00, -3.85218378200208633e+00, 3.46631167144697816e+00, -1.11134881830776777e+00, 1.00825093638115693e+00, -3.63716361589508708e-01], [1.81590557848873702e+00, 2.28720660921004759e+00, -6.66864481273273135e-01, 2.01325494425147161e-01, -4.08377423925027472e+00, -2.17651678392460601e+00], 5.00000000000000010e-04, 0),
    ([1.27139791909034061e+00, 1.31113042456264584e+00, 1.93751117434272069e+00], [-1.33783994172328669e+00, -1.38511073087870873e+00, 2.65027420307789274e+00, -3.90285573895254956e+00, -3.95792253413447082e+00, -3.54306692898222853e+00], [2.30011516427580398e+00, 8.41530730610627664e-01, 1.39094146503580562e+00, 2.79263604561213041e+00, -1.16437953671955507e+00, 1.73418421729091565e+00], 2.00000000000000010e-04, 0),
    ([1.10650089069290236e+00, 1.79716431487778538e+00, 2.20267555088182476e+00], [2.37807039840408851e+00, -2.01511131755336859e+00, -1.83807119440726474e-01, -3.63527411032158110e+00, -1.69697172698220111e+00, 3.84774124775242310e+00], [2.92889658155975185e+00, 5.93514664428716365e-01, 1.27531334090080795e+00, 1.79042289477488392e+00, -6.10693422434179900e-01, 1.74515207045784249e+00], 1.00000000000000008e-05, 0),
    ([8.97241054399424254e-01, 1.82518566649791847e+00, 1.98410664523131230e+00], [1.54267663777589537e+00, 1.36762342571474260e+00, -3.66843880353449059e+00, 1.11235579817663499e+00, 3.37299975930325324e+00, 1.36084973737310122e+00], [4.23712899502120166e+00, 1.02605121494756368e-01, 1.35967978305933546e+00, -1.57672747365224386e+00, 1.59483884239052176e+00, -2.56284453814920932e-01], 5.00000000000000010e-04, 0),
    ([2.69167333397243391e+00, 1.30207117744447931e+00, 2.06080074904473154e+00], [2.77274823830603445e+00, 1.69085545084245226e+00, 9.59431519975046498e-01, 3.76517561427028813e+00, -3.29602649748663712e+00, -3.18403776696020735e+00], [1.53303318351621365e-01, 1.35139763619668862e+00, 1.14572844081091563e+00, -1.16955555829002522e+00, -1.52283497373967469e+00, -2.00752190147420878e-01], 2.00000000000000004e-03, -1),
    ([1.67068730282489564e+00, 2.04833839752214697e+00, 1.18200759423961754e+00], [3.04314499140651940e+00, -1.69336987134554828e+00, -3.14340944461934768e+00, -1.20833203054093818e+00, 3.48573342555786425e+00, -2.72461098714488870e+00], [-1.06452612323643403e+00, 2.04263271922344236e+00, -5.93333029793968403e-01, 1.33999201814845814e+00, -3.07658333933563233e+00, 2.62932083585991982e+00], 5.00000000000000010e-04, 1),
    ([1.09319553437603756e+00, 9.51416328402529654e-01, 2.35560869350484303e+00], [3.91569050229142679e+00, 3.67425922835870455e+00, 2.01850773897432045e+00, 3.31425884761755984e+00, 3.89649097047971260e+00, 4.22203275594170258e-01], [-2.60072787028988683e+00, 2.55680635946747481e+00, -8.75474788909997303e-01, -1.57530895041795271e+00, 3.09161107855181498e-03, -1.57965380014527312e+00], 2.00000000000000010e-04, 0),
    ([1.16078792581637757e+00, 1.98968876872892286e+00, 2.37353740262985990e+00], [3.65753918305446835e+00, 2.54512568343036527e+00, 3.52890256432592597e+00, 1.09071707926202510e+00, 1.87120660206140688e+00, 6.98949917653514774e-01], [-2.88471088892563898e+00, -3.01930211665610537e+00, 1.04595046701144279e+00, -1.19678497150078123e+00, 9.76835509297739990e-01, -1.21147062515036999e-01], 1.00000000000000008e-05, 0),
    ([5.17367585063153612e-01, 1.90382572732150379e+00, 1.60759485922952572e+00], [-3.40892641521308715e+00, -1.91140329991778124e+00, -1.54685146742004420e+00, -1.23455741694569365e+00, 1.65463359149587852e+00, -2.66181167062552504e+00], [-8.76865258633852007e-01, 1.65192484268691109e+00, 1.87554064649047358e+00, -1.61784721496588957e+00, -2.09346437666647267e+00, -1.53415143892705053e+00], 2.00000000000000004e-03, -1),
    ([9.06920269133183909e-01, 6.89990469926896388e-01, 1.02386369548757528e+00], [-2.21308231349285744e+00, -1.63607749428885008e-01, 1.18188421910354080e+00, -1.08905065327345074e+00, -2.65457647235188610e+00, -3.24262336752542790e+00], [-2.08991629367633225e+00, -2.90375377675716040e+00, -2.19881435835389938e-01, 4.30255801711313079e+00, 1.69124610019792576e+00, -2.31879110146838219e+00], 2.00000000000000010e-04, 0),
];

/// **Transcription check against the attached reference**, on the criterion rather than on a
/// trajectory.
///
/// Comparing trajectories would compare integrators: the reference is a fixed-step leapfrog and
/// this is Aarseth-Zare, and they diverge legitimately. What was ported is the *criterion*, so
/// that is what is checked -- the same `(m, r, v)` and the same closure value into both, and the
/// fired body compared.
///
/// **20 of the 40 rows fire, and 14 of those return a body that is NOT the one outside the
/// tightest pair.** That is what makes this a test of the ordering: the reference's
/// `b = np.argmax(fire, -1)` returns the lowest firing index, and the tightest-pair-first
/// ordering the `Distance` rule uses would fail on those 14.
#[test]
fn the_criterion_transcribes_the_reference_including_its_body_ordering() {
    let rule = outcome::EscapeRule::Closure(1e-3);
    let mut fired = 0usize;
    for (i, (m, r, v, dn, want)) in GOLDEN.iter().enumerate() {
        let mut s = Cart::<f64>::default();
        for k in 0..3 {
            s.r[k] = Vec2::new(r[2 * k], r[2 * k + 1]);
            s.v[k] = Vec2::new(v[2 * k], v[2 * k + 1]);
        }
        let got = outcome::escape_candidate_rule(&s, m, rule, 1.0, Some(*dn))
            .map_or(-1i8, |b| b as i8);
        assert_eq!(got, *want, "row {i}: reference says {want}, this port says {got}");
        if *want >= 0 {
            fired += 1;
        }
    }
    assert_eq!(fired, 20,
               "the golden set must keep firing on half its rows -- a table that fires on none \
                would pass under an implementation that never returns Some");
}
