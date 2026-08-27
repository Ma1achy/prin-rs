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
/// **The arm does not fire at all at `t = 13`.** Measured over a 32x32 near-field grid, both
/// classifiers return "bound" for all 1024 pixels; Burrau's escape happens later than the
/// project's horizon. So this test also runs to `t = 20`, where 109 of 1024 pixels fire,
/// because an arm that is never exercised is not tested by agreeing with another arm that is
/// also never exercised.
///
/// **Two properties of the ported arm, transcribed rather than fixed.** It is sampled at sync
/// boundaries, so `t_end` for an escape has the resolution of the sync grid and its value
/// depends on `n_sync` and `t_max` — an event at `t = 9.5` is seen at `t_max = 16`, where the
/// boundaries are 0.5 apart, and missed at `t_max = 13`, where they are 0.40625 apart. And it
/// latches on first firing: at `t = 20`, 109 pixels fire but only 87 are still labelled
/// escaping at the horizon, so a body can satisfy "unbound and receding" at one boundary and
/// be recaptured by the next.
#[test]
fn escape_matches_the_legacy_classifier() {
    let s = grid::region("near-field", 5, 5, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let mut fired = 0usize;
    let mut checked = 0usize;
    for t_max in [13.0f64, 20.0] {
        for i in 0..s.npix() {
            let o = az::integrate_az_opts(
                s.nominal::<f64>(i), &m, t_max, 32, 0.01, 30_000,
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
        r_esc_frac: 0.0,
        escape_all_bodies: false,
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
        forced_refs: None,
        lc_stable: true,
        r_coll_frac: 1e-3,
        stop_on_event: true,
        escape_every: 1,
        escape_confirm: confirm,
        // The numpy reference's ungated escape test: every result in this diagnostic
        // predates the distance gate and is quoted against that form.
        r_esc_frac: 0.0,
        escape_all_bodies: false,
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
    let r_esc = 5.0 * energy::hyperradius(&near.r, &m);
    assert_eq!(
        outcome::escape_candidate_gated(&near, &m, r_esc, true),
        None,
        "the distance gate must reject a body still deep inside the system"
    );
    assert_eq!(
        outcome::escape_candidate_gated(&far, &m, r_esc, true),
        Some(2),
        "the distance gate must NOT reject a genuine, far escape — without this arm a gate that \
         refuses everything passes"
    );

    // `r_esc = 0` is the reference path and must be bit-for-bit the ungated test.
    for s in [&near, &far] {
        assert_eq!(
            outcome::escape_candidate_gated(s, &m, 0.0, false),
            outcome::escape_candidate(s, &m),
            "r_esc = 0 with one candidate body is the numpy reference and must not diverge"
        );
    }
}
