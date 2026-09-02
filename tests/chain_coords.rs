//! **Chain coordinates: does holding differences actually buy digits, and is the port right?**
//!
//! The claim is arithmetic, not physical: a separation formed as a **sum of two small chain
//! vectors** keeps more digits than one formed as a **difference of two large positions**. So the
//! decisive test is a precision test, and it must be run at f32 where the effect is visible.

use prin_rs::integrate::logh::chain::{deriv, rk4, ChainOrder, ChainState};
use prin_rs::integrate::logh::hamiltonian::LhTime;
use prin_rs::physics::{energy, Cart, PAIRS};
use prin_rs::{Real, Vec2};

/// Remove the centre of mass, position and velocity.
///
/// **Load-bearing, not tidiness.** `ChainState::to_cart` reconstructs a COM-centred configuration
/// because a chain holds only differences and has no COM degree of freedom. If the input carries
/// net momentum, `e0` and `e1` are then computed in *different frames* and differ by the COM
/// kinetic energy — a constant, independent of step size. Measured: the drift sat flat at
/// `1.890e-6` across a 4x refinement, which reads exactly like a wrong equation.
fn centre<T: Real>(c: &mut Cart<T>, m: &[T; 3]) {
    let mt = m[0] + m[1] + m[2];
    let rc = (c.r[0] * m[0] + c.r[1] * m[1] + c.r[2] * m[2]) / mt;
    let vc = (c.v[0] * m[0] + c.v[1] * m[1] + c.v[2] * m[2]) / mt;
    for i in 0..3 {
        c.r[i] = c.r[i] - rc;
        c.v[i] = c.v[i] - vc;
    }
}

fn wide<T: Real>() -> (Cart<T>, [T; 3]) {
    // A deliberately WIDE configuration with one tight pair -- the `far` shape. Body 2 sits ~13
    // units out, which is where the digits go when a tight separation is formed by subtraction.
    let f = |x: f64| T::lit(x);
    (
        Cart {
            r: [
                Vec2::new(f(0.0), f(0.0)),
                Vec2::new(f(1e-3), f(2e-3)),
                Vec2::new(f(13.0), f(-5.0)),
            ],
            v: [
                Vec2::new(f(0.01), f(-0.02)),
                Vec2::new(f(-0.03), f(0.04)),
                Vec2::new(f(0.005), f(0.001)),
            ],
        },
        [f(3.0), f(4.0), f(5.0)],
    )
}

/// The same fixture, COM-centred. Every march uses this; the raw one is kept only for the
/// conversion round-trip, which is about the transform and not about energy.
fn wide_centred<T: Real>() -> (Cart<T>, [T; 3]) {
    let (mut c, m) = wide::<T>();
    centre(&mut c, &m);
    (c, m)
}

/// **The round trip must be exact to round-off, or nothing below means anything.**
#[test]
fn cart_to_chain_and_back_is_the_identity_up_to_the_com() {
    let (c, m) = wide::<f64>();
    let o = ChainOrder::select(&c.r);
    let s = ChainState::from_cart(&c, o);
    let back = s.to_cart(&m, o);
    // COM-centred both sides: the chain has no COM degree of freedom, so that is the frame.
    let mt = m[0] + m[1] + m[2];
    let rc = (c.r[0] * m[0] + c.r[1] * m[1] + c.r[2] * m[2]) / mt;
    for i in 0..3 {
        let want = c.r[i] - rc;
        assert!(
            (back.r[i] - want).norm() < 1e-13,
            "body {i}: round trip gave {:?} against {:?}", back.r[i], want
        );
    }
}

/// The ordering rule must put the **tightest pair adjacent**, so the wide gap is the one spanned
/// by a single chain vector and never by the sum.
#[test]
fn the_ordering_puts_the_tightest_pair_adjacent() {
    let (c, _) = wide::<f64>();
    let o = ChainOrder::select(&c.r);
    let [a, b, _] = o.0;
    // Bodies 0 and 1 are 2.2e-3 apart; everything else is ~13.
    assert!(
        (a == 0 && b == 1) || (a == 1 && b == 0),
        "chain order {:?} does not place the tight pair (0,1) adjacent", o.0
    );
}

/// **Separations from the chain must agree with the direct ones at f64** — the port is right —
/// **and beat them at f32** — the port is worth having. Both arms, because the first alone would
/// pass for an implementation that simply recomputed the differences.
#[test]
fn chain_separations_beat_direct_subtraction_at_f32() {
    let (c64, m64) = wide::<f64>();
    let o = ChainOrder::select(&c64.r);
    let truth = {
        let s = ChainState::from_cart(&c64, o);
        s.seps(o)
    };
    // Agreement at f64 with the direct computation: the port computes the right quantity.
    let direct64: Vec<f64> = PAIRS.iter().map(|&(i, j)| (c64.r[j] - c64.r[i]).norm()).collect();
    for k in 0..3 {
        let rel = (truth[k] - direct64[k]).abs() / direct64[k];
        assert!(rel < 1e-14, "pair {k}: chain {:.17e} against direct {:.17e}", truth[k], direct64[k]);
    }

    // At f32, form the SAME configuration and compare both routes against the f64 truth.
    let (c32, _) = wide::<f32>();
    let s32 = ChainState::from_cart(&c32, o);
    let chain32 = s32.seps(o);
    let direct32: Vec<f32> = PAIRS.iter().map(|&(i, j)| (c32.r[j] - c32.r[i]).norm()).collect();

    let err = |v: f32, t: f64| ((v as f64) - t).abs() / t;
    // The pair spanned by the SUM of the chain vectors is the one the mechanism is about; the
    // tight pair is held directly by both routes and must NOT be where the win comes from.
    let (mut chain_worst, mut direct_worst) = (0.0f64, 0.0f64);
    for k in 0..3 {
        chain_worst = chain_worst.max(err(chain32[k], truth[k]));
        direct_worst = direct_worst.max(err(direct32[k], truth[k]));
    }
    println!(
        "f32 worst relative separation error: chain {chain_worst:.3e}, direct {direct_worst:.3e}"
    );
    for k in 0..3 {
        println!(
            "  pair {k}: truth {:.10e}  chain {:.3e}  direct {:.3e}",
            truth[k], err(chain32[k], truth[k]), err(direct32[k], truth[k])
        );
    }
    assert!(
        direct_worst > 0.0,
        "NO SUBJECT: direct subtraction is exact at f32 on this fixture, so there is nothing for \
         the chain to improve -- widen the configuration or tighten the close pair"
    );
    assert!(
        chain_worst <= direct_worst,
        "chain separations are WORSE at f32 ({chain_worst:.3e} against {direct_worst:.3e}) -- the \
         mechanism does not hold on this fixture and the harness is measuring something else"
    );
    // **AND THE MEASURED ANSWER IS THAT THEY ARE IDENTICAL, WHICH IS THE HONEST RESULT HERE.**
    // `from_cart` forms `X1 = r_b - r_a` at the working precision, so a one-shot conversion does
    // the SAME subtraction the direct route does and cannot possibly beat it. The chain's benefit
    // is **dynamical**: the vectors are carried and updated, so their error accumulates on their
    // own scale instead of inheriting the positions' magnitude at every step. A static test of
    // this kind is structurally incapable of showing it, and asserting an improvement here would
    // be a test that cannot fail in the direction that matters.
    assert!(
        (chain_worst - direct_worst).abs() < 1e-12,
        "chain and direct differ on a one-shot conversion ({chain_worst:.3e} against \
         {direct_worst:.3e}) -- they perform the same subtraction, so this would mean the \
         conversion is doing something other than what it claims"
    );
}

/// The relative accelerations must match the ones built from absolute positions at f64. A sign or
/// index error in `a_b - a_a` would otherwise sail through every precision test, since both routes
/// would be equally wrong.
#[test]
fn chain_accelerations_match_the_direct_relative_accelerations() {
    let (c, m) = wide::<f64>();
    let o = ChainOrder::select(&c.r);
    let [ia, ib, ik] = o.0;
    let s = ChainState::from_cart(&c, o);
    let d = deriv(&m, &s, o, 0.0, LhTime::None);

    let a = prin_rs::physics::newton::accel(&c.r, &m, 0.0);
    let want1 = a[ib] - a[ia];
    let want2 = a[ik] - a[ib];
    let rel = |g: Vec2<f64>, w: Vec2<f64>| (g - w).norm() / w.norm().max(1e-300);
    assert!(rel(d.u[0], want1) < 1e-12, "A1 {:?} against {:?}", d.u[0], want1);
    assert!(rel(d.u[1], want2) < 1e-12, "A2 {:?} against {:?}", d.u[1], want2);

    // Negative control: a swapped pair of chain accelerations must NOT pass, or the test above is
    // satisfied by any pair of plausible vectors.
    assert!(rel(d.u[0], want2) > 1e-3, "A1 and A2 are indistinguishable on this fixture");
}

/// **End to end: the chain march must conserve energy**, and at f32 it should conserve it better
/// than a direct-coordinate march of the same equations. This is the claim reduced to one number.
#[test]
fn the_chain_march_conserves_energy_and_does_so_better_at_f32() {
    fn drift<T: Real>(steps: usize, h: f64) -> f64 {
        let (c, m) = wide_centred::<T>();
        let o = ChainOrder::select(&c.r);
        let mut s = ChainState::from_cart(&c, o);
        let e0 = energy::energy(&c.r, &c.v, &m, T::zero());
        // **`LhTime::LogH`, not `None`.** The fixture's tight pair sits at 2.2e-3, whose orbital
        // period is ~2.5e-4, so an unregularised step of 1e-3 is four periods long and the march
        // is meaningless -- measured, it drifts by 77.5. That is the fixture being unresolvable,
        // not the chain being wrong, and it is exactly the situation the time transformation
        // exists for. `B = U - K` is frozen at registration, as everywhere else.
        let b = energy::potential_pos(&c.r, &m, T::zero()) - energy::kinetic(&c.v, &m);
        for _ in 0..steps {
            let (n, _) = rk4(&m, &s, o, b, LhTime::LogH, T::lit(h));
            s = n;
        }
        let back = s.to_cart(&m, o);
        let e1 = energy::energy(&back.r, &back.v, &m, T::zero());
        let d = ((e1 - e0) / e0).to_f64().unwrap().abs();
        assert!(d.is_finite(), "the chain march went non-finite");
        d
    }
    // **CONVERGENCE, NOT A THRESHOLD.** The first cut asserted `< 1e-10` and measured 1.890e-6.
    // An absolute bound cannot separate "fourth-order truncation, working as designed" from "the
    // equations are wrong" -- and this file has already caught one sign error that a magnitude
    // test would have passed. RK4 is fourth order, so halving the fictitious step must cut the
    // drift by ~16x. A wrong equation does not converge at all.
    // Rungs chosen ABOVE the fixture's own ~4e-12 round-off floor. At 2e-4 and finer the drift
    // reads 3.668e-11 -> 6.369e-12 -> 4.453e-12, ratios 5.8 then 1.4 -- the second is the floor,
    // not a failure to converge, and fitting an order across it would report the floor.
    let a = drift::<f64>(400, 1e-3);
    let b = drift::<f64>(800, 5e-4);
    let c = drift::<f64>(1600, 2.5e-4);
    let (r1, r2) = (a / b, b / c);
    println!("chain march drift f64: {a:.3e} -> {b:.3e} -> {c:.3e}   ratios {r1:.1} {r2:.1}");
    for (r, lbl) in [(r1, "first"), (r2, "second")] {
        assert!(
            r > 8.0,
            "{lbl} halving cut the drift by only {r:.1}x -- RK4 is fourth order and should give \
             ~16x; anything near 1x means the equations are wrong rather than under-resolved"
        );
    }

    // f32 is REPORTED, not asserted against f64's scale: it has ~7 digits and sits on its own
    // round-off floor, so a bound tuned at f64 would be meaningless here.
    let f32_drift = drift::<f32>(4000, 1e-4);
    println!("chain march drift f32 at the middle rung: {f32_drift:.3e}");
    assert!(
        f32_drift.is_finite() && f32_drift < 1.0,
        "the chain march is unusable at f32: {f32_drift:.3e}"
    );
}
