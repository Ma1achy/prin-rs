//! Finite-difference `Lambda` against the analytic equations of motion, for logH.
//!
//! # Why the off-shell draw is the whole design of this file
//!
//! `Lambda = ln(K + B) - ln(U)`, and on the solution path `K + B == U`. So **the most plausible
//! transcription error — swapping which denominator the drift and the kick use — changes nothing
//! on shell.** A test written on physical initial conditions would pass with the two halves
//! exchanged and report a clean derivative.
//!
//! Every state here is therefore drawn with `K + B = lam * U` for a random `lam != 1`, so `lam`
//! *is* the off-shellness, and `lam = 1` is available as the negative control that demonstrates
//! the swap really is invisible there. This is the same move as the random off-shell `h` in
//! `tests/heggie_hamiltonian_fd.rs`, for the same reason and against a sharper hazard.
//!
//! # And the `h^2` law has a subject here, unlike `Gamma*`
//!
//! Heggie's `Gamma*` is at most **quadratic in each single component**, so its third derivative
//! vanishes identically, central differencing is exact, and the measured error *rises* as `h`
//! falls at the `1/h` roundoff law. `Lambda` is a logarithm of a quadratic; third derivatives
//! are abundant. So the ratio test applies here and the ladder is read before it is asserted on.
//!
//! Assertions, in order: the sign identity `potential_pos == -potential`; agreement at a small
//! step, on both arms; the `h^2` law; and four mutation arms, two of which are *negative*
//! controls asserting a mutation is invisible where it must be.

use prin_rs::integrate::logh::hamiltonian::{self, Dens, LhTime};
use prin_rs::integrate::logh::{LhState, LhSystem};
use prin_rs::physics::energy;
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

const N: usize = 256;

fn masses() -> [f64; 3] {
    [3.0, 4.0, 5.0]
}

/// `dLambda/ds[k]` for the twelve non-time components, from the analytic equations of motion.
///
/// `to_array14`'s layout is `r0 v0 r1 v1 r2 v2 t`, and `deriv_with` returns
/// `r[i] = dLambda/dp_i` and `v[i] = -(1/m_i) dLambda/dr_i`. Momenta are `p_i = m_i v_i`, so the
/// mass factors go back in here rather than in the physics — the march never needs them.
fn analytic_grad(sys: &LhSystem<f64>, s: &LhState<f64>, d: Dens<f64>) -> [f64; 12] {
    let k = hamiltonian::deriv_with(sys, s, d);
    let m = sys.masses;
    let mut out = [0.0; 12];
    for i in 0..3 {
        out[4 * i] = -m[i] * k.v[i].x; // dLambda/dr_i
        out[4 * i + 1] = -m[i] * k.v[i].y;
        out[4 * i + 2] = m[i] * k.r[i].x; // dLambda/dv_i
        out[4 * i + 3] = m[i] * k.r[i].y;
    }
    out
}

/// Per-component step, scaled by **the scale on which each denominator actually varies**.
///
/// Not `h * max(|s_k|, 1)`, which is the Heggie file's rule and is right for `Gamma*` — a
/// polynomial with no small scale to violate. `Lambda` is `ln(K + B) - ln(U)`, and its curvature
/// in a coordinate is set by how far that coordinate has to move to change the corresponding
/// **denominator** by order one:
///
/// ```text
///   position slot i:   U / (m_i |a_i|)          the distance over which U changes by ~100%
///   velocity slot i:   (K + B) / (m_i |v_i|)    the velocity change that moves K + B by ~100%
/// ```
///
/// capped by the state's own `d_min` and `|v|_max` so a near-zero gradient cannot inflate the
/// step. **The velocity rule is the load-bearing one and it is specific to an off-shell draw.**
/// Off shell `K + B = lam U` can be arbitrarily small next to `K` itself: the measured worst
/// case had `K = 2.0e3` against `K + B = 1.6`, so `ln(K+B)`'s third derivative in `v` runs as
/// `(m v / (K+B))^3 ~ 4.6e5` and a step scaled by `|v|` truncates at `~5e-5` — which is exactly
/// the tail that appeared, `3.4e-5` against a median of `8.6e-11`.
///
/// Two guards were tried and did not find it: **the median moved 5x and hid it**, and **every
/// ladder ratio still read a clean 0.2500** — because the ladder is a median too. The tail was
/// only ever visible as `worst`, and only diagnosable by printing the state behind it. This is
/// *what order is the function in this variable* asked about magnitude instead of degree, and
/// the answer is again that a smaller absolute step is not the conservative choice.
///
/// It also means the two requirements interact: the off-shell draw the denominator-swap mutation
/// needs is precisely what makes `K + B` small and the naive step rule fail.
fn fd_grad<F>(g: F, s: &LhState<f64>, b: f64, h: f64) -> [f64; 12]
where
    F: Fn(&LhState<f64>) -> f64,
{
    let m = masses();
    let base = s.to_array14();
    let d = prin_rs::physics::newton::pair_dists(&s.r);
    let d_min = d.iter().fold(f64::INFINITY, |a, &x| a.min(x));
    let v_max = s.v.iter().fold(0.0f64, |a, x| a.max(x.norm()));
    let u = energy::potential_pos(&s.r, &m, 0.0);
    let kb = energy::kinetic(&s.v, &m) + b;
    let a = prin_rs::physics::newton::accel(&s.r, &m, 0.0);

    let mut out = [0.0; 12];
    for k in 0..12 {
        let i = k / 4;
        // Layout is r0 v0 r1 v1 r2 v2: slots 0,1 of each group of four are positions.
        let scale = if k % 4 < 2 {
            (u / (m[i] * a[i].norm()).max(1e-300)).min(d_min)
        } else {
            (kb.abs() / (m[i] * s.v[i].norm()).max(1e-300)).min(v_max)
        }
        .max(1e-300);
        let hk = h * scale;
        let (mut hi, mut lo) = (base, base);
        hi[k] += hk;
        lo[k] -= hk;
        out[k] = (g(&LhState::from_array14(hi)) - g(&LhState::from_array14(lo))) / (2.0 * hk);
    }
    out
}

/// Random Cartesian states, with **each of the three pairs in turn made the closest**, and `B`
/// set so that `K + B = lam * U`.
///
/// `lam` is the off-shellness: `lam = 1` is the physical shell, where the two denominators
/// coincide and a swap is undetectable. Returns `(state, B, lam)`.
fn random_states(n: usize, seed: u64, on_shell: bool) -> Vec<(LhState<f64>, f64, f64)> {
    let m = masses();
    let mut rng = SplitMix64::new(seed);
    let mut out = Vec::with_capacity(n);
    let mut i = 0usize;
    while out.len() < n {
        let tight = i % 4;
        i += 1;
        let sr = 10f64.powf(rng.range(-1.0, 1.0));
        let sv = 10f64.powf(rng.range(-1.0, 1.0));
        let mut s = LhState { r: [Vec2::zero(); 3], v: [Vec2::zero(); 3], t: 0.0, w: 0.0 };
        for a in 0..3 {
            // Shrink one body toward the origin so the pairs it forms are the close ones.
            let shrink = if a == tight { 1e-2 } else { 1.0 };
            s.r[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * (sr * shrink);
            s.v[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * sv;
        }
        let u = energy::potential_pos(&s.r, &m, 0.0);
        let k = energy::kinetic(&s.v, &m);
        if !(u.is_finite() && u > 1e-12) {
            continue;
        }
        // `lam` bounded away from 1 on the off-shell draw, so "off shell" is not a rounding
        // accident on some sample; and bounded away from 0 so `ln(K+B)` stays defined under the
        // finite-difference perturbations too.
        let lam = if on_shell {
            1.0
        } else if rng.range(0.0, 1.0) < 0.5 {
            rng.range(0.3, 0.7)
        } else {
            rng.range(1.4, 3.0)
        };
        out.push((s, lam * u - k, lam));
    }
    out
}

/// Error normalised by the **magnitude of the gradient**, not by the individual component: a
/// component passing through zero makes a per-component relative error blow up while the
/// absolute discrepancy is negligible.
fn worst_rel(
    sys: &LhSystem<f64>,
    s: &LhState<f64>,
    b: f64,
    time: LhTime,
    step: f64,
    dens: Option<Dens<f64>>,
) -> f64 {
    let d = dens.unwrap_or_else(|| hamiltonian::denominators(sys, s, b, time));
    let an = analytic_grad(sys, s, d);
    let fd = fd_grad(|st| hamiltonian::lambda(sys, st, b, time), s, b, step);
    let scale = (0..12).map(|k| an[k].abs().max(fd[k].abs())).fold(0.0f64, f64::max).max(1e-300);
    (0..12).map(|k| (an[k] - fd[k]).abs() / scale).fold(0.0, f64::max)
}

fn stats(mut v: Vec<f64>) -> (f64, f64, f64) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (v[v.len() / 2], v[(v.len() * 99) / 100], *v.last().unwrap())
}

/// The identity that gives the sign convention a name.
///
/// `dU/dr_i = m_i a_i` is what makes the kick `a_i / U`, and it is true only for the **positive**
/// potential. Bitwise, because it is a negation and nothing else — an approximate assertion here
/// would be hiding a difference that cannot exist.
#[test]
fn potential_pos_is_exactly_the_negated_potential() {
    let m = masses();
    for (s, _, _) in random_states(64, 0x5164, false) {
        let a = energy::potential_pos(&s.r, &m, 0.0);
        let b = -energy::potential(&s.r, &m, 0.0);
        assert_eq!(a.to_bits(), b.to_bits(), "potential_pos != -potential: {a:e} vs {b:e}");
        assert!(a > 0.0, "the positive potential came out non-positive: {a:e}");
    }
}

#[test]
fn deriv_matches_finite_differenced_lambda() {
    let sys = LhSystem::new(masses());
    let mut seen: Vec<(f64, f64)> = Vec::new();
    for (time, label) in [(LhTime::LogH, "LogH"), (LhTime::None, "None (control)")] {
        let st = random_states(N, 0x1067, false);
        let v: Vec<f64> =
            st.iter().map(|(s, b, _)| worst_rel(&sys, s, *b, time, 1e-5, None)).collect();
        let (bi, bw) = v.iter().enumerate().fold((0, 0.0f64), |a, (i, &x)| if x > a.1 { (i, x) } else { a });
        {
            let (s, b, lam) = &st[bi];
            let d = prin_rs::physics::newton::pair_dists(&s.r);
            let u = energy::potential_pos(&s.r, &masses(), 0.0);
            let k = energy::kinetic(&s.v, &masses());
            println!(
                "  worst state [{label}] = {bw:.3e}: d = ({:.2e},{:.2e},{:.2e})  lam = {lam:.3}  \
                 U = {u:.3e}  K = {k:.3e}  K+B = {:.3e}  |r|max = {:.2e}  |v|max = {:.2e}",
                d[0], d[1], d[2], k + b,
                s.r.iter().fold(0.0f64, |a, x| a.max(x.norm())),
                s.v.iter().fold(0.0f64, |a, x| a.max(x.norm())),
            );
        }
        let (med, p99, worst) = stats(v);
        println!(
            "FD vs analytic, {label:<15} h = 1e-5 * (denominator scale), {N} OFF-SHELL states:\n  \
             median = {med:.3e}   p99 = {p99:.3e}   worst = {worst:.3e}"
        );
        assert!(med < 1e-9, "{label}: median = {med:e}");
        assert!(worst < 1e-7, "{label}: worst = {worst:e}");
        seen.push((med, worst));
    }
    // **The two arms must not agree bitwise.** They read the same states through the same
    // harness, so a `LhTime` that never reached `lambda` would make this file a test of the
    // `None` arm printed twice under two headings — a difference can be small because both sides
    // are right or because one side is dead, and four matching digits is not the same as
    // identity.
    assert!(
        seen[0] != seen[1],
        "the LogH and None arms are bitwise identical, so LhTime is inert in `lambda`: {seen:?}"
    );
    println!(
        "\n  The `None` arm is the control on the HARNESS: it differentiates `K - U`, whose\n  \
         gradient is elementary. If it were to fail, the mass factors in `analytic_grad` would\n  \
         be wrong and the LogH number would mean nothing."
    );
}

/// Median FD error across a ladder of step sizes, and the ratio between neighbouring rungs.
#[test]
fn finite_difference_error_falls_as_h_squared() {
    let sys = LhSystem::new(masses());
    let states = random_states(N, 0x21EF, false);
    println!("median FD error against step, LogH, off shell:");
    let mut rows = Vec::new();
    for e in 0..14 {
        let step = 1e-2 / 4f64.powi(e as i32).sqrt(); // 1e-2, 5e-3, 2.5e-3, ...
        let v: Vec<f64> =
            states.iter().map(|(s, b, _)| worst_rel(&sys, s, *b, LhTime::LogH, step, None)).collect();
        let (med, _, _) = stats(v);
        rows.push((step, med));
    }
    for (i, (step, med)) in rows.iter().enumerate() {
        let r = if i == 0 { f64::NAN } else { med / rows[i - 1].1 };
        println!("  h = {step:.3e}   median = {med:.4e}   ratio = {r:.4}");
    }
    // Read the halvings that are still truncation-dominated. Unlike `Gamma*`, which is exactly
    // differenced and whose error RISES as h falls, this one has third derivatives to truncate.
    let r1 = rows[1].1 / rows[0].1;
    let r2 = rows[2].1 / rows[1].1;
    println!("\n  second-order prediction is 0.25 per halving; measured {r1:.4} and {r2:.4}");
    assert!(r1 < 0.4, "the first halving did not fall as h^2: ratio {r1:.4}");
    assert!(r2 < 0.4, "the second halving did not fall as h^2: ratio {r2:.4}");
}

/// The four mutations, two positive arms and **two negative controls**.
///
/// The negative controls are the load-bearing half. *A test that cannot fail is
/// indistinguishable from a test that passes*, and its converse holds too: an assertion that a
/// mutation fires means nothing unless the conditions under which it must **not** fire are also
/// pinned. Here the swap is genuinely undetectable on shell, and saying so is the finding.
#[test]
fn the_fd_test_detects_a_denominator_swap_only_off_shell() {
    let sys = LhSystem::new(masses());

    let off = random_states(N, 0x9C0F, false);
    let on = random_states(N, 0x9C0F, true);

    let swapped = |sys: &LhSystem<f64>, s: &LhState<f64>, b: f64| {
        let d = hamiltonian::denominators(sys, s, b, LhTime::LogH);
        Dens { drift: d.kick, kick: d.drift }
    };

    // --- positive arm: off shell, the swap must be caught ---
    let v: Vec<f64> = off
        .iter()
        .map(|(s, b, _)| worst_rel(&sys, s, *b, LhTime::LogH, 1e-4, Some(swapped(&sys, s, *b))))
        .collect();
    let (med_off, _, _) = stats(v);

    // --- negative control: on shell, the swap is invisible and must be ---
    let v: Vec<f64> = on
        .iter()
        .map(|(s, b, _)| worst_rel(&sys, s, *b, LhTime::LogH, 1e-4, Some(swapped(&sys, s, *b))))
        .collect();
    let (med_on, _, _) = stats(v);

    // --- and the unmutated on-shell arm, so `med_on` being small is not just a dead harness ---
    let v: Vec<f64> =
        on.iter().map(|(s, b, _)| worst_rel(&sys, s, *b, LhTime::LogH, 1e-4, None)).collect();
    let (clean_on, _, _) = stats(v);

    println!("denominator swap, median relative FD disagreement:");
    println!("  OFF shell (lam != 1) : {med_off:.3e}    <- must be caught");
    println!("  ON  shell (lam == 1) : {med_on:.3e}    <- must NOT be, and is not");
    println!("  ON  shell, unmutated : {clean_on:.3e}");
    println!(
        "\n  **This is why every state in this file is drawn off shell.** A test written on\n  \
         physical initial conditions passes with the two halves exchanged."
    );
    assert!(med_off > 1e-3, "the swap was NOT caught off shell: {med_off:e}");
    // Not "both are small" -- **the same number**. The swapped and unswapped derivatives are the
    // same function on shell, so what is left in either is the harness's own truncation, and the
    // two must agree to it. A pair of loose bounds would pass on two different small numbers.
    assert!(
        (med_on - clean_on).abs() <= 1e-12 * clean_on.max(1e-300),
        "the swap CHANGED something on shell (mutated {med_on:e} against clean {clean_on:e}), so \
         `K + B == U` does not hold there and the whole off-shell argument is wrong"
    );
}

#[test]
fn the_fd_test_detects_a_kick_sign_flip_and_a_wrong_b_sign() {
    let sys = LhSystem::new(masses());
    let states = random_states(N, 0x3B71, false);

    // Kick sign flip: `dv/ds = -a/U`. Attraction becomes repulsion.
    let v: Vec<f64> = states
        .iter()
        .map(|(s, b, _)| {
            let d = hamiltonian::denominators(&sys, s, *b, LhTime::LogH);
            worst_rel(&sys, s, *b, LhTime::LogH, 1e-4, Some(Dens { drift: d.drift, kick: -d.kick }))
        })
        .collect();
    let (med_kick, _, _) = stats(v);

    // `B` given the wrong sign: `B = K - U` rather than `U - K`. This is `+E`, not `-E`, and it
    // is the error that would make `K + B` equal `2K - U` instead of `U`.
    let v: Vec<f64> = states
        .iter()
        .map(|(s, b, _)| {
            let bad = -*b;
            let d = hamiltonian::denominators(&sys, s, bad, LhTime::LogH);
            // FD the TRUE Lambda while the analytic side uses the wrong B, so the disagreement
            // is the sign error and not a consistent pair of wrong things.
            let an = analytic_grad(&sys, s, d);
            let fd = fd_grad(|st| hamiltonian::lambda(&sys, st, *b, LhTime::LogH), s, *b, 1e-4);
            let scale =
                (0..12).map(|k| an[k].abs().max(fd[k].abs())).fold(0.0f64, f64::max).max(1e-300);
            (0..12).map(|k| (an[k] - fd[k]).abs() / scale).fold(0.0, f64::max)
        })
        .collect();
    let (med_b, _, _) = stats(v);

    println!("median relative FD disagreement under mutation:");
    println!("  kick sign flipped    : {med_kick:.3e}");
    println!("  B sign flipped       : {med_b:.3e}");
    assert!(med_kick > 1e-3, "a flipped kick sign was not caught: {med_kick:e}");
    assert!(med_b > 1e-3, "a flipped B sign was not caught: {med_b:e}");
}
