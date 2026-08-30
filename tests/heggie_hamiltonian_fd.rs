//! Finite-difference `Gamma*` against the analytic `deriv`, for Heggie's global regularisation.
//!
//! The AZ analogue of this test caught two sign errors that were otherwise invisible — wrong AZ
//! algebra produces trajectories that look like physics. Heggie's `Gamma*` is a larger expression
//! with three coupled vectors instead of two, so the same hazard is larger here, not smaller.
//!
//! It means something only because `tests/heggie_identities.rs` anchors `gamma` independently
//! all the way down to the Cartesian energy. A sign error present in **both** `gamma` and `deriv`
//! would pass this file silently.
//!
//! Three assertions, not one:
//!   - agreement at a small step;
//!   - **the error falls as `h^2`** — an FD error insensitive to `h` is a wrong derivative, not a
//!     truncation problem;
//!   - two **mutation arms**, proving the test would catch the specific errors this algebra is
//!     prone to.

use prin_rs::integrate::az::lc;
use prin_rs::integrate::heggie::system::cyc;
use prin_rs::integrate::heggie::{hamiltonian, HgState, HgSystem};
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

const N: usize = 192;

fn masses() -> [f64; 3] {
    [3.0, 4.0, 5.0]
}

/// `dGamma*/ds[k]` for the twelve non-time components, from the analytic derivative.
///
/// `to_array13`'s layout is `u0 p0 u1 p1 u2 p2 t`, and `deriv` returns
/// `u[i] = +dGamma*/dP_i`, `p[i] = -dGamma*/dQ_i`. So the signs alternate by block, as in AZ.
fn analytic_grad(sys: &HgSystem<f64>, s: &HgState<f64>, h: f64) -> [f64; 12] {
    let d = hamiltonian::deriv(sys, s, h);
    let mut out = [0.0; 12];
    for i in 0..3 {
        out[4 * i] = -d.p[i].x; // dGamma*/dQ_i
        out[4 * i + 1] = -d.p[i].y;
        out[4 * i + 2] = d.u[i].x; // dGamma*/dP_i
        out[4 * i + 3] = d.u[i].y;
    }
    out
}

/// Per-component scaled step `h_k = h * max(|s_k|, 1)`. A fixed absolute step is wrong here: the
/// states span orders of magnitude in `|Q|` and `|P|`, so one absolute `h` is simultaneously too
/// coarse for the small components and too fine for the large ones.
fn fd_grad<F>(g: F, s: &HgState<f64>, h: f64) -> [f64; 12]
where
    F: Fn(&HgState<f64>) -> f64,
{
    let base = s.to_array13();
    let mut out = [0.0; 12];
    for k in 0..12 {
        let hk = h * base[k].abs().max(1.0);
        let (mut hi, mut lo) = (base, base);
        hi[k] += hk;
        lo[k] -= hk;
        out[k] = (g(&HgState::from_array13(hi)) - g(&HgState::from_array13(lo))) / (2.0 * hk);
    }
    out
}

/// Random states spanning orders of magnitude, with **each of the three pairs in turn made the
/// closest**.
///
/// AZ's states are implicitly conditioned by the reference-body choice — `R3` is always the
/// longest side there. Heggie has no such ordering, and sampling the same way would leave the
/// globality untested. The `tight == 3` cell is the unconditioned control.
fn random_states(n: usize, seed: u64) -> Vec<(HgState<f64>, f64)> {
    let mut rng = SplitMix64::new(seed);
    (0..n)
        .map(|i| {
            let tight = i % 4;
            let su = 10f64.powf(rng.range(-1.0, 1.0));
            let sp = 10f64.powf(rng.range(-1.0, 1.0));
            let mut s = HgState { u: [Vec2::zero(); 3], p: [Vec2::zero(); 3], t: 0.0 };
            for a in 0..3 {
                let shrink = if a == tight { 1e-2 } else { 1.0 };
                s.u[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * (su * shrink);
                s.p[a] = Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)) * sp;
            }
            (s, rng.range(-30.0, 5.0))
        })
        .collect()
}

/// Error normalised by the **magnitude of the gradient**, not by the individual component.
///
/// A gradient component passing through zero makes a per-component relative error blow up while
/// the absolute discrepancy is negligible. The object under test is the gradient as a vector, so
/// that is what sets the scale.
fn worst_rel(sys: &HgSystem<f64>, s: &HgState<f64>, h: f64, step: f64) -> f64 {
    worst_rel_with(sys, s, h, step, 0.0)
}

/// As above, with an optional `c * (Q_0.x)^3` added to `Gamma*` and to its analytic gradient.
///
/// `c = 0` is the shipped function. `c != 0` is the control for the exactness test below: the
/// shipped `Gamma*` has no third derivative for an `h^2` law to live in, and a claim of "no
/// truncation" is worth nothing unless the harness is shown able to detect truncation when it is
/// there.
fn worst_rel_with(sys: &HgSystem<f64>, s: &HgState<f64>, h: f64, step: f64, c: f64) -> f64 {
    let mut an = analytic_grad(sys, s, h);
    an[0] += 3.0 * c * s.u[0].x * s.u[0].x;
    let fd = fd_grad(|st| hamiltonian::gamma(sys, st, h) + c * st.u[0].x.powi(3), s, step);
    let scale = (0..12).map(|k| an[k].abs().max(fd[k].abs())).fold(0.0f64, f64::max).max(1e-300);
    (0..12).map(|k| (an[k] - fd[k]).abs() / scale).fold(0.0, f64::max)
}

#[test]
fn deriv_matches_finite_differenced_gamma() {
    let sys = HgSystem::new(masses());
    let mut all = Vec::new();
    let mut worst = 0.0f64;
    let mut worst_at = String::new();
    for (s, h) in random_states(N * 4, 0x4E66) {
        let w = worst_rel(&sys, &s, h, 1e-5);
        all.push(w);
        if w > worst {
            worst = w;
            worst_at =
                format!("R = ({:.2e}, {:.2e}, {:.2e})  h = {:.2e}", s.r(0), s.r(1), s.r(2), h);
        }
    }
    all.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let med = all[all.len() / 2];
    let p99 = all[(all.len() * 99) / 100];
    println!("FD vs analytic at AZ's step, h = 1e-5 * max(|s_k|,1), {} states:", all.len());
    println!("  median = {med:.3e}   p99 = {p99:.3e}   worst = {worst:.3e}");
    println!("  worst at {worst_at}");
    println!("  This step is AZ's, kept so the two are comparable. It is NOT where this");
    println!("  derivative is sharpest -- see the exactness test below, where the same");
    println!("  quantity reads 3.3e-16 at h = 1. All of the error here is roundoff.");
    assert!(med < 1e-9, "median = {med:e}");
    assert!(worst < 1e-6, "worst = {worst:e} at {worst_at}");
}

/// Median FD error across a ladder of step sizes.
fn ladder(sys: &HgSystem<f64>, states: &[(HgState<f64>, f64)], c: f64) -> Vec<(f64, f64)> {
    (0..8)
        .map(|e| {
            let step = 1.0 / 2f64.powi(e);
            let mut v: Vec<f64> =
                states.iter().map(|(s, h)| worst_rel_with(sys, s, *h, step, c)).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            (step, v[v.len() / 2])
        })
        .collect()
}

fn print_ladder(name: &str, l: &[(f64, f64)]) {
    println!("{name}");
    println!("  h        median FD error   ratio to previous");
    for (i, (h, m)) in l.iter().enumerate() {
        match i {
            0 => println!("  {h:.4}   {m:.4e}          --"),
            _ => println!("  {h:.4}   {m:.4e}        {:.4}", m / l[i - 1].1),
        }
    }
}

/// **`Gamma*` is at most QUADRATIC in every single component, so central differencing is exact
/// and there is no `h^2` regime to measure.**
///
/// The AZ analogue asserts the FD error falls as `h^2`. That test cannot be written here, and
/// finding out why is the result rather than an obstacle. `R_i = Q_ix^2 + Q_iy^2` enters every
/// term at most linearly, and `W_i = L(Q_i) P_i` is linear in each of its two arguments — so
/// although `Gamma*` is degree six *jointly*, its third derivative with respect to any one
/// coordinate is identically zero, and that is exactly what central differencing truncates at.
///
/// Measured: the median error is **3.3e-16 at `h = 1.0`**, a hundred-per-cent perturbation of
/// every component, and it **rises** as `h` falls at the `1/h` roundoff law. A floor at every
/// reachable step, not a slope. AZ has no such property: its `Gamma` carries `A B m_b m_c/|R3|`,
/// which has third derivatives in abundance.
///
/// **Asserting "no truncation" alone would be asserting that the harness produces small numbers**,
/// which a broken harness does equally well. The control arm adds a deliberately cubic term and
/// requires the same harness, over the same states, to recover the `0.25` that second-order
/// central differencing predicts. Without it this test could not fail.
#[test]
fn central_differencing_is_exact_because_gamma_is_quadratic_in_each_component() {
    let sys = HgSystem::new(masses());
    let states = random_states(N * 4, 0x1234);

    let flat = ladder(&sys, &states, 0.0);
    let cubic = ladder(&sys, &states, 1.0);
    print_ladder("Gamma* as shipped:", &flat);
    print_ladder("control -- Gamma* + (Q_0.x)^3, same harness, same states:", &cubic);

    assert!(
        flat[0].1 < 1e-14,
        "the FD error at h = 1 is {:e}. Gamma* is then not exactly differenced, so it is not \
         quadratic in each component and the reasoning this test rests on is wrong",
        flat[0].1
    );
    let rise = flat[7].1 / flat[0].1;
    println!("  shipped: error rises {rise:.2}x as h falls 128x -- the 1/h roundoff law");
    assert!(rise > 4.0, "the error did not rise as h fell (x{rise:.2}); that is not a floor");

    let mid = cubic[4].1 / cubic[3].1;
    println!("  control ratio at h = 0.0625: {mid:.4}, against the predicted 0.25");
    assert!(
        mid < 0.4,
        "the harness could not detect an h^2 law even with a cubic term present (ratio {mid}); \
         the exactness above would then be a property of the harness, not of Gamma*"
    );
}

/// The test must have teeth, and the two errors this algebra is actually prone to are not AZ's.
///
/// Heggie has no unregularised `R3` term to flip, so the analogues are:
///
///   1. **a sign flip on one coupling sub-term** in `dGamma*/dQ_i` — the piece where the argument
///      moves into the matrix slot via `L(u)w = L(w)u`, which is the step that looks like a typo;
///   2. **the cyclic `mu` pairing rotated by one**. `hamiltonian.rs` is written cyclically
///      precisely so this cannot happen in one branch and not another, and a claim that a whole
///      bug class is designed out is worth nothing unless the detector for it is shown to work.
#[test]
fn the_fd_test_detects_a_coupling_sign_error_and_a_rotated_mu() {
    let sys = HgSystem::new(masses());
    let mut worst_correct = 0.0f64;
    let mut worst_coupling = 0.0f64;
    let mut worst_rotated = 0.0f64;

    for (s, h) in random_states(N * 2, 0x9999) {
        let an = analytic_grad(&sys, &s, h);
        let fd = fd_grad(|st| hamiltonian::gamma(&sys, st, h), &s, 1e-5);
        let scale =
            (0..12).map(|k| an[k].abs().max(fd[k].abs())).fold(0.0f64, f64::max).max(1e-300);

        let w: [Vec2<f64>; 3] = std::array::from_fn(|i| HgSystem::w(&s, i));
        let r = [s.r(0), s.r(1), s.r(2)];

        for i in 0..3 {
            let (j, k) = cyc(i);

            // Arm 1: flip (R_j/m_j) L(P_i)^T W_k inside dGamma*/dQ_i. The correct term enters
            // with -1/4, so flipping it changes the gradient by twice that.
            let sub = lc::lt_apply(s.p[i], w[k]) * (r[j] * sys.inv_m[j]) / 4.0;
            let bad_q = [an[4 * i] + 2.0 * sub.x, an[4 * i + 1] + 2.0 * sub.y];

            // Arm 2: mu rotated by one in the leading term of dGamma*/dP_i.
            let delta = s.p[i] * (r[j] * r[k] / 4.0 * (1.0 / sys.mu[j] - 1.0 / sys.mu[i]));
            let bad_p = [an[4 * i + 2] + delta.x, an[4 * i + 3] + delta.y];

            for c in 0..2 {
                let (kq, kp) = (4 * i + c, 4 * i + 2 + c);
                worst_correct = worst_correct
                    .max((an[kq] - fd[kq]).abs() / scale)
                    .max((an[kp] - fd[kp]).abs() / scale);
                worst_coupling = worst_coupling.max((bad_q[c] - fd[kq]).abs() / scale);
                worst_rotated = worst_rotated.max((bad_p[c] - fd[kp]).abs() / scale);
            }
        }
    }
    println!("correct gradient:              {worst_correct:.3e}");
    println!("coupling sub-term sign flipped: {worst_coupling:.3e}");
    println!("mu pairing rotated by one:      {worst_rotated:.3e}");
    assert!(worst_correct < 1e-7, "the correct derivative did not agree: {worst_correct:e}");
    assert!(worst_coupling > 1e-3, "a flipped coupling sign still agreed at {worst_coupling:e}");
    assert!(worst_rotated > 1e-3, "a rotated mu still agreed at {worst_rotated:e}");
}

/// Heggie's Eq. (23)/(24) is a **rescaling** of Eq. (21) plus one control term, so it is checked
/// as an algebraic relation rather than finite-differenced: `Gamma*` is not the generator under
/// that time transformation and an FD test of it would be measuring the wrong object.
///
/// The `keep_gamma_term` arm must **differ** off shell and **agree** on shell — Eq. (25) is
/// legitimate precisely because `Gamma* = 0` on the solution path. A test that only checked the
/// agreement would pass for an arm that was inert everywhere.
#[test]
fn the_sum_power_time_transformation_rescales_and_adds_the_control_term() {
    use hamiltonian::HgTime;
    let sys = HgSystem::new(masses());
    let mut worst_scale = 0.0f64;
    let mut worst_onshell = 0.0f64;
    let mut biggest_offshell = 0.0f64;

    for (s, h) in random_states(N * 2, 0x2224) {
        let base = hamiltonian::deriv(&sys, &s, h);
        let f = s.s().powf(-1.5);

        // Eq. (23): the coordinate equations are a pure rescaling.
        let scaled = hamiltonian::deriv_time(&sys, &s, h, HgTime::SumPow32 { keep_gamma_term: true });
        for i in 0..3 {
            let want = base.u[i] * f;
            let sc = want.norm().max(scaled.u[i].norm()).max(1e-300);
            worst_scale = worst_scale.max((scaled.u[i] - want).norm() / sc);
        }
        // dt/dtau = R1 R2 R3 / S^{3/2}.
        let want_t = base.t * f;
        worst_scale = worst_scale.max((scaled.t - want_t).abs() / want_t.abs().max(1e-300));

        // Eq. (25) is Eq. (24) with the control term dropped: they must differ off shell.
        let dropped =
            hamiltonian::deriv_time(&sys, &s, h, HgTime::SumPow32 { keep_gamma_term: false });
        for i in 0..3 {
            let sc = scaled.p[i].norm().max(dropped.p[i].norm()).max(1e-300);
            biggest_offshell = biggest_offshell.max((scaled.p[i] - dropped.p[i]).norm() / sc);
        }

        // ...and agree on shell, where Gamma* = 0.
        let h_on = sys.energy_of(&s);
        let a = hamiltonian::deriv_time(&sys, &s, h_on, HgTime::SumPow32 { keep_gamma_term: true });
        let b = hamiltonian::deriv_time(&sys, &s, h_on, HgTime::SumPow32 { keep_gamma_term: false });
        for i in 0..3 {
            let sc = a.p[i].norm().max(b.p[i].norm()).max(1e-300);
            worst_onshell = worst_onshell.max((a.p[i] - b.p[i]).norm() / sc);
        }
    }
    println!("Eq. (23) is a pure S^(-3/2) rescaling of Eq. (21): {worst_scale:.3e}");
    println!("Eq. (24) vs Eq. (25) OFF shell (must differ):      {biggest_offshell:.3e}");
    println!("Eq. (24) vs Eq. (25) ON shell  (must agree):       {worst_onshell:.3e}");
    assert!(worst_scale < 1e-13, "Eq. (23) is not a pure rescaling: {worst_scale:e}");
    assert!(
        biggest_offshell > 1e-3,
        "the control term is inert off shell at {biggest_offshell:e} — this arm has no teeth"
    );
    assert!(worst_onshell < 1e-9, "Eqs. (24) and (25) disagree on shell: {worst_onshell:e}");
}
