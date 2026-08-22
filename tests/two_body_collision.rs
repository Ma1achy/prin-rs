//! Gate (b): the two-body radial collision.
//!
//! Equal masses released from rest with the third body far away. The pair falls radially
//! into an *exact* collision. Regularisation is what makes this survivable: in the
//! transformed Hamiltonian the singular `-G m_i m_j / |rho|` term becomes a constant, the
//! regularised two-body problem is a harmonic oscillator, and the trajectory passes through
//! collision at machine precision.
//!
//! Acceptance: `d_min < 1e-10` with `|dE/E| < 1e-12`.

use prin_rs::integrate::az;
use prin_rs::integrate::az::reference_body::choose_reference;
use prin_rs::physics::Cart;
use prin_rs::Vec2;

fn setup() -> (Cart<f64>, [f64; 3]) {
    // Bodies 0 and 1: equal masses, unit separation, at rest.
    // Body 2: far away, and deliberately off-axis so the two long sides are not an exact
    // tie — an argmax tie would make the reference-body choice depend on tie-breaking
    // rather than on geometry.
    let m = [1.0, 1.0, 1.0];
    let s = Cart::new(
        [
            Vec2::new(-0.5, 0.0),
            Vec2::new(0.5, 0.0),
            Vec2::new(0.1, 1000.0),
        ],
        [Vec2::zero(); 3],
    );
    (s, m)
}

#[test]
fn the_close_pair_is_the_regularised_one() {
    // Worth asserting rather than assuming. AZ regularises the two pairs sharing the
    // reference body; if the geometry put the colliding pair on the *unregularised* side,
    // this test would be measuring nothing.
    let (s, _) = setup();
    let a = choose_reference(&s.r);
    let (ra, rb, rc) = az::reference_body::triple(a);
    println!("reference body = {a}, regularised pairs = ({ra},{rb}) and ({ra},{rc})");
    let pair_is_regularised = (ra == 0 || rb == 0 || rc == 0) && (ra == 1 || rb == 1 || rc == 1);
    assert!(pair_is_regularised, "the colliding pair (0,1) is not regularised");
    // Specifically, (a,b) must be the close pair.
    assert!(
        (ra == 1 && rb == 0) || (ra == 0 && rb == 1),
        "the close pair is not the first regularised pair: a={ra} b={rb} c={rc}"
    );
}

#[test]
fn radial_collision_passes_through_at_machine_precision() {
    let (s0, m) = setup();

    // Free-fall time from rest at separation d for total mass M is
    // (pi/2) sqrt(d^3 / (2 G M)); with d = 1, M = 2 that is ~0.785. Integrate past it so
    // the trajectory goes through collision rather than stopping at it.
    let t_max = 1.0;
    // eta is pinned small deliberately — see the scan test below. d_min at a genuine
    // collision is sampling-limited, not physical, so this threshold is a statement about
    // resolution as much as about correctness.
    let eta = 1e-4;
    let n_sync = 1;

    let out = az::integrate_az(s0, &m, t_max, n_sync, eta, 20_000_000, None);

    println!("two-body radial collision, t_max = {t_max}, n_sync = {n_sync}, eta = {eta}");
    println!("  NOTE: eta = 1e-4 is 100x finer than the production eta = 1e-2. This gate is");
    println!("  not a production guarantee — see d_min_threshold_is_unreachable_at_production_eta.");
    println!("  d_min (regularised pairs) = {:.6e}", out.d_min_ref);
    println!("  d_min (all three pairs)   = {:.6e}", out.d_min_true);
    println!("  |dE/E|                    = {:.6e}", out.drift);
    println!("  max |Gamma| / |largest term| = {:.6e}", out.gamma_max);
    println!("  steps = {}, reference switches = {}", out.steps, out.switches);

    assert!(out.finite, "trajectory went non-finite");
    assert!(!out.budget_exhausted, "step budget exhausted");
    assert!(out.d_min_ref < 1e-10, "d_min = {:e}, want < 1e-10", out.d_min_ref);
    assert!(out.drift < 1e-12, "|dE/E| = {:e}, want < 1e-12", out.drift);
}

/// **d_min at a near-collision is a phase-limited sampling artefact.**
///
/// `u` crosses zero roughly linearly in fictitious time while `|R1| = |u1|^2`, so the closest
/// *sampled* separation is `(|u'| * dphase)^2` where `dphase` is the distance from the
/// crossing to the nearest sample — essentially a uniform draw in `[0, dtau/2]`.
///
/// The consequence matters more than the mechanism: the **scale** falls as `eta^2`, but any
/// single realisation is dominated by where the sample happens to land and scatters over
/// four or five orders. A scaling law read off one trajectory is therefore meaningless — the
/// realisation scatter is as wide as the shift between decades. It needs an ensemble, which
/// is what this test uses.
#[test]
fn d_min_at_collision_is_phase_limited() {
    fn perturbed(d: f64) -> (Cart<f64>, [f64; 3]) {
        (
            Cart::new(
                [Vec2::new(-0.5 * d, 0.0), Vec2::new(0.5 * d, 0.0), Vec2::new(0.0, 1000.0)],
                [Vec2::zero(); 3],
            ),
            [1.0, 1.0, 1.0],
        )
    }
    fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[((v.len() - 1) as f64 * q).round() as usize]
    }

    const N: usize = 64;
    println!("d_min over {N} perturbed collision times (initial separation 1 + k*1e-3)");
    println!("{:>8}{:>14}{:>14}{:>14}{:>16}", "eta", "p10", "median", "p90", "median ratio");
    let mut prev: Option<f64> = None;
    let mut ratios = Vec::new();
    for eta in [1e-2f64, 1e-3, 1e-4, 1e-5] {
        let mut ds = Vec::with_capacity(N);
        for k in 0..N {
            let (s, m) = perturbed(1.0 + (k as f64) * 1e-3);
            ds.push(az::integrate_az(s, &m, 1.0, 1, eta, 40_000_000, None).d_min_ref);
        }
        let (p10, med, p90) = (quantile(&mut ds, 0.10), quantile(&mut ds, 0.50), quantile(&mut ds, 0.90));
        let r = prev.map(|p| med / p);
        if let Some(r) = r {
            ratios.push(r);
        }
        let rs = r.map(|x| format!("{x:>16.3e}")).unwrap_or_else(|| format!("{:>16}", "-"));
        println!("{eta:>8.0e}{p10:>14.3e}{med:>14.3e}{p90:>14.3e}{rs}");
        prev = Some(med);
    }
    println!();
    println!("eta^2 predicts a median ratio of 1e-2 per decade. The p10-p90 spread within a");
    println!("single eta is wider than the shift between decades, which is exactly why a");
    println!("single trajectory cannot establish the scaling.");

    for r in &ratios {
        assert!(
            *r > 3e-3 && *r < 3e-2,
            "median d_min ratio {r:e} is not consistent with eta^2 (predicts 1e-2)"
        );
    }
}

/// **The mechanism behind the bit-identical d_min at eta=1e-4 and eta=1e-5: grid nesting.**
///
/// `dtau = eta * dt_left / (A0*B0)` is proportional to `eta` exactly, and every sub-interval
/// starts at `tau = 0`. So reducing `eta` by an integer factor produces a **strict refinement**
/// — every coarse sample point is still a sample point on the fine grid. `d_min` can
/// therefore never increase under refinement, and it only decreases when a genuinely new
/// point lands closer to the crossing.
///
/// For a regularised two-body collision the crossing is at a quarter period,
/// `tau_c = (pi/2) * dt_left/(A0*B0)`, so `tau_c/dtau = pi/(2 eta)` exactly. Write that as
/// `N + delta`. Refining by 10 gives `10N + 10delta`, and while `|10 delta| < 0.5` the
/// nearest sample is still the same *physical* point: the offset in `tau` is `delta*dtau`
/// either way. Identical offset, identical `d_min`.
///
/// So `d_min(eta)` is a **step function**, not a power law. Its envelope falls as `eta^2`;
/// individual decades can be flat.
///
/// This also eliminates the three alternative explanations:
///   - no clamping — at `eta = 1e-6` `d_min` reaches 1.7e-13, below `DIST_FLOOR = 1e-12`;
///   - it is recorded per RK4 step, not per sync boundary — the minimum lands at step 15708
///     of 24839 within a single sub-interval (`n_sync = 1`, so there is one registration);
///   - registration does not reset it — `d_min` accumulates across the whole run.
#[test]
fn d_min_is_a_step_function_because_the_sample_grids_nest() {
    let (_, m) = setup();
    println!("{:>8}{:>16}{:>18}{:>14}", "eta", "d_min", "tau_c/dtau = pi/2eta", "offset*dtau");
    for eta in [1e-4f64, 1e-5, 1e-6] {
        let (s, _) = setup();
        let o = az::integrate_az(s, &m, 1.0, 1, eta, 40_000_000, None);
        let x = std::f64::consts::PI / (2.0 * eta);
        let delta = x - x.round();
        // dtau in units where dt_left/(A0*B0) is absorbed; only the product matters here.
        println!("{eta:>8.0e}{:>16.7e}{x:>18.4}{:>14.3e}", o.d_min_ref, delta.abs() * eta);
    }
    println!();
    println!("The last column is the distance from the nearest sample to the crossing, in");
    println!("units of dt_left/(A0*B0). It is identical for 1e-4 and 1e-5 - which is why");
    println!("d_min is identical - and changes at 1e-6, which is why d_min drops there.");

    // The load-bearing claim: refinement never makes d_min worse.
    let mut prev = f64::INFINITY;
    for eta in [1e-2f64, 1e-3, 1e-4, 1e-5, 1e-6] {
        let (s, _) = setup();
        let d = az::integrate_az(s, &m, 1.0, 1, eta, 40_000_000, None).d_min_ref;
        assert!(d <= prev * (1.0 + 1e-9), "d_min rose under refinement: {d:e} > {prev:e}");
        prev = d;
    }
}

/// At the **production** step size, the §5 threshold is out of reach by two orders. Not a
/// defect — the trajectory simply is not sampled that finely.
#[test]
fn d_min_threshold_is_unreachable_at_production_eta() {
    let (s, m) = setup();
    let o = az::integrate_az(s, &m, 1.0, 1, 1e-2, 40_000_000, None);
    println!("at production eta = 1e-2: d_min = {:.3e} (the §5 threshold is 1e-10)", o.d_min_ref);
    println!("the gate above passes only at eta = 1e-4, 100x finer than production");
    assert!(o.d_min_ref > 1e-10, "unexpectedly reached the threshold at production eta");
}
