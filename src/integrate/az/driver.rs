//! The sync loop: choose a reference body, register, RK4-march to the boundary, return to
//! Cartesian, repeat.
//!
//! Transcribed from `reference/tb_az.py:integrate_az`.

use crate::outcome::{collision_pairs_from, Events};
use crate::physics::{energy, Cart};
use crate::Real;

use super::hamiltonian::gamma_residual;
use super::reference_body::{choose_reference, triple, RefPolicy};
use super::rk4;
use super::system::AzSystem;

#[derive(Clone, Debug)]
pub struct AzOut<T> {
    pub state: Cart<T>,
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|`, Cartesian and unsoftened, exactly as the reference.
    pub drift: T,
    /// `min(|R1|, |R2|)` — the two **regularised** pairs only. This is what the reference
    /// tracks, and it is what the cross-check compares.
    pub d_min_ref: T,
    /// `min` over all three pairs, including the unregularised side. BRIEF §4 defines
    /// `d_min` this way; the reference does not compute it. The gap between the two
    /// measures how well the reference-switching cadence tracks encounters (NOTES §2.1).
    pub d_min_true: T,
    pub switches: u32,
    /// `max |H - E| / |E|` over the run — the free residual from `Gamma == 0` along the
    /// trajectory.
    pub gamma_max: T,
    /// Reference body chosen at each sync boundary. Compared as a column in the
    /// cross-check: a near-tie in `argmax` broken differently by the two implementations
    /// fails the comparison while looking exactly like a transcription error.
    pub refs: Vec<u8>,
    pub steps: usize,
    pub finite: bool,
    pub budget_exhausted: bool,
    /// Termination events, per BRIEF §2.4. Collisions are sampled **inside** the RK4 loop;
    /// escape at sync boundaries, which is where the reference samples it.
    pub events: Events<T>,
    /// Identity of the **currently tightest pair** at each completed sync boundary, as an
    /// index into [`crate::physics::PAIRS`].
    ///
    /// This is what `spread_event` is built on, and it is deliberately *not* the terminal
    /// outcome. A terminal label is terminal-grain: early in the march nothing has
    /// terminated, so every copy agrees, so a spread built on it reports maximum confidence
    /// at exactly the playhead where least is known. The tightest-pair identity is defined at
    /// every playhead, needs no gate, and moves earlier. See NOTES §2.8.
    pub tight: Vec<u8>,
    /// Time of the first terminating event, or the time reached if none fired. Distinct from
    /// `t` only when `stop_on_event` is off, where the run continues past the event.
    pub t_end: T,
}

/// The options that grew past what a positional argument list can carry legibly.
#[derive(Clone, Copy, Debug)]
pub struct AzOpts<'a, T> {
    pub forced_refs: Option<&'a [u8]>,
    /// Conditioned inverse LC branch. `false` reproduces the reference's original.
    pub lc_stable: bool,
    /// **Canonical**: a fraction of the initial hyperradius `R`, evaluated once at `t = 0`
    /// from this trajectory's own initial condition and never updated. Zero disables
    /// collision detection.
    ///
    /// A fraction, not a length. An absolute `r_coll` would break the scale invariance the
    /// project quotients out — measured, the same physical system gave answers differing by
    /// 1.66x purely from an arbitrary overall size. A *co-moving* one would make the
    /// Hamiltonian time-dependent and destroy energy conservation.
    pub r_coll_frac: T,
    /// Stop at the first terminating event, per §2.4. With it off the run continues to
    /// `t_max` and the event is recorded but not acted on — which is what the reference does,
    /// and what keeps every copy's continuous fields evaluated at a common playhead.
    pub stop_on_event: bool,
}

impl<T: Real> Default for AzOpts<'_, T> {
    fn default() -> Self {
        Self { forced_refs: None, lc_stable: true, r_coll_frac: T::zero(), stop_on_event: true }
    }
}

/// `forced_refs`, when supplied, overrides the per-sync choice — this is how the shared
/// reference policy is applied, by handing every copy the nominal copy's choices.
pub fn integrate_az<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    forced_refs: Option<&[u8]>,
) -> AzOut<T> {
    integrate_az_lc(s0, m, t_max, n_sync, eta, max_steps, forced_refs, true)
}

/// As [`integrate_az`], with an explicit choice of inverse LC branch.
///
/// `lc_stable = false` reproduces the reference's branch exactly and is what the cross-check
/// uses; `true` is the production kernel. The two produce genuinely different trajectories —
/// the stable branch is more accurate, so it necessarily stops agreeing bit-for-bit with the
/// reference. Keeping both makes that a measured difference rather than a lost gate.
#[allow(clippy::too_many_arguments)]
pub fn integrate_az_lc<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    forced_refs: Option<&[u8]>,
    lc_stable: bool,
) -> AzOut<T> {
    // No r_coll and no early stop: these two entry points are the reference-matching path,
    // and the reference has no termination logic at all. Events are still *recorded* — that
    // costs nothing and changes no arithmetic — but nothing acts on them.
    integrate_az_opts(
        s0, m, t_max, n_sync, eta, max_steps,
        &AzOpts { forced_refs, lc_stable, r_coll_frac: T::zero(), stop_on_event: false },
    )
}

/// The full entry point.
#[allow(clippy::too_many_arguments)]
pub fn integrate_az_opts<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &AzOpts<T>,
) -> AzOut<T> {
    let (forced_refs, lc_stable) = (opts.forced_refs, opts.lc_stable);
    let mut cart = s0;
    let e0 = energy::energy(&s0.r, &s0.v, m, T::zero());
    let mut t = T::zero();
    let mut d_min_ref = T::infinity();
    let mut d_min_true = T::infinity();
    let mut gamma_max = T::zero();
    let mut switches = 0u32;
    let mut prev_ref: Option<usize> = None;
    let mut refs = Vec::with_capacity(n_sync);
    let mut tight = Vec::with_capacity(n_sync);
    let mut total_steps = 0usize;
    let mut finite = true;
    let mut budget_exhausted = false;
    let mut events = Events::default();
    let mut t_end: Option<T> = None;

    // Canonical and fixed at t=0: a fraction of *this* trajectory's initial hyperradius,
    // evaluated once, before anything moves. Never recomputed from the instantaneous
    // configuration — a co-moving length makes the Hamiltonian time-dependent.
    let r_coll = opts.r_coll_frac * energy::hyperradius(&s0.r, m);

    for kk in 0..n_sync {
        let t_target = T::lit((kk + 1) as f64) * t_max / T::lit(n_sync as f64);

        // `f.get(kk)`, not `f[kk]`. Since Step 5b the nominal copy can terminate early, so its
        // `refs` record is shorter than `n_sync` and the shared policy has no opinion past
        // that point. Falling back to the per-copy choice is the only defensible reading:
        // sharing applies where the nominal has a choice to share. (This indexed out of
        // bounds the first time the shared policy met a terminating run.)
        let a = match forced_refs.and_then(|f| f.get(kk)) {
            Some(&f) => f as usize,
            None => choose_reference(&cart.r),
        };
        refs.push(a as u8);
        if let Some(p) = prev_ref {
            if p != a {
                switches += 1;
            }
        }
        prev_ref = Some(a);

        // Matches the reference's row selection `t < t_target - 1e-15`.
        if t >= t_target - T::SYNC_EPS {
            continue;
        }

        let (ab, bb, cb) = triple(a);
        let sys = if lc_stable {
            AzSystem::new(ab, bb, cb, *m)
        } else {
            AzSystem::new(ab, bb, cb, *m).with_reference_lc()
        };
        let (mut s, e) = sys.to_reg(&cart);
        let dt_left = t_target - t;

        // dt = A*B*dtau, so the A*B factor ALREADY shrinks the physical step at close
        // approach — that is what regularisation buys. dtau must therefore NOT also shrink
        // with the separation: doing so drove dt -> 1e-13 and exhausted the step budget,
        // producing a false "this region is intractable". Size dtau so the FIRST physical
        // step is a fixed fraction of the interval; the rest follows.
        let a0 = s.a().max(T::TINY);
        let b0 = s.b().max(T::TINY);
        let dtau = eta * dt_left / (a0 * b0);

        let mut steps = 0usize;
        loop {
            // NaN >= x is false, so a non-finite trajectory never satisfies `done` and the
            // loop would burn the whole budget (measured 354 s against 3 s nominal).
            // Test is_finite explicitly.
            let bad = !s.is_finite();
            if bad {
                finite = false;
            }
            if s.t >= dt_left || bad {
                break;
            }
            if steps >= max_steps {
                budget_exhausted = true;
                break;
            }

            s = rk4::step(&sys, &s, e, dtau);
            steps += 1;

            let (r1, r2, _, _) = sys.phys_from_state(&s);
            let d1 = r1.norm();
            let d2 = r2.norm();
            let d3 = (r2 - r1).norm();
            // Only fold finite separations in. `min(x, NaN)` is NaN in the reference, which
            // silently poisons d_min for every consumer downstream (NOTES §2.3). The copy
            // is still kept — `finite` records the outcome.
            let m_ref = d1.min(d2);
            if m_ref.is_finite() && m_ref < d_min_ref {
                d_min_ref = m_ref;
            }
            let m_true = m_ref.min(d3);
            if m_true.is_finite() && m_true < d_min_true {
                d_min_true = m_true;
            }

            // Sampled here, inside the RK4 loop, not at sync boundaries: with n_sync = 32 and
            // t_max = 13 the boundaries are 0.4 apart, and a close encounter passes between
            // two of them entirely unseen.
            if events.collision.is_none() && r_coll > T::zero() {
                let mask = collision_pairs_from(ab, bb, cb, d1, d2, d3, r_coll);
                if mask != 0 {
                    let tc = t + s.t;
                    events.collision = Some((mask, tc));
                    t_end.get_or_insert(tc);
                    if opts.stop_on_event {
                        break;
                    }
                }
            }

            let g = gamma_residual(&sys, &s, e);
            if g.is_finite() && g > gamma_max {
                gamma_max = g;
            }
        }
        total_steps += steps;

        cart = sys.to_cartesian(&s);
        // The overshoot past the boundary is clipped in the time bookkeeping only; the
        // state written back is the overshot one. Sub-step interpolation is not done. This
        // is the reference's behaviour and part of what the cross-check measures.
        t += s.t.min(dt_left);

        tight.push(crate::outcome::binary_id(&cart));

        // The escape test runs at the sync boundary, where the state is Cartesian and every
        // trajectory shares a playhead — the reference's cadence, transcribed.
        if events.escape.is_none() {
            if let Some(b) = crate::outcome::escape_candidate(&cart, m) {
                events.escape = Some((b, t));
                t_end.get_or_insert(t);
            }
        }

        if budget_exhausted || (opts.stop_on_event && (events.collision.is_some() || events.escape.is_some())) {
            break;
        }
    }

    let e1 = energy::energy(&cart.r, &cart.v, m, T::zero());
    let drift = ((e1 - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs();

    AzOut {
        state: cart,
        t,
        drift,
        d_min_ref,
        d_min_true,
        switches,
        gamma_max,
        refs,
        tight,
        steps: total_steps,
        finite,
        budget_exhausted,
        events,
        t_end: t_end.unwrap_or(t),
    }
}

/// Convenience: run with the shared-reference policy by first computing the nominal copy's
/// choices, then forcing them on the copy.
pub fn integrate_with_policy<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    policy: RefPolicy,
    nominal_refs: Option<&[u8]>,
) -> AzOut<T> {
    let forced = match (policy, nominal_refs) {
        (RefPolicy::Shared, Some(f)) => Some(f),
        _ => None,
    };
    integrate_az(s0, m, t_max, n_sync, eta, max_steps, forced)
}
