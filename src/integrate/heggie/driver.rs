//! The Heggie march.
//!
//! **What this file does not contain is the point.** `az::driver::integrate_az_opts` re-chooses a
//! reference body at every sync boundary, rebuilds `AzSystem`, round-trips the state through
//! Cartesian and re-registers. Measured, doubling that count at fixed step size moves the drift
//! field by **0.444 decades** — 6000x the effect of changing *which* reference is chosen.
//!
//! Here there is no reference body. `HgSystem` depends on nothing but the masses, `to_reg` is
//! called **once at `t = 0`**, and the sync loop is a *sampling* cadence: it records the residuals
//! and the closest approach, and it does not touch the state. The regularised march runs
//! uninterrupted from start to finish.
//!
//! Two free residuals come out of it, and one of them has no AZ analogue:
//!   - `|Gamma*|` against its largest term — zero along the exact trajectory, as AZ's is;
//!   - `|sum q_i| / max|q_i|` — Heggie's Eq. (9), an integral of the enlarged motion. AZ has
//!     nothing like it because AZ has no enlarged phase space.

use std::collections::VecDeque;

use crate::outcome::{self, Events};
use crate::physics::{energy, Cart};
use crate::Real;

use super::hamiltonian::{gamma_residual, HgTime};
use super::rk4;
use super::state::HgState;
use super::system::{cyc, HgSystem};

/// How `dtau` is sized within a sync interval. The same three arms as
/// [`DtauMode`](crate::integrate::az::DtauMode), named identically so a comparison between the
/// two integrators can hold this constant.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HgDtauMode {
    /// Sized once at interval entry and held.
    FixedPerInterval,
    /// `eta * (dt_left - s.t) / (dt/dtau)`. **Zeno by arithmetic** — the remaining time in the
    /// numerator gives `rem_{n+1} = rem_n (1 - eta)`, so the interval is approached geometrically
    /// and completed only through the relative landing tolerance. Named and kept as an axis,
    /// never a candidate; carried here so the AZ result reproduces on this integrator too.
    PerStepRemaining,
    /// Recomputed per step from the current `dt/dtau` with `dt_left` **held fixed**, capped at
    /// the interval's entry value. AZ's default and this one's.
    #[default]
    PerStepInterval,
}

#[derive(Clone, Copy, Debug)]
pub struct HgOpts<T> {
    /// Which of Heggie's time transformations. Default is Eqs. (22)-(24).
    pub time: HgTime,
    pub dtau_mode: HgDtauMode,
    /// Land the final step of each interval **on** the boundary rather than past it. A
    /// correctness property, not a preference: on AZ it is worth 1.06 -> 2.08 in measured
    /// convergence order.
    pub clamp_final_step: bool,
    /// **Secant-correct the step that lands on a sync boundary.** See `AzOpts::land_iterate` —
    /// same knob, same default, deliberately. Heggie's measured figure-eight order is 2.40
    /// clamped against 1.03 unclamped, and the `O(h^2)` landing residual is what holds it there.
    pub land_iterate: bool,
    /// Cap on secant iterations per landing step. Each is a real step and is counted.
    pub land_max_iters: usize,
    /// `f` for the predictive per-step limit, `dtau <= f d_min / (|dq/dt|_max dt/dtau)`. Zero
    /// disables it. AZ ships `0.02`.
    pub step_limit_f: T,
    /// Collision radius as a **fraction of the initial hyperradius**, canonical and fixed at
    /// `t = 0`. Zero disables the test. Same semantics as `AzOpts::r_coll_frac`, so a config can
    /// hand the same number to either integrator.
    pub r_coll_frac: T,
    /// Stop the march at the first collision.
    pub stop_on_event: bool,
    /// Stop the march at the first confirmed escape. **Default off**, as in AZ: closure certifies
    /// what escaped and is silent about whether the displayed shape has settled, and stopping
    /// moves 37% of pixels by up to 0.6 on a sphere of diameter 2.
    pub stop_on_escape: bool,
    pub escape_rule: outcome::EscapeRule<T>,
    /// Closure window, in sync intervals. It is a **time**, so this must scale with `n_sync`.
    pub closure_k: usize,
    /// Test escape every `n` RK4 steps as well as at boundaries. Zero is the reference cadence.
    /// **Vacuous under `Closure`**, which is only defined on the boundary series.
    pub escape_every: usize,
    /// Hold an in-loop escape provisional until the next boundary. Vacuous when
    /// `escape_every == 0`, which is the shipping path.
    pub escape_confirm: bool,
    pub keep_boundary_shapes: bool,
    pub keep_drift_hist: bool,
    /// Re-derive the frozen `h` from the state at each sync boundary.
    ///
    /// **Default off, and it exists to be measured rather than used.** Refreshing `h` is the one
    /// thing this driver could do that would reintroduce a boundary-dependent quantity into an
    /// otherwise uninterrupted march — so it is the control arm for the claim that removing
    /// re-registration is what matters. A control that changes nothing proves nothing, so its
    /// effect is reported, not assumed.
    pub refresh_h_at_boundary: bool,
    pub lc_stable: bool,
}

impl<T: Real> Default for HgOpts<T> {
    fn default() -> Self {
        Self {
            time: HgTime::default(),
            dtau_mode: HgDtauMode::default(),
            clamp_final_step: true,
            land_iterate: true,
            land_max_iters: 4,
            step_limit_f: T::lit(0.02),
            r_coll_frac: T::zero(),
            stop_on_event: true,
            stop_on_escape: false,
            escape_rule: outcome::EscapeRule::Reference,
            closure_k: 1,
            escape_every: 0,
            escape_confirm: true,
            keep_boundary_shapes: false,
            keep_drift_hist: false,
            refresh_h_at_boundary: false,
            lc_stable: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct HgOut<T> {
    pub state: Cart<T>,
    /// Physical time actually reached.
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|`, from the **returned Cartesian state**, exactly as AZ reports it.
    ///
    /// This is the honest one and the one to quote: it is the energy of the state every
    /// downstream consumer actually reads — `shape_vec`, `classify`, `error_ratio` — so a drift
    /// measured anywhere else is a drift for a state nobody sees.
    pub drift: T,
    /// The same quantity measured in the **enlarged variables**, before reconstruction.
    ///
    /// Carried because the two can differ by orders after a deep collision and the gap is the
    /// reconstruction's round-off, not the integration's error. An earlier cut of this driver
    /// reported only this one and read **4.4e-15 where the Cartesian state was at 1.2e-12** — a
    /// 280x under-report, and flattering in exactly the direction that would not have been
    /// questioned.
    pub drift_reg: T,
    /// Closest approach over all three pairs, sampled at every RK4 step.
    ///
    /// **All three, symmetrically.** AZ's `d_min_ref` is over its two regularised pairs and needs
    /// `d_min_true` beside it because the third side is treated differently. Here the three `q_i`
    /// are the three pairs and there is only one number to report.
    pub d_min: T,
    pub steps: usize,
    /// Extra steps spent on secant landing corrections. **Zero is ambiguous** — either the knob is
    /// off, or every landing hit tolerance first. Print the flag beside it.
    pub land_iters: u64,
    pub finite: bool,
    pub budget_exhausted: bool,
    /// Running max of `|Gamma*|` against its largest term.
    pub gamma_max: T,
    /// Running max of `|sum q_i| / max|q_i|`, Heggie's Eq. (9). **No AZ analogue.**
    pub sum_q_max: T,
    /// Largest physical step actually taken, as an `s.t` difference across one step.
    ///
    /// The tripwire. On AZ a single step advanced the clock by `2.209e128` against an interval of
    /// `0.4` and the march recorded a clean landing, because the clamp corrects the *clock* and
    /// cannot un-take the step. Nothing recorded it until this field existed.
    pub dt_max: T,
    /// Steps whose `s.t` exceeded twice their own interval. Must be zero.
    pub n_overshoot: u32,
    /// `R1 R2 R3` hit the `T::TINY` floor, so `dtau` divided by a fabricated denominator and the
    /// step advanced anyway. An advance-anyway site with no telemetry is indistinguishable from
    /// one that never fires.
    pub r_floored: bool,
    /// The terminating events, in the same shape AZ reports them, so `outcome::classify` can be
    /// called on either integrator's output without a second code path.
    pub events: Events<T>,
    /// Time of the first terminating event, or `t_max`.
    ///
    /// **No replay refinement.** AZ refines an escape's `t_end` by replaying the interval it
    /// fired in; measured, `at entry` is 1.0000 everywhere, so it never finds a crossing and
    /// returns the boundary time regardless. This returns the boundary time directly rather than
    /// duplicating the replay to reach the same answer — stated, because "not implemented" and
    /// "implemented and inert" are different claims and only the second is measured.
    pub t_end: T,
    /// Which pair is the tightest binary, at each sync boundary. `spread_event`'s input.
    pub tight: Vec<u8>,
    /// Shape vector at each boundary, when asked for.
    pub boundary_shapes: Vec<[T; 3]>,
    pub drift_hist: Vec<T>,
    /// Closure value at each boundary, `NaN` until the window fills.
    pub closure_hist: Vec<T>,
    /// Per-body unbound flags at each boundary. Rule-blind, unlike `escape_flags`.
    pub unbound_flags: Vec<[bool; 3]>,
    /// Whether the escape rule was satisfied at each boundary, recorded whether or not it has
    /// already fired — this is the history a persistence guard reads.
    pub escape_flags: Vec<bool>,
    /// `d[1]/d[0]` at each boundary: how close the two shortest sides are to a tie.
    ///
    /// Pure geometry, so it means the same thing here as in AZ. Its **sibling** `ref_tie`
    /// (`d[1]/d[2]`, a tie for the longest side) is deliberately absent: it exists in AZ to say
    /// how near the reference-body argmax came to flipping, and there is no argmax here for it
    /// to be about. Reporting it would be reporting a quantity with no referent.
    pub tie_ratio: Vec<T>,
    /// Physical `d|q_i|/dt` at the end, all three pairs. Kept because it is free — the march
    /// needs it for the step limit anyway.
    pub q_dot_max: T,
}

/// `dt/dtau` under the time transformation in force: `R1 R2 R3` for Eq. (20), divided by
/// `S^{3/2}` for Eq. (22).
#[inline]
fn dt_dtau<T: Real>(s: &HgState<T>, time: HgTime) -> T {
    let rp = s.r_prod();
    match time {
        HgTime::Product => rp,
        HgTime::SumPow32 { .. } => {
            let ss = s.s().max(T::TINY);
            rp / (ss * ss.sqrt())
        }
    }
}

/// `dq_i/dt` from Heggie's Eq. (7): `q_i' = p_i/mu_i - p_k/m_j - p_j/m_k`.
///
/// The index pattern is the one place this is easy to get wrong — the momentum with the
/// **further** cyclic index carries the reciprocal of the **nearer** mass, and vice versa.
/// `tests/heggie_march.rs` finite-differences it against Eq. (6) rather than trusting it.
#[inline]
pub fn q_dot<T: Real>(sys: &HgSystem<T>, p: &[crate::Vec2<T>; 3]) -> [crate::Vec2<T>; 3] {
    std::array::from_fn(|i| {
        let (j, k) = cyc(i);
        p[i] / sys.mu[i] - p[k] * sys.inv_m[j] - p[j] * sys.inv_m[k]
    })
}

/// The predictive step limit, in `dtau` units.
///
/// `dt = (dt/dtau) dtau`, so a physical bound `dt <= f d_min / |dq/dt|_max` is
/// `dtau <= f d_min / (|dq/dt|_max (dt/dtau))`. Every input is already in hand: the three `q_i`
/// **are** the three pair separations, with no reference body making one of them special.
///
/// `+inf` when nothing is in force — the absence of a bound, not a step of zero.
#[inline]
fn predictive_dtau<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, f: T, time: HgTime) -> T {
    if f <= T::zero() {
        return T::infinity();
    }
    let (q, p) = sys.phys_from_state(s);
    let v = q_dot(sys, &p);
    let d_min = q.iter().fold(T::infinity(), |a, x| a.min(x.norm()));
    let v_max = v.iter().fold(T::zero(), |a, x| a.max(x.norm()));
    let denom = v_max * dt_dtau(s, time);
    if denom > T::zero() && d_min.is_finite() && denom.is_finite() {
        f * d_min / denom
    } else {
        T::infinity()
    }
}

/// The full entry point.
///
/// `n_sync` is a **sampling and event** cadence here, and nothing more: the loop records
/// residuals, samples the escape rule and reads the tightest pair, and it does **not**
/// re-register. That is the whole structural difference from `integrate_az_opts`.
///
/// The outcome machinery is not reimplemented — `outcome::collision_pairs_from`,
/// `escape_candidate_rule`, `closure`, `unbound` and `binary_id` are free functions on `Cart`
/// that the AZ driver also calls, so both integrators share one definition of what an event is.
pub fn integrate_hg<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &HgOpts<T>,
) -> HgOut<T> {
    let mut sys = HgSystem::new(*m);
    sys.lc_stable = opts.lc_stable;

    // Registration. Once, at t = 0, and never again.
    let (mut s, h0) = sys.to_reg(&s0);
    let mut h = h0;
    let e0 = energy::energy(&s0.r, &s0.v, m, T::zero());

    let n = n_sync.max(1);
    let dt_sync = t_max / T::lit(n as f64);

    // Canonical and fixed at t = 0: a fraction of *this* trajectory's initial hyperradius,
    // evaluated before anything moves. A co-moving length makes the Hamiltonian time-dependent.
    let r0 = energy::hyperradius(&s0.r, m);
    let r_coll = opts.r_coll_frac * r0;
    let rule = opts.escape_rule;
    let esc = |c: &Cart<T>, cl: Option<T>| outcome::escape_candidate_rule(c, m, rule, r0, cl);

    // The closure window reads only its two ENDS, so a ring buffer of `k + 1` is all it needs.
    let kw = opts.closure_k.max(1);
    let mut nbuf: VecDeque<[T; 3]> = VecDeque::with_capacity(kw + 1);

    let mut out = HgOut {
        land_iters: 0,
        state: s0,
        t: T::zero(),
        drift: T::zero(),
        drift_reg: T::zero(),
        d_min: T::infinity(),
        steps: 0,
        finite: true,
        budget_exhausted: false,
        gamma_max: T::zero(),
        sum_q_max: T::zero(),
        dt_max: T::zero(),
        n_overshoot: 0,
        r_floored: false,
        q_dot_max: T::zero(),
        events: Events::default(),
        t_end: t_max,
        tight: Vec::with_capacity(n),
        boundary_shapes: Vec::with_capacity(if opts.keep_boundary_shapes { n } else { 0 }),
        drift_hist: Vec::with_capacity(if opts.keep_drift_hist { n } else { 0 }),
        closure_hist: Vec::with_capacity(n),
        unbound_flags: Vec::with_capacity(n),
        escape_flags: Vec::with_capacity(n),
        tie_ratio: Vec::with_capacity(n),
    };

    let mut t_now = T::zero();
    let mut t_end: Option<T> = None;
    let mut pending_escape: Option<(u8, T)> = None;
    let mut cart = s0;

    'outer: for _ in 0..n {
        let dt_left = dt_sync;
        let tol = if opts.clamp_final_step { dt_left * T::LAND_EPS_REL } else { T::zero() };
        s.t = T::zero();

        let entry = dt_dtau(&s, opts.time);
        if entry <= T::zero() || !entry.is_finite() {
            out.finite = false;
            break;
        }
        let dtau_entry = eta * dt_left / entry;

        loop {
            // `is_finite` explicitly: `NaN >= x` is false, so a diverged trajectory never
            // satisfies the loop exit and burns its entire step budget.
            if !s.is_finite() {
                out.finite = false;
                break 'outer;
            }
            if s.t >= dt_left - tol {
                break;
            }
            if out.steps >= max_steps {
                out.budget_exhausted = true;
                out.finite = false;
                break 'outer;
            }

            if s.r(0) < T::TINY || s.r(1) < T::TINY || s.r(2) < T::TINY {
                out.r_floored = true;
            }
            let d = dt_dtau(&s, opts.time).max(T::TINY);

            let mut dtau = match opts.dtau_mode {
                HgDtauMode::FixedPerInterval => dtau_entry,
                HgDtauMode::PerStepRemaining => eta * (dt_left - s.t).max(T::zero()) / d,
                HgDtauMode::PerStepInterval => (eta * dt_left / d).min(dtau_entry),
            };
            let lim = predictive_dtau(&sys, &s, opts.step_limit_f, opts.time);
            if lim < dtau {
                dtau = lim;
            }
            // Land ON the boundary. A one-sided reduction of the final step only, using the same
            // `d` the sizing used — recomputing it here would let the clamp and the step disagree.
            if opts.clamp_final_step {
                dtau = dtau.min((dt_left - s.t).max(T::zero()) / d);
            }

            let before = s.t;
            let s_before = s;
            s = rk4::step(&sys, &s, h, opts.time, dtau);
            out.steps += 1;

            // **THE LANDING CORRECTION.** The clamp above sizes the final step with `d`, the
            // `dt/dtau` read *before* the step — a first-order predictor whose residual is
            // `O(h^2)`. The clock is then set to the boundary and the state is not, which is what
            // holds the measured order at 2.40 rather than at RK4's four. The condition is the
            // loop's own exit test on the state just reached, so a mid-interval step is never
            // touched.
            if opts.clamp_final_step && opts.land_iterate && s.t >= dt_left - tol {
                let want = dt_left - before;
                let mut cur = dtau;
                for _ in 0..opts.land_max_iters {
                    let got = s.t - before;
                    if !got.is_finite() || got <= T::zero() || (got - want).abs() <= tol {
                        break;
                    }
                    // Secant on `t(dtau)`, with the free exact point `t(0) = 0`.
                    cur = cur * want / got;
                    s = rk4::step(&sys, &s_before, h, opts.time, cur);
                    out.steps += 1;
                    out.land_iters += 1;
                }
            }

            let took = s.t - before;
            if took > out.dt_max {
                out.dt_max = took;
            }
            if opts.clamp_final_step && s.t > dt_left * T::lit(2.0) {
                out.n_overshoot += 1;
                debug_assert!(
                    false,
                    "step overshot its interval: s.t = {} against dt_left = {}",
                    s.t, dt_left
                );
            }

            let (q, p) = sys.phys_from_state(&s);
            let dm = q.iter().fold(T::infinity(), |a, x| a.min(x.norm()));
            if dm < out.d_min {
                out.d_min = dm;
            }
            let v = q_dot(&sys, &p);
            out.q_dot_max = out.q_dot_max.max(v.iter().fold(T::zero(), |a, x| a.max(x.norm())));

            // Collision, sampled **inside** the loop rather than at boundaries: with `n_sync = 32`
            // and `t_max = 13` the boundaries are 0.4 apart and a close encounter passes between
            // two of them unseen. `q[i]` is the separation of the pair `q_i` joins, and unlike AZ
            // all three are on an equal footing, so there is no reference-dependent index map.
            if r_coll > T::zero() && out.events.collision.is_none() {
                let mut mask = 0u8;
                for (i, qi) in q.iter().enumerate() {
                    // `q[i]` joins bodies `i+1` and `i+2`; `outcome::pair_index` is keyed on
                    // `PAIRS`, which is `(0,1), (0,2), (1,2)`.
                    if qi.norm() < r_coll {
                        let (a, b) = cyc(i);
                        mask |= 1 << outcome::pair_index(a.min(b), a.max(b));
                    }
                }
                if mask != 0 {
                    let tc = t_now + s.t;
                    out.events.collision = Some((mask, tc));
                    t_end.get_or_insert(tc);
                    if opts.stop_on_event {
                        break;
                    }
                }
            }

            // The escape test at RK4-step resolution, when asked for. Off by default: this is the
            // reference's cadence and changing it changes results. Vacuous under `Closure`, which
            // is defined from the boundary series and has no finer resolution.
            if opts.escape_every > 0
                && !matches!(rule, outcome::EscapeRule::Closure(_))
                && out.events.escape.is_none()
                && pending_escape.is_none()
                && out.steps % opts.escape_every == 0
            {
                let c = sys.to_cartesian(&s);
                if let Some(b) = esc(&c, None) {
                    let te = t_now + s.t;
                    if opts.escape_confirm {
                        // Provisional. Do NOT break: the trajectory has to reach the next boundary
                        // for the condition to be re-tested. Breaking here is what turned 895
                        // transients into terminal escapes on AZ.
                        pending_escape = Some((b, te));
                    } else {
                        out.events.escape = Some((b, te));
                        t_end.get_or_insert(te);
                        if opts.stop_on_escape {
                            break;
                        }
                    }
                }
            }

            let g = gamma_residual(&sys, &s, h);
            if g.is_finite() && g > out.gamma_max {
                out.gamma_max = g;
            }
        }

        t_now += s.t.min(dt_left);
        cart = sys.to_cartesian(&s);

        // ---- the boundary: sample, classify, and do NOT re-register ----
        out.sum_q_max = out.sum_q_max.max(sys.sum_q_residual(&s));
        out.tight.push(outcome::binary_id(&cart));
        let n_now = crate::physics::shape::shape_vec(&cart.r, m);
        if opts.keep_boundary_shapes {
            out.boundary_shapes.push(n_now);
        }
        if opts.keep_drift_hist {
            let ek = energy::energy(&cart.r, &cart.v, m, T::zero());
            out.drift_hist.push(((ek - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs());
        }
        nbuf.push_back(n_now);
        if nbuf.len() > kw + 1 {
            nbuf.pop_front();
        }
        let cl = if nbuf.len() == kw + 1 { Some(outcome::closure(&n_now, &nbuf[0])) } else { None };
        out.closure_hist.push(cl.unwrap_or(T::nan()));
        out.unbound_flags.push([
            outcome::unbound(&cart, m, 0),
            outcome::unbound(&cart, m, 1),
            outcome::unbound(&cart, m, 2),
        ]);
        {
            let mut d = crate::physics::newton::pair_dists(&cart.r);
            d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            out.tie_ratio.push(d[1] / d[0].max(T::TINY));
        }

        // Instantaneous candidacy, recorded whether or not it has already fired — this is the
        // history a persistence guard reads, and it must not stop being written.
        let candidate_now = esc(&cart, cl);
        out.escape_flags.push(candidate_now.is_some());

        // Confirm or discard a provisional in-loop escape. The committed time is the FIRST
        // crossing: the guard decides whether the event was real, not when it happened.
        if let Some((b, te)) = pending_escape.take() {
            if candidate_now.is_some() {
                out.events.escape = Some((b, te));
                t_end.get_or_insert(te);
            }
        }
        if out.events.escape.is_none() {
            if let Some(b) = candidate_now {
                out.events.escape = Some((b, t_now));
                t_end.get_or_insert(t_now);
            }
        }

        if opts.refresh_h_at_boundary {
            h = sys.energy_of(&s);
        }

        if out.budget_exhausted
            || (opts.stop_on_event && out.events.collision.is_some())
            || (opts.stop_on_escape && out.events.escape.is_some())
        {
            break;
        }
    }

    out.t = t_now;
    out.t_end = t_end.unwrap_or(t_max);
    if out.finite {
        out.state = cart;
        let e1 = energy::energy(&cart.r, &cart.v, m, T::zero());
        out.drift = ((e1 - e0) / e0.abs().max(T::DRIFT_FLOOR)).abs();
        // From the regularised state directly, never through the Cartesian reconstruction.
        let e1r = sys.energy_of(&s);
        out.drift_reg = ((e1r - h0) / h0.abs().max(T::DRIFT_FLOOR)).abs();
        out.finite = out.state.is_finite() && out.drift.is_finite();
    }
    out
}
