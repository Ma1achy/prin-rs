//! The logH march.
//!
//! Structurally the same loop as `heggie::driver`, with one thing removed and one added. Removed:
//! the regularised state and every reconstruction of it — the state here *is* Cartesian, so a
//! boundary sample is a read. Added: nothing, which is the point.
//!
//! # `n_sync` is a sampling cadence and nothing else
//!
//! As in Heggie. The loop records residuals, samples the escape rule and reads the tightest pair
//! at each boundary, and it does not re-register, because there is nothing to register into.
//! Where AZ's boundary does a round trip Cartesian -> LC -> Cartesian and re-freezes the energy,
//! this one does neither.
//!
//! The outcome machinery is not reimplemented: `outcome::collision_pairs_from`,
//! `escape_candidate_rule`, `closure`, `unbound` and `binary_id` are the same free functions the
//! other two drivers call, so all three share one definition of what an event is.
//!
//! # The predictive step limit DEFEATS this integrator at a collision, and it defaults off
//!
//! It was expected to be inert. It is not, and it fails in the direction the project already has
//! on record: *do not shrink the fictitious step at close approach — the time transformation
//! already shrinks the physical step, and shrinking `dtau` too drove `dt -> 1e-13` and produced
//! a false "this region is intractable".* Measured on the two-body radial collision:
//!
//! ```text
//!            limit f = 0.02 (the AZ/Heggie production value)   limit OFF
//!   RK4      40e6 steps, budget exhausted                      40e6 steps, budget exhausted
//!   KDK      40e6 steps, budget exhausted                      34034 steps, drift 1.1e-9
//! ```
//!
//! The arithmetic is direct. `dt/ds = 1/U`, so the limit `ds <= f d_min/(|v_rel| dt/ds)` is
//! `ds <= f d_min U/|v_rel|`. In free fall `U ~ 1/d` and `|v_rel| ~ 1/sqrt(d)`, so the bound
//! tends to `f sqrt(d)` while the unbounded step wants to *grow* as `1/d`. The limit therefore
//! forces `dt ~ d^{3/2}` — the free-fall time — and the step count to cross the encounter
//! diverges. Heggie escapes this because its `dt/dtau = R1 R2 R3/S^{3/2}` collapses much faster,
//! so the same bound never binds.
//!
//! **It is not simply harmful, and that is why it is a knob rather than a deletion.** On Burrau
//! at `t = 13`, `eta = 1e-2`, it moves drift `2.3e-6 -> 1.2e-9` for 36% more steps under RK4. So
//! it is a robustness/cost trade that **inverts** between an ordinary encounter and an exact
//! collision.
//!
//! # And the default stays at `0.02`, matching the other two, on purpose
//!
//! The tempting fix is to default this to `0.0`. It was written that way first and reverted,
//! because `EnsembleCfg::production()` asks for `StepLimit::Predictive` at `f = 0.02` and
//! `pixel.rs` maps that onto every occupant — so a `0.0` default here would make
//! `LhOpts::default()` disagree with what the ensemble path builds for the same config.
//! **That asymmetry is already on this project's record**: `AzOpts::default()` carried
//! `StepLimit::None` while `HgOpts::default()` carried the predictive limit, and a whole Phase 4
//! table compared a Heggie paying for a limit that AZ was not.
//!
//! The two failures are not equally bad. A uniform `0.02` makes a collision-heavy field exhaust
//! its budget, which lands in the `budget` column of every table and in `LhOut::finite` —
//! **visible**. A split default hides a knob difference inside two constructors of the same
//! options — **silent**. Prefer the visible one.
//!
//! So both arms are carried, the comparison harness runs both — *a knob held fixed for fairness
//! is a knob whose effect is unattributed*, and forcing a limit onto a method it defeats is not
//! fairness either — and every test whose subject is a collision pins `step_limit_f: 0.0` by
//! name rather than inheriting it.

use std::collections::VecDeque;

use crate::outcome::{self, Events};
use crate::physics::{energy, newton, Cart};
use crate::Real;

use super::hamiltonian::{denominators, residual, LhTime};
use super::state::LhState;
use super::step::{self, Stepper};
use super::system::LhSystem;

/// How the fictitious step `ds` is sized within a sync interval.
///
/// Named identically to [`az::DtauMode`](crate::integrate::az::DtauMode) and
/// [`HgDtauMode`](crate::integrate::heggie::HgDtauMode) so a three-way comparison can hold it
/// constant, and carrying the same warning: `PerStepRemaining` is Zeno by arithmetic and is a
/// measurement axis, never a candidate default.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LhDsMode {
    /// `ds = eta dt_left (K+B)`, sized once at interval entry and held.
    FixedPerInterval,
    /// `ds = eta (dt_left - s.t) (K+B)`, recomputed per step. Approaches the boundary
    /// geometrically and never reaches it. Kept as an axis because it looks right.
    PerStepRemaining,
    /// `ds = min(eta dt_left (K+B), ds_entry)`: `dt_left` held fixed, only `K+B` recomputed.
    #[default]
    PerStepInterval,
}

#[derive(Clone, Copy, Debug)]
pub struct LhOpts<T> {
    /// Whether the time transformation is in force. [`LhTime::None`] is the unregularised
    /// control, through this same code path.
    pub time: LhTime,
    /// RK4 (comparable to AZ and Heggie) or KDK (the method as designed).
    pub stepper: Stepper,
    pub ds_mode: LhDsMode,
    /// Land the final step of each interval **on** the boundary. A correctness property, not a
    /// preference: on AZ it is worth 1.06 -> 2.08 in measured convergence order, and on Heggie
    /// 1.03 -> 2.40.
    pub clamp_final_step: bool,
    /// `f` for the predictive limit `ds <= f d_min (K+B) / |v_rel|_max`. Zero disables it.
    ///
    /// Defaults to `0.02` as AZ and Heggie do, and **it is fatal at an exact collision** — see
    /// the module doc for why it is left uniform anyway, and why every test about a collision
    /// pins `0.0` by name.
    pub step_limit_f: T,
    /// Collision radius as a fraction of the initial hyperradius, fixed at `t = 0`.
    pub r_coll_frac: T,
    pub stop_on_event: bool,
    pub stop_on_escape: bool,
    pub escape_rule: outcome::EscapeRule<T>,
    /// The closure window in sync intervals. **It is a time**, so it must scale with `n_sync`.
    pub closure_k: usize,
    pub keep_boundary_shapes: bool,
    pub keep_drift_hist: bool,
}

impl<T: Real> Default for LhOpts<T> {
    fn default() -> Self {
        Self {
            time: LhTime::LogH,
            stepper: Stepper::Rk4,
            ds_mode: LhDsMode::PerStepInterval,
            clamp_final_step: true,
            step_limit_f: T::lit(0.02),
            r_coll_frac: T::zero(),
            stop_on_event: true,
            stop_on_escape: false,
            escape_rule: outcome::EscapeRule::Reference,
            closure_k: 1,
            keep_boundary_shapes: false,
            keep_drift_hist: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct LhOut<T> {
    pub state: Cart<T>,
    pub t: T,
    /// `|E(t) - E(0)| / |E(0)|` from the returned Cartesian state — the same quantity AZ and
    /// Heggie report, and here there is no second candidate, because there is no regularised
    /// energy to report instead. The AZ/Heggie pair were compared on two different measurements
    /// for a whole table once; this integrator cannot repeat that.
    pub drift: T,
    pub d_min: T,
    pub steps: usize,
    /// Force evaluations. **Not `steps * k`** — RK4 spends four per step and KDK one, so a
    /// derived count stops being true the moment two steppers appear in one table.
    pub force_evals: usize,
    pub finite: bool,
    pub budget_exhausted: bool,
    /// `max |K + B - U| / U` over the march. The energy defect normalised by the potential; see
    /// the module doc for why that is not an independent instrument.
    pub gamma_max: T,
    /// Largest physical step actually taken, as an `s.t` difference across one step. The
    /// tripwire that caught AZ advancing `2.209e128` against a sync interval of `0.4`.
    pub dt_max: T,
    pub n_overshoot: u32,
    /// A denominator went non-positive or non-finite before flooring. `K + B` is `U > 0` on
    /// shell, so this fires only when the trajectory has left the shell far enough that the
    /// transformation itself has failed — an *advance-anyway* site, recorded because otherwise
    /// it is invisible.
    pub den_degenerate: bool,
    pub events: Events<T>,
    pub t_end: T,
    pub tight: Vec<u8>,
    pub boundary_shapes: Vec<[T; 3]>,
    pub drift_hist: Vec<T>,
    pub closure_hist: Vec<T>,
    pub unbound_flags: Vec<[bool; 3]>,
    pub escape_flags: Vec<bool>,
    /// `d[1]/d[0]` over the sorted pair separations. The sibling `ref_tie` is deliberately
    /// absent: there is no argmax here to be near a tie of.
    pub tie_ratio: Vec<T>,
}

/// `dt/ds`, the physical time bought by a unit of fictitious time.
#[inline]
fn dt_ds<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T, time: LhTime) -> T {
    T::one() / denominators(sys, s, b, time).drift.max(T::TINY)
}

/// The predictive step limit, in `ds` units.
///
/// `dt = (dt/ds) ds`, so a physical bound `dt <= f d_min / |v_rel|_max` is
/// `ds <= f d_min / (|v_rel|_max (dt/ds))`. The three pair separations and the three relative
/// velocities are read straight off the state — there is no chart to reconstruct them from.
///
/// `+inf` when nothing is in force: the absence of a bound, not a step of zero.
#[inline]
fn predictive_ds<T: Real>(sys: &LhSystem<T>, s: &LhState<T>, b: T, time: LhTime, f: T) -> T {
    if f <= T::zero() {
        return T::infinity();
    }
    let d = newton::pair_dists(&s.r);
    let d_min = d.iter().fold(T::infinity(), |a, &x| a.min(x));
    let mut v_max = T::zero();
    for &(i, j) in crate::physics::PAIRS.iter() {
        v_max = v_max.max((s.v[j] - s.v[i]).norm());
    }
    let denom = v_max * dt_ds(sys, s, b, time);
    if denom > T::zero() && d_min.is_finite() && denom.is_finite() {
        f * d_min / denom
    } else {
        T::infinity()
    }
}

/// The full entry point. Signature matches `integrate_hg` field for field.
pub fn integrate_lh<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    n_sync: usize,
    eta: T,
    max_steps: usize,
    opts: &LhOpts<T>,
) -> LhOut<T> {
    let sys = LhSystem::new(*m);

    // The only thing frozen at t = 0, and it stays frozen: `dB/ds` is identically zero for
    // autonomous Newtonian gravity.
    let b = sys.b_of(&s0);
    let mut s = LhState::from_cart(&s0);
    let e0 = energy::energy(&s0.r, &s0.v, m, T::zero());

    let n = n_sync.max(1);
    let dt_sync = t_max / T::lit(n as f64);

    let r0 = energy::hyperradius(&s0.r, m);
    let r_coll = opts.r_coll_frac * r0;
    let rule = opts.escape_rule;
    let esc = |c: &Cart<T>, cl: Option<T>| outcome::escape_candidate_rule(c, m, rule, r0, cl);

    let kw = opts.closure_k.max(1);
    let mut nbuf: VecDeque<[T; 3]> = VecDeque::with_capacity(kw + 1);

    let mut out = LhOut {
        state: s0,
        t: T::zero(),
        drift: T::zero(),
        d_min: T::infinity(),
        steps: 0,
        force_evals: 0,
        finite: true,
        budget_exhausted: false,
        gamma_max: T::zero(),
        dt_max: T::zero(),
        n_overshoot: 0,
        den_degenerate: false,
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
    let mut cart = s0;

    'outer: for _ in 0..n {
        let dt_left = dt_sync;
        let tol = if opts.clamp_final_step { dt_left * T::LAND_EPS_REL } else { T::zero() };
        s.t = T::zero();

        let entry = dt_ds(&sys, &s, b, opts.time);
        if entry <= T::zero() || !entry.is_finite() {
            out.finite = false;
            break;
        }
        // `ds = eta dt_left / (dt/ds)`, the same shape the other two drivers use.
        let ds_entry = eta * dt_left / entry;

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

            let dens = denominators(&sys, &s, b, opts.time);
            if dens.degenerate() {
                out.den_degenerate = true;
            }
            let dtds = T::one() / dens.drift.max(T::TINY);

            let mut ds = match opts.ds_mode {
                LhDsMode::FixedPerInterval => ds_entry,
                LhDsMode::PerStepRemaining => eta * (dt_left - s.t).max(T::zero()) / dtds,
                LhDsMode::PerStepInterval => (eta * dt_left / dtds).min(ds_entry),
            };
            let lim = predictive_ds(&sys, &s, b, opts.time, opts.step_limit_f);
            if lim < ds {
                ds = lim;
            }
            // Land ON the boundary. A one-sided reduction of the final step only, using the same
            // `dt/ds` the sizing used — recomputing it here would let clamp and step disagree.
            if opts.clamp_final_step {
                ds = ds.min((dt_left - s.t).max(T::zero()) / dtds);
            }

            let before = s.t;
            let (next, evals) = step::step(&sys, &s, b, opts.time, opts.stepper, ds);
            s = next;
            out.steps += 1;
            out.force_evals += evals;
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

            let d = newton::pair_dists(&s.r);
            let dm = d.iter().fold(T::infinity(), |a, &x| a.min(x));
            if dm.is_finite() && dm < out.d_min {
                out.d_min = dm;
            }

            // Collision, sampled **inside** the loop: at `n_sync = 32` and `t_max = 13` the
            // boundaries are 0.4 apart and a close encounter passes between two of them unseen.
            // `pair_dists` is already in `PAIRS` order, so the free function takes it directly
            // and there is no reference-dependent index map to get wrong.
            if r_coll > T::zero() && out.events.collision.is_none() {
                let mask = outcome::collision_pairs_from(0, 1, 2, d[0], d[1], d[2], r_coll);
                if mask != 0 {
                    let tc = t_now + s.t;
                    out.events.collision = Some((mask, tc));
                    t_end.get_or_insert(tc);
                    if opts.stop_on_event {
                        break;
                    }
                }
            }

            let g = residual(&sys, &s, b);
            if g.is_finite() && g > out.gamma_max {
                out.gamma_max = g;
            }
        }

        t_now += s.t.min(dt_left);
        cart = s.to_cart();

        // ---- the boundary: sample and classify. No registration, because there is no chart. ----
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
            let mut dd = newton::pair_dists(&cart.r);
            dd.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            out.tie_ratio.push(dd[1] / dd[0].max(T::TINY));
        }

        // Instantaneous candidacy, recorded whether or not it has already fired — this is the
        // history a persistence guard reads, and it must not stop being written.
        let candidate_now = esc(&cart, cl);
        out.escape_flags.push(candidate_now.is_some());
        if out.events.escape.is_none() {
            if let Some(bd) = candidate_now {
                out.events.escape = Some((bd, t_now));
                t_end.get_or_insert(t_now);
            }
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
        out.finite = out.state.is_finite() && out.drift.is_finite();
    }
    out
}
