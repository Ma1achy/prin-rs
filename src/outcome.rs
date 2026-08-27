//! Outcome classification.
//!
//! Step 5a needs only the **legacy** classifier: `spread_event` is defined over outcome
//! classes, and `tb.classify` is the one class labelling with a reference to check against.
//!
//! BRIEF §2.4's 3-bit `state` / 2-bit `detail` encoding follows below. `classify_legacy` and
//! `classify` are kept apart on purpose: the first has a reference and the second largely does
//! not, and the ported part should be reviewable separately from the invented part.

use crate::physics::{newton, Cart, G, PAIRS, THIRD};
use crate::Vec2;
use crate::Real;

/// `tb.classify`, transcribed.
///
/// Returns the index of the escaping body (`0`, `1`, `2`) or `3` for still-bound. The
/// candidate escaper is the body not in the *tightest* pair; it is labelled escaping only if
/// its specific orbital energy relative to the other two is positive **and** it is receding.
/// Bodies not selected as the candidate can never be labelled escaping — a property of the
/// reference, transcribed rather than fixed.
pub fn classify_legacy<T: Real>(s: &Cart<T>, m: &[T; 3]) -> u8 {
    let d = newton::pair_dists(&s.r);
    let mut tight = 0usize;
    for k in 1..3 {
        if d[k] < d[tight] {
            tight = k;
        }
    }
    let b = THIRD[tight];
    let (o1, o2) = {
        let others: Vec<usize> = (0..3).filter(|&k| k != b).collect();
        (others[0], others[1])
    };

    let m_bin = m[o1] + m[o2];
    let rc = (s.r[o1] * m[o1] + s.r[o2] * m[o2]) / m_bin;
    let vc = (s.v[o1] * m[o1] + s.v[o2] * m[o2]) / m_bin;
    let dr = s.r[b] - rc;
    let dv = s.v[b] - vc;
    let dist = dr.norm();

    let half = T::lit(0.5);
    let spec = half * dv.norm_sq() - T::lit(G) * m_bin / dist.max(T::CLASSIFY_FLOOR);
    let receding = dr.dot(dv) > T::zero();

    if spec > T::zero() && receding {
        b as u8
    } else {
        3
    }
}

/// Which pair is tightest at the final state, as an index into [`PAIRS`]. Cheap, and useful
/// for interpreting `spread_event` — it says *which* binary formed, not merely that the
/// copies disagreed.
pub fn binary_id<T: Real>(s: &Cart<T>) -> u8 {
    let d = newton::pair_dists(&s.r);
    let mut best = 0usize;
    for k in 1..3 {
        if d[k] < d[best] {
            best = k;
        }
    }
    best as u8
}

// ---------------------------------------------------------------------------
// Step 5b: BRIEF §2.4's encoding.
//
// PORTED vs INVENTED, because the two carry different weight in review:
//
//   PORTED    the escape arm — `escape_candidate` transcribes the sync-boundary test from
//             `reference/tb_all_az.py:59-75`, including its `mb` (not `G*mb`, since G = 1),
//             its `1e-12` distance floor, and its restriction to the body outside the
//             *tightest* pair. Same test, same cadence.
//
//             **And the numpy tree is not the only reference.** The GLSL
//             (`frag.glsl:104`) tests a third condition the numpy form does not have —
//             `dist > r_esc` — on all three bodies. `escape_candidate_gated` carries both
//             arms; which one runs is a configuration, not a transcription.
//
//   INVENTED  everything else. `r_coll` appears nowhere in the numpy tree, so the collision
//             arm, the >=2-pair triple rule, `triple_ejection`, and the packing are new
//             construction with no oracle. They are validated by property tests only.

use crate::physics::energy;

/// BRIEF §2.4's 3-bit `state`.
///
/// **A framing note, reported rather than papered over.** §2.4's table gives conditions for
/// escape, collision, triple collision, triple ejection and `running` ("still bound at
/// `t_max`"), while the state list also contains `bounded`. As written, `bounded` and
/// `running` describe the same situation and one of the six is unreachable.
///
/// Resolved here by giving each a distinct job rather than leaving a dead state: [`Bounded`]
/// is reaching `t_max` with nothing having fired, and [`Running`] is *not* reaching it — the
/// step budget ran out, so the trajectory is genuinely still running and its final state is
/// not a terminal answer. That reading keeps all six reachable and keeps the difference
/// between "we integrated to the horizon" and "we stopped early" visible in the dump, which
/// otherwise has nowhere to record it. If the intended reading was the literal one, this is
/// the line to change.
///
/// [`Bounded`]: State::Bounded
/// [`Running`]: State::Running
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
pub enum State {
    Escape = 0,
    Bounded = 1,
    Collision = 2,
    Running = 3,
    SimFailed = 4,
    DecodeFailed = 5,
}

impl State {
    pub fn name(self) -> &'static str {
        match self {
            State::Escape => "escape",
            State::Bounded => "bounded",
            State::Collision => "collision",
            State::Running => "running",
            State::SimFailed => "sim_failed",
            State::DecodeFailed => "decode_failed",
        }
    }

    pub fn from_bits(b: u8) -> Option<State> {
        Some(match b {
            0 => State::Escape,
            1 => State::Bounded,
            2 => State::Collision,
            3 => State::Running,
            4 => State::SimFailed,
            5 => State::DecodeFailed,
            _ => return None,
        })
    }
}

/// `state` in the high 3 bits, `detail` in the low 2.
///
/// `detail = 3` means **"all three"** for both arms — triple collision under [`State::Collision`]
/// and triple ejection under [`State::Escape`]. One rule, both arms, per §2.4.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Outcome {
    pub state: State,
    pub detail: u8,
}

impl Outcome {
    pub fn new(state: State, detail: u8) -> Self {
        debug_assert!(detail < 4, "detail is 2 bits");
        Self { state, detail }
    }

    pub fn pack(self) -> u8 {
        ((self.state as u8) << 2) | (self.detail & 3)
    }

    pub fn unpack(b: u8) -> Option<Outcome> {
        State::from_bits(b >> 2).map(|state| Outcome { state, detail: b & 3 })
    }

    /// True for the two "all three" outcomes.
    pub fn is_triple(self) -> bool {
        self.detail == 3 && matches!(self.state, State::Collision | State::Escape)
    }
}

/// Bitmask over [`PAIRS`] of pairs separated by less than `r_coll`.
pub fn collision_pairs<T: Real>(r: &[Vec2<T>; 3], r_coll: T) -> u8 {
    if !(r_coll > T::zero()) {
        return 0;
    }
    let d = newton::pair_dists(r);
    let mut mask = 0u8;
    for (k, dk) in d.iter().enumerate() {
        if *dk < r_coll {
            mask |= 1 << k;
        }
    }
    mask
}

/// `detail` for a collision mask: the pair index for one pair, `3` for **two or more**.
///
/// **The >=2 rule, not "all three".** By the triangle inequality `|AB| < r_coll` and
/// `|AC| < r_coll` force `|BC| < 2 r_coll`, so the two-pair state is reachable and is already
/// a near-triple. Requiring all three would label it an ordinary binary collision, which is
/// the misclassification the rule exists to prevent.
pub fn collision_detail(mask: u8) -> u8 {
    match mask.count_ones() {
        0 => unreachable!("collision_detail on an empty mask"),
        1 => mask.trailing_zeros() as u8,
        _ => 3,
    }
}

/// Which escape condition the run is using. **Three references, and they do not agree.**
///
/// | rule | conditions | bodies | source |
/// |---|---|---|---|
/// | [`Reference`](EscapeRule::Reference) | unbound && receding | tightest-pair only | `reference/tb_all_az.py:73-74` |
/// | [`Distance`](EscapeRule::Distance) | + `dist > r_esc` | all three | `frag.glsl:104` |
/// | [`Closure`](EscapeRule::Closure) | `|dn| < tau` && unbound | all three | `reference/escape_criterion.py` |
///
/// The energy arm is identical in all three (`0.5 dv.dv - mb/max(d, floor)`, `G = 1`); only the
/// *gating* differs. `Reference` is what the cross-check measures and it is byte-preserved.
///
/// # Why `Closure` drops two arms rather than adding a third
///
/// `spec > 0 && receding` is **not absorbing** — during a close encounter it is transiently true
/// while the body is still deep inside the system, and 0 of 895 such firings in `deep interior`
/// were still unbound one boundary later. `Distance` guards that geometrically. `Closure` guards
/// it by measuring the thing that actually distinguishes an escape: on a hierarchical escape the
/// binary separation stays bounded while `lambda` grows linearly, so `tan(alpha) ~ t`,
/// `n_0 = cos(2 alpha) -> -1` with error `~1/t^2`, and `|dn/dt| ~ 1/t^3`. **The shape vector
/// converges to a pole.** Escape is a limit being approached, not a threshold being crossed.
///
/// Once closure and energy both hold the body *is* receding and far away by construction, so
/// `receding` and `dist > r_esc` add nothing — measured identical to the digit. Three tuned
/// constants become one, and `tau` is set from a measured gap rather than picked.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EscapeRule<T> {
    /// The numpy reference: unbound and receding, tightest-pair candidate only.
    Reference,
    /// The GLSL: adds a distance gate, **canonical** — a fraction of the initial hyperradius `R`
    /// fixed at `t = 0`, which the driver multiplies out. An absolute length would break the
    /// scale invariance the project quotients out (BRIEF §2.5).
    Distance(T),
    /// Closure and energy. `tau` is a chord on the **unit sphere**, so it is dimensionless and
    /// needs no canonical conversion — unlike `r_coll`, `epsilon` and `r_esc`. The *window* does,
    /// and it lives in the driver as a count of sync boundaries.
    Closure(T),
}

impl EscapeRule<f64> {
    /// Lift the configured (always `f64`) rule into the integrator's precision.
    ///
    /// `tau` is a chord on the unit sphere and `Distance`'s payload is a fraction of `R`; both
    /// are dimensionless, so this is a cast and not a conversion.
    pub fn lift<T: Real>(self) -> EscapeRule<T> {
        match self {
            EscapeRule::Reference => EscapeRule::Reference,
            EscapeRule::Distance(x) => EscapeRule::Distance(T::lit(x)),
            EscapeRule::Closure(x) => EscapeRule::Closure(T::lit(x)),
        }
    }
}

/// The default closure threshold.
///
/// **Not a picked number: the geometric midpoint of a measured gap.** Closure separates escapers
/// from bound trajectories by 383x on the config chart — `7.04e-05` against `2.70e-02`, stable
/// across `t = 25-30` — and any value in the middle two orders gives the same answer. The
/// reference carries `tau = 1e-3` for that reason and this matches it.
///
/// Setting it from the *distribution* rather than by eye is the whole point: an absolute cutoff
/// of `2e-3` picked by eye is what dismissed closure the first time, because it sits **inside**
/// the bound population's range. `examples/escape_closure.rs` re-measures the gap on this
/// implementation's own regions and reports the separation rather than adopting the number.
pub const CLOSURE_TAU: f64 = 1e-3;

/// `|n_now - n_past|`, the chord between the two **ends** of the closure window.
///
/// The reference buffers `nbuf` samples but reads only `buf[-1]` and `buf[0]`
/// (`escape_criterion.py`), so the window interior is never used. That is why sampling at sync
/// boundaries is a transcription rather than an approximation: at `t_max = 13, n_sync = 32` the
/// realised window is 0.406 against the reference's 0.400.
///
/// Non-finite in, non-finite out — a triple collision gives a NaN `shape_vec` and `shape.rs`
/// deliberately does not floor it. A NaN closure fails `< tau` and cannot fire, which is the
/// right answer: an undetermined shape has not settled.
pub fn closure<T: Real>(n_now: &[T; 3], n_past: &[T; 3]) -> T {
    let mut s = T::zero();
    for k in 0..3 {
        let d = n_now[k] - n_past[k];
        s = s + d * d;
    }
    s.sqrt()
}

/// Body `b`'s specific energy relative to the barycentre of the other two, and its separation.
///
/// `0.5 dv.dv - mb/max(d, floor)`. The reference floors the distance at `1e-12` and writes `mb`,
/// not `G*mb`; `G` is 1. Shared by all three rules — only the gating differs.
fn rel_two_body<T: Real>(s: &Cart<T>, m: &[T; 3], b: usize) -> (T, Vec2<T>, Vec2<T>) {
    let o: Vec<usize> = (0..3).filter(|&k| k != b).collect();
    let mb = m[o[0]] + m[o[1]];
    let rc = (s.r[o[0]] * m[o[0]] + s.r[o[1]] * m[o[1]]) / mb;
    let vc = (s.v[o[0]] * m[o[0]] + s.v[o[1]] * m[o[1]]) / mb;
    let dr = s.r[b] - rc;
    let dv = s.v[b] - vc;
    let spec = T::lit(0.5) * dv.norm_sq() - T::lit(G) * mb / dr.norm().max(T::DIST_FLOOR);
    (spec, dr, dv)
}

/// Is body `b` unbound from the other two? The one arm every rule shares.
pub fn unbound<T: Real>(s: &Cart<T>, m: &[T; 3], b: usize) -> bool {
    rel_two_body(s, m, b).0 > T::zero()
}

/// Does body `b` satisfy the **gated** escape test — far, receding, unbound?
///
/// Three conditions, in the order the GLSL writes them. `r_esc` is an absolute length (the caller
/// multiplies the canonical fraction by `R`); pass zero and the distance arm is vacuous, which
/// recovers the numpy reference exactly.
fn escapes<T: Real>(s: &Cart<T>, m: &[T; 3], b: usize, r_esc: T) -> bool {
    let (spec, dr, dv) = rel_two_body(s, m, b);
    dr.norm() > r_esc && dr.dot(dv) > T::zero() && spec > T::zero()
}

/// The escaping body, or `None`.
///
/// # Body ordering, and it differs by rule
///
/// `Reference` can label only the body outside the **tightest** pair — a property of the numpy
/// reference, transcribed rather than fixed. `Distance` tries that candidate first and then the
/// others, so its single-candidate answer is returned unchanged whenever it fires. `Closure`
/// scans **0, 1, 2 in index order**, because the reference's `b = np.argmax(fire, -1)` returns
/// the lowest firing index. That is a real divergence from the tightest-pair ordering and it is
/// transcribed, not chosen.
///
/// # `closure_now`
///
/// `None` means the window is not yet full, and under `Closure` **nothing can fire** — correct,
/// since nothing has settled at `t ~ 0`. It is ignored by the other two rules.
///
/// # Why closure gates once and energy gates per body
///
/// `|dn|` is a property of the whole configuration; `spec` is per body. The reference broadcasts
/// `dn[..., None]` against an `E` of shape `(..., 3)` — one scalar gate, then a per-body scan.
pub fn escape_candidate_rule<T: Real>(
    s: &Cart<T>,
    m: &[T; 3],
    rule: EscapeRule<T>,
    r_scale: T,
    closure_now: Option<T>,
) -> Option<u8> {
    let d = newton::pair_dists(&s.r);
    let mut tight = 0usize;
    for k in 1..3 {
        if d[k] < d[tight] {
            tight = k;
        }
    }
    let first = THIRD[tight];
    match rule {
        EscapeRule::Reference => escapes(s, m, first, T::zero()).then_some(first as u8),
        EscapeRule::Distance(frac) => {
            let r_esc = frac * r_scale;
            if escapes(s, m, first, r_esc) {
                return Some(first as u8);
            }
            (0..3)
                .find(|&b| b != first && escapes(s, m, b, r_esc))
                .map(|b| b as u8)
        }
        EscapeRule::Closure(tau) => {
            // NaN fails `<` and so cannot fire: an undetermined shape has not settled.
            if !closure_now.is_some_and(|c| c < tau) {
                return None;
            }
            (0..3).find(|&b| unbound(s, m, b)).map(|b| b as u8)
        }
    }
}

/// The numpy reference's form: no distance gate, tightest-pair candidate only.
///
/// Kept as the reference-matching path — `integrate_az_lc` hardcodes it and the cross-check
/// measures it. See [`EscapeRule`] for what it is missing and why.
pub fn escape_candidate<T: Real>(s: &Cart<T>, m: &[T; 3]) -> Option<u8> {
    escape_candidate_rule(s, m, EscapeRule::Reference, T::zero(), None)
}

/// All three pairs mutually unbound and receding, with total energy positive. **Invented** —
/// §2.4 names the outcome and gives the `E > 0` precondition but no pairwise test.
///
/// The pairwise test used is the two-body one: relative kinetic energy in the pair's reduced
/// mass exceeds its mutual potential, and the separation is growing. The `E > 0` check is not
/// redundant with it — three pairwise-unbound instants can occur transiently in a bound
/// system, and the total energy is the invariant that cannot.
pub fn triple_ejection<T: Real>(s: &Cart<T>, m: &[T; 3]) -> bool {
    if energy::energy(&s.r, &s.v, m, T::zero()) <= T::zero() {
        return false;
    }
    PAIRS.iter().all(|&(i, j)| {
        let dr = s.r[j] - s.r[i];
        let dv = s.v[j] - s.v[i];
        let mu = m[i] * m[j] / (m[i] + m[j]);
        let e = T::lit(0.5) * mu * dv.norm_sq() - T::lit(G) * m[i] * m[j] / dr.norm().max(T::DIST_FLOOR);
        e > T::zero() && dr.dot(dv) > T::zero()
    })
}

/// What the integrator saw, in the order it saw it. Times are physical, not fictitious.
#[derive(Clone, Copy, Debug, Default)]
pub struct Events<T> {
    /// First moment a pair fell below `r_coll`, with the pair bitmask **at that moment**.
    /// Sampled inside the RK4 loop, not at sync boundaries: with `n_sync = 32` and
    /// `t_max = 13` the boundaries are 0.4 apart and a close encounter passes between two of
    /// them unseen.
    pub collision: Option<(u8, T)>,
    /// First sync boundary at which the escape test fired, and which body. Sampled at
    /// boundaries because that is where the reference samples it and where the state is
    /// Cartesian.
    pub escape: Option<(u8, T)>,
}

/// BRIEF §2.4's encoding, from the events and the final state.
///
/// Precedence: a non-finite trajectory outranks everything, because its final state carries
/// no information to classify. **Then whichever terminating event happened FIRST, by time.**
/// Then budget exhaustion, then reaching `t_max`.
///
/// # The ordering bug this replaced
///
/// This used to rank collision above escape unconditionally, discarding both times, and
/// justified it as *"collision is sampled continuously, so it is the earliest thing that can
/// fire"*. Continuous sampling makes collision the earliest **detected**, not the earliest
/// **occurring**. Escape is tested only where the state is Cartesian, so an escape that truly
/// happened at `t = 5.0` may not be noticed until `t = 5.28`, and a collision at `t = 5.1` was
/// then reported as the outcome of a trajectory that had already terminated.
///
/// Deciding by `min(t)` removes the dependence on *when each arm happens to be sampled* rather
/// than reducing it. A tie goes to collision: at the same instant it is the more specific
/// event, and the two arms cannot be separated by anything this function can see.
///
/// `t_end` is set the same way, in the driver, so the state and the time it is quoted with
/// cannot disagree.
pub fn classify<T: Real>(
    ev: &Events<T>,
    final_state: &Cart<T>,
    m: &[T; 3],
    finite: bool,
    budget_exhausted: bool,
) -> Outcome {
    if !finite || !final_state.is_finite() {
        return Outcome::new(State::SimFailed, 0);
    }
    // Whichever fired first. `triple_ejection` sits between them as before: it is a *detail*
    // refinement of an escape read off the final state, not a fourth event with a time.
    match (ev.collision, ev.escape) {
        (Some((mask, tc)), Some((_, te))) if tc <= te => {
            return Outcome::new(State::Collision, collision_detail(mask));
        }
        (Some((mask, _)), None) => {
            return Outcome::new(State::Collision, collision_detail(mask));
        }
        _ => {}
    }
    if triple_ejection(final_state, m) {
        return Outcome::new(State::Escape, 3);
    }
    if let Some((b, _)) = ev.escape {
        return Outcome::new(State::Escape, b);
    }
    if let Some((mask, _)) = ev.collision {
        return Outcome::new(State::Collision, collision_detail(mask));
    }
    if budget_exhausted {
        return Outcome::new(State::Running, 0);
    }
    Outcome::new(State::Bounded, 0)
}

/// Index into [`PAIRS`] of the unordered pair `(i, j)`.
///
/// `PAIRS` ordering is load-bearing throughout the crate; this is the one place the mapping
/// from a body pair back to that index is written down.
pub fn pair_index(i: usize, j: usize) -> usize {
    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
    PAIRS
        .iter()
        .position(|&p| p == (lo, hi))
        .expect("every unordered pair of three bodies is in PAIRS")
}

/// Collision mask from the three AZ separations, which are labelled by the *reference triple*
/// rather than by pair index.
///
/// `R1 = r_b - r_a`, `R2 = r_c - r_a`, `R3 = r_c - r_b`. The reference body changes between
/// sync boundaries, so the mapping cannot be hoisted — mixing this up would attribute a
/// collision to the wrong pair while every magnitude stayed correct, which is exactly the
/// class of silent error this project keeps finding.
pub fn collision_pairs_from<T: Real>(
    a: usize,
    b: usize,
    c: usize,
    d1: T,
    d2: T,
    d3: T,
    r_coll: T,
) -> u8 {
    let mut mask = 0u8;
    for (d, (i, j)) in [(d1, (a, b)), (d2, (a, c)), (d3, (b, c))] {
        if d < r_coll {
            mask |= 1 << pair_index(i, j);
        }
    }
    mask
}
