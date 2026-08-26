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

/// The escaping body, or `None`. **Ported** from `tb_all_az.py:59-75`.
///
/// The candidate is the body outside the *tightest* pair; it escapes if its specific orbital
/// energy relative to that pair's barycentre is positive **and** it is receding. A body that
/// is not the candidate can never be labelled escaping — a property of the reference, kept.
pub fn escape_candidate<T: Real>(s: &Cart<T>, m: &[T; 3]) -> Option<u8> {
    let d = newton::pair_dists(&s.r);
    let mut tight = 0usize;
    for k in 1..3 {
        if d[k] < d[tight] {
            tight = k;
        }
    }
    let b = THIRD[tight];
    let o: Vec<usize> = (0..3).filter(|&k| k != b).collect();
    let mb = m[o[0]] + m[o[1]];
    let rc = (s.r[o[0]] * m[o[0]] + s.r[o[1]] * m[o[1]]) / mb;
    let vc = (s.v[o[0]] * m[o[0]] + s.v[o[1]] * m[o[1]]) / mb;
    let dr = s.r[b] - rc;
    let dv = s.v[b] - vc;
    // The reference floors the distance at 1e-12 and writes `mb`, not `G*mb`; G is 1.
    let spec = T::lit(0.5) * dv.norm_sq() - T::lit(G) * mb / dr.norm().max(T::DIST_FLOOR);
    if spec > T::zero() && dr.dot(dv) > T::zero() {
        Some(b as u8)
    } else {
        None
    }
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
