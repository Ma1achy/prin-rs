//! Chazy final-motion classification, **as of time `t`**.
//!
//! Chazy (1922) classifies the planar three-body problem by its motion as `t -> infinity`:
//! H (all three separate hyperbolically), HP, P, HE (a bound pair and a hyperbolic escaper), PE,
//! B (bounded) and OS (oscillatory: `limsup` of the largest separation infinite, `liminf` finite).
//! A finite march can only ever say what the motion *looks like at the playhead*, so every label
//! here is an approximation that holds **up to the current time** and is recomputed at every
//! snapshot. Watching it evolve is the point.
//!
//! Everything [`classify`] reads is either the live state or a running latch ([`Latch`]): the
//! largest moment of inertia seen so far and the previous snapshot's shape vector. So it is
//! live-compatible in the CLAUDE.md sense -- no re-run from `t = 0`, nothing after `t`.
//!
//! Two physical facts are used as certificates rather than heuristics:
//!
//! * **Lagrange–Jacobi.** `I'' = 4E - 2V > 4E`. With `E > 0` the moment of inertia grows without
//!   bound, so B and OS are impossible and an interacting `E > 0` system is [`Motion::Transient`],
//!   never "bounded so far". With `E < 0`, H, HP and P are impossible (all three cannot separate).
//! * **Parabolic is a band, not a class.** PE, HP and P need an exact zero of a two-body energy and
//!   have zero measure. A kinetic-to-potential ratio within `eps` of 1 is reported as the parabolic
//!   label and scored with its hyperbolic neighbour (see [`Motion::coarse`]).
//!
//! [`fewbody_resolved`] is the literature's standard for "this scattering is over", Fregeau et al.
//! (2004, MNRAS 352, 1), §3.3.5, specialised to three bodies. It is here as the **reference** the
//! proposal is measured against (`examples/chazy_probe.rs`, `results/chazy/README.md`), and its
//! tidal isolation test is one of the two gates [`classify`] can run.

use crate::physics::{Cart, PAIRS};
use crate::Vec2;

/// The motion class at the playhead.
#[derive(Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Motion {
    /// All three pairs receding with every two-body ratio hyperbolic.
    H = 0,
    /// All three receding, one or two pairs inside the parabolic band.
    HP = 1,
    /// All three receding, every pair inside the parabolic band.
    P = 2,
    /// A bound pair and a receding hyperbolic escaper.
    HE = 3,
    /// A bound pair and a receding escaper inside the parabolic band.
    PE = 4,
    /// Bounded and Mardling–Aarseth (2001) stable as a hierarchical triple.
    BStable = 5,
    /// `E < 0`, interacting: bounded so far.
    BSoFar = 6,
    /// `E > 0`, interacting: Lagrange–Jacobi says it cannot stay bounded.
    Transient = 7,
    /// `E < 0`, interacting, and has been far out (`I_max / I > os_ratio`) and come back.
    OsCandidate = 8,
    /// One body far and receding but still bound to its pair: it will fall back.
    Excursion = 9,
    /// Only the reference classifier uses this: not resolved by its test.
    Unresolved = 10,
}

impl Motion {
    pub const ALL: [Motion; 11] = [
        Motion::H, Motion::HP, Motion::P, Motion::HE, Motion::PE, Motion::BStable,
        Motion::BSoFar, Motion::Transient, Motion::OsCandidate, Motion::Excursion, Motion::Unresolved,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Motion::H => "H", Motion::HP => "HP", Motion::P => "P", Motion::HE => "HE",
            Motion::PE => "PE", Motion::BStable => "B-stable", Motion::BSoFar => "B-so-far",
            Motion::Transient => "transient", Motion::OsCandidate => "OS-candidate",
            Motion::Excursion => "excursion", Motion::Unresolved => "unresolved",
        }
    }

    /// A final-motion label, as opposed to "not decided yet".
    pub fn committed(self) -> bool {
        matches!(self, Motion::H | Motion::HP | Motion::P | Motion::HE | Motion::PE | Motion::BStable)
    }

    /// Parabolic labels scored with their hyperbolic neighbour: the band has width, the class
    /// has none, so agreement on the band is not a thing a finite run can be held to.
    pub fn coarse(self) -> Motion {
        match self {
            Motion::HP | Motion::P => Motion::H,
            Motion::PE => Motion::HE,
            c => c,
        }
    }
}

/// Which test says a runaway is asymptotic enough to label.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Gate {
    /// Hierarchy ratio (escaper distance over its pair's separation) above `h_far`, and the shape
    /// settled: the escape rule's closure chord over one snapshot window below `tau`.
    Hierarchy,
    /// Fewbody's isolation: the tidal perturbation on the pair from the escaper below `delta`;
    /// three singles need only recede with `E > 0`.
    Tidal,
}

#[derive(Clone, Copy, Debug)]
pub struct ChazyCfg {
    pub gate: Gate,
    /// Half-width of the parabolic band on the kinetic/potential ratio.
    pub eps: f64,
    /// Hierarchy-ratio gate, and the three-singles separation gate in units of `R0`.
    pub h_far: f64,
    /// Shape-settling chord over one snapshot window. `INFINITY` switches settling off.
    pub tau: f64,
    /// `I_max / I` above which a bounded-so-far system is an OS candidate.
    pub os_ratio: f64,
    /// Tidal tolerance for [`Gate::Tidal`].
    pub delta: f64,
    /// Label a Mardling–Aarseth-stable bound hierarchy [`Motion::BStable`].
    pub mardling_aarseth: bool,
}

impl Default for ChazyCfg {
    /// The measured choice: see `results/chazy/README.md`.
    fn default() -> Self {
        Self {
            gate: Gate::Hierarchy,
            eps: 0.05,
            h_far: 5.0,
            tau: f64::INFINITY,
            os_ratio: 25.0,
            delta: 1e-5,
            mardling_aarseth: true,
        }
    }
}

/// The two running quantities [`classify`] needs beyond the live state. Both are O(1) and update
/// once per snapshot.
#[derive(Clone, Copy, Debug)]
pub struct Latch {
    pub i_max: f64,
    pub n_prev: Option<[f64; 3]>,
}

impl Latch {
    pub fn new(i0: f64) -> Self {
        Self { i_max: i0, n_prev: None }
    }
}

/// `(i, j)`, the pair not containing `b`.
pub fn other(b: usize) -> (usize, usize) {
    match b {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

fn pair_index(i: usize, j: usize) -> usize {
    PAIRS.iter().position(|&(a, b)| (a, b) == (i.min(j), i.max(j))).unwrap()
}

/// Relative two-body `(kinetic, |potential|, separation, radial velocity)` of pair `(i, j)`.
pub fn pair_kup(s: &Cart<f64>, m: &[f64; 3], i: usize, j: usize) -> (f64, f64, f64, f64) {
    let mu = m[i] * m[j] / (m[i] + m[j]);
    let dr = s.r[j] - s.r[i];
    let dv = s.v[j] - s.v[i];
    let d = dr.norm();
    (0.5 * mu * dv.norm_sq(), m[i] * m[j] / d, d, dr.dot(dv) / d)
}

/// Body `b` against the centre of mass of the other two, same tuple as [`pair_kup`].
pub fn out_kup(s: &Cart<f64>, m: &[f64; 3], b: usize) -> (f64, f64, f64, f64) {
    let (i, j) = other(b);
    let mp = m[i] + m[j];
    let rc = (s.r[i] * m[i] + s.r[j] * m[j]) / mp;
    let vc = (s.v[i] * m[i] + s.v[j] * m[j]) / mp;
    let dr = s.r[b] - rc;
    let dv = s.v[b] - vc;
    let mu = m[b] * mp / (m[b] + mp);
    let d = dr.norm();
    (0.5 * mu * dv.norm_sq(), m[b] * mp / d, d, dr.dot(dv) / d)
}

/// Kepler `(a, e, h)` of a relative orbit with total mass `mtot` (`G = 1`); `a < 0` if unbound.
pub fn kepler(dr: Vec2<f64>, dv: Vec2<f64>, mtot: f64) -> (f64, f64, f64) {
    let r = dr.norm();
    let spec = 0.5 * dv.norm_sq() - mtot / r;
    let h = dr.x * dv.y - dr.y * dv.x;
    let a = -mtot / (2.0 * spec);
    let e = (1.0 + 2.0 * spec * h * h / (mtot * mtot)).max(0.0).sqrt();
    (a, e, h)
}

/// Mardling & Aarseth (2001) for a planar hierarchy, inner pair `(i, j)`, outer body `k`:
/// `R_p^out / a_in > 2.8 [(1 + q_out)(1 + e_out)/sqrt(1 - e_out)]^{2/5} (1 - 0.3 i/180 deg)`,
/// with `i` 0 (prograde) or 180 (retrograde). `None` if either level is unbound.
pub fn ma_stable(s: &Cart<f64>, m: &[f64; 3], i: usize, j: usize, k: usize) -> Option<bool> {
    let mb = m[i] + m[j];
    let (ai, _, hi) = kepler(s.r[j] - s.r[i], s.v[j] - s.v[i], mb);
    let rc = (s.r[i] * m[i] + s.r[j] * m[j]) / mb;
    let vc = (s.v[i] * m[i] + s.v[j] * m[j]) / mb;
    let (ao, eo, ho) = kepler(s.r[k] - rc, s.v[k] - vc, mb + m[k]);
    if !(ai > 0.0 && ao > 0.0 && eo < 1.0) {
        return None;
    }
    let q = m[k] / mb;
    let incl = if hi * ho >= 0.0 { 1.0 } else { 0.7 };
    let crit = 2.8 * ((1.0 + q) * (1.0 + eo) / (1.0 - eo).sqrt()).powf(0.4) * incl;
    Some(ao * (1.0 - eo) / ai > crit)
}

fn band(q: f64, eps: f64) -> i8 {
    if q > 1.0 + eps {
        1
    } else if q < 1.0 - eps {
        -1
    } else {
        0
    }
}

/// The motion class at the playhead. `e` is the total energy, `r0` the initial hyperradius,
/// `i_now` the current moment of inertia and `n_now` the current shape vector; `latch` holds the
/// previous snapshot's shape and `I_max` **before** this snapshot updates them.
pub fn classify(
    s: &Cart<f64>,
    m: &[f64; 3],
    e: f64,
    r0: f64,
    i_now: f64,
    n_now: &[f64; 3],
    latch: &Latch,
    cfg: &ChazyCfg,
) -> Motion {
    let settled = cfg.tau.is_infinite()
        || latch.n_prev.is_some_and(|np| {
            (0..3).map(|k| (n_now[k] - np[k]).powi(2)).sum::<f64>().sqrt() < cfg.tau
        });
    let pairs: [(f64, f64, f64, f64); 3] = [0, 1, 2].map(|k| {
        let (i, j) = PAIRS[k];
        pair_kup(s, m, i, j)
    });

    // All three apart.
    let dmin = pairs.iter().map(|x| x.2).fold(f64::INFINITY, f64::min);
    let apart = match cfg.gate {
        Gate::Hierarchy => dmin > cfg.h_far * r0 && settled,
        Gate::Tidal => e > 0.0,
    };
    if apart && pairs.iter().all(|x| x.3 > 0.0) {
        let bands = pairs.map(|x| band(x.0 / x.1, cfg.eps));
        if bands.iter().all(|&b| b >= 0) {
            return match bands.iter().filter(|&&b| b == 0).count() {
                0 => Motion::H,
                1 | 2 => Motion::HP,
                _ => Motion::P,
            };
        }
    }

    // Hierarchical: the body most separated from its pair.
    let (b, ratio) = (0..3)
        .map(|b| {
            let (i, j) = other(b);
            (b, out_kup(s, m, b).2 / pairs[pair_index(i, j)].2)
        })
        .fold((0, 0.0), |acc, x| if x.1 > acc.1 { x } else { acc });
    let (i, j) = other(b);
    let k = pair_index(i, j);
    let (ko, uo, dk, vr) = out_kup(s, m, b);
    // The pair's binding is a SIGN: a bound eccentric pair at pericentre has K/U = (1+e)/2,
    // which a parabolic band around 1 would misread as marginal.
    let qp = if pairs[k].0 < pairs[k].1 { -1 } else { band(pairs[k].0 / pairs[k].1, cfg.eps) };
    let gate = match cfg.gate {
        Gate::Hierarchy => ratio > cfg.h_far && settled,
        Gate::Tidal => {
            let (ai, ei, _) = kepler(s.r[j] - s.r[i], s.v[j] - s.v[i], m[i] + m[j]);
            let r_pair = if ai > 0.0 { ai * (1.0 + ei) } else { pairs[k].2 };
            2.0 * m[b] / (m[i] + m[j]) * (r_pair / dk).powi(3) < cfg.delta
        }
    };
    if gate && vr > 0.0 {
        match (qp, band(ko / uo, cfg.eps)) {
            (-1, 1) => return Motion::HE,
            (-1, 0) => return Motion::PE,
            (-1, _) => return Motion::Excursion,
            // Escaper gone with its former partners unbound to each other: a two-stage break-up.
            (_, qb) if qb >= 0 => return Motion::H,
            _ => {}
        }
    }

    if e > 0.0 {
        return Motion::Transient;
    }
    if cfg.mardling_aarseth && qp == -1 && ko < uo && ma_stable(s, m, i, j, b) == Some(true) {
        return Motion::BStable;
    }
    if latch.i_max / i_now > cfg.os_ratio {
        Motion::OsCandidate
    } else {
        Motion::BSoFar
    }
}

/// Fregeau et al. (2004) §3.3.5's termination test for three bodies, the reference.
///
/// Pairs the bound top-level nodes of smallest semimajor axis until none is bound to another,
/// then returns a resolved class iff every pair of top-level nodes recedes, the tidal
/// perturbation on every binary is below `delta`, every tree is Mardling–Aarseth stable, and the
/// top-level system has positive energy; otherwise [`Motion::Unresolved`].
pub fn fewbody_resolved(s: &Cart<f64>, m: &[f64; 3], e: f64, delta: f64) -> Motion {
    let mut best: Option<(usize, usize, f64)> = None;
    for &(i, j) in PAIRS.iter() {
        let (a, _, _) = kepler(s.r[j] - s.r[i], s.v[j] - s.v[i], m[i] + m[j]);
        if a > 0.0 && best.map_or(true, |b| a < b.2) {
            best = Some((i, j, a));
        }
    }
    match best {
        None => {
            let rec = PAIRS.iter().all(|&(i, j)| pair_kup(s, m, i, j).3 > 0.0);
            if rec && e > 0.0 { Motion::H } else { Motion::Unresolved }
        }
        Some((i, j, ai)) => {
            let k = 3 - i - j;
            let (ko, uo, dk, vr) = out_kup(s, m, k);
            if ko < uo {
                return match ma_stable(s, m, i, j, k) {
                    Some(true) => Motion::BStable,
                    _ => Motion::Unresolved,
                };
            }
            let (_, ei, _) = kepler(s.r[j] - s.r[i], s.v[j] - s.v[i], m[i] + m[j]);
            let tid = 2.0 * m[k] / (m[i] + m[j]) * (ai * (1.0 + ei) / dk).powi(3);
            if vr > 0.0 && tid < delta { Motion::HE } else { Motion::Unresolved }
        }
    }
}
