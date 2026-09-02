//! `Gamma*` and its derivatives, kept **adjacent on purpose**, as in the AZ module.
//!
//! Heggie's Eq. (20) sets `dtau = dt / (R1 R2 R3)`, and `Gamma* = (H - h) R1 R2 R3`. Written in
//! the planar variables with `W_i = L(Q_i) P_i`, `R_i = |Q_i|^2`:
//!
//! ```text
//! Gamma* = sum_i (1/8) (R_j R_k / mu_i) |P_i|^2
//!        - sum_i (1/4) (R_i / m_i) W_j . W_k
//!        - sum_i m_j m_k R_j R_k
//!        - h R1 R2 R3
//! ```
//!
//! **There is no `1/r` term anywhere.** A polynomial of degree six in the coordinates and two in
//! the momenta, and that absence *is* the globality: AZ's `Gamma` still carries
//! `-A B m_b m_c / |R3|`, the unregularised third side, which is exactly why AZ has to choose a
//! reference body so that side is the longest.
//!
//! **Everything here is written cyclically over `i`, never as three expanded branches.** AZ's
//! `g1`/`g2` are hand-written and the documented hazard is that `g2` differs from `g1` in two
//! easily-missed places — a cross mass pair and a sign. A cyclic form cannot carry a sign error
//! in one branch and not another, which removes that bug class by construction rather than by
//! vigilance. `tests/heggie_hamiltonian_fd.rs` still proves the test would catch one.
//!
//! As in AZ, an FD test alone is **not sufficient**: a sign error present in both `gamma` and
//! `deriv` passes it silently. `gamma` is independently anchored by
//! `Gamma* == (energy_enlarged - h) * R1 R2 R3`, and `energy_enlarged` in turn by the Cartesian
//! energy through Heggie's Eqs. (10) and (12). `tests/heggie_identities.rs` holds that chain.

use crate::{Real, Vec2};

use crate::integrate::az::lc;

use super::state::HgState;
use super::system::{cyc, HgSystem};

/// Which of Heggie's two time transformations is in force.
///
/// The choice is not cosmetic. His §3 reports that near a triple collision Eq. (20)/(21) **and**
/// the set (22),(23),(25) both carry modes growing as `R_i^{-1}`, and that **only (22)-(24) has
/// no unstable mode** — he likens the retained term to Baumgarte's control term. So the default
/// keeps it. `Product` and `keep_gamma_term: false` are named and retained because they are the
/// arms that claim is about.
///
/// A second argument he does not make: under `r -> alpha r`, `t -> alpha^{3/2} t`, Eq. (22) gives
/// `dtau ~ alpha^0`, so `tau` is **scale-invariant** and a fixed `dtau` introduces no length or
/// time scale — which is precisely what BRIEF §2.5 requires. Eq. (20) gives `dtau ~ alpha^{-3/2}`
/// and does not have this property.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HgTime {
    /// Eq. (20): `dtau = dt / (R1 R2 R3)`, Hamilton's equations directly (Eq. 21).
    Product,
    /// Eq. (22) at `n = 3/2`: `dtau = S^{3/2} / (R1 R2 R3) dt`, with Eq. (23) and Eq. (24).
    ///
    /// `keep_gamma_term: false` is Eq. (25), legitimate on the solution path because `Gamma* = 0`
    /// there, and **numerically different**: it is the arm Heggie found unstable.
    SumPow32 { keep_gamma_term: bool },
}

impl Default for HgTime {
    fn default() -> Self {
        HgTime::SumPow32 { keep_gamma_term: true }
    }
}

/// The ten terms of `Gamma*`, in a fixed order: three kinetic, three coupling, three potential,
/// one energy. Returned rather than summed so the residual can be scaled by the largest.
#[inline]
fn terms<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, h: T) -> [T; 10] {
    let eight = T::lit(8.0);
    let four = T::lit(4.0);
    let r = [s.r(0), s.r(1), s.r(2)];
    let w = [HgSystem::w(s, 0), HgSystem::w(s, 1), HgSystem::w(s, 2)];
    let mut out = [T::zero(); 10];
    for i in 0..3 {
        let (j, k) = cyc(i);
        out[i] = r[j] * r[k] * s.p[i].norm_sq() / (eight * sys.mu[i]);
        out[3 + i] = -(r[i] * sys.inv_m[i] * w[j].dot(w[k]) / four);
        out[6 + i] = -(sys.mm[i] * r[j] * r[k]);
    }
    out[9] = -(h * r[0] * r[1] * r[2]);
    out
}

/// `Gamma* = (H - h) R1 R2 R3`, Heggie's §2 final form.
pub fn gamma<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, h: T) -> T {
    terms(sys, s, h).into_iter().fold(T::zero(), |a, b| a + b)
}

/// Hamilton's equations under Eq. (20)/(21): `u[i] = dGamma*/dP_i`, `p[i] = -dGamma*/dQ_i`,
/// `t = dt/dtau = R1 R2 R3`.
///
/// The two bilinear slots are the place to look if this is ever wrong. `W_i = L(Q_i) P_i` is
/// bilinear and `L(u)w = L(w)u`, so differentiating `W_i . X` gives `L(Q_i)^T X` with respect to
/// `P_i` and `L(P_i)^T X` with respect to `Q_i` — the argument moves into the **matrix** slot.
/// It looks like a typo. It is the same identity `lc.rs` documents for AZ's cross term.
pub fn deriv<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, h: T) -> HgState<T> {
    let four = T::lit(4.0);
    let two = T::lit(2.0);
    let r = [s.r(0), s.r(1), s.r(2)];
    let w = [HgSystem::w(s, 0), HgSystem::w(s, 1), HgSystem::w(s, 2)];

    let mut d = HgState { u: [Vec2::zero(); 3], p: [Vec2::zero(); 3], t: r[0] * r[1] * r[2] };
    for i in 0..3 {
        let (j, k) = cyc(i);

        // dGamma*/dP_i: own kinetic term, plus the two coupling terms W_i appears in.
        d.u[i] = s.p[i] * (r[j] * r[k] / (four * sys.mu[i]))
            - (lc::lt_apply(s.u[i], w[k]) * (r[j] * sys.inv_m[j])
                + lc::lt_apply(s.u[i], w[j]) * (r[k] * sys.inv_m[k]))
                / four;

        // dGamma*/dQ_i. Five groups, every one of them carrying R_i or W_i:
        //   1. the two OTHER kinetic terms, which contain R_i
        //   2. own coupling term, through R_i
        //   3. the two other coupling terms, through W_i — argument in the matrix slot
        //   4. the two potential terms containing R_i
        //   5. the energy term
        let g = s.u[i]
            * ((r[k] * s.p[j].norm_sq() / sys.mu[j] + r[j] * s.p[k].norm_sq() / sys.mu[k]) / four)
            - (s.u[i] * (two * sys.inv_m[i] * w[j].dot(w[k]))
                + lc::lt_apply(s.p[i], w[k]) * (r[j] * sys.inv_m[j])
                + lc::lt_apply(s.p[i], w[j]) * (r[k] * sys.inv_m[k]))
                / four
            - s.u[i] * (two * (sys.mm[j] * r[k] + sys.mm[k] * r[j]))
            - s.u[i] * (two * h * r[j] * r[k]);

        d.p[i] = -g;
    }
    d
}

/// The equations of motion under a chosen time transformation.
///
/// Under [`HgTime::Product`] this is exactly [`deriv`]. Under [`HgTime::SumPow32`] it is
/// Heggie's Eqs. (23) and (24): everything scales by `S^{-3/2}`, and the momentum equation picks
/// up `+ (3/2) Gamma* d ln S / dQ_i`, which is `3 Gamma* Q_i / S`.
///
/// `Gamma*` is evaluated, not assumed zero. It **is** zero on the exact solution path, which is
/// what makes Eq. (25) formally legitimate — and its numerical value away from that path is
/// precisely the stabilising control term.
pub fn deriv_time<T: Real>(
    sys: &HgSystem<T>,
    s: &HgState<T>,
    h: T,
    time: HgTime,
) -> HgState<T> {
    let mut d = deriv(sys, s, h);
    let keep = match time {
        HgTime::Product => return d,
        HgTime::SumPow32 { keep_gamma_term } => keep_gamma_term,
    };
    let s_sum = s.s().max(T::TINY);
    let f = T::one() / (s_sum * s_sum.sqrt());
    let g = if keep { gamma(sys, s, h) } else { T::zero() };
    let three = T::lit(3.0);
    for i in 0..3 {
        d.u[i] = d.u[i] * f;
        // d.p[i] is already -dGamma*/dQ_i, so the control term enters with a leading minus in
        // Eq. (24) applied to a bracket that itself carries a minus: -(dG/dQ - (3/2) G dlnS/dQ).
        d.p[i] = (d.p[i] + s.u[i] * (three * g / s_sum)) * f;
    }
    d.t = d.t * f;
    d
}

/// `Gamma*` together with the magnitude of its largest term.
pub fn gamma_scaled<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, h: T) -> (T, T) {
    let t = terms(sys, s, h);
    let sum = t.into_iter().fold(T::zero(), |a, b| a + b);
    let scale = t.into_iter().fold(T::zero(), |a, b| a.max(b.abs()));
    (sum, scale)
}

/// A free integration-quality residual: `|Gamma*|` against the magnitude of its largest term.
///
/// **Not divided by `R1 R2 R3`.** That is the obvious thing and it is wrong for the same reason
/// dividing AZ's by `A*B` is: the factor vanishes at a close approach, so the quotient blows up
/// exactly where the integrator is working hardest and the "residual" then measures how closely
/// the sampling approached collision rather than how well energy was conserved.
pub fn gamma_residual<T: Real>(sys: &HgSystem<T>, s: &HgState<T>, h: T) -> T {
    let (g, scale) = gamma_scaled(sys, s, h);
    if scale <= T::zero() || !scale.is_finite() {
        return T::infinity();
    }
    (g / scale).abs()
}
