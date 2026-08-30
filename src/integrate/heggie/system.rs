//! The Heggie system: masses and the derived constants, plus the maps in and out.
//!
//! Unlike [`AzSystem`](crate::integrate::az::AzSystem) there is **no reference body**. The three
//! relative vectors sit on an equal footing and the system depends on nothing but the masses, so
//! it is built once per trajectory and never rebuilt. That absence is the whole point of the
//! method and the reason this file has no `choose_reference` analogue.
//!
//! Transcribed from Heggie (1974) §2, Eqs. (4)-(19). Indices here are 0-based; his are 1-based.

use crate::physics::{energy, Cart};
use crate::{Real, Vec2};

use crate::integrate::az::lc;

use super::state::HgState;

/// `(j, k)` for index `i`, i.e. Heggie's `(i+1, i+2)`.
#[inline(always)]
pub fn cyc(i: usize) -> (usize, usize) {
    ((i + 1) % 3, (i + 2) % 3)
}

#[derive(Clone, Copy, Debug)]
pub struct HgSystem<T> {
    pub masses: [T; 3],
    pub mtot: T,
    /// `mu[i]` is Heggie's `mu_{jk} = m_j m_k / (m_j + m_k)`, the reduced mass **of the pair
    /// `q_i` joins**. It pairs with `|P_i|^2`, not with `m_i`. Getting this pairing wrong is
    /// silent, which is why it is a named field rather than an inline expression.
    pub mu: [T; 3],
    /// `mm[i] = m_j m_k`, the product in the potential term `-m_j m_k / |q_i|`.
    pub mm: [T; 3],
    /// `inv_m[i] = 1 / m_i`, the coefficient of the coupling term `p_j . p_k / m_i`.
    ///
    /// Heggie §4 notes that Eq. (21) is inapplicable if any mass vanishes; these reciprocals are
    /// where that shows. They are stored so a caller can inspect them.
    pub inv_m: [T; 3],
    /// Use the numerically stable inverse LC map. Default true.
    pub lc_stable: bool,
}

impl<T: Real> HgSystem<T> {
    pub fn new(masses: [T; 3]) -> Self {
        let mtot = masses[0] + masses[1] + masses[2];
        let mut mu = [T::zero(); 3];
        let mut mm = [T::zero(); 3];
        let mut inv_m = [T::zero(); 3];
        for i in 0..3 {
            let (j, k) = cyc(i);
            mu[i] = masses[j] * masses[k] / (masses[j] + masses[k]);
            mm[i] = masses[j] * masses[k];
            inv_m[i] = T::one() / masses[i];
        }
        Self { masses, mtot, mu, mm, inv_m, lc_stable: true }
    }

    pub fn with_reference_lc(mut self) -> Self {
        self.lc_stable = false;
        self
    }

    /// Heggie's `A_i` (Eq. 18) reduced to the plane is `2 L(Q_i)^T`, so `A_i^T w = 2 L(Q_i) w`.
    ///
    /// Named rather than inlined because every appearance of `A_i` in the paper routes through
    /// here, and `tests/heggie_identities.rs::L0` checks this reduction against Eq. (18) written
    /// out literally as a 4x3 matrix. It is the one step of the transcription the paper does not
    /// state, so it is measured rather than asserted.
    #[inline(always)]
    pub fn a_transpose_apply(u: Vec2<T>, w: Vec2<T>) -> Vec2<T> {
        lc::l_apply(u, w) * T::lit(2.0)
    }

    /// `W_i = L(Q_i) P_i`, the quantity every coupling term in `Gamma*` is built from.
    ///
    /// `P_i^T A_i A_j^T P_j = 4 W_i . W_j`, which is why the coupling reads as a plain dot
    /// product of two vectors rather than a quadratic form.
    #[inline(always)]
    pub fn w(s: &HgState<T>, i: usize) -> Vec2<T> {
        lc::l_apply(s.u[i], s.p[i])
    }

    /// Cartesian -> the enlarged relative variables, Heggie Eq. (4) and Eq. (8b).
    ///
    /// `q_i(0) = r_j - r_k` and, **in the centre-of-mass frame**, `p_i(0) = (m_j v_j - m_k v_k)/3`.
    /// Eq. (8b)'s general form carries `-((m_j - m_k)/M) sum p'` which vanishes there.
    ///
    /// The COM frame is assumed and **asserted**, not hoped for: a construction that assumes a
    /// COM-centred input returns a drifting system without one, and this project has that on
    /// record from `momenta_for`. `com_defect` reports it.
    pub fn enlarged_from_cart(&self, s: &Cart<T>) -> ([Vec2<T>; 3], [Vec2<T>; 3]) {
        let three = T::lit(3.0);
        let mut q = [Vec2::zero(); 3];
        let mut p = [Vec2::zero(); 3];
        for i in 0..3 {
            let (j, k) = cyc(i);
            q[i] = s.r[j] - s.r[k];
            p[i] = (s.v[j] * self.masses[j] - s.v[k] * self.masses[k]) / three;
        }
        (q, p)
    }

    /// `|sum m_i r_i| + |sum m_i v_i|`, both relative to their own scale.
    ///
    /// Two terms summed but the two moments are asserted **separately** by the tests: zero
    /// momentum does not imply a zero first moment, and this project has a mass audit on record
    /// that only holds because the two are checked apart.
    pub fn com_defect(&self, s: &Cart<T>) -> (T, T) {
        let mut mr = Vec2::zero();
        let mut mv = Vec2::zero();
        let mut sr = T::zero();
        let mut sv = T::zero();
        for i in 0..3 {
            mr += s.r[i] * self.masses[i];
            mv += s.v[i] * self.masses[i];
            sr = sr.max(s.r[i].norm() * self.masses[i]);
            sv = sv.max(s.v[i].norm() * self.masses[i]);
        }
        (mr.norm() / sr.max(T::TINY), mv.norm() / sv.max(T::TINY))
    }

    /// Recover the physical Cartesian state, Heggie Eqs. (10) and (12), in the COM frame.
    ///
    /// `q_i* = (m_j q_k - m_k q_j)/M` and `p_i* = -p_j + p_k`. The first is the crossed-mass
    /// shape this project has recorded as invisible to a `sum p_i = 0` check; it is anchored by
    /// the round-trip test at unequal masses instead.
    pub fn cart_from_enlarged(&self, q: &[Vec2<T>; 3], p: &[Vec2<T>; 3]) -> Cart<T> {
        let mut r = [Vec2::zero(); 3];
        let mut v = [Vec2::zero(); 3];
        for i in 0..3 {
            let (j, k) = cyc(i);
            r[i] = (q[k] * self.masses[j] - q[j] * self.masses[k]) / self.mtot;
            v[i] = (p[k] - p[j]) * self.inv_m[i];
        }
        Cart::new(r, v)
    }

    /// Cartesian -> regularised. Returns the state and the frozen energy `h`.
    ///
    /// `Q_i = u_of_rho(q_i)` and `P_i = 2 L(Q_i)^T p_i`, the latter being Eq. (19) inverted using
    /// `L(u)^T L(u) = R I`.
    pub fn to_reg(&self, s: &Cart<T>) -> (HgState<T>, T) {
        let (q, p) = self.enlarged_from_cart(s);
        let inv = if self.lc_stable { lc::u_of_rho } else { lc::u_of_rho_reference };
        let two = T::lit(2.0);
        let mut st = HgState { u: [Vec2::zero(); 3], p: [Vec2::zero(); 3], t: T::zero() };
        for i in 0..3 {
            st.u[i] = inv(q[i]);
            st.p[i] = lc::lt_apply(st.u[i], p[i]) * two;
        }
        (st, self.energy_enlarged(&q, &p))
    }

    /// Regularised -> the enlarged relative variables. `q_i = rho_of_u(Q_i)` is Heggie's
    /// identity `q_i = (1/2) A_i^T Q_i`; `p_i = L(Q_i) P_i / (2 R_i)` is Eq. (19).
    ///
    /// `R_i` **is** floored here, and is not in `HgState::r`. Same asymmetry as AZ's, and for
    /// the same reason: the floor belongs where the division happens.
    pub fn phys_from_state(&self, s: &HgState<T>) -> ([Vec2<T>; 3], [Vec2<T>; 3]) {
        let two = T::lit(2.0);
        let mut q = [Vec2::zero(); 3];
        let mut p = [Vec2::zero(); 3];
        for i in 0..3 {
            let r = s.r(i).max(T::TINY);
            q[i] = lc::rho_of_u(s.u[i]);
            p[i] = lc::l_apply(s.u[i], s.p[i]) / (two * r);
        }
        (q, p)
    }

    /// Heggie's Eq. (6) in the centre-of-mass frame, where the `|<p>|^2 / 2M` term drops.
    ///
    /// Three pieces, and the pairing of each is the part that fails silently:
    ///   - `|p_i|^2 / (2 mu_i)`, where `mu_i` is the reduced mass of the pair `q_i` **joins**;
    ///   - `- p_j . p_k / m_i`, the coupling, linear in each momentum;
    ///   - `- m_j m_k / |q_i|`.
    pub fn energy_enlarged(&self, q: &[Vec2<T>; 3], p: &[Vec2<T>; 3]) -> T {
        let two = T::lit(2.0);
        let mut e = T::zero();
        for i in 0..3 {
            let (j, k) = cyc(i);
            e += p[i].norm_sq() / (two * self.mu[i]);
            e -= p[j].dot(p[k]) * self.inv_m[i];
            e -= self.mm[i] / q[i].norm().max(T::TINY);
        }
        e
    }

    /// Energy directly from a regularised state.
    pub fn energy_of(&self, s: &HgState<T>) -> T {
        let (q, p) = self.phys_from_state(s);
        self.energy_enlarged(&q, &p)
    }

    /// Heggie's Eq. (1), the ordinary Cartesian energy. Delegates to [`energy::energy`] at
    /// `eps2 = 0` rather than restating it — the anchor chain is only worth anything if its far
    /// end is the same function the rest of the project validates against.
    pub fn energy_cartesian(&self, s: &Cart<T>) -> T {
        energy::energy(&s.r, &s.v, &self.masses, T::zero())
    }

    /// Regularised -> Cartesian.
    pub fn to_cartesian(&self, s: &HgState<T>) -> Cart<T> {
        let (q, p) = self.phys_from_state(s);
        self.cart_from_enlarged(&q, &p)
    }

    /// `|sum q_i| / max |q_i|`, Heggie's Eq. (9) as a running residual.
    ///
    /// An integral of the enlarged motion, so it is zero along the exact flow and drifts under
    /// RK4. It exists only in this method — AZ has no analogue — and it is free.
    ///
    /// Normalised by the **largest** `|q_i|` and never by the smallest or by `R1 R2 R3`: a
    /// denominator that vanishes at a close approach turns a residual into a measure of how
    /// closely the sampling approached collision, which is the defect already on record for the
    /// `Gamma` residual.
    pub fn sum_q_residual(&self, s: &HgState<T>) -> T {
        let (q, _) = self.phys_from_state(s);
        let mut sum = Vec2::zero();
        let mut scale = T::zero();
        for v in q {
            sum += v;
            scale = scale.max(v.norm());
        }
        if scale <= T::zero() || !scale.is_finite() {
            return T::infinity();
        }
        sum.norm() / scale
    }
}
