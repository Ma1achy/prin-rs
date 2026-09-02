//! Heggie's global regularisation of the three-body problem.
//!
//! D. C. Heggie, *A global regularisation of the gravitational N-body problem*, Celestial
//! Mechanics **10** (1974) 217-241. Equation numbers throughout this module are his.
//!
//! # Why this exists
//!
//! AZ must choose a reference body — the one not in the longest side — and re-choose it at every
//! sync boundary. Measured on `config_stability` at 256^2, holding the step size fixed by scaling
//! `eta` with `n_sync`:
//!
//! ```text
//!   LC branch unconditioned            2.5e-6 decades
//!   hysteresis, switches 17.8 -> 6.7   7.5e-5 decades
//!   re-registration x2, fixed step     4.4e-1 decades    <- 6000x
//! ```
//!
//! It is not *which* chart is chosen, it is *how often the state is passed through one*. Every
//! remedy tried against AZ smooths the choice, and the measurement says the choice is the small
//! term. Heggie has **no reference body and therefore no re-registration at all**: three relative
//! vectors on an equal footing, one transformation each, one symmetric time transformation.
//!
//! Heggie's own §3 verdict is that his method is "significantly the weaker of the two" against AZ
//! on his 'realistic' problem and less efficient on the Pythagorean one, at ~1.6x the cost per
//! step. That is **accuracy per step on one trajectory**, which is the wrong axis for a defect
//! that is a discontinuity across *neighbouring* initial conditions. The remark that bears on
//! this project is his other one: the close-triple-encounter time reversals "[do] not depend on
//! any judicious choice for the initial labelling of the bodies, which is the case with the
//! method of Aarseth and Zare."
//!
//! # The planar reduction
//!
//! Heggie works in three dimensions and applies KS per vector, giving 4-vectors `Q_i` and the 4x3
//! matrix `A_i` of Eq. (18). This project is planar throughout (BRIEF §2.3). Setting
//! `Q_3 = Q_4 = 0` in Eq. (17) leaves `q_i = (Q_1^2 - Q_2^2, 2 Q_1 Q_2)` — the Levi-Civita map,
//! already implemented as [`lc::rho_of_u`](crate::integrate::az::lc::rho_of_u) — and Eq. (18)
//! reduces to the 2x2 block
//!
//! ```text
//! A_i = 2 L(Q_i)^T ,   L(u) = [[u.x, -u.y], [u.y, u.x]]
//! ```
//!
//! so `R_i = Q_i^T Q_i = |Q_i|^2`, Heggie's identity `q_i = (1/2) A_i^T Q_i` **is** `rho_of_u`,
//! and the coupling `P_j^T A_j A_k^T P_k` becomes `4 (L(Q_j)P_j) . (L(Q_k)P_k)` — the same shape
//! as AZ's cross term, through the same `l_apply`. The state is thirteen numbers, not
//! twenty-five.
//!
//! **This reduction is the one step of the transcription the paper does not state**, so it is
//! verified numerically against Eq. (18) written out literally as a 4x3 matrix, in
//! `tests/heggie_identities.rs`, rather than asserted here.
//!
//! # Index convention
//!
//! Heggie is 1-based; this module is 0-based. His `q_1 = q_2' - q_3'` is `q[0] = r[1] - r[2]`.
//! Every cyclic formula uses `(j, k) = cyc(i)`, which is his `(i+1, i+2)`.

pub mod driver;
pub mod hamiltonian;
pub mod rk4;
pub mod state;
pub mod system;

pub use driver::{integrate_hg, HgDtauMode, HgOpts, HgOut};
pub use hamiltonian::HgTime;
pub use state::HgState;
pub use system::HgSystem;
