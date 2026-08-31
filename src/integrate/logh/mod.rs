//! The logarithmic Hamiltonian — algorithmic regularisation, with no coordinate transformation.
//!
//! Mikkola & Tanikawa, MNRAS **310** (1999) 745, and independently Preto & Tremaine, AJ **118**
//! (1999) 2532. Implementation form after Mikkola & Merritt, AJ **135** (2008) 2398.
//!
//! # Why this exists: it is a falsification test, not a third opinion
//!
//! The standing claim, measured on `config_stability` at 256^2 with the step size held fixed by
//! scaling `eta` with `n_sync`:
//!
//! ```text
//!   LC branch unconditioned            2.5e-6 decades
//!   hysteresis, switches 17.8 -> 6.7   7.5e-5 decades
//!   re-registration x2, fixed step     4.4e-1 decades    <- 6000x
//! ```
//!
//! Heggie wins 31 of 32 gallery cases and has no reference body, which is *consistent* with that
//! mechanism and does not establish it — Heggie changes the chart **and** removes the
//! re-registration in one move.
//!
//! logH separates them. It has **no coordinate transformation at all**: a time transformation and
//! a good integrator, nothing else. That is a strictly stronger form of the property the Heggie
//! win is attributed to.
//!
//! - If the mechanism is real, logH matches or beats Heggie.
//! - **If logH loses to Heggie, the mechanism is wrong** and the win comes from something else.
//!   `INVESTIGATION.md` §5 names the next candidates: the KS square-root's own round-off, the
//!   per-boundary energy re-freeze, and the landing residual at each boundary.
//!
//! Both outcomes are informative, and the second is the one worth guarding against wishing away.
//!
//! # The method
//!
//! ```text
//!   U = +sum_pairs G m_i m_j / |r_i - r_j| > 0        K = 0.5 sum_i m_i |v_i|^2
//!   B = U - K = -E, frozen at t = 0 and CONSTANT
//!
//!   Lambda = ln(K + B) - ln(U)
//!
//!   drift(h):  dt = h/(K+B);  r_i += v_i dt;  t += dt      0 force evaluations
//!   kick(h):   dt = h/U;      v_i += a_i dt                1 force evaluation
//!   step:      drift(h/2) kick(h) drift(h/2)
//! ```
//!
//! On shell `K + B == U`, so the two denominators are the same number there and differ only off
//! it. That difference is where the method lives — and it is also why the most plausible
//! transcription error, swapping which denominator each half uses, is **invisible on shell**.
//! [`hamiltonian::Dens`] exists as a named pair so the swap can be constructed and asserted to
//! fire, and `tests/logh_hamiltonian_fd.rs` draws `B` independently of the state for exactly the
//! reason the `Gamma*` test draws a random off-shell `h`.
//!
//! # Two steppers, and one control
//!
//! [`step::Stepper`] selects RK4 (comparable to AZ and Heggie, four evaluations per step) or KDK
//! (the method as designed, one). [`hamiltonian::LhTime::None`] switches the transformation off
//! entirely, giving unregularised Cartesian dynamics **through the same code path** — so the
//! control is literally the regularisation removed and nothing else, and it inherits the event
//! sampling, escape rule and `t_end` that make its labels comparable.
//!
//! # No registration, and therefore nothing to count
//!
//! There is no `to_reg`, no `phys_from_state`, no reference body and no per-boundary rebuild.
//! A boundary sample is a read of the state, not a round trip through a chart. That absence is
//! the measurement.

pub mod driver;
pub mod gbs;
pub mod hamiltonian;
pub mod state;
pub mod step;
pub mod system;

pub use driver::{integrate_lh, LhDsMode, LhOpts, LhOut};
pub use hamiltonian::{Dens, LhTime};
pub use state::LhState;
pub use gbs::GbsOut;
pub use step::Stepper;
pub use system::LhSystem;
