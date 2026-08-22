//! Aarseth-Zare global regularisation (planar three-body).
//!
//! Levi-Civita regularises **one** pair, so a close approach of either other pair still
//! enters through an unregularised term — measured: near-field drift 5.0e-04 unregularised
//! became 1.7e+01 with LC on the wrong pair. AZ regularises **two** pairs at once.
//!
//! The reference body `a` is the one **not** in the longest side, so both regularised pairs
//! share it and the unregularised side `(b,c)` is the longest. Then
//! `|R3| >= max(|R1|,|R2|)`, so `R3 -> 0` only in a genuine triple collision, which is
//! provably non-regularisable anyway.
//!
//! Ported from `reference/tb_az.py`. **The algebra is not re-derived** — it is error-prone
//! and fails silently; two sign errors in the reference were invisible until someone
//! finite-differenced the Hamiltonian.

pub mod driver;
pub mod hamiltonian;
pub mod lc;
pub mod reference_body;
pub mod rk4;
pub mod state;
pub mod system;

pub use driver::{integrate_az, integrate_with_policy, AzOut};
pub use reference_body::RefPolicy;
pub use state::AzState;
pub use system::AzSystem;
