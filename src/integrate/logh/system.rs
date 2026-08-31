//! The constants of a logH march. There are almost none, and that is the finding.
//!
//! `AzSystem` carries a reference body, two reduced masses and the pair assignment that follows
//! from the choice; `HgSystem` carries three reduced masses, three inverse masses and the mass
//! products of Heggie's enlarged space. This carries the masses. Nothing here depends on the
//! configuration, so nothing here is ever rebuilt — there is no `to_reg`, no re-registration, and
//! no per-boundary refresh, because there is no chart to register into.
//!
//! It exists as a type rather than a bare `&[T; 3]` so the three drivers have the same signature
//! shape and a comparison harness can hold one against another without an argument-order
//! difference standing in for a physics difference.

use crate::Real;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LhSystem<T> {
    pub masses: [T; 3],
    pub mtot: T,
}

impl<T: Real> LhSystem<T> {
    pub fn new(masses: [T; 3]) -> Self {
        Self { masses, mtot: masses[0] + masses[1] + masses[2] }
    }

    /// `B = U - K`, evaluated once at registration and never again.
    ///
    /// This is `-E`, so it is the physical energy with a sign, frozen exactly the way AZ's `E`
    /// and Heggie's `h` are. Under autonomous Newtonian gravity it is **constant along the
    /// march** — Mikkola's general `dB = (ds/U) dU/dt` term vanishes when the force law has no
    /// explicit time dependence and no velocity dependence. It is non-zero in AR-CHAIN because
    /// of the post-Newtonian terms, and if those are ever ported this is the line that changes.
    pub fn b_of(&self, c: &crate::physics::Cart<T>) -> T {
        crate::physics::energy::potential_pos(&c.r, &self.masses, T::zero())
            - crate::physics::energy::kinetic(&c.v, &self.masses)
    }
}
