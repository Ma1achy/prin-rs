//! The Burrau (Pythagorean) configuration: masses 3, 4, 5 at rest at the vertices of a
//! 3-4-5 triangle. Released from rest, so `L_z = 0` for every trajectory — which is why
//! there is no `L_z` analogue of `error_ratio` (the ratio would be 0/0).

use crate::physics::Cart;
use crate::{Real, Vec2};

/// Reference masses, in `PAIRS`-compatible body order.
pub const MASSES: [f64; 3] = [3.0, 4.0, 5.0];

/// Reference positions.
pub const R0: [[f64; 2]; 3] = [[1.0, 3.0], [-2.0, -1.0], [1.0, -1.0]];

pub fn masses<T: Real>() -> [T; 3] {
    [T::lit(MASSES[0]), T::lit(MASSES[1]), T::lit(MASSES[2])]
}

/// The nominal Burrau state: at rest, so all velocities are zero.
pub fn state<T: Real>() -> Cart<T> {
    Cart {
        r: [
            Vec2::new(T::lit(R0[0][0]), T::lit(R0[0][1])),
            Vec2::new(T::lit(R0[1][0]), T::lit(R0[1][1])),
            Vec2::new(T::lit(R0[2][0]), T::lit(R0[2][1])),
        ],
        v: [Vec2::zero(); 3],
    }
}
