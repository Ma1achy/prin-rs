//! Newtonian three-body physics in the plane. `G = 1` throughout.

pub mod burrau;
pub mod decoder;
pub mod energy;
pub mod ftle;
pub mod newton;
pub mod shape;

use crate::{Real, Vec2};

/// Gravitational constant. The whole project works in `G = 1` units.
pub const G: f64 = 1.0;

/// Canonical pair ordering. **Load-bearing** — the third-body lookup `THIRD[k]` and every
/// pair-indexed field in the raw dump depend on it.
pub const PAIRS: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

/// `THIRD[k]` is the body not in `PAIRS[k]`.
pub const THIRD: [usize; 3] = [2, 1, 0];

/// A full planar three-body state: 12 numbers.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cart<T> {
    pub r: [Vec2<T>; 3],
    pub v: [Vec2<T>; 3],
}

/// A decoded initial condition: **masses and state together**.
///
/// Masses used to be the compile-time constant [`burrau::MASSES`], read at three separate sites.
/// Four of the five chart families in the chart reference vary them — the latent chart carries
/// `(z_mu1, z_mu2)` as two of its eight coordinates, and the mass simplex is nothing else — so a
/// mass cannot be a global. It is a property of the point being decoded, and travels with it.
///
/// **The guard that makes this safe:** every chart that existed before returns exactly
/// `burrau::MASSES`, so `BodyPlane` stays bitwise and the Python cross-check is untouched.
/// `tests/charts.rs` asserts it rather than trusting it.
///
/// [`burrau::MASSES`]: crate::physics::burrau::MASSES
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ic<T> {
    pub m: [T; 3],
    pub s: Cart<T>,
}

impl<T: Real> Ic<T> {
    pub fn cast<U: Real>(&self) -> Ic<U> {
        Ic {
            m: [
                U::lit(self.m[0].to_f64().unwrap()),
                U::lit(self.m[1].to_f64().unwrap()),
                U::lit(self.m[2].to_f64().unwrap()),
            ],
            s: self.s.cast(),
        }
    }
    pub fn is_finite(&self) -> bool {
        self.m.iter().all(|x| x.is_finite()) && self.s.is_finite()
    }
}

impl<T: Real> Cart<T> {
    pub fn new(r: [Vec2<T>; 3], v: [Vec2<T>; 3]) -> Self {
        Self { r, v }
    }

    pub fn is_finite(&self) -> bool {
        self.r.iter().chain(self.v.iter()).all(|p| p.is_finite())
    }

    /// Cast to another precision. Initial conditions are generated once in f64 and cast
    /// down, never generated separately per precision — otherwise an IC difference is
    /// indistinguishable from a genuine f32 arithmetic effect.
    pub fn cast<U: Real>(&self) -> Cart<U> {
        Cart {
            r: [self.r[0].cast(), self.r[1].cast(), self.r[2].cast()],
            v: [self.v[0].cast(), self.v[1].cast(), self.v[2].cast()],
        }
    }
}
