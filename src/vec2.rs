//! Planar vectors, with operators, so the AZ algebra reads as vector algebra rather than
//! index soup. This is the single biggest reduction in transcription risk in the port.

use crate::Real;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec2<T> {
    pub x: T,
    pub y: T,
}

impl<T: Real> Vec2<T> {
    #[inline(always)]
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    #[inline(always)]
    pub fn zero() -> Self {
        Self { x: T::zero(), y: T::zero() }
    }

    #[inline(always)]
    pub fn dot(self, o: Self) -> T {
        self.x * o.x + self.y * o.y
    }

    /// The 2D cross product (a scalar): `self.x * o.y - self.y * o.x`.
    #[inline(always)]
    pub fn cross(self, o: Self) -> T {
        self.x * o.y - self.y * o.x
    }

    #[inline(always)]
    pub fn norm_sq(self) -> T {
        self.dot(self)
    }

    #[inline(always)]
    pub fn norm(self) -> T {
        self.norm_sq().sqrt()
    }

    #[inline(always)]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// Cast to another precision. Used to generate initial conditions once in f64 and
    /// cast down, so an f32/f64 comparison isolates arithmetic from IC differences.
    #[inline(always)]
    pub fn cast<U: Real>(self) -> Vec2<U> {
        Vec2 {
            x: U::lit(self.x.to_f64().unwrap()),
            y: U::lit(self.y.to_f64().unwrap()),
        }
    }
}

impl<T: Real> Add for Vec2<T> {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y)
    }
}

impl<T: Real> Sub for Vec2<T> {
    type Output = Self;
    #[inline(always)]
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y)
    }
}

impl<T: Real> Neg for Vec2<T> {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y)
    }
}

/// Scalar multiply. Written `v * s` (not `s * v`) — a blanket `impl Mul<Vec2<T>> for T`
/// is not permitted for a generic `T`.
impl<T: Real> Mul<T> for Vec2<T> {
    type Output = Self;
    #[inline(always)]
    fn mul(self, s: T) -> Self {
        Self::new(self.x * s, self.y * s)
    }
}

impl<T: Real> Div<T> for Vec2<T> {
    type Output = Self;
    #[inline(always)]
    fn div(self, s: T) -> Self {
        Self::new(self.x / s, self.y / s)
    }
}

impl<T: Real> AddAssign for Vec2<T> {
    #[inline(always)]
    fn add_assign(&mut self, o: Self) {
        self.x += o.x;
        self.y += o.y;
    }
}

impl<T: Real> SubAssign for Vec2<T> {
    #[inline(always)]
    fn sub_assign(&mut self, o: Self) {
        self.x -= o.x;
        self.y -= o.y;
    }
}
