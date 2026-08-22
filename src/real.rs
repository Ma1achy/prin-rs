//! The precision chokepoint.
//!
//! One supertrait, so bounds are written once per `impl` block and never per function.
//!
//! The associated floors are **not cosmetic**. The NumPy reference guards divisions with
//! literals chosen for f64, and two of them are not representable at f32:
//!
//! - `1e-300` underflows to zero at f32 (min normal ~1.18e-38), so the guard stops guarding.
//! - `1e-15` is far below f32 epsilon (~1.19e-7), so the sync test `t < t_target - SYNC_EPS`
//!   degenerates to `t < t_target`.
//!
//! f64 therefore carries the reference's literal values, so the Step-4 cross-check is a clean
//! equality; f32 carries precision-appropriate ones. This is the one place where f32 and f64
//! are not running the same algorithm, and it lives exactly where degenerate geometry lives.
//! It must be stated up front in the f32 PR.

use num_traits::Float;
use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

pub trait Real:
    Float
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + Sum
    + Send
    + Sync
    + Debug
    + Display
    + Default
    + 'static
{
    /// Norm / denominator guard. Reference: `1e-300`.
    const TINY: Self;
    /// Denominator floor for the relative energy drift. Reference: `1e-30`.
    const DRIFT_FLOOR: Self;
    /// Slack in the sync-boundary time test. Reference: `1e-15`.
    const SYNC_EPS: Self;
    /// Distance floor in the escape test. Reference: `1e-12`.
    const DIST_FLOOR: Self;
    /// Distance floor in the legacy classifier. Reference: `1e-9`.
    const CLASSIFY_FLOOR: Self;

    /// A literal, without `T::from(x).unwrap()` at every call site. Constant-folds.
    fn lit(x: f64) -> Self;
}

impl Real for f64 {
    const TINY: Self = 1e-300;
    const DRIFT_FLOOR: Self = 1e-30;
    const SYNC_EPS: Self = 1e-15;
    const DIST_FLOOR: Self = 1e-12;
    const CLASSIFY_FLOOR: Self = 1e-9;

    #[inline(always)]
    fn lit(x: f64) -> Self {
        x
    }
}

impl Real for f32 {
    // 1e-300 is zero at f32. Smallest normal is ~1.18e-38; leave headroom so squaring a
    // guarded value does not itself underflow.
    const TINY: Self = 1e-37;
    const DRIFT_FLOOR: Self = 1e-30;
    // f32 eps is ~1.19e-7, and at t ~ 13 the ulp is ~1e-6. A 1e-15 slack is no slack at all.
    const SYNC_EPS: Self = 1e-6;
    const DIST_FLOOR: Self = 1e-12;
    const CLASSIFY_FLOOR: Self = 1e-9;

    #[inline(always)]
    fn lit(x: f64) -> Self {
        x as f32
    }
}

/// `T::lit(0.5)`, shorter.
#[macro_export]
macro_rules! r {
    ($x:expr) => {
        <T as $crate::real::Real>::lit($x)
    };
}
