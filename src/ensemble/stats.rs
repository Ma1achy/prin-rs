//! Robust ensemble statistics.
//!
//! `error_ratio` uses MAD internally rather than a standard deviation, because a std returns
//! NaN the moment one copy is non-finite — precisely the pathological pixel the statistic
//! exists to flag. Over pixels it aggregates by **max**, not median: max tracks damage at
//! Spearman +0.956 against +0.599 for median. Its magnitude is unstable, so it is a boolean
//! flag, not a measurement.

use crate::Real;

/// Sort key placing non-finite values at the top.
///
/// A non-finite energy is an *infinitely bad* outcome, not missing data, so the copy keeps
/// its slot in the ordering and pushes the median. This is how "never discard a copy"
/// survives contact with a median: the copy is counted, it is simply extreme. As long as
/// fewer than half the copies are non-finite the median is still meaningful.
fn ordered<T: Real>(v: &[T]) -> Vec<T> {
    let mut s: Vec<T> = v
        .iter()
        .map(|&x| if x.is_finite() { x } else { T::infinity() })
        .collect();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    s
}

/// Median, numpy convention for even counts: the mean of the two central order statistics.
pub fn median<T: Real>(v: &[T]) -> T {
    if v.is_empty() {
        return T::nan();
    }
    let s = ordered(v);
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) * T::lit(0.5)
    }
}

/// `1.4826 * median(|x - median(x)|)` — the consistent estimator of the standard deviation
/// for normally distributed data.
pub fn mad<T: Real>(v: &[T]) -> T {
    if v.len() < 2 {
        return T::zero();
    }
    let m = median(v);
    let dev: Vec<T> = v
        .iter()
        .map(|&x| if x.is_finite() { (x - m).abs() } else { T::infinity() })
        .collect();
    median(&dev) * T::lit(1.4826)
}

/// `sigma_E(t) / sigma_E(0)`, exactly 1.0 under exact dynamics.
///
/// Each trajectory conserves its own energy, so the ensemble's *spread* of energies is fixed
/// at `t = 0` and any growth is pure integration error — no threshold, no tuned constant.
///
/// **`sigma_E(0)` is proportional to the jitter and therefore to the cell width.** As
/// resolution rises `sigma_E(0)` shrinks while integration error does not, so the ratio
/// inflates for a purely trivial reason. Both `sigma_E` values are returned so the confound
/// is visible and correctable rather than baked into a single number.
pub fn error_ratio<T: Real>(e0: &[T], et: &[T]) -> (T, T, T) {
    let s0 = mad(e0);
    let st = mad(et);
    let ratio = if s0 > T::zero() { st / s0 } else { T::nan() };
    (ratio, s0, st)
}

/// Maximum absolute deviation from the median, with non-finite treated as infinite.
///
/// The **non-robust** companion to [`mad`], and it exists because robustness cuts both ways.
/// MAD is specified in BRIEF §4 so the statistic survives a non-finite copy — but a spread
/// estimator that a single wild copy cannot move is also one that cannot *see* a single wild
/// copy, which is precisely the damage `error_ratio` exists to flag. Measured: a pixel whose
/// worst copy drifted by `1.2e+02` reported a MAD-based `error_ratio` of 1.1369.
///
/// Reported alongside rather than instead: the two answer different questions, and which one
/// the acceptance test should use is a decision for data, not for taste.
pub fn max_dev<T: Real>(v: &[T]) -> T {
    if v.len() < 2 {
        return T::zero();
    }
    let m = median(v);
    let mut worst = T::zero();
    for &x in v {
        let d = if x.is_finite() { (x - m).abs() } else { T::infinity() };
        if d > worst {
            worst = d;
        }
    }
    worst
}

/// `error_ratio` built on [`max_dev`] instead of [`mad`]. Sensitive to a single damaged copy.
pub fn error_ratio_range<T: Real>(e0: &[T], et: &[T]) -> T {
    let s0 = max_dev(e0);
    let st = max_dev(et);
    if s0 > T::zero() {
        st / s0
    } else {
        T::nan()
    }
}

/// Fraction of copies not sharing the modal class, normalised by `1 - 1/(E+1)` so a
/// maximally split ensemble reads 1.0.
pub fn spread_event<T: Real>(classes: &[u8]) -> T {
    let n = classes.len();
    if n < 2 {
        return T::zero();
    }
    let mut best = 0usize;
    for c in classes {
        let k = classes.iter().filter(|x| *x == c).count();
        if k > best {
            best = k;
        }
    }
    let frac = T::one() - T::lit(best as f64) / T::lit(n as f64);
    frac / (T::one() - T::one() / T::lit(n as f64))
}
