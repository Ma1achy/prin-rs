//! Ensemble statistics.
//!
//! **Two independent decisions both land on the word "max", for different reasons. Do not
//! collapse them.**
//!
//! *Within* a footprint, `error_ratio` uses the **maximum deviation** from the median. BRIEF §4
//! requires the inner statistic survive a non-finite copy, because a standard deviation returns
//! NaN on exactly the pathological pixel the field exists to flag. MAD satisfies that but
//! overshoots: an estimator a single wild copy cannot move is one that cannot *see* a single
//! wild copy, and with eight copies one bad value sits above the median of eight deviations and
//! is arithmetically invisible. Measured on 23 damaged against 1001 healthy pixels, the
//! damaged/healthy separation is **1.06 with MAD and 59.51 with max deviation**; a pixel whose
//! worst copy drifted 120x the total energy reported a MAD ratio of 1.1369, inside the healthy
//! p99 of 1.0756. Max deviation meets the original requirement *better*: a non-finite copy
//! gives an infinite deviation, which is the correct answer where a std gives NaN — and NaN is
//! indistinguishable from "could not compute".
//!
//! *Across* footprints, [`crate::render`] aggregates by max rather than median, at Spearman
//! +0.956 against +0.599. That is an aggregation-level result about `error_ratio_max` versus
//! `error_ratio_median`, correlated against per-quad exponent damage. It is not a statement
//! about which estimator goes inside a footprint, and the two must not be conflated when
//! either is re-measured.
//!
//! The magnitude of the ratio is unstable either way: it is a boolean flag, not a measurement.

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

/// `sigma_E(t) / sigma_E(0)`, exactly 1.0 under exact dynamics. **This is the field.**
///
/// Each trajectory conserves its own energy, so the ensemble's *spread* of energies is fixed
/// at `t = 0` and any growth is pure integration error — no threshold, no tuned constant.
///
/// `sigma_E` here is [`max_dev`]; see the module docs for why not MAD. [`error_ratio_mad`] is
/// kept and dumped alongside so the change stays measurable rather than asserted.
///
/// **`sigma_E(0)` is proportional to the jitter and therefore to the cell width.** As
/// resolution rises `sigma_E(0)` shrinks while integration error does not, so the ratio
/// inflates for a purely trivial reason. Both `sigma_E` values are returned so the confound
/// is visible and correctable rather than baked into a single number.
pub fn error_ratio<T: Real>(e0: &[T], et: &[T]) -> (T, T, T) {
    let s0 = max_dev(e0);
    let st = max_dev(et);
    let ratio = if s0 > T::zero() { st / s0 } else { T::nan() };
    (ratio, s0, st)
}

/// Maximum absolute deviation from the median, with non-finite treated as infinite.
///
/// The deliberately **non-robust** spread estimator, and the one [`error_ratio`] is built on.
/// Robustness cuts both ways: see the module docs. A non-finite copy yields an infinite
/// deviation — the correct answer, and never a NaN.
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

/// `error_ratio` built on [`mad`] instead of [`max_dev`] — BRIEF §4's original wording.
///
/// Retained and dumped so the switch to [`error_ratio`] stays a measured difference rather
/// than a silent replacement. Insensitive to a single damaged copy, which is why it is no
/// longer the field.
pub fn error_ratio_mad<T: Real>(e0: &[T], et: &[T]) -> T {
    let s0 = mad(e0);
    let st = mad(et);
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
