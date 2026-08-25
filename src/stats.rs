//! Rank statistics, promoted out of the examples.
//!
//! [`spearman`] was written inline in `examples/open_items.rs` and, tie-unaware, again in
//! `examples/worst_pixels.rs`. The `+0.956 / +0.599` aggregation figures quoted in
//! `src/render.rs` and `src/ensemble/stats.rs` come from those. It is promoted here because the
//! criterion work reads it in the library and in several new examples, and a third copy would
//! be a third chance for the tie handling to drift.
//!
//! # Why a correlation is not on its own admissible
//!
//! Two standing rules bear directly on how these are used, and both are enforced by what this
//! module *returns* rather than by convention:
//!
//! - **"Never conclude 'no effect' from an aggregate without the per-pixel distribution."** A
//!   `rho` can only say the ordering did not move; it cannot say the items did not.
//!   [`rank_displacement`] returns the per-item movement so the distribution is always
//!   available beside the scalar.
//! - **"A difference can be small because both sides are right or because both are dead."**
//!   [`spearman`] returns `NaN` — never 0, never 1 — when either input has no variance at all,
//!   because a degenerate input is a measurement outcome and not a correlation of zero.

/// Ranks of `v`, ties by average rank, 1-based.
pub fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for &k in &idx[i..=j] {
            r[k] = avg;
        }
        i = j + 1;
    }
    r
}

/// Spearman rank correlation, ties handled by average rank.
///
/// `NaN` below three points or when either input is constant. Transcribed from
/// `examples/open_items.rs` without change, so the figures already on record still reproduce.
pub fn spearman(x: &[f64], y: &[f64]) -> f64 {
    if x.len() < 3 || x.len() != y.len() {
        return f64::NAN;
    }
    let (rx, ry) = (ranks(x), ranks(y));
    let n = x.len() as f64;
    let (mx, my) = (rx.iter().sum::<f64>() / n, ry.iter().sum::<f64>() / n);
    let num: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let dx: f64 = rx.iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let dy: f64 = ry.iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
    if dx == 0.0 || dy == 0.0 {
        f64::NAN
    } else {
        num / (dx * dy)
    }
}

/// `|rank_x(i) - rank_y(i)|` per item, normalised by `len - 1` so it reads as a fraction of
/// the list.
///
/// **This is the distribution behind the `rho`.** Two orderings can correlate at +0.99 while a
/// tail of items moves most of the way across the list, and that tail is exactly the population
/// a scheduler spends its budget on. Report the interdecile of this, not its variance — the
/// spread quantities in this project carry excess kurtosis 110 and the variance lives in the
/// tail.
pub fn rank_displacement(x: &[f64], y: &[f64]) -> Vec<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return Vec::new();
    }
    let (rx, ry) = (ranks(x), ranks(y));
    let d = (x.len() - 1) as f64;
    rx.iter().zip(&ry).map(|(a, b)| (a - b).abs() / d).collect()
}

/// Quantile by nearest rank on an already-collected sample. Mutates `v` by sorting.
pub fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    crate::quad::quantile(v, q)
}

/// `(p10, p50, p90)` and the interdecile width `p90 - p10`.
///
/// The interdecile rather than the variance, per the standing rule: on `alpha_shape` the excess
/// kurtosis is 110, `interdecile/sd` is 0.866 against a normal's 2.563, and the Halton switch
/// cut the variance 267,000x while moving the interdecile not at all.
pub fn interdecile(v: &[f64]) -> (f64, f64, f64, f64) {
    let mut s: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if s.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    let p10 = quantile(&mut s.clone(), 0.10);
    let p50 = quantile(&mut s.clone(), 0.50);
    let p90 = quantile(&mut s, 0.90);
    (p10, p50, p90, p90 - p10)
}
