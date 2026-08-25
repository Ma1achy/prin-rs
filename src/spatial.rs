//! Where the hot footprints sit, not just how many.
//!
//! A quad computes `N x N` footprint values and keeps one scalar. **The layout of the hot set
//! is free information** — O(N^2) over data already in hand, against 512 trajectories per quad —
//! and it distinguishes the two cases that want opposite decisions:
//!
//! - a **boundary** shows as a *connected, thin* structure -> split
//! - **chaos** shows as *scattered* hot footprints -> floor
//!
//! Every current aggregate conflates them: a quad with 5% hot footprints has a filament, one
//! with a high median is uniformly blurred, and mean/median/p90 cannot tell those apart.
//!
//! # The perimeter convention, stated because it changes what the number means
//!
//! `perimeter` counts **internal edges only** — an edge between a hot cell and a cold cell
//! inside the grid. Edges on the grid border are not counted.
//!
//! The alternative (treating outside as cold) makes every quad's border contribute perimeter
//! whether or not there is structure there, so a featureless fully-hot quad would read as
//! having a boundary. Under this convention it reads 0, which is right. Measured against the
//! shapes that matter: an isolated cell gives ratio 4.0, a one-cell-wide filament ~2.0, a
//! compact blob of area `A` ~`4/sqrt(A)`, a half-plane `2/N`, and fully hot exactly 0. Thin is
//! high, blobby is low, and the separation is wide.

/// What the hot set looks like. All counts are of footprints, not of pixels.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Layout {
    pub n_hot: u32,
    /// Connected components of the hot set under **4-connectivity**. Diagonal-only contact is
    /// not connection: a checkerboard is scattered, which is the reading that matters here.
    pub n_components: u32,
    pub largest_component: u32,
    /// `perimeter / n_hot`. See the module docs for the convention and the reference shapes.
    /// `NaN` when nothing is hot — **not 0**, which would be indistinguishable from a
    /// featureless fully-hot quad.
    pub perimeter_ratio: f64,
}

impl Layout {
    /// Fraction of the `N^2` footprints above the threshold.
    pub fn frac_hot(&self, n: usize) -> f64 {
        self.n_hot as f64 / (n * n) as f64
    }

    /// Is the hot set a *connected thin structure* rather than scatter?
    ///
    /// Reported alongside the raw fields, never instead of them: it is one reading of four
    /// numbers and the thresholds are stated here rather than buried in a decision.
    pub fn looks_like_boundary(&self, n: usize, thin: f64) -> bool {
        self.n_hot > 0
            && self.largest_component as f64 >= 0.5 * self.n_hot as f64
            && self.largest_component >= (n as u32) / 2
            && self.perimeter_ratio >= thin
    }
}

/// Connected components and perimeter of `mask`, an `n x n` grid in row-major order
/// (`index = jy*n + jx`, x fastest — the same ordering as `Slice::decode_pos`).
pub fn layout(mask: &[bool], n: usize) -> Layout {
    assert_eq!(mask.len(), n * n, "mask must be n*n");
    let idx = |jx: usize, jy: usize| jy * n + jx;

    let n_hot = mask.iter().filter(|&&h| h).count() as u32;
    if n_hot == 0 {
        return Layout { n_hot: 0, n_components: 0, largest_component: 0, perimeter_ratio: f64::NAN };
    }

    // Internal edges between a hot cell and a cold one. Each edge is visited once by scanning
    // only the +x and +y neighbours.
    let mut perimeter = 0u32;
    for jy in 0..n {
        for jx in 0..n {
            let h = mask[idx(jx, jy)];
            if jx + 1 < n && h != mask[idx(jx + 1, jy)] {
                perimeter += 1;
            }
            if jy + 1 < n && h != mask[idx(jx, jy + 1)] {
                perimeter += 1;
            }
        }
    }

    // Flood fill, 4-connectivity. An explicit stack rather than recursion: at N = 8 the depth
    // is trivial, but this is called once per quad on every descent and a stack overflow in a
    // scheduler is not a failure mode worth leaving available.
    let mut seen = vec![false; n * n];
    let mut stack: Vec<usize> = Vec::new();
    let (mut n_components, mut largest) = (0u32, 0u32);
    for start in 0..n * n {
        if !mask[start] || seen[start] {
            continue;
        }
        n_components += 1;
        let mut size = 0u32;
        seen[start] = true;
        stack.push(start);
        while let Some(c) = stack.pop() {
            size += 1;
            let (jx, jy) = (c % n, c / n);
            let push = |k: usize, seen: &mut Vec<bool>, stack: &mut Vec<usize>| {
                if mask[k] && !seen[k] {
                    seen[k] = true;
                    stack.push(k);
                }
            };
            if jx > 0 {
                push(idx(jx - 1, jy), &mut seen, &mut stack);
            }
            if jx + 1 < n {
                push(idx(jx + 1, jy), &mut seen, &mut stack);
            }
            if jy > 0 {
                push(idx(jx, jy - 1), &mut seen, &mut stack);
            }
            if jy + 1 < n {
                push(idx(jx, jy + 1), &mut seen, &mut stack);
            }
        }
        largest = largest.max(size);
    }

    Layout {
        n_hot,
        n_components,
        largest_component: largest,
        perimeter_ratio: perimeter as f64 / n_hot as f64,
    }
}

/// How a footprint is called **hot**.
///
/// # Both rules are computed on every quad, and that is not redundancy
///
/// The obvious reading of "make the threshold relative" is to replace the absolute one. It is
/// wrong twice.
///
/// **`n_hot` stops being a signal under any quantile rule.** On a field with distinct values the
/// count above the cut is set by the rule, not the field — 31 of 64 at `N = 8, q = 0.5`, given
/// nearest-rank and a strict comparison. So `frac_hot` carries essentially no information once
/// the mask is relative. Under [`HotRule::Quantile`] the whole signal is the *shape* of the mask:
/// `n_components`, `largest_component`, `perimeter_ratio`.
///
/// The one exception, measured rather than assumed: on a **tied** field the count is set by the
/// tie structure. A two-valued field reads the same count at `q = 0.5, 0.75, 0.9` alike — which
/// is the case that occurs when the event arm, with five distinct values, dominates a footprint
/// field.
///
/// **And `frac_hot_between/median` is the best criterion measured on this project** — the only one
/// beating the random band in both measurable regions. Replacing the absolute mask would have
/// deleted the best-performing signal in the system and read as an improvement.
///
/// So the absolute rule keeps `frac_above_tau_*` exactly as it was, and the relative rule is added
/// beside it to desaturate the shape statistics. Measured cause of the saturation: with
/// `tau_display = 1e-4` sitting at the **0.4th percentile** of the observed spread distribution,
/// `n_hot_within == N^2` in **98.8%** of the 75,359 committed leaves and `n_components == 1` in
/// **99.5%** of them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HotRule {
    /// Above a fixed level. The shipped rule, and the cause of the saturation above.
    AbsTau(f64),
    /// Above the quad's **own** `q`-quantile. A shape statistic rather than a magnitude one, so
    /// it does not drift as the signal rises globally with `t`.
    Quantile(f64),
}

impl Default for HotRule {
    fn default() -> Self {
        HotRule::Quantile(0.5)
    }
}

impl HotRule {
    pub fn name(self) -> String {
        match self {
            HotRule::AbsTau(t) => format!("abs[{t:.3e}]"),
            HotRule::Quantile(q) => format!("q[{q:.2}]"),
        }
    }
}

/// The hot mask for `vals` under `rule`.
///
/// **Non-finite is hot**, under both rules and for the same reason the absolute form already had
/// it: a footprint that could not be determined is not evidence of calm, and treating it as cold
/// hides the pathological case from exactly the statistic built to find structure.
///
/// A quad with fewer than two finite values has no distribution to take a quantile of. It yields
/// an **all-hot** mask — undetermined, not resolved — rather than an empty one, which would read
/// as a calm quad. `Layout::n_hot == vals.len()` is how that case is countable downstream.
pub fn hot_mask(vals: &[f64], rule: HotRule) -> Vec<bool> {
    let cut = match rule {
        HotRule::AbsTau(t) => t,
        HotRule::Quantile(q) => {
            let mut finite: Vec<f64> = vals.iter().cloned().filter(|x| x.is_finite()).collect();
            if finite.len() < 2 {
                return vec![true; vals.len()];
            }
            crate::quad::quantile(&mut finite, q)
        }
    };
    vals.iter().map(|&v| !v.is_finite() || v > cut).collect()
}

/// RMS of the forward-difference gradient across the `n x n` footprint grid.
///
/// The magnitude companion to [`layout`]: `layout` says where the hot set sits, this says how
/// fast the field moves, and neither needs a threshold to say it.
///
/// **`NaN` when no adjacent pair is finite**, never 0 — the same convention as
/// `scheduler::termination_gradient`, and for the same reason. A zero here would be a null
/// presented as a measurement about the field, and indistinguishable from a genuinely flat quad.
pub fn grad_rms(vals: &[f64], n: usize) -> f64 {
    assert_eq!(vals.len(), n * n, "field must be n*n");
    let idx = |jx: usize, jy: usize| jy * n + jx;
    let (mut acc, mut pairs) = (0.0f64, 0usize);
    let take = |a: f64, b: f64, acc: &mut f64, pairs: &mut usize| {
        if a.is_finite() && b.is_finite() {
            let d = a - b;
            *acc += d * d;
            *pairs += 1;
        }
    };
    for jy in 0..n {
        for jx in 0..n {
            let a = vals[idx(jx, jy)];
            if jx + 1 < n {
                take(a, vals[idx(jx + 1, jy)], &mut acc, &mut pairs);
            }
            if jy + 1 < n {
                take(a, vals[idx(jx, jy + 1)], &mut acc, &mut pairs);
            }
        }
    }
    if pairs == 0 {
        f64::NAN
    } else {
        (acc / pairs as f64).sqrt()
    }
}
