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
/// `n_hot_within == N^2` in **98.8%** of the 75,359 committed `charts/` leaves (**87.1%** over
/// the whole 92,880-leaf corpus) and `n_components == 1` in **99.6%** (**92.6%**). Two scopes,
/// both stated: the chart dumps are the saturated end and the zoom ladders the unsaturated one,
/// and quoting one under the other's name is how this got written up wrong the first time.
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

// -------------------------------------------------------------------------------------------
// Straightness of a connected structure.
// -------------------------------------------------------------------------------------------

/// One connected component's **straightness**: `sqrt(lambda_2 / lambda_1)` of its coordinate
/// covariance. `0.0` is a perfect line, `1.0` is isotropic.
///
/// This is total least squares in closed form. `lambda_2` is the mean squared **perpendicular**
/// distance to the best-fit line and `lambda_1` the mean squared distance along it, so the ratio
/// is the RMS residual divided by the extent — the same quantity as *"RMS 4.28 px over 621
/// rows"*, made scale-free and computable without a hand-drawn mask.
///
/// # Why this discriminates, and why density could not
///
/// A decision boundary and a fractal boundary both raise the density of sharp neighbour steps;
/// only the first is *straight*. The wedge edge in `config_stability` was measured at RMS 4.28 px
/// over 621 rows against a chaotic ribbon's 42.05 px over 186 in the same image — a tenfold
/// separation — but by hand, on one edge, under one integrator. This is that measurement as a
/// function.
///
/// # What it is not
///
/// **Curvature is indistinguishable from scatter here.** A circular arc and a cloud of the same
/// covariance score alike; the eigenvalue ratio knows nothing about ordering along the curve.
/// That is a real limit and it is the reason [`straightness`] is reported beside the component
/// size rather than alone: a large component that is straight is a boundary, a large component
/// that is not may be an arc *or* a blob, and this number cannot say which.
///
/// Returns `NaN` for fewer than three pixels — two points are collinear by construction, so a
/// number there would be an artefact of the count rather than a measurement of shape.
pub fn straightness(pts: &[(usize, usize)]) -> f64 {
    if pts.len() < 3 {
        return f64::NAN;
    }
    let n = pts.len() as f64;
    let (mut mx, mut my) = (0.0f64, 0.0f64);
    for &(x, y) in pts {
        mx += x as f64;
        my += y as f64;
    }
    mx /= n;
    my /= n;
    let (mut sxx, mut syy, mut sxy) = (0.0f64, 0.0f64, 0.0f64);
    for &(x, y) in pts {
        let (dx, dy) = (x as f64 - mx, y as f64 - my);
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    sxx /= n;
    syy /= n;
    sxy /= n;
    // Eigenvalues of the symmetric 2x2 [[sxx, sxy], [sxy, syy]].
    let tr = sxx + syy;
    let det = sxx * syy - sxy * sxy;
    let disc = (tr * tr / 4.0 - det).max(0.0).sqrt();
    let (l1, l2) = (tr / 2.0 + disc, tr / 2.0 - disc);
    if l1 <= 0.0 {
        return f64::NAN;
    }
    (l2.max(0.0) / l1).sqrt()
}

/// Every 4-connected component of `mask`, as pixel coordinate lists, largest first.
///
/// [`layout`] returns component *counts*; this returns the components themselves, which is what
/// a shape statistic needs. Same 4-connectivity convention: diagonal contact is not connection.
pub fn components(mask: &[bool], n: usize) -> Vec<Vec<(usize, usize)>> {
    assert_eq!(mask.len(), n * n, "mask must be n*n");
    let mut seen = vec![false; n * n];
    let mut out: Vec<Vec<(usize, usize)>> = Vec::new();
    for s in 0..n * n {
        if !mask[s] || seen[s] {
            continue;
        }
        let (mut stack, mut comp) = (vec![s], Vec::new());
        seen[s] = true;
        while let Some(i) = stack.pop() {
            let (x, y) = (i % n, i / n);
            comp.push((x, y));
            for (dx, dy) in [(1i32, 0i32), (0, 1), (-1, 0), (0, -1)] {
                let (u, v) = (x as i32 + dx, y as i32 + dy);
                if u < 0 || v < 0 || u >= n as i32 || v >= n as i32 {
                    continue;
                }
                let j = v as usize * n + u as usize;
                if mask[j] && !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        out.push(comp);
    }
    out.sort_by_key(|c| std::cmp::Reverse(c.len()));
    out
}

/// Boundary pixels of `mask`: in the set, with a 4-neighbour outside it or off the grid.
pub fn boundary(mask: &[bool], n: usize) -> Vec<(usize, usize)> {
    assert_eq!(mask.len(), n * n, "mask must be n*n");
    let mut out = Vec::new();
    for y in 0..n {
        for x in 0..n {
            if !mask[y * n + x] {
                continue;
            }
            let edge = [(1i32, 0i32), (0, 1), (-1, 0), (0, -1)].iter().any(|&(dx, dy)| {
                let (u, v) = (x as i32 + dx, y as i32 + dy);
                u < 0 || v < 0 || u >= n as i32 || v >= n as i32 || !mask[v as usize * n + u as usize]
            });
            if edge {
                out.push((x, y));
            }
        }
    }
    out
}

/// **Median LOCAL straightness of a region's boundary**, at scale `radius`.
///
/// The global [`straightness`] of a closed boundary is ~1 whatever its shape — a square outline
/// and a blob outline are both isotropic — so a region's edge has to be fitted *locally*. For
/// every boundary pixel, the boundary pixels within `radius` are fitted by total least squares
/// and the median residual ratio is returned. A straight edge is straight in every window; a
/// fractal edge wiggles inside each one.
///
/// This is the hand measurement generalised: *"pale edge RMS 4.28 px over 621 rows against the
/// red band's 42.05 px over 186"* was one long fit on one edge, chosen by eye. This is the same
/// quantity over every boundary pixel of every component, with no mask drawn by hand.
///
/// # The limit, stated rather than discovered later
///
/// **Local straightness at `radius` cannot tell a straight line from a curve whose radius of
/// curvature is much larger than `radius`.** A big circle reads straight at small `radius` and
/// that is not a defect — it is what "local" means. Pick `radius` from the structure being
/// tested and quote it with the number.
///
/// `NaN` when fewer than `min_pts` boundary pixels have a full enough neighbourhood, because a
/// median over a handful of windows describes the sample and not the shape.
pub fn boundary_straightness(mask: &[bool], n: usize, radius: usize, min_pts: usize) -> f64 {
    let b = boundary(mask, n);
    if b.len() < min_pts {
        return f64::NAN;
    }
    let r2 = (radius * radius) as i64;
    let mut vals: Vec<f64> = Vec::with_capacity(b.len());
    for &(x, y) in &b {
        let local: Vec<(usize, usize)> = b
            .iter()
            .copied()
            .filter(|&(u, v)| {
                let (dx, dy) = (u as i64 - x as i64, v as i64 - y as i64);
                dx * dx + dy * dy <= r2
            })
            .collect();
        // A window needs enough points to have a shape; `radius` pixels is a bare line's worth.
        if local.len() >= radius.max(3) {
            let s = straightness(&local);
            if s.is_finite() {
                vals.push(s);
            }
        }
    }
    if vals.len() < min_pts {
        return f64::NAN;
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vals[vals.len() / 2]
}
