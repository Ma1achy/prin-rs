//! §2-§5 — audit every available signal against the DP's own labels.
//!
//! # Why this exists
//!
//! PR #22's honest reading is **"no criterion reliably clears breadth-first"**. That has three
//! ways of being wrong, and [`prin_rs::metric::Cache::dp_optimal`] is what makes them testable:
//! it gives, for the first time, a **ground-truth label per quad** -- did the exact optimum split
//! this one -- which no previous exclusion had.
//!
//! 1. **Wrong axis.** `error(B)` has `B` in quads, and a cost-aware criterion could lose per quad
//!    and win per second.
//! 2. **Wrong normalisation.** `ensemble_spread = max(...)`, so a field living an order of
//!    magnitude lower can carry real information and never win the `max`. That is how
//!    **diffusion** was excluded -- a statement about *scale* recorded as a statement about
//!    *information* -- and diffusion and FTLE have **never been scored by `error(B)` at all**.
//! 3. **Wrong signal.** A combination could beat every single one.
//!
//! If all three come back negative, *"on this manifold adaptive refinement does not beat uniform
//! sampling at equal cost"* is a real result and a better outcome than a criterion that looks
//! clever and is not.
//!
//! # One premise corrected before it is built on
//!
//! The cost figure this audit was commissioned against -- quads varying **100x** in
//! `total_substeps`, p1 `2.04e4` / p50 `4.09e4` / p99 `2.05e6`, bimodal -- is a
//! **between-configuration** spread pooled over 68,685 leaves from 60 dumps at different regions,
//! playheads and configurations, quoted as a within-region one. The committed within-region
//! measurement (`results/output/cost_and_anisotropy.txt`, level 5) reads `max/p50` of **1.00** in
//! `far`, **1.03** in `near-field`, **3.02** in `deep interior`.
//!
//! **Where cost is constant, `error(C)` is `error(B)` rescaled and cannot reorder anything.** So
//! stage 1 measures the distribution and stops. And it measures it **per level**: every quad
//! carries the same `N^2 x (E+1)` trajectories whatever its level, so a level-7 quad lying wholly
//! inside a collision zone is all-hot where a level-0 quad averages the region. The spread should
//! widen with depth by that mechanism, and a pooled p1/p99 across levels would hide exactly what
//! going to level 7 is meant to reveal.
//!
//! # Writes
//!
//! `<root>/output/signal_audit.txt` is this file's stdout; `<root>/audit/*.tsv` is the per-quad
//! table. **The root is an argument** -- an output path that cannot be redirected is the same
//! defect as an argument hardcoded past, and it has already cost this project two round trips.

use std::collections::HashMap;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Cache, Key, Rank};
use prin_rs::physics::ftle::FtleOpts;
use prin_rs::quad::{Agg, Criterion};
use prin_rs::stats;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const SHIPPING: metric::Colouring =
    metric::Colouring::Bivariate(prin_rs::output::colour::Scalar::ShapeSpread);

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { f64::NAN } else { stats::quantile(v, p) }
}

/// Quantiles of a finite-filtered sample, or all-`NaN` if nothing is finite.
fn quant(v: &[f64], ps: &[f64]) -> Vec<f64> {
    let mut f: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    ps.iter().map(|p| q(&mut f, *p)).collect()
}

// ---------------------------------------------------------------------------------------------
// The signal table
// ---------------------------------------------------------------------------------------------

/// One named column of the per-quad table.
struct Col {
    name: &'static str,
    /// What the column means, for the report -- not decoration: a column of `NaN` and a column
    /// that is structurally undefined are different findings.
    v: Vec<f64>,
}

/// The base columns: everything on `QuadReduction`, plus the two that have never been scored.
///
/// **No pre-filtering.** The last filter excluded a signal on a scale artefact and it took three
/// months to notice, so every scalar goes in and the ones at Spearman ~0 are reported as findings
/// rather than omitted for tidiness.
fn base_columns(cache: &Cache, px_of: &HashMap<Key, Vec<PixelOut>>, keys: &[Key]) -> Vec<Col> {
    // `diffusion` and `ftle` are per-footprint scalars of the NOMINAL copy only, and are NaN
    // unless `EnsembleCfg::ftle` is set. They are deliberately NOT added to `QuadReduction`:
    // the corpus is already mixed-version v1/v2 on PRNQ and a third version for an audit's
    // convenience is not worth the hazard. Aggregated here, from the footprints the build
    // already returns.
    //
    // **The aggregation is swept for these two specifically**, because §7.13's whole diffusion
    // verdict turned on a `max`. The best aggregation is reported with its name.
    let agg_of = |f: &dyn Fn(&PixelOut) -> f64, p: f64| -> Vec<f64> {
        keys.iter()
            .map(|k| {
                let mut v: Vec<f64> =
                    px_of[k].iter().map(f).filter(|x: &f64| x.is_finite()).collect();
                q(&mut v, p)
            })
            .collect()
    };
    let mean_of = |f: &dyn Fn(&PixelOut) -> f64| -> Vec<f64> {
        keys.iter()
            .map(|k| {
                let v: Vec<f64> =
                    px_of[k].iter().map(f).filter(|x: &f64| x.is_finite()).collect();
                if v.is_empty() {
                    f64::NAN
                } else {
                    v.iter().sum::<f64>() / v.len() as f64
                }
            })
            .collect()
    };
    let nan_frac = |f: &dyn Fn(&PixelOut) -> f64| -> Vec<f64> {
        keys.iter()
            .map(|k| {
                let n = px_of[k].len().max(1) as f64;
                px_of[k].iter().map(f).filter(|x: &f64| !x.is_finite()).count() as f64 / n
            })
            .collect()
    };

    let diff = |p: &PixelOut| p.diffusion;
    let ftle = |p: &PixelOut| p.ftle;

    // The region-wide median of each, so `frac_above` is a real signal rather than a restatement
    // of the quad's own median. A quantile rule keyed on the quad itself makes the count a fact
    // about the rule and not about the field -- measured, and the reason both masks are kept.
    let region_median = |f: &dyn Fn(&PixelOut) -> f64| -> f64 {
        let mut all: Vec<f64> =
            keys.iter().flat_map(|k| px_of[k].iter().map(f)).filter(|x| x.is_finite()).collect();
        q(&mut all, 0.5)
    };
    let frac_above = |f: &dyn Fn(&PixelOut) -> f64, thr: f64| -> Vec<f64> {
        keys.iter()
            .map(|k| {
                let v: Vec<f64> =
                    px_of[k].iter().map(f).filter(|x: &f64| x.is_finite()).collect();
                if v.is_empty() {
                    f64::NAN
                } else {
                    v.iter().filter(|x| **x > thr).count() as f64 / v.len() as f64
                }
            })
            .collect()
    };
    let d_thr = region_median(&diff);
    let f_thr = region_median(&ftle);

    let red = |f: &dyn Fn(&prin_rs::quad::QuadReduction) -> f64| -> Vec<f64> {
        keys.iter().map(|k| f(&cache.get(*k).red)).collect()
    };

    // `alpha` is `None` on every quad of a metric cache -- the scheduler fills it, this path does
    // not. It is still *computable*, because in a complete tree every quad's parent is present:
    // `alpha = log2(spread_parent / spread_child)`. Named `_from_parent` so it is never read as
    // the scheduler's own field.
    let sig = |k: Key| cache.get(k).red.signal(Criterion::Within, Agg::Median);
    let parent = |k: Key| -> Option<Key> {
        if k.0 == 0 { None } else { Some((k.0 - 1, k.1 / 2, k.2 / 2)) }
    };
    let alpha_of = |k: Key| -> f64 {
        match parent(k) {
            None => f64::NAN,
            Some(p) => (sig(p) / sig(k)).log2(),
        }
    };
    let alpha_from_parent: Vec<f64> = keys.iter().map(|k| alpha_of(*k)).collect();
    let alpha_sibling_range: Vec<f64> = keys
        .iter()
        .map(|k| {
            if k.0 >= cache.levels {
                return f64::NAN;
            }
            let a: Vec<f64> = Cache::children(*k).iter().map(|c| alpha_of(*c)).collect();
            let lo = a.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = a.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            hi - lo
        })
        .collect();

    let mut c: Vec<Col> = Vec::new();
    let mut push = |name: &'static str, v: Vec<f64>| c.push(Col { name, v });

    push("spread_mean", red(&|r| r.spread_mean));
    push("spread_median", red(&|r| r.spread_median));
    push("spread_p90", red(&|r| r.spread_p90));
    push("spread_shape_median", red(&|r| r.spread_shape_median));
    push("spread_event_median", red(&|r| r.spread_event_median));

    push("between_shape", red(&|r| r.between_shape));
    push("between_event", red(&|r| r.between_event));
    push("between_spread", red(&|r| r.between_spread));
    push("between_matched", red(&|r| r.between_matched));
    push("within_pooled", red(&|r| r.within_pooled));

    push("lay_w_n_hot", red(&|r| r.layout_within.n_hot as f64));
    push("lay_w_n_components", red(&|r| r.layout_within.n_components as f64));
    push("lay_w_largest", red(&|r| r.layout_within.largest_component as f64));
    push("lay_w_perimeter", red(&|r| r.layout_within.perimeter_ratio));
    push("lay_b_n_hot", red(&|r| r.layout_between.n_hot as f64));
    push("lay_b_n_components", red(&|r| r.layout_between.n_components as f64));
    push("lay_b_largest", red(&|r| r.layout_between.largest_component as f64));
    push("lay_b_perimeter", red(&|r| r.layout_between.perimeter_ratio));

    push("layrel_w_n_hot", red(&|r| r.layout_rel_within.n_hot as f64));
    push("layrel_w_n_components", red(&|r| r.layout_rel_within.n_components as f64));
    push("layrel_w_largest", red(&|r| r.layout_rel_within.largest_component as f64));
    push("layrel_w_perimeter", red(&|r| r.layout_rel_within.perimeter_ratio));
    push("layrel_b_n_hot", red(&|r| r.layout_rel_between.n_hot as f64));
    push("layrel_b_n_components", red(&|r| r.layout_rel_between.n_components as f64));
    push("layrel_b_largest", red(&|r| r.layout_rel_between.largest_component as f64));
    push("layrel_b_perimeter", red(&|r| r.layout_rel_between.perimeter_ratio));

    push("frac_above_tau_within", red(&|r| r.frac_above_tau_within));
    push("frac_above_tau_between", red(&|r| r.frac_above_tau_between));
    push("grad_rms_within", red(&|r| r.grad_rms_within));
    push("grad_rms_between", red(&|r| r.grad_rms_between));

    push("running_max_divergence", red(&|r| r.running_max_divergence_median));
    push("divergence_trend", red(&|r| r.divergence_trend_median));
    push("first_divergence", red(&|r| r.first_divergence_median));
    push("frac_diverged", red(&|r| r.frac_diverged));

    push("terminated_fraction", red(&|r| r.terminated_fraction));
    push("escape_fraction", red(&|r| r.escape_fraction));
    push("t_end_gradient", red(&|r| r.t_end_gradient));

    push("error_ratio_max", red(&|r| r.error_ratio_max));
    push("worst_energy_drift", red(&|r| r.worst_energy_drift));
    push("n_nonfinite", red(&|r| r.n_nonfinite as f64));
    push("total_substeps", red(&|r| r.total_substeps as f64));

    push("level", keys.iter().map(|k| k.0 as f64).collect());
    push(
        "cell_width",
        keys.iter().map(|k| 2.0 * cache.half / (1u64 << k.0) as f64).collect(),
    );
    push("alpha_from_parent", alpha_from_parent);
    push("alpha_sibling_range", alpha_sibling_range);

    // ---- the two that have never been scored ----
    push("diffusion_mean", mean_of(&diff));
    push("diffusion_median", agg_of(&diff, 0.5));
    push("diffusion_p90", agg_of(&diff, 0.9));
    push("diffusion_frac_above", frac_above(&diff, d_thr));
    push("ftle_mean", mean_of(&ftle));
    push("ftle_median", agg_of(&ftle, 0.5));
    push("ftle_p90", agg_of(&ftle, 0.9));
    push("ftle_frac_above", frac_above(&ftle, f_thr));
    push("ftle_nan_frac", nan_frac(&ftle));

    c
}

/// `max` over the four edge-neighbours at the same level of `|v_self - v_neighbour|`.
///
/// Over a complete cache this is better defined than over a live tree: every quad exists at every
/// level, so the neighbour is always read at the same level rather than falling back up-tree.
fn contrast_of(v: &[f64], keys: &[Key], idx: &HashMap<Key, usize>) -> Vec<f64> {
    keys.iter()
        .enumerate()
        .map(|(i, k)| {
            let (l, ix, iy) = *k;
            let w = 1u32 << l;
            let mut best = f64::NAN;
            for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                let (nx, ny) = (ix as i64 + dx, iy as i64 + dy);
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= w as i64 {
                    continue;
                }
                let j = idx[&(l, nx as u32, ny as u32)];
                let d = (v[i] - v[j]).abs();
                if d.is_finite() && (!best.is_finite() || d > best) {
                    best = d;
                }
            }
            best
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// Logistic regression, IRLS with L2
// ---------------------------------------------------------------------------------------------

/// Fit `p(y=1) = sigmoid(x.w + b)` by iteratively reweighted least squares with ridge `lambda`.
///
/// Ridge is not optional here: with ~55 features and a collinear design -- most of them are
/// functions of the same footprint spreads -- the unpenalised normal equations are singular in
/// practice, and a fit that fails to converge would be read as a null.
fn logistic_fit(x: &[Vec<f64>], y: &[f64], lambda: f64, iters: usize) -> Vec<f64> {
    let n = x.len();
    let p = if n == 0 { 0 } else { x[0].len() };
    let mut w = vec![0.0f64; p + 1]; // last is the intercept
    if n == 0 {
        return w;
    }
    for _ in 0..iters {
        // Gauss-Newton step: solve (X' S X + lambda I) d = X' (y - mu)
        let mut h = vec![0.0f64; (p + 1) * (p + 1)];
        let mut g = vec![0.0f64; p + 1];
        for (row, yi) in x.iter().zip(y) {
            let mut z = w[p];
            for (j, v) in row.iter().enumerate() {
                z += w[j] * v;
            }
            let mu = 1.0 / (1.0 + (-z).exp());
            let s = (mu * (1.0 - mu)).max(1e-6);
            let r = yi - mu;
            for a in 0..=p {
                let xa = if a == p { 1.0 } else { row[a] };
                g[a] += r * xa;
                for b in 0..=p {
                    let xb = if b == p { 1.0 } else { row[b] };
                    h[a * (p + 1) + b] += s * xa * xb;
                }
            }
        }
        for a in 0..p {
            h[a * (p + 1) + a] += lambda;
        }
        let d = solve(&mut h, &mut g, p + 1);
        let mut moved = 0.0;
        for a in 0..=p {
            w[a] += d[a];
            moved += d[a].abs();
        }
        if moved < 1e-10 {
            break;
        }
    }
    w
}

/// Gaussian elimination with partial pivoting. Returns zeros on a singular system rather than
/// `NaN`, so a degenerate fold reports "no fit" instead of poisoning every downstream number.
fn solve(a: &mut [f64], b: &mut [f64], n: usize) -> Vec<f64> {
    for c in 0..n {
        let mut piv = c;
        for r in c + 1..n {
            if a[r * n + c].abs() > a[piv * n + c].abs() {
                piv = r;
            }
        }
        if a[piv * n + c].abs() < 1e-12 {
            return vec![0.0; n];
        }
        if piv != c {
            for k in 0..n {
                a.swap(c * n + k, piv * n + k);
            }
            b.swap(c, piv);
        }
        for r in c + 1..n {
            let f = a[r * n + c] / a[c * n + c];
            if f == 0.0 {
                continue;
            }
            for k in c..n {
                a[r * n + k] -= f * a[c * n + k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    x
}

fn predict(w: &[f64], row: &[f64]) -> f64 {
    let p = row.len();
    let mut z = w[p];
    for (j, v) in row.iter().enumerate() {
        z += w[j] * v;
    }
    z
}

/// Rank-based AUC. `NaN` when one class is absent -- **a fold with one class has no AUC, and
/// saying so beats reporting 0.5**, which reads as "no information" rather than "no test".
fn auc(score: &[f64], y: &[f64]) -> f64 {
    let pos = y.iter().filter(|v| **v > 0.5).count();
    let neg = y.len() - pos;
    if pos == 0 || neg == 0 {
        return f64::NAN;
    }
    let r = stats::ranks(score);
    let sum_pos: f64 = r.iter().zip(y).filter(|(_, v)| **v > 0.5).map(|(a, _)| *a).sum();
    (sum_pos - pos as f64 * (pos as f64 + 1.0) / 2.0) / (pos as f64 * neg as f64)
}

fn log_loss(score: &[f64], y: &[f64]) -> f64 {
    if score.is_empty() {
        return f64::NAN;
    }
    let mut s = 0.0;
    for (z, yi) in score.iter().zip(y) {
        let mu = (1.0 / (1.0 + (-z).exp())).clamp(1e-12, 1.0 - 1e-12);
        s += -(yi * mu.ln() + (1.0 - yi) * (1.0 - mu).ln());
    }
    s / score.len() as f64
}

// ---------------------------------------------------------------------------------------------

/// Everything one region contributes to the cross-region stages.
struct Region {
    name: String,
    /// Kept so a fitted score can be put through `replay_scored` on the held-out region --
    /// **the operational half of stage 4.** An AUC says the fit ranks quads; only `error(B)`
    /// says whether that buys a better image, and every other row in this report is an
    /// `error(B)`. ~13 MB per region at `levels = 7`.
    cache: Cache,
    keys: Vec<Key>,
    /// `[column][quad]`, standardisation applied later so the fit set decides the scale.
    cols: Vec<Col>,
    /// `label[budget_index][quad]` — `NaN` where the quad is not in that budget's optimal tree.
    labels: Vec<Vec<f64>>,
}

fn main() {
    let levels: u32 = arg(1, 7);
    let n: usize = arg(2, 8);
    let tau: f64 = arg(3, 1e-4);
    let t_max: f64 = arg(4, 13.0);
    let root: String = std::env::args().nth(5).unwrap_or_else(|| "results".into());
    let only: String = std::env::args().nth(6).unwrap_or_else(|| "all".into());
    let res = (1usize << levels) * n;

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        keep_boundary_shapes: true,
        // **The named gap.** `diffusion` and `ftle` are NaN without this, which is why they have
        // never been scored -- and §1 of the audit exists because diffusion was excluded on a
        // scale artefact and never given the good test. Running without it would repeat the
        // error with better tooling.
        ftle: Some(FtleOpts::default()),
        ..Default::default()
    };

    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;
    let max_splits = (full - 1) / 4;

    // The budget ladder for the labels, as **fractions of the full tree** rather than absolute
    // budgets. A fixed ladder collapses to one rung at small `levels`, and a stability statistic
    // with one rung reads 0 flips and looks like perfect stability -- *a test that cannot fail is
    // indistinguishable from a test that passes*, and the smoke pass at `levels = 3` produced
    // exactly that before this was relative.
    //
    // **Small budgets are excluded on purpose**: the optimal tree at `B = 5` has five nodes, and
    // a Spearman over five points is not a measurement. Population is printed with every row.
    let ladder: Vec<usize> = {
        let mut v: Vec<usize> = [128usize, 32, 8, 2]
            .into_iter()
            .map(|d| {
                let s = (full / d).saturating_sub(1) / 4;
                1 + 4 * s.max(5)
            })
            .filter(|b| *b <= full)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    // A separate, wider ladder for `error(B)`, which has no population problem.
    let err_budgets: Vec<usize> = {
        let mut b = vec![23usize];
        while *b.last().unwrap() * 4 < full {
            b.push(b.last().unwrap() * 4 + 3);
        }
        b.push(full);
        b
    };

    println!(
        "complete tree to level {levels}, N={n}, E+1={}, res {res}^2, tau={tau:e}, \
         t={t_max} n_sync={n_sync}, f64, FTLE dt={:e} ({} leapfrog steps/footprint)\n\
         {full} quads per region, {max_splits} splits at the full tree\n\
         label ladder B = {ladder:?}   (labels are the DP optimum's own split decisions)\n",
        ens.n_extra + 1,
        ens.ftle_dt,
        (t_max / ens.ftle_dt).round() as u64
    );

    let _ = std::fs::create_dir_all(format!("{root}/audit"));

    let mut regions: Vec<Region> = Vec::new();

    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if only != "all" && only != region {
            continue;
        }
        if !["far", "near-field", "deep interior"].contains(&region) {
            continue;
        }

        let t0 = std::time::Instant::now();
        let (mut caches, px_of) = metric::build_multi_with_footprints(
            region,
            cx,
            cy,
            0.05,
            body,
            Chart::BodyPlane,
            levels,
            n,
            res,
            tau,
            &ens,
            &[SHIPPING],
        );
        let cache = caches.pop().unwrap();
        let build_s = t0.elapsed().as_secs_f64();
        let e_root = cache.error_of(&[(0, 0, 0)]);
        println!("=== {region} === built in {build_s:.1}s; error(root)={e_root:.5}");

        let stem = region.replace(' ', "_");
        // **No `.fcache` here, deliberately.** At `levels = 7` it is 1.4M footprints per region,
        // a few hundred MB, and it would buy nothing: there is no `Cache` reader, so a later
        // process cannot rebuild the tree from it anyway. The TSV below is the durable artefact
        // and every stage after this one is a re-read of it.

        let keys: Vec<Key> = {
            let mut v: Vec<Key> = cache.quads.keys().copied().collect();
            v.sort_unstable();
            v
        };
        let idx: HashMap<Key, usize> =
            keys.iter().enumerate().map(|(i, k)| (*k, i)).collect();

        // ---- STAGE 1: the cost axis, measured per level, and then stopped ----
        println!(
            "\n  [1] total_substeps per quad, BY LEVEL. Every quad carries the same N^2 x (E+1)\n\
             \x20     trajectories whatever its level, so a deep quad lying wholly inside a\n\
             \x20     collision zone is all-hot where the root averages the region -- the spread\n\
             \x20     should widen with depth, and a pooled row would hide it."
        );
        println!(
            "{:>7} {:>8} {:>12} {:>12} {:>12} {:>12} {:>10}",
            "level", "quads", "p1", "p50", "p99", "max", "p99/p1"
        );
        for l in 0..=levels {
            let w = 1u32 << l;
            let v: Vec<f64> = (0..w)
                .flat_map(|iy| (0..w).map(move |ix| (l, ix, iy)))
                .map(|k| cache.get(k).red.total_substeps as f64)
                .collect();
            let s = quant(&v, &[0.01, 0.5, 0.99, 1.0]);
            println!(
                "{l:>7} {:>8} {:>12.4e} {:>12.4e} {:>12.4e} {:>12.4e} {:>10.4}",
                v.len(),
                s[0],
                s[1],
                s[2],
                s[3],
                s[2] / s[0]
            );
        }
        {
            let v: Vec<f64> =
                keys.iter().map(|k| cache.get(*k).red.total_substeps as f64).collect();
            let s = quant(&v, &[0.01, 0.5, 0.99, 1.0]);
            println!(
                "{:>7} {:>8} {:>12.4e} {:>12.4e} {:>12.4e} {:>12.4e} {:>10.4}   <- POOLED over \
                 levels; the two framings answer different questions",
                "pooled",
                v.len(),
                s[0],
                s[1],
                s[2],
                s[3],
                s[2] / s[0]
            );
        }

        // ---- the labels ----
        let dp = cache.dp_optimal(max_splits);
        let mut labels: Vec<Vec<f64>> = Vec::new();
        let mut label_maps: Vec<HashMap<Key, bool>> = Vec::new();
        for b in &ladder {
            let m = dp.labels((b.saturating_sub(1) / 4).min(dp.max_splits));
            labels.push(
                keys.iter()
                    .map(|k| match m.get(k) {
                        Some(true) => 1.0,
                        Some(false) => 0.0,
                        None => f64::NAN,
                    })
                    .collect(),
            );
            label_maps.push(m);
        }

        // ---- STAGE 2: label stability across budgets ----
        println!(
            "\n  [2] LABEL STABILITY. A quad not in the optimal tree was never decided, so it is\n\
             \x20     absent rather than labelled `keep`. Churn is over SHARED quads only -- a quad\n\
             \x20     present at one budget and not the other has not changed its decision, and\n\
             \x20     counting it folds the tree's growth into a statistic about its stability."
        );
        println!(
            "{:>8} {:>10} {:>8} {:>8} {:>9} {:>9} {:>9} {:>9}",
            "B", "in tree", "split", "keep", "shared", "flips", "0->1", "1->0"
        );
        for (i, b) in ladder.iter().enumerate() {
            let m = &label_maps[i];
            let split = m.values().filter(|v| **v).count();
            let (mut shared, mut flip, mut up, mut down) = (0usize, 0usize, 0usize, 0usize);
            if i > 0 {
                let prev = &label_maps[i - 1];
                for (k, v) in m {
                    if let Some(pv) = prev.get(k) {
                        shared += 1;
                        if pv != v {
                            flip += 1;
                            if *v {
                                up += 1;
                            } else {
                                down += 1;
                            }
                        }
                    }
                }
            }
            println!(
                "{b:>8} {:>10} {split:>8} {:>8} {shared:>9} {flip:>9} {up:>9} {down:>9}",
                m.len(),
                m.len() - split
            );
        }
        // The per-level view, because the PR #22 result being explained is per level.
        println!("\n      flips between the first and last ladder rungs, by level:");
        println!("{:>7} {:>9} {:>9} {:>9}", "level", "shared", "0->1", "1->0");
        {
            let (a, b) = (&label_maps[0], &label_maps[label_maps.len() - 1]);
            for l in 0..=levels {
                let (mut sh, mut up, mut dn) = (0usize, 0usize, 0usize);
                for (k, v) in b {
                    if k.0 != l {
                        continue;
                    }
                    if let Some(pv) = a.get(k) {
                        sh += 1;
                        if pv != v {
                            if *v {
                                up += 1
                            } else {
                                dn += 1
                            }
                        }
                    }
                }
                if sh > 0 {
                    println!("{l:>7} {sh:>9} {up:>9} {dn:>9}");
                }
            }
        }

        // ---- the columns ----
        let mut cols = base_columns(&cache, &px_of, &keys);
        drop(px_of);

        // Derived: contrast, and per unit cost. Both are in the audit's brief and both are cheap.
        let n_base = cols.len();
        let subs: Vec<f64> =
            keys.iter().map(|k| cache.get(*k).red.total_substeps.max(1) as f64).collect();
        let mut derived: Vec<Col> = Vec::new();
        for c in cols.iter().take(n_base) {
            derived.push(Col {
                name: Box::leak(format!("contrast:{}", c.name).into_boxed_str()),
                v: contrast_of(&c.v, &keys, &idx),
            });
            derived.push(Col {
                name: Box::leak(format!("{}/cost", c.name).into_boxed_str()),
                v: c.v.iter().zip(&subs).map(|(a, b)| a / b).collect(),
            });
        }

        // ---- the table ----
        {
            use std::io::Write;
            let path = format!("{root}/audit/signal_audit_{stem}.tsv");
            if let Ok(f) = std::fs::File::create(&path) {
                let mut w = std::io::BufWriter::new(f);
                let _ = write!(w, "level\tix\tiy\terr_sum\tgain");
                for c in &cols {
                    let _ = write!(w, "\t{}", c.name);
                }
                for b in &ladder {
                    let _ = write!(w, "\tlabel_B{b}");
                }
                let _ = writeln!(w);
                for (i, k) in keys.iter().enumerate() {
                    let _ = write!(
                        w,
                        "{}\t{}\t{}\t{:.9e}\t{:.9e}",
                        k.0,
                        k.1,
                        k.2,
                        cache.get(*k).err_sum,
                        cache.gain(*k)
                    );
                    for c in &cols {
                        let _ = write!(w, "\t{:.9e}", c.v[i]);
                    }
                    for lab in &labels {
                        let _ = write!(w, "\t{}", lab[i]);
                    }
                    let _ = writeln!(w);
                }
                println!("\n      table -> {path} ({} rows x {} signal columns)", keys.len(), cols.len());
            }
        }

        // ---- STAGE 3: Spearman against the label, plus error(B) beside it ----
        println!(
            "\n  [3] SPEARMAN against the DP label, per budget, with error(B) beside it.\n\
             \x20     High rho with poor error(B) says the signal is worth RESCALING; low both\n\
             \x20     says worth DROPPING. Nothing in the corpus could tell those apart.\n\
             \x20     `distinct` and `modal%` decide what a flat error row means -- a bad ordering\n\
             \x20     and NO ordering give the same curve. `nan%` is a property to read, not a\n\
             \x20     defect to hide.\n\
             \x20     `blk` is BLOCKED BY LEVEL and `pool` is not. Read `blk`. `level` itself\n\
             \x20     scores |rho| = 0.99 pooled -- the optimum splits shallow quads -- so the\n\
             \x20     pooled column is largely a picture of depth and every signal tracking cell\n\
             \x20     width inherits it. Rows are sorted by |blk|."
        );
        let mut all: Vec<&Col> = cols.iter().chain(derived.iter()).collect();
        all.sort_by_key(|c| c.name);
        let mut rows: Vec<(f64, String, &str, Vec<f64>)> = Vec::new();
        for c in &all {
            // Spearman over the quads the optimum actually decided, per budget.
            let rhos: Vec<f64> = labels
                .iter()
                .map(|lab| {
                    let (x, y): (Vec<f64>, Vec<f64>) = keys
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| lab[*i].is_finite() && c.v[*i].is_finite())
                        .map(|(i, _)| (c.v[i], lab[i]))
                        .unzip();
                    if x.len() < 8 { f64::NAN } else { stats::spearman(&x, &y) }
                })
                .collect();
            // **BLOCKED BY LEVEL, and this is not optional.** `level` itself reads |rho| = 0.993
            // against the label -- the optimum splits shallow quads and keeps deep ones, so a
            // pooled Spearman is largely a picture of depth and every signal that tracks cell
            // width inherits it. The form with the confound removed asks: among the quads AT
            // level L, did the ones the optimum split have the higher signal? Same repair as
            // `rho(depth, spread)`, which was confounded twice and fixed by blocking.
            //
            // Averaged over levels that have both outcomes; a level with one outcome has no
            // correlation and contributes NaN rather than 0.
            let blocked: Vec<f64> = labels
                .iter()
                .map(|lab| {
                    let mut acc = Vec::new();
                    for l in 0..=cache.levels {
                        let (x, y): (Vec<f64>, Vec<f64>) = keys
                            .iter()
                            .enumerate()
                            .filter(|(i, k)| {
                                k.0 == l && lab[*i].is_finite() && c.v[*i].is_finite()
                            })
                            .map(|(i, _)| (c.v[i], lab[i]))
                            .unzip();
                        if x.len() >= 8 && y.iter().any(|v| *v > 0.5) && y.iter().any(|v| *v < 0.5)
                        {
                            let r = stats::spearman(&x, &y);
                            if r.is_finite() {
                                acc.push(r);
                            }
                        }
                    }
                    if acc.is_empty() {
                        f64::NAN
                    } else {
                        acc.iter().sum::<f64>() / acc.len() as f64
                    }
                })
                .collect();
            let nanp = c.v.iter().filter(|x| !x.is_finite()).count() as f64
                / c.v.len() as f64
                * 100.0;
            let (distinct, modal) = {
                let mut b: Vec<u64> = c.v.iter().map(|x| x.to_bits()).collect();
                b.sort_unstable();
                let mut d = 0usize;
                let mut best = 0usize;
                let mut i = 0;
                while i < b.len() {
                    let mut j = i;
                    while j < b.len() && b[j] == b[i] {
                        j += 1;
                    }
                    d += 1;
                    best = best.max(j - i);
                    i = j;
                }
                (d, best as f64 / b.len() as f64 * 100.0)
            };
            let score: HashMap<Key, f64> =
                keys.iter().enumerate().map(|(i, k)| (*k, c.v[i])).collect();
            let (pts, _) = metric::replay_scored(&cache, &score, full);
            let curve = metric::curve_at(&pts, &err_budgets);
            // **Ranked by the BLOCKED rho**, not the pooled one, because the pooled ordering is
            // an ordering by depth.
            let best_rho = blocked.iter().copied().filter(|x| x.is_finite()).fold(
                0.0f64,
                |a, b| if b.abs() > a.abs() { b } else { a },
            );
            let mut line = format!(
                "{:>26} {distinct:>7} {modal:>7.1} {nanp:>6.1}",
                c.name
            );
            for r in &blocked {
                line += &format!(" {r:>7.3}");
            }
            line += " |";
            for r in &rhos {
                line += &format!(" {r:>7.3}");
            }
            line += "  |";
            for e in &curve {
                line += &format!(" {e:>8.5}");
            }
            rows.push((-best_rho.abs(), line, c.name, curve.clone()));
        }
        // **Sorted by the OPERATIONAL metric, not the scale-free one.** The first cut of this
        // table sorted by |blocked rho| and printed the top 60, which is a table sorted by one
        // statistic while a different statistic decides -- and it hid rows. Every row is printed
        // now, ordered by `error(B)` at the mid rung, with the rho columns beside them.
        let mid_e = err_budgets.len() / 2;
        rows.sort_by(|a, b| {
            let k = |x: f64| if x.is_finite() { x } else { f64::INFINITY };
            k(a.3[mid_e]).partial_cmp(&k(b.3[mid_e])).unwrap()
        });
        print!("{:>26} {:>7} {:>7} {:>6}", "signal", "distinct", "modal%", "nan%");
        for b in &ladder {
            print!(" {:>7}", format!("blk{b}"));
        }
        print!(" |");
        for b in &ladder {
            print!(" {:>7}", format!("pool{b}"));
        }
        print!("  |");
        for b in &err_budgets {
            print!(" {b:>8}");
        }
        println!();
        for r in rows.iter() {
            println!("{}", r.1);
        }
        println!("      ({} rows, sorted by error(B={}))", rows.len(), err_budgets[mid_e]);

        // The references, on the same error ladder.
        for r in [Rank::Uniform, Rank::GreedyLookahead1] {
            let pts = metric::replay(&cache, r, full);
            let curve = metric::curve_at(&pts, &err_budgets);
            print!("{:>26} {:>7} {:>7} {:>6}", r.name(), "-", "-", "-");
            for _ in &ladder {
                print!(" {:>7}", "-");
            }
            print!(" |");
            for _ in &ladder {
                print!(" {:>7}", "-");
            }
            print!("  |");
            for e in &curve {
                print!(" {e:>8.5}");
            }
            println!();
        }
        print!("{:>26} {:>7} {:>7} {:>6}", "dp_optimal (CEILING)", "-", "-", "-");
        for _ in &ladder {
            print!(" {:>7}", "-");
        }
        print!(" |");
        for _ in &ladder {
            print!(" {:>7}", "-");
        }
        print!("  |");
        for b in &err_budgets {
            print!(" {:>8.5}", dp.at_budget(*b));
        }
        println!();

        // ---- WHO BEATS BREADTH-FIRST, counted rather than eyeballed ----
        //
        // This is the whole question. PR #22 found no criterion clearing `uniform` in
        // `near-field` at any budget; that was over the eleven `Criterion` variants, and this
        // table is 162 signals wide. A row beating uniform is the finding; a table where the
        // reader has to spot it is not a report.
        {
            let upts = metric::replay(&cache, Rank::Uniform, full);
            let uni = metric::curve_at(&upts, &err_budgets);
            println!(
                "
      SIGNALS THAT BEAT BREADTH-FIRST, per budget. `captured` is
                       (uniform - row)/(uniform - dp): 1.0 is the exact optimum, 0.0 is no
                       better than refining uniformly."
            );
            println!(
                "{:>10} {:>10} {:>10} {:>8} {:>28} {:>10} {:>10}",
                "B", "uniform", "dp", "n better", "best signal", "error", "captured"
            );
            for (j, b) in err_budgets.iter().enumerate() {
                let d = dp.at_budget(*b);
                let better = rows.iter().filter(|r| r.3[j] < uni[j] - 1e-12).count();
                let best = rows
                    .iter()
                    .filter(|r| r.3[j].is_finite())
                    .min_by(|x, y| x.3[j].partial_cmp(&y.3[j]).unwrap());
                match best {
                    None => println!("{b:>10} {:>10.5} {d:>10.5} {better:>8}", uni[j]),
                    Some(r) => {
                        let den = uni[j] - d;
                        let cap = if den.abs() <= 1e-12 * uni[j].abs().max(f64::MIN_POSITIVE) {
                            "uni==dp".to_string()
                        } else {
                            format!("{:.4}", (uni[j] - r.3[j]) / den)
                        };
                        println!(
                            "{b:>10} {:>10.5} {d:>10.5} {better:>8} {:>28} {:>10.5} {cap:>10}",
                            uni[j], r.2, r.3[j]
                        );
                    }
                }
            }
        }

        // The FTLE caveat, printed where it is read rather than buried in a footer.
        println!(
            "\n      NOTE: the FTLE march is the UNREGULARISED fixed-step leapfrog, so near a\n\
             \x20     close approach it is not trustworthy. An ftle_* result in `deep interior`\n\
             \x20     carries that; a spread_* result does not."
        );

        regions.push(Region {
            name: region.into(),
            cache,
            keys,
            cols: std::mem::take(&mut cols),
            labels,
        });
        println!();
    }

    // ---- STAGE 4: does more information help ----
    if regions.len() > 1 {
        stage4(&regions, &ladder, &err_budgets, full);
    } else {
        println!("[4] skipped: the held-out designs need more than one region.");
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         [1] is a MEASUREMENT AND A STOP. If p99/p1 sits near 1.0 at every level, cost-aware\n\
         priority is arithmetically incapable of reordering anything at this configuration and\n\
         the cost axis is written off with the numbers rather than argued about. `far` being\n\
         exactly constant is the mechanism showing: smooth region, no encounters, every\n\
         trajectory the same work.\n\n\
         [2] conditions everything after it. If the optimum's label for a quad FLIPS between\n\
         budgets, then no budget-independent signal can be optimal at both, and the per-level\n\
         `captured` result is a symptom of that rather than of a bad signal. Read the shared\n\
         count before the churn.\n\n\
         [3] separates `does not know` from `knows but is drowned out`. Spearman is scale-free;\n\
         error(B) is operational. A signal high on the first and poor on the second is being\n\
         wasted by the aggregation, which is exactly the failure that excluded diffusion.\n\n\
         [4] says whether the signals are REDUNDANT. Read the held-out columns against the\n\
         in-sample one: 55 features on tens of thousands of quads will fit near-perfectly and\n\
         prove nothing. And three heterogeneous regions is NOT a validation set -- the regions\n\
         differ in kind, so the folds are never averaged and a null from them is weak evidence\n\
         either way."
    );
}

/// Leave-one-region-out and low-budget-to-high-budget, both, never averaged.
fn stage4(regions: &[Region], ladder: &[usize], err_budgets: &[usize], full: usize) {
    let names: Vec<&str> = regions[0].cols.iter().map(|c| c.name).collect();
    // **Two feature sets, because `level` alone nearly solves the task.** The optimum splits
    // shallow quads and keeps deep ones, so a fit with `level` and `cell_width` in it is largely a
    // depth model wearing 55 names, and its held-out AUC would be read as "the signals carry
    // information" when what carries it is the tree geometry. The second set removes them and asks
    // the question the audit was commissioned for: does what we MEASURE about a quad add anything
    // beyond where it sits.
    let geom: Vec<usize> = names
        .iter()
        .enumerate()
        .filter(|(_, n)| **n == "level" || **n == "cell_width")
        .map(|(i, _)| i)
        .collect();
    for (label, drop_geom) in [("ALL signals", false), ("WITHOUT level/cell_width", true)] {
        let use_j: Vec<usize> = (0..names.len())
            .filter(|j| !(drop_geom && geom.contains(j)))
            .collect();
        stage4_one(regions, ladder, &names, &use_j, label, err_budgets, full);
    }
}

fn stage4_one(
    regions: &[Region],
    ladder: &[usize],
    names: &[&str],
    use_j: &[usize],
    set_label: &str,
    err_budgets: &[usize],
    full: usize,
) {
    let p = use_j.len();

    // Standardise on the FIT set only, and impute a non-finite feature at the fit set's median.
    // Imputation is stated rather than hidden: a `NaN` column is a measurement outcome, and
    // dropping its rows would bias the sample toward the quads every signal happened to score.
    let design = |r: &Region, bi: usize| -> (Vec<Vec<f64>>, Vec<f64>) {
        let lab = &r.labels[bi];
        let mut x = Vec::new();
        let mut y = Vec::new();
        for i in 0..r.keys.len() {
            if !lab[i].is_finite() {
                continue;
            }
            x.push(use_j.iter().map(|j| r.cols[*j].v[i]).collect::<Vec<f64>>());
            y.push(lab[i]);
        }
        (x, y)
    };
    // Returns the fit set's `(mean, sd, median)` per column, so a later score of the WHOLE
    // held-out region uses the same transform. Recomputing it from the scored rows would let the
    // held-out data set its own scale, which is leakage wearing the name of standardisation.
    let standardise = |fit: &mut Vec<Vec<f64>>, score: &mut Vec<Vec<f64>>| -> Vec<(f64, f64, f64)> {
        let mut stat = Vec::with_capacity(p);
        for j in 0..p {
            let mut v: Vec<f64> =
                fit.iter().map(|r| r[j]).filter(|x| x.is_finite()).collect();
            let med = if v.is_empty() { 0.0 } else { stats::quantile(&mut v, 0.5) };
            let mut w: Vec<f64> = fit.iter().map(|r| r[j]).map(|x| if x.is_finite() { x } else { med }).collect();
            let n = w.len().max(1) as f64;
            let mu = w.iter().sum::<f64>() / n;
            let sd = (w.iter().map(|x| (x - mu) * (x - mu)).sum::<f64>() / n).sqrt();
            let sd = if sd > 0.0 && sd.is_finite() { sd } else { 1.0 };
            for (r, x) in fit.iter_mut().zip(w.drain(..)) {
                r[j] = (x - mu) / sd;
            }
            for r in score.iter_mut() {
                let x = if r[j].is_finite() { r[j] } else { med };
                r[j] = (x - mu) / sd;
            }
            stat.push((mu, sd, med));
        }
        stat
    };

    let bi = ladder.len() / 2; // the mid rung, where the population is largest and the anomaly reads
    println!(
        "[4] {set_label} -- logistic on {p} standardised signals against the DP\n\
    \x20   label at B = {}. Ridge L2 = 1.0. Non-finite features imputed at the FIT set's median,\n\
    \x20   never by dropping rows: a NaN is a measurement outcome, and dropping biases the sample\n\
    \x20   toward the quads every signal happened to score.\n",
        ladder[bi]
    );

    println!("  (a) LEAVE-ONE-REGION-OUT, three folds, REPORTED SEPARATELY. The regions differ in\n\
    \x20     KIND -- smooth / localised / everywhere -- so a mean over them is a mean over three\n\
    \x20     different questions.");
    println!(
        "{:>16} {:>9} {:>9} {:>11} {:>11} {:>11}",
        "held out", "n fit", "n score", "AUC(held)", "AUC(in)", "logloss(h)"
    );
    for (h, hr) in regions.iter().enumerate() {
        let mut xf: Vec<Vec<f64>> = Vec::new();
        let mut yf: Vec<f64> = Vec::new();
        for (i, r) in regions.iter().enumerate() {
            if i == h {
                continue;
            }
            let (x, y) = design(r, bi);
            xf.extend(x);
            yf.extend(y);
        }
        let (mut xs, ys) = design(hr, bi);
        let mut xf2 = xf.clone();
        let stat = standardise(&mut xf2, &mut xs);
        let w = logistic_fit(&xf2, &yf, 1.0, 25);
        let sh: Vec<f64> = xs.iter().map(|r| predict(&w, r)).collect();
        let si: Vec<f64> = xf2.iter().map(|r| predict(&w, r)).collect();
        println!(
            "{:>16} {:>9} {:>9} {:>11.4} {:>11.4} {:>11.4}",
            hr.name,
            xf2.len(),
            xs.len(),
            auc(&sh, &ys),
            auc(&si, &yf),
            log_loss(&sh, &ys)
        );
        // **The operational half.** An AUC says the fit orders the labelled quads; error(B) says
        // whether that orders the IMAGE better, which is what every other row in this report is
        // measured by. Scored on the held-out region only, and the score is computed for EVERY
        // quad -- not only the labelled ones -- because a replay descends through quads the
        // optimal tree never reached.
        let scores: HashMap<Key, f64> = hr
            .keys
            .iter()
            .enumerate()
            .map(|(i, k)| {
                let mut row: Vec<f64> =
                    use_j.iter().map(|j| hr.cols[*j].v[i]).collect();
                for (jj, v) in row.iter_mut().enumerate() {
                    let (mu, sd, med) = stat[jj];
                    let x = if v.is_finite() { *v } else { med };
                    *v = (x - mu) / sd;
                }
                (*k, predict(&w, &row))
            })
            .collect();
        let (pts, _) = metric::replay_scored(&hr.cache, &scores, full);
        let curve = metric::curve_at(&pts, err_budgets);
        print!("{:>16} {:>31}", "", "error(B) of the fit:");
        for e in &curve {
            print!(" {e:>8.5}");
        }
        println!();
        // The rows it has to beat, on the same ladder, in the same region.
        for r in [Rank::Uniform, Rank::Signal(Criterion::FracHotBetween, Agg::Median)] {
            let pts = metric::replay(&hr.cache, r, full);
            let c = metric::curve_at(&pts, err_budgets);
            print!("{:>16} {:>31}", "", format!("{}:", r.name()));
            for e in &c {
                print!(" {e:>8.5}");
            }
            println!();
        }
    }

    println!(
        "\n  (b) BUDGET HOLD-OUT -- fit on the low rungs, score on the high ones. This is the\n\
    \x20     better of the two: it is stage [2]'s question asked of the fit. A fit that survives\n\
    \x20     low->high transfer says the signal set is budget-stable; one that does not localises\n\
    \x20     the failure exactly."
    );
    println!(
        "{:>16} {:>10} {:>10} {:>9} {:>9} {:>11} {:>11}",
        "region", "fit B", "score B", "n fit", "n score", "AUC(held)", "AUC(in)"
    );
    let (lo, hi) = (0usize, ladder.len() - 1);
    if lo == hi {
        println!(
            "      SKIPPED: the ladder has one rung, so `fit low, score high` would be fitting\n\
        \x20     and scoring the same labels and would report a perfect transfer it never tested."
        );
    }
    for r in regions.iter().filter(|_| lo != hi) {
        let (mut xf, yf) = design(r, lo);
        let (mut xs, ys) = design(r, hi);
        let _ = standardise(&mut xf, &mut xs);
        let w = logistic_fit(&xf, &yf, 1.0, 25);
        let sh: Vec<f64> = xs.iter().map(|q| predict(&w, q)).collect();
        let si: Vec<f64> = xf.iter().map(|q| predict(&w, q)).collect();
        println!(
            "{:>16} {:>10} {:>10} {:>9} {:>9} {:>11.4} {:>11.4}",
            r.name,
            ladder[lo],
            ladder[hi],
            xf.len(),
            xs.len(),
            auc(&sh, &ys),
            auc(&si, &yf)
        );
    }

    // ---- collinearity, which can explain a null before it is mistaken for a ceiling ----
    println!(
        "\n  (c) COLLINEARITY. Most of these are functions of the same footprint spreads, so if\n\
    \x20     everything is collinear that explains a null BEFORE it is read as `the criterion\n\
    \x20     cannot be improved`. Multiple R^2 of each signal on the other {} , pooled over\n\
    \x20     regions at B = {}:",
        p - 1,
        ladder[bi]
    );
    {
        let mut x: Vec<Vec<f64>> = Vec::new();
        for r in regions {
            let (xi, _) = design(r, bi);
            x.extend(xi);
        }
        let mut dummy: Vec<Vec<f64>> = Vec::new();
        let mut xs = x.clone();
        let _ = standardise(&mut xs, &mut dummy);
        let mut r2: Vec<(f64, &str)> = Vec::new();
        for j in 0..p {
            let y: Vec<f64> = xs.iter().map(|r| r[j]).collect();
            let others: Vec<Vec<f64>> = xs
                .iter()
                .map(|r| (0..p).filter(|k| *k != j).map(|k| r[k]).collect())
                .collect();
            r2.push((multiple_r2(&others, &y), names[use_j[j]]));
        }
        // NaN-safe: a structurally undefined R^2 sorts last rather than panicking. A level with
        // one outcome has no correlation, and saying so is the point.
        let key = |x: f64| if x.is_finite() { x } else { f64::NEG_INFINITY };
        r2.sort_by(|a, b| key(b.0).partial_cmp(&key(a.0)).unwrap());
        let vals: Vec<f64> = r2.iter().map(|v| v.0).filter(|v| v.is_finite()).collect();
        println!(
            "      median R^2 = {:.4}; {} of {p} above 0.99",
            q(&mut vals.clone(), 0.5),
            vals.iter().filter(|v| **v > 0.99).count()
        );
        println!("      most collinear:  {}", fmt5(&r2[..5.min(r2.len())]));
        println!("      least collinear: {}", fmt5(&r2[r2.len().saturating_sub(5)..]));
    }
}

fn fmt5(v: &[(f64, &str)]) -> String {
    v.iter().map(|(r, n)| format!("{n} {r:.4}")).collect::<Vec<_>>().join("  ")
}

/// `R^2` of `y` regressed on `x`, by ridge-stabilised normal equations.
fn multiple_r2(x: &[Vec<f64>], y: &[f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    let p = x[0].len();
    let mut a = vec![0.0f64; p * p];
    let mut b = vec![0.0f64; p];
    for (row, yi) in x.iter().zip(y) {
        for i in 0..p {
            b[i] += row[i] * yi;
            for j in 0..p {
                a[i * p + j] += row[i] * row[j];
            }
        }
    }
    for i in 0..p {
        a[i * p + i] += 1e-6 * n as f64;
    }
    let w = solve(&mut a, &mut b, p);
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    let mu = y.iter().sum::<f64>() / n as f64;
    for (row, yi) in x.iter().zip(y) {
        let pred: f64 = row.iter().zip(&w).map(|(a, b)| a * b).sum();
        ss_res += (yi - pred) * (yi - pred);
        ss_tot += (yi - mu) * (yi - mu);
    }
    if ss_tot <= 0.0 { f64::NAN } else { 1.0 - ss_res / ss_tot }
}
