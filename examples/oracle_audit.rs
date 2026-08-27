//! Clear `error(B)` before any result read from it is trusted.
//!
//! # The anomaly
//!
//! `results/output/criterion_metric.txt`, `far`, `B = 1535`: `greedy_lookahead_1` reads **0.54760**
//! against a random band of **0.48550-0.52047** and every criterion at **0.36557**. Greedy worse
//! than random. Myopia puts greedy between random and optimal; it does not put it below random.
//!
//! # What was proposed to catch it, and why neither could
//!
//! The accounting identity `error_of(leaves) == err_sum(root) - sum(gains)` telescopes directly
//! from the definitions of [`prin_rs::metric::Cache::error_of`] (a sum of `err_sum` over leaves)
//! and `Cache::gain` (parent minus children). It holds for any ranking, any sequence, and any
//! values `err_sum` happens to hold -- random numbers included. And a choice check re-runs
//! `replay_with_leaves`'s own argmax over a pure, static `gain`. **Both report PASS and would have
//! been read as clearing the metric.**
//!
//! # The bound that can fail
//!
//! [`prin_rs::metric::Cache::dp_optimal`] is the exact minimum over **all** tree-shaped leaf sets
//! at a given budget. *No ranking may beat it.* If one does, the harness is wrong and every
//! `error(B)` number in the corpus is suspect. If none does, greedy's collapse is genuine myopia.
//!
//! **The prediction, written down before the run:** `far` is smooth, a quad's spread over a
//! gradient `g` is `~ g*w` and so tracks cell width, and argmax-on-spread therefore *is*
//! breadth-first -- which on a smooth field is near optimal. So `dp` should land near **0.366**.
//! Near **0.548**, at greedy, and something is wrong.
//!
//! # Reading the far block
//!
//! `far`'s thirteen non-greedy rows are **one scan order wearing thirteen names**, by two routes:
//! the non-constant signals because spread tracks cell size, the constant ones because they fall
//! through to the level-first tie-break. The distinct-value count and the leaf-level histogram are
//! printed beside every row so this cannot be read as thirteen criteria independently agreeing.
//!
//! And the consequence generalises: **a criterion can only differ from breadth-first where the
//! field is not smooth**, and smoothness is where there is nothing to find. `far` degenerating is
//! correct behaviour, not a defect -- it is the control that shows what a featureless field looks
//! like.
//!
//! # Writes
//!
//! stdout only. **No validation run writes into `results/`** -- a `criterion_metric` pass once
//! overwrote committed 512^2 artefacts with 128x64 ones, and a small raster reads as a rendering
//! fault rather than a stale file.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Key, Rank};
use prin_rs::quad::{Agg, Criterion};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const SHIPPING: metric::Colouring =
    metric::Colouring::Bivariate(prin_rs::output::colour::Scalar::ShapeSpread);

/// Pearson `r` over paired samples. `NaN` on fewer than two pairs or a degenerate arm -- a
/// constant field has no correlation, and saying so beats reporting 0.
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 2 {
        return f64::NAN;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for (x, y) in a.iter().zip(b) {
        sxy += (x - ma) * (y - mb);
        sxx += (x - ma) * (x - ma);
        syy += (y - mb) * (y - mb);
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return f64::NAN;
    }
    sxy / (sxx * syy).sqrt()
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { f64::NAN } else { prin_rs::stats::quantile(v, p) }
}

/// `1 + 3s` leaves binned by level, as a compact `l:count` string.
fn hist(leaves: &[Key], levels: u32) -> String {
    let mut h = vec![0usize; levels as usize + 1];
    for k in leaves {
        h[k.0 as usize] += 1;
    }
    h.iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(l, c)| format!("{l}:{c}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    let levels: u32 = arg(1, 7);
    let n: usize = arg(2, 8);
    let tau: f64 = arg(3, 1e-4);
    let t_max: f64 = arg(4, 13.0);
    let only: String = std::env::args().nth(5).unwrap_or_else(|| "all".into());
    let res = (1usize << levels) * n;

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        // The numpy reference's ungated escape test, with escape terminal: every result in
        // this diagnostic predates both the distance gate and the closure criterion, and is
        // quoted against that form.
        escape_rule: prin_rs::outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_escape: true,
        keep_boundary_shapes: true,
        ..Default::default()
    };
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;
    let max_splits = (full - 1) / 4;

    let budgets: Vec<usize> = {
        let mut b = vec![5usize];
        while *b.last().unwrap() * 2 < full {
            b.push(b.last().unwrap() * 2 + 1);
        }
        b.push(full);
        b
    };

    println!(
        "complete tree to level {levels}, N={n}, E+1={}, res {res}^2, tau={tau:e}, \
         t={t_max} n_sync={n_sync}, f64\n\
         {full} quads per region, {max_splits} splits at the full tree\n\
         floor = random band | reference = greedy_lookahead_1 | ceiling = dp_optimal\n",
        ens.n_extra + 1
    );

    // **`far` first**: it is where the anomaly lives and the cheapest region to build.
    let order = ["far", "near-field", "deep interior"];
    for name in order {
        if only != "all" && only != name {
            continue;
        }
        let &(region, cx, cy, body) =
            grid::REGIONS.iter().find(|r| r.0 == name).expect("region");

        let t0 = std::time::Instant::now();
        let cache = metric::build(
            region, cx, cy, 0.05, body, Chart::BodyPlane, levels, n, res, tau, &ens, SHIPPING,
        );
        let build_s = t0.elapsed().as_secs_f64();

        let e_root = cache.error_of(&[(0, 0, 0)]);
        println!("=== {region} === built in {build_s:.1}s; error(root)={e_root:.5}");

        // ---- the ceiling ----
        let dp = cache.dp_optimal(max_splits);
        println!(
            "  dp_optimal: {} splits (B up to {}) in {:.2}s -- MEASURED, not the estimate. \
             The 4-way merge is three 2-way convolutions, O(cap^2) not O(cap^4), and each node's \
             cap is bounded by its own subtree; that is what makes it affordable.",
            dp.max_splits,
            1 + 4 * dp.max_splits,
            dp.elapsed_s
        );
        // Where the prefix-min binds, a split made the image WORSE -- the measurable consequence
        // of a parent's N x N grid and its children's being different approximation families.
        if dp.prefix_min_binds.is_empty() {
            println!(
                "  prefix-min binds NOWHERE: f_root(s) is monotone, so no split of the optimal \
                 tree ever increased the error. Gains are non-negative along the optimal path."
            );
        } else {
            let b = &dp.prefix_min_binds;
            println!(
                "  prefix-min binds at {} of {} split counts (first {}, last {}): spending more \
                 made the image worse there. NEGATIVE GAIN, measured.",
                b.len(),
                dp.max_splits + 1,
                b[0],
                b[b.len() - 1]
            );
        }

        // ---- the curves ----
        let runs: Vec<Rank> = vec![
            Rank::Uniform,
            Rank::GreedyLookahead1,
            Rank::GreedyLookahead1PerCost,
            Rank::Signal(Criterion::Within, Agg::Median),
            Rank::Signal(Criterion::Between, Agg::Median),
            Rank::Signal(Criterion::FracHotWithin, Agg::Median),
            Rank::Signal(Criterion::FracHotBetween, Agg::Median),
            Rank::Signal(Criterion::Layout, Agg::Median),
            Rank::Signal(Criterion::GradRms, Agg::Median),
            Rank::Signal(Criterion::TerminationGradient, Agg::Median),
            Rank::Signal(Criterion::RunningMax, Agg::Median),
            Rank::Signal(Criterion::FirstDivergence, Agg::Median),
            Rank::Random(1),
            Rank::Random(2),
            Rank::Random(3),
            Rank::Random(4),
            Rank::Random(5),
        ];
        let mut rows: Vec<(String, Vec<f64>, usize)> = Vec::new();
        for r in &runs {
            let pts = metric::replay(&cache, *r, full);
            // Distinct values of the signal itself: a flat curve has two causes -- a bad
            // ordering and NO ordering -- and error(B) alone cannot tell them apart.
            let distinct = match r {
                Rank::Random(_) => 0,
                _ => {
                    let mut bits: Vec<u64> = cache
                        .quads
                        .keys()
                        .map(|k| metric::score(&cache, *k, *r).to_bits())
                        .collect();
                    bits.sort_unstable();
                    bits.dedup();
                    bits.len()
                }
            };
            rows.push((r.name(), metric::curve_at(&pts, &budgets), distinct));
        }

        println!();
        print!("{:>24} {:>8}", "B =", "distinct");
        for b in &budgets {
            print!(" {b:>9}");
        }
        println!();

        // The ceiling, first, so every row below is read against it.
        print!("{:>24} {:>8}", "dp_optimal (CEILING)", "-");
        for b in &budgets {
            print!(" {:>9.5}", dp.at_budget(*b));
        }
        println!();

        for (nm, curve, distinct) in &rows {
            if nm.starts_with("random") {
                continue;
            }
            print!("{nm:>24} {distinct:>8}");
            for e in curve {
                print!(" {e:>9.5}");
            }
            println!();
        }
        let rnd: Vec<&Vec<f64>> =
            rows.iter().filter(|r| r.0.starts_with("random")).map(|r| &r.1).collect();
        for (lbl, f) in [("random lo", f64::min as fn(f64, f64) -> f64), ("random hi", f64::max)] {
            print!("{lbl:>24} {:>8}", "-");
            for j in 0..budgets.len() {
                let v = rnd.iter().skip(1).fold(rnd[0][j], |a, r| f(a, r[j]));
                print!(" {v:>9.5}");
            }
            println!();
        }

        // ---- THE BOUND, checked and reported as a margin ----
        //
        // A rounded table cell cannot distinguish "equal" from "below by 1e-6", and below is a
        // harness bug. Report the worst signed margin `row - dp` in absolute units, and which row
        // and budget it fell at.
        {
            let mut worst = (f64::INFINITY, String::new(), 0usize);
            for (nm, curve, _) in &rows {
                for (&b, &e) in budgets.iter().zip(curve) {
                    let m = e - dp.at_budget(b);
                    if m < worst.0 {
                        worst = (m, nm.clone(), b);
                    }
                }
            }
            println!(
                "  bound check: worst margin (row - dp) = {:+.3e} at {} B={}  --  {}",
                worst.0,
                worst.1,
                worst.2,
                if worst.0 >= -1e-12 {
                    "HOLDS. No ranking beats the exact optimum."
                } else {
                    "**VIOLATED**. A ranking beat the exact optimum, so the harness is wrong and \
                     every error(B) number in the corpus is suspect."
                }
            );
        }

        // ---- the gap that is the value of lookahead ----
        //
        // `(row - dp) / (root - dp)`: the share of ACHIEVABLE improvement a ranking leaves on the
        // table. Against `root - row` it would flatter every row equally; against the ceiling it
        // says whether there is headroom worth chasing.
        println!(
            "\n  share of achievable improvement left on the table, (row - dp)/(root - dp):"
        );
        print!("{:>24}", "B =");
        for b in &budgets {
            print!(" {b:>9}");
        }
        println!();
        let mut best: Option<(String, Vec<f64>)> = None;
        for (nm, curve, _) in &rows {
            if nm.starts_with("random") {
                continue;
            }
            let g: Vec<f64> = budgets
                .iter()
                .zip(curve)
                .map(|(b, e)| {
                    let d = dp.at_budget(*b);
                    let denom = e_root - d;
                    if denom.abs() < 1e-15 { f64::NAN } else { (e - d) / denom }
                })
                .collect();
            if nm.starts_with("greedy_lookahead_1") && !nm.contains("cost") {
                print!("{nm:>24}");
                for v in &g {
                    print!(" {v:>9.4}");
                }
                println!();
            }
            // "best" by the mid-ladder point, which is where the anomaly reads largest.
            // **Greedy is excluded**: it is the reference being measured against, and letting it
            // win its own comparison would print the reference twice under two labels.
            if nm.starts_with("greedy_lookahead_1") {
                continue;
            }
            let mid = budgets.len() / 2;
            if best.as_ref().is_none_or(|(_, bg)| g[mid] < bg[mid]) {
                best = Some((nm.clone(), g));
            }
        }
        if let Some((nm, g)) = &best {
            print!("{:>24}", format!("best: {nm}"));
            for v in g {
                print!(" {v:>9.4}");
            }
            println!();
        }

        // ---- gain by level: the consequence the tautologies were reaching for ----
        println!(
            "\n  gain(k) by level -- a parent's N x N grid and its children's are different\n\
             \x20 approximation families, so gain need not be positive:"
        );
        println!(
            "{:>7} {:>8} {:>9} {:>11} {:>11} {:>11} {:>11}",
            "level", "quads", "gain<0", "min", "median", "p90", "max"
        );
        for l in 0..levels {
            let mut g: Vec<f64> = {
                let w = 1u32 << l;
                (0..w).flat_map(|iy| (0..w).map(move |ix| (l, ix, iy))).map(|k| cache.gain(k)).collect()
            };
            let neg = g.iter().filter(|x| **x < 0.0).count();
            let lo = g.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = g.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{l:>7} {:>8} {neg:>9} {lo:>11.3e} {:>11.3e} {:>11.3e} {hi:>11.3e}",
                g.len(),
                q(&mut g, 0.5),
                q(&mut g, 0.9)
            );
        }

        // ---- err_sum by level: does it fall gradually, or collapse only at the bottom? ----
        //
        // Flat across 0..levels-1 and zero at `levels` means essentially ALL gain sits at the last
        // split, greedy's shallow gains are noise, and concentrating budget is the worst available
        // allocation. Normalised per pixel so levels are comparable.
        println!(
            "\n  err_sum by level, per pixel of the quad (so levels are comparable), \
             normalised by the root's:"
        );
        println!("{:>7} {:>12} {:>12} {:>12}", "level", "median", "p10", "p90");
        let root_pp = cache.get((0, 0, 0)).err_sum / (res * res) as f64;
        for l in 0..=levels {
            let w = 1u32 << l;
            let span = res >> l;
            let mut e: Vec<f64> = (0..w)
                .flat_map(|iy| (0..w).map(move |ix| (l, ix, iy)))
                .map(|k| cache.get(k).err_sum / (span * span) as f64 / root_pp)
                .collect();
            println!(
                "{l:>7} {:>12.5} {:>12.5} {:>12.5}",
                q(&mut e, 0.5),
                q(&mut e, 0.1),
                q(&mut e, 0.9)
            );
        }

        // ---- CRITERION VS UNIFORM, PER LEVEL ----
        //
        // "Beats random" says nothing where the optimum IS breadth-first, and on a smooth field
        // it is. The measurement that matters is criterion vs **uniform**, and it is read per
        // level rather than aggregated: if a criterion only beats breadth-first below some
        // depth, that is evidence for a depth-dependent strategy with the depth MEASURED. If the
        // gap is flat across levels, the rank formulation is already self-adapting and nothing
        // needs adding.
        //
        // Evaluated at the budgets where uniform **exactly completes** a level,
        // `B_d = 1 + 4*(4^d - 1)/3`. At any other budget uniform sits mid-level and the
        // comparison would be scoring where its partial row happened to stop.
        //
        // `captured` = `(uniform - row) / (uniform - dp)`: the share of the improvement available
        // over breadth-first that the criterion takes. 1.0 is the optimum, 0.0 is no better than
        // refining uniformly, negative is worse than doing nothing clever.
        {
            let uni = rows.iter().find(|r| r.0 == "uniform").expect("uniform row");
            println!(
                "\n  criterion vs UNIFORM, per level -- at the budgets where uniform completes a\n\
                 \x20 level exactly. captured = (uniform - row)/(uniform - dp); 1.0 = optimum,\n\
                 \x20 0.0 = no better than breadth-first, <0 = worse:"
            );
            println!(
                "{:>7} {:>9} {:>10} {:>10} {:>10} {:>26} {:>10}",
                "level", "B", "uniform", "dp", "uni-dp", "best criterion", "captured"
            );
            for d in 1..=levels {
                let s_d = ((1usize << (2 * d)) - 1) / 3;
                let b_d = 1 + 4 * s_d;
                if b_d > full {
                    break;
                }
                // The printed ladder need not contain `b_d`, so replay each row at it directly.
                let u = cache.error_of(&cache.leaves_at(Rank::Uniform, b_d));
                let dpv = dp.at_budget(b_d);
                let mut bestc = (f64::INFINITY, String::new());
                for (nm, _, _) in &rows {
                    if nm.starts_with("random")
                        || nm.starts_with("greedy_lookahead_1")
                        || nm == "uniform"
                    {
                        continue;
                    }
                    let r = *runs.iter().find(|r| &r.name() == nm).unwrap();
                    let e = cache.error_of(&cache.leaves_at(r, b_d));
                    if e < bestc.0 {
                        bestc = (e, nm.clone());
                    }
                }
                // **A denominator at the arithmetic floor is the finding, not a ratio.** Where
                // uniform IS the optimum -- which is what a smooth field looks like -- `u - dpv`
                // is a rounding epsilon and can be either sign, so the quotient reads as a large
                // positive or negative number describing nothing. Print the identity instead.
                let denom = u - dpv;
                let cap = if denom.abs() <= 1e-12 * u.abs().max(f64::MIN_POSITIVE) {
                    "uni==dp".to_string()
                } else {
                    format!("{:.4}", (u - bestc.0) / denom)
                };
                println!(
                    "{d:>7} {b_d:>9} {u:>10.5} {dpv:>10.5} {denom:>10.2e} {:>26} {cap:>10}",
                    format!("{} {:.5}", bestc.1, bestc.0)
                );
            }
            let _ = uni;
        }

        // ---- IS THE FIELD SMOOTH, OR IS IT AMPLIFIED NOISE? ----
        //
        // The lightness ramp is each region's own p1-p99, so a field with no dynamic range has
        // its NOISE stretched to full scale. `criterion_metric`'s guard did not fire on `far`,
        // and `far` being auto-ranged noise is a standing finding -- so it matters whether
        // `far`'s flat `err_sum` is physics being absent or noise failing to resolve.
        //
        // **Energy drift is the wrong floor to compare against**, and that is why the absolute
        // arm cleared: a tame region has a tiny drift AND a tiny spread, so the floor falls with
        // the field it is meant to bound and the arm is a ratio in disguise.
        //
        // The discriminator that is not: **spatial coherence**. Noise is incoherent between
        // neighbouring quads; a smooth field is coherent by definition. Lag-1 neighbour
        // correlation of the ramped scalar, per level.
        {
            let mut d: Vec<f64> = cache
                .quads
                .values()
                .map(|q| q.red.worst_energy_drift)
                .filter(|x| x.is_finite() && *x > 0.0)
                .collect();
            let drift = q(&mut d, 0.5);
            println!(
                "\n  ramp window (p1,p99) = ({:.3e}, {:.3e}), span x{:.3}; median energy drift \
                 {:.3e}, guard floor 100x = {:.3e}",
                cache.ramp.0,
                cache.ramp.1,
                cache.ramp.1 / cache.ramp.0.max(f64::MIN_POSITIVE),
                drift,
                100.0 * drift
            );
            println!(
                "  lag-1 neighbour correlation of the ramped scalar (spread_shape median), by \
                 level.\n\x20 ~0 is incoherent = noise; ~1 is a smooth field the ramp is \
                 magnifying:"
            );
            println!("{:>7} {:>8} {:>12} {:>12} {:>12}", "level", "quads", "rho(lag-1)", "p1", "p99");
            for l in 1..=levels {
                let w = 1u32 << l;
                let f = |ix: u32, iy: u32| -> f64 {
                    cache.get((l, ix, iy)).red.signal(Criterion::Within, Agg::Median)
                };
                let mut v: Vec<f64> = Vec::new();
                let (mut a, mut b): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
                for iy in 0..w {
                    for ix in 0..w {
                        let x = f(ix, iy);
                        if x.is_finite() {
                            v.push(x);
                        }
                        for (jx, jy) in [(ix + 1, iy), (ix, iy + 1)] {
                            if jx < w && jy < w {
                                let (y, z) = (x, f(jx, jy));
                                if y.is_finite() && z.is_finite() {
                                    a.push(y);
                                    b.push(z);
                                }
                            }
                        }
                    }
                }
                let rho = pearson(&a, &b);
                println!(
                    "{l:>7} {:>8} {rho:>12.4} {:>12.3e} {:>12.3e}",
                    v.len(),
                    q(&mut v.clone(), 0.01),
                    q(&mut v, 0.99)
                );
            }
        }

        // ---- who spends the budget where ----
        //
        // The direct answer to "why the corner". If the optimum is uniform-depth and greedy's is
        // concentrated, greedy's collapse is an allocation failure and nothing to do with the
        // metric.
        let probe = *budgets.iter().find(|&&b| b >= 1535).unwrap_or(budgets.last().unwrap());
        println!("\n  leaf levels at B = {probe} -- where each strategy spent the budget:");
        let s_probe = (probe - 1) / 4;
        println!("{:>24}  {}", "dp_optimal", hist(&dp.leaves(s_probe.min(dp.max_splits)), levels));
        for r in [
            Rank::GreedyLookahead1,
            Rank::Signal(Criterion::Within, Agg::Median),
            Rank::Signal(Criterion::FracHotBetween, Agg::Median),
            Rank::Random(1),
        ] {
            println!("{:>24}  {}", r.name(), hist(&cache.leaves_at(r, probe), levels));
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         `dp_optimal` is the exact minimum over ALL tree-shaped leaf sets at each budget. It is a\n\
         BOUND: no row may sit below it. A row that does is a harness bug, and every error(B)\n\
         number in the corpus is suspect until it is found.\n\n\
         `greedy_lookahead_1` is greedy on immediate delta-error. It is a reference and NOT a\n\
         bound -- it has been measured below the random band on `far`. Read the three roles from\n\
         the header: floor, reference, ceiling.\n\n\
         The `distinct` column decides what a flat row means. A flat curve on thousands of\n\
         distinct values is an ordering that is actively bad; a flat curve on one distinct value\n\
         is no ordering at all, and the row is the tie-break's scan order. They are different\n\
         faults with different fixes and error(B) alone cannot tell them apart.\n\n\
         On a SMOOTH field the two converge on the same answer: a quad's spread over a gradient\n\
         is proportional to its cell width, so argmax-on-spread IS breadth-first, and a constant\n\
         signal falls through to a level-first tie-break which is also breadth-first. That is why\n\
         `far`'s rows agree to five digits across every budget. A criterion can only differ from\n\
         breadth-first where the field is not smooth -- and smoothness is exactly where there is\n\
         nothing to find. `far` degenerating is correct behaviour, not a defect."
    );
}
