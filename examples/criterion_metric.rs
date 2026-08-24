//! §2 — image change per quad spent, with both controls.
//!
//! # What this measures
//!
//! Render the region at the screen floor -> REFERENCE. Run a criterion under budget `B`, render
//! adaptively -> `IMAGE(B)`. `error(B)` is the mean per-pixel OKLab distance between them. A
//! criterion is better if `error(B)` falls faster. **The whole curve is the result**, not one
//! number: a criterion can win at small budget and lose at large.
//!
//! # The two controls, without which the metric says nothing
//!
//! - **random** — several seeds, reported as a band. A single random trace is a draw.
//! - **greedy_oracle** — greedy on immediate `Δerror`. **A strong reference, not a ceiling.**
//!   Greedy is optimal only when gains are independent and immediately available; here a quad
//!   whose own split gains little may unlock children with large gains two levels down, and
//!   greedy declines it. **A criterion beating it indicates lookahead value, not an error.**
//!
//! # How to misread this table
//!
//! **Check the oracle-to-random separation first.** If they are close, the metric is not
//! discriminating in that region and no criterion result read from it means anything there.
//!
//! **`error = 0` means "matches this sampling", not "correct".** The reference is the
//! fully-refined tree at one sample per pixel — a specific finite sampling. At the screen floor
//! sub-pixel structure is sampled arbitrarily: which side of a filament a pixel lands on is an
//! accident of where its sample fell. It is the right common target for comparing criteria and
//! it is not a statement about image quality.
//!
//! **Criteria enter as ORDERINGS, never against `tau`.** That is §2's reframe, and it also
//! disposes of a confound: the between arm runs 1.17x the within arm in near-field and 9.56x in
//! `far`, so a threshold comparison would have scored that rescaling instead of the signal.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Rank};
use prin_rs::quad::{Agg, Criterion};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let levels: u32 = arg(1, 6);
    let n: usize = arg(2, 8);
    let tau: f64 = arg(3, 1e-4);
    let t_max: f64 = arg(4, 13.0);
    let res = (1usize << levels) * n;

    // **`n_sync` scales with `t_max`.** `dtau = eta*dt_left/(A0*B0)`, so holding `n_sync` fixed
    // while `t_max` moves changes the step size, and the two runs are different discretisations
    // rather than one trajectory at two playheads.
    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    // The temporal accumulators need each copy's per-boundary shape vector. Enabled here
    // because §5's question -- do they add anything at t = 13? -- is a §2 question, and the
    // only honest way to answer it is to put them through the same curve as everything else.
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        keep_boundary_shapes: true,
        ..Default::default()
    };
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;

    println!(
        "complete tree to level {levels}, N={n}, E+1={}, res {res}^2, tau={tau:e}, \
         t={t_max} n_sync={n_sync}, f64\n\
         {full} quads per region, {} trajectories\n",
        ens.n_extra + 1,
        full * n * n * (ens.n_extra + 1)
    );

    let budgets: Vec<usize> = {
        let mut b = vec![5usize];
        while *b.last().unwrap() * 2 < full {
            b.push(b.last().unwrap() * 2 + 1);
        }
        b.push(full);
        b
    };

    for &(region, cx, cy, body) in grid::REGIONS
        .iter()
        .filter(|r| matches!(r.0, "far" | "near-field" | "deep interior"))
    {
        let t0 = std::time::Instant::now();
        let cache = metric::build(
            region, cx, cy, 0.05, body, Chart::BodyPlane, levels, n, res, tau, &ens, metric::Colouring::Outcome,
        );
        let build_s = t0.elapsed().as_secs_f64();

        // Reference sanity, before any curve is read: the deepest-level tree must give exactly
        // zero, and the root alone must give something clearly non-zero. If the root error is
        // near zero the region has no structure at this scale and no ranking can be judged on
        // it.
        let deepest: Vec<metric::Key> = {
            let w = 1u32 << levels;
            (0..w).flat_map(|iy| (0..w).map(move |ix| (levels, ix, iy))).collect()
        };
        let e_full = cache.error_of(&deepest);
        let e_root = cache.error_of(&[(0, 0, 0)]);

        println!(
            "--- {region} --- built in {build_s:.1}s, {} trajectories; \
             error(root)={e_root:.5} error(full)={e_full:.1e}",
            cache.trajectories
        );

        // ---- degeneracy, read BEFORE the curves ----
        //
        // A ranking is only a ranking if the signal takes different values. In a saturated
        // region a criterion can be constant across every quad, at which point the "ranking"
        // is whatever the tie-break happens to be -- a fixed scan order wearing a criterion's
        // name. That is a different failure from ranking badly and the curves cannot tell them
        // apart, so it is measured directly.
        // Termination and escape fractions, printed before the term_grad row is readable at
        // all. **These are different quantities and the difference matters**: `t_end` is set by
        // whichever terminating event came first, so in `deep interior` `term` reads ~0.99
        // while `escape` is ~0 -- those are collisions. Quoting the first as an escape fraction
        // would contradict the standing result that zero of 1024 near-field pixels escape at
        // t = 13, while appearing to agree with it.
        let tm: Vec<f64> = cache.quads.values().map(|q| q.red.terminated_fraction).collect();
        let ec: Vec<f64> = cache.quads.values().map(|q| q.red.escape_fraction).collect();
        println!(
            "  terminated: mean {:.4}, {} of {} quads with any  |  escaped: mean {:.4}, {} quads",
            tm.iter().sum::<f64>() / tm.len() as f64,
            tm.iter().filter(|&&x| x > 0.0).count(),
            tm.len(),
            ec.iter().sum::<f64>() / ec.len() as f64,
            ec.iter().filter(|&&x| x > 0.0).count(),
        );

        println!(
            "{:>22} {:>8} {:>8} {:>7} {:>9}",
            "signal", "distinct", "modal%", "nan%", "spread"
        );
        for (c, a) in [
            (Criterion::Within, Agg::Median),
            (Criterion::Within, Agg::Mean),
            (Criterion::Within, Agg::P90),
            (Criterion::Between, Agg::Median),
            (Criterion::MaxOfBoth, Agg::Median),
            (Criterion::FracHotWithin, Agg::Median),
            (Criterion::FracHotBetween, Agg::Median),
            (Criterion::Layout, Agg::Median),
            (Criterion::TerminationGradient, Agg::Median),
            (Criterion::RunningMax, Agg::Median),
            (Criterion::FirstDivergence, Agg::Median),
        ] {
            let vals: Vec<f64> = cache.quads.values().map(|q| q.red.signal(c, a)).collect();
            let nanf = vals.iter().filter(|x| !x.is_finite()).count() as f64 / vals.len() as f64;
            let mut bits: Vec<u64> = vals.iter().map(|v| v.to_bits()).collect();
            bits.sort_unstable();
            let distinct = { let mut b = bits.clone(); b.dedup(); b.len() };
            let modal = {
                let (mut best, mut run, mut prev) = (0usize, 0usize, u64::MAX);
                for b in &bits {
                    if *b == prev { run += 1 } else { run = 1; prev = *b }
                    best = best.max(run);
                }
                best as f64 / bits.len() as f64
            };
            let lo = vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{:>22} {distinct:>8} {:>7.1}% {:>6.1}% {:>9.3e}",
                format!("{}/{}", c.name(), a.name()),
                100.0 * modal,
                100.0 * nanf,
                hi - lo
            );
        }
        println!();

        let mut rows: Vec<(String, Vec<f64>)> = Vec::new();
        let runs: Vec<Rank> = vec![
            Rank::GreedyOracle,
            Rank::Signal(Criterion::Within, Agg::Median),
            Rank::Signal(Criterion::Within, Agg::Mean),
            Rank::Signal(Criterion::Within, Agg::P90),
            Rank::Signal(Criterion::Between, Agg::Median),
            Rank::Signal(Criterion::MaxOfBoth, Agg::Median),
            Rank::Signal(Criterion::FracHotWithin, Agg::Median),
            Rank::Signal(Criterion::FracHotBetween, Agg::Median),
            Rank::Signal(Criterion::Layout, Agg::Median),
            Rank::Signal(Criterion::TerminationGradient, Agg::Median),
            Rank::Signal(Criterion::RunningMax, Agg::Median),
            Rank::Signal(Criterion::FirstDivergence, Agg::Median),
            Rank::Contrast(Criterion::Within, Agg::Median),
            Rank::Contrast(Criterion::Between, Agg::Median),
            Rank::GreedyOraclePerCost,
            Rank::Random(1),
            Rank::Random(2),
            Rank::Random(3),
            Rank::Random(4),
            Rank::Random(5),
        ];
        for r in runs {
            let pts = metric::replay(&cache, r, full);
            rows.push((r.name(), metric::curve_at(&pts, &budgets)));
        }

        print!("{:>22}", "B =");
        for b in &budgets {
            print!(" {b:>9}");
        }
        println!();
        for (name, curve) in &rows {
            if name.starts_with("random") {
                continue;
            }
            print!("{name:>22}");
            for e in curve {
                print!(" {e:>9.5}");
            }
            println!();
        }
        // The random band: min and max across seeds at each budget, never one trace.
        let rnd: Vec<&Vec<f64>> = rows
            .iter()
            .filter(|(n, _)| n.starts_with("random"))
            .map(|(_, c)| c)
            .collect();
        for (label, pick) in [("random lo", true), ("random hi", false)] {
            print!("{label:>22}");
            for (j, _) in budgets.iter().enumerate() {
                let vals: Vec<f64> = rnd.iter().map(|c| c[j]).collect();
                let v = if pick {
                    vals.iter().cloned().fold(f64::INFINITY, f64::min)
                } else {
                    vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                };
                print!(" {v:>9.5}");
            }
            println!();
        }
        // ---- images: the reference, and the best and worst criteria at one budget ----
        //
        // Drawn at TRUE per-quad texel sizes, so a coarse leaf is visibly coarse. A uniform
        // render of the same tree would show where boundaries fell rather than what the system
        // displays, which is the instrument this project has already been caught using once.
        let stem = region.replace(' ', "_");
        let _ = prin_rs::output::adaptive::save(
            &format!("results/criterion/{stem}_reference.png"),
            res,
            &cache.reference,
        );
        let mid = full / 8;
        for r in [
            Rank::GreedyOracle,
            Rank::Signal(Criterion::Within, Agg::Median),
            Rank::Signal(Criterion::Between, Agg::Median),
            Rank::Signal(Criterion::FracHotBetween, Agg::Median),
            Rank::Random(1),
        ] {
            let leaves = cache.leaves_at(r, mid);
            let img = cache.render(&leaves);
            let tag = r.name().replace('/', "_").replace(['[', ']'], "");
            let _ = prin_rs::output::adaptive::save(
                &format!("results/criterion/{stem}_B{mid}_{tag}.png"),
                res,
                &img,
            );
        }
        println!("  images at B={mid} written to results/criterion/{stem}_*.png\n");
    }

    println!(
        "Read the oracle-to-random gap FIRST, per region. If greedy_oracle and the random band\n\
         overlap, the metric does not discriminate there and nothing below it is readable.\n\
         \n\
         `greedy_oracle` is a reference, not a bound: greedy declines a low-gain split that\n\
         unlocks large gains two levels down, so a criterion beating it at some budget indicates\n\
         LOOKAHEAD VALUE and is not a bug. There is no assertion anywhere that it dominates.\n\
         \n\
         `nan%` is how much of the region a signal declines to score. A signal that is NaN\n\
         nearly everywhere is not ranking; NaN never wins a comparison, so its curve is the\n\
         tie-break's scan order. Read `escaped:` above for why escape_grad is that at t=13.\n\
         \n\
         Read `distinct`/`modal%` before the curves. A criterion whose signal takes ONE value\n\
         across the region is not ranking at all -- its curve is the tie-break's scan order, and\n\
         a flat curve there means `no ordering`, not `a bad ordering`. Those are different\n\
         faults with different fixes and the error curve alone cannot separate them.\n\
         \n\
         error(full) is exactly 0 by construction -- the reference IS the fully-refined tree.\n\
         That makes the zero locatable, not meaningful as image quality: at the screen floor\n\
         sub-pixel structure is sampled arbitrarily, so `error = 0` reads as `matches this\n\
         sampling`."
    );
}
