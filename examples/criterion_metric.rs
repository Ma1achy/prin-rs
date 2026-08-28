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
//! - **greedy_lookahead_1** — greedy on immediate `Δerror`. **A strong reference, not a ceiling.**
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

/// The colouring the criterion is scored under.
///
/// **This is a criterion parameter, not a presentation one.** `error(B)` measures image change,
/// so what is displayed decides which quads matter -- the standing rule is *choose a criterion
/// under the colouring that will ship*. PR #13 scored every criterion under `Outcome`, which at
/// `t = 13` is a saturated categorical label: near-field's reference image is two flat colours.
const SHIPPING: metric::Colouring =
    metric::Colouring::Bivariate(prin_rs::output::colour::Scalar::ShapeSpread);

fn main() {
    let levels: u32 = arg(1, 6);
    let n: usize = arg(2, 8);
    let tau: f64 = arg(3, 1e-4);
    let t_max: f64 = arg(4, 13.0);
    // **Where the artefacts land, and it is an argument because it has to be.** This example
    // writes images, caches and an APNG under `<root>/criterion` and `<root>/animated`, and the
    // committed ones there are 512^2. A validation pass at reduced `levels` -- the whole point of
    // which is to fire the `dp_optimal` bound assertion cheaply -- would overwrite them with a
    // small raster, and a small raster reads as a rendering fault rather than a stale file. That
    // has cost this project two round trips. Point it at a scratch directory instead.
    let root: String = std::env::args().nth(5).unwrap_or_else(|| "results".into());
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
        // The numpy reference's ungated escape test, with escape terminal: every result in
        // this diagnostic predates both the distance gate and the closure criterion, and is
        // quoted against that form.
        escape_rule: prin_rs::outcome::EscapeRule::Reference,
        closure_k: 1,
        stop_on_escape: true,
        keep_boundary_shapes: true,
        keep_drift_hist: false,
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
        // The production colouring, and the outcome control beside it, from ONE integration
        // pass. `outcome` is categorical and saturated at t = 13 -- near-field's reference image
        // is two flat colours -- so it is kept as a control rather than as the target.
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
            &[SHIPPING, metric::Colouring::Outcome],
        );
        let control = caches.pop().unwrap();
        let cache = caches.pop().unwrap();

        // The footprints, so no future colouring change costs another integration. PRQC stores
        // a BAKED err_sum and cannot be replayed under a new colouring; this can.
        let stem0 = region.replace(' ', "_");
        if let Ok(f) = std::fs::File::create(format!("{root}/criterion/{stem0}_t{t_max}.fcache")) {
            let mut w = std::io::BufWriter::new(f);
            let fp = cache.footprints_from(&px_of, t_max);
            let _ = prin_rs::output::fcache::write(&mut w, &fp);
        }
        // A cheap assertion that the replay path is live on real data, not only in the unit
        // test: recolouring to the control must reproduce the control bitwise.
        {
            let fp = cache.footprints_from(&px_of, t_max);
            match cache.recolour(&fp, metric::Colouring::Outcome) {
                Ok(r) => println!(
                    "  replay check: recolour to `outcome` reproduces the control reference {}",
                    if r.reference == control.reference { "BITWISE" } else { "**DIFFERENTLY**" }
                ),
                Err(e) => println!("  replay check FAILED: {e}"),
            }
        }
        drop(px_of);
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
             error(root)={e_root:.5} error(full)={e_full:.1e}\n\
             {:>14} ramp [{:.4e}, {:.4e}]  span x{:.3}{}",
            cache.trajectories,
            "",
            cache.ramp.0,
            cache.ramp.1,
            cache.ramp.1 / cache.ramp.0.max(f64::MIN_POSITIVE),
            // Two arms, because the ratio alone is not enough. `far` reads a span of x8 -- above
            // any sensible ratio threshold -- over a window of (1.3e-9, 1.1e-8). `spread_shape`
            // is a mean chord distance on the unit sphere and is dimensionless, so a p99 of 1e-8
            // means the copies agree to eight digits: the field is at the level of the
            // integrator's own arithmetic, not of the physics.
            //
            // The absolute arm compares against a MEASURED floor rather than a chosen constant:
            // the region's own median energy drift. A field whose whole range sits within two
            // orders of that is not distinguishable from integration noise.
            {
                let mut d: Vec<f64> = cache
                    .quads
                    .values()
                    .map(|q| q.red.worst_energy_drift)
                    .filter(|x| x.is_finite() && *x > 0.0)
                    .collect();
                let floor = if d.is_empty() {
                    0.0
                } else {
                    100.0 * prin_rs::stats::quantile(&mut d, 0.5)
                };
                // **THE THIRD ARM, AND THE ONE THAT DISCRIMINATES.** The two above read
                // AMPLITUDE, and amplitude cannot separate a small real signal from noise --
                // which is exactly why neither fired on `far`. Worse, the absolute arm's floor
                // is the region's own median energy drift, so in a tame region the floor falls
                // with the field it is meant to bound: a ratio in disguise. Measured on `far`,
                // `ramp.1 = 1.064e-8` against a floor of `4.478e-9` -- clears by 2.4x, and the
                // field really is at the eighth digit.
                //
                // Noise is spatially INCOHERENT between neighbouring quads; a smooth field is
                // coherent by definition. Lag-1 neighbour correlation of the ramped scalar over
                // the level-3 grid separates them at any amplitude. `far` reads **0.9984** there
                // and its p1/p99 halve exactly per level, which is `spread ~ g*w` measured: a
                // real gradient of tiny magnitude, not an amplified noise floor.
                let rho = {
                    let l = 3u32.min(levels);
                    let w = 1u32 << l;
                    let f = |ix: u32, iy: u32| {
                        cache.get((l, ix, iy)).red.signal(Criterion::Within, Agg::Median)
                    };
                    let (mut a, mut b): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
                    for iy in 0..w {
                        for ix in 0..w {
                            for (jx, jy) in [(ix + 1, iy), (ix, iy + 1)] {
                                if jx < w && jy < w {
                                    let (y, z) = (f(ix, iy), f(jx, jy));
                                    if y.is_finite() && z.is_finite() {
                                        a.push(y);
                                        b.push(z);
                                    }
                                }
                            }
                        }
                    }
                    let n = a.len() as f64;
                    if a.len() < 2 {
                        f64::NAN
                    } else {
                        let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
                        let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
                        for (x, y) in a.iter().zip(&b) {
                            sxy += (x - ma) * (y - mb);
                            sxx += (x - ma) * (x - ma);
                            syy += (y - mb) * (y - mb);
                        }
                        if sxx <= 0.0 || syy <= 0.0 { f64::NAN } else { sxy / (sxx * syy).sqrt() }
                    }
                };
                println!(
                    "  ramp (p1,p99) = ({:.3e}, {:.3e}) span x{:.3}; noise floor (100x median                      drift) {:.3e}; lag-1 coherence at level 3 rho={:.4}",
                    cache.ramp.0,
                    cache.ramp.1,
                    cache.ramp.1 / cache.ramp.0.max(f64::MIN_POSITIVE),
                    floor,
                    rho
                );
                if rho.is_finite() && rho < 0.5 {
                    "  <-- AUTO-RANGED OVER NOISE: the ramped scalar is spatially INCOHERENT between neighbouring quads, so the p1-p99 window is stretched over something with no structure in it. This arm reads coherence rather than amplitude, which is what the two below could not do."
                } else if cache.ramp.1 / cache.ramp.0.max(f64::MIN_POSITIVE) < 2.0
                    || cache.ramp.1 < floor
                {
                    "  <-- AUTO-RANGED OVER NOISE: the ramp is normalised to this region's own p1-p99, so a field with no dynamic range -- or one whose whole range sits at the integrator's own arithmetic floor -- is stretched to full scale and error(B) becomes nonzero for a region with nothing in it. Read this before the curve."
                } else {
                    ""
                }
            }
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
            // **The baseline, and it was missing from every table in the corpus.** Random is a
            // floor no strategy should sit below; breadth-first is the bar a criterion has to
            // clear, and measured against it `frac_hot_between/median` never wins in
            // `near-field` at any budget while sitting well clear of the random band.
            Rank::Uniform,
            Rank::GreedyLookahead1,
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
            Rank::GreedyLookahead1PerCost,
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

        // ---- the ceiling ----
        //
        // The exact minimum over ALL tree-shaped leaf sets at each budget. `greedy_lookahead_1`
        // was read as a bound for two PRs and is not one -- it sits BELOW the random band on
        // `far`. Printed first so every row beneath is read against it, and asserted rather than
        // trusted: a row below the ceiling is a harness bug.
        let dp = cache.dp_optimal((full - 1) / 4);
        println!(
            "  dp_optimal: {} splits in {:.2}s; prefix-min binds at {} split counts \
             (where it binds, a split made the image WORSE)",
            dp.max_splits,
            dp.elapsed_s,
            dp.prefix_min_binds.len()
        );

        print!("{:>22}", "B =");
        for b in &budgets {
            print!(" {b:>9}");
        }
        println!();
        {
            let mut worst = f64::INFINITY;
            for (_, curve) in &rows {
                for (&b, &e) in budgets.iter().zip(curve) {
                    worst = worst.min(e - dp.at_budget(b));
                }
            }
            assert!(
                worst >= -1e-9,
                "a ranking beat the exact tree optimum by {worst:e} -- the harness is wrong and \
                 every error(B) number in this table is suspect"
            );
            print!("{:>22}", "dp_optimal");
            for b in &budgets {
                print!(" {:>9.5}", dp.at_budget(*b));
            }
            println!("   <- CEILING (worst margin {worst:+.1e})");
        }
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
        // ---- the raw dump ----
        //
        // Without this the complete tree lives only in RAM for one process, and reproducing any
        // table above means paying the 2.8M-trajectory integration again. Every criterion's
        // scalar is dumped whatever this run ranked on, which is what makes offline comparison
        // real rather than aspirational.
        if let Ok(f) = std::fs::File::create(format!("{root}/criterion/{stem0}_t{t_max}.qcache")) {
            let mut w = std::io::BufWriter::new(f);
            let _ = prin_rs::output::qcache::write(&mut w, &cache, &ens, tau);
        }

        // ---- the error(B) curves as a figure ----
        //
        // Log y, because the curves span four decades and the interesting part is the bottom
        // one. Exact zeros are NOT snapped to the floor: they get their own band below the log
        // panel, one row per series, because reaching zero IS the result on several of these
        // curves and the previous figure drew 15 of 17 `far` series on top of each other.
        {
            use prin_rs::output::plot::{palette, Figure, Series};
            let live: Vec<&(String, Vec<f64>)> = rows
                .iter()
                .filter(|(name, _)| !name.starts_with("random") || name == "random[1]")
                .collect();
            let pal = palette(live.len());
            let mut ser: Vec<Series> = Vec::new();
            for (i, (name, curve)) in live.iter().enumerate() {
                let control = name.starts_with("random") || name.starts_with("greedy");
                let rgb = if name.starts_with("greedy_lookahead_1") && !name.contains("cost") {
                    (255, 255, 255)
                } else if name.starts_with("random") {
                    (150, 150, 160)
                } else {
                    pal[i]
                };
                let pts = budgets.iter().zip(curve.iter()).map(|(&b, &e)| (b as f64, e)).collect();
                let sr = Series::new(name.clone(), pts, rgb);
                ser.push(if control { sr.dashed() } else { sr });
            }
            let fig = Figure {
                title: format!("{region} — error(B) against the fully-refined tree"),
                x_label: "budget B (quads computed)".into(),
                y_label: "mean per-pixel OKLab distance".into(),
                series: ser,
                y_lo: 1e-6,
                y_hi: 0.2,
                notes: vec![
                    format!(
                        "t_max = {t_max}, n_sync = {n_sync}, N = {n}, E+1 = {}, eta = {}, \
                         levels = {levels}, {res}x{res}, f64, criterion = display colouring",
                        ens.n_extra + 1,
                        ens.eta
                    ),
                    "error = 0 means MATCHES THIS SAMPLING, not correct: the reference is the \
                     fully-refined tree at one sample per pixel, and at the screen floor \
                     sub-pixel structure is sampled arbitrarily."
                        .into(),
                    "THREE ROLES: floor = random band (several seeds); reference = \
                     greedy_lookahead_1, greedy on immediate delta-error, neither optimal nor a \
                     bound; ceiling = dp_optimal, the exact minimum over all tree-shaped leaf \
                     sets. Greedy has been measured BELOW the random band on far."
                        .into(),
                    "A label reading (k/n) means n-k points were non-finite and were DROPPED. \
                     A high NaN fraction is a property to read, not a defect: term_grad is NaN \
                     on 97.1% of near-field and still reaches zero by B = 383."
                        .into(),
                ],
            };
            if let Err(e) = fig.save(&format!("{root}/criterion/curve_{stem0}_t{t_max}")) {
                eprintln!("figure failed: {e}");
            }
        }

        // ---- images: the reference, and the best and worst criteria at one budget ----
        //
        // Drawn at TRUE per-quad texel sizes, so a coarse leaf is visibly coarse. A uniform
        // render of the same tree would show where boundaries fell rather than what the system
        // displays, which is the instrument this project has already been caught using once.
        let stem = region.replace(' ', "_");
        let _ = prin_rs::output::adaptive::save(
            &format!("{root}/criterion/{stem}_reference.png"),
            res,
            &cache.reference,
        );
        // The reference's own tree is the complete one at the screen floor; its wire is the
        // uniform mesh every other tree is being compared against.
        let deepest: Vec<metric::Key> = {
            let w = 1u32 << levels;
            (0..w).flat_map(|iy| (0..w).map(move |ix| (levels, ix, iy))).collect()
        };
        let _ = prin_rs::output::adaptive::save(
            &format!("{root}/criterion/{stem}_reference_wire.png"),
            res,
            &cache.render_wire(&deepest),
        );
        let mid = full / 8;
        for r in [
            Rank::Uniform,
            Rank::GreedyLookahead1,
            Rank::Signal(Criterion::Within, Agg::Median),
            Rank::Signal(Criterion::Between, Agg::Median),
            Rank::Signal(Criterion::FracHotBetween, Agg::Median),
            Rank::Random(1),
        ] {
            let leaves = cache.leaves_at(r, mid);
            let tag = r.name().replace('/', "_").replace(['[', ']'], "");
            // Both, always. The texel size says what is displayed; the wire says where the tree
            // cut. Neither substitutes for the other, and PR #11's failure was reading one as
            // though it were the other.
            let _ = prin_rs::output::adaptive::save(
                &format!("{root}/criterion/{stem}_B{mid}_{tag}.png"),
                res,
                &cache.render(&leaves),
            );
            let _ = prin_rs::output::adaptive::save(
                &format!("{root}/criterion/{stem}_B{mid}_{tag}_wire.png"),
                res,
                &cache.render_wire(&leaves),
            );
        }
        // ---- the animation: the budget being spent, oracle against the shipped default ----
        //
        // Side by side in one frame on purpose. Two separate animations would make a reader
        // hold one in memory while watching the other, which is exactly the comparison the
        // picture exists to remove. Drawn at TRUE per-quad texel sizes, so a coarse leaf is
        // visibly coarse and the shipped default's misspent budget is legible rather than
        // inferred from a table.
        {
            let mut frames: Vec<Vec<u8>> = Vec::new();
            let mut wire_frames: Vec<Vec<u8>> = Vec::new();
            let mut b = 5usize;
            while b <= full {
                let la = cache.leaves_at(Rank::GreedyLookahead1, b);
                let lc = cache.leaves_at(Rank::Signal(Criterion::Within, Agg::Median), b);
                let mut fr =
                    prin_rs::output::apng::side_by_side(&cache.render(&la), &cache.render(&lc), res, res);
                prin_rs::output::apng::divide(&mut fr, res, res, [230, 230, 240]);
                frames.push(fr);
                let mut fw = prin_rs::output::apng::side_by_side(
                    &cache.render_wire(&la),
                    &cache.render_wire(&lc),
                    res,
                    res,
                );
                prin_rs::output::apng::divide(&mut fw, res, res, [230, 230, 240]);
                wire_frames.push(fw);
                b = (b * 2 + 1).min(full).max(b + 1);
                if b == full {
                    let la = cache.leaves_at(Rank::GreedyLookahead1, full);
                    let lc = cache.leaves_at(Rank::Signal(Criterion::Within, Agg::Median), full);
                    let mut fr = prin_rs::output::apng::side_by_side(
                        &cache.render(&la), &cache.render(&lc), res, res);
                    prin_rs::output::apng::divide(&mut fr, res, res, [230, 230, 240]);
                    frames.push(fr);
                    let mut fw = prin_rs::output::apng::side_by_side(
                        &cache.render_wire(&la), &cache.render_wire(&lc), res, res);
                    prin_rs::output::apng::divide(&mut fw, res, res, [230, 230, 240]);
                    wire_frames.push(fw);
                    break;
                }
            }
            let _ = prin_rs::output::apng::write(
                &format!("{root}/animated/budget_{stem0}_t{t_max}_animated.png"),
                res * 2,
                res,
                &frames,
                1,
                2,
            );
            let _ = prin_rs::output::apng::write(
                &format!("{root}/animated/budget_{stem0}_t{t_max}_wire_animated.png"),
                res * 2,
                res,
                &wire_frames,
                1,
                2,
            );
            // Three representative frames as ordinary PNGs -- first, middle and last -- so
            // nothing here depends on APNG support. Not all of them: at 1024^2 the side-by-side
            // is 2048x1024 and thirteen of them per region per horizon came to 154 MB of
            // duplicated content, since the animation already carries every frame.
            let picks = [0usize, frames.len() / 2, frames.len().saturating_sub(1)];
            for &i in picks.iter() {
                if let Some(fr) = frames.get(i) {
                    let _ = prin_rs::output::adaptive::save_rect(
                        &format!("{root}/criterion/budget_{stem0}_t{t_max}_{i:02}.png"),
                        res * 2,
                        res,
                        fr,
                    );
                }
                if let Some(fr) = wire_frames.get(i) {
                    let _ = prin_rs::output::adaptive::save_rect(
                        &format!("{root}/criterion/budget_{stem0}_t{t_max}_wire_{i:02}.png"),
                        res * 2,
                        res,
                        fr,
                    );
                }
            }
            println!(
                "  {} animation frames (oracle | within/median) at true texel sizes",
                frames.len()
            );
        }

        println!("  images at B={mid} written to {root}/criterion/{stem}_*.png\n");
    }

    println!(
        "Read the oracle-to-random gap FIRST, per region. If greedy_lookahead_1 and the random band\n\
         overlap, the metric does not discriminate there and nothing below it is readable.\n\
         \n\
         THE TABLE HAS THREE ROLES AND THEY ARE NOT INTERCHANGEABLE.\n\
           floor      random lo/hi          several seeds, read as a band, never one trace\n\
           reference  greedy_lookahead_1    greedy on immediate delta-error -- NEITHER OPTIMAL\n\
                                            NOR A BOUND. Measured BELOW the random band on\n\
                                            `far`: 0.54760 against 0.48550-0.52047 at B = 1535.\n\
           BASELINE   uniform               breadth-first -- the bar a criterion must clear.\n\
           ceiling    dp_optimal            the exact minimum over ALL tree-shaped leaf sets.\n\
                                            No row may sit below it; one that does is a harness\n\
                                            bug, and it is asserted, not trusted.\n\
         \n\
         A criterion beating `greedy_lookahead_1` indicates LOOKAHEAD VALUE and is not a bug:\n\
         greedy declines a low-gain split that unlocks large gains two levels down. There is no\n\
         assertion anywhere that it dominates, and there must not be.\n\
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
