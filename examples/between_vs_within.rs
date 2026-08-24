//! §1 — does the criterion read the wrong disagreement?
//!
//! # The claim under test
//!
//! `ensemble_spread` is a statistic over the `E+1` copies of one footprint, and `reduce`
//! aggregates those `N^2` numbers. The brief reads that as a category error: refinement buys
//! more footprints, so only *between*-footprint variation is reducible, and a
//! *within*-footprint statistic is measuring something no resolution changes.
//!
//! # The prediction on record, made before this run
//!
//! **The conclusion is right and the stated mechanism is wrong.** The premise — "the ICs there
//! are identical up to perturbation" — does not describe this implementation. `jitter_frac` is
//! 0.5 and `halton_offset` returns `[-1, 1)^2` scaled by cell width, so the copies span the
//! **whole cell, edge to edge**: a quasi-random sample of exactly the area the footprint stands
//! for. The corroboration was already on record before this example existed — the Halton
//! control's true `alpha` is exactly **1.0**, and an irreducible within-point statistic would
//! have `alpha == 0` by construction, because splitting would not shrink it.
//!
//! So both arms are between-point statistics. What differs is **scale** (cell against quad),
//! **sample count** (`E+1` against `N^2`), and the **aggregation**, which is the one that
//! bites: a median over footprints discards where the hot ones are, and a clean boundary
//! touches few of them.
//!
//! If that is right, `rho` is **high overall and lower on boundary-containing quads**, and the
//! scale ratio is near 1 while the aggregation flips decisions. If `rho` is low everywhere, the
//! prediction is wrong and the brief is right.
//!
//! # How to misread this table
//!
//! **Do not read the `all` column.** A region-wide `rho` is dominated by tame quads, where both
//! arms read near zero and agree trivially; it would read high whatever happens at the
//! boundaries, which is the entire population a scheduler spends its budget on. The `bnd`
//! column — quads whose hot set looks like a boundary — is the measurement.
//!
//! **Do not read `rho` without `d90`.** Two orderings can correlate at +0.99 while a tail
//! crosses the whole list. `d90`/`dmax` are the per-quad rank displacements behind the scalar.
//!
//! **Check `dead` and `coll` first.** A difference can be small because both sides are right or
//! because both are dead. `coll` counts quads whose decode collapsed — fewer distinct ICs than
//! footprints, where every spread is a spread over repeats and reads as perfect resolution.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::quad::{Agg, Criterion};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::stats;
use prin_rs::{grid, quad};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn med(v: &[f64]) -> f64 {
    let mut s: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if s.is_empty() {
        return f64::NAN;
    }
    quad::quantile(&mut s, 0.5)
}

fn main() {
    let budget: usize = arg(1, 2000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let viewport: usize = arg(4, 512);

    println!(
        "budget {budget} quads, tau={tau:e}, alpha_hi={alpha_hi}, N=8, E+1=8, t=13, f64, \
         viewport {viewport}^2 (veto ON)\n"
    );

    // The copies are needed for `within_pooled` — the within arm at the between arm's sample
    // count, which is what separates a scale effect from small-sample bias.
    let ens = EnsembleCfg { keep_copy_shapes: true, refine_flagged: false, ..Default::default() };

    println!(
        "{:>14} {:>6} {:>5} {:>5} {:>9} {:>9} {:>9} {:>7} {:>7} {:>8} {:>6} {:>6}",
        "region", "quads", "coll", "dead", "rho all", "rho mix", "rho bnd",
        "d90", "dmax", "med hot", "n mix", "n bnd"
    );

    let mut carry: Vec<(String, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> = Vec::new();

    for &(region, _, _, _) in grid::REGIONS.iter().filter(|r| {
        matches!(r.0, "far" | "near-field" | "deep interior")
    }) {
        let s = grid::region(region, 2, 2, 0.05).unwrap();
        let cfg = SchedCfg {
            keep_pixels: true,
            budget,
            tau_display: tau,
            alpha_hi,
            alpha_lo: alpha_hi,
            agg: Agg::Median,
            criterion: Criterion::Within,
            camera: Some(Camera::framing(s.cx, s.cy, 0.05, viewport)),
            ..Default::default()
        };
        let (t, st) = scheduler::descend(s.cx, s.cy, 0.05, s.body, &cfg, &ens, Precision::F64);

        // The adaptive render of the tree this criterion built, and the same tree's spread
        // overlay -- the diagnostic that exposed `deep interior`'s failure in the first place.
        {
            let cam = Camera::framing(s.cx, s.cy, 0.05, 512);
            let (img, _t) = prin_rs::output::adaptive::render(
                &t, &st.pixels, &cam, 512,
                prin_rs::output::adaptive::TexelMode::Adaptive,
                prin_rs::output::png::outcome_rgb,
            );
            let _ = prin_rs::output::adaptive::save(
                &format!("results/criterion/tree_{}.png", region.replace(' ', "_")),
                512,
                &img,
            );
            let mut wimg = img.clone();
            let boxes = prin_rs::output::wire::boxes_from_tree(&t, &cam, 512);
            let deepest = boxes.iter().map(|b| b.level).max().unwrap_or(1);
            prin_rs::output::wire::draw(&mut wimg, 512, 512, &boxes, deepest.max(1));
            let _ = prin_rs::output::adaptive::save(
                &format!("results/criterion/tree_{}_wire.png", region.replace(' ', "_")),
                512,
                &wimg,
            );
        }

        // The full v2 quad dump, so every column in this table can be recomputed offline and
        // every column NOT in it is still available.
        let stem = format!("results/criterion/between_{}.prnq", region.replace(' ', "_"));
        if let Ok(f) = std::fs::File::create(&stem) {
            let mut w = std::io::BufWriter::new(f);
            let _ = prin_rs::output::tree::write(&mut w, &t, &cfg, &ens, &st, region, "f64");
        }

        // Every computed quad, leaf or not. The root has no reduction worth reading until it
        // is computed, and `descend` computes every node it creates.
        let all: Vec<usize> = (0..t.nodes.len()).filter(|&i| t.nodes[i].red.n_footprints > 0).collect();

        let coll = all.iter().filter(|&&i| t.nodes[i].red.between_collapsed()).count();
        // Alive on both sides: not collapsed, and neither arm identically zero. A quad where
        // one arm is dead contributes an arbitrary rank and would be read as agreement.
        let live: Vec<usize> = all
            .iter()
            .cloned()
            .filter(|&i| {
                let r = &t.nodes[i].red;
                !r.between_collapsed()
                    && r.spread_median.is_finite()
                    && r.between_spread.is_finite()
                    && (r.spread_median > 0.0 || r.between_spread > 0.0)
            })
            .collect();
        let dead = all.len() - live.len() - coll;

        let w: Vec<f64> = live.iter().map(|&i| t.nodes[i].red.spread_median).collect();
        let b: Vec<f64> = live.iter().map(|&i| t.nodes[i].red.between_spread).collect();
        let rho_all = stats::spearman(&w, &b);

        // Boundary-containing, by the hot-set layout rather than by a spread threshold: a
        // connected thin structure is the case the two arms are predicted to disagree on.
        let bnd: Vec<usize> = live
            .iter()
            .cloned()
            .filter(|&i| t.nodes[i].red.layout_within.looks_like_boundary(t.n, 1.5))
            .collect();
        let wb: Vec<f64> = bnd.iter().map(|&i| t.nodes[i].red.spread_median).collect();
        let bb: Vec<f64> = bnd.iter().map(|&i| t.nodes[i].red.between_spread).collect();
        let rho_bnd = stats::spearman(&wb, &bb);

        // **Why `n bnd` is small has to be visible.** A quad that is *uniformly* hot has no
        // internal hot/cold edge, so `perimeter_ratio` is 0 and it is correctly not a
        // boundary — it is saturated. In a chaotic region that is most quads, which leaves the
        // boundary stratum with almost no population and a rho that means nothing.
        //
        // `mixed` is the weaker, better-populated stratum: any quad whose hot set is a proper
        // subset, i.e. one containing a transition of some kind, filamentary or not.
        let frac_hot: Vec<f64> = live
            .iter()
            .map(|&i| t.nodes[i].red.frac_above_tau_within)
            .collect();
        let mix: Vec<usize> = live
            .iter()
            .cloned()
            .filter(|&i| {
                let f = t.nodes[i].red.frac_above_tau_within;
                f > 0.02 && f < 0.98
            })
            .collect();
        let wm: Vec<f64> = mix.iter().map(|&i| t.nodes[i].red.spread_median).collect();
        let bm: Vec<f64> = mix.iter().map(|&i| t.nodes[i].red.between_spread).collect();
        let rho_mix = stats::spearman(&wm, &bm);

        let disp = stats::rank_displacement(&w, &b);
        let (_, _, d90, _) = stats::interdecile(&disp);
        let dmax = disp.iter().cloned().fold(0.0f64, f64::max);

        println!(
            "{:>14} {:>6} {:>5} {:>5} {:>9.4} {:>9.4} {:>9.4} {:>7.4} {:>7.4} {:>8.3} {:>6} {:>6}",
            region, all.len(), coll, dead, rho_all, rho_mix, rho_bnd,
            d90, dmax, med(&frac_hot), mix.len(), bnd.len()
        );

        let scale: Vec<f64> = live
            .iter()
            .map(|&i| {
                let r = &t.nodes[i].red;
                r.between_matched / r.spread_median
            })
            .collect();
        let count: Vec<f64> = live
            .iter()
            .map(|&i| {
                let r = &t.nodes[i].red;
                r.within_pooled / r.between_shape
            })
            .collect();
        carry.push((region.to_string(), w, b, scale, count));
    }

    // ---------------------------------------------------------------------------------
    println!(
        "\n{:>14} {:>12} {:>12} {:>12} {:>12}",
        "region", "med within", "med between", "scale", "count"
    );
    println!(
        "{:>14} {:>12} {:>12} {:>12} {:>12}",
        "", "(cell,E+1)", "(quad,N^2)", "matched/w", "pooled/b"
    );
    for (region, w, b, scale, count) in &carry {
        println!(
            "{:>14} {:>12.4e} {:>12.4e} {:>12.4} {:>12.4}",
            region,
            med(w),
            med(b),
            med(scale),
            med(count)
        );
    }

    println!(
        "\nRead `coll` and `dead` before anything else: a rho over collapsed or one-sided quads\n\
         is a correlation between two things that are not being measured.\n\
         \n\
         Then read `rho bnd`, never `rho all`. The all-quads column is dominated by tame quads\n\
         where both arms read near zero and agree trivially; it is high by construction and\n\
         says nothing about the population the budget is spent on.\n\
         \n\
         `scale` is between_matched/within at EQUAL sample count, so it isolates cell-extent\n\
         against quad-extent. `count` is within_pooled/between at EQUAL extent, so it isolates\n\
         the small-sample bias. Neither is meaningful without the other: comparing the raw\n\
         within and between numbers moves both at once, and a spread estimator's expectation\n\
         depends on its sample count (E+1=2 reports 0.539 of E+1=32's value in near-field).\n\
         \n\
         A rho over fewer than ~20 quads is a draw, not a measurement. Read `n mix`/`n bnd`\n\
         beside every stratified rho and discard the ones with no population -- `med hot` says\n\
         why: where it sits at 1.000 the quad is SATURATED, every footprint hot, no internal\n\
         edge and correctly not a boundary. That is a real reading, not a failure of the mask.\n\
         \n\
         If `rho mix` is high and `scale` is near 1, the two arms are the same statistic at two\n\
         extents and §1's mechanism is wrong -- the fault is the aggregation, and the fix is\n\
         §3.1/§3.2 rather than a new arm. If `rho mix` is low, §1 is right and this note is\n\
         wrong."
    );
}
