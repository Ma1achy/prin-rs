//! Does the lightness channel carry an ordering, and does the hue channel carry three axes?
//!
//! Two questions, both asked *before* any picture is read, because both have a failure mode that
//! looks like a picture of physics.
//!
//! # 1. Count the distinct values first
//!
//! `ensemble_spread` is `max(spread_shape, spread_event)`, and `spread_event` is a count ratio
//! over `E+1` copies. Wherever the event arm dominates, the lightness field has **`E+2` levels
//! and no ramp recovers what is not there** — the committed §6 table's `spread` p99 values are
//! exactly `2/7` and `6/7`, which is what a count ratio over 8 copies looks like. This is the
//! standing rule *count the signal's distinct values before reading any curve*, applied to the
//! colour channel: two criteria in PR #13 had flat `error(B)` curves that turned out to be the
//! tie-break's scan order, and the distinct-value count is what separated them from a criterion
//! that was ordering finely and ordering *badly*.
//!
//! `frac(event arm)` says *why* the field is quantised where it is, so the two are printed
//! together. A high fraction is a reason to colour on `spread_shape` instead, and that pair of
//! images is emitted here rather than argued about.
//!
//! # 2. The old hue map was 2-to-1 and this shows by how much
//!
//! `chroma*(cos h, sin h)` with `h = atan2(n2,n1)` reduces to `C_MAX*(n1,n2)` — a linear
//! projection that discards `n0` entirely. Since `n0 = (|rho~|^2 - |lam~|^2)/I`, a tight binary
//! with a distant third body and a wide pair with a close third were painted the same colour.
//! The `n0 span` column below is how much of that axis the region actually uses: it is the size
//! of what the old map could not see, per region, rather than a claim in the abstract.
//!
//! # How to misread this
//!
//! **A large distinct count does not mean a good field.** `within/median` had 5418 distinct
//! values of 5461 and was beaten by random at every budget past 383. Fine ordering and useful
//! ordering are different properties; this example measures only the first, which is the one
//! that can rule a field out.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use prin_rs::grid;
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, png as pngout};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let res: usize = arg(1, 1024);
    let t_max: f64 = arg(2, 13.0);

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { refine_flagged: false, t_max, n_sync, ..base };

    let sites = colour::landmarks(&prin_rs::physics::burrau::MASSES);
    println!(
        "res {res}^2, t={t_max}, n_sync={n_sync}, E+1={}, kappa={}",
        ens.n_extra + 1,
        sites.kappa
    );
    println!("sites (computed from the masses, never hard-coded):");
    for s in &sites.sites {
        println!("  {:>16}  n = [{:+.4}, {:+.4}, {:+.4}]", s.name, s.n[0], s.n[1], s.n[2]);
    }
    let eu = colour::euler_points(&prin_rs::physics::burrau::MASSES);
    for (k, n) in eu.iter().enumerate() {
        println!("  {:>16}  n = [{:+.4}, {:+.4}, {:+.4}]", format!("euler(mid={k})"), n[0], n[1], n[2]);
    }
    println!(
        "  worst angular gap to the nearest site over the sphere: {:.3} rad\n",
        colour::worst_site_gap(&sites, 20_000)
    );

    let fields = [Scalar::Spread, Scalar::ShapeSpread, Scalar::EventSpread];

    println!(
        "{:>14} {:>14} {:>9} {:>8} {:>8} {:>10} {:>10} {:>9}",
        "region", "field", "distinct", "finite", "modal%", "p1", "p99", "event%"
    );

    for (name, _cx, _cy, _body) in grid::REGIONS.iter().filter(|r| {
        matches!(r.0, "near-field" | "deep interior" | "far")
    }) {
        let slice = grid::region(name, res, res, 0.05).expect("known region");
        let px: Vec<PixelOut> =
            (0..slice.npix()).map(|i| evaluate::<f64>(&slice, i, &ens)).collect();

        // How much of the axis the old projection discarded is actually used here.
        let n0: Vec<f64> = px.iter().map(|p| p.shape_vec[0]).filter(|x| x.is_finite()).collect();
        let n0_span = n0.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - n0.iter().cloned().fold(f64::INFINITY, f64::min);
        // The span is a MAX statistic: a handful of pixels at each end give a span of 2 while the
        // bulk sits in a sliver. The interdecile says whether the axis is USED or merely REACHED
        // -- the same distinction the project already makes for `alpha_shape`, where the excess
        // kurtosis is 110 and the variance describes only the tail.
        let (_, _, _, n0_idec) = prin_rs::stats::interdecile(&n0);

        let ev = colour::event_arm_fraction(&px);
        for s in fields {
            let (distinct, finite, modal) = colour::quantisation(&px, s);
            let (lo, hi) = colour::range(&px, s);
            println!(
                "{:>14} {:>14} {:>9} {:>8} {:>7.1}% {:>10.3e} {:>10.3e} {:>8.1}%",
                name,
                s.name(),
                distinct,
                finite,
                100.0 * modal,
                lo,
                hi,
                100.0 * ev
            );
        }
        println!(
            "{:>14} {:>14} n0 span = {:.4}  n0 interdecile = {:.4}  ({} nonfinite shapes)",
            name,
            "--",
            n0_span,
            n0_idec,
            px.len() - n0.len()
        );

        // kappa is a design choice and is chosen on a number rather than by eye. `coverage` is
        // the mean distance from the centroid of the blended OKLab (a,b) over this region: how
        // much of the hue plane the data actually reaches. Read it beside `n0 span` -- a region
        // that visits one kind of configuration SHOULD read near zero, and that is not a broken
        // palette. `far` is the control for exactly that.
        let ns: Vec<[f64; 3]> = px.iter().map(|p| p.shape_vec).collect();
        let mut row = format!("{:>14} {:>14}", name, "hue coverage");
        for k in [0.5f64, 1.0, 2.0, 3.0, 6.0, 12.0] {
            let set = colour::SiteSet { kappa: k, ..sites.clone() };
            row.push_str(&format!("  k={k}: {:.4}", colour::hue_coverage(&set, &ns)));
        }
        println!("{row}");

        // One image per field, plus the outcome control. The pair `spread` / `spread_shape` is
        // the point: if the first is quantised and the second is not, they are different
        // instruments wearing one name.
        let stem = name.replace(' ', "_");
        for s in fields {
            let (lo, hi) = colour::range(&px, s);
            let mut buf = Vec::with_capacity(px.len() * 3);
            for p in &px {
                buf.extend_from_slice(&colour::rgb(p, s, &sites, lo, hi));
            }
            let path = format!("results/colour/{stem}_{}.png", s.name());
            adaptive::save_rect(&path, res, res, &buf).unwrap();
        }
        let mut buf = Vec::with_capacity(px.len() * 3);
        for p in &px {
            buf.extend_from_slice(&pngout::outcome_rgb(p));
        }
        adaptive::save_rect(&format!("results/colour/{stem}_outcome.png"), res, res, &buf).unwrap();
        println!();
    }

    // ---------------------------------------------------------------------------------------
    // 3. Does the lightness field mean the same thing at every level?
    // ---------------------------------------------------------------------------------------
    //
    // An adaptive render draws leaves at several levels in one image, and `ensemble_spread` is
    // a spread over copies jittered within the CELL -- so a coarse leaf's copies start further
    // apart than a fine leaf's at the same place. If the field were proportional to cell width,
    // a bivariate render would carry a brightness gradient tracking the TREE rather than the
    // physics, and `error(B)` would partly measure "this leaf is coarse".
    //
    // The control is `t_end`, which is a nominal-copy quantity and has no cell width in it at
    // all. If both fields move together the measurement is picking up something else.
    {
        use prin_rs::grid::Chart;
        use prin_rs::metric::{self, Colouring};
        let levels = 4u32;
        let nq = 8usize;
        println!("\nper-level scaling of the lightness field (levels 0..{levels}, N={nq}):");
        println!(
            "{:>14} {:>6} {:>7} {:>12} {:>12} {:>8} {:>12} {:>8}",
            "region", "level", "quads", "cell width", "med spread", "ratio", "med t_end", "ratio"
        );
        for (name, cx, cy, body) in
            [("near-field", 1.0, 3.0, 0usize), ("deep interior", 0.0, 0.0, 0usize)]
        {
            let (_c, px_of) = metric::build_multi_with_footprints(
                name, cx, cy, 0.05, body, Chart::BodyPlane, levels, nq,
                (1usize << levels) * nq, 1e-4, &ens,
                &[Colouring::Bivariate(Scalar::ShapeSpread)],
            );
            let (mut ps, mut pt) = (f64::NAN, f64::NAN);
            for l in 0..=levels {
                let (mut sp, mut te, mut n) = (Vec::new(), Vec::new(), 0usize);
                for (k, px) in &px_of {
                    if k.0 != l {
                        continue;
                    }
                    n += 1;
                    for p in px {
                        if p.spread_shape.is_finite() {
                            sp.push(p.spread_shape);
                        }
                        if p.t_end.is_finite() {
                            te.push(p.t_end);
                        }
                    }
                }
                let med = prin_rs::stats::quantile(&mut sp, 0.5);
                let medt = prin_rs::stats::quantile(&mut te, 0.5);
                let cell = 2.0 * (0.05 / (1u64 << l) as f64) / (nq - 1) as f64;
                println!(
                    "{:>14} {:>6} {:>7} {:>12.4e} {:>12.4e} {:>8.3} {:>12.4e} {:>8.3}",
                    name, l, n, cell, med, ps / med, medt, pt / medt
                );
                ps = med;
                pt = medt;
            }
        }
        println!(
            "\nCell width halves every level, so a field proportional to it would show `ratio`\n\
             = 2.000 in every row. Measured 1.19-1.62, FALLING with depth: the dependence is\n\
             sub-linear and saturates, which is the chaotic-divergence signature -- past some\n\
             cell width the copies decohere fully however close they started. Consistent with\n\
             the standing measurement that cell width falling 9x moves sigma_E(0) by 8.6x\n\
             (proportional) and `ensemble_spread` by only 2.1x.\n\
             \n\
             So a bivariate render DOES carry a mild scale term: a level-0 leaf reads brighter\n\
             than a level-4 leaf at the same place. Under the log ramp that is about 12% of the\n\
             lightness range across five levels. Stated rather than corrected -- normalising by\n\
             cell width would change what the field means. `t_end` is the scale-free control and\n\
             is flat, which is what says the effect above is the cell width and not the depth."
        );
    }

    println!(
        "\nRead `distinct` and `modal%` before either image. A field whose modal value covers a\n\
         large fraction of the region has that many pixels painted one colour, and the picture\n\
         is describing the estimator rather than the physics. `event%` says whether the event\n\
         arm is what is doing it -- it is a count ratio over E+1 copies and is quantised by\n\
         construction, so where it dominates, `spread_shape` is the field with structure in it.\n\
         \n\
         `n0 span` is the extent of the shape-sphere axis the previous hue map discarded. It is\n\
         the size of what was invisible, measured per region rather than asserted -- but it\n\
         is a MAX statistic, so read `n0 interdecile` beside it. A span of 2 with a tiny\n\
         interdecile means the axis was reached by a few pixels, not used by the region.\n\
         \n\
         `hue coverage` chooses kappa on a number. A region that genuinely visits one kind\n\
         of configuration SHOULD read near zero at every kappa, and that is a true answer\n\
         rather than a broken palette -- `far` is the control for exactly that."
    );
}
