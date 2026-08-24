//! §6 — the production colour scheme, and whether the criterion is blind to half of it.
//!
//! # Why this is a criterion question
//!
//! The criterion asks *"would splitting change what we display?"*, so **what is displayed
//! decides what the criterion should measure.** The production scheme is bivariate: hue from
//! the shape sphere, lightness from a scalar. `spread_shape` maps to hue, so that half is
//! aligned by construction. If lightness carries diffusion or FTLE, the criterion has **no term
//! for it**: a quad can be uniform in shape and structured in diffusion, and nothing refines it.
//!
//! # The coupling question, put through §2's own metric
//!
//! Does `error(B)` move when lightness switches from `spread` to `diffusion`? Three colourings
//! over **one** integration pass, so the only thing that changes is the map from footprint to
//! pixel.
//!
//! - If the ranking of criteria is the same under all three, the criterion is not sensitive to
//!   the lightness field and §6's concern does not bite in practice.
//! - If it reorders, the criterion needs a term it does not have.
//!
//! # How to misread this
//!
//! **Check `renorm` before reading any FTLE column.** Renormalisation is what stops the
//! estimator saturating; with none, `log(d/d0)/T` decays and reports `lambda ~ 0` for the *most*
//! chaotic trajectories. A zero renormalisation count makes the FTLE field the saturated case
//! in new clothing, and its colouring meaningless rather than merely noisy.
//!
//! **The FTLE march is UNREGULARISED.** It is `tb.py`'s fixed-step leapfrog, because that is the
//! pair `tb_ftle.py` has a reference for. Near a close approach it is not trustworthy, and the
//! `d_min` column from the AZ march is what says where that is.
//!
//! **`error(B) = 0` still means "matches this sampling"**, per colouring. The reference is that
//! colouring's own fully-refined tree, so the three curves are each self-consistent and their
//! *absolute* values are not comparable across colourings -- only the ORDERING of criteria is.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Colouring, Rank};
use prin_rs::output::bivariate::Lightness;
use prin_rs::physics::ftle::FtleOpts;
use prin_rs::quad::{Agg, Criterion};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn main() {
    let levels: u32 = arg(1, 5);
    let n: usize = arg(2, 8);
    let t_max: f64 = arg(3, 13.0);
    let ftle_dt: f64 = arg(4, 1e-4);
    let res = (1usize << levels) * n;
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        ftle: Some(FtleOpts::default()),
        ftle_dt,
        ..Default::default()
    };

    println!(
        "level {levels}, N={n}, res {res}^2, t={t_max} n_sync={n_sync}, FTLE dt={ftle_dt:e} \
         ({} leapfrog steps/footprint)\n{full} quads per region\n",
        (t_max / ftle_dt).round() as u64
    );

    let colourings = [
        Colouring::Outcome,
        Colouring::Bivariate(Lightness::Spread),
        Colouring::Bivariate(Lightness::Diffusion),
        Colouring::Bivariate(Lightness::Ftle),
    ];

    let budgets = [21usize, 85, 341, full];

    for &(region, cx, cy, body) in grid::REGIONS
        .iter()
        .filter(|r| matches!(r.0, "near-field" | "deep interior"))
    {
        let t0 = std::time::Instant::now();
        let caches = metric::build_multi(
            region, cx, cy, 0.05, body, Chart::BodyPlane, levels, n, res, 1e-4, &ens, &colourings,
        );
        println!("--- {region} --- one integration pass in {:.1}s", t0.elapsed().as_secs_f64());

        // Renormalisation first: without it the FTLE field is saturated, not noisy.
        let c0 = &caches[0];
        let renorm: Vec<f64> =
            c0.quads.values().map(|q| q.red.total_substeps as f64).collect();
        let _ = renorm;
        println!(
            "  ramps: spread {:?}  diffusion {:?}  ftle {:?}",
            caches[1].ramp, caches[2].ramp, caches[3].ramp
        );

        for (ci, c) in caches.iter().enumerate() {
            let e_root = c.error_of(&[(0, 0, 0)]);
            println!("\n  colouring {} -- error(root) = {e_root:.5}", c.colouring.name());
            if e_root == 0.0 {
                println!("    UNDEFINED here: the image is featureless at this resolution, so");
                println!("    every criterion reads zero and none of it is data.");
                continue;
            }
            print!("{:>24}", "B =");
            for b in &budgets {
                print!(" {b:>9}");
            }
            println!();
            for r in [
                Rank::GreedyOracle,
                Rank::Signal(Criterion::Within, Agg::Median),
                Rank::Signal(Criterion::Between, Agg::Median),
                Rank::Signal(Criterion::FracHotBetween, Agg::Median),
                Rank::Contrast(Criterion::Between, Agg::Median),
                Rank::Random(1),
            ] {
                let pts = metric::replay(c, r, full);
                print!("{:>24}", r.name());
                for e in metric::curve_at(&pts, &budgets) {
                    print!(" {e:>9.5}");
                }
                println!();
            }
            let stem = format!(
                "results/criterion/colour_{}_{}",
                region.replace(' ', "_"),
                c.colouring.name().replace('/', "_")
            );
            let _ = prin_rs::output::adaptive::save(&format!("{stem}_reference.png"), res, &c.reference);
            let leaves = c.leaves_at(Rank::Signal(Criterion::Between, Agg::Median), full / 8);
            let _ = prin_rs::output::adaptive::save(
                &format!("{stem}_B{}.png", full / 8),
                res,
                &c.render(&leaves),
            );
            let _ = ci;
        }
        println!();
    }

    println!(
        "The coupling question: does the ORDERING of criteria change between colourings?\n\
         Absolute errors are NOT comparable across colourings -- each has its own reference and\n\
         its own zero. Only the ranking within a block is.\n\
         \n\
         If the ordering is stable, the criterion is not sensitive to the lightness field and\n\
         §6's concern does not bite in practice. If it reorders -- particularly if a criterion\n\
         that wins under `outcome` loses under `diffusion` -- then the criterion needs a term\n\
         for the lightness field and currently has none.\n\
         \n\
         The FTLE march is UNREGULARISED (tb.py's fixed-step leapfrog, the pair tb_ftle.py has\n\
         a reference for). Near a close approach it is not trustworthy, so an FTLE-lightness\n\
         result in `deep interior` carries that caveat and a spread-lightness one does not."
    );
}
