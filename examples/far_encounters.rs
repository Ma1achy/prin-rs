//! **Does `far` contain close approaches?** The standing text says it does not. It does, and they
//! are the deepest of any Burrau region measured here.
//!
//! The claim on record, from the first clean AZ-vs-Heggie gallery:
//!
//! > `far` is the only AZ win... That **fits** the mechanism rather than contradicting it —
//! > `far` is smooth and wide, **no pair ever approaches**, so there is no unregularised third
//! > side to remove and the longer reconstruction chain costs a constant with nothing bought
//! > back.
//!
//! The mechanism half of that was already withdrawn and replaced with a scale guess, labelled a
//! guess. The premise underneath it was never checked, and it is false.
//!
//! # Why it survived
//!
//! "far" names **where the sampled body sits**, not how the system behaves. `grid::region("far")`
//! puts body 0 at `(1.0, 13.0)` and varies it over `half = 0.05`; bodies 1 and 2 keep their
//! Burrau positions **identically across the whole slice**. So the two of them fall together and
//! pass very close on every pixel, thirteen units away from anything the slice is varying. A
//! region named for a large separation is not a region without an encounter.
//!
//! # The discriminator, and it refuted what this comment first predicted
//!
//! A uniform `d_min` is consistent with "the encounter is between the two bodies the slice does
//! not move", and consistency is not a measurement. The second arm moves body 0 further out.
//! **The prediction written here first was "`d_min` must barely move".** It moves enormously:
//!
//! ```text
//!   body 0 at        d_min p50      ratio    spread p99/p01
//!   (1.0, 13.0)      7.1324e-7     1.0000        1.18
//!   (2.0, 26.0)      1.0874e-8     0.0152        1.10
//!   (4.0, 52.0)      1.6335e-10    0.0002        1.05
//!   (10.0, 130.0)    6.4587e-13    0.0000        1.02
//! ```
//!
//! The framing was binary — *does the encounter involve body 0?* — and the answer is not binary.
//! Body 0 does not **participate** in the encounter; it **controls its depth**. It is the
//! perturber that keeps bodies 1 and 2 from falling straight into each other, and taking it away
//! makes the collision more nearly exact. Both readings agree the encounter is the `(1,2)` pair;
//! the depth dependence is a further fact that the wrong prediction is what exposed.
//!
//! And it is quantitative: the ratios per doubling are **65.6, 66.5** and the last rung is
//! `252.9` against `(10/4)^6 = 244`, so
//!
//! ```text
//!   d_min ~ r0^-6      over three decades of perturber distance
//! ```
//!
//! which is the tidal signature: the induced impact parameter goes as the tidal acceleration
//! `~1/r0^3` integrated over the infall, and the closest approach goes as its square. The spread
//! tightening monotonically to 1.02 says the same thing from the other side — the further the
//! perturber, the more identical every pixel's encounter becomes.
//!
//! # And the committed gallery already said so, in the column beside the one being read
//!
//! `results/output/integrator_gallery.txt` carries an `escape`/`coll` pair per row. At 1024^2:
//!
//! ```text
//!   far             escape 0        coll 1048576   <- 100% of the frame collides
//!   deep_interior   escape 0        coll 1033184      98.5%
//!   preset_shape    escape 188      coll  850590      81.1%
//!   near-field      escape 0        coll   24649       2.4%
//! ```
//!
//! **`far` collides on every pixel and `near-field` on one in forty.** So "no pair ever
//! approaches" was contradicted by the project's own committed output the whole time, one column
//! to the right of `escape`, which reads `0` for `far` and is the number that was being read.
//! Nothing new had to be measured to falsify it — only the next column looked at. The `d_min`
//! table above is what says *how* deep and the discriminator is what says *why*, but the
//! refutation was free.
//!
//! Args: `res`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::energy;

fn q(v: &[f64], p: f64) -> f64 {
    let mut w = v.to_vec();
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    w[(((w.len() - 1) as f64) * p).round() as usize]
}

/// The DIAGNOSTIC configuration: no collision radius and no termination, so `d_min` is the closest
/// approach the trajectory actually reaches rather than the radius it was stopped at.
fn cfg() -> EnsembleCfg {
    EnsembleCfg::production().with_overrides(&[
        // **Pinned, and it must be.** This harness measures `d_min` in `far` to explain why AZ
        // wins there, so the trajectories have to be AZ's. It took the integrator from the
        // default, which was `Az` when this was written and is `Heggie` from 2026-09-02 --
        // a silent switch that would have measured the wrong integrator's close approaches
        // while the header still said `far` is the only AZ win.
        Override::Integrator(prin_rs::integrate::Integrator::Az),
        Override::TMax(13.0),
        Override::NSync(32),
        Override::RCollFrac(0.0),
        Override::StopOnEvent(false),
        Override::EscapeRule(EscapeRule::Closure(CLOSURE_TAU)),
        Override::RefineFlagged(false),
        Override::MaxSteps(400_000),
    ])
}

fn run(sl: &grid::Slice) -> Vec<f64> {
    let px: Vec<PixelOut> =
        (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(sl, k, &cfg())).collect();
    px.iter().map(|p| p.d_min_true).filter(|x| x.is_finite()).collect()
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(64);
    println!("{res}^2, t_max = 13, r_coll = 0 (no termination), refine_flagged = false\n");
    println!(
        "  {:>15} {:>10} {:>10} {:>10} {:>10} {:>9} {:>11} {:>11}",
        "region", "d_min p01", "p10", "p50", "p99", "R", "d_min/R p50", "spread p99/p01"
    );

    for case in ["far", "mid-field", "near-field", "body2 core", "deep interior"] {
        let Some(r) = grid::region(case, 4, 4, 0.05) else { continue };
        let sl = grid::Slice::body_plane(res, res, r.cx, r.cy, r.half, r.body).with_chart(r.chart);
        let d = run(&sl);
        let ic = grid::decode_state(&sl.chart, sl.body, sl.cx, sl.cy);
        let big_r = energy::hyperradius(&ic.s.r, &ic.m);
        println!(
            "  {case:>15} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>9.4} {:>11.3e} {:>11.2}",
            q(&d, 0.01), q(&d, 0.10), q(&d, 0.50), q(&d, 0.99), big_r,
            q(&d, 0.50) / big_r, q(&d, 0.99) / q(&d, 0.01)
        );
    }

    // ---- the discriminator ----
    println!(
        "\nDISCRIMINATOR: move the sampled body further out.\n\n\
         **The prediction written here first was `d_min must barely move`, and it is wrong.**\n\
         Body 0 does not participate in the encounter -- it CONTROLS ITS DEPTH, as the perturber\n\
         that keeps bodies 1 and 2 from falling straight into each other. `exponent` is the\n\
         local slope of `d_min` against perturber distance; the tidal prediction is -6, the\n\
         induced impact parameter going as the tidal acceleration `~1/r0^3` integrated over the\n\
         infall and the closest approach as its square.\n"
    );
    let r = grid::region("far", 4, 4, 0.05).unwrap();
    println!(
        "  {:>22} {:>12} {:>12} {:>10} {:>12}",
        "body 0 at", "d_min p50", "ratio", "exponent", "spread p99/p01"
    );
    let mut base: Option<f64> = None;
    let mut prev: Option<(f64, f64)> = None;
    for scale in [1.0f64, 2.0, 4.0, 10.0] {
        let sl = grid::Slice::body_plane(res, res, r.cx * scale, r.cy * scale, r.half, r.body)
            .with_chart(r.chart);
        let d = run(&sl);
        let p50 = q(&d, 0.50);
        let b = *base.get_or_insert(p50);
        let ex = prev.map(|(ps, pd): (f64, f64)| (p50 / pd).ln() / (scale / ps).ln());
        println!(
            "  {:>22} {p50:>12.4e} {:>12.4} {:>10} {:>12.2}",
            format!("({:.1}, {:.1})", r.cx * scale, r.cy * scale),
            p50 / b,
            ex.map(|x| format!("{x:.2}")).unwrap_or_else(|| "-".into()),
            q(&d, 0.99) / q(&d, 0.01)
        );
        prev = Some((scale, p50));
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **A region named for a large separation is not a region without an encounter.** `far`\n\
         names where the SAMPLED body sits. Bodies 1 and 2 hold their Burrau positions on every\n\
         pixel of the slice, so their mutual encounter is the same one everywhere -- which is why\n\
         `d_min` is both the deepest here and the most uniform.\n\n\
         This matters for every comparison on `far`: an integrator is being graded on a region\n\
         whose whole frame contains a near-collision, not on a smooth field. The unregularised\n\
         control failing on all of it is then the expected result rather than a surprise, and\n\
         logH doing worst there is the collision finding at field scale -- it slows the clock but\n\
         does not remove the singularity, so the encounter still has to be RESOLVED."
    );
}
