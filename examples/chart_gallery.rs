//! Every chart family, rendered under the shipping colouring, with its wireframe and its dumps.
//!
//! # What each artefact answers
//!
//! The plain render says **what is displayed** — texels at true per-quad sizes, so a coarse leaf
//! is visibly coarse. The `_wire` twin says **where the tree cut**, brightness graded by level.
//! They answer different questions and neither substitutes for the other: PR #11 drew boundaries
//! over a *uniform* base, which conflated them, and `deep interior`'s bad tree survived a whole
//! build unnoticed. The level ladder (`_animated`) says **how the tree got there** — the same
//! descent truncated at each depth, so it is one tree seen at several playheads rather than
//! several unrelated trees.
//!
//! # How to misread this table
//!
//! **Leaf counts are slice-conditional to 4.3x. Compare within a chart, never across.** The
//! `alpha` distribution is the cross-chart quantity, and the interdecile rather than the
//! variance: excess kurtosis on `alpha_shape` is 110, so the variance is a statement about the
//! tail and the interdecile describes the bulk.
//!
//! **A chart that produces a prettier picture is not a better chart.** The measurement is
//! whether the criterion behaves consistently across charts, not which one looks best — and that
//! temptation is stronger now that the renders have structure in them. The `error(B)` curve is
//! the result; this gallery is a diagnostic.
//!
//! **Read the control line first.** `plane_00deg` must be bitwise `body_plane` — it is the same
//! chart written a second way. If it is not, the bases are wrong and every other row is
//! comparing different physics rather than different slices. It is an assertion here, not a
//! printed remark: PR #13's version only printed it.
//!
//! **The hue sites are computed from each chart's own nominal masses.** On the mass simplex the
//! landmarks move across the slice, so a single site set is a choice: the palette describes the
//! centre configuration and is held fixed across the image, because a per-pixel palette is not a
//! picture of anything. Every dump's `chart_params` records the chart it was built from.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Domain};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, png, wire};
use prin_rs::physics::decoder::Latent;
use prin_rs::quad::{Agg, Decision, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::{decode, stats};
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The leaf set of the same tree **truncated at `cap`**: nodes at or above `cap` that are either
/// leaves already or sit exactly at the cap.
///
/// One descent, six pictures. A budget ladder would need a fresh descent per frame, which is a
/// different tree each time and would make the animation a sequence of unrelated runs rather
/// than one refinement seen at several depths.
fn leaves_capped(t: &QuadTree, cap: u32) -> Vec<usize> {
    (0..t.nodes.len())
        .filter(|&i| {
            let q = &t.nodes[i];
            q.level <= cap && (q.children.is_none() || q.level == cap)
        })
        .collect()
}

fn render_leaves(
    t: &QuadTree,
    pixels: &[Vec<PixelOut>],
    cam: &Camera,
    res: usize,
    leaves: &[usize],
    rgb: &dyn Fn(&PixelOut) -> [u8; 3],
) -> Vec<u8> {
    // `adaptive::render` walks the tree's own leaves. To draw a truncated tree, build a shadow
    // tree whose leaf set is the truncated one -- cheaper and less error-prone than duplicating
    // the rasteriser, and it keeps the one endpoint-inclusive overhang rule in one place.
    let mut shadow = t.clone();
    let keep: std::collections::HashSet<usize> = leaves.iter().cloned().collect();
    for i in 0..shadow.nodes.len() {
        if keep.contains(&i) {
            shadow.nodes[i].children = None;
        }
    }
    adaptive::render(&shadow, pixels, cam, res, adaptive::TexelMode::Adaptive, |p| rgb(p)).0
}

fn main() {
    // Budget high enough that the descent stops on the CRITERION rather than on the budget.
    // At 4000 every chart but `body_plane` and `shape_sphere` hit the cap, which made their
    // leaf counts a fact about the budget and left large areas as coarse leaves -- and a coarse
    // leaf is drawn as one flat tile, because the render never interpolates. That reads as blur
    // and is not: it is an honest picture of an unrefined tree.
    let budget: usize = arg(1, 40000);
    let tau: f64 = arg(2, 1e-4);
    let alpha_hi: f64 = arg(3, 0.2);
    let res: usize = arg(4, 1024);

    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    // A base latent point. Deliberately not the origin: at z = 0 every sigmoid sits at 0.5 and
    // several coordinates would be at a symmetry point, which is exactly where a sign error
    // hides.
    let z0 = Latent {
        z_alpha: 0.35,
        z_beta: -0.45,
        z_q: [0.25, -0.15, 0.40, 0.05],
        z_mu: [0.20, -0.30],
    };

    // Twelve instances across all five families, **each carrying its own centre**.
    //
    // The centre is a property of the chart INSTANCE, not of its variant, and inferring it from
    // the variant is how this example first failed: `Chart::plane_for_body` zeroes the varying
    // body in `origin` and carries the whole position in `(u,v)`, so it is centred like
    // `BodyPlane` at (1,3) -- while `slice_gallery`'s oblique planes carry the configuration in
    // `origin` and are centred at zero. Centring the first at (0,0) sampled a box three units
    // away and gave 29 quads against 549. The control caught it; the fix is to stop inferring.
    let cases: Vec<(&str, Chart, f64, f64, f64)> = vec![
        ("body_plane", Chart::BodyPlane, 1.0, 3.0, 0.05),
        ("plane_00deg", Chart::plane_for_body(0), 1.0, 3.0, 0.05),
        ("shape_sphere", Chart::shape_at_burrau(0.0), 0.0, 0.0, 0.05),
        // The five named axis-aligned latent planes from the reference's table.
        ("latent_shape", Chart::latent_axes(z0, 0, 1), 0.0, 0.0, 1.5),
        ("latent_inner_p", Chart::latent_axes(z0, 2, 3), 0.0, 0.0, 1.5),
        ("latent_outer_p", Chart::latent_axes(z0, 4, 5), 0.0, 0.0, 1.5),
        ("latent_mass", Chart::latent_axes(z0, 6, 7), 0.0, 0.0, 1.5),
        ("latent_mixed", Chart::latent_axes(z0, 0, 4), 0.0, 0.0, 1.5),
        // Two oblique planes. The bases are recorded by `Chart::params` in every dump.
        (
            "latent_oblique_a",
            Chart::latent_oblique(
                z0,
                [0.3, -1.2, 0.5, 0.9, -0.4, 0.1, 0.7, -0.6],
                [1.1, 0.2, -0.8, 0.3, 0.6, -0.9, 0.15, 0.4],
            ),
            0.0,
            0.0,
            1.5,
        ),
        (
            "latent_oblique_b",
            Chart::latent_oblique(
                z0,
                [-0.7, 0.4, 1.0, -0.2, 0.55, 0.8, -0.35, 0.6],
                [0.25, 1.3, -0.15, 0.7, -0.9, 0.3, 0.85, -0.4],
            ),
            0.0,
            0.0,
            1.5,
        ),
        ("burrau_nu_k", Chart::BurrauFamily { nu_lo: 0.05, nu_hi: 0.95, k_max: 4.0, gamma_k: 1.5 }, 0.5, 0.5, 0.45),
        (
            "invariant_lz_k",
            Chart::Invariant { base: z0, k_max: 4.0, gamma_k: 1.5, report_e: false },
            0.5,
            0.5,
            0.45,
        ),
        (
            "mass_simplex",
            Chart::MassSimplex { z_alpha: 0.35, z_beta: -0.45, z_q: [0.25, -0.15, 0.40, 0.05], margin: 0.02 },
            0.5,
            0.5,
            0.45,
        ),
    ];

    println!(
        "budget {budget}, tau={tau:e}, alpha_hi={alpha_hi}, N=8, E+1={}, t={}, f64, {res}^2, \
         screen floor ON.\n\
         Colouring: hue = shape sphere by vMF site-blend, lightness = spread_shape on a log ramp\n\
         over each chart's own p1-p99. The window is printed because a false-colour image without\n\
         its scale is decoration.\n",
        ens.n_extra + 1,
        ens.t_max
    );
    println!(
        "{:>18} {:>14} {:>6} {:>7} {:>7} {:>6} {:>7} {:>9} {:>10} {:>10} {:>9} {:>6}",
        "case", "chart", "domain", "quads", "leaves", "depth", "screen", "distinct", "alpha med",
        "alpha idec", "ramp span", "bound"
    );

    let mut frames: Vec<Vec<u8>> = Vec::new();
    let mut wire_frames: Vec<Vec<u8>> = Vec::new();
    let mut control_ics: Option<Vec<prin_rs::physics::Cart<f64>>> = None;

    for (name, chart, cx, cy, half) in &cases {
        let (cx, cy, half) = (*cx, *cy, *half);
        if let Err(e) = chart.validate(0.0, cx, cy, half) {
            println!("{name:>18}  REFUSED: {e}");
            continue;
        }

        let cam = Camera::framing(cx, cy, half, res);
        let cfg = SchedCfg {
            budget,
            tau_display: tau,
            alpha_hi,
            alpha_lo: alpha_hi,
            agg: Agg::Median,
            chart: *chart,
            camera: Some(cam),
            keep_pixels: true,
            ..Default::default()
        };
        let (t, st) = scheduler::descend(cx, cy, half, 0, &cfg, &ens, Precision::F64);

        let leaves: Vec<usize> = t.leaves().collect();
        let depth = leaves.iter().map(|&i| t.nodes[i].level).max().unwrap_or(0);
        let screen = leaves
            .iter()
            .filter(|&&i| t.nodes[i].decision == Decision::ScreenFloor)
            .count();
        let alphas: Vec<f64> = leaves.iter().filter_map(|&i| t.nodes[i].alpha).collect();
        let (_, amed, _, aidec) = stats::interdecile(&alphas);

        // The colouring: one ramp and one site set for the whole chart, built from the pixels
        // this tree actually produced. Per-quad normalisation would make a quad's colour depend
        // on which quads happen to be leaves.
        let all_px: Vec<PixelOut> =
            leaves.iter().flat_map(|&i| st.pixels[i].iter().cloned()).collect();
        let (lo, hi) = colour::range(&all_px, Scalar::ShapeSpread);
        let (distinct, _, _) = colour::quantisation(&all_px, Scalar::ShapeSpread);
        let m_here = grid::decode_state(chart, 0, cx, cy).m;
        let sites = colour::landmarks(&m_here);
        let rgb = move |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);

        println!(
            "{:>18} {:>14} {:>6} {:>7} {:>7} {:>6} {:>7} {:>9} {:>10.4} {:>10.4} {:>9.3} {:>6}",
            name,
            chart.name(),
            if chart.domain() == Domain::Unit { "unit" } else { "free" },
            t.nodes.len(),
            leaves.len(),
            depth,
            screen,
            distinct,
            amed,
            aidec,
            hi / lo.max(f64::MIN_POSITIVE),
            // Whether the descent stopped because the CRITERION was satisfied or because the
            // BUDGET ran out. A budget-bound row's leaf count is a fact about the budget.
            if t.nodes.len() + 4 > budget { "BUDGET" } else { "crit" }
        );

        let stem = format!("results/charts/{name}");
        let img = render_leaves(&t, &st.pixels, &cam, res, &leaves, &rgb);
        let mut wimg = img.clone();
        let boxes = wire::boxes_from_tree(&t, &cam, res);
        let deepest = boxes.iter().map(|b| b.level).max().unwrap_or(1);
        wire::draw(&mut wimg, res, res, &boxes, deepest.max(1));
        let _ = adaptive::save(&format!("{stem}.png"), res, &img);
        let _ = adaptive::save(&format!("{stem}_wire.png"), res, &wimg);

        // **The chart itself, at one sample per pixel on a uniform grid.**
        //
        // The adaptive render above is a picture of the SCHEDULER: near-field's `alpha` median
        // is 0.14 against `alpha_hi = 0.2`, so the criterion says refinement does not pay and
        // keeps coarse leaves, and a coarse leaf is drawn as one flat tile because the render
        // never interpolates. That reads as blur and is an honest picture of an unrefined tree
        // -- but it means showing only the adaptive render is never showing the chart.
        //
        // So both. `_uniform` is what the chart looks like; the adaptive one and its wire twin
        // are what the scheduler made of it. Reading either alone is how a criterion's failure
        // gets mistaken for a rendering artefact, or a rendering choice for a finding.
        {
            // Full resolution: this is the sharpest artefact and the only one that shows the
            // chart rather than the tree, so it is the one worth paying for. One sample per
            // pixel, no interpolation anywhere.
            let ures = res;
            let usl = grid::Slice::body_plane(ures, ures, cx, cy, half, 0).with_chart(*chart);
            let upx: Vec<PixelOut> = (0..usl.npix())
                .into_par_iter()
                .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&usl, k, &ens))
                .collect();
            let (ulo, uhi) = colour::range(&upx, Scalar::ShapeSpread);
            let usites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(upx.len() * 3);
            for p in &upx {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &usites, ulo, uhi));
            }
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), ures, ures, &buf);
            let mut obuf = Vec::with_capacity(upx.len() * 3);
            for p in &upx {
                obuf.extend_from_slice(&png::outcome_rgb(p));
            }
            let _ =
                adaptive::save_rect(&format!("{stem}_uniform_outcome.png"), ures, ures, &obuf);
        }

        // The outcome control on the ADAPTIVE tree, so the pair says whether a feature is in the
        // physics or in the colouring. At t = 13 the outcome label is saturated, which is the
        // point.
        let (oimg, _) = adaptive::render(
            &t, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, png::outcome_rgb,
        );
        let _ = adaptive::save(&format!("{stem}_outcome.png"), res, &oimg);

        if let Ok(f) = std::fs::File::create(format!("{stem}.prnq")) {
            let mut w = std::io::BufWriter::new(f);
            let _ = prin_rs::output::tree::write(&mut w, &t, &cfg, &ens, &st, name, "f64");
        }

        // The level ladder: ONE descent, truncated at each depth.
        //
        // Rendered at the same resolution as the stills. Re-rasterised rather than
        // re-integrated, so it costs nothing but the raster.
        let ares = res;
        let acam = Camera::framing(cx, cy, half, ares);
        let aboxes = wire::boxes_from_tree(&t, &acam, ares);
        let mut ladder: Vec<Vec<u8>> = Vec::new();
        let mut wladder: Vec<Vec<u8>> = Vec::new();
        for cap in 0..=depth {
            let lv = leaves_capped(&t, cap);
            let f = render_leaves(&t, &st.pixels, &acam, ares, &lv, &rgb);
            let mut wf = f.clone();
            let b: Vec<wire::Box2> = aboxes.iter().cloned().filter(|b| b.level <= cap).collect();
            wire::draw(&mut wf, ares, ares, &b, cap.max(1));
            ladder.push(f);
            wladder.push(wf);
        }
        let _ = apng::write(&format!("{stem}_levels.png"), ares, ares, &ladder, 1, 2);
        let _ = apng::write(&format!("{stem}_levels_wire.png"), ares, ares, &wladder, 1, 2);

        // The control: `plane_00deg` is `body_plane` written a second way. Compared on INITIAL
        // CONDITIONS, which is exact -- comparing images conflates "same chart" with "the
        // rasteriser rounds the same way at O(1) and O(0) coordinate magnitudes".
        let sl = grid::Slice::body_plane(64, 64, cx, cy, half, 0).with_chart(*chart);
        let ics: Vec<prin_rs::physics::Cart<f64>> =
            (0..sl.npix()).map(|k| sl.nominal::<f64>(k)).collect();
        if *name == "body_plane" {
            control_ics = Some(ics);
        } else if *name == "plane_00deg" {
            let a = control_ics.as_ref().expect("body_plane must run first");
            let d = a
                .iter()
                .zip(ics.iter())
                .map(|(x, y)| decode::max_abs_diff(x, y))
                .fold(0.0f64, f64::max);
            assert_eq!(
                d, 0.0,
                "CONTROL FAILED: plane_00deg is not bitwise body_plane (max |dIC| = {d:e}). \
                 The bases are wrong and every other row compares different physics."
            );
            println!(
                "{:>18}  [control] plane_00deg vs body_plane: max |dIC| = {d:e} -- the same \
                 chart, asserted",
                ""
            );
        }

        frames.push(img);
        wire_frames.push(wimg);
    }

    let _ = apng::write("results/charts/gallery.png", res, res, &frames, 1, 1);
    let _ = apng::write("results/charts/gallery_wire.png", res, res, &wire_frames, 1, 1);
    println!(
        "\n{} charts: still + wire twin + outcome control + level ladder (both) + .prnq each,\n\
         plus the two gallery APNGs. Everything at {res}^2.",
        frames.len(),
    );

    println!(
        "\n\
         `distinct` is how many distinct values the lightness field takes over the chart. Read it\n\
         before the picture: a field with few distinct values has that many colours in it, and no\n\
         ramp recovers what is not there.\n\
         \n\
         `ramp span` is hi/lo of the p1-p99 window. The ramp is auto-ranged per chart, so a chart\n\
         with no dynamic range has its NOISE stretched to full scale and reads as structure. A\n\
         span near 1 means the picture is of the estimator, not the physics.\n\
         \n\
         `screen` is leaves stopped by the screen floor rather than by the criterion. It is a\n\
         VETO on scale and never a trigger, and it is view-relative -- a floored quad refines\n\
         again when zoomed into.\n\
         \n\
         Leaf counts are slice-conditional to 4.3x and compare WITHIN a chart, never across.\n\
         `alpha med` and `alpha idec` are the cross-chart quantities. The interdecile rather than\n\
         the variance: excess kurtosis on alpha_shape is 110, so the variance describes the tail.\n\
         \n\
         READ THE `bound` COLUMN BEFORE ANY LEAF COUNT. A row marked BUDGET stopped because the\n\
         budget ran out, not because the criterion was satisfied -- its quads, leaves and depth\n\
         are facts about the budget and are the same for every such row by construction. Only a\n\
         `crit` row's tree shape is a statement about the chart.\n\
         \n\
         And the substantive finding this table carries: the reference's chart families, centred\n\
         where they are centred here, are TAME. `alpha med` sits at 0.99-1.01 on every latent,\n\
         Burrau and simplex row against 0.14 for `body_plane` and 0.19 for the shape sphere.\n\
         alpha near 1 means splitting halves the spread -- refinement pays, so the scheduler\n\
         refines everywhere and runs to the budget. That is correct behaviour on a tame region,\n\
         not a scheduler fault.\n\
         \n\
         But it means these charts are not exercising the criterion where it is hard. Tameness is\n\
         a property of WHERE a chart is centred, and the base latent point here was chosen to\n\
         avoid the z = 0 symmetry rather than to find chaos. Moving it until the picture gets\n\
         interesting would be tuning; the honest report is that the interesting question -- does\n\
         the criterion behave consistently across charts -- is not answered by a set of charts\n\
         that are all tame."
    );
}
