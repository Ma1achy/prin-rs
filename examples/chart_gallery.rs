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
use prin_rs::grid::{self, Domain};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, png, wire};
use prin_rs::quad::{Agg, Criterion, Decision, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::{decode, logln, stats};
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
    // The leaf set is the argument. This used to build a shadow tree AND empty the samples of
    // every node outside the set, because the render keyed painting on "has samples" -- and
    // that key was also the coarse-ancestor fill, so no truncated frame could show one. Both
    // live in `adaptive::render_leaves` now: the set and its ancestors, nothing else.
    adaptive::render_leaves(t, pixels, cam, res, adaptive::TexelMode::Adaptive, |p| rgb(p), leaves).0
}

fn main() {
    // Budget high enough that the descent stops on the CRITERION rather than on the budget.
    // At 4000 every chart but `body_plane` and `shape_sphere` hit the cap, which made their
    // leaf counts a fact about the budget and left large areas as coarse leaves -- and a coarse
    // leaf is drawn as one flat tile, because the render never interpolates. That reads as blur
    // and is not: it is an honest picture of an unrefined tree.
    let budget: usize = arg(1, 40000);
    // **The struct's default is the one default.** These read `1e-4` and `0.2` here while
    // `SchedCfg::default()` said `1e-2` and `0.5` -- two defaults for one knob, and every
    // committed tree was cut at the argument's value. The committed corpus names its arguments
    // in `results/charts/README.md`; a bare run now means the shipped configuration.
    let tau: f64 = arg(2, SchedCfg::default().tau_display);
    let alpha_hi: f64 = arg(3, SchedCfg::default().alpha_hi);
    // **`alpha_lo` was TIED to `alpha_hi` here, which is a `Policy::Alpha`-era coupling.**
    //
    // Under that policy they were the two ends of one band -- split above `alpha_hi`, floor below
    // `alpha_lo`, keep between -- and setting them equal collapsed the band to a single threshold,
    // which is what the corpus was measured at. Under `Policy::Tolerance` they are **different
    // mechanisms**: `alpha_hi` is inert (the split test is the tolerance, not an exponent) and
    // `alpha_lo` is the **area floor's dimension threshold**, `alpha_area = 2 - d`. Carrying the
    // coupling into a tolerance run silently sets that floor to `alpha_hi` -- the documented
    // command's `0.2`, or the argument default's **0.5** -- against its own measured default of
    // `0.005`, where `0.2` alone costs 11% of `config_stability`'s resolvable pixels and puts its
    // tree ABOVE uniform at its own error.
    //
    // So it is argument 11, defaulting to the struct's value. The `Policy::Alpha` reproduction
    // passes it explicitly, exactly as `k_frac = 1.0` does for the unranked control.
    let alpha_lo: f64 = arg(11, SchedCfg::default().alpha_lo);
    let res: usize = arg(4, 1024);
    // **The knob that made the whole committed gallery a uniform-mode render.** `k_frac = 1`
    // takes the top 100% of the frontier, so the ranking runs and changes nothing. It was the
    // silent default when the first `results/charts` was made; it is now an argument with the
    // shipped value as its default, and passing `1.0` writes to `charts_unranked` instead --
    // so the control cannot land on top of the corpus.
    let k_frac: f64 = arg(5, scheduler::K_FRAC_RANKED);
    let crit = std::env::args()
        .nth(6)
        .map(|c| Criterion::parse(&c).expect("criterion"))
        .unwrap_or(Criterion::Within);
    // **The naming is inverted from what it was, because the split it encoded is gone.**
    // `charts_ranked` was the *after* of the `k_frac` change; `k_frac = 0.25` has been the
    // shipped default since PR #21, so the after IS the corpus and the canonical name should
    // hold the canonical run. The unranked arm keeps a name that says what it is.
    //
    // The old `results/charts` and `results/charts_ranked` were both written 25-26 August and
    // are superseded by every integrator fix from 27 August on; they are recoverable at
    // `9d48510` and are not preserved under a third name here. `results/README.md` says so.
    let ranked = k_frac < scheduler::K_FRAC_UNRANKED;
    // **An output root is an argument, not a constant.** `criterion_metric` was fixed for exactly
    // this and `chart_gallery` was not: with the root hardcoded, a reduced-`res` validation pass
    // -- the only way to check whether a flag is inert before spending hours on the real run --
    // overwrites the committed 1024^2 corpus with a small raster, and *softness in an image is a
    // raster size* reads it back as a rendering fault rather than a stale file. Third site.
    let root: String = std::env::args().nth(7).unwrap_or_else(|| "results".into());
    let dir = if ranked { format!("{root}/charts") } else { format!("{root}/charts_unranked") };
    let adir = if ranked { format!("{root}/animated") } else { format!("{root}/animated_unranked") };
    let dir = dir.as_str();
    let adir = adir.as_str();
    let _ = std::fs::create_dir_all(dir);
    let _ = std::fs::create_dir_all(adir);

    // **`refine_flagged` is argument 8, and its default is production's.** The hardcoded `false`
    // that stood here is the line the record names as spread-by-copy out of the experiment
    // harnesses it was correct in, into render harnesses it was never argued for -- while
    // `results/README.md` asserted renders had it on. It is a real choice with two defensible
    // readings (the repair pass has no live-playhead analogue, so a scheduler corpus arguably
    // wants it off; the standing invariant says renders have it on), so it is a *named argument*
    // that the provenance sidecar records, not a constant nothing prints.
    let refine: bool = std::env::args()
        .nth(8)
        .map(|v| v == "1" || v == "true")
        .unwrap_or(EnsembleCfg::production().refine_flagged);
    let ens = EnsembleCfg { refine_flagged: refine, ..EnsembleCfg::production() };
    // **The uniform panels are argument 9, default off.** `<case>_uniform*.png` is the chart at
    // one sample per pixel -- 8.4M trajectories per chart at 1024^2, about 95% of a run -- and
    // it used to be skipped whenever the tree was ranked, on the argument that a scheduler
    // change cannot move it. True, and the PHYSICS moved: the committed `_uniform*` panels were
    // 25 August beside adaptive twins from 3 September, mirror-imaged and on a different colour
    // window. When asked for, the grid is evaluated FIRST and its window colours both panels.
    let uniform: bool =
        std::env::args().nth(9).map(|v| v == "1" || v == "true").unwrap_or(false);
    // **Argument 10: which charts, comma-separated; `all` or absent runs the gallery.** So a
    // regeneration can be staged a chart at a time -- under the tolerance policy a chart can
    // cost minutes to hours -- rather than committed to as twenty-six at once.
    let only: Option<Vec<String>> = std::env::args()
        .nth(10)
        .filter(|s| s != "all")
        .map(|s| s.split(',').map(|x| x.trim().to_string()).collect());
    // A tree under `results/` on any kernel but production's is the superseded corpus again.
    scheduler::assert_production_kernel(&ens, dir);
    let log = prin_rs::output::Log::tee(&format!("{root}/output/chart_gallery.txt"));
    let log = &log;
    // **The column, not the instance.** Nine harnesses feeding the refinement work printed no
    // provenance at all -- the `refine_flagged` failure exactly: *the failure was never the
    // choice, it is that nothing recorded the choice.*
    logln!(log, "  config: {}", ens.provenance());
    logln!(log, "  uniform panels: {}", if uniform { "ON (8.4M trajectories per chart at 1024^2)" }
                                          else { "off -- pass 1 as argument 9 to regenerate them" });


    // A base latent point. Deliberately not the origin: at z = 0 every sigmoid sits at 0.5 and
    // several coordinates would be at a symmetry point, which is exactly where a sign error
    // hides.
    let cases = grid::gallery_cases();


    logln!(log, 
        "budget {budget}, tau={tau:e}, alpha_hi={alpha_hi}, N=8, E+1={}, t={}, f64, {res}^2, \
         screen floor ON.\n\
         Colouring: hue = shape sphere by vMF site-blend, lightness = spread_shape on a log ramp\n\
         over each chart's own p1-p99. The window is printed because a false-colour image without\n\
         its scale is decoration.\n",
        ens.n_extra + 1,
        ens.t_max
    );
    logln!(log, 
        "{:>18} {:>14} {:>6} {:>7} {:>7} {:>6} {:>7} {:>9} {:>10} {:>10} {:>9} {:>9}",
        "case", "chart", "domain", "quads", "leaves", "depth", "screen", "distinct", "alpha med",
        "alpha idec", "ramp span", "bound"
    );

    let mut frames: Vec<Vec<u8>> = Vec::new();
    let mut wire_frames: Vec<Vec<u8>> = Vec::new();
    let mut control_ics: Option<Vec<prin_rs::physics::Cart<f64>>> = None;

    for (name, chart, cx, cy, half) in &cases {
        let (cx, cy, half) = (*cx, *cy, *half);
        if let Some(list) = &only {
            if !list.iter().any(|n| n == name) {
                continue;
            }
        }
        if let Err(e) = chart.validate(0.0, cx, cy, half) {
            logln!(log, "{name:>18}  REFUSED: {e}");
            continue;
        }

        let cam = Camera::framing(cx, cy, half, res);
        let cfg = SchedCfg {
            budget,
            tau_display: tau,
            alpha_hi,
            alpha_lo,
            agg: Agg::Median,
            chart: *chart,
            camera: Some(cam),
            keep_pixels: true,
            criterion: crit,
            k_frac,
            ..Default::default()
        };
        // Refuses only when this run would write a headline artefact from the uniform-mode
        // control. `ranked` is the deliberate-control flag: at `k_frac = 1` the destination is
        // `results/charts`, which IS the before, and the guard must not block reproducing it.
        scheduler::assert_not_uniform_in_disguise(&cfg, dir, !ranked);
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
        // The uniform grid, when asked for, is evaluated here so that ONE window colours the
        // adaptive panel and the uniform panel beside it. Two auto-ranges made the pair
        // incomparable pixel for pixel.
        let upx: Option<Vec<PixelOut>> = if uniform {
            let usl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(*chart);
            Some(
                (0..usl.npix())
                    .into_par_iter()
                    .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&usl, k, &ens))
                    .collect(),
            )
        } else {
            None
        };
        let (lo, hi) = match &upx {
            Some(u) => colour::range(u, Scalar::ShapeSpread),
            None => colour::range(&all_px, Scalar::ShapeSpread),
        };
        let (distinct, _, _) = colour::quantisation(&all_px, Scalar::ShapeSpread);
        let m_here = grid::decode_state(chart, 0, cx, cy).m;
        let sites = colour::landmarks(&m_here);
        let rgb = move |p: &PixelOut| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi);

        logln!(log, 
            "{:>18} {:>14} {:>6} {:>7} {:>7} {:>6} {:>7} {:>9} {:>10.4} {:>10.4} {:>9.3} {:>9}",
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
            // **What actually stopped the descent.** This read `crit` unless the BUDGET ran
            // out, which is wrong and was hiding the largest fact about this table: on most
            // charts the tree is set by `Decision::MaxRelDepth`, a CAMERA VETO, and the leaf
            // count is a fact about the cap rather than about the criterion. A row stopped by a
            // veto has not exercised the criterion at all, and its `alpha` describes quads the
            // cap forced rather than quads the criterion chose. Same lesson as the screen floor.
            {
                let code_of = |d: Decision| {
                    leaves.iter().filter(|&&i| t.nodes[i].decision == d).count()
                };
                let (cap, floor_n, keep) = (
                    code_of(Decision::MaxRelDepth) + code_of(Decision::ScreenFloor),
                    code_of(Decision::Floor),
                    code_of(Decision::Keep),
                );
                let n = leaves.len().max(1);
                if t.nodes.len() + 4 > budget {
                    "BUDGET".to_string()
                } else if 100 * cap / n >= 50 {
                    format!("VETO {}%", 100 * cap / n)
                } else {
                    format!("crit {}%", 100 * (floor_n + keep) / n)
                }
            }
        );

        // ---- The mechanism test: leaf depth against `terminated_fraction` ------------------
        //
        // The finding that refinement goes to smooth regions and ignores the filaments was read
        // off a **wireframe at the wrong window**. A wireframe is an appearance; this tests the
        // proposed *cause* directly, and it either shows or it does not.
        //
        // The mechanism: terminated regions are absorbing, so nearby copies share an outcome,
        // `spread_event` collapses to zero and the criterion sees a resolved quad. Still-running
        // regions keep diverging and hold high spread forever. If that is what drives the tree,
        // leaf depth should be **anti-correlated** with `terminated_fraction` -- refinement
        // chasing non-convergence rather than structure.
        //
        // Both quantities are already in the `PRNQ` dump, so this costs a plot and no
        // integration. Read the per-depth median and interdecile, never the variance: the
        // standing rule on `alpha` (excess kurtosis 110) is about this same family of per-quad
        // statistics, and a mean here would be a statement about the tail.
        {
            let pts: Vec<(f64, f64)> = leaves
                .iter()
                .map(|&i| (t.nodes[i].level as f64, t.nodes[i].red.terminated_fraction))
                .collect();
            let (xs, ys): (Vec<f64>, Vec<f64>) = pts.iter().cloned().unzip();
            let rho = stats::spearman(&xs, &ys);
            let mut med: Vec<(f64, f64)> = Vec::new();
            let mut rows: Vec<String> = Vec::new();
            for lv in 0..=depth {
                let v: Vec<f64> = pts
                    .iter()
                    .filter(|&&(x, y)| x as u32 == lv && y.is_finite())
                    .map(|&(_, y)| y)
                    .collect();
                if v.is_empty() {
                    continue;
                }
                let (p10, m, p90, _) = stats::interdecile(&v);
                med.push((lv as f64, m));
                rows.push(format!("L{lv}:n={} med={m:.3} [{p10:.3},{p90:.3}]", v.len()));
            }
            // `escape_fraction` separately, because `t_end` termination is NOT escape and
            // conflating them contradicts a standing result while appearing to agree with it.
            let esc: f64 = leaves
                .iter()
                .map(|&i| t.nodes[i].red.escape_fraction)
                .sum::<f64>()
                / leaves.len().max(1) as f64;
            logln!(log, 
                "{:>18}  depth~terminated_fraction: spearman = {rho:+.4}, mean escape_fraction = \
                 {esc:.4}\n{:>20}{}",
                "",
                "",
                rows.join("  ")
            );
            let sc = prin_rs::output::plot::Scatter {
                title: format!("{name}: leaf depth against terminated_fraction"),
                x_label: "leaf depth".into(),
                y_label: "terminated_fraction".into(),
                notes: vec![
                    format!(
                        "{} leaves, spearman(depth, terminated) = {rho:+.4}, \
                         mean escape_fraction = {esc:.4}",
                        leaves.len()
                    ),
                    "Anti-correlation is the prediction: the criterion chases non-convergence, \
                     not structure. Terminated regions are absorbing, so copies agree, \
                     spread_event collapses and the quad reads resolved."
                        .into(),
                    "Orange is the per-depth MEDIAN. Read the cloud's spread, not a variance."
                        .into(),
                ],
                points: pts,
                overlay: med,
            };
            let _ = sc.save(&format!("{dir}/{name}_termdepth"));
        }

        let stem = format!("{dir}/{name}");
        let img = render_leaves(&t, &st.pixels, &cam, res, &leaves, &rgb);
        let mut wimg = img.clone();
        let boxes = wire::boxes_from_tree(&t, &cam, res);
        let deepest = boxes.iter().map(|b| b.level).max().unwrap_or(1);
        wire::draw(&mut wimg, res, res, &boxes, deepest.max(1));
        let _ = adaptive::save(&format!("{stem}.png"), res, &img);
        let _ = adaptive::save(&format!("{stem}_wire.png"), res, &wimg);
        // **The PNGs were the blind spot.** The `.prnq` carries a settings header; the panels
        // carried nothing, so a picture could not say which integrator drew it. One sidecar per
        // chart rather than per frame: the frames of a chart share a config by construction, and
        // 208 near-identical files would be noise rather than provenance.
        let _ = prin_rs::output::provenance_sidecar(
            &format!("{stem}.png"),
            &ens,
            &format!(
                "chart={} leaves={} depth={} stop={} scalar=ShapeSpread window=({lo:.4e},{hi:.4e}) \
                 window_from={} res={res} viewport={res} budget={budget} tau_display={tau:e} \
                 alpha_hi={alpha_hi} criterion={} k_frac={k_frac}\n",
                chart.name(),
                leaves.len(),
                depth,
                t.stop_breakdown(),
                if upx.is_some() { "uniform_grid" } else { "tree_leaves" },
                crit.name(),
            ),
        );

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
        //
        // **And it is the one block that CANNOT depend on `k_frac`.** It builds its own
        // `res x res` slice and evaluates it directly; no quad, no tree, no decision enters it.
        // So a scheduler change cannot move `results/charts/*_uniform*.png` -- but a PHYSICS
        // change does, and this block was skipped whenever the tree was ranked on the strength
        // of the first fact alone. It costs `res^2 * (E+1)` trajectories, 8.4M per chart at
        // 1024, so it is argument 9 and off by default; the grid itself was evaluated above,
        // before the adaptive render, so the two panels share a window.
        if let Some(upx) = upx.as_deref() {
            // Full resolution: this is the sharpest artefact and the only one that shows the
            // chart rather than the tree, so it is the one worth paying for. One sample per
            // pixel, no interpolation anywhere.
            let ures = res;
            let usites = colour::landmarks(&m_here);
            let mut buf = Vec::with_capacity(upx.len() * 3);
            for p in upx {
                buf.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &usites, lo, hi));
            }
            let _ = prin_rs::output::provenance_sidecar(
                &format!("{stem}_uniform.png"),
                &ens,
                &format!(
                    "chart={} panel=uniform scalar=ShapeSpread window=({lo:.4e},{hi:.4e}) \
                     window_from=uniform_grid res={ures} one_sample_per_pixel=true\n",
                    chart.name()
                ),
            );
            let _ = adaptive::save_rect(&format!("{stem}_uniform.png"), ures, ures, &buf);
            let mut obuf = Vec::with_capacity(upx.len() * 3);
            for p in upx {
                obuf.extend_from_slice(&png::outcome_rgb(p));
            }
            let _ =
                adaptive::save_rect(&format!("{stem}_uniform_outcome.png"), ures, ures, &obuf);

            // **The event-class panel, on viridis** -- the mode the reference's WebGPU panel
            // renders (`Colour mode: Event class, Palette: viridis`). A reference comparison
            // has to be made under a matched mode: a continuous field and a categorical map
            // cannot look alike even when both are correct, and comparing across modes is how a
            // rendering choice gets mistaken for a physics bug. That is most of what went wrong
            // with this port.
            //
            // The outcome panel above **stays as the control**. The pair says whether a feature
            // is in the class definition or in the physics, and the two classes differ: the
            // event class is the currently-tightest pair joined with the terminal outcome and
            // is defined at every playhead, where the outcome label at t = 13 is saturated.
            let mut ebuf = Vec::with_capacity(upx.len() * 3);
            for p in upx {
                ebuf.extend_from_slice(&png::event_class_rgb(p));
            }
            let _ = adaptive::save_rect(&format!("{stem}_uniform_event.png"), ures, ures, &ebuf);

            // The histogram, before the image. 27 slots on one ramp means adjacent classes are
            // close in colour by construction, so the legend and the counts are the instrument.
            // A class that never fires is a fact about the slice and reads as a zero here.
            let (rows, undet) = png::event_class_histogram(upx);
            let live: Vec<String> = rows
                .iter()
                .filter(|&&(_, n)| n > 0)
                .map(|&(c, n)| format!("{}={n}", png::event_class_name(c)))
                .collect();
            logln!(log, 
                "{:>18}  event classes ({} of {} fire, {undet} undetermined): {}",
                "",
                live.len(),
                png::N_EVENT_CLASSES,
                if live.is_empty() { "none".into() } else { live.join(", ") }
            );
        }

        // The outcome control on the ADAPTIVE tree, so the pair says whether a feature is in the
        // physics or in the colouring. At t = 13 the outcome label is saturated, which is the
        // point.
        let (oimg, _) = adaptive::render(
            &t, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, png::outcome_rgb,
        );
        let _ = adaptive::save(&format!("{stem}_outcome.png"), res, &oimg);

        // The same matched-mode panel on the adaptive tree.
        let (eimg, _) = adaptive::render(
            &t, &st.pixels, &cam, res, adaptive::TexelMode::Adaptive, png::event_class_rgb,
        );
        let _ = adaptive::save(&format!("{stem}_event.png"), res, &eimg);

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

        let mut ladder: Vec<Vec<u8>> = Vec::new();
        let mut wladder: Vec<Vec<u8>> = Vec::new();
        for cap in 0..=depth {
            let lv = leaves_capped(&t, cap);
            let f = render_leaves(&t, &st.pixels, &acam, ares, &lv, &rgb);
            let mut wf = f.clone();
            // From the capped leaf set, not from a level filter over the finished tree's boxes.
            // Filtering by level keeps every deep quad whose level happens to be <= cap while
            // dropping nothing that the cap actually removed, so the wire drifted out of step
            // with the colour frame beside it.
            let b = wire::boxes_from_leaves(&t, &acam, ares, &lv);
            // The FINISHED tree's depth grades every frame. `cap.max(1)` regraded each frame,
            // so a level-3 box was bright in frame 3 and dim in frame 6 -- the ramp moving
            // rather than the tree, the fault the colour frames avoid by holding one window.
            wire::draw(&mut wf, ares, ares, &b, depth.max(1));
            ladder.push(f);
            wladder.push(wf);
        }
        // **Animations live in `results/animated/`, not beside the stills.** They are the only
        // artefacts here that show *how the tree got there* rather than what it settled on, and
        // 72 of them scattered through three directories were effectively unfindable.
        let anim = format!("{adir}/{name}");
        let _ = apng::write(&format!("{anim}_levels.png"), ares, ares, &ladder, 1, 2);
        let _ = apng::write(&format!("{anim}_levels_wire.png"), ares, ares, &wladder, 1, 2);
        // The duplicate-frame count, printed and kept: a ladder whose frames repeat is a still.
        // And a sidecar per animation, because a frame carries no header of its own.
        let (dup, wdup) =
            (apng::adjacent_duplicates(&ladder), apng::adjacent_duplicates(&wladder));
        for (suffix, d) in [("levels", dup), ("levels_wire", wdup)] {
            let _ = prin_rs::output::provenance_sidecar(
                &format!("{anim}_{suffix}.png"),
                &ens,
                &format!(
                    "chart={} animation={suffix} frames={} adjacent_duplicates={d} \
                     scalar=ShapeSpread window=({lo:.4e},{hi:.4e}) res={ares} viewport={ares} \
                     budget={budget} tau_display={tau:e} alpha_hi={alpha_hi} criterion={} \
                     k_frac={k_frac} stop={}\n",
                    chart.name(), ladder.len(), crit.name(), t.stop_breakdown()
                ),
            );
        }
        if dup > 0 || wdup > 0 {
            logln!(log, "{:>18}  levels ladder: {dup}/{wdup} identical adjacent frame pairs of {} \
                          (colour/wire)", "", ladder.len() - 1);
        }

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
            logln!(log, 
                "{:>18}  [control] plane_00deg vs body_plane: max |dIC| = {d:e} -- the same \
                 chart, asserted",
                ""
            );
        }

        frames.push(img);
        wire_frames.push(wimg);
    }

    let _ = apng::write(&format!("{adir}/gallery.png"), res, res, &frames, 1, 1);
    let _ = apng::write(&format!("{adir}/gallery_wire.png"), res, res, &wire_frames, 1, 1);
    for (suffix, fr) in [("gallery", &frames), ("gallery_wire", &wire_frames)] {
        let _ = prin_rs::output::provenance_sidecar(
            &format!("{adir}/{suffix}.png"),
            &ens,
            &format!(
                "animation={suffix} one_frame_per_chart frames={} adjacent_duplicates={} \
                 res={res} budget={budget} tau_display={tau:e} alpha_hi={alpha_hi} \
                 criterion={} k_frac={k_frac}\n",
                fr.len(), apng::adjacent_duplicates(fr), crit.name()
            ),
        );
    }
    logln!(log, 
        "\n{} charts: still + wire twin + outcome control + level ladder (both) + .prnq each,\n\
         plus the two gallery APNGs. Everything at {res}^2.",
        frames.len(),
    );

    logln!(log, 
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
         where they are centred here, are TAME in `alpha` -- but the numbers this paragraph\n\
         used to quote were measured on the PRE-FIX kernel and no longer describe the table\n\
         above it. `body_plane` read 0.14 and the shape sphere 0.19; they now read 1.02 and\n\
         1.26, and no row runs to the budget at all. Read the printed columns, not this text.\n\
         alpha near 1 means splitting halves the spread -- refinement pays, so the scheduler\n\
         refines everywhere. On the fixed kernel it does NOT run to the budget: the spread fell\n\
         with the integration failure that inflated it, so at a fixed `tau` the quads now read\n\
         resolved and stop as `Keep`. That is correct behaviour on a tame region,\n\
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
