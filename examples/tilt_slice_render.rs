//! **`tilt_plambda` — the resolution ladder and the two time animations.**
//!
//! The first slice in the project with a non-zero tilt. `tests/tilt_slice.rs` holds what is
//! checkable without integrating: the two dead slots at exactly `0.000e+00`, the unequal masses
//! `(0.35333, 0.27523, 0.37144)` against the presets' `1/3` each, the orthonormal frame turned by
//! exactly `amt + gamma`, and the fact that this tilt is **in-plane** — the same 2-plane as
//! `preset_plambda` with the frame rotated, which is a different *slice* and not a different
//! *plane*.
//!
//! Run: `cargo run --release --example tilt_slice_render [root] [anim_res] [t_max] [stride]`
//!
//! # What it writes, under `<root>/tilt/`
//!
//! - `tilt_plambda_<res>_uniform.png`, `_outcome.png`, `_event.png` at 64, 256, 512 and 1024.
//! - `tilt_plambda_uniform_time.png` — an APNG, one frame per recorded sync boundary, each a
//!   **full uniform render** of the grid as it stood at that playhead. No tree.
//!
//! The refinement-over-time animation is `examples/live_animation.rs`, which already does exactly
//! that and reaches this slice through `grid::named_slice`. It is not duplicated here.
//!
//! # One colour window for the whole ladder
//!
//! Every panel takes the p1–p99 window of the **1024²** pass. *An auto-ranged ramp cannot tell
//! "no signal" from "signal"*, and a ramp re-ranged per panel would stretch each raster's own
//! quantiles to full scale — on a question about what changes with resolution that manufactures
//! the answer. The window is printed and lands in every sidecar. The animation carries its own
//! fixed window for the same reason, computed once over **every** frame rather than per frame.
//!
//! # `Veto::None`
//!
//! A presentation render does not colour a pixel by a debug flag. `n_nonfinite` counts copies the
//! *driver* could not use; the nominal `shape_vec` is finite on the flagged set and
//! `spread_shape` is an ordinary number over the copies that ran. The count is printed and named
//! in every sidecar instead.
//!
//! Memory: the 1024² pass is 1,048,576 footprints. They are held as [`PixelSlim`] (48 bytes
//! against `PixelOut`'s 656) so the ladder fits, which is the same route `chart_gallery` takes
//! and is pinned bitwise by `tests/pixel_slim.rs`.

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut, PixelSlim};
use prin_rs::grid::{self, Chart};
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::{adaptive, apng, png, Log};
use prin_rs::scheduler;
use prin_rs::logln;
use rayon::prelude::*;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// Evaluate a uniform grid strip by strip, keeping only the six fields the panels read.
fn uniform_slim(sl: &grid::Slice, ens: &EnsembleCfg, strip: usize) -> Vec<PixelSlim> {
    let mut out = Vec::with_capacity(sl.npix());
    let mut row = 0usize;
    while row < sl.ny {
        let hi = (row + strip).min(sl.ny);
        let mut band: Vec<PixelSlim> = (row * sl.nx..hi * sl.nx)
            .into_par_iter()
            .map(|k| PixelSlim::thin(&pixel::evaluate::<f64>(sl, k, ens)))
            .collect();
        out.append(&mut band);
        row = hi;
    }
    out
}

fn main() {
    let root: String = arg(1, "results".into());
    let anim_res: usize = arg(2, 256);
    let t_max: f64 = arg(3, EnsembleCfg::production().t_max);
    let stride: usize = arg(4, 1);

    let dir = format!("{root}/tilt");
    let _ = std::fs::create_dir_all(&dir);
    let _ = std::fs::create_dir_all(format!("{root}/output"));
    let log = Log::tee(&format!("{root}/output/tilt_slice_render.txt"));
    let log = &log;

    let (chart, cx, cy, half) = Chart::tilt_plambda();
    let base = EnsembleCfg::production();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;

    // **The stills take the production kernel unmodified**, so `assert_production_kernel` is
    // free to refuse a `results/` path if it ever stops being production.
    let ens = EnsembleCfg { t_max, n_sync, ..base.clone() };
    scheduler::assert_production_kernel(&ens, &dir);

    // Argument five overrides the ladder for a smoke pass; the default is the shipped one.
    let ladder: Vec<usize> = std::env::args().nth(5)
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![1024usize, 512, 256, 64]);
    assert!(!ladder.is_empty(), "empty resolution ladder");

    logln!(log, "tilt_plambda — the first slice with a non-zero tilt.");
    logln!(log, "config: {}", ens.provenance());
    logln!(log, "reproduce: cargo run --release --example tilt_slice_render {root} {anim_res} {t_max} {stride} {}",
           ladder.iter().map(|r| r.to_string()).collect::<Vec<_>>().join(","));
    logln!(log, "window: cx {cx:.17} cy {cy:.17} half {half:.17}   (pan is the window CENTRE: cx = 2*pan-1, half = zoom)");
    if let Chart::Latent { z0, q1, q2 } = chart {
        let (m, _) = prin_rs::physics::decoder::masses(z0.z_mu);
        logln!(log, "masses: {:.5} {:.5} {:.5}   <- UNEQUAL. 1/3 each here means the decode is being overridden.",
               m[0], m[1], m[2]);
        logln!(log, "q1: {q1:?}");
        logln!(log, "q2: {q2:?}");
        let leak: f64 = (0..8).filter(|k| *k != 4 && *k != 5)
            .map(|k| q1[k].abs().max(q2[k].abs())).fold(0.0, f64::max);
        logln!(log, "component outside span(e4,e5): {leak:.3e}  -> the tilt is IN-PLANE: \
                     preset_plambda's plane, frame turned by {:.4} deg.", q1[5].atan2(q1[4]).to_degrees());
    }
    logln!(log, "");

    // ---- the resolution ladder -------------------------------------------------------------
    //
    // Descending, so the finest pass runs first and its window governs the whole ladder.
    //
    // **One pass resident at a time.** The first cut computed all four into a `Vec` and painted
    // at the end; it was **OOM-killed** inside the 1024² pass, holding that raster plus a strip
    // of full `PixelOut` plus everything already finished. Each pass is now painted and dropped
    // before the next begins, so the peak is the finest raster alone and the only thing carried
    // between passes is the ramp window, which is two numbers. The strip is 16 rows rather than
    // 64 for the same reason: a `PixelOut` is 656 bytes against `PixelSlim`'s 48, so the
    // transient dominates the resident set at 1024 wide.
    let finest = ladder[0];
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);
    let mut window: Option<(f64, f64)> = None;

    for &res in ladder.iter() {
        let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
        let t0 = std::time::Instant::now();
        let px = uniform_slim(&sl, &ens, 16);
        logln!(log, "{res:>5}^2 uniform: {} footprints in {:.1}s", px.len(), t0.elapsed().as_secs_f64());
        let (lo, hi) = *window
            .get_or_insert_with(|| colour::range_q_of(px.iter().map(|p| p.spread_shape), 0.01, 0.99));
        if res == finest {
            logln!(log, "shared ramp window from the {finest}^2 pass (the finest): ({lo:.4e}, {hi:.4e})");
        }
        let px = &px;
        let res = &res;
        let vetoed = px
            .iter()
            .filter(|p| colour::vetoed(&p.fat(), Scalar::ShapeSpread, &sites, lo, hi))
            .count();
        let stem = format!("{dir}/tilt_plambda_{res}");
        let mut u = Vec::with_capacity(px.len() * 3);
        let mut o = Vec::with_capacity(px.len() * 3);
        let mut e = Vec::with_capacity(px.len() * 3);
        for p in px.iter().map(|p| p.fat()) {
            u.extend_from_slice(&colour::rgb_veto(&p, Scalar::ShapeSpread, &sites, lo, hi, colour::Veto::None));
            o.extend_from_slice(&png::outcome_rgb_veto(&p, colour::Veto::None));
            e.extend_from_slice(&png::event_class_rgb_veto(&p, colour::Veto::None));
        }
        for (suffix, buf, what) in [
            ("uniform", &u, "scalar=ShapeSpread"),
            ("uniform_outcome", &o, "colouring=outcome_class"),
            ("uniform_event", &e, "colouring=event_class_viridis"),
        ] {
            let path = format!("{stem}_{suffix}.png");
            let _ = adaptive::save_rect(&path, *res, *res, buf);
            let _ = prin_rs::output::provenance_sidecar(
                &path,
                &ens,
                &format!(
                    "chart=latent_tilt_plambda panel={suffix} {what} \
                     window=({lo:.4e},{hi:.4e}) window_from=uniform_{finest} res={res} \
                     one_sample_per_pixel=true veto=none flagged_footprints={vetoed} of={}\n",
                    px.len()
                ),
            );
        }
        logln!(log, "{res:>5}^2 panels written; {vetoed}/{} ({:.4}%) flagged undetermined, drawn normally",
               px.len(), 100.0 * vetoed as f64 / px.len().max(1) as f64);
    }

    // ---- uniform over time -------------------------------------------------------------------
    //
    // One integration of a full grid with the live series kept, then one frame per recorded
    // boundary through `scheduler::project_at`. No tree anywhere: this is the control the
    // refinement animation is read against, same physics and same playhead.
    logln!(log, "\nuniform-over-time animation at {anim_res}^2, live_stride {stride}:");
    let ens_live = EnsembleCfg { keep_live_series: true, live_stride: stride, ..ens.clone() };
    let sl = grid::Slice::body_plane(anim_res, anim_res, cx, cy, half, 0).with_chart(chart);
    let t0 = std::time::Instant::now();

    // **Streamed strip by strip, and projected as it goes.** The first cut held
    // `Vec<PixelOut>` for the whole grid with the live series attached — about 2.3 KB a
    // footprint against `PixelSlim`'s 48 bytes, which is 2.4 GB at 1024² and is why the
    // animation was capped at 256. Each strip is projected into `frames_px[j]` and dropped, so
    // the resident set is the projected store alone.
    //
    // **`PixelSlim` and not a narrower record.** The panel is `Scalar::ShapeSpread` under
    // `Veto::None`, which reads only `shape_vec` and `spread_shape` — but the flagged count is
    // taken at `Veto::Debug`, which reads `n_nonfinite` and `state`. A store holding only the
    // two colouring fields would have made that count silently zero.
    let npix = sl.npix();
    let mut nb = 0usize;
    let mut store: Vec<PixelSlim> = Vec::new();
    let mut row = 0usize;
    while row < sl.ny {
        let hi = (row + 16).min(sl.ny);
        let band: Vec<PixelOut> = (row * sl.nx..hi * sl.nx)
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &ens_live))
            .collect();
        let here = band.iter().map(|p| p.live_t.len()).min().unwrap_or(0);
        if store.is_empty() {
            nb = here;
            assert!(nb > 1, "only {nb} live boundaries: keep_live_series produced nothing");
            store = vec![PixelSlim::thin(&band[0]); nb * npix];
        }
        // A strip may record fewer boundaries than the first one did; the original took the min
        // over the whole grid, so the frame count falls to the min and the store keeps its
        // allocation. It never rises above the first strip's count, which is the same rule.
        nb = nb.min(here);
        for (i, p) in band.iter().enumerate() {
            let k = row * sl.nx + i;
            for j in 0..nb {
                store[j * npix + k] = PixelSlim::thin(&scheduler::project_at(p, j));
            }
        }
        row = hi;
    }
    logln!(log, "  {npix} footprints in {:.1}s, {nb} boundaries recorded",
           t0.elapsed().as_secs_f64());

    // **One window over EVERY frame, not per frame.** A per-frame range would rescale the ramp as
    // the field develops and turn a growing signal into a constant one -- the auto-ranged-ramp
    // failure, on the axis the animation exists to show.
    let (alo, ahi) =
        colour::range_q_of(store[..nb * npix].iter().map(|p| p.spread_shape), 0.01, 0.99);
    logln!(log, "  fixed ramp window over all {nb} frames: ({alo:.4e}, {ahi:.4e})");

    let mut frames: Vec<Vec<u8>> = Vec::with_capacity(nb);
    let mut vetoed_any = 0usize;
    for j in 0..nb {
        let mut buf = Vec::with_capacity(npix * 3);
        for p in &store[j * npix..(j + 1) * npix] {
            let q = p.fat();
            if colour::vetoed(&q, Scalar::ShapeSpread, &sites, alo, ahi) {
                vetoed_any += 1;
            }
            buf.extend_from_slice(&colour::rgb_veto(
                &q, Scalar::ShapeSpread, &sites, alo, ahi, colour::Veto::None,
            ));
        }
        frames.push(buf);
    }
    // A duplicate-frame count is the arm that says the playhead moved. Nine identical frames and
    // a frozen playhead are the same picture, and only this number separates them.
    // **The control for the streaming rewrite, run wherever it is affordable.** At `anim_res <=
    // 64` the whole grid is evaluated the direct way — one `Vec<PixelOut>`, projected per frame,
    // exactly what this harness did before — and the frames are compared **bitwise**. A
    // restructure justified by memory has to show it changed nothing else, and the assertion
    // fires on every smoke pass rather than being a sentence in a commit message.
    if anim_res <= 64 {
        let direct: Vec<PixelOut> = (0..npix)
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, &ens_live))
            .collect();
        let dnb = direct.iter().map(|p| p.live_t.len()).min().unwrap_or(0);
        assert_eq!(dnb, nb, "the direct path records {dnb} boundaries against the streamed {nb}");
        for (j, f) in frames.iter().enumerate() {
            let mut want = Vec::with_capacity(npix * 3);
            for p in &direct {
                want.extend_from_slice(&colour::rgb_veto(
                    &scheduler::project_at(p, j), Scalar::ShapeSpread, &sites, alo, ahi,
                    colour::Veto::None,
                ));
            }
            assert_eq!(*f, want, "streamed frame {j} differs from the direct path");
        }
        logln!(log, "  control: {nb} frames bitwise identical to the direct whole-grid path");
    }

    let dup = apng::adjacent_duplicates(&frames);
    logln!(log, "  {nb} frames, {dup} adjacent duplicates, {vetoed_any} flagged footprint-frames");
    assert!(dup + 1 < nb, "every frame is a duplicate of its neighbour: the playhead never moved");
    let path = format!("{dir}/tilt_plambda_uniform_time.png");
    let _ = apng::write(&path, anim_res, anim_res, &frames, 1, 3);
    let _ = prin_rs::output::provenance_sidecar(
        &path,
        &ens_live,
        &format!(
            "chart=latent_tilt_plambda panel=uniform_over_time scalar=ShapeSpread \
             window=({alo:.4e},{ahi:.4e}) window_from=all_frames res={anim_res} frames={nb} \
             adjacent_duplicates={dup} live_stride={stride} tree=none veto=none\n"
        ),
    );
    logln!(log, "  wrote {path}");
}
