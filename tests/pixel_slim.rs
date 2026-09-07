//! **The thinned pixel colours identically, on real pixels, with both populations present.**
//!
//! [`PixelSlim`] exists so a full-resolution uniform grid need not be held as `PixelOut` -- 656
//! bytes against 40, which is 688 MB against 42 MB at `1024^2` and the difference between a panel
//! that can be regenerated on this machine and one that cannot. It is only sound if widening it
//! back reproduces every panel byte for byte.
//!
//! **The arm with teeth is the population check.** `fat()` writes `Default` into every field it
//! does not carry, so a fixture with no flagged footprint and no unusual terminal state would
//! agree trivially -- the round trip would be exercised on the easy half only. The test asserts
//! the fixture contains **both** a flagged and an unflagged footprint and more than one terminal
//! state before it reads the agreement, and says so in the failure message.

use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg, PixelSlim};
use prin_rs::grid;
use prin_rs::output::colour::{self, Scalar, Veto};
use prin_rs::output::png;

#[test]
fn widening_a_slim_pixel_reproduces_all_three_panels() {
    let (_, chart, cx, cy, half) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == "latent_mixed_h3")
        .expect("latent_mixed_h3 is in the gallery");
    let n = 12;
    let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
    // The old budget, deliberately: it truncates ~12% of this chart, which is what puts flagged
    // footprints in the fixture at all. Under the shipped budget none of them are flagged and the
    // population arm below would refuse the run.
    let ens = EnsembleCfg { max_steps: 30_000, ..EnsembleCfg::production() };
    let px: Vec<_> = (0..sl.npix()).map(|k| evaluate::<f64>(&sl, k, &ens)).collect();

    let flagged = px.iter().filter(|p| p.n_nonfinite > 0).count();
    let states: std::collections::BTreeSet<u8> = px.iter().map(|p| p.state).collect();
    println!(
        "latent_mixed_h3 {n}^2: {} pixels, {flagged} flagged, {} distinct terminal states",
        px.len(),
        states.len()
    );
    assert!(
        flagged > 0 && flagged < px.len(),
        "the fixture must hold both flagged and unflagged footprints or the round trip is \
         exercised on one population only (flagged {flagged} of {})",
        px.len()
    );
    assert!(states.len() > 1, "one terminal state in the fixture: the outcome panel is constant");

    let m = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m);
    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
    // A window that is not degenerate, or the ramp maps everything to one colour and agrees.
    assert!(hi > lo, "degenerate ramp window ({lo:.4e}, {hi:.4e})");

    for (k, p) in px.iter().enumerate() {
        let q = PixelSlim::thin(p).fat();
        for v in [Veto::None, Veto::Quiet, Veto::Debug] {
            assert_eq!(
                colour::rgb_veto(p, Scalar::ShapeSpread, &sites, lo, hi, v),
                colour::rgb_veto(&q, Scalar::ShapeSpread, &sites, lo, hi, v),
                "spread panel disagrees at pixel {k} under {v:?}"
            );
            assert_eq!(
                png::outcome_rgb_veto(p, v),
                png::outcome_rgb_veto(&q, v),
                "outcome panel disagrees at pixel {k} under {v:?}"
            );
            assert_eq!(
                png::event_class_rgb_veto(p, v),
                png::event_class_rgb_veto(&q, v),
                "event panel disagrees at pixel {k} under {v:?}"
            );
        }
        assert_eq!(
            colour::vetoed(p, Scalar::ShapeSpread, &sites, lo, hi),
            colour::vetoed(&q, Scalar::ShapeSpread, &sites, lo, hi),
            "veto set disagrees at pixel {k}"
        );
    }

    // The whole-slice statistics the panel's header carries.
    let slim: Vec<PixelSlim> = px.iter().map(PixelSlim::thin).collect();
    assert_eq!(
        colour::range(&px, Scalar::ShapeSpread),
        colour::range_q_of(slim.iter().map(|p| p.spread_shape), 0.01, 0.99)
    );
    assert_eq!(
        colour::quantisation(&px, Scalar::ShapeSpread),
        colour::quantisation_of(slim.iter().map(|p| p.spread_shape))
    );
    assert_eq!(
        png::event_class_histogram(&px),
        png::event_class_histogram_of(
            slim.iter().map(|p| (p.n_nonfinite, p.state, p.event_class))
        )
    );

    // **And the negative control: the thinning is not vacuous.** A `PixelOut` differing only in a
    // field `PixelSlim` drops must still be a different `PixelOut` -- if this passed, `fat()`
    // would be reproducing the input rather than a projection of it.
    let mut mutated = px[0].clone();
    mutated.energy_drift_max = 1.0;
    assert_ne!(mutated.energy_drift_max, PixelSlim::thin(&mutated).fat().energy_drift_max);
}

/// **The strip loop covers the grid exactly, at a strip height that does not divide it.**
///
/// The memory saving is worthless if the loop drops or repeats a row, and a strip height chosen
/// to divide the raster would never show it -- `1024 / 64` is exact, so the shipped configuration
/// is the one case that cannot fail. This runs `STRIP = 5` over a 12-row grid.
#[test]
fn the_strip_pass_reproduces_the_whole_grid_pass() {
    let (_, chart, cx, cy, half) = grid::gallery_cases()
        .into_iter()
        .find(|c| c.0 == "latent_shape")
        .expect("latent_shape is in the gallery");
    let n = 12;
    let sl = grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart);
    let ens = EnsembleCfg::production();
    let flat: Vec<PixelSlim> =
        (0..sl.npix()).map(|k| PixelSlim::thin(&evaluate::<f64>(&sl, k, &ens))).collect();

    const STRIP: usize = 5;
    assert_ne!(n % STRIP, 0, "a strip height that divides the grid cannot catch a bound error");
    let mut striped: Vec<PixelSlim> = Vec::with_capacity(sl.npix());
    let mut row = 0usize;
    while row < n {
        let hi_row = (row + STRIP).min(n);
        striped.extend(
            (row * n..hi_row * n).map(|k| PixelSlim::thin(&evaluate::<f64>(&sl, k, &ens))),
        );
        row = hi_row;
    }
    assert_eq!(striped.len(), flat.len(), "the strip loop covered {} of {} pixels", striped.len(), flat.len());
    assert_eq!(striped, flat, "the strip pass differs from the whole-grid pass");
}
