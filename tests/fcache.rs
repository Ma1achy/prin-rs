//! Tests for the per-footprint cache.
//!
//! The claim this file exists to support is exact: **recolouring from a footprint file gives
//! bitwise the same cache as integrating again under that colouring.** If it were merely close,
//! the replay would produce `error(B)` curves that looked like measurements and were artefacts
//! of the replay path, and nothing downstream could tell the difference.

use std::collections::HashMap;

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::Chart;
use prin_rs::metric::{self, Colouring};
use prin_rs::output::colour::Scalar;
use prin_rs::output::fcache::{self, Row};

fn small() -> (Vec<metric::Cache>, HashMap<metric::Key, Vec<prin_rs::ensemble::pixel::PixelOut>>) {
    // Level 1, N = 4: 5 quads, 80 footprints. Small enough to be a test and large enough that a
    // ramp, a reference image and an err_sum all have something to do.
    let ens = EnsembleCfg { refine_flagged: false, t_max: 2.0, n_sync: 5, ..Default::default() };
    metric::build_multi_with_footprints(
        "near-field",
        1.0,
        3.0,
        0.05,
        0,
        Chart::BodyPlane,
        1,
        4,
        8,
        1e-4,
        &ens,
        &[Colouring::Outcome, Colouring::Bivariate(Scalar::Spread)],
    )
}

#[test]
fn a_recolour_is_bitwise_identical_to_a_fresh_build() {
    let (caches, px_of) = small();
    let outcome = &caches[0];
    let fresh_bivariate = &caches[1];

    let fp = outcome.footprints_from(&px_of, "body_plane", 2.0);
    let replayed = outcome.recolour(&fp, Colouring::Bivariate(Scalar::Spread)).unwrap();

    // The control first: the two colourings must actually DIFFER, or "identical" below is
    // satisfied by everything being the same and the test proves nothing.
    assert_ne!(
        outcome.reference, fresh_bivariate.reference,
        "the two colourings render identically, so this test cannot fail"
    );

    assert_eq!(replayed.ramp, fresh_bivariate.ramp, "the region-wide ramp did not survive");
    assert_eq!(
        replayed.reference, fresh_bivariate.reference,
        "the reference image differs between a replay and a fresh build"
    );
    for (k, q) in &fresh_bivariate.quads {
        let r = replayed.quads.get(k).expect("replay lost a quad");
        assert_eq!(&r.rgb, &q.rgb, "quad {k:?} coloured differently on replay");
        assert_eq!(
            r.err_sum.to_bits(),
            q.err_sum.to_bits(),
            "quad {k:?} err_sum differs: replay {} vs fresh {}",
            r.err_sum,
            q.err_sum
        );
    }
}

#[test]
fn the_reductions_carry_over_untouched_by_a_recolour() {
    // A QuadReduction is a property of the physics and does not know what colour anything is
    // drawn. If a recolour moved one, the criterion would be scoring a different tree.
    let (caches, px_of) = small();
    let c = &caches[0];
    let fp = c.footprints_from(&px_of, "body_plane", 2.0);
    let r = c.recolour(&fp, Colouring::Bivariate(Scalar::Ftle)).unwrap();
    for (k, q) in &c.quads {
        let a = &r.quads[k].red;
        let b = &q.red;
        assert_eq!(a.spread_median.to_bits(), b.spread_median.to_bits());
        assert_eq!(a.between_shape.to_bits(), b.between_shape.to_bits());
        assert_eq!(a.n_distinct_ic, b.n_distinct_ic);
    }
}

#[test]
fn a_footprint_file_round_trips_through_disk() {
    let (caches, px_of) = small();
    let fp = caches[0].footprints_from(&px_of, "body_plane", 2.0);

    let mut buf: Vec<u8> = Vec::new();
    fcache::write(&mut buf, &fp).unwrap();
    let back = fcache::read(&mut std::io::Cursor::new(&buf)).unwrap();

    assert_eq!(back.region, fp.region);
    assert_eq!(back.chart, fp.chart);
    assert_eq!(back.levels, fp.levels);
    assert_eq!(back.n, fp.n);
    assert_eq!(back.res, fp.res);
    assert_eq!(back.cx.to_bits(), fp.cx.to_bits());
    assert_eq!(back.quads.len(), fp.quads.len());
    for (k, rows) in &fp.quads {
        let got = back.quads.get(k).expect("a quad did not survive the round trip");
        assert_eq!(got.len(), rows.len());
        for (a, b) in got.iter().zip(rows.iter()) {
            // Bitwise, including NaN: a NaN ftle is a measurement outcome and must come back a
            // NaN rather than a zero, which would read as "no divergence".
            assert_eq!(a.shape[0].to_bits(), b.shape[0].to_bits());
            assert_eq!(a.packed, b.packed);
            assert_eq!(a.n_nonfinite, b.n_nonfinite);
            assert_eq!(a.ensemble_spread.to_bits(), b.ensemble_spread.to_bits());
            assert_eq!(a.ftle.to_bits(), b.ftle.to_bits());
            assert_eq!(a.diffusion.to_bits(), b.diffusion.to_bits());
        }
    }

    // And a recolour from the DISK copy matches one from the in-memory copy.
    let from_disk = caches[0].recolour(&back, Colouring::Bivariate(Scalar::Spread)).unwrap();
    let from_mem = caches[0].recolour(&fp, Colouring::Bivariate(Scalar::Spread)).unwrap();
    assert_eq!(from_disk.reference, from_mem.reference);
}

#[test]
fn a_nan_survives_the_round_trip_as_a_nan() {
    // `ftle` is NaN when n_renorm == 0 and `shape_vec` is NaN at a triple collision. Both are
    // measurement outcomes. A format that wrote them as zero would turn "could not be
    // determined" into "determined to be the quietest value in the region".
    let r = Row {
        shape: [f64::NAN, 0.5, -0.5],
        ftle: f64::NAN,
        diffusion: f64::NEG_INFINITY,
        packed: (5 << 2) | 1,
        n_nonfinite: 3,
        ..Default::default()
    };
    let fp = fcache::Footprints {
        region: "t".into(),
        chart: "body_plane".into(),
        cx: 0.0,
        cy: 0.0,
        half: 1.0,
        body: 0,
        levels: 0,
        n: 1,
        res: 1,
        t_max: 1.0,
        quads: [((0u32, 0u32, 0u32), vec![r])].into_iter().collect(),
    };
    let mut buf: Vec<u8> = Vec::new();
    fcache::write(&mut buf, &fp).unwrap();
    let back = fcache::read(&mut std::io::Cursor::new(&buf)).unwrap();
    let got = back.quads[&(0, 0, 0)][0];
    assert!(got.shape[0].is_nan(), "a NaN shape came back as {}", got.shape[0]);
    assert!(got.ftle.is_nan());
    assert_eq!(got.diffusion, f64::NEG_INFINITY);
    assert_eq!(got.packed, r.packed);
    assert_eq!(got.n_nonfinite, 3);

    // And the projection preserves the state/detail split the palette reads.
    let p = got.to_pixel();
    assert_eq!(p.state, 5, "DecodeFailed must survive as DecodeFailed");
    assert_eq!(p.detail, 1);
    assert_eq!(p.n_nonfinite, 3);
}

#[test]
fn a_mismatched_footprint_file_is_refused_rather_than_recoloured() {
    // A file from a different region would recolour without complaint and produce an error(B)
    // curve for a tree that was never integrated. A self-describing header only prevents that
    // if something reads it.
    let (caches, px_of) = small();
    let mut fp = caches[0].footprints_from(&px_of, "body_plane", 2.0);

    fp.region = "far".into();
    assert!(caches[0].recolour(&fp, Colouring::Outcome).is_err(), "a wrong region was accepted");

    fp.region = "near-field".into();
    fp.res = 16;
    assert!(caches[0].recolour(&fp, Colouring::Outcome).is_err(), "a wrong resolution was accepted");

    fp.res = 8;
    fp.cx = 1.5;
    assert!(caches[0].recolour(&fp, Colouring::Outcome).is_err(), "a wrong centre was accepted");

    // The control: put it back and it is accepted, so the arms above are not failing on
    // something unrelated.
    fp.cx = 1.0;
    assert!(caches[0].recolour(&fp, Colouring::Outcome).is_ok());
}
