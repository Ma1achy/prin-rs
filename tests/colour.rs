//! Tests for the production colouring.
//!
//! Every one was chosen by asking what would have to be true for it to fire. Three of them fire
//! on faults that were actually shipped: [`the_hue_map_separates_shapes_the_old_projection_
//! merged`] fires on the `atan2` map's blindness to `n0`,
//! [`a_log_curve_spreads_a_decade_spanning_field_and_linear_does_not`] fires on the linear ramp
//! that produced a flat navy image, and
//! [`every_undetermined_case_is_debug_nan_and_never_the_background`] fires on the `NaN -> u8`
//! cast that made an undetermined pixel bitwise identical to un-rendered background.
//!
//! One test here was written to fire on the `atan2` map's *seam* and did not, because that map
//! has no seam: it reduces algebraically to `(a,b) = C_MAX*(n1,n2)`, a linear projection. The
//! test was replaced rather than loosened. See [`the_hue_map_is_continuous_on_the_sphere`].

use prin_rs::ensemble::pixel::PixelOut;
use prin_rs::outcome::State;
use prin_rs::output::colour::*;
use prin_rs::output::oklab;
use prin_rs::physics::burrau;
use prin_rs::physics::shape::shape_vec;
use prin_rs::Vec2;

fn sites() -> SiteSet {
    landmarks(&burrau::MASSES)
}

/// A footprint that is entirely determined, so a test that expects a *valid* colour is testing
/// the map rather than the null path.
fn healthy(n: [f64; 3], v: f64) -> PixelOut {
    PixelOut {
        shape_vec: n,
        ensemble_spread: v,
        spread_shape: v,
        state: State::Collision as u8,
        detail: 0,
        n_nonfinite: 0,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------------------------
// The vMF blend
// ---------------------------------------------------------------------------------------------

#[test]
fn vmf_weights_sum_to_one_and_are_invariant_to_a_constant_shift() {
    let d = [0.9, -0.2, 0.35, -0.87, 0.01];
    for &kappa in &[0.5, 3.0, 12.0] {
        let w = vmf_weights(&d, kappa).unwrap();
        let s: f64 = w.iter().sum();
        assert!((s - 1.0).abs() < 1e-15, "weights sum to {s}, not 1");

        // Subtracting the maximum is conditioning, not a parameter. If this fails, every colour
        // in every image moves by a small amount and nothing looks wrong.
        for &c in &[-40.0, -1.0, 0.25, 40.0] {
            let shifted: Vec<f64> = d.iter().map(|x| x + c).collect();
            let w2 = vmf_weights(&shifted, kappa).unwrap();
            for (a, b) in w.iter().zip(w2.iter()) {
                // Exact in exact arithmetic; in f64 the shifted subtraction loses low bits, so
                // the claim is a relative one. Measured worst drift at kappa=12, c=-40: 7.5e-15
                // relative. Stated as the tolerance rather than hidden behind an absolute one.
                let rel = (a - b).abs() / a.max(1e-300);
                assert!(rel < 1e-13, "shift by {c} moved a weight by {rel:e} relative: {a} vs {b}");
            }
        }
    }
}

#[test]
fn vmf_weights_decline_a_non_finite_input_rather_than_returning_a_number() {
    assert!(vmf_weights(&[0.1, f64::NAN], 3.0).is_none());
    assert!(vmf_weights(&[], 3.0).is_none());
}

#[test]
fn kappa_interpolates_between_a_smooth_blend_and_a_hard_voronoi() {
    // At high kappa the nearest site dominates; at low kappa the blend is near-uniform. If this
    // were flat in kappa the concentration would be an inert knob dressed as a design choice.
    let d = [1.0, 0.6, 0.2];
    let lo = vmf_weights(&d, 0.5).unwrap();
    let hi = vmf_weights(&d, 12.0).unwrap();
    assert!(hi[0] > lo[0] + 0.3, "kappa did not concentrate: {:?} -> {:?}", lo, hi);
    assert!(lo[2] > hi[2], "kappa did not suppress the far site");
}

#[test]
fn the_hue_map_is_continuous_on_the_sphere() {
    // Sweep a great circle and measure the largest step the OKLab (a,b) path takes between
    // adjacent samples. A continuous map's largest step falls with the sample count; a map with
    // a branch cut has one step that does not.
    //
    // NOTE this was originally written with the shipped `atan2` map as the negative control, on
    // the belief that it seams. It does not: `chroma*(cos h, sin h)` with `h = atan2(n2,n1)` and
    // `chroma = C_MAX*hypot(n1,n2)` is identically `C_MAX*(n1,n2)` (agreement 4.2e-17 over a
    // sphere sweep), which is linear and therefore continuous everywhere. The control is now a
    // map that genuinely wraps, so this test has teeth against the failure it names.
    let set = sites();
    let n_of = |t: f64| {
        let (z, r) = (0.35, (1.0f64 - 0.35 * 0.35).sqrt());
        [z, r * t.cos(), r * t.sin()]
    };

    let step = |samples: usize, f: &dyn Fn([f64; 3]) -> (f64, f64)| -> f64 {
        let mut worst: f64 = 0.0;
        let mut prev = f(n_of(0.0));
        for i in 1..=samples {
            let t = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
            let c = f(n_of(t));
            worst = worst.max(((c.0 - prev.0).powi(2) + (c.1 - prev.1).powi(2)).sqrt());
            prev = c;
        }
        worst
    };

    let blend = |n: [f64; 3]| hue_ab(&set, n).unwrap();
    let coarse = step(360, &blend);
    let fine = step(3600, &blend);
    assert!(fine < 0.002, "site-blend step at 3600 samples is {fine}");
    assert!(fine < coarse * 0.5, "step did not fall with sampling: {coarse} -> {fine}");

    // The negative control: an angle taken through a colour wheel, which is the construction the
    // colour reference warns about. Its step across the cut is O(1) and does not fall.
    let wheel = |n: [f64; 3]| {
        let h = n[2].atan2(n[1]);
        let t = (h + std::f64::consts::PI) / (2.0 * std::f64::consts::PI); // in [0,1], wraps
        let lab = oklab::srgb_to_oklab([(255.0 * t) as u8, 40, (255.0 * (1.0 - t)) as u8]);
        (lab[1], lab[2])
    };
    let (wc, wf) = (step(360, &wheel), step(3600, &wheel));
    assert!(
        wf > wc * 0.5,
        "the wheel control converged, so this test has no teeth: {wc} -> {wf}"
    );
}

#[test]
fn the_hue_map_separates_shapes_the_old_projection_merged() {
    // The real fault in the shipped map, as a test. `(a,b) = C_MAX*(n1,n2)` discards n0, so it
    // is exactly 2-to-1: `n` and its `n0 -> -n0` partner get bitwise identical colours. n0 is
    // (|rho~|^2 - |lam~|^2)/I, so those two are a tight binary with a distant third body and a
    // wide pair with a close third -- the hierarchical/anti-hierarchical distinction, painted
    // the same colour.
    let set = sites();
    let old = |n: [f64; 3]| {
        let h = n[2].atan2(n[1]);
        let c = C_MAX * (n[1] * n[1] + n[2] * n[2]).sqrt();
        (c * h.cos(), c * h.sin())
    };

    let mut worst_new = f64::INFINITY;
    for &n0 in &[0.2f64, 0.45, 0.7, 0.9] {
        let r = (1.0 - n0 * n0).sqrt();
        let p = [n0, r * 0.6, r * 0.8];
        let q = [-n0, r * 0.6, r * 0.8];

        // The control: the old map merges them exactly.
        let (op, oq) = (old(p), old(q));
        assert_eq!(op, oq, "the old map did not merge n0=+-{n0}, so this test proves nothing");

        let (a, b) = hue_ab(&set, p).unwrap();
        let (c, d) = hue_ab(&set, q).unwrap();
        worst_new = worst_new.min(((a - c).powi(2) + (b - d).powi(2)).sqrt());
    }
    assert!(
        worst_new > 0.01,
        "the site blend also merges the n0 pairs (worst separation {worst_new}) -- it reads all \
         three components only if the sites are not all in one plane"
    );
}

#[test]
fn chroma_shrinks_between_two_sites() {
    // The desaturation feature, asserted rather than eyeballed: a direction sitting between two
    // sites reads greyer, so "uncertain which regime this is" is visible as greyness.
    let set = sites();
    let a = set.sites[0].n;
    let b = set.sites[1].n;
    let mid = {
        let m = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        let l = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
        [m[0] / l, m[1] / l, m[2] / l]
    };
    let c = |n: [f64; 3]| {
        let (x, y) = hue_ab(&set, n).unwrap();
        (x * x + y * y).sqrt()
    };
    let (ca, cb, cm) = (c(a), c(b), c(mid));
    assert!(cm < ca && cm < cb, "chroma did not shrink at the midpoint: {ca}, {cm}, {cb}");
}

// ---------------------------------------------------------------------------------------------
// Landmarks
// ---------------------------------------------------------------------------------------------

#[test]
fn landmark_sites_are_unit_vectors_and_mutually_distinct() {
    let set = sites();
    assert_eq!(set.sites.len(), 6);
    for s in &set.sites {
        let l = (s.n[0] * s.n[0] + s.n[1] * s.n[1] + s.n[2] * s.n[2]).sqrt();
        assert!((l - 1.0).abs() < 1e-12, "site {} is not a unit vector: {l}", s.name);
    }
    for i in 0..set.sites.len() {
        for j in (i + 1)..set.sites.len() {
            let (a, b) = (set.sites[i].n, set.sites[j].n);
            let d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            assert!(
                d < 0.999,
                "sites {} and {} coincide (dot {d}) -- the palette has a hole",
                set.sites[i].name,
                set.sites[j].name
            );
        }
    }
}

#[test]
fn the_collision_sites_lie_on_the_collinear_circle_and_lagrange_does_not() {
    // A binary collision is a degenerate collinear configuration, so `q = 0` and `n[2] = 0`.
    // An equilateral triangle is as far from collinear as a configuration gets. If both landed
    // on the same circle the site set would be describing one feature, not two.
    let set = sites();
    for s in set.sites.iter().take(4) {
        assert!(s.n[2].abs() < 1e-12, "{} is not collinear: n2 = {}", s.name, s.n[2]);
    }
    for s in set.sites.iter().skip(4) {
        assert!(s.n[2].abs() > 0.5, "{} is nearly collinear: n2 = {}", s.name, s.n[2]);
    }
}

#[test]
fn the_landmarks_move_with_the_masses() {
    // The reason they are computed rather than hard-coded. If this were flat, a hard-coded set
    // would have been fine and the mass-varying charts would have shipped against a palette
    // that silently did not describe them.
    let a = landmarks(&burrau::MASSES);
    let b = landmarks(&[1.0, 1.0, 10.0]);
    let mut moved = 0;
    for (x, y) in a.sites.iter().zip(b.sites.iter()) {
        let d = x.n[0] * y.n[0] + x.n[1] * y.n[1] + x.n[2] * y.n[2];
        if d < 0.999 {
            moved += 1;
        }
    }
    assert!(moved >= 3, "only {moved} of 6 landmarks moved under a mass change");
}

#[test]
fn euler_points_are_central_configurations() {
    // The bisection is a measurement, so it gets checked against the defining property rather
    // than trusted. A central configuration has a_i = -lambda (r_i - R_com) with ONE lambda.
    // This is the check that a transcribed Euler quintic would have needed and would not have
    // got, which is why the quintic is not transcribed.
    use prin_rs::physics::newton::accel;

    let m = burrau::MASSES;
    let pts = euler_points(&m);
    for (k, n) in pts.iter().enumerate() {
        assert!(n.iter().all(|x| x.is_finite()), "euler point {k} did not bracket");
        assert!(n[2].abs() < 1e-9, "euler point {k} is not collinear: n2 = {}", n[2]);
    }

    // Re-derive the configuration for body 1 in the middle and check lambda directly.
    let (i, j, k) = (0usize, 2usize, 1usize);
    let build = |x: f64| {
        let mut r = [Vec2::zero(); 3];
        r[i] = Vec2::new(-1.0, 0.0);
        r[j] = Vec2::new(1.0, 0.0);
        r[k] = Vec2::new(x, 0.0);
        r
    };
    // Recover x by matching the shape vector the solver returned.
    let mut best = (f64::INFINITY, 0.0);
    for s in 0..200_001 {
        let x = -1.0 + 2.0 * s as f64 / 200_000.0;
        if x.abs() >= 1.0 - 1e-9 {
            continue;
        }
        let n = shape_vec(&build(x), &m);
        let d = (0..3).map(|c| (n[c] - pts[k][c]).powi(2)).sum::<f64>();
        if d < best.0 {
            best = (d, x);
        }
    }
    let r = build(best.1);
    let a = accel(&r, &m, 0.0);
    let mtot: f64 = m.iter().sum();
    let c = (0..3).map(|q| r[q].x * m[q]).sum::<f64>() / mtot;
    let lam: Vec<f64> = (0..3).map(|q| -a[q].x / (r[q].x - c)).collect();
    let spread = lam.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - lam.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(
        spread / lam[0].abs() < 1e-3,
        "lambda is not shared across bodies: {lam:?} (spread {spread})"
    );
    assert!(lam[0] > 0.0, "lambda must be positive for a bound central configuration");
}

// ---------------------------------------------------------------------------------------------
// Lightness
// ---------------------------------------------------------------------------------------------

#[test]
fn polarity_is_monotone_in_the_declared_direction_for_every_field() {
    // Driven off `Scalar::direction()`, so the test cannot pass by agreeing with a hard-coded
    // expectation that was itself copied from the implementation.
    let fields = [
        Scalar::Spread,
        Scalar::ShapeSpread,
        Scalar::EventSpread,
        Scalar::Ftle,
        Scalar::Diffusion,
        Scalar::ErrorRatio,
        Scalar::TEnd,
        Scalar::DMin,
    ];
    for s in fields {
        let (lo, hi) = match s.curve() {
            Curve::Log => (1e-6, 1.0),
            Curve::SymLog => (-1.0, 1.0),
            _ => (0.0, 1.0),
        };
        let mut prev = range_norm(s, lo, lo, hi).unwrap();
        let mut saw_change = false;
        for i in 1..=64 {
            let v = lo + (hi - lo) * i as f64 / 64.0;
            let t = range_norm(s, v, lo, hi).unwrap();
            match s.direction() {
                Direction::HighIsUnstable => assert!(
                    t >= prev - 1e-12,
                    "{} declared HighIsUnstable but L fell at v={v}",
                    s.name()
                ),
                Direction::HighIsSettled => assert!(
                    t <= prev + 1e-12,
                    "{} declared HighIsSettled but L rose at v={v}",
                    s.name()
                ),
            }
            if (t - prev).abs() > 1e-9 {
                saw_change = true;
            }
            prev = t;
        }
        // A ramp that never moves is monotone in both senses and would pass either arm above.
        assert!(saw_change, "{} produced a constant ramp -- the test could not fail", s.name());
    }
}

#[test]
fn a_log_curve_spreads_a_decade_spanning_field_and_linear_does_not() {
    // This is the shipped fault, as a test. `ensemble_spread` in near-field ran
    // (4.29e-5, 0.286) with a median near 1e-3; under a linear ramp the whole image sat at
    // L_MIN. The synthetic field here is log-uniform over the same span.
    let (lo, hi) = (4.29e-5f64, 0.286f64);
    let vals: Vec<f64> = (0..1000)
        .map(|i| lo * (hi / lo).powf(i as f64 / 999.0))
        .collect();

    assert_eq!(Scalar::Spread.curve(), Curve::Log, "the spread field must be log-ramped");

    let log_t: Vec<f64> =
        vals.iter().map(|&v| range_norm(Scalar::Spread, v, lo, hi).unwrap()).collect();
    let lin_t: Vec<f64> = vals.iter().map(|&v| ((v - lo) / (hi - lo)).clamp(0.0, 1.0)).collect();

    let bottom = |t: &[f64]| t.iter().filter(|&&x| x < 0.1).count() as f64 / t.len() as f64;
    let lin_bottom = bottom(&lin_t);
    let log_bottom = bottom(&log_t);
    assert!(
        lin_bottom > 0.6,
        "the linear control did not collapse ({lin_bottom:.3} in the bottom decile), so this \
         test does not demonstrate the fault it was written for"
    );
    assert!(
        log_bottom < 0.15,
        "the log ramp still piles up at the bottom: {log_bottom:.3}"
    );
}

#[test]
fn lightness_never_reaches_pure_black_or_pure_white() {
    // Hue lives in (a,b) but is invisible at either extreme of L, so an image that saturated
    // would be univariate while looking bivariate.
    for i in 0..=100 {
        let l = lightness(i as f64 / 100.0);
        assert!(l >= L_MIN - 1e-12 && l <= L_MAX + 1e-12);
        assert!(l > 0.0 && l < 1.0);
    }
    assert!(lightness(-5.0) >= L_MIN - 1e-12);
    assert!(lightness(5.0) <= L_MAX + 1e-12);
}

#[test]
fn range_norm_declines_rather_than_clamping_a_non_finite_value() {
    // The old map clamped a NaN scalar to t = 0, which rendered it at L_MIN -- the same colour
    // as the quietest genuine pixel in the region.
    assert!(range_norm(Scalar::Spread, f64::NAN, 1e-6, 1.0).is_none());
    assert!(range_norm(Scalar::Ftle, f64::INFINITY, 0.0, 1.0).is_none());
    // A degenerate window is a different failure from a degenerate value.
    assert!(range_norm(Scalar::Ftle, 0.5, 1.0, 1.0).is_none());
}

// ---------------------------------------------------------------------------------------------
// The reserved null
// ---------------------------------------------------------------------------------------------

#[test]
fn every_undetermined_case_is_debug_nan_and_never_the_background() {
    let set = sites();
    let good = shape_vec(&[Vec2::new(0.0, 0.0), Vec2::new(1.0, 0.0), Vec2::new(0.3, 0.9)], &burrau::MASSES);

    // The control: a healthy footprint must NOT be DEBUG_NAN, or every arm below passes for
    // the wrong reason.
    let ok = rgb(&healthy(good, 1e-3), Scalar::Spread, &set, 1e-6, 1.0);
    assert_ne!(ok, DEBUG_NAN, "a healthy footprint rendered as undetermined");
    assert_ne!(ok, BACKGROUND);

    let mut nan_shape = healthy(good, 1e-3);
    nan_shape.shape_vec = [f64::NAN; 3];
    assert_eq!(rgb(&nan_shape, Scalar::Spread, &set, 1e-6, 1.0), DEBUG_NAN);

    let mut nan_scalar = healthy(good, 1e-3);
    nan_scalar.ensemble_spread = f64::NAN;
    assert_eq!(rgb(&nan_scalar, Scalar::Spread, &set, 1e-6, 1.0), DEBUG_NAN);

    let mut nonfinite_copy = healthy(good, 1e-3);
    nonfinite_copy.n_nonfinite = 1;
    assert_eq!(rgb(&nonfinite_copy, Scalar::Spread, &set, 1e-6, 1.0), DEBUG_NAN);

    for st in [State::SimFailed, State::DecodeFailed] {
        let mut p = healthy(good, 1e-3);
        p.state = st as u8;
        assert_eq!(
            rgb(&p, Scalar::Spread, &set, 1e-6, 1.0),
            DEBUG_NAN,
            "{} did not render as undetermined",
            st.name()
        );
    }

    // An invalid state byte is also undetermined, not a dark grey pixel.
    let mut bad = healthy(good, 1e-3);
    bad.state = 7;
    assert_eq!(rgb(&bad, Scalar::Spread, &set, 1e-6, 1.0), DEBUG_NAN);

    assert_ne!(DEBUG_NAN, BACKGROUND, "the null and the background must be distinguishable");
    assert_ne!(DEBUG_NAN, [0, 0, 0]);
    assert_ne!(BACKGROUND, [0, 0, 0], "a NaN cast lands on black, so background must not be");
}

#[test]
fn the_outcome_palette_flags_a_failed_decode() {
    // `State::DecodeFailed` previously fell to the catch-all grey, indistinguishable from an
    // invalid state byte -- so a pixel whose IC could not be formed read as ordinary data.
    let mut p = PixelOut { state: State::DecodeFailed as u8, ..Default::default() };
    assert_eq!(prin_rs::output::png::outcome_rgb(&p), DEBUG_NAN);
    p.state = State::SimFailed as u8;
    assert_eq!(prin_rs::output::png::outcome_rgb(&p), DEBUG_NAN);
    // And the control: an ordinary state is not flagged.
    p.state = State::Collision as u8;
    assert_ne!(prin_rs::output::png::outcome_rgb(&p), DEBUG_NAN);
}

// ---------------------------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------------------------

#[test]
fn quantisation_counts_exact_distinct_values() {
    // The diagnostic that says whether a lightness field carries an ordering at all. A count
    // ratio over E+1 copies has E+2 levels and no ramp recovers what is not there.
    let px: Vec<PixelOut> = (0..800)
        .map(|i| {
            let mut p = PixelOut::default();
            p.spread_event = (i % 8) as f64 / 7.0;
            p.spread_shape = 0.0;
            p.ensemble_spread = p.spread_event;
            p
        })
        .collect();
    let (distinct, finite, modal) = quantisation(&px, Scalar::Spread);
    assert_eq!(distinct, 8, "a count ratio over 8 copies has 8 levels");
    assert_eq!(finite, 800);
    assert!((modal - 0.125).abs() < 1e-9);

    assert!((event_arm_fraction(&px) - 1.0).abs() < 1e-12);

    // The contrast: a continuous field resolves finely.
    let cont: Vec<PixelOut> = (0..800)
        .map(|i| {
            let mut p = PixelOut::default();
            p.ensemble_spread = i as f64 * 1e-4;
            p
        })
        .collect();
    assert_eq!(quantisation(&cont, Scalar::Spread).0, 800);
}

#[test]
fn the_site_set_leaves_no_large_hole_on_the_sphere() {
    // A hole means a whole neighbourhood of shapes shares one flat blended colour.
    let gap = worst_site_gap(&sites(), 4000);
    assert!(gap < 1.15, "worst angular gap to the nearest site is {gap} rad");
}

#[test]
fn oklab_round_trips_and_matches_ottosson() {
    // Kept from the previous build: the transcription check for the constants everything above
    // rests on.
    for c in [[0u8, 0, 0], [255, 255, 255], [255, 0, 0], [18, 200, 77], [130, 40, 210]] {
        let back = oklab::oklab_to_srgb(oklab::srgb_to_oklab(c));
        for k in 0..3 {
            assert!(
                (back[k] as i32 - c[k] as i32).abs() <= 1,
                "round trip moved {c:?} to {back:?}"
            );
        }
    }
    let w = oklab::linear_to_oklab(1.0, 1.0, 1.0);
    assert!((w[0] - 1.0).abs() < 1e-4, "white is not L=1: {w:?}");
    assert!(w[1].abs() < 1e-4 && w[2].abs() < 1e-4, "white is not neutral: {w:?}");
}
