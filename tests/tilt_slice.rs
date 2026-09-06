//! **`tilt_plambda` — the first slice in the project with a non-zero tilt.**
//!
//! Supplied as a reference-UI config. Three things about it are checkable without integrating
//! anything, and all three are the kind that fail silently:
//!
//! 1. the ten-slot `z0` carries **two dead dimensions** the canonical frame consumes;
//! 2. it decodes to **unequal masses**, where every `z0 = 0` preset decodes to `1/3` each;
//! 3. the tilt is **in-plane**, so this is a rotated frame and not a new 2-plane — the thing
//!    `shape_pl` is on record for getting wrong in the other direction.
//!
//! Each is asserted with the arm that says the test could have failed, because *a test that
//! cannot fail is indistinguishable from a test that passes*.

use prin_rs::decode::{self, Path};
use prin_rs::grid::{self, Chart};
use prin_rs::physics::{decoder, Cart};

/// The supplied ten-slot config, verbatim.
const Z10: [f64; 10] = [2.23, -0.56, -0.05, 0.0, 0.05, -0.04, -0.02, 0.12, -0.1, 0.02];

fn ic_at(chart: Chart, cx: f64, cy: f64, half: f64, du: f64, dv: f64) -> Cart<f64> {
    let lin = decode::linearise(&chart, 0, cx, cy, half);
    decode::sample(Path::DirectF64, &chart, 0, cx, cy, half, du, dv, &lin)
}

fn max_dic(a: &Cart<f64>, b: &Cart<f64>) -> f64 {
    (0..3)
        .map(|i| {
            (a.r[i].x - b.r[i].x)
                .abs()
                .max((a.r[i].y - b.r[i].y).abs())
                .max((a.v[i].x - b.v[i].x).abs())
                .max((a.v[i].y - b.v[i].y).abs())
        })
        .fold(0.0, f64::max)
}

/// **Slots 2 and 3 are dead: `decodeIC` never reads them.** The supplier verified this by
/// perturbation; it is committed here rather than remembered, because "two of these ten numbers
/// do nothing" is exactly the fact that rots.
///
/// The control is the other eight: perturbing any of them **must** move the IC, or the probe is
/// measuring a decoder that ignores its input and the zeros above mean nothing.
#[test]
fn two_of_the_ten_slots_are_dead_and_the_other_eight_are_not() {
    let base = Chart::latent_ui_slice(Z10, 4, 5, &[], 0.0, 0.5, (0.5, 0.5)).0;
    let a = ic_at(base, 0.0, 0.0, 0.5, 0.3, -0.2);

    for slot in 0..10 {
        let mut z = Z10;
        z[slot] += 0.7;
        let c = Chart::latent_ui_slice(z, 4, 5, &[], 0.0, 0.5, (0.5, 0.5)).0;
        let d = max_dic(&a, &ic_at(c, 0.0, 0.0, 0.5, 0.3, -0.2));
        if slot == 2 || slot == 3 {
            assert_eq!(d, 0.0, "slot {slot} is supposed to be dead but moved the IC by {d:e}");
        } else {
            assert!(d > 1e-6, "slot {slot} is live but moved the IC by only {d:e}");
        }
    }
}

/// **Unequal masses, and they are the supplier's own numbers.** A run reporting `1/3` each here
/// is not decoding this chart — that is the check asked for, and it needs the preset control
/// beside it or "unequal" is a claim about nothing.
#[test]
fn the_slice_decodes_to_unequal_masses_where_the_presets_do_not() {
    let (chart, ..) = Chart::tilt_plambda();
    let Chart::Latent { z0, .. } = chart else { panic!("tilt_plambda is not a Latent chart") };
    let (m, flag) = decoder::masses(z0.z_mu);
    assert_eq!(flag, None);
    for (got, want) in m.iter().zip([0.353_33, 0.275_23, 0.371_44]) {
        assert!((got - want).abs() < 5e-6, "masses {m:?} against the supplied {want}");
    }
    assert!((m.iter().sum::<f64>() - 1.0).abs() < 1e-15);

    // The control: every `z0 = 0` preset is equal-mass, so an equal-mass reading is reachable
    // and this test is discriminating rather than describing a decoder that cannot produce it.
    let Chart::Latent { z0: zp, .. } = Chart::preset_plambda() else { unreachable!() };
    let (mp, _) = decoder::masses(zp.z_mu);
    for x in mp {
        assert!((x - 1.0 / 3.0).abs() < 1e-15, "the preset control is not equal-mass: {mp:?}");
    }
}

/// **The basis is orthonormal, and the frame angle is exactly `amt + gamma`.**
///
/// The angle arm is what pins the ORDER of orthonormalisation and rotation. Rotating a
/// non-orthonormal pair by a matrix is not a rotation: applied before Gram-Schmidt this slice's
/// `gamma = 4.5 deg` turns the frame by 2.45 deg, so `gammaDeg` would not be degrees. The
/// negative control asserts that wrong order is measurably different rather than merely stating
/// it.
#[test]
fn the_basis_is_orthonormal_and_turned_by_exactly_amt_plus_gamma() {
    let (chart, ..) = Chart::tilt_plambda();
    let Chart::Latent { q1, q2, .. } = chart else { unreachable!() };
    let dot = |a: &[f64; 8], b: &[f64; 8]| (0..8).map(|k| a[k] * b[k]).sum::<f64>();
    assert!((dot(&q1, &q1) - 1.0).abs() < 1e-14);
    assert!((dot(&q2, &q2) - 1.0).abs() < 1e-14);
    assert!(dot(&q1, &q2).abs() < 1e-14);

    // The shipped slice runs `gamma = 2.0`; the supplied config said 4.5 and is carried below,
    // so a reduction of the shipped value cannot quietly retire the arithmetic.
    let want = -1.04 + 2.0f64.to_radians();
    let got = q1[5].atan2(q1[4]);
    assert!((got - want).abs() < 1e-12, "frame at {got} rad, expected amt+gamma = {want}");

    let supplied = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 5, -1.04)], 4.5, 0.5, (0.5, 0.5)).0;
    let Chart::Latent { q1: s1, .. } = supplied else { unreachable!() };
    let want45 = -1.04 + 4.5f64.to_radians();
    assert!((s1[5].atan2(s1[4]) - want45).abs() < 1e-12, "the supplied 4.5 arm no longer holds");
    assert!((want45 - want).abs() > 0.04, "the two gammas are not distinguishable by angle");

    // The wrong order, computed here rather than asserted about: mix first, normalise after.
    let (c, s) = ((-1.04f64).cos(), (-1.04f64).sin());
    let (g, h) = (4.5f64.to_radians().cos(), 4.5f64.to_radians().sin());
    let (x, y) = (g * c, g * s + h); // q1 = cos(g)*tilted + sin(g)*e5, un-normalised
    let wrong = y.atan2(x);
    assert!(
        (wrong - want45).abs() > 0.03,
        "the two orders agree to {:.4} rad, so the angle arm cannot tell them apart",
        (wrong - want45).abs()
    );
}

/// **The tilt is IN-PLANE: this is `preset_plambda`'s 2-plane with a rotated frame.**
///
/// `dim 5` is `q2` itself, so the rotation never leaves `span{e4, e5}`. Stated as a measurement
/// because the opposite mistake is on this project's record — `shape_pl`'s crossed basis looked
/// like a reorientation and was a genuinely different 2-plane, told apart by `max |dIC|`.
///
/// The negative control is a tilt into a genuinely hidden **live** dimension, which must leave
/// the span. Without it, "the leak is zero" would pass on a constructor that ignored `tilts`
/// entirely.
#[test]
fn the_tilt_rotates_the_frame_without_leaving_the_plane_and_a_hidden_dim_does_leave_it() {
    let leak = |c: Chart| {
        let Chart::Latent { q1, q2, .. } = c else { unreachable!() };
        (0..8)
            .filter(|k| *k != 4 && *k != 5)
            .map(|k| q1[k].abs().max(q2[k].abs()))
            .fold(0.0, f64::max)
    };
    assert_eq!(leak(Chart::tilt_plambda().0), 0.0, "the shipped tilt left span(e4,e5)");

    // Into `z_mu2` (live-8D dim 7), a hidden and live coordinate: the plane must move.
    let hidden = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 7, -1.04)], 4.5, 0.5, (0.5, 0.5)).0;
    assert!(leak(hidden) > 0.5, "a tilt into a hidden dim did not leave the plane: the \
                                 constructor is ignoring `tilts` and the zero above is vacuous");

    // And the in-plane tilt is not inert either — the FRAME turns, so a pixel maps elsewhere.
    // Two ways to be indistinguishable from `preset_plambda`, and this rules out the second.
    let (t, cx, cy, half) = Chart::tilt_plambda();
    let flat = Chart::latent_ui_slice(Z10, 4, 5, &[], 0.0, half, (0.5, 0.5)).0;
    let d = max_dic(&ic_at(t, cx, cy, half, 0.6, -0.4), &ic_at(flat, cx, cy, half, 0.6, -0.4));
    assert!(d > 1e-3, "the rotated frame maps the same pixel to the same IC ({d:e}); the tilt \
                       changes nothing observable");
}

/// **`pan` is the window CENTRE, and the arm that decides it is the gamma sweep.**
///
/// The first cut used `config_slice`'s `2*pan - 1 + zoom`, which places the lower-left *corner* at
/// `pan` — half a window out in each axis, which is large and does not look like an error. Two
/// independently supplied figures refuse it, and both are asserted here rather than the bare
/// arithmetic: the camera sits **0.0419 uv** from `z0`, and because gamma pivots about `z0` a
/// rotation of 4.5 deg sweeps the image by **20 px** at 1024. Under the corner reading those are
/// 0.1585 and 75.9. Asserting `cx == 2*pan - 1` alone would pass on any convention that happened
/// to be typed in twice; these two are properties of the picture.
#[test]
fn pan_is_the_window_centre_and_gamma_sweeps_the_camera_by_twenty_pixels() {
    let (_, cx, cy, half) = Chart::tilt_plambda();
    let (zoom, px, py) =
        (0.167_798_447_231_782_42, 0.516_965_993_956_604_7, 0.538_322_670_350_332_3);
    assert_eq!(half, zoom);
    assert_eq!(cx, 2.0 * px - 1.0);
    assert_eq!(cy, 2.0 * py - 1.0);

    // `z0` sits at signed (0,0), so the camera's uv offset from it is |pan - (0.5, 0.5)|.
    let off_uv = ((cx + 1.0) / 2.0 - 0.5).hypot((cy + 1.0) / 2.0 - 0.5);
    assert!((off_uv - 0.0419).abs() < 5e-5, "camera {off_uv:.4} uv off-centre, expected 0.0419");

    // Gamma pivots about `z0`, so an off-centre camera sweeps. 1024 px across a window of 2*half.
    // **The 20 px figure is the SUPPLIED 4.5 deg**, which is the number that decided the
    // convention, so it is pinned at the value it was supplied at and not at the shipped one.
    let sweep = |deg: f64| 2.0 * cx.hypot(cy) * (deg.to_radians() / 2.0).sin() / (2.0 * zoom) * 1024.0;
    let shift_px = sweep(4.5);
    assert!((shift_px - 20.0).abs() < 0.5, "gamma sweeps {shift_px:.1} px, expected 20");

    // The shipped slice runs 2.0 deg, which is the same pivot at a smaller angle: linear in
    // `sin(gamma/2)`, so a little under half. Asserted so a change of the shipped value has to
    // come past this line.
    let shipped = sweep(2.0);
    assert!((shipped - 8.9).abs() < 0.5, "the shipped 2 deg sweeps {shipped:.1} px, expected ~8.9");

    // The negative control: pivoting about the CAMERA instead would move nothing at all, and
    // the corner convention would sweep nearly four times as far. Both are excluded by number.
    let (bx, by) = (2.0 * px - 1.0 + zoom, 2.0 * py - 1.0 + zoom);
    let g = 4.5f64.to_radians();
    let corner_px = 2.0 * bx.hypot(by) * (g / 2.0).sin() / (2.0 * zoom) * 1024.0;
    assert!(corner_px > 70.0, "the corner convention sweeps {corner_px:.1} px; the two \
                               conventions are not distinguishable by this test");
}

/// **Gamma stays in the plane; a tilt into a hidden dimension does not.** The two are different
/// operations and the test says so in the quantity that separates them: gamma is a rotation about
/// the normal through `z0`, so it can never change which slice is seen — only which direction is
/// "right" on screen.
#[test]
fn gamma_never_leaves_the_plane_and_a_hidden_tilt_always_does() {
    let span_leak = |c: Chart, keep: [usize; 2]| {
        let Chart::Latent { q1, q2, .. } = c else { unreachable!() };
        (0..8)
            .filter(|k| !keep.contains(k))
            .map(|k| q1[k].abs().max(q2[k].abs()))
            .fold(0.0, f64::max)
    };
    // Gamma alone, swept hard: the plane is untouched at every angle.
    for deg in [0.0, 4.5, 37.0, 90.0, 180.0, -73.5] {
        let c = Chart::latent_ui_slice(Z10, 4, 5, &[], deg, 0.5, (0.5, 0.5)).0;
        assert_eq!(span_leak(c, [4, 5]), 0.0, "gamma {deg} deg left span(e4,e5)");
    }
    // A tilt into a hidden live dim leaves it at every non-zero amount.
    for amt in [0.3, -1.04, 1.2] {
        let c = Chart::latent_ui_slice(Z10, 4, 5, &[(0, 7, amt)], 0.0, 0.5, (0.5, 0.5)).0;
        assert!(span_leak(c, [4, 5]) > 0.2, "a tilt of {amt} rad did not leave the plane");
    }
}

/// Resolvable by name from the one table, so every harness that takes a chart string reaches it.
#[test]
fn the_slice_resolves_by_name() {
    let (c, cx, cy, half) = grid::named_slice("tilt_plambda").expect("not in named_slice");
    let (c2, cx2, cy2, half2) = Chart::tilt_plambda();
    assert_eq!(format!("{c:?}"), format!("{c2:?}"));
    assert_eq!((cx, cy, half), (cx2, cy2, half2));
    assert!(grid::named_slice("config_stability").is_some());
    assert!(grid::named_slice("no_such_slice").is_none());
}
