//! **The sampling-measure weight, against an independent estimate of the same area.**
//!
//! The scheduler contract's Part 2 is the reason this exists: *refinement density is not
//! probability density.* The quadtree concentrates compute at boundaries because they are
//! interesting, not because those ICs are more probable, so leaf density must never feed a
//! quantitative claim. What may is the local area scale factor of the chart→IC map — and the CPU
//! already computes the Jacobian it comes from, so it is free.
//!
//! **The spec calls it `|det J_D|` and that name is wrong here.** `J_D` maps a 2-plane into a
//! 12-dimensional state, so it is `12 × 2` and has no determinant. The area scale factor of a
//! 2-form is the square root of the Gram determinant, which reduces to `|det J|` exactly when the
//! target is 2-D — the case the phrase is written for.

use prin_rs::decode::linearise;
use prin_rs::grid::Chart;

/// **Against a finite-difference estimate of the same quantity**, computed a different way.
///
/// The weight is `sqrt(|ju|² |jv|² − (ju·jv)²)` from the *secant* columns. The independent estimate
/// is the area of the parallelogram spanned by two decoded *displacements*, built by differencing
/// decoded states directly rather than through `Lin`. On an affine chart the two are the same
/// number by construction and agreeing to round-off is the check that the Gram algebra is right;
/// on a curved one they agree to `O(h²)`, and the ladder below is what says which regime a chart
/// is in rather than assuming.
#[test]
fn the_weight_matches_an_independent_area_estimate() {
    let chart = Chart::BodyPlane;
    let (cu, cv) = (1.0, 3.0);

    // The independent path: decode four points and form the parallelogram area in the full state
    // space, with no reference to `Lin`.
    let area_fd = |h: f64| {
        let d = |a: (f64, f64), b: (f64, f64)| {
            let (sa, sb) = (
                prin_rs::grid::decode_state(&chart, 0, a.0, a.1).s,
                prin_rs::grid::decode_state(&chart, 0, b.0, b.1).s,
            );
            let mut o = [0.0f64; 12];
            for k in 0..3 {
                o[k * 4] = (sa.r[k].x - sb.r[k].x) * 0.5;
                o[k * 4 + 1] = (sa.r[k].y - sb.r[k].y) * 0.5;
                o[k * 4 + 2] = (sa.v[k].x - sb.v[k].x) * 0.5;
                o[k * 4 + 3] = (sa.v[k].y - sb.v[k].y) * 0.5;
            }
            o
        };
        let du = d((cu + h, cv), (cu - h, cv));
        let dv = d((cu, cv + h), (cu, cv - h));
        let dot = |a: &[f64; 12], b: &[f64; 12]| (0..12).map(|i| a[i] * b[i]).sum::<f64>();
        (dot(&du, &du) * dot(&dv, &dv) - dot(&du, &dv).powi(2)).max(0.0).sqrt()
    };

    let mut checked = 0;
    for h in [1e-2, 1e-3, 1e-4] {
        let w = linearise(&chart, 0, cu, cv, h).measure_weight();
        let fd = area_fd(h);
        assert!(w > 0.0, "the weight vanished at h={h}, so this chart carries no area");
        let rel = (w - fd).abs() / fd.max(1e-300);
        assert!(rel < 1e-9, "h={h}: weight {w:e} against the independent area {fd:e}, rel {rel:e}");
        checked += 1;
    }
    assert_eq!(checked, 3);
}

/// **The weight scales as `h²`**, because it is an area and the columns each carry the half-width.
///
/// The arm that makes this non-trivial: a weight that was secretly a *length* would scale as `h`,
/// and a constant would not scale at all. Both would pass an "it is positive" check.
#[test]
fn the_weight_is_an_area_and_scales_as_h_squared() {
    let chart = Chart::BodyPlane;
    let w = |h: f64| linearise(&chart, 0, 1.0, 3.0, h).measure_weight();

    let (a, b) = (w(1e-2), w(1e-3));
    let ratio = a / b;
    assert!(
        (ratio - 100.0).abs() / 100.0 < 1e-6,
        "a tenfold h must give a hundredfold area, got {ratio}"
    );
    // Not a length, and not a constant.
    assert!((ratio - 10.0).abs() > 1.0, "the weight scales as h, so it is a length not an area");
    assert!((ratio - 1.0).abs() > 1.0, "the weight does not scale, so it is not a measure at all");
}

/// **A degenerate plane carries no area, and the weight says so rather than being floored.**
///
/// If the two columns are parallel the plane has collapsed to a curve. Zero is the true answer and
/// a floored small number would be a quantitative claim about a measure that is not there — which
/// is exactly the failure the whole weight exists to prevent.
#[test]
fn a_collapsed_plane_weighs_nothing() {
    use prin_rs::decode::Lin;
    use prin_rs::physics::Cart;

    let mut col = Cart::<f64>::default();
    col.r[0].x = 1.0;
    let parallel = Lin { x0: Cart::default(), ju: col, jv: col };
    assert_eq!(parallel.measure_weight(), 0.0, "parallel columns span no area");

    let mut other = Cart::<f64>::default();
    other.r[0].y = 1.0;
    let square = Lin { x0: Cart::default(), ju: col, jv: other };
    assert!((square.measure_weight() - 1.0).abs() < 1e-15, "orthogonal unit columns span area 1");
}
