//! `shape_vec` quotients out translation, rotation and scale. Those three invariances are
//! the foundation of the Step 5a gauge-invariance gate, so they are asserted directly here
//! rather than only end-to-end.

use prin_rs::physics::{burrau, shape};
use prin_rs::Vec2;

fn rotate(r: &[Vec2<f64>; 3], th: f64) -> [Vec2<f64>; 3] {
    let (c, s) = (th.cos(), th.sin());
    let f = |p: Vec2<f64>| Vec2::new(c * p.x - s * p.y, s * p.x + c * p.y);
    [f(r[0]), f(r[1]), f(r[2])]
}

fn max_dev(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max)
}

#[test]
fn shape_vec_is_a_unit_vector() {
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();
    let n = shape::shape_vec(&s.r, &m);
    let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    println!("n = {n:?}, |n| = {norm}");
    assert!((norm - 1.0).abs() < 1e-15, "|n| = {norm}");
}

#[test]
fn shape_vec_quotients_translation_rotation_scale() {
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();
    let n0 = shape::shape_vec(&s.r, &m);

    let shift = Vec2::new(3.7, -1.2);
    let t = [s.r[0] + shift, s.r[1] + shift, s.r[2] + shift];
    let dt = max_dev(&shape::shape_vec(&t, &m), &n0);
    println!("translation: max dev = {dt:e}");
    assert!(dt < 1e-14, "translation moved shape_vec by {dt:e}");

    for th in [0.3f64, 1.9, -2.4] {
        let d = max_dev(&shape::shape_vec(&rotate(&s.r, th), &m), &n0);
        println!("rotation {th}: max dev = {d:e}");
        assert!(d < 1e-14, "rotation {th} moved shape_vec by {d:e}");
    }

    for alpha in [0.25f64, 4.0, 16.0] {
        let sc = [s.r[0] * alpha, s.r[1] * alpha, s.r[2] * alpha];
        let d = max_dev(&shape::shape_vec(&sc, &m), &n0);
        println!("scale {alpha}: max dev = {d:e}");
        assert!(d < 1e-14, "scale {alpha} moved shape_vec by {d:e}");
    }
}

#[test]
fn spread_statistics_are_zero_on_identical_copies_and_bounded() {
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();
    let n = shape::shape_vec(&s.r, &m);

    let same = vec![n; 8];
    assert!(shape::svar(&same).abs() < 1e-15);
    assert!(shape::spread_shape(&same).abs() < 1e-15);

    // Antipodal split: spread_shape divides by the chord bound 2, so a maximally split
    // ensemble reads exactly 0.5 here (mean distance from centroid is 1). The point is that
    // the normalisation is a *geometric* constant, carrying no dependence on sigma_E(0) and
    // hence none on cell width — which is what keeps ensemble_spread free of the resolution
    // confound that afflicts error_ratio.
    let anti: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    let ss = shape::spread_shape(&anti);
    let sv = shape::svar(&anti);
    println!("antipodal: spread_shape = {ss}, svar = {sv}");
    assert!((ss - 0.5).abs() < 1e-15, "spread_shape = {ss}");
    assert!((sv - 1.0).abs() < 1e-15, "svar = {sv}");
}
