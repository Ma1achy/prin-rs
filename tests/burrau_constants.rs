//! BRIEF §2.1 reference constants. Exact, cheap, and the anchor for everything downstream:
//! the AZ validation chain in Step 4 bottoms out here, via `energy_phys == energy`.

use prin_rs::physics::{burrau, energy};

const TOL: f64 = 5e-5; // the brief quotes these to 4-5 significant figures

#[test]
fn burrau_reference_constants() {
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();

    let mtot: f64 = m.iter().sum();
    let e = energy::energy(&s.r, &s.v, &m, 0.0);
    let rg = energy::hyperradius(&s.r, &m);
    let tc = energy::crossing_time(&s.r, &m);

    println!("M  = {mtot}");
    println!("R  = {rg}");
    println!("E  = {e}");
    println!("tc = {tc}");

    assert!((mtot - 12.0).abs() < 1e-15, "M = {mtot}, want 12");
    assert!((rg - 2.2361).abs() < TOL, "R = {rg}, want 2.2361");
    assert!((e - (-12.8167)).abs() < TOL, "E = {e}, want -12.8167");
    assert!((tc - 0.9652).abs() < TOL, "crossing time = {tc}, want 0.9652");
}

#[test]
fn burrau_is_at_rest_with_com_at_origin() {
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();

    // Released from rest: v = 0 for every body, so L_z = 0 for every copy. This is why
    // BRIEF §4 forbids an L_z analogue of error_ratio — sigma_Lz(0) = 0 makes it 0/0.
    for k in 0..3 {
        assert_eq!(s.v[k].norm_sq(), 0.0);
    }

    let c = energy::com(&s.r, &m);
    assert!(c.x.abs() < 1e-15 && c.y.abs() < 1e-15, "COM = {c:?}, want origin");
}

#[test]
fn hyperradius_is_scale_covariant() {
    // R -> alpha R exactly. This is the invariance the whole project quotients out; if it
    // fails here, canonical units are meaningless downstream.
    let m = burrau::masses::<f64>();
    let s = burrau::state::<f64>();
    let r0 = energy::hyperradius(&s.r, &m);

    for alpha in [0.25f64, 1.0, 4.0] {
        let scaled = [s.r[0] * alpha, s.r[1] * alpha, s.r[2] * alpha];
        let r = energy::hyperradius(&scaled, &m);
        let rel = ((r / alpha) - r0).abs() / r0;
        println!("alpha = {alpha}: R/alpha = {}, rel err = {rel:e}", r / alpha);
        assert!(rel < 1e-14, "alpha = {alpha}: R/alpha = {}, want {r0}", r / alpha);
    }
}
