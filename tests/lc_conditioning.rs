//! Conditioning of the Levi-Civita inverse.
//!
//! `u_of_rho` computes `u0 = sqrt((|rho| + rho.x)/2)`. When `rho` points along **negative
//! x**, `|rho| + rho.x` is a difference of near-equal numbers: catastrophic cancellation.
//! `u1 = rho.y / (2 u0)` then divides by that damaged value and amplifies it.
//!
//! The reference guards only with `u0 > 1e-300`, which prevents division by zero and does
//! nothing for conditioning. This is inherent to the branch, not a transcription error — but
//! it sets a floor on registration accuracy, and registration happens at every one of the
//! `n_sync` boundaries.
//!
//! Reported rather than asserted tight: the point is to know the number.

use prin_rs::integrate::az::lc;
use prin_rs::Vec2;

#[test]
fn lc_inverse_conditioning_versus_angle() {
    println!("{:>12}{:>18}{:>18}{:>12}", "angle (deg)", "reference branch", "stable branch", "u0/|u|");
    let mut worst_ref = 0.0f64;
    let mut worst_stable = 0.0f64;
    for deg in [0.0f64, 45.0, 90.0, 135.0, 170.0, 179.0, 179.9, 179.99, 180.0] {
        let th = deg.to_radians();
        let rho = Vec2::new(th.cos(), th.sin());

        let ur = lc::u_of_rho_reference(rho);
        let rel_r = (lc::rho_of_u(ur) - rho).norm() / rho.norm();
        let us = lc::u_of_rho(rho);
        let rel_s = (lc::rho_of_u(us) - rho).norm() / rho.norm();

        println!("{deg:>12}{rel_r:>18.3e}{rel_s:>18.3e}{:>12.3e}", ur.x / ur.norm().max(1e-300));
        worst_ref = worst_ref.max(rel_r);
        worst_stable = worst_stable.max(rel_s);
    }
    println!("\nreference worst = {worst_ref:.3e}, stable worst = {worst_stable:.3e}");
    println!("The reference always computes u0 first; that sum cancels when rho points along");
    println!("-x. The stable branch computes whichever component is larger and derives the");
    println!("other, which removes the loss entirely rather than merely bounding it.");

    assert!(worst_ref > 1e-10, "the reference branch was expected to lose accuracy here");
    assert!(worst_stable < 1e-15, "the stable branch should be at roundoff: {worst_stable:e}");
}

/// The defect is not only precision: the reference's branch cut is fixed along negative x,
/// so its accuracy depends on the **absolute orientation** of a configuration. Rotating a
/// configuration is a symmetry of the physics and must not change the numerics.
#[test]
fn the_stable_branch_is_rotationally_uniform() {
    let mut worst_ref = 0.0f64;
    let mut worst_stable = 0.0f64;
    for k in 0..3600 {
        let th = (k as f64) * std::f64::consts::TAU / 3600.0;
        let rho = Vec2::new(2.5 * th.cos(), 2.5 * th.sin());
        worst_ref = worst_ref
            .max((lc::rho_of_u(lc::u_of_rho_reference(rho)) - rho).norm() / rho.norm());
        worst_stable = worst_stable
            .max((lc::rho_of_u(lc::u_of_rho(rho)) - rho).norm() / rho.norm());
    }
    println!("worst round-trip over 3600 orientations:");
    println!("  reference branch = {worst_ref:.3e}");
    println!("  stable branch    = {worst_stable:.3e}");
    assert!(worst_stable <= 1e-15, "stable branch is orientation-dependent: {worst_stable:e}");
}
