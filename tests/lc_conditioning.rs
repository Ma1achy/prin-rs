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
    println!("{:>12}{:>16}{:>16}", "angle (deg)", "round-trip rel", "u0/|u|");
    let mut worst = 0.0f64;
    for deg in [0.0f64, 45.0, 90.0, 135.0, 170.0, 179.0, 179.9, 179.99, 180.0] {
        let th = deg.to_radians();
        let rho = Vec2::new(th.cos(), th.sin());
        let u = lc::u_of_rho(rho);
        let back = lc::rho_of_u(u);
        let rel = (back - rho).norm() / rho.norm();
        let frac = u.x / u.norm().max(1e-300);
        println!("{deg:>12}{rel:>16.3e}{frac:>16.3e}");
        worst = worst.max(rel);
    }
    println!("\nThe loss is concentrated where rho points along -x: |rho| + rho.x cancels.");
    println!("This is a property of the LC branch, present identically in the reference.");
    // Loose: this documents a conditioning floor, it is not a correctness gate.
    assert!(worst < 1e-7, "LC round trip lost more than expected: {worst:e}");
}
