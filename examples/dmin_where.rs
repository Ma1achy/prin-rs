//! Where, exactly, is d_min attained? Mirrors the driver's inner loop and records the
//! sequence of |R1| around the minimum, so the mechanism can be read off rather than guessed.
use prin_rs::integrate::az::rk4;
use prin_rs::integrate::az::{reference_body::triple, AzSystem};
use prin_rs::physics::Cart;
use prin_rs::Vec2;

fn setup() -> (Cart<f64>, [f64; 3]) {
    (
        Cart::new(
            [Vec2::new(-0.5, 0.0), Vec2::new(0.5, 0.0), Vec2::new(0.1, 1000.0)],
            [Vec2::zero(); 3],
        ),
        [1.0, 1.0, 1.0],
    )
}

fn probe(eta: f64) {
    let (cart, m) = setup();
    let (a, b, c) = triple(1); // reference body 1, as the gate asserts
    let sys = AzSystem::new(a, b, c, m).with_reference_lc();
    let (mut s, e) = sys.to_reg(&cart);

    let dt_left = 1.0f64;
    let dtau = eta * dt_left / (s.a().max(f64::MIN_POSITIVE) * s.b().max(f64::MIN_POSITIVE));

    let mut best = f64::INFINITY;
    let mut best_i = 0usize;
    let mut best_tau = 0.0;
    let mut best_t = 0.0;
    let mut hist: Vec<(usize, f64, f64, f64)> = Vec::new();
    let mut i = 0usize;
    let mut tau = 0.0f64;

    while s.t < dt_left && s.is_finite() && i < 40_000_000 {
        s = rk4::step(&sys, &s, e, dtau);
        i += 1;
        tau += dtau;
        let (r1, _, _, _) = sys.phys_from_state(&s);
        let d1 = r1.norm();
        hist.push((i, tau, s.t, d1));
        if d1 < best {
            best = d1;
            best_i = i;
            best_tau = tau;
            best_t = s.t;
        }
    }

    println!("eta = {eta:.0e}:  dtau = {dtau:.6e}, steps = {i}");
    println!("  d_min = {best:.9e} at step {best_i}, tau = {best_tau:.9e}, t = {best_t:.9}");
    println!("  neighbouring samples of |R1| (step, tau, |R1|):");
    for k in best_i.saturating_sub(2)..=(best_i + 2).min(hist.len()) {
        if k == 0 { continue; }
        let (si, st, _, d) = hist[k - 1];
        let mark = if si == best_i { "  <- min" } else { "" };
        println!("    {si:>10}  {st:.9e}  {d:.9e}{mark}");
    }
    println!();
}

fn main() {
    println!("Instrumented inner loop, n_sync=1, reference LC branch (as gate (b) runs it).\n");
    for eta in [1e-4f64, 1e-5] {
        probe(eta);
    }
}
