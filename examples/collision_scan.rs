//! Diagnostic: d_min is sampling-limited near an exact collision. It scales as eta^2
//! because u passes through zero linearly in fictitious time and |R1| = |u1|^2, so the
//! closest *sampled* separation goes as (step in u)^2. The physical minimum is zero.
use prin_rs::integrate::az;
use prin_rs::physics::Cart;
use prin_rs::Vec2;

fn setup(far: f64) -> (Cart<f64>, [f64; 3]) {
    (
        Cart::new(
            [Vec2::new(-0.5, 0.0), Vec2::new(0.5, 0.0), Vec2::new(0.1, far)],
            [Vec2::zero(); 3],
        ),
        [1.0, 1.0, 1.0],
    )
}

fn main() {
    println!("far = 100, n_sync = 1 (one registration, no re-registration loss)");
    println!("{:>10}{:>14}{:>14}{:>14}{:>10}", "eta", "d_min", "|dE/E|", "gamma_max", "steps");
    let mut prev: Option<(f64, f64)> = None;
    for eta in [1e-2f64, 1e-3, 1e-4, 1e-5, 3e-6, 1e-6] {
        let (s, m) = setup(100.0);
        let o = az::integrate_az(s, &m, 1.0, 1, eta, 20_000_000, None);
        let scal = prev.map(|(pe, pd): (f64, f64)| {
            let want = (eta / pe).powi(2);
            format!("  d_min ratio {:.2e}, eta^2 predicts {:.2e}", o.d_min_ref / pd, want)
        }).unwrap_or_default();
        println!(
            "{eta:>10.0e}{:>14.3e}{:>14.3e}{:>14.3e}{:>10}{scal}",
            o.d_min_ref, o.drift, o.gamma_max, o.steps
        );
        prev = Some((eta, o.d_min_ref));
    }
}
