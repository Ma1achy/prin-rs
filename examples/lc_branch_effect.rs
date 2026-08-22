//! What the stable inverse LC branch buys the kernel.
//!
//! The cross-check pins the reference branch, so this is the only place the change is
//! measured. The comparison is on integration quality — drift and the Gamma residual — not
//! on agreement with the reference, which the change deliberately gives up.
use prin_rs::grid;
use prin_rs::integrate::az;
use prin_rs::physics::burrau;

fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * q).round() as usize]
}

/// Generic over precision. Initial conditions are always built in f64 and cast down, so an
/// f32/f64 difference here is arithmetic and never the initial conditions.
fn run_t<T: prin_rs::Real>(t_max: f64, stable: bool) -> (f64, f64, f64) {
    let m = burrau::masses::<T>();
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    let n_sync = ((t_max / (13.0 / 32.0)).round() as usize).max(1);
    let mut drift = Vec::new();
    let mut gam = Vec::new();
    for i in 0..s.npix() {
        let o = az::integrate_az_lc(
            s.nominal::<T>(i), &m, T::lit(t_max), n_sync, T::lit(0.01), 30_000, None, stable,
        );
        drift.push(o.drift.to_f64().unwrap());
        gam.push(o.gamma_max.to_f64().unwrap());
    }
    let med = quantile(&mut drift, 0.5);
    let mx = quantile(&mut drift, 1.0);
    let g = quantile(&mut gam, 0.5);
    (med, mx, g)
}

fn run(t_max: f64, stable: bool) -> (f64, f64, f64) {
    run_t::<f64>(t_max, stable)
}

fn main() {
    println!("near-field 8x8 nominal, eta=0.01, n_sync scaled to a fixed 13/32 sub-interval");
    println!();
    println!(
        "{:>7} | {:>11} {:>11} {:>7} | {:>11} {:>11} {:>7}",
        "t_max", "med ref", "med stable", "gain", "max ref", "max stable", "gain"
    );
    println!("{}", "-".repeat(76));
    for t_max in [0.5f64, 1.0, 2.0, 4.0, 8.0, 13.0] {
        let (mr, xr, _) = run(t_max, false);
        let (ms, xs, _) = run(t_max, true);
        println!(
            "{t_max:>7} | {mr:>11.3e} {ms:>11.3e} {:>6.1}x | {xr:>11.3e} {xs:>11.3e} {:>6.1}x",
            mr / ms,
            xr / xs
        );
    }
    println!();
    println!("=== f32 preview (Step 6 proper will do this correctly) ===");
    println!("The f64 gain vanishes past t=2 because larger errors mask it. At f32 the same");
    println!("cancellation is ~5e8 times worse, so the question is whether it is masked there");
    println!("too. Initial conditions are built in f64 and cast down, so this isolates");
    println!("arithmetic from IC differences.");
    println!();
    println!("{:>7} | {:>13} {:>13} {:>8}", "t_max", "f32 reference", "f32 stable", "gain");
    println!("{}", "-".repeat(46));
    for t_max in [0.5f64, 1.0, 2.0, 4.0, 13.0] {
        let (mr, _, _) = run_t::<f32>(t_max, false);
        let (ms, _, _) = run_t::<f32>(t_max, true);
        println!("{t_max:>7} | {mr:>13.3e} {ms:>13.3e} {:>7.1}x", mr / ms);
    }
    println!();
    println!("Gamma residual (median over the grid):");
    println!("{:>7} | {:>13} {:>13} {:>8}", "t_max", "reference", "stable", "gain");
    println!("{}", "-".repeat(46));
    for t_max in [1.0f64, 4.0, 13.0] {
        let (_, _, gr) = run(t_max, false);
        let (_, _, gs) = run(t_max, true);
        println!("{t_max:>7} | {gr:>13.3e} {gs:>13.3e} {:>7.1}x", gr / gs);
    }
}
