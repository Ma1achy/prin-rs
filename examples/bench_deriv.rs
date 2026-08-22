//! Throughput of the AZ derivative, and an extrapolation to the §9 target render.
use prin_rs::grid;
use prin_rs::integrate::az;
use prin_rs::physics::burrau;
use std::time::Instant;

fn main() {
    let m = burrau::masses::<f64>();
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    let n = s.npix();

    let t0 = Instant::now();
    let mut steps = 0usize;
    for i in 0..n {
        let o = az::integrate_az(s.nominal::<f64>(i), &m, 13.0, 32, 0.01, 30_000, None);
        steps += o.steps;
    }
    let el = t0.elapsed().as_secs_f64();

    println!("near-field {n} nominal trajectories to t=13, n_sync=32, eta=0.01");
    println!("  wall clock      : {el:.3} s");
    println!("  total RK4 steps : {steps}");
    println!("  per trajectory  : {:.2} ms, {} steps", 1e3 * el / n as f64, steps / n);
    println!("  RK4 steps/s     : {:.2e}  (4 deriv evaluations each)", steps as f64 / el);

    let target = 1024.0 * 1024.0 * 8.0;
    let secs = el / n as f64 * target;
    println!();
    println!("Extrapolated to 1024x1024 with E+1 = 8 copies ({:.1e} trajectories):", target);
    println!("  single core : {:.0} s  ({:.1} h)", secs, secs / 3600.0);
    let cores = std::thread::available_parallelism().map(|c| c.get()).unwrap_or(1);
    println!("  {cores} cores    : {:.0} s  ({:.2} h)  at perfect scaling", secs / cores as f64,
             secs / cores as f64 / 3600.0);
    println!();
    println!("Perf levers, in order: threading (rayon over pixels), f32, n_sync.");
    println!("The algorithm is not a lever — it is the thing being measured.");
}
