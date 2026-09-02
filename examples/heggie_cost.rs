//! **What does a Heggie step actually cost, on this machine, in this port?**
//!
//! `examples/heggie_machinery.rs` reports `steps`, which is the machine-independent cost proxy
//! this project insists on — *read `steps`, not `secs`*, because under load a run doing 1.7% more
//! work timed faster than the baseline. But `steps` is only half the accounting: it says how many
//! steps were taken and nothing about what one costs, and Heggie's `Gamma*` has ten terms against
//! AZ's seven, three coupled vectors against two, and thirteen state components against nine.
//!
//! Heggie's §3 measures "a factor of about 1.6" in computing time per step against Aarseth-Zare.
//! That is his number, on his machine, in 1974, on a Runge-Kutta-Fehlberg 7(8) implementation in
//! three dimensions. Carrying it as though it described this port would be quoting where a
//! measurement is available.
//!
//! # Why this is allowed to be a timing measurement
//!
//! The standing rule is about comparing *trajectories*, where step counts differ, thread
//! contention varies and wall clock scores the machine rather than the work. Here the two sides
//! do **exactly the same number of calls** on states of the same magnitude, single-threaded, with
//! nothing else running. There is no work-count difference for the clock to confound, which is
//! the one case where seconds are the honest unit.
//!
//! `deriv` is the right subject: RK4 calls it four times per step and does almost nothing else,
//! so the per-call ratio is the per-step ratio to within the accumulation.
//!
//! Args: `iters`.

use std::hint::black_box;
use std::time::Instant;

use prin_rs::integrate::az::{hamiltonian as az_h, reference_body::triple, AzState, AzSystem};
use prin_rs::integrate::heggie::{hamiltonian as hg_h, HgState, HgSystem, HgTime};
use prin_rs::physics::burrau;
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

fn main() {
    let iters: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20_000_000);

    let m = burrau::masses::<f64>();
    let (a, b, c) = triple(0);
    let az_sys = AzSystem::new(a, b, c, m);
    let hg_sys = HgSystem::new(m);

    // A pool of states rather than one, so the loop cannot be hoisted and the branch predictor
    // cannot learn a single path. Same generator and same magnitudes for both sides.
    let mut rng = SplitMix64::new(0xC057);
    const POOL: usize = 64;
    let az_states: Vec<AzState<f64>> = (0..POOL)
        .map(|_| AzState {
            u1: Vec2::new(rng.range(0.2, 2.0), rng.range(-2.0, 2.0)),
            p1: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
            u2: Vec2::new(rng.range(0.2, 2.0), rng.range(-2.0, 2.0)),
            p2: Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)),
            t: 0.0,
        })
        .collect();
    let hg_states: Vec<HgState<f64>> = (0..POOL)
        .map(|_| HgState {
            u: std::array::from_fn(|_| Vec2::new(rng.range(0.2, 2.0), rng.range(-2.0, 2.0))),
            p: std::array::from_fn(|_| Vec2::new(rng.range(-2.0, 2.0), rng.range(-2.0, 2.0))),
            t: 0.0,
        })
        .collect();
    let e = -12.8167;

    // Warm up, so neither side pays for a cold cache or a first-touch page fault.
    for i in 0..POOL {
        black_box(az_h::deriv(&az_sys, &az_states[i], e));
        black_box(hg_h::deriv(&hg_sys, &hg_states[i], e));
        black_box(hg_h::deriv_time(&hg_sys, &hg_states[i], e, HgTime::default()));
    }

    let row = |name: &str, ns: f64| {
        println!("  {name:>34} {ns:>10.2} ns/call");
        ns
    };

    // Interleaved repeats, so a drift in machine state over the run cannot land on one side.
    let (mut t_az, mut t_p, mut t_s) = (f64::MAX, f64::MAX, f64::MAX);
    for _ in 0..3 {
        let t0 = Instant::now();
        for i in 0..iters {
            black_box(az_h::deriv(&az_sys, &az_states[i % POOL], e));
        }
        t_az = t_az.min(t0.elapsed().as_secs_f64() / iters as f64 * 1e9);

        let t0 = Instant::now();
        for i in 0..iters {
            black_box(hg_h::deriv(&hg_sys, &hg_states[i % POOL], e));
        }
        t_p = t_p.min(t0.elapsed().as_secs_f64() / iters as f64 * 1e9);

        let t0 = Instant::now();
        for i in 0..iters {
            black_box(hg_h::deriv_time(&hg_sys, &hg_states[i % POOL], e, HgTime::default()));
        }
        t_s = t_s.min(t0.elapsed().as_secs_f64() / iters as f64 * 1e9);
    }

    println!("{iters} calls per timing, best of 3 interleaved repeats, single-threaded.\n");
    let az = row("AZ  deriv (Gamma)", t_az);
    let hp = row("HG  deriv, Eq. (20)/(21)", t_p);
    let hs = row("HG  deriv_time, Eqs. (22)-(24)", t_s);

    // **The control for the comparison above, and it is not optional.**
    //
    // AZ's `deriv` computes `r3.powf(3.0)` rather than `r3*r3*r3`, deliberately: `powf` routes to
    // the same libm call NumPy's `**3` does, removing one ulp-level divergence source from the
    // cross-check. That is a **fidelity** choice, not algebra, and a libm `powf` is far more
    // expensive than two multiplies. Reporting a Heggie cost win without separating it would be
    // scoring AZ's cross-check discipline as though it were the method.
    //
    // Timed in isolation and reported as an amount rather than folded in, because the compiler
    // may schedule it differently inside `deriv` than in a bare loop -- so it bounds the effect,
    // it does not subtract it exactly.
    let xs: Vec<f64> = (0..POOL).map(|i| 0.1 + i as f64 * 0.05).collect();
    let (mut t_pw, mut t_mul) = (f64::MAX, f64::MAX);
    for _ in 0..3 {
        let t0 = Instant::now();
        for i in 0..iters {
            black_box(black_box(xs[i % POOL]).powf(black_box(3.0)));
        }
        t_pw = t_pw.min(t0.elapsed().as_secs_f64() / iters as f64 * 1e9);
        let t0 = Instant::now();
        for i in 0..iters {
            let x = black_box(xs[i % POOL]);
            black_box(x * x * x);
        }
        t_mul = t_mul.min(t0.elapsed().as_secs_f64() / iters as f64 * 1e9);
    }
    println!();
    row("powf(3.0), as AZ's deriv calls it", t_pw);
    row("x*x*x, the same value", t_mul);
    println!(
        "  -> up to {:.2} ns/call of AZ's {az:.2} is the ulp-fidelity choice, not the algebra",
        t_pw - t_mul
    );

    println!("\n  Heggie Eq. (20)/(21) against AZ : {:.2}x", hp / az);
    println!("  Heggie Eqs. (22)-(24) against AZ: {:.2}x", hs / az);
    println!("  the extra cost of the control term and S^-3/2: {:.2}x", hs / hp);

    println!(
        "\nAgainst Heggie's own §3: he measures ~1.6x per step for Eqs. (20)/(21) and a further\n\
         ~15% for (22)-(24), on RKF7(8) in three dimensions in 1974. The planar reduction is why\n\
         this port can differ -- his 4-vectors collapse to 2-vectors and his 4x3 A_i to a 2x2\n\
         block, so a chunk of the arithmetic he paid for does not exist here.\n\n\
         **Read the `powf` row before the ratios.** AZ pays a libm call per `deriv` for bitwise\n\
         agreement with NumPy; that is cross-check discipline and not the method, so a Heggie\n\
         win smaller than that row is not a win at all.\n\n\
         **Multiply this by the step ratio for the total.** `heggie_machinery` measures 1.23x the\n\
         steps on `config_stability`; the compute ratio is that times the number above, and\n\
         neither factor alone is the cost."
    );
}
