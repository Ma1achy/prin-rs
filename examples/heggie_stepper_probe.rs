//! Which steppers can march Heggie's regularised Hamiltonian?
//!
//! The Principia corpus names "Heggie with KDK leapfrog" as the GUI's default occupant, and
//! `src/integrate/mod.rs` says there is no Heggie + leapfrog because `Gamma*` couples position and
//! momentum. This probe measures what that combination actually does, against RK4 (the stepper
//! every Heggie result here was taken under) and against the two genuinely symplectic options for
//! a non-separable Hamiltonian, both implicit.
//!
//! Arms, all in fictitious time at a FIXED step `h` (no step control), so the stepper alone is
//! measured:
//!
//!   rk4        classic RK4, 4 deriv calls per step
//!   kdk_naive  the explicit drift-kick-drift a reader of "Heggie with KDK" would write:
//!              kick p at (u, p), drift u and t at (u, p_half), kick p at (u1, p_half). 3 calls.
//!   gen_lf     the generalised (Stormer-Verlet) leapfrog for a non-separable H: implicit in
//!              p_half and in u1, solved by fixed-point iteration. Symplectic, symmetric, order 2.
//!   imp_mid    the implicit midpoint rule, fixed-point iteration. Symplectic, symmetric, order 2.
//!
//! Two measurements:
//!   1. ORDER: state error at a fixed tau against a fine RK4 reference, over a halving ladder.
//!   2. ENERGY over ~100 periods: the regularised energy's error, max and final. A symplectic
//!      method's stays bounded; a non-symplectic one's drifts.
//!
//! Measured (2026-10-21): rk4 order 3.99, bounded energy; kdk_naive order 0.8-1.2, diverges within 100 orbits;
//! imp_mid order 2.00, bounded, 9-10 calls/step; gen_lf order 2 over one orbit but **unbounded over 100 and worse
//! at smaller steps** -- that smell usually means this implementation, not the method, and it is NOT explained.
//! `DBG=1` prints gen_lf's fixed-point cap hits and first bad step.
//!
//! Fixture: the Chenciner-Montgomery figure-eight (equal masses, exactly periodic, no close
//! approach), under both time transformations.

use prin_rs::integrate::heggie::hamiltonian::deriv_time;
use prin_rs::integrate::heggie::{HgState, HgSystem, HgTime};
use prin_rs::physics::Cart;
use prin_rs::Vec2;

fn figure_eight() -> Cart<f64> {
    let (x, y) = (0.97000436, -0.24308753);
    let (vx, vy) = (-0.93240737, -0.86473146);
    Cart::new(
        [Vec2::new(x, y), Vec2::new(0.0, 0.0), Vec2::new(-x, -y)],
        [
            Vec2::new(-vx / 2.0, -vy / 2.0),
            Vec2::new(vx, vy),
            Vec2::new(-vx / 2.0, -vy / 2.0),
        ],
    )
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Arm {
    Rk4,
    KdkNaive,
    GenLf,
    ImpMid,
}

struct Ctx<'a> {
    sys: &'a HgSystem<f64>,
    h0: f64,
    time: HgTime,
    calls: u64,
    max_iters: u32,
    caps: u64,
}

impl Ctx<'_> {
    fn f(&mut self, s: &HgState<f64>) -> HgState<f64> {
        self.calls += 1;
        deriv_time(self.sys, s, self.h0, self.time)
    }
}

fn mid(a: &HgState<f64>, b: &HgState<f64>) -> HgState<f64> {
    let (x, y) = (a.to_array13(), b.to_array13());
    let mut m = [0.0; 13];
    for k in 0..13 {
        m[k] = 0.5 * (x[k] + y[k]);
    }
    HgState::from_array13(m)
}

fn diff(a: &HgState<f64>, b: &HgState<f64>) -> f64 {
    let (x, y) = (a.to_array13(), b.to_array13());
    (0..13).map(|k| (x[k] - y[k]).abs()).fold(0.0, f64::max)
}

/// Replace the p-part of `s` with `s.p + c * d.p`.
fn kick(s: &HgState<f64>, d: &HgState<f64>, c: f64) -> HgState<f64> {
    let mut o = *s;
    for i in 0..3 {
        o.p[i] = s.p[i] + d.p[i] * c;
    }
    o
}

/// Replace the u- and t-parts of `s` with `s.u + c * d.u`, `s.t + c * d.t`.
fn drift(s: &HgState<f64>, d: &HgState<f64>, c: f64) -> HgState<f64> {
    let mut o = *s;
    for i in 0..3 {
        o.u[i] = s.u[i] + d.u[i] * c;
    }
    o.t = s.t + d.t * c;
    o
}

const TOL: f64 = 1e-15;

fn step(ctx: &mut Ctx, arm: Arm, s: &HgState<f64>, h: f64) -> HgState<f64> {
    match arm {
        Arm::Rk4 => {
            let k1 = ctx.f(s);
            let k2 = ctx.f(&s.axpy(0.5 * h, &k1));
            let k3 = ctx.f(&s.axpy(0.5 * h, &k2));
            let k4 = ctx.f(&s.axpy(h, &k3));
            let mut acc = HgState { u: [Vec2::new(0.0, 0.0); 3], p: [Vec2::new(0.0, 0.0); 3], t: 0.0 };
            for i in 0..3 {
                acc.u[i] = k1.u[i] + k2.u[i] * 2.0 + k3.u[i] * 2.0 + k4.u[i];
                acc.p[i] = k1.p[i] + k2.p[i] * 2.0 + k3.p[i] * 2.0 + k4.p[i];
            }
            acc.t = k1.t + k2.t * 2.0 + k3.t * 2.0 + k4.t;
            s.axpy(h / 6.0, &acc)
        }
        Arm::KdkNaive => {
            let d1 = ctx.f(s);
            let a = kick(s, &d1, 0.5 * h);
            let d2 = ctx.f(&a);
            let b = drift(&a, &d2, h);
            let d3 = ctx.f(&b);
            kick(&b, &d3, 0.5 * h)
        }
        Arm::GenLf => {
            // p_half = p + h/2 * pdot(u, p_half)            (implicit in p_half)
            let mut a = kick(s, &ctx.f(s), 0.5 * h);
            let mut ok = false;
            for _ in 0..ctx.max_iters {
                let n = kick(s, &ctx.f(&a), 0.5 * h);
                let done = diff(&n, &a) <= TOL;
                a = n;
                if done {
                    ok = true;
                    break;
                }
            }
            if !ok { ctx.caps += 1; }
            // u1 = u + h/2 * (udot(u, p_half) + udot(u1, p_half))   (implicit in u1)
            let da = ctx.f(&a);
            let mut b = drift(&a, &da, h);
            for _ in 0..ctx.max_iters {
                let db = ctx.f(&b);
                let mut avg = da;
                for i in 0..3 {
                    avg.u[i] = (da.u[i] + db.u[i]) * 0.5;
                }
                avg.t = 0.5 * (da.t + db.t);
                let n = drift(&a, &avg, h);
                let done = diff(&n, &b) <= TOL;
                b = n;
                if done {
                    ok = true;
                    break;
                }
                ok = false;
            }
            if !ok { ctx.caps += 1; }
            // p1 = p_half + h/2 * pdot(u1, p_half)          (explicit)
            let db = ctx.f(&b);
            kick(&b, &db, 0.5 * h)
        }
        Arm::ImpMid => {
            let mut y = s.axpy(h, &ctx.f(s));
            for _ in 0..ctx.max_iters {
                let n = s.axpy(h, &ctx.f(&mid(s, &y)));
                let done = diff(&n, &y) <= TOL;
                y = n;
                if done {
                    break;
                }
            }
            y
        }
    }
}

fn march(sys: &HgSystem<f64>, s0: &HgState<f64>, h0: f64, time: HgTime, arm: Arm, h: f64, tau: f64)
    -> (HgState<f64>, u64, f64) {
    let n = (tau / h).round() as usize;
    let mut ctx = Ctx { sys, h0, time, calls: 0, max_iters: 60, caps: 0 };
    let mut s = *s0;
    let mut emax: f64 = 0.0;
    let mut first_bad = None;
    for k in 0..n {
        s = step(&mut ctx, arm, &s, h);
        let e = ((sys.energy_of(&s) - h0) / h0.abs()).abs();
        emax = emax.max(e);
        if first_bad.is_none() && !(e < 1e-2) { first_bad = Some((k, ctx.caps, s.r(0).min(s.r(1)).min(s.r(2)))); }
    }
    if arm == Arm::GenLf && std::env::var("DBG").is_ok() { eprintln!("    gen_lf h={h:.4}: cap hits {}, first |dE/E|>1e-2 at step/caps/minR {:?}", ctx.caps, first_bad); }
    (s, ctx.calls, emax)
}

fn main() {
    let cart = figure_eight();
    let m = [1.0, 1.0, 1.0];
    let sys = HgSystem::new(m);
    let (s0, h0) = sys.to_reg(&cart);
    let arms = [Arm::Rk4, Arm::KdkNaive, Arm::GenLf, Arm::ImpMid];

    for (tname, time) in [("Eq.20 product", HgTime::Product), ("Eq.22 n=3/2 (production)", HgTime::default())] {
        // tau for one period: march a fine RK4 until t passes T.
        let period = 6.325_913_98;
        let mut ctx = Ctx { sys: &sys, h0, time, calls: 0, max_iters: 0, caps: 0 };
        let (mut s, mut tau, hh) = (s0, 0.0, 1e-4);
        while s.t < period {
            s = step(&mut ctx, Arm::Rk4, &s, hh);
            tau += hh;
        }
        let tau_end = (tau * 64.0).round() / 64.0;
        println!("\n=== time transformation: {tname}; one period is tau = {tau:.4}, using tau_end = {tau_end}");

        // 1. order
        let hs = [tau_end / 64.0, tau_end / 128.0, tau_end / 256.0, tau_end / 512.0, tau_end / 1024.0];
        let (sref, _, _) = march(&sys, &s0, h0, time, Arm::Rk4, tau_end / 16384.0, tau_end);
        println!("  ORDER (state error at tau_end against RK4 at h/16 of the finest)");
        println!("  {:>10} {:>12} {:>12} {:>12} {:>12}", "h", "rk4", "kdk_naive", "gen_lf", "imp_mid");
        let mut errs = vec![vec![0.0; hs.len()]; arms.len()];
        let mut calls = vec![vec![0u64; hs.len()]; arms.len()];
        for (j, &h) in hs.iter().enumerate() {
            for (a, &arm) in arms.iter().enumerate() {
                let (s1, c, _) = march(&sys, &s0, h0, time, arm, h, tau_end);
                errs[a][j] = diff(&s1, &sref);
                calls[a][j] = c;
            }
            println!("  {:>10.5} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e}", h, errs[0][j], errs[1][j], errs[2][j], errs[3][j]);
        }
        let ord = |e: &Vec<f64>| (e[0] / e[e.len() - 1]).log2() / ((e.len() - 1) as f64);
        println!("  {:>10} {:>12.2} {:>12.2} {:>12.2} {:>12.2}", "order", ord(&errs[0]), ord(&errs[1]), ord(&errs[2]), ord(&errs[3]));
        println!("  {:>10} {:>12.1} {:>12.1} {:>12.1} {:>12.1}", "calls/step",
            calls[0][2] as f64 / (tau_end / hs[2]), calls[1][2] as f64 / (tau_end / hs[2]),
            calls[2][2] as f64 / (tau_end / hs[2]), calls[3][2] as f64 / (tau_end / hs[2]));

        // 2b. is gen_lf's long-run growth a step-size instability? the same 100 periods, finer
        for div in [256.0, 512.0, 1024.0] {
            let (s1, c, emax) = march(&sys, &s0, h0, time, Arm::GenLf, tau_end / div, 100.0 * tau_end);
            println!("  gen_lf at h = tau_end/{div}: max |dE/E| {:.3e}, calls/step {:.1}, t reached {:.2}", emax, c as f64 / (100.0 * div), s1.t);
        }
        // 2. energy over ~100 periods at a moderate step
        let h = tau_end / 128.0;
        println!("  ENERGY over 100 periods at h = tau_end/128 (regularised energy, |dE/E|)");
        println!("  {:>10} {:>12} {:>12} {:>12} {:>12}", "arm", "max", "final", "calls", "t reached");
        for (a, &arm) in arms.iter().enumerate() {
            let (s1, c, emax) = march(&sys, &s0, h0, time, arm, h, 100.0 * tau_end);
            let ef = ((sys.energy_of(&s1) - h0) / h0.abs()).abs();
            println!("  {:>10} {:>12.3e} {:>12.3e} {:>12} {:>12.2}", ["rk4", "kdk_naive", "gen_lf", "imp_mid"][a], emax, ef, c, s1.t);
        }
    }
}
