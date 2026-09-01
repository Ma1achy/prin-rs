//! **STAGE 2 — do chain coordinates buy digits on `far`, and only at f32?**
//!
//! Stage 0 measured that `far`'s AZ advantage **inverts** between precisions: AZ better by 0.890
//! decades at f64 on every pixel, Heggie better by 0.768 at f32 on every pixel, saturation guard
//! LIVE at both. That says the `far` result is precision-limited, and chain coordinates are the
//! published repair — hold the configuration as inter-particle vectors so a separation is a sum of
//! small quantities rather than a difference of large ones.
//!
//! **PREDICTION, recorded in `results/overnight/PREDICTIONS.md` before Stage 0 ran: chain helps at
//! f32 and is near-invisible at f64**, because it is a round-off fix and f64 has digits to spare.
//!
//! # The comparison is controlled to one thing
//!
//! Both arms run the **same** time transformation (`LogH`), the **same** stepper (RK4), the same
//! fixed fictitious step and the same number of steps, on the same initial conditions. The only
//! difference is whether the state is held as three positions or as two chain vectors. No sync
//! boundaries, no events, no step control — those would each put a second thing in the measurement.
//!
//! The chain ordering is selected **once** and frozen. Re-selecting is a re-registration, which is
//! the mechanism this whole investigation exists to isolate.
//!
//! # Guards
//!
//! - **`nonfin` per arm.** A win over a dead arm is not a win.
//! - **The f64 row is the control.** If chain moves f64 substantially, it is not behaving as a
//!   round-off fix and the f32 row cannot be attributed to precision.
//! - Both arms are compared against the **same f64 direct reference trajectory**, so "better"
//!   means closer to the truth rather than merely smaller.
//!
//! Args: `res span steps`.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::integrate::logh::chain::{rk4 as chain_rk4, ChainOrder, ChainState};
use prin_rs::integrate::logh::hamiltonian::LhTime;
use prin_rs::integrate::logh::state::LhState;
use prin_rs::integrate::logh::step::rk4 as direct_rk4;
use prin_rs::integrate::logh::system::LhSystem;
use prin_rs::physics::{energy, Cart, Ic};
use prin_rs::Real;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &[f64], p: f64) -> f64 {
    let mut w: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if w.is_empty() {
        return f64::NAN;
    }
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    w[(((w.len() - 1) as f64) * p).round() as usize]
}

/// Energy drift of a march held in **chain** coordinates.
fn drift_chain<T: Real>(c: &Cart<T>, m: &[T; 3], steps: usize, h: f64) -> f64 {
    let o = ChainOrder::select(&c.r);
    let mut s = ChainState::from_cart(c, o);
    let b = energy::potential_pos(&c.r, m, T::zero()) - energy::kinetic(&c.v, m);
    let e0 = energy::energy(&c.r, &c.v, m, T::zero());
    for _ in 0..steps {
        let (n, _) = chain_rk4(m, &s, o, b, LhTime::LogH, T::lit(h));
        s = n;
        if !s.is_finite() {
            return f64::INFINITY;
        }
    }
    let back = s.to_cart(m, o);
    let e1 = energy::energy(&back.r, &back.v, m, T::zero());
    ((e1 - e0) / e0).to_f64().unwrap().abs()
}

/// Energy drift of the same march held as **three positions**.
fn drift_direct<T: Real>(c: &Cart<T>, m: &[T; 3], steps: usize, h: f64) -> f64 {
    let sys = LhSystem::new(*m);
    let b = sys.b_of(c);
    let mut s = LhState::from_cart(c);
    let e0 = energy::energy(&c.r, &c.v, m, T::zero());
    for _ in 0..steps {
        let (n, _) = direct_rk4(&sys, &s, b, LhTime::LogH, T::lit(h));
        s = n;
        if !s.is_finite() {
            return f64::INFINITY;
        }
    }
    let e1 = energy::energy(&s.r, &s.v, m, T::zero());
    ((e1 - e0) / e0).to_f64().unwrap().abs()
}

fn centred<T: Real>(ic: &Ic<T>) -> Cart<T> {
    let m = ic.m;
    let mt = m[0] + m[1] + m[2];
    let mut c = ic.s;
    let rc = (c.r[0] * m[0] + c.r[1] * m[1] + c.r[2] * m[2]) / mt;
    let vc = (c.v[0] * m[0] + c.v[1] * m[1] + c.v[2] * m[2]) / mt;
    for i in 0..3 {
        c.r[i] = c.r[i] - rc;
        c.v[i] = c.v[i] - vc;
    }
    c
}

fn main() {
    let res: usize = arg(1, 64);
    let steps: usize = arg(2, 4000);
    let h: f64 = arg(3, 1e-3);
    let cfg = EnsembleCfg::production();

    println!("STAGE 2: chain vs direct coordinates on `far`, {res}^2, {steps} fixed steps of");
    println!("fictitious size {h:.1e}, LogH transformation, RK4, both precisions.\n");
    println!(
        "  **PREDICTION (recorded before Stage 0 ran): chain helps at f32 and is near-invisible\n  \
         at f64.** The f64 row is the control -- a large move there means it is not behaving as a\n  \
         round-off fix and the f32 row cannot be attributed to precision.\n"
    );
    println!(
        "  Controlled to ONE thing: same transformation, stepper, step size, step count and ICs.\n  \
         Only the coordinates differ. Chain ordering selected once and FROZEN.\n"
    );

    let sl = grid::region("far", res, res, 0.05).expect("far");
    let ics: Vec<Ic<f64>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, 0, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )[0]
        })
        .collect();

    println!(
        "  {:>5} {:>8} {:>11} {:>11} {:>11} {:>7}",
        "prec", "coords", "drift p50", "drift p90", "drift p99", "nonfin"
    );

    let mut row: Vec<(String, String, f64)> = Vec::new();
    for prec in ["f64", "f32"] {
        for coords in ["direct", "chain"] {
            let d: Vec<f64> = ics
                .par_iter()
                .map(|ic| {
                    let c64 = centred::<f64>(ic);
                    match (prec, coords) {
                        ("f64", "direct") => drift_direct::<f64>(&c64, &ic.m, steps, h),
                        ("f64", _) => drift_chain::<f64>(&c64, &ic.m, steps, h),
                        (_, "direct") => {
                            let ic32 = ic.cast::<f32>();
                            let c = centred::<f32>(&ic32);
                            drift_direct::<f32>(&c, &ic32.m, steps, h)
                        }
                        _ => {
                            let ic32 = ic.cast::<f32>();
                            let c = centred::<f32>(&ic32);
                            drift_chain::<f32>(&c, &ic32.m, steps, h)
                        }
                    }
                })
                .collect();
            let nonfin = d.iter().filter(|x| !x.is_finite()).count();
            println!(
                "  {prec:>5} {coords:>8} {:>11.3e} {:>11.3e} {:>11.3e} {nonfin:>7}",
                q(&d, 0.50), q(&d, 0.90), q(&d, 0.99)
            );
            row.push((prec.into(), coords.into(), q(&d, 0.50)));
        }
    }

    println!("\n  CHAIN GAIN = log10(direct / chain), positive means chain is better");
    for prec in ["f64", "f32"] {
        let dir = row.iter().find(|r| r.0 == prec && r.1 == "direct").unwrap().2;
        let ch = row.iter().find(|r| r.0 == prec && r.1 == "chain").unwrap().2;
        println!("    {prec}: {:+.3} decades   (direct {dir:.3e}, chain {ch:.3e})", (dir / ch).log10());
    }

    println!(
        "\nHOW TO READ THIS\n\n\
         **The f64 row is the control.** Chain is a round-off fix; if it moves f64 by more than a\n\
         fraction of a decade it is doing something other than preserving digits, and the f32\n\
         number cannot then be attributed to precision.\n\n\
         **A null at f32 is a real answer.** It would say the `far` inversion measured in Stage 0\n\
         is not a separation-differencing effect, and that the mechanism behind it is still\n\
         unidentified -- which is where this investigation already stands, since the conditioning\n\
         story proposed for `far` was refuted by its own prediction."
    );
}
