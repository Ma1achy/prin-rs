//! **Two questions left over from the 2x2, and each one has a way to come out the wrong way.**
//!
//! # §1 Does the cheap configuration keep the cadence insensitivity?
//!
//! `HG eta x8` costs 0.98x AZ's work for 21x less drift, which looks like a free performance win.
//! But the reason to want Heggie is **not** drift, it is that its field barely moves when the sync
//! cadence changes — chord `4.83e-2` against AZ's `5.16e-1`. A coarse run might be spending
//! exactly that.
//!
//! So the controlled cadence pair is re-run at a coarse base. The protocol is the standing one:
//! double `n_sync` **and** `eta` together so the step size is held fixed and "more boundaries" is
//! not "finer stepping" in disguise. `steps p50` is the check that it worked.
//!
//! **Compared within the pair, never against the fine baseline.** A coarse run differs from a fine
//! one for the ordinary reason that it is coarse; the question is whether *doubling the cadence*
//! moves it, and that is a chord between two coarse rows.
//!
//! Run at 256^2 deliberately. **No chord ratio may be quoted from a coarse grid** — this project
//! has one overstated 26x and another understated 8x by doing so — and §1 is nothing but a chord
//! ratio.
//!
//! # §2 Is Heggie's 19x the unregularised third side?
//!
//! With the step limit off, Heggie's drift is 19x AZ's, and the standing explanation is that AZ's
//! `Gamma` still carries `-A B m_b m_c/|R3|` where `Gamma*` has no `1/r` term anywhere. **That is
//! a story until it is split.**
//!
//! `AzOut` already carries the discriminator: `d_min_ref` is the min over the two **regularised**
//! pairs and `d_min_true` the min over all three, so `d_min_true < d_min_ref` says the closest
//! approach of the whole run landed on the side AZ does not regularise. The prediction is that
//! AZ's disadvantage is **concentrated** there.
//!
//! **The null is informative.** If the ratio is flat across the split, the 19x is not the third
//! side and the mechanism on record is wrong. And the split is on AZ's own geometry, which is the
//! only honest way to define it — Heggie has no third side to condition on.
//!
//! The limit is **off** in both integrators here, because §2 is about the regularisation and the
//! 2x2 showed the limit contributes an independent 29x that would otherwise be mixed in.
//!
//! # The two sections need OPPOSITE settings for the step limit
//!
//! §1 is about the **shipping** configuration, so the limit is **on** there: the 4.83e-2 the
//! performance win would be spending was measured with it on, and a chord measured without it is
//! a chord for a configuration nobody runs. §2 is about the **regularisation**, so the limit is
//! **off** there: the 2x2 showed it contributes an independent 29x that would otherwise be mixed
//! into the answer.
//!
//! The first cut of this harness ran one setting for both and §1 read `2.56e-1` where the shipping
//! number is `4.83e-2` — a fivefold error, from one flag serving two sections with opposite
//! requirements. That is the `prin --size` shape at a new site, and the fix is the same: carry the
//! setting per arm and print it.
//!
//! Args: `res root max_steps`. Resumable; re-run until it prints.

use rayon::prelude::*;

use prin_rs::ensemble::jitter;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::output::ckpt::Ckpt;
use prin_rs::physics::Ic;

const WINDOW: f64 = 0.4;
const T_MAX: f64 = 50.0;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// Per pixel: max drift over copies, summed steps, and the two `d_min` measures.
///
/// `d_min_ref` / `d_min_true` are AZ's and are `NaN` for Heggie, which has no unregularised side
/// and therefore no such distinction to report. Carrying them as `NaN` rather than as zero keeps
/// that an absence rather than a value.
#[derive(Clone, Copy)]
struct Rec {
    drift: f64,
    steps: f64,
    d_ref: f64,
    d_true: f64,
}

#[derive(Clone, Copy, PartialEq)]
enum Arm {
    Az,
    Hg,
}

fn pixel(cs: &[Ic<f64>], arm: Arm, n_sync: usize, eta: f64, max_steps: usize, limit: bool) -> Rec {
    let f = if limit { 0.02 } else { 0.0 };
    let (mut drift, mut steps) = (0.0f64, 0.0f64);
    let (mut d_ref, mut d_true) = (f64::INFINITY, f64::INFINITY);
    for c in cs {
        match arm {
            Arm::Az => {
                let o = az::integrate_az_opts(
                    c.s, &c.m, T_MAX, n_sync, eta, max_steps,
                    &AzOpts {
                        stop_on_event: false,
                        r_coll_frac: 0.0,
                        step_limit: if limit {
                            az::StepLimit::Predictive
                        } else {
                            az::StepLimit::None
                        },
                        step_limit_f: f,
                        ..Default::default()
                    },
                );
                steps += o.steps as f64;
                drift = drift.max(if o.finite { o.drift } else { f64::INFINITY });
                d_ref = d_ref.min(o.d_min_ref);
                d_true = d_true.min(o.d_min_true);
            }
            Arm::Hg => {
                let o = integrate_hg(
                    c.s, &c.m, T_MAX, n_sync, eta, max_steps,
                    &HgOpts { step_limit_f: f, ..Default::default() },
                );
                steps += o.steps as f64;
                drift = drift.max(if o.finite { o.drift } else { f64::INFINITY });
                d_true = d_true.min(o.d_min);
            }
        }
    }
    Rec { drift, steps, d_ref, d_true }
}

fn encode(v: &[Rec]) -> Vec<u8> {
    let mut b = Vec::with_capacity(v.len() * 32);
    for r in v {
        for x in [r.drift, r.steps, r.d_ref, r.d_true] {
            b.extend_from_slice(&x.to_le_bytes());
        }
    }
    b
}

fn decode(b: &[u8]) -> Vec<Rec> {
    b.chunks_exact(32)
        .map(|c| {
            let f = |i: usize| f64::from_le_bytes(c[i * 8..i * 8 + 8].try_into().unwrap());
            Rec { drift: f(0), steps: f(1), d_ref: f(2), d_true: f(3) }
        })
        .collect()
}

fn lg(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 { x.log10() } else { 2.0 }
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let max_steps: usize = arg(3, 400_000);
    let dir = format!("{root}/step_control/heggie_attribution");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let n0 = (T_MAX / WINDOW).round() as usize;
    let cfg = EnsembleCfg::production();
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let e0 = cfg.eta;

    println!("config_stability {res}^2, termination OFF. Step limit per arm -- see §1/§2 above.");
    println!("t_max={T_MAX} max_steps={max_steps} copies={}\n", cfg.n_extra + 1);

    let copies: Vec<Vec<Ic<f64>>> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            jitter::copies_with_path::<f64>(
                &sl, k, cfg.n_extra, cfg.jitter_frac, cfg.seed, cfg.jitter_scheme, cfg.decode_path,
            )
        })
        .collect();

    // (label, arm, n_sync, eta, step limit)
    let arms: Vec<(&str, Arm, usize, f64, bool)> = vec![
        // §1 -- the SHIPPING configuration, limit on.
        ("HG fine   n=125 eta x1  limit ON ", Arm::Hg, n0, e0, true),
        ("HG fine   n=250 eta x2  limit ON ", Arm::Hg, n0 * 2, e0 * 2.0, true),
        ("HG coarse n=125 eta x4  limit ON ", Arm::Hg, n0, e0 * 4.0, true),
        ("HG coarse n=250 eta x8  limit ON ", Arm::Hg, n0 * 2, e0 * 8.0, true),
        // §2 -- the REGULARISATION alone, limit off.
        ("AZ        n=125 eta x1  limit OFF", Arm::Az, n0, e0, false),
        ("HG        n=125 eta x1  limit OFF", Arm::Hg, n0, e0, false),
    ];

    let key = format!(
        "heggie_attribution v1 res={res} t_max={T_MAX} eta={e0} max_steps={max_steps} n0={n0}\n{}",
        cfg.provenance()
    );
    let (mut ck, have) = Ckpt::open(&format!("{dir}/arms.ckpt"), &key).expect("checkpoint");
    if !have.is_empty() {
        println!("  resuming: {} of {} arms already computed\n", have.len(), arms.len());
    }

    let mut out: Vec<Vec<Rec>> = Vec::new();
    for (i, (label, arm, ns, eta, lim)) in arms.iter().enumerate() {
        let v = match have.get(&(i as u64)) {
            Some(b) => decode(b),
            None => {
                let v: Vec<Rec> =
                    copies.par_iter().map(|cs| pixel(cs, *arm, *ns, *eta, max_steps, *lim)).collect();
                ck.put(i as u64, &encode(&v)).expect("checkpoint write");
                v
            }
        };
        let mut st: Vec<f64> = v.iter().map(|r| r.steps).collect();
        let mut dr: Vec<f64> = v.iter().map(|r| r.drift).filter(|x| x.is_finite()).collect();
        println!(
            "  {label}   steps p50 {:>10.3e}   drift p50 {:>10.3e}   nonfin {:>5}",
            q(&mut st, 0.5),
            q(&mut dr, 0.5),
            v.len() - dr.len()
        );
        out.push(v);
    }

    // ---------------------------------------------------------------------------------------
    println!("\n== §1  DOES THE CHEAP CONFIGURATION KEEP THE CADENCE INSENSITIVITY? ==");
    println!("  Predictive step limit ON, the shipping configuration. Chord is |delta log10 drift|
  WITHIN each pair. Doubling n_sync and eta together");
    println!("  holds the step size fixed, so `steps p50` must be roughly flat inside a pair.\n");
    println!("  {:>18} {:>12} {:>12} {:>12}", "pair", "steps ratio", "chord p50", "chord p90");
    for (name, a, b) in [("fine  x1 -> x2", 0usize, 1usize), ("coarse x4 -> x8", 2, 3)] {
        let mut sa: Vec<f64> = out[a].iter().map(|r| r.steps).collect();
        let mut sb: Vec<f64> = out[b].iter().map(|r| r.steps).collect();
        let mut c: Vec<f64> = (0..out[a].len())
            .map(|i| (lg(out[b][i].drift) - lg(out[a][i].drift)).abs())
            .filter(|x| x.is_finite())
            .collect();
        println!(
            "  {name:>18} {:>12.3} {:>12.3e} {:>12.3e}",
            q(&mut sb, 0.5) / q(&mut sa, 0.5),
            q(&mut c.clone(), 0.5),
            q(&mut c, 0.9)
        );
    }
    println!("\n  **If the coarse chord is comparable to the fine one, the performance win is");
    println!("  real and costs nothing that matters. If it is an order larger, the cheap");
    println!("  configuration is spending the exact property the port exists for.**");

    // ---------------------------------------------------------------------------------------
    println!("\n== §2  IS HEGGIE'S ADVANTAGE THE UNREGULARISED THIRD SIDE? ==");
    println!("  Split on AZ's OWN geometry: `d_min_true < d_min_ref` means the closest approach");
    println!("  of the run landed on the pair AZ does not regularise. Heggie has no third side,");
    println!("  so there is no Heggie-side split to make and this is the only honest definition.\n");

    let (az, hg) = (&out[4], &out[5]);
    let n = az.len();
    // A strict `<` would count float noise as a hit. The unregularised side has to be closer by a
    // clear margin for the trajectory to have spent time there.
    let unreg: Vec<bool> = (0..n).map(|i| az[i].d_true < az[i].d_ref * 0.999).collect();
    let n_un = unreg.iter().filter(|x| **x).count();
    println!("  pixels whose closest approach was on the UNREGULARISED side: {n_un} of {n}");
    if n_un == 0 || n_un == n {
        println!("  **SATURATED -- the split has one class and cannot discriminate.**");
    }
    println!(
        "\n  {:>26} {:>8} {:>12} {:>12} {:>10}",
        "population", "n", "AZ drift p50", "HG drift p50", "AZ/HG"
    );
    for (name, want) in [("closest on UNREG side", true), ("closest on a reg pair", false)] {
        let idx: Vec<usize> = (0..n).filter(|&i| unreg[i] == want).collect();
        let mut a: Vec<f64> = idx.iter().map(|&i| az[i].drift).filter(|x| x.is_finite()).collect();
        let mut h: Vec<f64> = idx.iter().map(|&i| hg[i].drift).filter(|x| x.is_finite()).collect();
        let (qa, qh) = (q(&mut a, 0.5), q(&mut h, 0.5));
        println!("  {name:>26} {:>8} {qa:>12.3e} {qh:>12.3e} {:>10.1}", idx.len(), qa / qh);
    }
    println!("\n  Also, non-finite counts by population, because a divergence is the extreme of");
    println!("  the same failure and a median cannot see it.\n");
    println!("  {:>26} {:>8} {:>12} {:>12}", "population", "n", "AZ nonfin", "HG nonfin");
    for (name, want) in [("closest on UNREG side", true), ("closest on a reg pair", false)] {
        let idx: Vec<usize> = (0..n).filter(|&i| unreg[i] == want).collect();
        let na = idx.iter().filter(|&&i| !az[i].drift.is_finite()).count();
        let nh = idx.iter().filter(|&&i| !hg[i].drift.is_finite()).count();
        println!("  {name:>26} {:>8} {na:>12} {nh:>12}", idx.len());
    }
    println!("\n  **A ratio that is flat across the split refutes the mechanism on record.**");
}
