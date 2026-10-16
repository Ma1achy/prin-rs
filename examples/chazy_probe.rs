//! Chazy final-motion classification "as of time t": a design probe.
//!
//! Classifiers run on the same snapshots of the same trajectories:
//!
//! * `chazy`   — the proposal: per-step, from the live state only (plus two running latches,
//!               `I_max` and the previous shape). Hierarchy ratio, outward motion, shape settling
//!               (the escape rule's closure chord over one window), kinetic/potential ratios with a
//!               parabolic band, the Lagrange–Jacobi certificate for E > 0, and I_max/I for OS
//!               candidates. Optionally (`ma`) the Mardling–Aarseth stability test for a bound
//!               hierarchy, borrowed from the gold standard.
//! * `fewbody` — the literature gold standard for deciding that a few-body scattering is over:
//!               Fregeau et al. 2004 (MNRAS 352, 1), §3.3. Hierarchy by repeatedly pairing the
//!               bound top-level nodes of smallest semimajor axis; resolved when every pair of
//!               top-level nodes recedes, the tidal perturbation on every binary is below `delta`,
//!               every tree is Mardling–Aarseth (2001) stable, and the top-level system has
//!               positive energy.
//! * `energy`  — naive snapshot: tightest pair is the binary, the sign of the third body's
//!               two-body energy decides.
//! * `dist`    — naive snapshot: farthest body beyond 10 R0 has escaped.
//!
//! Truth is the same trajectory marched to `T_LONG` and classified by Chazy's own definition:
//! the power-law growth exponent of each pair separation over the second half of the run
//! (1 hyperbolic, 2/3 parabolic, 0 bounded), at hierarchy ratio > 100. A trajectory with no
//! runaway at `T_LONG` is `unresolved`, and is `B-stable` if Mardling–Aarseth stable there.
//!
//! Run: cargo run --release --example chazy_probe -- [n_ics] [t_eval] [eps] [os_ratio] [unused] [unused] [t_truth] [seed0]

use prin_rs::chazy::{self, ChazyCfg, Gate, Latch, Motion as Cls};
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::physics::{decoder, energy, shape, Cart, Ic, PAIRS};
use rayon::prelude::*;
use std::collections::BTreeMap;

const SNAP: f64 = 0.4; // snapshot spacing = the escape rule's closure window (R-29)
const HORIZONS: [f64; 8] = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0];

use chazy::{kepler, ma_stable, other, out_kup, pair_kup};

fn pair_index(i: usize, j: usize) -> usize {
    PAIRS.iter().position(|&(a, b)| (a, b) == (i.min(j), i.max(j))).unwrap()
}

fn energy_naive(s: &Cart<f64>, m: &[f64; 3]) -> Cls {
    let d: Vec<f64> = PAIRS.iter().map(|&(i, j)| (s.r[j] - s.r[i]).norm()).collect();
    let tight = (0..3).min_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap()).unwrap();
    let (i, j) = PAIRS[tight];
    let b = 3 - i - j;
    let (kp, up, _, _) = pair_kup(s, m, i, j);
    let (ko, uo, _, _) = out_kup(s, m, b);
    match (kp < up, ko > uo) {
        (true, true) => Cls::HE,
        (false, true) => Cls::H,
        _ => Cls::BSoFar,
    }
}

fn dist_naive(s: &Cart<f64>, m: &[f64; 3], r0: f64) -> Cls {
    let (b, d) = (0..3).map(|b| (b, out_kup(s, m, b).2)).fold((0, 0.0), |a, x| if x.1 > a.1 { x } else { a });
    if d > 10.0 * r0 {
        let (i, j) = other(b);
        if (s.r[j] - s.r[i]).norm() > 10.0 * r0 { Cls::H } else { Cls::HE }
    } else {
        Cls::BSoFar
    }
}

/// Chazy's definition, applied at T_LONG: growth exponents of the pair separations.
fn truth(series: &[(f64, [f64; 3], f64, usize)], s: &Cart<f64>, m: &[f64; 3]) -> Cls {
    let n = series.len();
    let (t2, d2, hr, b) = series[n - 1];
    let (t1, d1, _, _) = *series.iter().find(|x| x.0 >= 0.5 * t2).unwrap();
    let beta: Vec<f64> = (0..3).map(|k| (d2[k] / d1[k]).ln() / (t2 / t1).ln()).collect();
    let dmin = d2.iter().cloned().fold(f64::INFINITY, f64::min);
    if dmin > 50.0 && beta.iter().all(|&x| x > 0.5) {
        let npar = beta.iter().filter(|&&x| x < 0.85).count();
        return match npar { 0 => Cls::H, 1 | 2 => Cls::HP, _ => Cls::P };
    }
    if hr > 100.0 {
        // The escaper's growth exponent from its two separations. The pair's own separation is
        // NOT a usable exponent: an eccentric binary's separation at two arbitrary phases differs
        // by up to (1+e)/(1-e), which reads as growth. At hierarchy > 100 the tidal perturbation
        // is < 2e-6 x mass ratio, so the pair's Kepler energy is the rigorous binding test.
        let (i, j) = other(b);
        let k = pair_index(i, j);
        let g = (0..3).filter(|&q| q != k).map(|q| beta[q]).fold(f64::INFINITY, f64::min);
        let (a_pair, _, _) = kepler(s.r[j] - s.r[i], s.v[j] - s.v[i], m[i] + m[j]);
        if a_pair > 0.0 && g > 0.5 {
            return if g > 0.85 { Cls::HE } else { Cls::PE };
        }
    }
    // No runaway: bounded at T_LONG. Stable hierarchy?
    let d: Vec<f64> = PAIRS.iter().map(|&(i, j)| (s.r[j] - s.r[i]).norm()).collect();
    let tight = (0..3).min_by(|&a, &c| d[a].partial_cmp(&d[c]).unwrap()).unwrap();
    let (i, j) = PAIRS[tight];
    if ma_stable(s, m, i, j, 3 - i - j) == Some(true) { Cls::BStable } else { Cls::Unresolved }
}

/// One classifier configuration. All of them read the same snapshots of the same trajectory.
#[derive(Clone, Copy)]
enum Method {
    Chazy(ChazyCfg),
    Fewbody(f64),
    Energy,
    Dist,
}

fn methods(base: &ChazyCfg) -> Vec<(String, Method)> {
    let mut v = Vec::new();
    let ch = |g: Gate, h_far: f64, tau: f64, delta: f64| Method::Chazy(ChazyCfg { gate: g, h_far, tau, delta, ..*base });
    v.push(("chazy h5 tau.01".to_string(), ch(Gate::Hierarchy, 5.0, 1e-2, base.delta)));
    v.push(("chazy h5 tau.05".to_string(), ch(Gate::Hierarchy, 5.0, 5e-2, base.delta)));
    v.push(("chazy h5 no-settle".to_string(), ch(Gate::Hierarchy, 5.0, f64::INFINITY, base.delta)));
    v.push(("chazy h3 no-settle".to_string(), ch(Gate::Hierarchy, 3.0, f64::INFINITY, base.delta)));
    v.push(("chazy h10 tau.01".to_string(), ch(Gate::Hierarchy, 10.0, 1e-2, base.delta)));
    v.push(("chazy-tidal d1e-5".to_string(), ch(Gate::Tidal, base.h_far, base.tau, 1e-5)));
    v.push(("chazy-tidal d1e-3".to_string(), ch(Gate::Tidal, base.h_far, base.tau, 1e-3)));
    v.push(("fewbody d1e-5".to_string(), Method::Fewbody(1e-5)));
    v.push(("fewbody d1e-3".to_string(), Method::Fewbody(1e-3)));
    v.push(("naive energy".to_string(), Method::Energy));
    v.push(("naive distance".to_string(), Method::Dist));
    v
}

struct Run {
    e: f64,
    truth: Cls,
    at: Vec<Vec<Cls>>,                // [horizon][method]
    first_commit: Vec<Option<f64>>,  // first time the method commits
    first_correct: Vec<Option<f64>>,
    flips: Vec<usize>,               // coarse-label changes after first commit
    fail: bool,
    t_truth_end: f64,
}

fn march(ic: &Ic<f64>, t_eval: f64, t_truth: f64, ms: &[(String, Method)]) -> Run {
    let m = ic.m;
    let mut s = ic.s;
    let r0 = energy::hyperradius(&s.r, &m);
    let e = energy::energy(&s.r, &s.v, &m, 0.0);
    let opts = HgOpts::<f64> { stop_on_event: false, r_coll_frac: 0.0, ..Default::default() };
    let mut h = Latch::new(energy::inertia(&s.r, &m));
    let mut hz = 0usize;
    let mut at = Vec::new();
    let mut series = Vec::new();
    let mut cls_hist: Vec<Vec<Cls>> = Vec::new();
    let mut fail = false;
    let record = |s: &Cart<f64>, t: f64, series: &mut Vec<(f64, [f64; 3], f64, usize)>| {
        let d: [f64; 3] = [0, 1, 2].map(|k| { let (i, j) = PAIRS[k]; (s.r[j] - s.r[i]).norm() });
        let (b, r) = (0..3)
            .map(|b| { let (i, j) = other(b); (b, out_kup(s, &m, b).2 / d[pair_index(i, j)]) })
            .fold((0, 0.0), |a, x| if x.1 > a.1 { x } else { a });
        series.push((t, d, r, b));
        (d.iter().cloned().fold(f64::INFINITY, f64::min), r)
    };
    // Phase 1: every classifier at every snapshot, to t_eval.
    let n = (t_eval / SNAP).round() as usize;
    let mut t = 0.0;
    for step in 1..=n {
        let o = integrate_hg(s, &m, SNAP, 2, 1e-2, 2_000_000, &opts);
        if !o.finite || o.budget_exhausted {
            fail = true;
            break;
        }
        s = o.state;
        t = step as f64 * SNAP;
        let i_now = energy::inertia(&s.r, &m);
        h.i_max = h.i_max.max(i_now);
        let nv = shape::shape_vec(&s.r, &m);
        let c: Vec<Cls> = ms.iter().map(|(_, mm)| match mm {
            Method::Chazy(pp) => chazy::classify(&s, &m, e, r0, i_now, &nv, &h, pp),
            Method::Fewbody(d) => chazy::fewbody_resolved(&s, &m, e, *d),
            Method::Energy => energy_naive(&s, &m),
            Method::Dist => dist_naive(&s, &m, r0),
        }).collect();
        cls_hist.push(c.clone());
        h.n_prev = Some(nv);
        while hz < HORIZONS.len() && t >= HORIZONS[hz] - 1e-9 {
            at.push(c.clone());
            hz += 1;
        }
        if step % 5 == 0 { record(&s, t, &mut series); }
    }
    // Phase 2: the truth run continues alone, in coarse chunks, until the motion is
    // unambiguously asymptotic (hierarchy ratio or minimum separation over 1e3) or t_truth.
    let chunk = 10.0;
    while !fail && t < t_truth - 1e-9 {
        let o = integrate_hg(s, &m, chunk, 4, 1e-2, 20_000_000, &opts);
        if !o.finite || o.budget_exhausted {
            fail = true;
            break;
        }
        s = o.state;
        t += chunk;
        let (dmin, hr) = record(&s, t, &mut series);
        if (hr > 1e3 || dmin > 1e3 * r0) && t > 2.0 * t_eval.min(200.0) { break; }
    }
    let truth = if fail || series.len() < 4 { Cls::Unresolved } else { truth(&series, &s, &m) };
    let mut first_commit = vec![None; ms.len()];
    let mut first_correct = vec![None; ms.len()];
    let mut flips = vec![0usize; ms.len()];
    for (q, c) in cls_hist.iter().enumerate() {
        let t = (q + 1) as f64 * SNAP;
        for k in 0..ms.len() {
            if c[k].committed() {
                if first_commit[k].is_none() { first_commit[k] = Some(t); }
                if first_correct[k].is_none() && c[k].coarse() == truth.coarse() { first_correct[k] = Some(t); }
            }
        }
    }
    for k in 0..ms.len() {
        let mut last: Option<Cls> = None;
        for c in &cls_hist {
            if c[k].committed() {
                if let Some(l) = last { if l.coarse() != c[k].coarse() { flips[k] += 1; } }
                last = Some(c[k]);
            }
        }
    }
    Run { e, truth, at, first_commit, first_correct, flips, fail, t_truth_end: t }
}

/// A tiny xorshift, so the probe needs no new dependency.
fn rng(seed: u64) -> impl FnMut() -> f64 {
    let mut x = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
    move || {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
}
fn normal(r: &mut impl FnMut() -> f64) -> f64 {
    let (u, v) = (r().max(1e-300), r());
    (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
}

fn pct(a: usize, b: usize) -> f64 { if b == 0 { f64::NAN } else { 100.0 * a as f64 / b as f64 } }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let arg = |i: usize, d: f64| a.get(i).and_then(|x| x.parse().ok()).unwrap_or(d);
    let n_ics = arg(1, 300.0) as usize;
    let t_long = arg(2, 1000.0); // evaluation horizon
    let base = ChazyCfg { eps: arg(3, 0.05), os_ratio: arg(4, 25.0), ..ChazyCfg::default() };
    let t_truth = arg(7, 10000.0);
    let seed0 = arg(8, 1000.0) as u64;
    let ms = methods(&base);

    // Mix: half free fall (E < 0), a quarter small momenta (mostly E < 0), a quarter O(1) momenta (mostly E > 0).
    let ics: Vec<Ic<f64>> = (0..n_ics)
        .map(|k| {
            let mut r = rng(seed0 + k as u64);
            let mut z = decoder::Latent::default();
            z.z_alpha = 1.5 * normal(&mut r);
            z.z_beta = 1.5 * normal(&mut r);
            z.z_mu = [normal(&mut r), normal(&mut r)];
            let sig = [0.0, 0.0, 0.3, 1.0][k % 4];
            for q in 0..4 { z.z_q[q] = sig * normal(&mut r); }
            decoder::decode(&z).ic
        })
        .collect();
    let t0 = std::time::Instant::now();
    let runs: Vec<Run> = ics.par_iter().map(|ic| march(ic, t_long, t_truth, &ms)).collect();
    let ok: Vec<&Run> = runs.iter().filter(|r| !r.fail).collect();
    println!("# chazy_probe n_ics={n_ics} t_eval={t_long} t_truth<={t_truth} eps={} os_ratio={} seed0={seed0} snap={SNAP}",
             base.eps, base.os_ratio);
    println!("# integrator heggie eta=1e-2, r_coll=0 (regularised through binary collisions), stop_on_event=false");
    println!("# {} finished, {} failed (non-finite or budget), {:.1}s", ok.len(), runs.len() - ok.len(),
             t0.elapsed().as_secs_f64());

    let mut tc: BTreeMap<(&str, &str), usize> = BTreeMap::new();
    for r in &ok { *tc.entry((r.truth.name(), if r.e > 0.0 { "E>0" } else { "E<0" })).or_insert(0) += 1; }
    println!("\n## truth: Chazy's definition (separation growth exponents), run continued until asymptotic or t = {t_truth}");
    for ((c, s), n) in &tc { println!("  {c:<11} {s}  {n}"); }
    let mut te: Vec<f64> = ok.iter().map(|r| r.t_truth_end).collect();
    te.sort_by(|x, y| x.partial_cmp(y).unwrap());
    println!("  truth run end time p50 {:.0} p90 {:.0} max {:.0}", te[te.len() / 2], te[te.len() * 9 / 10], te[te.len() - 1]);

    let hs: Vec<(usize, f64)> = HORIZONS.iter().cloned().enumerate().filter(|x| x.1 <= t_long).collect();
    let resolved: Vec<&&Run> = ok.iter().filter(|r| r.truth != Cls::Unresolved).collect();
    for (label, sel) in [("E<0", -1i8), ("E>0", 1), ("all", 0)] {
        let sub: Vec<&&&Run> = resolved.iter().filter(|r| sel == 0 || (sel < 0) == (r.e < 0.0)).collect();
        println!("\n## [{label}] correct% / wrong% at horizon T, over the {} trajectories with a resolved truth", sub.len());
        print!("{:<20}", "method");
        for (_, th) in &hs { print!(" {:>13}", format!("T={th}")); }
        println!();
        for (k, (name, _)) in ms.iter().enumerate() {
            print!("{name:<20}");
            for &(hi, _) in &hs {
                let (mut ok_, mut w) = (0, 0);
                for r in &sub {
                    let x = r.at[hi][k];
                    if x.committed() { if x.coarse() == r.truth.coarse() { ok_ += 1 } else { w += 1 } }
                }
                print!(" {:>6.1}/{:<6.2}", pct(ok_, sub.len()), pct(w, sub.len()));
            }
            println!();
        }
    }

    let unres: Vec<&&Run> = ok.iter().filter(|r| r.truth == Cls::Unresolved).collect();
    println!("\n## commits on the {} trajectories with NO resolved truth (count; not necessarily wrong)", unres.len());
    print!("{:<20}", "method");
    for (_, th) in &hs { print!(" {:>7}", format!("T={th}")); }
    println!();
    for (k, (name, _)) in ms.iter().enumerate() {
        print!("{name:<20}");
        for &(hi, _) in &hs { print!(" {:>7}", unres.iter().filter(|r| r.at[hi][k].committed()).count()); }
        println!();
    }

    println!("\n## [E<0, resolved] time to first correct commit (p50, p90), never-correct count, flips after commit");
    let neg: Vec<&&&Run> = resolved.iter().filter(|r| r.e < 0.0).collect();
    for (k, (name, _)) in ms.iter().enumerate() {
        let mut fr: Vec<f64> = neg.iter().filter_map(|r| r.first_correct[k]).collect();
        fr.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let q = |f: f64| if fr.is_empty() { f64::NAN } else { fr[((fr.len() - 1) as f64 * f) as usize] };
        let flips: usize = ok.iter().map(|r| r.flips[k]).sum();
        let flipped = ok.iter().filter(|r| r.flips[k] > 0).count();
        let premature = neg.iter().filter(|r| match (r.first_commit[k], r.first_correct[k]) {
            (Some(a), Some(b)) => a < b,
            (Some(_), None) => true,
            _ => false,
        }).count();
        println!("  {name:<20} p50 {:>7.1}  p90 {:>7.1}  never {:>4}  wrong-first {:>3}  flips {:>5} on {:>4} trajectories",
                 q(0.5), q(0.9), neg.len() - fr.len(), premature, flips, flipped);
    }

    // Labels at the last horizon against truth, for the proposal and the gold standard.
    let hl = hs.last().unwrap().0;
    for k in [0usize, 2, 5, 7] {
        let mut cm: BTreeMap<(&str, &str), usize> = BTreeMap::new();
        for r in &ok { *cm.entry((r.at[hl][k].name(), r.truth.name())).or_insert(0) += 1; }
        println!("\n## {} at T = {} (row) -> truth", ms[k].0, HORIZONS[hl]);
        for ((x, y), n) in &cm { println!("  {x:<11} -> {y:<11} {n}"); }
    }

    println!("\n## class mix over time under '{}', % of all finished trajectories", ms[2].0);
    println!("   (bounded-side labels are approximations that hold up to T; see the truth census above)");
    let names = [Cls::H, Cls::HP, Cls::P, Cls::HE, Cls::PE, Cls::BStable, Cls::Transient, Cls::OsCandidate, Cls::Excursion, Cls::BSoFar];
    print!("{:>6}", "T");
    for c in names { print!(" {:>9}", c.name()); }
    println!();
    for &(hi, th) in &hs {
        print!("{:>6}", th);
        for c in names { print!(" {:>9.1}", pct(ok.iter().filter(|r| r.at[hi][2] == c).count(), ok.len())); }
        println!();
    }
    let e_pos_unres = ok.iter().filter(|r| r.e > 0.0 && r.truth == Cls::Unresolved).count();
    println!("\nE > 0 with no resolved truth (Lagrange–Jacobi says 0 in the limit): {e_pos_unres}");
}
