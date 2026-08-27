//! The `dtau` step-control measurement: what the fix costs, what it buys, and whether the
//! damage it removes is the damage seen in the image.
//!
//! # The mechanism
//!
//! `dt = A*B*dtau`, and `dtau` was sized **once** at each sync interval's entry. So the physical
//! step is `eta*dt_left` only while `A*B` stays near its entry value. A trajectory sitting at a
//! close encounter *at a sync boundary* has a tiny `A0*B0`, so `dtau` is enormous; as the bodies
//! separate through the interval `A*B` grows by orders and `dt` grows with it. Giant physical
//! steps immediately after an encounter.
//!
//! The population that suffers is **encounters coinciding with a boundary**, which is a thin set
//! — so the damage clusters spatially instead of tracking `d_min`, and `d_min` split by whether a
//! trajectory was affected is the control that says so.
//!
//! # The order matters
//!
//! §1 is the **gate**. The fix trades accuracy for step count, and if the budget is simply
//! exhausted elsewhere then non-finite pixels have been swapped for budget-exhausted ones — a
//! different failure wearing a different colour. It is reported as its own number and never
//! folded into the drift result.
//!
//! One nominal trajectory per pixel, no ensemble: the quantities here are per-trajectory
//! (`steps`, `ab_min`, `budget_exhausted`) and an ensemble reduction would hide them behind a
//! max. A grid is still needed for §3, which is a spatial statistic.

use rayon::prelude::*;

use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{integrate_az_opts, AzOpts, AzOut, DtauMode};
use prin_rs::physics::Cart;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const ETA: f64 = 1e-2;
const MAX_STEPS: usize = 30_000;
/// The drift above which a trajectory is called damaged. The user's NumPy measurement used the
/// same cut, so the reproduction in §2 is like-for-like.
const HOT: f64 = 1e-6;

const MODES: [(DtauMode, &str); 3] = [
    (DtauMode::FixedPerInterval, "fixed"),
    (DtauMode::PerStepRemaining, "per-step-rem"),
    (DtauMode::PerStepInterval, "per-step-int"),
];

struct Case {
    name: String,
    t_max: f64,
    n_sync: usize,
    r_coll: f64,
    n: usize,
    ics: Vec<(Cart<f64>, [f64; 3])>,
}

fn sample(chart: &Chart, body: usize, cx: f64, cy: f64, half: f64, n: usize)
    -> Vec<(Cart<f64>, [f64; 3])>
{
    let mut out = Vec::with_capacity(n * n);
    for j in 0..n {
        for i in 0..n {
            let u = cx - half + 2.0 * half * (i as f64 + 0.5) / n as f64;
            let v = cy - half + 2.0 * half * (j as f64 + 0.5) / n as f64;
            let ic = grid::decode_state(chart, body, u, v);
            out.push((ic.s, ic.m));
        }
    }
    out
}

fn run(c: &Case, mode: DtauMode) -> Vec<AzOut<f64>> {
    let o = AzOpts::<f64> {
        r_coll_frac: c.r_coll,
        // Nothing terminal. A run stopped early has a shorter step count and a smaller drift for
        // a reason that has nothing to do with the step control, and `d_min` over a truncated
        // run inherits the truncation.
        stop_on_event: false,
        stop_on_escape: false,
        dtau_mode: mode,
        ..Default::default()
    };
    c.ics
        .par_iter()
        .map(|(s, m)| integrate_az_opts(*s, m, c.t_max, c.n_sync, ETA, MAX_STEPS, &o))
        .collect()
}

fn q(v: &[f64], p: f64) -> f64 {
    let mut w: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    if w.is_empty() { f64::NAN } else { prin_rs::stats::quantile(&mut w, p) }
}

/// Fraction of hot pixels having a hot 4-neighbour, and that fraction over the base rate.
///
/// **A ratio of 1 is chance.** Under the defect the damaged set is the thin "encounter at a
/// boundary" set, which is spatially coherent, so the ratio runs above 1; if the fix works it
/// falls toward 1 as the survivors become scattered rather than clustered.
fn clustering(hot: &[bool], n: usize) -> (f64, f64, usize) {
    let nh = hot.iter().filter(|x| **x).count();
    if nh == 0 {
        return (f64::NAN, f64::NAN, 0);
    }
    let base = nh as f64 / hot.len() as f64;
    let mut with = 0usize;
    for j in 0..n {
        for i in 0..n {
            if !hot[j * n + i] {
                continue;
            }
            let mut any = false;
            for (di, dj) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                let (a, b) = (i as i64 + di, j as i64 + dj);
                if a >= 0 && b >= 0 && (a as usize) < n && (b as usize) < n && hot[b as usize * n + a as usize] {
                    any = true;
                }
            }
            with += any as usize;
        }
    }
    let f = with as f64 / nh as f64;
    // Chance: probability at least one of up to four independent neighbours is hot.
    (f, f / (1.0 - (1.0 - base).powi(4)), nh)
}

fn main() {
    let n: usize = arg(1, 48);

    let mut cases: Vec<Case> = Vec::new();
    let nsync13 = 33usize;
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            cases.push(Case { name: region.into(), t_max: 13.0, n_sync: nsync13, r_coll: 1e-3, n,
                              ics: sample(&Chart::BodyPlane, body, cx, cy, 0.05, n) });
        }
    }
    let (cs, sx, sy, sh) = Chart::config_stability();
    cases.push(Case { name: "config_stability".into(), t_max: 50.0, n_sync: 125, r_coll: 0.005, n,
                      ics: sample(&cs, 0, sx, sy, sh, n) });
    // The user's own NumPy settings, so §2 and §3 are a like-for-like reproduction rather than a
    // different measurement that happens to agree. Their numbers: 8 non-finite -> 2,
    // frac drift>1e-6 35.6% -> 18.6%, median on affected 4.32e-06 -> 2.96e-08, clustering 1.7x,
    // and 6 of 6 non-finite pixels with a high-drift neighbour.
    cases.push(Case { name: "config_stability@t6".into(), t_max: 6.0, n_sync: 12, r_coll: 0.005, n,
                      ics: sample(&cs, 0, sx, sy, sh, n) });

    println!("dtau_step -- {n}x{n} nominal trajectories per case, eta={ETA}, max_steps={MAX_STEPS}\n");
    println!("MODES: `fixed` is the shipped-until-now sizing and what EVERY committed Rust and");
    println!("NumPy number was taken under. `per-step-rem` is the obvious repair and is Zeno by");
    println!("arithmetic: dt ~ eta*rem gives rem_(n+1) = rem_n (1-eta), so the interval is");
    println!("approached geometrically and never completed. `per-step-int` holds dt_left fixed and");
    println!("recomputes only A*B, capped at the entry value.\n");

    let mut all: Vec<(String, DtauMode, Vec<AzOut<f64>>)> = Vec::new();
    for c in &cases {
        for (mode, _) in MODES {
            let t0 = std::time::Instant::now();
            let o = run(c, mode);
            eprintln!("{} {:?} {:.1}s", c.name, mode, t0.elapsed().as_secs_f64());
            all.push((c.name.clone(), mode, o));
        }
    }
    let get = |name: &str, mode: DtauMode| -> &Vec<AzOut<f64>> {
        &all.iter().find(|(n, m, _)| n == name && *m == mode).unwrap().2
    };

    // -------------------------------------------------------------------------------------
    println!("== 1. STEP-COUNT DISTRIBUTION -- THE GATE ==");
    println!("The fix trades accuracy for step count. If the budget is exhausted elsewhere, the");
    println!("non-finite pixels have been swapped for budget-exhausted ones -- a different failure");
    println!("wearing a different colour. `budget` is that count and it is not folded into §2.\n");
    println!("**`t/t_max` IS THE DISCRIMINATOR AND IT BELONGS HERE, NOT IN §2.** A mode that stalls");
    println!("has a beautiful drift -- it barely moved. Reading §2 without this column would score");
    println!("`per-step-rem` as five orders better than either real mode. A difference can be small");
    println!("because both sides are right or because one side is dead.\n");
    println!("{:<22} {:<14} {:>9} {:>9} {:>9} {:>8} {:>8} {:>8}",
             "case", "mode", "steps p50", "steps p99", "steps max", "budget", "frac", "t/t_max");
    for c in &cases {
        for (mode, mn) in MODES {
            let o = get(&c.name, mode);
            let st: Vec<f64> = o.iter().map(|x| x.steps as f64).collect();
            let be = o.iter().filter(|x| x.budget_exhausted).count();
            let reached: Vec<f64> = o.iter().map(|x| x.t / c.t_max).collect();
            println!("{:<22} {:<14} {:>9.0} {:>9.0} {:>9.0} {:>8} {:>8.4} {:>8.4}",
                     c.name, mn, q(&st, 0.5), q(&st, 0.99),
                     st.iter().cloned().fold(0.0, f64::max), be, be as f64 / o.len() as f64,
                     q(&reached, 0.5));
        }
        println!();
    }

    // -------------------------------------------------------------------------------------
    println!("== 2. DRIFT AND NON-FINITE ==");
    println!("`affected` is the union over modes of trajectories any mode drives above {HOT:e} --");
    println!("a FIXED population, so the medians compare the same trajectories rather than each");
    println!("mode's own tail. `d_min` is printed for affected and unaffected: if they agree, the");
    println!("damage is not simply `close encounters`, it is encounters coinciding with a");
    println!("boundary, and that thin set is why it clusters.\n");
    println!("{:<22} {:<14} {:>10} {:>10} {:>9} {:>7} {:>11} {:>11}",
             "case", "mode", "drift p50", "drift p99", "frac>hot", "nonfin", "med(affctd)", "max drift");
    for c in &cases {
        let np = c.ics.len();
        let mut affected = vec![false; np];
        for (mode, _) in MODES {
            for (k, o) in get(&c.name, mode).iter().enumerate() {
                if !o.drift.is_finite() || o.drift > HOT {
                    affected[k] = true;
                }
            }
        }
        for (mode, mn) in MODES {
            let o = get(&c.name, mode);
            let d: Vec<f64> = o.iter().map(|x| x.drift).collect();
            let nf = o.iter().filter(|x| !x.finite).count();
            let hot = d.iter().filter(|x| !x.is_finite() || **x > HOT).count();
            let aff: Vec<f64> = d.iter().enumerate().filter(|(k, _)| affected[*k])
                                 .map(|(_, v)| *v).collect();
            println!("{:<22} {:<14} {:>10.3e} {:>10.3e} {:>9.4} {:>7} {:>11.3e} {:>11.3e}",
                     c.name, mn, q(&d, 0.5), q(&d, 0.99), hot as f64 / np as f64, nf,
                     q(&aff, 0.5), d.iter().cloned().filter(|x| x.is_finite()).fold(0.0, f64::max));
        }
        let dm = |sel: bool| -> f64 {
            let v: Vec<f64> = get(&c.name, DtauMode::PerStepInterval).iter().enumerate()
                .filter(|(k, _)| affected[*k] == sel).map(|(_, o)| o.d_min_true).collect();
            q(&v, 0.5)
        };
        println!("{:<22} d_min p50  affected {:.3e}   unaffected {:.3e}   (n_affected {})",
                 "", dm(true), dm(false), affected.iter().filter(|x| **x).count());
        println!();
    }

    // -------------------------------------------------------------------------------------
    println!("== 3. THE SPATIAL TEST -- what ties the numbers to the image ==");
    println!("`clust` is the fraction of hot pixels with a hot 4-neighbour divided by the same");
    println!("under a random field of the same density. **1.0 is chance.** `nf w/ hot nbr` is the");
    println!("fraction of non-finite pixels having a high-drift neighbour -- magenta cores inside");
    println!("speckled halos. The user measured 1.7x and 6/6 under `fixed`.\n");
    println!("Chance is `1 - (1-base)^4` for base density `n_hot/N^2`, so **the ratio is not");
    println!("comparable across densities** -- it rises simply because the fix thins the mask.");
    println!("Read `n_hot` as the measurement and the ratio as a shape statistic. And where the");
    println!("mask SATURATES (deep interior, 92% hot) chance is ~1 and the ratio says nothing:");
    println!("that is the standing regional mask-saturation result, at a third statistic.\n");
    println!("{:<22} {:<14} {:>7} {:>8} {:>7} {:>8} {:>14}",
             "case", "mode", "n_hot", "nbr frac", "clust", "n_nonfin", "nf w/ hot nbr");
    for c in &cases {
        for (mode, mn) in MODES {
            let o = get(&c.name, mode);
            let hot: Vec<bool> = o.iter().map(|x| !x.drift.is_finite() || x.drift > HOT).collect();
            let (f, r, nh) = clustering(&hot, c.n);
            let nfv: Vec<usize> = (0..o.len()).filter(|&k| !o[k].finite).collect();
            let mut nf_nbr = 0usize;
            for &k in &nfv {
                let (i, j) = (k % c.n, k / c.n);
                for (di, dj) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                    let (a, b) = (i as i64 + di, j as i64 + dj);
                    if a >= 0 && b >= 0 && (a as usize) < c.n && (b as usize) < c.n
                        && hot[b as usize * c.n + a as usize] && o[b as usize * c.n + a as usize].finite
                    {
                        nf_nbr += 1;
                        break;
                    }
                }
            }
            println!("{:<22} {:<14} {:>7} {:>8.4} {:>7.3} {:>8} {:>7}/{:<6}",
                     c.name, mn, nh, f, r, nfv.len(), nf_nbr, nfv.len());
        }
        println!();
    }

    // -------------------------------------------------------------------------------------
    println!("== 4. THE TRAJECTORIES THAT GOT WORSE ==");
    println!("Max drift moved the WRONG way in the NumPy smoke test (3.43e-06 -> 6.93e-06). A");
    println!("handful of trajectories do worse under the fix; averaging that away is how a second");
    println!("mechanism stays hidden. The worst regressions are named with the quantities that");
    println!("would explain them.\n");
    println!("**Ranked by the ABSOLUTE rise, not the ratio.** A ratio ranking finds a small");
    println!("denominator, not a bad outcome: near-field's largest ratio is a pixel going");
    println!("9.2e-10 -> 4.5e-4, still among the region's best, on a region whose `drift max`");
    println!("IMPROVED 36x overall.\n");
    println!("{:<22} {:>6} {:>11} {:>11} {:>9} {:>9} {:>10} {:>10}",
             "case", "k", "fixed", "per-step-int", "steps fix", "steps new", "d_min", "ab_min");
    for c in &cases {
        let a = get(&c.name, DtauMode::FixedPerInterval);
        let b = get(&c.name, DtauMode::PerStepInterval);
        let mut worse: Vec<usize> = (0..a.len())
            .filter(|&k| a[k].drift.is_finite() && b[k].drift.is_finite() && b[k].drift > a[k].drift)
            .collect();
        worse.sort_by(|&x, &y| (b[y].drift - a[y].drift)
            .partial_cmp(&(b[x].drift - a[x].drift))
            .unwrap_or(std::cmp::Ordering::Equal));
        println!("{:<22} {} of {} regressed", c.name, worse.len(), a.len());
        for &k in worse.iter().take(3) {
            println!("{:<22} {:>6} {:>11.3e} {:>11.3e} {:>9} {:>9} {:>10.3e} {:>10.3e}",
                     "", k, a[k].drift, b[k].drift, a[k].steps, b[k].steps,
                     b[k].d_min_true, b[k].ab_min);
        }
        println!();
    }

    // -------------------------------------------------------------------------------------
    println!("== 5. T::TINY -- IS THE FLOOR EVER THE GUARD? ==");
    println!("`dtau` divides by `A*B`, floored at TINY (1e-300 at f64, 1e-37 at f32). At f32");
    println!("`TINY*TINY` UNDERFLOWS, so a doubly-degenerate state gives `dtau = inf` -- caught by");
    println!("the explicit `is_finite` test and NOT by the floor it is named for. `ab_min` is the");
    println!("raw product before the floor, so `floored` says whether the clamp ever bound.\n");
    println!("{:<22} {:<14} {:>12} {:>12} {:>9}", "case", "mode", "ab_min p1", "ab_min min", "floored");
    for c in &cases {
        for (mode, mn) in MODES {
            let o = get(&c.name, mode);
            let v: Vec<f64> = o.iter().map(|x| x.ab_min).collect();
            println!("{:<22} {:<14} {:>12.3e} {:>12.3e} {:>9}",
                     c.name, mn, q(&v, 0.01),
                     v.iter().cloned().filter(|x| x.is_finite()).fold(f64::INFINITY, f64::min),
                     o.iter().filter(|x| x.ab_floored).count());
        }
        println!();
    }
    println!("DONE");
}
