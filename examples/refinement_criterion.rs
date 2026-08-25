//! **Experiment A (BRIEF §8): the refinement criterion, tested without a scheduler.**
//!
//! The criterion compares a parent quad against its children. A fine uniform grid already
//! *contains* every coarser scale by aggregation: pool the 2x2 block's copies to synthesise
//! the parent, compare against the children, and the whole exponent machinery is testable with
//! no quadtree at all. The absence of one is a feature of this test — nothing here can be an
//! artefact of a scheduler, because there is no scheduler.
//!
//! `alpha = log2(spread_parent / spread_child)`.
//!
//! **`alpha` for `sigma_E(0)` is the control and its true value is exactly 1.0.** `sigma_E(0)`
//! is proportional to the jitter and therefore to the cell width, so doubling the cell must
//! double it. Any estimator that cannot return 1.0 there is broken, and the deviation is a
//! direct error measure — which is how the sample-size bias below was found rather than
//! assumed.

use rayon::prelude::*;

use prin_rs::ensemble::jitter::Scheme;
use prin_rs::ensemble::{jitter, stats};
use prin_rs::grid::{self, Slice};
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::physics::{energy, shape};

const FINE: usize = 64;
const N_EXTRA: usize = 7;

/// Everything one pixel's ensemble contributes to a block statistic.
struct Copies {
    e0: Vec<f64>,
    et: Vec<f64>,
    shapes: Vec<[f64; 3]>,
}

fn copies_of(s: &Slice, idx: usize) -> Copies {
    copies_of_with(s, idx, N_EXTRA, Scheme::default())
}

fn copies_of_with(s: &Slice, idx: usize, n_extra: usize, scheme: Scheme) -> Copies {
    let c = jitter::copies_with::<f64>(s, idx, n_extra, 0.5, 0, scheme);
    let mut e0 = Vec::new();
    let mut et = Vec::new();
    let mut shapes = Vec::new();
    for x in &c {
        e0.push(energy::energy(&x.s.r, &x.s.v, &x.m, 0.0));
        let o = az::integrate_az_opts(
            x.s, &x.m, 13.0, 32, 0.01, 30_000,
            &AzOpts { r_coll_frac: 1e-3, stop_on_event: true, ..Default::default() },
        );
        et.push(energy::energy(&o.state.r, &o.state.v, &x.m, 0.0));
        shapes.push(shape::shape_vec(&o.state.r, &x.m));
    }
    Copies { e0, et, shapes }
}

/// Root-mean-square deviation from the median. Unlike the max, its expectation does not grow
/// with the sample size — which matters here, because the parent pools 4x as many copies as a
/// child and any order statistic would rise for that reason alone.
fn rms_dev(v: &[f64]) -> f64 {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    if x.len() < 2 {
        return f64::NAN;
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if x.len() % 2 == 1 { x[x.len() / 2] } else { 0.5 * (x[x.len() / 2 - 1] + x[x.len() / 2]) };
    (x.iter().map(|q| (q - med).powi(2)).sum::<f64>() / x.len() as f64).sqrt()
}

fn spread_shape_of(sh: &[[f64; 3]]) -> f64 {
    shape::spread_shape(sh)
}

fn qs(v: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mut x: Vec<f64> = v.iter().cloned().filter(|q| q.is_finite()).collect();
    if x.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| x[(((x.len() - 1) as f64) * f).round() as usize];
    (q(0.0), q(0.1), q(0.5), q(0.9), q(1.0))
}

fn run_region(region: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let s = grid::region(region, FINE, FINE, 0.05).unwrap();
    let all: Vec<Copies> = (0..s.npix()).into_par_iter().map(|i| copies_of(&s, i)).collect();

    let mut a_e0_rms = Vec::new();
    let mut a_e0_max = Vec::new();
    let mut a_et = Vec::new();
    let mut a_shape = Vec::new();

    // Parents are 2x2 blocks of the fine grid. Row-major, x fastest: idx = jy*nx + jx.
    for py in 0..FINE / 2 {
        for px in 0..FINE / 2 {
            let kids = [
                2 * py * FINE + 2 * px,
                2 * py * FINE + 2 * px + 1,
                (2 * py + 1) * FINE + 2 * px,
                (2 * py + 1) * FINE + 2 * px + 1,
            ];
            let pool = |f: fn(&Copies) -> &Vec<f64>| -> Vec<f64> {
                kids.iter().flat_map(|&k| f(&all[k]).iter().cloned()).collect()
            };
            let pe0 = pool(|c| &c.e0);
            let pet = pool(|c| &c.et);
            let pshape: Vec<[f64; 3]> = kids.iter().flat_map(|&k| all[k].shapes.iter().cloned()).collect();

            // The child value is the median over the four children, so one wild child does not
            // set the exponent for the quad.
            let child = |f: &dyn Fn(usize) -> f64| -> f64 {
                let mut v: Vec<f64> = kids.iter().map(|&k| f(k)).filter(|x| x.is_finite()).collect();
                if v.is_empty() {
                    return f64::NAN;
                }
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                v[v.len() / 2]
            };

            let lg = |p: f64, c: f64| if p > 0.0 && c > 0.0 { (p / c).log2() } else { f64::NAN };
            a_e0_rms.push(lg(rms_dev(&pe0), child(&|k| rms_dev(&all[k].e0))));
            a_e0_max.push(lg(
                stats::max_dev(&pe0),
                child(&|k| stats::max_dev(&all[k].e0)),
            ));
            a_et.push(lg(rms_dev(&pet), child(&|k| rms_dev(&all[k].et))));
            a_shape.push(lg(
                spread_shape_of(&pshape),
                child(&|k| spread_shape_of(&all[k].shapes)),
            ));
        }
    }
    (a_e0_rms, a_e0_max, a_et, a_shape)
}

/// The control at `t = 0` needs no integration at all — `sigma_E(0)` is a function of the
/// initial conditions — so its noise floor can be measured against ensemble size cheaply.
///
/// Two things to separate: the **scatter** (how far a single quad's alpha can be from 1.0 by
/// chance) and the **bias** (the parent pools `4(E+1)` samples against a child's `E+1`, and a
/// spread estimator's expectation depends on sample size). The subsampled column removes the
/// bias by drawing `E+1` of the parent's pooled copies, leaving only the scatter.
fn control_noise_floor() {
    for scheme in [Scheme::Halton, Scheme::Pcg] {
        control_noise_floor_for(scheme);
    }
}

fn control_noise_floor_for(scheme: Scheme) {
    use prin_rs::rng::SplitMix64;
    println!("=== the control's noise floor against ensemble size — {scheme:?} ===");
    println!("alpha for sigma_E(0). True value is exactly 1.0: sigma_E(0) is proportional to");
    println!("the jitter and so to cell width, so doubling the cell doubles it.");
    println!();
    println!("{:>8}{:>12}{:>10}{:>10}{:>10}{:>14}{:>10}",
             "E+1", "estimator", "p10", "median", "p90", "p90-p10", "subsamp med");

    for n_copies in [4usize, 8, 16, 32, 64] {
        let s = grid::region("near-field", FINE, FINE, 0.05).unwrap();
        let e0: Vec<Vec<f64>> = (0..s.npix())
            .into_par_iter()
            .map(|i| {
                jitter::copies_with::<f64>(&s, i, n_copies - 1, 0.5, 0, scheme)
                    .iter()
                    .map(|x| energy::energy(&x.s.r, &x.s.v, &x.m, 0.0))
                    .collect()
            })
            .collect();

        for (name, est) in [
            ("rms", rms_dev as fn(&[f64]) -> f64),
            ("max_dev", |v: &[f64]| stats::max_dev(v)),
        ] {
            let mut a = Vec::new();
            let mut a_sub = Vec::new();
            let mut rng = SplitMix64::new(0xA5A5_1234);
            for py in 0..FINE / 2 {
                for px in 0..FINE / 2 {
                    let kids = [
                        2 * py * FINE + 2 * px,
                        2 * py * FINE + 2 * px + 1,
                        (2 * py + 1) * FINE + 2 * px,
                        (2 * py + 1) * FINE + 2 * px + 1,
                    ];
                    let pool: Vec<f64> = kids.iter().flat_map(|&k| e0[k].iter().cloned()).collect();
                    let mut cv: Vec<f64> = kids.iter().map(|&k| est(&e0[k])).collect();
                    cv.sort_by(|x, y| x.partial_cmp(y).unwrap());
                    let child = cv[cv.len() / 2];
                    a.push((est(&pool) / child).log2());
                    // Same estimator, same sample count on both sides.
                    let mut sub = pool.clone();
                    for j in (1..sub.len()).rev() {
                        let r = (rng.next_u64() % (j as u64 + 1)) as usize;
                        sub.swap(j, r);
                    }
                    a_sub.push((est(&sub[..n_copies]) / child).log2());
                }
            }
            let (_, p10, med, p90, _) = qs(&a);
            let (_, _, smed, _, _) = qs(&a_sub);
            println!("{n_copies:>8}{name:>12}{p10:>10.4}{med:>10.4}{p90:>10.4}{:>14.4}{smed:>10.4}",
                     p90 - p10);
        }
    }
    println!();
    println!("p90-p10 is the noise floor: the width a single quad's alpha scatters over when");
    println!("the true value is exactly 1.0. Any per-quad refinement decision finer than that");
    println!("width is reading noise. The subsampled column matches sample counts on both");
    println!("sides, so the gap between it and the median column is pure sample-size bias.");
    println!();
}

fn main() {
    control_noise_floor();
    println!("Experiment A: the refinement criterion by aggregation, no quadtree.");
    println!("{FINE}x{FINE} fine grid -> {}x{} parents, E+1=8 copies, t=13, eta=0.01, f64",
             FINE / 2, FINE / 2);
    println!("alpha = log2(spread_parent / spread_child); parent pools the 2x2 block's copies,");
    println!("child is the median over the four.");
    println!();

    for region in ["near-field", "mid-field", "body2 core", "far"] {
        let (e0r, e0m, et, sh) = run_region(region);
        println!("=== {region} ===");
        println!("{:>34}{:>10}{:>10}{:>10}{:>10}{:>10}", "quantity", "min", "p10", "median", "p90", "max");
        for (name, v) in [
            ("alpha sigma_E(0), rms  [truth 1.0]", &e0r),
            ("alpha sigma_E(0), max_dev [control]", &e0m),
            ("alpha sigma_E(t), rms", &et),
            ("alpha spread_shape", &sh),
        ] {
            let (a, b, c, d, e) = qs(v);
            println!("{name:>34}{a:>10.4}{b:>10.4}{c:>10.4}{d:>10.4}{e:>10.4}");
        }
        println!();
    }
}
