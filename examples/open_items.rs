//! **The two open items carried forward from the PR #11 review.**
//!
//! ## 1. Does `floored` correlate with `worst_energy_drift` in `deep interior`?
//!
//! 40.9% floored is the highest of the three regions, and it is also where the integrator works
//! hardest. If the floor fires because `alpha` is corrupted by *integration error* rather than
//! because the physics is irreducible, that is a different bug with a different fix — and the
//! floor fraction alone cannot tell them apart.
//!
//! **Reported with and without the veto, with `n` and the level distribution of the floored
//! quads under each.** If the floored population was mostly quads descending past level 6, the
//! veto removes the very sample the correlation is measured on, and a weak Spearman would mean
//! *"no data"*, not *"no relationship"*. Without those columns the null is unreadable.
//!
//! ## 2. Does p90 aggregation fix `deep interior`'s tree?
//!
//! The median-blindness attribution for q3 is plausible and unproven. **If p90 descends where
//! median does not, it is aggregation. If p90 also stalls, it is not.**
//!
//! Run: `cargo run --release --example open_items [budget]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::quad::{Agg, Decision as D};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};

/// Spearman rank correlation, ties handled by average rank.
fn spearman(x: &[f64], y: &[f64]) -> f64 {
    fn ranks(v: &[f64]) -> Vec<f64> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
        let mut r = vec![0.0; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0 + 1.0;
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    }
    if x.len() < 3 {
        return f64::NAN;
    }
    let (rx, ry) = (ranks(x), ranks(y));
    let n = x.len() as f64;
    let (mx, my) = (rx.iter().sum::<f64>() / n, ry.iter().sum::<f64>() / n);
    let num: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let dx: f64 = rx.iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let dy: f64 = ry.iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
    if dx == 0.0 || dy == 0.0 {
        f64::NAN
    } else {
        num / (dx * dy)
    }
}

fn main() {
    let budget: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(6000);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("=== item 1: is the floor firing on integration error? ===");
    println!("budget {budget}, tau=1e-4, alpha_hi=0.2, N=8, E+1=8, t=13, f64.\n");
    println!("{:>14} {:>6} {:>7} {:>8} {:>10} {:>12} {:>12} {:>20}",
             "region", "veto", "leaves", "floored", "spearman", "drift med F", "drift med K",
             "floored levels");

    for region in ["deep interior", "near-field", "far"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        for veto in [false, true] {
            let cfg = SchedCfg {
                budget, tau_display: 1e-4, alpha_hi: 0.2, alpha_lo: 0.08,
                camera: veto.then(|| Camera::framing(root.cx, root.cy, 0.05, 512)),
                ..Default::default()
            };
            let (t, _) = scheduler::descend(
                root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
            let leaves: Vec<usize> = t.leaves().collect();
            // Correlation over ALL leaves: floored as 1/0 against the worst drift in the quad.
            let (mut f, mut d) = (Vec::new(), Vec::new());
            for &i in &leaves {
                let q = &t.nodes[i];
                if q.red.worst_energy_drift.is_finite() {
                    f.push(if q.decision == D::Floor { 1.0 } else { 0.0 });
                    d.push(q.red.worst_energy_drift);
                }
            }
            let med = |mut v: Vec<f64>| -> f64 {
                if v.is_empty() { return f64::NAN; }
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                v[v.len() / 2]
            };
            let df: Vec<f64> = leaves.iter().filter(|&&i| t.nodes[i].decision == D::Floor)
                .map(|&i| t.nodes[i].red.worst_energy_drift).filter(|x| x.is_finite()).collect();
            let dk: Vec<f64> = leaves.iter().filter(|&&i| t.nodes[i].decision == D::Keep)
                .map(|&i| t.nodes[i].red.worst_energy_drift).filter(|x| x.is_finite()).collect();
            let mut lv: Vec<u32> = leaves.iter().filter(|&&i| t.nodes[i].decision == D::Floor)
                .map(|&i| t.nodes[i].level).collect();
            lv.sort_unstable();
            let mut hist = String::new();
            let mut k = 0;
            while k < lv.len() {
                let l = lv[k];
                let c = lv.iter().filter(|&&x| x == l).count();
                hist.push_str(&format!("{l}:{c} "));
                k += c;
            }
            println!("{:>14} {:>6} {:>7} {:>8} {:>10.4} {:>12.3e} {:>12.3e} {:>20}",
                     region, if veto { "on" } else { "OFF" }, leaves.len(), df.len(),
                     spearman(&f, &d), med(df), med(dk),
                     if hist.is_empty() { "-".into() } else { hist.trim().to_string() });
        }
    }
    println!();
    println!("Read `floored` and `floored levels` before the Spearman. If the veto removes the");
    println!("levels the floored quads lived at, a weak correlation means NO DATA, not no");
    println!("relationship — the population it was measured on is gone.");

    println!("\n=== item 2: does p90 aggregation fix deep interior's tree? ===");
    println!("PR #11 q3: the tree leaves the largest high-spread structures at level 2. The");
    println!("median-blindness attribution is plausible and unproven.\n");
    println!("{:>14} {:>6} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>22}",
             "region", "veto", "agg", "quads", "leaves", "depth", "floor", "keep", "depth histogram");
    for region in ["deep interior", "near-field"] {
        let root = grid::region(region, 2, 2, 0.05).unwrap();
        for veto in [false, true] {
            for agg in [Agg::Median, Agg::Mean, Agg::P90] {
                let cfg = SchedCfg {
                    budget, tau_display: 1e-4, alpha_hi: 0.2, alpha_lo: 0.08, agg,
                    camera: veto.then(|| Camera::framing(root.cx, root.cy, 0.05, 512)),
                    ..Default::default()
                };
                let (t, st) = scheduler::descend(
                    root.cx, root.cy, 0.05, root.body, &cfg, &ens, Precision::F64);
                let leaves: Vec<usize> = t.leaves().collect();
                let c = |d: D| leaves.iter().filter(|&&i| t.nodes[i].decision == d).count();
                let hist = t.depth_histogram();
                let hs: String = hist.iter().enumerate().filter(|(_, &n)| n > 0)
                    .map(|(l, n)| format!("{l}:{n} ")).collect();
                println!("{:>14} {:>6} {:>8} {:>7} {:>7} {:>6} {:>6} {:>6} {:>22}",
                         region, if veto { "on" } else { "OFF" }, agg.name(),
                         st.quads_computed, leaves.len(), hist.len().saturating_sub(1),
                         c(D::Floor), c(D::Keep), hs.trim());
            }
        }
        println!();
    }
    println!("If p90 descends where median does not, q3's failure is AGGREGATION. If p90 also");
    println!("stalls at level 2-4, it is not, and the attribution has to be withdrawn.");
}
