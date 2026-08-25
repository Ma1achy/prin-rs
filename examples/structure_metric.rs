//! **§2.2 settled by measurement: does structure REPLACE the signal, or MULTIPLY it?**
//!
//! The recommendation on record is multiply — *uncertain AND structured*, which floors the
//! uniform sea by construction while keeping the determinacy question `ensemble_spread` answers.
//! The recommendation does not decide. `error(B)` does.
//!
//! # Three targets, and why each is here
//!
//! - **`near-field`** — the region every prior criterion result was measured in.
//! - **`deep interior`** — a change that only improves near-field is tuning. This is the arm that
//!   makes the result mean something.
//! - **`preset_shape`** — the recognisable slice, and the **only tree in the corpus whose leaves
//!   are entirely its own decisions**: 0% camera veto, 8 `floor` + 8 `keep`. Everywhere else the
//!   veto stops 95%+ of leaves, so a criterion comparison there is partly a comparison of what
//!   the veto left over.
//!
//! # Two controls, and the second one is the point
//!
//! `structure_only` ranks on the structure term with **no signal in it at all**. Without it,
//! `multiply` beating `within` says nothing about structure: it could be the rescaling alone. And
//! `off` is the identity row — `multiply` must be read against it, not against the field.
//!
//! # How to misread this table
//!
//! **Check the oracle-to-random separation first.** Where they are close the metric is not
//! discriminating in that region, and no criterion result read from it means anything there.
//! `far` is deliberately absent for exactly that reason: its reference has a p1-p99 window of
//! `(1.3e-9, 1.1e-8)`, which is the integrator's arithmetic and not physics.
//!
//! **`error = 0` means "matches this sampling", not "correct".**
//!
//! **A criterion beating `greedy_oracle` indicates lookahead value, not a bug.** Greedy on
//! immediate delta-error is optimal only when gains are independent and immediately available,
//! and on a tree they are neither. Nothing here asserts it dominates.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Rank};
use prin_rs::output::colour::Scalar;
use prin_rs::output::plot::{palette, Figure, Series};
use prin_rs::quad::{Agg, Criterion, StructureMode};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

/// The colouring the criterion is scored under. **A criterion parameter, not a presentation one**
/// — `error(B)` measures image change, so what is displayed decides which quads matter.
const SHIPPING: metric::Colouring = metric::Colouring::Bivariate(Scalar::ShapeSpread);

fn main() {
    let levels: u32 = arg(1, 5);
    let n: usize = arg(2, 8);
    let tau: f64 = arg(3, 1e-4);
    let t_max: f64 = arg(4, 13.0);
    let res = (1usize << levels) * n;

    // `n_sync` scales with `t_max` or the rows are different discretisations.
    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let ens = EnsembleCfg { refine_flagged: false, t_max, n_sync, ..Default::default() };
    let full = ((1usize << (2 * (levels + 1))) - 1) / 3;

    let budgets: Vec<usize> = {
        let mut b = vec![5usize];
        while *b.last().unwrap() * 2 < full {
            b.push(b.last().unwrap() * 2 + 1);
        }
        b.push(full);
        b
    };

    println!("structure modes by error(B). levels {levels}, N={n}, E+1={}, res {res}^2, \
              tau={tau:.0e}, t={t_max}, n_sync={n_sync}",
             ens.n_extra + 1);
    println!("colouring: {} -- the one that ships. {full} quads per target.", SHIPPING.name());
    println!("hot rule for the structure term: the quad's own median (the relative mask).");
    println!();

    // (label, cx, cy, half, body, chart)
    let mut targets: Vec<(String, f64, f64, f64, usize, Chart)> = Vec::new();
    for &(region, cx, cy, body) in
        grid::REGIONS.iter().filter(|r| matches!(r.0, "near-field" | "deep interior"))
    {
        targets.push((region.to_string(), cx, cy, 0.05, body, Chart::BodyPlane));
    }
    let ps = Chart::preset_shape();
    targets.push(("preset_shape".into(), 0.0, 0.0, ps.default_half(), 0, ps));

    let runs: Vec<Rank> = vec![
        Rank::GreedyOracle,
        // The identity row. `multiply` is read against THIS, not against the field.
        Rank::Structured(StructureMode::Off, Criterion::Within, Agg::Median),
        Rank::Structured(StructureMode::Multiply, Criterion::Within, Agg::Median),
        // `Replace` discards the criterion, so this row is identically `structure_only`. Kept
        // once, as the arm named in the brief, and NOT repeated per criterion -- two rows
        // agreeing to five digits because they are the same expression is not evidence.
        Rank::Structured(StructureMode::Replace, Criterion::Within, Agg::Median),
        Rank::Structured(StructureMode::Off, Criterion::Between, Agg::Median),
        Rank::Structured(StructureMode::Multiply, Criterion::Between, Agg::Median),
        // The best criterion measured so far, with and without the term.
        Rank::Structured(StructureMode::Off, Criterion::FracHotBetween, Agg::Median),
        Rank::Structured(StructureMode::Multiply, Criterion::FracHotBetween, Agg::Median),
        // The threshold-free control on the whole mask family.
        Rank::Signal(Criterion::GradRms, Agg::Median),
        Rank::Signal(Criterion::LayoutRel, Agg::Median),
        // Structure with no signal in it: says whether `multiply` buys structure or rescaling.
        Rank::StructureOnly,
        Rank::Random(1),
        Rank::Random(2),
        Rank::Random(3),
        Rank::Random(4),
        Rank::Random(5),
    ];

    for (label, cx, cy, half, body, chart) in targets {
        let t0 = std::time::Instant::now();
        let cache = metric::build(
            &label, cx, cy, half, body, chart, levels, n, res, tau, &ens, SHIPPING,
        );
        println!("=== {label} ===  chart {}, half {half}, built in {:.1} s",
                 chart.name(), t0.elapsed().as_secs_f64());

        // **Distinct values before any curve.** Two different faults give the same flat error
        // curve -- a BAD ordering and NO ordering -- and error(B) alone cannot tell them apart.
        let keys: Vec<_> = cache.quads.keys().cloned().collect();
        println!("{:>26} {:>9} {:>8} {:>7}", "ranking", "distinct", "modal%", "nan%");
        for r in &runs {
            if matches!(r, Rank::Random(_) | Rank::GreedyOracle) {
                continue;
            }
            let vals: Vec<f64> = keys.iter().map(|&k| metric::score(&cache, k, *r)).collect();
            let nanf = vals.iter().filter(|v| !v.is_finite()).count() as f64 / vals.len() as f64;
            let mut bits: Vec<u64> = vals.iter().map(|v| v.to_bits()).collect();
            bits.sort_unstable();
            let distinct = { let mut b = bits.clone(); b.dedup(); b.len() };
            let modal = {
                let (mut best, mut run, mut i) = (1usize, 1usize, 1usize);
                while i < bits.len() {
                    run = if bits[i] == bits[i - 1] { run + 1 } else { 1 };
                    best = best.max(run);
                    i += 1;
                }
                best as f64 / bits.len() as f64
            };
            println!("{:>26} {distinct:>9} {:>7.1}% {:>6.1}%",
                     r.name(), 100.0 * modal, 100.0 * nanf);
        }
        println!();

        let mut rows: Vec<(String, Vec<f64>)> = Vec::new();
        for r in &runs {
            let pts = metric::replay(&cache, *r, full);
            rows.push((r.name(), metric::curve_at(&pts, &budgets)));
        }

        print!("{:>26}", "B =");
        for b in &budgets {
            print!(" {b:>9}");
        }
        println!();
        for (name, curve) in &rows {
            if name.starts_with("random") {
                continue;
            }
            print!("{name:>26}");
            for e in curve {
                print!(" {e:>9.5}");
            }
            println!();
        }
        let rnd: Vec<&Vec<f64>> =
            rows.iter().filter(|(n, _)| n.starts_with("random")).map(|(_, c)| c).collect();
        for (lbl, lo) in [("random lo", true), ("random hi", false)] {
            print!("{lbl:>26}");
            for j in 0..budgets.len() {
                let v: Vec<f64> = rnd.iter().map(|c| c[j]).collect();
                let x = if lo {
                    v.iter().cloned().fold(f64::INFINITY, f64::min)
                } else {
                    v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                };
                print!(" {x:>9.5}");
            }
            println!();
        }

        // Is the metric discriminating here at all? Read this before anything above.
        let mid = budgets.len() / 2;
        let oracle = rows[0].1[mid];
        let rlo = rnd.iter().map(|c| c[mid]).fold(f64::INFINITY, f64::min);
        println!();
        println!("  at B = {}: greedy_oracle {oracle:.5}, best random {rlo:.5}, separation {:.5}",
                 budgets[mid], rlo - oracle);
        if rlo - oracle < 1e-4 {
            println!("  ^ THE METRIC IS NOT DISCRIMINATING HERE. No row above means anything.");
        }
        println!();

        let cols = palette(rows.len());
        let series: Vec<Series> = rows
            .iter()
            .enumerate()
            .map(|(i, (nm, c))| {
                let pts: Vec<(f64, f64)> =
                    budgets.iter().zip(c).map(|(&b, &e)| (b as f64, e)).collect();
                let s = Series::new(nm.clone(), pts, cols[i]);
                if nm.starts_with("random") || nm.starts_with("greedy") { s.dashed() } else { s }
            })
            .collect();
        let stem = label.replace(' ', "_");
        Figure {
            title: format!("error(B) by structure mode -- {label}"),
            x_label: "budget B (quads)".into(),
            y_label: "mean per-pixel OKLab error".into(),
            series,
            notes: vec![
                format!("colouring {} -- the criterion is scored under what ships", SHIPPING.name()),
                "dashed: the controls. greedy_oracle is a strong reference, NOT a ceiling."
                    .into(),
                "`off` is the identity row; `multiply` must be read against it, and against \
                 `structure_only`, which has no signal in it at all."
                    .into(),
            ],
            y_lo: 1e-5,
            y_hi: 1.0,
        }
        .save(&format!("results/criterion/structure_{stem}_t{t_max}"))
        .unwrap();
    }

    println!("Read `off` against `multiply` against `structure_only`. If multiply beats off but");
    println!("structure_only is no better than random, the term is re-weighting the spread rather");
    println!("than finding structure -- which is a real outcome and not a failed measurement.");
}
