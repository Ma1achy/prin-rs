//! **Is the ranked tree SELECTIVE, or merely SPARSE?**
//!
//! The widened `tau x k_frac` sweep says `k_frac` is the knob: depth variance runs
//! `0.577 -> 2.109` from `k = 1` to `k = 0.25` while a whole decade of `tau` moves it
//! `1.900 -> 1.866`. But it also says the tree at `k = 0.25` holds **64 leaves against 1755**,
//! and a small tree is not automatically a good one. Under-refining everywhere produces exactly
//! the same headline: fewer leaves, higher depth variance, a tidier picture.
//!
//! **Only an equal-budget comparison separates them**, and `src/metric.rs` is already the
//! machine for it. Build the fully-refined reference once, then score every tree by
//! `Cache::error_of` -- the summed OKLab distance from what it displays to what the complete
//! tree displays. Three things are scored at the SAME leaf count `B`:
//!
//! - **the ranked descent's own leaf set**, at whatever `B` it settled on;
//! - **`greedy_oracle`** at that `B` -- a strong reference and deliberately not a ceiling, since
//!   a criterion beating it indicates lookahead value rather than a bug;
//! - **`random` over several seeds**, read as a band. A single random trace is a draw.
//!
//! The reading is the sign of the gap. Below the random band at its own `B`, the tree spent a
//! small budget in the right places and is selective. Inside or above it, the tree is sparse and
//! the depth variance is a picture of under-refinement.
//!
//! **What this cannot say.** `error = 0` means "matches the reference sampling", not "correct":
//! the reference is the complete tree at one sample per pixel, and at the screen floor which
//! side of a filament a pixel lands on is an accident of where its sample fell. The zero is
//! exactly locatable, which is what makes it good for *comparing* orderings, and it is not a
//! statement about image quality.
//!
//! Run: `cargo run --release --example equal_budget [levels] [n] [seeds]`

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Colouring, Key, Rank};
use prin_rs::quad::{Agg, Criterion};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Mode, SchedCfg, K_FRAC_RANKED, K_FRAC_UNRANKED};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

struct Target {
    name: &'static str,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
}

fn targets() -> Vec<Target> {
    let mut v = Vec::new();
    for &(region, cx, cy, body) in
        grid::REGIONS.iter().filter(|r| matches!(r.0, "near-field" | "deep interior"))
    {
        v.push(Target {
            name: if region == "near-field" { "near-field" } else { "deep_interior" },
            chart: Chart::BodyPlane, cx, cy, half: 0.05, body,
        });
    }
    let ps = Chart::preset_shape();
    v.push(Target { name: "preset_shape", chart: ps, cx: 0.0, cy: 0.0, half: ps.default_half(), body: 0 });
    v
}

/// Map a live quad's box onto the cache's uniform `(level, ix, iy)` index.
///
/// The two trees are built by different code paths, so this is the one joint where they meet and
/// it is checked rather than assumed: an index outside the level's range means the descent left
/// the cache's root box and the whole comparison is void.
fn key_of(c: &metric::Cache, cx: f64, cy: f64, level: u32) -> Option<Key> {
    let h = c.half / (1u64 << level) as f64;
    let ix = ((cx - (c.cx - c.half)) / (2.0 * h)).round() as i64;
    let iy = ((cy - (c.cy - c.half)) / (2.0 * h)).round() as i64;
    let lim = 1i64 << level;
    // `.round()` lands on the cell index only if the centre really is at `(2i+1)h`; a half-cell
    // offset would round to a neighbour silently, so the reconstruction is checked back.
    let (bx, by) = (c.cx - c.half + (2 * ix + 1) as f64 * h, c.cy - c.half + (2 * iy + 1) as f64 * h);
    if ix < 0 || iy < 0 || ix >= lim || iy >= lim
        || (bx - cx).abs() > h * 1e-6 || (by - cy).abs() > h * 1e-6 {
        return None;
    }
    Some((level, ix as u32, iy as u32))
}

fn main() {
    let levels: u32 = arg(1, 6);
    let n: usize = arg(2, 8);
    let seeds: u64 = arg(3, 5);
    let res = (1usize << levels) * n;
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };

    println!("equal-budget comparison. levels={levels}, N={n}, res={res}^2, E+1={}, t={}, f64.",
             ens.n_extra + 1, ens.t_max);
    println!("Colouring: the SHIPPING one (hue = shape sphere, lightness = spread_shape). A");
    println!("criterion must be scored under the colouring that will ship -- measured, the gap");
    println!("between the best criterion and random runs TOTAL under `outcome` and 6.6x under");
    println!("lightness=spread, and the best criterion changes identity between them.\n");

    // **Restricted to one target by default, and the sweep is the reason rather than the cost.**
    // `deep_interior` returns 29 quads / 22 leaves / depth variance 0.614 at EVERY `k_frac`, and
    // `preset_shape` returns 21/16/0.000 at every `k_frac`: `k` truncates the set that already
    // decided to split, and neither produces enough splits per round for a fraction to bite.
    // Identical trees give identical leaf sets and therefore identical `error_of`, so scoring
    // them would print the same number five times and read as a result. Pass a name to run one
    // of them anyway.
    let only: String = std::env::args().nth(4).unwrap_or_else(|| "near-field".into());
    for t in targets().iter().filter(|t| only == "all" || t.name == only) {
        println!("=== {} ===", t.name);
        let cache = metric::build(t.name, t.cx, t.cy, t.half, t.body, t.chart, levels, n, res,
                                  1e-4, &ens, Colouring::Bivariate(prin_rs::output::colour::Scalar::ShapeSpread));

        // The random band first. Where the band is narrow against the oracle gap the metric is
        // not discriminating and nothing read off it means anything -- so it is printed before
        // any tree is scored, not after a result is in hand.
        println!("{:>7} {:>10} {:>11} {:>11} {:>11} {:>9}",
                 "k_frac", "leaves B", "tree err", "greedy@B", "random@B", "verdict");

        let mut ks: Vec<f64> = vec![K_FRAC_UNRANKED, 0.5, K_FRAC_RANKED, 0.1, 0.05];
        ks.dedup();
        for k in ks {
            let cfg = SchedCfg {
                n, budget: 40000, tau_display: 1e-4, alpha_hi: 0.2, alpha_lo: 0.2,
                criterion: Criterion::Within, agg: Agg::Median, mode: Mode::Balanced, k_frac: k,
                camera: Some(Camera::framing(t.cx, t.cy, t.half, res)), chart: t.chart,
                ..Default::default()
            };
            let (tree, _) = scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, &ens,
                                               Precision::F64);
            let mut keys: Vec<Key> = Vec::new();
            let mut missed = 0usize;
            for i in tree.leaves() {
                let q = &tree.nodes[i];
                match key_of(&cache, q.cx, q.cy, q.level) {
                    Some(kk) if q.level <= levels && cache.quads.contains_key(&kk) => keys.push(kk),
                    _ => missed += 1,
                }
            }
            if missed > 0 {
                // Deeper than the cache, or off its box. Reported, never dropped quietly: a
                // partial leaf set does not tile the root and `error_of` over it is meaningless.
                println!("{k:>7.2} {:>10} {:>11} {:>11} {:>11} {:>9}",
                         keys.len() + missed, "-", "-", "-", "NOT SCORED");
                println!("        ^ {missed} leaves fall outside the cache (deeper than level \
                          {levels} or off the root box), so the leaf set does not tile and \
                          error_of is undefined.");
                continue;
            }
            let b = keys.len();
            let e = cache.error_of(&keys);

            let g = *metric::curve_at(&metric::replay(&cache, Rank::GreedyOracle, b + 1), &[b])
                .first().unwrap_or(&f64::NAN);
            let mut rs: Vec<f64> = (0..seeds)
                .map(|s| *metric::curve_at(&metric::replay(&cache, Rank::Random(s), b + 1), &[b])
                    .first().unwrap_or(&f64::NAN))
                .collect();
            rs.sort_by(|a, c| a.partial_cmp(c).unwrap());
            let (rlo, rhi) = (rs[0], rs[rs.len() - 1]);

            let verdict = if !e.is_finite() { "-" }
                else if e < rlo { "SELECTIVE" }
                else if e > rhi { "SPARSE" }
                else { "in band" };
            println!("{k:>7.2} {b:>10} {e:>11.5} {g:>11.5} {:>11} {verdict:>9}",
                     format!("{rlo:.5}-{rhi:.5}"));
        }
        println!();
    }

    println!("SELECTIVE means the tree's error at its own leaf count is BELOW the random band at");
    println!("that same count: a small budget spent in the right places. SPARSE means above it --");
    println!("the same leaf count spent at random would have displayed more. `in band` is the");
    println!("honest third answer and is not a failure to measure; it says the ordering is not");
    println!("distinguishable from arbitrary at that budget.");
    println!();
    println!("greedy@B is a REFERENCE, not a ceiling. It is greedy on immediate delta-error, and");
    println!("on a tree gains are neither independent nor immediately available -- measured at");
    println!("t = 20 in near-field, greedy plateaus at 0.00048 while first_div reaches 0.00000.");
    println!("A tree beating it indicates lookahead value. Nothing here asserts it dominates.");
}
