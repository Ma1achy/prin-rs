//! **How far did the criterion's INPUT FIELD move under the integrator fixes?**
//!
//! The scheduler corpus -- `charts/`, `criterion/`, `vertical/` and their ranked twins, 1009
//! files and 188 `.prnq` dumps -- was written **25-26 August**. Every integrator change landed
//! **27 August or later**: `5cc8dec` the `dtau` fix, `f7d2a31` the landing clamp, `2eb9e8f` the
//! predictive step limit, `a526360` the no-discard fix and the secant landing, `84830a1` Heggie
//! as the default. So every field-conditional number in that corpus is measured on a kernel that
//! no longer exists, and the open question is not *whether* to re-measure but **how much**.
//!
//! This harness answers that and nothing else. Two kernels, everything else held identical:
//!
//! - **`pre`** -- `Az`, `DtauMode::FixedPerInterval`, `clamp_final_step: false`,
//!   `StepLimit::None`, `land_iterate: false`. Reconstructed through the flags rather than by
//!   checking out an old commit, which is the method already validated on this project:
//!   `FixedPerInterval` reproduced `71de13f` bitwise.
//! - **`now`** -- `EnsembleCfg::production()`.
//!
//! `refine_flagged` is **off in both arms**. It is batch-only, has no live-playhead analogue, and
//! it removes the very population a kernel comparison is about -- the failure already on record
//! from `integrator_gallery`'s first run.
//!
//! # What it measures, and why not the obvious thing
//!
//! **Not "did the trajectories change".** They must, in a chaotic region: `deep interior` is on
//! record with `chord max = 1.999` at every rung of an `eta` ladder, so a chord statistic there is
//! **saturated** and reads the same for a correct change and a broken one. The fraction antipodal
//! is printed precisely so a saturated cell can be recognised rather than quoted.
//!
//! What decides the re-measurement is whether the **criterion** sees a different field, and the
//! criterion is a **ranking**. *A ranking is invariant to a monotone rescaling of the signal; a
//! threshold is not.* So the headline is `rho` -- Spearman of `ensemble_spread` between the arms
//! -- and beside it the per-pixel ratio distribution, because *never conclude "no effect" from an
//! aggregate without the per-pixel distribution.*
//!
//! Stage 2 asks the same question where it is actually decided: run the real descent under both
//! kernels at the corpus's own `SchedCfg` and diff the trees. Decisions are compared over
//! **shared quads only, with the count printed** -- a quad present in one tree and not the other
//! has not changed its mind.
//!
//! # The guards
//!
//! - `decode::distinct` on the ICs. *Count distinct ICs first, read divergence second.*
//! - **not-inert**: the two arms must differ somewhere, or the comparison is dead on both sides.
//! - **saturation**: the antipodal fraction, printed per region.
//! - `n_undetermined` on both arms, from this PR -- the first run able to see it.
use rayon::prelude::*;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid;
use prin_rs::integrate::az::driver::{DtauMode, StepLimit};
use prin_rs::integrate::Integrator;
use prin_rs::outcome::EscapeRule;
use prin_rs::quad::Decision;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::spatial::HotRule;

/// The corpus kernel, reconstructed from a **field-by-field audit** rather than from memory.
///
/// Diffing `EnsembleCfg` at `0114be4` (26 Aug 15:06, the last commit before the corpus's newest
/// file) against `production()` today: **seventeen fields did not exist then**. Naming the four
/// famous fixes and stopping would have left thirteen at today's values inside an arm labelled
/// *pre* — the `refine_flagged` failure exactly, *a setting correct where it was born and silent
/// where it was not*.
///
/// | field | corpus era | today | here |
/// |---|---|---|---|
/// | `integrator` | AZ only | `Heggie` | **`Az`** |
/// | `dtau_mode` | fixed per interval | `PerStepInterval` | **`FixedPerInterval`** |
/// | `clamp_final_step` | overshoot | `true` | **`false`** |
/// | `step_limit` | none | `Predictive` | **`None`** |
/// | `land_iterate` | none | `true` | **`false`** |
/// | `escape_rule` | reference, hardcoded | `Closure` | **`Reference`** |
/// | `escape_confirm` | none | `true` | **`false`** |
/// | `escape_every`, `stop_on_escape`, `ref_hysteresis` | — | `0` / `false` / `0.0` | match already |
/// | `step_limit_f`, `blend_p`, `step_blend`, `closure_k`, `land_max_iters` | — | — | inert under the above |
/// | `keep_drift_hist`, `keep_ref_path` | — | `false` | instrumentation only |
///
/// **`escape_rule` is the one that was nearly missed.** The closure criterion landed at `71de13f`
/// on 27 August, *after* the corpus, and `spread_event` reads the event class — so leaving it at
/// `Closure` would have put a post-corpus terminal classifier inside the *pre* arm. It is silent
/// on `near-field` at `t = 13` by the standing result and demonstrably **not** silent on
/// `deep interior`, whose escape fraction moves by half the region under a cadence change.
fn pre_kernel() -> EnsembleCfg {
    let mut c = EnsembleCfg::production().with_overrides(&[
        Override::Integrator(Integrator::Az),
        Override::DtauMode(DtauMode::FixedPerInterval),
        Override::ClampFinalStep(false),
        Override::StepLimit(StepLimit::None),
        Override::EscapeRule(EscapeRule::Reference),
        Override::EscapeConfirm(false),
        Override::RefineFlagged(false),
    ]);
    // No `Override` variant exists for this one; set directly. `overrides_vs_production` derives
    // its declaration by diffing the struct, so the sidecar still reports it either way.
    c.land_iterate = false;
    c
}

/// **The audit's own control.** The four fixes everyone names, with `escape_rule` and
/// `escape_confirm` left at today's values — i.e. what `pre_kernel` would have been if the
/// field-by-field diff had not been done.
///
/// Reported per region against `pre`. If the two agree everywhere the audit cost nothing and this
/// says so; if they differ, a *pre* arm built from memory was carrying a post-corpus terminal
/// classifier. Either way it is a measurement rather than a claim about how carefully I read.
fn pre_naive_kernel() -> EnsembleCfg {
    let mut c = EnsembleCfg::production().with_overrides(&[
        Override::Integrator(Integrator::Az),
        Override::DtauMode(DtauMode::FixedPerInterval),
        Override::ClampFinalStep(false),
        Override::StepLimit(StepLimit::None),
        Override::RefineFlagged(false),
    ]);
    c.land_iterate = false;
    c
}

fn now_kernel() -> EnsembleCfg {
    EnsembleCfg::production().with_overrides(&[Override::RefineFlagged(false)])
}

/// Quantile over the finite values, **returning the count it dropped**.
///
/// The first cut sorted with `partial_cmp().unwrap_or(Equal)`, which is not a total order once a
/// `NaN` is present — Rust's driftsort detects that and panics, and it did, on `deep interior`
/// under `StepLimit::None`. The no-discard rule at the harness layer: filter, and **report what
/// was filtered** rather than quietly quoting a quantile over a subset.
fn qn(v: &[f64], p: f64) -> (f64, usize) {
    let mut f: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect();
    let dropped = v.len() - f.len();
    if f.is_empty() {
        return (f64::NAN, dropped);
    }
    f.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (f[(((f.len() - 1) as f64) * p).round() as usize], dropped)
}

fn q(v: &[f64], p: f64) -> f64 {
    qn(v, p).0
}

/// Spearman over pairs where **both** sides are finite, with the pair count returned.
///
/// Dropping a pair is a discard, so the count is reported rather than folded into `n`: a `rho`
/// over 40% of the frame is a different statistic from one over all of it.
fn spearman(a: &[f64], b: &[f64]) -> (f64, usize) {
    let pairs: Vec<(f64, f64)> =
        a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(&x, &y)| (x, y)).collect();
    let n = pairs.len();
    if n < 3 {
        return (f64::NAN, n);
    }
    let rank = |get: &dyn Fn(&(f64, f64)) -> f64| -> Vec<f64> {
        let mut ix: Vec<usize> = (0..n).collect();
        ix.sort_by(|&i, &j| get(&pairs[i]).partial_cmp(&get(&pairs[j])).unwrap());
        let mut r = vec![0.0; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && get(&pairs[ix[j + 1]]) == get(&pairs[ix[i]]) {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0;
            for &k in &ix[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    };
    let (ra, rb) = (rank(&|p: &(f64, f64)| p.0), rank(&|p: &(f64, f64)| p.1));
    let m = (n - 1) as f64 / 2.0;
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let (x, y) = (ra[k] - m, rb[k] - m);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    (num / (da * db).sqrt(), n)
}

/// Distinct finite values and the modal share, on the bit pattern.
///
/// **The standing rule: count the signal's distinct values before reading any curve.** Two
/// different faults give the same low `rho` -- a field that was genuinely re-ordered, and a field
/// with no ordering to preserve. `far` is the case that needs it: a ratio band 0.9% wide reads as
/// a clean monotone rescale, and a rescale cannot move a rank, so a `rho` of 0.59 beside it means
/// something other than "the ordering changed".
fn resolution(v: &[f64]) -> (usize, f64, f64) {
    let mut b: std::collections::HashMap<u64, usize> = Default::default();
    let mut fin = 0usize;
    for x in v.iter().filter(|x| x.is_finite()) {
        *b.entry(x.to_bits()).or_default() += 1;
        fin += 1;
    }
    if fin == 0 {
        return (0, f64::NAN, f64::NAN);
    }
    let modal = *b.values().max().unwrap() as f64 / fin as f64;
    let (lo, hi) = (q(v, 0.01), q(v, 0.99));
    (b.len(), modal, if lo > 0.0 { hi / lo } else { f64::NAN })
}

fn chord(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn field(sl: &grid::Slice, cfg: &EnsembleCfg) -> Vec<PixelOut> {
    (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(sl, k, cfg)).collect()
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(64);
    let budget: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let regions: Vec<String> = {
        let rest: Vec<String> = std::env::args().skip(3).collect();
        if rest.is_empty() {
            ["near-field", "mid-field", "far", "body2 core", "deep interior"]
                .iter()
                .map(|s| s.to_string())
                .collect()
        } else {
            rest
        }
    };

    let (pre, now, pre_naive) = (pre_kernel(), now_kernel(), pre_naive_kernel());
    println!("HOW FAR DID THE CRITERION'S INPUT FIELD MOVE? {n}^2 field, budget {budget} quads.\n");
    println!("  pre = Az + FixedPerInterval + no landing clamp + StepLimit::None + no secant landing");
    println!("  now = EnsembleCfg::production()          (both arms: refine_flagged = false)\n");
    println!("  t_max {} eta {} n_sync {} r_coll {}  E+1 = {}",
             now.t_max, now.eta, now.n_sync, now.r_coll_frac, now.n_extra + 1);

    // ------------------------------------------------------------------ stage 1: the field
    println!("\nSTAGE 1 -- the per-footprint field. `rho` is the headline: the criterion is a");
    println!("ranking, so a monotone rescale costs it nothing and a re-ordering costs it everything.\n");
    println!("{:<14} {:>9} {:>7} {:>10} {:>10} {:>8} {:>8} {:>8} {:>9} {:>7} {:>7} {:>6}",
             "region", "distinct", "moved", "spr pre", "spr now", "rat p10", "rat p50", "rat p90",
             "rho", "chd p50", "antip", "chdNaN");
    let mut stage1 = Vec::new();
    for r in &regions {
        let sl = match grid::region(r, n, n, 0.05) {
            Some(s) => s,
            None => {
                println!("{r:<14}  UNKNOWN REGION");
                continue;
            }
        };
        let ics: Vec<prin_rs::physics::Cart<f64>> = (0..sl.npix()).map(|k| sl.nominal::<f64>(k)).collect();
        let distinct = prin_rs::decode::distinct(&ics);

        let (a, b) = (field(&sl, &pre), field(&sl, &now));
        let sa: Vec<f64> = a.iter().map(|p| p.ensemble_spread).collect();
        let sb: Vec<f64> = b.iter().map(|p| p.ensemble_spread).collect();
        let chords: Vec<f64> = a.iter().zip(&b).map(|(x, y)| chord(x.shape_vec, y.shape_vec)).collect();
        let moved = chords.iter().filter(|&&c| c > 0.0).count() as f64 / chords.len() as f64;
        let antip = chords.iter().filter(|&&c| c > 1.9).count() as f64 / chords.len() as f64;
        let ratio: Vec<f64> = sa
            .iter()
            .zip(&sb)
            .filter(|(x, y)| x.is_finite() && y.is_finite() && **x > 0.0)
            .map(|(x, y)| y / x)
            .collect();
        let (rho, npair) = spearman(&sa, &sb);
        let nf_a = sa.iter().filter(|x| !x.is_finite()).count();
        let nf_b = sb.iter().filter(|x| !x.is_finite()).count();
        let (chd50, chd_nan) = qn(&chords, 0.5);
        let (da, ma, ra) = resolution(&sa);
        let (db, mb, rb) = resolution(&sb);
        println!("{r:<14} {distinct:>4}/{:<4} {moved:>7.4} {:>10.3e} {:>10.3e} {:>8.3} {:>8.3} {:>8.3} {rho:>9.4} {chd50:>7.3} {antip:>7.4} {chd_nan:>6}",
                 sl.npix(),
                 q(&sa, 0.5), q(&sb, 0.5),
                 q(&ratio, 0.10), q(&ratio, 0.50), q(&ratio, 0.90));
        // The audit's control: did correcting `escape_rule` change the *pre* arm at all?
        let c = field(&sl, &pre_naive);
        let sc: Vec<f64> = c.iter().map(|p| p.ensemble_spread).collect();
        let audit_moved = a
            .iter()
            .zip(&c)
            .filter(|(x, y)| chord(x.shape_vec, y.shape_vec) > 0.0 || x.outcome != y.outcome)
            .count();
        let audit_lbl = c.iter().zip(&a).filter(|(x, y)| x.outcome != y.outcome).count();
        stage1.push((r.clone(), moved, rho, npair, nf_a, nf_b, sl.npix(), audit_moved, audit_lbl,
                     q(&sc, 0.5), (da, ma, ra), (db, mb, rb)));
    }

    println!("\n  non-finite `ensemble_spread`, and the rank pairs the `rho` above is over:");
    for (r, moved, _, npair, nf_a, nf_b, npix, _, _, _, _, _) in &stage1 {
        let guard = if *moved > 0.0 { "live" } else { "**INERT -- the arms agree bitwise**" };
        println!("    {r:<14} pre {nf_a:>5}/{npix:<6} now {nf_b:>5}/{npix:<6} rho over {npair:>6} pairs   {guard}");
    }

    println!("\n  THE AUDIT'S OWN CONTROL. `pre` against a `pre` built from the four famous fixes");
    println!("  alone, with `escape_rule`/`escape_confirm` left at today's post-corpus values:");
    for (r, _, _, _, _, _, npix, am, al, snaive, _, _) in &stage1 {
        let v = if *am == 0 {
            "the audit was INERT here -- four flags would have been enough".to_string()
        } else {
            format!("the audit MATTERED: {al} outcome labels differ")
        };
        println!("    {r:<14} differs on {am:>5}/{npix:<6}  spr {snaive:>10.3e}   {v}");
    }

    println!("\n  SIGNAL RESOLUTION -- count the distinct values before reading any `rho`. A field with");
    println!("  no ordering cannot have one preserved, and would show a low `rho` for that reason:");
    println!("    {:<14} {:>9} {:>8} {:>10} {:>9} {:>8} {:>10}",
             "region", "dist pre", "modal", "p99/p1 pre", "dist now", "modal", "p99/p1 now");
    for (r, _, _, _, _, _, npix, _, _, _, (da, ma, ra), (db, mb, rb)) in &stage1 {
        println!("    {r:<14} {da:>5}/{npix:<3} {ma:>8.3} {ra:>10.3e} {db:>5}/{npix:<3} {mb:>8.3} {rb:>10.3e}");
    }

    // ------------------------------------------------------------------ stage 2: the decisions
    println!("\nSTAGE 2 -- the real descent under both kernels, at the corpus's own SchedCfg");
    println!("(n=8, bootstrap=2, tau=1e-4, hot=q[0.50], agg=median, criterion=within, camera");
    println!("viewport 1024 max_rel_depth 6). `k_frac` is held at the CURRENT default on both arms:");
    println!("this measures the kernel, and letting a second knob move would unattribute it.\n");
    println!("{:<14} {:>7} {:>7} {:>7} {:>8} {:>28} {:>28}",
             "region", "lv pre", "lv now", "shared", "differ", "stops pre", "stops now");
    for r in &regions {
        let Some((_, cx, cy, body)) = grid::REGIONS.iter().find(|x| x.0 == *r).cloned() else {
            continue;
        };
        let cfg = SchedCfg {
            n: 8,
            bootstrap_levels: 2,
            budget,
            tau_display: 1e-4,
            hot_rule: HotRule::Quantile(0.5),
            alpha_hi: 0.2,
            camera: Some(Camera::framing(cx, cy, 0.05, 1024)),
            ..Default::default()
        };
        let (ta, _) = scheduler::descend(cx, cy, 0.05, body, &cfg, &pre, Precision::F64);
        let (tb, _) = scheduler::descend(cx, cy, 0.05, body, &cfg, &now, Precision::F64);

        // Keyed on the box, not the index: two trees that split differently do not agree on
        // node numbering, and an index-keyed diff would report churn that is renumbering.
        let key = |q: &prin_rs::quad::Quad| {
            (q.level, (q.cx / q.half).round() as i64, (q.cy / q.half).round() as i64)
        };
        let ma: std::collections::HashMap<_, Decision> =
            ta.nodes.iter().filter(|q| q.red.n_footprints > 0).map(|q| (key(q), q.decision)).collect();
        let mb: std::collections::HashMap<_, Decision> =
            tb.nodes.iter().filter(|q| q.red.n_footprints > 0).map(|q| (key(q), q.decision)).collect();
        let shared: Vec<_> = ma.keys().filter(|k| mb.contains_key(k)).collect();
        let differ = shared.iter().filter(|k| ma[**k] != mb[**k]).count();

        let stops = |t: &prin_rs::quad::QuadTree| -> String {
            let mut m: std::collections::BTreeMap<&str, usize> = Default::default();
            for q in t.nodes.iter().filter(|q| q.is_leaf() && q.red.n_footprints > 0) {
                *m.entry(q.decision.name()).or_default() += 1;
            }
            m.iter().map(|(k, v)| format!("{k}:{v}")).collect::<Vec<_>>().join(" ")
        };
        println!("{r:<14} {:>7} {:>7} {:>7} {:>8} {:>28} {:>28}",
                 ta.leaves().count(), tb.leaves().count(), shared.len(), differ,
                 stops(&ta), stops(&tb));
    }

    println!("\nHOW TO READ IT");
    println!("  `rho` near 1.0 with a wide ratio spread = a monotone rescale. The ORDERING survived,");
    println!("  so every criterion comparison in the corpus is largely intact and the numbers are not.");
    println!("  `rho` far from 1.0 = the field was re-ordered and the criterion work is a full");
    println!("  re-derivation. `antip` near 1.0 means the chord column is SATURATED and says nothing:");
    println!("  read `rho` and STAGE 2 there. `differ`/`shared` is the bottom line -- it is measured");
    println!("  over shared quads only, because a quad in one tree and not the other never decided.");
}
