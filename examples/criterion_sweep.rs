//! **The sweep the criterion work never had.** Nothing is overwritten; everything lands in
//! `results/sweep/`.
//!
//! # Why this exists
//!
//! PR #18 built the machinery and then ran nothing with it enabled. All 69 committed dumps carry
//! `tau_display=1e-4  structure=off  k_frac=1  criterion=within` — the pre-fix configuration with
//! new columns attached. Three things reproduce the old behaviour exactly:
//!
//! - **`tau = 1e-4` still gates the split.** The relative hot rule feeds the *shape* statistics,
//!   which is the right design, but the decision never consults it. `tau` sits at the 0.4th
//!   percentile of the observed spread distribution, so the predicate is true for ~99.6% of quads.
//! - **`k_frac = 1` takes the top 100% of the frontier.** The ranking runs and changes nothing.
//! - **`structure = off`, `criterion = within`** — neither new signal is in play.
//!
//! So the before/after the plan promised cannot be read off the corpus: there is no after.
//!
//! **One correction to that diagnosis, because it matters for reading the control row.**
//! `mode=balanced, k_frac=1` is *not* `Mode::Uniform`. Uniform returns `Split` unconditionally
//! and bypasses the `tau` and `alpha` gates entirely; balanced at `k_frac = 1` still applies both
//! and only declines to truncate the frontier. The two produce different trees. What is true is
//! that the **rank truncation** never engaged, which is the mechanism that was supposed to be new.
//!
//! # Sweep, do not pick
//!
//! One setting is not a result. Stage 1 sweeps `tau × k_frac` at `structure=off`, because those
//! two alone decide whether the tree is ever selective. Stage 2 fixes the best pair and sweeps
//! `structure × criterion`.
//!
//! `tau = 1e-4` and `k_frac = 1.0` are retained as **labelled degenerate controls**: the sweep
//! should show the failure it is correcting rather than merely avoiding it.
//!
//! # Three targets, and why each
//!
//! `near-field` is where every prior result was measured. `deep interior` because a change that
//! only improves near-field is tuning. **`preset_shape` because it is the only chart where the
//! camera veto does not bind** — 0% of leaves stopped by `MaxRelDepth` — so it is the only one
//! where the criterion's own decisions are visible rather than a cap's.
//!
//! # Read the stop-reason breakdown, not the leaf count
//!
//! A leaf count says nothing without knowing what stopped it. Outside `preset_shape` most leaf
//! counts are largely a fact about `MaxRelDepth`.

use std::io::BufWriter;

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart};
use prin_rs::output::tree as treeout;
use prin_rs::quad::{Criterion, Decision, StructureMode};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Mode, SchedCfg};

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
            chart: Chart::BodyPlane,
            cx,
            cy,
            half: 0.05,
            body,
        });
    }
    let ps = Chart::preset_shape();
    v.push(Target { name: "preset_shape", chart: ps, cx: 0.0, cy: 0.0, half: ps.default_half(), body: 0 });
    v
}

/// A filename that carries its own settings, so a directory listing is a settings table.
fn stem(t: &Target, tau: f64, k: f64, st: StructureMode, cr: Criterion) -> String {
    format!(
        "results/sweep/{}__tau{:.0e}__k{:.2}__struct-{}__crit-{}",
        t.name, tau, k, st.name(), cr.name()
    )
}

#[allow(clippy::too_many_arguments)]
fn run(
    t: &Target,
    tau: f64,
    k: f64,
    st: StructureMode,
    cr: Criterion,
    alpha_hi: f64,
    budget: usize,
    res: usize,
    ens: &EnsembleCfg,
    write: bool,
) -> (usize, f64, f64, usize, [usize; 11], f64) {
    let cfg = SchedCfg {
        budget,
        tau_display: tau,
        alpha_hi,
        alpha_lo: alpha_hi,
        criterion: cr,
        structure: st,
        mode: Mode::Balanced,
        k_frac: k,
        camera: Some(Camera::framing(t.cx, t.cy, t.half, res)),
        chart: t.chart,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    let (tree, sst) = scheduler::descend(t.cx, t.cy, t.half, t.body, &cfg, ens, Precision::F64);
    let wall = t0.elapsed().as_secs_f64();

    let leaves: Vec<usize> = tree.leaves().collect();
    let lv: Vec<f64> = leaves.iter().map(|&i| tree.nodes[i].level as f64).collect();
    let m = lv.iter().sum::<f64>() / lv.len().max(1) as f64;
    let var = lv.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / lv.len().max(1) as f64;
    let mx = lv.iter().cloned().fold(f64::MIN, f64::max);
    let pmax = lv.iter().filter(|&&x| x == mx).count() as f64 / lv.len().max(1) as f64;
    let mut levels: Vec<i64> = lv.iter().map(|&x| x as i64).collect();
    levels.sort_unstable();
    levels.dedup();

    let mut dec = [0usize; 11];
    for &i in &leaves {
        let c = tree.nodes[i].decision.code() as usize;
        if c < 11 {
            dec[c] += 1;
        }
    }
    // `Split` never lands on a leaf -- deciding to split makes it internal -- so it is counted
    // over all quads or the column reads zero and looks like the criterion never fired.
    dec[Decision::Split.code() as usize] =
        tree.nodes.iter().filter(|q| q.decision == Decision::Split).count();

    if write {
        let _ = std::fs::create_dir_all("results/sweep");
        if let Ok(f) = std::fs::File::create(format!("{}.prnq", stem(t, tau, k, st, cr))) {
            let mut w = BufWriter::new(f);
            let _ = treeout::write(&mut w, &tree, &cfg, ens, &sst, t.name, "f64");
        }
    }
    (leaves.len(), 100.0 * pmax, var, levels.len(), dec, wall)
}

fn header() {
    println!("{:>14} {:>8} {:>6} {:>8} {:>6} {:>7} {:>5} {:>7} {:>6} {:>6} {:>7} {:>7} {:>7}",
             "target", "tau", "k", "struct", "crit", "leaves", "lvls", "depthvar", "%max",
             "split", "floor", "keep", "veto");
}

#[allow(clippy::too_many_arguments)]
fn row(t: &Target, tau: f64, k: f64, st: StructureMode, cr: Criterion,
       r: &(usize, f64, f64, usize, [usize; 11], f64)) {
    let (leaves, pmax, var, lvls, dec, _) = r;
    let veto = dec[Decision::ScreenFloor.code() as usize] + dec[Decision::MaxRelDepth.code() as usize];
    println!("{:>14} {:>8.0e} {:>6.2} {:>8} {:>6} {:>7} {:>5} {:>8.3} {:>5.0}% {:>6} {:>7} {:>7} {:>7}",
             t.name, tau, k, st.name(),
             match cr { Criterion::Within => "within", Criterion::FracHotBetween => "fhb",
                        Criterion::LayoutRel => "lrel", Criterion::GradRms => "grms", _ => cr.name() },
             leaves, lvls, var, pmax,
             dec[Decision::Split.code() as usize],
             dec[Decision::Floor.code() as usize],
             dec[Decision::Keep.code() as usize], veto);
}

fn main() {
    let stage: usize = arg(1, 0);
    let budget: usize = arg(2, 40000);
    let alpha_hi: f64 = arg(3, 0.2);
    let res: usize = arg(4, 1024);
    let ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let tg = targets();

    const TAUS: [f64; 5] = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2];
    const KS: [f64; 4] = [1.0, 0.5, 0.25, 0.1];

    println!("budget {budget}, alpha_hi={alpha_hi}, N=8, E+1={}, t={}, f64, viewport {res}^2",
             ens.n_extra + 1, ens.t_max);
    println!("tau = 1e-4 and k_frac = 1.0 are the LABELLED DEGENERATE CONTROLS: they reproduce");
    println!("the configuration every committed dump was made with.\n");

    if stage == 0 {
        // ---- the wiring check ------------------------------------------------------------
        //
        // A sweep over a parameter that is not plumbed through produces identical trees at every
        // setting and reads as "the criterion cannot be fixed" -- which is exactly the wrong
        // conclusion to draw from a wiring bug. Two settings per knob, diffed, BEFORE the run.
        println!("=== WIRING CHECK: do the knobs reach the tree at all? ===\n");
        header();
        let t = &tg[0];
        let base = run(t, 1e-4, 1.0, StructureMode::Off, Criterion::Within, alpha_hi, budget, res, &ens, false);
        row(t, 1e-4, 1.0, StructureMode::Off, Criterion::Within, &base);
        let hi_tau = run(t, 3e-3, 1.0, StructureMode::Off, Criterion::Within, alpha_hi, budget, res, &ens, false);
        row(t, 3e-3, 1.0, StructureMode::Off, Criterion::Within, &hi_tau);
        let lo_k = run(t, 1e-4, 0.25, StructureMode::Off, Criterion::Within, alpha_hi, budget, res, &ens, false);
        row(t, 1e-4, 0.25, StructureMode::Off, Criterion::Within, &lo_k);
        let mult = run(t, 1e-4, 1.0, StructureMode::Multiply, Criterion::Within, alpha_hi, budget, res, &ens, false);
        row(t, 1e-4, 1.0, StructureMode::Multiply, Criterion::Within, &mult);
        let fhb = run(t, 1e-4, 1.0, StructureMode::Off, Criterion::FracHotBetween, alpha_hi, budget, res, &ens, false);
        row(t, 1e-4, 1.0, StructureMode::Off, Criterion::FracHotBetween, &fhb);

        println!();
        let mut ok = true;
        for (name, r) in [("tau", &hi_tau), ("k_frac", &lo_k), ("structure", &mult), ("criterion", &fhb)] {
            let same = r.0 == base.0 && r.2 == base.2;
            println!("  {name:>10}: {} the baseline tree ({} leaves against {})",
                     if same { "**DOES NOT CHANGE**" } else { "changes" }, r.0, base.0);
            ok &= !same;
        }
        println!();
        if ok {
            println!("All four knobs reach the tree. Run stage 1.");
        } else {
            println!("AT LEAST ONE KNOB IS NOT PLUMBED THROUGH. Do not run the sweep: identical");
            println!("trees at every setting would read as `the criterion cannot be fixed`, which");
            println!("is the wrong conclusion to draw from a wiring bug.");
        }
        return;
    }

    if stage == 1 {
        println!("=== STAGE 1: tau x k_frac at structure=off, criterion=within ===\n");
        header();
        for t in &tg {
            for &tau in &TAUS {
                for &k in &KS {
                    let r = run(t, tau, k, StructureMode::Off, Criterion::Within,
                                alpha_hi, budget, res, &ens, true);
                    row(t, tau, k, StructureMode::Off, Criterion::Within, &r);
                }
            }
            println!();
        }
        println!("DEPTH VARIANCE is the headline. If no configuration produces a tree with");
        println!("meaningful depth variance, that is the finding and it is bigger than any");
        println!("individual setting.");
        return;
    }

    if stage == 3 {
        // ---- stage 3: alpha, because stage 1 says it is the only knob that binds ---------
        //
        // **Not in the brief, and stage 1 is why.** `tau x k_frac` is inert on two of the three
        // targets: `deep interior` returns 22/3/0.614 at every `k` and `preset_shape` returns
        // 16/1/0.000 in all twenty cells. Neither knob can reach them, and the reason is
        // structural rather than a matter of degree:
        //
        // - **`tau` cannot gate `preset_shape`.** Its leaf-spread median is `2.86e-1`, 3400x
        //   above the largest `tau` swept, so every quad clears the spread gate at every setting.
        // - **`k_frac` has nothing to rank.** It truncates the set that already decided to
        //   *split*; `preset_shape` produces **zero** splits past the bootstrap, and `deep
        //   interior`'s frontier is 1-2 quads a round, where `ceil(1 * 0.1) = 1` truncates
        //   nothing.
        //
        // What decides both is `alpha`: `preset_shape`'s sixteen level-2 quads are 8 `Floor`
        // (alpha below `alpha_lo`) and 8 `Keep` (no computable alpha at all). So the sweep as
        // specified could not have moved them, and running only it would have reported "no
        // configuration produces depth variance" while leaving the binding constraint untouched.
        let tau: f64 = arg(5, 1e-4);
        let k: f64 = arg(6, 0.25);
        println!("=== STAGE 3: alpha_hi at tau={tau:.0e}, k_frac={k}, structure=off ===");
        println!("Stage 1 showed tau and k_frac are inert on deep_interior and preset_shape.");
        println!("alpha is what binds them, and it is the one knob the brief did not sweep.\n");
        header();
        for t in &tg {
            for &a in &[0.5f64, 0.2, 0.1, 0.05, 0.0, -1.0] {
                let r = run(t, tau, k, StructureMode::Off, Criterion::Within, a, budget, res,
                            &ens, true);
                println!("{:>14} {:>8.0e} {:>6.2} {:>8} {:>6} {:>7} {:>5} {:>8.3} {:>5.0}%                           {:>6} {:>7} {:>7} {:>7}   alpha_hi={a}",
                         t.name, tau, k, "off", "within", r.0, r.3, r.2, r.1,
                         r.4[1], r.4[2], r.4[3], r.4[7] + r.4[8]);
            }
            println!();
        }
        println!("alpha_hi = -1.0 is the DEGENERATE CONTROL: every finite alpha clears it, so the");
        println!("alpha gate is effectively off and the tree runs to the veto or the budget. If a");
        println!("target is still flat there, nothing in the criterion can move it.");
        return;
    }

    // ---- stage 2: structure x criterion at a fixed (tau, k) --------------------------
    let tau: f64 = arg(5, 1e-4);
    let k: f64 = arg(6, 0.25);
    println!("=== STAGE 2: structure x criterion at tau={tau:.0e}, k_frac={k} ===\n");
    header();
    for t in &tg {
        for st in StructureMode::ALL {
            for cr in [Criterion::Within, Criterion::FracHotBetween, Criterion::LayoutRel,
                       Criterion::GradRms] {
                let r = run(t, tau, k, st, cr, alpha_hi, budget, res, &ens, true);
                row(t, tau, k, st, cr, &r);
            }
        }
        println!();
    }
}
