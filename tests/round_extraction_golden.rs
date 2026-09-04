//! **The round extraction moved nothing.**
//!
//! `descend_with`'s loop body became `scheduler::round`, so the one-shot descent and the
//! per-frame step are one implementation rather than two that agree today. A refactor claiming
//! bitwise identity is not tested by the claim, and this project has already reported "reproduces
//! bitwise" from eleven dumps while nineteen had moved — so the check is a hash of the **whole
//! record table** across a matrix of configurations, over the analytic fields in `prin_rs::testing`
//! at **zero trajectories**, which is what lets it be a committed test rather than a manual
//! regeneration.
//!
//! **"Bitwise" needs restating, and that is a finding rather than a weakening.** `tree::write` puts
//! `wall_seconds` in the `.prnq` header, so file bytes were never reproducible and never could be.
//! What is preservable — and what a corpus diff must compare — is the record block plus the named
//! header tokens, with the `decision` column diffed specifically: a tree with the same leaf count
//! and a moved stop-reason column is the documented failure mode of a parameter change and has
//! cost this project a round trip twice.
//!
//! The hashes below were taken **before** the extraction, on `4afef35`. If one moves, the round
//! body is not a move.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use prin_rs::grid::Slice;
use prin_rs::quad::{Agg, Criterion, QuadTree, StructureMode};
use prin_rs::scheduler::{self, Mode, Policy, SchedCfg};

/// Hash every field of every node, in node order, plus the stats that are not a clock.
///
/// Deliberately not a leaf count: a leaf count is exactly what stays put while the `decision`
/// column moves underneath it.
fn fingerprint(tree: &QuadTree, st: &scheduler::SchedStats) -> u64 {
    let mut h = DefaultHasher::new();
    tree.nodes.len().hash(&mut h);
    for q in &tree.nodes {
        q.level.hash(&mut h);
        q.cx.to_bits().hash(&mut h);
        q.cy.to_bits().hash(&mut h);
        q.half.to_bits().hash(&mut h);
        q.parent.hash(&mut h);
        q.children.hash(&mut h);
        q.sib_index.hash(&mut h);
        q.iteration.hash(&mut h);
        q.decision.code().hash(&mut h);
        q.merged.hash(&mut h);
        for v in [q.alpha, q.alpha_mean, q.alpha_p90, q.alpha_sibling_spread, q.alpha_area,
                  q.alpha_spread_set, q.no_gain_weight] {
            v.map(f64::to_bits).hash(&mut h);
        }
        // The reduction, through the one accessor every consumer reads it by.
        q.red.n_footprints.hash(&mut h);
        q.red.n_unresolved.hash(&mut h);
        q.red.n_nonfinite.hash(&mut h);
        q.red.n_distinct_ic.hash(&mut h);
        q.red.total_substeps.hash(&mut h);
        for a in [Agg::Mean, Agg::Median, Agg::P90] {
            q.red.spread(a).to_bits().hash(&mut h);
        }
    }
    // `wall_seconds` is excluded on purpose: it is a clock, and including it would make this test
    // fail for the one reason that says nothing.
    st.quads_computed.hash(&mut h);
    st.footprints.hash(&mut h);
    st.iterations.hash(&mut h);
    st.budget_exhausted.hash(&mut h);
    st.balance_forced.hash(&mut h);
    st.leaves_per_iteration.hash(&mut h);
    h.finish()
}

fn cell(policy: Policy, mode: Mode, k_frac: f64, balance: bool, budget: usize) -> SchedCfg {
    SchedCfg {
        // **`N = 8`, the production value, because the fields and the agreement arm were
        // calibrated there.** At `n = 4` the whole matrix went blind: every structured field
        // floored at level 2 (21 quads), `smooth` ran to a complete 1024, so no tree carried any
        // depth contrast and `balance` moved nothing in 160 cells -- and `filament_through_sea`
        // produced a fingerprint BITWISE EQUAL to `sea`, the filament invisible. A fixture chosen
        // for speed that cannot see its own subject.
        n: 8,
        budget,
        tau_display: 1e-2,
        policy,
        mode,
        k_frac,
        balance,
        criterion: Criterion::Within,
        agg: Agg::Median,
        structure: StructureMode::Off,
        max_level: Some(5),
        camera: None,
        ..Default::default()
    }
}

/// Every field, every policy, every mode, both `k_frac` regimes, both balance arms.
#[test]
fn the_round_extraction_is_a_move() {
    let fields: Vec<(&str, Box<dyn Fn(&Slice, usize) -> prin_rs::ensemble::pixel::PixelOut + Sync>)> = vec![
        ("smooth", Box::new(prin_rs::testing::smooth(1.0, 13.0))),
        // **The feature must not sit on a quad boundary.** At `x0 = 0.0` on a root spanning
        // [-1, 1] the step falls exactly on the midline at every level, so no quad ever straddles
        // it, every footprint agrees, and the field is featureless: the tree stopped at the
        // complete bootstrap (21 quads) and `balance_forced` was 0 in all 160 cells. An offset
        // that is not a dyadic fraction cuts through quads at every depth.
        ("step", Box::new(prin_rs::testing::step(0.137, 13.0))),
        ("sea", Box::new(prin_rs::testing::sea(7, 13.0))),
        ("filament_in_sea", Box::new(prin_rs::testing::filament_in_sea(0.137, 7, 13.0))),
        ("filament_through_sea", Box::new(prin_rs::testing::filament_through_sea(0.137, 7, 13.0))),
    ];

    let mut seen: Vec<String> = Vec::new();
    let mut fps: std::collections::HashMap<(usize, usize, usize, usize, usize, usize), u64> =
        std::collections::HashMap::new();
    let policies = [Policy::Tolerance, Policy::Alpha];
    let modes = [Mode::Balanced, Mode::Uniform];
    let kfs = [0.25f64, 1.0];
    let bals = [false, true];
    let budgets = [2000usize, 40];

    for (fi, (fname, f)) in fields.iter().enumerate() {
        for (pi, &policy) in policies.iter().enumerate() {
            for (mi, &mode) in modes.iter().enumerate() {
                for (ki, &k_frac) in kfs.iter().enumerate() {
                    for (bi, &balance) in bals.iter().enumerate() {
                        for (gi, &budget) in budgets.iter().enumerate() {
                            let cfg = cell(policy, mode, k_frac, balance, budget);
                            let (tree, st) =
                                scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg, 13.0, &**f);
                            let fp = fingerprint(&tree, &st);
                            fps.insert((fi, pi, mi, ki, bi, gi), fp);
                            seen.push(format!(
                                "{fname} {policy:?} {mode:?} k{k_frac} b{balance} B{budget} \
                                 -> {} quads {} leaves f{} [{}] {fp:016x}",
                                st.quads_computed,
                                tree.leaves().count(),
                                st.balance_forced,
                                tree.stop_breakdown()
                            ));
                        }
                    }
                }
            }
        }
    }

    for line in &seen {
        println!("{line}");
    }
    let distinct: std::collections::HashSet<u64> = fps.values().cloned().collect();

    // **The arm that says the matrix exercises anything, per AXIS.** A bare distinct count is the
    // wrong shape: many cells coincide for good reasons — a smooth field keeps everything whatever
    // the policy — and the first cut of this asserted `distinct > cells/4` and fired at 18 of 160
    // on correct code. What matters is that **each knob moves at least one cell**, or that knob is
    // untested and a refactor could break it silently.
    let axis = |name: &str, n: usize, key: &dyn Fn(&(usize, usize, usize, usize, usize, usize), usize)
        -> (usize, usize, usize, usize, usize, usize)| {
        let moved = fps
            .iter()
            .filter(|(k, v)| {
                (0..n).any(|j| fps.get(&key(k, j)).is_some_and(|w| w != *v))
            })
            .count();
        println!("axis {name:<8} moves {moved} of {} cells", fps.len());
        assert!(moved > 0, "axis `{name}` changes nothing -- it is untested by this matrix");
    };
    axis("field", fields.len(), &|k, j| (j, k.1, k.2, k.3, k.4, k.5));
    axis("policy", policies.len(), &|k, j| (k.0, j, k.2, k.3, k.4, k.5));
    axis("mode", modes.len(), &|k, j| (k.0, k.1, j, k.3, k.4, k.5));
    axis("k_frac", kfs.len(), &|k, j| (k.0, k.1, k.2, j, k.4, k.5));
    axis("balance", bals.len(), &|k, j| (k.0, k.1, k.2, k.3, j, k.5));
    axis("budget", budgets.len(), &|k, j| (k.0, k.1, k.2, k.3, k.4, j));

    let mut all = DefaultHasher::new();
    for line in &seen {
        line.hash(&mut all);
    }
    let combined = all.finish();
    println!("\ncombined {combined:016x} over {} cells, {} distinct", seen.len(), distinct.len());

    // Pinned from the pre-extraction build. A move does not change it.
    assert_eq!(
        combined, GOLDEN,
        "the round extraction changed a tree. Diff the per-cell lines above against the \
         pre-extraction run; the `decision` column is where a parameter change hides."
    );
}

/// Taken on `34c44cc`, with `src/scheduler.rs` and `src/camera.rs` at that commit — i.e. **before**
/// `descend_with`'s body moved into `scheduler::round` and before cursor bias existed. Both changes
/// must be moves, and this is what says so.
const GOLDEN: u64 = 0xf7e1_414f_38b6_1261;
