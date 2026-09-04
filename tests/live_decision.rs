//! Phase 2 of the refinement rebuild: the tolerance policy, pinned against the legacy one.
//!
//! Every test here carries the arm that shows it can fire. The first is the pin: under
//! `Policy::Alpha` the near-field fixture must reproduce the tree the shipped code produced
//! before the frontier change, level histogram and stop breakdown alike, so the legacy policy
//! stays a faithful reconstruction of every committed tree.

use prin_rs::camera::Camera;
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Policy, SchedCfg};

/// The `tests/colour.rs` fixture: near-field, `N = 4`, budget 400, `tau = 1e-4`,
/// `alpha_hi = alpha_lo = 0.2`, a 64-pixel camera, `t = 13`.
fn alpha_fixture() -> (prin_rs::quad::QuadTree, scheduler::SchedStats) {
    let ens = EnsembleCfg { t_max: 13.0, refine_flagged: false, ..Default::default() };
    let cam = Camera::framing(1.0, 3.0, 0.05, 64);
    let cfg = SchedCfg {
        n: 4,
        budget: 400,
        tau_display: 1e-4,
        alpha_hi: 0.2,
        alpha_lo: 0.2,
        policy: Policy::Alpha,
        camera: Some(cam),
        ..Default::default()
    };
    scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64)
}

fn histogram(t: &prin_rs::quad::QuadTree) -> Vec<(u32, usize)> {
    let mut m = std::collections::BTreeMap::new();
    for i in t.leaves() {
        *m.entry(t.nodes[i].level).or_insert(0usize) += 1;
    }
    m.into_iter().collect()
}

/// A hash of the leaf geometry and decisions, order-independent: FNV over
/// `(level, cx bits, cy bits, decision code)` sorted. Pins the tree, not only its counts.
fn tree_hash(t: &prin_rs::quad::QuadTree) -> u64 {
    let mut rows: Vec<(u32, u64, u64, u8)> = t
        .leaves()
        .map(|i| {
            let q = &t.nodes[i];
            (q.level, q.cx.to_bits(), q.cy.to_bits(), q.decision.code())
        })
        .collect();
    rows.sort_unstable();
    let mut h: u64 = 0xcbf29ce484222325;
    for (l, x, y, d) in rows {
        for b in l.to_le_bytes().iter().chain(x.to_le_bytes().iter()).chain(y.to_le_bytes().iter()).chain([d].iter()) {
            h ^= *b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}

/// **T10 — the legacy policy reproduces the pre-change tree.** Recorded from the shipped code
/// at `042fb73` before the frontier and policy changes: 40 leaves, levels `2:12 3:12 4:16`, and
/// the stop breakdown printed below. A change in either means `Policy::Alpha` is no longer the
/// thing every committed `.prnq` was cut with.
#[test]
fn the_legacy_alpha_policy_reproduces_the_pre_change_tree() {
    let (t, st) = alpha_fixture();
    let h = histogram(&t);
    let stop = t.stop_breakdown();
    let hash = tree_hash(&t);
    println!("alpha fixture: {} quads, {} leaves, levels {h:?}, stop [{stop}], hash {hash:#018x}",
             st.quads_computed, t.leaves().count());
    assert_eq!(st.quads_computed, 53);
    assert_eq!(h, vec![(2, 12), (3, 12), (4, 16)], "level histogram moved under Policy::Alpha");
    assert_eq!(t.leaves().count(), 40);
    assert_eq!(stop, GOLDEN_STOP, "stop breakdown moved under Policy::Alpha");
    if GOLDEN_HASH != 0 {
        assert_eq!(hash, GOLDEN_HASH, "leaf geometry or decisions moved under Policy::Alpha");
    }
}

/// Recorded at `042fb73`, before the policy change; see the test above.
const GOLDEN_STOP: &str = "keep:24 screen_floor:16";
const GOLDEN_HASH: u64 = 0x0df085b5c7922bd3;

// ---------------------------------------------------------------------------------------
// The tolerance policy on analytic fields, where the right tree is known in advance.
// ---------------------------------------------------------------------------------------

use prin_rs::quad::{Decision, QuadTree};
use prin_rs::scheduler::Mode;

const T: f64 = 13.0;

fn cfg(max_level: u32, budget: usize) -> SchedCfg {
    SchedCfg {
        n: 8,
        budget,
        tau_display: 0.01,
        policy: Policy::Tolerance,
        max_level: Some(max_level),
        camera: None,
        k_frac: prin_rs::scheduler::K_FRAC_UNRANKED,
        ..Default::default()
    }
}

fn count(t: &QuadTree, d: Decision) -> usize {
    t.leaves().filter(|&i| t.nodes[i].decision == d).count()
}

/// **T1 — a smooth field is `Keep` everywhere above the bootstrap.** The control: the same
/// field at a tolerance below its own cell spread splits to the cap.
#[test]
fn a_smooth_field_is_resolved_at_the_bootstrap_and_splits_only_when_the_tolerance_is_below_it() {
    let field = prin_rs::testing::smooth(0.2, T);
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(4, 400), T, &field);
    println!("smooth: {} quads, {} leaves, stop [{}]", st.quads_computed, t.leaves().count(), t.stop_breakdown());
    assert_eq!(t.leaves().count(), 16, "a smooth field must stop at the bootstrap");
    assert_eq!(count(&t, Decision::Keep), 16);
    assert_eq!(count(&t, Decision::Floor), 0);
    assert_eq!(count(&t, Decision::Split), 0);

    let tight = SchedCfg { tau_display: 1e-4, ..cfg(4, 400) };
    let (t2, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &tight, T, &field);
    assert_eq!(t2.leaves().count(), 256, "below its own spread the field must split to the cap");
    assert_eq!(count(&t2, Decision::MaxLevel), 256);
}

/// **T2 — a step is refined along the step and nowhere else.** Leaves at the cap all straddle
/// the step; every leaf off it is `Keep` at the bootstrap; the count is `O(2^L)` against `4^L`.
/// Controls: `Mode::Uniform` fills the cap, and `Policy::Alpha` on the same field never reaches
/// it — the one hot column of eight footprints in sixty-four is outvoted by the median, which is
/// the defect this policy replaces.
#[test]
fn a_step_is_refined_along_the_step_and_nowhere_else() {
    let x0 = 0.123;
    let field = prin_rs::testing::step(x0, T);
    let levels = 5u32;
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(levels, 2000), T, &field);
    let leaves: Vec<usize> = t.leaves().collect();
    println!("step: {} quads, {} leaves, stop [{}]", st.quads_computed, leaves.len(), t.stop_breakdown());
    assert!(leaves.len() < 4usize.pow(levels - 1), "{} leaves is not O(2^L)", leaves.len());
    // A footprint's copies span half a cell past the quad edge (`Slice::axis` puts a sample ON
    // the edge), so the quad beside the step can straddle it through its edge column. The
    // geometric claim is therefore: an unresolved quad is within half a footprint cell of the
    // step, a resolved one is `Keep` wherever it sits, and only unresolved quads reach the cap.
    for &i in &leaves {
        let q = &t.nodes[i];
        let hx = 2.0 * q.half / (t.n - 1) as f64;
        let near = (x0 - q.cx).abs() <= q.half + 0.5 * hx;
        match q.decision {
            Decision::MaxLevel => {
                assert_eq!(q.level, levels);
                assert!(q.red.n_unresolved > 0, "a cap leaf with nothing unresolved");
                assert!(near, "a cap leaf far from the step: cx {} half {}", q.cx, q.half);
            }
            Decision::Keep => {
                assert_eq!(q.red.n_unresolved, 0, "a Keep leaf with something unresolved");
                assert!(q.level >= 2);
            }
            d => panic!("unexpected stop {d:?} on a step field"),
        }
    }
    let at_cap = count(&t, Decision::MaxLevel);
    assert!(at_cap > 0, "nothing reached the cap");
    assert!(at_cap <= 4 * (1usize << levels), "{at_cap} cap leaves is more than four columns");
    assert_eq!(count(&t, Decision::Floor), 0);

    let uni = SchedCfg { mode: Mode::Uniform, ..cfg(levels, 2000) };
    let (tu, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &uni, T, &field);
    assert_eq!(tu.leaves().count(), 4usize.pow(levels), "uniform mode must fill the cap");

    let alpha = SchedCfg { policy: Policy::Alpha, alpha_hi: 0.2, alpha_lo: 0.2, ..cfg(levels, 2000) };
    let (ta, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &alpha, T, &field);
    let deepest = ta.leaves().map(|i| ta.nodes[i].level).max().unwrap();
    println!("step under Policy::Alpha: {} leaves, deepest level {deepest}, stop [{}]",
             ta.leaves().count(), ta.stop_breakdown());
    assert!(deepest < levels, "the legacy policy reached the cap on a step, which it never did");
}

/// **T3 — a sea refines uniformly to the cap** under the tolerance policy without a stationarity
/// stop: the sea cost, made visible. The control is `Mode::Uniform` giving the identical tree.
#[test]
fn a_sea_refines_uniformly_to_the_cap_when_full_depth_is_allowed() {
    let field = prin_rs::testing::sea(7, T);
    // The stop OFF and the area floor OFF (`alpha_lo = 0`, the opt-in): the sea cost, visible.
    let off = SchedCfg { stationary: false, alpha_lo: 0.0, ..cfg(4, 400) };
    let (t, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &off, T, &field);
    assert_eq!(t.leaves().count(), 256);
    assert_eq!(count(&t, Decision::MaxLevel), 256);
    let uni = SchedCfg { mode: Mode::Uniform, ..cfg(4, 400) };
    let (tu, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &uni, T, &field);
    let a: Vec<(u32, u64, u64)> = t.leaves().map(|i| (t.nodes[i].level, t.nodes[i].cx.to_bits(), t.nodes[i].cy.to_bits())).collect();
    let b: Vec<(u32, u64, u64)> = tu.leaves().map(|i| (tu.nodes[i].level, tu.nodes[i].cx.to_bits(), tu.nodes[i].cy.to_bits())).collect();
    assert_eq!(a, b, "on a sea the tolerance policy and uniform mode must build the same tree");
}

/// **T5 — `Deferred` and `Keep` are exclusive on what they say about the quad.** Under a ranked
/// frontier a `Deferred` leaf is always unresolved and a `Keep` leaf never is, and both
/// populations exist.
#[test]
fn deferred_is_never_resolved_and_keep_is_never_unresolved() {
    let field = prin_rs::testing::step(0.123, T);
    let ranked = SchedCfg { k_frac: 0.25, budget: 120, ..cfg(6, 120) };
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &ranked, T, &field);
    let (mut n_def, mut n_keep) = (0, 0);
    for i in t.leaves() {
        let q = &t.nodes[i];
        match q.decision {
            Decision::Deferred => {
                n_def += 1;
                assert!(q.red.n_unresolved > 0, "a Deferred leaf with nothing unresolved");
            }
            Decision::Keep => {
                n_keep += 1;
                assert_eq!(q.red.n_unresolved, 0, "a Keep leaf with something unresolved");
            }
            _ => {}
        }
    }
    println!("ranked step: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    assert!(n_def > 0 && n_keep > 0, "both populations must exist: deferred {n_def}, keep {n_keep}");
}

/// **T6 — leaf depth follows the between-footprint variation, not the other way round.** On the
/// step, Spearman(depth, between_shape) over the leaves is positive; under `Policy::Alpha` the
/// tree has one level and the correlation is undefined, which is the control.
#[test]
fn leaf_depth_tracks_between_shape_under_the_tolerance_policy() {
    let field = prin_rs::testing::step(0.123, T);
    let (t, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(5, 2000), T, &field);
    let (xs, ys): (Vec<f64>, Vec<f64>) = t
        .leaves()
        .map(|i| (t.nodes[i].level as f64, t.nodes[i].red.between_shape))
        .unzip();
    let rho = prin_rs::stats::spearman(&xs, &ys);
    println!("spearman(depth, between_shape) = {rho:+.4} over {} leaves", xs.len());
    assert!(rho > 0.0, "depth does not follow between_shape: rho {rho}");
    // The per-level medians, not only the pooled correlation (the standing rule): the median
    // between-footprint variation must not fall with depth.
    let mut by_level: std::collections::BTreeMap<u32, Vec<f64>> = Default::default();
    for (x, y) in xs.iter().zip(ys.iter()) {
        by_level.entry(*x as u32).or_default().push(*y);
    }
    let meds: Vec<(u32, f64)> = by_level
        .into_iter()
        .map(|(l, mut v)| {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            (l, v[v.len() / 2])
        })
        .collect();
    println!("per-level median between_shape: {meds:?}");
    assert!(meds.len() >= 2, "one level only; nothing to compare");
    assert!(meds.last().unwrap().1 >= meds.first().unwrap().1,
            "median between_shape falls with depth: {meds:?}");

    let alpha = SchedCfg { policy: Policy::Alpha, alpha_hi: 0.2, alpha_lo: 0.2, ..cfg(5, 2000) };
    let (ta, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &alpha, T, &field);
    let (xa, ya): (Vec<f64>, Vec<f64>) = ta
        .leaves()
        .map(|i| (ta.nodes[i].level as f64, ta.nodes[i].red.between_shape))
        .unzip();
    let ra = prin_rs::stats::spearman(&xa, &ya);
    println!("under Policy::Alpha: rho = {ra:+.4} over {} leaves", xa.len());
    assert!(!(ra > 0.0), "the control shows the same correlation, so the test cannot discriminate");
}

// ---------------------------------------------------------------------------------------
// Phase 2b: the stationarity stop.
// ---------------------------------------------------------------------------------------

/// **T3′ — a sea stops as `Stationary` at the bootstrap when the stop is on**, and refines to
/// the cap when it is off (T3 is the control). The sea is white at the footprint scale, its
/// mixture is the same at every scale, and its spread does not fall: all three arms pass.
#[test]
fn a_sea_stops_as_stationary_when_the_stop_is_on() {
    let field = prin_rs::testing::sea(7, T);
    let on = SchedCfg { stationary: true, ..cfg(4, 400) };
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &on, T, &field);
    println!("sea, stationary on: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    let cohs: Vec<f64> = t.leaves().map(|i| t.nodes[i].red.coherence()).collect();
    let tvq: Vec<f64> = t.leaves().map(|i| t.nodes[i].red.mix_tv_quadrants).collect();
    let tvp: Vec<f64> = t.leaves().map(|i| t.nodes[i].red.mix_tv_parent).collect();
    println!("  coherence max {:.3}, quadrant tv max {:.3}, parent tv max {:.3}",
             cohs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
             tvq.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
             tvp.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    assert_eq!(count(&t, Decision::Stationary), 16, "the sea must stop at the bootstrap as Stationary");
    assert_eq!(t.leaves().count(), 16);
    let off = SchedCfg { stationary: false, alpha_lo: 0.0, ..cfg(4, 400) };
    let (t2, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &off, T, &field);
    assert_eq!(t2.leaves().count(), 256, "with the stop off the sea must refine to the cap");
}

/// **T4 — a filament through a sea: the sea stops, the filament refines, the basin keeps.**
/// The arm that shows the stop can fire wrongly: the quads on the filament must read coherent
/// or mixed (they fail arm 1 or arm 2), and none of them may be `Stationary`.
#[test]
fn a_filament_through_a_sea_is_refined_while_the_sea_stops() {
    let x0 = 0.123;
    let field = prin_rs::testing::filament_in_sea(x0, 11, T);
    let levels = 5u32;
    let on = SchedCfg { stationary: true, ..cfg(levels, 2000) };
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &on, T, &field);
    println!("filament in sea: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    let (mut n_stat, mut n_cap, mut n_keep) = (0, 0, 0);
    for i in t.leaves() {
        let q = &t.nodes[i];
        let hx = 2.0 * q.half / (t.n - 1) as f64;
        let on_filament = (x0 - q.cx).abs() <= q.half + 0.5 * hx;
        let in_sea = q.cx + q.half < x0;
        match q.decision {
            Decision::Stationary | Decision::Floor => {
                // Floor: the area floor reads a sea quad with no coherent population as noise
                // and stops it too, where the mixture arm's sampling noise let it past the stop.
                n_stat += 1;
                assert!(in_sea && !on_filament, "a {:?} leaf off the sea: cx {} half {}", q.decision, q.cx, q.half);
            }
            Decision::MaxLevel => {
                n_cap += 1;
                assert!(on_filament, "a cap leaf off the filament: cx {} half {}", q.cx, q.half);
            }
            Decision::Keep => {
                n_keep += 1;
                assert!(q.red.n_unresolved == 0);
            }
            d => panic!("unexpected stop {d:?}"),
        }
    }
    println!("  stationary {n_stat}, at cap {n_cap}, keep {n_keep}");
    assert!(n_stat > 0, "the sea never read as stationary");
    assert!(n_cap > 0, "the filament never reached the cap");
    assert!(n_keep > 0, "the basin never read as resolved");
    // The filament quads at the cap fail the stationarity arms by measurement, not by luck.
    for i in t.leaves() {
        let q = &t.nodes[i];
        if q.decision == Decision::MaxLevel {
            let coh = q.red.coherence();
            let mixed = q.red.mix_tv_quadrants;
            assert!(coh.is_nan() || coh >= 0.3 || mixed >= 0.25,
                    "a filament quad passes the stationarity arms: coh {coh:.3} tv {mixed:.3}");
        }
    }
}

// ---------------------------------------------------------------------------------------
// Phase 2c: the live descent.
// ---------------------------------------------------------------------------------------

/// **T12 — the projection cannot read the future.** A footprint that collides at `t = 5`,
/// viewed at a boundary before that, is bounded at that boundary with every accumulator that
/// reads past it masked; viewed after, it has collided. The control is the terminal view.
#[test]
fn the_live_view_masks_everything_after_the_boundary() {
    use prin_rs::ensemble::pixel::PixelOut;
    use prin_rs::outcome::State;
    let mut p = PixelOut::default();
    p.state = State::Collision as u8;
    p.outcome = (State::Collision as u8) << 2;
    p.t_end = 5.0;
    p.censored = false;
    p.error_ratio = 3.0;
    p.energy_drift_max = 1e-6;
    p.running_max_divergence = 0.4;
    p.first_divergence_t = 4.0;
    p.total_substeps = 1000;
    p.live_t = vec![3.0, 6.0, 13.0];
    p.live_spread_shape = vec![0.001, 0.2, 0.3];
    p.live_spread_event = vec![0.0, 0.143, 0.143];
    p.live_shape = vec![[1.0, 0.0, 0.0]; 3];
    p.live_class = vec![0, 9, 9];

    let early = scheduler::project_at(&p, 0);
    assert_eq!(early.state, State::Bounded as u8, "before it collided it was bounded");
    assert!(early.censored);
    assert_eq!(early.t_end, 3.0);
    assert_eq!(early.spread_shape, 0.001);
    assert_eq!(early.event_class, 0);
    assert!(early.error_ratio.is_nan() && early.energy_drift_max.is_nan());
    assert!(early.running_max_divergence.is_nan() && early.first_divergence_t.is_nan());
    assert!(early.total_substeps < p.total_substeps);

    let late = scheduler::project_at(&p, 1);
    assert_eq!(late.state, State::Collision as u8, "by t = 6 it had collided");
    assert_eq!(late.t_end, 5.0);
    assert!(!late.censored);
    assert_eq!(late.event_class, 9);

    let terminal = scheduler::project_at(&p, 2);
    assert_eq!(terminal.spread_shape, 0.3);
    assert_eq!(terminal.total_substeps, p.total_substeps);
}

/// **T11 — the live tree only grows, and it is the running union of the static trees.**
/// On the pulse the unresolved band rises to mid-march and collapses by the horizon. The live
/// leaf set at every boundary refines the one before; the static tree at the horizon sees only
/// the step, and every one of its leaves is a live leaf or has live leaves beneath it, with
/// at least one strictly refined — the band the static tree never saw.
#[test]
fn the_live_tree_only_grows_and_is_the_running_union_with_merging_off() {
    let field = prin_rs::testing::pulse(0.123, 0.35, 8, T);
    let levels = 4u32;
    let no_merge = SchedCfg { merge: false, ..cfg(levels, 2000) };
    let (t, st) = scheduler::descend_live_with(0.0, 0.0, 1.0, 0, &no_merge, T, &field);
    println!("live pulse: {} quads, {} boundaries, catchup {} substeps, stop [{}]",
             st.quads_computed, st.live.len(), st.catchup_substeps, t.stop_breakdown());
    for p in &st.live {
        println!("  j={} t={:.2} computed={} leaves={} split={} keep={} stationary={} deferred={}",
                 p.j, p.t, p.computed, p.leaves, p.split, p.keep, p.stationary, p.deferred);
    }
    assert_eq!(st.live_leaves.len(), 8);
    // Refinement: every leaf at j+1 is a leaf at j or descends from one.
    let is_desc = |t: &QuadTree, mut i: usize, anc: usize| -> bool {
        loop {
            if i == anc { return true; }
            match t.nodes[i].parent { Some(p) => i = p, None => return false }
        }
    };
    for j in 1..st.live_leaves.len() {
        let prev: std::collections::HashSet<usize> = st.live_leaves[j - 1].iter().cloned().collect();
        for &l in &st.live_leaves[j] {
            assert!(prev.iter().any(|&a| is_desc(&t, l, a)),
                    "boundary {j}: leaf {l} does not refine the previous leaf set");
        }
        assert!(st.live[j].computed >= st.live[j - 1].computed);
    }
    let grew = st.live.windows(2).any(|w| w[1].computed > w[0].computed);
    assert!(grew, "the live tree never grew past the bootstrap");

    // The control: the static tree at the horizon is coarser, and nowhere finer.
    let (ts, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(levels, 2000), T, &field);
    let live_leaves: Vec<usize> = t.leaves().collect();
    let mut strictly_refined = 0usize;
    for l in ts.leaves() {
        let q = &ts.nodes[l];
        // Find the live leaves inside this static leaf's box.
        let inside: Vec<usize> = live_leaves
            .iter()
            .cloned()
            .filter(|&i| {
                let r = &t.nodes[i];
                (r.cx - q.cx).abs() <= q.half - r.half + 1e-12 && (r.cy - q.cy).abs() <= q.half - r.half + 1e-12
            })
            .collect();
        assert!(!inside.is_empty(), "a static leaf has no live leaf under it");
        if inside.len() > 1 {
            strictly_refined += 1;
        }
        assert!(inside.iter().all(|&i| t.nodes[i].level >= q.level), "a live leaf is coarser than the static one");
    }
    println!("static leaves {}, live leaves {}, static leaves strictly refined live: {strictly_refined}",
             ts.leaves().count(), live_leaves.len());
    assert!(strictly_refined > 0, "the live tree should carry the band the static tree collapsed");
    assert!(st.catchup_substeps > 0, "a late split must cost catch-up");
}

// ---------------------------------------------------------------------------------------
// Phase 2d: the area floor and merging.
// ---------------------------------------------------------------------------------------

/// **T14 — the area floor stops a sea at the bootstrap and never a step.** A sea is unresolved
/// over its whole area at every level, so every split reads `alpha_area = 0` and the level-2
/// children are `Floor` under the default `alpha_lo`; T3 is the control, where `alpha_lo = 0`
/// runs the same field to the cap. A step's unresolved area halves per level -- the edge
/// columns weigh half, so a step on a quad edge is not counted twice -- and it reads
/// `alpha_area = 1` at every split.
#[test]
fn the_area_floor_stops_a_sea_and_never_a_step() {
    let sea = prin_rs::testing::sea(7, T);
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(4, 400), T, &sea);
    println!("sea under the floor: {} quads, {} leaves, stop [{}]", st.quads_computed, t.leaves().count(), t.stop_breakdown());
    assert_eq!(t.leaves().count(), 16);
    assert_eq!(count(&t, Decision::Floor), 16);
    let exps: Vec<f64> = t.nodes.iter().filter_map(|q| q.alpha_area).collect();
    assert_eq!(exps.len(), 5, "the root and its four children are the splits that happened");
    for a in &exps {
        assert!(a.abs() < 1e-9, "a sea split read alpha_area {a}");
    }

    let step = prin_rs::testing::step(0.123, T);
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(5, 2000), T, &step);
    println!("step under the floor: {} quads, {} leaves, stop [{}]", st.quads_computed, t.leaves().count(), t.stop_breakdown());
    assert_eq!(count(&t, Decision::Floor), 0);
    assert!(count(&t, Decision::MaxLevel) > 0, "the step never reached the cap");
    let exps: Vec<f64> = t.nodes.iter().filter_map(|q| q.alpha_area).collect();
    assert!(!exps.is_empty());
    let lo = exps.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = exps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("step alpha_area over {} splits: min {lo:.4} max {hi:.4}", exps.len());
    // Over two levels the exponent skips the parent's own grid; the straddles at the coarse
    // and fine ends can each move it by half a level, and the sea's 0 is far below either.
    assert!(lo >= 0.5 && hi <= 1.5, "a line's unresolved area halves per level; read {lo}..{hi}");
}

/// **T15 -- a filament THROUGH a sea is invisible to the area floor and found by the agreement
/// arm; a shore is found by the floor alone.** With sea on both sides every footprint is
/// unresolved, so no split buys less unresolved area and by area alone the level-2 children
/// floor, filament and sea alike -- the limit, stated. With the stationarity stop on the sea
/// children read white and count for nothing, the filament children carry the whole remaining
/// area, and the column refines to the cap while the sea stops. `filament_in_sea` is a shore:
/// sea on one side, basin on the other, and the sea's edge is a line the area floor follows
/// by itself -- the first cut of this test assumed it could not, and 32 cap leaves said otherwise.
#[test]
fn a_filament_through_a_sea_needs_the_whiteness_arm_and_a_shore_does_not() {
    let x0 = 0.123;
    let levels = 5u32;
    let on_filament = |t: &QuadTree, i: usize| {
        let q = &t.nodes[i];
        let hx = 2.0 * q.half / (t.n - 1) as f64;
        (x0 - q.cx).abs() <= q.half + 0.5 * hx
    };
    // `off`: neither the stationarity stop nor the agreement arm -- the area floor alone.
    let off = SchedCfg { stationary: false, agreement: false, ..cfg(levels, 2000) };
    let on = SchedCfg { stationary: true, ..cfg(levels, 2000) };

    let through = prin_rs::testing::filament_through_sea(x0, 11, T);
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &off, T, &through);
    println!("filament through sea, floor only: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    assert_eq!(count(&t, Decision::MaxLevel), 0, "by area alone the filament through the sea cannot clear the floor");
    assert_eq!(count(&t, Decision::Floor), 16, "everything floors at the bootstrap");

    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &on, T, &through);
    println!("filament through sea, both arms: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    assert!(count(&t, Decision::MaxLevel) > 0, "with the whiteness arm the filament should reach the cap");
    for i in t.leaves() {
        let q = &t.nodes[i];
        match q.decision {
            Decision::MaxLevel => assert!(on_filament(&t, i), "a cap leaf off the filament: cx {} half {}", q.cx, q.half),
            Decision::Floor | Decision::Stationary => {
                assert!(!on_filament(&t, i), "a filament quad stopped as {:?}: cx {} half {}", q.decision, q.cx, q.half)
            }
            _ => {}
        }
    }

    let shore = prin_rs::testing::filament_in_sea(x0, 11, T);
    let (t, st) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &off, T, &shore);
    println!("shore, floor only: {} quads, stop [{}]", st.quads_computed, t.stop_breakdown());
    assert!(count(&t, Decision::MaxLevel) > 0, "the shore is a line and the floor alone should follow it");
    assert!(count(&t, Decision::Floor) > 0, "the sea behind the shore should floor");
    for i in t.leaves() {
        let q = &t.nodes[i];
        match q.decision {
            Decision::MaxLevel => assert!(on_filament(&t, i), "a cap leaf off the shore: cx {} half {}", q.cx, q.half),
            // A floored leaf sits on the sea side: a quad whose right edge just touches the
            // filament holds sea and filament and no basin, gains no area by splitting, and
            // floors, while its basin-side sibling carries the shore to the cap.
            Decision::Floor => assert!(q.cx < x0, "a floored leaf off the sea side: cx {} half {}", q.cx, q.half),
            _ => {}
        }
    }
}

/// **T16 — the live tree merges back once the band has collapsed.** The pulse's band widens
/// and then collapses to a step; with merging on, the parents of the band's leaves become
/// resolved at a late boundary and their children are released, so the resident count peaks and
/// falls and the tree at the end **is** the static tree at the horizon. T11 is the control:
/// with merging off the same march ends strictly finer than the static tree.
#[test]
fn the_live_tree_merges_back_after_the_band_collapses() {
    let field = prin_rs::testing::pulse(0.123, 0.35, 8, T);
    let levels = 4u32;
    let (t, st) = scheduler::descend_live_with(0.0, 0.0, 1.0, 0, &cfg(levels, 2000), T, &field);
    println!("live pulse with merging: {} quads computed, merged {}, resident peak {} final {}, stop [{}]",
             st.quads_computed, st.merged, st.resident_peak, st.resident_final, t.stop_breakdown());
    for p in &st.live {
        println!("  j={} t={:.2} computed={} leaves={} resident={} split={} keep={} merged={}",
                 p.j, p.t, p.computed, p.leaves, p.resident, p.split, p.keep, p.merged);
    }
    assert!(st.merged > 0, "nothing merged");
    // The resident count falls after its FIRST peak: the memory a live design gives back. The
    // first peak, because the tree may legitimately regrow to the same size later -- at the
    // horizon the collapsed band is a step column those quads re-split to follow, and under the
    // `alpha_lo = 0.005` default the final tree is exactly as large as the band was at its widest
    // (149 quads, after a trough of 117). Reading the last maximum found nothing after it.
    let peak = st.live.iter().map(|p| p.resident).max().unwrap();
    let peak_at = st.live.iter().position(|p| p.resident == peak).unwrap();
    let trough = st.live[peak_at..].iter().map(|p| p.resident).min().unwrap();
    println!("  resident first peak {peak} at j={}, trough after it {trough}", st.live[peak_at].j);
    assert!(trough < peak, "the resident count never fell after its first peak");
    assert_eq!(st.resident_final, t.resident());
    // No merged quad is a leaf, and every leaf tiles the root: the leaf areas sum to the root's.
    let area: f64 = t.leaves().map(|i| 4.0 * t.nodes[i].half * t.nodes[i].half).sum();
    assert!((area - 4.0).abs() < 1e-9, "the leaves do not tile the root: area {area}");
    assert!(t.leaves().all(|i| !t.nodes[i].merged));
    // The tree at the end is the static tree at the horizon.
    let (ts, _) = scheduler::descend_with(0.0, 0.0, 1.0, 0, &cfg(levels, 2000), T, &field);
    let key = |t: &QuadTree, i: usize| (t.nodes[i].level, t.nodes[i].cx.to_bits(), t.nodes[i].cy.to_bits());
    let mut a: Vec<_> = t.leaves().map(|i| key(&t, i)).collect();
    let mut b: Vec<_> = ts.leaves().map(|i| key(&ts, i)).collect();
    a.sort_unstable();
    b.sort_unstable();
    let hist = |t: &QuadTree| { let mut h = std::collections::BTreeMap::new(); for i in t.leaves() { *h.entry(t.nodes[i].level).or_insert(0) += 1; } h };
    println!("  live leaves by level {:?}, static {:?}; live-only {}, static-only {}",
             hist(&t), hist(&ts), a.iter().filter(|k| !b.contains(k)).count(), b.iter().filter(|k| !a.contains(k)).count());
    println!("  static stop [{}]; live stop [{}]", ts.stop_breakdown(), t.stop_breakdown());
    // Under `alpha_lo = 0.005` this is the assertion that caught the cap gap: with capped leaves
    // terminal, 24 resolved parents held 96 capped children and the live tree ended at 149 quads
    // against the static 69. A capped leaf is re-decided every boundary and reads `Keep` once its
    // region resolves, and the merge pass takes a capped child as settled.
    assert_eq!(a, b, "the merged live tree is not the static tree at the horizon");
    // Before the first merge the tree only grew.
    let first_merge = st.live.iter().position(|p| p.merged > 0).expect("a merge happened");
    assert!(first_merge > 0);
    for j in 1..=first_merge {
        assert!(st.live[j].computed >= st.live[j - 1].computed);
    }
}
