//! Quadtree invariants — the geometry `alpha` depends on, and the guards that keep the descent
//! honest. These run before any measurement, because every question the scheduler exists to
//! answer is meaningless if the tree is not what it claims to be.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::quad::{Decision, QuadTree, MIN_SAMPLES_PER_AXIS, PRECISION_MARGIN};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, Order, Policy, SchedCfg};

const N: usize = 8;

/// **The load-bearing geometry.** `cell_width = 2h/(N-1)`, and a child's is *exactly* half its
/// parent's. That factor of two is what `alpha = log2(sp_parent/sp_child)` is measured against; if
/// it drifts, every exponent in the tree is wrong by an unknown amount.
#[test]
fn a_childs_cell_width_is_exactly_half_its_parents() {
    let mut t = QuadTree::new(1.0, 3.0, 0.05, N, 0);
    println!("{:>6}{:>16}{:>16}{:>10}", "level", "half", "cell width", "ratio");
    let mut cur = 0usize;
    let mut prev = t.nodes[0].cell_width(N);
    println!("{:>6}{:>16.6e}{:>16.6e}{:>10}", 0, t.nodes[0].half, prev, "-");
    for l in 1..=10 {
        let kids = t.split(cur, l as u32);
        cur = kids[0];
        let w = t.nodes[cur].cell_width(N);
        // Bitwise: halving is exact in binary, so anything else is a formula error.
        assert_eq!(w, prev / 2.0, "level {l}: cell width {w} is not exactly half of {prev}");
        assert_eq!(t.nodes[cur].level, l as u32);
        println!("{:>6}{:>16.6e}{:>16.6e}{:>10.4}", l, t.nodes[cur].half, w, prev / w);
        prev = w;
    }
    println!();
    println!("Exact in binary at every level, so alpha's factor of two holds by construction.");
}

/// `Slice::axis` corner-anchors at `n <= 1` and `cell_widths` clamps to `nx.max(2)-1`, so a
/// one-sample quad would sit at its lower-left corner and be jittered by the whole box. Rejected,
/// not left as a trap.
#[test]
#[should_panic(expected = "samples per quad axis")]
fn a_single_sample_quad_is_rejected() {
    let t = QuadTree::new(1.0, 3.0, 0.05, N, 0);
    let _ = t.nodes[0].slice(MIN_SAMPLES_PER_AXIS - 1, 0);
}

/// The quad `Slice` must sample the box it claims to: corners inclusive, spacing uniform.
#[test]
fn a_quads_slice_spans_its_own_box() {
    let t = QuadTree::new(1.0, 3.0, 0.05, N, 0);
    let q = &t.nodes[0];
    let s = q.slice(N, 0);
    let (lo, _) = s.decode_pos(0);
    let (hi, _) = s.decode_pos(N - 1);
    assert_eq!(lo, q.cx - q.half, "first sample is not the lower edge");
    assert_eq!(hi, q.cx + q.half, "last sample is not the upper edge");
    let (hx, hy) = s.cell_widths();
    assert_eq!(hx, q.cell_width(N));
    assert_eq!(hy, q.cell_width(N));
    println!("N={N}: samples span [{lo}, {hi}] with cell width {hx:.6e}");
}

/// Leaves **tile** the root exactly — no overlap, no gap. Areas summing is necessary; checking
/// that no two leaves overlap is what makes it sufficient.
#[test]
fn leaves_tile_the_root_exactly() {
    let mut t = QuadTree::new(1.0, 3.0, 0.05, N, 0);
    // An irregular tree, so this is not just testing a uniform grid.
    let a = t.split(0, 1);
    let b = t.split(a[0], 2);
    t.split(a[3], 2);
    t.split(b[2], 3);

    let leaves: Vec<usize> = t.leaves().collect();
    let root_area = (2.0 * t.nodes[0].half).powi(2);
    let sum: f64 = leaves.iter().map(|&i| (2.0 * t.nodes[i].half).powi(2)).sum();
    assert!((sum - root_area).abs() < 1e-15 * root_area,
            "leaf areas sum to {sum}, root is {root_area}");

    for (u, &i) in leaves.iter().enumerate() {
        for &j in leaves.iter().skip(u + 1) {
            let (p, q) = (&t.nodes[i], &t.nodes[j]);
            let sep_x = (p.cx - q.cx).abs();
            let sep_y = (p.cy - q.cy).abs();
            let touch_x = p.half + q.half;
            let overlap = sep_x < touch_x * (1.0 - 1e-12) && sep_y < touch_x * (1.0 - 1e-12);
            assert!(!overlap, "leaves {i} and {j} overlap");
        }
    }
    println!("{} leaves tile the root: areas sum exactly, no pair overlaps", leaves.len());
}

/// **Never pooled.** A parent's reduction must come from its own footprints at its own cell width,
/// not from its children — the surrogate error measured for pooling is +38.6%, flat in E.
///
/// Checked by construction *and* by consequence: the parent is computed at an iteration strictly
/// before its children exist, and its `n_footprints` is exactly `N²` rather than `4N²`.
#[test]
fn a_parent_is_never_synthesised_from_its_children() {
    let cfg = SchedCfg { n: N, budget: 40, bootstrap_levels: 2, ..Default::default() };
    let ens = EnsembleCfg { t_max: 2.0, n_sync: 8, refine_flagged: false, ..Default::default() };
    let (t, _) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64);

    for i in 0..t.nodes.len() {
        let q = &t.nodes[i];
        if q.red.n_footprints == 0 {
            continue;
        }
        assert_eq!(q.red.n_footprints as usize, N * N,
                   "quad {i} reduced {} footprints, not N^2 = {}", q.red.n_footprints, N * N);
        if let Some(kids) = q.children {
            for k in kids {
                assert!(t.nodes[k].iteration > q.iteration,
                        "child {k} was computed at iteration {}, not after its parent's {}",
                        t.nodes[k].iteration, q.iteration);
            }
        }
    }
    println!("every computed quad reduced exactly {} footprints of its own", N * N);
}

/// The precision guard must fire before the cell width reaches `f64::EPSILON`, and it must be
/// distinguishable from a physical stop — a descent that hits it has hit a *numerical* floor.
#[test]
fn the_precision_floor_fires_above_machine_epsilon() {
    let mut t = QuadTree::new(1.0, 3.0, 0.05, N, 0);
    let mut cur = 0usize;
    let mut fired = None;
    for l in 1..60u32 {
        let kids = t.split(cur, l);
        cur = kids[0];
        if t.nodes[cur].below_precision_floor(N) {
            fired = Some(l);
            break;
        }
    }
    let l = fired.expect("the guard never fired in 60 levels");
    let w = t.nodes[cur].cell_width(N);
    let scale = t.nodes[cur].cx.abs().max(t.nodes[cur].cy.abs()).max(1.0);
    println!("guard fires at level {l}: cell width {w:.4e}, {:.1}x eps*scale",
             w / (f64::EPSILON * scale));
    println!("machine epsilon would be reached at level ~{:.1}",
             (0.1f64 / ((N - 1) as f64 * f64::EPSILON)).log2());
    assert!(w > f64::EPSILON * scale, "the guard fired below epsilon, which is too late");
    assert!(w < PRECISION_MARGIN * 2.0 * f64::EPSILON * scale, "the guard fired far too early");
    assert!((30..42).contains(&l), "guard fired at level {l}, expected ~36");
}

/// The budget is a cap on **quads**, respected exactly, and quads that wanted to split but could
/// not are **reported** rather than silently dropped.
#[test]
fn the_budget_is_respected_and_exhaustion_is_reported() {
    let ens = EnsembleCfg { t_max: 2.0, n_sync: 8, refine_flagged: false, ..Default::default() };
    for budget in [5usize, 21, 60] {
        let cfg = SchedCfg { n: 4, budget, bootstrap_levels: 6, ..Default::default() };
        let (t, st) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64);
        assert!(st.quads_computed <= budget,
                "computed {} quads against a budget of {budget}", st.quads_computed);
        let reported = t.nodes.iter().filter(|q| q.decision == Decision::BudgetExhausted).count();
        println!("budget {budget:>3}: computed {:>3}, exhausted {}, {} quads flagged",
                 st.quads_computed, st.budget_exhausted, reported);
        if st.budget_exhausted {
            assert!(reported > 0, "budget ran out but nothing was flagged");
        }
    }
}

/// Both policies must run and produce a tree. Whether the sibling one is *better* is a
/// measurement, not a test.
#[test]
fn both_policies_descend() {
    let ens = EnsembleCfg { t_max: 2.0, n_sync: 8, refine_flagged: false, ..Default::default() };
    for policy in [Policy::Alpha, Policy::Sibling] {
        for order in [Order::Spread, Order::SpreadArea, Order::Shuffled] {
            let cfg = SchedCfg { n: 4, budget: 60, policy, order, ..Default::default() };
            let (t, st) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64);
            assert!(st.quads_computed > 0);
            println!("{:>8} / {:>12}: {:>3} quads, {:>3} leaves, depth {}",
                     policy.name(), order.name(), st.quads_computed,
                     t.leaves().count(), t.depth_histogram().len() - 1);
        }
    }
}

/// The sibling signal is the **range** of four exponents, set on the parent once all four children
/// exist — not before, and not on a child.
#[test]
fn the_sibling_range_is_set_on_the_parent_from_four_children() {
    let cfg = SchedCfg { n: 4, budget: 200, bootstrap_levels: 3, ..Default::default() };
    let ens = EnsembleCfg { t_max: 2.0, n_sync: 8, refine_flagged: false, ..Default::default() };
    let (t, _) = scheduler::descend(1.0, 3.0, 0.05, 0, &cfg, &ens, Precision::F64);

    let mut checked = 0usize;
    for i in 0..t.nodes.len() {
        let q = &t.nodes[i];
        match (q.children, q.alpha_sibling_spread) {
            (None, Some(s)) => panic!("leaf {i} carries a sibling range of {s}"),
            (Some(kids), Some(s)) => {
                let a: Vec<f64> = kids.iter().filter_map(|&k| t.nodes[k].alpha).collect();
                assert_eq!(a.len(), 4, "quad {i} has a range but not four child exponents");
                let want = a.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                    - a.iter().cloned().fold(f64::INFINITY, f64::min);
                assert_eq!(s, want, "quad {i}: range {s} is not max-min of {a:?}");
                checked += 1;
            }
            _ => {}
        }
    }
    assert!(checked > 0, "no parent carried a sibling range, so this test proves nothing");
    println!("{checked} parents carry the range of their four children's exponents");
}
