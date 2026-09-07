//! The physics-space metric, pinned: `error(full)` is exactly zero, the tolerance saturates at
//! both ends, the three forms order as they must, the DP ceiling holds under it, and the v2
//! footprint file carries the event class while a v1 file refuses to pretend it does.

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::Chart;
use prin_rs::metric::{self, Cache, ClassArm, Colouring, Metric, PayloadForm, Rank};
use prin_rs::output::fcache::{self, EVENT_CLASS_ABSENT};
use prin_rs::quad::{Agg, Criterion, StructureMode};

const LEVELS: u32 = 2;
const N: usize = 2;

fn ens() -> EnsembleCfg {
    EnsembleCfg { n_extra: 1, t_max: 2.0, n_sync: 4, refine_flagged: false, ..Default::default() }
}

fn payload(eps: f64, form: PayloadForm) -> Metric {
    Metric::Payload { eps, class: ClassArm::EventClass, form }
}

fn build(metrics: &[Metric]) -> (Vec<Cache>, std::collections::HashMap<metric::Key, Vec<prin_rs::ensemble::pixel::PixelOut>>) {
    let res = (1usize << LEVELS) * N;
    metric::build_metrics_with_footprints(
        "deep interior", 0.0, 0.0, 0.05, 0, Chart::BodyPlane, LEVELS, N, res, 1e-2, &ens(), metrics,
    )
}

fn deepest(c: &Cache) -> Vec<metric::Key> {
    let w = 1u32 << c.levels;
    (0..w).flat_map(|iy| (0..w).map(move |ix| (c.levels, ix, iy))).collect()
}

/// **`error(full) == 0` exactly, and the tolerance saturates at both ends.** At `eps = 0` any
/// difference counts and the root must read above zero (the control that the field is not
/// constant); at `eps = 1` nothing exceeds the antipodal bound and the root must read exactly
/// zero. A metric that cannot reach both ends is not a tolerance.
#[test]
fn error_full_is_zero_and_the_tolerance_saturates_at_both_ends() {
    let (caches, _) = build(&[
        payload(0.0, PayloadForm::Indicator),
        payload(1.0, PayloadForm::Indicator),
        payload(0.0, PayloadForm::Hinge),
        payload(0.0, PayloadForm::Resolvable),
    ]);
    for c in &caches {
        assert_eq!(c.error_of(&deepest(c)), 0.0, "{}: error(full) must be exactly zero", c.metric.name());
        let root = c.error_of(&[(0, 0, 0)]);
        assert!((0.0..=1.0).contains(&root), "{}: root error {root} outside [0, 1]", c.metric.name());
    }
    let root0 = caches[0].error_of(&[(0, 0, 0)]);
    let root1 = caches[1].error_of(&[(0, 0, 0)]);
    assert!(root0 > 0.0, "at eps = 0 the root must be unresolved somewhere (field is not constant)");
    assert_eq!(root1, 0.0, "at eps = 1 nothing can exceed the antipodal chord; read {root1}");
    println!("root error: eps=0 {root0:.4}, eps=1 {root1:.4}");
}

/// Hinge and resolvable never exceed the indicator, pixel for pixel, so quad for quad.
#[test]
fn hinge_and_resolvable_are_bounded_by_the_indicator() {
    let eps = 0.05;
    let (caches, _) = build(&[
        payload(eps, PayloadForm::Indicator),
        payload(eps, PayloadForm::Hinge),
        payload(eps, PayloadForm::Resolvable),
    ]);
    let (ind, hin, res) = (&caches[0], &caches[1], &caches[2]);
    let mut any_gap = false;
    for (k, q) in &ind.quads {
        let h = hin.quads[k].err_sum;
        let r = res.quads[k].err_sum;
        assert!(h <= q.err_sum + 1e-12, "quad {k:?}: hinge {h} > indicator {}", q.err_sum);
        assert!(r <= q.err_sum + 1e-12, "quad {k:?}: resolvable {r} > indicator {}", q.err_sum);
        if r < q.err_sum {
            any_gap = true;
        }
    }
    let sea = ind.sea_fraction(eps);
    println!("sea_fraction({eps}) = {sea:.4}; indicator-resolvable gap seen: {any_gap}");
    // The sea fraction is a fraction, and it falls as the tolerance widens.
    assert!((0.0..=1.0).contains(&sea));
    let mut last = f64::INFINITY;
    for e in [0.0, 1e-3, 1e-2, 1e-1, 1.0] {
        let s = ind.sea_fraction(e);
        assert!(s <= last + 1e-15, "sea_fraction is not monotone: {s} after {last} at eps {e}");
        last = s;
    }
    assert_eq!(ind.sea_fraction(1.0), 0.0, "no spread exceeds 1.0");
    assert!(Colouring::Outcome.name().len() > 0);
}

/// The exact ceiling is a ceiling under the payload metric too: no ranking beats it, at any
/// budget, and both ends are pinned.
#[test]
fn no_ranking_beats_the_exact_tree_optimum_under_the_payload_metric() {
    let (caches, _) = build(&[payload(0.01, PayloadForm::Indicator)]);
    let cache = &caches[0];
    let full = cache.quads.len();
    let max_splits = (full - 1) / 4;
    let dp = cache.dp_optimal(max_splits);
    assert_eq!(dp.curve[0], cache.error_of(&[(0, 0, 0)]));
    assert_eq!(dp.curve[max_splits], 0.0);
    for s in 1..=max_splits {
        assert!(dp.curve[s] <= dp.curve[s - 1] + 1e-15);
    }
    let mut ranks: Vec<Rank> = vec![
        Rank::Uniform,
        Rank::GreedyLookahead1,
        Rank::GreedyLookahead1PerCost,
        Rank::StructureOnly,
        Rank::Random(1),
        Rank::Random(7),
    ];
    for c in Criterion::ALL {
        for a in [Agg::Mean, Agg::Median, Agg::P90] {
            ranks.push(Rank::Signal(c, a));
            ranks.push(Rank::Contrast(c, a));
            ranks.push(Rank::Structured(StructureMode::Multiply, c, a));
        }
    }
    let mut worst = 0.0f64;
    for r in ranks {
        for p in metric::replay(cache, r, full) {
            let m = p.error - dp.at_budget(p.budget);
            worst = worst.min(m);
            assert!(m >= -1e-12, "{} beats the optimum at B={} by {m:e}", r.name(), p.budget);
        }
    }
    // The memory number is well-defined at both ends.
    assert_eq!(dp.budget_needed(dp.curve[0]), Some(1), "the root reaches its own error at B = 1");
    assert_eq!(dp.budget_needed(0.0), Some(1 + 4 * max_splits).min(dp.budget_needed(0.0)),
               "zero error is reached no later than the full tree");
    println!("worst margin over every ranking: {worst:+.1e}");
}

/// **v2 carries the event class; v1 reads it back as ABSENT and refuses the event-class arm.**
///
/// The v1 stream is written by hand here, byte for byte the old layout, so the reader's
/// backward path is exercised by a file that never had a 15th column rather than by a v2 file
/// with its version number changed.
#[test]
fn a_v1_footprint_file_refuses_the_event_class_arm_and_replays_the_outcome_arm() {
    let (caches, px_of) = build(&[
        payload(0.01, PayloadForm::Indicator),
        Metric::Payload { eps: 0.01, class: ClassArm::Outcome, form: PayloadForm::Indicator },
    ]);
    let (ev, oc) = (&caches[0], &caches[1]);
    let fp = ev.footprints_from(&px_of, 2.0);
    assert_eq!(fp.version, fcache::VERSION);
    assert!(fp.has_event_class());
    assert!(fp.quads.values().flatten().any(|r| r.event_class != EVENT_CLASS_ABSENT),
            "a v2 footprint set must carry real event classes");

    // v2 round trip keeps the class bitwise, and a remeasure from disk is bitwise the build.
    let mut buf = Vec::new();
    fcache::write(&mut buf, &fp).unwrap();
    let back = fcache::read(&mut std::io::Cursor::new(&buf)).unwrap();
    assert_eq!(back.version, 2);
    for (k, rows) in &fp.quads {
        for (a, b) in back.quads[k].iter().zip(rows) {
            assert_eq!(a.event_class, b.event_class);
        }
    }
    let re = ev.remeasure(&back, ev.metric).unwrap();
    for (k, q) in &ev.quads {
        assert_eq!(re.quads[k].err_sum.to_bits(), q.err_sum.to_bits(), "quad {k:?} moved on replay");
    }

    // The v1 stream, by hand: same header, version 1, 14 columns.
    let mut v1 = Vec::new();
    v1.extend_from_slice(fcache::MAGIC);
    v1.extend_from_slice(&1u32.to_le_bytes());
    let header = format!(
        "region={} body={} cx={:?} cy={:?} half={:?} levels={} n={} res={} t_max={}\nchart={}\nfields=v1\n",
        fp.region, fp.body, fp.cx, fp.cy, fp.half, fp.levels, fp.n, fp.res, fp.t_max, fp.chart
    );
    v1.extend_from_slice(&(header.len() as u32).to_le_bytes());
    v1.extend_from_slice(header.as_bytes());
    let mut keys: Vec<metric::Key> = fp.quads.keys().cloned().collect();
    keys.sort_by_key(|&(l, ix, iy)| (l, iy, ix));
    v1.extend_from_slice(&(keys.len() as u64).to_le_bytes());
    v1.extend_from_slice(&14u32.to_le_bytes());
    for k in &keys {
        let rows = &fp.quads[k];
        v1.extend_from_slice(&k.0.to_le_bytes());
        v1.extend_from_slice(&k.1.to_le_bytes());
        v1.extend_from_slice(&k.2.to_le_bytes());
        v1.extend_from_slice(&(rows.len() as u32).to_le_bytes());
        for r in rows {
            for v in [r.shape[0], r.shape[1], r.shape[2], r.packed as f64, r.n_nonfinite as f64,
                      r.ensemble_spread, r.spread_shape, r.spread_event, r.ftle, r.diffusion,
                      r.t_end, r.error_ratio, r.d_min_true, r.energy_drift_max] {
                v1.extend_from_slice(&v.to_le_bytes());
            }
        }
    }
    let old = fcache::read(&mut std::io::Cursor::new(&v1)).unwrap();
    assert_eq!(old.version, 1);
    assert!(!old.has_event_class());
    assert!(old.quads.values().flatten().all(|r| r.event_class == EVENT_CLASS_ABSENT));

    let err = ev.remeasure(&old, ev.metric).err().expect("a v1 file must refuse the event-class arm");
    assert!(err.contains("v1") && err.contains("event class"), "the refusal must say why: {err}");
    let err = Cache::from_footprints(&old, ev.metric).err().expect("from_footprints must refuse too");
    assert!(err.contains("v1"), "{err}");

    // And the outcome arm replays from the v1 stream, bitwise against the fresh outcome build.
    let re = oc.remeasure(&old, oc.metric).unwrap();
    for (k, q) in &oc.quads {
        assert_eq!(re.quads[k].err_sum.to_bits(), q.err_sum.to_bits(), "outcome-arm quad {k:?} moved on a v1 replay");
    }
    // `from_footprints` reaches the same err_sums with no reductions at all.
    let bare = Cache::from_footprints(&old, oc.metric).unwrap();
    for (k, q) in &oc.quads {
        assert_eq!(bare.quads[k].err_sum.to_bits(), q.err_sum.to_bits());
    }
    let d = bare.dp_optimal((bare.quads.len() - 1) / 4);
    assert_eq!(d.curve[0], bare.error_of(&[(0, 0, 0)]));
}
