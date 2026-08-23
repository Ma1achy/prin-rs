//! The vertical slice: chart, camera, adaptive render, SSAA, decode paths.
//!
//! Every assertion here was chosen by asking what would have to be true for it to **fire**.
//! Three of them exist specifically because the natural version could not: a curvature term on
//! an affine chart, a texel-scaling check that a uniform render would pass, and a
//! divergence-agreement check that two collapsed decode paths would pass together.

use prin_rs::camera::Camera;
use prin_rs::decode::{self, Path};
use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::output::adaptive::{self, TexelMode};
use prin_rs::physics::{burrau, shape, Cart};
use prin_rs::quad::{Decision, Quad, QuadTree};
use prin_rs::render::Precision;
use prin_rs::scheduler::{self, SchedCfg};
use prin_rs::Vec2;

// ---------------------------------------------------------------------------------------
// A. the chart
// ---------------------------------------------------------------------------------------

/// `BodyPlane` must be **bitwise** what it was before the chart existed. Written against the
/// semantics (numpy's linspace, the other two bodies untouched, released from rest) rather
/// than against a recorded output, so it is not circular.
#[test]
fn body_plane_decode_is_bitwise_unchanged() {
    let s = grid::region("near-field", 5, 3, 0.05).unwrap();
    let (a, b) = (1.0 - 0.05, 1.0 + 0.05);
    for idx in 0..s.npix() {
        let (jx, jy) = (idx % 5, idx / 5);
        let ex = if jx == 4 { b } else { a + jx as f64 * ((b - a) / 4.0) };
        let ey = if jy == 2 { 3.05 } else { 2.95 + jy as f64 * (0.1 / 2.0) };
        let (x, y) = s.decode_pos(idx);
        assert_eq!(x.to_bits(), ex.to_bits(), "x at {idx}");
        assert_eq!(y.to_bits(), ey.to_bits(), "y at {idx}");

        let c: Cart<f64> = s.nominal(idx);
        assert_eq!(c.r[0].x.to_bits(), ex.to_bits());
        assert_eq!(c.r[0].y.to_bits(), ey.to_bits());
        assert_eq!(c.r[1], Vec2::new(-2.0, -1.0));
        assert_eq!(c.r[2], Vec2::new(1.0, -1.0));
        assert!(c.v.iter().all(|v| *v == Vec2::zero()));
    }
}

/// The oblique chart's degenerate case must reproduce the axis-aligned one **exactly**, or the
/// two families are not comparable and §3.5's "does slice type matter" has a confound.
#[test]
fn plane_chart_reproduces_body_plane_exactly() {
    for body in 0..3 {
        let base = Slice::body_plane(7, 7, 0.4, -1.1, 0.3, body);
        let obl = base.with_chart(Chart::plane_for_body(body));
        for idx in 0..base.npix() {
            let (a, b): (Cart<f64>, Cart<f64>) = (base.nominal(idx), obl.nominal(idx));
            assert_eq!(decode::bits(&a), decode::bits(&b), "body {body} idx {idx}");
        }
    }
}

/// The Hopf inverse round-trips. **This can fail**: `shape_vec` computes
/// `q = rt.y*lt.x - rt.x*lt.y`, the negative of the standard cross product, so the inverse
/// needs `atan2(-q, p)`. The reflected sign is checked to be caught, not merely absent.
#[test]
fn shape_chart_round_trips_and_the_sign_matters() {
    let m = burrau::MASSES;
    let mut worst = 0.0f64;
    for a in 1..16 {
        for b in 0..16 {
            let (th, ph) = (
                std::f64::consts::PI * a as f64 / 16.0,
                2.0 * std::f64::consts::PI * b as f64 / 16.0,
            );
            let n = [th.cos(), th.sin() * ph.cos(), th.sin() * ph.sin()];
            let r = shape::from_shape(n, 12.0, 0.4, &m);
            let back = shape::shape_vec(&r, &m);
            worst = worst.max((0..3).map(|k| (back[k] - n[k]).abs()).fold(0.0, f64::max));

            let com = (r[0] * m[0] + r[1] * m[1] + r[2] * m[2]) / (m[0] + m[1] + m[2]);
            assert!(com.x.abs() < 1e-12 && com.y.abs() < 1e-12, "com not at origin");
        }
    }
    assert!(worst < 1e-13, "round-trip worst |dn| = {worst:.3e}");

    // The teeth: reflecting n0 must break the round-trip, or the test is decoration.
    let n = [0.3, 0.5, (1.0f64 - 0.09 - 0.25).sqrt()];
    let bad = shape::shape_vec(&shape::from_shape([-n[0], n[1], n[2]], 5.0, 0.7, &m), &m);
    let d = (0..3).map(|k| (bad[k] - n[k]).abs()).fold(0.0, f64::max);
    assert!(d > 0.1, "a reflected shape point must not round-trip, got {d:.3e}");
}

/// The whole reason the shape chart exists: its decode is **not** affine, so the linearised
/// path is a genuine approximation. On `BodyPlane` the same quantity is exactly zero — and
/// that zero is structural, never a measurement.
#[test]
fn shape_chart_is_nonlinear_and_body_plane_is_not() {
    let sh = Chart::shape_at_burrau(0.4);
    let lin = decode::linearise(&sh, 0, 0.0, 0.0, 0.25);
    let a = decode::sample(Path::DirectF64, &sh, 0, 0.0, 0.0, 0.25, 0.5, -0.5, &lin);
    let b = decode::sample(Path::LinSplitF64, &sh, 0, 0.0, 0.0, 0.25, 0.5, -0.5, &lin);
    assert!(decode::max_abs_diff(&a, &b) > 1e-6, "shape chart must carry curvature");
    assert!(!sh.is_affine());

    let bp = Chart::BodyPlane;
    let lin = decode::linearise(&bp, 0, 1.0, 3.0, 0.25);
    let a = decode::sample(Path::DirectF64, &bp, 0, 1.0, 3.0, 0.25, 0.5, -0.5, &lin);
    let b = decode::sample(Path::LinSplitF64, &bp, 0, 1.0, 3.0, 0.25, 0.5, -0.5, &lin);
    assert_eq!(decode::max_abs_diff(&a, &b), 0.0, "affine chart: curvature is structurally 0");
    assert!(bp.is_affine());
}

// ---------------------------------------------------------------------------------------
// B. the screen floor is a veto
// ---------------------------------------------------------------------------------------

fn quad_at(level: u32, half: f64) -> Quad {
    let mut t = QuadTree::new(1.0, 3.0, 0.05, 8, 0);
    t.nodes[0].half = half;
    t.nodes[0].level = level;
    t.nodes[0].clone()
}

/// **Structural**: the camera's entire interface with the scheduler cannot return `Split`.
/// The contract strikes "below screen resolution -> refine" explicitly, and this is what
/// stops it being reintroduced.
#[test]
fn the_camera_can_only_ever_stop_a_descent() {
    let cam = Camera::framing(1.0, 3.0, 0.05, 512);
    for level in 0..24u32 {
        let q = quad_at(level, 0.05 / (2f64).powi(level as i32));
        match cam.veto(&q, 8, 0.05) {
            None | Some(Decision::ScreenFloor) | Some(Decision::MaxRelDepth) => {}
            other => panic!("the camera returned {other:?}; it is a veto, not a trigger"),
        }
    }
}

/// The arithmetic in the brief, asserted: `N = 8`, 512² viewport, root framed — samples stop
/// being displayable at **level 6**, and PR #11 descended to 12.
#[test]
fn screen_floor_lands_at_level_six_at_n8_on_512() {
    let cam = Camera::framing(1.0, 3.0, 0.05, 512);
    for level in 0..10u32 {
        let q = quad_at(level, 0.05 / (2f64).powi(level as i32));
        assert_eq!(cam.screen_floor(&q, 8), level >= 6, "level {level}");
    }
}

/// **View-relative and uncached.** The same quad, floored at one zoom, must be refinable at
/// the next — otherwise zooming in would show a permanently coarse patch.
#[test]
fn a_floored_quad_refines_again_when_zoomed_into() {
    let q = quad_at(6, 0.05 / 64.0);
    let framed = Camera::framing(1.0, 3.0, 0.05, 512);
    assert!(framed.screen_floor(&q, 8), "must be floored at the framing zoom");
    let zoomed = Camera::at_depth(q.cx, q.cy, 0.05, 512, 3);
    assert!(!zoomed.screen_floor(&q, 8), "must refine again once zoomed in");
    // And `max_rel_depth` moves with the view rather than capping it absolutely.
    assert!(zoomed.veto(&q, 8, 0.05).is_none());
}

// ---------------------------------------------------------------------------------------
// C. the adaptive render
// ---------------------------------------------------------------------------------------

fn toy_tree() -> (QuadTree, Vec<Vec<prin_rs::ensemble::pixel::PixelOut>>) {
    let mut t = QuadTree::new(1.0, 3.0, 0.05, 4, 0);
    let kids = t.split(0, 1);
    let deeper = t.split(kids[0], 2);
    let mut px = vec![Vec::new(); t.nodes.len()];
    for i in 0..t.nodes.len() {
        px[i] = vec![prin_rs::ensemble::pixel::PixelOut::default(); 16];
    }
    let _ = deeper;
    (t, px)
}

/// **The acceptance test, and the negative case is the point.** A level-3 leaf must be drawn
/// with 4x the texel size of a level-5 leaf. PR #11's uniform render draws them the same, and
/// the *same* assertion has to reject it, or the failure passes unnoticed.
#[test]
fn texel_size_varies_as_two_to_the_minus_level_and_uniform_is_rejected() {
    let (t, px) = toy_tree();
    let cam = Camera::framing(1.0, 3.0, 0.05, 128);

    let (_, ad) = adaptive::render(&t, &px, &cam, 128, TexelMode::Adaptive, |_| [0, 0, 0]);
    let slope = adaptive::texel_scaling(&ad).expect("more than one leaf level");
    assert!((slope + 1.0).abs() < 1e-12, "adaptive slope {slope}, want -1");

    let (_, un) = adaptive::render(&t, &px, &cam, 128, TexelMode::Uniform, |_| [0, 0, 0]);
    let slope_u = adaptive::texel_scaling(&un).expect("more than one leaf level");
    assert!(
        (slope_u + 1.0).abs() > 0.5,
        "a uniform render must FAIL the texel assertion; slope was {slope_u}"
    );
}

/// The SSAA footprint has to scale with the texel, or a coarse quad anti-aliases over the
/// wrong area. It does so by construction — the copies are jittered by `cell_width` — and this
/// asserts the construction rather than trusting it.
#[test]
fn ssaa_copy_footprint_scales_with_the_texel() {
    let spread = |half: f64| -> f64 {
        let s = Slice::body_plane(4, 4, 1.0, 3.0, half, 0);
        let c: Vec<Cart<f64>> = jitter::copies_with(&s, 5, 7, 0.5, 0, Scheme::Halton);
        let (lo, hi) = c.iter().fold((f64::MAX, f64::MIN), |(l, h), x| {
            (l.min(x.r[0].x), h.max(x.r[0].x))
        });
        hi - lo
    };
    let (a, b) = (spread(0.05), spread(0.05 / 2.0));
    assert!((a / b - 2.0).abs() < 1e-9, "footprint ratio {} want 2", a / b);
}

// ---------------------------------------------------------------------------------------
// D. decode: distinctness, not agreement
// ---------------------------------------------------------------------------------------

fn grid_samples(path: Path, chart: &Chart, cu: f64, cv: f64, half: f64, n: usize) -> Vec<Cart<f64>> {
    let lin = decode::linearise(chart, 0, cu, cv, half);
    let d = decode::deltas(n);
    let mut v = Vec::with_capacity(n * n);
    for &dv in &d {
        for &du in &d {
            v.push(decode::sample(path, chart, 0, cu, cv, half, du, dv, &lin));
        }
    }
    v
}

/// **Distinctness, not agreement.** Asserting that the linearised path *agrees* with f64 at
/// depth would pass when both paths have collapsed to a single initial condition — they agree
/// perfectly, and the agreement means nothing. This asserts the linearised path still resolves
/// its samples where the direct one no longer does.
#[test]
fn lin_split_f32_keeps_samples_distinct_where_direct_f32_has_collapsed() {
    let chart = Chart::BodyPlane;
    let (cu, cv, n) = (1.0, 3.0, 8);
    let half = 0.05 / (2f64).powi(30); // ~5e-11: well past f32's reach on an O(1) coordinate

    let d32 = decode::distinct(&grid_samples(Path::DirectF32, &chart, cu, cv, half, n));
    let split = decode::distinct(&grid_samples(Path::LinSplitF32, &chart, cu, cv, half, n));
    assert_eq!(d32, 1, "direct f32 should have collapsed here, got {d32} distinct");
    assert_eq!(split, n * n, "lin-split must still resolve all {} samples, got {split}", n * n);
}

/// The control has teeth, and it caught the design claim being weaker than advertised.
///
/// `L-naive` collapses on **exactly the same curve as plain `direct_f32`** — 56, 18, 2, 1 at
/// depths 16, 18, 20, 22. The literal formula buys nothing at all: adding a ~1e-13 term to an
/// O(1) f32 quantity is the same operation as forming the coordinate in f32 in the first place.
/// A change that made this path work would fire this test rather than pass quietly.
#[test]
fn lin_naive_f32_collapses_on_the_same_curve_as_plain_f32() {
    let chart = Chart::BodyPlane;
    for depth in [16, 18, 20, 22, 26, 40] {
        let half = 0.05 / (2f64).powi(depth);
        let naive = decode::distinct(&grid_samples(Path::LinNaiveF32, &chart, 1.0, 3.0, half, 8));
        let plain = decode::distinct(&grid_samples(Path::DirectF32, &chart, 1.0, 3.0, half, 8));
        assert_eq!(naive, plain, "at depth {depth}: naive {naive}, plain f32 {plain}");
    }
    let half = 0.05 / (2f64).powi(22);
    assert_eq!(decode::distinct(&grid_samples(Path::LinNaiveF32, &chart, 1.0, 3.0, half, 8)), 1);
}

/// **The measured limit of the linearisation, asserted as the negative it is.**
///
/// `L-split_f32` tracks `direct_f64` *exactly* — both hold all 64 samples to depth 44 and both
/// reach 1 by depth 50. So the split form buys ~24 levels over f32 (22 -> 46) and **exactly
/// zero over f64**. The contract's "~23 to ~50+" is real only in the sense that it lets an f32
/// GPU reach the f64 CPU's floor; it does not push past it, because the initial conditions must
/// still be formed as absolute O(1) numbers before integration.
///
/// This asserts the equality, so a future change claiming to beat the f64 floor has to break a
/// test to do it.
#[test]
fn lin_split_f32_reaches_the_f64_floor_and_no_further() {
    let chart = Chart::BodyPlane;
    for depth in [0, 20, 40, 44, 45, 46, 48, 50, 52] {
        let half = 0.05 / (2f64).powi(depth);
        let split = decode::distinct(&grid_samples(Path::LinSplitF32, &chart, 1.0, 3.0, half, 8));
        let f64d = decode::distinct(&grid_samples(Path::DirectF64, &chart, 1.0, 3.0, half, 8));
        assert_eq!(split, f64d, "at depth {depth}: split {split}, direct f64 {f64d}");
    }
    let deep = 0.05 / (2f64).powi(50);
    assert_eq!(decode::distinct(&grid_samples(Path::LinSplitF32, &chart, 1.0, 3.0, deep, 8)), 1,
               "the split form has a floor too, and it is f64's");
}

// ---------------------------------------------------------------------------------------
// E. MAX_REL_DEPTH is scheduler state, not payload
// ---------------------------------------------------------------------------------------

fn cheap_ens() -> EnsembleCfg {
    EnsembleCfg { n_extra: 1, t_max: 1.0, n_sync: 4, eta: 0.05, refine_flagged: false, ..Default::default() }
}

/// Lowering the relative-depth cap must stop scheduling deeper and change **nothing** about
/// the quads that survive. If it moved a payload, a zoom-out would invalidate cached work.
#[test]
fn lowering_max_rel_depth_changes_no_payload() {
    let ens = cheap_ens();
    let run = |m: u32| {
        let mut cam = Camera::framing(1.0, 13.0, 0.05, 512);
        cam.max_rel_depth = Some(m);
        let cfg = SchedCfg {
            n: 4, budget: 400, camera: Some(cam), tau_display: 1e-12, alpha_hi: -10.0,
            ..Default::default()
        };
        scheduler::descend(1.0, 13.0, 0.05, 0, &cfg, &ens, Precision::F64)
    };
    let (deep, _) = run(4);
    let (shallow, _) = run(2);
    assert!(shallow.nodes.len() < deep.nodes.len(), "the cap must actually bind");
    for (a, b) in shallow.nodes.iter().zip(deep.nodes.iter()) {
        assert_eq!(a.level, b.level);
        assert_eq!(a.cx.to_bits(), b.cx.to_bits());
        assert_eq!(
            a.red.spread_median.to_bits(),
            b.red.spread_median.to_bits(),
            "a surviving quad's payload moved when only the cap changed"
        );
    }
}
