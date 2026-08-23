//! Deep zoom: the direct decode against the linearised one, and **the reason distinctness is
//! measured before divergence**.
//!
//! The canonical spec: *"CPU (f64) owns global/nonlinear precision (quad centre/half-width,
//! `x0`, `J_D`); the GPU (f32) does only quad-local relative arithmetic (`x = x0 + J_D . delta`)."*
//! The claimed benefit is extending usable zoom from depth ~23 to ~50+.
//!
//! # The formula does not work as written, and the failure is invisible
//!
//! At depth 40 a quad spans ~1e-13. `x0` is O(1) and `J_D . delta` is ~1e-13, so an **f32 sum
//! returns `x0` for every delta**: all `N²` samples collapse to one initial condition — exactly
//! what the direct path does at that depth. A divergence ladder then compares two *collapsed*
//! sets, finds them in agreement, and reports the linearised path tracking f64 beautifully,
//! from a path that had lost every sample.
//!
//! So the primary measurement here is **distinctness**, not divergence: [`distinct`] counts how
//! many of a set of decoded states are actually different, and a divergence figure is only
//! admissible at a depth where both paths are fully distinct. That measurement can fail; the
//! divergence ladder cannot tell success from mutual collapse.
//!
//! # What the linearisation can and cannot buy
//!
//! **The initial conditions must be formed as absolute O(1) numbers before integration.** The
//! three-body separations are O(1) and no nonlinear integrator can carry `(x0, delta)`
//! separately through the march. So the linearised path can escape a floor set by the *chart
//! coordinate*; it cannot escape one set by the *IC magnitude*. That is a weaker claim than
//! "no floor", and whether the two differ at all is what [`Path::LinSplitF32`] measures.
//!
//! # What was measured
//!
//! On `body_plane` at a chart centre of magnitude ~3, with 64 samples per quad:
//!
//! | path | all 64 distinct to | collapsed to 1 by |
//! |---|---|---|
//! | [`Path::DirectF32`] | depth 14 | depth 22 |
//! | [`Path::LinNaiveF32`] | depth 14 | depth 22 |
//! | [`Path::LinSplitF32`] | depth 44 | depth 50 |
//! | [`Path::DirectF64`] | depth 44 | depth 50 |
//!
//! **The literal formula buys nothing** — it collapses on exactly the same curve as forming the
//! chart coordinate in f32 in the first place. **The split form reaches f64's floor and stops
//! there**, so the gain is ~24 levels *for an f32 consumer* and zero over f64. The contract's
//! "~50+" is f64's floor, not something the linearisation creates. That is the bound above,
//! confirmed rather than escaped.
//!
//! And the floor itself is conditional: the **same box at a chart centre of zero has no
//! cell-width floor at all** in the tested range, on either precision, because there is no O(1)
//! neighbour for the increment to be absorbed into. Quote the coordinate magnitude with any floor
//! depth.
//!
//! # The linearisation is a secant, not a derivative
//!
//! `J_D` is taken as `(D(c + half) - D(c - half)) / 2` per axis — two f64 decodes per axis, as
//! the contract specifies, and the natural *per-quad* linearisation. It makes the quad's edge
//! midpoints exact and the centre worst. On an affine chart it equals the derivative exactly,
//! so the choice is only visible on [`crate::grid::Chart::Shape`].

use crate::grid::{decode_state, Chart};
use crate::physics::Cart;
use crate::Vec2;

/// Which arithmetic forms a sample's initial condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Path {
    /// `u = cu + du*half` in f64, decoded in f64. Floors when the cell width crosses
    /// `f64::EPSILON` relative to the chart coordinate — PR #11's level **45.87**.
    DirectF64,
    /// The same, with the chart coordinate formed in **f32**. No `delta`, so the coordinate
    /// crosses f32 resolution early. Isolates the coordinate-formation floor; a fully-f32
    /// decode cannot floor any later than this.
    DirectF32,
    /// The literal formula: `x0`, `J_D . delta` and the sum all in f32. **Expected to collapse
    /// around depth 20**, and included so that it is *seen* to — it is the control that makes
    /// the distinctness measurement load-bearing.
    LinNaiveF32,
    /// `x0` in f64 on the CPU; `delta` and `J_D . delta` in f32 — the GPU's only quantities,
    /// both well inside f32 range and neither shrinking toward an O(1) neighbour — promoted and
    /// summed in **f64**. The variant the design claim rests on.
    LinSplitF32,
    /// The same split, all in f64. Its difference from [`Path::DirectF64`] is the **curvature**
    /// the linearisation discards — identically zero on any affine chart, and stated as
    /// structural rather than measured there.
    LinSplitF64,
}

impl Path {
    pub fn name(self) -> &'static str {
        match self {
            Path::DirectF64 => "direct_f64",
            Path::DirectF32 => "direct_f32",
            Path::LinNaiveF32 => "lin_naive_f32",
            Path::LinSplitF32 => "lin_split_f32",
            Path::LinSplitF64 => "lin_split_f64",
        }
    }
    pub fn is_linearised(self) -> bool {
        matches!(self, Path::LinNaiveF32 | Path::LinSplitF32 | Path::LinSplitF64)
    }
}

/// A quad's linearisation. Built once per quad on the CPU, in f64.
#[derive(Clone, Copy, Debug)]
pub struct Lin {
    pub x0: Cart<f64>,
    /// d(state)/d(delta_u) — already carries the quad half-width, so `delta` is in `[-1, 1]`.
    pub ju: Cart<f64>,
    pub jv: Cart<f64>,
}

/// **Jacobian cost: four f64 decodes per quad**, two per axis, plus one for the centre.
///
/// Against 512 trajectories per quad at `N = 8, E+1 = 8` this is negligible; the caching
/// contract records it piling up badly enough to hitch a gesture, so the count is reported
/// rather than assumed small.
pub const DECODES_PER_LINEARISATION: usize = 5;

fn axpy(a: &Cart<f64>, s: f64, b: &Cart<f64>) -> Cart<f64> {
    let mut o = *a;
    for k in 0..3 {
        o.r[k] = o.r[k] + b.r[k] * s;
        o.v[k] = o.v[k] + b.v[k] * s;
    }
    o
}

fn sub_scaled(p: &Cart<f64>, m: &Cart<f64>, s: f64) -> Cart<f64> {
    let mut o = Cart::<f64>::default();
    for k in 0..3 {
        o.r[k] = (p.r[k] - m.r[k]) * s;
        o.v[k] = (p.v[k] - m.v[k]) * s;
    }
    o
}

pub fn linearise(chart: &Chart, body: usize, cu: f64, cv: f64, half: f64) -> Lin {
    let x0 = decode_state(chart, body, cu, cv);
    let up = decode_state(chart, body, cu + half, cv);
    let um = decode_state(chart, body, cu - half, cv);
    let vp = decode_state(chart, body, cu, cv + half);
    let vm = decode_state(chart, body, cu, cv - half);
    Lin { x0, ju: sub_scaled(&up, &um, 0.5), jv: sub_scaled(&vp, &vm, 0.5) }
}

/// One sample's initial condition, by path. `du`, `dv` are quad-local in `[-1, 1]`.
pub fn sample(
    path: Path,
    chart: &Chart,
    body: usize,
    cu: f64,
    cv: f64,
    half: f64,
    du: f64,
    dv: f64,
    lin: &Lin,
) -> Cart<f64> {
    match path {
        Path::DirectF64 => decode_state(chart, body, cu + du * half, cv + dv * half),
        Path::DirectF32 => {
            // The chart coordinate itself formed in f32 — where this path's floor lives.
            let u = (cu as f32 + (du as f32) * (half as f32)) as f64;
            let v = (cv as f32 + (dv as f32) * (half as f32)) as f64;
            decode_state(chart, body, u, v)
        }
        Path::LinSplitF64 => axpy(&axpy(&lin.x0, du, &lin.ju), dv, &lin.jv),
        Path::LinSplitF32 => {
            // delta and J.delta in f32; the sum with the O(1) centre in f64.
            let mut o = lin.x0;
            for k in 0..3 {
                let dr = Vec2::new(
                    (lin.ju.r[k].x as f32 * du as f32 + lin.jv.r[k].x as f32 * dv as f32) as f64,
                    (lin.ju.r[k].y as f32 * du as f32 + lin.jv.r[k].y as f32 * dv as f32) as f64,
                );
                let dv_ = Vec2::new(
                    (lin.ju.v[k].x as f32 * du as f32 + lin.jv.v[k].x as f32 * dv as f32) as f64,
                    (lin.ju.v[k].y as f32 * du as f32 + lin.jv.v[k].y as f32 * dv as f32) as f64,
                );
                o.r[k] = o.r[k] + dr;
                o.v[k] = o.v[k] + dv_;
            }
            o
        }
        Path::LinNaiveF32 => {
            // The literal formula. The sum happens in f32 against an O(1) x0.
            let mut o = Cart::<f64>::default();
            for k in 0..3 {
                o.r[k] = Vec2::new(
                    (lin.x0.r[k].x as f32
                        + lin.ju.r[k].x as f32 * du as f32
                        + lin.jv.r[k].x as f32 * dv as f32) as f64,
                    (lin.x0.r[k].y as f32
                        + lin.ju.r[k].y as f32 * du as f32
                        + lin.jv.r[k].y as f32 * dv as f32) as f64,
                );
                o.v[k] = Vec2::new(
                    (lin.x0.v[k].x as f32
                        + lin.ju.v[k].x as f32 * du as f32
                        + lin.jv.v[k].x as f32 * dv as f32) as f64,
                    (lin.x0.v[k].y as f32
                        + lin.ju.v[k].y as f32 * du as f32
                        + lin.jv.v[k].y as f32 * dv as f32) as f64,
                );
            }
            o
        }
    }
}

/// Bit pattern of a state: 12 doubles, compared exactly.
pub fn bits(c: &Cart<f64>) -> [u64; 12] {
    let mut b = [0u64; 12];
    for k in 0..3 {
        b[4 * k] = c.r[k].x.to_bits();
        b[4 * k + 1] = c.r[k].y.to_bits();
        b[4 * k + 2] = c.v[k].x.to_bits();
        b[4 * k + 3] = c.v[k].y.to_bits();
    }
    b
}

/// **The primary measurement.** How many of these states are actually distinct?
///
/// A path that has collapsed reports 1. Read this before any divergence figure: two collapsed
/// sets agree perfectly, and their agreement means nothing.
pub fn distinct(states: &[Cart<f64>]) -> usize {
    let mut v: Vec<[u64; 12]> = states.iter().map(bits).collect();
    v.sort_unstable();
    v.dedup();
    v.len()
}

/// Largest absolute component difference between two states.
pub fn max_abs_diff(a: &Cart<f64>, b: &Cart<f64>) -> f64 {
    let mut m = 0.0f64;
    for k in 0..3 {
        m = m
            .max((a.r[k].x - b.r[k].x).abs())
            .max((a.r[k].y - b.r[k].y).abs())
            .max((a.v[k].x - b.v[k].x).abs())
            .max((a.v[k].y - b.v[k].y).abs());
    }
    m
}

/// The `N` values of `delta` on one axis, matching `Slice::axis`'s endpoint-inclusive grid.
pub fn deltas(n: usize) -> Vec<f64> {
    (0..n).map(|i| if n <= 1 { -1.0 } else { -1.0 + 2.0 * i as f64 / (n - 1) as f64 }).collect()
}
