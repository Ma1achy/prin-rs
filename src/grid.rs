//! The initial-condition slice: a 2D box of chart coordinates, decoded to full states.
//!
//! Row order matches `tb.burrau_grid`: `meshgrid(indexing='xy')` then C-order `ravel()`,
//! i.e. **`index = jy*nx + jx`, x fastest**. The cross-check compares row by row, so this
//! ordering is load-bearing and is asserted in a test against a hardcoded case.
//!
//! Initial conditions are always built in `f64` and cast down. Generating them separately
//! per precision would make an IC difference indistinguishable from a genuine f32
//! arithmetic effect — the decomposition that made the earlier f32 investigation
//! interpretable at all.
//!
//! **The chart is a parameter.** Every experiment before the vertical slice varied one
//! body's position in the plane, which is an *affine* decode: `J_D` is constant, so the
//! linearised path `x = x0 + J_D.delta` is exact rather than approximate and "where does
//! the linearisation start to matter" answers "never" at every depth. [`Chart::Shape`]
//! exists so that question has something to measure. [`Chart::BodyPlane`] is the old
//! behaviour, preserved **bitwise** and asserted so in `tests/vertical_slice.rs`.

use crate::physics::{burrau, shape, Cart, Ic};
use crate::{Real, Vec2};

/// How a chart coordinate `(u, v)` becomes a three-body state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chart {
    /// One body's position over a 2D box; the other two fixed at their Burrau values.
    /// **Affine** — this is every result before the vertical slice.
    BodyPlane,
    /// An oblique 2-plane in the 6D position space: `origin + u*U + v*V`, where `U` and `V`
    /// may move more than one body. Still affine, but no longer axis-aligned, which is what
    /// §3.5 asks about. [`Chart::BodyPlane`] is the special case
    /// `U = e(body, x)`, `V = e(body, y)`.
    Plane {
        origin: Cart<f64>,
        u: [Vec2<f64>; 3],
        v: [Vec2<f64>; 3],
    },
    /// The shape sphere: `(u, v)` is a tangent step from `n0`, mapped by the exponential
    /// map and inverted through the Hopf map at fixed scale and fibre phase.
    ///
    /// **Nonlinear.** The exp map's `cos|t|`, `sin|t|/|t|` is where the curvature lives, and
    /// it is the only chart here on which the linearised decoder is an approximation at all.
    Shape {
        n0: [f64; 3],
        e1: [f64; 3],
        e2: [f64; 3],
        inertia: f64,
        phase: f64,
    },
}

impl Default for Chart {
    fn default() -> Self {
        Chart::BodyPlane
    }
}

impl Chart {
    pub fn name(&self) -> &'static str {
        match self {
            Chart::BodyPlane => "body_plane",
            Chart::Plane { .. } => "plane",
            Chart::Shape { .. } => "shape",
        }
    }

    /// Is the decode affine? If so the linearised path is exact and its curvature term is
    /// **structurally zero** — a fact to report as structural, never as a measurement.
    pub fn is_affine(&self) -> bool {
        !matches!(self, Chart::Shape { .. })
    }

    /// The `Plane` that reproduces `BodyPlane` for `body`, used to assert the two agree.
    pub fn plane_for_body(body: usize) -> Chart {
        let mut u = [Vec2::zero(); 3];
        let mut v = [Vec2::zero(); 3];
        u[body] = Vec2::new(1.0, 0.0);
        v[body] = Vec2::new(0.0, 1.0);
        // The origin's varying body sits at zero; `u`, `v` carry the whole position.
        let mut origin = burrau::state::<f64>();
        origin.r[body] = Vec2::zero();
        Chart::Plane { origin, u, v }
    }

    /// The shape chart through Burrau's own configuration: same point on the sphere, same
    /// scale, deterministic tangent frame. A slice of it is a slice of the shape sphere
    /// around the reference triangle.
    pub fn shape_at_burrau(phase: f64) -> Chart {
        let m = burrau::MASSES;
        let r: [Vec2<f64>; 3] = [
            Vec2::new(burrau::R0[0][0], burrau::R0[0][1]),
            Vec2::new(burrau::R0[1][0], burrau::R0[1][1]),
            Vec2::new(burrau::R0[2][0], burrau::R0[2][1]),
        ];
        let n0 = shape::shape_vec(&r, &m);
        let (e1, e2) = shape::tangent_frame(n0);
        Chart::Shape { n0, e1, e2, inertia: shape::inertia(&r, &m), phase }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Slice {
    pub nx: usize,
    pub ny: usize,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    /// Which body's position varies under [`Chart::BodyPlane`]: 0, 1 or 2. Retained on every
    /// chart because it names the slice family in dumps and headers.
    pub body: usize,
    pub chart: Chart,
}

impl Slice {
    /// The pre-vertical-slice constructor. Every existing result is a `BodyPlane` slice and
    /// this is the only way one is built, so nothing can acquire a different chart silently.
    pub fn body_plane(nx: usize, ny: usize, cx: f64, cy: f64, half: f64, body: usize) -> Self {
        Slice { nx, ny, cx, cy, half, body, chart: Chart::BodyPlane }
    }

    pub fn with_chart(mut self, chart: Chart) -> Self {
        self.chart = chart;
        self
    }

    pub fn npix(&self) -> usize {
        self.nx * self.ny
    }

    /// Cell widths, **per axis**.
    ///
    /// The reference computes only `hx = 2*half/max(nx-1,1)` and uses it for *both* axes.
    /// Every experiment run against it used a square grid, so it never bit; on any
    /// non-square grid the y-jitter is silently mis-scaled. Fixed here.
    pub fn cell_widths(&self) -> (f64, f64) {
        let hx = 2.0 * self.half / (self.nx.max(2) - 1) as f64;
        let hy = 2.0 * self.half / (self.ny.max(2) - 1) as f64;
        (hx, hy)
    }

    /// One axis of the grid, matching `numpy.linspace` including its exact endpoint.
    fn axis(c: f64, half: f64, n: usize, i: usize) -> f64 {
        let (a, b) = (c - half, c + half);
        if n <= 1 {
            return a;
        }
        if i == n - 1 {
            return b; // numpy sets the last sample exactly, rather than accumulating
        }
        a + (i as f64) * ((b - a) / (n - 1) as f64)
    }

    /// Chart position of pixel `idx`, with `idx = jy*nx + jx`.
    pub fn decode_pos(&self, idx: usize) -> (f64, f64) {
        let jx = idx % self.nx;
        let jy = idx / self.nx;
        (
            Self::axis(self.cx, self.half, self.nx, jx),
            Self::axis(self.cy, self.half, self.ny, jy),
        )
    }

    /// **The decode**: chart coordinate `(u, v)` to a full state, in `f64`.
    ///
    /// Public because the deep-zoom ladder decodes directly, without a grid.
    pub fn decode_state(&self, u: f64, v: f64) -> Ic<f64> {
        decode_state(&self.chart, self.body, u, v)
    }

    /// Nominal (un-jittered) initial condition for pixel `idx`.
    ///
    /// This is copy 0 of the reference's ensemble — `mask[::reps] = False` leaves it
    /// un-jittered and therefore completely seed-independent. That is what makes a
    /// nominal-only cross-check possible with no RNG on either side.
    pub fn nominal<T: Real>(&self, idx: usize) -> Cart<T> {
        let (x, y) = self.decode_pos(idx);
        self.decode_state(x, y).s.cast::<T>()
    }

    /// Nominal initial condition **with its masses**.
    ///
    /// [`Self::nominal`] drops them, which is correct at the many call sites that predate the
    /// chart families and are all Burrau. Anything that decodes a chart which can vary mass must
    /// use this instead, or it integrates the right configuration with the wrong bodies.
    pub fn nominal_ic<T: Real>(&self, idx: usize) -> Ic<T> {
        let (x, y) = self.decode_pos(idx);
        self.decode_state(x, y).cast::<T>()
    }
}

/// The decode, free of any grid. `body` is read only by [`Chart::BodyPlane`].
pub fn decode_state(chart: &Chart, body: usize, u: f64, v: f64) -> Ic<f64> {
    match *chart {
        // Bitwise what this function did before the chart existed. Asserted in a test.
        Chart::BodyPlane => {
            let mut s = burrau::state::<f64>();
            s.r[body] = Vec2::new(u, v);
            Ic { m: burrau::MASSES, s }
        }
        Chart::Plane { origin, u: uu, v: vv } => {
            let mut s = origin;
            for k in 0..3 {
                s.r[k] = s.r[k] + uu[k] * u + vv[k] * v;
            }
            Ic { m: burrau::MASSES, s }
        }
        Chart::Shape { n0, e1, e2, inertia, phase } => {
            let t = [
                u * e1[0] + v * e2[0],
                u * e1[1] + v * e2[1],
                u * e1[2] + v * e2[2],
            ];
            let n = shape::exp_map(n0, t);
            let r = shape::from_shape(n, inertia, phase, &burrau::MASSES);
            // Released from rest, like every configuration in this project.
            Ic { m: burrau::MASSES, s: Cart { r, v: [Vec2::zero(); 3] } }
        }
    }
}

/// BRIEF §2.2's named regions, all Burrau: `(name, cx, cy, body)`.
pub const REGIONS: [(&str, f64, f64, usize); 8] = [
    ("near-field", 1.0, 3.0, 0),
    ("mid-field", 1.0, 6.0, 0),
    ("far", 1.0, 13.0, 0),
    ("body2 core", 1.0, -1.0, 2),
    ("body2 mid", 1.0, -5.0, 2),
    ("body1 slice", -2.0, -1.0, 1),
    ("body1 far", -2.0, -7.0, 1),
    // Pathological: drives all three bodies together. Not regularisable; expected to hit
    // the triple-collision outcome. That is correct behaviour, not a bug (BRIEF §2.6).
    ("deep interior", 0.0, 0.0, 0),
];

pub fn region(name: &str, nx: usize, ny: usize, half: f64) -> Option<Slice> {
    REGIONS
        .iter()
        .find(|r| r.0 == name)
        .map(|&(_, cx, cy, body)| Slice::body_plane(nx, ny, cx, cy, half, body))
}
