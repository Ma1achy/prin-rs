//! The initial-condition slice: one body's position varies over a 2D box, the other two
//! stay fixed (BRIEF §2.2).
//!
//! Row order matches `tb.burrau_grid`: `meshgrid(indexing='xy')` then C-order `ravel()`,
//! i.e. **`index = jy*nx + jx`, x fastest**. The cross-check compares row by row, so this
//! ordering is load-bearing and is asserted in a test against a hardcoded case.
//!
//! Initial conditions are always built in `f64` and cast down. Generating them separately
//! per precision would make an IC difference indistinguishable from a genuine f32
//! arithmetic effect — the decomposition that made the earlier f32 investigation
//! interpretable at all.

use crate::physics::{burrau, Cart};
use crate::{Real, Vec2};

#[derive(Clone, Copy, Debug)]
pub struct Slice {
    pub nx: usize,
    pub ny: usize,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    /// Which body's position varies: 0, 1 or 2.
    pub body: usize,
}

impl Slice {
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

    /// Grid position of pixel `idx`, with `idx = jy*nx + jx`.
    pub fn decode_pos(&self, idx: usize) -> (f64, f64) {
        let jx = idx % self.nx;
        let jy = idx / self.nx;
        (
            Self::axis(self.cx, self.half, self.nx, jx),
            Self::axis(self.cy, self.half, self.ny, jy),
        )
    }

    /// Nominal (un-jittered) initial condition for pixel `idx`.
    ///
    /// This is copy 0 of the reference's ensemble — `mask[::reps] = False` leaves it
    /// un-jittered and therefore completely seed-independent. That is what makes a
    /// nominal-only cross-check possible with no RNG on either side.
    pub fn nominal<T: Real>(&self, idx: usize) -> Cart<T> {
        let (x, y) = self.decode_pos(idx);
        let mut s = burrau::state::<f64>();
        s.r[self.body] = Vec2::new(x, y);
        s.cast::<T>()
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
    REGIONS.iter().find(|r| r.0 == name).map(|&(_, cx, cy, body)| Slice {
        nx,
        ny,
        cx,
        cy,
        half,
        body,
    })
}
