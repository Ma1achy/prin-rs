//! **The missing middle space** — §12.
//!
//! The deep-zoom spec nests three coordinate systems and this codebase has only had two.
//!
//! | space | what it is | who reads it |
//! |---|---|---|
//! | **screen** | pixels, `[0, res)²` | [`crate::camera::Camera::to_px`], the renderers |
//! | **UV** | quad addressing, `[0, 1]²` against a **fixed** frame | this module |
//! | **chart** | initial-condition coordinates, whatever the chart's units are | [`crate::grid`] |
//!
//! `Slice::axis` is a `linspace` over `[cx - half, cx + half]` in **chart** units and the tree's
//! root half is `0.05` or `3.0`. So the spec's *"`h` is an exact power of two at every depth"* is
//! not a claim about this repo's arithmetic — it is a claim about a space this repo did not have.
//! Chart-space cell widths satisfy the weaker `w == prev / 2.0` exactly (halving is exact) and are
//! not powers of two; UV widths are exactly `2^-(d+1)`, which is what makes an integer address
//! meaningful. `tests/uv.rs` asserts both, side by side, because the distinction is the whole
//! reason this file exists.
//!
//! # The frame is FIXED, and that is what makes an address absolute
//!
//! The caching contract's `QuadID (z, tx, ty)` is *"a Mandelbrot-style absolute address"*. Absolute
//! against what? Not the tree root — [`crate::quad::QuadTree::grow_root`] doubles the root box, and
//! an address relative to it would renumber every quad in the store on a zoom-out, which is exactly
//! the cost re-rooting exists to avoid. So the frame is fixed at session start and never moves, and
//! a grown root simply addresses **negative depth and negative indices**. Refusing to name that
//! area would make the address relative after all.
//!
//! # The half-cell trap, which has already cost this project a measurement
//!
//! Recovering `(depth, tx, ty)` from a box is a division and a round, and
//! `metric::Cache::key_of` carries a paragraph about getting it wrong: a centre sits at
//! `(2i + 1) h`, so dividing by the **cell width** `2h` gives `i + 0.5` and `.round()` lands on
//! `i + 1` — mapping every quad to its right/upper neighbour, and scoring a perfectly coherent
//! leaf set belonging to a shifted tree. [`UvFrame::id_of`] therefore checks the recovered address
//! back against the box it came from and returns `None` rather than a plausible wrong answer.

/// The fixed reference box, in **chart** coordinates. UV is this box mapped to `[0, 1]²`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UvFrame {
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
}

/// An absolute quad address in UV space.
///
/// The cell is `[tx * 2^-d, (tx + 1) * 2^-d] × [ty * 2^-d, (ty + 1) * 2^-d]`.
///
/// **`depth` is signed and the indices are signed**, for the regrow case above: depth `-1` is a box
/// twice the frame, and `tx = -1` is the cell immediately left of it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct QuadId {
    pub depth: i32,
    pub tx: i64,
    pub ty: i64,
}

/// Half-width of a depth-`d` cell in UV, **exactly** `2^-(d+1)`.
///
/// `exp2` of an integer is exact at f64 for `|d| <= 1023`, so this is a power of two in the strict
/// sense: `uv_half(d).log2()` is an integer and the mantissa is 1.
pub fn uv_half(depth: i32) -> f64 {
    f64::from(-(depth + 1)).exp2()
}

impl QuadId {
    /// Centre in UV: `(2*tx + 1) * 2^-(d+1)`. Exact while `|2 tx + 1| < 2^53`.
    pub fn centre_uv(&self) -> (f64, f64) {
        let h = uv_half(self.depth);
        (((2 * self.tx + 1) as f64) * h, ((2 * self.ty + 1) as f64) * h)
    }

    /// The four children, in [`crate::quad::Quad::child_boxes`] order:
    /// lower-left, lower-right, upper-left, upper-right.
    pub fn children(&self) -> [QuadId; 4] {
        let (d, x, y) = (self.depth + 1, self.tx * 2, self.ty * 2);
        [
            QuadId { depth: d, tx: x, ty: y },
            QuadId { depth: d, tx: x + 1, ty: y },
            QuadId { depth: d, tx: x, ty: y + 1 },
            QuadId { depth: d, tx: x + 1, ty: y + 1 },
        ]
    }

    /// The parent. Uses a floor division, so it is right for negative indices too — `-1 / 2` is
    /// `0` in Rust and the cell left of the frame would claim the frame's own parent.
    pub fn parent(&self) -> QuadId {
        QuadId {
            depth: self.depth - 1,
            tx: self.tx.div_euclid(2),
            ty: self.ty.div_euclid(2),
        }
    }

    /// A stable per-quad seed. **Depends on the address and nothing else** — not on the tree, not
    /// on the root box, not on iteration order — so a quad re-reached after a re-root, an eviction
    /// or a fresh session draws the same numbers. That is the property `Scheme::Halton`'s fixed
    /// prefix already has per copy index, at the level of the quad.
    pub fn seed(&self, salt: u64) -> u64 {
        // SplitMix64's finaliser over the three components, mixed one at a time.
        let mut h = salt;
        for v in [self.depth as i64 as u64, self.tx as u64, self.ty as u64] {
            h ^= v.wrapping_add(0x9E37_79B9_7F4A_7C15).wrapping_add(h << 6).wrapping_add(h >> 2);
            h ^= h >> 30;
            h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
            h ^= h >> 27;
            h = h.wrapping_mul(0x94D0_49BB_1331_11EB);
            h ^= h >> 31;
        }
        h
    }
}

impl UvFrame {
    pub fn new(cx: f64, cy: f64, half: f64) -> UvFrame {
        assert!(half > 0.0, "a UV frame needs a positive half-width, got {half}");
        UvFrame { cx, cy, half }
    }

    /// Chart → UV. The frame's box maps to `[0, 1]²`; outside it the values simply leave `[0, 1]`,
    /// which is what a grown root needs.
    pub fn to_uv(&self, x: f64, y: f64) -> (f64, f64) {
        let w = 2.0 * self.half;
        ((x - self.cx) / w + 0.5, (y - self.cy) / w + 0.5)
    }

    /// UV → chart.
    pub fn to_chart(&self, u: f64, v: f64) -> (f64, f64) {
        let w = 2.0 * self.half;
        (self.cx + (u - 0.5) * w, self.cy + (v - 0.5) * w)
    }

    /// The chart box of an address: `(cx, cy, half)`, matching [`crate::quad::Quad`]'s fields.
    pub fn box_of(&self, id: QuadId) -> (f64, f64, f64) {
        let h = uv_half(id.depth);
        let (u, v) = id.centre_uv();
        let (x, y) = self.to_chart(u, v);
        (x, y, h * 2.0 * self.half)
    }

    /// **The address of a chart box, or `None` if it is not on the lattice.**
    ///
    /// Two ways to be off-lattice and both are refusals rather than roundings: a half-width that is
    /// not `frame.half * 2^-d` for integer `d`, and a centre that does not sit at an odd multiple
    /// of that width. The second is the half-cell trap — the recovered index is checked back
    /// against the centre it came from, at a tolerance of a thousandth of a cell, because a
    /// silently-shifted address is a *coherent wrong answer* and this project has scored one.
    pub fn id_of(&self, cx: f64, cy: f64, half: f64) -> Option<QuadId> {
        if !(half > 0.0) || !cx.is_finite() || !cy.is_finite() {
            return None;
        }
        // depth from the width. `log2` of an exact ratio of powers of two is exact.
        let ratio = half / self.half; // = 2^-d
        let d = -ratio.log2();
        let depth = d.round();
        if (d - depth).abs() > 1e-9 || !(-1023.0..=1023.0).contains(&depth) {
            return None;
        }
        let depth = depth as i32;
        let h = uv_half(depth);
        let (u, v) = self.to_uv(cx, cy);
        // The centre is `(2t + 1) h`, so `u / h` is `2t + 1`: **an ODD integer**. Solving for `t`
        // through `u / (2h) - 0.5` is the form that lands on `i + 0.5` and rounds the wrong way.
        let (fx, fy) = (u / h, v / h);
        let (ox, oy) = (fx.round(), fy.round());
        if (fx - ox).abs() > 1e-6 || (fy - oy).abs() > 1e-6 {
            return None;
        }
        if (ox as i64).rem_euclid(2) != 1 || (oy as i64).rem_euclid(2) != 1 {
            return None;
        }
        let id = QuadId { depth, tx: (ox as i64 - 1) / 2, ty: (oy as i64 - 1) / 2 };
        // The guard that earns its place: recover the box and require it back.
        let (bx, by, bh) = self.box_of(id);
        let tol = bh * 1e-3;
        if (bx - cx).abs() > tol || (by - cy).abs() > tol || (bh - half).abs() > tol {
            return None;
        }
        Some(id)
    }
}

/// **Quad-local sample coordinates: `du, dv ∈ [-1, 1]`, formed DIRECTLY.**
///
/// The defect §12 names is not that the codebase lacks these — [`crate::decode::sample`] already
/// takes them — it is that `ensemble::jitter` **recovers** `du` as `(u - cx) / half` from a `u` the
/// global `linspace` had just formed. Precision is spent building the offset and then the offset is
/// subtracted back off, which is free at f64 on an O(1) chart coordinate and total at f32 or at
/// depth 40.
///
/// `axis_local(n, i)` is `-1 + 2i/(n-1)` with the last sample pinned to `+1`, mirroring
/// `Slice::axis`'s endpoint handling. **It is not bitwise `(axis(c,h,n,i) - c)/h`** and does not
/// claim to be: the two orders of operation round differently, which is precisely why one of them
/// carries information the other has already thrown away. `examples/sample_space.rs` measures the
/// gap rather than asserting it.
pub fn axis_local(n: usize, i: usize) -> f64 {
    if n <= 1 {
        return -1.0;
    }
    if i == n - 1 {
        return 1.0;
    }
    if i == 0 {
        return -1.0;
    }
    -1.0 + 2.0 * (i as f64) / ((n - 1) as f64)
}
