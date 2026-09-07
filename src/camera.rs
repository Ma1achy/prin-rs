//! The camera, and the screen floor — **the everyday refinement stop**.
//!
//! The scheduler contract: *"Once a quad's tiles have shrunk to pixel size, splitting further
//! produces sub-pixel samples that cannot be displayed distinctly — wasted compute by
//! definition. This is the everyday refinement stop: in normal exploration you hit it far
//! shallower than any precision floor."*
//!
//! PR #11 had no camera, so it could not apply this, and its descent reached **level 12** in
//! near-field. The arithmetic: one sample, one tile, no interpolation, so a fully-refined tree
//! at level `L` holds `4^L * N^2` samples; at `N = 8` on a 512² viewport that reaches 262144
//! samples at **L = 6**. Everything past level 6 was 4096x beyond anything displayable. So
//! PR #11's q1, q2, q3 and q7 measured the criterion *minus its principal stop condition*.
//!
//! **Three properties, all load-bearing:**
//!
//! 1. **A veto, never a trigger.** The tile-to-pixel ratio is never itself a reason to refine.
//!    The contract strikes the "below screen resolution -> refine" rule explicitly. Enforced
//!    structurally here: [`Camera::veto`] returns only `None` or a *stopping* decision, and it
//!    is the only place the scheduler reads a camera. There is no path from tile size to
//!    `Decision::Split`.
//! 2. **View-relative, and never cached as a quad fact.** A quad floored at one zoom must
//!    refine again when zoomed into, with real new samples. Nothing here is stored on a `Quad`.
//! 3. **`MAX_REL_DEPTH` replaces absolute `max_level`.** The absolute form caps infinite zoom
//!    at ~14; the real predicate is `level < camera_depth + MAX_REL_DEPTH`. It is view-relative
//!    *scheduler* state and never on the sim key, so lowering it while zoomed invalidates no
//!    payload — it just stops scheduling deeper.

use crate::quad::{Decision, Quad};

/// Sensible range is 4–8 below the view (scheduler contract). 6 matches the screen floor at
/// `N = 8` on a 512² viewport, which is the configuration everything here is measured in.
pub const MAX_REL_DEPTH: u32 = 6;

/// Gaussian width of the fovea, in units of the viewport half-width. A quarter of the viewport,
/// so the falloff is gentle across the visible field rather than a spot.
pub const FOVEA_SIGMA: f64 = 0.25;

/// **Where the pointer is, and how settled it is.** Lives beside the camera — never on a `Quad`,
/// for the same reason camera state does not.
///
/// `dwell` is the low-pass: `0.0` while the pointer is moving fast, rising toward `1.0` as it
/// settles. It is supplied rather than computed here because the smoothing window is a property of
/// the input loop, and a `Camera` has no clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cursor {
    pub cx: f64,
    pub cy: f64,
    pub dwell: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    pub cx: f64,
    pub cy: f64,
    /// Half-width of the visible world box. Zoom is `root_half / half_world`.
    pub half_world: f64,
    /// Viewport edge, in pixels. Square.
    pub viewport: usize,
    /// `None` disables the relative-depth cap; the screen floor still applies.
    pub max_rel_depth: Option<u32>,
}

impl Camera {
    /// The whole root box, filling the viewport. `camera_depth = 0`.
    pub fn framing(cx: f64, cy: f64, root_half: f64, viewport: usize) -> Self {
        Camera { cx, cy, half_world: root_half, viewport, max_rel_depth: Some(MAX_REL_DEPTH) }
    }

    /// Zoomed in by `2^depth` about the same centre, i.e. framing a level-`depth` quad.
    pub fn at_depth(cx: f64, cy: f64, root_half: f64, viewport: usize, depth: u32) -> Self {
        Camera {
            cx,
            cy,
            half_world: root_half / (2f64).powi(depth as i32),
            viewport,
            max_rel_depth: Some(MAX_REL_DEPTH),
        }
    }

    /// World units per screen pixel.
    pub fn pixel_size(&self) -> f64 {
        2.0 * self.half_world / self.viewport as f64
    }

    /// **The one world-to-pixel projection**, for a raster of `res` pixels centred on the camera.
    ///
    /// Row 0 is the **minimum** `y`: pixel `y` grows with world `y`, exactly as `Slice` index
    /// order does, so an adaptive render, a wireframe and a `Slice` buffer written through
    /// `save_rect` agree row for row. The adaptive render and the wire used to carry three
    /// private copies of this closure that flipped `y`, and every uniform panel did not; the
    /// panels sat beside each other as mirror images. `res` is the raster and may differ from
    /// `viewport`, which sets the scale — the screen floor is a property of the viewport, the
    /// image size is a property of the file.
    pub fn to_px(&self, res: usize, x: f64, y: f64) -> (f64, f64) {
        let px = self.pixel_size();
        ((x - self.cx) / px + res as f64 / 2.0, (y - self.cy) / px + res as f64 / 2.0)
    }

    pub fn zoom(&self, root_half: f64) -> f64 {
        root_half / self.half_world
    }

    /// The level whose quad width matches the viewport width. Fractional on purpose — the
    /// integer form is only taken where a level comparison needs it.
    pub fn depth(&self, root_half: f64) -> f64 {
        (root_half / self.half_world).log2()
    }

    /// One tile is one sample. `tile_size(quad, zoom) = quad_width * zoom / N`, in pixels.
    ///
    /// **Nominal, and 14.3% smaller than what is painted at `N = 8`.** `Slice::axis` is
    /// endpoint-inclusive, so the samples sit `2h/(N-1)` apart and the adaptive render paints
    /// cells of that width; this reads `2h/N`. Kept, because moving the floor to the painted
    /// width would push the everyday stop one level deeper at the standard viewport (level 6 to
    /// 7 at `N = 8` on 512²) — a regime change to every committed tree, not a rendering fix —
    /// and because a floor that fires slightly early is the conservative direction for a veto.
    pub fn tile_size_px(&self, q: &Quad, n: usize) -> f64 {
        (2.0 * q.half / n as f64) / self.pixel_size()
    }

    /// Have this quad's tiles shrunk to pixel size?
    pub fn screen_floor(&self, q: &Quad, n: usize) -> bool {
        self.tile_size_px(q, n) <= 1.0
    }

    /// **The whole camera-scheduler interface.** Returns a *stopping* decision or nothing.
    ///
    /// That this returns `Option<Decision>` and never `Decision::Split` is the structural
    /// enforcement of "the screen floor is a veto, complexity the sole trigger". A test asserts
    /// the return value is one of the two stopping variants.
    /// Does this quad's box intersect the viewport?
    ///
    /// **A measurement, not a policy.** It is deliberately *not* consulted by [`Self::veto`]:
    /// adding view culling to the scheduler would be a caching/eviction decision, and this
    /// build has neither. It exists so §7 can report what *would* be evictable under a pan
    /// without anything being evicted.
    ///
    /// It also names a real property of the current camera: `veto` reads `tile_size_px`, which
    /// depends on the quad's width and the camera's `half_world` and `viewport` — **and not on
    /// `cx`/`cy` at all**. So panning changes no scheduling decision whatsoever today, and a
    /// pan study that measured tree persistence without saying so would be reporting "the tree
    /// persists perfectly" as a result when it is an identity.
    pub fn covers(&self, cx: f64, cy: f64, half: f64) -> bool {
        (cx - self.cx).abs() <= self.half_world + half
            && (cy - self.cy).abs() <= self.half_world + half
    }

    /// **How much of this quad the viewer can actually see**, in `[0, 1]`: the fraction of its
    /// area inside the viewport box, widened by `margin` quad-widths.
    ///
    /// # This is a PRIORITY term and never a veto term, and the distinction is load-bearing
    ///
    /// [`Self::veto`] reads `tile_size_px`, `half_world` and `viewport` — and **not `cx`/`cy`**.
    /// That is deliberate: a quad's *decision* must not depend on where the camera points, or a
    /// pan would silently invalidate the tree, which is what "never cached as a quad fact"
    /// exists to prevent. Ranking may depend on it precisely because nothing about it is stored:
    /// it is recomputed at query time on every frame, and no `Quad` gains a camera field.
    ///
    /// It also names the reason the committed pan sequence measured nothing. Nine steps moving
    /// `cx` 0.95 -> 1.05 produced a byte-identical tree, which was read as "the tree persists
    /// perfectly" — an identity, since no scheduling term read `cx` at all. Setting
    /// `max_rel_depth` would not have changed that; **this** is the term that makes a pan mean
    /// something.
    ///
    /// # The destination model
    ///
    /// `margin` is the honest baseline of §4.3's three options: widen the viewport and drop
    /// prediction entirely. Velocity extrapolation fails on flick-and-stop — it prefetches past
    /// where the user lands — and the swept-path variant must beat this before its complexity is
    /// justified. Shipping the baseline first is what makes that comparison possible.
    pub fn relevance(&self, cx: f64, cy: f64, half: f64, margin: f64) -> f64 {
        let w = self.half_world + margin * 2.0 * half;
        let overlap = |c: f64, cc: f64| {
            let (lo, hi) = ((c - half).max(cc - w), (c + half).min(cc + w));
            (hi - lo).max(0.0)
        };
        let (ox, oy) = (overlap(cx, self.cx), overlap(cy, self.cy));
        let area = (2.0 * half) * (2.0 * half);
        if area <= 0.0 {
            0.0
        } else {
            (ox * oy / area).clamp(0.0, 1.0)
        }
    }

    /// **Cursor bias (§18): where attention is, modulating relevance — never a third factor.**
    ///
    /// Returns a factor in `[1/cap, 1]` that multiplies [`Self::relevance`]. It is not added
    /// beside it and it is not a separate priority term: §18 is explicit that the pointer
    /// *modulates* camera relevance, so everything §4.3 says still holds, including that this
    /// lives in **priority and never in veto** and that no cursor field goes on a `Quad`.
    ///
    /// **The fallback is the default path.** No cursor gives `1.0` everywhere and the whole
    /// mechanism vanishes, which is exactly §2's uniform-over-viewport behaviour. Keyboard
    /// navigation, touch after the finger lifts, an unfocused window and every headless render
    /// take that path — so it is the primary case, and foveation is the modulation on top. The
    /// reverse arrangement makes every non-mouse path a special case.
    ///
    /// **A smooth falloff, never a hard radius.** A hard edge is a disc of sharpness that moves
    /// with the mouse, which reads worse than no foveation at all. The kernel is Gaussian in
    /// screen-space distance from the cursor, in units of the viewport half-width.
    ///
    /// **Weighted by dwell.** A cursor still for a moment is a far stronger signal than one flying
    /// across the canvas, and a fovea that *chases* a moving pointer spends the whole budget on
    /// regions the user has already left. At `dwell = 0` this returns `1.0` — the mechanism is off
    /// during motion, not merely weakened.
    ///
    /// **The periphery is slowed and never starved.** `cap` bounds the ratio: at `cap = 4` the far
    /// edge keeps a quarter of its priority. This is *attention* bias and not acuity exploitation
    /// — unlike VR the viewer can look away without moving the mouse, and a fovea tight enough
    /// that the edges never resolve makes the image look broken the moment they do.
    pub fn foveation(&self, cx: f64, cy: f64, cur: &Cursor, cap: f64) -> f64 {
        if cap <= 1.0 || cur.dwell <= 0.0 || self.half_world <= 0.0 {
            return 1.0;
        }
        // Distance in units of the viewport half-width, so the kernel is resolution-independent.
        let (dx, dy) = ((cx - cur.cx) / self.half_world, (cy - cur.cy) / self.half_world);
        let d2 = dx * dx + dy * dy;
        let w = (-d2 / (2.0 * FOVEA_SIGMA * FOVEA_SIGMA)).exp();
        let floor = 1.0 / cap;
        // dwell 0 -> 1.0 everywhere; dwell 1 -> 1.0 at the cursor falling to `floor` far away.
        1.0 - cur.dwell.clamp(0.0, 1.0) * (1.0 - floor) * (1.0 - w)
    }

    pub fn veto(&self, q: &Quad, n: usize, root_half: f64) -> Option<Decision> {
        if self.screen_floor(q, n) {
            return Some(Decision::ScreenFloor);
        }
        if let Some(m) = self.max_rel_depth {
            // `camera_depth` floors: a quad is capped relative to the coarsest level the view
            // fully contains, so zooming in always buys depth rather than losing it.
            let cam = self.depth(root_half).max(0.0).floor() as u32;
            if q.level >= cam + m {
                return Some(Decision::MaxRelDepth);
            }
        }
        None
    }
}
