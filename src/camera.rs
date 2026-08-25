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

    pub fn zoom(&self, root_half: f64) -> f64 {
        root_half / self.half_world
    }

    /// The level whose quad width matches the viewport width. Fractional on purpose — the
    /// integer form is only taken where a level comparison needs it.
    pub fn depth(&self, root_half: f64) -> f64 {
        (root_half / self.half_world).log2()
    }

    /// One tile is one sample. `tile_size(quad, zoom) = quad_width * zoom / N`, in pixels.
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
