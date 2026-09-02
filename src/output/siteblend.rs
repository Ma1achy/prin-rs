//! **Mode 2 — the shape-sphere site blend.** `blend(colours, weights(kernel, {n_hat . p_i}), space)`
//!
//! The specced Family A map, factored so each axis is an explicit choice rather than an
//! implementation accident. `colour::hue_ab` is the special case `Vmf(KAPPA)` over
//! `Physics(Mixed)` in `Oklab`; this is the general form it was a hard-coded instance of.
//!
//! # Kernel = support x temperature, one primitive at different temperatures
//!
//! - `Vmf(kappa)` — all sites, `w_i = exp(kappa*(d_i - d_max))`, `kappa in [0.5, 12]`
//! - `TopK(k, ks)` — nearest `k` only, softmax over the `k` largest `d_i`, `ks in [1, 20]`
//! - `Nearest` — `w_i = [i = argmax d_i]`
//!
//! **`Nearest` is the `kappa -> inf` limit and is a DISCRETE code path, not a large `kappa`.**
//! `exp` overflows long before the limit is reached, so approaching it numerically produces
//! `NaN`, not a Voronoi map. `d_max` is subtracted for stability throughout, which shifts the
//! weights and leaves their ratios exactly as they were.
//!
//! The UI shape follows: **one blend-sharpness dial with a hard detent at the top edge** that
//! snaps to `Nearest`, and `support` — all-site vMF against top-k Voronoi — as the separate
//! discrete choice.
//!
//! # Blend space is a parameter, because the two say different things
//!
//! - `Oklab` — blend `(a, b)` at the sites' own `(L, C)`. A weighted mean of chroma **shrinks
//!   toward site boundaries**, so **uncertainty reads as desaturation**. That is a feature and
//!   the perceptual default.
//! - `Rgb` — for site colours that already came from a perceptual LUT, where a second perceptual
//!   blend would double-correct.
//!
//! # Sites are two kinds and the distinction is in the TYPE
//!
//! [`StaticGen`] is uniform-hoistable by construction. [`PhysicsGen`] depends on the masses and
//! **moves with them**: when a slice axis or tilt touches a mass dimension the masses are
//! per-pixel state and there is no per-slice constant to bake *even in principle*.
//!
//! [`Sites::hoistable_over`] is the optimisation: when neither basis axis touches a mass
//! dimension, physics sites are constant over the slice and lift to a uniform. **Semantics
//! per-pixel, cost per-slice when possible** — and the check is exact rather than heuristic,
//! because `decoder::Latent` puts `z_mu` at indices **6 and 7**.
//!
//! # Site colours are always swatches
//!
//! Hue tables are dissolved: Full-OKLAB and Okabe-Ito are preset **swatch sets**, not a separate
//! type, so a caller cannot reach a hue-table code path that a swatch set cannot express.
use crate::output::oklab;
use crate::physics::shape::shape_vec;
use crate::Vec2;

/// Support x temperature. One primitive; `Nearest` is its limit as a discrete path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kernel {
    /// All sites, `w_i = exp(kappa*(d_i - d_max))`. `kappa in [0.5, 12]`.
    Vmf(f64),
    /// Nearest `k` sites, softmax at sharpness `ks in [1, 20]` over their `d_i`.
    TopK(usize, f64),
    /// `w_i = [i = argmax d_i]`. The `kappa -> inf` limit, taken discretely because `exp`
    /// overflows before it.
    Nearest,
}

impl Kernel {
    /// The dial position, `[0, 1]`, with `1.0` the hard detent that snaps to [`Kernel::Nearest`].
    pub fn from_dial(t: f64, top_k: Option<usize>) -> Kernel {
        if t >= 1.0 {
            return Kernel::Nearest;
        }
        let t = t.clamp(0.0, 1.0);
        match top_k {
            Some(k) => Kernel::TopK(k, 1.0 + 19.0 * t),
            None => Kernel::Vmf(0.5 + 11.5 * t),
        }
    }

    /// Weights over `d`, normalised to sum 1. `None` if `d` is empty or non-finite throughout —
    /// never a uniform fallback, which would render an undetermined pixel as a valid blend.
    pub fn weights(self, d: &[f64]) -> Option<Vec<f64>> {
        if d.is_empty() || d.iter().any(|x| !x.is_finite()) {
            return None;
        }
        let dmax = d.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if !dmax.is_finite() {
            return None;
        }
        match self {
            Kernel::Nearest => {
                // Discrete. No `exp` is evaluated at all, which is the whole point of the variant.
                let i = d.iter().position(|&x| x == dmax)?;
                let mut w = vec![0.0; d.len()];
                w[i] = 1.0;
                Some(w)
            }
            Kernel::Vmf(kappa) => {
                let mut w: Vec<f64> = d.iter().map(|&x| (kappa * (x - dmax)).exp()).collect();
                let s: f64 = w.iter().sum();
                if !(s > 0.0) {
                    return None;
                }
                for x in w.iter_mut() {
                    *x /= s;
                }
                Some(w)
            }
            Kernel::TopK(k, ks) => {
                let k = k.max(1).min(d.len());
                let mut idx: Vec<usize> = (0..d.len()).collect();
                idx.sort_by(|&a, &b| d[b].partial_cmp(&d[a]).unwrap());
                let keep = &idx[..k];
                let mut w = vec![0.0; d.len()];
                let mut s = 0.0;
                for &i in keep {
                    let e = (ks * (d[i] - dmax)).exp();
                    w[i] = e;
                    s += e;
                }
                if !(s > 0.0) {
                    return None;
                }
                for x in w.iter_mut() {
                    *x /= s;
                }
                Some(w)
            }
        }
    }
}

/// Where the weighted mean is taken. **An explicit parameter, not an implementation accident.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlendSpace {
    /// OKLab `(a, b)`. Weighted-mean chroma shrinks toward site boundaries, so **uncertainty
    /// reads as desaturation**. The perceptual default.
    Oklab,
    /// Linear in sRGB, for site colours already drawn from a perceptual LUT.
    Rgb,
}

/// Static site generators. Uniform-hoistable by construction: they do not depend on the masses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StaticGen {
    Axes6,
    Corner8,
    Ico12,
    Fib(usize),
    /// `N` points on a great circle, tilted by `tilt` and rotated by `rot`, both radians.
    Ring(usize, f64, f64),
}

/// Physics site generators. **Move with the masses**, so they are per-pixel state whenever a
/// slice axis touches a mass dimension.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhysicsGen {
    /// The three binary-collision singularities.
    Bc,
    /// The three Euler (collinear) configurations.
    Euler,
    /// The two Lagrange (equilateral) configurations.
    Lagrange,
}

/// The two kinds, distinguished in the type so the hoist decision is a match and not a guess.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Sites {
    Static(StaticGen),
    Physics(PhysicsGen),
}

impl Sites {
    /// Unit directions on the shape sphere. `m` is ignored by [`Sites::Static`].
    pub fn points(self, m: &[f64; 3]) -> Vec<[f64; 3]> {
        match self {
            Sites::Static(g) => static_points(g),
            Sites::Physics(g) => physics_points(g, m),
        }
    }

    /// **The hoist.** `true` when the sites are constant over the slice, so they can be computed
    /// once instead of per pixel. Static sites always; physics sites only when no basis axis
    /// touches a mass dimension.
    ///
    /// `q1`/`q2` are the slice's basis axes in `decoder::Latent` coordinates, where `z_mu` sits
    /// at indices **6 and 7**. Exact, not heuristic. A caller with a chart that varies masses by
    /// some other route must pass `false`.
    pub fn hoistable_over(self, q1: &[f64; 8], q2: &[f64; 8]) -> bool {
        match self {
            Sites::Static(_) => true,
            Sites::Physics(_) => {
                const MASS_DIMS: [usize; 2] = [6, 7];
                MASS_DIMS.iter().all(|&i| q1[i] == 0.0 && q2[i] == 0.0)
            }
        }
    }
}

fn norm(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n > 0.0 { [v[0] / n, v[1] / n, v[2] / n] } else { [f64::NAN; 3] }
}

fn static_points(g: StaticGen) -> Vec<[f64; 3]> {
    match g {
        StaticGen::Axes6 => vec![
            [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0], [0.0, 0.0, -1.0],
        ],
        StaticGen::Corner8 => {
            let mut v = Vec::with_capacity(8);
            for sx in [-1.0f64, 1.0] {
                for sy in [-1.0f64, 1.0] {
                    for sz in [-1.0f64, 1.0] {
                        v.push(norm([sx, sy, sz]));
                    }
                }
            }
            v
        }
        StaticGen::Ico12 => {
            let p = (1.0 + 5f64.sqrt()) / 2.0;
            let mut v = Vec::with_capacity(12);
            for s1 in [-1.0f64, 1.0] {
                for s2 in [-1.0f64, 1.0] {
                    v.push(norm([0.0, s1, s2 * p]));
                    v.push(norm([s1, s2 * p, 0.0]));
                    v.push(norm([s2 * p, 0.0, s1]));
                }
            }
            v
        }
        StaticGen::Fib(n) => {
            let n = n.max(1);
            let ga = std::f64::consts::PI * (3.0 - 5f64.sqrt());
            (0..n)
                .map(|i| {
                    let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
                    let r = (1.0 - z * z).max(0.0).sqrt();
                    let th = ga * i as f64;
                    [r * th.cos(), r * th.sin(), z]
                })
                .collect()
        }
        StaticGen::Ring(n, tilt, rot) => {
            let n = n.max(1);
            (0..n)
                .map(|i| {
                    let a = std::f64::consts::TAU * i as f64 / n as f64 + rot;
                    // A great circle in the xy-plane, tilted about x.
                    let (x, y, z) = (a.cos(), a.sin(), 0.0f64);
                    let (c, s) = (tilt.cos(), tilt.sin());
                    [x, y * c - z * s, y * s + z * c]
                })
                .collect()
        }
    }
}

fn physics_points(g: PhysicsGen, m: &[f64; 3]) -> Vec<[f64; 3]> {
    let at = |r: [Vec2<f64>; 3]| shape_vec(&r, m);
    let o = Vec2::new(0.0, 0.0);
    let e = Vec2::new(1.0, 0.0);
    match g {
        PhysicsGen::Bc => vec![at([o, o, e]), at([o, e, o]), at([e, o, o])],
        PhysicsGen::Lagrange => {
            let s3 = 3f64.sqrt() / 2.0;
            vec![at([o, e, Vec2::new(0.5, s3)]), at([o, e, Vec2::new(0.5, -s3)])]
        }
        PhysicsGen::Euler => crate::output::colour::euler_points(m).to_vec(),
    }
}

/// The full Mode 2 map. Every axis explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct SiteBlend {
    pub sites: Sites,
    pub kernel: Kernel,
    /// One OKLab swatch per site. **Always a swatch** — hue tables are dissolved into preset
    /// swatch sets rather than surviving as a distinct type.
    pub colours: Vec<[f64; 3]>,
    pub space: BlendSpace,
}

impl SiteBlend {
    /// `n_hat` -> OKLab. `None` when the direction or the weights are undetermined, so the caller
    /// renders the explicit invalid colour rather than inheriting a valid one.
    pub fn blend(&self, n_hat: [f64; 3], m: &[f64; 3]) -> Option<[f64; 3]> {
        if !n_hat.iter().all(|x| x.is_finite()) {
            return None;
        }
        let pts = self.sites.points(m);
        if pts.is_empty() || pts.len() != self.colours.len() {
            return None;
        }
        let d: Vec<f64> = pts
            .iter()
            .map(|p| n_hat[0] * p[0] + n_hat[1] * p[1] + n_hat[2] * p[2])
            .collect();
        let w = self.kernel.weights(&d)?;
        match self.space {
            BlendSpace::Oklab => {
                let mut out = [0.0f64; 3];
                for (wi, c) in w.iter().zip(&self.colours) {
                    for k in 0..3 {
                        out[k] += wi * c[k];
                    }
                }
                Some(out)
            }
            BlendSpace::Rgb => {
                let mut acc = [0.0f64; 3];
                for (wi, c) in w.iter().zip(&self.colours) {
                    let s = oklab::oklab_to_srgb(*c);
                    for k in 0..3 {
                        acc[k] += wi * s[k] as f64;
                    }
                }
                Some(oklab::srgb_to_oklab([
                    acc[0].round().clamp(0.0, 255.0) as u8,
                    acc[1].round().clamp(0.0, 255.0) as u8,
                    acc[2].round().clamp(0.0, 255.0) as u8,
                ]))
            }
        }
    }

    /// The provenance line. Emitted with every render, per the spec's requirement that both modes
    /// declare mode, kernel, temperature, blend space, site set, brightness field and polarity,
    /// and ramp window. The brightness half is appended by the caller, which owns that choice.
    pub fn provenance(&self) -> String {
        let k = match self.kernel {
            Kernel::Vmf(x) => format!("vmf(kappa={x})"),
            Kernel::TopK(k, ks) => format!("topk(k={k}, ks={ks})"),
            Kernel::Nearest => "nearest".into(),
        };
        let s = match self.sites {
            Sites::Static(g) => format!("static({g:?})"),
            Sites::Physics(g) => format!("physics({g:?})"),
        };
        format!("mode=siteblend kernel={k} space={:?} sites={s} n_sites={}", self.space, self.colours.len())
    }
}

// -------------------------------------------------------------------------------------------
// The brightness channel
// -------------------------------------------------------------------------------------------

use crate::ensemble::pixel::PixelOut;
use crate::output::colour::{self, Scalar};

/// **The brightness slot, as a first-class occupant.** Any §1.2 scalar field can fill it.
///
/// `combine(colour = site blend, brightness = <field> greyscale, Replace-L)`. The three named
/// defaults are `Ftle` (white = high, chaos pops), `Diffusion` (white = high, spreading pops) and
/// `TEnd` (white = low / early, quick-resolving pops and bounded stays black) -- but they are
/// **not hardcoded**: the slot takes a [`Scalar`], so a field added later occupies it with no
/// change here.
///
/// **Magnitude fields default to greyscale rather than a sequential LUT**, because a greyscale
/// magnitude field *is* a lightness and therefore doubles as the brightness occupant with no
/// conversion. A sequential LUT would have to be inverted back out to be used here.
///
/// **Polarity is per-field and deliberately not uniform.** It lives on the [`Scalar`] variant as
/// `Direction`, so it cannot be got wrong through an argument: `Ftle` is `HighIsUnstable`, `TEnd`
/// is `HighIsSettled`, and the salient end is bright in both.
///
/// **The window is fixed and shareable.** Auto-ranging per panel manufactures or hides the very
/// difference it is meant to show -- measured on this project, not hypothesised, where `far` read
/// `error(root) = 0.60` under an auto-ranged ramp whose window was `(1.3e-9, 1.1e-8)`. Auto-range
/// is opt-in via [`Brightness::auto`] and **prints its window**.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Brightness {
    pub field: Scalar,
    pub lo: f64,
    pub hi: f64,
    /// True when `lo`/`hi` came from the data rather than from a constant. Carried so the
    /// provenance line can say so.
    pub auto_ranged: bool,
}

impl Brightness {
    /// A fixed, shareable window. The default: two renders made this way are comparable.
    pub fn fixed(field: Scalar, lo: f64, hi: f64) -> Self {
        Self { field, lo, hi, auto_ranged: false }
    }

    /// Opt-in auto-range over this panel's own p1..p99. Comparable with nothing else, which is
    /// why it is not the default and why `provenance` announces it.
    pub fn auto(field: Scalar, px: &[PixelOut]) -> Self {
        let (lo, hi) = colour::range(px, field);
        Self { field, lo, hi, auto_ranged: true }
    }

    /// OKLab `L`, or `None` for an undetermined value.
    pub fn l(&self, p: &PixelOut) -> Option<f64> {
        colour::range_norm(self.field, self.field.value(p), self.lo, self.hi)
            .map(colour::lightness)
    }

    pub fn provenance(&self) -> String {
        format!(
            "brightness={:?} polarity={:?} curve={:?} window=({:e},{:e}){}",
            self.field,
            self.field.direction(),
            self.field.curve(),
            self.lo,
            self.hi,
            if self.auto_ranged { " AUTO-RANGED (comparable with nothing else)" } else { "" }
        )
    }
}

/// Mode 2's colour for one footprint: site blend for hue and chroma, `b` for lightness.
///
/// **Replace-L**, so the two channels stay independent: the sites' own `L` is discarded and the
/// brightness field's substituted. Modulate-L would let a site's palette bleed into the scalar.
///
/// Returns `inv` -- the **explicit invalid colour** -- for a non-finite shape direction, an
/// undetermined brightness value, any non-finite copy in the ensemble, and the two failure
/// states. *Do not let an invalid pixel inherit a valid colour by accident.*
pub fn rgb(
    p: &PixelOut,
    blend: &SiteBlend,
    b: &Brightness,
    m: &[f64; 3],
    inv: [u8; 3],
) -> [u8; 3] {
    use crate::outcome::State;
    if p.n_nonfinite > 0 {
        return inv;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => return inv,
        _ => {}
    }
    let lab = match blend.blend(p.shape_vec, m) {
        Some(x) => x,
        None => return inv,
    };
    let l = match b.l(p) {
        Some(x) => x,
        None => return inv,
    };
    oklab::oklab_to_srgb([l, lab[1], lab[2]])
}

/// [`rgb`] **resolved over the ensemble** -- supersampling. Every sub-sample is a full
/// simulation, blended independently; the pixel is the mean of their colours.
///
/// The mean is taken **after** the site blend, never on the sphere before it: the blend is
/// nonlinear in `n_hat`, so the two orders differ, and averaging directions first would collapse
/// diverging copies toward the origin and render them pale -- manufacturing the appearance the
/// wedge investigation was about.
///
/// Requires `keep_copy_shapes`; without it this returns [`rgb`] unchanged.
pub fn rgb_resolved(
    p: &PixelOut,
    blend: &SiteBlend,
    b: &Brightness,
    m: &[f64; 3],
    inv: [u8; 3],
) -> [u8; 3] {
    use crate::outcome::State;
    use crate::output::compose;
    if p.copy_shapes.len() < 2 {
        return rgb(p, blend, b, m, inv);
    }
    if p.n_nonfinite > 0 {
        return inv;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => return inv,
        _ => {}
    }
    // The brightness field is **per-pixel**, not per-copy: the scalars are footprint reductions
    // and there are no per-copy values to average. So the colour channel is supersampled and the
    // lightness is not, and that asymmetry is stated rather than hidden.
    let l = match b.l(p) {
        Some(x) => x,
        None => return inv,
    };
    let lab = compose::resolve(p.copy_shapes.len(), |i| {
        blend.blend(p.copy_shapes[i], m).map(|c| [l, c[1], c[2]])
    });
    compose::finish(lab, inv)
}
