//! The production colouring: **hue from the shape sphere by vMF site-blend, lightness from a
//! scalar with a declared polarity, and a reserved colour for everything undetermined.**
//!
//! This replaces `bivariate.rs`, which had three faults that its own outputs demonstrate.
//!
//! # Fault one: the old hue map is 2-to-1, and the seam it was defended against never existed
//!
//! The old map was `hue = atan2(n[2], n[1])`, `chroma = C_MAX * sqrt(n1^2 + n2^2)`, fed into
//! OKLab as `(a, b) = chroma * (cos hue, sin hue)`. **That composition is algebraically just
//! `(a, b) = C_MAX * (n1, n2)`** — measured agreement `4.2e-17` over a sphere sweep. So it is a
//! linear orthogonal projection of `S^2` onto a disc: perfectly continuous, and the seam its doc
//! comment defended against was never there. That defence was transcribed from a warning about a
//! *different* construction — mapping a hue angle through a colour wheel, which does wrap.
//!
//! The real fault is a different one. The projection **discards `n0` entirely**, so it is exactly
//! 2-to-1: `n` and its `n0 -> -n0` partner render **bitwise identically**. And `n0` is
//! `(a - b) / I` with `a = |rho~|^2`, `b = |lam~|^2`, so flipping it exchanges the inner-pair
//! separation with the outer one. A tight binary with a distant third body and a wide pair with a
//! close third body were painted the same colour.
//!
//! **How much that cost, measured rather than asserted — and it is less than the above implies.**
//! `n0` is *reached* end to end (span `1.9946` in `near-field`, `1.9994` in `deep interior`,
//! against a maximum of 2) but its interdecile is only `0.0665` and `0.1684`. Span is a max
//! statistic; the interdecile says the bulk sits in a sliver. So the merge bit in the **tail**,
//! not the bulk, and **the flat images were the ramp, not the hue map.** Two independent faults,
//! and it is worth being clear which one produced the picture.
//!
//! The construction on record instead places the shape sphere's landmarks as **von Mises–Fisher
//! poles** and blends between them. It reads all three components, so it separates those pairs;
//! `tests/colour.rs` asserts the separation with the old map as the negative control.
//!
//! ```text
//! d_i = n . p_i                       cosine similarity to site i
//! w_i = exp(kappa * (d_i - d_max))    d_max subtracted for stability, NOT a parameter
//! w   = w / sum(w)
//! (a, b) = sum_i w_i * (a_i, b_i)     weighted mean in OKLab
//! ```
//!
//! Blending in OKLab rather than sRGB is what makes the desaturation *mean* something: a
//! weighted mean of two chromas shrinks toward the midpoint, so a direction sitting between two
//! sites reads as greyer. **Uncertainty about which regime you are in reads as greyness**, which
//! is the property the space was chosen for, not a side effect.
//!
//! [`hue_ab`] is continuous everywhere on `S^2` because it is a smooth function of the dot
//! products, and `tests/colour.rs` sweeps a great circle asserting the OKLab path is continuous.
//! That arm carries a genuinely discontinuous negative control (an angle through a hue wheel) so
//! it is a test rather than a formality — the `atan2` map would have passed it.
//!
//! # Fault two: a linear ramp over four decades
//!
//! `ensemble_spread` in `near-field` spans `(4.19e-5, 0.286)` with a median near `1e-3`. Under
//! `(v - lo)/(hi - lo)` essentially every pixel lands at `L_MIN`, and the committed
//! `colour_near-field_bivariate_spread_reference.png` is a flat navy field. This is the failure
//! already recorded in `RESULTS.md` — *"a linear ramp hid the structure entirely"* — reintroduced
//! in code that took the percentile half of the lesson and dropped the log half. Here the curve
//! is a property of the field ([`Scalar::curve`]), not a call-site choice.
//!
//! **And the window was stretched by a different estimator.** `ensemble_spread` is
//! `max(spread_shape, spread_event)`. The event arm is a count ratio over `E+1` copies and takes
//! **5 distinct values** in `near-field` (modal `98.2%`), but it dominates only `1.7%` of
//! footprints — all of them in the top tail. So it sets the p99 and nothing else: the window ran
//! to `2.857e-1` (exactly `2/7`) where the continuous arm's own p99 is `2.244e-2`, **12.7x
//! narrower**. A linear ramp over a window an order of magnitude too wide, set by a staircase
//! that describes 1.7% of the region. That is why [`Scalar::ShapeSpread`] exists as a separate
//! field and why [`quantisation`] and [`event_arm_fraction`] are printed before any image.
//!
//! # Fault three: undetermined pixels rendered as plausible physics
//!
//! A NaN `shape_vec` — a triple collision, which `shape.rs` deliberately does not floor — drove
//! `oklab_to_srgb` through `NaN.round() as u8`, which saturates to `[0, 0, 0]`, and
//! `metric::Cache::render` filled background with `0u8`. **An undetermined pixel was bitwise
//! identical to un-rendered background, and the metric scored it as a perfect match.** Every
//! such case now returns [`DEBUG_NAN`], and [`BACKGROUND`] is deliberately not black.
//!
//! # Polarity is a property of the field, not an argument
//!
//! [`Scalar::direction`] is a method with no parameter, so there is no call site at which the
//! sense of the ramp can be got backwards. Get it wrong once and the image is a photographic
//! negative of the physics while still looking entirely plausible. The rule, applied everywhere:
//!
//! ```text
//! HIGH uncertainty / divergence / chaos  ->  BRIGHT
//! LOW  uncertainty / settled  / regular  ->  DARK
//! ```
//!
//! The eye goes to the bright regions, and the regions that need attention are the uncertain
//! ones. `t_end`, `d_min` and `terminated_fraction` run the other way — a *longer*-lived or
//! *wider*-passing trajectory is the *more* settled one — so they carry
//! [`Direction::HighIsSettled`] and are inverted before the curve is applied.

use crate::ensemble::pixel::PixelOut;
use crate::outcome::State;
use crate::output::oklab;
use crate::physics::newton::accel;
use crate::{Vec2};

/// Undetermined: non-finite shape, non-finite scalar, a non-finite copy, a failed decode, or a
/// collapsed footprint. **Never black, never white, never interpolated.** A pixel with no value
/// must be visibly null rather than plausibly real.
pub const DEBUG_NAN: [u8; 3] = [255, 0, 255];

/// Un-rendered background. Deliberately not `[0,0,0]` and not [`DEBUG_NAN`]: "nothing was drawn
/// here", "this was drawn and is undetermined" and "this was drawn and is dark" are three
/// different statements and must be three different colours.
pub const BACKGROUND: [u8; 3] = [18, 18, 22];

/// Maximum chroma carried by a site colour.
pub const C_MAX: f64 = 0.13;
/// Lightness range. Not `[0, 1]`: pure black and pure white carry no chroma, so hue would vanish
/// at both ends of the ramp and the bivariate map would silently become univariate.
pub const L_MIN: f64 = 0.30;
pub const L_MAX: f64 = 0.92;

/// Default vMF concentration. Low blends smoothly, high tends to a hard Voronoi partition of the
/// sphere. A design choice to sweep and state, never to tune until a picture looks right.
pub const KAPPA: f64 = 3.0;

// ---------------------------------------------------------------------------------------------
// Sites
// ---------------------------------------------------------------------------------------------

/// One vMF pole: a direction on the shape sphere and the OKLab colour placed there.
#[derive(Clone, Copy, Debug)]
pub struct Site {
    pub n: [f64; 3],
    pub lab: [f64; 3],
    pub name: &'static str,
}

/// An ordered set of sites and the concentration to blend them with.
#[derive(Clone, Debug)]
pub struct SiteSet {
    pub sites: Vec<Site>,
    pub kappa: f64,
}

/// Normalised vMF weights for a vector of cosine similarities.
///
/// Split out as a pure function of `d` so the invariance claim is **testable**: subtracting the
/// maximum is numerical conditioning, not a free parameter, and `tests/colour.rs` asserts that
/// adding any constant to every `d_i` leaves the weights bitwise unchanged. Asserted rather than
/// stated, because a shift that did move the weights would move every colour in the image by a
/// small amount and look like nothing at all.
pub fn vmf_weights(d: &[f64], kappa: f64) -> Option<Vec<f64>> {
    if d.is_empty() || !d.iter().all(|x| x.is_finite()) {
        return None;
    }
    let dmax = d.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = d.iter().map(|x| (kappa * (x - dmax)).exp()).collect();
    let z: f64 = w.iter().sum();
    if !(z > 0.0) || !z.is_finite() {
        return None;
    }
    Some(w.into_iter().map(|x| x / z).collect())
}

/// OKLab `(a, b)` for a direction on the shape sphere.
///
/// Returns `None` for a non-finite direction, so the caller renders [`DEBUG_NAN`] rather than
/// letting NaN propagate into a `u8` cast and come out as black.
pub fn hue_ab(set: &SiteSet, n: [f64; 3]) -> Option<(f64, f64)> {
    if !n.iter().all(|x| x.is_finite()) || set.sites.is_empty() {
        return None;
    }
    let d: Vec<f64> = set
        .sites
        .iter()
        .map(|s| n[0] * s.n[0] + n[1] * s.n[1] + n[2] * s.n[2])
        .collect();
    let w = vmf_weights(&d, set.kappa)?;
    let mut a = 0.0;
    let mut b = 0.0;
    for (wi, s) in w.iter().zip(set.sites.iter()) {
        a += wi * s.lab[1];
        b += wi * s.lab[2];
    }
    Some((a, b))
}

fn site_lab(hue_deg: f64, chroma: f64) -> [f64; 3] {
    let h = hue_deg.to_radians();
    // L is carried only so the triple is a well-formed OKLab colour; the production combiner is
    // Replace-L, so the site's own lightness is discarded and the scalar's substituted.
    [0.65, chroma * h.cos(), chroma * h.sin()]
}

/// The distinguished points of the shape sphere, **computed for the given masses**.
///
/// Six sites: the three binary-collision singularities, the two Lagrange (equilateral) central
/// configurations, and the `lambda = 0` degeneracy. They are derived by running
/// [`crate::physics::shape::shape_vec`] on configurations built here, never hard-coded — their
/// coordinates depend on the mass ratios, and masses become a chart coordinate in the latent and
/// mass-simplex charts, so a hard-coded site set would be silently wrong on exactly the charts
/// that were added to vary them.
///
/// **The sixth site was added on a measurement, not for symmetry.** With five, the worst angular
/// gap to the nearest site was `1.193 rad` and it sat at `n0 = +1`: `b = |lam~|^2 -> 0`, the third
/// body at the inner pair's barycentre. That is the *antipode* of the `(0,1)` collision — the
/// `n0` axis runs from "bodies 0 and 1 coincident" to "body 2 at their barycentre" — so leaving
/// it uncovered left a whole pole of the sphere blending to one wash. `near-field` and
/// `deep interior` use that axis end to end (measured `n0` span `1.9946` and `1.9994` against a
/// maximum of 2), so the hole was in the part of the sphere the data actually visits.
///
/// The three collisions get a maximally separated hue triad: they are the singular points and
/// the ones a reader most needs to locate. The Lagrange configurations are regular, and get
/// reduced chroma so a regular region does not shout.
///
/// The [`euler_points`] configurations also lie on the collinear great circle, between adjacent
/// collision points. They are computed and available as overlays and as a test, but are **not**
/// sites: adding a site whose colour is the interpolation the blend already produces localises
/// nothing and only makes the palette harder to read.
pub fn landmarks(m: &[f64; 3]) -> SiteSet {
    use crate::physics::shape::shape_vec;

    let at = |r: [Vec2<f64>; 3]| shape_vec(&r, m);
    let o = Vec2::new(0.0, 0.0);
    let e = Vec2::new(1.0, 0.0);

    let coll01 = at([o, o, e]);
    let coll02 = at([o, e, o]);
    let coll12 = at([e, o, o]);

    let s3 = 3f64.sqrt() / 2.0;
    let lag_p = at([o, e, Vec2::new(0.5, s3)]);
    let lag_m = at([o, e, Vec2::new(0.5, -s3)]);

    // `lambda = 0`: body 2 at the inner pair's barycentre. Built so `com01` is the origin for
    // these masses, so the site moves with them like every other one.
    let lam0 = at([Vec2::new(-m[1], 0.0), Vec2::new(m[0], 0.0), o]);

    SiteSet {
        sites: vec![
            Site { n: coll01, lab: site_lab(30.0, C_MAX), name: "collision(0,1)" },
            Site { n: coll02, lab: site_lab(150.0, C_MAX), name: "collision(0,2)" },
            Site { n: coll12, lab: site_lab(270.0, C_MAX), name: "collision(1,2)" },
            Site { n: lam0, lab: site_lab(210.0, C_MAX), name: "lambda=0" },
            Site { n: lag_p, lab: site_lab(90.0, C_MAX * 0.45), name: "lagrange(+)" },
            Site { n: lag_m, lab: site_lab(330.0, C_MAX * 0.45), name: "lagrange(-)" },
        ],
        kappa: KAPPA,
    }
}

/// Spread of the blended OKLab `(a, b)` over a set of shape directions: how much work the hue
/// channel is actually doing on this data.
///
/// Returned as the mean distance from the centroid in `(a, b)`, the same shape of statistic
/// `spread_shape` uses. A near-zero value means every pixel got the same hue — which is a
/// truthful answer when a region visits one kind of configuration, and a broken palette when it
/// does not. Read it beside the region's `n0` span, never alone.
pub fn hue_coverage(set: &SiteSet, ns: &[[f64; 3]]) -> f64 {
    let pts: Vec<(f64, f64)> = ns.iter().filter_map(|&n| hue_ab(set, n)).collect();
    if pts.len() < 2 {
        return f64::NAN;
    }
    let k = pts.len() as f64;
    let (ca, cb) = (
        pts.iter().map(|p| p.0).sum::<f64>() / k,
        pts.iter().map(|p| p.1).sum::<f64>() / k,
    );
    pts.iter().map(|p| ((p.0 - ca).powi(2) + (p.1 - cb).powi(2)).sqrt()).sum::<f64>() / k
}

/// The three Euler collinear central configurations, as shape-sphere directions.
///
/// `euler_points(m)[k]` is the configuration with **body `k` between the other two**.
///
/// Found by bisection on the middle body's position against the central-configuration condition
/// `a_i = -lambda (r_i - R_com)` with one `lambda` shared by all three, evaluated through the
/// already-validated [`accel`]. Solving it this way rather than transcribing the Euler quintic
/// is deliberate: the quintic is exactly the kind of algebra that fails silently, and this form
/// is checkable against a function the project already trusts. `tests/colour.rs` asserts the
/// residual, which is what makes the bisection a measurement rather than an assertion of faith.
pub fn euler_points(m: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut out = [[f64::NAN; 3]; 3];
    for k in 0..3 {
        let (i, j) = match k {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        // Ends pinned at -1 and +1; the middle body slides between them.
        let build = |x: f64| {
            let mut r = [Vec2::zero(); 3];
            r[i] = Vec2::new(-1.0, 0.0);
            r[j] = Vec2::new(1.0, 0.0);
            r[k] = Vec2::new(x, 0.0);
            r
        };
        let resid = |x: f64| {
            let r = build(x);
            let a = accel(&r, m, 0.0);
            let mtot = m[0] + m[1] + m[2];
            let c = (r[0].x * m[0] + r[1].x * m[1] + r[2].x * m[2]) / mtot;
            a[i].x / (r[i].x - c) - a[j].x / (r[j].x - c)
        };
        let (mut lo, mut hi) = (-1.0 + 1e-9, 1.0 - 1e-9);
        let (flo, fhi) = (resid(lo), resid(hi));
        if !(flo.is_finite() && fhi.is_finite()) || flo * fhi > 0.0 {
            // No bracket: report an undetermined landmark rather than a plausible number.
            continue;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if resid(lo) * resid(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        out[k] = crate::physics::shape::shape_vec(&build(0.5 * (lo + hi)), m);
    }
    out
}

/// Largest angular gap between the shape direction and its nearest site, over a sphere sweep.
///
/// A site set that leaves a large hole gives a whole neighbourhood one flat blended colour. Used
/// by the diagnostics rather than by the render.
pub fn worst_site_gap(set: &SiteSet, samples: usize) -> f64 {
    let mut worst: f64 = 0.0;
    for i in 0..samples {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / samples as f64;
        let r = (1.0 - z * z).max(0.0).sqrt();
        // Golden-angle spiral: a deterministic near-uniform sphere sample, no RNG.
        let th = std::f64::consts::PI * (1.0 + 5f64.sqrt()) * i as f64;
        let n = [r * th.cos(), r * th.sin(), z];
        let best = set
            .sites
            .iter()
            .map(|s| n[0] * s.n[0] + n[1] * s.n[1] + n[2] * s.n[2])
            .fold(f64::NEG_INFINITY, f64::max);
        worst = worst.max(best.clamp(-1.0, 1.0).acos());
    }
    worst
}

// ---------------------------------------------------------------------------------------------
// Lightness
// ---------------------------------------------------------------------------------------------

/// Which way the field runs. Declared per field, never per call site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Direction {
    /// Higher means less settled: render brighter.
    HighIsUnstable,
    /// Higher means more settled: **invert** before the curve, then render brighter for less
    /// settled. `t_end`, `d_min` and `terminated_fraction` are the three of these.
    HighIsSettled,
}

/// How the normalised value is warped before it becomes lightness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Linear,
    /// Log between `lo` and `hi`. The right default for anything spanning decades.
    Log,
    /// Signed log, for fields that take both signs. `diffusion` is the one.
    SymLog,
    Gamma(f64),
    Sqrt,
}

/// The scalar driving lightness. Its direction and curve are properties of the variant, so
/// there is no argument through which either can be got wrong.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scalar {
    /// `max(spread_shape, spread_event)` — what the criterion currently reads.
    Spread,
    /// The continuous arm alone. Worth colouring separately: `spread_event` is a count ratio
    /// over `E+1` copies, so wherever it dominates `Spread` the field is a staircase with
    /// `E+2` levels and no ramp can recover what was never there. See [`quantisation`].
    ShapeSpread,
    EventSpread,
    Ftle,
    Diffusion,
    ErrorRatio,
    TEnd,
    DMin,
    /// `energy_drift_max` -- the **diagnostic** field, not a science field.
    ///
    /// Already in the payload per footprint, so this is a colouring and not a computation. It
    /// is what made the `dtau` blow-up visible: mapping drift directly showed coherent ARCS of
    /// high drift with the non-finite pixels sitting inside them, where the science fields only
    /// show the artefact after it has propagated into an outcome or a spread.
    Drift,
}

impl Scalar {
    pub fn name(self) -> &'static str {
        match self {
            Scalar::Spread => "spread",
            Scalar::ShapeSpread => "spread_shape",
            Scalar::EventSpread => "spread_event",
            Scalar::Ftle => "ftle",
            Scalar::Diffusion => "diffusion",
            Scalar::ErrorRatio => "error_ratio",
            Scalar::TEnd => "t_end",
            Scalar::DMin => "d_min",
            Scalar::Drift => "energy_drift_max",
        }
    }

    pub fn value(self, p: &PixelOut) -> f64 {
        match self {
            Scalar::Spread => p.ensemble_spread,
            Scalar::ShapeSpread => p.spread_shape,
            Scalar::EventSpread => p.spread_event,
            Scalar::Ftle => p.ftle,
            Scalar::Diffusion => p.diffusion,
            Scalar::ErrorRatio => p.error_ratio,
            Scalar::TEnd => p.t_end,
            Scalar::DMin => p.d_min_true,
            Scalar::Drift => p.energy_drift_max,
        }
    }

    /// Which end of the field is the unsettled one. See the module note on polarity.
    pub fn direction(self) -> Direction {
        match self {
            Scalar::TEnd | Scalar::DMin => Direction::HighIsSettled,
            _ => Direction::HighIsUnstable,
        }
    }

    pub fn curve(self) -> Curve {
        match self {
            // Spans decades with a median near 1e-3. This is the fix for the flat images.
            // Drift spans ~fifteen decades between a clean trajectory and a blown-up one.
            Scalar::Spread | Scalar::ShapeSpread | Scalar::DMin | Scalar::Drift => Curve::Log,
            // Signed: the measured ramps start negative.
            Scalar::Diffusion => Curve::SymLog,
            // A count ratio in [0,1]; a log would only stretch the quantisation.
            Scalar::EventSpread => Curve::Linear,
            Scalar::Ftle | Scalar::ErrorRatio | Scalar::TEnd => Curve::Linear,
        }
    }
}

/// `v` into `[0, 1]` against `[lo, hi]`, with the field's own curve and direction applied.
///
/// Returns `None` for a non-finite input, so the caller renders [`DEBUG_NAN`] rather than
/// silently clamping an undetermined value to one end of the ramp — which is what the old code
/// did, and which made "this could not be determined" indistinguishable from "this is the
/// quietest pixel in the region".
pub fn range_norm(s: Scalar, v: f64, lo: f64, hi: f64) -> Option<f64> {
    if !v.is_finite() || !lo.is_finite() || !hi.is_finite() || !(hi > lo) {
        return None;
    }
    let t = match s.curve() {
        Curve::Linear => (v - lo) / (hi - lo),
        Curve::Sqrt => ((v - lo) / (hi - lo)).clamp(0.0, 1.0).sqrt(),
        Curve::Gamma(g) => ((v - lo) / (hi - lo)).clamp(0.0, 1.0).powf(g),
        Curve::Log => {
            // A p1 of exactly 0 is common (spread is 0 wherever the copies agree exactly).
            // Floor it rather than returning None: the value is determined, the *window* is
            // degenerate, and those are different failures.
            let l = lo.max(1e-12);
            let h = hi.max(l * 10.0);
            ((v.max(l).ln() - l.ln()) / (h.ln() - l.ln())).clamp(0.0, 1.0)
        }
        Curve::SymLog => {
            let scale = lo.abs().max(hi.abs()).max(1e-12);
            let lin = scale * 1e-3;
            let f = |x: f64| (x / lin).abs().ln_1p() * x.signum();
            let (a, b) = (f(lo), f(hi));
            if (b - a).abs() < 1e-300 {
                return None;
            }
            (f(v) - a) / (b - a)
        }
    };
    let t = t.clamp(0.0, 1.0);
    Some(match s.direction() {
        Direction::HighIsUnstable => t,
        Direction::HighIsSettled => 1.0 - t,
    })
}

/// Normalised value to OKLab `L`, in `[L_MIN, L_MAX]`.
pub fn lightness(t: f64) -> f64 {
    L_MIN + (L_MAX - L_MIN) * t.clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------------------------
// The production map
// ---------------------------------------------------------------------------------------------

/// The colour of one footprint. Hue and chroma from the shape sphere, lightness from `s`.
///
/// Returns [`DEBUG_NAN`] for every case in which the pixel has no value: a non-finite
/// `shape_vec` (triple collision), a non-finite scalar, any non-finite copy in the ensemble
/// (`n_nonfinite > 0`), and the two failure states. Each of those was previously rendered as a
/// valid colour, three of them as the *quietest* colour on the ramp.
pub fn rgb(p: &PixelOut, s: Scalar, set: &SiteSet, lo: f64, hi: f64) -> [u8; 3] {
    if p.n_nonfinite > 0 {
        return DEBUG_NAN;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => return DEBUG_NAN,
        _ => {}
    }
    let (a, b) = match hue_ab(set, p.shape_vec) {
        Some(x) => x,
        None => return DEBUG_NAN,
    };
    let t = match range_norm(s, s.value(p), lo, hi) {
        Some(t) => t,
        None => return DEBUG_NAN,
    };
    // Replace-L: the sites' own lightness is discarded and the scalar's substituted, so the two
    // channels stay independent. Modulate-L would let a site's palette bleed into the scalar.
    oklab::oklab_to_srgb([lightness(t), a, b])
}

/// [`rgb`] **resolved over the ensemble** -- supersampling, not anti-aliasing.
///
/// The right name is **supersampling**: every sub-sample is a full simulation, coloured
/// independently, and the pixel is their mean. Anti-aliasing is what that *buys*, not what it
/// *is*. `src/output/ssaa.rs` already does exactly this for the categorical outcome palette --
/// its module doc calls it *resolve*, "an average that drives display", against
/// `ensemble_spread`'s *disagreement*, which drives scheduling. **This is the shape-sphere arm of
/// the same operation**, and it is the arm that had never been written.
///
/// **The nominal sample is included.** `copy_shapes` is built from all `outs`, so
/// `copy_shapes[0]` *is* `shape_vec` -- the nominal is one sub-sample among `E+1`, not a
/// privileged one, which is what makes this a mean over the footprint rather than a correction
/// applied to a centre point.
///
/// # The samples already exist and the production map discards them
///
/// Every footprint carries `E+1` copies jittered across the **whole cell**, edge to edge --
/// `jitter_frac = 0.5` with the fixed Halton (2,3) prefix. `spread_shape` reduces over all of
/// them and sets the **lightness**; the **hue** reads `p.shape_vec`, which is `shapes[0]`, the
/// nominal copy alone. So one channel is ensemble-averaged and the other is a single point sample
/// at the cell centre, on a field whose structure runs below pixel scale. That is the textbook
/// setup for aliasing, with 8 free samples per pixel already computed and thrown away.
///
/// # It must average in COLOUR space, not on the sphere
///
/// `hue_ab` is a **von Mises-Fisher weighted blend** of landmark colours, which is nonlinear in
/// `n`. So `hue_ab(mean(shapes)) != mean(hue_ab(shapes))` and only the second is the box filter.
/// Averaging the shape vectors first would also be actively wrong where the copies diverge: their
/// mean collapses toward the origin, chroma drops, and the pixel renders **pale** -- manufacturing
/// the very appearance the wedge investigation was about. This averages `(a, b)` in OKLab, after
/// the map, which is what supersampling means.
///
/// (`CLAUDE.md` records the shipped hue map as *linear* and identically `C_MAX*(n1,n2)`. That
/// describes the earlier `chroma*(cos h, sin h)` form, not this one. Under a linear map the two
/// orders agree; under this one they do not, and the note is stale rather than wrong.)
///
/// # What it costs and what it does not fix
///
/// `E+1` extra `hue_ab` evaluations **at render time**. No extra integration -- the copies are
/// already marched. It is live-playhead compatible: it reduces over samples that exist at the
/// current playhead time, with no lookahead and no re-integration.
///
/// **It is a sampling fix, not a physics fix.** If the banding is a *cadence* artefact -- `t_end`
/// quantised to sync boundaries -- every copy in a footprint snaps to the same boundary and
/// averaging them changes nothing. That is the discriminating prediction and it is why this ships
/// as a toggle beside `rgb` rather than replacing it.
///
/// **The sample count is not chosen for rendering.** `E+1` is set by the refinement criterion, so
/// the anti-aliasing rate is whatever that decided. Eight is a reasonable rate; one is not, and
/// with `n_extra = 0` this falls back to [`rgb`] exactly.
///
/// Requires [`crate::ensemble::pixel::EnsembleCfg::keep_copy_shapes`]. With it unset
/// `copy_shapes` is empty and this **returns [`rgb`] unchanged** rather than a silently different
/// picture.
pub fn rgb_resolved(p: &PixelOut, s: Scalar, set: &SiteSet, lo: f64, hi: f64) -> [u8; 3] {
    use crate::output::compose;
    if p.copy_shapes.len() < 2 {
        return rgb(p, s, set, lo, hi);
    }
    if p.n_nonfinite > 0 {
        return DEBUG_NAN;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => return DEBUG_NAN,
        _ => {}
    }
    let l = match range_norm(s, s.value(p), lo, hi) {
        Some(t) => lightness(t),
        None => return DEBUG_NAN,
    };
    // **One supersampler.** `compose::resolve` knows nothing about this map; it is handed a
    // per-sub-sample colour and averages. That independence is the whole point -- a `*_resolved`
    // twin per mode is two copies of one operation waiting to drift.
    let lab = compose::resolve(p.copy_shapes.len(), |i| {
        hue_ab(set, p.copy_shapes[i]).map(|(a, b)| [l, a, b])
    });
    compose::finish(lab, DEBUG_NAN)
}

/// Robust `[p1, p99]` of a scalar over a set of footprints.
///
/// Percentiles rather than min/max: one undetermined footprint at `1e12` would compress every
/// other pixel into the bottom of the range and the image would read as featureless — the same
/// failure as reading a variance where the excess kurtosis is 110.
pub fn range(px: &[PixelOut], s: Scalar) -> (f64, f64) {
    range_q(px, s, 0.01, 0.99)
}

/// [`range`] at stated quantiles. The diagnostic drift field is auto-ranged over p2-p98.
///
/// Percentiles rather than min/max, for the same reason: one undetermined footprint at `1e12`
/// would compress every other pixel into the bottom of the range.
pub fn range_q(px: &[PixelOut], s: Scalar, lo_q: f64, hi_q: f64) -> (f64, f64) {
    let mut v: Vec<f64> = px.iter().map(|p| s.value(p)).filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return (0.0, 1.0);
    }
    let lo = crate::quad::quantile(&mut v.clone(), lo_q);
    let hi = crate::quad::quantile(&mut v, hi_q);
    (lo, hi)
}

/// Inferno control points, sampled at `t = 0, 1/8, ..., 1`. Interpolated linearly in sRGB —
/// close enough for a diagnostic ramp, and it needs no colour-space machinery.
const INFERNO: [[f64; 3]; 9] = [
    [0.0, 0.0, 3.0],
    [22.0, 11.0, 57.0],
    [66.0, 10.0, 104.0],
    [106.0, 23.0, 110.0],
    [147.0, 38.0, 103.0],
    [188.0, 55.0, 84.0],
    [221.0, 81.0, 58.0],
    [243.0, 136.0, 20.0],
    [252.0, 255.0, 164.0],
];

/// The **diagnostic** colouring: `energy_drift_max` on an inferno ramp, magenta where there is
/// no value.
///
/// Univariate on purpose. [`rgb`] is bivariate (hue from the shape sphere, lightness from a
/// scalar) and a diagnostic asked under that colouring inherits the shape field's structure,
/// which is exactly what obscures a numerical defect: the science fields only show it once it
/// has propagated. **When a numerical defect is suspected, render this, not the science field.**
///
/// Returns [`DEBUG_NAN`] for the same veto set [`rgb`] applies — a non-finite value, any
/// non-finite copy in the ensemble, and the two failure states — reused rather than re-derived,
/// so "no value" means the same thing in both panels.
pub fn drift_rgb(p: &PixelOut, lo: f64, hi: f64) -> [u8; 3] {
    if p.n_nonfinite > 0 {
        return DEBUG_NAN;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => return DEBUG_NAN,
        _ => {}
    }
    let Some(t) = range_norm(Scalar::Drift, Scalar::Drift.value(p), lo, hi) else {
        return DEBUG_NAN;
    };
    let x = (t.clamp(0.0, 1.0) * 8.0).min(7.999_999);
    let i = x as usize;
    let f = x - i as f64;
    let mut out = [0u8; 3];
    for k in 0..3 {
        let c = INFERNO[i][k] * (1.0 - f) + INFERNO[i + 1][k] * f;
        out[k] = c.round().clamp(0.0, 255.0) as u8;
    }
    out
}

/// How much ordering the lightness channel actually carries: `(distinct, finite, modal_frac)`.
///
/// **Count the distinct values before reading any picture.** The committed `spread` p99 values
/// are exactly `2/7` and `6/7` — `ensemble_spread` is `max(spread_shape, spread_event)` and
/// `spread_event` is a count ratio over `E+1 = 8` copies, so wherever the event arm dominates
/// the lightness field is an eight-level staircase and no ramp recovers what is not there. The
/// same rule already caught two criteria whose flat `error(B)` curves were the tie-break's scan
/// order rather than a signal.
///
/// Values are compared by bit pattern, so this counts *exact* distinct values and never merges
/// two that a tolerance would.
pub fn quantisation(px: &[PixelOut], s: Scalar) -> (usize, usize, f64) {
    use std::collections::HashMap;
    let mut counts: HashMap<u64, usize> = HashMap::new();
    let mut finite = 0usize;
    for p in px {
        let v = s.value(p);
        if v.is_finite() {
            finite += 1;
            *counts.entry(v.to_bits()).or_insert(0) += 1;
        }
    }
    let modal = counts.values().cloned().max().unwrap_or(0);
    let frac = if finite == 0 { f64::NAN } else { modal as f64 / finite as f64 };
    (counts.len(), finite, frac)
}

/// Fraction of footprints where `ensemble_spread` is the **event** arm rather than the shape arm.
///
/// The companion to [`quantisation`]: it says *why* the field is quantised where it is. Reported
/// beside any image coloured on [`Scalar::Spread`].
pub fn event_arm_fraction(px: &[PixelOut]) -> f64 {
    let mut n = 0usize;
    let mut k = 0usize;
    for p in px {
        if p.ensemble_spread.is_finite() {
            n += 1;
            if p.spread_event >= p.spread_shape {
                k += 1;
            }
        }
    }
    if n == 0 {
        f64::NAN
    } else {
        k as f64 / n as f64
    }
}
