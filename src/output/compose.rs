//! **The backbone.** `combine(colour: Option, brightness: Option)` + supersampling, §4.1.
//!
//! Two-tier mutability by design: the backbone is a **fixed typed topology with switchable
//! occupants**. Occupants are switched off, never structurally deleted, so the topology stays
//! legible and `Option<Occupant>` carries removal as data.
//!
//! # `None` is the identity element of `combine`
//!
//! | colour | brightness | Replace-L | Multiply | meaning |
//! |---|---|---|---|---|
//! | C | B | `OKLab(L=B, Ca, Cb)` | `C · B` | normal |
//! | C | None | `C` (its own L kept) | `C · 1` | just the colour map (default) |
//! | None | B | `OKLab(B, 0, 0)` | `white · B` | greyscale of the brightness field |
//! | None | None | flat mid-grey `OKLab(0.6, 0, 0)` | flat mid-grey | well-defined and visible |
//!
//! The `None`/`B` cell is why magnitude fields default to greyscale rather than a sequential LUT:
//! a greyscale magnitude field **is** a lightness, so the same field drops into either slot with
//! no conversion.
//!
//! # Supersampling is independent of the colouring
//!
//! [`resolve`] takes a **per-sub-sample colour function** and averages. It knows nothing about
//! site blends, palettes, ramps or combinators, and it is the only place the averaging exists —
//! there is exactly one supersampler and every map gets it by passing a closure. Writing a
//! `*_resolved` twin per mode, which is what this replaces, is how two copies of one operation
//! drift apart.
//!
//! *Adding a new map is wiring, not a new pixel function*, and that has to hold for the resolve
//! path too or the claim is only about the forward path.
use crate::output::oklab;

/// Mid-grey, the `None`/`None` cell. Well-defined, harmless, and instantly visible as "nothing is
/// wired here" rather than looking like data.
pub const NEUTRAL: [f64; 3] = [0.6, 0.0, 0.0];

/// The default invalid colour: a conspicuous, out-of-gamut-adjacent magenta, overridable per node.
///
/// **Validity is not optional.** Every field carries a validity lane and every ramp an explicit
/// invalid colour; without it a debug view lies at exactly the pixels it exists to expose, because
/// a NaN would ramp to *some* colour and look like data.
pub const INVALID: [u8; 3] = [0xFF, 0x00, 0xFF];

/// How the two occupants are combined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Combiner {
    /// Substitute the brightness as OKLab `L`, keeping the colour's `(a, b)`. Keeps the two
    /// channels independent; `Multiply` lets a site's palette bleed into the scalar.
    ReplaceL,
    Multiply,
}

/// `combine`, with `None` as the identity element in both slots.
///
/// `colour` is OKLab; `brightness` is a lightness already through its compaction.
pub fn combine(colour: Option<[f64; 3]>, brightness: Option<f64>, c: Combiner) -> [f64; 3] {
    match (colour, brightness, c) {
        (Some(col), Some(l), Combiner::ReplaceL) => [l, col[1], col[2]],
        (Some(col), None, _) => col,
        (None, Some(l), Combiner::ReplaceL) => [l, 0.0, 0.0],
        (Some(col), Some(l), Combiner::Multiply) => [col[0] * l, col[1] * l, col[2] * l],
        // `white * B`: white is L = 1 with no chroma, so the product is the greyscale again.
        (None, Some(l), Combiner::Multiply) => [l, 0.0, 0.0],
        (None, None, _) => NEUTRAL,
    }
}

/// **Supersample.** `n` sub-sample colours in OKLab, one pixel colour, and nothing about how they
/// were produced.
///
/// Every sub-sample is a full simulation coloured independently; the pixel is their mean. This is
/// **supersampling** — anti-aliasing is what it buys, not what it is.
///
/// **One undetermined sub-sample makes the pixel undetermined**, and that is deliberate. Dropping
/// it and averaging the survivors is the no-discard violation this project has now fixed at three
/// separate sites: a copy fails because its integration was hard, integration is hard at a close
/// encounter, and close encounters are what the instrument exists to measure — so discarding
/// biases the picture toward the tame exactly where it matters.
///
/// **The mean is taken after the map, never before it.** The maps are nonlinear in `n_hat`, so the
/// two orders differ; and averaging directions first collapses diverging copies toward the origin,
/// dropping their chroma and rendering them pale — manufacturing the very appearance the wedge
/// investigation was about.
pub fn resolve<F>(n: usize, sample: F) -> Option<[f64; 3]>
where
    F: Fn(usize) -> Option<[f64; 3]>,
{
    if n == 0 {
        return None;
    }
    let mut acc = [0.0f64; 3];
    for i in 0..n {
        let c = sample(i)?;
        for k in 0..3 {
            acc[k] += c[k];
        }
    }
    let d = n as f64;
    Some([acc[0] / d, acc[1] / d, acc[2] / d])
}

/// Finish: OKLab to sRGB, or the explicit invalid colour.
pub fn finish(lab: Option<[f64; 3]>, invalid: [u8; 3]) -> [u8; 3] {
    match lab {
        Some(c) => oklab::oklab_to_srgb(c),
        None => invalid,
    }
}

// -------------------------------------------------------------------------------------------
// Display stage (§4.3) — terminal, settings not nodes
// -------------------------------------------------------------------------------------------

/// **Sub-LSB dither before 8-bit quantisation.** The cure for contour banding on a smooth field.
///
/// A gradient whose float value moves less than `1/255` per pixel quantises to **plateaus**, and
/// the plateau edges read as contour lines. Measured on `config_stability`, one row of a red
/// ribbon: **18 of 254 adjacent pairs render to identical bytes while all 18 have a different
/// underlying float** — the field moves and the output does not. In flatter regions the plateaus
/// run 5 px and wider.
///
/// **This is NOT the ribbon banding, and saying it was is a claim this project went on to
/// refute.** The bands in `config_stability`'s ribbons are the bound pair's orbital phase winding
/// through IC space — measured on the FLOAT field, invariant to a 3.07x step change and to every
/// sampler variant, and present in seven single-trajectory observables at one orientation
/// (`results/osc/README.md`). Quantisation plateaus are a real and separate display artefact, and
/// this function is for those. *A remedy aimed at the wrong cause is still a remedy for
/// something*, but the write-up has to say which.
///
/// # This is not the dither that would have been cheating
///
/// A per-pixel dither of the **step phase** would decorrelate integration error — converting a
/// coherent artefact into incoherent noise of the same amplitude, removing the evidence rather
/// than the error. This is the opposite case and the distinction is the whole point: the error
/// being removed here is **purely a display encoding artefact**. The float is exact, the 8-bit
/// grid is what cannot represent it, and dithering trades an artificial contour for sub-LSB noise
/// that carries the true value in its local mean. Nothing about the data is hidden, because
/// nothing about the data was wrong.
///
/// # Deterministic, so renders stay reproducible
///
/// The offset is a hash of `(x, y)`, not an RNG: two runs of the same render are byte-identical,
/// which this project requires of every committed artefact. Triangular PDF over `[-1, 1]` LSB —
/// the standard choice, because a uniform dither leaves the quantisation error correlated with
/// the signal and a triangular one does not.
///
/// Applied at the display stage, so it is a **setting and not a node**: the codegen never sees
/// it and it is not part of any occupant or preset.
pub fn dither_lsb(lab: [f64; 3], x: usize, y: usize, amount: f64) -> [f64; 3] {
    // Two decorrelated hashes -> triangular PDF as the difference of two uniforms.
    let h = |a: usize, b: usize, s: u64| {
        let mut v = (a as u64).wrapping_mul(0x9E3779B97F4A7C15)
            ^ (b as u64).wrapping_mul(0xC2B2AE3D27D4EB4F)
            ^ s;
        v ^= v >> 29;
        v = v.wrapping_mul(0xBF58476D1CE4E5B9);
        v ^= v >> 32;
        (v >> 11) as f64 / (1u64 << 53) as f64
    };
    let tri = h(x, y, 0x1234_5678) - h(x, y, 0x9ABC_DEF0);
    // OKLab L is [0,1]; one 8-bit step of sRGB is ~1/255 there to within the transfer curve, and
    // `amount` is in those units. Only L is dithered: the measured banding is pure-lightness --
    // the RGB deltas across a contour are equal in all three channels.
    [lab[0] + tri * amount / 255.0, lab[1], lab[2]]
}
