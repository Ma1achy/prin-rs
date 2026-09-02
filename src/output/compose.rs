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
