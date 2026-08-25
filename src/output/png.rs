//! Two images, per BRIEF §7: outcome, and ensemble spread.
//!
//! Both are diagnostics. Anything read off them should be confirmed against the raw dump.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use crate::ensemble::pixel::PixelOut;
use crate::grid::Slice;
use crate::outcome::State;

/// Colour by `state` with `detail` shading it, per BRIEF §7. A pixel with any non-finite copy
/// is flagged separately regardless of the nominal copy's label, because "undetermined" is a
/// distinct answer from any of the others and must not be painted as one.
///
/// `detail = 3` — the two "all three" outcomes — gets the brightest shade of its family, so a
/// triple reads at a glance rather than blending into ordinary collisions or escapes.
pub fn outcome_rgb(p: &PixelOut) -> [u8; 3] {
    if p.n_nonfinite > 0 {
        return crate::output::colour::DEBUG_NAN; // deliberately loud
    }
    let base = match State::from_bits(p.state) {
        Some(State::Escape) => [220, 80, 60],
        Some(State::Collision) => [110, 190, 110],
        Some(State::Bounded) => [70, 150, 220],
        Some(State::Running) => [200, 190, 90],
        // Both failure states are undetermined and take the reserved colour. `DecodeFailed`
        // previously fell to the catch-all grey, where it was indistinguishable from an
        // invalid state byte -- a pixel whose IC could not be formed read as ordinary data.
        Some(State::SimFailed) | Some(State::DecodeFailed) => {
            return crate::output::colour::DEBUG_NAN
        }
        _ => [40, 40, 48],
    };
    let k = 0.55 + 0.15 * p.detail as f64;
    [
        (base[0] as f64 * k).min(255.0) as u8,
        (base[1] as f64 * k).min(255.0) as u8,
        (base[2] as f64 * k).min(255.0) as u8,
    ]
}

/// The number of event classes in the alphabet: three tightest-pair identities, then every
/// `(state, detail)` a terminated copy can carry, offset past them by
/// [`crate::ensemble::stats::TERMINAL_TAG`].
///
/// `pack()` is `(state << 2) | detail` with six states and four details, so the terminal arm is
/// 24 wide. **Fixed, not derived from the data**: a data-dependent alphabet would give the same
/// class a different colour in two slices, and two images that cannot be compared are the fault
/// this whole colouring exists to remove.
pub const N_EVENT_CLASSES: usize = 3 + 24;

/// The dense ordinal of an event class, or `None` if the byte is outside the alphabet.
///
/// The two alphabets are deliberately disjoint in [`crate::ensemble::stats::event_class_at`] --
/// pair indices `0..2`, terminal classes `TERMINAL_TAG..` -- and the gap between them is *not*
/// part of either. A byte landing in it is a bug, and returns `None` so it paints as
/// undetermined rather than as a plausible class.
pub fn event_class_ordinal(c: u8) -> Option<usize> {
    const TAG: u8 = crate::ensemble::stats::TERMINAL_TAG;
    if c < 3 {
        Some(c as usize)
    } else if c >= TAG && ((c - TAG) as usize) < N_EVENT_CLASSES - 3 {
        Some(3 + (c - TAG) as usize)
    } else {
        None
    }
}

/// A human-readable name for an event class, for the legend.
pub fn event_class_name(c: u8) -> String {
    const TAG: u8 = crate::ensemble::stats::TERMINAL_TAG;
    const PAIRS: [&str; 3] = ["tightest (0,1)", "tightest (0,2)", "tightest (1,2)"];
    if let Some(n) = PAIRS.get(c as usize) {
        return (*n).to_string();
    }
    match crate::outcome::Outcome::unpack(c.wrapping_sub(TAG)) {
        Some(o) if c >= TAG => format!("terminated: {} d{}", o.state.name(), o.detail),
        _ => format!("invalid ({c})"),
    }
}

/// **Colour by event class, on viridis** — the categorical mode the `Ma1achy/principia-ii`
/// WebGPU panel renders, and the one a reference comparison must be made under.
///
/// The class is the identity of the **currently tightest pair**, joined with the terminal
/// `(state, detail)` once a copy has terminated — [`crate::ensemble::stats::event_class_at`].
/// That is the quantity `spread_event` is built on, and it is defined at every playhead, unlike
/// the terminal outcome which is terminal-grain and reads pure under lockstep.
///
/// Undetermined takes [`crate::output::colour::DEBUG_NAN`] and **never a colourmap entry**: a
/// pixel with no value must be visibly null rather than plausibly a class.
///
/// Adjacent ordinals are close in colour by construction — 27 slots on one ramp — so **the
/// legend and the per-class histogram are the instrument, not the image**. Both are printed
/// beside every render.
pub fn event_class_rgb(p: &PixelOut) -> [u8; 3] {
    if p.n_nonfinite > 0 {
        return crate::output::colour::DEBUG_NAN;
    }
    match State::from_bits(p.state) {
        Some(State::SimFailed) | Some(State::DecodeFailed) | None => {
            return crate::output::colour::DEBUG_NAN
        }
        _ => {}
    }
    match event_class_ordinal(p.event_class) {
        Some(k) => {
            crate::output::viridis::viridis(k as f64 / (N_EVENT_CLASSES - 1) as f64)
        }
        None => crate::output::colour::DEBUG_NAN,
    }
}

/// Per-class pixel counts over the alphabet, in ordinal order, plus the undetermined count.
///
/// Printed before any event-class image. **A class that never fires is a fact about the slice**
/// and reads as a zero here; without the histogram an image with three colours in it and one
/// with twenty are indistinguishable at a glance.
pub fn event_class_histogram(px: &[PixelOut]) -> (Vec<(u8, usize)>, usize) {
    let mut counts = vec![0usize; N_EVENT_CLASSES];
    let mut undetermined = 0usize;
    for p in px {
        let bad = p.n_nonfinite > 0
            || matches!(
                State::from_bits(p.state),
                Some(State::SimFailed) | Some(State::DecodeFailed) | None
            );
        match (bad, event_class_ordinal(p.event_class)) {
            (false, Some(k)) => counts[k] += 1,
            _ => undetermined += 1,
        }
    }
    const TAG: u8 = crate::ensemble::stats::TERMINAL_TAG;
    let rows = counts
        .iter()
        .enumerate()
        .map(|(k, &n)| (if k < 3 { k as u8 } else { TAG + (k - 3) as u8 }, n))
        .collect();
    (rows, undetermined)
}

/// Perceptually monotone ramp for a value in `[0, 1]`. Not a scientific colourmap; it is
/// only here to make structure visible at a glance.
fn ramp(x: f64) -> [u8; 3] {
    let t = x.clamp(0.0, 1.0);
    let r = (255.0 * t.powf(0.6)) as u8;
    let g = (255.0 * (t * (1.0 - t) * 4.0).powf(0.8)) as u8;
    let b = (255.0 * (1.0 - t).powf(0.6)) as u8;
    [r, g, b]
}

fn save(path: &Path, w: u32, h: u32, data: &[u8]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), w, h);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    writer.write_image_data(data)?;
    Ok(())
}

/// Returns the `[lo, hi]` window the spread image was scaled over, so the caller can print it.
/// A false-colour image without its scale is decoration.
pub fn write_pair(stem: &str, slice: &Slice, pixels: &[PixelOut]) -> std::io::Result<(f64, f64)> {
    let (w, h) = (slice.nx as u32, slice.ny as u32);

    let mut a = Vec::with_capacity(pixels.len() * 3);
    for p in pixels {
        a.extend_from_slice(&outcome_rgb(p));
    }
    save(Path::new(&format!("{stem}_outcome.png")), w, h, &a)?;

    // `ensemble_spread` spans several decades and its median sits near 1e-3, so a linear
    // [0,1] ramp paints the whole grid at the bottom of the scale and the structure — thin
    // filaments where the copies decohere — disappears into flat background. Mapped on a log
    // scale between the grid's own p1 and p99 instead, which is a *diagnostic* choice: the
    // image is not the product, the raw dump is, and the image only has to make structure
    // visible. The window is printed alongside so the picture is not read as absolute.
    let mut fin: Vec<f64> = pixels
        .iter()
        .map(|p| p.ensemble_spread)
        .filter(|x| x.is_finite() && *x > 0.0)
        .collect();
    fin.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, hi) = if fin.len() < 2 {
        (1e-6, 1.0)
    } else {
        let q = |f: f64| fin[(((fin.len() - 1) as f64) * f).round() as usize];
        (q(0.01).max(1e-12), q(0.99).max(q(0.01) * 10.0))
    };
    let (ll, lh) = (lo.ln(), hi.ln());

    let mut b = Vec::with_capacity(pixels.len() * 3);
    for p in pixels {
        let v = p.ensemble_spread;
        // Non-finite is painted at full scale: undetermined is the loudest thing a pixel can
        // be, and must not be painted as quiet.
        let t = if !v.is_finite() {
            1.0
        } else if v <= 0.0 {
            0.0
        } else {
            ((v.ln() - ll) / (lh - ll)).clamp(0.0, 1.0)
        };
        b.extend_from_slice(&ramp(t));
    }
    save(Path::new(&format!("{stem}_spread.png")), w, h, &b)?;
    Ok((lo, hi))
}
