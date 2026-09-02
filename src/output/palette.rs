//! **Mode 1 — the outcome/event-class categorical map.** The specced default, replacing ad-hoc.
//!
//! # The assignment is a mnemonic and that is the point
//!
//! ```text
//!   COLLISIONS = ADDITIVE primaries, keyed by the colliding PAIR
//!     collision 1-2   red      #DE2D2D
//!     collision 1-3   green    #2EBC4E
//!     collision 2-3   blue     #3462E0
//!
//!   ESCAPES    = SUBTRACTIVE primaries, keyed by the escaping BODY
//!     body 1 escape   yellow   #F0DE32
//!     body 2 escape   magenta  #E034C6
//!     body 3 escape   cyan     #30C8DC
//!
//!   NON-GENERIC
//!     bounded          black   #141418
//!     collision @ t=0  orange  #F29620
//!     degenerate       white   #ECECF0
//! ```
//!
//! The two event **families** separate at a glance while pair/body identity stays legible. A
//! generated palette -- a golden-angle cycle, viridis over 27 ordinals -- cannot encode that,
//! which is why this one is not substitutable. `png::event_class_rgb`'s 27 fixed viridis slots
//! are the ad-hoc map this supersedes.
//!
//! # It reads `state` plus the `detail` union, and there is no separate escaper field
//!
//! BRIEF §2.4 packs `state` in 3 bits and `detail` in 2. Collision -> pair index, escape -> body
//! index, so R/G/B and Y/M/C both land on `detail`. **The escaping body IS `detail | state =
//! escape`** -- looking for an `escaper` field is looking for something that does not exist.
//! `PAIRS = [(0,1),(0,2),(1,2)]`, so `detail` 0/1/2 is 1-2 / 1-3 / 2-3 in one-based naming.
//!
//! # Brightness carries the event time, and the polarity is deliberately inverted
//!
//! Hue is identity, lightness is `t_end` through the standard compaction, with **white = LOW /
//! EARLY**: quick-resolving pixels pop and late ones darken, which keeps `bounded` (never
//! resolving) consistent with its black swatch. That is the **opposite** of the FTLE/diffusion
//! convention, where white = high. The inconsistency is deliberate and specced, and the legend
//! says so rather than leaving a reader to infer it.
//!
//! `Scalar::TEnd` already carries `Direction::HighIsSettled`, so the existing compaction gives
//! this polarity with no special case.
//!
//! # Two things this codebase has that the nine classes do not name
//!
//! **`detail == 3` is "all three"** for both arms -- triple collision and triple ejection, the
//! `>=2-pair` rule. It is a determinate physical outcome, not an invalid pixel, and it is not one
//! of the nine. It gets its own editable slot defaulting to the degenerate white, and is called
//! out here rather than folded silently into a neighbour.
//!
//! **`Running`, `SimFailed`, `DecodeFailed`** are not outcomes. They take the explicit invalid
//! colour. *Do not let an invalid pixel inherit a valid colour by accident -- that is how the
//! wedges hid.*
use crate::ensemble::pixel::PixelOut;
use crate::outcome::State;

/// One slot of the categorical map. The ordering is the legend order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum EventClass {
    Collision12,
    Collision13,
    Collision23,
    CollisionTriple,
    Escape1,
    Escape2,
    Escape3,
    EscapeTriple,
    Bounded,
    CollisionAtZero,
    Degenerate,
    /// Not an outcome: `Running`, `SimFailed`, `DecodeFailed`, or an unrecognised bit pattern.
    Invalid,
}

impl EventClass {
    pub fn name(self) -> &'static str {
        match self {
            EventClass::Collision12 => "collision 1-2",
            EventClass::Collision13 => "collision 1-3",
            EventClass::Collision23 => "collision 2-3",
            EventClass::CollisionTriple => "collision triple",
            EventClass::Escape1 => "body 1 escape",
            EventClass::Escape2 => "body 2 escape",
            EventClass::Escape3 => "body 3 escape",
            EventClass::EscapeTriple => "triple ejection",
            EventClass::Bounded => "bounded",
            EventClass::CollisionAtZero => "collision at t=0",
            EventClass::Degenerate => "degenerate",
            EventClass::Invalid => "invalid",
        }
    }

    pub const ALL: [EventClass; 12] = [
        EventClass::Collision12,
        EventClass::Collision13,
        EventClass::Collision23,
        EventClass::CollisionTriple,
        EventClass::Escape1,
        EventClass::Escape2,
        EventClass::Escape3,
        EventClass::EscapeTriple,
        EventClass::Bounded,
        EventClass::CollisionAtZero,
        EventClass::Degenerate,
        EventClass::Invalid,
    ];
}

/// `state` + `detail` -> class. **`detail` carries the pair for a collision and the body for an
/// escape**; there is no second field.
///
/// `t_end <= 0` under `Collision` is the at-`t=0` case, which is a property of the initial
/// conditions rather than of the march and is why it has its own swatch.
pub fn classify(p: &PixelOut) -> EventClass {
    let st = match State::from_bits(p.state) {
        Some(s) => s,
        None => return EventClass::Invalid,
    };
    let d = p.detail & 3;
    match st {
        State::Collision => {
            if p.t_end.is_finite() && p.t_end <= 0.0 {
                return EventClass::CollisionAtZero;
            }
            match d {
                0 => EventClass::Collision12,
                1 => EventClass::Collision13,
                2 => EventClass::Collision23,
                _ => EventClass::CollisionTriple,
            }
        }
        State::Escape => match d {
            0 => EventClass::Escape1,
            1 => EventClass::Escape2,
            2 => EventClass::Escape3,
            _ => EventClass::EscapeTriple,
        },
        State::Bounded => EventClass::Bounded,
        State::Running | State::SimFailed | State::DecodeFailed => EventClass::Invalid,
    }
}

/// The class -> swatch table. **User-editable: this is the default, not a fixed mapping.**
#[derive(Clone, Debug, PartialEq)]
pub struct Swatches {
    /// Parallel to [`EventClass::ALL`].
    pub rgb: [[u8; 3]; 12],
}

impl Default for Swatches {
    fn default() -> Self {
        Self {
            rgb: [
                [0xDE, 0x2D, 0x2D], // collision 1-2   red
                [0x2E, 0xBC, 0x4E], // collision 1-3   green
                [0x34, 0x62, 0xE0], // collision 2-3   blue
                [0xEC, 0xEC, 0xF0], // collision triple -- NOT one of the nine; see module docs
                [0xF0, 0xDE, 0x32], // body 1 escape   yellow
                [0xE0, 0x34, 0xC6], // body 2 escape   magenta
                [0x30, 0xC8, 0xDC], // body 3 escape   cyan
                [0xEC, 0xEC, 0xF0], // triple ejection -- NOT one of the nine; see module docs
                [0x14, 0x14, 0x18], // bounded         black
                [0xF2, 0x96, 0x20], // collision @ t=0 orange
                [0xEC, 0xEC, 0xF0], // degenerate      white
                [0xFF, 0x00, 0xFF], // invalid -- the explicit invalid colour, never inherited
            ],
        }
    }
}

impl Swatches {
    pub fn get(&self, c: EventClass) -> [u8; 3] {
        self.rgb[EventClass::ALL.iter().position(|&x| x == c).unwrap()]
    }

    /// Replace one class's swatch. The table is a default, not a fixture.
    pub fn set(&mut self, c: EventClass, rgb: [u8; 3]) {
        let i = EventClass::ALL.iter().position(|&x| x == c).unwrap();
        self.rgb[i] = rgb;
    }
}

/// A categorical filter: show these classes, mute the rest.
///
/// **A general operation on any categorical mode, not a separate render mode.** *"Just
/// collisions"*, *"just body-2 escape"* are this with a different `show` set; a standalone
/// which-body-escaped view is this map filtered, and building it as its own mode would duplicate
/// the palette and let the two drift.
#[derive(Clone, Debug, PartialEq)]
pub struct Filter {
    pub show: Vec<EventClass>,
    /// What a muted class renders as. Deliberately not the invalid colour: muted means *not
    /// selected*, and undetermined means *not known*, and a reader must be able to tell them
    /// apart.
    pub muted: [u8; 3],
}

impl Default for Filter {
    fn default() -> Self {
        Self { show: EventClass::ALL.to_vec(), muted: [0x2A, 0x2A, 0x2E] }
    }
}

impl Filter {
    pub fn only(classes: &[EventClass]) -> Self {
        Self { show: classes.to_vec(), muted: Self::default().muted }
    }
    pub fn collisions() -> Self {
        Self::only(&[
            EventClass::Collision12,
            EventClass::Collision13,
            EventClass::Collision23,
            EventClass::CollisionTriple,
        ])
    }
    pub fn escapes() -> Self {
        Self::only(&[
            EventClass::Escape1,
            EventClass::Escape2,
            EventClass::Escape3,
            EventClass::EscapeTriple,
        ])
    }
    pub fn passes(&self, c: EventClass) -> bool {
        self.show.contains(&c)
    }
}

// -------------------------------------------------------------------------------------------
// The render
// -------------------------------------------------------------------------------------------

use crate::output::colour::{self, Scalar};
use crate::output::oklab;

/// Mode 1's colour for one footprint: **hue = class identity, lightness = event time.**
///
/// Replace-L, the same combiner the shape mode uses: the swatch's `(a, b)` are kept and its own
/// `L` is discarded in favour of the compacted `t_end`. So the class stays identifiable at every
/// brightness and the two channels do not contaminate each other.
///
/// **`Invalid` is returned unmodulated**, and so is a muted class. An undetermined pixel must not
/// be lightened into looking like an early-resolving one, and a *muted* pixel means "not
/// selected" rather than "not known" -- three states, three appearances.
///
/// **Bounded falls out black without a special case.** `Scalar::TEnd` carries
/// `Direction::HighIsSettled`, so a footprint that never resolved sits at `t_end = t_max`,
/// normalises to 0, and takes `L_MIN`. That is why the inverted polarity was chosen rather than
/// worked around.
pub fn event_rgb(p: &PixelOut, sw: &Swatches, f: &Filter, lo: f64, hi: f64) -> [u8; 3] {
    let c = classify(p);
    if c == EventClass::Invalid {
        return sw.get(EventClass::Invalid);
    }
    if !f.passes(c) {
        return f.muted;
    }
    let base = sw.get(c);
    match colour::range_norm(Scalar::TEnd, p.t_end, lo, hi) {
        Some(t) => {
            let lab = oklab::srgb_to_oklab(base);
            oklab::oklab_to_srgb([colour::lightness(t), lab[1], lab[2]])
        }
        // A determinate class with an unreadable time: keep the class, do not invent a time.
        None => base,
    }
}

/// [`event_rgb`] **resolved over the ensemble** -- the categorical arm of the same supersampling
/// `colour::rgb_resolved` does for the shape sphere, and `ssaa::resolve_rgb` already did for the
/// older outcome palette.
///
/// Every sub-sample is a full simulation with its own class; the pixel is the mean of their
/// colours. A footprint split 4/4 between two outcomes resolves to a blend, which is what it
/// *looks like*, while `ensemble_spread` separately says "refine here". Substituting either for
/// the other loses a real distinction.
///
/// Requires `keep_copy_outcomes`; without it `copy_outcomes` is empty and this returns
/// [`event_rgb`] unchanged rather than a silently different picture.
pub fn event_rgb_resolved(
    p: &PixelOut,
    sw: &Swatches,
    f: &Filter,
    lo: f64,
    hi: f64,
) -> [u8; 3] {
    use crate::output::compose;
    if p.copy_outcomes.len() < 2 {
        return event_rgb(p, sw, f, lo, hi);
    }
    if State::from_bits(p.state).is_none() {
        return sw.get(EventClass::Invalid);
    }
    // Lightness is the footprint's own `t_end`: per-copy event times are not retained, so the
    // colour channel resolves and the lightness does not. Stated, not implied.
    let l = colour::range_norm(Scalar::TEnd, p.t_end, lo, hi)
        .map(colour::lightness)
        .unwrap_or_else(|| oklab::srgb_to_oklab(sw.get(classify(p)))[0]);
    // The same supersampler the continuous maps use, given a categorical sub-sample colour.
    let lab = compose::resolve(p.copy_outcomes.len(), |i| {
        let packed = p.copy_outcomes[i];
        let st = State::from_bits(packed >> 2)?;
        let probe = PixelOut { state: st as u8, detail: packed & 3, ..p.clone() };
        let c = classify(&probe);
        let rgb = if f.passes(c) { sw.get(c) } else { f.muted };
        let lab = oklab::srgb_to_oklab(rgb);
        Some([l, lab[1], lab[2]])
    });
    compose::finish(lab, sw.get(EventClass::Invalid))
}

/// The legend, including the polarity note. Printed beside every Mode 1 render.
pub fn legend(sw: &Swatches, f: &Filter) -> String {
    let mut out = String::from(
        "MODE 1 -- outcome/event class. Hue = identity, lightness = t_end.\n\
         **WHITE = LOW / EARLY**: quick-resolving pixels pop, late ones darken, and `bounded`\n\
         (never resolving) sits at black consistent with its swatch. This is the OPPOSITE of the\n\
         FTLE/diffusion convention, where white = high. The inconsistency is deliberate.\n\n",
    );
    for c in EventClass::ALL {
        let [r, g, b] = sw.get(c);
        let mark = if f.passes(c) { ' ' } else { '-' };
        out.push_str(&format!("  {mark} #{r:02X}{g:02X}{b:02X}  {}\n", c.name()));
    }
    out.push_str("\n  '-' is muted by the active filter: NOT SELECTED, which is a different\n");
    out.push_str("  thing from `invalid` (not known) and renders differently.\n");
    out
}
