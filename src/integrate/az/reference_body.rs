//! Choosing the reference body, and the policy for whether ensemble copies share it.

use crate::physics::{newton, THIRD};
use crate::{Real, Vec2};

/// Whether the ensemble copies of a pixel share the nominal copy's reference body.
///
/// **The flag governs cross-copy sharing only, never freezing across time.** Sharing (across
/// copies) and switching (across time) are separate knobs; freezing the reference across
/// time would break AZ outright, since the whole point of re-choosing is to keep the
/// reference out of the longest side as the triangle deforms.
///
/// Default is `PerCopy`. Forcing a copy onto another copy's reference can put its reference
/// body *inside* its longest side, at which point `|R3| >= max(|R1|,|R2|)` no longer holds
/// and an unregularised pair can close — the AZ guarantee is not degraded but void. See
/// `NOTES.md` §1.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RefPolicy {
    #[default]
    PerCopy,
    Shared,
}

/// The body **not** in the longest side, so the unregularised side `(b,c)` is the longest.
///
/// Tie-breaking matches `numpy.argmax`, which returns the *first* maximal index. A near-tie
/// broken the other way would fail the cross-check while looking exactly like a
/// transcription error, so the reference body is logged per sync and compared as a column.
pub fn choose_reference<T: Real>(r: &[Vec2<T>; 3]) -> usize {
    let d = newton::pair_dists(r);
    let mut longest = 0usize;
    // Strict `>` only, so the first maximum wins — numpy's convention.
    for k in 1..3 {
        if d[k] > d[longest] {
            longest = k;
        }
    }
    THIRD[longest]
}

/// `(a, b, c)` for a given reference body, matching `tb_az.TRIPLES`.
pub fn triple(a: usize) -> (usize, usize, usize) {
    match a {
        0 => (0, 1, 2),
        1 => (1, 0, 2),
        _ => (2, 0, 1),
    }
}
