//! Outcome classification.
//!
//! Step 5a needs only the **legacy** classifier: `spread_event` is defined over outcome
//! classes, and `tb.classify` is the one class labelling with a reference to check against.
//!
//! BRIEF §2.4's 3-bit `state` / 2-bit `detail` encoding, `r_coll`, and the >=2-pair triple
//! rule are Step 5b. They have no reference at all — `r_coll` appears nowhere in the numpy
//! tree — so they are separated deliberately: the ported part and the invented part should
//! be reviewable apart from each other.

use crate::physics::{newton, Cart, G, PAIRS, THIRD};
use crate::Real;

/// `tb.classify`, transcribed.
///
/// Returns the index of the escaping body (`0`, `1`, `2`) or `3` for still-bound. The
/// candidate escaper is the body not in the *tightest* pair; it is labelled escaping only if
/// its specific orbital energy relative to the other two is positive **and** it is receding.
/// Bodies not selected as the candidate can never be labelled escaping — a property of the
/// reference, transcribed rather than fixed.
pub fn classify_legacy<T: Real>(s: &Cart<T>, m: &[T; 3]) -> u8 {
    let d = newton::pair_dists(&s.r);
    let mut tight = 0usize;
    for k in 1..3 {
        if d[k] < d[tight] {
            tight = k;
        }
    }
    let b = THIRD[tight];
    let (o1, o2) = {
        let others: Vec<usize> = (0..3).filter(|&k| k != b).collect();
        (others[0], others[1])
    };

    let m_bin = m[o1] + m[o2];
    let rc = (s.r[o1] * m[o1] + s.r[o2] * m[o2]) / m_bin;
    let vc = (s.v[o1] * m[o1] + s.v[o2] * m[o2]) / m_bin;
    let dr = s.r[b] - rc;
    let dv = s.v[b] - vc;
    let dist = dr.norm();

    let half = T::lit(0.5);
    let spec = half * dv.norm_sq() - T::lit(G) * m_bin / dist.max(T::CLASSIFY_FLOOR);
    let receding = dr.dot(dv) > T::zero();

    if spec > T::zero() && receding {
        b as u8
    } else {
        3
    }
}

/// Which pair is tightest at the final state, as an index into [`PAIRS`]. Cheap, and useful
/// for interpreting `spread_event` — it says *which* binary formed, not merely that the
/// copies disagreed.
pub fn binary_id<T: Real>(s: &Cart<T>) -> u8 {
    let d = newton::pair_dists(&s.r);
    let mut best = 0usize;
    for k in 1..3 {
        if d[k] < d[best] {
            best = k;
        }
    }
    let _ = PAIRS;
    best as u8
}
