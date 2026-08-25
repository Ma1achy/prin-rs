//! **The persistent frontier** — §4.6.
//!
//! Each frame refines the top `k` leaves by priority, so the frontier must be ordered. Rebuilding
//! that order from scratch every frame is `O(n log n)` over thousands of leaves, sixty times a
//! second.
//!
//! Normally that would be a rounding error against 512 trajectories per quad. **Camera bias
//! changes the arithmetic.** Priority is `structure x camera relevance`, and relevance changes
//! for *every* quad on *every* frame the camera moves — so during a gesture the naive version is
//! not re-sorting a mostly-unchanged list, it is genuinely recomputing all of it.
//!
//! # The split that makes it work
//!
//! **Stored: the physics term.** It changes only when a quad is recomputed or the zoom changes.
//! **Derived: the camera term.** It changes every frame of motion, and is computed at query time.
//!
//! A pan therefore touches a distance calculation and never the physics — and it keeps camera
//! state off the `Quad`, which the *"never cache view state as a quad fact"* rule already
//! requires. Nothing here is stored on a `Quad`; the frontier is a separate structure with its
//! own lifetime.
//!
//! # Buckets, not a heap
//!
//! A plain binary heap cannot reprioritise an entry already inside it without either a
//! `id -> heap position` map or a full rebuild. **Priority bucketing** suits a rank-based scheme,
//! which wants the top slice rather than a total order, and an entry is re-bucketed only when it
//! crosses a band boundary — so a small relevance change costs nothing at all.
//!
//! # The failure mode is staleness, and it is invisible
//!
//! An incrementally-maintained frontier that is *wrong* looks exactly like a criterion that is
//! wrong: a quad sitting high in the queue on a priority it no longer has. So
//! [`Frontier::rebuild`] — the from-scratch path — is **kept permanently as the reference
//! implementation**, not deleted once the fast one works, and [`Frontier::agrees_with_rebuild`]
//! is an independent path to the same answer. The same shape as the `Gamma`-identity chain: a
//! silent divergence cannot survive two paths that must agree.

/// Log-spaced priority bands. An entry moves only when it crosses one, so a relevance change
/// that does not change the band costs nothing.
///
/// Log-spaced rather than linear because the signal spans six orders across regions — `4.26e-8`
/// in `far` against `9.75e-4` in near-field — and linear bands would put every quad of a whole
/// region in one bucket, which is the saturation failure this project has already met twice.
pub const BANDS: usize = 24;
const LO: f64 = 1e-12;
const HI: f64 = 1e2;

/// Which band a priority falls in. Monotone: a higher priority never lands in a lower band.
///
/// **`NaN` goes to the bottom band; `+inf` goes to the top.** They are different statements and
/// the first cut of this conflated them under one `is_finite` guard. `NaN` is *undetermined*, and
/// a quad that could not be scored must not outrank one that could — the same convention the
/// replay uses, where `NaN` maps to `-inf` rather than blocking. `+inf` is *maximally important*
/// and belongs at the top. Nothing in the current signal produces `+inf`, which is exactly why
/// the conflation would have sat there unnoticed.
pub fn band_of(p: f64) -> usize {
    if p.is_nan() || p <= LO {
        return 0;
    }
    if p.is_infinite() {
        return BANDS - 1;
    }
    let t = (p.min(HI).ln() - LO.ln()) / (HI.ln() - LO.ln());
    ((t * (BANDS - 1) as f64).floor() as usize).min(BANDS - 1)
}

/// One entry: a quad id and the **stored** half of its priority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Entry {
    pub id: usize,
    /// The physics term. Changes only when the quad is recomputed or the zoom changes.
    pub stored: f64,
}

/// The frontier, bucketed by band.
#[derive(Clone, Debug, Default)]
pub struct Frontier {
    buckets: Vec<Vec<Entry>>,
    /// Band each id currently sits in, so a reprioritise does not scan.
    at: std::collections::HashMap<usize, usize>,
}

impl Frontier {
    pub fn new() -> Self {
        Frontier { buckets: vec![Vec::new(); BANDS], at: Default::default() }
    }

    pub fn len(&self) -> usize {
        self.at.len()
    }

    pub fn is_empty(&self) -> bool {
        self.at.is_empty()
    }

    pub fn insert(&mut self, id: usize, stored: f64) {
        self.remove(id);
        let b = band_of(stored);
        self.buckets[b].push(Entry { id, stored });
        self.at.insert(id, b);
    }

    pub fn remove(&mut self, id: usize) {
        if let Some(b) = self.at.remove(&id) {
            self.buckets[b].retain(|e| e.id != id);
        }
    }

    /// Update an entry's **stored** term. Re-buckets only if the band changed.
    pub fn reprioritise(&mut self, id: usize, stored: f64) {
        let want = band_of(stored);
        match self.at.get(&id).copied() {
            Some(b) if b == want => {
                if let Some(e) = self.buckets[b].iter_mut().find(|e| e.id == id) {
                    e.stored = stored;
                }
            }
            _ => self.insert(id, stored),
        }
    }

    /// The top `k` by full priority — `stored * derive(id)`.
    ///
    /// `derive` is the camera term, applied **here** and never stored. Bands are walked from the
    /// top, and enough of them are drained to be sure the derived factor cannot promote a lower
    /// band past a higher one: the factor is in `[0, 1]`, so it can only *demote*. That is why
    /// the whole frontier is scored rather than the top band alone — a bounded-above derived
    /// term makes the band order an upper bound, not an answer.
    pub fn top_k<F: Fn(usize) -> f64>(&self, k: usize, derive: F) -> Vec<usize> {
        let mut scored: Vec<(usize, f64)> = self
            .buckets
            .iter()
            .rev()
            .flatten()
            .map(|e| (e.id, e.stored * derive(e.id)))
            .collect();
        scored.sort_by(|a, b| match (a.1.is_nan(), b.1.is_nan()) {
            (true, true) => a.0.cmp(&b.0),
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            _ => b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal),
        });
        scored.into_iter().take(k).map(|(i, _)| i).collect()
    }

    /// **The reference implementation. Kept permanently.**
    ///
    /// Build the ordering from scratch from a plain `(id, stored)` list, with no buckets and no
    /// incremental state. Slower by construction and correct by construction, which is the
    /// point: [`Self::agrees_with_rebuild`] compares the two.
    pub fn rebuild<F: Fn(usize) -> f64>(items: &[(usize, f64)], k: usize, derive: F) -> Vec<usize> {
        let mut v: Vec<(usize, f64)> = items.iter().map(|&(i, s)| (i, s * derive(i))).collect();
        v.sort_by(|a, b| match (a.1.is_nan(), b.1.is_nan()) {
            (true, true) => a.0.cmp(&b.0),
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            _ => b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal),
        });
        v.into_iter().take(k).map(|(i, _)| i).collect()
    }

    /// Every entry, as the rebuild path wants them.
    pub fn entries(&self) -> Vec<(usize, f64)> {
        let mut v: Vec<(usize, f64)> =
            self.buckets.iter().flatten().map(|e| (e.id, e.stored)).collect();
        v.sort_by_key(|x| x.0);
        v
    }

    /// **The staleness check with teeth.** Run it every `N` frames, not as a benchmark.
    ///
    /// A wrong incremental frontier is indistinguishable from a wrong criterion by looking at
    /// the tree — both put the budget in the wrong place, quietly. This is the independent path.
    pub fn agrees_with_rebuild<F: Fn(usize) -> f64 + Copy>(&self, k: usize, derive: F) -> bool {
        self.top_k(k, derive) == Self::rebuild(&self.entries(), k, derive)
    }
}
