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

/// **The one ordering, and it is TOTAL.** Descending priority, `NaN` last, and **ties broken by
/// id ascending**.
///
/// The tie-break is the load-bearing part and it was missing. `top_k` flattens buckets top-down
/// and `rebuild` reads an id-sorted list, so with a merely *stable* sort two entries of equal
/// priority in different bands come out in opposite orders — measured: priorities `0.2 * 6/7` and
/// `0.4 * 3/7` are bitwise equal, and the two paths disagreed on which came first. That makes
/// [`Frontier::agrees_with_rebuild`] — the staleness check with teeth — report a **false**
/// disagreement on any tied data, which is the check firing on the wrong thing rather than
/// failing to fire. A total order removes the dependence on input order entirely, so a
/// disagreement between the paths can only be real.
fn by_priority(a: &(usize, f64), b: &(usize, f64)) -> std::cmp::Ordering {
    match (a.1.is_nan(), b.1.is_nan()) {
        (true, true) => a.0.cmp(&b.0),
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal).then(a.0.cmp(&b.0)),
    }
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
        scored.sort_by(by_priority);
        scored.into_iter().take(k).map(|(i, _)| i).collect()
    }

    /// Can any entry in band `b` or below still beat a `k`-th score of `kth`?
    ///
    /// **Expressed through [`band_of`] rather than by inverting it.** The obvious form computes
    /// band `b`'s upper bound analytically and compares — and the log round-trip does not land on
    /// the boundary: `exp(ln LO + 1/23 * d)` for band 0 returns `4.0616e-12`, which `band_of`
    /// places back in **band 0**. A bound that is too small stops the walk early with a contender
    /// unseen, which is unsound in the silent direction. Every entry of band `b` has
    /// `stored < hi_b` and `derive <= 1`, so its priority is below `hi_b` too and therefore lands
    /// in band `b` or lower; so `band_of(kth) > b` is exactly the condition, with no inverse.
    ///
    /// `NaN` maps to band 0 and so never stops the walk — conservative, and the right direction
    /// for a frontier that could not be scored.
    fn nothing_below_can_contend(kth: f64, b: usize) -> bool {
        band_of(kth) > b
    }

    /// The top `k`, walking bands from the top and stopping once no lower band can contend.
    ///
    /// **Exact**, on [`Self::top_k`]'s own argument taken one step further: the derived factor is
    /// in `[0, 1]`, so `stored * derive <= stored`, and every entry of band `b` has
    /// `stored < band_upper(b)`. Once `k` candidates are held whose `k`-th score is at least
    /// `band_upper(b)`, nothing in band `b` or below can enter, and the walk stops.
    ///
    /// **Returns the number of entries scored alongside the ids**, and that is the point. Whether
    /// bucketing buys anything over `top_k`'s full sort is an *empirical* question about band
    /// occupancy — if the signal piles into two or three bands the walk degenerates to the full
    /// scan and the frontier is a `HashMap` with extra steps. A version that could not report its
    /// own scan count could not settle that, and this project has shipped two mechanisms that
    /// computed, sorted and changed nothing.
    ///
    /// `NaN` scores never satisfy the stopping test, so a frontier of undetermined quads walks to
    /// the bottom — conservative, and the right direction.
    pub fn top_k_bounded<F: Fn(usize) -> f64>(&self, k: usize, derive: F) -> (Vec<usize>, usize) {
        if k == 0 {
            return (Vec::new(), 0);
        }
        let mut best: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        let mut scanned = 0usize;
        for b in (0..BANDS).rev() {
            if best.len() >= k && Self::nothing_below_can_contend(best[k - 1].1, b) {
                break;
            }
            if self.buckets[b].is_empty() {
                continue;
            }
            for e in &self.buckets[b] {
                scanned += 1;
                best.push((e.id, e.stored * derive(e.id)));
            }
            // Stable, so a tie keeps the entry from the higher band — which is the order
            // `top_k`'s single full sort produces, since it flattens bands top-down.
            best.sort_by(by_priority);
            best.truncate(k);
        }
        (best.into_iter().map(|(i, _)| i).collect(), scanned)
    }

    /// Band occupancy, for the measurement that decides whether the buckets earn their place.
    pub fn band_histogram(&self) -> Vec<usize> {
        self.buckets.iter().map(|b| b.len()).collect()
    }

    /// **The reference implementation. Kept permanently.**
    ///
    /// Build the ordering from scratch from a plain `(id, stored)` list, with no buckets and no
    /// incremental state. Slower by construction and correct by construction, which is the
    /// point: [`Self::agrees_with_rebuild`] compares the two.
    pub fn rebuild<F: Fn(usize) -> f64>(items: &[(usize, f64)], k: usize, derive: F) -> Vec<usize> {
        let mut v: Vec<(usize, f64)> = items.iter().map(|&(i, s)| (i, s * derive(i))).collect();
        v.sort_by(by_priority);
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
