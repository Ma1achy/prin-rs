//! **Payload residency: which quads still hold their samples, and what it cost to lose them.**
//!
//! The caching contract's Part 7 is the shape of this: what a session keeps is **current state
//! only — one timestep per quad, never a trajectory, never a time-series.** A revisited quad is a
//! *resume point*, marched `t_cached -> playhead` rather than re-booted from zero, and fixed-`dt`
//! determinism makes the resumed state *equal* the never-evicted one. So eviction is safe by
//! construction, and the only question it raises is cost.
//!
//! **Three states, never two.** "Never computed" and "computed and released" look identical to a
//! renderer and must not look identical to telemetry: a hole from an eviction is a budget fact and
//! a hole from an unreached quad is a coverage fact. Pooling them is the stop-reason conflation one
//! level down, which this project has already paid for.
//!
//! **The tree never learns that eviction exists.** `Quad` keeps its reduction, its exponents and
//! its decision; only the `Vec<PixelOut>` goes. That is the caching contract's *"refinement latches
//! are not cache entries"* — the verdict that a quad proved interesting is CPU metadata that
//! survives its state being dropped. It is also how *"camera state is never a `Quad` fact"*
//! survives: `Residency` and the LRU clock are camera-derived and live here, not on the node.
//!
//! **Eviction degrades exactly into the coarse-ancestor fill, which is already built.**
//! `adaptive::render_leaves` skips a leaf whose samples are empty and its nearest painted ancestor
//! shows through — §4.5's option 1, *the only one that never lies*. No renderer change is needed
//! for eviction to be correct; the only change is being able to measure it.

use crate::ensemble::pixel::PixelOut;

/// What a quad's `N x N` payload is doing.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Residency {
    /// Never computed. A hole here is a **coverage** fact.
    #[default]
    Absent,
    /// Computed and held.
    Resident,
    /// Computed, reduced, released. `Quad::red` is still a real measurement, and `substeps` says
    /// what getting it back would cost **without** getting it back. A hole here is a **budget**
    /// fact.
    Evicted { frame: u64, substeps: u64 },
}

impl Residency {
    pub fn name(self) -> &'static str {
        match self {
            Residency::Absent => "absent",
            Residency::Resident => "resident",
            Residency::Evicted { .. } => "evicted",
        }
    }
}

/// How much a store keeps.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Retain {
    /// Records nothing and allocates nothing — today's `keep_pixels = false`, and the default,
    /// so the ~100 MB the batch descent does not want stays unallocated.
    #[default]
    Never,
    /// Keeps everything. The pre-session behaviour under `keep_pixels = true`.
    All,
    /// Keeps at most `cap` payloads resident, releasing the least relevant.
    Capped(usize),
}

/// Payloads by node index, with residency and an LRU clock beside them.
#[derive(Clone, Debug, Default)]
pub struct PixelStore {
    slots: Vec<Vec<PixelOut>>,
    state: Vec<Residency>,
    /// Frame each slot was last written or read. **Not on the `Quad`.**
    touched: Vec<u64>,
    retain: Retain,
    resident: usize,
}

impl PixelStore {
    pub fn new(retain: Retain) -> PixelStore {
        PixelStore { retain, ..Default::default() }
    }

    pub fn retain(&self) -> Retain {
        self.retain
    }

    pub fn resident(&self) -> usize {
        self.resident
    }

    fn grow(&mut self, i: usize) {
        if self.slots.len() <= i {
            self.slots.resize(i + 1, Vec::new());
            self.state.resize(i + 1, Residency::Absent);
            self.touched.resize(i + 1, 0);
        }
    }

    /// Record a freshly computed payload. Under [`Retain::Never`] this is a no-op that allocates
    /// nothing, which is what keeps the one-shot descent's memory exactly where it was.
    pub fn put(&mut self, i: usize, px: Vec<PixelOut>, frame: u64) {
        if self.retain == Retain::Never {
            return;
        }
        self.grow(i);
        if self.state[i] != Residency::Resident {
            self.resident += 1;
        }
        self.slots[i] = px;
        self.state[i] = Residency::Resident;
        self.touched[i] = frame;
    }

    pub fn residency(&self, i: usize) -> Residency {
        self.state.get(i).copied().unwrap_or_default()
    }

    /// The payload, or an empty slice. **An evicted quad and an absent one both read empty here**
    /// — that is deliberate, because the renderer must treat them alike (both fall back to the
    /// ancestor). Telemetry asks [`Self::residency`], which distinguishes them.
    pub fn get(&self, i: usize) -> &[PixelOut] {
        self.slots.get(i).map_or(&[], |v| v.as_slice())
    }

    /// Mark a slot used this frame, for the LRU tie-break.
    pub fn touch(&mut self, i: usize, frame: u64) {
        if i < self.touched.len() {
            self.touched[i] = frame;
        }
    }

    pub fn last_touched(&self, i: usize) -> u64 {
        self.touched.get(i).copied().unwrap_or(0)
    }

    /// Release payloads until at most `cap` remain resident, lowest score first.
    ///
    /// **Deterministic**: the comparison is total — score, then LRU, then node index — so a frame's
    /// evictions are a function of `(tree, camera, frame)` and nothing else. A non-total order here
    /// would make the session's memory behaviour depend on `HashMap` iteration order, which is the
    /// scheduler-firewall violation one level down.
    ///
    /// `keep` is the exempt class: the caching contract pins the baseline cover, the visible
    /// ancestor chain and the backdrop's leaf cover, because *breaking the visible fallback chain
    /// is the one eviction that visibly hurts — the difference between "blurry for a moment" and
    /// "blank"*.
    ///
    /// Returns the indices released.
    pub fn evict_to(
        &mut self,
        cap: usize,
        frame: u64,
        keep: &dyn Fn(usize) -> bool,
        score: &dyn Fn(usize) -> f64,
        substeps: &dyn Fn(usize) -> u64,
    ) -> Vec<usize> {
        if self.resident <= cap {
            return Vec::new();
        }
        let mut cands: Vec<usize> = (0..self.state.len())
            .filter(|&i| self.state[i] == Residency::Resident && !keep(i))
            .collect();
        cands.sort_by(|&a, &b| {
            let (sa, sb) = (score(a), score(b));
            sa.partial_cmp(&sb)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(self.touched[a].cmp(&self.touched[b]))
                .then(a.cmp(&b))
        });
        let want = self.resident.saturating_sub(cap);
        let mut freed = Vec::new();
        for &i in cands.iter().take(want) {
            self.slots[i] = Vec::new();
            self.state[i] = Residency::Evicted { frame, substeps: substeps(i) };
            self.resident -= 1;
            freed.push(i);
        }
        freed
    }

    /// The dense `Vec<Vec<PixelOut>>` the one-shot descent has always returned.
    ///
    /// Must reproduce `st.pixels.resize(i + 1, Vec::new())` exactly — same length, same empty
    /// slots — or `SchedStats::pixels` changes shape and every consumer that indexes it moves.
    pub fn into_dense(self) -> Vec<Vec<PixelOut>> {
        self.slots
    }
}
