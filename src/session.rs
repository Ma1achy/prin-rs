//! **A quadtree that survives the camera: the frame loop, persistence, and bounded eviction.**
//!
//! Every scheduler experiment before this spends a fixed quad count **once**. A slippy map has
//! ~16.7 ms per frame, **forever** — which turns *"which quads"* into *"which quads this frame"*,
//! and is what makes breadth-first, the persistent frontier and the coarse-ancestor fill
//! requirements rather than preferences (§4.2).
//!
//! # What this is not
//!
//! Not a second scheduler. `Session::step` calls [`scheduler::round`] — the same eight steps the
//! one-shot descent runs, extracted rather than rewritten, and pinned by
//! `tests/round_extraction_golden.rs` against a hash taken before the extraction existed.
//!
//! Not a wall-clock deadline. The quota is a **count**, because a deadline makes the tree a
//! function of machine load, and every result in `results/` is reproducible precisely because
//! nothing in the descent reads a clock. The harness measures what the quota *cost* in
//! milliseconds and reports whether it fits 16.7 and 41.7; the tree does not depend on the answer.
//!
//! # The three things a session adds
//!
//! **A frame quota.** Bounds what one frame computes. Under the scheduler contract's lockstep loop
//! the playhead advances one fixed `dt` per frame regardless, and the quota governs the
//! *catching-up* tier — so a frame that spends none of its quota still advances time.
//!
//! **Persistence across camera moves.** The caching contract is explicit: pan, zoom, tilt, slice
//! and chart switch change *which identities are requested*, never the validity of anything
//! computed. **Navigation re-addresses; it does not invalidate.** What a zoom *does* invalidate is
//! the **priority**, because structure is pixel-relative (§1.2) — no physics moves, and every
//! stored ranking term does.
//!
//! **Bounded eviction**, in `crate::store`. Payloads go; `Quad` nodes, reductions and decisions
//! stay. That is the caching contract's *"refinement latches are not cache entries"*, and it is
//! what makes eviction unable to change a decision under a frozen playhead: `decide` does not take
//! pixels and nothing it calls does, so after `reduce` there is no path from a payload to a
//! decision at all.
//!
//! # The guard
//!
//! `payload_cap: None` with `regrow: Off` is a `Session` in which nothing a session adds can
//! engage — a `descend_with` with logging, wearing the session's name. That is `k_frac = 1.0`
//! again, at a new site, and [`assert_session_engages`] refuses it rather than trusting a
//! convention.

use crate::camera::Camera;
use crate::frontier::Frontier;
use crate::scheduler::{self, DescentState, SchedCfg, Spend, Stop};
use crate::store::{PixelStore, Retain};

/// Deterministic per-frame quota. **Counts, never a deadline.**
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameQuota {
    /// Quads this frame may compute. The primary knob.
    pub quads: usize,
    /// Substeps this frame may integrate, or `None`.
    ///
    /// **Enforced post hoc and it says so**: a quad's substep cost is known only *after* it is
    /// computed, so this stops the frame at the first round that crosses the line and can overshoot
    /// by one round's work. A soft cap that reported itself as hard would be the more dangerous
    /// arrangement.
    pub substeps: Option<u64>,
    /// Rounds this frame may run. Needed because a small `quads` still pays the decide/order pass
    /// once per round, so without it a tight frame is all overhead.
    pub rounds: u32,
}

impl Default for FrameQuota {
    fn default() -> Self {
        FrameQuota { quads: 64, substeps: None, rounds: 8 }
    }
}

/// Why a frame stopped. Printed, never inferred — the stop-reason discipline one level up.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum QuotaHit {
    /// The frontier emptied: everything the criterion wanted is done.
    #[default]
    Drained,
    Quads,
    Substeps,
    Rounds,
    /// `cfg.budget`, the run-wide cap. Distinct from the three above: this is terminal.
    TotalBudget,
}

impl QuotaHit {
    pub fn name(self) -> &'static str {
        match self {
            QuotaHit::Drained => "drained",
            QuotaHit::Quads => "quads",
            QuotaHit::Substeps => "substeps",
            QuotaHit::Rounds => "rounds",
            QuotaHit::TotalBudget => "total_budget",
        }
    }
}

/// How the camera moved since the last frame. **Logged separately from the playhead delta**,
/// because §4.5 is explicit that every measurement taken during a gesture has two causes and a
/// churn spike that cannot be attributed to one of them is not evidence.
#[derive(Clone, Copy, Debug, Default)]
pub struct CameraDelta {
    /// Centre motion in world units.
    pub d_centre: f64,
    /// The comparable one: centre motion in units of the **new** pixel size.
    pub d_centre_px: f64,
    /// `log2(half_old / half_new)`. Positive is a zoom **in**.
    pub d_zoom_octaves: f64,
    /// Frontier entries whose **stored** priority was recomputed. **Zero on a pure pan** — the
    /// stored term is position-free — and every entry on a zoom, because structure is
    /// pixel-relative. That asymmetry is the measurement §4.5 asks for.
    pub restored: usize,
    /// Levels the root grew to contain the camera. Zero unless `Regrow::Upto(k)` is set **and**
    /// the camera actually left the root box — a zoom-in can never trigger it, which is what
    /// makes a non-zero value here evidence rather than bookkeeping.
    pub regrown: u32,
}

impl CameraDelta {
    pub fn moved(&self) -> bool {
        self.d_centre != 0.0 || self.d_zoom_octaves != 0.0
    }
}

/// Whether the root box may grow to contain a camera that has left it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Regrow {
    /// The root is fixed; a camera outside it simply sees nothing new.
    #[default]
    Off,
    /// Grow the root upward at most `k` levels, keeping the old subtree.
    Upto(u32),
}

pub struct SessionCfg {
    pub sched: SchedCfg,
    pub quota: FrameQuota,
    /// Resident payloads. `None` keeps everything — the pre-session behaviour.
    pub payload_cap: Option<usize>,
    /// Viewport margin, in quad-widths, for the **eviction** relevance.
    ///
    /// Deliberately separate from `SchedCfg::camera_bias`, which is the **priority** margin:
    /// sharing one number would make a change to the eviction policy silently change the ranking.
    pub evict_margin: f64,
    pub regrow: Regrow,
    /// Frames between `Frontier::agrees_with_rebuild`. `0` means never, and is named as such
    /// rather than left to mean "every frame" by accident.
    pub audit_every: u32,
    /// Rank `pending` by frame priority before the quota truncates it.
    ///
    /// **`false` is the named control and it is the pre-wiring behaviour**: the quota keeps
    /// whichever quads happen to be first in the `Vec`, which is split order. Kept reachable
    /// because the question the wiring has to answer is whether ranking moves a tree at all --
    /// and the answer is not obvious. Children of one parent all inherit that parent's stored
    /// term, so a round whose `pending` comes from a single split is a **total tie**, ties break
    /// by id ascending, and id order *is* split order. The ranking can only bite where `pending`
    /// spans parents that disagree. Measured in `results/session/`, never assumed.
    pub rank_frame: bool,
}

impl Default for SessionCfg {
    fn default() -> Self {
        SessionCfg {
            sched: SchedCfg::default(),
            quota: FrameQuota::default(),
            payload_cap: None,
            evict_margin: 0.5,
            regrow: Regrow::Off,
            audit_every: 16,
            rank_frame: true,
        }
    }
}

/// **Refuse to write session telemetry from a session that is a one-shot descent with logging.**
///
/// `payload_cap: None` and `regrow: Off` together mean nothing is ever evicted and the root never
/// grows: every mechanism a session adds is inert, and a dump written from that cell is
/// `descend_with` wearing the session's name. That is exactly `k_frac = 1.0` shipping as the
/// default — a configuration that silently reproduces the old behaviour — and the project's answer
/// to that class is a guard, not a convention.
///
/// `allow` is explicit so the degenerate cell stays reachable as a **named control**, which it has
/// to be: it is the arm that says what a session buys.
pub fn assert_session_engages(cfg: &SessionCfg, path: &str, allow: bool) {
    let inert = cfg.payload_cap.is_none() && cfg.regrow == Regrow::Off;
    if inert && !allow && path.contains("results/") {
        panic!(
            "session telemetry to `{path}` from the degenerate cell: payload_cap None, regrow Off. \
             Nothing a session adds can engage here, so this is a one-shot descent with a frame \
             counter. Pass `allow = true` to keep it as the named control."
        );
    }
}

/// A quadtree, a camera, and a store, across frames.
pub struct Session {
    cfg: SessionCfg,
    ds: DescentState,
    cam: Camera,
    t_max: f64,
    frame: u64,
    /// The recorded boundary the tree is currently decided at.
    playhead: usize,
    /// The zoom the stored priority terms were computed at. A change here invalidates every
    /// ranking and no physics.
    stored_zoom: f64,
    /// The persistent frontier over `ds.pending`, maintained incrementally across frames.
    fr: Frontier,
    /// Entries the last frame's banded walk actually scored, and how many it held. `scan/len` is
    /// the statistic that decides whether the bucketing earns its place, and a frontier that
    /// could not report it could not settle that — this project has shipped two mechanisms that
    /// computed, sorted and changed nothing.
    last_scan: usize,
    last_len: usize,
    /// `1.0` / `0.0` from the last audit, or **`NaN` on a frame where it did not run**. A check
    /// that did not run must not report a pass.
    last_agrees: f64,
    /// Levels the root has grown, against `Regrow::Upto(k)`.
    grown: u32,
}

impl Session {
    pub fn new(cx: f64, cy: f64, half: f64, body: usize, cam: Camera, cfg: SessionCfg, t_max: f64) -> Session {
        let retain = match cfg.payload_cap {
            Some(c) => Retain::Capped(c),
            None if cfg.sched.keep_pixels => Retain::All,
            None => Retain::Never,
        };
        let mut ds = DescentState::new(cx, cy, half, body, &cfg.sched);
        ds.px = PixelStore::new(retain);
        let stored_zoom = cam.half_world;
        Session {
            cfg,
            ds,
            cam,
            t_max,
            frame: 0,
            playhead: 0,
            stored_zoom,
            fr: Frontier::new(),
            last_scan: 0,
            last_len: 0,
            last_agrees: f64::NAN,
            grown: 0,
        }
    }

    pub fn tree(&self) -> &crate::quad::QuadTree {
        &self.ds.tree
    }

    pub fn stats(&self) -> &scheduler::SchedStats {
        &self.ds.stats
    }

    pub fn store(&self) -> &PixelStore {
        &self.ds.px
    }

    /// The scheduler config the session is running, **including the camera it keeps in sync**.
    /// `set_camera` writes `sched.camera`, so `priority` and the frontier read the same camera the
    /// veto does — two cameras drifting apart would be an invisible staleness of exactly the kind
    /// the frontier audit exists to catch, one level up.
    pub fn sched_cfg(&self) -> &SchedCfg {
        &self.cfg.sched
    }

    /// Quads queued for the next round. Uncomputed by construction — their reductions are default
    /// until a round reaches them.
    pub fn pending(&self) -> &[usize] {
        &self.ds.pending
    }

    pub fn frame(&self) -> u64 {
        self.frame
    }

    /// Move the camera. **Cheap and total: never computes, never evicts, never splits.**
    ///
    /// It classifies the move and re-derives exactly what the move invalidated. A **pan** touches
    /// nothing stored — the stored term is position-free, which is the whole reason the frontier
    /// splits stored from derived. A **zoom** invalidates every stored ranking term and no physics,
    /// because structure is pixel-relative. `restored` reports which happened, and reporting `0` on
    /// a pan is the §4.5 measurement rather than an omission.
    pub fn set_camera(&mut self, cam: Camera) -> CameraDelta {
        let old = self.cam;
        let d_centre = ((cam.cx - old.cx).powi(2) + (cam.cy - old.cy).powi(2)).sqrt();
        let d_zoom_octaves = if cam.half_world > 0.0 && old.half_world > 0.0 {
            (old.half_world / cam.half_world).log2()
        } else {
            0.0
        };
        self.cam = cam;
        self.cfg.sched.camera = Some(cam);
        let zoomed = cam.half_world != self.stored_zoom;
        let restored = if zoomed {
            self.stored_zoom = cam.half_world;
            // Every leaf's ranking is stale; no payload is. The frontier is re-scored here rather
            // than counted and forgotten: `reprioritise` re-buckets only entries that crossed a
            // band, which is the whole reason a bucketed frontier beats a re-sort.
            let ids: Vec<usize> = self.fr.entries().into_iter().map(|(i, _)| i).collect();
            for i in ids {
                let v = scheduler::stored_term(&self.ds.tree, i, &self.cfg.sched);
                self.fr.reprioritise(i, v);
            }
            self.ds.tree.leaves().count()
        } else {
            0
        };
        let regrown = self.regrow_for_camera();
        CameraDelta {
            d_centre,
            d_centre_px: if cam.pixel_size() > 0.0 { d_centre / cam.pixel_size() } else { 0.0 },
            d_zoom_octaves,
            restored,
            regrown,
        }
    }

    /// **Grow the root upward until it contains the camera, or the bound is reached.**
    ///
    /// A zoom-out past the root box needs area the tree does not hold. Discarding and re-rooting
    /// would throw away exactly the quads the zoom-out is about to display — turning the one
    /// gesture the caching contract calls *nearly free* into the most expensive in the system. So
    /// the root grows: `grow_root` pushes the new root at the **end**, so no index moves and
    /// neither the store nor the frontier is disturbed.
    ///
    /// The three new siblings are queued, not computed — `set_camera` stays free of integration,
    /// and `step` spends the frame's quota on them like any other pending quad. They may well be
    /// vetoed instead, which is correct: a quad below the screen floor is not worth a trajectory
    /// however new it is.
    ///
    /// **Past the bound the camera is clamped in the record, not in the camera.** `regrown` stops
    /// rising while `Camera::half_world` keeps growing, so a harness can see the clamp; silently
    /// pinning the camera would make a zoom-out stop working with nothing saying why.
    fn regrow_for_camera(&mut self) -> u32 {
        let Regrow::Upto(max) = self.cfg.regrow else { return 0 };
        let mut made = 0u32;
        while self.grown < max {
            let r = self.ds.tree.root_node();
            let (rx, ry, rh) = (r.cx, r.cy, r.half);
            // Contained? The camera's own box, not its centre: a camera whose centre is inside
            // but whose viewport spills over still has nothing to show at the edge.
            let (c, h) = (self.cam, self.cam.half_world);
            if (c.cx - h) >= (rx - rh) && (c.cx + h) <= (rx + rh)
                && (c.cy - h) >= (ry - rh) && (c.cy + h) <= (ry + rh)
            {
                break;
            }
            // Grow TOWARD the camera: the old root becomes the child on the far side.
            // `child_boxes` is (lower-left, lower-right, upper-left, upper-right), and
            // `grow_root(q)` puts the old root at `q`, so the new root's centre moves opposite.
            let quadrant = usize::from(c.cx < rx) + 2 * usize::from(c.cy < ry);
            // Index: 0 = old root lower-left  => grows up-right   (camera right and up)
            //        1 = old root lower-right => grows up-left    (camera left  and up)
            //        2 = old root upper-left  => grows down-right (camera right and down)
            //        3 = old root upper-right => grows down-left  (camera left  and down)
            let kids = self.ds.tree.grow_root(quadrant, self.ds.iteration);
            self.ds.pending.extend_from_slice(&kids);
            self.grown += 1;
            made += 1;
        }
        made
    }

    /// Levels the root has grown over the session's life.
    pub fn grown(&self) -> u32 {
        self.grown
    }

    pub fn camera(&self) -> Camera {
        self.cam
    }

    /// What a frame did to the store.
    ///
    /// `recomputed` is kept apart from `computed` on purpose: one is progress and the other is the
    /// price of the cap, and a session reporting their sum reports a budget being spent well while
    /// it thrashes.
    pub fn evict_after_frame(&mut self, keep: &dyn Fn(usize) -> bool) -> Vec<usize> {
        let Some(cap) = self.cfg.payload_cap else { return Vec::new() };
        let (cam, margin) = (self.cam, self.cfg.evict_margin);
        let boxes: Vec<(f64, f64, f64)> =
            self.ds.tree.nodes.iter().map(|q| (q.cx, q.cy, q.half)).collect();
        let steps: Vec<u64> =
            self.ds.tree.nodes.iter().map(|q| q.red.total_substeps as u64).collect();
        let frame = self.frame;
        self.ds.px.evict_to(
            cap,
            frame,
            keep,
            &|i| boxes.get(i).map_or(0.0, |&(x, y, h)| cam.relevance(x, y, h, margin)),
            &|i| steps.get(i).copied().unwrap_or(0),
        )
    }

    /// **Advance the playhead, re-deciding every quad from the boundary it has now reached.**
    ///
    /// Under the scheduler contract's lockstep loop the playhead moves one fixed `dt` per frame
    /// whether or not the frame spends its quota, and every already-computed quad must be re-read
    /// at the new time — not left holding the reduction it had at an earlier boundary. Without
    /// this a session is a **camera-only** march: the field never changes, the tree converges once
    /// and then sits, and the depth-variance curve is flat because **nothing moves**. That is
    /// *frozen*, not *balanced*, and §3.2 is explicit that the two are indistinguishable in a
    /// variance plot alone — the churn column is what tells them apart, and it read 0.0000 for nine
    /// consecutive frames before this existed.
    ///
    /// Requires `keep_live_series`. Returns the number of quads re-reduced, so a caller can see
    /// the arm is live rather than assume it.
    pub fn set_playhead(&mut self, j: usize) -> usize {
        let (n, tau, hot, t_max) =
            (self.cfg.sched.n, self.cfg.sched.tau_display, self.cfg.sched.hot_rule, self.t_max);
        let mut moved = 0usize;
        for i in 0..self.ds.tree.nodes.len() {
            let px = self.ds.px.get(i);
            // A quad whose payload is absent or evicted keeps the reduction it has: the caching
            // contract's resume point, and re-reducing from nothing would be worse than stale.
            if px.is_empty() || j >= px[0].live_t.len() {
                continue;
            }
            let projected: Vec<crate::ensemble::pixel::PixelOut> =
                px.iter().map(|p| scheduler::project_at(p, j)).collect();
            self.ds.tree.nodes[i].red = scheduler::reduce(&projected, n, tau, hot, t_max);
            moved += 1;
        }
        if moved > 0 {
            // The frontier must be re-decided at the new time, so leaves that were settled are
            // eligible again. `decide` is pure on the reduction, so this costs nothing.
            // **Excluding anything already pending.** `round`'s frontier is `pending` chained
            // with `deferred`, so a quad in both is decided twice, enters `want` twice and is
            // split twice -- "quad 26 already split". Before this, `deferred` only ever held quads
            // the ranking had dropped, which are by construction not pending; re-deciding *every*
            // leaf breaks that invariant unless the overlap is removed here.
            let pending: std::collections::HashSet<usize> =
                self.ds.pending.iter().cloned().collect();
            let leaves: Vec<usize> =
                self.ds.tree.leaves().filter(|i| !pending.contains(i)).collect();
            self.ds.deferred.extend(leaves);
            self.ds.deferred.sort_unstable();
            self.ds.deferred.dedup();
        }
        self.playhead = j;
        moved
    }

    pub fn playhead(&self) -> usize {
        self.playhead
    }

    /// **Bring the frontier level with `pending`, then move the top `k` to the front.**
    ///
    /// This is the site the frontier measurement pointed at. `order_queue` ranks the quads that
    /// *want to split*; the **frame quota** truncates `pending`, the quads about to be *computed*,
    /// and it truncated them by position — which is split order, i.e. raster order within a
    /// parent. So a frame spent its budget on whichever children happened to be first in a `Vec`.
    ///
    /// Ranking here is what the persistent frontier is for, and `k` is the frame quota: the
    /// small-`k` regime where the banded walk beats a sort by 83–94%, against 0–67% at
    /// `k_frac = 0.25`. `descend_with` was measured to be the wrong site for exactly that reason
    /// and is left alone.
    ///
    /// **Inert when the quota does not bind**, and it returns early to say so rather than
    /// reordering a list nothing will truncate — a ranking applied where nothing is dropped is
    /// `k_frac = 1.0` again.
    fn rank_pending(&mut self, room: usize) -> usize {
        if !self.cfg.rank_frame {
            self.last_scan = 0;
            self.last_len = 0;
            self.last_agrees = f64::NAN;
            return 0;
        }
        // Sync incrementally: the frontier persists across frames precisely because held-over
        // pending does. Removing computed ids and inserting new ones is O(n) with no sort.
        let held: std::collections::HashSet<usize> =
            self.fr.entries().into_iter().map(|(i, _)| i).collect();
        let now: std::collections::HashSet<usize> = self.ds.pending.iter().cloned().collect();
        for i in held.difference(&now) {
            self.fr.remove(*i);
        }
        for &i in &self.ds.pending {
            if !held.contains(&i) {
                let v = scheduler::stored_term(&self.ds.tree, i, &self.cfg.sched);
                self.fr.insert(i, v);
            }
        }
        self.last_len = self.fr.len();

        // **The audit, on the population the frame actually ranks.** A staleness check run on a
        // different set from the one selected would be answering a different question.
        if self.cfg.audit_every > 0 && self.frame % u64::from(self.cfg.audit_every) == 0 {
            let (tree, cfg) = (&self.ds.tree, &self.cfg.sched);
            let d = |i: usize| scheduler::derived_term(tree, i, cfg);
            self.last_agrees = f64::from(u8::from(self.fr.agrees_with_rebuild(room.max(1), d)));
        } else {
            self.last_agrees = f64::NAN;
        }

        if self.ds.pending.len() <= room {
            self.last_scan = 0;
            return 0;
        }
        let (top, scanned) = {
            let (tree, cfg) = (&self.ds.tree, &self.cfg.sched);
            self.fr.top_k_bounded(room, |i| scheduler::derived_term(tree, i, cfg))
        };
        let set: std::collections::HashSet<usize> = top.iter().cloned().collect();
        // The remainder keeps its relative order, so a frame that holds work over hands the next
        // one a stable tail rather than a reshuffled one.
        let rest: Vec<usize> =
            self.ds.pending.iter().cloned().filter(|i| !set.contains(i)).collect();
        self.ds.pending = top.into_iter().chain(rest).collect();
        self.last_scan = scanned;
        scanned
    }

    /// Entries the last frame's banded walk scored, the frontier's length, and whether the audit
    /// agreed (`NaN` if it did not run).
    pub fn frontier_telemetry(&self) -> (usize, usize, f64) {
        (self.last_scan, self.last_len, self.last_agrees)
    }

    /// **One frame.** Run rounds until the quota binds, then evict to the cap.
    ///
    /// Returns what stopped it and what it spent. The tree is *not* rebuilt: this continues the
    /// descent the session has been running since `new`, which is what makes a quad computed on
    /// frame 3 still there on frame 300.
    pub fn step(&mut self, sampler: scheduler::Sampler<'_>) -> (QuotaHit, Spend) {
        self.frame += 1;
        let mut spend = Spend::default();
        let stop = Stop::Frame { quads: self.cfg.quota.quads, substeps: self.cfg.quota.substeps };
        let hit;

        loop {
            if self.ds.pending.is_empty() && self.ds.deferred.is_empty() {
                hit = QuotaHit::Drained;
                break;
            }
            if spend.rounds >= self.cfg.quota.rounds {
                hit = QuotaHit::Rounds;
                break;
            }
            if spend.quads >= self.cfg.quota.quads {
                hit = QuotaHit::Quads;
                break;
            }
            if self.cfg.quota.substeps.is_some_and(|s| spend.substeps >= s) {
                hit = QuotaHit::Substeps;
                break;
            }
            if self.ds.stats.budget_exhausted {
                hit = QuotaHit::TotalBudget;
                break;
            }
            // **Before the round, because the round's quota truncation is what it feeds.**
            self.rank_pending(self.cfg.quota.quads.saturating_sub(spend.quads));
            let before = spend.quads;
            scheduler::round(&mut self.ds, &self.cfg.sched, self.t_max, sampler, stop, &mut spend);
            // A round that computed nothing and left nothing pending has drained; without this a
            // frame with an empty frontier would spin to its round cap.
            if spend.quads == before && self.ds.pending.is_empty() {
                hit = QuotaHit::Drained;
                break;
            }
        }
        (hit, spend)
    }
}

