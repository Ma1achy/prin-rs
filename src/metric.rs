//! §2 — a metric to judge refinement criteria by, and the cache that makes it affordable.
//!
//! Nothing before this measured whether a criterion is *good*. Leaf counts, depths and `alpha`
//! distributions describe what a criterion **did**, not whether it was right.
//!
//! # The reframe, and why this is a ranking problem
//!
//! The screen floor stops 61.2% of near-field's leaves, so in normal use the criterion is not
//! deciding *when to stop* — it is deciding **which quads get the budget first, within what is
//! displayable.** That is a ranking problem, and ranking tolerates noise that would wreck a
//! threshold.
//!
//! It also disposes of a confound found while measuring §1: the between-footprint arm runs
//! 1.17x the within arm in `near-field` and **9.56x** in `far`, so swapping criterion at a fixed
//! `tau` silently changes the effective threshold by up to 8x per region. **A ranking is
//! invariant under any monotone rescaling of the signal**, so comparing criteria by rank costs
//! nothing to that factor, where comparing them by threshold would have scored the rescaling.
//!
//! # Precompute once, replay many
//!
//! One integration pass per region builds a **complete uniform tree** to the screen floor. After
//! that every criterion, both controls, the whole `error(B)` curve and the panning study are
//! traversals of the cache with **no re-integration**.
//!
//! Two facts make that exact rather than an approximation:
//!
//! - **Quads are disjoint**, so refining a leaf changes only the pixels it covers. Each quad's
//!   contribution to the image error is therefore a *constant*, [`CachedQuad::err_sum`],
//!   independent of what the rest of the tree does. The greedy replay becomes a static priority
//!   queue rather than a re-evaluation per step.
//! - **The reference colouring is `E`-independent.** It is the nominal copy's outcome, and copy
//!   0 is never jittered, so the reference image does not move with the ensemble size. Only the
//!   *signals* depend on `E`. That removes the brief's reference-resolution caveat for the base
//!   metric outright; it returns the moment the colouring is an SSAA resolve, and is stated
//!   again there.
//!
//! # What `error = 0` means, and what it does not
//!
//! The reference is the fully-refined tree at **one sample per pixel** — a specific finite
//! sampling, not the true image. At the screen floor sub-pixel structure is sampled arbitrarily:
//! which side of a filament a pixel lands on is an accident of where its sample fell. So
//! `error = 0` means **"matches this sampling"**, not "correct".
//!
//! That is the right target for *comparing* criteria, which need a common yardstick rather than
//! truth, and the exactly-locatable zero is a real virtue. It is not a statement about image
//! quality, and every table that quotes the curve says so.

use rayon::prelude::*;
use std::collections::HashMap;

use crate::ensemble::pixel::{evaluate, EnsembleCfg, PixelOut};
use crate::grid::{Chart, Slice};
use crate::output::oklab;
use crate::output::png::outcome_rgb;
use crate::quad::{Agg, Criterion, StructureMode, QuadReduction};
use crate::rng::SplitMix64;
use crate::scheduler::reduce;
use crate::spatial::HotRule;

/// How footprints become pixels.
///
/// **This is a criterion parameter, not a presentation one.** `error(B)` measures image change,
/// so changing what is displayed changes which quads matter. §6's coupling question is exactly
/// whether the curve moves when lightness switches from spread to diffusion: if it does, the
/// criterion needs a term for the lightness field and has none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Colouring {
    /// The nominal copy's outcome class. **`E`-independent**, because copy 0 is never
    /// jittered — which is what removes the brief's reference-resolution caveat for the base
    /// metric. The caveat returns for any colouring that reads the ensemble.
    Outcome,
    /// Hue from the shape sphere by vMF site-blend, lightness from a scalar. The production
    /// scheme. See [`crate::output::colour`].
    Bivariate(crate::output::colour::Scalar),
    /// **The event class on viridis** — the categorical mode the reference's WebGPU panel
    /// renders, and the one a reference comparison must be made under. Reads the ensemble
    /// (the class is joined with each copy's terminal outcome), so the brief's
    /// reference-resolution caveat applies to it and not to [`Colouring::Outcome`].
    EventClass,
}

impl Colouring {
    pub fn name(self) -> String {
        match self {
            Colouring::Outcome => "outcome".into(),
            Colouring::Bivariate(l) => format!("bivariate/{}", l.name()),
            Colouring::EventClass => "event_class/viridis".into(),
        }
    }
}

/// A quad's address in the complete tree: level and integer position within it.
pub type Key = (u32, u32, u32);

/// One quad of the precomputed tree. **The footprints are not kept** — only what a criterion or
/// the metric can read, which is why a full tree costs megabytes rather than gigabytes.
#[derive(Clone, Debug)]
pub struct CachedQuad {
    pub key: Key,
    pub red: QuadReduction,
    /// The `N x N` sample colours, in the grid's own order.
    pub rgb: Vec<[u8; 3]>,
    /// **Sum** of per-pixel OKLab distance to the reference over this quad's screen footprint,
    /// were it drawn as a leaf. A constant of the quad, which is what makes the greedy replay
    /// a static priority queue.
    pub err_sum: f64,
}

#[derive(Clone)]
pub struct Cache {
    pub region: String,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    pub body: usize,
    /// The chart this tree was integrated on. Held so a dump can record its **parameters**, not
    /// only its name: a `Plane`'s basis and a `Latent`'s `(z0, q1, q2)` are free, so two dumps
    /// with the same chart name can be different configurations.
    pub chart: Chart,
    /// Deepest level present. `2^levels * n == res`, so the leaves are exactly one sample per
    /// pixel — asserted at construction, not assumed.
    pub levels: u32,
    pub n: usize,
    pub res: usize,
    pub quads: HashMap<Key, CachedQuad>,
    /// The reference image, RGB8.
    pub reference: Vec<u8>,
    pub trajectories: u64,
    pub colouring: Colouring,
    /// The `[p1, p99]` the lightness ramp was normalised against, region-wide.
    pub ramp: (f64, f64),
    /// The hue sites, **fixed for the whole cache**. A colouring has to be one map: if the site
    /// set were rebuilt per footprint from that footprint's own masses, two pixels with
    /// different masses would be read against different palettes and the image would not be a
    /// picture of anything. Built from the chart's nominal masses.
    pub sites: crate::output::colour::SiteSet,
}

impl Cache {
    pub fn get(&self, k: Key) -> &CachedQuad {
        &self.quads[&k]
    }

    /// The four children of `k`, in `(jy, jx)` order to match `Quad::child_boxes`.
    pub fn children(k: Key) -> [Key; 4] {
        let (l, ix, iy) = k;
        [
            (l + 1, 2 * ix, 2 * iy),
            (l + 1, 2 * ix + 1, 2 * iy),
            (l + 1, 2 * ix, 2 * iy + 1),
            (l + 1, 2 * ix + 1, 2 * iy + 1),
        ]
    }

    /// What refining `k` buys: the drop in total image error, in OKLab-distance units summed
    /// over pixels. Zero at the deepest level, where there is nothing to refine into.
    pub fn gain(&self, k: Key) -> f64 {
        if k.0 >= self.levels {
            return 0.0;
        }
        let kids: f64 = Self::children(k).iter().map(|c| self.get(*c).err_sum).sum();
        self.get(k).err_sum - kids
    }

    /// Same-level edge neighbour, or `None` outside the root box.
    pub fn neighbour(&self, k: Key, dir: crate::quad::Dir) -> Option<Key> {
        let (l, ix, iy) = k;
        let w = 1u32 << l;
        let (nx, ny) = match dir {
            crate::quad::Dir::NegX => (ix.checked_sub(1)?, iy),
            crate::quad::Dir::PosX => (ix + 1, iy),
            crate::quad::Dir::NegY => (ix, iy.checked_sub(1)?),
            crate::quad::Dir::PosY => (ix, iy + 1),
        };
        (nx < w && ny < w).then_some((l, nx, ny))
    }

    /// §3.3's contrast, at matched level. `NaN` on a quad with no in-box neighbour.
    pub fn contrast(&self, k: Key, criterion: Criterion, agg: Agg) -> f64 {
        let me = self.get(k).red.signal(criterion, agg);
        if !me.is_finite() {
            return f64::NAN;
        }
        let (mut best, mut seen) = (0.0f64, 0u8);
        for d in crate::quad::Dir::ALL {
            if let Some(nk) = self.neighbour(k, d) {
                let v = self.get(nk).red.signal(criterion, agg);
                if v.is_finite() {
                    best = best.max((me - v).abs());
                    seen += 1;
                }
            }
        }
        if seen == 0 {
            f64::NAN
        } else {
            best
        }
    }

    /// Mean per-pixel error of a tree given by its leaf set.
    pub fn error_of(&self, leaves: &[Key]) -> f64 {
        let s: f64 = leaves.iter().map(|k| self.get(*k).err_sum).sum();
        s / (self.res * self.res) as f64
    }
}

impl Cache {
    /// Rasterise a tree given by its leaf set, at true per-quad texel sizes.
    ///
    /// One sample, one tile, no interpolation — a level-3 leaf's texels are 4x the linear size
    /// of a level-5 leaf's. Never upsampled smoothly, which would fabricate structure that was
    /// not sampled.
    pub fn render(&self, leaves: &[Key]) -> Vec<u8> {
        // Background is `BACKGROUND`, never black. A NaN shape renders `DEBUG_NAN`, and both
        // must differ from "no leaf covered this pixel" -- previously a NaN pixel came out
        // `[0,0,0]` from the `u8` cast and was bitwise identical to the background, so
        // `err_sum` scored an undetermined pixel as a perfect match.
        let mut img: Vec<u8> = crate::output::colour::BACKGROUND
            .iter()
            .cloned()
            .cycle()
            .take(self.res * self.res * 3)
            .collect();
        for &k in leaves {
            let (l, ix, iy) = k;
            let span = self.res >> l;
            let tile = span / self.n;
            let (px0, py0) = (ix as usize * span, iy as usize * span);
            let q = self.get(k);
            for dy in 0..span {
                for dx in 0..span {
                    let c = q.rgb[(dy / tile) * self.n + (dx / tile)];
                    let o = ((py0 + dy) * self.res + px0 + dx) * 3;
                    img[o] = c[0];
                    img[o + 1] = c[1];
                    img[o + 2] = c[2];
                }
            }
        }
        img
    }

    /// The same render with the leaf wireframe drawn over it.
    ///
    /// Emitted **beside** the plain render, never instead of it: the texel size says what is
    /// displayed and the wire says where the tree cut, and neither substitutes for the other.
    pub fn render_wire(&self, leaves: &[Key]) -> Vec<u8> {
        let mut img = self.render(leaves);
        let boxes: Vec<crate::output::wire::Box2> = leaves
            .iter()
            .map(|&(l, ix, iy)| {
                let span = (self.res >> l) as f64;
                crate::output::wire::Box2 {
                    x0: ix as f64 * span,
                    y0: iy as f64 * span,
                    x1: (ix + 1) as f64 * span,
                    y1: (iy + 1) as f64 * span,
                    level: l,
                }
            })
            .collect();
        let deepest = leaves.iter().map(|k| k.0).max().unwrap_or(0);
        crate::output::wire::draw(&mut img, self.res, self.res, &boxes, deepest.max(1));
        img
    }

    /// The leaf set a replay holds at a given budget — what `render` needs to draw it.
    pub fn leaves_at(&self, rank: Rank, budget: usize) -> Vec<Key> {
        let pts = replay_with_leaves(self, rank, budget);
        pts.1
    }
}

/// Quad geometry from a key.
fn box_of(c: &Cache, k: Key) -> (f64, f64, f64) {
    let (l, ix, iy) = k;
    let h = c.half / (1u64 << l) as f64;
    (
        c.cx - c.half + (2 * ix + 1) as f64 * h,
        c.cy - c.half + (2 * iy + 1) as f64 * h,
        h,
    )
}

/// Build the complete tree to `levels`, colouring by the nominal copy's outcome.
///
/// **Cost is the whole of this stage**: `(4^(levels+1) - 1)/3` quads, each `N^2 * (E+1)`
/// trajectories. Every replay afterwards is free.
pub fn build(
    region: &str,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    chart: Chart,
    levels: u32,
    n: usize,
    res: usize,
    tau: f64,
    ens: &EnsembleCfg,
    colouring: Colouring,
) -> Cache {
    build_multi(region, cx, cy, half, body, chart, levels, n, res, tau, ens, &[colouring])
        .pop()
        .unwrap()
}

/// As [`build_multi_with_footprints`], discarding the footprints.
#[allow(clippy::too_many_arguments)]
pub fn build_multi(
    region: &str,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    chart: Chart,
    levels: u32,
    n: usize,
    res: usize,
    tau: f64,
    ens: &EnsembleCfg,
    colourings: &[Colouring],
) -> Vec<Cache> {
    build_multi_with_footprints(region, cx, cy, half, body, chart, levels, n, res, tau, ens, colourings).0
}

/// As [`build`], for several colourings over **one** integration pass.
///
/// The footprints are integrated once and coloured several ways. That matters for §6: the
/// lightness field costs a second, fixed-step, unregularised march per footprint, and paying it
/// three times to answer whether the *colouring* changes `error(B)` would be paying for the
/// thing being held fixed.
#[allow(clippy::too_many_arguments)]
pub fn build_multi_with_footprints(
    region: &str,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
    chart: Chart,
    levels: u32,
    n: usize,
    res: usize,
    tau: f64,
    ens: &EnsembleCfg,
    colourings: &[Colouring],
) -> (Vec<Cache>, HashMap<Key, Vec<PixelOut>>) {
    assert!(!colourings.is_empty());
    assert_eq!(
        (1usize << levels) * n,
        res,
        "the deepest level must be exactly one sample per pixel: 2^{levels} * {n} != {res}"
    );

    let mut c = Cache {
        region: region.to_string(),
        cx,
        cy,
        half,
        body,
        chart,
        levels,
        n,
        res,
        quads: HashMap::new(),
        reference: crate::output::colour::BACKGROUND
            .iter()
            .cloned()
            .cycle()
            .take(res * res * 3)
            .collect(),
        trajectories: 0,
        colouring: colourings[0],
        ramp: (0.0, 1.0),
        sites: crate::output::colour::landmarks(&crate::physics::burrau::MASSES),
    };

    // ---- integrate every quad at every level ----
    let mut all: Vec<Key> = Vec::new();
    for l in 0..=levels {
        let w = 1u32 << l;
        for iy in 0..w {
            for ix in 0..w {
                all.push((l, ix, iy));
            }
        }
    }

    let computed: Vec<(Key, QuadReduction, Vec<PixelOut>)> = all
        .par_iter()
        .map(|&k| {
            let (qx, qy, qh) = box_of(&c, k);
            let slice = Slice::body_plane(n, n, qx, qy, qh, body).with_chart(chart);
            let px: Vec<PixelOut> =
                (0..slice.npix()).map(|i| evaluate::<f64>(&slice, i, ens)).collect();
            let mut red = reduce(&px, n, tau, HotRule::default(), ens.t_max);
            let ics: Vec<crate::physics::Cart<f64>> =
                (0..slice.npix()).map(|i| slice.nominal::<f64>(i)).collect();
            red.n_distinct_ic = crate::decode::distinct(&ics) as u32;
            (k, red, px)
        })
        .collect();

    c.trajectories = computed.len() as u64 * (n * n) as u64 * (ens.n_extra + 1) as u64;

    let base_trajectories = computed.len() as u64 * (n * n) as u64 * (ens.n_extra + 1) as u64;

    let px_of: HashMap<Key, Vec<PixelOut>> =
        computed.iter().map(|(k, _, px)| (*k, px.clone())).collect();

    let mut out: Vec<Cache> = Vec::with_capacity(colourings.len());
    for &colouring in colourings {
        let mut c = Cache { colouring, trajectories: base_trajectories, ..c.clone() };
        for (k, red, _px) in &computed {
            c.quads.insert(*k, CachedQuad { key: *k, red: *red, rgb: Vec::new(), err_sum: 0.0 });
        }
        repaint(&mut c, &px_of);
        out.push(c);
    }

    (out, px_of)
}

/// Colour every quad, build the reference image, and compute each quad's `err_sum`.
///
/// Shared by [`build_multi`] and [`Cache::recolour`] so the two cannot drift: a replay that
/// coloured by a slightly different path would produce an `error(B)` curve that looked like a
/// measurement and was an artefact of the replay. `c.quads` must already hold every key with
/// its `red`; only `rgb`, `reference` and `err_sum` are written.
fn repaint(c: &mut Cache, px_of: &HashMap<Key, Vec<PixelOut>>) {
    // **The ramp is normalised over the whole region, once.** Per-quad normalisation would make
    // a quad's colour depend on which quads happen to be leaves, so refining one quad would
    // change the colour of another and `err_sum` would stop being a constant of the quad --
    // which is the property the greedy replay rests on.
    let ramp = match c.colouring {
        Colouring::Outcome | Colouring::EventClass => (0.0, 1.0),
        Colouring::Bivariate(sc) => {
            let all_px: Vec<PixelOut> = px_of.values().flat_map(|v| v.iter().cloned()).collect();
            crate::output::colour::range(&all_px, sc)
        }
    };
    c.ramp = ramp;
    let sites = c.sites.clone();
    let colouring = c.colouring;
    for (k, px) in px_of {
        let rgb: Vec<[u8; 3]> = match colouring {
            Colouring::Outcome => px.iter().map(outcome_rgb).collect(),
            Colouring::EventClass => {
                px.iter().map(crate::output::png::event_class_rgb).collect()
            }
            Colouring::Bivariate(sc) => px
                .iter()
                .map(|p| crate::output::colour::rgb(p, sc, &sites, ramp.0, ramp.1))
                .collect(),
        };
        if let Some(q) = c.quads.get_mut(k) {
            q.rgb = rgb;
        }
    }

    // ---- the reference image, from the deepest level ----
    let (levels, n, res) = (c.levels, c.n, c.res);
    let w = 1u32 << levels;
    for iy in 0..w {
        for ix in 0..w {
            let q = &c.quads[&(levels, ix, iy)];
            for sy in 0..n {
                for sx in 0..n {
                    let px = ix as usize * n + sx;
                    let py = iy as usize * n + sy;
                    let rgb = q.rgb[sy * n + sx];
                    let o = (py * res + px) * 3;
                    c.reference[o] = rgb[0];
                    c.reference[o + 1] = rgb[1];
                    c.reference[o + 2] = rgb[2];
                }
            }
        }
    }

    // ---- each quad's error contribution, were it drawn as a leaf ----
    let keys: Vec<Key> = c.quads.keys().cloned().collect();
    let sums: Vec<(Key, f64)> = keys
        .par_iter()
        .map(|&k| {
            let (l, ix, iy) = k;
            let span = res >> l; // pixels across this quad
            let tile = span / n; // pixels per sample
            let (px0, py0) = (ix as usize * span, iy as usize * span);
            let q = &c.quads[&k];
            let mut acc = 0.0;
            for dy in 0..span {
                for dx in 0..span {
                    let s = q.rgb[(dy / tile) * n + (dx / tile)];
                    let o = ((py0 + dy) * res + px0 + dx) * 3;
                    acc += oklab::delta(s, [c.reference[o], c.reference[o + 1], c.reference[o + 2]]);
                }
            }
            (k, acc)
        })
        .collect();
    for (k, s) in sums {
        c.quads.get_mut(&k).unwrap().err_sum = s;
    }
}

impl Cache {
    /// Rebuild this cache under a different colouring, from a footprint file.
    ///
    /// **The point of `PRQF`.** `err_sum` is a function of the colouring, so PR #13's curves
    /// could not be recomputed under a new one without re-integrating 2.8 million trajectories
    /// per region. With the footprints on disk, `error(B)` under any colouring is a replay.
    ///
    /// The reductions carry over unchanged: a `QuadReduction` is a property of the physics and
    /// does not know what colour anything is drawn.
    pub fn recolour(
        &self,
        fp: &crate::output::fcache::Footprints,
        colouring: Colouring,
    ) -> Result<Cache, String> {
        fp.agrees_with(self)?;
        let mut c = Cache { colouring, ..self.clone() };
        let px_of: HashMap<Key, Vec<PixelOut>> = fp
            .quads
            .iter()
            .map(|(k, rows)| (*k, rows.iter().map(|r| r.to_pixel()).collect()))
            .collect();
        for k in c.quads.keys().cloned().collect::<Vec<_>>() {
            if !px_of.contains_key(&k) {
                return Err(format!("footprint file is missing quad {k:?}"));
            }
        }
        repaint(&mut c, &px_of);
        Ok(c)
    }

    /// Every footprint of this cache's tree, ready to write with [`crate::output::fcache::write`].
    ///
    /// Only available from [`build_multi_with_footprints`], because a `Cache` deliberately does
    /// not retain its footprints — that is what keeps a complete tree in megabytes.
    pub fn footprints_from(
        &self,
        px_of: &HashMap<Key, Vec<PixelOut>>,
        t_max: f64,
    ) -> crate::output::fcache::Footprints {
        crate::output::fcache::Footprints {
            region: self.region.clone(),
            chart: format!("{} {}", self.chart.name(), self.chart.params()),
            cx: self.cx,
            cy: self.cy,
            half: self.half,
            body: self.body,
            levels: self.levels,
            n: self.n,
            res: self.res,
            t_max,
            quads: px_of
                .iter()
                .map(|(k, px)| {
                    (*k, px.iter().map(crate::output::fcache::Row::of).collect())
                })
                .collect(),
        }
    }
}

/// How a replay chooses which leaf to refine next.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rank {
    /// A criterion's own signal, read as a **priority order** and never against a threshold.
    Signal(Criterion, Agg),
    /// Uniformly random. The floor: any criterion must beat this. Run several seeds and read
    /// the band, never one trace — a single random run is a draw.
    Random(u64),
    /// Breadth-first: refine the shallowest quad available. **The honest baseline**, and the
    /// one that was missing while every table read "beats random" as though random were the
    /// thing to beat.
    ///
    /// On a smooth field this is near-optimal, so `random` is the wrong comparison there — a
    /// criterion can be far above random and still buy nothing over refining uniformly.
    /// Measured on `far`, `dp_optimal` and every non-greedy row produce the *identical* leaf
    /// set, which is this.
    ///
    /// It is also what the rest of the enum **degenerates to** when its signal goes flat: the
    /// tie-break at the argmax is lexicographic on `(level, ix, iy)`, level first. Scoring
    /// `-level` here makes that explicit rather than emergent, so the baseline can be read as a
    /// row instead of inferred from a coincidence.
    Uniform,
    /// Greedy on immediate `Δerror`. **Neither optimal nor a bound**, and named for what it is
    /// after being read as a ceiling in every table for two PRs.
    ///
    /// Greedy is optimal only when gains are independent and immediately available, and here
    /// they are neither: a quad whose own split gains little may unlock children with large
    /// gains two levels down, and greedy declines it. That is the classic failure of greedy on
    /// a sequential tree problem. **A criterion beating this indicates lookahead value, not an
    /// error**, and there is deliberately no assertion anywhere that it dominates.
    ///
    /// **It can lose to random, and to breadth-first, and has been measured doing both.** On
    /// `far` at `B = 1535` it reads **0.54760** against a random band of **0.48550-0.52047** and
    /// every criterion at **0.36557** — the worst strategy in the table. `far` is smooth, so a
    /// quad's spread tracks its cell width and argmax-on-spread *is* breadth-first, which is
    /// near-optimal there; greedy chases fluctuations in an immediate `Δerror` that is noise at
    /// every level above the last, and concentrates the budget in a flat corner.
    ///
    /// The ceiling this was mistaken for is [`Cache::dp_optimal`], which is exact.
    GreedyLookahead1,
    /// Greedy on `Δerror / cost`, where cost is the quad's measured substeps. §8.
    GreedyLookahead1PerCost,
    /// §3.3 — `max` over the four edge-neighbours of `|signal_self - signal_neighbour|`.
    ///
    /// Interesting regions are where the signal **changes**, not where it is high: a uniformly
    /// chaotic quad and a uniformly smooth one are both featureless, and the boundary between
    /// them is the structure. It also **sidesteps `tau` entirely**, which matters because the
    /// vertical slice promoted `tau` to the dominant knob under the screen floor (64x on `far`).
    ///
    /// Over the cache this is better defined than over a live tree: every quad exists at every
    /// level, so the neighbour is read **at the same level** always, rather than falling back
    /// up-tree when it happens not to have been refined yet.
    Contrast(Criterion, Agg),
    /// §2.2 — the signal with the spatial-structure term applied, `Replace` or `Multiply`.
    ///
    /// **Enters as an ordering, like everything else here.** That matters more for this one than
    /// for most: `structure` is bounded in `[0, 1]` while `ensemble_spread` is not, so
    /// `Multiply` rescales the signal by up to its whole range. Compared against a threshold
    /// that rescaling would be scored instead of the structure; compared as an ordering it
    /// costs nothing, which is the standing rule about the between arm at 9.56x.
    Structured(StructureMode, Criterion, Agg),
    /// The structure term **alone**, with no signal in it at all — the control that says whether
    /// `Structured` is buying structure or just re-weighting the spread.
    StructureOnly,
}

impl Rank {
    pub fn name(self) -> String {
        match self {
            Rank::Signal(c, a) => format!("{}/{}", c.name(), a.name()),
            Rank::Random(s) => format!("random[{s}]"),
            Rank::Uniform => "uniform".into(),
            Rank::GreedyLookahead1 => "greedy_lookahead_1".into(),
            Rank::GreedyLookahead1PerCost => "greedy_lookahead_1/cost".into(),
            Rank::Contrast(c, a) => format!("contrast:{}/{}", c.name(), a.name()),
            Rank::Structured(m, c, a) => format!("{}x{}/{}", m.name(), c.name(), a.name()),
            Rank::StructureOnly => "structure_only".into(),
        }
    }
}

/// One point on an `error(B)` curve.
#[derive(Clone, Copy, Debug)]
pub struct Point {
    /// Quads computed so far, root included.
    pub budget: usize,
    pub leaves: usize,
    pub error: f64,
}

/// Replay a ranking over the cache, recording `error(B)` after every split.
///
/// The tree starts at the root and refines the top-ranked leaf until the budget runs out or
/// every leaf is at the deepest level. **No threshold is consulted** — a criterion enters here
/// purely as an ordering, which is what §2's reframe asks for and what makes the comparison
/// invariant to the region-dependent scale factor between the two arms.
pub fn replay(cache: &Cache, rank: Rank, budget: usize) -> Vec<Point> {
    replay_with_leaves(cache, rank, budget).0
}

/// As [`replay`], also returning the final leaf set so the tree can be drawn.
pub fn replay_with_leaves(cache: &Cache, rank: Rank, budget: usize) -> (Vec<Point>, Vec<Key>) {
    replay_ordered(cache, Order::Ranked(rank), budget)
}

/// What the replay ranks by: one of the enumerated [`Rank`]s, or an arbitrary per-quad score.
///
/// The second exists so a **fitted** signal -- a logistic combination of many reduction fields,
/// which cannot be a `Rank` variant because `Rank` is `Copy` and carries no weights -- goes
/// through the *same* argmax, the same non-finite convention and the same level-first tie-break
/// as every criterion it is being compared against. A separate loop would have made the fitted
/// row incomparable to the rows above it while looking like a fair comparison.
enum Order<'a> {
    Ranked(Rank),
    Scored(&'a HashMap<Key, f64>),
}

/// Replay an arbitrary per-quad score, for a signal that is not a [`Rank`].
///
/// A quad absent from the map scores `-inf`, the same convention a non-finite signal gets:
/// **undetermined never wins, and never blocks.**
pub fn replay_scored(
    cache: &Cache,
    score: &HashMap<Key, f64>,
    budget: usize,
) -> (Vec<Point>, Vec<Key>) {
    replay_ordered(cache, Order::Scored(score), budget)
}

fn replay_ordered(cache: &Cache, order: Order, budget: usize) -> (Vec<Point>, Vec<Key>) {
    let mut leaves: Vec<Key> = vec![(0, 0, 0)];
    let mut rng = match order {
        Order::Ranked(Rank::Random(s)) => Some(SplitMix64::new(s)),
        _ => None,
    };
    let mut spent = 1usize;
    let mut out = vec![Point { budget: spent, leaves: 1, error: cache.error_of(&leaves) }];

    loop {
        // Only leaves that can still be refined.
        let cand: Vec<usize> = (0..leaves.len()).filter(|&i| leaves[i].0 < cache.levels).collect();
        if cand.is_empty() || spent + 4 > budget {
            break;
        }
        let pick = match order {
            Order::Ranked(Rank::Random(_)) => {
                let r = rng.as_mut().unwrap().next_u64() as usize;
                cand[r % cand.len()]
            }
            _ => {
                // **Non-finite maps to -inf, not to itself.** `v > bv` is false whenever `bv`
                // is NaN, so a NaN at the first candidate would block every finite score
                // behind it and collapse the ranking to scan order — silently, and worst on
                // exactly the signals that decline to score part of the region (`term_grad`
                // is NaN on 97.1% of near-field). "Undetermined never wins" is the intended
                // semantics; "undetermined blocks everything" was the bug.
                let score = |i: usize| -> f64 {
                    let v = match order {
                        Order::Ranked(r) => raw_score(cache, leaves[i], r),
                        Order::Scored(m) => m.get(&leaves[i]).copied().unwrap_or(f64::NAN),
                    };
                    if v.is_finite() {
                        v
                    } else {
                        f64::NEG_INFINITY
                    }
                };
                // Ties broken by the earliest candidate, deterministically. A random tie-break
                // would make a criterion's curve depend on a seed it never sees.
                let mut best = cand[0];
                let mut bv = score(best);
                for &i in &cand[1..] {
                    let v = score(i);
                    if v > bv || (v == bv && leaves[i] < leaves[best]) {
                        best = i;
                        bv = v;
                    }
                }
                best
            }
        };

        let k = leaves.swap_remove(pick);
        leaves.extend_from_slice(&Cache::children(k));
        spent += 4;
        out.push(Point { budget: spent, leaves: leaves.len(), error: cache.error_of(&leaves) });
    }
    (out, leaves)
}

/// A ranking's raw score for one quad, before the non-finite convention is applied.
///
/// Public because **the distinct-value count of a ranking has to be readable before its curve
/// is**: a flat `error(B)` has two causes -- a bad ordering and no ordering -- and the curve
/// alone cannot tell them apart. `Random` scores `NaN` by construction and is not meaningful
/// here.
pub fn score(cache: &Cache, k: Key, rank: Rank) -> f64 {
    raw_score(cache, k, rank)
}

fn raw_score(cache: &Cache, k: Key, rank: Rank) -> f64 {
    match rank {
        Rank::Signal(c, a) => cache.get(k).red.signal(c, a),
        // Shallowest first. Ties (every quad at one level) fall to the same level-first
        // lexicographic tie-break, so the order is total and deterministic.
        Rank::Uniform => -(k.0 as f64),
        Rank::GreedyLookahead1 => cache.gain(k),
        Rank::GreedyLookahead1PerCost => {
            cache.gain(k) / cache.get(k).red.total_substeps.max(1) as f64
        }
        Rank::Contrast(c, a) => cache.contrast(k, c, a),
        Rank::Structured(m, c, a) => cache.get(k).red.signal_with(c, a, m),
        Rank::StructureOnly => cache.get(k).red.structure(true),
        Rank::Random(_) => f64::NAN,
    }
}

/// The error each ranking reaches at a set of budgets — the shape a table wants.
pub fn curve_at(points: &[Point], budgets: &[usize]) -> Vec<f64> {
    budgets
        .iter()
        .map(|&b| {
            points
                .iter()
                .take_while(|p| p.budget <= b)
                .last()
                .map(|p| p.error)
                .unwrap_or(f64::NAN)
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// The exact optimum, by tree DP.
// ---------------------------------------------------------------------------------------------

/// The pairwise argmins of one node's 4-way merge, kept so the optimal leaf set can be recovered.
#[derive(Clone, Debug)]
struct Back {
    /// `b01[t]` — the split count given to child 0 in `min_{a+b=t} f_c0[a] + f_c1[b]`.
    b01: Vec<u32>,
    /// `b012[t]` — the count given to the `(c0,c1)` pair in that merge against child 2.
    b012: Vec<u32>,
    /// `b0123[t]` — the count given to the `(c0,c1,c2)` triple in the merge against child 3.
    b0123: Vec<u32>,
}

/// The exact minimum `error(B)` over **all** tree-shaped leaf sets — a true ceiling.
///
/// [`Rank::GreedyLookahead1`] is greedy on immediate `Δerror` and is neither optimal nor a bound;
/// it has been measured *below the random band* on `far`. This is the bound that table was
/// pretending to have, and the assertion it supports — *no ranking may beat it at any budget* —
/// is one that can actually fail.
///
/// **Not a [`Rank`]**, deliberately. The DP is not a ranking, and putting it through
/// [`replay_with_leaves`] would re-impose a greedy order on it.
pub struct Dp {
    /// `raw[s]` — the mean per-pixel error of the best tree using **exactly** `s` splits.
    pub raw: Vec<f64>,
    /// `curve[s] = min_{s' <= s} raw[s']` — the ceiling at budget `1 + 4s`.
    ///
    /// **Prefix-minimised rather than read off `raw[s]` directly**, because more splits only help
    /// if every gain is non-negative and that is exactly the open question: a parent's `N x N`
    /// sample grid and its children's are different approximation families, so a split can make
    /// the image worse. Where the two differ, it did.
    pub curve: Vec<f64>,
    /// The `s` where `curve[s] < raw[s]` — a direct measurement of negative gain.
    pub prefix_min_binds: Vec<usize>,
    pub max_splits: usize,
    pub elapsed_s: f64,
    back: HashMap<Key, Back>,
    levels: u32,
    res: usize,
}

/// `min_{a+b=t} x[a] + y[b]`, with the chosen `a` recorded. `INFINITY` where unreachable.
fn merge(x: &[f64], y: &[f64], cap: usize) -> (Vec<f64>, Vec<u32>) {
    let mut v = vec![f64::INFINITY; cap + 1];
    let mut b = vec![0u32; cap + 1];
    for (a, &xa) in x.iter().enumerate() {
        if !xa.is_finite() || a > cap {
            continue;
        }
        for (bb, &yb) in y.iter().enumerate() {
            let t = a + bb;
            if t > cap {
                break;
            }
            let s = xa + yb;
            if s < v[t] {
                v[t] = s;
                b[t] = a as u32;
            }
        }
    }
    (v, b)
}

impl Cache {
    /// Splits available inside the subtree rooted at level `l`, capped at `max_splits`.
    ///
    /// **This cap is what makes the DP affordable**, and the naive reading hides it. A 4-way merge
    /// looks like `O(cap^4)` and is done as three successive 2-way convolutions, `O(cap^2)`; then
    /// the per-node cap is bounded by that node's own subtree, so only the top two levels ever see
    /// the full budget. At `levels = 7` the whole DP is ~120M f64 min-adds at the *complete* tree.
    fn split_cap(&self, l: u32, max_splits: usize) -> usize {
        let nodes = ((1usize << (2 * (self.levels - l + 1))) - 1) / 3;
        max_splits.min((nodes - 1) / 4)
    }

    /// Solve the tree DP up to `max_splits`.
    ///
    /// `f_k(0) = err_sum(k)`; `f_k(s) = min_{s0+s1+s2+s3 = s-1} sum_i f_ci(si)` for `s >= 1`; only
    /// `f_k(0)` exists at the deepest level. Budget and splits are locked to [`replay`]'s own
    /// accounting: `spent` starts at 1 and each split adds 4, so `B = 1 + 4s`.
    pub fn dp_optimal(&self, max_splits: usize) -> Dp {
        let t0 = std::time::Instant::now();
        let mut back: HashMap<Key, Back> = HashMap::new();
        let mut prev: HashMap<Key, Vec<f64>> = HashMap::new();

        for l in (0..=self.levels).rev() {
            let cap = self.split_cap(l, max_splits);
            let w = 1u32 << l;
            let keys: Vec<Key> =
                (0..w).flat_map(|iy| (0..w).map(move |ix| (l, ix, iy))).collect();

            let done: Vec<(Key, Vec<f64>, Option<Back>)> = keys
                .par_iter()
                .map(|&k| {
                    let mut f = vec![f64::INFINITY; cap + 1];
                    f[0] = self.get(k).err_sum;
                    if cap == 0 || l >= self.levels {
                        return (k, f, None);
                    }
                    let c = Self::children(k);
                    // Three 2-way convolutions, not one 4-way loop. Each is capped at `cap - 1`
                    // because the split of `k` itself has already been paid for.
                    let inner = cap - 1;
                    let (g01, b01) = merge(&prev[&c[0]], &prev[&c[1]], inner);
                    let (g012, b012) = merge(&g01, &prev[&c[2]], inner);
                    let (g0123, b0123) = merge(&g012, &prev[&c[3]], inner);
                    for s in 1..=cap {
                        f[s] = g0123[s - 1];
                    }
                    (k, f, Some(Back { b01, b012, b0123 }))
                })
                .collect();

            let mut cur: HashMap<Key, Vec<f64>> = HashMap::with_capacity(done.len());
            for (k, f, b) in done {
                if let Some(b) = b {
                    back.insert(k, b);
                }
                cur.insert(k, f);
            }
            // The level below is fully consumed; drop it rather than hold the whole tree.
            prev = cur;
        }

        let px = (self.res * self.res) as f64;
        let raw: Vec<f64> = prev[&(0, 0, 0)].iter().map(|v| v / px).collect();
        let mut curve = raw.clone();
        let mut binds = Vec::new();
        for s in 1..curve.len() {
            if curve[s - 1] < curve[s] {
                curve[s] = curve[s - 1];
            }
            if curve[s] < raw[s] {
                binds.push(s);
            }
        }

        Dp {
            raw,
            curve,
            prefix_min_binds: binds,
            max_splits,
            elapsed_s: t0.elapsed().as_secs_f64(),
            back,
            levels: self.levels,
            res: self.res,
        }
    }
}

impl Dp {
    /// The ceiling at budget `b`, using `replay`'s accounting `b = 1 + 4s`. Saturates at
    /// `max_splits` — a caller past that is reading a curve that was not computed, so it is
    /// clamped and [`Dp::covers`] says where the curve honestly stops.
    pub fn at_budget(&self, b: usize) -> f64 {
        let s = b.saturating_sub(1) / 4;
        self.curve[s.min(self.max_splits)]
    }

    /// Whether `at_budget(b)` is a computed value rather than the clamp. **Print this beside any
    /// truncated curve**; a curve that stops is fine, a curve silently extrapolated is not.
    pub fn covers(&self, b: usize) -> bool {
        b.saturating_sub(1) / 4 <= self.max_splits
    }

    /// The optimal leaf set at exactly `splits` splits, for the level histogram.
    ///
    /// Walks the recorded pairwise argmins back down the tree. The allocation recovered here is
    /// what says whether the optimum is uniform-depth or concentrated — which is the direct
    /// answer to why greedy prefers a flat corner.
    pub fn leaves(&self, splits: usize) -> Vec<Key> {
        let mut out = Vec::new();
        let mut stack = vec![((0u32, 0u32, 0u32), splits)];
        while let Some((k, s)) = stack.pop() {
            if s == 0 || k.0 >= self.levels {
                out.push(k);
                continue;
            }
            let b = &self.back[&k];
            let t = s - 1;
            let a012 = b.b0123[t] as usize;
            let s3 = t - a012;
            let a01 = b.b012[a012] as usize;
            let s2 = a012 - a01;
            let s0 = b.b01[a01] as usize;
            let s1 = a01 - s0;
            let c = Cache::children(k);
            for (ci, si) in c.iter().zip([s0, s1, s2, s3]) {
                stack.push((*ci, si));
            }
        }
        out
    }

    /// The optimum's **decision** for every quad it decided: `true` = split, `false` = keep.
    ///
    /// This is the ground-truth label the whole signal audit is scored against, and its
    /// population matters as much as its values. A quad that is **not in the optimal tree was
    /// never decided** -- the optimum never reached it -- so it is absent from the map rather
    /// than labelled `false`. Treating an undecided quad as a `keep` would invent a label, and it
    /// would invent tens of thousands of them: at `s = 383` the tree holds ~1500 nodes of 21845.
    ///
    /// So `internal -> true`, `leaf -> false`, `absent -> not in the map`. **Report the
    /// population size with any statistic taken over it.**
    pub fn labels(&self, splits: usize) -> HashMap<Key, bool> {
        let leaves = self.leaves(splits);
        let mut m: HashMap<Key, bool> = HashMap::with_capacity(leaves.len() * 2);
        for k in &leaves {
            m.insert(*k, false);
            // Every strict ancestor of a leaf is an internal node of the same tree. Walking up
            // from the leaves reaches exactly the internal set and nothing else, because the
            // leaves tile the root.
            let (mut l, mut ix, mut iy) = *k;
            while l > 0 {
                l -= 1;
                ix /= 2;
                iy /= 2;
                m.insert((l, ix, iy), true);
            }
        }
        m
    }

    /// Total pixels, so a caller can check a leaf set tiles the root without reaching for `Cache`.
    pub fn res(&self) -> usize {
        self.res
    }
}
