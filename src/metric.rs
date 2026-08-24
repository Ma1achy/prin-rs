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
use crate::quad::{Agg, Criterion, QuadReduction};
use crate::rng::SplitMix64;
use crate::scheduler::reduce;

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

pub struct Cache {
    pub region: String,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    pub body: usize,
    /// Deepest level present. `2^levels * n == res`, so the leaves are exactly one sample per
    /// pixel — asserted at construction, not assumed.
    pub levels: u32,
    pub n: usize,
    pub res: usize,
    pub quads: HashMap<Key, CachedQuad>,
    /// The reference image, RGB8.
    pub reference: Vec<u8>,
    pub trajectories: u64,
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
        let mut img = vec![0u8; self.res * self.res * 3];
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
) -> Cache {
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
        levels,
        n,
        res,
        quads: HashMap::new(),
        reference: vec![0u8; res * res * 3],
        trajectories: 0,
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

    let computed: Vec<(Key, QuadReduction, Vec<[u8; 3]>)> = all
        .par_iter()
        .map(|&k| {
            let (qx, qy, qh) = box_of(&c, k);
            let slice = Slice::body_plane(n, n, qx, qy, qh, body).with_chart(chart);
            let px: Vec<PixelOut> =
                (0..slice.npix()).map(|i| evaluate::<f64>(&slice, i, ens)).collect();
            let mut red = reduce(&px, n, tau, ens.t_max);
            let ics: Vec<crate::physics::Cart<f64>> =
                (0..slice.npix()).map(|i| slice.nominal::<f64>(i)).collect();
            red.n_distinct_ic = crate::decode::distinct(&ics) as u32;
            let rgb: Vec<[u8; 3]> = px.iter().map(outcome_rgb).collect();
            (k, red, rgb)
        })
        .collect();

    c.trajectories = computed.len() as u64 * (n * n) as u64 * (ens.n_extra + 1) as u64;
    for (k, red, rgb) in computed {
        c.quads.insert(k, CachedQuad { key: k, red, rgb, err_sum: 0.0 });
    }

    // ---- the reference image, from the deepest level ----
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

    c
}

/// How a replay chooses which leaf to refine next.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rank {
    /// A criterion's own signal, read as a **priority order** and never against a threshold.
    Signal(Criterion, Agg),
    /// Uniformly random. The floor: any criterion must beat this. Run several seeds and read
    /// the band, never one trace — a single random run is a draw.
    Random(u64),
    /// Greedy on immediate `Δerror`. **A strong reference, not a ceiling.**
    ///
    /// Greedy is optimal only when gains are independent and immediately available, and here
    /// they are neither: a quad whose own split gains little may unlock children with large
    /// gains two levels down, and greedy declines it. That is the classic failure of greedy on
    /// a sequential tree problem. **A criterion beating this indicates lookahead value, not an
    /// error**, and there is deliberately no assertion anywhere that it dominates.
    GreedyOracle,
    /// Greedy on `Δerror / cost`, where cost is the quad's measured substeps. §8.
    GreedyOraclePerCost,
}

impl Rank {
    pub fn name(self) -> String {
        match self {
            Rank::Signal(c, a) => format!("{}/{}", c.name(), a.name()),
            Rank::Random(s) => format!("random[{s}]"),
            Rank::GreedyOracle => "greedy_oracle".into(),
            Rank::GreedyOraclePerCost => "greedy_oracle/cost".into(),
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
    let mut leaves: Vec<Key> = vec![(0, 0, 0)];
    let mut rng = match rank {
        Rank::Random(s) => Some(SplitMix64::new(s)),
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
        let pick = match rank {
            Rank::Random(_) => {
                let r = rng.as_mut().unwrap().next_u64() as usize;
                cand[r % cand.len()]
            }
            _ => {
                let score = |i: usize| -> f64 {
                    let k = leaves[i];
                    match rank {
                        Rank::Signal(c, a) => cache.get(k).red.signal(c, a),
                        Rank::GreedyOracle => cache.gain(k),
                        Rank::GreedyOraclePerCost => {
                            let c = cache.get(k).red.total_substeps.max(1) as f64;
                            cache.gain(k) / c
                        }
                        Rank::Random(_) => unreachable!(),
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
