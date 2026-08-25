# Build brief: structure-preferring refinement, and the slippy map

**Prerequisite:** the four `principia-ii` latent presets must render first. Everything here is
judged by eye against slices the user recognises, and until those exist nothing below is
evaluable — which is the state the vertical slice is currently stuck in.

**Deliverable:** a refinement criterion that prefers **spatial structure** over **uncertainty**, a
rank-based priority replacing threshold-based splitting, and the camera behaviour a slippy map
needs.

---

---

> **Four corrections to this brief, found by reading the repo before building.** §2's
> "no connected-components code exists", §2.1's neighbour contrast, §4.4's "no neighbour lookup"
> and §6's "`error(B)` may not be built yet" all describe code that **already exists and is
> tested**. Each is marked inline below. The genuinely missing pieces are the hot **rule**,
> `grad_rms`, rank-based priority *inside the live descent*, the two modes, the 2:1 balance pass,
> camera-biased priority and the persistent frontier.
>
> **And §5's headline acceptance test cannot fail.** *"Assert `n_hot_within < N^2` for a stated
> majority"* passes trivially and unconditionally under any quantile rule, because `n_hot` is set
> by the rule and not by the field. See §2.1's note below and `RESULTS.md` §14.5-14.6.

## 0. MEASURED: the criterion is currently always-split, and the cause is arithmetic

**Read this before §1.** All 60 committed `.prnq` tree dumps in `results/` were parsed. The
refinement failure has a single root cause and it is not conceptual.

### 0.1 `tau_display` sits at the 0.4th percentile of the spread distribution

`tau_display = 1e-4` in every committed run. The **median quad spread across all 17 charts is
6.24e-4** — so `tau` is nearly four times *below* the median, and:

| `tau` | % of quads exceeding it |
|---|---|
| **1e-4 (current)** | **99.6%** |
| 1e-3 | 22.4% |
| 3e-3 | 2.1% |
| 1e-2 | 1.0% |
| 3e-2 | 0.3% |

**The split predicate is true for essentially every quad.** The criterion is not selecting; it is
an always-split rule wearing a threshold's clothing.

### 0.2 Consequence: the trees are uniform

On **13 of 17 charts, 99.4–100% of leaves sit at max depth.** Three (`latent_shape`,
`latent_mass`, `mass_simplex`) are *exactly* 4096 = 4^6 — completely uniform to the screen floor.
Depth variance runs 0.002–0.031. Only `body_plane`/`plane_00deg` (61% at max, var 1.015),
`preset_shape` (78%, 0.691) and `shape_sphere` (79%, 0.448) show any selectivity at all.

**Balanced mode does not currently exist.** It degenerates to uniform on almost every chart, which
is the §3.2 failure already occurring rather than a risk to guard against.

### 0.3 Consequence: the structure fields are saturated and carry no information

`n_hot_within = 64` — **every** footprint hot — in 99–100% of quads on 15 of 17 charts, and
`n_components_within = 1` universally. The hot mask is a single blob covering the whole quad,
everywhere.

**So §2's spatial measures cannot discriminate as currently thresholded.** This is the same
instrument-cannot-see-its-own-signal failure as the linear-ramp spread image, and it is why §2.1's
**relative** threshold is load-bearing evidence rather than a precaution.

### 0.4 The termination mechanism, confirmed directly

In `preset_shape`, comparing shallow leaves against max-depth leaves:

| | shallow (129) | at max depth (448) |
|---|---|---|
| `escape_fraction` | **0.094** | **1.000** |
| `spread_median` | 1.0e-4 | 1.3e-2 |

**Refinement is following escape, not structure.** Non-escaping regions have collapsed spread and
are abandoned; fully-escaped regions keep high spread and get refined to the floor. Terminal states
are absorbing, so nearby copies share one, `spread_event` collapses — and the criterion reads that
as "resolved".

> **CORRECTED at the corrected window.** `preset_shape` re-rendered at `half = 3.0` is **16
> leaves**, not 577, and its leaf-spread median is `2.86e-1` — **3400x above `tau`**. It clears the
> spread gate on *every* leaf and is stopped by `alpha` (8 `floor` + 8 `keep`, zero spread-gate
> failures). It is the only tree in the corpus exercising the alpha gate, so the shallow/max-depth
> comparison above is not available on it any more. Read the decision column; a mechanism read off
> a tree shape is a guess.

`preset_shape_pl` shows the same thing in the wireframe: 2974 leaves, of which **14 are shallow**
(12 `Floor`, 2 `Keep`). Those fourteen are the coarse blocks sitting over the structured centre.
Everything else went to max depth. **Both observations are true at once: it refines nearly
everything, and the handful it declines are the interesting parts.**

### 0.5 What this fixes in the plan below

1. **§2.1's relative threshold is mandatory**, and §0.3 is the evidence. An absolute cutoff is what
   saturated the mask.
2. **Sweep `tau` over 1e-4 … 1e-1**, not around 1e-4. All the interesting behaviour is above 1e-3;
   below it the predicate is constant.
3. **§3's rank-over-threshold argument is now decisive rather than preferred.** A threshold at the
   0.4th percentile is exactly the failure rank-based priority is immune to — a ranking does not
   care where the absolute level sits.
4. **`n_components` cannot help until the mask desaturates.** Fix the threshold first, then measure
   the spatial fields; measuring them now would report a null that was guaranteed.

### 0.6 One caveat on this analysis

These trees were rendered at `half = 1.0` where the reference uses **3.0** — a 3x crop. §0.1's
threshold arithmetic and §0.3's saturation are **framing-independent** (they are statements about
the spread distribution, not about which patch is shown). But **§0.4's spatial claim — that the
coarse quads sit over the structured centre — must be re-taken at the corrected window** before it
is acted on. The mechanism is physical and should survive; that is a prediction, not a result.

**Test it directly and cheaply:** per quad, scatter **leaf depth against `terminated_fraction`**.
The prediction is a clear anti-correlation. One plot, no extra compute, and it tests the *cause*
rather than the appearance.

---

## 1. The diagnosis: the criterion measures the wrong thing

`ensemble_spread` measures **uncertainty** — how much a footprint's copies disagree. What is wanted
is **spatial structure** — how much the field varies *across* the quad. These come apart exactly
where the criterion fails:

| region | uncertainty | spatial frequency | refining reveals |
|---|---|---|---|
| uniformly chaotic sea | **high** | **low** | more mush |
| filament / basin boundary | **high** | **high** | structure |
| smooth region | low | low | nothing |

The current criterion cannot separate the first two, so it spends budget on the sea. **That is the
`deep interior` failure** — 29 quads against `near-field`'s 4617, in the richer region, with the
large high-spread structures left at level 2.

Note this is a *better* diagnosis than the earlier median-blindness attribution: changing the
aggregate from median to p90 changes *which footprints count*, but **neither measures variation
across the quad at all.**

### 1.1 It reconciles shader-agnostic with visually good

The criterion must not depend on the colour map — changing a shader must not force recomputation.
But it must prefer refining things that *look* interesting.

**These reconcile, because any reasonable colouring of a high-spatial-frequency physics field is
itself high-frequency.** Refine on structure in the physics and the image follows.

**But agnostic to the colour map is not agnostic to the field.** A filament in `ensemble_spread`
and a filament in `FTLE` are not the same set. So the criterion reads a **fixed named field**, not
whatever the active graph happens to display — the same argument that made `ensemble_spread` a
single named scalar rather than a per-graph quantity.

**The field is `ensemble_spread`.** It is the only quantity in the system with a settled definition
(`max(spread_shape, spread_event)`, both bounded by their achievable maxima) and it is already
computed per footprint. **Do not invent a set.** If a second field is ever added it is a spec change
with its own justification, not an implementation choice.

### 1.2 Structure is only meaningful relative to pixel size

A filament narrower than one pixel is not structure — it is noise, and refining to resolve it
produces detail nobody can see. So the structure measure is **view-relative**, exactly like the
screen floor, and for the same reason.

**Consequence: zoom invalidates the criterion, not the data.** Zooming changes no quad's physics
but does change whether its structure is visible. A quad that floored at one zoom **must
re-evaluate** at the next. Never cache a structure verdict as a quad fact.

---

## 2. The measure

`ensemble_spread` is already computed per footprint, `N x N` per quad. **The spatial layout of
those values is free information currently discarded.** Compute, on that field:

| field | definition |
|---|---|
| `n_hot` | count of footprints above the quad's **own median** — see below |
| `n_components` | connected components of the hot mask, 4-connectivity |
| `largest_component` | size of the biggest |
| `perimeter_ratio` | perimeter / area of the hot mask — **thin and connected -> filament** |
| `grad_rms` | RMS of the finite-difference gradient across the footprint grid |

> **CORRECTED, and this is the one that would have cost the most.** Under **any** quantile rule
> `n_hot` is fixed by the rule and not by the field — 31 of 64 at `N = 8, q = 0.5`. So `frac_hot`
> carries essentially no information once the mask is relative, and
> **`frac_hot_between/median` is the best criterion measured on this project**. Replacing the
> absolute mask, as written below, would have silently deleted the best-performing signal in the
> system and read as an improvement.
>
> **Both masks are now computed and dumped.** `spatial::HotRule` selects which one the *shape*
> criteria read; `frac_above_tau_*` is untouched. Two further measurements: the saturation is
> **regional** — in `far` the absolute mask is *empty*, not full, and every mask-derived criterion
> takes one distinct value there — and desaturating **coarsens** the ordering, `LayoutRel` running
> 78 -> 26 -> 17 -> 9 distinct values across `abs, q[0.50], q[0.75], q[0.90]` against `Layout`'s
> steady 58.
>
> **The hot threshold must be RELATIVE, not absolute — and this is measured, not argued.** With the
> current absolute cutoff, `n_hot_within = 64` (every footprint hot) in 99–100% of quads and
> `n_components = 1` universally (§0.3). The mask is saturated and the spatial fields carry no
> information at all. An absolute cutoff also reintroduces exactly the tunable constant §3 exists to
> remove, and it would be chosen by eye. Use the quad's **own median**
> (or a fixed quantile of its own distribution). That makes `n_hot` and its companions **shape**
> statistics rather than **magnitude** ones — which is what §2 is for, and it means they do not drift
> as the signal rises globally with `t`. Report the quantile used; sweep it if it looks load-bearing.

**Thin and connected -> boundary -> split. Scattered or filling -> chaos -> floor.** That answers
the split/floor question directly from one quad, without needing a second level — which matters
given the vertical slice found only four discretionary levels beneath the screen floor.

> **CORRECTED — this was written without checking the repo.** Connected components **already
> exist**: `src/spatial.rs:59` `layout()` is a 4-connectivity flood fill returning `n_hot`,
> `n_components`, `largest_component` and `perimeter_ratio`, with the internal-edges perimeter
> convention documented in the module header. The hand-built-mask tests asked for below are
> already in `tests/criterion.rs` — single blob, two blobs, a one-cell diagonal filament, full,
> empty, and the checkerboard. All four fields are in the `.prnq` dump and always were.
>
> **`grad_rms` was the only one genuinely missing**, and it is now `src/spatial.rs::grad_rms`.
> **Only the hot RULE was wrong** — an absolute cut at `tau`, which is what saturated the mask.

### 2.1 Neighbour contrast is REQUIRED, not optional — the edge-filament blind spot

A filament running along a quad **boundary** shows as low internal variation in *both* quads. Every
within-quad structure measure is blind to it — and boundaries are precisely what the criterion
exists to find. **The blind spot is systematically aligned with the target.**

Neighbour contrast fixes it: `contrast = max over the 4 neighbours of |signal_self −
signal_neighbour|`.

That makes the neighbour lookup do **three independent jobs**: the 2:1 balance constraint (§4.4),
steady-state stability under a rising signal (§3), and edge-filament detection here. **Three
independent arguments for one mechanism is a strong sign it is the right mechanism.** Build it
once, use it three ways.

> **CORRECTED — it is already built, twice over.** `QuadTree::contrast(i, criterion, agg)`
> (`src/quad.rs:592`) is exactly `max` over the four edges of `|signal_self - signal_neighbour|`,
> and it returns the **edge count** beside it so a root-border quad's low bias is visible rather
> than absorbed. It is already scored in the `error(B)` harness as `metric::Rank::Contrast`.
>
> **The gap is that `scheduler::decide` does not read it.** That is a wiring change, not a build.

### 2.2 The open question, and the recommendation

Does structure **replace** `ensemble_spread` in the decision, or **multiply** it?

- **Replace** says structure is the whole answer.
- **Multiply** says *"uncertain AND structured"*, which floors the uniform sea by construction.

**Implement both as selectable variants and measure.** The recommendation is multiply — it keeps
the determinacy question that `ensemble_spread` answers while adding the structure question it
cannot — but this is exactly the sort of thing that has been settled by measurement rather than
argument throughout this project.

---

## 3. Rank, not threshold — and why this is forced

**This is not a prediction — it has already happened.** `tau_display = 1e-4` sits at the 0.4th
percentile of the observed spread distribution, so 99.6% of quads exceed it and 13 of 17 charts
refined to 99.4–100% max depth (§0.1–0.2).

Any criterion comparing a signal against a **fixed threshold** must degenerate to uniform depth,
because spread grows with `t` everywhere. That is the treadmill; it is not avoidable by choosing a
better `tau`. **Only a comparison invariant to a global rise can hold a steady state.**

Two are:

- **Neighbour contrast** — if everything rises together, contrast does not move.
- **Rank within the current frontier** — always refine the top `k`% by structure, so the tree
  *redistributes* rather than deepening.

**Rank is the recommendation**, for three reasons: it makes the **budget** the constraint rather
than a threshold, which is what a frame-budgeted slippy map needs anyway; it removes a tunable
constant, which has been the recurring defect throughout; and it makes the two modes one mechanism.

### 3.1 The two modes are one mechanism with different budgets

- **Uniform mode** — refine everything to the screen floor. The criterion is **off**.
- **Balanced mode** — the same descent, but the frontier is **ranked** and only the top `k` gets
  budget per frame. The criterion is a **priority ordering**.

`k` is the balanced-mode budget and the whole mode turns on it. **It has no recommended value —
sweep it**, and note that `tau` (if retained at all) must be swept over **1e-4 … 1e-1**: below 1e-3
the predicate is effectively constant (§0.1) and report how depth variance and churn (§3.2) move with it. Do not pick one because a
tree looked right; that is the constant-tuning defect in its most tempting form.

This is the reframe that matters: **the screen floor stops things; the criterion decides what gets
attention first.** The criterion was never a stop condition.

**A quad can be demoted.** Its neighbours grow more structured, it falls down the ranking, it stops
being refined further. That is the tree "updating over time" — and it needs **no merging and no
eviction**, only not spending on it.

### 3.2 The acceptance test for balanced mode — it must be able to fail

Balanced mode is a **steady state, not a frozen state**: the tree keeps updating; it must not drift
toward uniform depth as `t` advances.

```
march the playhead forward; plot DEPTH VARIANCE across leaves against t
  degenerating : variance -> 0   (everything converges to max depth)
  balanced     : variance roughly constant, while individual quads still move
```

**Report both**: the variance curve *and* the per-quad churn (how many leaves change decision per
step). A tree that is stable because nothing moves is not balanced, it is frozen — and the two look
identical in a variance plot alone.

---

## 4. The slippy map

### 4.1 Breadth-first, and it is forced

Depth-first drives one region deep while the rest stays coarse, which is wrong for progressive
display: the user sees one sharp patch in a blurry field. **Breadth-first improves everything
visible together.** Quads are disjoint, so it parallelises naturally — refine the whole frontier
slice concurrently.

### 4.2 Frame budget, not total budget

Every scheduler experiment so far spends a fixed quad count **once**. A slippy map has ~16 ms per
frame, **forever**. That changes the question from *"which quads"* to *"which quads this frame"*,
and is what makes breadth-first a requirement rather than a preference.

Report **quads per frame achieved** against the budget, and what happens when the budget is missed.

### 4.3 Camera-biased priority

*"Bias refining to where the camera is going, then apply the standard criterion. Do not lazily
refine from the middle outwards — that is incorrect."*

So priority is a **product of two terms**: camera relevance and structure. Neither alone.

**The destination needs defining, and the obvious choice is wrong.** Velocity extrapolation fails
on flick-and-stop — it prefetches past where the user lands. Three options; **implement one, state
which, and say why**:

1. extrapolate velocity forward by a fixed horizon
2. refine the whole **swept path** between current and predicted position
3. simply widen the viewport margin and drop prediction entirely

Option 3 is the honest baseline and the others must beat it to justify their complexity.

### 4.4 The 2:1 balance constraint — needed regardless

No two adjacent leaves may differ by more than one level, or the adaptive render has cracks. This
is **separate from** the neighbour-contrast idea and is required for rendering.

**It forces splits the criterion did not ask for**, which interacts with the budget: report what
fraction of splits are balance-forced rather than criterion-driven. If it is large, the budget is
being spent on geometry rather than physics.

> **CORRECTED.** `QuadTree::neighbour(i, dir)` **already exists** (`src/quad.rs:556`), built
> exactly this way — root descent toward a probe point just outside the edge, stopping at the
> deepest computed node, so it returns the same-or-coarser neighbour. `O(level)`, nothing cached,
> computed at decision time. `tests/criterion.rs:439` already checks it against an independent
> geometric box-touch predicate over a whole tree.
>
> **What is missing is the 2:1 balance pass itself**, and the `Decision` variant that makes
> balance-forced splits countable in the dump rather than indistinguishable from criterion-driven
> ones.

---

## 4.5 What the slippy map does to the tree

**The playhead does NOT freeze during interaction.** A frozen sim during navigation means the thing
being explored stops being alive exactly when it is being explored. **Consequence: every measurement
taken during a gesture has two causes** — view change and time change. Log the camera delta and the
playhead delta per frame so a churn spike can be attributed to one or the other.

**Pan and zoom invalidate different things.**

| gesture | physics | structure verdict | tree |
|---|---|---|---|
| **pan** | unchanged | unchanged | valid; just need *more* of it at the same level |
| **zoom** | unchanged | **stale for every visible quad at once** | valid, but the whole frontier's ranking must be recomputed |

Because structure is pixel-relative (§1.2), a zoom step invalidates the entire frontier's priority
while invalidating no physics. **Measure the two gestures separately** — how much of a frame's
budget goes to re-ranking after a zoom versus after a pan. If they differ markedly the two gestures
will *feel* different, and that is a design fact to state rather than a bug to hide.

**Zoom-out should be nearly free, and there is a test for it.** Zooming *in* reveals deeper quads
that may not exist. Zooming *out* reveals shallower ones whose parents are already in the tree with
valid reductions. **Assert: the count of newly-computed quads after a zoom-out is ≈ 0.** If it is
not, something is recomputing what the tree already holds.

**The tree lags the camera, and the gap is filled by drawing the coarse ancestor.** At ~16 ms per
frame a newly-revealed region cannot be refined before it is drawn. Three options were considered;
**take the first**:

1. **draw the coarse ancestor and let it sharpen** — honest, visibly blocky during motion
2. hold the previous frame until ready — shows *stale* data as current
3. draw coarse with a "still resolving" cue

Option 1 is the only one that never lies, and it is what *"never mixed-time, never frozen, never
lying"* already commits to. **Big texels during motion are a deliberate choice**, and if it looks
bad the remedy is more budget or better ranking, never concealment.

**Rank on the VISIBLE part of a quad, not the whole quad.** A quad half off-screen has structure the
user cannot see, and ranking on it spends budget on nothing. This makes the ranking shift as the
camera pans without any physics changing — which is the *same* class of view-dependence as §1.2, not
a new one.

### 4.6 The persistent frontier

Each frame refines the top `k` leaves by priority, so the frontier must be ordered. Rebuilding that
order from scratch every frame is O(n log n) over thousands of leaves, 60 times a second.

Normally that would be a rounding error against 512 trajectories per quad — **but camera bias
changes the arithmetic.** Priority is structure × camera relevance, and **camera relevance changes
for every quad on every frame the camera moves.** So during a gesture the naive version is not
re-sorting a mostly-unchanged list; it is genuinely recomputing all of it, every frame. **Build the
persistent frontier.**

Three things to get right:

**Split the priority into a stored term and a derived term.** Structure changes only when a quad is
recomputed or the zoom changes; camera relevance changes every frame of motion. **Store structure on
the quad; compute camera relevance at query time.** A pan then touches a distance calculation and
never the physics term — and it keeps camera state off the quad, which the *"never cache view state
as a quad fact"* rule already requires.

**A plain binary heap will not do the update needed.** Reprioritising an entry already inside needs
either a `quad_id → heap_position` map so a changed entry can be sifted, or **priority bucketing**,
re-bucketing only when an entry crosses a band boundary. Bucketing suits a rank-based scheme, which
needs the top slice rather than a total order.

**The failure mode is staleness, not slowness — and it is invisible.** An incrementally-maintained
frontier that is *wrong* looks exactly like a criterion that is wrong: a quad sitting high in the
queue on a priority it no longer has. **Test: rebuild from scratch every N frames and assert the
incremental frontier matches.** Not a benchmark — a correctness check with teeth, the same shape as
the `Gamma`-identity chain: an independent path to the same answer, so a silent divergence cannot
survive.

**Keep the from-scratch path permanently** as the reference implementation. Do not delete it once
the fast one works.

---

## 5. Tests that can fail

- connected components on hand-built masks: one blob, two blobs, a one-cell diagonal filament,
  full, empty
- `perimeter_ratio` is **high for a filament and low for a blob** of equal area — asserted as an
  inequality, not eyeballed
- **balanced mode does not degenerate**: depth variance stays bounded away from zero over a march,
  *and* per-quad churn stays non-zero (§3.2)
- **uniform mode does degenerate**, by construction — the control that proves the test can tell
  them apart
- **the hot mask is NOT saturated**: assert `n_hot_within < N^2` for a stated majority of quads.
  At present it is exactly `N^2` in 99–100% of them (§0.3), so this test fails on today's code —
  which is what makes it worth having
- **depth anti-correlates with `terminated_fraction`** on the current criterion (§0.6), and the
  correlation weakens or reverses under the structure criterion. The scatter is the evidence that
  the mechanism was correctly identified
- a quad's structure verdict **changes with zoom** at fixed physics (§1.2) — the view-relative
  property, and it fails if the verdict is cached
- `neighbour(i, dir)` agrees with the geometric predicate over a whole tree
- the 2:1 constraint holds over every produced tree
- structure is **not** a `Quad` field — a compile-level guarantee, as the screen floor's veto has
- **neighbour contrast detects an edge filament** that within-quad structure misses: a synthetic
  field with a filament laid exactly along a quad boundary must be flagged by contrast and missed by
  every within-quad measure — the test that proves §2.1's blind spot is real and closed
- **zoom-out recomputes ≈ 0 quads** (§4.5)
- **the incremental frontier matches a from-scratch rebuild** every N frames (§4.6) — the staleness
  check, and the from-scratch path is retained as the reference
- camera relevance is computed at query time, never stored on a `Quad` (§4.6)
- the four `principia-ii` presets render, and `z = 0` gives the equilateral Lagrange configuration:
  masses `(1/3, 1/3, 1/3)`, `I = 1.0`, `sum p = 0`, positions
  `(-0.866, -0.5), (0.866, -0.5), (0, 1)`

---

## 6. Cautions

**Judge on the recognisable slices.** The four presets are the point — they are the only images
whose correctness can be assessed by eye. A criterion that improves an unrecognisable slice has
improved nothing checkable.

**Validate on `deep interior`.** A change that only improves `near-field` is tuning.

**Do not tune to a nice-looking tree.** The temptation is strongest here, because the render finally
looks like the product and the criterion is explicitly about looking good. **The depth-variance test
and — if it exists — the `error(B)` curve are the results; the picture is a diagnostic.**

> **CORRECTED — `error(B)` is built and mature.** `src/metric.rs` carries `Cache`, `err_sum`,
> `gain`, `replay`, `curve_at` and a `Rank` enum with random *bands* (never one trace) and both
> greedy oracles; `examples/criterion_metric.rs` runs twenty rankings and writes the curves. So
> §2.2's replace-vs-multiply comparison **does** have a quantitative measure and no reordering is
> needed.
>
> Two properties of it that constrain how this brief's changes are scored: a criterion enters the
> replay **as an ordering and never against `tau`** (`tests/criterion.rs:408` enforces it), and
> `error = 0` means *matches this sampling*, not *correct*.

**Report negative and messy results.** Every PR in this project has corrected something stated with
more confidence than it deserved. If structure-preference does not fix `deep interior`, that is the
headline.

---

## 6.1 Ordering

Each stage depends on being able to judge the previous one by eye. **Land them in this order:**

1. **The four `principia-ii` presets.** Nothing below is evaluable without slices the user
   recognises — that is the state the vertical slice is currently stuck in.
2. **Structure measures** (§2), both variants, dumped and rendered. Judge on the presets.
3. **Rank-based priority and the two modes** (§3), with the depth-variance test.
4. **The slippy map** (§4.4–§4.6) — neighbour lookup, 2:1 balance, camera bias, persistent frontier.

Open a PR per stage. Stage 1 is small and unblocks everything; do not bundle it.

## 7. Definition of done

- Structure fields computed per quad, both replace and multiply variants selectable
- Rank-based priority, with the two modes as one mechanism at different budgets
- Balanced mode passes §3.2's depth-variance test, and uniform mode fails it as a control
- Camera-biased priority with the destination model stated and beaten against the margin baseline
- Persistent frontier with stored-structure / derived-camera split, and the from-scratch reference
  retained and cross-checked
- Coarse-ancestor fill during motion; ranking on the visible part of a quad
- Pan and zoom re-ranking costs measured separately; zoom-out asserted ≈ free
- 2:1 balance enforced, with the balance-forced split fraction reported
- All four presets rendered and committed
- Still no eviction, no async, no promotion — **if one appeared, that is a bug**
