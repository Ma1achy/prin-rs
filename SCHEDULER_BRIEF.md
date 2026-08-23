# Build brief: minimal refinement scheduler (prin-rs)

**Audience:** you have built the uniform kernel and know the codebase. This adds the one thing it was
deliberately built without — **the quadtree and the descent loop.**

**Deliverable:** a binary that starts from a coarse quad, refines adaptively using the criterion, and
dumps the resulting tree with every decision it made. It exists to answer questions that only appear
when the criterion runs **in a loop**, which nothing so far has tested.

**Scope discipline.** No eviction, no caching, no async, no promotion, no GUI, no camera, no
interaction. Batch job: descend, dump tree, exit. Everything cut is deliberate.

---

## 1. Why this exists

Every measurement so far is **one split, measured in isolation**. The criterion has never run
iteratively, and every remaining question is dynamic:

| question | why it matters |
|---|---|
| **Does the descent terminate?** | With Wada-like boundaries dense at every scale, spread may stay high however far you refine. This was flagged at the outset and never tested. |
| **Does the floor actually engage?** | Or does it refine chaotic regions indefinitely until a cap stops it? |
| **Under a budget, does it spend well?** | The real scheduler question. Nothing touches it. |
| **Does the tree look right?** | Dense at boundaries, sparse in the sea. Visible, checkable. |
| **Does per-quad noise cause thrash?** | Chaotic-region α scatters 1.1–1.3, so neighbouring quads get different decisions from the same underlying physics. |

---

## 2. The criterion, stated whole

### 2.1 The three levels, and what is measured at each

| level | what it is |
|---|---|
| **QUAD** | a square patch of IC space. The unit of scheduling and the thing that splits. Contains **N × N footprints**, `N = SAMPLES_PER_QUAD_AXIS`, 8–32 across quality tiers; **8 here**. |
| **FOOTPRINT** | one pixel = one nominal initial condition. Carries `E+1` **copies**. This is the unit `ensemble_spread` is defined on. |
| **COPY** | one simulation. Offsets are the fixed Halton (2,3) prefix indexed by `copy_index`, scaled to the **cell width** — the footprint spacing — **not** the quad width. |

At `N = 8`, `E+1 = 8`: **one quad = 64 footprints = 512 trajectories.** That is why §3.6 counts the
budget in quads and why the descent needs a cap.

Copy offsets are fixed, not per-pixel: copy *k* sits at the same relative offset in every footprint
at every level. This is what makes parent and child comparable — **common random numbers by
construction**.

**`ensemble_spread` is a per-footprint quantity.** A quad holds `N²` of them, and §3.4 is how they
become one number.

```
ensemble_spread = max( spread_shape , spread_event )      -- per footprint
```

**Why one footprint per quad would break the criterion**, not merely the cost model — stated so it
is not reintroduced later as an optimisation:

- There would be **no between-footprint variation**. `ensemble_spread` would be pure within-pixel
  disagreement, and the spatial structure the criterion exists to detect would be invisible. A quad
  has to sample an **area**.
- Splitting would mean "the same single point at a finer cell width", which measures nothing about
  whether the region is spatially resolvable.
- §3.4 would collapse — the mean/median/p90 question exists *because* a quad has many footprints,
  and with excess kurtosis 110 that choice materially changes decisions.
- It silently reintroduces the single-scale cancellation of §2.3: the only two scales left would be
  within-footprint, which measures the SSAA jitter fraction rather than the dynamics.

**`E+1` and `N` do different jobs and must not be traded off against each other casually.** `E+1`
controls how well a *footprint* knows its own value; `N` controls how well a *quad* knows its own
*area*. The criterion needs both, and **only `N` affects the between-footprint signal that drives
splitting**. Too low an `N` and a quad misclassifies itself as coherent by undersampling its own
area.

Both bounded by **achievable maxima** — chord ≤ 2 on the unit sphere; `1 − 1/(E+1)` for a
disagreement fraction. The normalisations are therefore *facts, not choices*, and no free constant
decides the answer. (This is why kinetic energy was dropped despite winning most footprints: it is
unbounded and no domain worked in more than one configuration.)

**Nothing is ever discarded.** Every footprint always carries `E+1` copies; a badly-integrated
trajectory is a *measurement outcome*, not missing data. Trust is measured alongside, never enforced:
`error_ratio` (which sees spread, so correlated drift is invisible to it) paired with
`worst_energy_drift`.

### 2.2 The exponent is a measurement, not a forecast

When a quad splits, the exponent falls out of two numbers the tree already holds:

```
alpha = log( spread_parent / spread_child ) / log 2
```

`alpha ~ 1` → refining halved the uncertainty, **continue**.
`alpha ~ 0` → refining changed nothing, **floor**.
spread ~ 0 at both scales → **keep coarse**.

**This is available only *after* splitting**, and that is the point. α forecasts realised gain poorly
(median log error 0.219), but prediction is unnecessary when you are measuring what the split bought.
**Refinement is descent with feedback**: split, observe, decide whether to continue.

**Use α ordinally** — it ranks, it does not forecast `2^-alpha`.

### 2.3 Two things that are settled and must not be re-derived

**Never a single-scale ratio.** Within one quad, sample spacing and jitter both amplify by the same
`e^(lambda t)`, so any within-quad ratio cancels to the SSAA jitter fraction — it measures the
anti-aliasing knob, not the dynamics. **Genuinely different cell widths are required**, which means
genuinely different levels.

**Never pool children to synthesise a parent.** With fixed offsets a pooled 2×2 block is four exact
repeats of one pattern at four cell centres, *not* a wider-footprint ensemble. A true parent carries
offsets scaled to **its** width. Measured surrogate error: **+38.6%**, flat in E. The uniform kernel
had to pool because it has no tree; **this build must not.**

### 2.4 The measured limits

| region type | α interdecile | reading |
|---|---|---|
| tame (mid-field, far) | **0.0004–0.001** | resolves per quad, trivially |
| chaotic (near-field, body2 core) | **1.1–1.3** | does not resolve per quad, and no ensemble scheme fixes it |

Separation between regions is **0.9862**, measured the same way as the scatter.

**93% of chaotic scatter is chaotic divergence, not sampling noise** — the Halton switch cut the
control's floor 480× and moved α's scatter by 0.2%. More copies do not help.

**"Not resolvable per quad" is the correct answer for a chaotic quad**, not a defect: the criterion
is asked "would splitting help?" and returns "no, and no measurement will tell you otherwise". The
noise sits in the branch where the decision is *floor* either way, and the failure direction is
conservative — contamination inflates spread, pushing toward *refine*, wasting budget rather than
losing structure.

**The distribution is heavily tailed** (excess kurtosis 110; interdecile/sd 0.866 against a normal
2.563). A scheduler decides per *typical* quad, so **the interdecile is the measure** — never quote a
variance reduction as the improvement.

### 2.5 Two shared-footprint effects the grid produces, measured

`Slice::axis` is endpoint-inclusive, and a child's sample grid is a **strict refinement** of its
parent's with a shared origin. Two consequences, both verified against the actual formula rather
than reasoned about:

**Siblings share an edge.** Adjacent children are at the same level, so identical cell width,
identical Halton offsets, **identical copies** — the shared column is duplicated work, not a
differing ensemble. One column of `N` footprints, so `1/N` of a quad's data: **25% at N=4, 12.5% at
N=8, 6.25% at N=16.**

This is a non-physical correlation between neighbours, and **§4 question 4 is exactly what it
corrupts**: shared footprints make adjacent quads more alike than the physics makes them, so thrash
is **under-reported**. Not fatal, but the thrash figure must be quoted with the overlap fraction
beside it — and the under-reporting is worst at small `N`, which is where the `N` sweep starts.

**Parent and child share footprints too, and this one is beneficial.** Every parent sample inside a
child's extent is *also* a child sample — structurally, for every `N ≥ 2`, because the grids share
an origin and the spacing halves exactly. Copy 0 is un-jittered, so at those footprints the parent
and child nominal trajectories are **identical**: common random numbers between the two scales
`alpha` is a ratio of, which should reduce noise in that ratio.

It is **not** parity-dependent. Verified overlap fractions: exactly **25.00%** at every even `N`
(4, 6, 8, 16), and **higher** at odd `N` — 36.00% at 5, **32.65% at 7**, 30.86% at 9, 28.03% at 17.
Odd `N` does not lose the CRN, it strengthens it.

That still makes an odd `N` worth including in the sweep, for the opposite reason to the obvious
one: it is the only available lever that **varies** CRN strength, which is what separates "coarse
`N` under-splits" from "parent–child CRN is doing the work". Sweep `N ∈ {4, 7, 8, 16}`.

### 2.6 The precision floor is real and must be detected, not hit

At level `l`, `half = half0 / 2^l`, so `cell = 2·half/(N-1)`. At `half0 = 0.05`, `N = 8` the cell
width crosses `f64::EPSILON` at **level 45.87** — below which the copies stop being distinct initial
conditions and the spread is pure noise.

Guard at `1e3 × eps`, which triggers at **level 35.90**: comfortably above any physically meaningful
descent, and principled rather than an arbitrary cap. **A descent that reaches level ~36 has hit a
numerical floor, not a physical one, and must not be reported as "did not terminate".**

---

## 3. What the scheduler must do

### 3.1 The descent loop

```
start:  one quad covering the slice, at level 0
loop:   integrate any un-computed quads to the playhead
        reduce each to its QuadReduction
        for each leaf with a parent: alpha = log(sp_parent/sp_child)/log 2
        decide split / floor / keep per §3.2
        split the chosen leaves
until:  no leaf wants splitting, or budget exhausted, or a floor is hit everywhere
```

**Level 0 has no parent, so no α exists.** Split unconditionally for the first level or two to
bootstrap — state which, and make it a parameter, since it is a policy choice not a physical one.

### 3.2 The decision

From the scheduler contract, adapted to what is actually implemented here:

- **Guards first.** Terminal or budget-exhausted → stop. **Default is keep.**
- **Split** if `ensemble_spread > tau_display` **and** α is high **and** `level < max_level`.
- **Floor** if α ≈ 0 with spread still high — splitting will not reduce it. **Flag it and do not
  re-queue.** This is the branch that must be shown to engage.
- **Keep coarse** if spread ≈ 0 at both scales.

`tau_display` and the α threshold are parameters. **Sweep them; do not tune them to a nice picture.**

### 3.3 The reliability signal — the one genuine addition to try

Separation in α's **value** is 0.9862 against a chaotic scatter of 1.1–1.3 — marginal. Separation in
α's **reliability** is 0.001 against 1.2 — **three orders of magnitude**. That signal is currently
discarded, and it is free:

**A split produces four children, therefore four α values.** Their spread is a per-quad reliability
estimate at no extra cost.

- four α tightly clustered → smooth quad, α trustworthy, **act on its value**
- four α scattered over ~1 → chaotic quad → **floor**, without needing a reliable α at all

The unreliability *is* the answer, which removes the awkwardness of thresholding a quantity whose
noise is comparable to its range.

**Implement it, dump `alpha_sibling_spread` per quad, and compare the two policies**: threshold on α
alone versus threshold on sibling spread. This is the most promising available improvement and it has
never been tested.

### 3.4 Aggregation within a quad

A quad holds `N²` footprints, each with its own `ensemble_spread` (§2.1). `QuadReduction` needs
**one** number per quad, so those `N²` values must be aggregated. With excess kurtosis 110, **a mean
is dominated by a single footprint.**

Check what the design docs specify; if unspecified, dump **mean, median and p90** and report which
the decisions are sensitive to. **Do not silently pick one** — and report the **decision-level**
disagreement (how many quads decide differently under each), not merely the three numbers.

### 3.5 Playhead

All quads at one playhead. When a quad splits, children are integrated to that same `t` before any
comparison. State is a pure function of `(IC, sim key, t)`, so reaching `t` is path-independent and
the comparison cannot be contaminated by a time difference.

**Do not implement:** live promotion, catch-up-while-marching, mixed-time display. Those are
production concerns; here everything is at one `t` for the whole run.

### 3.6 Budget

A cap on **total quads computed** — quads, not trajectories, not footprints. The real cost is
`N² × (E+1)` trajectories per quad, so at `N = 8`, `E+1 = 8` a budget of `B` quads is `512·B`
trajectories; at the measured ~0.73 ms per footprint that is **~47 ms per quad**, and a 50 000-quad
descent is ~39 minutes. State the cost alongside the budget so it is visible.

Every quad costs the same regardless of level — depth is free, breadth is not.

When exhausted, stop and report what was left queued. Order the queue by priority (spread, or
spread × area) and **report whether the order mattered** — run the same budget with a shuffled queue
and compare the trees.

---

## 4. What to measure

The tree and the decisions are the output. Dump per quad: `level`, bounds, `ensemble_spread`
(mean/median/p90), `alpha`, `alpha_sibling_spread`, `error_ratio`, `worst_energy_drift`, the decision
taken, and the reason.

**Answer these:**

1. **Does it terminate?** Run to a large budget with no `max_level`. Report leaf count against
   iteration and the final depth histogram. If it does not terminate, **that is the headline result** —
   report where it kept splitting and what α was doing there.
2. **Does the floor engage?** What fraction of leaves floor rather than hitting the budget or
   `max_level`? If ~0, the floor branch is not working.
3. **Is the tree sensible?** Render leaf boundaries over the outcome image. Dense at fractal
   boundaries, sparse in smooth regions — visible, and the honest check.
4. **Does per-quad noise cause thrash?** Do neighbouring quads in a chaotic region get opposite
   decisions? Quantify: fraction of adjacent leaf pairs at different levels with similar spread.
5. **Does the sibling-spread policy beat the α policy?** Same budget, both policies, compare trees.
6. **Does budget order matter?** Priority-ordered vs shuffled, same budget.
7. **Threshold sensitivity.** Sweep `tau_display` and the α threshold; report how leaf count and
   depth distribution move. A criterion whose output is dominated by an arbitrary threshold is not a
   criterion — that finding would matter.

---

## 5. Cautions

**Do not tune to a nice-looking tree.** The picture is a diagnostic; the sweep is the result. A
threshold chosen because the image looked right is an arbitrary constant, which is the defect that
has already disqualified two candidate designs on this project.

**Do not smooth α over neighbouring quads.** It is the obvious variance reduction and it is wrong
here: α varies smoothly *except at boundaries*, and boundaries are exactly what refinement is
deciding about. It would blur the signal being detected.

**Watch for a criterion that never says stop.** The conservative failure direction means
contamination pushes toward *refine*. Under a budget that is not benign — it spends everything on the
chaotic sea and starves the structure. **If that happens it is a finding, not a bug to tune away.**

**Report negative and messy results.** Every PR in this project has corrected something stated with
more confidence than it deserved; that is the process working.

---

## 6. Definition of done

- Descends from one quad to a tree, adaptively, at a fixed playhead
- Never pools children to synthesise a parent
- Both policies implemented (α alone; sibling spread) and compared at equal budget
- All seven questions in §4 answered with raw output committed
- Threshold sweeps reported, not tuned
- No eviction, no caching, no async, no camera — **if it grew one, that is a bug**
