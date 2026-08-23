# Build brief: the vertical slice

**Deliverable:** the whole system end to end — physics → ensemble → quadtree → screen floor →
adaptive render — with **texels at their true per-quad sizes**. Plus deep-zoom relative coordinates
and the linearised decoder, and two parameter sweeps that are user-facing.

**Why now.** Every build so far was deliberately isolated, and **each isolation hid something the
next one found**: pooling hid a +38.6% surrogate error; the missing quadtree hid every dynamic
question; the missing screen floor hid that PR #11's termination result describes a regime the real
system never enters. **The remaining errors are in the interactions**, which is what this build
exposes.

**Scope.** Still no interaction, no eviction, no async, no GUI. But the camera, the screen floor,
the relative-depth predicate, SSAA resolve and adaptive rendering are all **in**, because without
them nothing measured is representative.

---

## 0. Three corrections made to this brief during the build

Recorded here rather than silently applied, because two of them are the difference between a
measurement and a decoration.

**0.1 — §3.4's curvature term is identically zero on the chart this repo had.** `Slice::decode_pos`
is a linspace and `Slice::nominal` writes `(x, y)` into `r[body]`, so the decode was **affine**:
`J_D` is constant and `x = x0 + J_D . delta` is *exact*, not an approximation. "Where does the
linearisation start to matter" would have answered "never", at every depth, in exact arithmetic.
Resolved by adding a **shape-sphere chart** — the closed-form inverse of the Hopf map that
`physics::shape::shape_vec` already computes forward, so it is an inverse rather than an
invention, and `shape_vec(decode(u, v)) == n(u, v)` is a round-trip test that can fail.

**0.2 — §3.5's "latent-aligned slices (two `z` axes varied)" has no referent here.** There is no 8D
`z` chart in this repo; the concept comes from the full Principia design. Resolved by generalising
`Slice` with a `Chart` enum: axis-aligned, oblique 2-plane in the 6D position space, and shape
sphere, with the historical behaviour preserved **bitwise** and asserted so.

**0.3 — §3.4's formula does not work as written, and the failure is invisible.** At depth 40 a quad
spans ~1e-13; `x0` is O(1), so an f32 sum `x0 + J_D . delta` returns `x0` for **every** delta and
all `N²` samples collapse to one initial condition — exactly what the direct path does there. The
divergence ladder would then compare two *collapsed* sets, find them in agreement, and report the
linearised path tracking f64 beautifully from a path that had lost every sample. Resolved by
making **distinctness the primary measurement** and specifying three named arithmetic variants
rather than one formula. See §3.4 as amended.

A consequence that bounds what §3.4 can conclude, stated before the run: **the initial conditions
must be formed as absolute O(1) numbers before integration**, because the three-body separations
are O(1) and no nonlinear integrator can carry `(x0, delta)` separately through the march. So the
linearised path can escape a floor set by the *chart coordinate*; it cannot escape one set by the
*IC magnitude*. That is a weaker claim than "no floor".

---

## 1. What PR #11 could not see, and why it matters

### 1.1 The screen floor is the everyday stop and it was absent

The scheduler contract:

> *"Screen-space floor — `tile_size(quad, zoom) <= pixel_size` -> stop refining. Once a quad's tiles
> have shrunk to pixel size, splitting further produces **sub-pixel** samples that cannot be
> displayed distinctly — wasted compute by definition… **This is the everyday refinement stop: in
> normal exploration you hit it far shallower than any precision floor.**"*

And it is a **veto**, with complexity the sole trigger — *"the tile-to-pixel ratio is never itself a
reason to refine; the screen floor is a veto, complexity the sole trigger."*

**The arithmetic.** One sample, one tile, no interpolation. At `N=8` a quad is 64 samples, so a
fully-refined tree at level `L` holds `4^L x 64` samples. For a 512² viewport:
`4^L x 64 = 262144` -> **L = 6**. PR #11's descent reached **level 12** in near-field: **4096x past
the point where samples stop being displayable.**

So q1's "terminates at 4617 quads" and q7's "tau dominates" are findings about the criterion
**minus its principal stop condition**. They must be re-run.

### 1.2 The split predicate under test was the superseded one

The contract explicitly rewrites it, because the absolute form *"caps infinite zoom at ~14"*:

```
split(C) <=> S_quad > tau(l)  AND  l < camera_depth + MAX_REL_DEPTH
```

PR #11 used absolute `max_level`. It had no camera, so it could not do otherwise — but
`MAX_REL_DEPTH` is the real predicate and it is **view-relative scheduler state, never on the sim
key**. Lowering it while zoomed invalidates no payload; it just stops scheduling deeper. Sensible
range **4–8 below the view**.

### 1.3 The render showed boundaries, not the system

PR #11's overlay draws leaf boundaries on a **uniform** render, so **every texel is the same size**.
That is what the earlier brief asked for, and it is the wrong instrument: it shows where boundaries
fell, not what the system displays. **A leaf at level 3 must be drawn with 4x the texel size of a
leaf at level 5.** Until then the tree's quality cannot be judged by eye at all.

### 1.4 SSAA is untested

The ensemble is currently used only for spread. Its other job is **resolve**: many sub-pixel samples
-> one pixel colour. That path has never run.

---

## 2. The reframe — the criterion is a priority signal, not a stop condition

This follows from §1.1 and it changes what "improving the criterion" means.

If the screen floor is the everyday stop, then in normal use **the criterion is not deciding when to
stop.** It is deciding **which quads get the frame budget first, within what is displayable.**

- **Thresholding** needs alpha resolved per quad against a fixed cutoff — exactly what the measured
  chaotic scatter of 1.1–1.3 makes unreliable.
- **Ranking** needs alpha only to order quads correctly *on average*. Noise that flips a threshold
  decision often preserves a ranking, and errors wash out over thousands of decisions.

**The per-quad noise is far less damaging to a priority signal than to a stop condition.** PR #11
supports this: `spread` and `spread x area` gave **byte-identical** trees, so the ordering is robust
to formula changes even where alpha is not robust per quad.

**So the improvement question is "what is the best priority ordering", not "what threshold".** Report
priority-quality directly: does the budget get spent on the quads that most change the image?

---

## 3. What to build

### 3.1 Camera and screen floor

A camera stub: centre, zoom, viewport in pixels. No interaction; set it, render, dump.

```
tile_size(quad, zoom) = quad_width * zoom / N        -- N tiles across a quad
screen_floor(quad)    = tile_size <= pixel_size
camera_depth          = level whose quad width ~ viewport width
```

- **Screen floor is a VETO**: never a reason to split, always a reason not to. **View-relative, not
  terminal, not cached as a quad fact** — zoom in and the same patch covers more screen, tiles regrow
  above pixel size, refinement resumes with real new samples.
- **`MAX_REL_DEPTH` replaces absolute `max_level`**, default 6, `MAX_REL_DEPTH <= screen floor`
  always.

### 3.2 Adaptive render — texels at true size

Each leaf rasterises its `N x N` samples across **its own screen footprint**. A level-3 leaf's texels
are 4x the linear size of a level-5 leaf's. **One sample, one tile, no interpolation** — never
upsample a coarse quad to fill pixels smoothly, which would fabricate structure.

**Acceptance:** measure texel size per leaf and assert it varies with level as `2^-level`. A render
where all texels are equal is the PR #11 failure and must fail the test.

### 3.3 SSAA resolve

Resolve the `E+1` copies to one pixel colour. Note the copies serve **two** purposes — spread
(scheduling) and resolve (display) — and they must not be confused: spread is a *disagreement*
statistic, resolve is an *average*.

Report what resolve does to the image at each `E`, and whether a resolved render differs visibly from
a nominal-copy-only render.

### 3.4 Deep zoom: relative coordinates and the linearised decoder

From the canonical spec:

> *"CPU (f64) owns global/nonlinear precision (quad centre/half-width, `x0`, `J_D`); the GPU (f32)
> does only quad-local relative arithmetic (`x = x0 + J_D . delta`)."*

Implement both paths and compare:

- **direct**: decode absolute coordinates at full precision
- **linearised**: per quad, CPU computes `x0` (centre decode) and `J_D` (Jacobian, two f64 decodes
  per axis); sample positions come from `x = x0 + J_D . delta` with `delta` quad-local in `[-1,1]`

**What to measure — as amended, see §0.3.** **Distinctness first**: per quad, per path, per depth,
how many of the `N²` sample ICs and the `E+1` copies are actually distinct. A divergence figure is
admissible only where both sides are fully distinct; elsewhere the rung reports collapse and no
agreement number at all.

Then the divergence, over three named variants rather than one formula:

| variant | arithmetic |
|---|---|
| **direct-f64** | `u = cu + du*half` in f64, decoded in f64 |
| **L-naive** | `x0`, `J_D . delta` and the sum all in f32 — the literal formula |
| **L-split** | `x0` in f64 on the CPU; `delta` and `J_D . delta` in f32 (both O(1), neither shrinking); promoted and summed in **f64** |

The claimed benefit is extending usable zoom from depth ~23 to ~50+. Two distinct limits: the
linear-decoder `AT_F32_FLOOR`, and PR #11's plain-f64 cell-width floor at level 45.87. **Whether
the linearised path escapes the second is the thing to verify rather than assume** — and given
§0.3's bound, the honest prior is that it does not.

Also report **Jacobian cost**: two f64 decodes per axis per deep quad. The caching contract records
that this piled up badly enough to hitch a gesture, so the per-quad cost is worth having a number for
even without interaction.

**Note the linearisation is an approximation** — it discards curvature. Report where `|direct -
linearised|` exceeds the sample spacing, since that is where the approximation starts to matter.

### 3.5 Slice variety

Every experiment so far is three regions of one Burrau slice family. The chart is 8D and slices
through it are not equivalent.

Run at least: **latent-aligned** slices (two `z` axes varied, others fixed) and **arbitrary/oblique**
slices (a rotated 2-plane). Report whether tree shape, leaf count and the alpha distribution differ
by slice type. If oblique slices behave differently, every prior conclusion is slice-conditional.

---

## 4. The two sweeps

### 4.1 Ensemble copies `E` — user-facing, tier-gated

**Prediction, stated before the run:** low `E` biases toward *refine*, exactly as low `N` did. Same
mechanism — a noisy spread estimate inflates apparent disagreement, and the conservative failure
direction turns that into extra splits. PR #11 measured this for `N`: **N=4 spent 4x the quads of
N=16.**

**If it holds for `E`, the cheap tier spends more quads than the expensive one**, and the tier design
partly cancels itself.

Sweep `E+1 in {2, 4, 8, 16, 32}` at fixed everything else. Report leaf count, depth distribution,
and — the quantity that matters — **total trajectories = leaves x N^2 x (E+1)**. If leaf count rises
faster than `E` falls, the cheap tier is a false economy end to end.

Report the **spread estimator's own noise vs E** separately, so the mechanism is visible rather than
inferred.

### 4.2 Re-run PR #11's q1, q2, q3, q7 with the screen floor

These four were measured without the veto and are conditional on its absence. **q3 especially**: if
near-field descended 6 levels past displayable, its 4617 quads are not evidence the tree is right,
and deep interior's 29 may look entirely different when both are capped at what is visible.

Also carry forward the two open items from that review:

- **Does floored correlate with `worst_energy_drift` in deep interior?** 40.9% floored is the highest
  of the three regions, and it is also where the integrator works hardest. If the floor fires because
  alpha is corrupted by integration error rather than because the physics is irreducible, that is a
  different bug with a different fix — and the floor fraction alone cannot distinguish them.
- **Does p90 aggregation fix deep interior's tree?** The median-blindness attribution for q3 is
  plausible and unproven. If p90 descends where median does not, it is aggregation. If p90 also
  stalls, it is not.

---

## 5. Cautions

**Do not tune to a nice-looking tree.** Now that the render is adaptive the temptation is stronger,
because the picture finally looks like the product. The sweep is the result; the picture is a
diagnostic.

**The screen floor is a veto, never a trigger.** Sub-pixel tiles are a reason to stop, never a reason
to continue. The contract strikes the "below screen resolution -> refine" rule explicitly.

**Do not cache the screen floor as a quad fact.** It is view-relative and evaluated live. A quad
floored at one zoom must refine again when zoomed into.

**Watch for the interaction that is not visible in isolation.** Every isolated build hid something
the next found. Expect this one to as well, and look for it in the *seams* — screen floor against
alpha, SSAA against spread, linearised decode against the criterion — rather than in the components.

**Report negative and messy results.** Every PR in this project has corrected something stated with
more confidence than it deserved.

---

## 6. Definition of done

- Adaptive render with **texel size varying as `2^-level`**, asserted by test
- Screen floor implemented as a veto, view-relative, uncached
- `MAX_REL_DEPTH` replacing absolute depth in the split predicate
- SSAA resolve, with its effect on the image reported per `E`
- Both decode paths, with divergence-vs-depth and the two floors separated
- Latent-aligned and oblique slices compared
- `E` sweep reported in **total trajectories**, not leaf count alone
- q1, q2, q3, q7 re-run under the screen floor, with the two open items answered
- Still no interaction, no eviction, no async, no GUI — **if one appeared, that is a bug**
