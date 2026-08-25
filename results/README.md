# Committed results

Raw output and images for every experiment. [`../RESULTS.md`](../RESULTS.md) is the findings
document; this directory is the evidence behind it.

Everything here is regenerable — `output/*.txt` is the captured stdout of the correspondingly
named example, and the images and dumps come from `src/bin/prin.rs`. Committed anyway, because a
findings document whose numbers cannot be checked without a 20-minute re-run is a summary rather
than a record.

## Images

`<region>_outcome.png` and `<region>_spread.png`, 256×256, f64, `t = 13`, `E+1 = 8`,
`eta = 1e-2` with the flagged-pixel refinement pass on.

The **outcome** image is coloured by `state` with `detail` shading it (BRIEF §2.4): red escape,
green collision, blue bounded, yellow running, magenta for any pixel with a non-finite copy.
`detail = 3` — the "all three" outcomes — takes the brightest shade of its family.

The **spread** image is `ensemble_spread`, **log scaled between the grid's own p1 and p99**. A
linear ramp paints almost everything at the bottom of the scale, since the median sits near
`1e-3`, and the filament structure disappears. The window is printed in the matching
`output/render-*.txt`. Non-finite is painted at full scale: undetermined is the loudest thing a
pixel can be and must not be shown as quiet.

`near-field-f32_*` is the same slice at f32, for direct comparison against `near-field_*`.

## Raw dumps

**The dump is the product; the images are diagnostics.** Both resolutions are committed:

- `<region>.raw` — 256×256, the same runs as the committed images, 21 MB each.
- `raw/<region>-64.raw` — 64×64, ~1.3 MB each, for reading and testing a parser against.

Together they are about 155 MB, which is large for a repository and is a deliberate choice: the
findings in [`../RESULTS.md`](../RESULTS.md) are re-derivable from these files without a
re-run. If that ever needs undoing, `git lfs migrate` moves them out of the tree retroactively.

The format is self-describing: magic `PRIN`, a version, a length-prefixed text header carrying
every parameter and the field names, then 40 `f64` per pixel in row-major order
(`index = jy*nx + jx`, x fastest). Fields are always `f64` regardless of kernel precision, so an
f32 run and an f64 run produce byte-comparable dumps.

`near-field.raw` and `near-field-f32.raw` are the same slice at the two precisions, so a
field-by-field diff between them is the f32 question answered directly from the data.

Read one with:

```python
import struct
d = open("near-field.raw", "rb").read()   # or raw/near-field-64.raw
assert d[:4] == b"PRIN"
hl, = struct.unpack_from("<I", d, 8)
hdr = d[12:12+hl].decode()
off = 12 + hl
n,  = struct.unpack_from("<Q", d, off); off += 8
nf, = struct.unpack_from("<I", d, off); off += 4
fields = [l for l in hdr.splitlines() if l.startswith("fields=")][0][7:].split(",")
rows = [struct.unpack_from(f"<{nf}d", d, off + i*nf*8) for i in range(n)]
```

## Captured output

| file | what it measures |
|---|---|
| `output/pooled_vs_true_parent.txt` | the pooled block against a real two-resolution parent, and the corrected per-quad scatter |
| `output/halton_noise_floor.txt` | offset schemes: noise floor, parent/child correlation, the `alpha_E` control variate |
| `output/scheme_f32_effect.txt` | whether the fixed prefix changes the f32 answer (one pixel, not the scheme) |
| `output/refinement_criterion.txt` | Experiment A — the refinement exponent, its control, the noise floor against `E` (**pooled**; superseded by `pooled_vs_true_parent`) |
| `output/convergence.txt` | Experiment B — which conclusions survive at large `n` |
| `output/worst_128.txt` | the seven high-drift pixels at 128², and drift against `eta` |
| `output/refine_pass.txt` | the flag-then-re-integrate remedy, on against off |
| `output/latching_decision.txt` | near-tie ratios, persistence, and why the latch ships inert |
| `output/f32_report.txt` | all four {f32,f64} × shared-reference {on,off} combinations |
| `output/spread_branch_cut.txt` | whether `spread_shape` inherits the LC branch cut |
| `output/lc_cut_proximity.txt` | branch-cut conditioning against orientation |
| `output/lc_branch_effect.txt` | what conditioning the inverse LC map is worth |
| `output/deep_interior.txt` | BRIEF §2.6's pixel, probed pair by pair |
| `output/r_coll_sweep.txt` | `r_coll` sensitivity, and whether the branch cut reaches the labels |
| `output/spread_event_correction.txt` | event class against terminal outcome |
| `output/escape_check.txt` | when the escape arm first fires |
| `output/flag_effect.txt` | per-pixel effect of flags whose aggregates do not move |
| `output/resolution_confound.txt` | `sigma_E(0)` against cell width; `ref_disagree` against `error_ratio` |
| `output/collision_scan.txt` | the two-body collision gate against `eta` |
| `output/dmin_where.txt` | where `d_min` is set along a trajectory |
| `output/worst_pixels.txt` | the damaged tail, and MAD against max deviation |
| `output/bench_deriv.txt` | `deriv` throughput, extrapolated to 1024²×8 |
| `output/render-*.txt` | the 256² renders, one per region |

## The scheduler

`sched-*_tree.png` — leaf boundaries drawn over a 512² uniform render of the same box, one per
region. Deeper leaves are brighter, so the depth is visible without a legend, and the base is dimmed
so the boundaries read. **The picture is a diagnostic; the threshold sweep is the result** — a
threshold chosen because the image looked right is an arbitrary constant.

`sched-*.tree` — the tree dump, magic `PRNQ`, one record per quad with 24 fields: `level`, bounds,
`cell_width`, the three spread aggregations, `alpha` under each of them, `alpha_sibling_spread`,
`error_ratio_max`, `worst_energy_drift`, and the decision taken. Same self-describing header shape
as the pixel dump, carrying every threshold, the policy, the order, the aggregation, the budget and
the trajectory cost per quad.

| file | what it measures |
|---|---|
| `output/sched_sweep.txt` | §4 q7 — threshold sensitivity. **Runs first, because it is what sets `tau`** |
| `output/sched_terminate.txt` | §4 q1, q2 — does the descent terminate, does the floor engage |
| `output/sched_thrash.txt` | §4 q4 — adjacent leaf pairs at different levels with similar spread |
| `output/sched_policies.txt` | §4 q5 — sibling-spread policy against the alpha policy, equal budget |
| `output/sched_order.txt` | §4 q6 — priority against shuffled, same budget |
| `output/sched_n_sweep.txt` | the N sweep, including `N = 7` to vary parent–child CRN strength |
| `output/prinq-*.txt` | the per-region descents that produced the overlays |

## Test and cross-check output

`tests/*.txt` is the raw output of each test binary, run with `--nocapture --test-threads=1` so
the printed tables appear in order. These carry BRIEF §5's acceptance evidence, which otherwise
lives only in pull-request descriptions.

| file | what it carries |
|---|---|
| `tests/burrau_constants.txt` | `M=12`, `R=2.2361`, `E=-12.8167`, crossing time |
| `tests/az_identities.txt` | the validation chain: energy round-trip, `Gamma == A*B*(H-E)` |
| `tests/az_hamiltonian_fd.txt` | finite-difference `Gamma` against analytic `deriv`, h² convergence |
| `tests/two_body_collision.txt` | gate (b): `d_min < 1e-10` with `\|dE/E\| < 1e-12` |
| `tests/gauge_invariance.txt` | the `alpha` rescaling gate — catches an absolute-length leak |
| `tests/error_ratio_acceptance.txt` | the `error_ratio` gate, MAD against max deviation, `d_min_gap` |
| `tests/error_ratio.txt` | the five invariants, including step-size convergence |
| `tests/outcome_encoding.txt` | BRIEF §2.4's encoding, the >=2-pair rule, scale invariance of `t_end` |
| `tests/f32_precision.txt` | the floor divergence, and gate (b) parameterised by precision |
| `tests/quadtree.txt` | quadtree geometry: exact cell-width halving, tiling, no pooling, precision floor, budget |
| `tests/halton_offsets.txt` | the fixed prefix: radical inverse, fixedness across resolutions, discrepancy against PCG |
| `tests/lc_conditioning.txt` | inverse-LC branch conditioning |
| `tests/spread_branch_cut.txt` | whether `spread_shape` inherits the branch cut |
| `tests/xcheck.txt` | gate (c): per-column comparison against the NumPy reference |
| `tests/horizon.txt` | the divergence-vs-horizon table, both sides conditioned |
| `tests/horizon-lc-unstable.txt` | the same table on the reference's original LC branch |

The two horizon tables are the reason the branch-cut finding is legible: `az_t13` reads
`1.930e-10` on the original branch and `2.718e-13` with both sides conditioned. Keep them
together — the pair is the evidence, either alone is not.

`tests/xcheck.txt` needs Python and NumPy; regenerate it with
`cargo test --release --test xcheck -- --ignored --nocapture`, and the horizon tables with
`python3 tools/xcheck/horizon.py [--lc-unstable]`.

## A note on refinement

The experiment examples run with `refine_flagged: false`, deliberately. Experiments A and B
characterise the kernel whose behaviour motivated the second pass, and measuring the repaired
kernel would hide the thing being measured. Precision comparisons also run it off: the pass is
threshold-triggered on `error_ratio` and f32 and f64 flag different pixel sets, so with it on the
comparison would be of pipelines rather than of arithmetic.

The `render-*.txt` runs have it **on**, and report the before and after drift maxima on the same
line.

## `vertical/` — the vertical slice

Everything the adaptive render, SSAA and zoom-ladder runs produced. **These are the first images
in this project drawn at true per-quad texel sizes**; every image above and every overlay in
PR #11 is a *uniform* render, where a level-3 leaf and a level-5 leaf are drawn at the same size.

### Adaptive renders — `<region>_*.png`

Four panels per region, from `examples/adaptive_render.rs`, 512×512, camera framing the root box,
screen floor on:

| file | what it is |
|---|---|
| `_adaptive_spread.png` | `ensemble_spread`, texels at true size, **fixed** window `1e-8 .. 1e-1` log — comparable across regions |
| `_adaptive_spread_auto.png` | the same, on this tree's own p1..p99 window — legible within a region |
| `_uniform_spread.png` / `_uniform_spread_auto.png` | **the wrong instrument, kept deliberately.** Every texel the same size, which is what PR #11 drew. It exists so the acceptance test has a negative case, and it fits a texel-scaling exponent of `+0.000000` where the adaptive render fits `-1.000000` |
| `_adaptive_outcome.png` | the nominal copy's `(state, detail)` label, same palette as the images above |
| `_adaptive_resolved.png` | the **SSAA resolve** — the mean colour over the `E+1` copies. Compare against `_adaptive_outcome.png`: they differ on 0.6% of near-field's footprints and 14.1% of `deep interior`'s, and not at all in `far` |

`deep_interior_adaptive_spread.png` is worth looking at directly: it is §5's q3 failure made
visible. The largest high-spread structures — the red regions top-left and bottom-right — are
drawn in the **coarsest** texels in the image, because the tree left them at level 2. A uniform
render cannot show that, which is why it went unnoticed.

### SSAA panels — `ssaa_<region>_<E+1>_{nominal,resolved}.png`

256×256 uniform renders at `E+1 ∈ {1, 8, 32}`. The `E+1 = 1` pair is the control and the two
images are **identical by construction** — the resolve of one copy is that copy. A difference
there would be a bug in the resolve rather than a finding.

### Zoom ladders — `zoom_<region>_NN.png`, `zoom_<region>_animated.png`

Nine frames at 384×384, each re-descending from a root box of `half = 0.05 / 2^k` with the camera
framing it. **This is the only artefact that shows the screen floor is view-relative**, and a
still image cannot: the 252 quads floored in frame 0 are refined in frames 1–3 with genuinely new
samples, and a fresh population — 868 by frame 2 — is floored in their place. Nothing is cached
and nothing is upsampled.

The animation is an **APNG** — written by the `png` crate rather than by adding a GIF dependency,
and named `.png` because an APNG *is* a PNG and viewers that do not animate show the first frame.
Every frame is also on disk as an ordinary PNG, so nothing here depends on APNG support to be
readable.

### Tree dumps — `<region>.prnq`, `zoom_<region>_NN.prnq`

The same self-describing `PRNQ` format as `results/sched-*.tree`: a magic, a version, a text
header naming every parameter including `chart`, then one record of 24 `f64` fields per quad.
One dump per rendered tree and one per zoom frame, so every picture above can be checked against
the numbers that produced it.

---

## `results/criterion/` — improving the criterion

### Reference and criterion renders — `<region>_reference.png`, `<region>_B682_<rank>.png`

`<region>_reference.png` is the **fully-refined tree at 512², one sample per pixel**. It is what
`error(B)` is measured against, and it is a *specific finite sampling* rather than the true
image: at the screen floor, which side of a filament a pixel lands on is an accident of where
its sample fell. `error = 0` means "matches this sampling".

The `_B682_` panels are each ranking's tree at one eighth of the full budget, drawn at **true
per-quad texel sizes** so a coarse leaf is visibly coarse. Compare `within_median` — the shipped
default — against `greedy_oracle` and `frac_hot_between`: the shipped one spends its budget in
visibly the wrong places, which is the picture behind §10.3's table.

`far_*.png` are included and are the control: they are featureless, which is exactly the point.
`error(root) = 0.00000` there, so the metric is **undefined** on `far` and every criterion reads
zero. It is not that they agree.

### Colour-coupling renders — `colour_<region>_<colouring>_*.png`

The production bivariate scheme — hue from the shape sphere, lightness from a selectable scalar —
against the diagnostic outcome colouring. §6's coupling question is whether `error(B)` reorders
the criteria when lightness switches from `spread` to `diffusion`; if it does, the criterion
needs a term for the lightness field and has none.

The hue map is `atan2(n2, n1)` with chroma tied to `sqrt(n1² + n2²)`, in OKLCh. The azimuthal
discontinuity is invisible **by construction**: hue is undefined at the poles and chroma goes to
zero there, so the two colours either side of the cut converge on the same grey.

### Wire twins — `*_wire.png`

**Every image in this directory has a `_wire` twin, and neither replaces the other.** The plain
render says *what is displayed* — texels at true per-quad sizes, so a coarse leaf is visibly
coarse. The wire says *where the tree cut*, with brightness graded by level so a boundary can be
attributed to a depth.

They answer different questions. A coarse texel tells you a leaf is coarse; only the wire tells
you whether the structure around it was subdivided *around* it or straight *through* it. PR #11
drew boundaries over a **uniform** base, which conflated the two, and that is how `deep
interior`'s bad tree went unnoticed for a whole build.

### The budget animation — `budget_<region>_t<T>_animated.png`, `..._wire_animated.png`

The single most useful artefact here. Each frame is **`greedy_oracle` on the left,
`within/median` — the shipped default — on the right**, at the same budget, drawn at true texel
sizes. Side by side in one frame on purpose: two separate animations would make you hold one in
memory while watching the other, which is exactly the comparison the picture exists to remove.

Watch where the right-hand side spends its budget. §10.3's table says it is beaten by random past
`B = 383`; this is what that looks like.

Every frame is also on disk as `budget_<region>_t<T>_NN.png`, so nothing depends on APNG support.

### Slice gallery — `slice_<case>.png`, `slice_gallery_animated.png`

Ten charts through **one shared centre configuration**, so only the 2-plane changes: the
axis-aligned body plane, oblique planes at 15/30/45°, cross-body mixes, and the shape chart at
three fibre phases. Bases are orthonormal in the 6D position metric, so a unit of chart
coordinate moves the system equally far in each — otherwise a "different slice" would be a
different *scale*.

`slice_body_plane.png` and `slice_plane_00deg.png` are the **control pair**: the same chart
written two ways, asserted **bitwise identical**. If they differ, every other panel is comparing
different physics rather than different slices.

`slice_variety` measured tree size as slice-conditional to **4.3×** while the `alpha` distribution
stayed put (median 0.172–0.289). These are what that looks like.

### `error(B)` curves — `curve_<region>_t<T>.png`

Log y, because the curves span four decades and the interesting part is the bottom one. **An
exact zero cannot sit on a log axis, so it is drawn at the floor and ticked rather than dropped**
— on several of these curves reaching zero *is* the result. Controls (`greedy_oracle`, `random`)
are dashed so they read as references rather than candidates.

### The complete-tree cache — `<region>_t<T>.qcache`

`PRQC`: magic, version, a text header naming every parameter, then 50 `f64` per quad for all
5461 quads of the complete tree. **Without this the tree lived only in RAM for one process**, and
reproducing any §10 table meant paying the 2.8-million-trajectory integration again.

It carries `err_sum` — this quad's summed OKLab distance to the reference were it drawn as a leaf.
That is a **constant of the quad**, because quads are disjoint, which is what makes the replay
exact and the greedy priority queue static. It also carries **every criterion's scalar** whatever
the run ranked on, so criteria can be compared offline without re-integrating.

### Pan frames — `pan_<region>_NN.png`, `pan_<region>_animated.png`

Nine camera positions across the region. **The animation showing nothing change is the result**:
`Camera::veto` reads `tile_size_px`, which does not depend on `cx`/`cy`, so panning changes no
scheduling decision at all. The `_wire` twin makes that unmissable — the mesh is identical in
every frame.

### Quad dumps — `between_<region>.prnq`, `slice_<case>.prnq`, `pan_<region>_NN.prnq`

`PRNQ` **version 2**: the same self-describing format, now 48 `f64` fields per quad. v2 appends
the between-footprint arm, the two matched-count controls, the hot-set layout, the
termination-gradient pair, the cost column, the IC-distinctness count and the temporal
accumulators. A reader that indexes by name still works on v1; one that indexes by position does
not, which is why the version moved.

---

## `charts/` — every chart family, from `examples/chart_gallery.rs`

Thirteen chart instances across the reference's five families, all at **1024²**, budget 40000 so
every descent stops on the criterion rather than the cap. Per chart:

| file | what it answers |
|---|---|
| `<case>.png` | **what is displayed** — adaptive render, texels at true per-quad sizes |
| `<case>_wire.png` | **where the tree cut** — the same image with leaf boundaries, graded by level |
| `<case>_uniform.png` | **what the chart looks like** — one sample per pixel, no tree at all |
| `<case>_outcome.png`, `<case>_uniform_outcome.png` | the categorical control |
| `<case>_levels.png`, `<case>_levels_wire.png` | **how the tree got there** — one descent truncated at each depth |
| `<case>.prnq` | the quad dump, with `chart_params` carrying the full basis |

**Read all three of the first three, never one alone.** The adaptive render is a picture of the
*scheduler*: near-field's `alpha` median is 0.14 against `alpha_hi = 0.2`, so the criterion says
refinement does not pay and keeps coarse leaves — and a coarse leaf is one flat tile, because the
render never interpolates. That reads as blur and is an honest picture of an unrefined tree.
`_uniform` is the chart without the tree in the way. Reading either alone is how a criterion's
failure gets mistaken for a rendering artefact, or a rendering choice for a finding.

The level ladder is **one** descent truncated at each depth, not a fresh descent per frame — so
it is one refinement seen at several playheads rather than several unrelated trees.

## `colour/` — from `examples/colour_check.rs`

Three regions × three lightness fields at 1024², plus the outcome control. The point of the set
is the **pair** `<region>_spread.png` and `<region>_spread_shape.png`: `ensemble_spread` is
`max(spread_shape, spread_event)` and the event arm is a count ratio over `E+1` copies, so where
it dominates the field is a staircase. `results/output/colour_check.txt` prints the distinct-value
count and the event-arm fraction, and those are what to read before either image.

## Footprint caches — `<region>_t<T>.fcache`, **not committed**

`PRQF`: the colour-relevant projection of every footprint of a complete tree — enough to recolour
`error(B)` under any colouring without re-integrating. At level 7 that is 1.4M footprints × 14
`f64`, about 160 MB each and 940 MB for six, with no redundancy to remove at one sample per pixel.

They exist so a colouring change costs a **replay within a working session**. Committing a
gigabyte so that survives a clone is the wrong trade, so they are gitignored with the regeneration
command recorded beside them in `.gitignore`.

`PRQC` is still committed: it is the per-quad reductions and every criterion's scalar, which is
what every table is read from.

## A note on raster sizes

Everything generated in this build is **1024²** (figures are 1400×800; budget side-by-sides are
2048×1024). Anything smaller in `results/` predates it.

If an image looks blurry, **measure its pixel dimensions first**. Wireframe lines are written with
integer pixel `set` calls and adaptive texels are nearest-neighbour, so neither can be soft in the
file — softness is always a viewer upscaling a small raster. That diagnosis cost a round trip in
this build, after a `criterion_metric -- 3 8` validation run was allowed to overwrite committed
512² artefacts with 128×64 ones.

### Regenerating the top-level region images and raw dumps

`prin --size` drives **both** the PNG pair and the per-pixel `.raw` dump, and they want opposite
sizes: the images want 1024², while the dumps are 64×64 by design and would be **320 MB per
region** at 1024². So it is run twice per region, with the large run's dump kept out of
`results/`:

```sh
for r in near-field mid-field far "body2 core" "body1 slice" "deep interior"; do
  stem="${r// /-}"
  cargo run --release --bin prin -- --region "$r" --size 1024 --out "/tmp/big/$stem"
  cp "/tmp/big/${stem}_outcome.png" "/tmp/big/${stem}_spread.png" results/
  cargo run --release --bin prin -- --region "$r" --size 64 --out "results/$stem"
done
```
