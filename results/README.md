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

### Zoom ladders — `zoom_<region>_NN.png`, `zoom_<region>.apng`

Nine frames at 384×384, each re-descending from a root box of `half = 0.05 / 2^k` with the camera
framing it. **This is the only artefact that shows the screen floor is view-relative**, and a
still image cannot: the 252 quads floored in frame 0 are refined in frames 1–3 with genuinely new
samples, and a fresh population — 868 by frame 2 — is floored in their place. Nothing is cached
and nothing is upsampled.

The animation is an **APNG**, written by the `png` crate rather than by adding a GIF dependency.
Every frame is also on disk as an ordinary PNG, so nothing here depends on APNG support to be
readable.

### Tree dumps — `<region>.prnq`, `zoom_<region>_NN.prnq`

The same self-describing `PRNQ` format as `results/sched-*.tree`: a magic, a version, a text
header naming every parameter including `chart`, then one record of 24 `f64` fields per quad.
One dump per rendered tree and one per zoom frame, so every picture above can be checked against
the numbers that produced it.
