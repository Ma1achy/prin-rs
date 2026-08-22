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

`raw/<region>-64.raw` — 64×64 rather than 256×256, because a 256² dump is 21 MB and seven of them
do not belong in a repository. The format is self-describing: magic `PRIN`, a version, a
length-prefixed text header carrying every parameter and the field names, then 40 `f64` per pixel
in row-major order (`index = jy*nx + jx`, x fastest). Fields are always `f64` regardless of kernel
precision, so an f32 run and an f64 run produce byte-comparable dumps.

Read one with:

```python
import struct
d = open("raw/near-field-64.raw", "rb").read()
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
| `output/refinement_criterion.txt` | Experiment A — the refinement exponent, its control, the noise floor against `E` |
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

## A note on refinement

The experiment examples run with `refine_flagged: false`, deliberately. Experiments A and B
characterise the kernel whose behaviour motivated the second pass, and measuring the repaired
kernel would hide the thing being measured. Precision comparisons also run it off: the pass is
threshold-triggered on `error_ratio` and f32 and f64 flag different pixel sets, so with it on the
comparison would be of pipelines rather than of arithmetic.

The `render-*.txt` runs have it **on**, and report the before and after drift maxima on the same
line.
