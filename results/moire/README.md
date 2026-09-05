# `results/moire` — is the ribbon banding drawn by the render, the stepper, or the physics?

Harness `examples/moire.rs`. Committed artefacts: `moire.txt` and its panels are
`config_stability` at 256² (five arms, from an earlier cut of the harness that swept the cadence);
`config_basin/` and `config_stability_192/` are the 2026-09-05 pair, two arms each at a **matched**
raster, which is what makes them comparable.

```sh
cargo run --release --example moire -- 256 results/moire                          # moire.txt
cargo run --release --example moire -- 192 results/moire/config_basin config_basin
cargo run --release --example moire -- 192 results/moire/config_stability_192 config_stability
cargo run --release --example moire --  64 <scratch>/basin config_basin           # ramp_guard_64
cargo run --release --example moire --  64 <scratch>/ctl   config_stability
```

**These commands no longer reproduce the two 192² logs byte for byte, and the difference is named
rather than left to be discovered.** The ramp-guard `[ramp]` line was added to the harness *after*
those two runs, so re-running them now emits one extra line per arm that the committed files do not
carry. Nothing else moved: the table rows, the panels and every number quoted below are from the
same code path, and the guard is a separate print of values `colour::range` was already computing.
The 64² pair was run with the later binary and is where the guard output lives. *A documented
reproduction command can be wrong, and only running it finds out* — this one was checked against
what actually ran rather than reconstructed from memory.

## Read `l dst` and `l sd` first, then the ramp guard, then `prom`

A low prominence has two causes and no spectrum can separate them: a **structured** window with no
periodic beat in it, and a **featureless** window with nothing in it at all. This is the standing
*a difference can be small because both sides are right or because both are dead*, at a spectral
statistic, and the harness carries three arms against it.

**`l dst` / `l sd`** — the distinct 8-bit luminance levels the render paints, and their spread. A
window carrying structure has hundreds; a dead one has a handful.

**The ramp guard, both arms.** The lightness window is each case's own p1–p99, so a field with no
dynamic range has its **noise** stretched to full scale and every number below is computed on the
stretch. A span test alone is not enough — the record's `far` failure cleared a span of ×8 — so the
second arm compares the window's **floor** against the region's own median energy drift.

Measured at **64²** — a p1–p99 window is a *distribution* statistic and 4096 samples estimate it
well, unlike the chord and max statistics this project requires at the shipping resolution. The
spectrum table below is at 64² and 192², and says which.

| case | window (64²) | span | drift p50 | **lo/drift** |
|---|---|---|---|---|
| `config_basin` | (2.187e-3, 2.078e-2) | ×9.5 | 1.516e-9 | **1.44e6** |
| `config_stability` | (7.320e-5, 4.926e-1) | ×6730 | 3.352e-6 | **2.18e1** |

Both pass, and they pass differently: the basin's range is narrow but sits six orders clear of the
integrator's arithmetic, while `config_stability` has four decades of range whose *floor* is only
22× the drift. Neither is an amplified noise floor.

## The measurement, matched raster, on the FLOAT field

`lam:8bit` and `lam:f64` are the same statistic on the rendered luminance and on `spread_shape`
itself. **Read the float column** — the 8-bit one can be drawn by quantisation contours of a smooth
ramp, which is one of the six mechanisms the `osc/` work had to exclude.

| raster | case | `lam:f64` | **`prom`** | `t_end dst` | `on bnd` | `l dst` | `l sd` |
|---|---|---|---|---|---|---|---|
| 64² | `config_basin` | 59 | **8.20** | 1 | 1.0000 | 183 | 51.84 |
| 64² | `config_stability` | 45 | **1.60** | 918 | 0.7771 | 173 | 51.54 |
| 192² | `config_basin` | 57 | **8.14** | 1 | 1.0000 | 184 | 49.21 |
| 192² | `config_stability` | 64 | **2.11** | 8455 | 0.7708 | 197 | 54.81 |
| 256² | `config_stability` | 55 | 2.14 | 14916 | 0.7723 | — | — |

**`config_basin` bands at 3.9× the prominence of the field the banding was established on**, and its
value is raster-stable — 8.20 at 64², 8.14 at 192², across a 3× change — where the control's is
still converging (1.60 → 2.11 → 2.14).

## So the standing prediction is REFUTED, and the reason is arithmetic

`results/osc/README.md` records: *"`config_basin` has no beat at all … its window is `zoom =
0.009095`, **70× tighter** than `config_stability`'s 0.63763, so the pair period does not vary
measurably across it. **Prediction, untested: `config_basin` should show no ribbon banding.**"*

It bands, more prominently than the reference. The premise is where it goes wrong: **`zoom` is
compared across two charts with different magnification.** `Chart::config_slice` puts `mag` on both
basis vectors, and it is **4.0** for `config_basin` against **1.0** for `config_stability`.

| | `zoom` | `mag` | latent half-span |
|---|---|---|---|
| `config_stability`, full view | 0.63763 | 1.0 | **0.63763** |
| `config_basin`, full view | 0.009095 | **4.0** | **0.036381** |
| `config_stability`, the ribbon window this harness measures | 0.63763 × 0.18 | 1.0 | **0.11477** |

So the basin window is **17.5×** tighter than the full view and only **3.15×** tighter than the
window the banding was actually measured in — not 70×. The record's own
*a default that spans two coordinate systems silently means two different things* (`half = 0.05` as
a body position against a sigmoid pre-image), at `zoom` across two magnifications.

**And a lag of "exactly 0.0000" is a resolution statement, not a zero.** The survey's probes sit on
interior grids, so their separation scales with the window; a phase shift 3× smaller than the one
the estimator was calibrated on can read as zero without being zero. `band_guard.rs` separates lag-0
from two frozen inputs and two identical ones — all three arms pass — but it does not ask whether
the estimator *resolves* the shift it is looking for.

## What the basin run does confirm, directly rather than by correlation

**Nothing terminates.** `t_end dst = 1` at `on bnd = 1.0000` at both rasters: every trajectory in
the window ends at the same time, the horizon. That is the regular-island premise as a measurement.

**And a 4× refinement of `eta` moves not one pixel.** The two cases clear the stepper by opposite
routes, and the difference between the routes **is** the regularity:

| 192² | `steps p50` | × | pixels moved | worst Δ | `prom` |
|---|---|---|---|---|---|
| `config_stability` | 1.326e5 → 4.593e5 | 3.46× | **1496 / 36864 (4.06%)** | 3/255 | 2.11 → 2.11 |
| `config_basin` | 3.000e5 → 6.281e5 | 2.09× | **0 / 36864 (0.0000%)** | 0/255 | 8.14 → 8.14 |

The control is the arm that matters: its panels *do* move, so the harness can see a stepper change
in this field — and the spectral peak is unmoved anyway. That is *never conclude "no effect" from
an aggregate without the per-pixel distribution*, satisfied rather than assumed: 4.06% of pixels
moved and the conclusion survives because the movement is 3 levels of 255 and the peak does not
budge. The basin's zero is then a statement about the basin, not about the instrument.

At 8 bits that zero says *converged below the display precision* rather than exactly equal. It is
still the cleanest demonstration of regularity in this project, because it needs no threshold.

**And the committed file sizes say it a third way, for free.** At identical 192×192 the basin
panels are **8,558 bytes** against the control's **108,523** — 12.7×. A regular field is
low-entropy and a chaotic one is not, which is the same measurement that showed the stale
`_uniform*` gallery panels to be speckle rather than merely old. Check the dimensions before
reading a file size, always; here they match.

So the window **is** a regular island, and being one does not remove the beat. If anything it
cleans it: on the chaotic slice the beat is buried under divergence, which is why the `osc/` work
needed a high-pass and a 2D transform to see it at all, and why the 1-D form reads only ~2.1 there.

## A number withdrawn, and why it was wrong

A first pass applied a 2D detrended spectrum to the **8-bit panels** and read `prom = 80,459` for
the basin against `85.78` for `config_stability` — a 938× that was reported and then withdrawn.
8-bit quantisation of a smooth ramp produces exactly periodic contour bands, and that is the fourth
of the six mechanisms `osc/` excluded; nothing here excluded it. The tell was in the same output:
the harness's float column reads `lam = 57` where the quantised-panel tool read `lam = 34`. **Two
instruments disagreeing about the wavelength of the same field means one of them is measuring the
display.** The float column is the one quoted above.
