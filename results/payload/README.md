# `results/payload` — the physics-space ground truth, and the tolerance policy against it

Phases 1 and 2 of the refinement rebuild (plan: *the refinement mechanism, rebuilt:
physics-space, live-playhead, memory-first*), measured 2026-09-03 on the fixed kernel
(`EnsembleCfg::production()`: Heggie, `StepLimit::Predictive`, `PerStepInterval`, landing clamp).
Every run here carries its config line; `refine_flagged=false` throughout, because the repair pass
is batch-only and has no live-playhead analogue.

## What moved, in one sentence

Every `error(B)` curve before this scored a tree by OKLab distance under the shipping colouring,
whose lightness is auto-ranged to each region's own p1–p99 — so a smooth region's `1e-8` residual
was stretched to full contrast and counted as error at every depth, and breadth-first came out
near-optimal by construction of the metric. Scored on the **payload** — the nominal copy's point
on the shape sphere and its event class, against a fixed tolerance `eps = 0.01` (chord/2, the
units of `spread_shape`) — the same footprints say the opposite: `far` is resolved at the root,
`near-field` at 137 quads against uniform's 5449, `deep interior` at 329 against 5377, and only
`preset_shape_h1` is a sea.

## The files

| file | what |
|---|---|
| `<target>_t13_L6.fcache` | PRQF **v2** footprints (with `event_class`), levels 6, `N = 8`, `E+1 = 8`, 512², `t = 13`, `n_sync = 32`. 2,796,032 trajectories each. Any metric, present or future, is a replay of this file. **Not committed** (`results/**/*.fcache` is gitignored, 42 MB each); regenerate with the `build` line below, 3–12 minutes per target. |
| `<target>_t13_L6.qcache` | the quad cache of the headline metric (`payload/event_class/indicator/eps=1e-2`), with reductions, for the signal rankings. |
| `output/payload_metric_<target>_L6.txt` | the build log: `sea_fraction(eps)` ladder, `error(B)`, `captured`, `headroom/err`, `B_needed`, per form and for the OKLab comparison arm, with the static `Policy::Alpha` tree scored beside them. |
| `output/payload_live_<target>_t13_L6_<policy>.txt` | one static descent under the named policy, mapped onto the cache and scored: `mem_all`, `mem_leaf`, error, and the budget the ceiling and uniform need for that error. |
| `output/payload_march_<target>_t13_L6.txt` | the **live** descent: grown boundary by boundary from one march per quad, the growth curve, the catch-up cost, the final tree scored the same way. |

Reproduce: `payload_metric build <target> 6 8 0.01 13 <root>`, then `payload_metric live <fcache>
<policy> 0.01 0.25 <root> <stationary>` and `payload_metric march <fcache> 0.01 0.25 <root> 1 4`.
The four builds took 682 s, 409 s, 710 s, 194 s on an idle 11-core machine; a live or march
descent is minutes.

## The metric, and the three forms

Per pixel, nominal copy: `d = max(chord(shape_tree, shape_ref)/2, [class_tree != class_ref])`;
`1.0` when exactly one side is unusable, `0` when both are, so `error(full)` is exactly zero.
**`indicator`** is the fraction of the image with `d > eps` — the headline. **`hinge`** is
`max(0, d - eps)`. **`resolvable`** counts only pixels whose deepest-level footprint has
`ensemble_spread <= eps`, i.e. pixels a finer grid could resolve; the gap to `indicator` is the
**sea cost**, and `sea_fraction(eps)` is printed before any curve. The event-class arm and the
outcome arm agree to five digits on every target at `eps = 0.01`: any pixel whose tightest pair
differs from the reference's already differs on the sphere by more than the tolerance.

## The numbers, levels 6, `eps = 0.01`, 512² camera, `k_frac = 0.25`

`B` counts quads computed (kept parents included — the memory model the user chose). `dp` is the
exact tree optimum (`Cache::dp_optimal`), `uniform` is breadth-first.

| target | root error | sea fraction | dp reaches 0 at | uniform reaches 1% at |
|---|---|---|---|---|
| `far` | 0.00000 | 0.0000 | `B = 1` | `B = 1` |
| `near-field` | 0.12464 | 0.0005 | 191 | 21 |
| `deep interior` | 0.07045 | 0.0025 | 341 | 249 |
| `preset_shape_h1` | 0.55365 | **0.3415** | 5461 (full) | never below 20% before 3071 |

| target | policy | quads | leaves | error | ratio to dp | ratio to uniform | stop |
|---|---|---|---|---|---|---|---|
| near-field | tolerance | 153 | 115 | **0.00000** | 1.12× | **0.03×** | `keep:81 screen_floor:34` |
| near-field | tolerance, stationarity off | 153 | 115 | 0.00000 | 1.12× | 0.03× | identical |
| near-field | alpha (legacy) | 21 | 16 | 0.00489 | 1.24× | 1.00× | `keep:16` |
| deep interior | tolerance | 353 | 265 | **0.00000** | 1.07× | **0.07×** | `keep:147 screen_floor:117 stationary:1` |
| deep interior | tolerance, stationarity off | 357 | 268 | 0.00000 | 1.09× | 0.07× | `keep:151 screen_floor:117` |
| deep interior | alpha (legacy) | 21 | 16 | 0.01882 | 2.33× | 1.24× | `keep:16` |
| preset_shape_h1 | tolerance, stationarity on | 3397 | 2548 | 0.02720 | 1.21× | 0.67× | `keep:450 screen_floor:2064 stationary:34` |
| preset_shape_h1 | **tolerance** (stationarity off, the default) | 3585 | 2689 | **0.00021** | **1.03×** | 0.66× | `keep:452 screen_floor:2237` |
| preset_shape_h1 | alpha (legacy) | 29 | 22 | 0.47227 | 1.16× | 1.00× | `floor:8 keep:14` |
| far | all three | 21 | 16 | 0.00000 | 21× (the bootstrap; dp needs 1) | 21× | `keep:16` |

Under the `resolvable` form, which ignores the sea proper, `preset_shape_h1` reads **0.0094 with
the stop on against 0.0002 with it off** (1.60× against 1.11× the optimum): the 34 stationary
quads saved 5.5% of the quads by stopping structure a finer grid resolves. The stop's default is
therefore **off**; the arms are still computed and dumped, for the sweep.

Live marches (`descend_live`, `live_stride = 4`, eight boundaries plus post-horizon rounds):

| target | quads (static) | quads (live) | live/static | error | catch-up share | first split at |
|---|---|---|---|---|---|---|
| near-field | 153 | 153 | 1.00 | 0.00000 | 84.3% | boundary 5, `t = 9.75` |
| deep interior | 353 | 405 | 1.15 | 0.00000 | 90.6% | boundary 0, `t = 1.6` |
| preset_shape_h1 (stationarity on) | 3397 | 3565 | 1.05 | 0.00315 | 95.4% | boundary 0 |
| preset_shape_h1 (stationarity off, the default) | 3585 | 3585 | 1.00 | 0.00021 | 95.5% | boundary 0 |

near-field ends on the **same 153-quad tree** as the static descent, at zero error: nothing
splits until the sixth boundary, one level per boundary to the horizon (37 quads), then
seventeen post-horizon rounds. deep interior grows from the first boundary and ends 15% larger
than the static tree — the running-union price — at the same zero error. The catch-up share is
what a child requested late pays: near-field's structure appears late, so most of its tree was
requested at the horizon. On the sea chart the live march reaches 949 quads at the horizon and
3565 after the post-horizon rounds; with the stop on it ends at 0.3% unresolved against the
static stop-on descent's 2.7%, because a stationary quad is re-tested at every boundary and some
later split — the running union again.

## The `eps` sweep, and a calibration the metric's own tiling was hiding

The tolerance tree on `near-field` is **bitwise the same** at `eps` 0.003, 0.01 and 0.03 (153
quads): its footprints sit far below or far above all three. But at `eps = 0.003` that tree
leaves **5.0%** of pixels unresolved, and a descent cut three times tighter (`tau = 0.001`, 165
quads) still leaves **2.0%**, where the exact optimum reaches 2% at 37 quads. `deep interior` at
`eps = 0.03` is unchanged at 357 quads and zero. The within-footprint spread is a *mean* deviation
from the copies' centroid, halved, so a footprint reading `spread <= tau` can hold a pixel whose
chord to the texel exceeds `tau`; the geometry of copies spread over the cell says by a factor of
about two, and the probe says more. The difference is the metric's rasteriser: `err_sum` assigned
pixel `dx` to sample `dx / tile`, so a sample scored the pixels to its **right**, up to a full
cell away, while `Slice::axis` puts the samples at the corners and the render paints cells
centred on them. The metric and the render tiled differently, and the metric's was the older
one. Fixed to nearest-sample tiling; every `err_sum` moves slightly, `error(full)` stays exactly
zero. **Re-run under the corrected tiling, near-field at `eps = 0.003`: `tau = eps` leaves
0.095% unresolved (was 5.0%), `tau = eps/3` leaves 0.00000 (was 2.0%), 165 quads at 1.14× the
optimum.** At `eps = 0.01` the tolerance tree is unchanged at zero error and the ceiling moved
from 137 to 125 quads. So the descent's `tau` and the metric's `eps` agree to within a factor
of about three, and the metric's tiling was most of what looked like a larger one. The tables
above were taken under the old tiling; the trees are unchanged and the ceilings move by a few
percent.

## The wider chart set: the GLSL presets and `config_stability`

The first four targets were thin -- two tame Burrau regions, one control and one sea. The
standard GLSL presets and the chart the integrator work was measured on are built the same way
(levels 6, `N = 8`, `E+1 = 8`, 512², `t = 13`, `eps = 0.01`; logs in `output/`). The `alpha`
column is the legacy policy's tree from the build log; the tolerance descents against these
caches are in the sections below as they land.

| target | root error | sea fraction at `eps` 2e-3 / 1e-2 / 5e-2 | dp reaches 1% at | uniform reaches 1% at | alpha (legacy) |
|---|---|---|---|---|---|
| `preset_shape` | 0.282 | 0.226 / 0.193 / 0.133 | 1737 | 4145 | 21 quads, error 0.271, `floor:8 keep:8` |
| `preset_prho` | 0.293 | 0.089 / 0.042 / 0.014 | 793 | 3345 | 21, 0.206, `floor:2 keep:14` |
| `preset_plambda` | 0.189 | 0.070 / 0.033 / 0.014 | 557 | 3941 | 21, 0.114, `keep:16` |
| `config_stability` | 0.923 | 0.658 / 0.238 / 0.068 | 4873 | 5325 | 73, 0.763, `floor:5 keep:42 screen_floor:8` |
| `preset_shape_pl` | 0.202 | 0.061 / 0.034 / 0.015 | 557 | 3277 | 21, 0.111, `keep:16` |
| `mid-field` | 0.000 | 0 / 0 / 0 | 1 | 1 | 21, 0.000, `keep:16` (resolved at the root, like `far`) |
| `body2 core` | 0.397 | 0.009 / 0.002 / 0.002 | 113 | 1297 | 21, 0.061, `keep:16` |

Three things the wider set adds. The momentum presets, which are one triangle at different
initial velocities, are resolved by the optimum at a sixth to a seventh of uniform's memory, so
the mechanism has room on them, and `body2 core` at an eleventh (113 against 1297 for 1%);
`mid-field` is a second `far` -- resolved at the root, nothing for any policy to do; `preset_shape` at the corrected window carries a 19% sea and
still shows 2.4x. `config_stability` is the other kind of chart: 92% of its pixels differ from
their level-6 reference at the root, a quarter of the frame is sea at `eps = 0.01` and two
thirds at `eps = 0.002`, and the exact optimum needs 4873 of 5461 quads to reach 1% -- there is
almost nothing a ranking can save, because the structure is everywhere. That is the chart the
area floor is for, and where "do not degenerate into uniform" has to be paid for in error. The
legacy policy stops at the bootstrap on every one of them.

<!-- phase3:begin -->
## Phase 3: the area floor, the agreement arm, merging and the live march, on six charts

Measured 2026-09-04 on six targets -- `near-field`, `deep interior`, `preset_prho`,
`preset_shape`, `config_stability`, `preset_shape_h1` -- at levels 6, `N = 8`, `E+1 = 8`, the
512² camera, `t = 13`, `eps = tau = 0.01`, `k_frac = 0.25`, `refine_flagged = false`, all from
the `.fcache` caches (`payload_metric live` and `payload_metric march`; logs `output/phase3_*.txt`,
batch scripts `output/phase3_*.sh`, tables by `output/sweep_condense.py`). Three pinned binaries:
the first sweep predates `272527c` and `61ff00c`, so its `alpha_lo = 0` rows still floored and are
replaced here by the guarded rows, and its marches use the old expiry and are replaced on the
three charts re-run; the control, guarded and ladder rows are the `61ff00c` build; the
dimension-floor rows are this commit's. **The error columns are disagreement with the level-6
reference at one sample per pixel**, so a complete tree reads near zero on a chart that is a third
sea: the sea's cost shows in the quad count, never in that column. `resolvable left` counts only
pixels a finer grid resolves, so it is the structure the tree gave up.

#### The area floor, static, all six targets

| target | sea | full depth: quads | vs ref | floor 0.2: quads | vs ref | resolvable left | vs optimum | vs uniform | floors | arm off: quads | vs ref | resolvable left | vs uniform | floors |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| near-field | 0.0005 | 153 | 0.0000 | 125 | 0.0003 | 0.0002 | 1.62x | 0.03x | 3 | 153 | 0.0000 | 0.0000 | 0.03x | 0 |
| deep_interior | 0.0025 | 357 | 0.0000 | 225 | 0.0029 | 0.0022 | 3.46x | 0.19x | 8 | 357 | 0.0000 | 0.0000 | 0.07x | 0 |
| preset_prho | 0.0415 | 1061 | 0.0027 | 761 | 0.0396 | 0.0232 | 1.81x | 0.27x | 52 | 1045 | 0.0052 | 0.0029 | 0.30x | 24 |
| preset_shape | 0.1934 | 2225 | 0.0000 | 1645 | 0.0482 | 0.0103 | 1.35x | 0.44x | 262 | 1317 | 0.1025 | 0.0303 | 0.51x | 141 |
| config_stability | 0.2379 | 5089 | 0.0107 | 3421 | 0.2357 | 0.1237 | 1.55x | 1.15x | 229 | 3801 | 0.1967 | 0.1038 | 1.11x | 174 |
| preset_shape_h1 | 0.3415 | 3585 | 0.0001 | 2005 | 0.1894 | 0.0434 | 1.40x | 0.58x | 334 | 2609 | 0.1285 | 0.0295 | 0.66x | 274 |

#### The stationarity stop, static

| target | full depth: fires | vs ref off -> on | floor 0.2: fires | vs ref off -> on |
|---|---|---|---|---|
| near-field | 0 | 0.0000 -> 0.0000 | 0 | 0.0003 -> 0.0003 |
| deep_interior | 0 | 0.0000 -> 0.0000 | 0 | 0.0029 -> 0.0029 |
| preset_prho | 2 | 0.0027 -> 0.0638 | 2 | 0.0396 -> 0.0904 |
| preset_shape | 3 | 0.0000 -> 0.0006 | 6 | 0.0482 -> 0.0494 |
| config_stability | 51 | 0.0107 -> 0.0214 | 25 | 0.2357 -> 0.2413 |
| preset_shape_h1 | 10 | 0.0001 -> 0.0124 | 12 | 0.1894 -> 0.2140 |

#### The live march against the static tree, floor 0.2

| target | static: quads | vs ref | march: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up | pin | arm off march: computed | resident | merged | vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| near-field | 125 | 0.0003 | 125 | 121 | 4 | 0.0005 | 0.0005 | 0.10x | 81% | expiry | 153 | 153 | 0 | 0.0000 |
| deep_interior | 225 | 0.0029 | 293 | 261 | 32 | 0.0044 | 0.0038 | 0.92x | 85% | pre-expiry | 405 | 377 | 28 | 0.0000 |
| preset_prho | 761 | 0.0396 | 841 | 777 | 64 | 0.0404 | 0.0282 | 0.30x | 89% | pre-expiry | 1061 | 1021 | 40 | 0.0047 |
| preset_shape | 1645 | 0.0482 | 1645 | 1361 | 284 | 0.0805 | 0.0192 | 0.48x | 93% | pre-expiry | 1545 | 1369 | 176 | 0.0977 |
| config_stability | 3421 | 0.2357 | 3865 | 3649 | 216 | 0.1909 | 0.1055 | 1.12x | 95% | expiry | 4233 | 4073 | 160 | 0.1495 |
| preset_shape_h1 | 2005 | 0.1894 | 2105 | 1621 | 484 | 0.2351 | 0.0597 | 0.96x | 93% | expiry | 2649 | 2221 | 428 | 0.1766 |

#### The alpha_lo ladder, static, arm on

| target | alpha_lo | quads | vs ref | resolvable left | vs optimum | vs uniform | floors | quads saved |
|---|---|---|---|---|---|---|---|---|
| near-field | 0 | 153 | 0.0000 | 0.0000 | 1.22x | 0.03x | 0 | 0% |
| near-field | 0.005 | 125 | 0.0003 | 0.0002 | 1.62x | 0.03x | 3 | 18% |
| near-field | 0.2 | 125 | 0.0003 | 0.0002 | 1.62x | 0.03x | 3 | 18% |
| deep_interior | 0 | 357 | 0.0000 | 0.0000 | 1.10x | 0.07x | 0 | 0% |
| deep_interior | 0.005 | 225 | 0.0029 | 0.0022 | 3.46x | 0.19x | 8 | 37% |
| deep_interior | 0.2 | 225 | 0.0029 | 0.0022 | 3.46x | 0.19x | 8 | 37% |
| preset_prho | 0 | 1061 | 0.0027 | 0.0027 | 1.23x | 0.30x | 0 | 0% |
| preset_prho | 0.005 | 977 | 0.0103 | 0.0085 | 1.47x | 0.29x | 25 | 8% |
| preset_prho | 0.2 | 761 | 0.0396 | 0.0232 | 1.81x | 0.27x | 52 | 28% |
| preset_shape | 0 | 2225 | 0.0000 | 0.0000 | 1.02x | 0.52x | 0 | 0% |
| preset_shape | 0.001 | 1837 | 0.0322 | 0.0068 | 1.33x | 0.47x | 160 | 17% |
| preset_shape | 0.005 | 1837 | 0.0322 | 0.0068 | 1.33x | 0.47x | 168 | 17% |
| preset_shape | 0.02 | 1837 | 0.0322 | 0.0068 | 1.33x | 0.47x | 172 | 17% |
| preset_shape | 0.05 | 1801 | 0.0350 | 0.0075 | 1.33x | 0.46x | 189 | 19% |
| preset_shape | 0.1 | 1785 | 0.0358 | 0.0078 | 1.32x | 0.46x | 214 | 20% |
| preset_shape | 0.2 | 1645 | 0.0482 | 0.0103 | 1.35x | 0.44x | 262 | 26% |
| preset_shape | 0.3 | 1569 | 0.0558 | 0.0123 | 1.37x | 0.43x | 327 | 29% |
| config_stability | 0 | 5089 | 0.0107 | 0.0107 | 1.06x | 0.96x | 0 | 0% |
| config_stability | 0.001 | 4729 | 0.0582 | 0.0259 | 1.22x | 1.00x | 206 | 7% |
| config_stability | 0.005 | 4729 | 0.0582 | 0.0259 | 1.22x | 1.00x | 211 | 7% |
| config_stability | 0.02 | 4637 | 0.0709 | 0.0302 | 1.25x | 1.01x | 214 | 9% |
| config_stability | 0.05 | 4445 | 0.0988 | 0.0413 | 1.31x | 1.03x | 220 | 13% |
| config_stability | 0.1 | 4165 | 0.1374 | 0.0663 | 1.39x | 1.05x | 208 | 18% |
| config_stability | 0.2 | 3421 | 0.2357 | 0.1237 | 1.55x | 1.15x | 229 | 33% |
| config_stability | 0.3 | 3149 | 0.2652 | 0.1407 | 1.58x | 1.22x | 256 | 38% |
| preset_shape_h1 | 0 | 3585 | 0.0001 | 0.0001 | 1.04x | 0.66x | 0 | 0% |
| preset_shape_h1 | 0.001 | 2185 | 0.1628 | 0.0350 | 1.35x | 0.59x | 207 | 39% |
| preset_shape_h1 | 0.005 | 2169 | 0.1654 | 0.0353 | 1.36x | 0.59x | 203 | 39% |
| preset_shape_h1 | 0.02 | 2137 | 0.1705 | 0.0375 | 1.37x | 0.59x | 226 | 40% |
| preset_shape_h1 | 0.05 | 2089 | 0.1779 | 0.0394 | 1.38x | 0.59x | 240 | 42% |
| preset_shape_h1 | 0.1 | 2041 | 0.1849 | 0.0418 | 1.40x | 0.58x | 283 | 43% |
| preset_shape_h1 | 0.2 | 2005 | 0.1894 | 0.0434 | 1.40x | 0.58x | 334 | 44% |
| preset_shape_h1 | 0.3 | 1697 | 0.2213 | 0.0562 | 1.39x | 0.66x | 322 | 53% |

#### The live march at alpha_lo 0.005 against 0.2

| target | static 0.005: quads | vs ref | march 0.005: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up | march 0.2: computed | resident | merged | vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| preset_shape_h1 | 2169 | 0.1654 | 2225 | 1925 | 300 | 0.1936 | 0.0478 | 0.66x | 93% | 2105 | 1621 | 484 | 0.2351 |
| config_stability | 4729 | 0.0582 | 4813 | 4629 | 184 | 0.0638 | 0.0292 | 1.03x | 96% | 3865 | 3649 | 216 | 0.1909 |
| preset_shape | 1837 | 0.0322 | 1725 | 1517 | 208 | 0.0698 | 0.0158 | 0.48x | 92% | 1645 | 1361 | 284 | 0.0805 |

#### The dimension floor off, noise stop on, alpha_lo 0.2, arm on

| target | sea | full depth: quads | floor on: quads | vs ref | resolvable left | vs uniform | floors | dim off: quads | vs ref | resolvable left | vs optimum | vs uniform | floors | quads saved |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| near-field | 0.0005 | 153 | 125 | 0.0003 | 0.0002 | 0.03x | 3 | 153 | 0.0000 | 0.0000 | 1.22x | 0.03x | 2 | 0% |
| deep_interior | 0.0025 | 357 | 225 | 0.0029 | 0.0022 | 0.19x | 8 | 349 | 0.0000 | 0.0000 | 1.07x | 0.06x | 3 | 2% |
| preset_prho | 0.0415 | 1061 | 761 | 0.0396 | 0.0232 | 0.27x | 52 | 1045 | 0.0030 | 0.0029 | 1.23x | 0.30x | 11 | 2% |
| preset_shape | 0.1934 | 2225 | 1645 | 0.0482 | 0.0103 | 0.44x | 262 | 2181 | 0.0004 | 0.0001 | 1.06x | 0.51x | 13 | 2% |
| config_stability | 0.2379 | 5089 | 3421 | 0.2357 | 0.1237 | 1.15x | 229 | 5045 | 0.0149 | 0.0114 | 1.08x | 0.96x | 37 | 1% |
| preset_shape_h1 | 0.3415 | 3585 | 2005 | 0.1894 | 0.0434 | 0.58x | 334 | 3353 | 0.0300 | 0.0067 | 1.24x | 0.67x | 9 | 6% |

#### The live march with the dimension floor off

| target | static dim off: quads | vs ref | march: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up |
|---|---|---|---|---|---|---|---|---|---|
| preset_shape_h1 | 3353 | 0.0300 | 3353 | 3353 | 0 | 0.0300 | 0.0067 | 0.67x | 94% |
| config_stability | 5045 | 0.0149 | 5065 | 5065 | 0 | 0.0114 | 0.0113 | 0.96x | 96% |
| preset_shape | 2181 | 0.0004 | 2181 | 2181 | 0 | 0.0004 | 0.0001 | 0.51x | 91% |

#### The live march at 0.005 after the cap re-decision fix

| target | static 0.005: quads | vs ref | march before fix: computed | resident | merged | vs ref | march after fix: computed | resident | merged | vs ref | resolvable left | vs uniform | catch-up |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| preset_shape_h1 | 2169 | 0.1654 | 2225 | 1925 | 300 | 0.1936 | 2225 | 1909 | 316 | 0.1961 | 0.0483 | 0.67x | 92% |
| config_stability | 4729 | 0.0582 | 4813 | 4629 | 184 | 0.0638 | 4813 | 4601 | 212 | 0.0689 | 0.0293 | 1.04x | 95% |
| preset_shape | 1837 | 0.0322 | 1725 | 1517 | 208 | 0.0698 | 1725 | 1513 | 212 | 0.0702 | 0.0159 | 0.48x | 91% |

### What it says

1. **The area floor saves quads on every chart and floors resolvable structure on every chart.**
   Against the clean full-depth tree it saves 18% (`near-field`) to 44% (`preset_shape_h1`) and
   stays under uniform on five of six. What it gives up runs from two hundredths of a percent of
   the frame to **11% on `config_stability`**, where the floored tree is bigger than a uniform tree
   at its own error (1.15x). Even where the cost is negligible it is not zero: with the arm off
   `near-field` and `deep interior` never floor at all, so every floor there was the noise stop
   firing on structure.
2. **The threshold is not the lever.** Across `alpha_lo` 0.05-0.3 the sea charts give the same
   trade within a few percent and `config_stability` floors 208-256 boxes at every rung and never
   gets under uniform. The boxes it floors have an exponent under 0.05 at the coarse levels: their
   structured fill falls by under 7% per level, a set of box dimension about **1.94**. Under the
   floor's definition (`alpha_area = 2 - d`) that is noise; under the brief's words it is the
   fractal structure the mechanism exists for, and 16% of its area resolves at level 6. The
   saturation account -- the grandparent's four-by-four coarse end saturates, the exponent
   becomes minus half the log of the children's fill, and 0.2 floors any box over three quarters
   unresolved -- is correct arithmetic and **was refuted as the explanation**: at 0.05 the same
   boxes floor.
3. **Real seas are coherent, so the noise stop is nearly inert on real charts.** The agreement
   arm was calibrated on a synthetic white sea and decided the synthetic tests; on the six charts
   its noise stop floors **2-37 quads** where the dimension floor floors 262-334. Footprints whose
   copies diverge still agree with their neighbours on class and nominal shape. Same finding as
   the stationarity arm's coherent sponge, at a second arm. The arm's real work is in the
   exponent's weights: on `preset_shape` it takes the tree from 0.51x to **0.44x** uniform and the
   resolvable loss from 3% to 1%, because the sea drops out of both ends and the exponent reads
   the structure's own scaling; on `preset_prho` it costs 2% of the resolvable pixels for 27% of
   the quads; on the two seas it moves along the trade curve.
4. **The dimension floor off is the degeneration the brief forbids, and the fine ladder moves
   the default instead.** With only the noise stop, `preset_shape_h1` runs to 3353 of 3585 quads,
   `config_stability` to 5045 of 5089, `preset_shape` to 2181 of 2225 -- uniform depth on every
   sea at no resolvable loss -- so the floor stays on and `SchedCfg::dim_floor = false` is the
   named opt-out beside `alpha_lo = 0`. The fine rungs then say that exact saturation does not
   separate a sea from a sponge either: at `alpha_lo = 0.001` `config_stability` still floors 206
   boxes with 2.6% of the frame's resolvable pixels inside them, because a sponge that thins only
   below the coarse sampling scale is saturated seen from above, and no statistic of the levels
   computed can see the levels not computed. What the rung changes is the trade. **The default
   is now `alpha_lo = 0.005`**: a split is floored only where its children resolved nothing to
   within a noise margin, which is what "no gain" says. Against 0.2 it keeps the sea chart's
   saving (39% against 44%, cost 3.5% against 4.3%), halves `preset_prho`'s and `preset_shape`'s
   costs for a third less saving, and takes `config_stability` from 15% worse than uniform to
   break-even (1.00x, cost 2.6% against 12.4%). The marches at 0.005 merge less and trail their
   static trees less on all three charts (sea chart 0.1936 against 0.2351). 0.2 is the
   dimension-threshold rung, recorded, and `alpha_lo` is still that threshold for anyone who
   wants it: `alpha_area = 2 - d`, so 0.2 floors where the unresolved set is fatter than
   `d = 1.8`.
5. **Stationarity is worse on every chart where it fires.** Never on `near-field` or
   `deep interior`; on `preset_prho` two coarse stops take the full-depth tree from 0.3% to 6.4%
   unresolved and 3x uniform on the resolvable arm; `config_stability` doubles its error on 51.
   Off by default, now on six charts rather than one.
6. **The live march trails the static tree exactly through the no-gain merges.** With the
   dimension floor off the march reproduces the static tree quad for quad on all three charts run
   -- zero merges, resident equal to computed. With it on, the sea chart's march computes more
   quads than the static tree, holds a fifth fewer resident and displays worse (0.2351 against
   0.1894): a no-gain merge taken at an early boundary, on footprints that had not yet developed
   the structure, survives its appearance. The structured-weight expiry is **inert** on
   `near-field` and `config_stability` (rows identical to the old pin) and moves the sea chart
   from 0.2686 to 0.2351. Catch-up is 81-96% of the work on every chart. Two live-compatible
   repairs, unbuilt: a time-to-live on the memory, or an expiry keyed on the exponent's own
   inputs.
7. **A capped leaf was terminal in the live frontier, so a resolved parent could never merge it.**
   Found by the pulse test the hour the margin tightened: under 0.005 the live tree ended at 149
   quads against the static 69, finer, with 24 parents reading zero unresolved footprints at the
   horizon and still holding 96 capped children. A leaf stopped by a cap left the frontier for
   good, so it was never re-decided and kept its cap label after its region resolved, and the
   merge pass took only kept, floored, stationary or deferred children as settled. At 0.2 those
   parents had merged for no gain while the band was wide, which hid it. Caps are now re-tested
   every boundary and read `Keep` once the region resolves, any leaf that did not split this
   round is settled, and the pulse reads 69 against 69 with 128 merged. On the three real charts
   the route exists and gives back little at `t = 13`: 16 to 28 more children merged and the error
   within half a point (last table above), because their screen-floored regions have not resolved
   by the horizon. The fix is a correctness property of the merge, measured on the field that
   exercises it; on these charts at this horizon it is nearly inert, and that is recorded rather
   than inferred from the pulse.
8. **A build that prints nothing built nothing.** `cargo` dropped off the shell's PATH mid-session;
   two builds silently did not run, a bitwise pin check compared two stale binaries and read
   "identical", and a batch ran the old code under a new name. Caught because the third build's
   filtered output was empty where it should have carried a `Finished` line. Every build in a
   batch now prints that line unfiltered, and the binary is checked for a string the change adds.

<!-- phase3:end -->

## What the numbers say

1. **The room was an order of magnitude, not 3–12%.** On `near-field` and `deep interior` the
   tolerance descent reaches zero unresolved pixels within 1.1× of the exact optimum and at 3–7%
   of uniform's memory. The legacy policy stops at the bootstrap on both, which is why the
   refinement mechanism looked scatterbrained: the criterion had been inverted (split only where
   halving the cell had already halved the spread) and the metric it was scored under could not
   tell resolved from unresolved.
2. **`deep interior` is not "structure everywhere".** `sea_fraction(0.01) = 0.0025`. The standing
   reading was the OKLab metric's, which reads 0.379 at the root there because it stretches the
   spread field's texture to full contrast.
3. **`preset_shape_h1` is the sea, and the stationarity stop is worse than off on it.** It fires
   on 34 of 2064 floor leaves, and where it fires it stops structure a finer grid resolves. The
   real mixing region is a coherent sponge at the footprint scale — arcs, bands, perforation —
   and coherence reads it as structure, which it is. The stop is built for a field white at the
   footprint scale, which the synthetic sea is and this chart is not. The memory there is what
   the field costs at `eps = 0.01`: **1.03× the optimum, 0.66× uniform** with the stop off. The
   lever is `eps` — `sea_fraction` is 0.34 at 0.01 and 0.15 at 0.05 — and that is the Phase 3
   sweep. The user chose the stop as in scope from the start; it is built, tested and measured,
   and the measurement puts its default at off.
4. **The signal rankings under the payload metric line up with the ceiling** where the field is
   tame (five of six reach 1% at `B = 9` on near-field), and on the sea chart `within/median`, the
   old default, is the best single ranking (0.027 at `B = 3071` against dp 0.010 and uniform
   0.206), where under OKLab it was the worst. Same rankings, different metric, opposite verdicts.

## Calibration facts found on the way, each by a test that could fire

- The quadrant mixture arm is a **mean** over four quadrants: the max of four 16-footprint
  multinomial deviations reaches 0.33 by sampling noise on three classes, and a pure synthetic sea
  split on its own noise.
- The class coherence is **class-conditional**, the max over classes present: a global
  agreement-above-chance statistic moves by a few percent for one coherent column of eight
  footprints, while that column's own class clusters at 0.6 against a base rate of an eighth.
- A filament in a quad's **edge column** has the same mixture as its parent (an eighth of both),
  so the two-scale arm is blind to it; only class-conditional coherence catches it, and only when
  the filament's class differs from the sea's. Stated as a limit of `N = 8`.
- Under `Policy::Tolerance` a resolved or stationary quad is decided **ahead of the caps**, so
  `MaxLevel` and `ScreenFloor` count only quads that wanted to split.
- The live descent gains at most **one level per boundary** (children must catch up to the
  playhead before they can be decided) and continues in **post-horizon rounds** at the last
  boundary; without those it stopped with leaves still pending.
