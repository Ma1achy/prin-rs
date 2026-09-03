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

Three things the wider set adds. The momentum presets, which are one triangle at different
initial velocities, are resolved by the optimum at a sixth to a seventh of uniform's memory, so
the mechanism has room on them; `preset_shape` at the corrected window carries a 19% sea and
still shows 2.4x. `config_stability` is the other kind of chart: 92% of its pixels differ from
their level-6 reference at the root, a quarter of the frame is sea at `eps = 0.01` and two
thirds at `eps = 0.002`, and the exact optimum needs 4873 of 5461 quads to reach 1% -- there is
almost nothing a ranking can save, because the structure is everywhere. That is the chart the
area floor is for, and where "do not degenerate into uniform" has to be paid for in error. The
legacy policy stops at the bootstrap on every one of them.

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
