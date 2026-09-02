# The ribbon oscillation: what it is

**Question.** The ribbons in `config_stability` show alternating dark/light bands. Physics, or a
numerical artefact of the stepper, the substep count, the sampler, or the raster?

**Answer.** Physics. The bands are the **bound pair's orbital phase winding through
initial-condition space**, and the window measured is a *regular island* — neighbouring
trajectories do not diverge there. The moiré that appears at deep zoom is that structure crossing
the pixel grid, not a defect.

Harnesses: `examples/osc_zoom.rs`, `osc_source.rs`, `osc_ss.rs`, `band_phase.rs`,
`band_survey.rs`. Fields and panels in this directory. Predictions recorded before each run are
in `PREDICTION.md`, including the ones that were refuted.

---

## 1. Four candidates, three killed by measurement

All at the `config_stability` window `zf = 0.060` (`z1`), 256², `t_max = 50`, `n_sync = 125`,
`refine_flagged` off.

### Not the stepper

`eta` 1.0e-2 → 2.5e-3, **3.07× the steps** (`steps p50` 1.435e4 → 4.408e4). High-passed
correlation against the baseline, with a half-frame shifted control:

| σ | corr | shifted |
|---|---|---|
| 2 | 0.9585 | −0.0365 |
| 4 | 0.9777 | −0.1056 |
| 8 | 0.9852 | −0.1002 |
| 16 | 0.9904 | −0.1330 |

**The arm is not inert** — 56.3% of pixels moved, worst 10.6/255. The step size rewrites the
render and leaves the bands where they were. λ identical at 36.6 px under a row spectrum and
12.64 px under the 2D one.

### Not the substep count

`total_substeps` is banded — the sharpest beat in the whole table, prominence 152858 — which
made it the leading suspect. But profiled along the band normal it is **not a staircase**: 58–67%
of increments are exactly zero and the risers spread 2–29, and the change across one band is
**22 steps summed over 8 copies**, ~2.8 per trajectory. One extra revolution costs hundreds. The
step count is banded *because the trajectories are*, not the reverse.

### Not the pixel grid, and not the sampler

Pixel-locked test — `z1` against `z2` at 1:1, same grid, different chart region:
**−0.03 to −0.06**, a flat null.

Sampler test — the stencil changed completely, at the `z1` window:

| arm | offsets | E+1 | extent | λ px | angle | prom | vs E8 | shifted |
|---|---|---|---|---|---|---|---|---|
| `halton_E8` | Halton | 8 | 0.5 | **12.64** | **69.8°** | 10837 | 1.0000 | 0.0010 |
| `pcg_E8` | PCG, per-pixel | 8 | 0.5 | **12.64** | **69.8°** | 1958 | 0.6759 | 0.0013 |
| `halton_E4` | Halton | 4 | 0.5 | **12.64** | **69.8°** | 2013 | 0.7905 | 0.0006 |
| `halton_E16` | Halton | 16 | 0.5 | **12.64** | **69.8°** | 21922 | 0.9328 | 0.0012 |
| `jitter_q` | Halton | 8 | 0.25 | **12.64** | **69.8°** | 1512 | 0.6440 | 0.0007 |
| `jitter_f` | Halton | 8 | 1.0 | **12.64** | **69.8°** | 30817 | 0.5523 | 0.0014 |

λ and angle do not move at all. **The positive evidence is stronger than the null**: prominence
rises monotonically with sample count (2013 → 10837 → 21922 at E+1 = 4, 8, 16) and with extent
(1512 → 10837 → 30817). A sampling artefact weakens with more samples; this sharpens.

And it is in the **nominal** fields, where nothing is sampled at all — seven single-trajectory
observables, one orientation:

| field | ensemble? | λ px | angle | prom |
|---|---|---|---|---|
| `halton_E8_spread` | yes | 12.64 | 69.8° | 10837 |
| `nominal_shape0` | no | 10.59 | 65.6° | 2382 |
| `nominal_shape1` | no | 12.64 | 69.8° | 4451 |
| `nominal_shape2` | no | 14.20 | 70.6° | 38916 |
| `nominal_tend` | no | 11.39 | 69.1° | 14181 |
| `nominal_dmin` | no | 29.96 | 69.4° | 10007 |
| `nominal_drift` | no | 7.49 | 69.4° | 9569 |
| `nominal_steps` | no | 11.39 | 69.1° | 152858 |

### Not sync-cadence quantisation, and not a `t_end` staircase

`t_end` lands exactly on a sync boundary on 0.3048 of the frame — which reads alarming until the
horizon population is separated: **0.3036** of the frame sits at `t_end = t_max = 50`, itself a
multiple of `dt = 0.4`. Genuinely quantised: **0.0012**. This is the standing *"on a boundary
conflates quantised with finished"* result. Across the bands `t_end` runs 4.2847–4.5625,
continuous, 246 distinct values in 440 samples.

---

## 2. What it is

`examples/band_phase.rs`. Twenty-four ICs spanning exactly one band, integrated with
`n_sync = 20000` (`dt = 2.5e-3`), termination **off** so the series runs the full horizon.

**The trajectories are the same trajectory at different orbital phase.** Best-lag correlation
against sample 0:

| sample | lag (time) | corr |
|---|---|---|
| 0 | 0 | 1.0000 |
| 6 | −0.0500 | 1.0000 |
| 12 | −0.0975 | 1.0000 |
| 18 | −0.1475 | 0.9999 |
| 23 | −0.1900 | 0.9999 |

**The lag grows linearly in time, not exponentially.** Measured between the two ends of one band:

| window | lag | corr | × prev |
|---|---|---|---|
| t 5–10 | 2.500e-2 | 1.0000 | — |
| t 10–15 | 5.000e-2 | 0.9999 | 2.00 |
| t 15–20 | 7.500e-2 | 0.9999 | 1.50 |
| t 20–25 | 9.750e-2 | 0.9999 | 1.30 |
| t 25–30 | 1.200e-1 | 0.9999 | 1.23 |
| t 30–35 | 1.425e-1 | 0.9999 | 1.19 |
| t 35–40 | 1.650e-1 | 0.9999 | 1.16 |
| t 40–45 | 1.900e-1 | 0.9999 | 1.15 |

The ratios are exactly 2/1, 3/2, 4/3, 5/4, … — a constant increment. **This is a frequency beat,
not a divergence.**

**The mechanism.** Every IC here is a bound pair plus a third body. Neighbouring ICs give pairs
with slightly different periods; over the horizon the phase difference accumulates linearly; where
it completes one cycle every observable returns to where it was. That is one band.

**Arithmetic, from two independent directions.** The pair period is `T = 2.5` late in the run and
`3.125` early — it **hardens** as energy goes to the escaper. The shape signal runs at *twice* the
orbital frequency (the triangle repeats twice per orbit), so one shape cycle is `T = 1.25`. One
full cycle of accumulated phase is `12.64 × 1.25 / 0.198 =` **79.8 px**, against the raw
unfiltered field's coarse tier at **80.95 px**. A 1.6% match.

**The cascade is the harmonic ladder.** The shape waveform is a pericentre spike, not a sinusoid.
Measured on the true fundamental 0.400: 1.0, 0.033, 5.6e-3, 1.7e-3, 6e-4 — about 3× per harmonic.
Spatial tiers, all at one orientation:

| tier | λ px | ratio to coarse |
|---|---|---|
| raw field | 80.95 | 1 — one shape cycle |
| high-passed | 12.64 | 6.40 |
| tightest high-pass | 2.24 | 36.1 |

The visible striations sit near the sixth harmonic; the high-pass that made them visible is what
removed the fundamental above them. **Not pinned:** that 6.40 is exactly 6.

---

## 3. This window is a regular island

A correlation of **0.9999** between ICs a whole band apart, held to `t = 45`, means neighbouring
trajectories are **not diverging** here. The banding is the clean quasi-periodic structure such a
region has — the phenomenon Trani et al. 2024 (arXiv:2403.03247) call *"isles of regularity in a
sea of chaos"*.

Consequence for this codebase: **`spread_shape` here is not measuring chaotic divergence, it is
measuring a phase beat.** The refinement criterion reads it as structure worth resolving, and it
is structure — just not the kind the criterion was designed around.

---

## 4. The deep-zoom moiré, and what fixes it

At `z4` the dominant spatial peak is at **λ 2.44 px, angle 106°** — pinned at the Nyquist floor
with a folded angle, where every resolved tier reads 64–79°. That is what structure above Nyquist
looks like: not the structure, its reflection.

The magnification chain confirms it: `z1@2× vs z2` registers **0.78–0.88** against shifted
controls of −0.14 to −0.24, while `z2@2× vs z4` breaks to **−0.11**. Reading each band separately
rather than letting the loudest win shows why — the tiers *do* magnify (z1's 12.64 → z4's 62.09,
z1's 2.24 → z4's 7.86, both within a bin of 4×), but at z4 a new finer tier at the sampling limit
outranks the predicted one **54×**.

**Field supersampling reduces it and cannot remove it.** Marching a finer grid and box-averaging:

| f | samples/px | L hf rms | vs f=1 |
|---|---|---|---|
| 1 | 1 | 0.08514 | 1.000 |
| 2 | 4 | 0.07432 | 0.873 |
| 3 | 9 | 0.04810 | 0.565 |

Nine samples per pixel takes it to 56.5%; the broad coherent fringes break up and finer texture
replaces them. A harmonic ladder has a rung below whatever rate you sample at.

**`rgb_resolved` cannot help, and the reason is the window not the code.** With
`keep_copy_shapes` on, copies are kept on 36864 of 36864 pixels and the resolved render differs
from the plain one on **1 pixel** at f=1 and **0** at f=2 and f=3. It supersamples *hue*; the
fringes are in *lightness*, and `spread_shape` is a footprint statistic with no per-copy analogue,
so `l` is computed once and held constant across sub-samples. Separately, there is nothing for it
to average at that depth: OKLab `a` spans **0.00292 across the whole z4 frame** against ~0.0015
per 8-bit step, so a single cell spans ~1e-5, three orders below the quantiser. At `z1` the range
is 0.20071, 69× wider — and the control there confirms the mechanism is healthy: the resolved
render differs from the plain one on **1878 of 36864 pixels (5.09%)** against **1** at `z4`.

**`Scheme::Pcg` looks like a fix and is not.** Per-pixel random offsets convert the coherent
fringe into incoherent noise of the same amplitude. This project's own rule: *a remedy that only
changes the spatial correlation of an error is cosmetic.* It would also cost the fixed-stencil
property the refinement criterion rests on.

---

## 5. Instrument faults found on the way

Recorded because each was caught by measurement rather than by reading code.

- **The first registration number was 0.9844 and meant nothing** — it was the ribbon correlating
  with itself. High-passing is what made the statistic answer the question asked.
- **The turning-point counter read `px/half 2.7`** — Nyquist. It was counting single-pixel noise
  and could not see a 40-px band at all.
- **The shared colour window starved `z4`**: range 0.056, ~14 of 255 levels, registration 0.0256
  against a shifted control of 0.0693. A null meaning *this side resolves nothing*. Fixed by an
  explicit per-panel window arm — correct for period and registration, wrong for amplitude.
- **The `ssres` arm was inert and produced bitwise-identical files.** `keep_copy_shapes` is
  `false` in production so `rgb_resolved` returned `rgb` through its own guard. The harness now
  asserts the arms differ and prints the count.
- **The chroma channels appeared 7× noisier under supersampling.** Not 8-bit recovery — the
  distinct level counts *fall* (117 → 94 → 76) while the range grows. Cause: the downsampler
  averaged in **linear sRGB**, and two pixels of one hue at the fringes' lightness extremes land
  off that hue by `|Δa|` 0.0034–0.0112, matching the measured 0.011. The shipped
  `compose::resolve` averages in OKLab and is unaffected.
- **`band_phase` printed `t_end` span 0.0000** — termination was off, so `t_end` is the horizon on
  every sample and that arm was vacuous by construction.
- **A `jitter_frac = 0` guard read FAIL on `== 0.0`.** A variance taken as a difference of moments
  does not cancel bitwise on identical inputs; it reads `1.272e-16`, `2.14e10×` below the field it
  guards. The bar was wrong, not the code.

## 6. Predictions that were refuted

Recorded before their runs, in `PREDICTION.md`.

- **Sampler beat.** All four of its predictions failed: Pcg did not kill the coherence, E+1 of 4
  and 16 did not move the period, the sampling extent did not move it, and the nominal fields
  carried it.
- **z4 tier position.** Predicted the peak near λ 9 px at 63–72°; it is 2.48 px at 106°. The
  predicted tier exists at 7.86 px and is outranked 54× — right that it would be there, wrong
  that it would be the peak.
- **One-more-revolution at `d_min`.** Predicted `d(t_end)` per band ≈ the orbital period at
  closest approach: 3.104e-2 against 1.933e-4, a ratio of 160.6. Wrong orbit — `d_min` is the
  instant of collision, not the pair between encounters.

---

## 7. Survey: is the phase beat local, and does colour predict it?

`examples/band_survey.rs`, 96 probes on interior grids across six charts. Two ICs one coarse
pixel apart, integrated to `t = 50` with `n_sync = 20000`, termination off; the lag between their
shape series is measured early (`t ≈ 7`) and late (`t ≈ 42`).

The discriminator needs no spectrum:

| | lag growth | corr at `t ≈ 42` |
|---|---|---|
| **regular** — phase beat | linear in `t` | ~1.0 |
| **chaotic** | exponential | collapses |

**`corr` is the honest arm and the lag is the fragile one.** Once two trajectories decorrelate,
the "lag" is the argmax of noise — a number, always, and meaningless. Decorrelated rows print
`--` rather than a value.

### The guard that had to be added first

`config_basin` returned lag **exactly 0.0000** and correlation **exactly 1.0000** on all 16
probes, and several `preset_prho`/`preset_plambda` rows did the same. Perfect agreement is what a
regular island looks like, AND what two identical inputs look like, AND what two frozen outputs
look like — three states, one number, and `band_survey` could not tell them apart.

`examples/band_guard.rs` reports `|dr|`, `|dv|` and the late-window signal amplitude per probe.
**All 64 checked probes are live.** It also caught the keying trap already on record:
`preset_prho` and `preset_plambda` read `|dr|` **exactly 0** with `|dv|` of 8e-2 to 2e-1 — the
constant-configuration slices behaving as documented, where a position-keyed check would have
reported a collapsed decode.

### By chart

| case | n | regular | decorrelated | corr_l p50 | `sd_late` p50 |
|---|---|---|---|---|---|
| `config_stability` | 16 | 4 | 10 | 0.7990 | 8.1e-3 |
| `config_basin` | 16 | **16** | 0 | 1.0000 | 4.1e-2 |
| `preset_shape` | 16 | **0** | **16** | **0.1448** | — |
| `preset_prho` | 16 | **16** | 0 | 1.0000 | 5.9e-5 |
| `preset_plambda` | 16 | 14 | 2 | 1.0000 | 2.4e-5 |
| `preset_shape_pl` | 16 | 10 | 6 | 0.9966 | — |

**The split reproduces a standing result and supplies its mechanism.** The record says *a chart's
tameness is set by WHICH COORDINATES it varies, not by where it is centred* — `preset_shape` and
`preset_prho` share `z0 = 0` exactly and differ 5.7× in leaf count, because only the first sweeps
the configuration coordinates. That was measured through `alpha` interdecile and leaf counts.
Here it is measured on the trajectories directly: `preset_shape` is **chaotic at every probe**
(corr 0.1448), `preset_prho` **regular at every probe** (corr 1.0000).

**`config_basin` has no beat at all**, not merely a slow one: lag is exactly 0.0000 everywhere on
a live signal (`|dr|` 6.3e-4, `sd_late` 4.1e-2). Its window is `zoom = 0.009095`, **70× tighter**
than `config_stability`'s 0.63763, so the pair period does not vary measurably across it.
**Prediction, untested: `config_basin` should show no ribbon banding.**

**Amplitude caveat on the momentum presets.** `preset_prho` and `preset_plambda` are regular on a
signal of `sd_late` 2.4e-5 to 5.9e-5 — three orders below `config_basin`'s 4.1e-2. Real (six
orders above f64 round-off) but small: those slices are close to frozen, escape-dominated. Their
"regular" verdict is sound and their *amplitude* is not comparable to the others'.

### By colour, within `config_stability`

The question was whether the blue and green ribbons behave like the red one. Colour is
`hue_ab(shape_vec)` — the shape-sphere position, i.e. which configuration the system ended in.

| colour | n | regular | decorrelated | corr_l p50 | states |
|---|---|---|---|---|---|
| red/orange | 7 | 3 | 3 | 0.9665 | bounded |
| **blue** | 6 | **0** | **5** | **0.6539** | bounded, collision |
| green | 3 | 1 | 2 | 0.7273 | bounded, collision |

**Colour does not determine regularity, but it is not independent of it either.** Red/orange
splits 3–3, so hue alone decides nothing. Blue is **0 of 6 regular** and green 1 of 3 — the blue
ribbons look chaotic. **The sample is thin** (6 blue, 3 green) and this is the one claim here that
a larger probe grid could overturn; it is stated as a lead, not a result.

### And the deep-dive window was a regular patch of a mostly-chaotic slice

`config_stability` at the full frame is **4 regular against 10 decorrelated**. The window every
earlier section measured sits at frame fraction (0.10, 0.45), beside probe (0.12, 0.38), which
reads corr_l **0.9973** — regular. So §2's mechanism is established for a regular sub-region and
**is not a claim about the whole slice**. Most of `config_stability` is chaotic, and there the
`spread_shape` field is measuring divergence, exactly as the criterion assumes.

**That is the useful form of the finding:** one slice contains both regimes, the criterion cannot
tell them apart, and the discriminator that can is two trajectories and a correlation.

---

## Reproduction

Run from the repo root. These are the invocations that produced the committed artefacts, not a
reconstruction of them — a documented reproduction command has been wrong in this project before,
and only running it found out.

```sh
cargo run --release --example osc_zoom   -- 256 results/osc 0.10 0.45          # z*.png
cargo run --release --example osc_zoom   -- 256 results/osc 0.10 0.45 own      # z*_own.png
cargo run --release --example osc_source -- 256 results/osc/fields    0.060    # z1 fields
cargo run --release --example osc_source -- 256 results/osc/fields_z4 0.015    # z4 fields
cargo run --release --example osc_ss     -- 192 results/osc    0.015 1,2,3     # ss_*.png
cargo run --release --example osc_ss     -- 192 results/osc_z1 0.060 1         # z1 control
cargo run --release --example band_phase -- results/osc/phase 0.060 20000 24 12.64
cargo run --release --example band_survey -- results/osc 4                     # survey.txt
cargo run --release --example band_guard  -- 4                                 # guard.txt
```

**The raw `.f64` field dumps are NOT committed** — `fields/`, `fields_z4/` and `phase/` come to
25 MB and are regenerated exactly by the three commands above. The PNGs, the tables and this
record are committed. Said here rather than left as a silent omission.
