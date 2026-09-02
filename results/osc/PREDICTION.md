# Prediction, recorded 2026-09-02 before `osc_source` returned

## What is established

The ribbon oscillation is **not** a pixel-grid beat and **not** integer step-count level sets.
High-passed registration, shifted controls on every row:

| σ | z1@2× vs z2 (chart) | shifted | z1 vs z2 1:1 (pixel) | z1 vs z1_eta4 (step) | shifted |
|---|---|---|---|---|---|
| 4 | 0.8184 | −0.1464 | −0.0631 | 0.9777 | −0.1056 |
| 8 | 0.8713 | −0.1807 | −0.0205 | 0.9852 | −0.1002 |
| 16 | 0.8949 | −0.2452 | −0.0551 | 0.9904 | −0.1330 |

The step arm is **not inert**: `eta`/4 (3.07× the steps) moved 56.3% of pixels, worst 10.6/255.

## What is open, and where it actually lives

`z2@2× vs z4` reads **−0.1120 / −0.1127 / −0.1176** against shifted controls of
−0.0013 / −0.0037 / −0.0251 — the magnification chain **breaks** at the deepest step, on panels
whose own windows leave z4 with raw sd 0.1597 (so it is not the dead arm the shared window gave).
The figure says why: `z4` carries a fine diagonal **cross-hatch** that `z2` magnified does not
predict, and `z2` carries a fainter one that `z1` does not.

Split into OKLab, it is in the **lightness** channel — which under Replace-L *is* `spread_shape`,
the ensemble quantity. Chroma is `hue_ab(shape_vec)`, the nominal copy:

| panel | L hf rms | a hf rms | b hf rms | L hf/range |
|---|---|---|---|---|
| z1 | 0.02222 | 0.01303 | 0.00764 | 0.0364 |
| z2 | 0.02429 | 0.00283 | 0.00110 | 0.0424 |
| z4 | **0.11967** | 0.00050 | 0.00038 | **0.2908** |

240× more high-frequency content in the ensemble channel than the nominal one, and it **grows
with depth**.

## The hypothesis

`spread_shape` is a spread over `E+1 = 8` copies at a **fixed** Halton (2,3) stencil scaled to the
cell. Every cell uses the same 8 offsets. As the cell centre steps pixel by pixel across a field
with structure near the cell scale, the stencil slides through that structure and the estimate
beats against it. That is moiré in the estimator, not in the field — and it should **strengthen**
with zoom, because a smaller cell averages over less. It does: 0.0364 → 0.0424 → 0.2908.

## Predictions

1. **`pcg_E8`** — per-pixel random offsets, no fixed stencil. The cross-hatch's **coherence**
   collapses. Amplitude may not fall; Pcg is the noisier scheme on this project's own record
   (per-quad floor 0.4796 against 0.0010), so read coherence, not amplitude.
2. **`halton_E4`, `halton_E16`** — different stencils. Cross-hatch present with a **different
   period**. Same-period-different-count would refute the mechanism.
3. **`jitter_q` (0.25), `jitter_f` (1.0)** — same stencil, different extent. Period **moves with
   the extent**.
4. **Nominal fields** (`shape_vec`, `t_end`, `d_min`, `drift`, `steps`) — **no cross-hatch**.
   Nothing is sampled, so nothing can beat.
5. **`nominal_steps`** is the substep theory at field level. If its ripple tracks `spread`'s,
   step-count level sets are drawing it after all — the `eta` arm above says they are not.

## If it is a sampler beat

There is a real tension to state rather than fix quietly. The fixed stencil is a **non-negotiable**
— it is what makes copy `k` sit at the same offset in every footprint, giving common random numbers
and the parent/child correlation of 0.9998 the refinement criterion rests on. Randomising it to
kill the moiré would cost that. The resolution, if needed, is to decouple: the fixed stencil for
`spread` (scheduling) and a different estimator for display lightness. Not a change to make on one
render.

---

# Outcome, part 1 — prediction 4 REFUTED

The nominal fields carry the beat. All of them, at one orientation. `z1` window, band 2.2–40 px,
peak over the median bin:

| field | ensemble? | λ px | angle | prom | distinct |
|---|---|---|---|---|---|
| `halton_E8_spread` | **yes** | 12.64 | 69.8° | 10837 | 65536 |
| `nominal_shape0` | no | 10.59 | 65.6° | 2382 | 65535 |
| `nominal_shape1` | no | 12.64 | 69.8° | 4451 | 65536 |
| `nominal_shape2` | no | 14.20 | 70.6° | 38916 | 65536 |
| `nominal_tend` | no | 11.39 | 69.1° | 14181 | 45571 |
| `nominal_dmin` | no | 29.96 | 69.4° | 10007 | 65536 |
| `nominal_drift` | no | 7.49 | 69.4° | 9569 | 65536 |
| `nominal_steps` | no | 11.39 | 69.1° | **152858** | 19728 |

Nothing is sampled in seven of those eight rows, so nothing can beat against a stencil. The
striations are one family in the IC plane, read by seven different observables at 65–71°.

**And the channel reading that motivated the hypothesis was a scale error.** "240× more
high-frequency content in L than in chroma" compared absolute figures across channels whose
ranges differ by 180× (L 0.41157, a 0.00228 at `z4`). Against each channel's own range they are
**0.2908 / 0.2180 / 0.2188** — comparable. Both arms carry it; the ensemble arm is not special.

**`nominal_steps` is the sharpest beat in the table** — the integer substep count, prominence
152858, an order above anything else, at exactly `nominal_tend`'s λ and angle. The substep count
*is* strongly banded. The `eta` arm settles the direction: at 3.07× the steps the shape bands do
not move, so the trajectories are banded and the step count inherits it, not the reverse.

Arms 1 and 3 (offsets, sampling extent) are still running and are now expected to be **null** —
recorded here so a null is not read afterwards as a confirmation.

---

# What the striations are NOT — three mechanisms closed, one open

Measured on the `z1` fields, sampling **along the band normal** (−21°). A profile along the bands
reads flat and says nothing; two earlier row profiles looked clean for exactly that reason.

**Not sync-cadence quantisation.** `t_end` lands exactly on a sync boundary on 0.3048 of the
frame — which reads alarming until the horizon population is separated out: **0.3036** of the
frame is at `t_end = t_max = 50`, itself a multiple of `dt = 0.4`. The genuinely quantised
fraction is **0.0012**. This is the standing *"on a boundary conflates quantised with finished"*
result, and separating the two is what makes the level readable without a second cadence.

**Not a `t_end` staircase.** Across the bands `t_end` runs **4.2847 to 4.5625** with 246 distinct
values in 440 samples — continuous. 45571 distinct over the whole frame.

**Not one-more-revolution of the tight pair.** The standard mechanism for striations in a
three-body IC plane predicts `d(t_end)` across one band ≈ the inner-binary period. Measured:

| quantity | value |
|---|---|
| band wavelength | 12.64 px = 3.778e-3 chart units |
| `d(t_end)/dx` along the normal | −2.4554e-3 per px |
| `d(t_end)` across one band | 3.1037e-2 |
| `d_min` p50 on the cut | 9.8166e-4 |
| `2π√(d_min³/M)`, M = 1 | 1.9325e-4 |
| **ratio** | **160.6** |

Two orders out. Whatever sets the spacing, it is not the orbital period at closest approach.

**Open, and labelled as open.** The identity of the striations is not established. What is
established is that they are real structure in the IC plane at a single orientation (65–71° across
seven observables), invariant to a 3.07× change in step size, and chart-locked at the z1→z2
magnification step.

---

# Outcome, part 2 — arms 1 and 3 REFUTED. The sampler is not it.

`z1` window. Offsets, count and extent all changed; **λ and angle do not move at all.**

| arm | offsets | E+1 | extent | λ px | angle | prom | vs E8 | shifted |
|---|---|---|---|---|---|---|---|---|
| `halton_E8` | Halton | 8 | 0.5 | **12.64** | **69.8°** | 10837 | 1.0000 | 0.0010 |
| `pcg_E8` | **PCG, per-pixel** | 8 | 0.5 | **12.64** | **69.8°** | 1958 | 0.6759 | 0.0013 |
| `halton_E4` | Halton | **4** | 0.5 | **12.64** | **69.8°** | 2013 | 0.7905 | 0.0006 |
| `halton_E16` | Halton | **16** | 0.5 | **12.64** | **69.8°** | 21922 | 0.9328 | 0.0012 |
| `jitter_q` | Halton | 8 | **0.25** | **12.64** | **69.8°** | 1512 | 0.6440 | 0.0007 |
| `jitter_f` | Halton | 8 | **1.0** | **12.64** | **69.8°** | 30817 | 0.5523 | 0.0014 |

**The positive evidence is stronger than the null.** Prominence rises monotonically with sample
count (2013 → 10837 → 21922 at E+1 = 4, 8, 16) and with extent (1512 → 10837 → 30817 at 0.25,
0.5, 1.0). A sampling artefact weakens with more samples. This one **sharpens** — which is what
estimating a real thing better looks like.

**Guard correction: the bar was wrong, not the code.** `jitter0` asserted bitwise zero and read
max `1.272e-16`. A variance taken as a difference of moments does not cancel exactly on identical
inputs. Rebased on the field it guards: `2.14e10×` below it. PASS.

# A CASCADE AT ONE ORIENTATION, and it survives the filter

Sensitivity of the peak to the high-pass width and the search band:

| σ | 2.2–40 | 4–60 | 6–128 | 2.2–128 |
|---|---|---|---|---|
| 1.5 | 2.24 @ 64.6° | 12.64 @ 69.8° | 12.64 @ 69.8° | 2.24 @ 64.6° |
| 3.0 | **12.64 @ 69.8°** | **12.64** | **12.64** | **12.64** |
| 6.0 | **12.64 @ 69.8°** | **12.64** | **12.64** | **12.64** |
| 12.0 | 38.16 | 57.24 | 80.95 | 80.95 @ 71.6° |
| raw | — | — | — | 80.95 @ 71.6° |

Three tiers at **one orientation**: λ ≈ 81 px (the ribbon, present in the unfiltered field), 12.64
(the striations), 2.24 (near Nyquist). Ratios **6.40** and **5.64**.

This closes two loose ends rather than leaving them. `z1@2× vs z2` registered 0.78–0.88 and not
~1 because each magnification exposes a finer tier the coarser panel cannot carry; `z2@2× vs z4`
broke to −0.11 because by `z4` the dominant tier has reached the sampling limit (measured λ 2.44
px on that panel).

**Prediction for the `z4` arms, recorded before they return.** `z4` is 4× `z1`, so the tier at
λ 2.24 px in `z1` should appear near **λ 9 px** there, with a new tier near 2 px beneath it, all
at 63–72°, and all six sampler arms should again agree on λ and angle.
