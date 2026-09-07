# What 68,685 leaves say about the refinement criterion

**Source:** all 60 committed `.prnq` quadtree dumps in `prin-rs/results/` (17 `charts/`,
22 `criterion/`, 21 `vertical/`), parsed with `principia_prnq_parser.py`. Figures:
`principia_refinement_fig1.png`, `principia_refinement_fig2.png`.

**Caveat, stated first.** These trees were rendered at `half = 1.0` where the reference uses
**3.0** — a 3× crop. Findings below are marked **[framing-independent]** where they are statements
about a distribution rather than about which patch is shown, and **[re-take]** where they are not.

---

## 1. The criterion is an always-split rule wearing a threshold's clothing

`tau_display = 1e-4` in every committed run. The **median leaf spread is 6.44e-4** — the threshold
sits at the **0.4th percentile** of its own distribution. **[framing-independent]**

| `tau` | quads exceeding |
|---|---|
| **1e-4 (current)** | **99.6%** |
| 1e-3 | 22.4% |
| 3e-3 | 2.1% |
| 1e-2 | 1.0% |

**And the second condition does not save it.** The split predicate is
`spread > tau AND alpha > alpha_hi`, with `alpha_hi = 0.2`:

```
spread test passes : 94.4%
alpha  test passes : 85.8%
BOTH               : 80.3%      <- the effective split rate
```

`alpha`'s median is **0.996** and 76% of leaves exceed 0.8 — the field is in **linear response**
almost everywhere, so "refining halves the spread" is true nearly universally. Both gates are open.

**Consequence:** 13 of 17 charts have **99.4–100%** of leaves at max depth; three
(`latent_shape`, `latent_mass`, `mass_simplex`) are *exactly* 4096 = 4⁶, completely uniform.
Depth variance 0.002–0.031. **Balanced mode does not currently exist.**

---

## 2. The structure fields are saturated and carry no information **[framing-independent]**

`n_hot_within = 64` — every footprint hot — in **87% of all leaves** and 99–100% on 15 of 17
charts. `n_components_within = 1` universally. `perimeter_ratio_within` is **0 at the median**.

The hot mask is one blob covering the whole quad, everywhere. **The spatial measures cannot
discriminate as currently thresholded** — the same instrument-cannot-see-its-own-signal failure as
the linear-ramp spread image. This is why the relative threshold is mandatory rather than tidy.

---

## 3. The event arm is dead

`ensemble_spread = max(spread_shape, spread_event)`. Pooled over all 68,685 leaves:

```
shape arm wins : 99.8%
event arm wins :  0.2%
event arm EXACTLY ZERO : 99.8%
```

**One of the two contributors is identically zero almost everywhere.** Terminal states are
absorbing, so copies in a terminated quad share an outcome and the disagreement fraction collapses.
The two-field set is a one-field set in practice.

That is not an argument to remove it — it fires precisely where outcomes disagree, which is the
boundary case worth catching — but **any analysis treating `ensemble_spread` as a blend of two
signals is describing something that does not happen.**

---

## 4. Refinement follows escape, not structure **[re-take at half = 3.0]**

`preset_shape`, shallow leaves versus max-depth leaves:

| | shallow (129) | at max depth (448) |
|---|---|---|
| `escape_fraction` | **0.094** | **1.000** |
| `spread_median` | 1.0e-4 | 1.3e-2 |

Non-escaping regions have collapsed spread and are abandoned; fully-escaped regions keep high
spread and refine to the floor. `preset_shape_pl` shows the same in its wireframe: 2974 leaves, of
which **14 are shallow** (12 `Floor`, 2 `Keep`) — and those fourteen sit over the structured centre.

**Both complaints are true at once: it refines nearly everything, and the handful it declines are
the interesting parts.**

---

## 5. THE BIG ONE: selectivity is a property of the CHART, not the criterion

Per dump, spread dynamic range (`p99/p1`) against leaf depth variance, over 44 dumps:

**Spearman = +0.546** (and **−0.403** against % at max depth).

| chart | spread p99/p1 | % at max | depth var |
|---|---|---|---|
| `latent_mass` | **2×** | 100.0% | 0.000 |
| `latent_shape` | 3× | 100.0% | 0.000 |
| `preset_shape_pl` | 10× | 99.5% | 0.031 |
| `body_plane` | 73× | 61.2% | 1.015 |
| `shape_sphere` | 1962× | 79.2% | 0.448 |
| `preset_shape` | **10078×** | 77.6% | 0.691 |

**Charts whose spread is flat get uniform trees; charts with wide dynamic range get selective
ones.** `latent_mass` varies by a factor of **2 across the entire chart** — there is nothing there
for any criterion to find.

**This reframes the whole diagnosis.** The criterion is not choosing badly so much as being handed
inputs with no contrast. And it makes a prediction worth testing the moment the window fix lands:
**a 3× wider window should widen the dynamic range**, because more of the sphere is in frame — so
the framing fix may substantially improve the refinement *by itself*, without touching the
criterion. Measure the dynamic range before and after and report it as its own result.

---

## 6. The pan experiment could not have failed

Nine pan steps, camera moving `cx = 0.95 → 1.05`. **All nine produce a byte-identical tree** —
844 leaves, identical hash. `max_rel_depth: None`, so the camera-relative depth predicate was off.

**The camera is not wired into the scheduler in these runs**, so the pan sequence measured nothing.
Another member of the can't-fail family, and it means §4.5's pan/zoom asymmetry is untested rather
than tested.

**Zoom, by contrast, does work** — and interestingly:

| step | half_world | leaves | % at max | depth var |
|---|---|---|---|---|
| 0 | 0.05 | 412 | 61.2% | 1.015 |
| 2 | 0.0125 | 1045 | 83.1% | 0.395 |
| **4** | **0.003125** | **214** | **7.5%** | **0.727** |
| 6–8 | ≤7.8e-4 | 16 | 100.0% | 0.000 |

Step 4 is the **most selective tree in the entire corpus** — only 7.5% at max depth. And steps 6–8
collapse to 16 leaves: the screen floor has taken over completely. **There is a zoom band where the
criterion works well**, which is worth understanding rather than averaging away.

---

## 7. Cost is bimodal — cost-aware priority is worth building

`total_substeps` per quad: **p1 = 2.04e4, p50 = 4.09e4, p99 = 2.05e6 — a 100× range**, with a
clear two-population split (72.8% below 1e5, 27.2% above, ~50× apart). **[framing-independent]**

The criterion plan said *"report the cost distribution first — if it is narrow, `spread /
compute_cost` is moot and the idea stops there."* **It is not narrow.** A quad that costs 50× its
neighbour and yields the same information is exactly what a budget-limited scheduler should
deprioritise. Proceed.

---

## 8. Integration health, and two things that pass

`error_ratio_max`: **p50 = 1.000** (exact), p90 = 1.87, p99 = 152. `worst_energy_drift`:
p50 = 1.0e-11, p99 = 1.7e-2. Both are excellent in the bulk with a heavy tail in ~1–10% of quads —
the expected signature, and the reason `error_ratio` is a boolean flag rather than a magnitude.

`n_distinct_ic = 64` at **every percentile** — no collapsed decodes anywhere in the corpus. That
guard passes cleanly.

`t_end_gradient` is **0 up to the 90th percentile** — dead at `t_max = 13`, exactly as predicted for
an escape-derived signal at the project horizon.

---

## 9. What to do

1. **Fix `tau` first, and make it relative.** Everything else is measured through a predicate that
   is currently always true. A quantile of the observed distribution, not a constant.
2. **Measure the dynamic range before and after the window fix** (§5). If a 3× window widens it
   substantially, part of the refinement problem is the framing bug and should not be attributed to
   the criterion.
3. **Do not measure the spatial fields until the mask desaturates** (§2). Measuring now returns a
   guaranteed null.
4. **Wire the camera into the scheduler before re-running any pan experiment** (§6).
5. **Build cost-aware priority** (§7) — the distribution answers the question the plan asked.
6. **Investigate the zoom band around step 4** (§6) — the most selective tree in the corpus came
   from somewhere, and it is worth knowing where.
7. **Treat `ensemble_spread` as one signal, not two** (§3), in any analysis of its behaviour.
