# Results

What the uniform kernel measured, written to be read cold. Raw output for every number is in
[`results/output/`](results/output/); the mechanisms behind them are in
[`NOTES.md`](NOTES.md), and the working agreement they feed back into is
[`CLAUDE.md`](CLAUDE.md).

Everything below is Burrau (masses 3-4-5, released from rest), Aarseth–Zare with two-pair
Levi-Civita regularisation, `E+1 = 8` copies per pixel, `t_max = 13`, f64 unless stated.

Copy offsets are the spec's **fixed Halton (2,3) prefix indexed by copy index**. The port
originally inherited the reference's per-pixel PCG stream; §1 measures both, since every result
recorded before that correction was produced on the PCG path and `--jitter-pcg` reproduces it.

---

## 1. The refinement criterion

**The criterion works, and what limits it is not what an earlier version of this document
said.** It resolves individual quads exactly where the physics is smooth and not at all where it
is chaotic — and the second half is the correct answer rather than a limitation. The limit is
chaotic divergence, which no ensemble scheme and no amount of compute reaches.

The criterion compares a parent quad against its children. A fine uniform grid already contains
every coarser scale by aggregation, so the whole exponent machinery is testable with no
quadtree: pool a 2×2 block's copies to synthesise the parent and compare against the children.
Nothing here can be an artefact of a scheduler, because there is no scheduler.

`alpha = log2(spread_parent / spread_child)`, child taken as the median over the four.

### What it discriminates

In the tame regions the shape spread scales with cell width, as a smooth field must, so halving
the cell halves the spread and refinement pays. In the chaotic regions the parent spread is
barely above the child's, so refining buys almost nothing — those pixels are **undetermined, not
under-resolved**. Making that distinction is what the criterion is for, and it makes it cleanly.
The numbers are in the table two sections down, measured by rendering at two resolutions; the
figures this document previously carried were measured by pooling and are corrected there.

### The pooled parent is not a parent

Every `alpha` in earlier versions of this document was computed by **pooling** a 2×2 block's
copies to synthesise a parent. That is not a parent, and the error is systematic.

With the spec's fixed offsets, a pooled block is four **exact repeats** of one offset pattern at
four cell centres. A *true* parent at 2× cell width carries offsets scaled to **its own** width —
a wider footprint. They are not the same ensemble. `alpha` for `sigma_E(0)`, whose true value is
exactly 1.0, measured both ways on near-field (fine 64×64, coarse 32×32):

| E+1 | scheme | pooled median | true median | pooled \|err\| | **true \|err\|** |
|---|---|---|---|---|---|
| 4 | Halton | 1.6678 | 1.0231 | 0.6678 | **0.0231** |
| 8 | Halton | 1.3857 | 1.0227 | 0.3857 | **0.0227** |
| 16 | Halton | 1.1668 | 1.0227 | 0.1668 | **0.0227** |
| 32 | Halton | 1.0661 | 1.0227 | 0.0661 | **0.0227** |
| 4 | Pcg | 1.2798 | 0.8843 | 0.2798 | 0.1157 |
| 8 | Pcg | 1.0762 | 0.9447 | 0.0762 | 0.0553 |
| 16 | Pcg | 1.0137 | 0.9762 | 0.0137 | 0.0238 |
| 32 | Pcg | 0.9887 | 0.9872 | 0.0113 | 0.0128 |

The pooled error runs 0.67 → 0.07 as `E` grows; the true error is **flat at 0.0227**. The bias is
the surrogate, not the estimator, and it is removed by **rendering at two resolutions** rather
than by calibrating a correction factor. The residual 0.0227 is what the two-resolution method
itself costs and it does not shrink with `E`.

Under PCG the two errors partly cancel: pooling gives 4× the samples and per-footprint
randomisation blurs the surrogate mismatch, so the pooled figure looks better than the true one
at small `E`. Two wrongs, partially offsetting — which is why this went unnoticed until the
fixed scheme removed the randomisation.

**Any `alpha` measured by pooling carries this**, here and in the prior NumPy work.

### The criterion, re-measured properly

`alpha` for `spread_shape`, true two-resolution rendering, fixed Halton prefix, `t = 13`:

| region | p10 | median | p90 | **interdecile** |
|---|---|---|---|---|
| near-field | −0.6568 | 0.0368 | 0.6696 | **1.3264** |
| body2 core | −0.3673 | 0.1844 | 0.7405 | **1.1078** |
| mid-field | 1.0224 | 1.0229 | 1.0235 | **0.0010** |
| far | 1.0228 | 1.0230 | 1.0232 | **0.0004** |

Separation between region medians: **0.9862**.

This is a different picture from the pooled one, and a better one:

- **In the tame regions the exponent is essentially exact.** `alpha = 1.023` with an interdecile
  of `0.0004`–`0.001`. Per-quad decisions there are not merely possible, they are trivial. The
  pooled measurement reported ~0.63 scatter in these regions and **all of it was surrogate
  error**.
- **In the chaotic regions the scatter is 1.1–1.3** — twice what pooling suggested, and it is not
  sampling noise. It is chaotic divergence.
- The tame median `1.0229` matches the `0.0227` residual from the control above. Same constant:
  the two-resolution method carries a small systematic `+0.023`, and it is visible in two
  independent quantities.

**So "regions not quads" is too blunt.** The criterion resolves individual quads perfectly well
where the physics is smooth, and not at all where it is chaotic. That is the correct behaviour
rather than a limitation: **"not resolvable per quad" is the answer for a chaotic quad**, and the
scatter is the measurement, not an error bar around one.

**And no amount of compute changes the chaotic case.** Switching from the PCG stream to the
spec's fixed Halton prefix cut the control's sampling noise by a factor of 267,000
(`var(alpha_E)`: `3.75e-2 → 1.40e-7`) and moved `alpha_shape`'s interdecile not at all
(0.6313 → 0.6326 pooled; 1.3051 → 1.3264 true). Sampling error was ~7% of the variance; the rest
is the physics the instrument exists to measure. There is no ensemble size that buys per-quad
resolution in a chaotic region.

### Read the interdecile, not the variance

The two summaries disagree, and the disagreement is the point. For `alpha_shape` under Halton:

| | value |
|---|---|
| variance | 5.331e-1 |
| sd | 0.7302 |
| interdecile (p90−p10) | 0.6326 |
| interdecile / sd | **0.866** (normal: 2.563) |
| excess kurtosis | **110.0** (normal: 0) |

An interdecile three times *narrower* than a normal distribution's, with excess kurtosis of 110:
**the variance lives in the tails and the interdecile describes the bulk.** That resolves the
tension between "variance fell 6.9%" and "scatter unchanged" — both are true, of different parts
of the same distribution.

**A scheduler decides per typical quad, so the interdecile is the measure that matters, and it is
the one that did not move.** Do not quote the 6.9% variance reduction as the improvement; it is a
statement about the tail.

### The bias only a control could have caught

None of the above is visible without a quantity whose answer is known exactly. A smooth offset on
a statistic with no oracle does not look like an error, it looks like a result — `+38.6%` would
have read as "parent spread grows faster than linearly in cell width", a plausible physical
statement and completely wrong. It was visible only because `sigma_E(0)` has an *exactly* known
value, and it was *diagnosable* only because the fixed offset scheme removed the noise that had
been hiding it.

### Two footnotes

**The `alpha_E` control variate was dropped.** `rho` is −0.079 (Halton) and −0.042 (PCG), and the
regression form makes the floor slightly worse either way. Under Halton the fitted
`beta = −153.76` is a division by nothing — `var(alpha_E) = 1.4e-7`. Two variance reductions
targeted the same component and the cheaper one won: the offset scheme had already removed what
the control variate would have corrected.

**The fixed prefix's advantage grows with `E` rather than shrinking** — L2 star discrepancy ratios
against PCG are 0.748, 0.624, 0.489, 0.395 at `E+1 = 4, 8, 16, 32`. Low-discrepancy sequences are
usually most valuable at small `N`, so this suggests the raw unscrambled Halton prefix's early
terms are less well distributed than its later ones, which is a known property in low dimensions.
Reported, not chased; scrambling is the standard remedy if it ever matters, and nothing here
depends on it.

`alpha` was **not** smoothed over neighbouring quads. It is the obvious variance reduction and it
is wrong here: `alpha` varies smoothly except at boundaries, and boundaries are exactly what a
refinement decision is about.

---

## 2. Conclusions drawn at n = 64 do not survive

Subsampling one fixed 128×128 grid: interdecile spread of a statistic over 200 random draws of
`n` pixels, as a fraction of its full-grid value. The same physical quantity throughout, so this
isolates sampling error. Below ~0.1 a conclusion drawn from `n` pixels is stable.

| quantity | truth | n=16 | n=64 | n=256 | n=1024 | n=4096 |
|---|---|---|---|---|---|---|
| drift median | 2.1312e-9 | 1.544 | 1.070 | 0.851 | 0.445 | 0.227 |
| drift p99 | 5.0117e-3 | 7.716 | 2.021 | 3.051 | 2.056 | 1.146 |
| drift max | 1.4909e4 | **0.000** | **0.000** | **0.000** | 0.003 | 1.000 |
| error_ratio median | 1.0000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| error_ratio p99 | 1.3655e2 | 0.807 | 1.870 | 2.882 | 2.235 | 1.061 |
| spread_shape median | 1.3244e-3 | 0.514 | 0.243 | 0.154 | 0.065 | 0.033 |
| spread_shape p99 | 2.9705e-2 | 4.084 | 4.001 | 6.768 | 5.750 | 3.345 |
| d_min_true median | 8.4693e-3 | 0.370 | 0.234 | 0.140 | 0.060 | 0.034 |
| frac collision | 2.2095e-2 | 2.829 | **2.122** | 1.061 | 0.575 | 0.243 |
| frac er > 10 | 1.6235e-2 | 3.850 | **2.887** | 1.203 | 0.602 | 0.331 |
| frac drift > 1e-3 | 1.4099e-2 | 4.433 | **2.216** | 1.385 | 0.693 | 0.329 |
| frac spread_event > 0 | 2.3315e-2 | 2.681 | **2.010** | 1.005 | 0.545 | 0.251 |

**At `n = 64` every fraction has an interdecile scatter of 2.0 to 4.4 times the quantity
itself.** Two independent studies at that `n` would routinely disagree by more than a factor of
two — and a factor-of-two disagreement on a fraction is a qualitative disagreement, not a
quantitative one. That is the mechanism by which a 1.2× effect turns out to be 18.8×, quantified
rather than anecdotal.

**What this implies for any conclusion drawn at `n = 64`:**

- **A fraction, rate, or "how often" claim at `n = 64` carries no information about its own
  magnitude.** It can be off by more than 2×, in either direction, from sampling alone. Any such
  claim in prior work should be treated as an order-of-magnitude estimate at best, and any
  *comparison* between two such fractions should be treated as unresolved unless the effect
  exceeds about 4×.
- **A median of a well-behaved quantity is usable from a few hundred pixels.** `spread_shape` and
  `d_min_true` medians reach 0.15 at `n = 256`.
- **`drift median` is slower** — 0.85 at `n = 256` — because its distribution spans twelve
  orders. Medians of heavy-tailed quantities are not cheap.
- **A p99 of anything heavy-tailed is not estimable at any `n` tested.** `spread_shape p99` sits
  at 3.3–6.8 throughout and does not improve with `n`. Tail quantiles need either far more
  samples or a different estimator.

### The row that reads backwards

`drift max` shows 0.000 scatter at `n ≤ 256` — apparently the most stable quantity in the table —
and 1.000 at `n = 4096`. It is not stable. The tail is a **single pixel of 16384**, and small
samples essentially never draw it, so the statistic is perfectly reproducible at the wrong
answer. **A max statistic's apparent stability at small `n` is the statistic never seeing the
tail.**

This matters beyond the table: every "worst case observed" number in prior work at `n = 64` is
subject to it. Not noisy — *confidently wrong*, and stable across repetitions in a way that would
survive re-testing.

---

## 3. Seven pixels, a cliff, and the remedy

At 128×128, near-field carries **7 pixels of 16384 (0.043%) with `|dE/E| > 1`**, worst
`1.4909e4`. They have lost more than the total energy of the system; their trajectories mean
nothing. At 64×64 they are not grid points and the worst drift in the region is `5.1e-1`.

| pixel | (jx,jy) | drift_max | error_ratio | d_min_true | gamma_max |
|---|---|---|---|---|---|
| 16110 | (110,125) | 1.4909e4 | 3.9615e8 | 1.8381e-3 | 1.0911 |
| 15989 | (117,124) | 4.0856e1 | 7.3239e5 | 1.7260e-3 | 0.7969 |
| 16351 | (95,127) | 6.3484 | 1.9212e5 | 2.0439e-3 | 0.9829 |
| 16229 | (101,126) | 6.1218 | 1.7145e5 | 2.2125e-3 | 0.9041 |

All finite. All clustered in one corner. All at `d_min_true ≈ 2e-3` against
`r_coll = 1e-3 R = 2.214e-3` — a near-collision the run is not permitted to terminate on. All
with `gamma_max ≈ 1`, so the regularised Hamiltonian residual is order unity and the trajectory
is not being integrated so much as invented.

### `error_ratio` flagged 7 of 7, and MAD would have flagged none

This is the strongest evidence for the switch from MAD to the maximum deviation. The damaged /
healthy separation is **1.06 with MAD and 59.51 with max deviation**; a pixel whose worst copy
drifted 120× the total energy reported a MAD-based `error_ratio` of 1.1369, inside the healthy
p99 of 1.0756. In one 24×24 run, a pixel with 2.9% energy drift reported `error_ratio = 204.8`
under max deviation and **0.9018** under MAD — quieter than it was at `t = 0`.

MAD was specified for a sound reason: a standard deviation returns NaN the moment one copy is
non-finite, on exactly the pathological pixel the field exists to flag. But robustness to
outliers is the opposite of what a detector needs. With 8 copies, one wild value sits above the
median of eight deviations and is arithmetically invisible. The maximum deviation satisfies the
original requirement *better*: a non-finite copy gives an infinite deviation, which is the
correct answer where a std gives NaN.

### It is a cliff, not a slope — which rules out adaptive eta by region

Checked, because insensitivity to step size would mean a wrong equation rather than a resolution
problem:

| pixel | eta=1e-2 | eta=3e-3 | eta=1e-3 | eta=3e-4 |
|---|---|---|---|---|
| 16110 | 1.4909e4 | 5.9589e-9 | 3.1161e-11 | 1.4885e-11 |
| 15989 | 4.0856e1 | 5.3896e-9 | 3.6647e-11 | 1.3174e-11 |
| 16351 | 6.3484 | 7.9098e-9 | 3.9721e-11 | 1.7923e-11 |
| 16229 | 6.1218 | 9.8056e-9 | 3.5737e-11 | 1.5836e-11 |

**Thirteen orders of magnitude for a 3.3× change in `eta`.** The drift does fall with step size,
so the equations are right — but it falls off a cliff, not down a slope.

That geometry decides the remedy. **An adaptive-`eta`-by-region policy cannot work**, because
there is nothing to interpolate along: the error is `1e-9` on one side of the cliff and `1e4` on
the other, with no intermediate regime to tune against, and the cliff's location depends on where
a step happens to land relative to a close approach rather than on any property of the region.
A region-wide `eta` would either be far too fine everywhere or fail on the same pixels.

**Flag-then-re-integrate does work**, and it is what the kernel now does: run the grid, then
re-integrate pixels `error_ratio` flags at `eta/4`, up to three passes, recording the coarse
value, the refined value and the `eta` used for every pixel. Bounded, one extra evaluation of a
shrinking subset, no tree and no state carried between pixels.

Measured on 128×128 near-field:

| | no refinement | with refinement |
|---|---|---|
| drift max | 1.4909e4 | 5.0876e-4 |
| drift p99 | 5.0117e-3 | 3.2498e-5 |
| pixels `|dE/E| > 1` | 7 | **0** |
| pixels `|dE/E| > 1e-3` | 231 | **0** |
| pixels re-integrated | 0 | 266 (1.62%) |
| wall clock | 11.63 s | 12.00 s (+3%) |

The unrefined tail grows with resolution — worst drift `5.1e-1` at 64², `1.49e4` at 128²,
`1.38e7` at 256² — so this is not a fixed cost that can be ignored at scale.

### Two honest limits on the remedy

**One pass is not always enough.** `deep interior` needed three: one pass takes it from `1.10e12`
to `1.99e1`, still unusable, because 14% of that region is flagged and its close approaches are
far deeper. With three passes it reaches `1.15e-1` with **0 pixels still flagged**. The bound is
a parameter and a pixel still flagged after the last pass keeps its `error_ratio`, so an
unrepaired pixel is reported rather than silently accepted.

Across all six regions at 256×256, drift max before and after:

| region | pixels refined | before | after |
|---|---|---|---|
| far | 0 | 4.702e-11 | 4.702e-11 |
| mid-field | 0 | 4.129e-10 | 4.129e-10 |
| body1 slice | 0 | 9.542e-8 | 9.542e-8 |
| near-field | 890 | 1.377e7 | 3.120e-4 |
| body2 core | 736 | 3.110e10 | 6.139e-4 |
| deep interior | 9228 | 1.102e12 | 2.543e-1 |

Three of six regions need no refinement at all — `far`, `mid-field` and `body1 slice` flag zero
pixels. The cost is concentrated exactly where the physics is.

**And the pass budget is calibrated on f64.** The same near-field slice at f32 needs more passes
for the same grid: at 128×128, f64 clears with three and f32 needs about six (one pass leaves 421
of 16384 still flagged, six leaves none, drift max `4.25e1 → 5.37e-4`). At 256×256 with the
default three passes, f32 leaves **1578 of 65536 still flagged**. Finer `eta` means more steps,
and at f32 roundoff accumulation eats into what truncation error gives back, so convergence is
slower rather than absent. Raise `--refine-passes` for f32 runs, and read `n_still_flagged`.

**`error_ratio` detects spread, not drift.** After three passes `deep interior` at 128² has zero
pixels above the flag threshold and a worst drift of `1.146e-1` — 11% energy error (`2.543e-1` at
256²). A pixel whose
eight copies all drift *together* has a low ratio however badly they drift. This is the limitation
BRIEF §4 already names in a different form, and it means the refinement repairs what
`error_ratio` can see. `energy_drift_max` is dumped per pixel and is the quantity to threshold on
if absolute conservation is what matters.

---

## 4. The latching guard ships inert, and that is the result

`spread_event` measures disagreement over the **event class** — which pair is currently the
tightest binary, at each sync boundary, joined with the terminal outcome for copies that have
terminated.

Within a single `t = 13` run, **130 of the 165 pixels that ever disagree have re-agreed by the
horizon.** Prior numpy work saw the same shape in a continuous measure — ensemble spread falling
6× between `t = 6` and `t = 8` — which is why its divergence accumulator latches. The obvious
move is to latch this one too.

**A discrete label has a failure mode a continuous divergence measure does not.** If two pairs are
near-equal in separation, copies can disagree about which is *tightest* without their trajectories
having diverged at all. Measured — the second-tightest over tightest separation at the boundary
where the copies first disagree:

| population | min | p10 | median | p90 | max | n |
|---|---|---|---|---|---|---|
| all that ever disagree | 1.0000 | 1.0006 | 1.0040 | 1.0884 | 2.3587 | 165 |
| **re-agreed by the horizon** | 1.0000 | 1.0006 | **1.0030** | 1.0193 | 1.1098 | 130 |
| still disagreeing | 1.0004 | 1.0066 | 1.0797 | 1.1636 | 2.3587 | 35 |

**129 of the 130 that re-agreed were at a near-tie** — essentially exact. An unguarded running max
would light 165 pixels where 35 have genuinely diverged: **79% of the firing pixels lit
permanently for a labelling artefact**, and a latched artefact never clears.

The tie ratio cannot be the guard — the populations are shifted, not separated. **Persistence
can**: artefacts last one boundary (median run 1, max 2), genuine divergence persists (median run
10). A run of 3 admits 0 of 130 artefacts. `spread_event_latched` joins that guard with the
playhead value and is lit on 35/35 genuine and 0/130 artefact.

**And then it changes nothing.** Evaluated at every boundary of one run, the guarded latch tracks
the playhead value exactly:

| k | t | playhead | latched | unguarded max |
|---|---|---|---|---|
| 3 | 1.6250 | 0 | 0 | 0 |
| 7 | 3.2500 | 0 | 0 | 26 |
| 15 | 6.5000 | 0 | 0 | 26 |
| 23 | 9.7500 | 22 | 22 | 165 |
| 31 | 13.0000 | 35 | 35 | 165 |

Every genuine disagreement on this slice persists to the horizon, so latching adds nothing here.
It ships as cheap insurance for regions where one does re-agree, and it costs nothing to carry.
The unguarded version over-reports by 4.7×.

**The plain playhead `spread_event` stays the spec field and is what `ensemble_spread` uses.** All
three are dumped.

---

## 5. The scheduler — the criterion in a loop

Everything above measured the criterion on **one split in isolation**. `prinq` descends from one
quad, adaptively, at a fixed playhead. `N = 8` footprints per quad axis, `E+1 = 8` copies each, so
**one quad is 512 trajectories** and the budget is counted in quads.

### It terminates, and the floor engages

`alpha_hi = 0.2`, `tau = 1e-4`, **no depth cap**, budget 50 000 quads:

| region | quads used | leaves | depth | terminated | floored | hit the cap |
|---|---|---|---|---|---|---|
| far | **21** | 16 | 2 | 100% | 0.0% | no |
| near-field | **4617** | 3463 | 12 | 100% | 17.6% | no |
| deep interior | **29** | 22 | 4 | 100% | 40.9% | no |

The Wada-dense-boundary fear — that spread stays high however far you refine — does not
materialise here. No leaf hit the budget, the depth cap, or the precision floor. The floor branch
works and bites hardest where it should: 40.9% in the most chaotic region, 0% in the tame one,
which exits through *keep* because its spread is already below `tau`.

**But what terminates the descent is `tau`, not the floor.** In near-field's deepest two levels the
exponent has median **+3.945** while the spread is `1.4e-5`, below `tau`. Those quads split because
their *parents* were above `tau`, and the split collapsed the spread ~16×. 82.4% of leaves exit
through *keep*. So this is termination **at `tau = 1e-4`**, and at `tau = 1e-8` the same descent
exhausted a 2000-quad budget with 869 leaves still wanting to split. Bounded, not general.

### The threshold does more work than the criterion

near-field, budget 2000:

| tau | alpha_hi | quads | leaves | depth | budget-blocked |
|---|---|---|---|---|---|
| 1e-8 | **0.20** | 1997 | 1498 | 8 | 869 |
| 1e-8 | 0.50 | 25 | 19 | 3 | 0 |
| 1e-6 | **0.20** | 1997 | 1498 | 8 | 869 |
| 1e-4 | 0.20 | 1997 | 1498 | 8 | 343 |
| 1e-3 | 0.20 | 441 | 331 | 8 | 0 |
| 1e-2 | 0.20 | 21 | 16 | 2 | 0 |

`tau = 1e-8` and `1e-6` give **identical trees** — the spread never falls that low, so `tau` never
binds. Meanwhile `alpha_hi` from 0.20 to 0.50 collapses the tree **80×**. near-field's `alpha`
median is +0.389, so the threshold sits *inside* the distribution and a small move flips most
decisions. Which knob binds is region-dependent: at `tau = 1e-4` deep interior stops at 29 quads
while near-field runs to 4617.

### The aggregation is three schedulers wearing one name

A quad holds `N²` footprint spreads and needs one number. Budget 6000:

| region | agg | quads | leaves | depth | floor | leaf jaccard vs median | decisions differing |
|---|---|---|---|---|---|---|---|
| near-field | median | 4617 | 3463 | 12 | 609 | — | — |
| near-field | mean | 5997 | 4498 | 9 | 847 | 0.0963 | **54.1%** |
| near-field | p90 | 1141 | 856 | **14** | 472 | 0.0283 | **49.1%** |
| deep interior | median | 29 | 22 | 4 | 9 | — | — |
| deep interior | mean | 161 | 121 | 7 | 38 | 0.1260 | **34.5%** |
| deep interior | p90 | 81 | 61 | 7 | 30 | 0.2388 | **34.5%** |

**Half the shared decisions flip and the trees overlap by 3–13%.** median under-refines structure,
being blind to a thin filament crossing a quad; mean over-refines and is the only configuration
here to hit a cap; p90 refines deepest and narrowest and floors 55% of its leaves. Three different
intentions — resolve the typical footprint, the total, or the worst — and which is wanted is a
display question this cannot settle. What it settles is that the choice must be stated wherever a
tree is quoted.

### Coarse `N` over-refines, which is the opposite of the stated concern

Fixed budget 4000 quads. The worry on record was that a low `N` makes a quad call itself *coherent*
by undersampling its area:

| region | N | traj/quad | leaves | depth |
|---|---|---|---|---|
| near-field | 4 | 128 | **106** | 5 |
| near-field | 8 | 512 | 19 | 3 |
| near-field | 16 | 2048 | **16** | 2 |
| deep interior | 4 | 128 | **40** | 4 |
| deep interior | 16 | 2048 | **16** | 2 |

Leaf count falls monotonically with `N`. A coarse quad calls itself **uncertain**, not coherent —
a noisy spread estimate biases toward *refine*, the conservative failure direction — so `N = 4`
spends four times the quads of `N = 16` on the same region. The cheaper quad is a false economy.

The `N = 7` probe of parent–child common random numbers is **inconclusive**: it sits between 4 and
8 in near-field and matches 8 exactly in deep interior. The `N` trend dominates; no separable CRN
effect.

### The reliability signal works, and is 9× cheaper

Equal budget, near-field: the **sibling-spread** policy reaches depth 11 on **497 quads** against
the α policy's depth 12 on **4617**. It floors 63% of leaves against 18%. Where the four sibling
exponents scatter, the unreliability *is* the answer and no trustworthy α is needed.

Cheaper is not the same as better: its median quad spread is `7.97e-4` against `5.78e-5`, so it
leaves an order more uncertainty on the table. The two trees are not nested — 109 leaves are
sibling-only — so it is not a truncation of the other. `alpha_sibling_spread` is the **range of
four samples** and is itself noisy; this result makes the policy worth pursuing and that noise the
next thing to characterise.

### Thrash is real, and it is per-quad noise

Adjacent leaf pairs with spreads within 1.5× that sit at **different levels**. near-field:
**0.3392 at N=4, 0.2179 at N=8, 0.0733 at N=16.** More samples, less noise, fewer contradictory
neighbours — which is what identifies it as noise rather than structure.

Adjacent siblings share an edge column of footprints *identically* (`1/N` of a quad: 25% at N=4,
6.25% at N=16), which makes neighbours more alike and **suppresses** apparent thrash. The
most-suppressed row shows the most thrash, so the true N=4 figure is above 0.339 and the trend is
understated.

### Ordering matters; which ordering does not

At a budget that binds (1500 quads, near-field): shuffling changes **600 of 1123 leaves, 42% of the
tree**, at identical cost. But `spread` and `spread × area` produce **byte-identical** trees. The
priority function is load-bearing; the choice between these two candidates is free.

### The tree looks right in one region and wrong in another

Against the spread base, near-field refines in coherent bands — but leaves the brightest, thinnest
filaments in coarse quads. **`deep interior` fails outright**: it leaves the large high-spread wedge
and the bright diagonal bands at level 2 while spending its only fine refinement on an unremarkable
central patch. `far` is uniform at level 2, correctly.

The first overlay was drawn over the **outcome** image, which is 97.7% one colour in near-field, and
that version would have passed inspection. Both are now written; only the spread base is a check.

---

## 6. What this invalidates in the prior numpy investigation

Those documents live outside this repository. Each item below is a specific correction with the
measurement behind it.

### "Aarseth–Zare is unusable at f32" — **withdrawn. It was the branch cut, not arithmetic.**

The inverse Levi-Civita map computes `u0 = sqrt((|rho| + rho.x)/2)` first and derives `u1` from
it. That sum cancels catastrophically when `rho` points along negative x. **The Burrau default
sits exactly on the cut**: bodies 1 and 2 start at the same `y`, so their separation is `(3, 0)`
and, with reference body 2, registers at exactly 180° before anything moves.

Worst case over 3600 orientations: `6.206e-11` unstable against `4.108e-16` conditioned at f64;
**`2.2e-2` against `5.96e-8` at f32**. The defect is correctness, not precision — the cut is fixed
in the coordinate frame, so accuracy depends on a configuration's absolute orientation, while the
physics is rotationally invariant.

With the conditioned branch, f32 gives median drift `9.293e-6` against f64's `2.755e-9` — three
and a half orders, which is exactly what `eps ≈ 1.19e-7` compounded over ~5000 RK4 steps predicts
— and **outcome labels agree with f64 on 1022 of 1024 pixels**. That is a usable kernel.

With the unconditioned branch at f32, `spread_shape` inflates **32×** (median `6.1582e-2` against
the f64 truth `1.9095e-3`) while single-trajectory drift still looks reasonable. **That is the
exact shape of the original dispute**: the ensemble diagnostic breaks early while the trajectory
looks fine. The hypothesis on record was reference-body switching across copies. It is not that —
running the shared-reference flag both ways moves `spread_shape` by 1%, against the 18.8× that
motivated the concern.

**Two further consequences.** The unconditioned branch at f32 also flips **152 of 1024 outcome
labels** at the default `r_coll`, where at f64 it flips none — a discrete label changed by a
registration artefact, which no continuous-field check would see. And the f32 caveat that *does*
survive is different from the one on record: 2 pixels of 1024 carry `|dE/E| > 1` in the tail.
`error_ratio` flags them.

### "`deep interior` is a near-triple encounter" — **it is a binary collision**

The characterisation on record came from an unregularised run with drift ~37. Under Aarseth–Zare,
same initial condition, both implementations:

| | Rust | numpy reference |
|---|---|---|
| `d_min` | 2.2837794877e-5 | 2.2976014100e-5 |
| `|dE/E|` | 1.395286245e-7 | 1.393632170e-7 |
| reference switches | 2 | 2 |
| reaches `t = 13` | yes | yes |

Initial separations are 2.236, 1.414 and 3.0. Sweeping `r_coll` from `1e-4 R` up to `R` itself,
pairs (0,1) and (1,2) **never register at any threshold**; only (0,2) does, and already at
`1e-4 R`. There is no near-triple. It is a close binary approach with a distant third body —
which is precisely the case regularisation exists to handle, so the warning predates the method
that removes it.

The region is still the hardest one measured: at 256×256 it flags 9228 of 65536 pixels and needs
all three refinement passes. Hard is not the same as non-regularisable.

### Any quantity measured at `n = 64` — **read through §2's scatter**

Not a specific correction but a blanket one, and it is the most consequential item here. At
`n = 64` every fraction carries an interdecile scatter of 2.0–4.4× its own value. Fractions,
rates, and comparisons between them from that era should be re-derived, not adjusted. Medians of
well-behaved quantities are more defensible; tail quantiles and "worst observed" figures are not
defensible at all, and the `drift max` row shows they will have looked *stable* while being
wrong.

### "Event class fires ~4 time units earlier than terminal outcome" — **an artefact of the baseline**

The figure was measured against an **escape-based** terminal criterion, where the terminal label
genuinely lagged. Here the terminal arm is *collision at `r_coll`*, and a collision **is** the
tightest pair reaching threshold — so both fire at the same boundary by construction. Measured
lead on all 22 pixels where both fire: **exactly zero**.

The lead time was never a property of the event class. It was a property of what it was compared
against, and quoting it without the baseline makes it unreproducible.

**The justification that survives is stronger than the one it replaces.** The event class needs no
gate and is defined at every playhead: at `t = 13` it flags 165 pixels against the terminal
statistic's 22, strictly nested with none flagged by the terminal one alone; at `t_max = 8` it
flags 110 against **0**, because nothing has terminated yet. Coverage and horizon-independence,
not earliness.

---

## 8. The vertical slice — the screen floor, and what it changes

Every build before this one was isolated, and each isolation hid something the next found. This
one put the camera, the screen floor, the adaptive render, SSAA and the linearised decoder
together. **The screen floor is the change that matters**, and it revises §5.

### 8.1 The arithmetic, and why four of §5's answers had to be re-run

One sample, one tile, no interpolation, so a fully-refined tree at level `L` holds `4^L x N²`
samples. At `N = 8` on a 512² viewport that reaches 262144 at **level 6**. §5's descent reached
**level 12** — 4096x past the point where samples stop being displayable.

So §5's q1, q2, q3 and q7 measured the criterion **minus its principal stop condition**. They are
not wrong; they describe a regime the real system never enters.

**A structural cap, stated before the numbers:** under the veto a tree cannot exceed `4^6 = 4096`
leaves. near-field's 4617 cannot recur, and its absence is arithmetic, not improvement.

### 8.2 q1 and q2 re-run: in near-field the criterion decides 39% of leaves

`results/output/sched_screen.txt`, budget 50000 quads, `tau = 1e-4`, `alpha_hi = 0.2`, N = 8,
E+1 = 8, `t = 13`, f64, viewport 512², camera framing the root box.

| region | veto | quads | leaves | depth | floor | keep | screen | wall s |
|---|---|---|---|---|---|---|---|---|
| far | off | 21 | 16 | 2 | 0 | 16 | 0 | 0.2 |
| far | **on** | 21 | 16 | 2 | 0 | 16 | 0 | 0.2 |
| near-field | off | 4617 | 3463 | 12 | 609 | 2854 | 0 | 236.3 |
| near-field | **on** | 549 | 412 | 6 | 100 | 60 | **252** | 27.1 |
| deep interior | off | 29 | 22 | 4 | 9 | 13 | 0 | 0.4 |
| deep interior | **on** | 29 | 22 | 4 | 9 | 13 | 0 | 0.5 |

Three separate results:

**near-field: the view stops 252 of 412 leaves — 61.2%.** The criterion decides the other 39%.
§5's "terminates at 4617 quads" describes a descent the real system stops at 549, an **8.4x**
reduction in compute for a tree that is by construction all anyone can see. Depth 12 becomes 6.

**deep interior is byte-identical with the veto on and off.** Its tree never reaches level 6, so
the veto never fires. **The screen floor does not fix q3's failure** — the tree still leaves the
largest high-spread structures shallow, and the veto is not the explanation or the remedy.

**far is unchanged, which is the control working.** It stops at level 2 either way; a region that
never descends cannot be vetoed, and a difference there would have meant the veto was firing on
something other than tile size.

`relcap` is zero everywhere: `MAX_REL_DEPTH = 6` coincides with the screen floor at camera depth
0, and the floor is checked first. That is the contract's `MAX_REL_DEPTH <= screen floor`, visible
as a column rather than asserted.

### 8.3 q7 inverts: `tau` is not inert, and it is now the dominant knob

`results/output/sweep_screen.txt`, budget 2000 quads, against §5's `sched_sweep.txt`.

§5 concluded: *"`alpha_hi` from 0.20 to 0.50 collapses the tree 80x, while `tau` is inert over
four orders."* Under the veto, at `tau = 1e-4`:

| region | `alpha_hi` 0.20 -> 0.50 | §5 | `tau` sensitivity | §5 |
|---|---|---|---|---|
| far | **x1.00** (no effect at all) | x80 | 1e-8 -> 1e-6: **x64** | identical |
| near-field | **x21.7** | x80 | 1e-3 -> 1e-2: **x16** | identical |
| deep interior | x1.16 at 1e-4 | x80 | 1e-6 -> 1e-4: **x7** | identical |

Both halves of §5's answer move, and they move in opposite directions.

**The mechanism is depth, and it is not subtle.** `alpha` is a *rate* statistic —
`log2(spread_parent / spread_child)` — and it needs levels to express itself. With
`bootstrap_levels = 2` and a floor at level 6, the descent has **four discretionary levels**;
`alpha_hi`'s 80x collapse was measured over twelve. `tau` is a *level* statistic, compared against
the spread directly, and it has exactly the same room it always had. Truncating the tree therefore
demotes `alpha` and promotes `tau`.

**So "sweep both before quoting any tree" survives, but its emphasis inverts.** §5 said `alpha_hi`
does more work than the criterion and `tau` is often inert. Under the veto `tau` decides whether a
region descends at all — `far` goes from 1024 leaves at `tau = 1e-8` to 16 at `1e-6`, a 64x swing
on a region whose whole point is that it is tame — while `alpha_hi` does nothing there whatsoever.

### 8.4 The adaptive render, and the test that rejects the old instrument

§5's overlay drew leaf boundaries over a **uniform** render, so every texel was the same size. It
showed where boundaries fell, not what the system displays, and the tree could not be judged by
eye at all.

`output::adaptive` rasterises each leaf's `N²` samples across its own screen footprint: a level-3
leaf's texels are 4x the linear size of a level-5 leaf's, one sample one tile, no interpolation. A
coarse quad is never upsampled smoothly, because that fabricates structure — the one thing a chaos
instrument's picture must not do.

`texel_scaling` fits `log2(texel_px)` against level. **The assertion is that it is exactly -1, and
the same assertion must reject a uniform render** — which fits 0, because every texel is the same
size whatever the level. §5's failure now has a test that fires on it
(`tests/vertical_slice.rs::texel_size_varies_as_two_to_the_minus_level_and_uniform_is_rejected`).

A known geometric cost, stated rather than hidden: `Slice::axis` is endpoint-inclusive, so a
sample-centred tile overhangs its quad by half a cell. Leaves are painted coarsest-first so a
finer neighbour overwrites the overhang. It is the same endpoint-inclusive duplication already
recorded at sibling edges (1/N of a quad, 12.5% at N = 8), seen from the render side.

### 8.5 The `E` sweep refutes the prediction, and the veto would have hidden it

The prediction on record before the run: *low `E` biases toward refine, exactly as low `N` did* —
`N = 4` spent 4x the quads of `N = 16`. **It does the opposite.**

`results/output/e_sweep.txt`, budget 6000 quads, `tau = 1e-4`, `alpha_hi = 0.2`, N = 8, `t = 13`,
f64. near-field, veto **off**:

| E+1 | quads | leaves | depth | trajectories | sibling range | capped? |
|---|---|---|---|---|---|---|
| 2 | 989 | 742 | 10 | 94,976 | 0.8997 | no |
| 4 | 3617 | 2713 | 10 | 694,528 | 0.9852 | no |
| 8 | 4617 | 3463 | 12 | 1,773,056 | 1.2403 | no |
| 16 | 5997 | 4498 | 10 | 4,605,952 | 1.4259 | **yes** |
| 32 | 5997 | 4498 | 10 | 9,211,904 | 1.4085 | **yes** |

Leaf count **rises** with `E`, monotonically, over the three uncapped rows: 742 -> 2713 -> 3463.
**Low `E` under-refines.** The last two rows are budget-limited and their leaf counts are floors,
not measurements; they are excluded from the trend rather than read as a plateau.

**And the veto would have reported a null.** The same sweep with the screen floor on gives 535,
523, 412, 607, 520 — no trend, everything inside a factor of 1.5, because a veto-capped tree
saturates before over- or under-refinement can express itself. **Running only the veto-on rows
would have concluded that `E` does not matter.** That confound was written into the plan before
the run for exactly this reason.

`far` is flat at 16 leaves across every `E` and both settings — the control working: `E` cannot
matter in a region that never descends. `deep interior` is erratic by an order (16, 127, 22, 31,
187 with the veto) and is **not** a trend; its `E+1 = 32` veto-off row explodes to a
budget-capped 4498 leaves. Chaotic scatter, reported as scatter.

### 8.6 `N` and `E` fail in opposite directions, and it is a bias not a noise

`results/output/spread_bias_e.txt` measures the mechanism directly, on identical footprints with
only `E+1` varying. Because the Halton offsets are a fixed prefix, the `E+1 = 2` copies are a
**subset** of the `E+1 = 32` copies, so any movement is the estimator and not a different sample.

`spread_shape` is the mean distance of the copies' `shape_vec` from their centroid. With two
points the centroid sits exactly between them and the statistic measures half of one pair's
separation; it is systematically **smaller** than the same quantity over eight copies. A small
spread falls below `tau_display`, and the quad is **kept**.

That is the opposite failure direction from `N`, and the two are not interchangeable:

| knob | what it resolves | undersampling it | direction |
|---|---|---|---|
| `N` | how well a quad knows its own **area** | inflates the between-footprint variation that drives `alpha` | **over**-refines |
| `E` | how well a footprint knows its own **value** | deflates the within-footprint spread compared against `tau` | **under**-refines |

The `sibling range` column reads the same way: it **rises** with `E` (0.90 -> 0.99 -> 1.24 -> 1.43
without the veto). That is not noise falling, it is signal appearing — more copies find more
genuine disagreement inside the same cell.

**The consequence for the tier design is the reverse of the one predicted.** The cheap tier does
not cancel itself: `E+1 = 2` costs 94,976 trajectories against `E+1 = 8`'s 1,773,056, an **18.7x**
saving for a 4x reduction in copies. It buys that by **refining less** — by not seeing structure —
not by being efficient. The risk is silent under-resolution, not a blown budget.

### 8.7 The estimator bias, measured directly

`results/output/spread_bias_e.txt`. Uniform 48x48 renders, identical nominal footprints, only
`E+1` varying. The Halton offsets are a **fixed prefix**, so the `E+1 = 2` copies are a strict
subset of the `E+1 = 32` copies — any movement is the estimator and not a different sample.

Median `ensemble_spread`, as a fraction of the `E+1 = 32` value:

| region | E+1=2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| near-field | **0.539** | 0.806 | 0.882 | 0.965 | 1.000 |
| deep interior | **0.558** | 0.835 | 0.890 | 0.992 | 1.000 |
| far | **0.131** | 0.438 | 0.801 | 0.923 | 1.000 |

Monotone in every region, and the direction is unambiguous: **two copies report about half the
spread that thirty-two do, and in a tame region only 13% of it.** That is the under-refinement in
§8.5, arriving as a number rather than an inference.

The **p10** column moves further than the median — near-field 1.907e-4 -> 9.183e-4, a 4.8x rise
against the median's 1.9x. The bias is worst in the *low tail*, which is exactly the population
sitting near `tau_display` where the keep-or-split decision is made. So the effect on the tree is
larger than the median shift alone suggests.

`far` is the extreme case and it makes sense: where the copies barely diverge, a two-point spread
has almost nothing to measure. A cheap tier is therefore least trustworthy precisely in the
regions it would be assigned to.

### 8.8 SSAA resolve, and the zoom ladder

**Resolve** (`results/output/ssaa_resolve.txt`, 256² uniform renders). The `E+1 = 1` row is the
control and it is exact: `moved frac 0.0000`, because the resolve of one copy *is* the nominal
copy. A nonzero value there would have been a bug in the resolve rather than a finding.

| region | E+1=2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| near-field, fraction of pixels whose colour moved | 0.0026 | 0.0054 | 0.0064 | 0.0068 | 0.0075 |
| deep interior | 0.0514 | 0.0972 | 0.1141 | 0.1264 | — |

Resolve changes **0.75%** of near-field's pixels and **12.6%** of `deep interior`'s, with a worst
per-pixel `|dRGB|` of 172 and 195 respectively. It saturates by about `E+1 = 8`. In `far`
**nothing moves at all** — every copy lands on the same outcome label, so the resolved image is
the nominal image exactly.

So the ensemble's second job is real but small in the tame and mixed regions and substantial in
the chaotic one. Spread and resolve are measuring different things on the same data and neither
substitutes for the other: `far` has a nonzero spread (2.5e-8) and a *zero* resolve difference.

**The zoom ladder** (`results/output/zoom_near-field.txt`, and `zoom_near-field_animated.png` with nine
still frames beside it). Each frame re-descends from a root box of `half = 0.05 / 2^k` with the
camera framing it, so `camera_depth` is 0 throughout and the screen floor always sits six levels
below the view.

| frame | half | leaves | depth | screen-floored | spread median |
|---|---|---|---|---|---|
| 0 | 5.00e-2 | 412 | 6 | 252 | 9.747e-4 |
| 1 | 2.50e-2 | 586 | 6 | 468 | 3.023e-4 |
| 2 | 1.25e-2 | 1045 | 6 | 868 | 9.568e-5 |
| 3 | 6.25e-3 | 811 | 6 | 484 | 3.197e-5 |
| 4 | 3.13e-3 | 214 | 6 | 16 | 2.773e-5 |
| 5 | 1.56e-3 | 73 | 4 | 0 | 3.960e-5 |
| 8 | 1.95e-4 | 16 | 2 | 0 | 4.813e-6 |

**This is the view-relative property, and a still image cannot show it.** The 252 quads floored in
frame 0 are refined in frames 1-3 with genuinely new samples, and a fresh population is floored in
their place — 868 at frame 2. Nothing is cached, nothing is upsampled. Past frame 4 the
neighbourhood is smooth enough that the criterion stops before the view does, and the screen
column falls to zero: the veto is a veto, and where the criterion is satisfied first it never
fires.

### 8.9 The two open items from the PR #11 review

`results/output/open_items.txt`.

**Item 1 — does `floored` correlate with `worst_energy_drift`?** It depends on the region, and the
answer for `deep interior` is yes.

| region | veto | leaves | floored | Spearman | drift median, floored | drift median, kept |
|---|---|---|---|---|---|---|
| deep interior | off | 22 | 9 | **+0.2987** | **3.256e-1** | 2.971e-3 |
| deep interior | on | 22 | 9 | +0.2987 | 3.256e-1 | 2.971e-3 |
| near-field | off | 3463 | 609 | +0.1300 | 4.085e-9 | 3.266e-9 |
| near-field | on | 412 | 100 | +0.0716 | 4.124e-9 | 5.068e-9 |
| far | either | 16 | 0 | n/a | n/a | 4.501e-11 |

In `deep interior` the floored quads carry a median drift of **3.3e-1** against **3.0e-3** in the
kept ones — two orders, on quads whose energy is wrong by 33%. Those quads are not trustworthy at
all, so at least part of that region's floor is **integration error corrupting `alpha`**, which is
a different bug with a different fix: BRIEF §2.5's flag-and-re-integrate remedy, not a scheduler
change. In near-field the two populations are indistinguishable (4.1e-9 against 3.3e-9), so the
floor there is not an integration artefact.

**And the caveat the review asked for bites, exactly where it was predicted to.** near-field's
floored quads live at levels 2-12; the veto removes everything below level 6, which is **509 of
609** of them — 84% of the population. The weaker Spearman under the veto (+0.0716) therefore
means *"no data"*, not *"no relationship"*. The `floored levels` column is printed beside every
row so this is readable rather than inferable.

`deep interior` is unaffected: its floored quads live at levels 2-4 and the veto never fires
there, so its correlation is measured on the same population either way.

**Item 2 — does p90 aggregation fix `deep interior`'s tree?** Yes. It is aggregation. **And the
fix collides with the screen floor, mildly for p90 and severely for mean.**

Without a camera:

| agg | quads | leaves | depth | depth histogram |
|---|---|---|---|---|
| median | 29 | 22 | 4 | 2:15 3:3 4:4 |
| mean | 161 | 121 | **7** | 2:10 3:17 4:24 5:12 6:2 7:56 |
| p90 | 81 | 61 | **7** | 2:11 3:15 4:17 5:11 6:3 7:4 |

Median stalls at depth 4 with fifteen of its twenty-two leaves still at level 2. Both mean and p90
descend to **depth 7**. So §5's q3 failure is **median blindness** — a thin filament crossing a
quad does not move the median of 64 footprint spreads — and the attribution stated there as
"plausible and unproven" is now measured. It is *not* p90 specifically: mean descends further
still.

**But depth 7 is one level past what a 512² viewport can display**, so the fix and the veto want
different things in exactly the configuration production runs. `results/output/agg_vs_floor.txt`
measures the cost:

| viewport | `MAX_REL_DEPTH` | agg | leaves | depth | screen | relcap | histogram |
|---|---|---|---|---|---|---|---|
| 512² | 6 | mean | **79** | 6 | 16 | 0 | 2:10 3:17 4:24 5:12 6:16 |
| 512² | 6 | p90 | **58** | 6 | 4 | 0 | 2:11 3:15 4:17 5:11 6:4 |
| 1024² | **7** | mean | 121 | 7 | 0 | 0 | 2:10 3:17 4:24 5:12 6:2 7:56 |
| 1024² | **7** | p90 | 61 | 7 | 4 | 0 | 2:11 3:15 4:17 5:11 6:3 7:4 |
| 2048² | 8 | either | 121 / 61 | 7 | 0 | 0 | identical to no-camera |

**The truncation costs p90 three leaves of sixty-one (4.9%) and mean forty-two of a hundred and
twenty-one (34.7%).** Both keep an identical tree at levels 2–5; the difference is entirely in what
piles up at the cap. So p90's fix survives the production viewport nearly intact and mean's does
not — a ground for preferring p90 that has nothing to do with the statistic itself.

At **1024²** — where `4^7 x 64 = 1,048,576` makes level 7 displayable — the fix is fully realised
for mean and all but four leaves for p90. The collision is therefore a *resolution* limit, not a
design conflict: the aggregation the region needs is affordable one viewport step up.

**The two regions answer this oppositely, and that is the more useful result.** `deep interior`'s
descent is **criterion-bound**: the veto touches 4 of p90's leaves at 512² and none at 2048², so
raising the viewport hands the region straight back to the criterion. near-field's is
**view-bound at every viewport tested** — at 1024² with `MAX_REL_DEPTH = 7` it still floors 576 of
median's 844 leaves, 756 of mean's 988 and 88 of p90's 271, and at **2048² with `MAX_REL_DEPTH = 8`**
it still floors 2172 of mean's 2617 and 148 of p90's 382. Uncapped, p90 there reaches **depth 14**.
Its structure is dense at every scale, so more pixels buy more tree without ever reaching the point
where the criterion decides.

So "does the aggregation fix collide with the floor" has no single answer: in the region the fix
was *for*, it does not, at one viewport step up. In near-field the question does not arise,
because the criterion was never what was stopping the descent.

**And a trap the same table exposes.** At 1024² with `MAX_REL_DEPTH` left at its default 6, the
tree is **identical to the 512² tree** and the `screen` column reads **zero** — which looks
exactly like "the viewport made no difference". It is not: `MAX_REL_DEPTH` had taken over as the
binding cap, and it is a *policy default*, not arithmetic. The two coincide at 512² and diverge
above it. **Report both cap columns or the viewport looks inert**; the first version of this run
showed only `screen` and would have been read that way.

### 8.10 The deep-zoom decoder: two results, both narrower than the claim

`results/output/decode_ladder.txt` and `deep_zoom.txt`. **Distinctness is read before divergence**
throughout: two paths that have both collapsed agree perfectly, and their agreement means nothing.
`decode::distinct` compares all twelve state components bitwise, so a collapse is counted exactly.

**Result 1 — the floor is a property of *where* you zoom, not of the renderer.** §5 recorded a
plain-f64 cell-width floor at level 45.87. That is conditional on the chart coordinate being of
order 1, and the condition was never stated. With 64 samples per quad:

| chart | centre | `direct_f64` holds 64/64 to | `direct_f32` to |
|---|---|---|---|
| body_plane | \|c\| ~ 3 | depth 44 | depth 14 |
| body_plane | **0** | **depth 55+ — no floor in range** | **depth 55+** |
| shape | \|c\| ~ 3 | depth 35 | depth 14 |
| shape | 0 | depth 45 | depth 45 |

A quad centred at the chart origin has no O(1) neighbour for the increment to be absorbed into, so
it has no cell-width floor at all in the tested range, **on either precision**. Quote the
coordinate magnitude with any floor depth or the number means nothing.

**Result 2 — the linearised decoder buys ~24 levels over f32 and none over f64.** On body_plane at
\|c\| ~ 3:

| path | 64/64 distinct to | collapsed to 1 by |
|---|---|---|
| `direct_f32` | depth 14 | depth 22 |
| `L-naive_f32` — the literal formula | depth 14 | depth 22 |
| `L-split_f32` | depth 44 | depth 50 |
| `direct_f64` | depth 44 | depth 50 |

The literal formula — `x0`, `J_D.delta` and the sum all in f32 — collapses on **exactly the same
curve** as forming the chart coordinate in f32 in the first place (56, 18, 2, 1 at depths 16, 18,
20, 22 for both). It buys nothing. The split form — `x0` in f64 on the CPU, `delta` and
`J_D.delta` in f32, promoted and summed in f64 — tracks `direct_f64` rung for rung. So the gain
is real but it is *for an f32 consumer*: it lets an f32 GPU reach the f64 CPU's floor, and does
not push past it. The contract's "~50+" is f64's floor, not something the linearisation creates.

That follows from a bound stated **before** the run: the initial conditions must be formed as
absolute O(1) numbers before integration, because the three-body separations are O(1) and no
nonlinear integrator can carry `(x0, delta)` separately through the march. The linearisation
escapes the *chart-coordinate* floor and not the *IC-magnitude* one.

One modest exception: on the **nonlinear** chart at \|c\| ~ 3, `L-split` holds 64/64 to depth 45
where `direct_f64` has fallen to 10/64 — worth about six levels, from a single fused affine step
losing fewer bits than a decode through `cos`, `sin` and a renormalisation. Conditioning, not the
design's stated mechanism.

**Result 3 — the linearisation matters at the coarse end, not the deep end.** On the shape chart
`|direct - linearised|` against the sample spacing runs 0.39 at `half = 0.05`, 1.5e-3 at depth 8,
3.6e-7 at depth 20: it *falls* as the box shrinks, because the discarded term is `O(h²)` against a
spacing of `O(h)`. It exceeds one sample spacing only at `half >= 0.5`, larger than anything this
project renders. The approximation is worst where it is least needed and best where it is used —
the opposite of the intuition that a linearisation breaks down at depth. On `body_plane` the same
column is **structurally zero** and is reported as structural, never as a measurement.

**Result 4 — a collapsed decode makes the criterion maximally confident, and the collapse arrives
from the leaves upward.** `deep_zoom.txt` runs the descent under each path, counting collapsed
quads exactly:

| camera depth | path | root distinct | quads collapsed | spread of collapsed | tree |
|---|---|---|---|---|---|
| 14 | direct_f32 | **64/64** | **16 of 21** | 1.811e-7 | 21 quads, 16 leaves, depth 2 |
| 18 | direct_f32 | 18/64 | 21 of 21 | **5.551e-17** | 21 quads, 16 leaves, depth 2 |
| 22 | direct_f32 | 1/64 | 21 of 21 | 5.551e-17 | 21 quads, 16 leaves, depth 2 |
| 30 | lin_split_f32 | 64/64 | **0 of 21** | — | 21 quads, 16 leaves, depth 2 |

Two things there, and the first is the reason a root-level check is not enough. **At depth 14 the
root quad still resolves all 64 of its samples while 16 of its 21 descendants have collapsed** —
the children sit at half the cell width, so the failure begins at the leaves and works upward. A
distinctness check at the root would have reported everything fine.

And at depth 18 onward the tree is *21 quads, 16 leaves, depth 2* — the same shape `far` produces,
which is the tamest region in the study. The spread the criterion saw was **5.551e-17**, twelve
orders below `tau = 1e-4`. Nothing downstream can tell that apart from a perfectly resolved
region. This is the project's standing pattern from a new direction — *a statistic can report
maximum confidence precisely when it is least informed* — and it is why collapse is counted
exactly rather than thresholded.

The `1.811e-7` row is the dangerous middle: a partly-collapsed quad returns a small but not
absurd number, which no sanity check would flag.

**The precision floor is a separate limit and the dump keeps them apart.** At camera depth 40 the
root quad's own cell width is already below `PRECISION_MARGIN * f64::EPSILON`, so the descent
stops with `root decision = precision_floor` before any decode path is exercised. A numerical stop,
not a physical one, and not evidence about either.

### 8.11 Slice variety: are the prior conclusions slice-conditional?

`results/output/slice_variety.txt`. **The comparison holds the configuration fixed and rotates the
plane through it** — all cases share near-field's centre (body 0 at `(1, 3)`, released from rest)
with bases orthonormal in the 6D position metric, so a unit of chart coordinate moves the system
the same distance in every row.

A first version of this varied the chart at the same *coordinates* `(1.0, 3.0)`, which is not an
orientation test: an oblique plane evaluated there lands on a completely different configuration,
and the tamer spreads it reported were about the neighbourhood, not the angle. It was rewritten
rather than reinterpreted. The `plane 0deg` row is the control and the check is on the **initial
conditions**, not on the tree — the tree is downstream of a chaotic integration, so checking it
would be testing chaos rather than the charts. Measured: `max |dIC| = 0.000e0`, and the two trees
are identical quad for quad.

| case | leaves | depth | screen | spread median | alpha p10 | alpha median | alpha p90 |
|---|---|---|---|---|---|---|---|
| body_plane (control) | 412 | 6 | 252 | 9.747e-4 | -0.123 | 0.216 | 0.879 |
| plane 0deg (control) | 412 | 6 | 252 | 9.747e-4 | -0.123 | 0.216 | 0.879 |
| plane 15deg | 403 | 6 | 252 | 8.638e-4 | -0.132 | 0.234 | 1.102 |
| plane 30deg | 526 | 6 | 376 | 6.824e-4 | -0.125 | 0.289 | 1.257 |
| plane 45deg | 439 | 6 | 300 | 6.272e-4 | -0.092 | 0.287 | 1.213 |
| plane 45deg, half into body 2 | **226** | 6 | 128 | 1.280e-3 | -0.155 | 0.172 | 0.509 |
| plane 45deg, all into body 2 | 355 | 6 | 220 | 1.013e-3 | -0.109 | 0.211 | 0.629 |
| shape (nonlinear) | **970** | 6 | 768 | 1.700e-3 | -0.170 | 0.286 | 0.916 |

**Yes, the tree is slice-conditional — and no, the criterion is not, or much less so.**

Leaf count spans **226 to 970, a factor of 4.3**, at one fixed centre configuration with nothing
varying but the 2-plane through it. Median spread spans 6.3e-4 to 1.7e-3, a factor of 2.7. Among
the pure rotations alone (0 to 45 degrees, single body) it is 403 to 526, a factor of 1.3 — so
most of the variation comes from **which bodies the plane moves**, not from the angle within one
body's plane.

But the `alpha` distribution barely moves: median 0.172 to 0.289, p10 between -0.09 and -0.17, p90
between 0.51 and 1.26 across every case including the nonlinear chart. So the *exponent* the
criterion reads is a much more stable quantity than the tree it produces. That is a mildly
reassuring result for the criterion and a cautionary one for every leaf count quoted anywhere in
this document: **each is conditional on the slice, to about a factor of 4.**

**And an unplanned gauge check that could have failed.** The three `shape phase` rows — 0.0, 0.4,
1.3 — are **bitwise identical** in every column. The fibre phase is a global rotation of the
configuration and the three-body problem is rotation-invariant, so they must be. If the Hopf
inverse or the AZ port had broken rotational invariance, these rows would have separated. They
are the last row of `slice_variety.txt` and they cost nothing to keep.



## 10. Improving the refinement criterion

> **SUPERSEDED BY §11.** Every number in this section was measured under a colouring whose
> lightness ramp was linear over a window an order of magnitude too wide, and whose hue map was
> 2-to-1 in `n0`. The standing rule is *choose a criterion under the colouring that will ship*,
> so these tables score the criteria against a target that no longer exists. They are kept
> because what a criterion looks like through a broken instrument is itself a finding, and
> because §11's conclusions are stated as differences from these.
>
> The one conclusion that survives unchanged: `within/median` is beaten by random. §11 finds it
> beaten under the new colouring too, and finds `between/median` beaten as well — which is the
> part §10 got wrong.

Read §10.2 before anything else in this section. It is the only measurement here that says
whether a criterion is *good*, and every other result is judged against it.

### 10.1 The criterion does read a different quantity — but not for the reason proposed

`ensemble_spread` is a statistic over the `E+1` copies of one footprint; `reduce` aggregates
those `N²` numbers. The brief reads that as a category error, on the grounds that refinement
buys more footprints so only *between*-footprint variation is reducible.

**The premise does not describe this implementation.** `jitter_frac` is 0.5 and `halton_offset`
returns `[-1, 1)²` scaled by cell width, so the copies span the **whole cell, edge to edge** —
a quasi-random sample of exactly the area the footprint stands for. And the Halton control's
true `alpha` is exactly **1.0**, which an irreducible within-point statistic cannot be, since
splitting would not shrink it.

Measured, with sample counts matched so scale and count are separated:

```
        region  quads  coll  dead   rho all   rho mix   rho bnd     d90    dmax  med hot  n mix  n bnd
    near-field    549     0     0    0.7240    0.5818       NaN  0.3175  0.9872    1.000    192      2
           far     21     0     0    1.0000       NaN       NaN  0.0000  0.0000    0.000      0      0
 deep interior     29     0     0    0.6828    0.6361   -0.4000  0.4643  0.5714    0.438     27      4

        region   med within  med between        scale        count
                 (cell,E+1)   (quad,N^2)    matched/w     pooled/b
    near-field    1.0079e-3    2.0881e-3       1.1716       1.0127
           far    4.2680e-8    4.4921e-7       9.5550       1.0033
 deep interior    9.4462e-5    2.7824e-4       2.0617       1.0118
```

**`count` is 1.01 everywhere.** At equal extent and equal sample count the two arms are the
*same estimator* — not different quantities. **`scale` runs 1.17 to 9.56.** Widening the window
from a cell to a quad buys a factor of ~7 in a smooth region (the quad is `N-1` cells wide) and
almost nothing in a saturated one.

**But `rho mix` is 0.58–0.64** on the quads that contain a transition, so at their actual
settings the arms rank quads materially differently. §1's conclusion stands; its mechanism does
not, and the fix moves with it — the fault is the **aggregation**, which is what §10.4's
candidates address.

Read `rho mix`, never `rho all`: an all-quads correlation is dominated by tame quads where both
arms read near zero and agree trivially. `rho bnd` has a population of 2, 0 and 4 and is not
data — `med hot = 1.000` says why: in a chaotic region most quads are *uniformly* hot, so there
is no internal hot/cold edge and they are correctly not boundaries.

### 10.2 A metric to judge criteria by

```
reference = the fully-refined tree at the screen floor, one sample per pixel
IMAGE(B)  = the tree a ranking builds under budget B, rendered at true per-quad texel sizes
error(B)  = mean per-pixel OKLab distance between them
```

One integration pass per region (5461 quads, 2,796,032 trajectories, ~4 min) builds the complete
tree; every criterion, both controls and the whole curve are then replays over the cache with no
re-integration. Two facts make that exact: quads are disjoint, so each quad's error contribution
is a **constant** and the greedy replay is a static priority queue; and the reference colouring
is the nominal copy's outcome, which is **`E`-independent** because copy 0 is never jittered.

**`error = 0` means "matches this sampling", not "correct".** The reference is one finite
sampling; at the screen floor sub-pixel structure is sampled arbitrarily, and which side of a
filament a pixel lands on is an accident of where its sample fell. The exactly-locatable zero is
what makes the curve comparable between criteria; it is not image quality.

**`greedy_oracle` is a strong reference, not a ceiling.** Greedy is optimal only when gains are
independent and immediately available, and on a tree they are neither — a quad whose own split
gains little may unlock children two levels down, and greedy declines it. **A criterion beating
it indicates lookahead value, not an error**, and no test asserts it dominates.

**Criteria enter as orderings, never against `tau`.** §10.1's `scale` factor means a threshold
comparison would have scored the 1.17-vs-9.56 rescaling instead of the signal.

### 10.3 The shipped criterion is the worst one tested

`deep interior`, `t = 13`, N=8, E+1=8, budget in quads computed:

```
                   B =         5        11        23        47        95       191       383       767      1535      3071
         greedy_oracle   0.02715   0.02368   0.02038   0.01788   0.01626   0.01386   0.00874   0.00240   0.00004   0.00004
frac_hot_between/median   0.02715   0.02632   0.02155   0.01843   0.01696   0.01446   0.01199   0.00366   0.00002   0.00000
    running_max/median   0.02715   0.02368   0.02302   0.02274   0.02221   0.01786   0.01226   0.00425   0.00158   0.00080
        between/median   0.02715   0.02368   0.02282   0.02120   0.02037   0.01882   0.01478   0.00560   0.00003   0.00000
    max_of_both/median   0.02715   0.02368   0.02282   0.02120   0.02037   0.01885   0.01509   0.00560   0.00003   0.00000
contrast:between/median   0.02715   0.02715   0.02282   0.01973   0.01907   0.01812   0.01510   0.00975   0.00312   0.00000
           within/mean   0.02715   0.02632   0.02607   0.02532   0.02388   0.02118   0.01593   0.00509   0.00001   0.00000
            within/p90   0.02715   0.02632   0.02613   0.02536   0.02422   0.02143   0.01592   0.00489   0.00029   0.00006
         layout/median   0.02715   0.02368   0.02308   0.02281   0.02233   0.01966   0.01681   0.01505   0.00308   0.00000
      first_div/median   0.02715   0.02632   0.02610   0.02536   0.02384   0.02125   0.01593   0.00631   0.00467   0.00304
frac_hot_within/median   0.02715   0.02368   0.02308   0.02281   0.02233   0.01971   0.01681   0.01502   0.01076   0.00000
      term_grad/median   0.02715   0.02632   0.02603   0.02548   0.02401   0.02128   0.01596   0.00824   0.00521   0.00200
         within/median   0.02715   0.02368   0.02293   0.02242   0.01984   0.01861   0.01723   0.01509   0.01386   0.00032
             random lo   0.02715   0.02368   0.02114   0.02026   0.01881   0.01814   0.01580   0.01288   0.00932   0.00392
             random hi   0.02715   0.02715   0.02632   0.02632   0.02613   0.02452   0.02327   0.02132   0.02012   0.01579
```

**`within/median` — the shipped default — is beaten by random at every budget past 383**, in
both regions.

**And at `t = 20` a criterion beats the oracle, exactly as anticipated when it was named.** In
near-field `greedy_oracle` plateaus at **0.00048** from `B = 383` through `B = 3071`, while
`first_div` reaches **0.00000 at `B = 1535`**. Greedy declined a low-gain split that unlocked
children two levels down — the classic failure of greedy on a sequential tree problem. This is
why the control is called `greedy_oracle`, documented as a strong reference rather than a
ceiling, and why **no test asserts that it dominates**: such a test would have fired here, on
entirely correct behaviour. In `near-field` it is flat at 0.00394 to `B = 767` while `within/mean`,
`between/median` and `max_of_both` all reach the oracle's zero at `B = 191`.

**`far` cannot be measured**: `error(root) = 0.00000`. The outcome image is featureless at 512²,
so every criterion reads zero and none of it is data. Reported as undefined rather than as
agreement — and it reframes every earlier leaf-count comparison on `far`, where there was never
an image to get right.

### 10.4 A flat curve has two causes, and `error(B)` cannot tell them apart

I expected `within/median`'s flatness to be a degenerate ranking. It is not.

```
                signal distinct   modal%    nan%    spread      (near-field, of 5461 quads)
         within/median     5418     0.3%    0.0%  4.285e-1
           within/mean     5461     0.0%    0.0%  3.638e-1
        between/median     5130     1.3%    0.0%  5.079e-1
frac_hot_within/median       58    40.8%    0.0%   1.000e0
frac_hot_between/median      31    83.1%    0.0%  9.844e-1
         layout/median       78    40.8%    0.0%   1.000e0
      term_grad/median      159    97.1%   97.1%  1.113e-1
    running_max/median     5427     0.1%    0.0%  4.294e-1
```

`within/median` has **5418 distinct values of 5461**: a fine-grained ordering that is actively
bad. `frac_hot_within` and `layout` have 58 and 78: **no ordering at all**, and their flat curve
is the tie-break's scan order. Two different faults with different fixes, which is why the
distinct-value count is printed *above* the curves.

Two rows resist a tidy story. **`frac_hot_between` is the best criterion in `deep interior` on
65 distinct values**, beating a 4994-valued one — resolution is not what makes a ranking good.
And **`term_grad` is NaN on 97.1% of near-field** yet reaches the oracle's zero by `B = 383`,
because the 2.9% it scores are exactly the structured quads. A high `nan%` is a property to
read, not a defect to hide.

### 10.5 `alpha_sibling_spread`: usable as a signal, not at its shipped threshold

The obvious control could not fail. `sigma_E(0)` has true `alpha` exactly 1.0 and true sibling
range exactly 0 — and reads **0.003, flat in both `N` and `E+1`**. The flatness is the tell:
under the fixed Halton prefix the offsets *and* the footprints are fixed, so the whole quantity
is deterministic and there is no sampling noise in it. Kept as a geometric floor and labelled;
part 2 varies a `Pcg` seed for a real draw.

```
        region  seed parents     a p50    a idec   sib p50   sib p90  seed move
    near-field     0      21    0.2698    0.6605    0.4501    0.9466        NaN
    near-field     1      21    0.2602    0.8794    0.4193    0.9691     0.2641
    near-field     2      21    0.2991    0.7555    0.4826    0.9151     0.2091
 deep interior     0      21    0.1305    3.8809    0.7869   13.8295        NaN
 deep interior     1      21    0.0731    1.6555    0.9612   13.8885     0.2868
 deep interior     2      21    0.1462    2.7090    1.0514   13.8957     0.3636
```

Sampling noise alone (`seed move`, p90 of `|alpha(seed k) - alpha(seed 0)|`) is **0.21–0.36**
against `sib_tau = 0.5`. The sibling **median** is 0.45 in near-field and 0.79–1.05 in
`deep interior`: the shipped threshold sits inside the noise-broadened bulk in both, so a
typical quad's `Sibling` decision is close to a coin flip. It carries signal — `sib p90` reaches
13.9 in `deep interior`, far above the noise — but not at 0.5.

### 10.6 The FTLE port

`reference/tb_ftle.py` ported: Benettin shadow at `d0`, renormalised every 200 steps, `S`
accumulating `log(d/d0)`, `ftle = S/T`, plus the O(1) diffusion regression on `log(inertia)`.

Cross-checked against the live reference:

```
column               max abs       max rel  argmax row
ftle               8.882e-16     3.872e-16           0
diffusion          0.000e+00     0.000e+00          -1
dmin               0.000e+00     0.000e+00          -1
nre                0.000e+00     0.000e+00          -1
S                  1.776e-15     3.872e-16           0
rx0                0.000e+00     0.000e+00          -1
ry0                0.000e+00     0.000e+00          -1
PASS
```

Two limitations stated because they bound what follows. It sits on `tb.py`'s **unregularised**
fixed-step leapfrog, because that is the pair `tb_ftle.py` has a reference for — so an FTLE near
a close approach is not trustworthy, and `d_min` from the AZ march is the column that says
where. And the perturbation direction is pinned analytically on both sides: numpy's Ziggurat is
not ported, reproducing it is not required since the direction only seeds the shadow, and
**nothing about an RNG enters the validated path**.

`n_renorm` is asserted nonzero, not assumed. Renormalisation is what stops the estimator
saturating; without it `log(d/d0)/T` decays and reports `lambda ≈ 0` for the *most* chaotic
trajectories — the inversion this project has now met four times.

### 10.6b The coupling question: the criterion is largely blind to the lightness field

The production scheme is bivariate — hue from the shape sphere, lightness from a scalar.
`spread_shape` maps to hue, so that half is aligned by construction. §6 asks whether the
criterion is blind to the other half. Four colourings over **one** integration pass, so the only
thing that changes is the map from footprint to pixel.

near-field, `B = 341` of 1365, the gap between the **best criterion** and the **random** control:

```
colouring              random    best criterion            gap
outcome               0.00214   0.00000  between          -> total
bivariate/spread      0.01568   0.00238  between          -> 6.6x
bivariate/diffusion   0.02238   0.01449  frac_hot_between -> 1.5x
bivariate/ftle        0.03656   0.03044  frac_hot_between -> 1.2x
```

**The gap collapses as lightness moves away from what the criterion measures.** Under `outcome`
and under `spread` — the field the criterion actually reads — the best criterion beats random
outright. Under `diffusion` it is 1.5x better; under `ftle`, **1.2x**, which is barely better
than spending the budget at random.

The ordering changes too: `between/median` and `contrast:between` lead under `outcome` and
`spread`, and `frac_hot_between` takes over under `diffusion` and `ftle`. A criterion chosen on
the diagnostic colouring is not automatically the right one for the production colouring.

**§6's concern is real and this is the number for it.** If lightness ships carrying FTLE or
diffusion, the criterion needs a term for that field and has none — a quad can be uniform in
shape and structured in diffusion, and nothing would refine it.

One caveat bounding the `ftle` row specifically: that march is the **unregularised** fixed-step
leapfrog, because `tb_ftle.py` is built on `tb.py` and that is the pair with a reference. Near a
close approach it is not trustworthy, so the `ftle` row in `deep interior` carries that caveat
and the `spread` row does not.

### 10.7 Panning: three of §7's four questions are identities

`Camera::veto` reads `tile_size_px`, which depends on the quad's width and the camera's
`half_world` and `viewport` — **and not on `cx` or `cy`**. There is no view culling and no cache.
So panning changes no scheduling decision, and "the tree persists across a pan" is an identity
rather than a finding.

```
 step    cam cx   quads in view       new would-cache%  recompute%   floored
    0   0.95000     549     263       549         0.0%      100.0%       252
    1   0.96250     549     353         0       100.0%      100.0%       252
    4   1.00000     549     549         0       100.0%      100.0%       252
    8   1.05000     549     323         0       100.0%      100.0%       252

4392 quads recomputed after having been seen before.
0 came back with a DIFFERENT reduction. 0 with a different decision.
```

What is not an identity: **a quad recomputed after leaving view comes back bitwise identical**,
in reduction and in decision, across 4392 recomputes. That is the property a cache would need to
be sound, and it holds because the ensemble offsets are a fixed Halton prefix indexed by *copy* —
not by pixel, not by camera, not by time — so a quad's ensemble does not know the camera exists.

`in view` runs 263 → 549 → 323 as the camera crosses: at the extremes **half the tree is off
screen**, and that is what a cache would be sized against. `Camera::covers` computes it and the
scheduler never consults it. Making it consult it would be view culling, and putting a position
term into the screen floor would make a quad's decision depend on where the camera points —
which is what "never cached as a quad fact" exists to keep out.

### 10.8 Cost-aware priority is not worth building; anisotropy is narrower than it looks

```
region          p10       p50       p90       p99       max    p99/p50
near-field   1876649   1981003   2041477   2046765   2050044     1.03
far           321024    321024    321024    321024    321024     1.00
deep int      496253    645606    953387   1723569   1951532     2.67
```

The cost distribution is **narrow**: 1.03 in near-field, 1.00 in `far`, 2.67 in `deep interior`.
And ranking by `Δerror / cost` does not beat ranking by `Δerror` — it is identical in near-field
and marginally *worse* in `deep interior` (0.01483 against 0.01408 at `B = 85`). **Cost-aware
priority has nothing to move here.** Reported as a negative and not built.

Anisotropy, at three `tau` because the answer is a strong function of it:

```
                 >=3 children keep    all four keep    the anisotropy case
near-field  1e-4          0.9%              0.3%              0.6%
near-field  1e-3         33.1%             24.3%              8.8%
deep int    1e-4         66.3%             60.7%              5.6%
deep int    1e-3         81.2%             73.6%              7.6%
```

**The headline number is the wrong one.** 60.7% of splits in `deep interior` producing four
children that all immediately keep is an argument that *the split should not have happened* —
i.e. about the criterion, which §10.3 already says is the problem. The case anisotropic
splitting actually addresses is **exactly three keep**: one useful child, three wasted, where a
2-way split along the disagreement direction would have captured it. That band is **5.6–8.8%**.

So anisotropy is worth roughly a twentieth of the quads, not two thirds. Costing only; nothing
implemented.

### 10.9 Looking at the slices, and at the trees

`slice_variety` measured tree size as **slice-conditional to 4.3x** while the `alpha` distribution
stayed put. That is a table, and it cannot say whether an oblique plane cuts the same structure
at an angle or lands on different structure entirely. `slice_gallery` renders ten charts through
**one shared centre configuration** — the axis-aligned body plane, oblique planes at 15/30/45deg,
two cross-body mixes, and the shape chart at three fibre phases — with orthonormal bases in the
6D position metric, so a unit of chart coordinate moves the system equally far in each.

**The control pair caught a real error on the first run.** `body_plane` reads its centre from the
chart coordinate; `Plane` and `Shape` carry it in `origin` and must be centred at zero. Centring
them all at `(1, 3)` samples a box two units away from the shared configuration — a different
slice of different physics rather than a rotation of the same one. It reported 549 quads against
21 and was caught by the assertion that `plane_00deg` must be **bitwise** `body_plane`, which is
the same chart written twice. Without that pair the run would have produced a gallery of
plausible pictures of the wrong thing.

**Every image in this build now has a `_wire` twin.** The plain render says *what is displayed* —
texels at true per-quad sizes, so a coarse leaf is visibly coarse. The wire says *where the tree
cut*, brightness graded by level. They answer different questions: a coarse texel tells you a
leaf is coarse, and only the wire tells you whether the structure around it was subdivided
*around* it or straight *through* it. **PR #11 drew boundaries over a uniform base**, which
conflated the two, and that is how `deep interior`'s bad tree went unnoticed for a whole build.

The most useful artefact is `budget_<region>_animated.png`: `greedy_oracle` on the left and
`within/median` on the right, at the same budget, frame by frame.

**And the wire pair answers the question §10.3 could not.** At `B = 682` in near-field:

- `near-field_B682_within_median_wire.png` — the shipped default has spent almost its entire
  budget shredding the **top-left corner** into a fine mesh, while the collision region in the
  bottom-right sits untouched in **two enormous level-1 leaves**.
- `near-field_B682_greedy_oracle_wire.png` — the oracle refines along the collision region's
  boundary and the left edge, and leaves the uniform interior alone.

So `within/median` is not noisy, and it is not failing to order. **It is systematically
refining the wrong corner**, which is why it loses to random: random at least spreads its budget
evenly and hits the boundary by accident. That is the mechanism behind §10.4's finding that the
signal has 5418 distinct values and is still worse than a coin flip, and it is not visible in
any table — only in the wire.

## 11. The colouring, the charts, and the criterion re-measured

**§10 was measured under a colouring that could not see its own signal.** Its tables are kept
and marked, because what a criterion looks like through a broken instrument is itself a finding —
but every number in §10 scores the criteria against a target that no longer ships. The standing
rule is *choose a criterion under the colouring that will ship*.

### 11.1 The production colouring had three faults, and the one I named first was wrong

`bivariate::rgb` mapped hue as `chroma*(cos h, sin h)` with `h = atan2(n2,n1)` and
`chroma = C_MAX*hypot(n1,n2)`. I recorded that as a **seam**. It is not: that composition is
algebraically `C_MAX*(n1,n2)` — measured agreement `4.2e-17` over a sphere sweep — a linear
projection, continuous everywhere.

Its actual fault is that it is exactly **2-to-1**. It discards `n0 = (|rho~|^2 - |lam~|^2)/I`, so
`n` and its `n0 -> -n0` partner render **bitwise identically**: a tight binary with a distant
third body and a wide pair with a close third, painted the same colour.

**And that cost less than it sounds.** `n0` is *reached* end to end — span `1.9946` in
`near-field`, `1.9994` in `deep interior`, against a maximum of 2 — but its interdecile is only
`0.0665` and `0.1684`. Span is a max statistic. The merge bit in the tail, not the bulk.

**The flat images were the ramp.** Two separate causes, and only together do they explain the
picture:

| | value |
|---|---|
| `ensemble_spread` window, near-field | `(4.19e-5, 2.857e-1)` |
| the same field's *continuous arm* p99 | `2.244e-2` |
| ratio | **12.7x too wide** |
| `spread_event` distinct values | **5** (modal 98.2%) |
| footprints where the event arm dominates | **1.7%**, all in the top tail |

`ensemble_spread` is `max(spread_shape, spread_event)`. The event arm is a count ratio over
`E+1 = 8` copies, so it takes five values — and it sets the p99 while describing 1.7% of the
region. A **linear** ramp over a window an order of magnitude too wide, set by a staircase. The
committed `colour_near-field_bivariate_spread_reference.png` is a flat navy field.

The replacement (`src/output/colour.rs`) blends **six** vMF sites on the shape sphere, computed
from the run's own masses; the curve and the polarity are properties of the field
(`Scalar::curve`, `Scalar::direction`) rather than call-site arguments. The sixth site
(`lambda = 0`, body 2 at the inner barycentre) was added on a measurement, not for symmetry: with
five the worst angular gap was `1.193 rad` and it sat at `n0 = +1`. `kappa = 3` was chosen on
`hue_coverage` — near-field peaks there at `0.0058`, `deep interior` reaches 83% of its saturated
value.

### 11.2 The reserved null was reachable, and the metric was scoring it as correct

A NaN `shape_vec` — a triple collision, which `shape.rs` deliberately does not floor — went
through `oklab_to_srgb` into `NaN.round() as u8`, which saturates to `[0,0,0]`, and
`metric::Cache::render` filled background with `0u8`. **An undetermined pixel was bitwise
identical to un-rendered background, so `err_sum` scored it as a perfect match.**

Not theoretical: **18 of 65536** `deep interior` footprints at 256² have a non-finite shape.
`State::DecodeFailed` also had no palette arm and fell to the same grey as an invalid state byte.
Both now return `DEBUG_NAN`, and `BACKGROUND` is deliberately not black.

### 11.3 `far` is featureless in the data, and an auto-ranged ramp hides that

Predicted: a log ramp would rescue `far`, which §10 could not measure at all
(`error(root) = 0.00000`). **Refuted, then refuted again in the other direction.**

At 256²: `far` has `n0` span `0.0000`, `spread` p1/p99 differing by 2.6%, and hue coverage
**exactly `0.0000` at every kappa tested**. It is featureless in the data.

But under the shipping colouring at 1024² it reads `error(root) = 0.60219` — which *looks* like a
rescue. It is not. The ramp is auto-ranged to each region's own p1–p99, so a field with no
dynamic range has its **noise** stretched to full scale. `far`'s window is `(1.3e-9, 1.1e-8)`.
`spread_shape` is a dimensionless chord distance, so a p99 of `1e-8` means the copies agree to
eight digits.

**A ratio threshold missed it** — span `x8`, against `x1155` for near-field and `x220041` for
`deep interior`. The warning now carries a second arm comparing the window against the region's
**own median energy drift**, a measured floor rather than a chosen constant.

### 11.4 The lightness field carries a scale term

An adaptive render draws several levels in one image, and `ensemble_spread` is a spread over
copies jittered within the **cell**.

| level | cell width | median `spread_shape` | ratio | median `t_end` | ratio |
|---|---|---|---|---|---|
| 0 | 1.4286e-2 | 3.6675e-3 | — | 13.0 | — |
| 1 | 7.1429e-3 | 2.4591e-3 | 1.491 | 13.0 | 1.000 |
| 2 | 3.5714e-3 | 1.9301e-3 | 1.274 | 13.0 | 1.000 |
| 3 | 1.7857e-3 | 1.5990e-3 | 1.207 | 13.0 | 1.000 |
| 4 | 8.9286e-4 | 1.3423e-3 | 1.191 | 13.0 | 1.000 |

Cell width halves each level, so a proportional field would read `2.000` in every row. Measured
`1.19–1.62`, **falling with depth** — sub-linear and saturating, the chaotic-divergence
signature, and consistent with the standing measurement that a 9x fall in cell width moves
`sigma_E(0)` by 8.6x and `ensemble_spread` by only 2.1x. Under the log ramp it is about 12% of
the lightness range across five levels. `t_end` is the scale-free control and is flat.

Stated rather than corrected — normalising by cell width would change what the field means — but
`error(B)` under a spread colouring does include a small term for *this leaf is coarse*.

### 11.5 Two corrections to the chart reference

**`sum p_i = 0` does not catch the crossed-mass swap**, which is the one thing the reference says
it is for. Both forms give `p_lam*(1 - (m0+m1)/M01) = 0`.

| | crossed (correct) | uncrossed (the swap) |
|---|---|---|
| `\|sum p\|` | 7.9e-17 | 5.6e-17 |
| Jacobi round-trip error | 1.1e-16 | **6.8e-2** |
| kinetic-energy identity error | 4.4e-16 | **2.6e-1** |

Both catches are **empty at `m0 == m1`**, where the two forms are the same expression. Burrau has
`m0 = 3, m1 = 4`; the mass simplex passes through that line.

**`(Lz,E)` and `(Lz,K)` are one chart with two labels.** The reference lists them separately as
its most-machinery item, but its own warp parameterises both by `K(t) = K_max t^gamma` and then
reports `E = U + K(t)` — a relabelling, not a different sweep. Bitwise identical over the unit
square; only `gamma_k` makes them differ.

### 11.6 The chart families are tame where they are centred

Thirteen instances, budget 40000 so the descent stops on the criterion rather than the cap:

| case | quads | leaves | `alpha` med | `alpha` idec | ramp span |
|---|---|---|---|---|---|
| `body_plane` | 549 | 412 | 0.1402 | 1.1203 | 196 |
| `plane_00deg` *(control, bitwise)* | 549 | 412 | 0.1402 | 1.1203 | 196 |
| `shape_sphere` | 1293 | 970 | 0.1905 | 0.9709 | 2108 |
| `latent_shape` | 5441 | 4081 | 0.9995 | 0.0355 | 3.7 |
| `latent_inner_p` | 5461 | 4096 | 1.0082 | 0.1912 | 9.4 |
| `latent_outer_p` | 5101 | 3826 | 1.0045 | 0.1033 | 20.9 |
| `latent_mass` | 5461 | 4096 | 0.9956 | 0.0903 | 6.8 |
| `latent_mixed` | 4553 | 3415 | 0.9997 | 0.0486 | 11.3 |
| `latent_oblique_a` | 5293 | 3970 | 1.0005 | 0.0434 | 4.4 |
| `latent_oblique_b` | 4845 | 3634 | 1.0021 | 0.0791 | 3.6 |
| `burrau_nu_k` | 4753 | 3565 | 1.0038 | 0.0788 | 6.0 |
| `invariant_lz_k` | 5461 | 4096 | 1.0000 | 0.0817 | 7.3 |
| `mass_simplex` | 5461 | 4096 | 0.9996 | 0.1287 | 4.9 |

`alpha` sits at **0.99–1.01** on every new chart against **0.14** for `body_plane`. `alpha` near
1 means splitting halves the spread — refinement pays — so the scheduler refines everywhere.
That is correct behaviour on a tame region, not a scheduler fault.

**But it means these charts are not exercising the criterion where it is hard.** Tameness is a
property of *where* a chart is centred, and the base latent point was chosen to avoid the `z = 0`
symmetry rather than to find chaos. Moving it until the picture gets interesting would be tuning.
The honest report: the question these charts were added to answer — does the criterion behave
consistently across charts — is not answered by a set of charts that are all tame.

### 11.7 The criterion, re-measured under the colouring that ships

Complete tree to level 7, `N = 8`, `E+1 = 8`, 1024², `t = 13`, `f64`. 11,184,640 trajectories per
region. `error(B)` is the mean per-pixel OKLab distance to the fully-refined tree; `error = 0`
means **matches this sampling**, not correct.

At `B = 383`:

| ranking | near-field | deep interior |
|---|---|---|
| `greedy_oracle` | **0.12243** | **0.07387** |
| `frac_hot_between/median` | **0.12486** | **0.08075** |
| random band | 0.13059 – 0.14726 | 0.08334 – 0.10058 |
| `contrast:within/median` | 0.13668 | — |
| `layout/median` | 0.13968 | 0.10487 |
| `within/mean` | 0.15245 | 0.11460 |
| `between/median` | 0.15243 | 0.10205 |
| **`within/median` (shipped default)** | **0.15764** | **0.09977** |
| `running_max/median` | 0.15787 | 0.09095 |

**`frac_hot_between/median` is the only criterion that beats the random band, and it does so in
both measurable regions.** The shipped default is still beaten by random.

**This shifts §10.1.** Swapping the *within* arm for the *between* arm is not the fix —
`between/median` is beaten by random too. What fixes it is the **aggregation**: counting how many
footprints are hot rather than taking their median. That is what the build predicted before
running, and it now has a number under the colouring that ships.

### 11.8 Greedy is beaten by a scan order where the field is featureless

In `far` at `B = 12287`, `greedy_oracle` reads **0.26391** while every other criterion reads
**0.10885** — beaten **2.4x**. And every non-greedy row is identical to five digits, *including*
`frac_hot_within`, `frac_hot_between` and `first_div`, which each have **one distinct value** —
no ordering at all. So the ranking is irrelevant there, and greedy's loss is not to a better
ranking but to any ranking: immediate `delta-error` is noise in a featureless region and greedy
chases it.

This is the standing rule doing its work. A test asserting `greedy_oracle` dominates would have
fired here on correct behaviour.

Read §11.3 before quoting any `far` number: its ramp sits at the integrator's arithmetic.

### 11.9 What broke the images, and it was not the physics

`budget_far_t13_wire_01.png` was **128x64, 5 KB**. Wireframe lines are written with integer pixel
`set` calls and adaptive texels are nearest-neighbour, so neither can be soft in the file — a
blurry image is always a viewer upscaling a small raster.

The cause was a `criterion_metric -- 3 8` **validation run** allowed to write into `results/`,
overwriting committed 512² artefacts with 64²-derived ones. Everything is re-rendered at 1024².
`pan_sequence` also had `frame_res` hardcoded at 384 while its `viewport` argument set only the
camera.

`.fcache` files are gitignored: at level 7 they are 1.4M footprints x 14 `f64` each, 940 MB for
six, with no redundancy to remove at one sample per pixel. They make a colouring change a replay
**within a session**; committing a gigabyte so that survives a clone is the wrong trade. The
regeneration command is in `.gitignore` beside them.

## 12. The GLSL latent chart, and four slices you can recognise

Everything in §1–§11 was evaluated against nothing. The criterion was compared to random, the
colouring to its own percentiles, the tree to its own wireframe — but no render was of a
configuration anybody could name, so there was no independent handle on whether the physics or the
colouring was right. `Ma1achy/principia-ii` has one: a validated GLSL implementation with a pinned
decode and four named default slices.

`src/shaders/principia/frag.glsl:19-59` and `src/state.ts:71-76`. Nothing else from that repo was
needed.

### 12.1 The port was three constants

The reconstruction algebra in `src/physics/decoder.rs` already transcribed the GLSL correctly,
crossed mass factors included, and `Latent` was already 8 coordinates in the spec's order. Checked
line by line:

| GLSL | in tree | |
|---|---|---|
| `r01 = -m.z*lambda; r2 = M01*lambda` | `r01 = lam*(-m2/Mtot); r2 = r01 + lam` | same: `Mtot = 1`, and `r01+lam = (1-m2)lam` |
| `muLambda = m.z*M01` | `mu_lam = m2*M01/Mtot` | same |
| `r0 -= (m.y/M01)*rho; r1 += (m.x/M01)*rho` | identical | the crossing, verbatim |
| `p0 = -pRho - (m.x/M01)*pLambda` | identical | |
| `beta = PI*sigmoid(z)` | identical | half range — no wrap in this chart |

What was wrong was the saturation, in three places the LaTeX chart reference had guessed at:

| | was | GLSL |
|---|---|---|
| `MU_MAX` | `4.0` | `5.0` |
| `Q_MAX` | `1.0` | `2.0` |
| mass saturation | `MU_MAX*tanh(z)` | `MU_MAX*(2*sigmoid(z) - 1)` = `MU_MAX*tanh(z/2)` — **half the gain** |

Measured separation between the two saturation forms at `z_mu = (1.0, -0.5)`: worst `|dm| = 0.0899`
on a mass of ~0.2. Between `Q_MAX = 2` and `Q_MAX = 1` at `z_q = (0.7,-1.1,0.3,0.9)`:
`|dp_lambda| = 0.4474`. Both are far outside rounding, which is what makes the pinning test able
to fire.

The GLSL's `decodeIC` takes `z0..z9` and **never reads `z2` or `z3`** — dead slots from before the
chart was known to be 8D. Dropped. It also puts `beta` at index 0 where the spec names the chart
`(z_alpha, z_beta)`; this port follows the spec, so **the preset images are transposed relative to
the GLSL**. That is faithfulness, not a bug, and it is written into `decoder.rs`'s module header
next to the index table so the next reader is not caught by the same thing twice.

### 12.2 The landmark, and what it cannot see

At `z = 0` the decode gives the **equilateral Lagrange configuration** — masses `(1/3,1/3,1/3)`,
positions `[(-0.866025,-0.5), (0.866025,-0.5), (0,1)]`, `I = 1`, released from rest. Measured
separations:

```
1.732050807568878  1.732050807568877  1.732050807568877     (sqrt 3 = 1.732050807568877)
```

This is the strongest check in the build, because it is a named physical configuration and it can
be verified **by eye in the render** rather than only in a test. It is also, and this is the part
worth writing down, **blind to every constant §12.1 corrected**: at the origin the momentum
coordinates and the mass logits are all zero, so `MU_MAX`, `Q_MAX` and the choice of saturation
form each drop out of the arithmetic entirely.

`I = 1` and `COM = 0` are weaker still — algebraic identities of the canonical-frame decode
(`I = cos^2 a + sin^2 a`; `m0r0 + m1r1 = -M01 m2 lam` cancels `m2 r2`) that hold under **any** mass
factors whatever. All three are kept as wiring guards and labelled as such; the constants have
their own test.

Two further things the port does not get for free. The brief names `sum p = 0` as the test that
catches a crossed-mass swap; it cannot, for reasons already measured and recorded here
(`7.9e-17` crossed against `5.6e-17` uncrossed) — the position-side factors never enter the sum.
And the gauge: `decode` applies `canonicalise` and `scale_gauge` where `decodeIC` applies neither,
which should be inert on this chart since `rho~` already sits on `+x` and `lam~_y >= 0` for
`beta in [0,pi]`. It is, but not bitwise — `I = 1` only in exact algebra, so `scale_gauge` divides
by `sqrt(1 +- eps)`. Measured over four latent points: positions `7.161e-15`, momenta `3.140e-16`.

### 12.3 The four presets

`z0 = 0`, framed at `(0,0)` with `half = 3.0`. The reference's `z0 + (2u-1)q1 + (2v-1)q2` over
`[0,1]^2` is reproduced exactly by that framing, because `decode_state` is `z0 + u*q1 + v*q2` and
the slice already supplies a signed box — one fewer place for a factor of two to hide.

**The window shipped at `1.0` and was wrong.** The reference UI reads `Slice +/- 3.0e+0`:

```
half = 1.0  ->  alpha in [0.446, 1.125],  beta in [0.845, 2.297]  =  46% of the azimuth
half = 3.0  ->  alpha in [0.120, 1.451],  beta in [0.149, 2.993]  =  90% of the azimuth
```

So every first-cut preset image was a 3x zoom on the middle of the picture. In the GLSL the
fractal core is a small disk surrounded by large smooth regions; at `half = 1.0` it fills the
frame. Same structure, wrong crop — which is exactly why the port read as *similar but not the
same*. The number now comes from `Chart::default_half()`, which is chart-aware: a `BodyPlane`
coordinate is a body position in Burrau units and a `Latent` coordinate is a sigmoid pre-image, and
**one shared default silently meant two different things**, which is how this got through. The
`_h1` rows of §12.4 are the crop control — same chart, same basis, one number changed.

| case | GLSL id | `q1` | `q2` |
|---|---|---|---|
| `preset_shape` | `shape` | `e_alpha` | `e_beta` |
| `preset_prho` | `prho` | `e2` | `e3` |
| `preset_plambda` | `plambda` | `e4` | `e5` |
| `preset_shape_pl` | `shape_pl` | `e_alpha + e_pLambda_y` | `e_beta + e_pLambda_x` |

**`shape_pl` pairs by GLSL SLOT, and the first cut got it wrong.** The reference is
`q1 = e0 + e6`, `q2 = e1 + e7`; in *its* indexing (`z0 = beta`, `z1 = alpha`,
`z6/z7 = pLambda.x/y`) that pairs **beta with `pLambda.x`** and **alpha with `pLambda.y`**. This
port renumbers alpha and beta into the spec's order and must carry their momentum partners with
them. It did not, and paired alpha with `pLambda.x` — so the *pair assignment* transposes and each
pair stays intact.

**No transposition repairs it.** Swapping `q1` and `q2` gives `e_beta + e_pLy`,
`e_alpha + e_pLx` — still crossed. It is a genuinely different 2-plane through the 8D space, not a
reorientation of the same one, which is why it rendered as **twisted** rather than tilted: the
coupling sets how momentum co-varies with configuration across the slice, and the two pairings
give different shears. Measured at `(u,v) = (0.7,-1.3)`, `max |dIC|` against the correct plane:

```
shape_pl vs             crossed at (0.7,-1.3): max |dIC| = 5.4483e0
shape_pl vs crossed, transposed at (0.7,-1.3): max |dIC| = 1.2042e0
```

Both far outside rounding, which is what makes
`tests/charts.rs::the_shape_pl_preset_pairs_alpha_with_p_lambda_y` able to fire. An index
assertion alone would have passed on the transposition, and the whole finding is that transposing
does not work.

`shape_pl` is the **only** preset with a cross-coupling — the other three are pure-config or
pure-momentum — so it is the only one that can fail this way. That consistency is itself evidence
the diagnosis is right.

`preset_shape_pl` is constructed directly and **not** through `latent_oblique`: the reference's
basis is un-normalised (each direction has norm `sqrt 2`) and Gram–Schmidt would quietly render a
different slice while looking like a tidy-up. `tests/charts.rs` pins the norm so that fails loudly.
All four bases are built by `Chart::preset_*` rather than written at each call site: the `shape_pl`
literal appeared in three places and was wrong in all three.

These sit **beside** the existing `latent_*` rows, not in place of them. Those are deliberately at
an off-origin `z0` so no sigmoid rests at its symmetry point; the presets are at the origin for the
opposite reason — that point is the Lagrange configuration.

**`prho` and `plambda` are constant-configuration slices**, and this is the finding that cost a
test. Positions in this decode do not depend on the momentum coordinates at all, so every pixel of
those two slices is the **same triangle** released with a different initial velocity. A
distinctness check keyed on a body position reads `1 of 81` there:

```
     shape: 81 distinct ICs of 81, over 81 distinct configurations
      prho: 81 distinct ICs of 81, over  1 distinct configurations
   plambda: 81 distinct ICs of 81, over  1 distinct configurations
  shape_pl: 81 distinct ICs of 81, over 81 distinct configurations
```

Nothing has collapsed — the chart simply does not vary the quantity being counted. But a collapsed
decode gives `ensemble_spread` exactly zero, which reads as *perfectly resolved* and stops the
descent, so the two are worth being able to tell apart: the guard has to measure what the chart
actually moves.

**The consequence is that those two are a control.** Every pixel of `prho` and `plambda` starts at
the *same triangle*, so `spread_shape` at `t = 0` is identically zero across the whole slice. Any
structure in them is purely momentum-driven — which makes them the pair that separates
configuration effects from momentum effects, and it is why `preset_shape` (a configuration sweep)
and `preset_shape_pl` (a mixed basis) can be read against them at all.

### 12.4 The gallery, regenerated — and three corrections, one of them to this section

All 26 cases at 1024², budget 40000, `tau = 1e-4`, `alpha_hi = 0.2`, `N = 8`, `E+1 = 8`, `t = 13`,
f64, screen floor on. Thirteen pre-existing instances, the four presets at the corrected window,
their four `_h1` crop controls, and five `latent_*_h3` extent controls.

The control holds: `plane_00deg` against `body_plane`, `max |dIC| = 0e0`, asserted inside the run.
The thirteen pre-existing rows reproduce their committed values exactly; their `.prnq` dumps differ
from the committed ones in `wall_seconds` and in no other field.

| case | half | quads | leaves | alpha med | alpha idec | ramp span |
|---|---|---|---|---|---|---|
| `body_plane` | 0.05 | 549 | 412 | 0.1402 | 1.1203 | 196.2 |
| `plane_00deg` | 0.05 | 549 | 412 | 0.1402 | 1.1203 | 196.2 |
| `shape_sphere` | 0.05 | 1293 | 970 | 0.1905 | 0.9709 | 2107.6 |
| `latent_shape` | 1.5 | 5461 | 4096 | 1.0013 | 0.0383 | 5.5 |
| `latent_inner_p` | 1.5 | 5397 | 4048 | 1.0091 | 0.1170 | 18.5 |
| `latent_outer_p` | 1.5 | 5041 | 3781 | 1.0046 | 0.3119 | 28.9 |
| `latent_mass` | 1.5 | 5461 | 4096 | 0.9979 | 0.0413 | 3.8 |
| `latent_mixed` | 1.5 | 5117 | 3838 | 1.0001 | 0.0696 | 21.1 |
| `latent_oblique_a` | 1.5 | 5421 | 4066 | 1.0026 | 0.0745 | 5.8 |
| `latent_oblique_b` | 1.5 | 4697 | 3523 | 1.0024 | 0.2446 | 4.0 |
| `burrau_nu_k` | 0.45 | 4753 | 3565 | 1.0038 | 0.0788 | 6.0 |
| `invariant_lz_k` | 0.45 | 5429 | 4072 | 1.0001 | 0.0835 | 9.2 |
| `mass_simplex` | 0.45 | 5461 | 4096 | 1.0010 | 0.1319 | 5.2 |
| **`preset_shape`** | **3.0** | **21** | **16** | **1.2685** | **13.4654** | **18966.6** |
| `preset_prho` | 3.0 | 3593 | 2695 | 1.0001 | 0.7106 | 33.2 |
| `preset_plambda` | 3.0 | 3361 | 2521 | 0.9929 | 0.9491 | 44.7 |
| `preset_shape_pl` | 3.0 | 1729 | 1297 | 1.0311 | 0.5315 | 89.0 |
| `preset_shape_h1` | 1.0 | 769 | 577 | 0.6285 | 3.0724 | 10008.9 |
| `preset_prho_h1` | 1.0 | 4357 | 3268 | 1.0153 | 0.1198 | 19.7 |
| `preset_plambda_h1` | 1.0 | 4493 | 3370 | 0.9983 | 0.2893 | 14.7 |
| `preset_shape_pl_h1` | 1.0 | 4645 | 3484 | 0.9882 | 0.3613 | 13.3 |
| `latent_shape_h3` | 3.0 | 5269 | 3952 | 0.9989 | 0.5933 | 38.5 |
| `latent_inner_p_h3` | 3.0 | 3821 | 2866 | 1.0067 | 0.1618 | 77.2 |
| `latent_outer_p_h3` | 3.0 | 3877 | 2908 | 1.0023 | 0.6853 | 71.8 |
| `latent_mass_h3` | 3.0 | 4953 | 3715 | 0.9966 | 0.1216 | 18.7 |
| `latent_mixed_h3` | 3.0 | 4953 | 3715 | 0.9994 | 0.2110 | 66.8 |

#### 12.4a The previous section's headline was a fact about the crop

**The committed `preset_shape` row is reproduced exactly by `preset_shape_h1`** — 769 quads, 577
leaves, `alpha` med 0.6285, interdecile 3.0724, ramp span 10008.9, every figure identical. Same for
`preset_prho_h1` (4357/3268, 1.0153, 0.1198, 19.7) and `preset_plambda_h1` (4493/3370, 0.9983,
0.2893, 14.7). Three of the four committed preset rows were measured at `half = 1.0` and the
controls reproduce them to the digit.

So the claim **"`preset_shape` is the only chart family instance that is not tame"** does not
survive its own window. What replaces it is worse and more useful — see 12.4c.

#### 12.4b The `shape_pl` basis and the crop, separated

`preset_shape_pl` changed twice, so the `_h1` twin is what tells the two apart. It carries the
**corrected** basis at the **old** window, which isolates each cause to one comparison:

| | half | basis | leaves | alpha idec |
|---|---|---|---|---|
| committed | 1.0 | crossed | 2974 | 0.2344 |
| `preset_shape_pl_h1` | 1.0 | correct | 3484 | 0.3613 |
| `preset_shape_pl` | 3.0 | correct | 1297 | 0.5315 |

The basis alone, at a fixed window: **+17.2%** on leaves (2974 -> 3484). The window alone, at a
fixed basis: **2.69x** the other way (3484 -> 1297). Both real, and neither would have been
separable from the other without the control.

#### 12.4c **The trees are set by a camera veto, not by the criterion** — and the table said `crit`

This is the largest fact about the gallery and it was invisible in every previous version of this
table. `chart_gallery` printed a `bound` column reading `crit` unless the *budget* was exhausted.
It never asked what actually stopped the descent. The `.prnq` dumps carry `decision` per quad and
always did; read back, they say:

| case | leaves | veto % | floor % | keep % |
|---|---|---|---|---|
| `latent_shape`, `latent_mass`, `mass_simplex` | 4096 | **100.0** | 0.0 | 0.0 |
| `latent_oblique_a` | 4066 | 100.0 | 0.0 | 0.0 |
| `invariant_lz_k` | 4072 | 99.8 | 0.2 | 0.0 |
| `preset_shape_pl_h1` | 3484 | 99.4 | 0.5 | 0.1 |
| `preset_prho` | 2695 | 98.1 | 1.9 | 0.0 |
| `preset_shape_pl` | 1297 | 96.8 | 2.5 | 0.6 |
| `preset_plambda` | 2521 | 95.8 | 4.2 | 0.0 |
| `shape_sphere` | 970 | 79.2 | 18.1 | 2.7 |
| `preset_shape_h1` | 577 | 77.6 | 12.8 | 9.5 |
| `body_plane`, `plane_00deg` | 412 | 61.2 | 37.4 | 1.5 |
| **`preset_shape`** | **16** | **0.0** | **50.0** | **50.0** |

The veto is `Decision::MaxRelDepth` — `Camera::veto`, a cap at `camera_depth + max_rel_depth`
levels. **On 23 of 26 charts it stops 95% or more of the leaves.** On three of them it stops
*every* leaf, and every leaf sits at one depth: those are complete capped trees, and their leaf
counts are facts about the cap in exactly the way a budget-bound row's is a fact about the budget.

This reframes §11.6's standing result. *"The chart families sit at `alpha` 0.99–1.01 and do not
exercise the criterion"* is true, and now has a mechanism it never stated: **the criterion decides
under 1% of leaves on those rows.** The `alpha` near 1.0 describes quads a cap forced, not quads
the criterion chose. Same lesson as §8's screen floor, at a second stop condition, and it went
unnoticed for the same reason — nothing printed which one fired.

The `bound` column now reports `VETO n%` or `crit n%`. The full breakdown is
`results/output/gallery_table.txt`, re-derived from the dumps by `examples/gallery_table.rs`
rather than by a three-hour re-run for a label.

#### 12.4d At the corrected window, `preset_shape` is where the criterion fails outright

`preset_shape` is the **only** case in the set whose tree is entirely its own: 0% veto, 8 leaves
stopped by `Floor` (spread below `tau`) and 8 by `Keep` (`alpha` says refinement does not pay).
Sixteen leaves at depth 2, against a complete 4096.

It is not stopping on a featureless field. Its ramp span is **18966.6** — over four decades of
lightness range, the widest in the set — and its `alpha` interdecile is **13.47** against 0.04–0.99
for the tame rows. Widening the window brings in the large smooth surroundings that the GLSL shows
around the fractal core; their spread falls below `tau`, `Agg::Median` over each quad's `N x N`
footprints reads the quad as resolved, and the core — now a small disk inside a level-2 quad — gets
no refinement at all.

That is §5's *"median under-refines thin structure — blind to a filament crossing a quad"* at full
strength, on a chart with a named physical configuration at its centre. **The wrong window was
hiding it behind a plausible-looking 577-leaf tree.**

#### 12.4e The tameness result survives the extent axis

§11.6 and §12.3 held that a chart's tameness is set by which coordinates it varies, not by where it
is centred. The `latent_*_h3` twins give it an axis it had never been tested against, since every
`latent_*` row was measured at `half = 1.5`. It holds: all five stay at `alpha` med 0.9966–1.0067
across a 2x change in extent, and their veto fractions stay at 99.0–99.9%. Extent moves leaf counts
(`latent_inner_p` 4048 -> 2866) and it does not move tameness.

The one chart where extent matters enormously is the configuration sweep, and it matters by
breaking the criterion rather than by changing a distribution:

| | wide (3.0) | narrow (1.0) | ratio |
|---|---|---|---|
| `preset_shape` | 16 | 577 | **36x** |
| `preset_prho` | 2695 | 3268 | 1.21x |
| `preset_plambda` | 2521 | 3370 | 1.34x |

#### 12.4f The refinement-mechanism test is readable on 2 charts of 26

The proposal was that refinement chases non-convergence rather than structure: terminated regions
are absorbing, so nearby copies share an outcome, `spread_event` collapses to zero and the quad
reads *resolved*, while still-running regions keep diverging and hold high spread forever. That
predicts leaf depth **anti-correlated** with `terminated_fraction`, and both are already in the
`PRNQ` dump — one plot per chart and no extra integration.

**Before reading any correlation, count the distinct values.** There are three ways for this pair
to be uninformative and the Spearman alone cannot tell them apart:

- **x constant** — every leaf at one depth, so there is no depth axis. `preset_shape` (16 leaves,
  all at level 2), and the three complete capped trees.
- **y constant** — `terminated_fraction` takes exactly one value. Ten charts, all at 1.000.
- **y saturated** — one value holds over 90% of leaves, so the correlation is read off the thin
  remainder. Twelve charts.

That leaves **two readable charts**, and they disagree:

| | tf values | modal % | mean escape | spearman | per-depth medians |
|---|---|---|---|---|---|
| `shape_sphere` | 52 | 75.5 | 0.0993 | **-0.2245** | 0.766, 0.484, 0.078, 0.000, 0.000 |
| `preset_shape_h1` | 18 | 76.9 | 0.7053 | **+0.3756** | 0.812, 0.875, 0.984, 1.000, 1.000 |

`shape_sphere` matches the prediction and matches it more strongly than its Spearman suggests —
the medians fall by four steps while the rank correlation is diluted, because 908 of its 970 leaves
sit at the two deepest levels and their interdecile spans the full `[0,1]`. **Read the per-depth
medians, not the pooled correlation**; the depth distribution is far too unbalanced for one number.

`preset_shape_h1` runs the other way. Its `terminated_fraction` is climbing to 1.000 with depth and
its deepest levels sit at `1.000 [0.984, 1.000]`, which is close enough to saturation that it is
listed as readable on a 76.9% modal share and should be read with that in mind.

**So the mechanism is neither established nor refuted.** One chart for, one against, and 24 on
which the test cannot fire. What the run does establish is *why* it cannot fire: on 22 of 26 charts
`escape_fraction` is 0.9894–1.0000 — at `t = 13` on these charts essentially everything has escaped,
so `terminated_fraction` has no range to correlate against.

That is also a scope note on a standing result. **"The escape arm contributes nothing at `t = 13`"
is about Burrau's near-field body plane and does not generalise.** On the latent charts the escape
fraction is ~1.0. `preset_shape` is the counter-example within the latent family: escape 0.0547
with `terminated_fraction` median 0.984, which is to say its terminations are **collisions** — its
event histogram is dominated by `collision d0=361886, d1=234807, d2=234812`. Carrying
`terminated_fraction` and `escape_fraction` separately is what makes that legible.

## 14. The threshold, and the mask it saturated

The refinement criterion was diagnosed as *"an always-split rule wearing a threshold's clothing"*.
Re-derived here over **all 69 committed `.prnq` dumps, 92,880 leaves**, by
`examples/threshold_diagnosis.rs` — no integration, the dumps carry every quad's spread and
decision already.

### 14.1 `tau` sits in the bottom few percent of the distribution it is meant to cut

**Both scopes, because they differ and the first draft of this section mixed them.** It quoted the
`charts/`-only figures under a heading that said "all 69 dumps" — the mislabelled-denominator
fault this same section warns about, committed inside it.

| | `charts/` only (26 dumps, 75,359 leaves) | whole corpus (69 dumps, 92,880 leaves) |
|---|---|---|
| p1 / median / p99 | `1.47e-4` / `6.61e-4` / `7.87e-3` | `2.69e-5` / `6.73e-4` / `9.82e-3` |
| **`tau = 1e-4`** | **99.6% exceed, 0.4th pct** | **95.7% exceed, 4.3rd pct** |
| `3e-4` | 88.4%, 11.6 | 84.8%, 15.2 |
| `1e-3` | 29.1%, 70.9 | 30.2%, 69.8 |
| `3e-3` | 3.1%, 96.9 | 3.3%, 96.7 |
| `1e-2` | 0.9%, 99.1 | 1.0%, 99.0 |

The two agree from `1e-3` up and part company below it: the whole corpus carries the `vertical/`
zoom ladder and `far`, whose spreads reach `4.26e-8`, so its lower tail is three orders deeper.
**The conclusion is the same on either scope** — the predicate is true for 96–99.6% of quads and
`tau` sits in the bottom few percent of its own distribution. Which scope a figure comes from is
stated wherever one is quoted now.

The sweep ladder in `sched_sweep` and `sweep_screen` used to run `1e-8 … 1e-2`. It now runs
`1e-8 … 1e-1`: the top extended past the point where the predicate goes false everywhere, and
**both** low rungs kept as labelled degenerate controls — for *different regions*. Dropping one of
them was an error; see §14.10.

### 14.2 The failure has two sides, and that is the argument for rank

Both are in the corpus.

- **`tau` below the bulk** — everything splits, tree uniform at **max depth**. Most charts.
- **`tau` above the bulk** — everything keeps, tree uniform at **depth 2**: 16 leaves against a
  complete 4096. **16 of the 18 trees the camera veto does not bind** are stopped this way, with
  leaf-spread medians running `9.45e-5` down to `4.26e-8` against `tau = 1e-4`, and leaf decisions
  that are `keep` almost to the last one — `far`, `deep interior`, every deep zoom step.

Selectivity requires the threshold to **cut through the bulk**, and the chart's own dynamic range
decides whether any fixed value can. Measured across the 69 dumps, `spearman(p99/p1, depth
variance) = +0.727` and `spearman(p99/p1, % at max depth) = −0.702`. So a fixed threshold has a
narrow window of usefulness that varies per chart. **A ranking cannot land above or below a
distribution; it always cuts through it.** That is a stronger argument for rank than the treadmill
one, and it did not depend on predicting anything.

### 14.3 `preset_shape` is a third mode, and not the one it looks like

It was tempting to file `preset_shape` under the upper-side failure — 16 leaves, depth variance 0,
the widest dynamic range among the charts. **The decision column says otherwise.** Its leaf-spread
median is `2.86e-1`, **3400× above `tau`**: it clears the spread gate on every leaf and is stopped
by `alpha` — 8 `floor` (below `alpha_lo`) and 8 `keep` (between the thresholds). Not one leaf
failed the spread test.

So `preset_shape` is the only tree in the corpus where the **alpha gate** is what is exercised,
which makes it the cleanest instance of the standing result that `alpha_hi` does more work than
the criterion. Quoting it as a `tau` failure would have been a mechanism read off a shape.

### 14.4 The stop-reason column is the headline

`Decision::ScreenFloor` or `MaxRelDepth` — a **camera veto** — stops ≥95% of leaves on **21 of the
69 dumps**, and 100% on several. On those the criterion decides almost nothing.

That reframes *"13 of 17 charts at 99.4–100% max depth"*: it is not the criterion saying **split**
and being obeyed. It is the criterion never saying **stop**, with something else terminating the
descent. **The observed uniformity is what a permissive criterion looks like when a veto is doing
the stopping.** The eighteen trees where the veto binds on under 5% of leaves are the only ones
that are entirely their own decisions, and `preset_shape` is the only chart among them.

Never quote a leaf count without its stop-reason breakdown.

### 14.5 The hot mask: relative, and both rules kept

The instruction on record was to make the hot threshold relative — the quad's own median. Two
things that instruction missed, both now measured.

**`n_hot` stops being a signal under any quantile rule.** On a field with distinct values the count
above the cut is set by the rule, not the field — 31 of 64 at `N = 8, q = 0.5`, given nearest-rank
and a strict comparison. So `frac_hot` carries essentially no information once the mask is
relative. That matters because **`frac_hot_between/median` is the best criterion measured on this
project** (§11.7 — the only one beating the random band in both measurable regions). A relative
rule that *replaced* the absolute mask would have deleted the best-performing signal in the system
and read as an improvement.

**So both masks are computed and both are dumped.** `frac_above_tau_*` and the `frac_hot_*`
criteria are untouched; `spatial::HotRule` selects which mask the *shape* criteria read.

The one exception, measured rather than assumed away: on a **tied** field the count is set by the
tie structure. A two-valued field reads the same count at `q = 0.5, 0.75, 0.9` alike — the case
that occurs when the event arm, with its five distinct values, dominates a footprint field.

### 14.6 The desaturation, and what it costs

Paired against the absolute mask on **one descent** (`tests/criterion.rs`), because the form the
brief asked for — *"assert `n_hot < N²` for a stated majority"* — passes trivially and
unconditionally under a quantile rule and is decoration:

| | absolute, `tau = 1e-4` | relative, `q = 0.5` |
|---|---|---|
| `n_hot == N²` | **100.0%** | **0.0%** |
| `n_components ≤ 1` | **100.0%** | 59.7% |
| `n_components > 1` | 0.0% | **40.3%** |

The absolute arm is the control. Without it, a relative mask that happened to saturate too would
read the same as a working one.

Across regions (`examples/hot_rule_sweep.rs`, budget 600, `criterion=within` so the tree is
constant down each block and asserted so):

| region | rule | sat% | 1-comp% | median components | `d(Layout)` | `d(LayoutRel)` |
|---|---|---|---|---|---|---|
| `far` | abs 1e-4 | 0.0 | 100.0 | **0** | 1 | 1 |
| `far` | q 0.50 | 0.0 | 100.0 | 1 | 1 | 1 |
| near-field | abs 1e-4 | 48.0 | 82.1 | 1 | 58 | 78 |
| near-field | q 0.50 | 0.0 | **1.8** | **5** | 58 | 26 |
| near-field | q 0.90 | 0.0 | 2.7 | 4 | 58 | 9 |
| `deep interior` | abs 1e-4 | 0.0 | 50.0 | 2 | 17 | 12 |
| `deep interior` | q 0.50 | 0.0 | **0.0** | **7** | 17 | 16 |

**Two things here were not expected.**

**The saturation is not uniform across regions, and in `far` the absolute mask is EMPTY, not
full.** `far`'s leaf-spread median is `4.26e-8` against `tau = 1e-4`, so nothing clears the cut:
`n_hot == 0`, `perimeter_ratio` is `NaN` by the empty-set convention, and every criterion built on
it takes **one** distinct value over all 16 leaves. `deep interior` already resolves a median of 2
components under the absolute rule. It is near-field and the latent charts where the mask is full.
"Saturated everywhere" is the pooled number, not the regional one — and a full mask and an empty
mask are the same threshold landing on either side of the distribution, which is §14.2 again one
level down.

**The relative rule desaturates the mask and coarsens the ordering.** Near-field's median component
count runs 1 → 5 from absolute to `q[0.50]` — the mask finally describing something — while
`Criterion::LayoutRel`'s distinct-value count falls **78 → 26 → 17 → 9** across `abs, q[0.50],
q[0.75], q[0.90]`, against `Criterion::Layout` holding at 58. With `n_hot` pinned by the rule,
`largest/n_hot` can only take as many values as there are component sizes.

Reported, not hidden, and it settles nothing by itself: the standing result is that **signal
resolution is not what makes a ranking good** (§10.4 — `frac_hot_between` is the best criterion
measured here, on 65 distinct values, beating a 4994-valued one). `error(B)` decides. But a
criterion whose ordering coarsens as its input improves is worth watching.

### 14.7 `grad_rms` — the control on the whole mask family

The one structure measure with **no threshold in it**: RMS of the forward-difference gradient
across the `N×N` footprint grid, `NaN` (never 0) when no adjacent pair is finite. If a masked
signal cannot beat it, the mask is not earning its parameter.

### 14.9 The sweep, re-run — and a fourth stop reason

`sched_sweep` (no camera) and `sweep_screen` (camera framing the root, 512²), ladder
`1e-8 … 1e-1` × `alpha_hi ∈ {0.2, 0.5, 0.8, 1.0}`, budget 2000.

**`alpha_hi` dominates `tau` outright without the veto.** Near-field goes **1498 → 19 leaves**
between `alpha_hi` 0.20 and 0.50 — a **79× collapse** — and at `alpha_hi ≥ 0.5` **`tau` changes
nothing at any rung of the ladder**. `tau` is live in exactly one row of thirty-six. Under the
veto that ratio falls to **21.68×**, reproducing the standing result that the screen floor demotes
`alpha_hi` and promotes `tau`.

**Within the live row, three rungs give a bitwise identical tree.** `1e-8`, `1e-6`, `1e-4` and
`3e-4` all read 1997/1498/499 in near-field. The whole live range of `tau` is `3e-4 … 3e-3` — one
decade, spanning the 11.6th to the 96.9th percentile of §14.1's table. Everything outside it is a
constant predicate.

**The `tau` span over the whole ladder, at `alpha_hi = 0.20`, under the veto:**

| region | span | range |
|---|---|---|
| `far` | **×64.00** | 16 … 1024 |
| near-field | ×27.62 | 16 … 442 |
| `deep interior` | ×9.62 | 16 … 154 |

`far`'s **64×** is the number already on record, recovered exactly. That matters because of how it
was nearly lost — see §14.10.

**The fourth stop reason: budget.** The large low-`tau` trees are not criterion-bound either.
Near-field at `alpha_hi = 0.20` carries **869** `BudgetExhausted` leaves of 1498; `deep interior`
at `tau = 1e-8` carries **1357 of 1498**. The tree that looks selective there is a budget artefact.

So the criterion proper decides a minority of leaves in **every regime measured**:

| regime | what stops the descent |
|---|---|
| the chart gallery | `ScreenFloor` / `MaxRelDepth` — a camera veto, ≥95% of leaves on 21 of 69 dumps |
| low `tau`, no veto | `BudgetExhausted` — up to 91% of leaves |
| high `tau` | the spread gate — 16 leaves, all `keep` |
| `preset_shape` | the `alpha` gate — the only tree in the corpus that exercises it |

That is the shape of the problem, and it is a stronger statement than "the criterion refines too
much". **The criterion is rarely what decides anything.**

### 14.10 The ladder change was wrong once, and the region is why

The first cut of this dropped `1e-8` and `1e-6` as *"measuring the same always-split regime
twice"* — true in near-field, where `1e-8`, `1e-6` and `1e-4` are bitwise identical. It is **false
in `far`**, whose leaf-spread median is `4.26e-8`: `1e-8` is the only rung on the ladder below its
bulk. Without it `far` read 16 leaves in all 32 cells and the sweep said *"`tau` is inert here"* —
a statement about the ladder, not about the region, and it would have silently contradicted the
standing 64× result.

**Which rung is degenerate is a fact about the REGION, not about the ladder.** The regional spread
medians span six orders — `4.26e-8` in `far`, `9.45e-5` in `deep interior`, `9.75e-4` in near-field
— which is exactly what `sched_sweep`'s own module header has said since it was written. Both low
rungs stay, labelled as degenerate controls for different regions.

The same run also caught a hardcoded summary line: `sweep_screen` printed the `tau` span as the
ratio between two *named adjacent rungs*, `1e-8` and `1e-6`. With `1e-8` removed it read `×0.00`,
and with it present it would still have been reporting whichever pair happened to straddle the
bulk in one region. It now takes max/min over the whole ladder at the `alpha_hi` where `tau` is
live. **An argument hardcoded past is worse than an argument missing** — the same defect recorded
for `pan_sequence`'s viewport, at a different site.

### 14.8 The dump moved to PRNQ v3

Ten columns appended — the two relative layouts and the two `grad_rms` values — and the header now
carries `hot_rule=`. The new columns go **after** the existing 48, so a positional reader that
stops at 48 still reads every v2 field correctly. Both readers in this project parse the `fields=`
line by name and are unaffected. `.qcache` moves to PRQC v2 with the matching `sig_layout_rel`,
`sig_grad_rms` and their contrasts.

## 15. Rank, the two modes, and a premise that is wrong in sign

### 15.1 The queue never used the criterion it was configured with

`order_queue` sorted on `red.spread(agg)` and never read `cfg.criterion`. So every
`--order spread` run in the corpus ordered by the within arm whatever its header said, and the
budget-truncation point was decided by a different quantity than the one named. Fixed.

**It reproduces the corpus exactly**, and the reason is worth stating rather than asserting: every
committed run has `criterion=within`, and `signal(Within, agg)` *is* `spread(agg)` by definition.
The fix only changes runs that set a criterion the old code was ignoring. `tests/criterion.rs`
asserts the equality over the aggregations and over degenerate inputs.

### 15.2 The two modes are one mechanism

`Mode::Uniform` turns the criterion **off** — not "sets a permissive threshold", off — and splits
to the veto. `Mode::Balanced` ranks the frontier and gives the top `k_frac` its budget. A quad
that falls down the ranking is simply not spent on: the demotion §3.1 asks for, with no merging
and no eviction.

`k_frac = 1.0` refines the whole eligible frontier and reproduces the unranked descent exactly.
Deferred quads are marked **`Keep`, not `BudgetExhausted`** — they were outranked, not refused for
want of budget, and conflating them would hide the ranking inside the stop-reason column that
exists to expose it. Measured at `N = 4`, budget 300: `k_frac` 0.25 / 0.50 / 0.75 / 1.00 gives
13 / 46 / 121 / 223 leaves, with zero budget-exhausted below 1.0.

### 15.3 Balanced mode passes, and the control fails as it must

`examples/balanced_march.rs`, playhead `t ∈ {4, 6, 8, 10, 13, 16, 20}`, `n_sync` scaled with
`t_max`, viewport 64² so the uniform arm is stopped by the **veto** rather than the budget — a
budget-bound control is no control.

| | near-field | `deep interior` |
|---|---|---|
| balanced, depth variance | 0.004 – 0.574, no trend to zero | 0.402 – 0.740 |
| **uniform, depth variance** | **0.0000 at every `t`** | **0.0000 at every `t`** |
| balanced, churn | 0.000 – 0.429 | 0.000 – 0.194 |

The control is pinned at exactly zero across the whole march, which is what makes the balanced
row mean anything.

**Churn is reported over the SHARED quads only**, and flagged when that set is small. A quad
present at one playhead and not the other has not "changed decision", and counting it would fold
the tree's size change into a statistic about its stability. Near-field at `t = 16` shares only
14 quads, so its 0.4286 is 6 of 14 — printed as thin rather than quoted as a rate.

### 15.4 The treadmill does not happen. The opposite does.

§3's argument for rank is that *"spread grows with `t` everywhere"*, so any fixed threshold must
eventually fire on every quad and balanced mode must degenerate to uniform depth.

**Measured on the uniform arm** — a fixed tree, so this is the field and not the tree — the median
leaf spread does not grow:

| `t` | 4 | 6 | 8 | 10 | 13 | 16 | 20 |
|---|---|---|---|---|---|---|---|
| near-field | 1.62e-3 | 6.02e-4 | 1.98e-3 | **6.56e-3** | 1.87e-3 | 9.93e-5 | **8.09e-5** |
| `deep interior` | 1.63e-3 | 1.68e-4 | 1.45e-4 | 1.45e-4 | 5.36e-5 | 5.36e-5 | **5.31e-5** |

Near-field **peaks at `t = 10` and then falls 81×**, ending *below* `tau = 1e-4`. `deep interior`
falls **31× monotonically** and is below `tau` from `t = 13` on.

The mechanism is already on record one level down: **terminal states are absorbing.** As
termination saturates, copies share an outcome, `spread_event` collapses and `spread_shape` over
terminated trajectories stops growing. So at large `t` a fixed `tau` fires **nowhere**, the spread
gate keeps everything, and the tree **shrinks** — near-field 256 → 40 leaves, with the `keep`
count going 0 → 20 and `ScreenFloor` 164 → 12. That is the *upper-side* failure of §14.2 arriving
on the time axis.

**This strengthens the case for rank rather than weakening it.** The treadmill argument was that
no fixed `tau` can survive a monotone rise; the measured behaviour is a rise *and then a
collapse*, which is worse for a fixed threshold — there is no value that is correct at both ends,
and no monotone schedule that would track it either. A ranking is invariant to the whole curve.

### 15.5 The structure term needed a third factor, and a test found it

`QuadReduction::structure` is connectedness × thinness × extent, on the relative mask. The first
two were designed; the third was not.

| mask | structure |
|---|---|
| fully hot | **0.0000** — maximally connected, zero perimeter; thinness kills it |
| one-cell filament | **1.0000** — the target |
| checkerboard | **0.0039** — maximally thin, maximally scattered; connectedness kills it |
| single isolated cell | **0.1250** |

The isolated cell scored **1.0** on the first two factors alone: it *is* the largest component, so
connectedness is trivially 1, and `perimeter_ratio == 4` so thinness saturates. Maximum structure,
for one cell. **Extent** — `largest_component / N` — is the graded form of what
`Layout::looks_like_boundary` already encoded as `largest_component >= N/2`. Each of the three
factors catches a case the other two score at maximum.

`structure` is `NaN` on an empty mask, not 0. `far`'s absolute mask is empty on every leaf
(§14.6), and a 0 there would read as "no structure found" rather than "not measured".

## 16. §2.2 answered: structure neither replaces nor multiplies

`examples/structure_metric.rs`, levels 6 (5461 quads, 512² reference), `N = 8`, `E+1 = 8`,
`tau = 1e-4`, `t = 13`, under the **shipping** colouring. Three targets: `near-field`,
`deep interior` — because a change that only improves near-field is tuning — and `preset_shape`,
the **only tree in the corpus whose leaves are entirely its own decisions** (0% camera veto).
`far` is deliberately absent: its reference window is `(1.3e-9, 1.1e-8)`, the integrator's
arithmetic rather than physics.

Oracle-to-random separation, read first: **0.00597 / 0.00675 / 0.01768**. The metric discriminates
on all three, and by far the most on `preset_shape` — the one where the criterion is doing the
deciding.

### 16.1 The recommendation on record was multiply. The measurement says no.

`error(B)` at `B = 191`, `off` against `multiply` on the same arm:

| target | arm | `off` | `multiply` |
|---|---|---|---|
| near-field | within | 0.11622 | 0.11543 |
| near-field | between | **0.10110** | 0.10585 |
| near-field | frac_hot_between | **0.09313** | 0.11554 |
| `deep interior` | within | 0.07910 | 0.07898 |
| `deep interior` | between | **0.07958** | 0.08048 |
| `deep interior` | frac_hot_between | **0.06220** | 0.07003 |
| `preset_shape` | within | 0.12471 | 0.12471 |
| `preset_shape` | between | **0.12477** | 0.13218 |
| `preset_shape` | frac_hot_between | **0.07038** | 0.13133 |

**Multiply never helps and it wrecks the best criterion** — `frac_hot_between` on `preset_shape`
goes 0.07038 to 0.13133, from beating greedy's neighbourhood to worse than the random band's
*upper* edge. On the `within` arm it is a wash to five digits, which is the only place it does no
harm, and that is the arm that was already the worst criterion tested.

**`replace` is worse still, and it is not a second data point.** `signal_with(_, _, Replace)`
discards both arguments, so `replace × within`, `replace × between` and `structure_only` are the
same expression; their curves match to five digits because they *are* one row. Documented at the
enum. As `structure_only` it is the **worst row in the `preset_shape` table** — 0.13348 at
`B = 191` against a random *high* of 0.10814. Ranking on structure alone is worse than ranking at
random badly.

So §2.2's answer is **neither**. The structure term does not replace the signal and does not
usefully multiply it.

### 16.2 The winner is the criterion the instruction would have deleted

**`frac_hot_between/median`, with the structure term off**, beats the random band at nearly every
budget in all three targets:

| target | `B = 47` | `B = 191` | `B = 767` | `B = 1535` |
|---|---|---|---|---|
| near-field | **0.10457** / 0.10334 | **0.09313** / 0.09474 | **0.08046** / 0.08273 | 0.06914 / 0.06707 |
| `deep interior` | **0.06879** / 0.07140 | **0.06220** / 0.06258 | **0.04750** / 0.05173 | **0.02772** / 0.03983 |
| `preset_shape` | **0.09565** / 0.09970 | **0.07038** / 0.08649 | **0.04654** / 0.06676 | **0.03203** / 0.05044 |

(criterion / best random.) On `preset_shape` it is close to `greedy_oracle` itself — 0.07038
against 0.06881 at `B = 191` — which is the strongest showing any criterion has made on this
project, on the one tree the criterion actually controls.

**It does this on 31, 65 and 64 distinct values, with modal shares of 83.1%, 33.9% and 40.4%.**
Near-field's ordering has thirty-one distinct values over 5461 quads and an 83% mode, and it still
wins. That is the standing rule at full strength: *signal resolution is not what makes a ranking
good*. Meanwhile `within/median` carries 5418 distinct values in near-field and is beaten by
random at every budget — a fine-grained ordering that is actively bad.

**And this is the criterion that reads the ABSOLUTE mask** — the one the instruction to "make the
threshold relative" would have replaced. §14.5 caught that before it shipped. The relative mask
desaturated the spatial fields exactly as intended, and every criterion built on it still loses to
a saturated 31-valued count.

### 16.3 The threshold-free control does not rescue the family

`grad_rms` has **5461 distinct values — every quad distinct — and no threshold in it at all**. It
sits mid-pack: better than `within`, worse than `frac_hot_between`, worse than random. `layout_rel`
— the desaturated mask's own criterion — carries 46 / 174 / 214 distinct values and is the second
best non-`frac_hot` row on `preset_shape` (0.08970) while still losing to random there.

So the mask family's problem was never only the saturation. **Desaturating was necessary and is
not sufficient**, and nothing built on the spatial layout has yet beaten a plain count in the tail.

## 17. The slippy map

### 17.1 2:1 balance, and the share of the budget it costs

`Decision::BalanceForced` (code 10) is a separate decision so the geometry cost is **countable**;
a `Split` would be indistinguishable from a criterion-driven one. Measured on `deep interior`
under a camera: unbalanced 46 leaves with a worst adjacent gap of **2**, balanced 64 leaves with a
gap of **1**, and **14.1% of the quads computed** were balance-forced.

The unbalanced arm is the control, and it is not decorative: under the veto near-field reaches a
complete tree at one depth, where 2:1 holds trivially and the test would pass having measured
nothing.

### 17.2 The camera enters the priority and never the veto

`Camera::relevance(cx, cy, half, margin)` is the visible-area fraction, computed at query time.
`Camera::veto` stays position-free, which is what keeps a quad's *decision* independent of where
the camera points — the standing rule. Both halves are asserted: the veto returns the same
decision for a camera sitting on the quad and one ten units away, while relevance goes 1.0 → 0.0
and is strictly graded (0.5) for a quad straddling the edge.

**A pan now means something.** Without the bias, panning `cx` 1.00 → 1.04 leaves the tree
identical — the standing identity, now asserted as one rather than reported as a result. With it,
52 → 55 leaves. `margin` is §4.3's honest baseline: widen the viewport, drop prediction. Velocity
extrapolation fails on flick-and-stop, and the swept-path variant must beat this before its
complexity is justified.

### 17.3 The persistent frontier, and the reference kept beside it

`src/frontier.rs`: priority split into a **stored** physics term and a **derived** camera term,
bucketed into 24 log-spaced bands. Log-spaced because the signal spans six orders across regions;
linear bands would put a whole region in one bucket, which is the saturation failure this project
has already met twice.

`Frontier::rebuild` — the from-scratch path — is **kept permanently** and
`agrees_with_rebuild` compares them after inserts, cross-band reprioritisations and removals. The
failure mode is staleness, and it is invisible: a quad sitting high on a priority it no longer has
looks exactly like a bad criterion. Two paths that must agree is the only thing that catches it.

`band_of` conflated `NaN` with `+inf` under one `is_finite` guard — undetermined and
maximally-important sent to the same bucket. Nothing in the current signal produces `+inf`, which
is exactly why it would have sat there. `NaN` goes to the bottom, `+inf` to the top.

### 17.4 Zoom-out is free, but not in the form the brief states

§4.5 asks for *"the count of newly-computed quads after a zoom-out is ≈ 0"*. That presupposes a
tree persisting across frames, and this build deliberately has none — the scope discipline is *no
eviction, no caching, no async, no promotion*. Asserting it as written would require building the
thing the scope forbids; asserting it against a from-scratch descent would measure nothing.

What is available is the arithmetic underneath: measured, a zoomed-out descent computes **537
quads against 537 of the zoomed-in run's 597**, and **zero** of its boxes are absent from the
zoomed-in descent. So a persistent tree would have to compute none of them. Stated that way rather
than as the claim the build cannot support.

### 17.5 The coarse-ancestor fill was a missing filter, not a missing feature

`adaptive::render` drew only leaves, so a leaf without computed samples left raw background — a
hole, which reads as *"nothing here"* rather than *"not yet resolved"*, and is the one outcome
worse than a blocky texel. It now draws **every node with samples, coarsest first**, which is
§4.5's option 1 exactly: draw the coarse ancestor and let it sharpen.

**Wherever the tree is complete this is bitwise identical**, because leaves tile the root and a
parent is wholly overwritten by its children. It differs only where a leaf is missing — a
budget-exhausted quad truncated before compute, or a frontier the camera has just revealed. The
returned `LeafTexel` list stays leaves-only; including the fill would have doubled its rows and
halved the apparent texel size at every level.

## 18. The sweep — and PR #18 never ran anything with its own machinery enabled

All 69 dumps committed by PR #18 carry `tau_display=1e-4  hot_rule=q[0.50]  structure=off
k_frac=1  criterion=within`: the pre-fix configuration with new columns attached. `tau` at the
0.4th percentile still gated the split, `k_frac = 1` took the top 100% of the frontier so the
ranking changed nothing, and neither new signal was in play. **There was no "after" in the
corpus.** Everything below writes to `results/sweep/`; nothing existing is touched.

One correction to that diagnosis: `mode=balanced, k_frac=1` is *not* `Mode::Uniform`. Uniform
returns `Split` unconditionally and bypasses the `tau` and `alpha` gates; balanced at `k_frac = 1`
applies both and only declines to truncate. What never engaged was the **rank truncation**.

### 18.0 The wiring check found a bug in the thing being swept

Run before the sweep, exactly because a knob that is not plumbed through produces identical trees
at every setting and reads as *"the criterion cannot be fixed"*. All four knobs reached the tree —
and one reached too far.

**`k_frac` was truncating the bootstrap.** Levels below `bootstrap_levels` split unconditionally
because level 0 has no parent and therefore no `alpha`: there is no signal to rank them by. The
ranking demoted them to `Keep` anyway, so the tree never reached the depth where the criterion
could decide anything.

The tell: `near-field`, `deep interior` and `preset_shape` returned **byte-identical** leaf counts
and depth variances at every `k < 1` — `16/1/0.000`, `10/2/0.160`, `7/2/0.245`. Three unrelated
charts agreeing to the digit is chart-independent arithmetic, not physics. The `split` column
said so too: rows read 2 and 3 where the bootstrap alone requires 5. Fixed, and
`tests/criterion.rs::k_frac_never_truncates_the_bootstrap` asserts both the split count and that
no quad below the bootstrap is left a leaf.

### 18.1 Stage 1 — `tau × k_frac`. It works on one target of three

`near-field`, `structure=off`, `criterion=within`:

| `tau` | `k_frac` | leaves | levels | **depth var** | %max | veto |
|---|---|---|---|---|---|---|
| 1e-4 | **1.00** | 412 | 5 | **1.015** | 61% | **252** |
| 1e-4 | 0.50 | 103 | 5 | 1.900 | 35% | 36 |
| 1e-4 | **0.25** | 46 | 5 | **2.053** | 17% | **8** |
| 1e-4 | 0.10 | 31 | 5 | 2.046 | 13% | 4 |
| 1e-3 | 1.00 | 259 | 5 | 1.347 | 53% | 136 |
| **3e-3** | any | 16 | **1** | **0.000** | 100% | 0 |

Depth variance doubles and the veto share falls 61% → 13%: the tree stops being cap-decided and
becomes criterion-decided. `tau ≥ 3e-3` collapses it to one level at every `k` — §14.2's
upper-side failure, landing where the percentile table put it.

**The other two targets are inert across all twenty cells.** `deep_interior` returns `22/3/0.614`
at every `k`; `preset_shape` returns `16/1/0.000`. Neither knob can reach them, structurally:

- **`tau` cannot gate `preset_shape`** — its leaf-spread median is `2.86e-1`, 3400× above the
  largest `tau` swept, so every quad clears the gate everywhere.
- **`k_frac` has nothing to rank** — it truncates the set that already decided to *split*, and
  `preset_shape` produces **zero** splits past the bootstrap while `deep_interior`'s frontier is
  1–2 quads a round, where `ceil(1 × 0.1) = 1` truncates nothing.

### 18.2 Stage 3 — `alpha` is what binds them, and it still cannot move `preset_shape`

Swept because stage 1 showed the specified knobs could not reach two of three targets. Depth
variance at `tau = 1e-4, k = 0.25`:

| target | `alpha_hi` 0.5 | 0.2 | 0.1 | 0.0 | **−1.0 (gate off)** |
|---|---|---|---|---|---|
| near-field | 0.166 | 2.053 | 2.140 | 2.140 | 2.109 |
| `deep_interior` | 0.166 | 0.614 | 0.614 | **2.311** | **2.410** |
| **`preset_shape`** | 0.000 | 0.000 | 0.000 | 0.000 | **0.000** |

`alpha_hi = −1.0` is the **degenerate control**: every finite `alpha` clears it, so the gate is
effectively off. `deep_interior` needs `alpha_hi ≤ 0` to unlock at all — a 3.8× jump in depth
variance between 0.05 and 0.0.

**`preset_shape` is flat even with the gate off**: 16 leaves, one level, 5 splits (bootstrap
only), 8 `Floor` + 8 `Keep`. At `alpha_hi = −1` a `Floor` requires `alpha < −1`, i.e. the child's
spread is **more than twice the parent's**. So on half its quads refining makes the spread *grow*,
and on the other half `alpha` is not computable at all. No threshold on a convergence exponent can
help where there is no convergence.

### 18.3 Stage 2 — the criterion moves them, through the RANKING

I predicted stage 2 would be inert on those two for the same reason stage 1 was. **Wrong**:
`priority()` reads `signal_with(criterion, agg, structure)`, so with `k_frac < 1` the criterion
decides *which* quads get the budget, not only whether they pass a gate.

`preset_shape`, `tau = 1e-4`, `k = 0.25`:

| structure / criterion | leaves | levels | depth var |
|---|---|---|---|
| `off` / `within` | 16 | **1** | **0.000** |
| `off` / `frac_hot_between` | 31 | 4 | 1.193 |
| `off` / `layout_rel` | 28 | 4 | 1.167 |
| **`off` / `grad_rms`** | **31** | **5** | **2.046** |
| `multiply` / `within` | 16 | 1 | 0.000 |

**So a configuration does produce a selective tree on every target, and the knob that does it is
the criterion acting through the ranking** — not `tau`, not `alpha`. `grad_rms`, the
threshold-free control on the whole mask family, is what unlocks the one chart the others cannot.
No single criterion wins everywhere: `within` is best on near-field (2.053), `frac_hot_between`
and `layout_rel` on `deep_interior` (1.925, 2.046 under `multiply`), `grad_rms` on `preset_shape`.

Two controls firing correctly inside stage 2. **`replace` collapses to identical rows across all
four criteria** on every target — 40/5/1.569, 25/3/0.560, 28/4/1.167 — which is the documented
structural identity (`Replace` discards the criterion) confirming itself. And **`multiply/within`
drives both `preset_shape` and `deep_interior` back to `16/1/0.000`**, which is the `NaN` structure
term propagating through the product on an empty mask, exactly as the doc comment on
`signal_with` warns.

### 18.4 The dumps are self-describing

93 dumps in `results/sweep/`, named
`<target>__tau<t>__k<k>__struct-<s>__crit-<c>.prnq`, so a directory listing is a settings table
and the corpus can be re-derived by parsing filenames. The header carries the settings too.

## 19. The default was the control, and every headline render was made with it

`SchedCfg::default().k_frac` was **1.0** through PR #21. `Mode::Balanced` at `k_frac = 1` computes
the priority, sorts the frontier, and then refines all of it — the ranking runs and changes
nothing. §18 said this about the 69 dumps and then left the constant where it was, so the fix
landed in the sweep and nowhere else. Every dump in `results/charts`, `results/criterion` and
`results/vertical` carries `k_frac=1`, and so does every image derived from them.

The default is now `K_FRAC_RANKED = 0.25`, `K_FRAC_UNRANKED = 1.0` is a named constant rather than
a bare literal, and `scheduler::assert_not_uniform_in_disguise` refuses to let an example write
into `results/` from the degenerate cell. The old corpus is not touched: passing `1.0` reproduces
it bitwise and still lands in its own directory. New runs land in `results/charts_ranked`,
`results/criterion_ranked` and `results/animated_ranked`.

**One correction to the diagnosis as received.** `results/glsl/` was *not* made at `k_frac = 1` —
PR #21 (`0bc00ed`) re-rendered all sixteen files at `k_frac = 0.25, criterion = grad_rms`, taking
`shape` from 181 leaves to 31 and `prho` from 2695 to 49. The directory holds no `.prnq`, which is
why a parse of the dumps could not see it; it is now the one committed image set that was already
ranked. Everything else in the diagnosis stands.

### 19.1 The widened sweep — `k_frac` is the knob and `tau` is nearly inert

`tau ∈ {1e-4, 1e-3, 1e-2}` × `k_frac ∈ {1, 0.5, 0.25, 0.1, 0.05}`, three targets, `structure=off`,
`criterion=within`, budget 40000, viewport 1024², `alpha_hi = 0.2`. Near-field, where the knobs are
live:

| tau | k | quads | leaves | levels | depth var | %max | veto | rho_lvl |
|---|---|---|---|---|---|---|---|---|
| 1e-4 | **1.00** | 549 | 412 | 5 | **1.015** | 61% | 252 | **-0.295** |
| 1e-4 | 0.50 | 137 | 103 | 5 | 1.900 | 35% | 36 | -0.028 |
| 1e-4 | **0.25** | 61 | 46 | 5 | **2.053** | 17% | 8 | **+0.265** |
| 1e-4 | 0.10 | 41 | 31 | 5 | 2.046 | 13% | 4 | +0.137 |
| 1e-4 | 0.05 | 33 | 25 | **4** | 1.334 | 16% | 0 | +0.108 |
| 1e-3 | 1.00 | 345 | 259 | 5 | 1.347 | 53% | 136 | -0.320 |
| 1e-3 | 0.50 | 133 | 100 | 5 | 1.866 | 32% | 32 | -0.013 |
| 1e-3 | 0.25 | 57 | 43 | 5 | 1.691 | 9% | 4 | +0.264 |
| 1e-2 | any | 21 | 16 | 1 | 0.000 | 100% | 0 | NaN |

A whole decade of `tau` at fixed `k = 0.5` moves depth variance **1.900 → 1.866**. `k` from 1 to
0.25 moves it **1.015 → 2.053**. At `tau = 1e-2` every rung is `16/1/0.000` — the threshold has
gone above the bulk and keeps everything, the upper-side failure §14 named.

**There is an over-sparse end, and 0.25 is where the sweep peaks.** Depth variance falls at
`k = 0.05` and the tree loses a whole level (4 distinct against 5). The default is the peak of the
measured curve, not a value that made a picture look right.

`deep_interior` reads `29/22/0.614` at every `k` and `preset_shape` `21/16/0.000` at every `k`, as
§18 recorded: `k_frac` truncates the set that already decided to split, and neither produces enough
splits per round for a fraction to bite. `alpha` is what binds those two.

### 19.2 `rho(depth, spread)` is confounded twice, and the sign flip survives the form that is not

The naive statistic — leaf depth against the leaf's own spread — reads **-0.817 at `k = 1` and
+0.821 at `k = 0.25`**, which looks like the ranking reversing where the budget goes. It is not
readable: **refining a quad reduces the spread of the pieces it becomes**, so a deep leaf has a
small spread partly *because* it was refined.

Substituting the **parent's** spread removes that arm and leaves a second. `ensemble_spread` is a
spread over copies jittered within the **cell**; the cell halves every level and the measured
inter-level spread ratio runs 1.19–1.62 (§11). So a deep leaf's parent is a fine quad with a
systematically smaller spread than a shallow leaf's parent, and the correlation reads the
estimator's level-dependence. Measured, it is **negative at every `k`** — -0.419, -0.813, -0.472,
-0.348, -0.447 — including the rungs where the ranking demonstrably works.

**The form with neither confound is blocked by level.** Within one level every quad has the same
cell width, so the spreads are comparable, and the question is asked directly: among the quads at
level `L`, did the ones the descent split have the higher spread? Spearman of `was_split` against
`spread`, per level, pooled by quad count. That is the `rho_lvl` column above: **-0.295 → +0.265**
across `k = 1 → 0.25`, a third of the naive magnitude and the same sign change. It is `NaN` wherever
a level has one outcome only, which is every degenerate row — a level nothing was split at has no
correlation, and returning 0 there would read as "no relationship" where the truth is "one axis
does not vary".

### 19.2b The 26-chart gallery re-run: the criterion decides 1.5% of leaves before and 78% after

`results/charts_ranked/`, same command as the committed gallery with `k_frac = 0.25` in place of
`1.0`. Stop reasons read from the `decision` column of the `.prnq` dumps, not from the `bound`
column, which has been wrong before:

| corpus | charts | quads | leaves | Floor | Keep | MaxRelDepth | veto share |
|---|---|---|---|---|---|---|---|
| `charts/` (`k_frac = 1`) | 26 | 100,470 | 75,359 | 969 | 154 | 74,236 | **98.5%** |
| `charts_ranked/` (`k_frac = 0.25`) | 26 | 1,922 | 1,448 | 147 | 985 | 316 | **21.8%** |

**The standing result "the chart families do not exercise the criterion" was a fact about
`k_frac = 1`, not about the charts.** `MaxRelDepth` stopped 95%+ of leaves on 23 of 26 charts and
100% on three; at the ranked default it stops 21.8% and the criterion decides the rest. The
`screen` column is **0 on all 26 rows**.

`Keep` conflates two criterion decisions by design — the spread gate declining, and the ranking
deferring an outranked quad — because a deferred quad is deliberately `Keep` rather than
`BudgetExhausted`: it was outranked, not refused for want of budget, and conflating those two would
hide the ranking inside the column that exists to expose it. The 985 is not separable further from
the dump.

**And the `_uniform*` panels are not regenerated, for a reason worth stating.** That block builds
its own `res × res` slice and evaluates it directly — no quad, no tree, no decision enters it — so a
ranked run reproduces `results/charts/*_uniform*.png` bit for bit while costing `res² × (E+1)`
trajectories, 8.4M per chart at 1024 and about **95% of the run**. If that block's output moved with
`k_frac`, something would be very wrong. The committed ones stand for both corpora.

### 19.3 Selective, or merely sparse

64 leaves against 1755 is not a result on its own — under-refining everywhere produces the same
headline. `src/metric.rs` already answers it: build the fully-refined reference once, then score
each tree by `Cache::error_of` against `greedy_oracle` and a `random` band **at the tree's own leaf
count**. Below the band, the small budget went to the right places. Inside or above it, the tree is
sparse and the depth variance is a picture of under-refinement.

`greedy_oracle` is a reference and deliberately not a ceiling — measured at `t = 20` in near-field
it plateaus at 0.00048 while `first_div` reaches 0.00000, so a tree beating it indicates lookahead
value. Nothing asserts it dominates.

Near-field, reference tree complete to level 5, `N = 8`, `res = 256²`, `t = 13`, five random seeds:

| `k_frac` | leaves `B` | tree error | `greedy@B` | `random@B` | verdict |
|---|---|---|---|---|---|
| 1.00 | 223 | 0.05841 | 0.05670 | 0.06301–0.07380 | **SELECTIVE** |
| 0.50 | 76 | 0.07312 | 0.06635 | 0.07195–0.08633 | in band |
| 0.25 | 40 | 0.07539 | 0.07117 | 0.07502–0.09083 | in band |
| 0.10 | 28 | 0.07646 | 0.07441 | 0.07884–0.09196 | **SELECTIVE** |
| 0.05 | 25 | 0.07667 | 0.07441 | 0.07884–0.09196 | **SELECTIVE** |

**No rung is sparse.** Every ranked tree is at or below the random band at its own leaf count, and
the gap to `greedy_oracle` is small throughout — 0.07539 against 0.07117 at `B = 40`. So the small
tree is not under-refining everywhere; the budget it does spend goes where an ordering should send
it.

**And it is not a free lunch, which the depth-variance table alone would not say.** The tree error
rises monotonically as `k` falls — **0.05841 → 0.07667** — because the tree is displaying less. The
selectivity is in the *shape* of the tree (depth variance, criterion-bound stops, the level-blocked
`rho` turning positive), bought at a real cost in displayed error, because fewer quads are computed.
The correct reading of `k_frac` is a **budget-quality trade**, not an improvement at fixed cost.

Two things this cannot say. `error = 0` would mean "matches this sampling", not "correct" — the
reference is the complete tree at one sample per pixel, and at the floor which side of a filament a
pixel lands on is an accident of where its sample fell. And the descent is **capped at the
reference's own depth** (`max_level: Some(levels)`), because `Cache::error_of` is defined over a
leaf set that *tiles* the root: a leaf deeper than the cache has no entry and the number would be an
average over a hole. Without the cap the run correctly scored nothing — 170 of 223 leaves outside at
`k = 1` — rather than dropping them quietly.

**The mapping between the two trees is the one joint where they meet, and it was wrong.** A cell
centre sits at `(2i+1)h` from the low edge, so dividing by the cell width `2h` gives `i + 0.5` and
`.round()` lands on `i + 1` — every quad mapped to its right/upper neighbour. It was caught by the
reconstruction check that verifies the recovered index reproduces the centre, not by the numbers
looking wrong: without that check this would have scored a perfectly coherent leaf set belonging to
a shifted tree.

### 19.4 The §5 acceptance test, with a control that discriminates for the first time

`balanced_march` ran at the old default, so its "balanced" arm was uniform mode with the gates
still applied — **the treatment and the control were nearly the same tree**. Worse, the rank
truncation in `descend` ran regardless of mode, so passing `k_frac < 1` truncated the uniform arm
too: near-field at `t = 4` read **40 leaves and depth variance 0.6900 under both arms, to four
digits**, with the budget never exhausted. Two arms agreeing to the digit is the same tell as three
unrelated charts agreeing, one level up. `Mode::Uniform` is now exempt from the truncation, with a
two-armed test — the uniform tree identical across `k_frac`, the balanced tree different — because
a `k_frac` that reached nothing at all would pass the first assertion alone.

Near-field, budget 800, `N = 4`, viewport 64², `tau = 1e-4`, `k_frac = 0.25`:

| mode | t | leaves | depth var | churn | shared | screen |
|---|---|---|---|---|---|---|
| balanced | 4 | 40 | 0.6900 | – | 0 | 16 |
| balanced | 8 | 31 | 0.6514 | 0.3333 | 9 | 8 |
| balanced | 13 | 34 | 0.5744 | 0.2500 | 8 | 8 |
| balanced | 16 | 22 | 0.2314 | 0.5263 | 19 | 0 |
| balanced | 20 | 25 | 0.5600 | 0.0833 | 12 | 4 |
| **uniform** | 4 | **256** | **0.0000** | – | 0 | **256** |
| **uniform** | 8 | **256** | **0.0000** | 0.0000 | 256 | **256** |
| **uniform** | 13 | **256** | **0.0000** | 0.0000 | 256 | **256** |

`deep interior` reads the same shape: balanced 0.2496–0.6136 with churn to 0.3750, uniform
**0.0000 at every playhead** on 256 veto-stopped leaves.

The control reads **exactly 0.0000 at every playhead** and is stopped by the veto on all 256 leaves
rather than by the budget, which is the condition that makes it a control at all. The balanced arm
holds variance bounded away from zero **with churn nonzero** — a steady state rather than a frozen
one, which the variance plot alone cannot distinguish. Every churn row is annotated with its shared
count; at 8–19 shared quads these are thin and are labelled thin.

### 19.5 The magenta is four footprints, not 1426 pixels

`results/glsl/shape.png` carries 1046 `DEBUG_NAN` pixels and `plambda.png` 380, which reads as
scattered failure in the region of interest. The adaptive render is **nearest-neighbour**: one
footprint of a level-2 quad paints a texel roughly `res / (4N)` on a side, so at `res = 512, N = 8`
a single undetermined trajectory paints ~16×16 = 256 pixels. Measured on the committed frames, the
magenta is **3 axis-aligned blocks** (18×19, 18×18, 20×19) in `shape` and **1** (20×19) in
`plambda`. Four footprints. **A pixel count of a debug colour is a fact about the texel size.**

`colour::rgb` has four exits to `DEBUG_NAN` and they are different findings. The census
(`examples/nan_probe.rs`, uniform grid, `t = 13`, f64):

| chart | samples | non-finite copy | SimFailed | DecodeFailed | non-finite shape | total |
|---|---|---|---|---|---|---|
| `preset_shape` | 4096 | 4 | 0 | 0 | 6 | 0.244% |
| `preset_shape` | 16384 | 12 | 0 | 0 | 21 | **0.201%** |
| `preset_plambda` | 16384 | 0 | 0 | 0 | 0 | 0.000% |
| `preset_prho` | 16384 | 0 | 0 | 0 | 0 | 0.000% |
| `preset_shape_pl` | 16384 | 1 | 0 | 0 | 0 | 0.006% |

**Zero decode failures and zero sim failures at both resolutions.** All of it is non-finite copies
and non-finite `shape_vec` — a triple collision — and the rate is stable across a 4× change in
sampling, which is the tell that it is a property of the chart rather than of the grid.
`preset_shape` is the one chart in this set whose terminations are collisions (escape fraction
0.0547 against 0.9894–1.0000 for the momentum slices), so it is the one that passes through
collision-adjacent shapes. A non-finite copy is a **measurement outcome** and is never discarded; a
`DecodeFailed` would have been the chart handing back something that is not a three-body state, and
there are none.

## 20. The allocation is inverted, and the baseline was never random

### 20.1 The strongest number in this run is a histogram, not an error

`near-field` at `B = 1535`, leaf levels as `level:count`:

| | histogram |
|---|---|
| `dp_optimal` (the exact optimum) | `3:2 4:71 5:589 6:472 7:16` |
| `within/median` | `1:2 3:22 4:11 5:29 6:102 7:984` |

The error cell reads `0.14956` against `0.10664` -- somewhat behind, about 55% of the achievable
improvement forfeited. **The histogram says it is doing the opposite thing.** The optimum puts
nothing below level 3 and only **16** leaves at level 7; `within/median` leaves **2 leaves at
level 1** -- a quarter of the image at one texel each -- and drives **984** to the bottom. Its
allocation is not merely worse than the optimum, it is *inverted*.

No error cell showed that, and no error cell could: `error(B)` is one scalar per budget and the
failure is in the shape of the spend. **Read the leaf histogram beside every curve.** This is the
argument for changing the default criterion, and it is far stronger than any error digit.

### 20.2 And the honest baseline is breadth-first, not random

Every table in §10 through §19 read *beats random* as though random were the thing to beat. It is
not. **`Rank::Uniform` -- refine the shallowest quad available -- is a far stronger baseline, and
against it most of the corpus's conclusions change sign.**

`error(B)`, the best criterion against uniform, at each budget:

| region | B where a criterion beats uniform | by how much |
|---|---|---|
| `far` | **never** | uniform *is* the optimum; every criterion ties it or loses by `1e-5` |
| `near-field` | only `B >= 6143`, and only `term_grad/median` | `0.08439` against `0.08480`; `0.05255` against `0.05282` |
| `deep interior` | `B >= 767` for `frac_hot_between`, decisively at `B >= 6143` | `0.04035` against `0.04813` at `B = 6143` |

`frac_hot_between/median` -- the best criterion measured on this project -- **never beats uniform
in `near-field` at any budget**, at `0.11148` against `0.10984` at `B = 1535`. It was well clear of
the random band there. The random band was the wrong bar.

This is not a defect in the criterion so much as a missing row: on a field whose spread tracks cell
width, ranking by spread *is* breadth-first (§20.7), so a criterion can be far above random and
still buy nothing over refining uniformly. The comparison that decides anything is
**criterion against uniform**, and it is now in the table.

### 20.3 The gap rises with structure, and that is the result

Share of *achievable* improvement the best row leaves on the table, `(row - dp)/(root - dp)`, at
`B = 1535`:

| region | gap | what varies there |
|---|---|---|
| `far` | **0.0002** | nothing -- the field is one smooth gradient |
| `near-field` | **0.0336** | structure is **localised** |
| `deep interior` | **0.0999** | structure **everywhere** |

**A criterion can only earn its keep where structure is localised.** Where nothing varies there is
nothing to rank; where everything varies, the budget has to go everywhere and breadth-first is
close to right again. The headroom lives in between.

Which makes `far` degenerating **correct behaviour, not a defect**. It is the control that shows
what a featureless field looks like, and every criterion matching the optimum there is the right
answer -- not a region where the criteria mysteriously tie.

### 20.4 The gap is a curve, and a single budget hides its shape

`(row - dp)/(root - dp)` for the best row at every budget:

| region | 11 | 47 | 191 | 383 | 767 | 1535 | 3071 | 6143 | 12287 |
|---|---|---|---|---|---|---|---|---|---|
| `far` | 0.9001 | 1.0000 | 0.3294 | 0.0664 | 0.0216 | 0.0002 | 0.0013 | 0.0001 | 0.0007 |
| `near-field` | 0.0699 | 0.0744 | 0.0563 | 0.0250 | 0.0451 | 0.0336 | 0.0781 | **0.1106** | **0.1241** |
| `deep interior` | 0.1659 | 0.1352 | 0.0631 | 0.0772 | 0.1212 | 0.0999 | 0.1259 | 0.1157 | 0.0899 |

And greedy's, for the same budgets:

| region | 11 | 47 | 191 | 383 | 767 | 1535 | 3071 | 6143 | 12287 |
|---|---|---|---|---|---|---|---|---|---|
| `far` | 0.0002 | 0.3886 | 0.4201 | 0.5881 | 0.7623 | **0.7693** | 0.6744 | 0.5755 | 0.3147 |
| `near-field` | 0.0000 | 0.0000 | -0.0000 | 0.0000 | 0.0004 | 0.0004 | 0.0011 | 0.0038 | 0.0046 |
| `deep interior` | 0.0000 | 0.0000 | 0.0047 | 0.0045 | 0.0168 | 0.0301 | 0.0459 | 0.0752 | 0.0939 |

**The shape is region-dependent and the widening is a `near-field` phenomenon.** There the gap is
noisy through `B = 1535` and then rises steadily to **0.1241** -- the criterion is adequate at the
first few hundred splits and progressively worse at the later ones, which fits the inverted
histogram exactly: the failure is **late-stage allocation**. In `deep interior` the gap peaks at
`B = 3071` and *falls* to 0.0899; in `far` it collapses to zero by `B = 1535` and stays there. Any
one budget quoted alone would have supported a different story in each region.

`far`'s greedy row is worth reading on its own: **1.0000 at `B = 23`** -- greedy at that budget is
exactly as good as not refining at all -- rising again to 0.7693 at `B = 1535`. It is the only row
anywhere in the corpus that stays bad as the budget grows.

### 20.5 Criterion against uniform, per level

Read at the budgets where uniform **completes a level exactly**, `B_d = 1 + 4(4^d - 1)/3`. At any
other budget uniform sits mid-level and the comparison would be scoring where its partial row
happened to stop. `captured = (uniform - row)/(uniform - dp)`: **1.0 is the optimum, 0.0 is no
better than breadth-first, negative is worse than doing nothing clever.**

| level | B | `far` | `near-field` | `deep interior` |
|---|---|---|---|---|
| 2 | 21 | 0.0000 | `uni==dp` | **-9.38** |
| 3 | 85 | 0.0000 | 0.0000 | **-1.97** |
| 4 | 341 | 0.0000 | **-0.2988** | **-0.4497** |
| 5 | 1365 | `uni==dp` | **-0.7210** | -0.0394 |
| 6 | 5461 | `uni==dp` | 0.0019 | **+0.5296** |

**The gap is not flat across levels; it rises monotonically with depth**, and crosses zero only at
the deepest discretionary level, and only in `deep interior`. Shallow splits the criterion chooses
are actively worse than taking them in raster order.

That is the shape a depth-dependent strategy would exploit -- and **the same table is the argument
against adding one.** The crossing is at level 6 in `deep interior`, not reached at all in
`near-field`, and undefined in `far`. A depth parameter would have to be tuned per region, which is
the tunable-constant defect this project has been bitten by four times, and it is the same argument
that killed a fixed `tau`: a constant cannot land inside a distribution that moves between regions.

Two caveats, because the cell is generous to the criteria as printed. The `best criterion` column is
picked **per level in hindsight** -- `grad_rms` at level 6 in `deep interior`, `term_grad` at level
6 in `near-field`, `frac_hot_between` everywhere else -- so it is an oracle over criteria and still
loses to breadth-first below level 6. And where `uniform == dp` to machine precision the ratio has a
denominator at the arithmetic floor and describes nothing; it prints as `uni==dp` rather than as a
number, which is the difference between "both sides are right" and "both sides are dead".

### 20.6 The bound that can fail, and it holds

`Cache::dp_optimal` is the exact minimum over **all** tree-shaped leaf sets at a given budget, by
tree DP: `f_k(0) = err_sum(k)`, `f_k(s) = min over s0+s1+s2+s3 = s-1 of sum_i f_ci(si)`. Budget and
splits are locked to `replay`'s own accounting, `B = 1 + 4s`. *No ranking may beat it* -- an
assertion with a real failure mode, and `tests/criterion.rs` runs it over every `Rank` at every
budget.

`greedy_oracle` was read as an upper reference in every table since §10 and is not one. On `far` at
`B = 1535` it reads **0.54760** against a random band of **0.48550-0.52047** and every criterion at
**0.36557** -- the worst strategy in the table, under a name that says it cannot be. It is renamed
**`greedy_lookahead_1`** throughout the code; the numbers already recorded above keep the old name,
because renaming them would rewrite the record of what was measured.

Two invariants were proposed to find the fault and **neither could have fired.** `Cache::error_of`
is a sum of `err_sum` over the leaf set and `Cache::gain` is parent-minus-children, so splitting
`k` replaces `err_sum(k)` with `err_sum(k) - gain(k)` and the accounting identity
`error_of(leaves) == err_sum(root) - sum(gains)` telescopes for **any** ranking, any sequence, and
any values `err_sum` happens to hold -- random numbers included. And a choice check re-runs
`replay_with_leaves`'s own argmax over a pure, static `gain`. Both report PASS. *A test that cannot
fail is indistinguishable from a test that passes*, at the level of the metric this time.

**It reads at the complete tree, uncapped, in centiseconds.** The naive `O(quads x B^2)` reading is
wrong twice: the 4-way merge is three successive 2-way convolutions (`O(cap^2)`, not `O(cap^4)`),
and each node's split cap is bounded by its own subtree, so only the top two levels ever see the
full budget. Measured at `levels = 7`, 21845 quads, 5461 splits: **0.01 s**, against 349-1117 s to
build the cache. No cap is needed and none is applied. **Keep the ceiling in every `error(B)` table
permanently** -- at that price it converts "beats random" into "captures 95% of the achievable
improvement", which is a claim with a denominator.

**The bound holds everywhere.** Worst margin `row - dp` over every ranking and every budget:
`+0.000e0` (`far`), `-1.388e-16` (`near-field`), `-1.388e-17` (`deep interior`) -- summation order,
not violations. **The replay is sound and no `error(B)` number in the corpus is suspect.**

`criterion_metric` now computes the ceiling and **asserts** it. That assertion has been executed at
scale -- `4 8 1e-4 13` over all three regions, worst margins `-1.4e-17`, `-5.6e-17`, `+0.0e0` --
because an assert that has never run is the kind of gap that survives indefinitely. It writes to a
scratch root, which is now an **argument** rather than a hardcoded `results/`: a reduced-`levels`
validation pass would otherwise overwrite the committed 512^2 artefacts with a small raster, and a
small raster reads as a rendering fault rather than a stale file.

The committed `results/output/criterion_metric.txt` **predates both the ceiling row and the uniform
row**; it is unchanged, and §20 is where they are recorded until that table is next regenerated for
its own reasons. `results/output/oracle_audit.txt` carries every number in this section.

### 20.7 Why greedy fails on `far`: `err_sum` is flat until level 3

Per pixel, normalised by the root's, median over quads at each level:

| level | `far` | `near-field` | `deep interior` |
|---|---|---|---|
| 0 | 1.00000 | 1.00000 | 1.00000 |
| 1 | **1.00067** | 0.94264 | 0.76073 |
| 2 | **1.00016** | 0.81430 | 0.54474 |
| 3 | **0.99884** | 0.77844 | 0.43228 |
| 4 | 0.91387 | 0.69918 | 0.40615 |
| 5 | 0.61965 | 0.59543 | 0.38210 |
| 6 | 0.31066 | 0.46190 | 0.28214 |
| 7 | 0.00000 | 0.00000 | 0.00000 |

`far` is **flat through level 3** and its gains there are noise: `-3.022e-7` at the root,
`6.319e-7` at level 1, `-6.425e-9` at level 2 with **13 of 16 negative** -- against `8.634e2` at
level 3. So a level-2 split has to be paid for at zero or negative immediate gain to unlock a gain
five orders larger beneath it, and greedy declines. Once it has opened one subtree it descends
inside it to the bottom and never returns. That is the level-2 barrier, and it is the entire
failure.

`near-field` and `deep interior` fall gradually at every level, so gains are available immediately
everywhere and greedy has no barrier to decline. **The same statistic explains both the failure and
its absence**, which is what makes it the mechanism rather than a story.

### 20.8 `far` is a smooth field, not an amplified noise floor

`criterion_metric`'s AUTO-RANGED OVER NOISE guard did not fire on `far`, and `far` being
auto-ranged noise is a standing finding (§11). If the ramp were stretching a `x8` span of noise
across the full colour range, then the flat `err_sum` above would be **noise failing to resolve**
rather than physics being absent -- and the level-2 barrier would be a property of the metric on a
noise-dominated region rather than of the field. `far` is now load-bearing in §20.3, so it matters
which.

**It is a field.** Lag-1 neighbour correlation of the ramped scalar, by level:

| level | `far` | `near-field` | `deep interior` |
|---|---|---|---|
| 1 | 0.9434 | 0.3536 | -0.3332 |
| 2 | 0.9923 | 0.8392 | -0.0514 |
| 3 | 0.9984 | 0.8316 | 0.4781 |
| 4 | 0.9996 | -0.0253 | 0.7746 |
| 5 | 0.9999 | 0.6352 | 0.9054 |
| 6 | **1.0000** | 0.7652 | 0.9432 |

`far` is maximally coherent at every level. And its p1/p99 **halve exactly** with each level --
`8.453e-8, 4.211e-8, 2.105e-8, 1.051e-8, 5.255e-9, 2.627e-9, 1.314e-9`, a ratio of 2.000 -- which
is `spread ~ g*w` measured directly rather than argued. That is a real gradient of tiny magnitude,
not an amplified noise floor. **No conclusion from this run changes, and `far` stands as the
smooth-field control.**

It also says why the guard could not have fired. Both existing arms read **amplitude**, and
amplitude cannot separate a small real signal from noise. Worse, the absolute arm's floor is the
region's own median energy drift -- `4.478e-11` on `far`, so the floor is `4.478e-9` against a
`ramp.1` of `1.064e-8`. **The floor falls with the field it is meant to bound**: a ratio in
disguise, in a region tame enough that both go to zero together. A **third arm** now reads lag-1
coherence at level 3 and fires below `rho = 0.5`; `far` reads 0.9984 there and correctly does not
trip it.

The standing finding that inter-level spread ratios run `1.19-1.62, falling to 1.048` was measured
in a chaotic region and named 2.000 as the proportional control value. `far` reads exactly 2.000.
The two are consistent, and this is the control the earlier measurement was missing.

### 20.9 The leaf histograms in full

Leaf levels at `B = 1535`:

| | `far` | `near-field` | `deep interior` |
|---|---|---|---|
| `dp_optimal` | `5:982 6:168` | `3:2 4:71 5:589 6:472 7:16` | `3:15 4:99 5:200 6:724 7:112` |
| `greedy_lookahead_1` | `2:14 5:48 6:64 7:1024` | `3:3 4:70 5:577 6:468 7:32` | `2:1 3:18 4:93 5:134 6:584 7:320` |
| `within/median` | `5:982 6:168` | `1:2 3:22 4:11 5:29 6:102 7:984` | `1:3 2:1 3:3 4:3 5:41 6:119 7:980` |
| `frac_hot_between` | `5:982 6:168` | `3:2 4:14 5:870 6:264` | `1:1 2:4 3:10 4:15 5:40 6:984 7:96` |
| `random[1]` | `2:1 3:15 4:101 5:210 6:291 7:532` | (same seed, same tree shape) | (same) |

On `far`, greedy leaves **14 quads at level 2** -- most of the image -- while driving 1024 leaves
to the bottom. That is the uniform dark corner, counted.

### 20.10 On a smooth field, ranking by spread IS breadth-first

`far`'s thirteen non-greedy rows agree to five digits at every budget -- `within/median` on
**21845 distinct values** and `frac_hot_between` on **1**. They are not thirteen criteria
independently agreeing; they are **one allocation reached by two routes**, and the histogram proves
it rather than the digits: all three rows produce the *identical* leaf set `5:982 6:168`, which is
also `dp_optimal`'s and also `uniform`'s.

The routes are different. A quad of width `w` over a gradient `g` has internal variation `~ g*w`,
so **spread tracks cell size** -- measured in §20.8 at a ratio of exactly 2.000 per level -- and
argmax-on-spread picks the shallowest quad. A constant signal has no argmax at all and falls
through to the tie-break, which is lexicographic on `(level, ix, iy)` -- **level first**, so also
shallowest. Both are breadth-first, and on a smooth field breadth-first is within `2e-5` of exactly
optimal.

**So the criterion already self-adapts.** Where there is nothing to discriminate, the ranking
degenerates to uniform *on its own*. That is a property worth protecting rather than overriding,
and it is the second reason not to bolt a depth parameter onto it.

### 20.11 A split can make the image worse, and the prefix-min is why

`f_root(s)` is the best tree at **exactly** `s` splits; the ceiling at budget `B` is
`min over s <= (B-1)/4` of it, because a replay may stop early and because more splits only help if
every gain is non-negative. They are not. **The root's own gain on `far` is `-3.022e-7`** -- its
`err_sum` per pixel is `1.00000` against level 1's median `1.00067` -- so the prefix-min binds at
`s = 1`, measured, and nowhere else in any region.

This is the hypothesis that a parent's `N x N` sample grid and its children's are **different
approximation families**, confirmed directly. On a smooth field one 8x8 grid over the whole box is
marginally better than four 8x8 grids over quarters. Negative gains elsewhere: `far` 14 quads,
`near-field` 14, `deep interior` **102**. Reading `f_root(S)` instead of the prefix-min would have
quoted a ceiling above its own achievable minimum.

### 20.12 Four roles, stated on the table's face

The table had a floor and a reference and was read as though it had a ceiling. It now says which is
which, and it has gained the row that was missing:

```
floor      random lo/hi          several seeds, read as a band, never one trace
BASELINE   uniform               breadth-first -- the bar a criterion must actually clear
reference  greedy_lookahead_1    greedy on immediate delta-error -- NEITHER OPTIMAL NOR A BOUND
ceiling    dp_optimal            exact minimum over all tree-shaped leaf sets; asserted, not trusted
```


---

## 21. The concentric banding is a sync-cadence artefact, and it does not reach the near-field

Every committed latent chart carries fine dotted arcs across its smooth regions -- most clearly
`preset_plambda_uniform.png`, where they cover the whole field with spacing that widens outward.
That is the signature of contours of a smooth function at uniform value intervals, and it is a
**rendering** defect as well as the **measurement** defect already on record.

### 21.1 The mechanism, read out of the driver rather than guessed

**Collision is sampled inside the RK4 loop** (`driver.rs`, `tc = t + s.t`) and carries step
resolution. **Escape is sampled only at sync boundaries**, where the state is already Cartesian
and every trajectory shares a playhead -- the reference's cadence, transcribed. So `t_end` is
quantised to `n_sync` values **exactly where escape is the terminating event**, and is continuous
where collision is.

That predicts which images band, which is what makes it a mechanism rather than a story.

### 21.2 The decisive count, at `escape_every` 0 against 1

| target | escape frac | `t_end` distinct | on a sync boundary |
|---|---|---|---|
| `preset_plambda` | 0.9938 | **16 -> 2623** (164x) | **99.52% -> 0.26%** |
| `preset_shape_pl_h1` | 0.9721 | **41 -> 2316** (56x) | **98.60% -> 3.09%** |
| `deep interior` | 0.0945 | 2983 -> 3303 (1.1x) | 10.35% -> 3.25% |
| `near-field` | 0.0002 | **86 -> 86** | unchanged at every stride |

A quantised `t_end` does not merely take few values: **its values ARE the boundary times**, and
99.52% of `preset_plambda`'s land on one. That is the artefact demonstrated rather than inferred.

**`near-field` is bitwise unaffected.** `t_end`, `spread_event`, `spread_shape`,
`frac_hot_between`, `term_grad` and `within` are identical to the digit across strides
`0, 32, 4, 1`. Its escape arm is silent at `t = 13` and every termination is a collision, so
there is nothing for the cadence to quantise -- and the §20 and signal-audit results measured
there stand.

**One defect in the statistic itself, stated rather than corrected away.** `near-field` reads
97.85% "on a sync boundary" while being completely clean, because 97.8% of its footprints are
*Bounded* and sit at `t_end = t_max = 13`, which is itself a multiple of the sync interval. The
counter conflates "quantised by escape" with "reached the horizon". **Read the delta column, not
the level.**

### 21.3 The test originally proposed could not have fired

*"Recount `frac_hot_between`'s distinct values: 45 -> thousands means quantisation, 45 -> 45 means
Wada."* But `frac_hot_between` is `frac_above_tau_between`, a fraction over the quad's `N^2`
footprints, so at `N = 8` it can take **at most 65 distinct values by construction**. Its 45 is 45
of an arithmetic ceiling, and §16's own `31 / 65 / 64` is that ceiling showing. `45 -> thousands`
was unreachable under either hypothesis, so the test would have reported "it is Wada" whatever was
true upstream.

Measured anyway, and it moves the *wrong way*: under the finer cadence `frac_hot_between` gets
**more** saturated, not less -- `deep interior` 31 -> 12 distinct with the modal share 41.2% ->
**83.5%**, and the latent charts stay at 1-2 distinct. **The saturation is not quantisation.**

### 21.4 The real contamination is an outcome RE-LABELLING, and it is larger

`deep interior`'s `t_end` was already continuous. What the cadence changes there is the terminal
class: **escape 0.0945 -> 0.5494, collision 0.8965 -> 0.4482.** Half the region.

Under `stop_on_event`, a genuine early escape is not *noticed* until the next boundary; the run
keeps integrating, dips below `r_coll`, and a collision wrongly wins precedence. The finer test
corrects a **precedence error**, not a resolution one -- a different defect from the one the
banding pointed at, found by the same measurement.

### 21.5 The `d_min` discriminator was CONFOUNDED, and its answer was wrong

The first reading of whether a finer cadence is a fix or spurious mid-encounter firing used
`d_min_true` by terminal state: the re-labelled footprints carried **larger** separations
(escaped p50 `1.063e-3 -> 4.419e-3`), which was taken as evidence against firing during an
encounter.

**It is evidence of the truncation it was supposed to test.** `d_min_true` is a minimum over the
whole run, and a run stopped early by a spurious escape never reaches its close approach -- so its
`d_min` is larger *because* it terminated early. The statistic was confounded by the very effect
it was measuring. A direct test was needed and gives the opposite answer.

### 21.6 The escapes the finer stride adds are 100% transient

`escape_candidate` is relative energy `> 0` and receding, which during a close encounter is
transiently true. The test: take the trajectories that escape at `escape_every = 1` and **not** at
the reference cadence, and ask whether they are still unbound at later sync boundaries.

| region | escape only at stride 1 | still unbound at +1 / +2 / +3 / +4 / +8 boundaries |
|---|---|---|
| `deep interior` | **895** of 2304 | **0.000 / 0.000 / 0.000 / 0.000 / 0.000** |
| `near-field` | 1 | 0.000 at every window |
| `preset_plambda` | **0** | -- nothing to test |

Every one of the 895 has re-bound. Latching them took `deep interior`'s escape fraction from
0.0947 to **0.5494**, and all of that was invented. The `0.0945 -> 0.2153 -> 0.4423 -> 0.5494`
sequence was not converging -- its largest *relative* step was at the finest stride, which is the
wrong shape for a resolved quantity.

**And `preset_plambda` adds zero.** Where escape genuinely terminates, the finer stride changes no
labels at all and only sharpens the time. That split is what makes the guard safe.

**The first window sweep was confounded too, by this file's own instrument.** It re-ran to
`t_e + w` with `n_sync` rescaled per window, so every window was a different discretisation, and
it produced `0.162, 0.219, 0.011, 0.083, 0.335` -- read as "the escape condition flickers". It does
not. `AzOut::escape_flags` records candidacy at every boundary of **one** run at **one** step size,
and the answer is a flat zero. *`n_sync` fixed while `t_max` varies compares different
discretisations* -- this was the same defect inverted, in a diagnostic written to catch defects.

### 21.7 The guard, and what it turns the stride into

`AzOpts::escape_confirm` (default **on**) holds an **in-loop** escape provisional until the next
sync boundary and commits it only if the condition still holds -- with the **first crossing** as
the time, because the guard decides whether the event was real, not when it happened. Boundary
detections are the reference's own arm and are untouched. With `escape_every = 0` the guard is
inert, so nothing already measured moves.

Escape fraction across strides `0, 32, 4, 1`, guarded:

| region | 0 | 32 | 4 | 1 | `t_end` distinct, 0 -> 1 |
|---|---|---|---|---|---|
| `deep interior` | 0.0947 | 0.1268 | **0.1564** | **0.1564** | 2982 -> 3281 |
| `preset_shape_pl_h1` | 0.9721 | 0.9721 | 0.9721 | 0.9721 | **41 -> 2314** |
| `preset_plambda` | 0.9938 | 0.9956 | 0.9956 | 0.9956 | **16 -> 2099** |
| `near-field` | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 86 -> 86 |

**`preset_shape_pl_h1`'s labels are now stride-invariant while its `t_end` resolution improves
56x.** That is what a pure resolution gain looks like. `deep interior` converges at stride 4 and
4 -> 1 moves nothing, which is the convergence the unguarded sequence never showed.

**So the stride is a COST knob, not a correctness one** -- a coarse stride only delays detection,
it no longer changes which event wins or how many there are. That is the property that makes the
default safe to choose on cost grounds.

`deep interior`'s residual `0.0947 -> 0.1564` is the precedence repair arriving: those are escapes
that genuinely occurred mid-interval and persisted, but were pre-empted by a collision that fired
first under the coarse cadence. It is the same order as the directly counted precedence population
(5.21% of all trajectories), which is corroboration rather than a second effect.

### 21.8 Fix 1 measured alone: nearly inert in production, large on the reference path

`classify` used to rank collision above escape unconditionally, **discarding both times**, and
justified it as *"collision is sampled continuously, so it is the earliest thing that can fire"*.
Continuous sampling makes it the earliest **detected**, not the earliest **occurring**. Deciding by
`min(t)` -- with `t_end` set the same way, so the state and its time cannot disagree -- removes the
dependence on when each arm happens to be sampled.

Measured with the cadence untouched, it moves **one footprint of 5440**. Under `stop_on_event` the
loop breaks on the first *detected* event, so only one is ever recorded and `min(t)` has nothing to
compare. Its subject lives on the reference path, where both arms accumulate:

| region | fired both arms | escaped FIRST, labelled `collision` | median lead |
|---|---|---|---|
| `preset_plambda` | 996 | **990** (99.4%; **42.97% of all**) | 2.6201 = 6.45 intervals |
| `deep interior` | 1219 | 120 (9.8%; 5.21% of all) | 0.5084 = 1.25 intervals |
| `near-field` | 0 | 0 | -- |

**The ordering guarantee that makes the stride safe comes from sampling both arms at the same
cadence, not from `min(t)`.** `min(t)` is what stops the state and `t_end` disagreeing, and it is
what the reference path needed.

### 21.9 A spec change behind flags, with both defaults chosen from measurement

`AzOpts::escape_every` defaults to **0**, the reference's boundary-only cadence, and
`escape_confirm` defaults to **on** but is inert at that stride, so nothing already measured has
moved -- `cargo test --release` is **213 passed / 0 failed** with the Python cross-check green. Turning it on changes results, and the cross-check and the horizon table were both measured
coarse, so it is a spec change behind a flag rather than a tidy-up.

`tests/outcome_encoding.rs` asserts **both** arms, because a flag that does nothing passes as
easily as one that works: at the reference cadence an escape time must land exactly on a sync
boundary, and at the fine cadence some escaping trajectory's `t_end` must actually move off it.
It also fails if no trajectory in the test escaped at all -- Burrau's near-field could not have
exercised this, and a test whose subject never executes is decoration.

The guard carries **two** arms for the same reason. It must **cut** the escape count in
`deep interior`, where the transients are; it must **not** cut it on `preset_plambda`, where escape
genuinely terminates and the finer stride adds none. A guard that rejects everything passes the
first arm exactly as well as a correct one, and only the second tells them apart.

**Bisection within the firing interval is not implemented, and the counts say it is not needed at
these strides.** The stride alone takes `preset_shape_pl_h1` from 41 distinct `t_end` values to
2314 and drops the on-boundary fraction from 98.60% to 1.43% -- the step is already RK4-sized. The
confirming test is a render at 1024², which has not been run; the count evidence is not a
substitute for it and is not quoted as one.
### 21.10 The second diagnostic: the crisp edges are real, the banding is not

`preset_shape_pl_h1_uniform_outcome.png` is already committed. Under outcome-class colouring the
**banding vanishes entirely** -- consistent with it being a continuous-field artefact -- while the
**crisp polygonal edges survive and sharpen**: straight lines meeting at points, and a genuinely
*circular* boundary around the central fan. Those are outcome-class boundaries, so they are real
regime structure and not a colouring artefact.

A circle plus radiating wedges is **polar structure in the chart plane**: a radius threshold and an
angle threshold. Saturation is the candidate that would produce it. **Stated and stopped there** --
two image diagnoses on this project have already been settled by one targeted measurement, and
speculating past the render is how the earlier ones went wrong.


---

## 22. The closure criterion: check 1 passes outright, and the 383x gap does not reproduce

`examples/escape_closure.rs`, `examples/closure_render.rs`, n = 24 per side, `eta = 1e-2`,
`tau = 1e-3`, `closure_k = 1`. Raw output in `results/output/escape_closure.txt`.

The criterion, transcribed from `reference/escape_criterion.py`:

```text
ESCAPE  <=>  |dn| over a window < tau    AND    E_rel > 0
```

**Four details are read off the reference and three of them change the design.** `dn` is a chord
between the two **ends** of the window -- `nbuf` samples are buffered and only `buf[-1]` and
`buf[0]` are used -- so boundary sampling is a *transcription* rather than an approximation: at
`t_max = 13, n_sync = 32` the realised window is **0.406 against the reference's 0.400**. Closure
gates once per trajectory and energy per body (`dn[...,None]` against an `E` of shape `(...,3)`).
Body selection is `np.argmax(fire,-1)`, the **lowest firing index** -- not tightest-pair-first, which
is what the `Distance` rule does. And there is no latch beyond `esc_body < 0`: the window *is* the
persistence guard, and it cannot fire before `t = win` because the buffer is not full.

### 22.1 The window is a TIME, and it is 100-1000x the shortest inner timescale

`n_sync` is derived per case so every row realises the same ~0.4 window. Holding it fixed while
`t_max` varies would compare different discretisations and, here, different criteria -- at
`t_max = 50, n_sync = 32` the window is **1.5625**, 3.9x the reference's.

| case | n_sync | dt_sync | R p50 | w(k=1) | t_close p50 | nonfin | d_min p50 | bnd/n_sync |
|---|---|---|---|---|---|---|---|---|
| `near-field` | 33 | 0.3939 | 2.2361 | 0.3939 | 1.44e-3 | 0.0000 | 8.56e-3 | 1.0000 |
| `deep interior` | 33 | 0.3939 | 1.3694 | 0.3939 | 7.86e-5 | 0.0052 | 1.23e-3 | **0.3030** |
| `preset_plambda` | 33 | 0.3939 | 1.0000 | 0.3939 | 6.82e-4 | 0.0000 | 2.27e-3 | 1.0000 |
| `config_basin` | 125 | 0.4000 | 1.0000 | 0.4000 | 2.28e-2 | 0.0000 | 2.36e-2 | 1.0000 |
| `config_stability` | 125 | 0.4000 | 1.0000 | 0.4000 | 4.16e-3 | 0.0017 | 7.60e-3 | 1.0000 |

**THE WINDOW CANNOT RESOLVE INNER-BINARY PHASE ANYWHERE.** `t_close` -- the closest-approach
timescale `2 pi sqrt(d_min^3/M)`, a proxy for the shortest inner period -- runs **17x to 274x below
the window** in every region. A two-end chord cannot tell a full revolution from stationarity
(`tests/outcome_encoding.rs::a_full_revolution_aliases_to_zero_closure` holds that as a property,
not a bug), so the closure arm is structurally blind to a tight bound pair and **rejecting one rests
entirely on the energy arm**. That is not a defect of this port -- the reference uses the same 0.4 --
but it is the reason the two arms are not interchangeable and neither is redundant.

`bnd/n_sync = 0.3030` says the median `deep interior` trajectory completes **30% of its boundaries**
before collision stops it, so the criterion sees a third of that region's run.

### 22.2 The gap does not reproduce, and a wider window makes it worse

Closure at the final boundary of a run with **nothing terminal**, split by the geometric ground
truth. Collided trajectories are excluded rather than folded in: with `stop_on_event` on, a collided
run freezes, its shape stops changing, closure reads **exactly 0**, and it lands in the *bound*
population and destroys the gap it is supposed to measure.

| case | xT | k=1 | k=2 | k=3 | k=4 |
|---|---|---|---|---|---|
| `near-field` | 1 | 0.9 | 0.9 | 0.9 | 1.1 |
| `near-field` | 2 | 1.0 | 0.9 | 0.9 | 1.1 |
| `deep interior` | 1 | 1.5 | 1.4 | 1.2 | 1.1 |
| `deep interior` | 2 | **6.1** | **6.8** | 4.9 | 4.0 |
| `config_stability` | 1 | 0.6 | 0.5 | 0.5 | 0.4 |
| `config_stability` | 2 | 1.9 | 4.0 | 3.8 | 1.7 |

(`sep` = p50(bound)/p50(escaper); `preset_plambda` has **zero** bound trajectories to compare
against and `config_basin` has zero escapers, so neither yields a gap at all.)

**The reference's 383x is not reachable here. The best cell is 6.8x, and in `near-field` there is no
separation whatsoever** -- 0.9 to 1.1 at both horizons, which is the two populations sitting on top
of one another. Three things account for it and they are separable:

- **Maturity.** `|dn/dt| ~ 1/t^3`, so the gap is a function of how far the run has got. The
  reference quotes 383x at `t = 25-30`; this project ships at **13**. Every region that separates at
  all separates *better* at `xT = 2` -- `deep interior` 1.5 -> 6.1, `config_stability` 0.6 -> 1.9.
- **Population.** The ground truth here is geometric and counts a *triple dispersal* as an escape.
  The limit argument assumes a **hierarchical** escape -- binary bounded, `lambda` linear, so
  `alpha -> pi/2`. In a full dispersal nothing converges to a pole and closure stays large. The
  reference's ground truth was "unbound and receding at `t = 30`", which shares a term with the
  criterion; this one shares none, and the disagreement is informative rather than a discrepancy.
- **`near-field` has not escaped yet.** The standing result -- *the escape arm contributes nothing
  at `t = 13`* -- reappears exactly: **0 fires of 576**, against a ground truth that says 42% escape
  by `t = 39`. The criterion is not failing there; the population has not formed.

**`k` narrows the gap rather than widening it**, which is the opposite of what a longer window was
expected to buy. `k = 1` is the best or joint-best rung in every region that separates.

**A closure of exactly zero is not a settled trajectory**, and it is counted rather than absorbed
into a percentile: `deep interior` **10 of 63** bound, `config_stability` **4 of 248**, and **0 of
every escaper population**. Those zeros are what drive `bnd p1` to 0 and make the geometric-midpoint
`tau*` construction undefined in two regions. They are bitwise-identical consecutive shape vectors
on trajectories that have stopped moving numerically; the energy arm is what stops them firing, and
if one ever coincided with a positive relative energy it would read as maximally settled.

### 22.3 CHECK 1 PASSES OUTRIGHT -- 1.0000 against the old criterion's 0 of 895

Candidacy at every boundary of **one** unstopped run at **one** step size, out to `3 t_max` with
`n_sync` scaled so `dt_sync` is unchanged. The question is whether the fired body is still
**unbound** -- the energy arm alone.

| case | fired | u+1 | u+2 | u+4 | u+8 | u end | cand end |
|---|---|---|---|---|---|---|---|
| `deep interior` | 82 | **1.0000** | 1.0000 | 1.0000 | 1.0000 | **1.0000** | 0.9146 |
| `preset_plambda` | 314 | **1.0000** | 1.0000 | 1.0000 | 0.9968 | **1.0000** | 0.4777 |
| `config_stability` | 52 | **1.0000** | 0.9808 | 0.9808 | 0.9608 | **1.0000** | 0.6731 |

**Nothing re-binds.** The old criterion's number on the same question was **0 of 895**. This is the
check that killed the previous criterion and the strongest result in the run.

**And the two columns are not the same question.** `cand end` -- full criterion candidacy at the last
boundary -- reads 0.4777 on `preset_plambda` while `u end` reads 1.0000. Closure is a difference of
neighbouring samples, so it jitters above `tau` on a perfectly settled escape; reading persistence
off *candidacy* would have scored ordinary jitter as a re-binding and reported a correct criterion as
broken. Both are printed so the difference is visible rather than asserted.

**Precision must be read against `u end`, not alone.** `deep interior` reads precision 0.3214 and
`config_stability` 0.5385 -- but every one of those fires is still unbound at `3 t_max`. The ground
truth demands **3x separation growth by `3 t_max`**, so a slow genuine escape fails to be certified.
A precision shortfall with `u end` at 1.0000 is the **ground truth missing them**, not the criterion
inventing them.

Recall is low -- 0.4466 on `preset_plambda`, 0.0217 on `deep interior`, **0.0000** on `near-field` --
and the median firing time is **11.8 of 13**. The criterion is a late one by construction, and at
this horizon most of the escape population has not settled. That is the *right* failure direction:
firing late writes a late timestamp, firing early writes a wrong one permanently.

### 22.4 The `t_end` replay refinement is measurably decoration

When the conjunction first holds at boundary `k`, the interval `(t_{k-1}, t_k]` is replayed from the
saved entry state with the same stepper and the same `dtau`, and the first sub-step at which the
firing body is unbound is taken as `t_end`. Measured:

| case | escapes | at entry | on boundary | distinct t_end |
|---|---|---|---|---|
| `deep interior` | 28 | **1.0000** | 1.0000 | 16 |
| `preset_plambda` | 226 | **1.0000** | 1.0000 | 9 |
| `config_stability` | 26 | **1.0000** | 1.0000 | 14 |

**`at entry` is 1.0000 everywhere: the energy arm always already holds when closure settles, so there
is never a crossing inside the interval to find.** The replay runs and changes nothing. `t_end` under
this criterion is irreducibly quantised to the boundary cadence, and that is a property of the
criterion rather than a bug to route around -- closure is *defined* from the boundary series and has
no finer resolution. The mechanism is visible in §22.3: energy goes positive early and flickers,
closure settles at `t ~ 11.8`.

Kept, with the counter, because the counter is what says it is inert. `refine_escape_time` returns
the **boundary** time in that case, not the entry time -- reporting the entry time would claim an
escape at a playhead the criterion had not yet concluded one at.

### 22.5 The toggle moves 37% of pixels, and the median says zero

| case | stop | escape | collide | frozen | d med | d max | moved |
|---|---|---|---|---|---|---|---|
| `near-field` | false / true | 0.0000 | 0.0226 | 0.0226 / 0.0226 | 0.000e0 | 0.000e0 | 0 |
| `deep interior` | false / true | 0.1042 | 0.6372 | 0.6424 / 0.6875 | 0.000e0 | 6.353e-2 | 18 |
| `preset_plambda` | false / true | 0.4097 | 0.4306 | 0.4306 / **0.8021** | 0.000e0 | **7.365e-2** | **214** |
| `config_stability` | false / true | 0.0451 | 0.4323 | 0.4340 / 0.4792 | 0.000e0 | 3.059e-2 | 22 |

**THE MEDIAN IS EXACTLY ZERO AND 214 OF 576 PIXELS MOVED.** *Never conclude "no effect" from an
aggregate without the per-pixel distribution* -- the rule that has now caught this three times, and
it caught it again here on the run that was expected to confirm the prediction. `preset_plambda`'s
frozen fraction goes 0.4306 -> 0.8021, a **37.15%** increase, and **214/576 = 37.15%** of pixels
move: exactly the pixels that newly froze, by up to **7.4e-2** of chord on a sphere of diameter 2.

### 22.6 AT 1024^2 THE PREDICTION FAILS OUTRIGHT, AND THE RENDER IS THE EVIDENCE

`examples/closure_render.rs`, `preset_plambda`, 1024^2, one sample per pixel, `tau = 1e-3`,
`closure_k = 1`, `n_sync = 33`.

```text
preset_plambda stop=false  826.4s  escape 0.4128 collision 0.4345  frozen 0.8088
preset_plambda stop=true   572.2s  escape 0.4128 collision 0.4345  frozen 0.8088
                shape d: median 0.000e0  max 5.993e-1  pixels moved 392466
```

**392,466 of 1,048,576 pixels move -- 37.43% -- and the worst moves 5.993e-1**, a third of the
shape sphere's diameter. At `n = 24` the worst was 7.4e-2; **at full resolution it is eight times
that**, so the small grid understated the effect by nearly an order of magnitude. The median is
`0.000e0` in both.

**The two `_outcome.png` files are bitwise identical** (same md5). The toggle changes no label and no
event time -- events are recorded whether or not they terminate -- so the only difference is *when
`shape_vec` is read*, and the comparison is clean by construction.

**The images settle it.** Under `stop = false` the ribbons run continuously to every frame edge:
no domes, no tents, no seams. Under `stop = true`, with the same criterion and the same physics,
**two large smooth arcs sweep up from the lower corners and cut straight through the ribbon
structure**, with sharply truncated regions in both upper corners. That is the artefact family,
regrown by the single-variable change, under the *new* criterion.

So the prediction on record -- *"under the new criterion freezing barely matters, and the toggle
produces near-identical images"* -- **does not hold**, and neither does the mechanism stated with it.
The criterion fires at a median `t = 11.8` of 13 with persistence 1.0000, and the shape still moves
by up to 0.6 in the remaining 1.2 time units on more than a third of the frame.

**CLOSURE OVER A 0.4 WINDOW IS A LOCAL RATE TEST; IT DOES NOT CERTIFY STATIONARITY OVER THE REST OF
THE RUN.** `|dn/dt| ~ 1/t^3` bounds the *rate*, and a small rate integrated over 1.2 time units is
not a small displacement -- least of all on the pixels where the run has barely converged. The two
arms are doing different jobs than the design assumed: the criterion is right about **what escaped**
(check 1 at 1.0000 against 0 of 895) and silent about **whether the displayed quantity has settled**.

`stop_on_escape` stays **off**, which is the shipped continuous-ribbon image. Nothing is stacked on
top of this: the next candidate is the one already named -- not terminating on escape at all -- and
it is a decision to be taken on these numbers rather than a third fix.

### 22.7 THE BASIN CONTROL IS DEAD UNDER THIS RULE, AND ITS AGREEMENT PROVES NOTHING

`config_basin` at 1024^2, `t = 50`, `n_sync = 125`, `EscapeRule::Closure`:

```text
config_basin stop=false  4457.1s  escape 0.0000 collision 0.0000  frozen 0.0000
                                  t_end distinct 1  on bdry 100.00%
```

**Nothing terminates. Not one pixel of 1,048,576.** `t_end` takes a single value -- the horizon --
and `frozen = 0.0000`, so the toggle has nothing to act on and the two rows are identical *by
construction*. It reproduces the `n = 24` measurement exactly (0 escapers, 0 collisions, 576/576
bounded), so it is the region and not the resolution.

*A difference can be small because both sides are right or because both are dead.* This one is dead:
an agreement here would be the arithmetic of a field with one value in it. Asserted rather than
assumed -- `d_min p50` is `2.36e-2` against `r_coll = 2.0e-2`, so the collision arm is *just* out of
reach, and `t_close p50` is `2.28e-2` against a 0.4 window, the tightest window-to-inner-timescale
ratio in the set.

**And the control was pointed at the wrong rule.** The saved config is a **basin-mode** slice with
`rEsc = 5`, and basin mode in the reference colours by the terminal outcome under
`dist > r_esc && outward > 0 && E_out > 0`. Rendering it under `Closure` compares a different
criterion, so a disagreement with the reference here would not localise to the physics -- it would
just be the two rules differing, which is the thing being measured everywhere else in §22. **The
basin control has to be rendered under `EscapeRule::Distance(5.0)`**, the rule its own saved config
names, and that is a separate run rather than a conclusion drawn from this one.

### 22.8 What is settled and what is not

**Settled:** check 1 passes at 1.0000, on `deep interior` (check 3), against an independent
geometric ground truth (check 2). The criterion does not latch transients. The reference path is
untouched by construction -- `integrate_az_lc` hardcodes `EscapeRule::Reference` and the cross-check
reads 4/4 PASS.

**Not settled:** the 383x separation does not reproduce at this project's horizon, on this project's
regions, against a ground truth that shares no term with the criterion -- the best cell is 6.8x and
`near-field` shows none. Whether that is maturity, population, or the criterion is separated in
§22.2 but not decided.

**Refuted:** that this criterion makes freezing safe. §22.6 regrows the domes at 1024^2 with 37.43%
of pixels moving and a worst chord of 0.5993. `stop_on_escape` stays **off**.


---

## 23. The `dtau` step control: a blow-up after a boundary-coincident encounter

**Third bug in this sequence found by taking an odd-looking image seriously**, after the LC branch
cut and the escape-termination patchwork. The user saw clusters of magenta (non-finite) pixels on
`config_stability_stop0_*.png`, each ringed by a halo of speckle — red interrupting solid green
bands in outcome mode, salt-and-pepper over the ribbons in continuous mode — where the reference
GLSL renders the same slice solid and continuous.

### 23.1 The mechanism

`dt = A*B*dtau`, and `driver.rs` sized `dtau` **once** at each sync interval's entry:

```rust
let dtau = eta * dt_left / (a0 * b0);      // ONCE
loop { s = rk4::step(&sys, &s, e, dtau); } // never updated
```

`reference/tb_az.py:184` carried the identical line and the identical comment, so **every NumPy
number in the corpus inherited it**.

The physical step is `eta*dt_left` only while `A*B` stays near its entry value. A trajectory sitting
at a close encounter *at a sync boundary* has a tiny `A0*B0`, so `dtau` is enormous; as the bodies
separate through the interval, `A*B` grows by orders and `dt` grows with it. **Giant physical steps
immediately after an encounter.** The GLSL is clean because its step is fixed and cannot explode.

That is not "close encounters" — it is encounters **coinciding with a boundary**, a thin set, which
is why the damage clusters spatially rather than tracking `d_min`.

The comment above the line records how it got here: shrinking `dtau` with separation was a real
earlier bug (drove `dt -> 1e-13`, produced a false "intractable region"), and the correction removed
adaptivity **entirely**. Fixed-per-interval is only right if `A*B` is roughly constant across the
interval.

### 23.2 The obvious repair is Zeno by arithmetic, and the measurement says so

Putting the **remaining** time in the numerator — `dtau = eta*(dt_left - s.t)/(A*B)` — gives
`dt ~ eta*rem`, so `rem_{n+1} = rem_n (1 - eta)`. The interval is approached **geometrically and
never completed**. Measured, `DtauMode::PerStepRemaining` completes

| case | `t/t_max` | steps p50 | budget exhausted |
|---|---|---|---|
| config chart @ t6 | **0.0833** | 30000 | 2304 / 2304 |
| near-field | **0.0303** | 30000 | 2304 / 2304 |
| deep interior | **0.0303** | 30000 | 2304 / 2304 |
| config_stability | **0.0080** | 30000 | 2304 / 2304 |

**And its drift is the best in the whole table** — `1.3e-14` to `6.0e-13`, ten orders below
everything else — because it went nowhere. *A difference can be small because both sides are right
or because one side is dead*, arriving inside the diagnostic written to catch it. The `t/t_max`
column was added for exactly this and is printed first.

The form that keeps the intent holds `dt_left` fixed and recomputes only `A*B`:

```rust
let dtau_now = (eta * dt_left / ab).min(dtau_entry);
```

**The cap is one-sided in the right direction.** When `A*B` grows the recomputed value falls and the
blow-up is removed; when `A*B` *falls* at a close approach the cap holds `dtau` at nominal, so
`dt = A*B*dtau` shrinks with the separation — which is what regularisation buys and what the
original comment wanted. It removes the over-correction without reintroducing what was
over-corrected for.

### 23.3 The cost, which gates everything — about 10% more steps

`examples/dtau_step.rs`, 48x48 nominal trajectories per case, `max_steps = 30000`, nothing terminal.

| case | mode | steps p50 | steps p99 | budget |
|---|---|---|---|---|
| config chart @ t6 | fixed | 1464 | 2066 | 0 |
| | per-step-int | 1625 | 2260 | 0 |
| near-field | fixed | 3961 | 4109 | 0 |
| | per-step-int | 4311 | 4429 | 0 |
| deep interior | fixed | 4725 | 7034 | 3 |
| | per-step-int | 5427 | 7314 | **1** |
| config_stability | fixed | 13995 | 18914 | 0 |
| | per-step-int | 15308 | 21027 | 1 |

`t/t_max = 1.0000` on every row. The tail does **not** blow up and the budget count does not rise —
`deep interior` falls 3 to 1. **The fix does not swap non-finite pixels for budget-exhausted ones**,
which was the failure it had to be cleared of before any drift number could be read.

### 23.4 The NumPy numbers reproduce, and `med(affected)` falls 123x-2400x

Row one is the reproduction case at the settings the mechanism was found at (`t = 6`,
`n_sync = 12`, `eta = 0.01`).

| case | mode | nonfin | drift p50 | frac>1e-6 | med(affected) |
|---|---|---|---|---|---|
| config chart @ t6 | fixed | **11** | 1.721e-7 | **0.3555** | 2.818e-4 |
| | per-step-int | **2** | 2.248e-8 | **0.1888** | **1.192e-6** |
| near-field | fixed | 0 | 1.515e-9 | 0.0113 | 4.315e-6 |
| | per-step-int | 0 | 1.233e-9 | 0.0082 | **1.797e-9** |
| deep interior | fixed | **42** | 2.094e-1 | 0.9206 | 2.496e-1 |
| | per-step-int | **15** | 1.613e-4 | 0.7799 | **1.957e-4** |
| config_stability | fixed | **7** | 2.549e-6 | 0.5621 | 1.467e-4 |
| | per-step-int | **5** | 8.083e-8 | 0.3056 | **1.190e-6** |

The NumPy measurement gave non-finite 8 -> 2 and `frac>1e-6` 35.6% -> 18.6% on the same case; this
reads 11 -> 2 and 35.55% -> 18.88%. It reproduces.

**`med(affected)` is over the pixels hot under EITHER mode, the same set in every row.** The first
cut conditioned each row on its own hot set and read the improvement *backwards* (3.53e-4 ->
4.44e-4), because the selection moves with the thing being measured. Paired, it falls 236x, 2400x,
1275x and 123x.

### 23.5 The spatial test — the counts fall, and the clustering ratio RISES

`cluster` is the fraction of hot pixels having a hot 4-neighbour, divided by the same under a
random field of the same density — `1 - (1-base)^4`. **1.0 is chance.**

| case | mode | cluster | n hot | nf w/ hot nb | n nonfin |
|---|---|---|---|---|---|
| config chart @ t6 | fixed | 1.164 | 819 | 1.0000 | 11 |
| | per-step-int | **1.650** | **435** | 1.0000 | **2** |
| near-field | fixed | 12.133 | 26 | — | 0 |
| | per-step-int | **16.154** | **19** | — | 0 |
| deep interior | fixed | 1.000 | 2121 | 1.0000 | 42 |
| | per-step-int | **0.982** | **1797** | 1.0000 | **15** |
| config_stability | fixed | 1.000 | 1295 | 1.0000 | 7 |
| | per-step-int | **1.164** | **704** | 1.0000 | **5** |

**The prediction was that the clustering ratio falls toward chance. It rises, in three regions of
four.** Stated rather than dressed up. The counts fall — `n hot` by 47%, 27%, 15% and 46%,
non-finite by 82%, 64% and 29% — so what the fix removes is the *diffuse* high-drift population,
and what survives is the genuinely clustered core, which is a higher ratio on a smaller set.

**Read `n hot` as the measurement and the ratio as a shape statistic.** Chance depends on the
density, so the ratio is *not comparable across the two rows of a pair* — it rises partly because
the fix thins the mask. And in `deep interior` the mask **saturates** at 92% hot, where chance is
~1 and the ratio can say nothing at all; that is the standing regional mask-saturation result
landing on a third statistic. The ratio column is quoted here because it was the shape the bug was
found by, not because it decides anything.

**`nf w/ hot nb` is 1.0000 in every cell, before and after — which means it cannot fail here.**
Every non-finite pixel sits inside a high-drift neighbourhood under both modes, so the statistic
the user quoted as 6/6 is saturated in this configuration and is a description rather than a test.
It says the magenta that remains is magenta for the same reason it always was; it could not have
said otherwise.

### 23.6 A few trajectories get much worse, and they are near-singular

Ranked by the **absolute** rise, not the ratio. The first cut ranked by ratio and returned pixels
whose `fixed` drift was accidentally excellent — `9.2e-10 -> 4.5e-4` is a ratio of `4.9e5` on a pixel
that is still the region's best, while near-field's `drift max` *improved* 36x overall. A ratio
ranking finds a small denominator, not a bad outcome.

| case | pixel | fixed | per-step-int | `d_min` |
|---|---|---|---|---|
| config_stability | 1616 | 8.263e-6 | **3.617e183** | **2.102e-99** |
| deep interior | 151 | 1.642e-7 | **2.178e19** | 3.672e-9 |
| config chart @ t6 | 1439 | 4.637e1 | 1.640e7 | 2.069e-3 |

270 / 189 / 449 / 392 of ~2300 finite-both pixels get worse per case. The worst are at `d_min` of
`1e-99` and `1e-9` — genuine near-singular states, where finer stepping penetrates *further* toward
the singularity within the interval instead of blowing past it. Both modes are wrong there; the
fix's failure is louder. `error_ratio` and `n_nonfinite` are what flag those, and they are never
discarded.

### 23.7 `T::TINY` underflows at **f64**, not only f32

`ab_floored` fires on exactly **one** trajectory of 2304, in `config_stability` under
`per-step-int`, where the raw `A*B` reaches **exactly `0.000e0`**. The standing note has this as an
f32 property — `TINY*TINY` underflowing so a doubly-degenerate state gives `dtau = inf`. At f64
`TINY = 1e-300` and `TINY*TINY = 1e-600` **also underflows to zero**, so the same hole is open at
both precisions. Here the `.min(dtau_entry)` cap absorbs it (`inf.min(x) == x`); under
`FixedPerInterval` there is no cap and nothing to absorb it. Minimum raw `A*B` seen: `2.663e-217`
(config chart), `4.206e-3` (near-field), `3.295e-219` (deep interior), `1.114e-217`
(config_stability).

### 23.8 The 109 escaping near-field pixels at `t = 20` were the bug

This is the largest correction to the corpus in this section. The standing finding reads *"zero of
1024 near-field pixels fire at `t_max = 13`, 109 at `t_max = 20`"*. Under the corrected stepping,
**zero of 1024 fire at `t = 20`**.

That is not a chaotic reshuffle, and the discriminator is on the trajectories themselves:

| population (`fixed`, `t = 20`) | drift p50 | drift p90 | `d_min` p50 |
|---|---|---|---|
| the 109 that FIRE | **1.147** | 2.186e2 | 8.883e-3 |
| the 915 that stay silent | **6.233e-5** | 1.941e-2 | 6.876e-4 |
| the same 109, under `per-step-int` | **1.610e-3** | 2.182e-3 | 3.194e-4 |

**A median energy drift of 1.147 is 115% of the total energy.** Those trajectories are not physics.
A giant post-encounter step throws a body outward, it reads as unbound and receding at the next
boundary, and the arm latches. Under correct stepping the same pixels sit four orders lower and do
not escape.

**Genuine escape is not suppressed.** At `t = 40` the two modes give **280 and 308** of 1024, and at
`t = 30`, 2 and 2. Burrau's escape is simply later than `t = 20`; the earlier figure was measuring
the step control.

### 23.9 Two tests went vacuous, and both control arms caught it

Neither is a regression in the fix; both are guards firing because the thing under test stopped
being exercised.

- `escape_matches_the_legacy_classifier` asserted `fired > 0`, and at `t = 20` nothing fires any
  more. Re-pinned at `t = 40` (7 of 25 fire, worst drift 7.5e-6), and **`n_sync` is now scaled with
  `t_max`** — the old form ran `t = 13` at a 0.406 interval and `t = 20` at 0.625, which is this
  project's own discretisation trap inside a test.
- `the_two_to_one_constraint_holds_and_the_control_violates_it` asserted its unbalanced control
  actually violates 2:1. **The fixture region has now moved twice.** It was `deep interior`, moved
  to `near-field` when the escape distance gate flattened `deep interior`, and moves back now: under
  the corrected stepping near-field is gap 1 at **all twenty-four** cells of
  `alpha_hi x tau x n` swept, while `deep interior` recovers gap 2 at `n = 4`.

### 23.10 The cross-check, run three ways

| pair | result |
|---|---|
| new Rust vs new Python (`per-step-interval`) | **4/4 PASS** |
| old Rust vs old Python (`fixed`) | **4/4 PASS** |
| old Rust vs new Python (deliberate mismatch) | **refused** — header assertion |

Either of the first two alone would pass while the two transcriptions diverged, which is why both
are run. The mode is emitted into the TSV header and `compare.py` asserts headers match, so a
mismatched pair is a hard failure rather than a quiet disagreement; the third row confirms that
guard fires.

`fixed` reproduces the committed pre-fix `tb_az.py` **bitwise** on the smoke test —
`3.1966735584829495e-09` — so the mode switch is faithful. Under the fix the median is
`4.462793760861922e-10`. **Separately: the number `reference/README.md` used to quote,
`3.892633125701676e-09`, does not reproduce on the unmodified committed reference either.** It was
already wrong, and only running it found out.

BRIEF §5 acceptance is unchanged: `M = 12`, `R = 2.23606797749979`, `E = -12.816666666666666`; gauge
invariance `rel err = 0e0` at `alpha in {0.25, 1, 4}`; two-body radial collision
`d_min = 1.288146e-11` with `|dE/E| = 2.947299e-14`; `error_ratio` median `1.000000`. Suite: **223
passed, 0 failed**.

### 23.11 The modes have to be shown to differ, and where they must not

`tests/dtau_step_control.rs`. `dtau_mode` is threaded through four structs and a closure, and a
field that never reaches the stepper produces perfectly plausible numbers — so each test names the
configuration in which it fires and asserts a direction.

| test | what would have to be true for it to fail |
|---|---|
| `the_modes_disagree_where_ab_grows_across_an_interval` | `dtau_mode` not reaching `rk4::step`; asserts >half of 64 `deep interior` trajectories differ, and by more than round-off |
| `the_modes_agree_to_roundoff_where_ab_is_flat` | the difference above being two arbitrary steppers rather than the mechanism; over `t = 0.5` with no encounter the cap binds throughout and the modes agree to `< 1e-8` |
| `per_step_remaining_stalls_rather_than_integrating` | Zeno not happening; asserts the budget is exhausted, `t` stays inside the first interval, **and that the stalled mode's drift is smaller than the completed one's** — the trap, asserted rather than described |
| `ab_min_is_recorded_and_the_f64_floor_never_binds` | `ab_min` written from the entry state rather than per step; asserts it falls below `1e-3` on `deep interior` and that the floor binds nowhere there |

### 23.12 The renders — the prediction holds on all four clauses

`examples/dtau_render.rs`, 1024², one sample per pixel, `results/dtau_fix/`. Only `dtau_mode`
differs between the two rows; every other setting is `closure_render`'s, which is what the
committed "before" images used.

| case | mode | nonfin | simfail | hot | drift ramp p2–p98 | secs |
|---|---|---|---|---|---|---|
| config_stability | fixoff | **30109** | **0** | 0.9285 | (1.079e-8, 4.071e7) | 725.9 |
| config_stability | fixon | **178** | **0** | 0.8558 | (1.075e-8, 4.364e7) | 767.9 |
| preset_shape | fixon | 27 | 0 | 0.1767 | (5.664e-10, 2.105e3) | 318.2 |
| preset_prho | fixon | 1 | 0 | 0.0116 | (1.243e-13, 2.378e-8) | 760.3 |
| preset_plambda | fixon | 0 | 0 | 0.0033 | (3.943e-12, 9.306e-9) | 557.7 |
| preset_shape_pl | fixon | 1 | 0 | 0.0451 | (1.471e-11, 8.278e-5) | 447.1 |

**The magenta count falls 30109 → 178, a factor of 169** — from 2.87% of the frame to 0.017%. And
`simfail` is **0 on both rows**, which is the clause that had to be checked separately: a fix that
removed non-finite pixels by exhausting the budget instead would have moved the first column and
fixed nothing.

Against the prediction, clause by clause:

- **magenta clusters shrink sharply or vanish** — yes, 169x, and visibly nothing left in either
  `_outcome` or `_drift`.
- **the speckle halos go, and the outcome bands become solid** — yes. The blue *bounded* regions
  that were shot through with red-and-magenta stipple are now large and solid, and the green bands
  are unbroken. 348,314 of 1,048,576 outcome labels flip, concentrated exactly where the halos were.
- **the continuous ribbons stop being interrupted** — yes, and the filament structure is *more*
  legible afterwards, not less.
- **nothing new appears** — the large-scale geometry is the same slice, recognisably feature for
  feature.

**What remains is genuine.** The red regions still carry blue and green stipple, but that is
fractal mixing at the sampling limit, not the halo pattern: it has no magenta in it, it does not
ring a core, and it survives at every resolution. Distinguishing the two is what the `_drift` panel
is for, and there the "before" is grained with magenta through every arc while the "after" is
smooth.

`shape d` between the rows is median `1.113e-1`, max `2.000` — the full diameter of the shape
sphere — over 909,184 pixels moved. That is **not** a defect and is not evidence either way: the
two modes integrate genuinely different trajectories through a chaotic region, so they must
diverge. The question the images answer is whether the structure got cleaner, not whether the
pixels agree.

### 23.13 The diagnostic render mode

`Scalar::Drift` and `colour::drift_rgb` — `energy_drift_max` on an inferno ramp, `DEBUG_NAN` magenta
for the same veto set `colour::rgb` applies, auto-ranged over the field's own p2-p98 via
`colour::range_q`. Univariate on purpose: `rgb` is bivariate and a diagnostic asked under it
inherits the shape field's structure.

**When a numerical defect is suspected, render the diagnostic field, not the science field.** The
science fields show a defect only after it has propagated into a spread or a label; the drift map
shows it at source, as coherent arcs with the non-finite pixels inside them. `_drift.png` is written
for **both** modes — the "before" map is the artefact worth keeping, because it is what the
signature looks like.


---

## 24. The boundary overshoot: the `dtau` fix shipped without its partner

`§23` fixed how `dtau` is *sized* and left untouched how the interval *ends*. The two are not
separable, and shipping the first alone made the images worse while the drift numbers improved.

### 24.1 The mechanism, and why it interacts

The march exits a sync interval by **overshooting** it — the loop condition was `s.t >= dt_left`
— and only the *clock* was corrected:

```rust
// The overshoot past the boundary is clipped in the time bookkeeping only; the
// state written back is the overshot one. Sub-step interpolation is not done.
t += s.t.min(dt_left);
```

The Cartesian state handed to the next interval is the overshot one. That is a **first-order**
error injected at every one of `n_sync` boundaries, inside an RK4 march.

Under `FixedPerInterval` `dtau` is constant across the interval, so the overshoot is a fixed slice
of fictitious time: every trajectory in a neighbourhood overshoots by roughly the same amount. The
error is large but spatially **smooth** — it displaces the picture without breaking it.

Under `PerStepInterval` the final step's size is a function of the local `A*B`, so the overshoot
becomes **a function of local state** and neighbouring pixels overshoot by different amounts. A
spatially-varying error injected at every boundary is measurably worse than a smooth one, and
§24.4 and §24.6 measure exactly that.

**The nested-arc appearance that prompted this was NOT produced by it.** §24.8 renders all four
arms and the arcs are present in every one, including the arm that predates both changes. The
defect is real and independently measured; it is not the cause of that picture. The rest of this
section reports the defect. Read §24.8 before attributing any appearance to it.

### 24.2 The four arms

```text
  A  dtau fixed      + overshoot present   the original committed behaviour
  B  dtau per-step   + overshoot present   what §23 shipped
  C  dtau fixed      + overshoot clamped
  D  dtau per-step   + overshoot clamped   the default from here
```

Two knobs, four cells. Rendering only the diagonal would show a difference and say nothing about
which knob produced it, and the claim being made is specifically about the cross terms.

### 24.3 The fix

`AzOpts::clamp_final_step`, default **on**, applied after `dtau_for_step` returns so it composes
with every `DtauMode` rather than being a fourth mode:

```rust
let dtau = if clamp_final {
    dtau.min((dt_left - s.t).max(T::zero()) / ab)
} else { dtau };
```

`ab` is **the same floored product the mode already computed**. Recomputing it would let the clamp
and the step disagree.

**The landing tolerance is RELATIVE to `dt_left`, and an absolute one is a bug the suite caught.**
All times rescale by `alpha^{3/2}` under the project's scale gauge, so an absolute slack is a
different tolerance at every scale: at `alpha = 0.25` the same `1e-15` is eight times wider in
relative terms and a rescaled twin can land one step earlier. `gauge_invariance::
shape_spread_is_invariant_under_the_scale_symmetry` asserts **bitwise** equality and fired at
`4.24e-15`. `Real::LAND_EPS_REL` is `1e-14` at f64 and `1e-5` at f32.

### 24.4 The convergence order — the cheapest confirmation, and it is the ORDER that matters

Chenciner–Montgomery figure-eight, equal masses, integrated over exactly one period,
`n_sync = 32` fixed so the number of boundaries is the same at every `eta`. `closure` is the max
component difference between the state at `T` and the state at `0` — a pure error measure with no
reference trajectory and no chaotic amplification.

| arm | `eta = 0.02` | `0.01` | `0.005` | `0.002` | `0.001` | **order** |
|---|---|---|---|---|---|---|
| A `fixed` + overshoot | 8.665e-2 | 3.473e-2 | 1.975e-2 | 7.884e-3 | 2.976e-3 | **1.13** |
| B `per-step` + overshoot | 4.681e-2 | 2.182e-2 | 9.191e-3 | 4.082e-3 | 1.950e-3 | **1.06** |
| C `fixed` + clamp | 3.420e-5 | 6.735e-6 | 1.123e-6 | 3.241e-7 | 3.595e-9 | **3.06** |
| D `per-step` + clamp | 7.212e-5 | 1.805e-5 | 4.500e-6 | 6.885e-7 | 1.435e-7 | **2.08** |

**An error falls for many reasons; only the order says the leading term changed.** First-order to
roughly third, and the error at `eta = 0.001` falls **827,000x** (A → C).

**D lands at 2, not at C's 3, and that is stated rather than smoothed.** The clamp sizes the final
step from the *instantaneous* `A*B`, which predicts the time increment to first order, so the
landing residual is `O(h^2)` per boundary — and where that residual overshoots, the tolerance
accepts it rather than paying another step. Under `FixedPerInterval` the entry `dtau` is a better
predictor of the whole interval and the residual is smaller. Two is ample against one; three would
be better and is not what this implementation gives.

**Read the per-rung columns as noise.** Each is a two-point estimate over a factor of two and C's
run 2.34, 2.58, 1.36, 6.49. The endpoint slope across the decade is the measurement.

### 24.5 ENERGY DRIFT IS NEARLY BLIND TO THIS DEFECT

The diagnostic field that found §23's bug **cannot find this one.**

| case | arm | drift p50 | drift p99 | hot | nonfin | budget | steps p50 |
|---|---|---|---|---|---|---|---|
| near-field | A | 1.515e-9 | 2.428e-6 | 0.0113 | 0 | 0 | 3961 |
| near-field | B | 1.233e-9 | 2.427e-7 | 0.0082 | 0 | 0 | 4311 |
| near-field | C | 5.555e-8 | 1.768e-6 | 0.0135 | 0 | 0 | 4064 |
| near-field | D | **1.125e-9** | **2.168e-7** | **0.0061** | 0 | 0 | 4375 |
| deep interior | A | 2.094e-1 | 3.364e6 | 0.9206 | 42 | 3 | 4725 |
| deep interior | B | 1.613e-4 | 6.077e5 | 0.7799 | 15 | 1 | 5427 |
| deep interior | C | 1.824e-1 | 1.783e6 | 0.9193 | 45 | 2 | 4721 |
| deep interior | D | 1.960e-4 | **2.881e4** | 0.7808 | 17 | **0** | 5451 |
| config_stability | A | 2.549e-6 | 5.789e3 | 0.5621 | 7 | 0 | 13995 |
| config_stability | B | 8.083e-8 | 6.233e3 | 0.3056 | 5 | 1 | 15308 |
| config_stability | C | 2.065e-6 | 4.190e3 | 0.5473 | 4 | 0 | 14114 |
| config_stability | D | 8.688e-8 | **2.151e3** | 0.3090 | 7 | 1 | 15392 |

`48x48` nominal trajectories, `eta = 1e-2`, nothing terminal.

**A → C moves the median drift 37x the WRONG WAY in near-field** (1.5e-9 to 5.6e-8) while the same
change buys 24,000x on the figure-eight. The overshoot displaces the state in *time*, and the AZ
energy is nearly stationary along the flow, so a time displacement barely registers in `|dE/E|`.
The NumPy smoke test says the same: median drift 3.197e-9 → 4.047e-9 under the clamp at `fixed`.

**The generalisation, and it is the correction to §23's own lesson.** *Render the diagnostic field,
not the science field* is right, and incomplete: **a diagnostic field is specific to a class of
defect.** Energy drift finds a step that grew too large. It does not find a step that ended in the
wrong place. Ask what the diagnostic would say about the defect you are looking for before reading
it as clean.

What the drift table *does* confirm is the prediction's last clause: **D improves on B**, not merely
on A — `drift p99` 6.077e5 → 2.881e4 in `deep interior` (21x) and 6233 → 2151 in `config_stability`
(2.9x), with `budget` 1 → 0 and 1 → 1. The clamp is not paid for in step budget: `steps p50` rises
0.3–1.5% from B to D.

### 24.6 HOW MUCH OF THE FIELD MOVES — and THE TWO RESOLUTIONS DISAGREE BY 26x

Chord between final shape vectors, on a sphere of diameter 2. **At 1024², the shipping
resolution**, all four arms in one run so the pairs are exact:

| pair | moved | frac | chord p50 | chord max | labels flipped |
|---|---|---|---|---|---|
| A→B | 909184 | 0.8671 | **1.113e-1** | 2.000e0 | 348314 |
| C→D | 914569 | 0.8722 | **4.367e-2** | 2.000e0 | 301906 |
| B→D | 979644 | 0.9343 | 1.242e-2 | 2.000e0 | 145651 |
| A→C | 973512 | 0.9284 | 6.284e-2 | 2.000e0 | 291427 |
| A→D | 973739 | 0.9286 | 1.294e-1 | 2.000e0 | 353310 |

**A→B reproduces §23's render exactly** — 909184 moved, median 1.113e-1, 348314 labels — which is
what says this run and that one are the same measurement.

**THE INTERACTION IS REAL AND FAR SMALLER THAN THE 48x48 GRID SAID.** Under the clamp, switching
`dtau_mode` moves the field **2.5x** less (A→B 1.113e-1 against C→D 4.367e-2). On the 48x48 grid
the same ratio read **66x** in `config_stability` and **316x** in `near-field`:

| case | resolution | A→B chord p50 | C→D chord p50 | ratio |
|---|---|---|---|---|
| `config_stability` | 48x48 | 1.565e-2 | 2.374e-4 | **66x** |
| `config_stability` | **1024²** | 1.113e-1 | 4.367e-2 | **2.5x** |
| `near-field` | 48x48 | 1.854e-2 | 5.876e-5 | 316x |

**Quote the 1024² number.** The standing rule is *read the max and the moved count at the
resolution that ships*, and it applies to the median too — for the opposite reason to last time. In
§23 the coarse grid **understated** the effect eightfold; here it **overstates** the improvement
twenty-six-fold, because 2304 samples over the same window are dominated by the tame majority
while a million samples land squarely in the chaotic population, where any change to the step
control diverges regardless. Both directions are the same defect: a grid coarse enough to miss the
population that matters.

The 48x48 rows in §24.5 stand as they are — drift is a per-trajectory statistic and does not have
this problem — but **no chord ratio should be quoted from them.**

**`moved` ORDERS THE PAIRS BACKWARDS.** B→D moves the *most* pixels (0.9343) and displaces them the
*least* (1.242e-2); A→B moves the fewest (0.8671) and displaces them the most (1.113e-1). It counts
pixels differing in the last bit, which on a chaotic field is a fact about the field and not about
the change — the whole table sits in 0.867–0.934, and its ordering is the reverse of the magnitude
one. The 87% figure §23 reported is that statistic; it is a real and correct number that answers a
different question than the one it was read as answering. **`chord p50` is the discriminator.**

`chord max` is **2.000 — antipodal — in every pair**, which is the same statement one level up:
somewhere on a million-pixel chaotic slice, two step controls put a trajectory on opposite poles of
the shape sphere. That is not evidence about either change.

### 24.7 IS THE SLICE RESOLVED AT ALL? — arm D at `eta`, `eta/2`, `eta/4`

| case | pair | moved | chord p50 | chord max | drift p50 |
|---|---|---|---|---|---|
| near-field | 1.00e-2 → 5.00e-3 | 2304 | 2.116e-5 | 9.982e-4 | 7.001e-11 |
| near-field | 5.00e-3 → 2.50e-3 | 2304 | **5.272e-6** | 1.678e-5 | 4.358e-12 |
| deep interior | 1.00e-2 → 5.00e-3 | 2287 | 4.709e-2 | 1.999e0 | 6.839e-6 |
| deep interior | 5.00e-3 → 2.50e-3 | 2290 | **5.505e-3** | 1.999e0 | 2.924e-7 |
| config_stability | 1.00e-2 → 5.00e-3 | 2299 | 2.119e-5 | 1.996e0 | 4.584e-9 |
| config_stability | 5.00e-3 → 2.50e-3 | 2304 | **4.254e-6** | 1.990e0 | 2.640e-10 |

**`config_stability` IS resolved at the shipping settings, and the ratio says so quantitatively.**
The median displacement falls **4.98x for a 2x step reduction** — consistent with §24.4's order
2.08 — and its absolute level at `eta = 1e-2` is `2e-5` on a diameter-2 sphere, one part in
100,000. Horizon 50 was a live worry; it is not the answer here. `near-field` falls 4.01x, the same
story.

**`deep interior` does not converge in any useful sense and never will.** Its chord median falls
8.6x but from `4.7e-2`, and its max sits at 1.999 — antipodal — at every rung. That is chaotic
divergence over `t = 13`, not a discretisation artefact, and no step size buys it off. The
distinction matters: *a difference can be small because both sides are right or because both are
dead*, and here it is large because the physics is.

### 24.8 THE RENDERS — THE PREDICTION FAILS ON ITS MAIN CLAUSE

`examples/overshoot_render.rs`, 1024², one sample per pixel, `results/overshoot_fix/`. Panels are
`<case>_arm{A,B,C,D}_{uniform,outcome,drift}.png`. `results/dtau_fix/` is the "before" for this
comparison and is untouched.

| arm | escape | collision | **nonfin** | simfail | hot | spread ramp | drift ramp |
|---|---|---|---|---|---|---|---|
| A `fixed` + overshoot | 0.2618 | 0.3632 | **30109** | 0 | 0.9285 | (6.850e-5, 4.955e-1) | (1.079e-8, 4.071e7) |
| B `per-step` + overshoot | 0.2048 | 0.3477 | **178** | 0 | 0.8558 | (6.830e-5, 4.954e-1) | (1.075e-8, 4.364e7) |
| C `fixed` + clamp | 0.2415 | 0.3605 | **2071** | 0 | 0.9257 | (5.519e-5, 4.953e-1) | (1.082e-8, 4.828e7) |
| D `per-step` + clamp | 0.2017 | 0.3477 | **178** | 0 | 0.8541 | (5.705e-5, 4.955e-1) | (1.074e-8, 4.277e7) |

**Arms A and B reproduce `dtau_fix`'s `fixoff` and `fixon` rows exactly** — 30109 and 178
non-finite, ramps identical to four digits. The new flag does not leak into the old paths and the
comparison is clean by construction.

**The prediction, clause by clause:**

| clause | verdict |
|---|---|
| the stacked-crescent banding in the green and blue mid-field of `_uniform` **disappears in D** and is present in B | **FAILS** |
| arm C is smooth too, and probably smoother than A | holds — magenta 30109 → 2071, speckle largely gone |
| the white regions do not grow further, ideally shrinking toward A | holds — extent unchanged; A's pale regions are the *noisiest* |
| magenta stays at arm B's level, near zero | **holds exactly** — 178 = 178 |
| drift improves in D over B, not merely over A | holds — `drift p99` 6.077e5 → 2.881e4 in `deep interior` (§24.5) |

**THE CRESCENTS ARE PRESENT IN ALL FOUR ARMS, INCLUDING THE ONE THAT PREDATES BOTH CHANGES.** The
dark-green wedge left of centre carries the same stacked chevrons in A, B, C and D; the lower-right
blue band carries the same nested arcs in all four. Neither knob causes the banding and neither
removes it. **The mechanism proposed for it in §24.1 — a spatially-varying overshoot whose
accumulated level sets are nested arcs — is not what produces these.** That mechanism is real, and
§24.4 and §24.6 measure it; it is not the cause of this appearance.

**What the arcs are, as far as this measures.** Under outcome-class colouring arm D's arcs **vanish
entirely**, and the region boundaries survive and sharpen into solid blocks with genuinely
fractal-mixed zones between them. That is §21's standing result at a new site: *the banding is a
colouring artefact; the crisp edges are not.* Whether it is the **same** mechanism — closure's
`t_end` being irreducibly quantised to the boundary cadence — is **not tested here and is not
claimed.** What is established is narrower and firmer: the arcs live in the continuous field's
ramp, they predate both knobs, and no step-control change touches them.

**A → D is still the whole case for the pair.** A's outcome panel is shredded — magenta swathes
through the red bands and salt-and-pepper everywhere; D's is solid. And the attribution splits
cleanly: the `dtau` fix removes the magenta (30109 → 178, and 2071 → 178 on top of the clamp), the
clamp removes most of what is left on its own (30109 → 2071), and **neither touches the arcs**.

**The presets at arm D, against `dtau_fix`'s `fixon`:** non-finite `preset_prho` 1 → 0,
`preset_plambda` 0 → 0, `preset_shape_pl` 1 → 0, and **`preset_shape` 27 → 42**. The one that rises
is the chart whose undetermined pixels are already on record as *the chart, not the grid* — triple
collisions, stable across resolution, zero `SimFailed` at every grid size. 42 of 1,048,576 is
0.004%, and it is the instrument reporting rather than a fault. `simfail` is 0 on all four presets
under both arms, so nothing was traded.

*A finding read off a wireframe is a finding about an appearance* — the standing rule — and this is
the same at a render. The appearance was taken seriously, it produced a real and independently
measurable defect (§24.4, first order to third), and it turned out not to be the cause of the
appearance that prompted it. Both halves are worth recording.


---

## 13. Reproducing any of this

**Two of these commands were wrong, and only running them found it.** The `pan_sequence` line said
`9 2000 512` where the committed dumps were made at `9 20000 1024`, and the `slice_gallery` line
said `4000 … 512` against a committed `40000 … 1024`. Regenerating from the documented commands
produced nineteen dumps that did not reproduce: nine at a tenth of the budget, and ten with an
**identical tree and a different `decision` column** — 252 leaves moved from `MaxRelDepth` to
`ScreenFloor` purely by the viewport. Same leaf count, different stop reason, which is this
project's own standing lesson arriving through its documentation.

A sample of eleven dumps had already reported "reproduces bitwise". Checking all sixty-nine is what
caught it. **Verify a regeneration over the whole corpus, not a sample**, and diff the `decision`
column specifically — it is the one that moves when a parameter is wrong while the tree is not.

**And the default under them moved, which is a second way a documented command can be wrong.**
`SchedCfg::default().k_frac` was `1.0` when every artefact below `§18` was made and is now `0.25`.
So **every command in this table that does not name `k_frac` produces a different tree than the one
it reproduced before** — a ranked one rather than the uniform-mode control. That is the intended
direction and it is stated rather than quietly absorbed: where the committed artefact is the
*before*, the command carries `1.0` explicitly (see the `chart_gallery` and `balanced_march` rows).
Everywhere else the table now names the ranked run.

Every table above comes from a committed example. Raw output for all of them is in
[`results/output/`](results/output/), the acceptance-gate and cross-check output is in
[`results/tests/`](results/tests/), and the images and 64×64 raw dumps are in
[`results/`](results/). [`results/README.md`](results/README.md) indexes all of it.

| result | command |
|---|---|
| §1 refinement criterion | `cargo run --release --example refinement_criterion` |
| §2 statistical convergence | `cargo run --release --example convergence` |
| §3 the seven pixels | `cargo run --release --example worst_128` |
| §3 the remedy | `cargo run --release --example refine_pass` |
| §4 latching | `cargo run --release --example latching_decision` |
| §5 f32 | `cargo run --release --example f32_report` |
| §5 branch cut | `cargo run --release --example lc_cut_proximity` |
| §5 deep interior | `cargo run --release --example deep_interior` |
| §5 event class | `cargo run --release --example spread_event_correction` |
| §1 offset schemes | `cargo run --release --example halton_noise_floor` |
| §1 pooled vs true parent | `cargo run --release --example pooled_vs_true_parent` |
| the offset properties | `cargo test --release --test halton_offsets -- --nocapture` |
| a full slice | `cargo run --release --bin prin -- --region near-field --size 256 --out out` |
| §5 the descent | `cargo run --release --bin prinq -- --region near-field --budget 8000 --tau 1e-4 --alpha-hi 0.2 --alpha-band 0.4 --overlay 512 --out t` |
| §5 thresholds | `cargo run --release --example sched_sweep` |
| §5 termination | `cargo run --release --example sched_terminate -- 50000 1e-4 0.2` |
| §5 aggregation | `cargo run --release --example sched_agg -- 6000 1e-4 0.2` |
| §5 policies | `cargo run --release --example sched_policies -- 10000 1e-4 0.2` |
| §5 ordering | `cargo run --release --example sched_order -- 1500 1e-4 0.2` |
| §5 N sweep | `cargo run --release --example sched_n_sweep -- 4000 1e-4` |
| §5 thrash | `cargo run --release --example sched_thrash -- 1e-4 0.2` |
| §8 screen floor, q1/q2 | `cargo run --release --example sched_screen -- 50000 1e-4 0.2` |
| §8 q7 under the veto | `cargo run --release --example sweep_screen -- 2000` |
| §8 the E sweep | `cargo run --release --example e_sweep -- 6000 1e-4 0.2` |
| §8 the estimator bias in E | `cargo run --release --example spread_bias_e -- 48` |
| §8 slice variety | `cargo run --release --example slice_variety -- 4000 1e-4 0.2` |
| §8 the decode ladder | `cargo run --release --example decode_ladder` |
| §8 deep zoom in situ | `cargo run --release --example deep_zoom -- 400` |
| §8 SSAA resolve | `cargo run --release --example ssaa_resolve -- 256` |
| §8 the open items | `cargo run --release --example open_items -- 6000` |
| §8 aggregation vs the floor | `cargo run --release --example agg_vs_floor -- 6000 1e-4 0.2` |
| §8 adaptive render | `cargo run --release --example adaptive_render -- near-field 6000 1e-4 0.2` |
| §8 the zoom ladder (APNG) | `cargo run --release --example zoom_sequence -- near-field 9 2000` |
| §8 the vertical-slice tests | `cargo test --release --test vertical_slice -- --nocapture` |
| the acceptance gates | `cargo test --release -- --nocapture` |
| the NumPy cross-check | `cargo test --release --test xcheck -- --ignored --nocapture` |
| the horizon table | `python3 tools/xcheck/horizon.py [--lc-unstable]` |
| §10.1 between vs within | `cargo run --release --example between_vs_within -- 2000 1e-4 0.2 512` |
| §10.2-10.4 the metric | `cargo run --release --example criterion_metric -- 6 8 1e-4 13` (a fifth argument sets the output root; it defaults to `results` and any **validation** run must override it) |

| §10.3 at a longer horizon | `cargo run --release --example criterion_metric -- 6 8 1e-4 20` |
| §10.5 sibling noise | `cargo run --release --example sibling_noise -- 5 3` |
| §10.6 the FTLE cross-check | `cargo test --release --test xcheck -- --ignored ftle --nocapture` |
| §10.7 the bivariate colouring | `cargo run --release --example bivariate_colour -- 5 8 13 1e-4` |
| §10.8 panning | `cargo run --release --example pan_sequence -- 9 20000 1024 near-field` |
| §10.9 the two costings | `cargo run --release --example cost_and_anisotropy -- 5 8` |
| the criterion gates | `cargo test --release --test criterion -- --nocapture` |
| §10.9 the slice gallery | `cargo run --release --example slice_gallery -- 40000 1e-4 0.2 1024` |
| §12 the chart gallery as committed (k_frac = 1, the uniform-mode control) | `cargo run --release --example chart_gallery -- 40000 1e-4 0.2 1024 1.0` |
| §19 the chart gallery at the ranked default | `cargo run --release --example chart_gallery -- 40000 1e-4 0.2 1024 0.25` |
| §12 the decoder and preset gates | `cargo test --release --test charts -- --nocapture` |
| §14 the threshold diagnosis | `cargo run --release --example threshold_diagnosis` |
| §14.5-14.6 the hot rule swept | `cargo run --release --example hot_rule_sweep` |
| §15.3-15.4 the balanced march as committed (k_frac = 1, no truncation) | `cargo run --release --example balanced_march -- 800 4 1e-4 64 1.0` |
| §19.4 the balanced march at the ranked default | `cargo run --release --example balanced_march -- 800 4 1e-4 64` |
| §16 structure modes by error(B) | `cargo run --release --example structure_metric -- 6 8 1e-4 13` |
| §17 the slippy-map gates | `cargo test --release --test slippy -- --nocapture` |
| the refinement animations | `cargo run --release --example refinement_animation -- 40000 1e-4 0.2 512` |
| the four GLSL slices refining | `cargo run --release --example glsl_refinement -- 40000 1e-4 0.2 512 40` |
| §18.0 the wiring check | `cargo run --release --example criterion_sweep -- 0` |
| §18.1 stage 1, tau x k_frac (widened in §19.1) | `cargo run --release --example criterion_sweep -- 1` |
| §19.3 selective or sparse, by error(B) | `cargo run --release --example equal_budget -- 6 8 5` |
| §19.5 the undetermined-footprint census | `cargo run --release --example nan_probe -- 128` |
| §18.3 stage 2, structure x criterion | `cargo run --release --example criterion_sweep -- 2 40000 0.2 1024 1e-4 0.25` |
| §18.2 stage 3, alpha | `cargo run --release --example criterion_sweep -- 3 40000 0.2 1024 1e-4 0.25` |
| §20 the ceiling, uniform, and the audit | `cargo run --release --example oracle_audit -- 7 8 1e-4 13 all` |
| §21 the sync-cadence artefact | `cargo run --release --example sync_artefact -- 3 8 13 1 0` (unguarded) and `... 3 8 13 1 1` (guarded) |
| §22 the closure criterion | `cargo run --release --example escape_closure -- 24` |
| §22.5 the toggle renders | `cargo run --release --example closure_render -- 1024 results` — writes `results/closure/`, **never** over the committed "before" set |
| §23.3-23.7 the `dtau` measurement | `cargo run --release --example dtau_step -- 48` — stdout only, `results/output/dtau_step.txt` |
| §23.8 the escape-arm correction | inside `tests/outcome_encoding.rs::escape_matches_the_legacy_classifier`; the 109-vs-0 comparison is `t = 20`, `n_sync = 32`, near-field 32x32, both `DtauMode` arms |
| §23.11 the mode-discrimination tests | `cargo test --release --test dtau_step_control` — four tests, each naming the configuration in which it fires |
| §23.10 the cross-check, three ways | `cargo test --release --test xcheck -- --ignored` and again with `PRIN_DTAU_MODE=fixed`. The mode is in the TSV header and `compare.py` asserts headers match, so a mismatched pair is refused rather than compared |
| §23 renders | `cargo run --release --example dtau_render -- 1024 results all dtau_fix` — writes `results/dtau_fix/`, **never** over the committed "before" set. Args are `res root only sub`; `only` filters to one case |
| §24.4-24.7 the overshoot measurement | `cargo run --release --example overshoot -- 48` — stdout only, `results/output/overshoot.txt`. §1 is the figure-eight and needs no grid; it runs in under a second and is the cheapest confirmation the clamp works |
| §24 the clamp tests | `cargo test --release --test overshoot_clamp` — three tests, each carrying the control arm that says it could have failed |
| §24 renders | `cargo run --release --example overshoot_render -- 1024 results all overshoot_fix` — writes `results/overshoot_fix/`, **never** over `results/dtau_fix/`, which is this comparison's "before". Args are `res root only sub arms`; `arms` is a subset of `ABCD`, so an interrupted run resumes rather than repeating ~15 min per arm. **The pair table needs all four arms in ONE run** — `... 1024 results config_stability overshoot_fix ABCD` — and must be read at 1024², never at 48x48 (§24.6) |
| §24.3 the cross-check under both knobs | `PRIN_DTAU_MODE=per-step-interval PRIN_CLAMP_FINAL=1 cargo test --release --test xcheck -- --ignored`, and again with `fixed`/`0`. Both read **4/4 PASS**. A deliberately mismatched pair is refused by `compare.py`'s header assertion, which now carries `clamp_final` as well as `dtau_mode` |
| the reference smoke test | `reference/README.md`. it now runs **all four arms**. `fixed`+overshoot reproduces the committed pre-fix file bitwise at `3.1966735584829495e-09`; the four read `3.197e-9`, `4.463e-10`, `4.047e-9`, `4.153e-10` — and the clamp moves the median *up* under `fixed`, which is §24.5's point about the drift field being blind to this defect. The number the README used to quote, `3.892633125701676e-09`, does not reproduce on the unmodified committed reference either |
| §21.6 persistence and precedence | `cargo run --release --example escape_persistence -- 48 13 32` |
| the signal audit against the DP labels | `cargo run --release --example signal_audit -- 7 8 1e-4 13 <scratch>` (~2 h; it enables the FTLE march, which is ~2.5x a plain build. **Point it at a scratch root**) |
| §20.6 firing the bound assert cheaply | `cargo run --release --example criterion_metric -- 4 8 1e-4 13 /tmp/scratch` |
| §20.1 the bound over every ranking | `cargo test --release --test criterion no_ranking_beats -- --nocapture` |
