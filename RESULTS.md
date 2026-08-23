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

## 9. Reproducing any of this

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
| the acceptance gates | `cargo test --release -- --nocapture` |
| the NumPy cross-check | `cargo test --release --test xcheck -- --ignored --nocapture` |
| the horizon table | `python3 tools/xcheck/horizon.py [--lc-unstable]` |
