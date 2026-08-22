# Results

What the uniform kernel measured, written to be read cold. Raw output for every number is in
[`results/output/`](results/output/); the mechanisms behind them are in
[`NOTES.md`](NOTES.md), and the working agreement they feed back into is
[`CLAUDE.md`](CLAUDE.md).

Everything below is Burrau (masses 3-4-5, released from rest), Aarseth–Zare with two-pair
Levi-Civita regularisation, `E+1 = 8` copies per pixel, `t_max = 13`, f64 unless stated.

---

## 1. The refinement criterion resolves regions, not quads

**This is a scheduler design constraint, not a defect.** The criterion works. It works at a
resolution the ensemble size sets, and that resolution is measurable in advance.

The criterion compares a parent quad against its children. A fine uniform grid already contains
every coarser scale by aggregation, so the whole exponent machinery is testable with no
quadtree: pool a 2×2 block's copies to synthesise the parent and compare against the children.
Nothing here can be an artefact of a scheduler, because there is no scheduler.

`alpha = log2(spread_parent / spread_child)`, child taken as the median over the four.

### What it discriminates

`alpha` for `spread_shape`, 64×64 fine grid → 32×32 parents, `t = 13`:

| region | min | p10 | **median** | p90 | max |
|---|---|---|---|---|---|
| near-field | −0.6678 | −0.0989 | **0.1722** | 0.5324 | 10.7910 |
| body2 core | −0.5254 | 0.0551 | **0.3390** | 0.7312 | 10.7508 |
| mid-field | 0.7427 | 0.9546 | **1.1781** | 1.3930 | 1.7517 |
| far | 0.6000 | 0.9276 | **1.1716** | 1.4161 | 1.8814 |

In the tame regions `alpha ≈ 1.17`: the shape spread scales with cell width, as a smooth field
must, so halving the cell halves the spread and refinement pays. In the chaotic regions
`alpha ≈ 0.17–0.34`: the parent spread is barely above the child's, so refining buys almost
nothing — those pixels are **undetermined, not under-resolved**. Making that distinction is what
the criterion is for, and it makes it cleanly.

### The cost curve

`alpha` for `sigma_E(0)` is a control whose true value is **exactly 1.0**: `sigma_E(0)` is
proportional to the jitter and therefore to the cell width, so doubling the cell doubles it.
Measured, near-field:

| E+1 | estimator | p10 | median | p90 | **p90−p10** | matched-count median |
|---|---|---|---|---|---|---|
| 8 | rms | 0.8525 | 1.0762 | 1.3321 | **0.4796** | 1.0191 |
| 8 | max_dev | 0.8787 | 1.1248 | 1.3855 | 0.5067 | 0.9865 |
| 16 | rms | 0.8614 | 1.0137 | 1.1817 | **0.3203** | 0.9999 |
| 32 | rms | 0.8902 | 0.9887 | 1.1025 | **0.2123** | 0.9884 |
| 64 | rms | 0.9169 | 0.9862 | 1.0663 | **0.1493** | 0.9949 |

The interdecile width is the per-quad noise floor: how far a single quad's `alpha` scatters when
the true value is exactly 1. It falls as `1/sqrt(E)` — measured ratios 0.667, 0.663, 0.703
against the predicted 0.707.

**The design constraint, stated as a trade:**

- `E+1 = 8` buys a per-quad noise floor of **0.48** in `alpha`.
- The measured region separation is about **1.0** — roughly twice that floor.
- So `E+1 = 8` resolves **regions** comfortably and **individual quads** not at all.
- **Halving the floor costs 4× the compute.** `E+1 = 32` gives 0.21; `E+1 = 128` would give
  about 0.075 for sixteen times the work of 8.

A scheduler that thresholds `alpha` per quad at `E+1 = 8` is acting on noise. One that
thresholds on a region aggregate is not. If per-quad decisions are wanted, the ensemble size is
the price and the curve above says what it is.

### The bias only a control could have caught

There is also a **+7.6% median bias** at `E+1 = 8`. Its cause: a parent pools `4(E+1)` copies
against a child's `E+1`, and a spread estimator's expectation depends on sample size. Drawing
`E+1` of the parent's pooled copies puts matched counts on both sides, and the median goes to
1.0191 at `E+1 = 8` and 0.9999 at 16.

**No estimator without a known-true-value control could have found this.** The bias is a smooth
7.6% offset on a quantity with no independent oracle — it does not look like an error, it looks
like a result. It would have been read as "parent spread grows slightly faster than linearly in
cell width", which is a plausible physical statement and completely wrong. It was visible only
because `sigma_E(0)` has an *exactly* known answer.

The same reasoning chose the estimator: this experiment uses an rms deviation rather than
`error_ratio`'s max deviation, because an order statistic's sample-size bias is far larger — the
`max_dev` rows above carry it too, and worse.

**Match the counts, or the exponent is biased before any physics enters.**

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

## 5. What this invalidates in the prior numpy investigation

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

## 6. Reproducing any of this

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
| a full slice | `cargo run --release --bin prin -- --region near-field --size 256 --out out` |
| the acceptance gates | `cargo test --release -- --nocapture` |
| the NumPy cross-check | `cargo test --release --test xcheck -- --ignored --nocapture` |
| the horizon table | `python3 tools/xcheck/horizon.py [--lc-unstable]` |
