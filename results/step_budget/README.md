# The production step budget is sized for AZ, and the default integrator is Heggie

`examples/step_budget.rs` -> `output/step_budget.txt`. 35 targets, a fixed grid per rung,
`max_steps` the only knob. Nothing here is a change; it is the measurement asked for before one.

## What was being looked at

The chart gallery was rendering a scatter of magenta pixels, then -- after the debug flag was
taken out of the presentation path -- a scatter of dark ones in the same places. Neither is in the
Heggie regularisation panels. The difference is not the refinement work and not the colouring:

`EnsembleCfg::production()` carries `max_steps: 30_000`. `integrator_gallery` raises it to
**400_000 for both arms** and says why in its own header -- *"at the production
max_steps = 30_000 Heggie exhausts the budget on 8.6% of `config_stability` and its drift panel
comes back dominated by the magenta veto set"*. **27 diagnostic harnesses raise it; every render
and scheduler harness silently took 30_000.** So a footprint undetermined in a gallery panel is
determined in an integrator panel, and nothing took the difference as an argument. The default
integrator moved to Heggie on 2026-09-02 (`84830a1`) and its step budget did not follow -- the same
shape as `refine_flagged` spreading by copy and `k_frac = 1.0` shipping as the default.

## The design, and every choice here is a defect this project has already paid for

**A FIXED UNIFORM GRID, NOT A TREE.** The first pass read the flagged count off `chart_gallery`'s
live march, where the tree changes with the budget -- so the denominators ran 29632, 22912, 23296,
22144 and the "fraction" was over a moving population. Every rung here evaluates the identical
`48^2` grid, so a difference in the numerator is a difference in the physics.

**`steps TOTAL`, NOT `steps p50` -- AND THE MEDIAN IS THE TRAP.** The budget caps only the tail, so
the median trajectory is untouched **by construction**: `steps p50` is identical at every rung on
every one of the 35 targets. The first cut printed only that column, which would have been *a
statistic that cannot move being read as a statistic that did not move*. `steps TOTAL` is the cost.

**THE CONTROL IS A CHART THE BUDGET NEVER BINDS ON.** `latent_shape` is **bitwise constant** across
all five rungs -- flagged, budget, and total substeps to the digit. If it moved, every other row
would be suspect. The harness asserts it, and when a run selects only affected charts it reports
`NO CONTROL` rather than a clean sweep.

## The result

**Every target reaches `flagged = 0` at 480_000. Not one exception. And `flagged == budget`
exactly in every row of all 35 targets** -- 100% of it is budget exhaustion, with no second
mechanism anywhere in the sweep.

The ladder on the two worst, with the control:

```
                30k     60k    120k    240k    480k     steps TOTAL at clear
latent_shape      0       0       0       0       0        1.000x  (bitwise, all 5 rungs)
latent_mixed_h3 195      16       0       0       0        1.137x  at 120k
burrau_nu_k     521      65      65      57       0        1.351x  at 480k
```

`burrau_nu_k` has a long tail: the **same 65** footprints are capped at 60k and at 120k -- identical
rows but 319.7M against 341.1M substeps, so they burn the larger cap without completing -- 57 at
240k, and 0 by 480k.

12 of 35 targets flag at 30_000. **All eight Burrau regions read 0 at both rungs, bitwise**, which
is why this never appeared in the corpus or the scheduler tests: the budget does not bind there. It
binds on the latent charts.

```
                 flagged @30k       cost at 480k
burrau_nu_k         521  22.61%        1.351x
latent_mixed_h3     195   8.46%        1.137x
latent_oblique_a     19   0.83%        1.273x
latent_oblique_b     10   0.43%        1.143x
preset_shape_h1      10   0.43%        1.011x
config_stability      5   0.22%        1.000x
preset_shape          3   0.13%
latent_mass_h3        3   0.13%
mass_simplex          1   0.04%
preset_plambda_h1     1   0.04%
preset_shape_pl_h1    1   0.04%
latent_shape_h3       1   0.04%
```

Unaffected targets are **bitwise 1.000x** -- raising the budget cannot cost anything where it never
bound.

## AND IT IS NOT COSMETIC: THE BUDGET MOVES THE TREE

A truncated footprint reads as *undetermined*, so its quad can never resolve, so it splits. **A
too-small budget manufactures work rather than saving it.** Measured through `chart_gallery` at
192^2, static tree leaf counts:

```
burrau_nu_k       415 leaves at 30k, 373 at every rung from 60k up
latent_mixed_h3   214 leaves at 30k, 142 at 480k   -- and the whole chart ran 2.5x faster
```

That inflation lands in leaf counts and stop-reason breakdowns, which are the numbers the gallery
table exists to report.

## Read the drift columns with the selection in mind

`burrau_nu_k`'s `drift p50` runs `1.358e-9 -> 4.065e-9` across the ladder, which is not a
degradation. At 30_000 the 521 truncated footprints have non-finite drift and are excluded from the
quantile; at 480_000 they complete, contribute honest larger values, and the median rises. *A median
conditioned on each arm's own selection is a different statistic in each row* -- on this record
already, and live here.

## What is NOT settled

Whether `EnsembleCfg::production()`'s value should move. That is corpus-invalidating and is not
taken here. What the sweep supports: 30_000 is not defensible for the shipping integrator on the
latent charts, 480_000 clears every target measured, and the cost is 1.000x wherever the budget
never bound and at worst 1.351x where it did. `integrator_gallery`'s own 400_000 sits between two
tested rungs and is **not** measured here.

Also unmeasured: whether the residual on `burrau_nu_k` is a budget question at all rather than a
trajectory that needs a finer `eta`. `error_ratio` p99 there reads `1.24e2` at 480_000, so
something on that chart is not data whatever the budget.

**Corrected 2026-09-06 — this paragraph originally read "`1.24e2` at every rung", and both the
word "tail" and the word "every" were wrong.** The column is not flat: it runs
`1.127e2, 1.170e2, 1.170e2, 1.242e2, 1.244e2`, rising monotonically as truncation is removed. Two
readings of that were wrong on the way, and the second is the finding:

- **It is not a selection artefact.** The obvious explanation was that each row's quantile is taken
  over its own finite population — `.filter(|x| x.is_finite())` — and `stats::max_dev` returns
  `+inf` for a non-finite copy, so the 521 truncated footprints at 30_000 looked excluded from
  exactly the rung where they exist. Measured over the population finite at **every** rung:
  **2304 of 2304, 0.0% dropped.** Nothing was ever filtered. A starved footprint does not go
  non-finite — every copy stops at the same early point and agrees perfectly, which is the standing
  *a starved footprint reads exactly 1.0000* result, and 1.0000 is finite.
- **It is not a tail.** Fixing the population left the p99 where it was and exposed the **median**:
  `err p50` runs **`1.082e0 -> 3.189e1`**, a thirty-fold rise, while p99 moves 10%. At 30_000,
  22.6% of the frame was truncated and reading ~1.0, dragging the median beneath it. So this is not
  a handful of pixels that are not data — **more than half the chart sits at `error_ratio` around
  32**, against a `refine_threshold` of 10, and the p99 is the top of a broad bulk rather than an
  outlier.

`examples/step_budget.rs` now prints the fixed-population block with `n` beside every quantile;
`examples/err_ratio_residual.rs` carries the decomposition, and it settles the question.

## The residual is the `t_end` READOUT, not the integration

`results/output/err_ratio_residual.txt`, 48^2, production kernel, `burrau_nu_k` against
`near-field` as the control. Argument four is `stop_on_event`.

```
                        stop_on_event=1        stop_on_event=0
  burrau_nu_k
    ratio p50               3.1888e1               8.7151e-1
    ratio p99               1.2435e2               1.5659e0
    sigma_e_0 p50           1.565e-2               1.565e-2     <- bitwise, it is at t=0
    sigma_e_t p50           6.381e-1               1.270e-2
    re-integrated           1385 (60.1%)           0 (0.0%)
    terminated              1424 (61.8%)           1424 (61.8%)
  near-field (control)
    ratio p50               1.0000e0               1.0000e0
    terminated              64 (2.8%)              64 (2.8%)
```

**A trajectory stopped by an event is parked *at* a close approach, and its Cartesian energy there
is a cancellation of two enormous terms.** `error_ratio` is built from exactly that readout
(`pixel.rs:829`), so the whole residual collapses when termination is off: the ratio falls 37x, the
numerator falls 50x, and **nothing on the chart exceeds `refine_threshold` any more**. This is the
standing *`stop_on_event` contaminated the drift diagnostic and produced five wrong conclusions in
a row* result, at a sixth site and on a different statistic.

The denominator is **bitwise identical** across the arms, which is what says the mechanism is in
the numerator and not in the cell-width artefact `sigma_e_0` is dumped separately to catch.

**And it explains why this chart and not the control.** The parked readout can only bite where
trajectories terminate before the horizon: `burrau_nu_k` 61.8%, `near-field` **2.8%** — the
standing result that near-field at `t = 13` is 97.8% *Bounded* and sits at `t_end = t_max`, so
there is no parked state to read. The re-integrated count (1385) tracks the terminated count
(1424) to within 3%: **the flagged population is the parked population.**

### Four candidates ruled out, one of them by its own sign

- **Step size.** 60.1% of the chart was re-integrated, three passes at `eta/4` to `eta/64`, and the
  ratio moved **0.5%** (`|ln ratio|` p50 `1.8e-3`). *Insensitivity to step size is the tell.* The
  `eta` ladder this diagnosis was authorised to build is therefore **guaranteed flat and was not
  built.**
- **A small denominator.** `sigma_e_0` p50 is `1.565e-2` on `burrau_nu_k` against `1.050e-3` on
  `near-field` — fifteen times **larger**, the opposite of the artefact.
- **Copies disagreeing on outcome.** The physical explanation predicts the disagreeing pixels carry
  the large ratios. They read p50 **1.241** while the 2144 pixels that **agree** read **32.18**.
  Refuted by the sign, not merely unsupported.
- **The energy scale.** `|E|` p50 is `3.060e-1`, so `sigma_e_t = 6.381e-1` is **twice the total
  energy** while relative drift p50 is `4.065e-9`. Those cannot both describe a conserved
  trajectory, which is what said the readout and not the integration before the flag was tried.

### A real defect found on the way, which is NOT the cause

`pixel.rs:824-831` builds the two energy vectors the ratio is taken over:

```rust
let e0 = copies.iter().map(|c| energy(&c.s.r, &c.s.v, &c.m, 0));    //  own masses
let et = outs.iter().map(|o| energy(&o.state.r, &o.state.v, &m, 0)); //  NOMINAL masses
```

`m` is `copies[0].m`. **Eleven lines above, the comment introducing it states the invariant this
violates in as many words** — *"per-copy masses are used where a copy is integrated **or its
energy taken**, because on a mass chart two jittered copies are two different systems."* On
`BodyPlane` every copy shares the nominal's masses, so the two spellings are one expression and the
seam is invisible: the equal-mass blind spot already on record for the Jacobi round-trip.

Measured at `t = 0`, where no integration and no chaos intervene: mass spread across copies is
`3.231e-3` on `burrau_nu_k` and **exactly zero** on `near-field`, and the ratio of the two
spellings is p50 **0.857**, p99 1.242 — against **exactly 1.0000 at every quantile** on the
control. So the seam is real and its size is **~1.24 at the tail, not 32**. It is credited with
nothing.

What it *does* explain is the leftover: with termination off the ratio reads **0.8715** where a
clean chart reads 1.0000, and the independent `t = 0` measurement of the seam alone reads
**0.8573** — agreeing to 2%. Two routes to the same number, which is what makes it the seam rather
than a second unknown.

### Not taken

Per the decision this was run under: **diagnose and report, change no production setting.** Both
findings are corpus-invalidating — `stop_on_event` decides the terminal class, which *is* the
physics for `_uniform`/`_outcome`, and the mass fix moves every `error_ratio` on the four chart
families that vary masses. Neither is a bug in the integrator: the trajectories are fine and the
statistic reads them at a point where the coordinates are ill-conditioned.

The narrow consequence, stated so it is not lost: on a chart that terminates heavily,
`refine_flagged` re-integrates the **parked** population at `eta/64` and cannot help it, because
the flag is not measuring integration error there. On `burrau_nu_k` that is 60% of the frame.
