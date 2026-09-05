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

Also unmeasured: whether the tail that survives 240_000 on `burrau_nu_k` is a budget question at
all rather than a trajectory that needs a finer `eta`. `error_ratio` p99 there reads `1.24e2` at
every rung including 480_000, so something on that chart is not data whatever the budget.
