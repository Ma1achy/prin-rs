# Four step-control candidates — the numbers

Reproduce:

```sh
cargo run --release --example step_limit sample 128 48 results   #  ~6 min
cargo run --release --example step_limit frame  192 48 results   # ~25 min
cargo run --release --example warp_divergence   128    results   #  ~7 min
cargo run --release --example subsume           128 64 results   #  ~2 min
```

## 0. What is being fixed, and what is not a fix

`dt = A*B*dtau` is emergent — `dt/dtau = A*B` is integrated *by* the RK4 stepper — so the step
taken is not the step predicted. One step advanced the physical clock by **`2.209e128`** against a
sync interval of `0.4`, and the march recorded a clean landing. **Nothing asked whether the step
it just took was one it could afford.**

`eta/256` and a fourth refinement pass both repair this and **neither is a remedy**: the first
pays 256x everywhere for a local failure, the second is `refine_flagged` again — re-integration
from `t = 0`, which a live playhead cannot do. They are characterisation, and their value is that
*`eta/256` brings every flagged pixel to `error_ratio` 1.000* proves this is ordinary
under-resolution and not a wrong equation, a cap, or a threshold. That is why `eta/256` is the
**ground truth** here rather than a candidate.

## 1. Read `steps`, not `secs`

The plan called for `error_ratio` p99 against wall clock. **Wall clock on this machine is not
trustworthy** — load average ran 85–100 throughout, and `Predictive f=0.1` on `config_stability`
timed *faster than the baseline* (77.6 s against 136.9 s) while doing 1.7% **more** work. Every
cost figure below is `total_substeps`, which is machine-independent and was collected anyway.
Wall clock is printed and should not be quoted.

## 2. The decision — full frame, 192², `config_stability`

```text
   mode                steps p50      err p99   err>10   overshoot   retries
   None (baseline)       1.033e5      7.108e9   0.1110         634         0
   Predictive f=0.02     1.053e5        1.109   0.0000           0         0     <- +1.9%
   Reject     f=0.02     1.832e5        1.205   0.0000           0     5.5e9     <- +77%
   AbGrowth   f=2        1.033e5      7.108e9   0.1110         634         0     <- inert
   Global     f=0.25     4.134e5      9.066e7   0.0767         153         0     <- +300%
```

**B wins outright.** `error_ratio` p99 from `7.1e9` to `1.109`, the fraction above the flag
threshold from 0.1110 to **0.0000**, and the overshoot count from **634 to zero** — for **1.9%
more steps**. One divide per step, no trial step, no retry, no branch.

**And D — the dumb control — does not even fix it.** At four times the work `Global f=0.25` still
carries **153 overshoots** and `err>10 = 0.0767`. A uniform `eta` cut buys accuracy everywhere and
still fails to bound a step whose size is set by local geometry. That is the strongest single
argument for a per-step limit over a global one, and it comes from the control rather than from
the candidate.

## 3. The controls — it is not tuned to one slice

```text
   preset_shape          steps p50      err p99   err>10   overshoot
   None (baseline)         3.021e3     5.538e10   0.0824         584
   Predictive f=0.02       3.037e3        1.059   0.0000           0    <- +0.5%
   Reject     f=0.02       1.558e4      2.337e5   0.0259           0    <- +416%, PLATEAUS
   Global     f=0.25       1.184e4      4.580e8   0.0564         103    <- +292%

   deep interior         steps p50      err p99   err>10   overshoot
   None (baseline)         1.080e4      2.664e2   0.0381          84
   Predictive f=0.1        1.408e4        1.002   0.0000           0    <- +30%
   Reject     f=0.1        1.661e4        1.009   0.0000           0    <- +54%
   Global     f=0.25       4.293e4        1.047   0.0034           0    <- +298%
```

**`preset_shape` was supposed to be the clean control and is not.** Its baseline carries
`err>10 = 0.0824` and **584 overshoots** — the defect is there too, and "the chart families are
tame" was a statement about `alpha`, not about the step control. B fixes it for +0.5%.

**A plateaus above 1.0 on `preset_shape`.** Going from `f = 0.1` to `f = 0.02` moves `err p99`
from `2.087e5` to `2.337e5` — *up* — while costing 4x. The mechanism runs out of retries: 39 of 96
sampled pixels exhaust `MAX_RETRIES` and are marked **undetermined**. A halving ladder bounded at
8 cannot reach where a single well-chosen step goes directly. **A does not address the defect
here**, and that is the plan's own stated criterion for saying so.

## 4. A is not viable on a GPU — and the control is what makes the number mean anything

`config_stability` 128², one lane per pixel, warp 32:

```text
   mode              retry p90   retry max   div linear   div tiled   warps hit
   None (control)            0           0        1.577       1.432      0.0000
   Reject f=0.5            152       16013        1.611       1.446      0.9824
   Reject f=0.1          10619      394459        1.766       1.543      0.9980
   Reject f=0.02        368118     5205454        2.557       1.975      1.0000
```

The divergence factor is `mean(max per warp) / mean(per lane)`. **The absolute level is the
field's, not the mode's** — step counts vary lane to lane with no retries at all, which is why the
control row reads 1.577. The increase is A's cost: **+62% linear, +38% tiled** at `f = 0.02`.

**The killer is `warps hit`.** At the parameter A needs to work, **every warp contains a retrying
lane** — 1.0000, under both dispatch shapes. A warp executes in lockstep and pays its worst lane,
and the worst lane retried 5.2 million times. Rare retries that scatter perfectly are the bad
case, and this is that case at saturation. **CPU wall clock hides it completely**: A costs +8% of
wall clock at `f = 0.1` and would cost far more on the hardware it is meant for.

## 5. C was already shipped, under another name

`StepLimit::AbGrowth` is **bitwise identical to the baseline** at every parameter on every region.
The brief's formula assumes `dtau` is fixed across the interval, which is `FixedPerInterval`;
under the shipped `PerStepInterval`, `dtau = eta*dt_left/(A*B)` is recomputed every step so
`dt ~ eta*dt_left` *however much `A*B` grows*. **`DtauMode::PerStepInterval` IS an `A*B` growth
clamp at `C = 1`.** `tests/step_limit.rs` holds both halves — inert under `PerStepInterval`,
active under `FixedPerInterval` — so if the older mechanism ever changes, the thing that silently
starts mattering announces itself.

## 6. And the cap can now be removed — reported, not done

With B in force, `dtau_mode` barely matters:

```text
   region              dtau_mode          step_limit     err p99   err>10   steps p50   overshoot
   config_stability    PerStepInterval    Predictive       1.555   0.0000     8.988e4           0
   config_stability    FixedPerInterval   Predictive       1.686   0.0000     8.416e4           0
   deep interior       PerStepInterval    Predictive       1.000   0.0000     3.831e4           0
   deep interior       FixedPerInterval   Predictive       1.000   0.0000     3.831e4           0
```

`deep interior` is **identical to four digits including the step count**; `config_stability` moves
`err p99` `1.555 -> 1.686`, both far below the flag threshold of 10, on **6% fewer steps**. So the
cap is redundant once the per-step limit is in force and one of the three controllers could be
deleted. **Reported and not done** — that is a second corpus-invalidating change and it belongs in
its own measurement, next to the question of whether the GLSL app's `N_sub` bucket goes the same
way. It very likely does: the bucket, the cap and the limit are three controllers on a `dt` that
AZ already adapts.

`clamp_final_step` is **not** a candidate for removal. It is a correctness property, not a
step-size one: it takes the measured convergence order from 1.06 to 2.08 on the figure-eight.

## 7. The tripwire, and it is permanent

`AzOut::n_overshoot` counts steps after which the interval clock passed `2 * dt_left` —
`debug_assert`ed in debug, counted in release, aggregated to `PixelOut`, and printed in the `.raw`
header. `dt > dt_left` is a bug, not a condition to handle. It is conditioned on
`clamp_final_step`, because with the clamp off overshoot is the expected behaviour of a named
measurement axis and an assert that fires on a deliberate mode is a broken assert.

It reads **634 / 584 / 84** on the three baselines and **0** under B everywhere.

## 8. What is NOT established

- **The ground-truth chord and label comparison is saturated and says nothing.** `flips` reads
  1.0000 against the `eta/256` truth for a *correct* mode and a broken one alike: over `t = 50`
  any change of step size gives a different trajectory through a chaotic region. It is reported as
  `NaN`/saturated rather than quoted. What discriminates is `error_ratio`, which is normalised to
  1.0 under exact dynamics and therefore has an absolute meaning.
- **`f = 0.02` is the finest rung measured, not an optimum.** Between `0.1` and `0.02` the cost
  rises 1.7% and the p99 falls 635x on `config_stability`; nothing between was tried and no
  argument says 0.02 is where the curve turns.
- Wall clock, for the reason in §1.

## 9. Shipped, and what it costs the corpus

`EnsembleCfg::production()` now carries `step_limit: Predictive, step_limit_f: 0.02`.
`StepLimit::None` stays reachable and named, and `reference_opts` pins it so the NumPy
cross-check is unaffected — **4/4 PASS**.

**The committed corpus was taken under `None` and does not reproduce bitwise under this default.**
That is the cost and it is stated rather than discovered: `provenance()` names the setting in
every header and sidecar, so any future disagreement is a value in a log rather than six days of
copy-paste.

### Three tests failed, and every one of them failed correctly

- `refined_pixels_are_repaired` fell over on its **own** guard — *nothing was flagged, so this
  test has no subject*. The assertion written to catch a vacuous test caught one.
- `mad_based_error_ratio_cannot_separate_damaged_pixels` had no damaged population left to fail
  to separate.
- `error_ratio_minus_one_falls_with_step_size` went non-monotone because the residual is now
  **round-off, not truncation**: `3.1e-9` at `eta = 4e-2` falling only to `1.6e-9` across a
  decade.

All three are pinned to `StepLimit::None`, which is the kernel they are about. **That the limit
deletes the subject of three characterisation tests is the strongest corroboration in the suite**
— stronger than any single error digit — and the pins are recorded with that reason rather than as
a tolerance loosened to keep green. `error_ratio.rs` gains the arm the pin depends on:
`the_shipped_limit_is_already_at_the_floor_at_the_coarsest_step`, which asserts an eightfold step
cut moves the residual by **less** than an order — without it, the pin would read as a test
weakened to pass.

`cargo test --release`: **256 passed, 0 failed.**
