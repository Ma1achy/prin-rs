# The saturation check — what stops the march, and does it draw the artefact

Reproduce:

```sh
cargo run --release --example saturation_mask 512 results   # 1087 s, results/output/saturation_mask.txt
cargo run --release --example cliff_ladder 256 5 128 results # 336 s, results/output/cliff_ladder.txt
```

Settings are `closure_render`'s — `t_max = 50`, `n_sync = 125`, `r_coll = 0.005`,
`EscapeRule::Closure`, `stop_on_escape = false` — with **`refine_flagged` OFF**, deliberately:
that is the state the committed panel is in, and repairing it first would switch off the
diagnostic being measured. Every panel here carries a `.cfg.txt` sidecar naming its config.

## 1. The hypothesis, and its three concrete forms in this port

The proposal was that the cut-out regions are a *saturation* boundary: a substep cap is hit, the
wrapper stops refining and **advances anyway** with a step it knows is too coarse, and the
boundary where the cap engages is a sharp edge in IC space.

`principia_integrator_contract.md`'s `substep_bucket` / `N_sub` / `N_max` / descriptor bit 5 **do
not exist in this port** — that is the GLSL app's contract. Asked in this codebase's own terms —
*where does this integrator knowingly advance with a step it cannot resolve, and is that
recorded?* — there are three sites:

| site | behaviour | was it recorded? |
|---|---|---|
| `max_steps` per interval | **terminal** — the march breaks, the run ends | yes, as `PixelOut::state == 6` |
| `A`,`B` floored to `T::TINY` | **advance anyway**, fabricated denominator | `AzOut::ab_floored` — **computed on every march and read by nothing** |
| `.min(dtau_entry)` cap | **advance anyway**, refused the step it asked for | not recorded at all |

The middle row was the finding at the plumbing level: the sticky bit existed, was written, and
stopped one layer below `PixelOut`. *A sticky bit that nothing reads is indistinguishable from
one that never fires.* All three are now on `PixelOut`, along with `dt_max`.

## 2. The mask — refuted, in all three forms

```text
              mask   fraction     P(m|err>10)     P(m|err<=10)     lift
        ab_floored   0.000000        0.000000         0.000000    0.000
      n_cap_hits>0   1.000000        1.000000         1.000000    1.000
  budget_exhausted   0.000000        0.000000         0.000000    0.000
    error_ratio>10   0.111061              --               --       --
```

Two never fire on this slice. The third fires on **every pixel of 262144** — it is saturated, and
a saturated mask has a lift of 1.000 by arithmetic, not by measurement. The base-rate row is
printed first for exactly this reason. `capped` firing everywhere is not a fault: it fires
whenever `A*B` falls below its value at the interval's entry, which is what happens every time
bodies approach mid-interval. It is routine, and `tests/saturation_plumbing.rs` holds that it
fires under `PerStepInterval` and never under `FixedPerInterval`.

**The saturation hypothesis is refuted as a discriminator here.** No mask matches the artefact
because two are empty and one is full.

## 3. The cliff is a SLOPE, not a floor — the prediction inverts

`error_ratio > 10` on the flagged population, `eta` over four decades, 128 of 7289 flagged pixels
(evenly spaced, cap printed):

```text
         eta     err p50     err p90 sigE(t) p50   drift p50    slope   steps p50 dt_max p50
    1.000e-2     2.130e5     3.389e9     1.168e1     8.609e1       --     8.997e4   5.739e-3
    2.500e-3     6.139e3     1.235e8    2.380e-1     2.063e0    +2.56     3.609e5   1.054e-3
    6.250e-4     5.176e0     1.159e6    3.144e-4    3.471e-3    +5.11     1.613e6   2.517e-4
    1.563e-4     1.000e0     2.846e1    6.391e-5    3.638e-6    +1.19     5.690e6   6.259e-5
    3.906e-5     1.000e0     1.000e0    5.826e-5   3.922e-14    +0.00     2.400e5   1.562e-5
```

**This is characterisation, not a remedy.** A global `eta/256` pays 256x everywhere for a local
failure and does not survive a live playhead; it is used as the **ground truth** in
`results/step_control/` rather than as a candidate. What it establishes — and this is worth the
run on its own — is that the failure is ordinary under-resolution: not a wrong equation, not a
saturating cap, not a threshold.

Median `error_ratio` **2.13e5 -> 1.000** and median drift **8.6e1 -> 3.9e-14**, thirteen orders.
The p90 converges too. The control sits at 1.000 throughout. The `+0.00` at the last rung is not
a floor: 1.000 is `error_ratio`'s *converged* value, since `sigma_E(t) -> sigma_E(0)` under exact
dynamics. The prediction — *the cap is on `N_sub`, not on `eta`, so refining `eta` cannot raise
it* — is **wrong for this port**. Refining `eta` clears it completely.

## 4. And `error_ratio p99 = 35.6` is the pass count, not a mechanism

```text
  cleared at rung 1 (eta = 2.500e-3):   47   (0.3672 cumulative)
  cleared at rung 2 (eta = 6.250e-4):   21   (0.5312)
  cleared at rung 3 (eta = 1.563e-4):   37   (0.8203)   <- shipped refine_max_passes = 3
  cleared at rung 4 (eta = 3.906e-5):   23   (1.0000)
  NEVER cleared: 0 of 128 (0.0000)
```

A fourth pass is `refine_flagged` again and is a **batch** mechanism for the same reason, so this
row too is characterisation. **82.0% clear by the third pass**, so ~18% survive the shipped ladder — which is the p99 = 35.6
tail, exactly. The repair pass converges; it is stopped one rung early. `refine_max_passes = 5`
clears all 128 of this sample. That is a different claim from "the pass does not repair", and the
run was built to be able to make either.

## 5. What the plumbing DID find, and it is not a cap

`dt_max` is the largest **physical** step any copy took, as an actual `s.t` difference across one
RK4 step. Nominal is `eta * t_max / n_sync = 4.0e-3`.

```text
    err>10   dt_max p50 6.074e-3   p99 1.263e43   max 2.209e128
   err<=10   dt_max p50 4.183e-3   p99 1.874e-2   max 7.241e-2
```

**One RK4 step advanced the physical clock by up to `2.209e128` against a sync interval of
`0.4`.** The march then accepted it as a clean landing, and this is code, not inference:

- `bad = !s.is_finite()` — `1e128` is finite, so the divergence guard passes.
- `s.t >= dt_left - land_tol` — satisfied by 128 orders, so `landed = true`.
- `t += dt_left` under `clamp_final_step` — **the clock is corrected to the boundary while the
  state kept is the one reached at `s.t = 1e128`.**

The clamp corrects the clock and cannot un-take the step. `t` is clamped on both branches, so the
overshoot was invisible in **every** recorded quantity until `dt_max` was plumbed. This is a real
advance-anyway defect — the shape the hypothesis predicted, at a site it did not name, and it is
an **unbounded step with no acceptance test** rather than a cap.

It is also curable by `eta` (§3), which is what makes it under-resolution rather than a structural
ceiling. What it is not is *detected*: nothing in the march asks whether the step it just took was
one it could afford.

`field_dt_max.png` renders it — coherent bright ridges tracing the wedge structure, threshold-free,
and matching `mask_errhot.png`, which is the `error_ratio > 10` set and is the artefact: solid
wedges with straight edges, cut out of the ribbons, resuming outside.

## 6. Per box — `error_ratio` orders them, `dt_max` follows

```text
   box     verdict     err>10  dt_max p50        box     verdict     err>10  dt_max p50
    P5      BROKEN     0.9775    6.309e-3         P2       sound     0.0235    4.246e-3
   B10      BROKEN     0.8943    2.076e-2         P3       sound     0.1436    4.460e-3
    B9      BROKEN     0.9112    8.594e-3         P1       sound     0.2022    4.238e-3
    B2      BROKEN     0.7314    5.886e-3         B7       sound     0.3144    4.023e-3
    B1      BROKEN     0.6523    4.986e-3         B4  SOUND-CTRL     0.3163    4.548e-3
    B8      BROKEN     0.5964    4.888e-3      FRAME    baseline     0.1111    4.228e-3
```

The BROKEN set runs 0.60-0.98 against the sound set's 0.02-0.32 and a frame baseline of 0.11 —
`BUG_REPORT.md` §4's verdicts, recovered from an independent quantity. **The negative control
qualifies the result rather than confirming it**: `B4` reads 0.3163, three times the frame
baseline, so it is not pristine — it is simply far below the broken set. Its wedge survives
refinement and is real structure; that stands, and "sound" here means *not this fault*, not
*unflagged*.

## 7. Open, and not explained

- **Terminal labels move with `eta` and the step count collapses at the finest rung** (`steps p50`
  5.69e6 -> 2.40e5 for a 4x finer step). Trajectories are terminating differently, not merely more
  accurately. Consistent with the standing result that `hot` falls with `eta` while the labels do
  not converge; not investigated here.
- The ladder sampled **128 of 7289** flagged pixels at 256². The `2.209e128` is from the full
  512² frame and its pixel is not necessarily in the ladder's sample.
- No remedy is applied. The shape one would take — reject and retry a step whose taken increment
  exceeds its own remaining interval — is local and needs no re-integration, so unlike
  `refine_flagged` it has a live-playhead analogue. That is a design question, not a measurement,
  and it is not settled here.
