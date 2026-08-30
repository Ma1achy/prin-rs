# BUG REPORT — the render harnesses disable the repair pass, and the committed panels carry the fault

**Status:** confirmed, reproduced, cause located, **fixed** — see §11. The propagation
defect of §2–§7 is separate and is a *discrepancy*, not the defect — and see §0, which changes what
this report is about. The propagation defect is real and is a *discrepancy between renders and
production*. It is **not** the defect that draws the artefact.
**Found:** 2026-08-29, while diagnosing pale/magenta patches on `config_stability`.
**Severity:** every committed render made by an `examples/*_render.rs` harness is affected.
The production binary `prin` is **not** affected.

---

## 0. WHAT THIS REPORT IS NOT — read before §7

`refine_flagged` re-integrates a whole trajectory from `t = 0` at finer `eta` **after the fact**.
Principia marches a live playhead: there is nothing to re-integrate *from*, and a pixel bad at
`t = 30` cannot be repaired at `t = 30.1` without redoing thirty time units. **It is a batch-only
workaround with no live-path analogue.**

So the reading inverts. The render harnesses with `refine_flagged: false` were showing the
**unmasked kernel**; `prin.rs`, on the default, has been hiding the same defect behind the repair.
Closing the propagation gap **closes a discrepancy, not the defect** — and doing it first would
switch off the diagnostic that revealed it. §7 is rewritten accordingly and is still not applied.

**One qualification, and it does not contradict the above.** The renders are honest about the
*kernel* and remain invalid as *science images*: `spread_shape` saturates because the copies
diverged to garbage, and the terminal class is reclassified (escape `0.4645 -> 0.0021` in `B10`).
`_drift.png` is where the unmasked kernel is legible; `_uniform.png` is not. Both are true.

### What was measured instead — `results/saturation/README.md`

The proposal was a **saturation** boundary: a substep cap engages, the wrapper advances anyway
with a step it knows is too coarse, and the cap's boundary is a sharp edge in IC space.

- **`AzOut::ab_floored` and `ab_min` were computed on every march and read by nothing.** They
  stopped one layer below `PixelOut`, so no render, dump, criterion or test could see the floor
  fire. Now plumbed, with `dt_max` and `n_cap_hits`.
- **The saturation hypothesis is refuted in all three forms this port has**: `ab_floored`
  `0.000000`, `budget_exhausted` `0.000000`, and `n_cap_hits > 0` on **every pixel of 262144** —
  saturated, lift exactly `1.000` by arithmetic.
- **The cliff is a SLOPE.** Over four decades of `eta` the flagged population converges
  completely: median `error_ratio` `2.13e5 -> 1.000`, drift `8.6e1 -> 3.9e-14`, **0 of 128 fail to
  clear**. Refining `eta` does clear it, contrary to the prediction.
- **§8's `error_ratio p99 = 35.6` is the pass count, not a mechanism.** 82.0% of flagged pixels
  clear by rung 3 — the shipped `refine_max_passes` — so ~18% survive, which is that tail exactly.

**BOTH OF THOSE ARE CHARACTERISATION, NOT REMEDIES, AND MUST NOT BE READ AS THE FIX.** A global
`eta/256` pays **256x everywhere** for a failure that is local, and a fourth refinement pass is
`refine_flagged` again — re-integration from `t = 0`, which a live playhead cannot do. Their value
is diagnostic and it is large: *`eta/256` brings every flagged pixel to `error_ratio` 1.000* is
what proves this is ordinary under-resolution and not a wrong equation, a saturating cap, or a
threshold. That is why `eta/256` is used as the **ground truth** in §11's comparison rather than
as a candidate in it. **Only a per-step mechanism survives contact with marching.**
- **And the plumbing found a real advance-anyway defect, at a site nothing named.** One RK4 step
  advanced the physical clock by up to **`2.209e128`** against a sync interval of `0.4`. `1e128`
  is finite so the divergence guard passes; `s.t >= dt_left` is satisfied by 128 orders so the
  march records a clean landing; and `t += dt_left` then corrects the **clock** to the boundary
  while keeping the **state** from `1e128` time units later. The clamp cannot un-take the step,
  and `t` is clamped on both branches — so the overshoot was invisible in every recorded quantity
  until `dt_max` existed. An unbounded step with no acceptance test, not a cap.

---

## 1. Symptom

`results/closure/config_stability_stop0_uniform.png` — and the whole family of committed panels
— carries pale, magenta-speckled patches that interrupt otherwise continuous ribbon structure.
They read as a rendering or physics fault rather than as fractal structure, and they were reported
as such.

## 2. What is actually wrong

`EnsembleCfg::default()` sets **`refine_flagged: true`**. That pass re-integrates any pixel whose
`error_ratio` exceeds `refine_threshold = 10.0` at `eta * refine_eta_factor = eta/4`, up to
`refine_max_passes = 3` times. It is the remedy BRIEF §2.5 prescribes for a failure this project
has on record:

> `eta = 1e-2` is not sufficient above ~64x64, and the failure is a cliff. […] Re-integrate
> flagged pixels at finer `eta` rather than lowering `eta` globally — no scheduler needed.

**`examples/closure_render.rs:89` sets `refine_flagged: false`.** So do **62 files** under
`examples/`, including every render harness: `dtau_render`, `overshoot_render`,
`escape_gate_render`, `adaptive_render`, `slice_gallery`, `chart_gallery`, `glsl_refinement`,
`bivariate_colour`. `src/bin/prinq.rs:73` sets it false too.

`error_ratio` exists precisely to say *this pixel is not data*. Its healthy value is `1.0`. On the
committed panel its **p99 is 1.039e+10**, and on six of sixteen sampled regions its **median** is
between `9.8e5` and `4.9e7`. The flag has been firing the whole time and nothing acts on it.

## 3. Evidence

Full slice, same window, same settings, one field changed:

```text
                    refine OFF (committed)      refine ON
drift p50                  4.251e-07            2.560e-07
drift p99                  8.866e+06            1.635e-02      540,000,000x
drift max                  1.968e+12            6.738e-02
error_ratio p50                1.000               1.0000
error_ratio p99            1.039e+10               35.571
non-finite pixels             30,109                    0
escape fraction               0.0403               0.0067
pixels refined                     —               0.1114
```

**11.1% of the slice needs re-integration.** All 30,109 non-finite pixels are repairable. Median
energy drift is *not* the tell — it barely moves — which is why this survived: the damage is
entirely in the tail, and `error_ratio` is the statistic that sees it.

One box (`B10`, 256²) in isolation:

```text
                    refine OFF        refine ON
drift p50            2.094e+03        4.156e-04       5,000,000x
drift max            3.506e+15        5.965e-02
error_ratio p50      3.2943e+07          2.1340
spread_shape p50     3.905e-01        3.208e-04           1,218x
non-finite                  60                0
escape                  0.4645           0.0021
collision               0.4099           0.8596
wall clock               123 s           1216 s              10x
```

`spread_shape` saturating at `0.39` **is** the white: the ensemble copies have diverged to
garbage, the spread pins at the top of the lightness ramp, and the bivariate map paints it pale.

The outcome distribution is not merely recoloured — it is **reclassified**. Escape goes
`0.4645 -> 0.0021` and collision `0.4099 -> 0.8596` in that box, and `0.0403 -> 0.0067` over the
whole slice. Anything downstream that read the terminal class over these regions read the wrong
class.

## 4. Not every marked region is this bug

Sixteen regions were sampled at their own windows. Six are the fault; ten are sound.

```text
   box  err_ratio p50  drift p50   steps/time vs nominal  d_min p50  <r_coll   verdict
    P5      4.943e+07   4.30e+03        5.25             3.69e-03    0.758   BROKEN
   B10      4.541e+07   5.23e+03        8.06             3.09e-03    0.834   BROKEN
    B9      2.568e+06   3.84e+02        5.10             3.65e-03    0.795   BROKEN
    B2      2.523e+06   1.21e+02        5.39             4.23e-03    0.744   BROKEN
    B1      1.144e+06   5.56e+01        4.71             4.39e-03    0.647   BROKEN
    B8      9.811e+05   1.12e+02       11.25             3.59e-03    0.872   BROKEN
    B5           5616   6.83e-01        2.46             3.38e-03    0.839   marginal
    P6          109.8   1.15e-01        1.54             4.71e-03    0.609   marginal
    B4          1.051   4.14e-05        1.35             2.70e-02    0.080   sound
    P1          1.001   1.94e-05        1.38             6.81e-03    0.398   sound
 FRAME              1   4.25e-07        1.02             4.69e-03    0.594   baseline
```

The broken set runs at **4.7x-11.3x the nominal step rate** and still loses the energy — they are
working harder, not stepping coarser. Their `d_min` sits at `3.1e-3`-`4.4e-3` against
`r_coll = 5e-3`, with 65-87% below it: deep close approaches, exactly the cliff `refine_flagged`
exists for. `B4` and `P1` read `error_ratio = 1.001` — **their pale structure is real** and must
not be attributed to this bug.

## 5. Why the override exists — it was correct, and it spread

Introduced at **`c03fc85`, 2026-08-22 23:28**, in the same commit that *implemented* the
refinement pass. Its stated reason is still in `results/README.md:190` and is sound:

> The experiment examples run with `refine_flagged: false`, **deliberately**. Experiments A and B
> characterise the kernel whose behaviour motivated the second pass, and measuring the repaired
> kernel would hide the thing being measured. Precision comparisons also run it off: the pass is
> threshold-triggered on `error_ratio` and f32 and f64 flag different pixel sets, so with it on
> the comparison would be of pipelines rather than of arithmetic.
>
> **The `render-*.txt` runs have it ON**, and report the before and after drift maxima on the
> same line.

Both halves are right. Experiments must measure the unrepaired kernel; **renders must not.** That
last line is the invariant, and it is now false: `results/README.md` still asserts renders have it
on, while every render harness turns it off.

The spread is traceable. `refine_flagged: false` propagated from the experiment harnesses into the
render harnesses by copy, one file at a time, over six days:

```text
  c03fc85  08-22 23:28  introduced, correctly, in experiments and precision tests
  fa0c8c0  08-23 03:03  halton
  857dff4  08-23 03:16  pooled-parent
  badb308  ...          the refinement scheduler
  39f8c83  08-24 00:10  the vertical slice
  b869fd8  ...          the between-footprint arm
  9312223  ...          temporal accumulators
  71de13f  08-27 13:02  closure_render  <- the harness that made the panel in question
```

No commit message argues for it in a render. It was inherited, not decided.

**This is a convention that outlived its justification and was never re-examined at the boundary
it was supposed to stop at.** The rationale was written down; the invariant it implies was not
tested; nothing fires when a render harness copies the line.

## 6. What is not affected

- **`src/bin/prin.rs`** — the production renderer. It builds from `Config::default()`, whose
  `ens` is `EnsembleCfg::default()`, so refinement is **on**. It also prints
  `refined N of M pixels re-integrated` and `|dE/E| max … before refinement, … after`. The
  `--no-refine` flag exists to turn it off deliberately.
- **The experiment and precision-comparison harnesses.** Their override is correct and should
  stay.
- **`src/bin/xcheck.rs`** — the cross-check compares against the numpy reference, which has no
  refinement pass; turning it off there is required.

## 7. The remedy — a DISCREPANCY fix, not a defect fix, and still not applied

**Read §0 first.** None of the below repairs the artefact. It makes renders and production agree,
and it makes any future disagreement impossible to introduce silently. Those are worth having on
their own terms and they are not the same thing as a fix.

### 7.1 The architecture, which is the part that matters — DONE

`prin` and every harness constructed `EnsembleCfg` **independently**: production took
`::default()`, harnesses wrote their own struct literals. **There was no single source of truth**,
so "the production config" existed nowhere — it was whatever `default()` happened to return, and
**111 literal sites** could disagree with it. The failure was never that someone chose `false`; it
is that **nothing recorded the choice**, so it propagated by copy through five commits and six
days invisibly.

`src/ensemble/provenance.rs`:

- `EnsembleCfg::production()` is the one literal; `Default` delegates to it.
- `Override` is a **named** value per field, one variant per field, so
  `production().with_overrides(&[Override::RefineFlagged(false)])` declares what it changes.
- `overrides_vs_production()` **derives** the list by diffing, so a config declares itself
  *however it was built* — including all 111 existing literals, and including any future one whose
  author forgets. A hand-maintained list can itself go stale, which is the same class of failure
  this exists to stop; the derived diff cannot.
- Both the diff and `Override::apply` destructure and match **exhaustively, no `..`, no `_` arm**.
  Adding a field to `EnsembleCfg` breaks the build until it is handled. A mechanism reporting "no
  overrides" because it does not know about a field is the original defect, one level up.
- `provenance()` renders it for a header, and reads `production` — not an empty string — when
  there are none. A blank field and an absent field look the same in a log.
- `output::provenance_sidecar` writes `<stem>.cfg.txt` beside a panel. **This is where the six
  days lived**: the `.raw` and `.prnq` dumps have carried a full settings header since they were
  written; the PNGs carried nothing and the harnesses printed nothing.

It works as intended on first use — `cliff_ladder`'s own header now reads:

```text
config: production + 4 override(s): t_max=50.0 (production 13.0), n_sync=125 (production 32),
        r_coll_frac=0.005 (production 0.001), refine_flagged=false (production true)
```

`tests/provenance.rs` holds five properties, including the one that matters most: **a plain struct
literal still declares itself**, because the 111 legacy sites were never going to be rewritten in
one go and a mechanism that only declares configs someone remembered to annotate would have missed
every one of them.

### 7.2 Still to do

1. Render and gallery harnesses take `production()` with no `RefineFlagged` override. The
   experiment and precision harnesses (`c03fc85`'s original set) and `xcheck` **keep** theirs — it
   is correct there. `closure_render` and `saturation_mask` now print and sidecar their config;
   the other 25 PNG-writing harnesses in §8's table do not yet.
2. Re-render the committed corpus, or label it. Cost is ~3x over this slice, ~10x on a badly
   affected window.
3. Correct `results/README.md:190`, which still asserts renders have refinement on.
4. **`refine_max_passes = 3` stops one rung early** (§0). Raising it is a cost decision, not a
   correctness one, and it belongs to whoever owns the render budget.

## 8. What this does NOT explain, and what the reclassification touches

- **`error_ratio` p99 is 35.6 after three passes** — **explained** in §0. 82.0% of flagged pixels
  clear by rung 3, so ~18% survive the shipped ladder. Not a floor: the ladder converges
  completely by rung 4, and `NEVER cleared` is `0 of 128`.
- **Ten of the sixteen marked regions are sound and their structure is real. The pale wedge in
  `B4` survives refinement and is not this bug.** This must stay prominent or the next
  investigation starts from a false premise. §0's independent check qualifies it rather than
  confirming it: `B4` reads `err>10` at **0.3163** against a frame baseline of **0.1111** — three
  times the frame, but far below the broken set's 0.60-0.98. "Sound" means *not this fault*, not
  *unflagged*.
- **The escape fraction moving `0.0403 -> 0.0067`** is a six-fold change in terminal
  classification, not merely in colouring. Anything read off the outcome panel for this slice
  needs re-checking.

### The corpus in scope

**21 harnesses** carry `refine_flagged: false` **and** write PNGs; **19 of those read terminal
class, `t_end` or the event class**. `spread_event` reads the event class, so the criterion corpus
is in scope. This lists the harnesses, not the individual committed files — a file-level audit
needs a re-render to be worth anything.

| harness | writes into `results/` | reads terminal |
|---|---|---|
| `adaptive_render` | yes | class |
| `banding_render` | yes | class,t_end |
| `between_vs_within` | yes | class |
| `bivariate_colour` | yes | - |
| `box_panels` | no | class,t_end |
| `chart_gallery` | yes | class,t_end,event |
| `closure_render` | yes | class,t_end |
| `colour_check` | yes | class,t_end,event |
| `criterion_metric` | yes | class,t_end |
| `dtau_render` | yes | class |
| `escape_gate_render` | yes | class,t_end |
| `overshoot_render` | yes | class |
| `pan_sequence` | yes | class |
| `preset_control` | yes | class,event |
| `saturation_mask` | yes | class |
| `slice_gallery` | yes | class |
| `slice_refined` | yes | class |
| `ssaa_resolve` | yes | class |
| `switch_study` | no | class |
| `wedge_zoom` | no | class |
| `zoom_sequence` | yes | - |

Not affected: `src/bin/prin.rs` (takes the default), the experiment and precision harnesses (their
override is correct), and `src/bin/xcheck.rs` (the numpy reference has no refinement pass).

## 9. Reproduce

```sh
# the saturation check and the eta ladder -- results/saturation/README.md
cargo run --release --example saturation_mask 512 results     # 1087 s
cargo run --release --example cliff_ladder   256 5 128 results #  336 s

# the refinement A/B -- results/refine_bug/
cargo run --release --example slice_refined 384 on  <out>
cargo run --release --example slice_refined 384 off <out>
cargo run --release --example refine_ab 192 0.942 0.789 0.0522 B10 <out>
```

Every harness above prints its config as a `provenance` line and, where it writes a panel, a
`<stem>.cfg.txt` sidecar beside it.

## 10. Tests holding this

- `tests/provenance.rs` — five properties of the single-source-of-truth mechanism.
- `tests/saturation_plumbing.rs` — `dt_max`, `ab_min` and `n_cap_hits` reach the payload and carry
  values that could not be defaults (an unplumbed `dt_max` is exactly `0.0`, an unplumbed `ab_min`
  is `INFINITY`), plus a tame-region negative control.
- `src/integrate/az/driver.rs::step_control_tests` — **the `T::TINY` floor fires at f64 and the
  march advances anyway**, with a healthy-state negative control and a mode control for the cap.


---

## 11. THE FIX — a predictive per-step limit, `results/step_control/README.md`

Four candidates behind `StepLimit`, measured rather than argued. **B wins outright.**

```text
   config_stability, 192^2      steps p50      err p99   err>10   overshoot
   None (baseline)                1.033e5      7.108e9   0.1110         634
   Predictive f=0.02              1.053e5        1.109   0.0000           0   <- +1.9%
   Reject     f=0.02              1.832e5        1.205   0.0000           0   <- +77%
   AbGrowth   f=2                 1.033e5      7.108e9   0.1110         634   <- bitwise inert
   Global     f=0.25              4.134e5      9.066e7   0.0767         153   <- +300%, FAILS
```

`dtau <= f*d_min/(|v_rel|_max*A*B)` — one divide, no trial step, no retry, no branch, from values
`phys_from_state` already returns. Shipped as production at `f = 0.02`.

- **The dumb control does not fix it at four times the cost.** `Global` still leaves 153
  overshoots. A uniform `eta` cut cannot bound a step whose size is set by local geometry.
- **A is not GPU-viable**: at the parameter it needs, **every warp contains a retrying lane**
  (1.0000, both dispatch shapes; worst lane 5.2M retries), and it **plateaus above 1.0** on
  `preset_shape` where 39 of 96 pixels exhaust the retry budget.
- **C was already shipped** — `DtauMode::PerStepInterval` is an `A*B` growth clamp at `C = 1`.
- **The cap can now be removed** and that is reported, not done: a second corpus-invalidating
  change belongs in its own measurement.
- **Read `steps`, not `secs`** — under load 85–100 the winning row timed faster than the baseline
  while doing more work.

`cargo test --release` **256 passed, 0 failed**; xcheck **4/4** (`reference_opts` pins `None`).
Three tests failed when the default changed and **every one failed correctly** — the limit deletes
the damaged population those characterisation tests are about. They are pinned to `StepLimit::None`
with that reason recorded.

**The committed corpus was taken under `None` and does not reproduce bitwise under this default.**
Stated rather than discovered: `provenance()` names the setting in every header and sidecar.
