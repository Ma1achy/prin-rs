# BUG REPORT — the render harnesses disable the repair pass, and the committed panels carry the fault

**Status:** confirmed, reproduced, cause located, not yet fixed.
**Found:** 2026-08-29, while diagnosing pale/magenta patches on `config_stability`.
**Severity:** every committed render made by an `examples/*_render.rs` harness is affected.
The production binary `prin` is **not** affected.

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

## 7. Proposed fix — not applied

1. **Restore the invariant.** Every `examples/*render*.rs` and gallery harness takes
   `refine_flagged` from `EnsembleCfg::default()`, or exposes it as an argument defaulting to on.
2. **Guard it**, on the model of `scheduler::assert_not_uniform_in_disguise`: a render that writes
   into `results/` with `refine_flagged: false` and a non-trivial `error_ratio` tail should refuse.
   *A configuration that silently reproduces the old behaviour needs a guard, not a convention* is
   already a standing rule in this project — this is the same defect at a second site.
3. **Print the flag with every render**, as `prin` already does. None of the affected harnesses
   record `refine_flagged` in their stdout or in the PNG's companion text, which is why six days
   of renders carry it invisibly.
4. **Re-render the committed corpus**, or label it. The cost is real: ~10x on a badly affected
   window, ~3x over this slice.
5. **Correct `results/README.md:190`**, which currently asserts the opposite of the truth.

## 8. What this does NOT explain

- `error_ratio` p99 is **35.6** after refinement, not 1.0. A tail is still unresolved after three
  passes at `eta/4`. The cliff is steeper than the current `refine_max_passes` reaches.
- Ten of the sixteen marked regions are sound and their structure is real. The pale wedge in `B4`
  survives refinement and is not this bug.
- The escape fraction moving `0.0403 -> 0.0067` over the slice means the terminal classification
  was substantially wrong before, not merely the colouring. Any result read off the outcome panel
  for this slice needs re-checking.

## 9. Reproduce

```sh
cargo run --release --example slice_refined 384 on  <out>   # refinement on
cargo run --release --example slice_refined 384 off <out>   # as committed
cargo run --release --example refine_ab 256 0.942 0.789 0.0522 B10 <out>   # one box, both arms
cargo run --release --example box_report 128 all off        # the sixteen-region table
```
