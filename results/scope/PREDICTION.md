# Predictions, written before the full `corpus_scope` run

Recorded at `c63b536`, before `examples/corpus_scope.rs` was run at anything past a 16² / budget-100
smoke on `near-field`. The smoke is quoted where it already bears on a prediction, and labelled.

**The question.** The scheduler corpus is 25-26 August; every integrator fix is 27 August or later.
Is the re-measurement a full re-derivation of the criterion work, or a re-run of its numbers under
intact conclusions? The discriminator is `rho` — Spearman of `ensemble_spread` between the kernels
— because the criterion is a **ranking**, and *a ranking is invariant to a monotone rescaling of
the signal; a threshold is not.*

## 1. `rho` splits by region, and `far` is near 1.0

`far` is on record as featureless — thirteen non-greedy criteria produce the **identical leaf set**
there, `err_sum` is flat through level 3, and ranking by spread on a smooth field *is*
breadth-first. A smooth gradient's ordering is a property of the geometry, not of the integrator,
so any kernel should reproduce it. Predicted `rho > 0.95` in `far`.

`near-field`, `deep interior` and `body2 core` are chaotic and the ordering is set by where
trajectories land. Predicted `rho` **well below** 0.9. *(The 16² smoke reads 0.7305 on `near-field`
— consistent, and not yet at the resolution that ships.)*

## 2. `deep interior` saturates and its chord column says nothing

On record: `chord max = 1.999`, antipodal, at every rung of an `eta` ladder over `t = 13` — chaotic
divergence no step size buys off. Predicted `antip` near 1.0 there, and near 0 in `far`. If it
saturates, `rho` and stage 2 are the only columns that mean anything in that row, and the harness
prints `antip` so the cell can be recognised rather than quoted.

## 3. Stage 2 churn tracks `rho`, and `far` is near zero

Predicted `differ/shared` above 25% in the chaotic regions and under 5% in `far`.

## 4. `Floor` engagement falls under the new kernel

`Decision::Floor` is *"the branch that must be shown to engage"*. It fires when `alpha` says
refining will not reduce the spread — and a spread inflated by integration failure is exactly what
looks irreducible. So the corpus's `Floor` counts are predicted to be **partly a kernel artefact**.
*(The smoke reads `floor:17` under `pre` and `floor:0` under `now`, which is why this is written
down as a prediction rather than discovered afterwards.)*

**This is the one that would change a conclusion rather than a number**, because the `alpha` gate
is the branch `preset_shape` is on record as the only tree exercising.

## 5. The `pre` arm carries non-finite pixels the `now` arm does not

`StepLimit::None` is the arm that produced `err>10 = 0.1110` on `config_stability` and 634
overshoots. Predicted `nonfin pre > nonfin now` in at least the chaotic regions, and this is the
first run with `n_undetermined` available to see it at quad level.

## What would refute the whole framing

`rho > 0.95` in **every** region. Then the field was monotonically rescaled, every ordering-based
conclusion in the corpus stands, and the re-measurement is a re-run of numbers under intact
headings. I do not expect it, and it is written here so that outcome is reportable rather than
explained away.

## What this run does NOT decide

Whether the criterion is any *good* under the new kernel. It measures how far the input moved, not
what the right criterion is. `error(B)` against `dp_optimal` is a separate run and the corpus's
arithmetic findings survive to inform it.

---

# Outcomes, appended after the run

**Two confirmed, three refuted.** Full record in `README.md`; the short form:

1. **REFUTED.** `far` reads `rho = 0.5887`, the **lowest** in the table, not `> 0.95`. The reason
   is more useful than the prediction: `far`'s bulk is flat to 1.4% (p01-p90 `1.500e-8` to
   `1.521e-8`), so there is **no ordering to preserve** and `rho` measures noise reshuffling. The
   reasoning behind the prediction — a smooth field's ordering is geometric — was right about the
   geometry and wrong that there is an ordering.
2. **REFUTED.** `deep interior` reads `antip = 0.0034`, not near 1.0. The saturation result on
   record is over `t = 50` at an `eta` ladder and does not transfer to `t = 13` between kernels:
   `deep interior` terminates by **collision**, and a terminated trajectory's shape is frozen
   rather than diverging.
3. **Confirmed, with the caveat now measured.** Churn is 36% / 29% / 24% in the three regions that
   decide anything. `far`'s 0% is **not** stability — it and `mid-field` read 16 leaves, `keep:16`
   on both arms, the standing `tau`-below-the-bulk degeneracy. Nothing there decides.
4. **CONFIRMED, and it is the headline.** `Floor` collapses **17 -> 0, 16 -> 1, 9 -> 1**.
5. **REFUTED as stated, and this is the instructive one.** `nonfin` reads 0/4096 on both arms in
   every region. But the `pre` arm *does* carry undetermined footprints — `deep interior` has
   **11 with a `NaN` `spread_shape` against 0 under `now`** — and they are invisible in that column
   because `ensemble_spread` is `sp_shape.max(sp_event)` and **`f64::max` ignores `NaN`**. The
   physics was predicted correctly and the statistic chosen could not see it. *Ask what a
   diagnostic would say about the defect you are hunting before reading it as clean.*

**The refutation clause did not trigger.** `rho > 0.95` in every region would have meant the field
was monotonically rescaled and every ordering-based conclusion stood. `rho` is 0.59-0.84.

**Unpredicted, and the largest single correction:** `far` gains a **two-decade tail above p90**
(p99 `1.917e-8 -> 1.926e-6`, max `1.921e-8 -> 4.403e-6`) while its bulk stays flat. The control
region's "featureless" result was measured on a field with no tail at all.
