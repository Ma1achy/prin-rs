# How far did the criterion's input field move?

`examples/corpus_scope.rs` -> `corpus_scope.txt`. Predictions in `PREDICTION.md`, written before
the run. Reproduce with `cargo run --release --example corpus_scope -- 64 2000` (about a minute).

**The question.** The scheduler corpus is 25-26 August; every integrator fix is 27 August or later.
Is the re-measurement a **full re-derivation** of the criterion work, or a re-run of its numbers
under conclusions that still stand?

**The answer: a full re-derivation.** Four independent lines, below.

## The two kernels

`pre` reconstructs the corpus kernel through the flags rather than by checking out an old commit —
the method already validated here when `FixedPerInterval` reproduced `71de13f` bitwise.

Built from a **field-by-field diff** of `EnsembleCfg` at `0114be4` (26 Aug 15:06) against
`production()`: **seventeen fields did not exist then.** Naming the four famous fixes and stopping
would have left thirteen at today's values inside an arm labelled *pre*.

| field | corpus era | today | `pre` |
|---|---|---|---|
| `integrator` | AZ only | `Heggie` | `Az` |
| `dtau_mode` | fixed per interval | `PerStepInterval` | `FixedPerInterval` |
| `clamp_final_step` | overshoot | `true` | `false` |
| `step_limit` | none | `Predictive` | `None` |
| `land_iterate` | none | `true` | `false` |
| **`escape_rule`** | reference, hardcoded | `Closure` | **`Reference`** |
| `escape_confirm` | none | `true` | `false` |

`refine_flagged` is off in **both** arms: it is batch-only, has no live-playhead analogue, and it
removes the very population a kernel comparison is about.

**The audit has its own control, and it fired on exactly one region.** A third arm rebuilds `pre`
from the four famous fixes alone:

```
    near-field      differs on     0/4096    the audit was INERT -- four flags would have been enough
    mid-field       differs on     0/4096    INERT
    far             differs on     0/4096    INERT
    body2 core      differs on     0/4096    INERT
    deep interior   differs on   361/4096    the audit MATTERED: 361 outcome labels differ
```

`deep interior` is one of the three regions the corpus actually covers. A `pre` arm built from
memory would have carried a **post-corpus terminal classifier** into it, and nothing would have
said so.

## Stage 1 — the field

`rho` is the headline because the criterion is a **ranking**: *a ranking is invariant to a monotone
rescaling of the signal; a threshold is not.* A rescale would leave `rho` at 1.0.

```
region          distinct   moved    spr pre    spr now  rat p10  rat p50  rat p90       rho chd p50   antip
near-field     4096/4096  1.0000   1.525e-3   2.677e-4    0.099    0.179    0.354    0.6359   0.037  0.0000
mid-field      4096/4096  1.0000   9.399e-8   7.008e-8    0.736    0.746    0.756    0.8351   0.000  0.0000
far            4096/4096  1.0000   1.892e-8   1.512e-8    0.795    0.798    0.802    0.5887   0.000  0.0000
body2 core     4096/4096  1.0000   1.894e-3   5.998e-4    0.156    0.374    0.681    0.6456   0.038  0.0000
deep interior  4096/4096  0.9998   5.808e-5   8.478e-6    0.000    0.146    0.409    0.7281   0.000  0.0034
```

**`rho` runs 0.59-0.84 in every region.** Nowhere near 1.0, so this is not a rescale: the field was
re-ordered, and a re-ordering is what a criterion comparison is made of.

`antip` is 0.0000-0.0034 throughout, so the chord column is **not saturated** and can be read — the
standing `chord max = 1.999` result is over `t = 50` at an `eta` ladder, and does not transfer to
`t = 13` between kernels. `deep interior` is dominated by *collisions*, and a terminated trajectory's
shape is frozen rather than diverging.

## `far` is flat in the bulk under BOTH kernels, and the new one adds a two-decade tail

`far`'s cell is the surprise, and it needed the standing rule — *count the signal's distinct values
before reading any curve* — to be read at all. Its ratio band is **0.795-0.802**, a rescale to
within 0.9%, sitting beside `rho = 0.5887`, the lowest in the table. A near-perfect monotone map
cannot move a rank, so those two numbers contradict each other until the ladder is printed:

```
  ensemble_spread     p01        p10        p50        p90        p99        max
    far  pre     1.868e-8   1.876e-8   1.892e-8   1.908e-8   1.917e-8   1.921e-8
    far  now     1.500e-8   1.502e-8   1.512e-8   1.521e-8   1.926e-6   4.403e-6
```

The **bulk is flat under both** — p01 to p90 spans 2.8% (pre) and 1.4% (now). So within the bulk
there is **no ordering to preserve**, the rank is set by the last significant digits, and `rho`
there is measuring noise reshuffling rather than a change. *Two different faults give the same flat
curve: a bad ordering and no ordering* — and this is the second.

But **above p90 the new kernel carries a tail the old one does not**: p99 `1.917e-8 -> 1.926e-6`,
max `1.921e-8 -> 4.403e-6`. Two decades, in 5-10% of the region. `far` is the control region — the
one on record as *what a featureless field looks like*, where thirteen non-greedy criteria produce
the **identical leaf set**. That result was measured on a field with **no tail at all**, and a tail
is exactly what a criterion ranks on. Not overturned; measured on a different field, and it needs
re-taking.

## Stage 2 — the decisions, and `Floor` collapses

The real descent under both kernels, at the corpus's own `SchedCfg` (`n=8`, `bootstrap=2`,
`tau=1e-4`, `hot=q[0.50]`, `agg=median`, `criterion=within`, camera viewport 1024, `max_rel_depth`
6). Decisions compared over **shared quads only, with the count printed**.

```
region          lv pre  lv now  shared   differ        stops pre                 stops now
near-field          46      64      33       12   floor:17 keep:21 mrd:8     keep:48 mrd:16
mid-field           16      16      21        0             keep:16                  keep:16
far                 16      16      21        0             keep:16                  keep:16
body2 core          43      64      41       12   floor:16 keep:23 mrd:4    floor:1 keep:47 mrd:16
deep interior       22      16      21        5             floor:9 keep:13   floor:1 keep:15
```

**`Decision::Floor` collapses: 17 -> 0, 16 -> 1, 9 -> 1.** `Floor` is documented in `quad.rs` as
*"the branch that must be shown to engage"*, and it fires when `alpha < alpha_lo` — when refining
did **not** reduce the spread. A spread inflated by integration failure is precisely what looks
irreducible, because the failure is a property of the trajectory and not of the cell width, so
halving the cell does not halve it. Under the fixed kernel it essentially stops firing.

That is a **branch of the decision changing behaviour, not a value changing magnitude**, and it is
the single strongest reason the re-measurement is a re-derivation.

Churn on shared quads is **12/33, 12/41, 5/21** — 36%, 29%, 24%.

**`far` and `mid-field` decide nothing and their `differ = 0` is not evidence of stability.** Both
read 16 leaves, `keep:16`, on both arms: the standing `tau`-below-the-bulk degeneracy, *where the
bulk sits below `tau` everything keeps and the tree is uniform at depth 2 — 16 leaves against a
complete 4096*. Of five regions only three reach a decision at all, and all three moved.

## Predictions: two confirmed, three refuted

| # | prediction | outcome |
|---|---|---|
| 1 | `rho > 0.95` in `far`, low elsewhere | **REFUTED** — `far` is the *lowest* at 0.5887 |
| 2 | `deep interior` chord saturates | **REFUTED** — `antip` 0.0034, not 1.0 |
| 3 | churn > 25% chaotic, < 5% in `far` | **confirmed**, but `far`'s 0% is degeneracy not stability |
| 4 | `Floor` engagement falls | **CONFIRMED** — 17 -> 0, 16 -> 1, 9 -> 1 |
| 5 | `pre` carries non-finite where `now` does not | **REFUTED as stated** — 0/4096 both arms, everywhere |

**Prediction 1 was refuted for a reason that turned out to matter more than the prediction.** I
reasoned that a smooth field's ordering is geometric and any kernel reproduces it. The ordering
is geometric; there just is not one, because the bulk is flat to 1.4%.

**Prediction 5 is the instructive one: the physics was right and the statistic was blind to it.**
The `pre` arm *does* carry undetermined footprints — `deep interior` has **11 with a `NaN`
`spread_shape` against 0 under `now`**. They do not appear in the `nonfin` column because
`ensemble_spread` is `sp_shape.max(sp_event)` and **Rust's `f64::max` ignores `NaN`**, so all 11
report their *event* spread as an ordinary number. I chose the one field that cannot see the thing
I predicted. *A diagnostic field is specific to a class of defect* — ask what it would say about the
defect you are hunting **before** reading it as clean.

That finding closed a gap in this PR's own guard: `footprint_undetermined` caught all 11, but
through its `n_nonfinite` arm, by coincidence — every one of them also had an unusable copy. A
triple collision reaches the same state with every copy usable. The predicate now tests
`spread_shape` directly and `tests/no_discard.rs` asserts it on a footprint constructed to have
exactly that shape. **The `f64::max` swallowing is a defect in `pixel.rs` and is deliberately not
repaired here** — propagating the `NaN` changes `ensemble_spread` itself, which moves every tree
and every render, and it wants its own attribution.

## What this run does not decide

Whether the criterion is any **good** under the new kernel. It measures how far the input moved,
not what the right criterion is. `error(B)` against `dp_optimal` is the next run, and the corpus's
*arithmetic* findings survive to inform it.

It is also a **scoping instrument, not a one-off**: it costs about a minute and should be re-run
whenever the kernel moves again.
