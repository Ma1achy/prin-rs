# The criterion work, rebuilt: items A, B, D, E

`cargo run --release --example criterion_metric -- 6 8 1e-4 13 results production`
-> `results/output/criterion_metric.txt`, `results/criterion/*_t13.{qcache,fcache}`, panels.
Complete tree to level 6, `N = 8`, `E+1 = 8`, 512^2 reference, 5461 quads and 2,796,032
trajectories per region. `config: production + 2 override(s): refine_flagged=false,
keep_boundary_shapes=true` -- printed by the run, derived by diffing, not asserted by hand.

**Item C is not here because it already existed**: `signal_audit` ranks by the **blocked** rho and
prints the pooled one beside it. That was a methodology fix, not a measurement, and it survives.

---

## B. Does any criterion beat breadth-first? **It is a cliff, not a verdict.**

`near-field`, the share of achievable improvement and the scale it is normalised by:

```
B =                 11      23      47      95     191  |    383     767    1535    3071
headroom/err    3.8e-3  2.8e-2  5.1e-2  6.7e-2  9.1e-2  | 5.8e-4  9.5e-4  1.8e-3  3.5e-3
term_grad        1.000   0.958   0.610  -0.194  -1.352  |   ---  magnified noise  ---
contrast:between 0.003   0.666   0.544  -0.039  -1.146  |
perim_within     1.000   0.674   0.333  -0.206  -0.917  |
within/mean      1.000   0.800   0.334  -0.204  -1.181  |
frac_hot_between 0.000   0.000   0.000   0.000   0.000  |
greedy           1.000   1.000   0.977   0.909  -0.117  |
```

Headroom rises to **9.1% of the error at `B = 191`**, then **collapses 157-fold** at `B = 383` and
stays near 0.1%. Two regimes:

- **Below the cliff**, several criteria take most of what is available -- `term_grad/median`
  **0.958** at `B = 23` and **0.610** at 47, `contrast:between` 0.666/0.544, `perim_within`
  0.674/0.333.
- **Above it**, uniform is within **0.1%** of the exact optimum. There is nothing to beat, and
  every `captured` figure there (-89 to -621) is a magnified view of a small absolute gap.

**So the standing "nothing beats breadth-first" was true and uninformative.** It was quoted at
large budgets, where it holds because there is no headroom rather than because the criteria fail.
The interesting regime is the one nobody was reading.

**`headroom/err` is why this is legible, and it is new.** Without it the `B >= 383` block reads as
a catastrophic criterion failure instead of an absent denominator. The standing rule covers a
**zero** denominator; this is the finer case of a **small** one, and printing the scale beats
hiding the row.

## The three regions do not share a shape

```
                 headroom/err range      who wins
  near-field     3.8e-3 .. 9.1e-2        criteria below B=191, uniform above
  far            1.2e-13 .. 1.7e-2       uniform, always; greedy 2.1x WORSE
  deep interior  2.9e-3 .. 8.3e-2        greedy exactly optimal, every criterion far behind
```

**`deep interior`: `greedy_lookahead_1` reads captured 1.0000 at every budget**, and every
signal-based criterion is -25 to -236. The record's *"`deep interior` at `B >= 6143` shows a
criterion decisively ahead"* does not reproduce anywhere in the measurable range -- though that
figure was taken at `levels = 7`, where `B = 6143` exists and here the tree tops out at 5461, so
this narrows the claim rather than refuting it.

**`far` reproduces a standing result intact**: `greedy_lookahead_1` is **2.1x worse than uniform**
(0.37974 against 0.17968 at `B = 1535`), the *greedy is beaten by an arbitrary scan order where
the field is featureless* result, unchanged by the kernel. And `far`'s **AUTO-RANGED OVER NOISE**
guard fires (`ramp span x8.078`), so the whole region's curve sits on a stretched noise floor --
the standing caveat, printing itself.

## D. The cost axis, in substeps -- and it CHANGES THE VERDICT

`deep interior`, the same replays read against substeps instead of quads:

```
  substeps =    2.97e8    6.99e8    1.65e9    3.87e9
  uniform      0.14686   0.10919   0.07679   0.04191
  greedy       0.14623   0.10724   0.07081   0.04556
  greedy/cost  0.14231   0.10622   0.07010   0.03834   <- best at every rung
```

`greedy_lookahead_1/cost` optimises `Δerror / substeps` and **wins at every large rung on the axis
it optimises**. On the quad axis it reads an erratic 0.24 / 0.06 / 0.66 / -0.46 / 0.83 / -1.61 /
0.13 -- mediocre, and it had never been plotted any other way.

**And the axis flips plain greedy too**: captured **1.0000 at every quad budget**, yet *worse than
uniform* on substeps at the largest rung (0.04556 against 0.04191). A quad budget is only a cost
when quads cost the same, and they do not. This is not a refinement of the old axis; it is a
different question, and one of the two rankings built to answer it was being scored in the wrong
units.

## E. `lay_w_perimeter`, plotted for the first time

`perimeter_ratio` alone was never a `Criterion`, so `signal_audit`'s `lay_w_perimeter` -- 0.07605
against uniform 0.08480 and dp 0.07025 at `B = 6143`, the pre-fix figure that made it the lead --
could never go through the one measurement that decides a criterion here.

Now it can. `Criterion::PerimeterWithin` captures **0.674 at `B = 23`** and **0.333 at 47** in
`near-field`: respectable, mid-pack among the early winners, **behind `term_grad` at 0.958/0.610**.
Not the standout the pre-fix number implied. In `far` it reads **exactly 0.0000 at every budget**
-- its mask is empty on every leaf there, so it is `NaN` throughout, never wins a comparison, and
falls through to the level-first tie-break, which *is* breadth-first. Documented behaviour,
confirmed.

## The `frac_hot_between` result now has a MECHANISM, not a number

It reads captured **exactly 0.0000 at every budget in `near-field` (B = 23..767) and in `far`
(every rung)**: it reproduces the uniform leaf set **bit for bit**.

The cause is arithmetic. A fraction over `N^2 = 64` footprints has **at most 65 distinct values**,
and in `far` it has **1** (100% modal). Most quads tie, ties fall to the level-first tie-break, and
that tie-break is breadth-first. **The best criterion measured on this project matched
breadth-first by being it.** The distinct-value ceiling was already on record as a caution about
reading it as a signal property; this is the same fact showing up as behaviour.

`perim_within` and `frac_hot_within` do the same thing in `far`, for the same reason.

## The cliff was a LADDER ARTEFACT, and the real structure is a COMB

**Withdrawn: "the cliff at `B = 383`".** The budget ladder is `5, 11, 23, 47, ...` -- that is
`4k+1`, which lands one split **past** a complete tree every time and never **on** one. A complete
tree to level `L` holds `(4^(L+1)-1)/3` quads: **21, 85, 341, 1365, 5461**, and not one of them
was ever a rung. So when all three regions put their headroom minimum at `B = 383`, they were
agreeing about the ladder's nearest approach to **341**, not about the field. *Which rung of a
sweep is degenerate is a fact about the region* -- and which rung is **missing** is a fact about
the ladder.

With the complete-level counts merged in, `headroom/err`:

```
B =              1     5     11      21      23      47      85      95     191     341     383     767    1365    1535    3071
near-field       0     0  3.8e-3  2.8e-2  2.8e-2  5.1e-2  5.9e-2  6.7e-2  9.1e-2  4.7e-4  5.8e-4  9.5e-4 -6.9e-16 1.8e-3  3.5e-3
far              0     0 1.2e-13  1.7e-2  1.7e-2  1.3e-2  3.3e-3  2.7e-3  2.8e-3 -1.5e-16 4.3e-4  2.1e-3 -4.2e-15 4.9e-3  8.6e-3
deep interior    0     0     0    1.2e-16 1.2e-16 1.1e-2    0     2.9e-3  1.7e-2  2.2e-3  8.8e-3  5.2e-2  8.8e-5  1.2e-2  8.3e-2
```

**At a budget that buys a complete tree, breadth-first IS the exact optimum**, bit for bit:

```
                B =    21        85       341      1365
  near-field  dp    0.30018   0.28283   0.20081   0.10026
              uni   0.30018   0.28283   0.20090   0.10026
  deep int.   dp    0.23890   0.17346   0.11399   0.06129
              uni   0.23890   0.17346   0.11424   0.06129
```

That is **not** a tautology. `Dp::labels` optimises over *all* tree-shaped leaf sets at that count,
and an unbalanced tree of the same size was free to win. It does not: on this field, at a budget
affording a complete level, no allocation beats the flat one.

**So the structure is a comb, not a cliff.** Headroom returns to machine zero at every complete
level and rises between them, and the reframing of item B is exact:

> **A criterion can only ever win between complete levels.** How much is available is set by how
> far the budget sits from the next complete tree, not by the region being "hard".

The tooth heights are regional and their trend with depth differs: `near-field`'s shrink
(9.1e-2 at `B = 191`, then 9.5e-4 in the next interval), `deep interior`'s **grow**
(1.7e-2, then 5.2e-2, then 8.3e-2). That is the real regional difference, and the earlier reading
-- *"uniform is optimal above `B = 383` in near-field"* -- is wrong: uniform is optimal **at**
complete levels and near-optimal between them **in near-field only**, while `deep interior` opens
up with depth.

## What this does not answer

Whether the teeth would keep shrinking in `near-field` past level 5, which needs `levels = 7`.
**That run is blocked by memory and the number is worth recording**: `keep_boundary_shapes` holds
`quads x 64 x 8 x 33 x 3 f64` -- **2.2 GB at `levels = 6` and 8.9 GB at `levels = 7`** -- and it
OOMs on an 18 GB machine. Reaching level 7 means turning it off, which drops `running_max`,
`first_div` and `term_grad` from that table, including the best early performer. A deliberate
trade, not a retry.

And `deep interior`'s `.qcache` was, briefly, the only stale one in `results/criterion/` after an
interrupted run regenerated its two siblings. The whole run was relaunched rather than the gap
filled, so all three caches come from one run of one binary. **An interrupted regeneration
produces the mixed-version corpus live**, which is the failure the corpus write-up describes,
arriving by a route it did not anticipate.
