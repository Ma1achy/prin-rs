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

## THE COMB IS WITHDRAWN, AND THE ERROR WAS A TRANSCRIPTION SLIP OF MINE

*"At a budget that buys a complete tree, breadth-first IS the exact optimum, bit for bit"* rested
on a four-column table I copied out of the levels-6 output. The committed file says:

```
B =            1        5       11       21       23       47       85       95      191      341     1365
dp_optimal  0.30245  0.30130  0.30017  0.29178  0.29178  0.28183  0.26605  0.26290  0.22953  0.20081  0.10026
uniform     0.30245  0.30130  0.30131  0.30018  0.30018  0.29698  0.28283  0.28185  0.25246  0.20090  0.10026
```

What I wrote down as `dp` at `B = 21, 85` -- `0.30018` and `0.28283` -- are **uniform's** values.
Real dp is `0.29178` and `0.26605`, gaps of **2.8% and 5.9%**. The `headroom/err` row printed
`2.80e-2` and `5.93e-2` four lines below the row I was transcribing, so the harness stated the
contradiction and I did not cross-read it. *The measurement was correct and the column chosen
could not see it* has a blunter sibling: **the column was right there and I copied the wrong one.**

Scope: the slip is **near-field only**; `deep interior`'s half of that table was accurate, and it
really does have `dp == uniform` at `B = 21, 85, 1365`. But near-field was half the evidence, and
the generalisation across regions was drawn from the pair.

**What actually holds** is region-specific and much weaker: each region has *some* budgets where
uniform is exactly optimal -- `deep interior` at 21, 85, 1365; `far` at 341, 1365; near-field at
1365 alone -- and they are not the complete levels in general.

## AND THE MACHINE ZEROS MOVE WITH THE DEPTH CAP, SO THEY ARE NOT A FIELD PROPERTY

The ladder is `4k+1`, so quadrupling `levels` maps every levels-6 rung onto a levels-7 rung almost
exactly (`191 -> 764 ~ 767`). Laid out that way, near-field:

```
  lv6  B      11      21      23      47      85      95     191     341     383     767    1365    1535    3071
       h/e 3.8e-3 2.8e-2 2.8e-2 5.1e-2 5.9e-2 6.7e-2 9.1e-2 4.7e-4 5.8e-4 9.5e-4 -6.9e-16 1.8e-3 3.5e-3
  lv7  B      47      85      95     191     341     383     767    1365    1535    3071    5461    6143   12287
       h/e 2.6e-2 3.5e-2 3.8e-2 5.6e-2 7.0e-2 7.7e-2 9.7e-2 3.0e-3 3.2e-3 4.3e-3 -8.0e-16 1.5e-3 4.0e-3
```

Rise to a peak, a 30-190x collapse, a slow climb, machine zero, a slow climb -- the same shape at
the same rungs. Peak `9.08e-2 -> 9.74e-2` at `B/B_full` **0.0350 -> 0.0351**. `far` does it too,
with **both** its machine zeros moving one level (`341, 1365` becomes `1365, 5461`).

So in near-field and `far` the profile is a function of `B/B_full` and the zeros track the tree's
own depth cap. At levels 7, `B = 341` buys a complete level-4 tree and reads `7.04e-2`, the peak of
the row. **`deep interior` is the exception** -- its structure stays pinned to absolute budget
(`85, 341, 1365` are local minima at both depths) -- so this is an observation in two regions of
three, not a law.

**The question the run was launched to answer inverts.** I asked whether near-field's teeth keep
shrinking past level 5. They were never teeth: **the curve does not shrink, it rescales.**

## B, ANSWERED ON THE FIXED KERNEL: THE ROOM AND THE CAPTURE ARE IN DIFFERENT PLACES

Where the room is, at levels 7 (`headroom/err`, peak per region):

```
  near-field    peaks 9.74e-2 at B = 767  (B/B_full 0.035), then <= 4.3e-3 at every deeper rung
  far           peaks 2.17e-2 at B = 85,                    then <= 6.4e-3 from B = 767
  deep interior RISES monotonically to 1.00e-1 at B = 12287 -- the largest value in the table
```

Where the capture is (`captured`, near-field, against the headroom it is normalised by):

```
  B =              11      21      23      47      85      95     191     341     383     767    6143   12287
  headroom/err  4.8e-3  1.4e-2  1.4e-2  2.6e-2  3.5e-2  3.8e-2  5.6e-2  7.0e-2  7.7e-2  9.7e-2  1.5e-3  4.0e-3
  first_div     1.0000  1.0000  1.0000  0.7548  0.4637  0.4313  0.1831 -0.2651 -0.2858 -0.4828  0.6075  0.3774
  perim_within  1.0000  0.8411  0.8411  0.4825  0.4177  0.3888  0.1778 -0.3554 -0.3345 -0.5406  0.6075  0.3774
  frac_hot_btw  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 -279.37 -5.3601
```

**The wins sit where the headroom is small and vanish exactly as it grows.** `captured` falls
monotonically from 1.0000 at `B = 11` to negative at `B = 341-767`, which is precisely the peak.
The two rungs where a criterion is positive deep (`6143`, `12287`) carry headroom of `1.5e-3` and
`4.0e-3` -- the harness's own *"a MAGNIFIED view of a small absolute gap"*.

**`deep interior`: nothing beats breadth-first at any budget.** 0 of 13 for **every** criterion;
the least-bad is `first_div` at `-0.0060`. That is the region with the most headroom, and none of
it is taken. Raw errors, so this is not a magnification artefact -- at `B = 341`, uniform `0.17280`
against the best criterion's `0.32436`, **1.88x worse**.

## THE STANDING BEST CRITERION DOES NOT SURVIVE, AND IT IS BREADTH-FIRST IN DISGUISE

`frac_hot_between/median` -- on record as **the best criterion measured on this project** -- beats
uniform at **0 budgets in all three regions**. In `far` it reads `captured` exactly `0.0000` at
**every** budget; in near-field at ten consecutive rungs (`11` through `767`) before going
negative. It is reproducing the uniform leaf set bit for bit, which is the mechanism already on
record -- most quads tie, the tie-break is lexicographic on `(level, ix, iy)`, level first, and
that *is* breadth-first. Its 31-65 distinct values are the arithmetic ceiling of a fraction over
`N^2 = 64`, so ties are the normal case rather than the exception.

**Item E survives and is now joint-best.** `perim_within/median` -- `lay_w_perimeter`, the live
lead, never previously in the Rank list -- reaches `captured` 1.0000 and ties `first_div/median`
for the best row in near-field (positive at 9 of 15 budgets). The two converge to the **identical**
leaf set deep: errors `0.10652 / 0.09974 / 0.06037` at `B = 5461, 6143, 12287` on both.

## `greedy == dp_optimal` EXACTLY IN `deep interior`, AT BOTH DEPTHS

`greedy_lookahead_1` matches the exact optimum to five decimals at **all nineteen budgets** at
levels 7, and at all sixteen at levels 6. It is a field property, not a bug -- `deep interior` has
no level barrier in `err_sum`, so gains are independent and immediately available, which is the
condition under which greedy is optimal. The same table refutes any reading of it as a bound:
near-field `0.12865` against dp `0.06022` at `B = 12287`, **worse than uniform's `0.06046`**, and
`far` `0.26399` against `0.10837`. *Greedy is a strong reference, never a ceiling* holds, and the
one region where it is a ceiling says so by measurement rather than by name.

**And greedy's near-optimality in near-field does not survive the kernel fix**: the record quotes a
gap of `0.0004` there. It is `0.068` now, and negative `captured` from `B = 767` on.

## The confound, stated

Levels 6 and 7 score against **different references** -- `res = (1 << levels) * n`, so 512^2 and
1024^2. Raw `error` is therefore not comparable across the two runs; `headroom/err` is, because
both runs bottom out at exactly one sample per reference pixel and are self-similar by
construction. Every cross-depth claim above is made on `headroom/err` or on `captured`, never on
raw error.

## What this does not answer

Whether any criterion can be built that captures headroom where the headroom actually is. At
levels 7 every one tested captures most where there is least, and `deep interior` -- the region
with the most room -- is taken by none of them at any budget.

**WITHDRAWN: "that run is blocked by memory".** I computed `quads x 64 x 8 x 33 x 3 f64` =
**8.9 GB** for the boundary shapes, called it the cause of two kills, and wrote it into a commit
message as measured fact. It is arithmetic on a **false premise**: `boundary_shapes` lives on the
*march output* (`AzOut` / `HgOut` / `LhOut`), which is local to `evaluate` and dropped when it
returns. `PixelOut` never holds it -- the three temporal accumulators are extracted to scalars
inside `evaluate` and the series is discarded. Measured on the relaunched `levels = 7` run:
**RSS 0.01 GB**.

The likelier cause was in plain sight and I read past it: **four background tasks died
simultaneously**, two of them monitors that allocate nothing. That is a session-level cleanup, not
an out-of-memory kill. A shared failure time across unrelated processes is the signature, and no
`log show` entry supported the OOM reading either -- I checked, got nothing, and believed the
arithmetic anyway.

*The measurement was correct and the column chosen could not see it* has a sibling here: **the
measurement was never taken.** A number computed from a structure you have not read is not a
measurement, however carefully it is computed, and putting it in a commit message makes it
citable. `keep_boundary_shapes` stays **on** for `levels = 7`; nothing is traded away and the
table keeps `running_max`, `first_div` and `term_grad`.

**The finding this cost is still real and stands on its own**: the ladder never sampled a complete
level, and the comb is what was hiding behind that. It was found while looking for a cheaper route
around an obstacle that did not exist -- worth keeping, and worth saying plainly that the reason it
was looked for was wrong.

And `deep interior`'s `.qcache` was, briefly, the only stale one in `results/criterion/` after an
interrupted run regenerated its two siblings. The whole run was relaunched rather than the gap
filled, so all three caches come from one run of one binary. **An interrupted regeneration
produces the mixed-version corpus live**, which is the failure the corpus write-up describes,
arriving by a route it did not anticipate.
