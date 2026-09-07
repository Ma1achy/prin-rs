# Item C, predicted before the run

`level` scores `|rho| = 0.993` against the DP split label pooled, so any signal tracking cell
width inherits it. The blocked column asks what survives *within* a level.

## What item B just established, and the hypothesis it hands to C

1. Every criterion tested is beaten by breadth-first where the headroom is largest
   (`deep interior`: 0 of 13 budgets; near-field: negative across its whole `7.0e-2 - 9.7e-2`
   peak).
2. `frac_hot_between/median` reproduces the uniform leaf set **bit for bit** -- `captured`
   exactly `0.0000` at every budget in `far` -- because most quads tie and the tie-break is
   lexicographic on `(level, ix, iy)`, level first, which *is* breadth-first.
3. So the criteria that tie a lot ARE breadth-first, and the ones that order finely LOSE to it.

**Hypothesis: the signals are noisy depth proxies.** Pooled they score well because depth nearly
solves the labelling task; blocked within a level they carry little, and the noise they add on top
of depth is what makes them lose to the clean version of the same ordering.

## The prediction

**Blocked within level, every signal's correlation against the DP label collapses toward zero**,
and the surviving magnitudes are small enough that the pooled-vs-blocked gap accounts for most of
each signal's apparent skill.

- If it holds, item B's negative result has its mechanism, and **the search for a better criterion
  over this signal set is finished** -- not "no criterion found yet" but "these signals contain
  depth and little else". The next move would be a signal that is scale-free by construction.
- If some signal keeps a substantial blocked correlation, that signal is the lead, and the
  question becomes why its `error(B)` does not reflect what its correlation says.

## What would make this measurement vacuous, checked before reading it

- **A degenerate block.** A level where every quad carries the same label has no correlation;
  `NaN` is the honest value and must not fold to 0. `deep interior` at levels 6 already showed
  `rho` NaN on degenerate rows, and that behaviour is wanted.
- **Population.** `Dp::labels` maps internal -> split, leaf -> keep, **absent -> absent**. At
  `B = 2729` the tree holds 2729 nodes of 21845, so a statistic over "all quads" is mostly
  invented labels. The population must be reported with every number.
- **Budget dependence.** The DP's labels are monotone in budget (`1 -> 0` flips exactly zero at
  every rung), so the split sets are nested and a single fixed ordering is not structurally
  excluded. Blocked `rho` should therefore be read at several budgets, not one.

## The one prediction I expect to be wrong about

`perim_within/median` and `first_div/median` tie for the best `error(B)` row in near-field and
converge to an identical leaf set deep. If the collapse hypothesis is right they should *also*
collapse blocked -- but they are the two rows that actually beat uniform at shallow budgets, so a
signal set that is pure depth should not have produced them. **I do not have a reading that
accounts for both**, and that is stated here rather than resolved after the fact.
