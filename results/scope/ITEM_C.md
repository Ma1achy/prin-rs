# Item C — blocked, not pooled. The prediction splits: confirmed in form, refuted in strength.

`signal_audit 7 8 1e-4 13 <root> all`, levels 7, 21845 quads per region, DP label ladder
`B = [169, 681, 2729, 10921]`, 54 signals. Build times 5351 s / 2821 s / 4198 s. Full output in
`signal_audit_lv7.txt`; the per-quad signal matrices (21845 x 54) are in `audit/*.tsv.gz`.

The prediction was committed before the run in `PREDICTION_C.md`. Both halves are graded below.

## CONFIRMED: pooled correlation is largely depth, and blocked it does not merely fall — it INVERTS

`level` is the control and it settles the mechanism by itself: **blocked it is `NaN` at every
budget**, because level is constant within a level by construction, while pooled it scores
**-0.611, -0.977, -0.706, -0.690**. A signal with literally no within-level information scores
0.98 pooled. Every signal tracking cell width inherits exactly that, and `ensemble_spread` carries
a scale term by construction.

near-field, pooled against blocked at the same four budgets:

```
      signal                pooled                          blocked
  spread_median      0.486  0.701  0.601  0.572    -0.452   0.152  -0.260  -0.050
  spread_p90         0.557  0.727  0.607  0.579    -0.453   0.152  -0.260  -0.050
  between_shape      0.488  0.691  0.592  0.565    -0.038   0.152  -0.244  -0.016
  divergence_trend   0.563  0.740  0.607  0.571    -0.197  -0.116  -0.282  -0.083
  grad_rms_between   0.473  0.685  0.573  0.548    -0.047   0.152  -0.252  -0.016
  running_max_div    0.545  0.652  0.593  0.597     0.477  -0.066  -0.193   0.012
  cell_width         0.611  0.977  0.706  0.690       NaN     NaN     NaN     NaN
  level             -0.611 -0.977 -0.706 -0.690       NaN     NaN     NaN     NaN
```

**The sign flips.** Blocked, the spread signals are mostly *anti*-correlated with the DP label:
within a level, the optimum tends to split the quads with the **lower** spread. That is the
opposite of what the criterion does, and it reproduces the standing `spread_median` sign flip at
roughly twice the magnitude on the fixed kernel. `far` does the same —
`divergence_trend` pooled `-0.692, -0.744, -0.656, -0.632` against blocked
`-0.392, -0.148, -0.203, -0.043`.

## REFUTED: the signals are NOT depth in disguise, and the standing AUC number does not reproduce

The record says *"held-out AUC 0.88 -> 0.37 when `level` and `cell_width` are removed, so a fit
carrying them is a depth model wearing 55 names."* Measured on the fixed kernel, a ridge logistic
on the 54 standardised signals against the DP label at `B = 2729`:

```
                     leave-one-region-out AUC        budget hold-out AUC (fit 169 -> score 10921)
  held out        with depth   without   drop        with depth   without   drop
  near-field         0.9212     0.8735  -0.048          0.8811     0.8040  -0.077
  far                0.9204     0.6891  -0.231          0.7148     0.6935  -0.021
  deep interior      0.7678     0.6673  -0.101          0.8704     0.6436  -0.227
```

**0.88 -> 0.80, not 0.88 -> 0.37.** Removing every depth feature costs 0.02–0.23 AUC and leaves
0.64–0.80 held out. The signal set carries real information about the optimum's decisions that is
not depth. My prediction's second clause — *"the pooled-vs-blocked gap accounts for most of each
signal's apparent skill"* — is wrong.

## AND A MULTI-SIGNAL FIT BEATS BREADTH-FIRST WHERE NO SINGLE CRITERION DID

`error(B)` of the held-out fit against `uniform` and against the best single criterion:

```
  near-field   B =        23        95       383      1535      6143
    fit                0.30966   0.30565   0.28833   0.20435   0.11125
    uniform            0.30966   0.30796   0.28917   0.20434   0.09984
    frac_hot_between   0.30966   0.30796   0.28917   0.20479   0.14135
  far
    fit                0.60292   0.58973   0.54297   0.36602   0.17905
    uniform            0.60292   0.60128   0.54295   0.36598   0.17916
```

The fit is **below uniform at `B = 95` and `B = 383` in near-field** and at `B = 95` in `far`,
where item B found no single criterion ever beats it. Modest and real. It loses deep
(`0.11125` against `0.09984` at 6143), so this is the same *early-adequate, late-worse* shape
item B measured, not a reversal of it.

**`deep interior` is unbeaten by the fit too** — `0.34459 / 0.34276 / 0.33598` against uniform's
`0.29409 / 0.22886 / 0.17061`. Item B found 0 of 13 budgets for every single criterion; a 54-signal
fit does not rescue it either.

## THE MECHANISM FOR WHY SINGLE CRITERIA FAIL: 96% COLLINEARITY

Multiple `R^2` of each signal on the other 53, pooled at `B = 2729`: **median 0.9584, and 12 of 54
above 0.99**. `lay_b_n_hot`, `frac_above_tau_between`, `lay_w_n_hot` and `frac_above_tau_within`
all read **1.0000**. Most of these are functions of the same footprint spreads.

That is the honest reading of item B's null: it is not that the information is absent, it is that
54 names describe about the same thing, so picking a better *one* of them buys almost nothing —
while a fit that combines them does clear uniform at moderate budgets. **The harness's own header
says to read this before reading a null as "the criterion cannot be improved", and it was right
to.**

## The one thing the prediction flagged as unaccounted, resolved

`PREDICTION_C.md` recorded that `perim_within` and `first_div` beat uniform at shallow budgets,
which a pure-depth signal set should not have produced, and that *"I do not have a reading that
accounts for both"*. The audit supplies it: `lay_w_perimeter` is **`NaN` on 96.7%** of near-field
quads and on the 3.3% it scores it has a **positive blocked** correlation (`+0.269`, `+0.327`)
against a *negative* pooled one (`-0.012`, `-0.716`, `-0.297`, `-0.627`). It is the standing
`term_grad` shape — *a high `nan%` is a property to read, not a defect to hide* — and it is the
only family in the table whose blocked sign is opposite to its pooled sign in that direction.

## Constant signals, read before any curve

In near-field, `escape_fraction`, `n_nonfinite`, `spread_event_median`, `ftle_nan_frac`,
`layrel_w_n_hot` and `layrel_b_n_hot` all read **`distinct = 1`, `modal 100%`**. They have no
ordering at all there, and their flat `error(B)` rows are the tie-break's scan order, not a
verdict on the signal. `between_matched` and `within_pooled` are `nan% = 100`.

## Carried forward, confirmed on the fixed kernel

**The DP's labels are still monotone in budget**: `1 -> 0` flips are **exactly zero** at every rung
in near-field (48, 217, 979 flips, all `0 -> 1`), so the optimum's split sets remain strictly
nested and a single fixed priority order is not structurally excluded.

## Caveat carried from the harness, not smoothed over

The FTLE march is the **unregularised fixed-step leapfrog**, so near a close approach it is not
trustworthy — an `ftle_*` result in `deep interior` carries that and a `spread_*` result does not.
And three heterogeneous regions are not a validation set: the folds are never averaged, because
smooth / localised / everywhere are three different questions.
