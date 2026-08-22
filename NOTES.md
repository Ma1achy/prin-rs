# Notes from reading the brief and the reference

Working notes, not spec. `BRIEF.md` is authoritative; where this file disagrees with it, that
disagreement is the point and wants settling.

---

## 1. The shared-reference question — a recommendation

**BRIEF.md §3.** Should all `E+1` copies of a pixel share the nominal copy's AZ reference body?

**Recommendation: default UNSHARED, implement SHARED as the flag, and add a per-pixel
`ref_disagree` field so the question can be answered from one run rather than two.**

The reasoning, in three parts.

### The case for sharing is real

`error_ratio` and `spread_shape` are *cross-copy* statistics. If copies are integrated in different
regularisation charts, their integration errors are drawn from different distributions with
different structure. The energy spread then inflates for a reason that is not physics, and
`error_ratio` — whose whole value is that it has no threshold and no tuned constant — reads damage
that isn't there. Each trajectory's own `energy_drift` stays healthy throughout, so nothing else
flags it. That is exactly the failure mode the brief describes, and it is plausible.

### The case against sharing is also real, and is about correctness rather than diagnostics

AZ's guarantee rests on choosing the reference body **not in the longest side**, which makes the
unregularised separation `R3` the longest, so `R3 → 0` only in a genuine triple collision
(`reference/tb_az.py:145`). Force a copy onto another copy's reference and that copy may be
integrated with its reference body *inside* its longest side. Then `R3` is no longer the longest,
and an unregularised pair can close. The AZ guarantee is void for that copy — not degraded,
**void**.

So sharing trades a diagnostic artifact for a genuine integration failure. In general that is the
wrong direction.

### Why it is worse than a simple trade

The jitter is small (`0.5 ×` cell width), so copies usually agree on the reference anyway.
Disagreement concentrates where the triangle has two near-equal longest sides — where the
`argmax` in `choose_reference` is near-degenerate. Those are precisely the pixels where sharing is
most dangerous *and* where its diagnostic benefit would be largest. The two effects are not
independent; they are the same pixels. Neither default is safe there, which is why I would rather
**measure** than choose.

### The concrete suggestion

Record per pixel, per sync point, how many copies would have chosen a reference different from the
nominal copy's — accumulate to a `ref_disagree` count, alongside the existing `switches`. Then:

- Run unshared. Condition `error_ratio` on `ref_disagree`.
- If inflated `error_ratio` is uncorrelated with `ref_disagree`, the hypothesis is dead and the
  flag can default unshared permanently.
- If it correlates, the effect is real and quantified, and the shared run is worth its cost.

This answers BRIEF.md §8 experiment 2 with one run plus a confirmation, rather than two runs and a
difference of aggregates that carries no attribution. It costs one integer per pixel.

**One clarification for the flag's scope.** The brief says the reference "can change mid-run" and
raises sharing in the same breath. These are separate knobs: sharing is *across copies*, switching
is *across time*. The flag should govern sharing only. Freezing the reference across time would
break AZ outright — the whole point of re-choosing is to keep the reference out of the longest
side as the triangle deforms.

---

## 2. Three things in the brief I think are wrong or incomplete

Raised because the brief asks for it.

### 2.1 `d_min` as specified is not what the reference measures

BRIEF.md §4 defines `d_min` as "closest approach over the whole trajectory". The reference tracks
the minimum of `|R1|` and `|R2|` only (`reference/tb_az.py:201-203`) — the two *regularised* pairs.
`|R3|`, the unregularised side, never enters `d_min`.

At the moment the reference body is chosen, `|R3| ≥ max(|R1|,|R2|)`, so this looks safe. But the
reference is re-chosen only at sync points, and between them the triangle deforms freely. A close
approach on the `R3` pair inside a sync interval is invisible to `d_min`.

Whether this matters depends on `n_sync` and how fast the geometry turns over. It is at minimum a
discrepancy between the spec sentence and the validated code, and a Rust port matching the
reference to `1e-10` will inherit it. **Decide which behaviour is wanted before porting**, because
"closest approach over the whole trajectory" and "what `tb_az.py` returns" are not the same
quantity and the acceptance test cannot distinguish them.

### 2.2 The `error_ratio` acceptance test is in tension with its own specification

BRIEF.md §4 says `error_ratio` should be aggregated by `max`, treated as a **boolean flag**, and
that **its magnitude is unstable**. BRIEF.md §5 then requires it to equal `1.0000` — four decimal
places — as an acceptance test.

A max over ~10^6 pixels of a statistic whose magnitude is declared unstable is the least likely
quantity in the system to reproduce to four decimals. Either the test means the *median* or some
other robust aggregate (in which case §5 should say so), or it means max and the tolerance needs
stating as a bound (`max error_ratio < 1.01`, say) rather than an equality. As written a port can
fail this test while being correct, or pass it by aggregating the way that happens to be quiet.

### 2.3 Non-finite copies poison `d_min` silently

`reference/tb_az.py` correctly treats non-finite state as `done` so the step budget is not burned.
But `d_min` is updated with `np.minimum` *before* that check bites, and `np.minimum(x, nan)` is
`nan`. A copy that diverges leaves `d_min = NaN` for its slot.

This is consistent with "never discard a copy" — NaN *is* the measurement outcome — but only if
downstream reduction treats it as such. Any `min`/`max` over copies will propagate the NaN to the
pixel. The Rust port should decide this deliberately rather than inherit it: I would carry a
per-copy `finite` flag and let the pixel record "one copy undetermined" explicitly, which is the
same information without a value that silently contaminates every aggregate it touches.

---

## 3. Smaller observations

- `reference/tb_az.py` uses `eta=0.02` as its default; `tb_all_az` and the smoke test pass
  `eta=0.01`, which is what BRIEF.md §2.3 specifies. The bare default is the odd one out.
- The reference's `dtau` sizing (`eta*dt_left/(A0*B0)`) fixes the *first* physical step as a
  fraction of the remaining interval, so effective resolution depends on `n_sync` through
  `dt_left`. Worth pinning `n_sync` explicitly in the cross-check test rather than leaving it at a
  default, or the `1e-10` agreement target is comparing two different discretisations.
