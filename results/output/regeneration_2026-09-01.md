# Corpus regeneration, 2026-09-01

Record of the regeneration run after the no-discard fix and the secant landing. The
artefacts it produced are now in `results/integrator_gallery/` and `results/logh_arms/`;
this is the verification that put them there. The exact invocation is at the bottom.

other is not, and treating them alike would repeat a failure this project already has on record.

## What changed in the code since every committed panel was made

Three changes, all landed this session, all of which move numbers:

1. **The no-discard fix to `energy_drift_max` / `gamma_max`** (`src/ensemble/pixel.rs:959`).
   The old reduction was `if x.is_finite() { max }`, which *discarded* a diverged or truncated
   copy and left the pixel with a finite, healthy-looking maximum over the copies that survived.
   That is a §4.3 violation: **a badly-integrated trajectory is a measurement outcome, not
   missing data.** The guard is now `o.finite && !o.budget_exhausted` — tested separately,
   because the occupants do not agree on what `finite` means (`az::driver` leaves it `true` for
   a run it truncated; `logh::driver` does not).
2. **The secant landing** (`land_iterate`, default on) on AZ, Heggie **and** the NumPy
   reference. `clamp_final_step`'s first-order predictor left an `O(h^2)` residual that capped
   observable convergence order near 2; the secant iteration removes it.
3. `land_iterate` / `land_max_iters` routed through `EnsembleCfg` into all three drivers,
   replacing a hardcoded `land_iterate: false` in `LhOpts` whose justification had expired.

## Stage A — `integrator_gallery`, 32 cases, 256²

`integrator_gallery 256 results_regen all 400000 0` — 608 panels + 608 sidecars, 10649 s.

**The baseline it is checked against is not in the tree.** `results/output/integrator_gallery.txt`
is appended to by every run, and the 256² 32-case section — the artefact CLAUDE.md's strongest
Heggie claim rests on — was overwritten by the later 1024² run. It survives only at commit
`70cfbc4`, and that is where the comparison below was taken from:

```sh
git show 70cfbc4:results/output/integrator_gallery.txt
```

### Verification against `70cfbc4`

| claim | base `70cfbc4` | regen | verdict |
|---|---|---|---|
| Heggie wins on `drift p50` | 31 / 32 | **31 / 32** | holds |
| `err>10` total, AZ | 3916 | **3915** | holds |
| `err>10` total, Heggie | 73 | **74** | holds |
| the one Heggie-worse case | `burrau_nu_k` (10 v 9) | **`burrau_nu_k` (10 v 9)** | holds |
| shared chart-property rows | 4 | **4** (same 4) | holds |
| AZ failures cleared *completely* | 13 of 13 | **12 of 13** | **moved — see below** |

`drift p50` ratio, new over base, across all 64 rows: **min 0.987, median 1.000, max 1.004.**
The landing does not move the median drift. That is expected and is on record one level up —
*a diagnostic field is specific to a class of defect, and energy drift is blind to this one*:
the clamp bought 24,000x on the figure-eight while moving `near-field`'s median drift 37x the
wrong way.

Cost: **+0.7% of the force evaluations.** The wall clock read +13%, and that figure is the
machine, not the landing — see the correction under Stage B. Stage A's table has no `evals`
column, which is what made the wrong reading reachable here.

### The one claim that moved, and its cause is the landing, not the no-discard fix

CLAUDE.md: *"**Every integration failure clears completely**: 423, 438, 703, 2222, 33, and eight
smaller, all to zero."* That is now **12 of 13**. `deep_interior` reads AZ 434 / Heggie **1**,
where it read 438 / 0.

The cause is separable and was separated: `error_ratio` is computed at `pixel.rs:766` from the
energy arrays and **never consults `budget_exhausted`**, so the no-discard change cannot touch
it. The mover is the secant landing. One pixel of 65536, in the region already on record as
*not resolved at the shipping settings and never will be* — `chord max` 1.999, antipodal, at
every `eta` rung. Recorded rather than smoothed; it is not evidence of a regression and it is
not evidence of nothing.

### The no-discard fix fires exactly where the budget was exhausted, and nowhere else

Eight rows moved `nonfin`, every one upward, every one on a row already carrying
`budget_exhausted`. The remaining 56 rows are unmoved.

```
        config_stability      az   nonfin    0->1      budget    1->1
        config_stability  heggie   nonfin    0->1      budget    1->1
           deep_interior      az   nonfin    0->42     budget   41->42
           deep_interior  heggie   nonfin    0->199    budget  199->199    <- exact
            preset_shape      az   nonfin    0->4      budget    4->3
            preset_shape  heggie   nonfin    0->217    budget  219->219
         preset_shape_h1      az   nonfin    0->11     budget   17->8
         preset_shape_h1  heggie   nonfin    0->5      budget    9->9
            shape_sphere      az   nonfin    0->0      budget    1->0
```

`deep_interior heggie` is the exact match, and it is the row the code comment cites: 199 pixels
carried `budget_exhausted` while `nonfin` read 0. The others are near-matches rather than
identities, which is what the two counts mean — `nonfin` is per-pixel over copies, `budget`
counts pixels with any exhausted copy, and the shape path and the drift path have different
finiteness conditions. **A count that matched exactly on all nine rows would be the suspicious
result**, not this one.

The four rows whose `budget` count itself moved (`deep_interior az` 41→42, `preset_shape az`
4→3, `preset_shape_h1 az` 17→8, `shape_sphere az` 1→0) are the landing's extra steps pushing
pixels across the 400000 cap in both directions.

## Stage B — `logh_arms`, 6 cases x 6 arms, 256²

`logh_arms 256 results_regen all 400000 all`. Like-for-like: the committed `logh_arms/` is
**256² on both sides**, so unlike Stage A this one is a clean replacement.

### The unregularised arms are bitwise identical, and every regularised arm moved

Complete, 48 rows, eight arms. Identical means **every printed column** — `drift p50`, `p99`,
`err p99`, `err>10`, `steps`, `evals`, `nonfin`, `budget`, `over`, `escape`.

| arm | identical / total | has a time transformation |
|---|---|---|
| `az` | 0 / 6 | yes |
| `heggie` | 0 / 6 | yes |
| `logh_rk4` | 0 / 6 | yes |
| `logh_lf` | 0 / 6 | yes |
| `logh_rk4_nolim` | 0 / 6 | yes |
| `logh_lf_nolim` | 0 / 6 | yes |
| **`plain_rk4`** | **6 / 6** | **no — `dt/ds = 1`** |
| **`plain_lf`** | **6 / 6** | **no — `dt/ds = 1`** |

**12 / 12 and 0 / 36. The split falls exactly on the time transformation, with no exception on
either side.**

**That is the landing's mechanism confirmed by its own null, and the control was already in the
harness for a different reason.** The secant iteration acts on the residual left by landing a
step on a sync boundary, and that residual exists only where `dt` is a function of the state.
The `plain_*` arms run `dt/ds = 1`: the step lands on the boundary by construction, there is
nothing to iterate, and they come back bit for bit. The arms with a time transformation are
exactly the arms that moved. The `plain_*` pair was built as the stepper-only control — to
separate stepper from regularisation — and it happens to be precisely the null this change
needs.

### Correction: the cost is 0.7% of the force evaluations, not 13% of the wall clock

| arm | `evals` new/base | `secs` new/base |
|---|---|---|
| regularised (`az`, `heggie`, `logh_rk4`, `logh_lf`), total | **x1.0071** | x0.9948 |
| unregularised (`plain_rk4`, `plain_lf`), total | **x1.0000** | **x1.1350** |

**The plain arms did bitwise-identical work and timed 13.5% slower.** That is the whole of the
+13% quoted for Stage A: machine load, not the landing. The regularised arms did 0.7% *more*
work in 0.5% *less* time.

This is the project's own standing rule firing on the person quoting it — *read `steps`, not
`secs`; under load 85-100 the winning row timed faster than the baseline while doing 1.7% more
work.* Stage A's table carries no `evals` column, which is what made the wrong reading reachable
there and not here.

### Cross-harness agreement

`integrator_gallery` and `logh_arms` are separate binaries with separate render paths. On the
six cases they share, their `az` and `heggie` rows agree **to every printed digit** —
`config_stability az` 4.210e-8 / 424 / nonfin 1 / budget 1 in both; `deep_interior heggie`
9.243e-10 / 1 / 199 in both. Two independent harnesses, same numbers.

## Landing table — the two stages are NOT alike

| regen tree | committed tree | resolution | verdict |
|---|---|---|---|
| `results_regen/logh_arms/` | `results/logh_arms/` (270 panels) | **256² both sides** | like-for-like; clean replacement |
| `results_regen/integrator_gallery/` (608 panels, **32 cases, 256²**) | `results/integrator_gallery/` (123 panels, **7 cases, 1024²**) | **mismatch** | **must not overwrite** |

The committed gallery panels are **1024²**. This regeneration is **256²**. Copying one over the
other would replace committed high-resolution artefacts with a 16x smaller raster — the failure
already on record twice, where a `criterion_metric -- 3 8` validation pass overwrote committed
512² artefacts with 128x64 ones and the small raster read as a rendering fault rather than a
stale file. **The separation of generation from swap is what caught it here**; a runner that
wrote straight into `results/` would have destroyed them silently.

The two are not substitutes and neither supersedes the other:

- the **256² 32-case** set is what every gallery claim in `CLAUDE.md` is measured on, and it is
  complete;
- the **1024² 7-case** set is the high-resolution subset the wedge result was settled on by eye,
  and it was interrupted (`mid-field` has an `az` arm and no `heggie` arm).

Per-trajectory statistics — drift, `err>10`, force evaluations — are resolution-independent and
the project says so explicitly. Chord statistics are not, and **no chord ratio may be quoted
from a coarse grid**; none is quoted here.

The 1024² set is stale under the three changes above. It backs **no number in `CLAUDE.md`** —
every gallery claim there is measured at 256². It *is* read by `NOTES.md:2050`, which quotes
`coll` of 1048576/1048576 on `far`, 1033184 on `deep_interior` and 850590 on `preset_shape` to
argue that every region tested is collision-dominated. Those are 1024² counts and the 256²
table cannot substitute for them. Regenerating the 1024² set is a ~16 h job at the new cost and
is a separate decision.

## Found during verification, not fixed here: the same violation one layer up

`src/scheduler.rs:357-366` reduces `N x N` footprints to a quad with

```rust
error_ratio_max:    px.iter().map(|p| p.error_ratio)
                      .filter(|x| x.is_finite())      // <- discards the flagged pixel
                      .fold(0.0f64, f64::max),
worst_energy_drift: px.iter().map(|p| p.energy_drift_max)
                      .filter(|x| x.is_finite())      // <- discards the undetermined pixel
                      .fold(0.0f64, f64::max),
```

That is the **same §4.3 no-discard violation** the `pixel.rs` fix was for, at the aggregation
layer, on the same two quantities. `stats::max_dev` returns `T::infinity()` for a non-finite
copy — deliberately, its doc says so — and `error_ratio` is built on it, so both fields can be
`+inf` and both filters drop exactly the pixels the statistics exist to flag.

**The `pixel.rs` fix makes this site's blind spot strictly larger.** Before it, a
budget-truncated pixel contributed a finite, healthy-looking drift and at least reached the
`max`. After it, that pixel is `+inf` and the filter drops it entirely — so a quad whose pixels
are *all* undetermined folds to `0.0`, the identity, and reads as **perfectly clean**. That is
the project's own inversion signature at a new site: *a statistic can report maximum confidence
precisely when it is least informed.* Ask what it would say about a quad nothing is known
about; the answer is "zero drift".

Scope, measured rather than assumed:

- `worst_energy_drift` is **not** a shipping `Criterion` variant — `quad.rs:48` calls it
  *"trust, measured alongside and never enforced (§2.1)"*. **No tree in the corpus is wrong
  because of it.**
- It is dumped to every `.prnq` (`output/tree.rs:84`) and `.qcache` (`output/qcache.rs:77`), and
  read as a candidate signal by `examples/signal_audit.rs:234`, `oracle_audit.rs:485`,
  `criterion_metric.rs:187` and `open_items.rs`. So the damage is to **reported columns and
  signal-audit rows**, not to allocation.
- `src/render.rs:88` filters the same field and is **correct**: that one sets a colour-ramp
  window, and `colour::drift_rgb` paints the non-finite veto set magenta separately, so an
  undetermined pixel is visible rather than silently mixed into the ramp. Not every filter on
  this field is the bug — check what the reduction feeds before changing it.

**Not fixed in this pass, deliberately.** Stage B was in flight and rebuilding `target/release`
under a running three-hour job risks the run for no gain: neither `integrator_gallery` nor
`logh_arms` touches the scheduler (checked, not assumed), so the fix cannot change either
stage's output, and the scheduler corpus that *does* read these columns is not regenerated in
this pass. Fix first, regenerate `charts/`, `criterion/` and `vertical/` after — in that order,
or the corpus is mixed-version again.

`band_of`'s convention is already the answer for the propagated value: **`NaN` to the bottom,
`+inf` to the top** — undetermined and maximally-important are different, and the ranking layer
already distinguishes them.

---

## Two fixes made during verification

### 1. `scheduler.rs` — the no-discard violation at the aggregation layer

`error_ratio_max` and `worst_energy_drift` now reduce through `scheduler::max_no_discard`, which
distinguishes the three cases the filtered fold collapsed:

- **`+inf`** propagates — the quad is undetermined and says so;
- **`NaN`** is skipped — `error_ratio` is `0/0` when `sigma_E(0) == 0`, which is structurally
  undefined rather than damaged;
- **nothing determinable** returns `NaN`, **never `0.0`**.

The asymmetry is load-bearing and not tidiness: `f64::max` already ignores `NaN`, so deleting the
filter alone would have looked right while still folding an all-`NaN` quad to `0.0`.

`tests/no_discard.rs` gains two tests, and the third arm of the first is what makes it a test
rather than decoration: **the old filtered form is computed on the same data and asserted to
disagree.** It does — finite against `+inf`. Plus the usual pair: a *subject* arm (the starved
quad really does contain undetermined footprints) and a *control* arm (the healthy quad is
finite, so a reduction hard-wired to `inf` cannot pass).

No tree in the corpus changes: neither field is read by `signal()`, and `quad.rs:48` calls them
*"trust, measured alongside and never enforced"*. The damage was to `.prnq` / `.qcache` columns
and `signal_audit` rows.

### 2. `tools/contact_sheet.py` — the output root was a constant

Third site of *an output root is an argument, not a constant*, after `criterion_metric` and
`pan_sequence`. It hardcoded `results/logh_arms`, so running it against this tree would have
written the montages straight into the committed one — pairing new panels with an old sheet, or
the reverse. It now takes `--root`, and refuses a directory that does not exist rather than
creating one. The nine sheets here were generated with it, and `results/` stayed untouched,
which is the fix demonstrating itself.

## Swapping

`results_regen/swap.sh`, run from the repo root. It refuses to run unless `RUNLOG.txt` carries
`REGEN_COMPLETE`. What it does, and why the two stages are handled differently, is at the top of
the script and in the landing table above. It also rewrites the `image=` line in every sidecar,
which names the root the panel was generated under.

---

## The run, exactly as invoked

```sh
#!/bin/zsh
# Corpus regeneration. Writes to `results_regen/`, NEVER to `results/`.
#
# The swap into `results/` is deliberately NOT automatic: a run that dies halfway would otherwise
# leave a half-new, half-old directory, which is the mixed-version corpus this whole exercise
# exists to remove. Generation and swap are separate acts.
cd /Users/malachy/src/principia-rs-test
export PATH="$HOME/.cargo/bin:$PATH"
L=results_regen/RUNLOG.txt
log() { echo "[$(date -u +%H:%M:%S)] $*" >> $L; }

log "START -- regenerating under the post-session physics:"
log "  no-discard fix to energy_drift_max/gamma_max, secant landing on AZ+Heggie+reference,"
log "  land_iterate default on. Every committed panel predates these."

log "STAGE A: integrator_gallery, 32 cases, 256^2"
if ./target/release/examples/integrator_gallery 256 results_regen all 400000 0 \
     >> results_regen/integrator_gallery.txt 2>&1; then
  log "STAGE A done ($(ls results_regen/integrator_gallery 2>/dev/null | wc -l | tr -d ' ') files)"
else
  log "STAGE A FAILED (exit $?) -- stage B still runs"
fi

log "STAGE B: logh_arms, 6 cases x 6 arms, 256^2"
if ./target/release/examples/logh_arms 256 results_regen all 400000 all \
     >> results_regen/logh_arms.txt 2>&1; then
  log "STAGE B done ($(ls results_regen/logh_arms 2>/dev/null | wc -l | tr -d ' ') files)"
else
  log "STAGE B FAILED (exit $?)"
fi

log "END -- nothing has been swapped into results/. Verify counts, then swap."
echo REGEN_COMPLETE >> $L
```

## Timing

```
[13:29:32] START -- regenerating under the post-session physics:
[13:29:32]   no-discard fix to energy_drift_max/gamma_max, secant landing on AZ+Heggie+reference,
[13:29:32]   land_iterate default on. Every committed panel predates these.
[13:29:32] STAGE A: integrator_gallery, 32 cases, 256^2
[16:27:09] STAGE A done (1216 files)
[16:27:09] STAGE B: logh_arms, 6 cases x 6 arms, 256^2
[19:46:54] STAGE B done (540 files)
[19:46:54] END -- nothing has been swapped into results/. Verify counts, then swap.
REGEN_COMPLETE
```
