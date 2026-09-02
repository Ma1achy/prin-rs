# `integrator_gallery_1024/` — the high-resolution subset, and it is STALE

123 panels, **1024²**, seven cases, from `examples/integrator_gallery.rs`. Moved here from
`results/integrator_gallery/` when the complete 32-case set was regenerated at 256²; the
canonical path now holds that set, which is what every gallery claim in `CLAUDE.md` is measured
on. Its table is `results/output/integrator_gallery_1024.txt`.

**Neither set supersedes the other.** Per-trajectory statistics — drift, `err>10`, force
evaluations — are resolution-independent and this project says so explicitly, so the 256² table
is the right place to read a number. Chord statistics are *not*, and **no chord ratio may be
quoted from a coarse grid**. `NOTES.md` reads its collision fractions (1048576 on `far`, 1033184
on `deep_interior`, 850590 on `preset_shape`) from the 1024² table and the 256² one cannot
substitute.

## Stale, and by how much

Every panel here predates three changes that move numbers:

1. the no-discard fix to `energy_drift_max` / `gamma_max` (`pixel.rs:959`);
2. the secant landing on AZ, Heggie and `reference/tb_az.py`;
3. `land_iterate` routed through `EnsembleCfg` into all three drivers.

Measured on the 256² regeneration, the effect is small in the median and not nil: `drift p50`
ratios run 0.987–1.004 over 64 rows, `err>10` moves 3916 → 3915 (AZ) and 73 → 74 (Heggie), and
eight rows gain `nonfin` counts they should always have had. **The headline holds — Heggie still
wins 31 of 32.**

## Incomplete, and it says which arm is missing

The run was stopped to free the machine and never resumed. `mid-field` has an **`az` arm and no
`heggie` arm**; the other six cases have both. `far_heggie` has a `gain` panel and `far_az` does
not, which is by design — the gain map is Heggie-against-AZ and belongs to the second arm.

## Six of these panels are pre-fix renders committed in the fix's own commit

`far_az_norepair_{dmin,drift,outcome,spread,tend,uniform}.png` carry mtime **31 Aug 09:35:30**
and were committed by **`a526360` (1 Sep 12:51)** — the commit that landed the no-discard fix and
the secant landing. **They predate the code in their own commit by 27 hours.** A reader taking
the commit as provenance would read them as post-fix; they are not.

This is the failure already on record at `220d928`, where a pre-`dtau`-fix render sat in a
post-fix tree and anything using it as a *before* was measuring the previous investigation over
again. *Commit renders in the same commit as the code that made them, or name the commit in the
filename.* The rule was written down and hit anyway, one commit later, which is why it is
recorded here rather than quietly fixed by a re-render.

The other 117 panels come from `3c70aef` (08-30) and earlier and are honestly pre-fix by their
own commits.

## Regenerating

~16 h at the current cost — the 256² run of 32 cases took 10649 s and this is 16x the samples
over 7 cases. Not done, and it is a decision rather than an oversight: these panels back no
number in `CLAUDE.md`, and the argument for 1024² is the eye, which is how the wedge result was
actually settled.

```sh
cargo run --release --example integrator_gallery -- 1024 results all 400000 0
```

**Do not point that at `results/`** without first moving this directory aside — argument two is
the output root and it will write into `results/integrator_gallery/`, which now holds the 256²
set.
