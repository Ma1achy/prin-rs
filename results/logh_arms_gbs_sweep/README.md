# `logh_arms` GBS sweep — a SEPARATE run, deliberately not merged into `results/logh_arms/`

These panels come from `logh_arms 256 <case> 400000 **gbs**` — the four-arm set
(`heggie`, `logh_rk4`, `logh_gbs`, `logh_gbs_nolim`). The committed corpus in
`results/logh_arms/` comes from the six-arm set.

**They are not comparable by eye, and the reason is in the sidecars.** The shared shape window is
taken from the `az` arm, which the `gbs` set does not run, so the two sets normalise differently:

    committed  shape window = (6.275959216273652e-5, 4.91546835967575e-1)
    this run   shape window = (3.121206872935264e-5, 4.91601301118318e-1)

The sweep wrote into `results/logh_arms/` and overwrote 77 committed panels before this was
noticed. Those were reverted; these are the new ones, kept here rather than merged.

**Both sets predate three changes that move the field** — the no-discard fix to
`energy_drift_max`, the secant landing on AZ and Heggie, and the resulting default changes. Neither
directory is a current science image.

The sweep's actual result is the table in `results/output/logh_gbs_controlled.txt`, which carries a
full settings header per run and does not depend on a shared colour window.

*Commit renders in the same commit as the code that made them, or name the run in the path.* This
directory is the second option.
