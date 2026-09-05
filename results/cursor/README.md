# `results/cursor` — §18 foveation, measured against its own off-state

`Camera::foveation` modulates camera relevance. It is **not a third factor**: it lives in priority,
never in the veto, and no cursor field goes on a `Quad`. With no cursor, or a cursor still in motion
(`dwell = 0`), it returns exactly `1.0` and `priority` is §4.3's product unchanged — **the fallback
is the default path**, which every keyboard, touch, unfocused-window and headless route takes.

The plan set the rule before the numbers: *if it does not beat uniform on time-to-resolve-at-cursor
by a clear margin it is dropped.* Measured 2026-09-05. Reproduce:
`cargo run --release --example cursor_bias -- results 40 12`. Committed stdout in `output/`.

## The verdict: it stays OFF

It does **nothing at all** on a field with localised structure, and buys real frames only on a field
that is unresolved everywhere — and only at a cap that demotes the periphery sixteenfold. That is
not a clear margin; it is a narrow one with a named condition, and `SchedCfg::cursor` stays `None`.

```
== step                                         == filament_through_sea
       arm  t_cursor    t_edge                         arm  t_cursor    t_edge
       off        12        22                         off        16        14
    centre        12        22                      centre        16        14
    moving        12        22                      moving        16        14
  fovea x4        12        22                    fovea x4        14        16
 fovea x16        12        22                   fovea x16         8        16
 dwell 0.5        12        22                   dwell 0.5        14        16
 x16 @edge        12        22                   x16 @edge        16         8
```

**`step` is inert in every arm.** Only the quads straddling the discontinuity want to split, and
there are few enough that the quota never has to choose between distant ones — so there is no
ranking for the fovea to change. A foveation measured only on the sea chart would have looked like
a general result.

**On the sea, cap 4 is a wash and cap 16 is a gain.** `fovea x4` saves 2 frames at the cursor and
loses 2 at the edge: **the same budget moved**, which is what the edge column exists to catch.
`fovea x16` saves **8** and loses **2** — a net gain, on one field of two.

## Four controls, and one of them changed the reading

**The advantage follows the CURSOR, not the place.** The two probes are geometrically symmetric and
their baselines are not: the tie-break is lexicographic on `(level, ix, iy)`, so the lower-`y` probe
is reached first and the `off` arm resolves it two frames sooner. A gain measured at one probe could
therefore be the scan order. Moving the fovea onto the other probe takes **edge 14 → 8** while the
cursor probe stays at 16 — the effect is the fovea. Without this arm the 8 frames were not
attributable.

**`dwell = 0` reproduces `off` bitwise**, asserted rather than printed. That is the
fallback-is-the-default-path claim as a test.

**The `centre` arm is the named vacuous cell and behaves as named.** A cursor at the frame centre is
concentric with the viewport, so the fovea is nearly `relevance` itself; it moves nothing. A
measurement run only at the centre would have reported "no effect" and been reporting an identity.

**And the final tree is bitwise identical in every arm** — 196 leaves, 261 quads, 0 decisions moved.
Foveation changes *when* a region resolves and never *what* the tree becomes, which is the same
correctness property the frame frontier's drained arm carries.

## Two things the sweep says about the knobs themselves

**`fovea_cap` is not a dial over its whole range.** At this separation the Gaussian weight at the
edge is ~`1e-13`, so `foveation` is essentially `1/cap`. The first cut of this harness had both
probes off the structure and read x4 and x16 as *identical*, which was written up as saturation —
wrong, and caught by moving the probes: 2 frames against 8.

**`dwell` is continuous in the priority and discrete in the outcome.** `dwell = 0.5` at cap 4 gives
a peripheral factor of 0.625 against 1.0's 0.25 — a different number — and the **same tree
evolution**, because a ranking only moves when the demotion crosses another quad's priority.

## If it is ever reconsidered

The condition under which it pays is stated: a frame where the criterion wants far more than the
quota can supply, everywhere at once. That is the `uniform`-mode regime the frontier measurement
also picks out (`scan/len` 0.165 there against 0.66–0.74 elsewhere), and it is the regime a deep
zoom into a chaotic slice produces. It is not the regime a settled view is in.
