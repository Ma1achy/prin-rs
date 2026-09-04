# `results/session` — the frame loop marched, with the playhead and the camera both moving

§9's acceptance test. `session_frames` runs two arms and writes a `PRNF` record per frame beside a
provenance sidecar; `output/session_frames.txt` is the committed stdout. Reproduce:
`cargo run --release --example session_frames -- results near-field 8 96`.

## The test discriminates, and the control is what says so

| arm | leaves | final depth variance | frac over the 41.7 ms floor, during motion |
|---|---|---|---|
| `balanced` | 34 | **0.8789** | 0.2857 |
| `uniform` | 256 | **0.0000** | 1.0000 |

`Mode::Uniform` degenerates by construction — 256 leaves at one depth, variance exactly zero — and
that is the arm that proves the depth-variance test can tell the two apart. Without it a balanced
arm reading 0.88 would be an unfalsifiable number.

## Depth variance alone cannot tell BALANCED from FROZEN, and this run needed both extra columns

§3.2 is explicit that a tree stable because nothing moves and a tree stable because it is balanced
are **identical in a variance plot**. Two columns separate them here, and the first run needed both.

**`churn`, over shared quads with the count printed.** A quad present at one frame and not the other
has not changed its decision; counting it folds the tree's growth into a statistic about its
stability.

**`reproj` — the liveness arm for the playhead.** The first cut of this harness logged
`playhead_dt` while the sampler ran to `t_max` every frame, so the field never changed, the tree
converged once and sat: **churn 0.0000 for nine consecutive frames**. That is frozen, not balanced,
and the churn column caught it. `Session::set_playhead` now re-reduces every quad from the boundary
it has reached, and `reproj` reports how many — because a churn of zero means a *steady state* only
if the quads were actually re-read, and zero re-reads would be a playhead that is not wired. The
two are indistinguishable in the churn column alone.

Measured: `reproj` reads **45 on every frame from 3 onward** with churn 0.0000. All 45 quads are
genuinely re-read at each new boundary and every one decides the same thing. **A steady state, and
the arm that says so.**

## What the frame budget costs here

`p50` on the balanced arm is **0.00 ms** — most frames drain having computed nothing — against a
`max` of 2004 ms on the frames that do work. That spread is the honest shape of a quota measured on
a CPU integrator at 96²: a frame that computes 24 quads is ~2 s of trajectories, which is 48× the
41.7 ms floor. **The quota is a count and the milliseconds are a report**, not a target the
scheduler steers by — `tests/session.rs::the_quota_is_not_a_deadline` pins that a sleeping sampler
does not move a single decision.

So `frac_over_floor_moving` at 0.2857 and 1.0000 is a statement about *this integrator at this
resolution*, not about the scheduler. What it demonstrates is that the statistic exists, is measured
during motion (`camera_delta > 0` on every frame), and separates the arms.

## The frontier, now actually in the frame loop — `output/frame_rank.txt`

The frame quota truncates `pending`, the quads about to be *computed*, and before this it truncated
by **position** — which is split order, i.e. raster order within a parent. `Session::rank_pending`
puts the persistent frontier there. `rank_frame = false` is the named control and is the pre-wiring
behaviour.

**The physics ordering alone moves NOTHING, and the arithmetic says why.** A round whose `pending`
comes from one split is a **total tie**: all four children inherit their parent's stored term (they
have no reduction of their own — ranking on that would rank every child of every parent at exactly
zero, which is *no* ordering rather than a weak one), ties break by id ascending, and id order **is**
split order. Measured on two analytic fields at a binding quota on all 24 frames:

| case | arm | leaves | quads | dvar | scan/len | only-off | only-arm | moved |
|---|---|---|---|---|---|---|---|---|
| `step` | off | 145 | 189 | 1.6323 | — | – | – | – |
| | rank | 145 | 189 | 1.6323 | 0.6604 | 0 | 0 | 0 |
| | rank+cam | 145 | 189 | 1.6323 | 0.6604 | 0 | 0 | 0 |
| `filament_through_sea` | off | 145 | 189 | 1.6323 | — | – | – | – |
| | rank | 145 | 189 | 1.6323 | 0.6604 | 0 | 0 | 0 |
| | rank+cam | 145 | 189 | 1.6323 | 0.7385 | **47** | **47** | 0 |

**Only the camera term moves a tree, and only where the quota has a choice to make.** `step` reads
`keep:96` — most quads are resolved, few want to split, `pending` stays at 8.8 and the truncation
rarely bites. `filament_through_sea` reads `floor:96` — everything is unresolved, everything wants
to split, `pending` reaches 10.8, and 47 boxes enter and 47 leave. **Read the stop column**: the two
fields give the *same geometry* by opposite mechanisms, and without it the identical 145/189/1.6323
rows read as one field measured twice.

`moved` is 0 everywhere and that is expected rather than a null: a decision is a function of a
reduction, and reordering the schedule does not change any reduction. It can only move where a
differently-shaped tree gives a parent different children.

**And the drained arm is the correctness property, asserted rather than printed.** Run on to drain,
the ranking has only reordered work that all happened, so all three arms must reach the identical
tree — measured, the same 196 leaves in the same 34 frames on both fields. A ranking that changed
the *destination* would be changing the criterion rather than the schedule.

### The regime is set by the arithmetic, and the 83–94% headline does not transfer

`pending` is about `4s` times the previous round's quota, where `s` is the fraction that split, so
`k/n ~ 1/(4s)` sits in `[0.25, 1]` **by construction**. `results/frontier/`'s 83–94% was measured at
`k/n = 0.01` on the **visible frontier** — a different population — and quoting it here would be
quoting a number about another measurement. What this site actually reads is **0.66–0.74**, i.e.
26–34% saved, squarely inside that measurement's own 0–67% band for large `k`.

**Except where the criterion wants a lot.** In the `uniform` arm of `session_frames` the frontier
reaches 224–232 entries against a quota of 24 and `scan/len` falls to **0.165–0.211**, 79–83% saved
— then returns to 0.96 two frames later as the entries pile back into a few bands. So the regime is
a property of *how much the criterion wants*, not of the frame loop, and both ends are in the same
committed run.

## Re-rooting, driven by the camera

`Regrow::Upto(k)` grows the root when the camera's box leaves it, toward the camera, at most `k`
levels. `regrown` is a record column. The asymmetry is the test: a zoom-**in** must never trigger
it, or the containment test is reading something other than containment and the count would still
look plausible. Measured in `tests/session.rs`: zoom-in 0, zoom-out to 2.5× the root **2**, further
zoom-out clamped to the remaining 1, then 0.

**A regrow disturbs nothing.** `grow_root` pushes the new root at the end, so no index moves and
neither the store nor the frontier is touched — asserted through the session over every existing
box, decision and the resident payload count.

**And a grown root leaves the absolute UV lattice in three directions of four**, which is the seam
between this rooted tree and the caching contract's flat store. See `tests/uv.rs`.

## Notes carried in every record

`upload_ms` and `present_ms` are **`NaN`**: there is no GPU and no window in this build, and `0.0`
would read as *instant* where the truth is *absent*. `frontier_agrees` is the audit's real answer
now that the frontier is wired, and `NaN` on frames where it did not run — a check that did not run
must not report a pass. `frontier_scan`/`frontier_len` are **0 on a frame whose quota did not
bind**, where `scan/len` is undefined rather than perfect.

The record is **PRNF v2**: `frontier_scan`, `frontier_len` and `regrown` joined it. Version bumped
rather than appended silently — a reader that trusts a field *count* and not a version is the
mixed-version corpus defect one level down.
