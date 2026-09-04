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

## Notes carried in every record

`upload_ms` and `present_ms` are **`NaN`**: there is no GPU and no window in this build, and `0.0`
would read as *instant* where the truth is *absent*. `frontier_agrees` is `NaN` because the audit
did not run — a check that did not run must not report a pass.
