# `results/live` — the tree growing while the simulation runs

Phase 4 of the refinement rebuild: the **time-axis** animation, one `descend_live` per chart on the
committed mechanism (`Policy::Tolerance`, `alpha_lo = 0.005`, `dim_floor` on, the agreement arm
on, merging on, stationarity off), rendered 2026-09-04 at the commit after `62df9b2`. Each frame is
the leaf set at one sync boundary, drawn on the footprints **as they stood at that boundary**
(`scheduler::project_at`: the shape and spread the copies had reached by then, a copy's terminal
class only once it had terminated), coarsest ancestors filled underneath; then one frame per
post-horizon round at `t = t_max`, and a six-frame hold on the finished tree. `_live.png` is the
colour field (hue from the shape sphere, lightness from `spread_shape` on one window per chart,
taken from the finished tree so the ramp does not move between frames); `_live_wire.png` is the
same with the leaf boundaries over it. Every decision in every frame reads the trajectories only
up to that boundary; nothing after it.

**Diagnostic viewport.** These are 256², so a chart is a minute or two; the screen floor binds at
level 5 here where the 512² tables in `results/payload/README.md` reach level 6, and the leaf
counts below are facts about this viewport. Do not read a leaf count off a frame. Every panel's
`.cfg.txt` carries the full configuration, the growth curve (`boundary:t:quads_computed`), the
stop-reason breakdown and the catch-up share.

| chart | quads | leaves | depth | frames | catch-up | stop |
|---|---|---|---|---|---|---|
| near-field | 85 | 55 | 5 | 25 | 65.4% | floor:1 keep:42 screen_floor:12 |
| deep interior | 145 | 70 | 5 | 23 | 41.8% | floor:3 keep:32 screen_floor:35 |
| preset_shape_h1 | 833 | 439 | 5 | 27 | 51.9% | floor:59 keep:41 screen_floor:339 |
| config_stability | 1373 | 955 | 5 | 25 | 60.1% | floor:24 keep:85 screen_floor:846 |

What to look for. On `preset_shape_h1` the tree grows from 33 quads at the first boundary to 721
at the horizon and 833 after the post-horizon rounds, and the growth is where the mixing region
develops: the regular island stays coarse. On `config_stability` the ribbons sharpen boundary by
boundary and the tree follows them. Leaves that merge back are drawn as their parent again in the
next frame. `adjacent_duplicates=7` in each sidecar is the six-frame hold plus one post-horizon
round that changed nothing.

Reproduce (writes here; `refine_flagged` is off because the repair pass has no live analogue):

```
cargo run --release --example live_animation -- results "near-field,deep interior,preset_shape_h1,config_stability" 256 4000
```
