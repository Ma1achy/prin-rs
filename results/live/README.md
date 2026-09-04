# `results/live` — the tree growing while the simulation runs

Phase 4 of the refinement rebuild: the **time-axis** animation, one `descend_live` per slice on the
committed mechanism (`Policy::Tolerance`, `alpha_lo = 0.005`, the dimension floor and agreement arm
on, merging on, stationarity off), rendered 2026-09-04. Thirty-one slices: the eight Burrau regions,
`config_stability`, and all twenty-two gallery charts.

Each frame is the leaf set at one sync boundary, drawn on the footprints **as they stood at that
boundary** (`scheduler::project_at`: the shape and spread the copies had reached by then, a copy's
terminal class only once it had terminated, and the copies known unusable only once they had gone),
with the coarsest ancestors filled underneath. Then one frame per post-horizon round at `t = t_max`,
and a six-frame hold on the finished tree. `_live.png` is the colour field — hue from the shape
sphere, lightness from `spread_shape` on one window per slice, taken from the finished tree so the
ramp does not move between frames. `_live_wire.png` is the same with the leaf boundaries over it.
Every decision in every frame reads the trajectories only up to that boundary and nothing after it.

## The debug flag is off

`colour::DEBUG_NAN` — magenta — marks a footprint with no value so it cannot be read as a dark one.
That is a diagnostic device, and in an animation it is a speckle the eye reads before the field.
These panels use `colour::Veto::Quiet` instead: the nominal copy's own hue at the floor of the
lightness ramp, which is a real colour, because the vetoed footprints here are ones where *some
copy* exhausted `max_steps` while the **nominal** `shape_vec` stayed finite. Zero magenta pixels
survive in the 62 files.

**What that costs**: under this style an undetermined footprint is indistinguishable from a resolved
dark one, which is exactly what the flag exists to prevent. So the count moves into the record
rather than vanishing — it is the `vetoed` column below, and `veto=quiet vetoed_footprints=N of=M`
in every sidecar. `colour::rgb` is unchanged, so every diagnostic render keeps the flag.

## The slices

Leaf counts are facts about this viewport: 256², where the screen floor binds at level 5 against the
level 6 the 512² tables in `results/payload/README.md` reach. Do not read a leaf count off a frame.
Each panel's `.cfg.txt` carries the configuration, the growth curve (`boundary:t:quads_computed`),
the stop-reason breakdown and the catch-up share.

| slice | quads | leaves | depth | frames | catch-up | stop | vetoed |
|---|---|---|---|---|---|---|---|
| `near-field` | 85 | 55 | 5 | 25 | 65.4% | floor:1 keep:42 screen_floor:12 | 0 / 3520 |
| `mid-field` | 21 | 16 | 2 | 22 | 0.0% | keep:16 | 0 / 1024 |
| `far` | 21 | 16 | 2 | 22 | 0.0% | keep:16 | 0 / 1024 |
| `body2 core` | 261 | 157 | 5 | 26 | 85.6% | floor:2 keep:105 screen_floor:50 | 0 / 10048 |
| `body2 mid` | 29 | 16 | 2 | 22 | 0.0% | keep:16 | 0 / 1024 |
| `body1 slice` | 261 | 163 | 5 | 26 | 87.6% | keep:107 screen_floor:56 | 0 / 10432 |
| `body1 far` | 113 | 76 | 5 | 24 | 53.8% | floor:5 keep:52 screen_floor:19 | 0 / 4864 |
| `deep interior` | 145 | 70 | 5 | 23 | 41.8% | floor:3 keep:32 screen_floor:35 | 0 / 4480 |
| `config_stability` | 1373 | 955 | 5 | 25 | 60.0% | floor:24 keep:85 screen_floor:846 | 49 / 61120 |
| `body_plane` | 85 | 55 | 5 | 25 | 65.4% | floor:1 keep:42 screen_floor:12 | 0 / 3520 |
| `plane_00deg` | 85 | 55 | 5 | 25 | 65.4% | floor:1 keep:42 screen_floor:12 | 0 / 3520 |
| `shape_sphere` | 909 | 622 | 5 | 27 | 90.5% | floor:16 keep:268 screen_floor:338 | 0 / 39808 |
| `latent_shape` | 209 | 121 | 5 | 24 | 61.7% | keep:103 screen_floor:18 | 0 / 7744 |
| `latent_inner_p` | 469 | 352 | 5 | 24 | 45.1% | keep:202 screen_floor:150 | 0 / 22528 |
| `latent_outer_p` | 309 | 226 | 5 | 25 | 45.1% | floor:3 keep:137 screen_floor:86 | 0 / 14464 |
| `latent_mass` | 309 | 226 | 5 | 24 | 42.9% | keep:144 screen_floor:82 | 0 / 14464 |
| `latent_mixed` | 465 | 346 | 5 | 24 | 46.1% | floor:1 keep:152 screen_floor:193 | 0 / 22144 |
| `burrau_nu_k` | 521 | 376 | 5 | 26 | 51.8% | floor:8 keep:226 screen_floor:142 | 1 / 24064 |
| `preset_shape` | 653 | 382 | 5 | 25 | 49.8% | floor:33 keep:74 screen_floor:275 | 156 / 24448 |
| `preset_prho` | 429 | 313 | 5 | 26 | 58.6% | floor:6 keep:170 screen_floor:137 | 1 / 20032 |
| `preset_plambda` | 361 | 232 | 5 | 27 | 55.6% | floor:12 keep:91 screen_floor:129 | 0 / 14848 |
| `preset_shape_pl` | 377 | 262 | 5 | 26 | 63.8% | floor:15 keep:116 screen_floor:131 | 2 / 16768 |
| `preset_shape_h1` | 833 | 439 | 5 | 27 | 51.9% | floor:59 keep:41 screen_floor:339 | 120 / 28096 |
| `preset_prho_h1` | 573 | 418 | 5 | 26 | 54.7% | floor:6 keep:219 screen_floor:193 | 4 / 26752 |
| `preset_plambda_h1` | 529 | 361 | 5 | 27 | 57.4% | floor:10 keep:152 screen_floor:199 | 12 / 23104 |
| `preset_shape_pl_h1` | 493 | 352 | 5 | 27 | 59.3% | floor:8 keep:148 screen_floor:196 | 4 / 22528 |
| `latent_shape_h3` | 505 | 268 | 5 | 24 | 45.0% | floor:2 keep:171 screen_floor:95 | 4 / 17152 |
| `latent_inner_p_h3` | 485 | 358 | 5 | 25 | 50.1% | floor:1 keep:209 screen_floor:148 | 0 / 22912 |
| `latent_outer_p_h3` | 245 | 172 | 5 | 25 | 48.8% | floor:5 keep:108 screen_floor:59 | 0 / 11008 |
| `latent_mass_h3` | 337 | 253 | 5 | 26 | 48.7% | keep:173 screen_floor:80 | 4 / 16192 |
| `latent_mixed_h3` | 453 | 331 | 5 | 26 | 74.7% | floor:7 keep:152 screen_floor:164 undetermined:8 | 2238 / 21184 |

<!-- 31 slices -->

Three things the table says. **`mid-field`, `far` and `body2 mid` stop at the bootstrap** — 16 leaves,
`keep` on all of them, nothing to refine, and `dup 19` of 22 frames identical. **Catch-up runs
0–91%**: a child requested late pays nearly a full march to reach the playhead, and where it reads
0.0% the tree never split after the bootstrap. **`latent_mixed_h3` is the outlier** — 2238 vetoed
footprints of 21184 (10.6%) and the only slice whose stop breakdown carries `undetermined`, so it is
the one chart in the set where the step budget is genuinely biting.

## Reproduce

`refine_flagged` is off because the repair pass re-integrates from `t = 0` and has no live-playhead
analogue; `keep_live_series` is on and the stride is 2. Both are named in every sidecar.

```
cargo run --release --example live_animation -- results "$(cat <<'CHARTS'
near-field,mid-field,far,body2 core,body2 mid,body1 slice,body1 far,deep interior,config_stability,body_plane,plane_00deg,shape_sphere,latent_shape,latent_inner_p,latent_outer_p,latent_mass,latent_mixed,burrau_nu_k,preset_shape,preset_prho,preset_plambda,preset_shape_pl,preset_shape_h1,preset_prho_h1,preset_plambda_h1,preset_shape_pl_h1,latent_shape_h3,latent_inner_p_h3,latent_outer_p_h3,latent_mass_h3,latent_mixed_h3
CHARTS
)" 256 4000
```

`examples/live_magenta.rs` is the harness that attributes what a vetoed footprint is: per boundary,
the run-wide count against the count known *at* that boundary, and a census of the flagged
footprints against the healthy ones.
