# `config_stability_stop0_*.png` does not reproduce at its own commit — it is STALE

## The command

From `results/output/closure_render_1024.txt`'s own header:

```sh
cargo run --release --example closure_render 1024 results
```

`only` and `sub` default, so the targets are all four and the output lands in `results/closure`.
To regenerate beside the committed set rather than over it, pass them:

```sh
cargo run --release --example closure_render 1024 <root> config_stability <sub>
```

At `4b26466` the harness has **no `sub` argument** — `tau` is argument 4 — so there the
invocation is `1024 <root> config_stability` and it writes to `<root>/closure`.

## Step 3 FAILS, and step 4 says why

| run | escape | collision | frozen | t_end distinct | on bdry | ramp | vs committed PNG |
|---|---|---|---|---|---|---|---|
| **committed** | 0.0583 | 0.4257 | 0.4832 | 446650 | 57.42% | 5.810e-5 .. 4.911e-1 | — |
| re-run at **`220d928`** | 0.0425 | 0.4304 | 0.4724 | 451601 | 56.94% | 5.764e-5 .. 4.907e-1 | **84.12% of pixels differ** |
| re-run at **`4b26466`** | 0.0583 | 0.4257 | 0.4832 | 446650 | 57.42% | 5.810e-5 .. 4.911e-1 | **BITWISE IDENTICAL**, both panels |

**The build is deterministic.** Every number matches to the digit at `4b26466`, and the PNG
matches byte for byte — `_uniform` and `_outcome` both. What fails is not reproducibility; it is
that **the artefact was committed two commits after the tree that produced it**.

## The timeline, from the reflog

```text
  08-27 13:52:53   commit 4b26466                     <- HEAD from here
  08-27 16:38      config_stability_stop0_*.png rendered   (file mtime)
  08-27 21:19:50   commit 5cc8dec   the dtau fix      <- the tree changes under the artefact
  08-27 21:36:47   commit 220d928   the PNGs added    <- committed here
```

No stash, no uncommitted state, no reset in that window. The image is a **pre-`dtau`-fix render
committed into a post-fix tree**, and the only thing between the two is `5cc8dec`.

`closure_render.rs` itself changed in `220d928`, by 13 lines — `sub` becomes argument 4 and `tau`
moves to 5, plus `preset_shape_pl_h1` joins the targets. **None of that is in the render path**,
and the committed invocation passes neither argument, so the harness change is not the cause and
was checked before the tree was.

## What it means for anything compared against this image

The magenta fraction is **0.0029** committed against **0.0001** at `220d928` — a factor of 29 —
and the escape fraction moves 0.0583 -> 0.0425. Any "before/after" that used this file as the
*before* was comparing a pre-`dtau` render against a post-`dtau` one, which is the difference
`RESULTS.md` §23 already measures at a different `n_sync`. The artefact is not evidence of a
later regression.

## Not fixed here

Two options, and the choice is not mine to make: re-render at HEAD and replace the file, or keep
it and rename it to say which commit it is from. The second preserves the only committed record
of the pre-fix state at this discretisation, which nothing else in `results/` holds.

## The render path between `220d928` and HEAD, for completeness

Two commits touch it, and three files change:

```text
  f7d2a31  the overshoot clamp        driver.rs, pixel.rs, real.rs
  c7bdece  switch instrumentation     driver.rs, pixel.rs   (default off, asserted bitwise inert)
```

`outcome.rs`, `grid.rs`, `jitter.rs`, `decoder.rs`, `energy.rs`, `shape.rs`, `colour.rs` and
`png.rs` are the **same blob** at both ends. `f7d2a31`'s mechanism is three lines:

```rust
let dtau = if clamp_final { dtau.min((dt_left - s.t).max(T::zero()) / ab) } else { dtau };
...
t += if landed && opts.clamp_final_step { dt_left } else { s.t.min(dt_left) };
```

plus `Real::LAND_EPS_REL` (`1e-14` at f64, `1e-5` at f32) making the landing tolerance relative
to `dt_left` rather than absolute.
