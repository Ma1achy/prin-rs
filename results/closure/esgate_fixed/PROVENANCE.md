# `config_stability_stop0_uniform.png` — when it was made, and reproduced byte for byte

This is the image the whole bleaching thread is about. It is **256 x 256**, not 1024^2.

## Reproduced, bitwise, both panels

```sh
git checkout 4b26466
cargo run --release --example closure_render 256 <root> config_stability
# writes <root>/closure/config_stability_stop{0,1}_{uniform,outcome}.png
```

| file | committed | re-render at `4b26466` |
|---|---|---|
| `config_stability_stop0_uniform.png` | 169469 bytes | **byte-for-byte identical** |
| `config_stability_stop0_outcome.png` | 63946 bytes | **byte-for-byte identical** |

At `4b26466`, `closure_render` has no `sub` argument — `tau` is argument 4 — so the output lands
in `<root>/closure`. The `esgate_fixed` subdirectory name comes from the *later* harness; the
pixels do not depend on it.

**So there is no uncommitted working-tree difference to hunt for.** The code that produced this
image is exactly `4b26466`, and the build is deterministic across a fresh checkout and rebuild.

## When

```text
  08-27 13:52:53   commit 4b26466                       <- HEAD from here to 21:19
  08-27 14:45      esgate_fixed/preset_plambda_*        (file mtimes)
  08-27 15:20      esgate_fixed/config_stability_stop0  <- THIS FILE
  08-27 16:38      results/closure/ at 1024^2
  08-27 21:19:50   commit 5cc8dec   the dtau fix        <- the tree moves under it
  08-27 21:36:47   commit 220d928   the PNGs added      <- committed here
```

**It is a pre-`dtau`-fix render, committed into a post-fix tree.** Re-run at `220d928`, the commit
that adds it, the same command gives a *different* image — `config_stability_stop0_uniform_AT_220d928.png`
is committed beside it so the pair can be looked at rather than described.

```text
  at 4b26466 (the committed file)   escape 0.0588  collision 0.4251  frozen 0.4829  t_end distinct 27994
  at 220d928 (its own commit)       escape 0.0433  collision 0.4295  frozen 0.4722  t_end distinct 28280
```

## What is different in the code, then against now

Three commits touch the render path between `4b26466` and HEAD, and only two can move a pixel:

```text
  5cc8dec  the dtau step control       driver.rs, pixel.rs, az/mod.rs, colour.rs
  f7d2a31  the boundary-overshoot clamp driver.rs, pixel.rs, real.rs
  c7bdece  switch instrumentation       driver.rs, pixel.rs   (default off, asserted bitwise inert)
```

Everything else in the path is the **same blob** at both ends: `outcome.rs`, `grid.rs`,
`jitter.rs`, `stats.rs`, `reference_body.rs`, `decoder.rs`, `energy.rs`, `shape.rs`, `newton.rs`,
`png.rs`, `adaptive.rs`, `decode.rs`. `colour.rs` changes only additively — `Scalar::Drift`, the
inferno table, and `range` delegating to a new `range_q(px, s, 0.01, 0.99)` — so the uniform
panel's colouring is unchanged.

## The resolution, which is the first thing to read off it

**256 x 256.** One footprint is 4 x 4 screen pixels when this is viewed at the size the 1024^2
renders are, so a magenta cluster of a few footprints reads as a *blob*, and speckle reads as
texture. That is this project's standing result at a fourth site — *softness in an image is a
raster size, not a rendering fault; measure the dimensions first.* The magenta fraction is
**0.0030** here against **0.0029** in the 1024^2 render of the same configuration, so the physics
agrees and only the raster differs.
