# prin-rs

A Rust kernel that renders a uniform-resolution slice of three-body initial-condition space —
one simulation per pixel, integrated to a playhead, classified, coloured, written out.

The physics is the product. The image is a diagnostic.

## Start here

**[`principia_brief_kernel_build.md`](principia_brief_kernel_build.md)** is the entry point and the
specification. It assumes no prior context: the system, the slice, the regularised integrator, the
ensemble, every per-pixel field and why it exists, the acceptance tests, and the definition of done.

Read it before writing code. In particular §2.3 (integration) and §5 (verification).

## Layout

| path | contents |
|---|---|
| `principia_brief_kernel_build.md` | the build brief — spec, physics, acceptance tests |
| `reference/` | validated NumPy implementation. Port `tb_az.py`; do not re-derive the algebra |
| `src/`, `Cargo.toml` | the Rust kernel (not yet landed) |

`reference/README.md` carries the module map and the smoke test whose number a port must match.

## Scope

Uniform grid, one pass, write files, exit. No quadtree, no scheduler, no GUI, no streaming, no
interaction. Each omission is deliberate — see the brief. If it grows one, that is a bug.
