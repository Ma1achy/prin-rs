# prin-rs

A Rust kernel that renders a uniform-resolution slice of three-body initial-condition space.
Each pixel is one three-body simulation: integrate to a playhead `t`, classify the outcome, colour
it, write a PNG and a raw dump.

**The physics is the product. The image is a diagnostic.**

## The spec

**[`BRIEF.md`](BRIEF.md)** is authoritative. Read it in full before writing code — it assumes no
prior context and carries everything: the system, the slice, the regularised integrator, the
ensemble, every per-pixel field and why it exists, the failure signatures, and the definition of
done.

In particular §2.3 (integration), §4 (outputs), §5 (verification). [`CLAUDE.md`](CLAUDE.md) carries
the working agreement — scope, non-negotiables, build order.

**[`RESULTS.md`](RESULTS.md) is what the kernel measured**, written to be read cold: the
refinement criterion's resolution limit and its cost curve, which conclusions survive at large
`n`, the high-drift pixels and their remedy, and what all of it invalidates in the prior NumPy
work. [`NOTES.md`](NOTES.md) carries the mechanisms and the standing rules behind those findings.

## Layout

| path | contents |
|---|---|
| `BRIEF.md` | the build brief — authoritative spec |
| `CLAUDE.md` | working agreement: scope, invariants, review protocol |
| `RESULTS.md` | findings — what the kernel measured, and what it corrects |
| `NOTES.md` | mechanisms, corrections, and standing rules |
| `TODO.md` | the gated build checklist, kept current |
| `reference/` | validated NumPy implementation. Port `tb_az.py`; do not re-derive the algebra |
| `src/`, `Cargo.toml` | the Rust kernel |
| `examples/` | one per measurement in `RESULTS.md`; each prints its own raw table |
| `results/` | committed images, 64×64 raw dumps, and captured output for every example |
| `tools/xcheck/` | the cross-check harness against the NumPy reference |

## Running the kernel

```bash
cargo run --release --bin prin -- --region near-field --size 256 --out out
```

Writes `out.raw` (the product), `out_outcome.png` and `out_spread.png` (diagnostics). `--help`
lists the flags; `--precision f32`, `--shared-reference`, `--r-coll` and `--no-refine` are the
ones that change what is measured rather than how much.

```bash
cargo test --release                      # the invariant and acceptance suite
cargo test --release -- --ignored         # plus the NumPy cross-check
```

The cross-check needs Python and NumPy on `PATH`; nothing else does.

## Running the reference

Pure NumPy, no other dependencies. Reproduce the smoke test before porting anything, so you know
the reference runs on your machine:

```bash
cd reference
python3 -c "
import numpy as np, tb, tb_all_az as AA, warnings; warnings.filterwarnings('ignore')
r0,v0,gid,_,_ = tb.burrau_grid(3,3, 1.0,3.0, 0.05, ens=3, jitter_frac=0.5, seed=0)
res = AA.integrate_all_az(r0,v0, t_max=13.0, n_sync=32, eta=0.01)
print('median |dE/E| =', np.median(res['drift']))   # expect ~3.9e-09
"
```

See [`reference/README.md`](reference/README.md) for the module map.

## Building the kernel

Not yet landed. When it is, per BRIEF.md §7 it will be:

```bash
cargo build --release          # f64, CPU, rayon over pixels
cargo test                     # acceptance tests, BRIEF.md §5
cargo run --release -- <config>
```

Config is a small TOML or CLI: slice region, `W`, `H`, `t_max`, `E`, `eta`, `r_coll`, `epsilon`,
precision, shared-reference flag. Dependencies are `rayon`, a PNG writer, and `ndarray` or plain
`Vec` — nothing else.

## Scope

Uniform grid, one pass, write files, exit. No quadtree, no scheduler, no GUI, no streaming, no
interaction. Each omission is deliberate. **If the program grows a scheduler, that is a bug, not
progress.**
