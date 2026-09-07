# `results/uv` — the missing middle space, and what §12's defect actually costs

Two things, and the second is a negative result that took a control to establish.

Measured 2026-09-05. Reproduce: `cargo run --release --example sample_space -- results 8 60`.
Committed stdout is `output/sample_space.txt`. `src/uv.rs` and `tests/uv.rs` carry the addressing.

## The exactness claim was true, in a space this repo did not have

The deep-zoom spec's *"`h` is an exact power of two at every depth"* was objected to as false here.
The objection was in the wrong coordinate system. The spec nests **three** spaces — screen, UV quad
addressing, IC chart — and prin-rs had only the first and third: `Slice::axis` is a `linspace` in
**chart** units and the tree's root half is `0.05` or `3.0`.

In UV, where the root half is exactly `1/2`, `uv_half(d) = 2^-(d+1)` is a power of two in the strict
sense — integer exponent, unit mantissa. In chart space `0.05` halves **exactly** forever and never
becomes one. `tests/uv.rs` asserts both, side by side, because the pair is the finding rather than
either half.

**The frame is fixed, and that is what makes an address absolute.** Not the tree root — `grow_root`
doubles the root box, and an address relative to it would renumber the whole store on a zoom-out,
which is the cost re-rooting exists to avoid.

**And a grown root leaves the lattice in three directions of four.** A cell has exactly one parent,
so only growth toward it stays addressable. That is the seam between this **rooted** tree and the
caching contract's **flat** store, as a measurement rather than a paragraph: `id_of` returns `None`
for the other three, which is the same refusal the half-cell guard makes.

## `SampleSpace::QuadLocal` buys nothing, and a different global sum is why

§12's defect is real: `jitter` recovered `du` as `(u - cx)/half` from a `u` the `linspace` had just
built globally, so precision was spent on the offset and then the offset subtracted back off.
`SampleSpace::QuadLocal` forms `du` directly from `Slice::local_pos` and never builds the sum except
where a decode needs one. `Global` stays the default.

**At production settings it is a last-ulp reordering and nothing else.** `near-field` moves **0 of
64** — bitwise identical, so the plan's acceptance condition holds there. `config_stability` moves 8
of 64 at `6.1e-16` of a **cell width**, which is a reordering of a sum and not a different sample.

**And at depth all four f64 cells fail at exactly the same rung**, depth 48:

```
 depth        half  g/f64  l/f64  g/lin  l/lin  g/f32  g/L32  l/L32        |ju|
    44   2.842e-15     64     64     64     64      1     64     64    2.887e-15
    48   1.776e-16      4      4      4      4      1      4      4    2.220e-16
    52   1.110e-17      1      1      1      1      1      1      1      0.000e0
```

The `|ju|` column is the mechanism. `decode::linearise` differences the chart at `cu ± half` — **the
same global sum, one level up** — so the Jacobian's own secant meets the chart-coordinate floor
first: it tracks `half` exactly to depth 44, quantises to the ulp of `cx` at 48, and is **exactly
zero** at 52. A linearisation whose secant has collapsed carries no information for a carried `du`
to preserve, so the four f64 cells *must* fail together, and they do.

**The fix §12 actually needs is a Jacobian not built by a secant at the cell width.** Not built;
named, with the measurement that says so.

### The control that stops the f32 row being misread

`l/L32` holds 64 through depth 44 against `g/f32`'s collapse at 16 — 28 levels, and it would read as
this change buying them. `g/L32` is **identical at every rung**, so those levels are `LinSplitF32`'s
own: the standing *"the linearised decoder buys ~24 levels over f32 and none over f64"*, measured
under `Global` long before this flag existed. Without that column the table would have credited
`QuadLocal` with someone else's result.

### And the origin ladder is the control for the whole thing

At `cx = cy = 0` **nothing collapses at any rung in any cell** and `|ju|` tracks `half` to `4.3e-20`,
because there is no O(1) neighbour for the increment to be absorbed into. The standing *"the
deep-zoom floor is a property of where you zoom"*, reproduced here on a second construction. Quote
the coordinate magnitude with any floor depth.

## Why the flag stays

`Global` is every committed number on this project and moving it moves them all, for a measured
benefit of zero. `QuadLocal` is kept and named because it is the **correct** construction and
because the thing blocking it — `linearise`'s secant — is separable and fixable; when that is fixed,
the flag is what makes the comparison one line.
